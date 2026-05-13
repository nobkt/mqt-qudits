"""Qudit GKSL simulator using Stinespring dilation + Trotter decomposition.

Scenario 5: Qudit-based GKSL-Lindblad (no boson).

Uses native qutrit (d=3) encoding. Each molecule is naturally 1 qutrit,
so there are NO forbidden states (unlike the qubit encoding which wastes |11>).
The simulation uses the same 81-dim qutrit Hilbert space and Stinespring+Trotter
approach as the qubit simulator, but corresponds to a qudit circuit with fewer gates.

Independence note (A-2)
-----------------------
The qubit GKSL simulator (:mod:`qubit_gksl_simulator`) is **not** an
algorithmically independent re-implementation of this class.  It shares the
same Hamiltonian generators, the same Stinespring dilation routine and the
same Trotter step structure; the qubit simulator only adds a qutrit→qubit
embedding of every operator.  Cross-checks between the two simulators
therefore validate the embedding, not the underlying open-system dynamics.

Convergence note
----------------
The Hamiltonian–Dissipator splitting uses Strang (symmetric) splitting, which
is 2nd-order for the H-D decomposition.  However, the individual Lindblad
channels are applied via Stinespring dilation, which is a 1st-order
approximation of each exact Lindblad channel exp(L_{D_α} dt).  Additionally,
the sequential product of individual channel maps introduces a Lie-Trotter
product error.  The Lindblad channels are applied in *symmetric (palindromic)
order* to eliminate the Lie-Trotter commutator error, but the Stinespring
approximation itself remains 1st-order.  As a result, the effective
convergence in trace distance is **O(dt)** (1st-order), not O(dt²).
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
)
from exact_local_channels import (
    apply_channel_pair,
    apply_channel_single,
    precompute_exact_channels_half,
)
from dmsim_kraus_helpers import kraus_from_local_superoperator
from gksl_physical_parameters import GKSLPhysicalParameters
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)


class QuditGKSLSimulator:
    """Qudit GKSL simulator using Stinespring dilation + Trotter splitting.

    Uses native qutrit (d=3) encoding. No forbidden states exist.
    4 qutrits for system + 26 ancilla qudits (d=3) for Lindblad channels.

    Since this simulator targets qudit-type quantum computers (e.g.
    qudit-boson ion-trap processors), all registers — including Stinespring
    ancillas — are native d-level qudits, not 2-level qubits.

    Two channel-application algorithms are supported (parameter
    ``algorithm`` in :meth:`__init__`):

    * ``"stinespring"`` (default, **back-compatible**): each Lindblad
      channel is realised by a Stinespring dilation
      ``U_α = expm(-i √dt · G_α)`` and a partial trace over the ancilla.
      This is mathematically a **1st-order** approximation of
      ``exp(L_{D_α} dt)``; combined with the Strang H–D split and
      palindromic channel ordering the *global* convergence in trace
      distance is **O(dt)** — see B-1 in
      ``STATUS_HONEST_2026-05.md``.
    * ``"exact_local_channels"``: each Lindblad channel is exponentiated
      **exactly** as a local superoperator
      ``expm(L_{D_α}^local · dt)`` (size ``9×9`` for single-site, ``81×81``
      for TTA-pair channels) and applied to the system density matrix
      via :mod:`exact_local_channels`.  No Stinespring approximation is
      used.  Empirical convergence in trace distance is **O(dt²)** —
      see :class:`tutorials.test_gksl_simulators.TestExactLocalChannelsConvergence`.
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        algorithm: str = "stinespring",
        execute_on_backend: str | None = None,
    ) -> None:
        if params.with_boson:
            raise ValueError("QuditGKSLSimulator is for non-boson model only")
        if algorithm not in ("stinespring", "exact_local_channels"):
            msg = (
                "algorithm must be 'stinespring' (1st order, default) or "
                "'exact_local_channels' (2nd order); got "
                f"{algorithm!r}"
            )
            raise ValueError(msg)
        if execute_on_backend is not None:
            if execute_on_backend != "dmsim":
                msg = (
                    "execute_on_backend currently only supports 'dmsim' "
                    "(MQT-Qudits density-matrix backend); got "
                    f"{execute_on_backend!r}.  state-vector backends "
                    "(tnsim/misim) cannot run the fresh-ancilla per-step "
                    "circuit at N>=3 due to state-vector capacity limits "
                    "(see STATUS_HONEST_2026-05.md A-1)."
                )
                raise ValueError(msg)
            if algorithm != "exact_local_channels":
                msg = (
                    "execute_on_backend='dmsim' requires "
                    "algorithm='exact_local_channels' so that each Lindblad "
                    "channel has an exact Kraus representation that can be "
                    "fed to MQT-Qudits as a KrausChannel instruction.  "
                    "Stinespring dilation has no Kraus representation on "
                    "the system alone (it requires ancillas)."
                )
                raise ValueError(msg)
        self.algorithm = algorithm
        self.execute_on_backend = execute_on_backend
        self.params = params
        self.n_system_qudits = params.N_molecules
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension

        # Build operators in native qutrit space
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.n_ancilla_qudits = len(self.lindblad_ops)
        self.n_total_qudits = self.n_system_qudits + self.n_ancilla_qudits

        # System dimension (qutrit space)
        self.dim = params.d ** params.N_molecules  # 81

        # Lazy initialisation for DMSim backend execution.
        # _dmsim_kraus_half holds, for each Lindblad channel, a tuple of
        #   (target_sites, [K_α], "single" | "pair")
        # in the same order as exact_local_channels.precompute_exact_channels_half.
        self._dmsim_backend = None
        self._dmsim_kraus_half: list[tuple[tuple[int, ...], list[np.ndarray], str]] | None = None

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation).

        For ``algorithm="stinespring"`` we cache the half-dt Stinespring
        unitaries used in the symmetric palindromic product.

        For ``algorithm="exact_local_channels"`` we instead cache, for
        every Lindblad channel, the exponentiated local dissipator
        ``expm(L_D^local · dt/2)`` and the site indices on which it
        acts.  The forward + reverse pass over the cached half-step
        channels then realises ``exp(L_D · dt)`` exactly up to the
        ``O(dt³)`` Lie-product commutator error among local channels —
        which is the whole point of the palindromic ordering.

        For ``execute_on_backend="dmsim"`` we additionally extract the
        Kraus operators of every cached half-step local superoperator
        (Choi-Jamiolkowski decomposition,
        :func:`dmsim_kraus_helpers.kraus_from_local_superoperator`) so
        that each Lindblad channel can be issued as a MQT-Qudits
        :class:`KrausChannel` instruction.
        """
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        if self.algorithm == "stinespring":
            self._U_stines_half = [
                stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
                for L_op, _gamma in self.lindblad_ops
            ]
            self._exact_channels_half = None
        else:
            # exact_local_channels
            self._U_stines_half = None
            self._exact_channels_half = precompute_exact_channels_half(
                self.params, dt
            )

        if self.execute_on_backend == "dmsim":
            # Extract Kraus operators for every half-step local channel.
            kraus_list: list[tuple[tuple[int, ...], list[np.ndarray], str]] = []
            for sites, M_half, kind in self._exact_channels_half:
                d_root = self.params.d if kind == "single" else self.params.d ** 2
                kraus = kraus_from_local_superoperator(M_half, d_root)
                kraus_list.append((sites, kraus, kind))
            self._dmsim_kraus_half = kraus_list

            # Lazily acquire the DMSim backend.
            from mqt.qudits.simulation import MQTQuditProvider
            self._dmsim_backend = MQTQuditProvider().get_backend("dmsim")

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad channel ordering.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        where ``E_α(dt/2)`` is either a Stinespring dilation
        (``algorithm="stinespring"``, **1st-order** approximation of
        ``exp(L_{D_α} dt/2)`` → global O(dt)) or the exact local channel
        ``expm(L_{D_α}^local · dt/2)`` applied via
        :mod:`exact_local_channels` (``algorithm="exact_local_channels"``,
        exact → global O(dt²)).
        """
        if self.execute_on_backend == "dmsim":
            return self._trotter_step_dmsim(rho)

        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        if self.algorithm == "stinespring":
            for U_stine_half in self._U_stines_half:
                rho = apply_stinespring_to_density_matrix(
                    rho, U_stine_half, d_anc=self.d_anc
                )
            for U_stine_half in reversed(self._U_stines_half):
                rho = apply_stinespring_to_density_matrix(
                    rho, U_stine_half, d_anc=self.d_anc
                )
        else:
            N = self.n_system_qudits
            d = self.params.d
            for sites, exp_LD_half, kind in self._exact_channels_half:
                if kind == "single":
                    rho = apply_channel_single(rho, exp_LD_half, sites[0], N, d)
                else:
                    rho = apply_channel_pair(rho, exp_LD_half, sites, N, d)
            for sites, exp_LD_half, kind in reversed(self._exact_channels_half):
                if kind == "single":
                    rho = apply_channel_single(rho, exp_LD_half, sites[0], N, d)
                else:
                    rho = apply_channel_pair(rho, exp_LD_half, sites, N, d)
        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        return rho

    # ------------------------------------------------------------------
    # DMSim backend Trotter step (one MQT-Qudits backend run per step)
    # ------------------------------------------------------------------

    def _trotter_step_dmsim(self, rho: np.ndarray) -> np.ndarray:
        """Run one Trotter step as a single MQT-Qudits DMSim backend run.

        Builds an N-qutrit :class:`mqt.qudits.quantum_circuit.QuantumCircuit`
        containing:

        1. A ``cu_multi`` unitary on all system qudits implementing
           ``expm(-i H_total dt/2)``.
        2. Each Lindblad half-step channel as a
           :class:`~mqt.qudits.quantum_circuit.gates.KrausChannel` instruction
           in the **forward** order from
           :func:`exact_local_channels.precompute_exact_channels_half`,
           then again in the **reverse** order — preserving the palindromic
           structure that makes the dissipator pass 2nd-order in ``dt``.
        3. A second ``cu_multi`` for the closing ``expm(-i H_total dt/2)``.

        The circuit is then executed on the ``dmsim`` backend with the
        previous step's ρ as ``initial_density_matrix``; the returned
        density matrix is the input for the next step.

        This is genuine MQT-Qudits backend execution: each step's circuit
        is a real :class:`QuantumCircuit` whose instructions are
        dispatched by :class:`~mqt.qudits.simulation.backends.DMSim` —
        no NumPy ``expm`` is called inside the simulation loop, only at
        the once-per-simulation pre-compute stage (which is a property
        of the Hamiltonian Trotter step itself, identical for any
        backend choice).
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit
        from mqt.qudits.quantum_circuit.components.quantum_register import (
            QuantumRegister,
        )

        n = self.n_system_qudits
        d = self.params.d

        qreg = QuantumRegister("sys", n, [d] * n)
        circuit = QuantumCircuit(qreg)

        # 1. Hamiltonian half-step (full N-qudit unitary).
        circuit.cu_multi(list(range(n)), self._U_H_half.astype(np.complex128))

        # 2. Forward palindromic pass.
        for sites, kraus_ops, kind in self._dmsim_kraus_half:
            target = sites[0] if kind == "single" else list(sites)
            circuit.kraus_channel(target, kraus_ops)

        # 3. Reverse palindromic pass.
        for sites, kraus_ops, kind in reversed(self._dmsim_kraus_half):
            target = sites[0] if kind == "single" else list(sites)
            circuit.kraus_channel(target, kraus_ops)

        # 4. Closing Hamiltonian half-step.
        circuit.cu_multi(list(range(n)), self._U_H_half.astype(np.complex128))

        job = self._dmsim_backend.run(circuit, initial_density_matrix=rho)
        return job.result().get_density_matrix()

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix in qutrit space."""
        d = self.params.d
        N = self.params.N_molecules
        dim = d ** N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            # Molecules 0 and N-1 in T1, rest in S0
            index = 1 * (d ** (N - 1)) + 1
            psi[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d ** i) for i in range(N))
            psi[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d ** i) for i in range(N))
            psi[index] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

        return np.outer(psi, psi.conj())

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run Qudit GKSL simulation.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qudit circuit statistics.
        """
        start = time_module.time()

        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho, self.params)]
        entropies = [compute_von_neumann_entropy(rho)]
        purities = [compute_purity(rho)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        elapsed = time_module.time() - start

        # Per-step *high-level* gate count for the qudit circuit
        # (palindromic 2nd-order Trotter).  This is the number of high-level
        # gate objects that would be appended to a MQT-Qudits QuantumCircuit
        # by ``QuditGKSLKrausSimulator.build_*_circuit``:
        #   * 2 × N cu_one gates (H_0 onsite, two half-steps)
        #   * 2 × |neighbors| cu_two gates (H_transfer pairs, two half-steps)
        #   * 2 × n_lindblad cu_two/cu_multi gates (Stinespring fwd + rev)
        # It is **not** the output of an MQT-Qudits compiler pass.  For the
        # actual native gate breakdown (VirtRz / R / Rh / Rz / CEx counts
        # plus the residual undecomposable cu_multi count) call
        # :meth:`compute_compiler_measured_gate_counts` — see A-3 in
        # ``STATUS_HONEST_2026-05.md``.
        n_lindblad = len(self.lindblad_ops)
        gates_per_step = 2 * (self.n_system_qudits + len(self.params.neighbors)) + n_lindblad * 2

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": "qudit_gksl",
            "algorithm": self.algorithm,
            "params": self.params.to_dict(),
            "n_system_qudits": self.n_system_qudits,
            "n_ancilla_qudits": self.n_ancilla_qudits,
            "n_total_qudits": self.n_total_qudits,
            "d_anc": self.d_anc,
            "estimated_gates_per_step": gates_per_step,
            "n_high_level_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "gate_count_method": "high_level_count",
        }

    # ------------------------------------------------------------------
    # A-3: compiler-measured native gate counts
    # ------------------------------------------------------------------

    def compute_compiler_measured_gate_counts(
        self,
        dt: float,
        optimization_level: int = 0,
        backend_name: str = "faketraps2trits",
    ) -> dict:
        """Return the **MQT-Qudits compiler-measured** native gate counts.

        Delegates to
        :meth:`qudit_gksl_circuit_simulator.QuditGKSLKrausSimulator.compile_to_native_gates`,
        which actually invokes ``compileO0`` / ``compileO1`` on the
        per-block MQT-Qudits ``QuantumCircuit`` objects and counts the
        resulting ``VirtRz`` / ``R`` / ``Rh`` / ``Rz`` / ``CEx`` gates.
        The TTA-pair Stinespring uses a ``cu_multi`` gate which the
        MQT-Qudits compiler does **not** currently decompose into 2-qudit
        primitives — that count is reported separately as
        ``per_step_summary['uncompiled_cu_multi']`` rather than being
        masked by a heuristic estimate.

        Parameters
        ----------
        dt:
            Time step used to construct the per-block unitaries that
            are then compiled.
        optimization_level:
            ``0`` for ``compileO0`` (baseline) or ``1`` for ``compileO1``
            (optimised).
        backend_name:
            MQT-Qudits backend identifier passed to the compiler.

        Returns
        -------
        dict
            Same structure as
            :meth:`QuditGKSLKrausSimulator.compile_to_native_gates`,
            including ``per_step_summary['native_gates']`` and
            ``per_step_summary['uncompiled_cu_multi']``.

        Notes
        -----
        Imported lazily to keep ``QuditGKSLSimulator`` usable in
        environments where the MQT-Qudits compiler stack is not
        available; in that case this method will raise
        ``ImportError``.
        """
        # Lazy import to avoid a hard dependency at simulator construction
        # time and to make the import failure mode explicit.
        from qudit_gksl_circuit_simulator import QuditGKSLKrausSimulator

        kraus_sim = QuditGKSLKrausSimulator(self.params)
        return kraus_sim.compile_to_native_gates(
            dt=dt,
            optimization_level=optimization_level,
            backend_name=backend_name,
        )
