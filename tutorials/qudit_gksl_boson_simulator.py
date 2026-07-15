"""Qudit GKSL simulator with boson (phonon) interaction using Stinespring + Trotter.

Scenario 6: Qudit-based GKSL-Lindblad with boson.

Uses native qutrit (d=3) encoding for the electronic degrees of freedom and a
``(n_max+1)``-level qudit per phonon mode.  No forbidden states for either
subsystem.  Works in the extended electronic+phonon Hilbert space
``dim_total = d^N * (n_max+1)^N``, using Stinespring + 2nd-order Trotter
(``algorithm="stinespring"``, default, **back-compatible**) or
**exact local channels** (``algorithm="exact_local_channels"``).

DMSim backend execution
-----------------------
With ``algorithm="exact_local_channels"`` and ``execute_on_backend="dmsim"``
this simulator builds, for each Trotter step, a mixed-dimensional MQT-Qudits
``QuantumCircuit`` over ``[d]*N + [n_max+1]*N`` qudits and dispatches the
Hamiltonian half-steps (``cu_multi``) and Lindblad channels
(``KrausChannel`` on the electronic qudits only — the dissipator is purely
electronic and acts trivially on phonon) to
:class:`mqt.qudits.simulation.backends.DMSim`.  Each Kraus channel is
extracted exactly via the Choi-Jamiolkowski isomorphism from the
electronic local superoperator ``expm(L_D^local · dt/2)`` — no
heuristic, no truncation other than discarding strictly non-positive
Choi eigenvalues whose magnitude is indistinguishable from numerical
noise (see :mod:`dmsim_kraus_helpers`).
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import (
    build_H_total_boson,
    build_lindblad_operators,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    extend_lindblad_operators,
    partial_trace_phonon,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)
from exact_local_channels import (
    apply_channel_pair,
    apply_channel_single,
    precompute_exact_channels_half,
)
from dmsim_kraus_helpers import kraus_from_local_superoperator


class QuditGKSLBosonSimulator:
    """Qudit GKSL simulator with boson (phonon) interaction.

    Uses native qutrit (d=3) encoding for electronic states and qutrit encoding
    for phonon modes (when n_max=2). No forbidden states for either subsystem.
    Works in the extended space dim_total = d^N * (n_max+1)^N.

    Since this targets qudit quantum computers, all registers — including
    Stinespring ancillas — are native d-level qudits.
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        algorithm: str = "stinespring",
        execute_on_backend: str | None = None,
    ) -> None:
        if not params.with_boson:
            raise ValueError("QuditGKSLBosonSimulator requires with_boson=True")
        if algorithm not in ("stinespring", "exact_local_channels"):
            msg = (
                "algorithm must be 'stinespring' (1st order, default, "
                "back-compatible) or 'exact_local_channels' (2nd order); "
                f"got {algorithm!r}"
            )
            raise ValueError(msg)
        if execute_on_backend is not None:
            if execute_on_backend != "dmsim":
                msg = (
                    "execute_on_backend currently only supports 'dmsim' "
                    "for the boson model.  state-vector backends "
                    "(tnsim/misim) cannot run the fresh-ancilla per-step "
                    "circuit at these sizes (see STATUS_HONEST_2026-05.md "
                    f"A-1); got {execute_on_backend!r}."
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
        self.N = params.N_molecules
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension

        # Dimensions
        self.dim_el = params.d ** params.N_molecules
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_total = self.dim_el * self.dim_ph

        # Qudit counts
        self.n_system_qudits = params.N_molecules  # electronic qutrits
        # Number of d-level qudits needed per phonon mode: ceil(log_d(n_max+1))
        n_ph_per_mol = int(np.ceil(np.log(params.n_max + 1) / np.log(params.d)))
        self.n_phonon_qudits = n_ph_per_mol * params.N_molecules
        self.n_ancilla_qudits = 2 * len(params.neighbors) + 5 * params.N_molecules
        self.n_total_qudits = self.n_system_qudits + self.n_phonon_qudits + self.n_ancilla_qudits

        # Build extended Hamiltonian
        self.H_total = build_H_total_boson(params)

        # Build extended Lindblad operators (used by Stinespring path).
        lindblad_ops_el = build_lindblad_operators(params)
        self.lindblad_ops = extend_lindblad_operators(lindblad_ops_el, self.dim_ph)

        # ----------------------------------------------------------------
        # exact_local_channels / DMSim backend bookkeeping
        # ----------------------------------------------------------------
        # The full electronic Lindblad operators are local (single-site or
        # adjacent pair); on the joint el+phonon space they act trivially
        # on the phonon subsystem (``L_ext = L_el ⊗ I_ph``).  Therefore the
        # exact-local-channels machinery built for the non-boson case
        # applies *unchanged* to the electronic qudits — we just need a
        # twin parameter set with ``with_boson=False`` to drive
        # :func:`exact_local_channels.precompute_exact_channels_half`.
        self._params_el = GKSLPhysicalParameters(
            E_T=params.E_T,
            E_S=params.E_S,
            V=params.V,
            gamma_TTA=params.gamma_TTA,
            Gamma_fl=params.Gamma_fl,
            Gamma_ph=params.Gamma_ph,
            k_IC=params.k_IC,
            k_ISC_ST=params.k_ISC_ST,
            k_ISC_TS=params.k_ISC_TS,
            N_molecules=params.N_molecules,
            d=params.d,
            with_boson=False,
        )

        # DMSim backend — initialise eagerly when requested so that
        # provider/backend lookup failures surface at construction time
        # rather than at the first ``simulate`` call.
        self._dmsim_backend = None
        self._dmsim_kraus_half: list[
            tuple[tuple[int, ...], list[np.ndarray], str]
        ] | None = None
        if self.execute_on_backend == "dmsim":
            from mqt.qudits.simulation import MQTQuditProvider
            self._dmsim_backend = MQTQuditProvider().get_backend("dmsim")
        # Phonon dimensions per mode (for KrausChannel mixed-dim register).
        self._phonon_dim = params.n_max + 1
        # Number of phonon qudits per electronic site = 1 here (we keep the
        # phonon Fock space as one qudit per mode, dimension n_max+1).
        # Total qudits in the DMSim register: N electronic + N phonon.
        self._dmsim_dims = [params.d] * params.N_molecules + [self._phonon_dim] * params.N_molecules

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation).

        For ``algorithm="stinespring"``: half-dt Stinespring unitaries on
        the extended Lindblad operators.

        For ``algorithm="exact_local_channels"``: exponentiated electronic
        local dissipators ``expm(L_D^local · dt/2)`` (3×3 single, 9×9 pair),
        identical to the non-boson case — see the class docstring for why
        this is mathematically the same as the extended local channel.

        For ``execute_on_backend="dmsim"``: additionally extract the
        electronic Kraus operators of every cached half-step local
        superoperator (Choi-Jamiolkowski decomposition,
        :func:`dmsim_kraus_helpers.kraus_from_local_superoperator`) so
        that each Lindblad channel can be issued as a MQT-Qudits
        :class:`KrausChannel` instruction acting on the electronic
        qudits only.
        """
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        if self.algorithm == "stinespring":
            self._U_stines_half = [
                stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
                for L_op, _gamma in self.lindblad_ops
            ]
            self._exact_channels_half = None
        else:
            # exact_local_channels — electronic-only, on the d^N space.
            self._U_stines_half = None
            self._exact_channels_half = precompute_exact_channels_half(
                self._params_el, dt
            )

        if self.execute_on_backend == "dmsim":
            kraus_list: list[tuple[tuple[int, ...], list[np.ndarray], str]] = []
            for sites, M_half, kind in self._exact_channels_half:
                d_root = self.params.d if kind == "single" else self.params.d ** 2
                kraus = kraus_from_local_superoperator(M_half, d_root)
                kraus_list.append((sites, kraus, kind))
            self._dmsim_kraus_half = kraus_list

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad channel ordering.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        With ``algorithm="stinespring"`` the channels ``E_α(dt/2)`` are
        Stinespring dilations (1st-order → global O(dt)); with
        ``algorithm="exact_local_channels"`` they are exact local
        superoperator exponentials applied to the electronic subsystem
        (2nd-order → global O(dt²)).  With ``execute_on_backend="dmsim"``
        the entire step is dispatched to DMSim as one mixed-dimensional
        ``QuantumCircuit``.
        """
        if self.execute_on_backend == "dmsim":
            return self._trotter_step_dmsim(rho)

        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        if self.algorithm == "stinespring":
            for U_stine_half in self._U_stines_half:
                rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
            for U_stine_half in reversed(self._U_stines_half):
                rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
        else:
            # exact_local_channels: apply each electronic local channel
            # to the joint el+phonon density matrix.  Two palindromic
            # passes (forward + reverse) over the cached half-step
            # channels; each channel acts trivially on the phonon
            # subsystem.
            rho = self._apply_exact_channels_pass(rho, reverse=False)
            rho = self._apply_exact_channels_pass(rho, reverse=True)
        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        return rho

    def _apply_exact_channels_pass(
        self, rho: np.ndarray, *, reverse: bool
    ) -> np.ndarray:
        """One palindromic pass of exact electronic local channels on the joint state.

        The joint Hilbert space has dimension ``dim_el * dim_ph`` with the
        Kronecker ordering ``el ⊗ ph`` (matching
        :func:`gksl_math_utils.build_H_total_boson`).  Each electronic
        Lindblad channel acts non-trivially on one (or two) electronic
        qutrit and trivially on phonon, so we can apply
        :func:`exact_local_channels.apply_channel_single`/`pair` slice by
        slice over phonon-block indices ``(b_r, b_c)``.
        """
        d = self.params.d
        N = self.N
        dim_el = self.dim_el
        dim_ph = self.dim_ph

        # Reshape: ρ[(el_r, ph_r), (el_c, ph_c)] -> (el_r, ph_r, el_c, ph_c)
        rho_t = rho.reshape(dim_el, dim_ph, dim_el, dim_ph)
        # Transpose phonon legs to the back: (el_r, el_c, ph_r, ph_c)
        rho_t = rho_t.transpose(0, 2, 1, 3)
        rho_block = rho_t.reshape(dim_el, dim_el, dim_ph * dim_ph).copy()

        channels = (
            list(reversed(self._exact_channels_half))
            if reverse
            else list(self._exact_channels_half)
        )

        out_block = np.empty_like(rho_block)
        for b in range(dim_ph * dim_ph):
            rho_el = rho_block[:, :, b]
            for sites, M_half, kind in channels:
                if kind == "single":
                    rho_el = apply_channel_single(rho_el, M_half, sites[0], N, d)
                else:
                    rho_el = apply_channel_pair(rho_el, M_half, sites, N, d)
            out_block[:, :, b] = rho_el

        rho_t = out_block.reshape(dim_el, dim_el, dim_ph, dim_ph)
        rho_t = rho_t.transpose(0, 2, 1, 3)
        return rho_t.reshape(dim_el * dim_ph, dim_el * dim_ph)

    # ------------------------------------------------------------------
    # DMSim backend Trotter step (one MQT-Qudits backend run per step)
    # ------------------------------------------------------------------

    def _trotter_step_dmsim(self, rho: np.ndarray) -> np.ndarray:
        """Run one Trotter step as a single MQT-Qudits DMSim backend run.

        Builds a mixed-dimensional ``QuantumCircuit`` over
        ``[d]*N + [n_max+1]*N`` qudits (electronic qutrits first, phonon
        qudits second — matching the el⊗ph Kronecker ordering used by
        :func:`gksl_math_utils.build_H_total_boson`) containing:

        1. A ``cu_multi`` on **all** qudits implementing
           ``expm(-i H_total dt/2)`` — this mixes electronic and phonon.
        2. Each electronic Lindblad half-step channel as a
           :class:`~mqt.qudits.quantum_circuit.gates.KrausChannel`
           instruction in the **forward** order from
           :func:`exact_local_channels.precompute_exact_channels_half`,
           targeting only the electronic qudits (``sites`` indices).
        3. Same channels again in the **reverse** order — palindromic
           structure preserves 2nd-order convergence.
        4. A second ``cu_multi`` for the closing Hamiltonian half-step.

        The previous step's ρ is supplied as ``initial_density_matrix``
        and the returned density matrix is the input for the next step.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit
        from mqt.qudits.quantum_circuit.components.quantum_register import (
            QuantumRegister,
        )

        N = self.N
        d = self.params.d
        d_ph = self._phonon_dim
        n_total = 2 * N

        qreg = QuantumRegister("sys_ph", n_total, list(self._dmsim_dims))
        circuit = QuantumCircuit(qreg)

        all_qudits = list(range(n_total))

        # 1. Hamiltonian half-step on all qudits.
        circuit.cu_multi(all_qudits, self._U_H_half.astype(np.complex128))

        # 2. Forward palindromic pass over electronic Kraus channels.
        for sites, kraus_ops, kind in self._dmsim_kraus_half:
            if kind == "single":
                circuit.kraus_channel(sites[0], kraus_ops)
            else:
                circuit.kraus_channel(list(sites), kraus_ops)

        # 3. Reverse palindromic pass.
        for sites, kraus_ops, kind in reversed(self._dmsim_kraus_half):
            if kind == "single":
                circuit.kraus_channel(sites[0], kraus_ops)
            else:
                circuit.kraus_channel(list(sites), kraus_ops)

        # 4. Closing Hamiltonian half-step.
        circuit.cu_multi(all_qudits, self._U_H_half.astype(np.complex128))

        job = self._dmsim_backend.run(circuit, initial_density_matrix=rho)
        return job.result().get_density_matrix()

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic state ⊗ phonon vacuum |00...0⟩."""
        d = self.params.d
        N = self.N
        dim_el = d ** N

        # Electronic state
        psi_el = np.zeros(dim_el, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi_el[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d ** i) for i in range(N))
            psi_el[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d ** i) for i in range(N))
            psi_el[index] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

        # Phonon vacuum: |00...0⟩
        psi_ph = np.zeros(self.dim_ph, dtype=np.complex128)
        psi_ph[0] = 1.0

        # Full state: |psi_el⟩ ⊗ |0_ph⟩
        psi_total = np.kron(psi_el, psi_ph)
        return np.outer(psi_total, psi_total.conj())

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run Qudit GKSL+boson simulation.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qudit circuit statistics. Populations are computed by
        partial-tracing over phonon degrees of freedom.
        """
        start = time_module.time()

        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        # Partial trace for initial populations
        rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_el, self.params)]
        entropies = [compute_von_neumann_entropy(rho_el)]
        purities = [compute_purity(rho_el)]
        traces = [float(np.real(np.trace(rho_el)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_el))))
            populations.append(
                compute_populations_from_density_matrix(rho_el, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        elapsed = time_module.time() - start

        # Return electronic-only reduced density matrix for consistency
        # with ClassicalGKSLBosonSimulator and to match time-series observables
        rho_el_final = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

        # Gate count estimate for qudit circuit in extended space
        # (palindromic 2nd-order Trotter)
        # Qudit advantage: native d-level ops, no forbidden states
        # 2 × N cu_one gates (electronic on-site, two half-steps)
        # 2 × N cu_one gates (phonon on-site, two half-steps)
        # 2 × (N-1) cu_two gates (H_transfer, two half-steps)
        # 2 × N cu_two gates (H_eph electron-phonon coupling, two half-steps)
        # 2 × n_lindblad Stinespring gates (fwd + rev)
        gates_per_step = (
            2 * (self.n_system_qudits + self.n_phonon_qudits
                 + len(self.params.neighbors) + self.params.N_molecules)
            + self.n_ancilla_qudits * 2
        )

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_el_final,
            "elapsed_time": elapsed,
            "method": "qudit_gksl_boson",
            "algorithm": self.algorithm,
            "execute_on_backend": self.execute_on_backend,
            "params": self.params.to_dict(),
            "n_system_qudits": self.n_system_qudits,
            "n_phonon_qudits": self.n_phonon_qudits,
            "n_ancilla_qudits": self.n_ancilla_qudits,
            "d_anc": self.d_anc,
            "n_total_qudits": self.n_total_qudits,
            "dim_total": self.dim_total,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
