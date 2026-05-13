"""Backend-executed "shot-based" GKSL simulators for the comparison notebook.

Honest design statement (no heuristics, no lies)
================================================

The original :mod:`qudit_gksl_shot_simulator` and
:mod:`qubit_gksl_shot_simulator` modules implement the quantum-trajectory
unraveling of the GKSL dynamics using **NumPy + scipy.linalg.expm** only.
They do not invoke any quantum backend (this was verified by exhaustive
search for ``backend|QuantumCircuit|Provider|execute|Aer|MQTQuditProvider|
tnsim|misim|dmsim`` in those files: zero hits).

This module provides backend-executed *replacements* for the shot cells
of ``tutorials/quantum_dynamics_gksl_comparison.ipynb`` (cells 5b/5c/3b/3c).

What this module DOES
---------------------
* ρ(t) for the full Trotter trajectory is produced by a real backend:

  - Qudit path: MQT-Qudits :class:`~mqt.qudits.simulation.backends.DMSim`
    (``execute_on_backend="dmsim"``, ``algorithm="exact_local_channels"``).
    Each Trotter step is one MQT-Qudits :class:`QuantumCircuit` containing
    ``cu_multi`` (Hamiltonian half-step) plus
    :class:`~mqt.qudits.quantum_circuit.gates.KrausChannel` instructions
    (Lindblad channels, palindromic forward+reverse).
  - Qubit path: Qiskit Aer ``AerSimulator(method="density_matrix")``.
    Each Trotter step is one :class:`qiskit.QuantumCircuit` containing
    ``UnitaryGate`` (Hamiltonian half-step) plus
    :class:`qiskit.quantum_info.Kraus` instructions (Lindblad channels,
    palindromic forward+reverse).

* Hardware noise is injected as **exact CPTP Kraus channels** at the
  same gate locations as :class:`QuditGKSLNoisySimulator` /
  :class:`QubitGKSLNoisySimulator`:

    after every Hamiltonian half-step,
        for each nearest-neighbour pair (i, j),
            apply pair-depolarisation Kraus channel.

  The depolarisation Kraus operators are constructed analytically from
  the d-dimensional Weyl–Heisenberg twirling identity

      (1/d²) Σ_{a,b} W_{a,b} ρ W_{a,b}†  =  (Tr ρ) I/d                (1)

  giving the d² Kraus operators
      K₀          = sqrt(1 − p + p/d²) · I
      K_{(a,b)≠0} = sqrt(p/d²)        · W_{a,b}                       (2)
  for single-site depolarisation; the pair channel uses the d⁴ Kraus
  operators ``W_{a,b} ⊗ W_{c,d}`` with the same prefactors but on
  ``d²``-dim (i.e. d → d² in (2)).  No fudge factor, no truncation:
  ``Σ K†K = I`` is verified to machine precision and the construction
  raises if it fails.

* "Shot counts" are produced by **exact Born sampling** of the backend's
  final density matrix:

      counts ~ Multinomial(n_shots, diag(ρ_final))

  via :func:`numpy.random.multinomial`, seeded for reproducibility.

What this module DOES NOT do (and would be dishonest to claim)
--------------------------------------------------------------
* This is **not** a per-trajectory quantum-circuit simulation.  The
  backend is run **once** per Trotter step on the density-matrix
  representation, not ``n_shots`` independent state-vector trajectories
  with mid-circuit ancilla measurement.  The latter is currently
  impossible on the MQT-Qudits backends (the per-step "fresh ancilla +
  measure + reset" pattern that the trajectory unraveling needs has no
  state-vector backend that can run it at the notebook's N=4 due to
  state-vector capacity limits — see ``STATUS_HONEST_2026-05.md`` A-1)
  and is impractical on Qiskit Aer's state-vector method at N=4.

* The intermediate-time populations are the **exact** density-matrix
  values, not averages over n_shots realisations.  Shot-to-shot
  fluctuations of intermediate-time observables are therefore not
  produced.  Only the final-time ``counts`` carry genuine multinomial
  Born-sampling variance, and only for diagonal (computational-basis)
  observables.  For non-diagonal observables, sampling diag(ρ_final)
  is **not** equivalent to true trajectory variance.

* The hardware-noise model implemented here is the
  ``depol_pair_only=True, p_dephasing=0.0`` regime (the only regime the
  notebook actually uses for cells 5c/3c).  Other combinations are
  rejected at construction time rather than silently producing wrong
  noise; extending to single-site noise / dephasing requires inserting
  additional Kraus channels in the per-step circuit.

References
----------
* Notebook cells affected: 5b / 5c / 3b / 3c of
  ``tutorials/quantum_dynamics_gksl_comparison.ipynb``.
* Underlying density-matrix backend simulators reused without
  modification:

  - :class:`tutorials.qudit_gksl_simulator.QuditGKSLSimulator`
  - :class:`tutorials.qiskit_qubit_gksl_simulator.QiskitQubitGKSLSimulator`

* Trajectory references kept for cross-checks:

  - :class:`tutorials.qudit_gksl_shot_simulator.QuditGKSLShotSimulator`
    (NumPy ideal trajectory).
  - :class:`tutorials.qudit_gksl_shot_simulator.QuditGKSLNoisyShotSimulator`
    (NumPy stochastic-Pauli trajectory).
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_simulator import (  # noqa: E402  (sys.path bootstrap above)
    compute_forbidden_state_population,
)


# ---------------------------------------------------------------------------
# Weyl–Heisenberg / depolarisation Kraus operator construction (exact)
# ---------------------------------------------------------------------------

_KRAUS_CPTP_TOL = 1e-10
"""Tolerance for ``||Σ K†K − I||_F`` of constructed Kraus channels."""


def _weyl_heisenberg(d: int) -> list[np.ndarray]:
    """Return the d² Weyl–Heisenberg unitaries ``W_{a,b} = X^a Z^b``.

    ``X|k⟩ = |(k+1) mod d⟩`` and ``Z|k⟩ = ω^k |k⟩`` with ``ω = e^{2πi/d}``.
    These satisfy the twirling identity (1) in the module docstring.
    """
    if d < 2:
        raise ValueError(f"d must be >= 2, got {d}")
    omega = np.exp(2j * np.pi / d)
    # Z is diagonal: Z[k,k] = ω^k
    Z = np.diag([omega**k for k in range(d)]).astype(np.complex128)
    # X is the cyclic shift: X[i,j] = δ_{i,(j+1) mod d}
    X = np.zeros((d, d), dtype=np.complex128)
    for j in range(d):
        X[(j + 1) % d, j] = 1.0
    Xa = [np.linalg.matrix_power(X, a) for a in range(d)]
    Zb = [np.linalg.matrix_power(Z, b) for b in range(d)]
    return [Xa[a] @ Zb[b] for a in range(d) for b in range(d)]


def depolarisation_kraus(d: int, p: float) -> list[np.ndarray]:
    """Build exact Kraus operators of the d-dim depolarisation channel.

    ``E_p[ρ] = (1−p) ρ + p · I/d`` factored via the Weyl–Heisenberg
    twirling identity, giving exactly d² Kraus operators of the form
    given by equation (2) in the module docstring.

    Parameters
    ----------
    d:
        Local Hilbert space dimension.
    p:
        Depolarisation probability in ``[0, 1]``.

    Returns
    -------
    list of d² ``(d, d)`` complex matrices satisfying ``Σ K†K = I`` to
    machine precision.  Raises if closure deviates by more than
    ``_KRAUS_CPTP_TOL`` (no heuristic correction is applied).
    """
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"p must be in [0,1], got {p}")
    if p == 0.0:
        return [np.eye(d, dtype=np.complex128)]
    Ws = _weyl_heisenberg(d)
    coef0 = np.sqrt(1.0 - p + p / (d * d))
    coef_other = np.sqrt(p / (d * d))
    # Index 0 in Ws is W_{0,0} = I.
    kraus = [coef0 * Ws[0]] + [coef_other * W for W in Ws[1:]]
    _verify_cptp(kraus, d, channel_name=f"depolarisation(d={d}, p={p})")
    return kraus


def pair_depolarisation_kraus(d: int, p: float) -> list[np.ndarray]:
    """Pair-depolarisation Kraus operators on the ``(d²)``-dim local space.

    ``E_p[ρ] = (1−p) ρ + p · I/d²`` on the d²-dim joint space, factored
    as ``W_{a,b} ⊗ W_{c,d}`` over the d⁴ pairs.  This is exactly the
    same channel that
    :func:`qudit_gksl_noisy_simulator._apply_local_depolarization_pair`
    applies via partial trace; the Weyl–Heisenberg twirling identity
    (1) on each factor gives the operator-sum form used here.
    """
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"p must be in [0,1], got {p}")
    d2 = d * d
    if p == 0.0:
        return [np.eye(d2, dtype=np.complex128)]
    Ws = _weyl_heisenberg(d)
    coef0 = np.sqrt(1.0 - p + p / (d2 * d2))
    coef_other = np.sqrt(p / (d2 * d2))
    kraus: list[np.ndarray] = []
    for a, Wa in enumerate(Ws):
        for b, Wb in enumerate(Ws):
            K = np.kron(Wa, Wb)
            if a == 0 and b == 0:
                kraus.append(coef0 * K)
            else:
                kraus.append(coef_other * K)
    _verify_cptp(kraus, d2, channel_name=f"pair_depolarisation(d={d}, p={p})")
    return kraus


def _verify_cptp(kraus: list[np.ndarray], d_local: int, *, channel_name: str) -> None:
    """Assert ``||Σ K†K − I||_F < tol`` (no heuristic adjustment)."""
    closure = np.zeros((d_local, d_local), dtype=np.complex128)
    for K in kraus:
        closure += K.conj().T @ K
    defect = float(np.linalg.norm(closure - np.eye(d_local), ord="fro"))
    if defect > _KRAUS_CPTP_TOL:
        msg = (
            f"{channel_name}: ||Σ K†K - I||_F = {defect:.3e} > "
            f"{_KRAUS_CPTP_TOL:g}.  No heuristic correction is applied."
        )
        raise AssertionError(msg)


# ---------------------------------------------------------------------------
# Born-sampling helper (exact final-state multinomial)
# ---------------------------------------------------------------------------


def sample_counts_from_density_matrix(
    rho: np.ndarray, n_shots: int, seed: int | None = None
) -> dict[int, int]:
    """Draw ``n_shots`` computational-basis outcomes from ``diag(ρ)``.

    Parameters
    ----------
    rho:
        Density matrix in computational basis (dim × dim, hermitian,
        unit trace up to numerical noise).
    n_shots:
        Number of Born-rule samples to draw.
    seed:
        Seed for :class:`numpy.random.Generator`.

    Returns
    -------
    Mapping ``{basis_index: count}`` with ``Σ count == n_shots``.
    Zero-count outcomes are omitted.

    Notes
    -----
    This is the standard "measure ρ in the computational basis" of an
    exactly-known final density matrix — the same statistics a real
    quantum hardware run would produce **for diagonal observables**, in
    the limit where the device's per-shot state preparation and
    evolution match ρ exactly.  No heuristic re-normalisation: the
    diagonal is rescaled only to compensate for floating-point trace
    noise (and the original trace is reported by the caller).
    """
    diag = np.real(np.diag(rho))
    # Clip non-physical tiny negatives that come from finite-precision
    # arithmetic on a CPTP map; abort if the deficit is non-trivial
    # rather than silently masking it.
    min_val = float(diag.min())
    if min_val < -1e-10:
        msg = (
            f"diag(ρ) has a negative entry {min_val:.3e} below the "
            f"finite-precision floor 1e-10.  Refusing to renormalise."
        )
        raise ValueError(msg)
    diag = np.clip(diag, 0.0, None)
    total = float(diag.sum())
    if total <= 0.0:
        raise ValueError("diag(ρ) sums to zero; cannot sample.")
    probs = diag / total
    rng = np.random.default_rng(seed)
    outcomes = rng.multinomial(n_shots, probs)
    return {i: int(c) for i, c in enumerate(outcomes) if c > 0}


# ---------------------------------------------------------------------------
# Population accounting (qudit; matches qudit_gksl_shot_simulator output)
# ---------------------------------------------------------------------------


def _populations_from_diagonal(
    diag: np.ndarray, params: GKSLPhysicalParameters
) -> dict:
    N = params.N_molecules
    d = params.d
    dim = d**N

    N_S0 = 0.0
    N_T1 = 0.0
    N_S1 = 0.0
    per_molecule: list[dict[int, float]] = [{0: 0.0, 1: 0.0, 2: 0.0} for _ in range(N)]

    for idx in range(dim):
        p = float(diag[idx])
        remainder = idx
        for mol in range(N - 1, -1, -1):
            local_state = remainder % d
            remainder //= d
            if local_state == 0:
                N_S0 += p
            elif local_state == 1:
                N_T1 += p
            else:
                N_S1 += p
            per_molecule[mol][local_state] += p

    return {
        "N_S0": float(N_S0),
        "N_T1": float(N_T1),
        "N_S1": float(N_S1),
        "per_molecule_populations": per_molecule,
    }


# ---------------------------------------------------------------------------
# Qudit DMSim shot simulator (replaces 5b / 5c)
# ---------------------------------------------------------------------------


class QuditDMSimShotSimulator:
    """Backend-executed "shot-based" qudit GKSL simulator.

    Runs the full Trotter trajectory on the MQT-Qudits ``dmsim``
    backend (density-matrix method), then performs exact Born sampling
    of ``diag(ρ_final)`` to produce ``n_shots`` computational-basis
    outcomes.  Optionally injects pair-depolarisation Kraus channels
    after each Hamiltonian half-step (matching the gate-noise placement
    of :class:`QuditGKSLNoisySimulator` for ``depol_pair_only=True,
    p_dephasing=0.0``).

    Parameters
    ----------
    params:
        Physical parameters (must be ``with_boson=False``).
    p_depol:
        Depolarisation probability per pair gate.  ``0.0`` for the
        ideal (5b) case; ``0.001`` for the notebook's noisy (5c) case.
    p_dephasing:
        Currently must be ``0.0``.  Non-zero dephasing is not yet
        wired through the dmsim circuit — supporting it requires
        additional Kraus channels per single qudit; raise rather than
        silently ignore.
    depol_pair_only:
        Currently must be ``True``.  ``False`` would require additional
        single-site noise channels per Lindblad operator.

    Notes
    -----
    See the module docstring for the full honest scope statement
    (intermediate-time populations are exact, only final counts carry
    Born-sampling variance, etc.).
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.0,
        p_dephasing: float = 0.0,
        depol_pair_only: bool = True,
    ) -> None:
        if params.with_boson:
            raise ValueError(
                "QuditDMSimShotSimulator targets the non-boson model only"
            )
        if p_dephasing != 0.0:
            raise NotImplementedError(
                "p_dephasing != 0 is not implemented in the dmsim shot "
                "simulator (would require additional single-qudit Kraus "
                "channels). Refusing to silently drop it."
            )
        if not depol_pair_only:
            raise NotImplementedError(
                "depol_pair_only=False is not implemented in the dmsim "
                "shot simulator (would require single-site Kraus per "
                "Lindblad). Refusing to silently drop it."
            )
        if not 0.0 <= p_depol <= 1.0:
            raise ValueError(f"p_depol must be in [0,1], got {p_depol}")

        # Lazy import to avoid pulling MQT-Qudits at module import time
        # (mirrors the pattern used elsewhere in the tutorials package).
        from qudit_gksl_simulator import QuditGKSLSimulator

        self.params = params
        self.p_depol = float(p_depol)
        self.p_dephasing = 0.0
        self.depol_pair_only = True

        # Underlying density-matrix backend simulator (DMSim).
        # We always use exact_local_channels for backend execution —
        # Stinespring has no Kraus form on the system alone.
        self._dm_sim = QuditGKSLSimulator(
            params,
            algorithm="exact_local_channels",
            execute_on_backend="dmsim",
        )

        # Pre-build pair-depolarisation Kraus operators (independent of dt).
        if self.p_depol > 0.0:
            self._pair_depol_kraus = pair_depolarisation_kraus(
                params.d, self.p_depol
            )
        else:
            self._pair_depol_kraus = None

    # ------------------------------------------------------------------
    # Override the dmsim Trotter step to inject pair-depol Kraus.
    # ------------------------------------------------------------------

    def _trotter_step_dmsim_with_noise(self, rho: np.ndarray) -> np.ndarray:
        """Build one MQT-Qudits Trotter circuit including noise and run on DMSim.

        Mirrors :meth:`QuditGKSLSimulator._trotter_step_dmsim` but inserts
        pair-depolarisation Kraus channels after each Hamiltonian
        half-step at every nearest-neighbour pair.  This matches the
        noise placement of :class:`QuditGKSLNoisyShotSimulator` for the
        ``depol_pair_only=True, p_dephasing=0.0`` regime.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit
        from mqt.qudits.quantum_circuit.components.quantum_register import (
            QuantumRegister,
        )

        n = self._dm_sim.n_system_qudits
        d = self.params.d

        qreg = QuantumRegister("sys", n, [d] * n)
        circuit = QuantumCircuit(qreg)

        # 1. Hamiltonian half-step.
        circuit.cu_multi(
            list(range(n)),
            self._dm_sim._U_H_half.astype(np.complex128),
        )

        # 1b. Pair-depolarisation noise on every neighbour pair (if enabled).
        if self._pair_depol_kraus is not None:
            for i, j in self.params.neighbors:
                circuit.kraus_channel([i, j], self._pair_depol_kraus)

        # 2. Forward palindromic pass (Lindblad channels as Kraus).
        for sites, kraus_ops, kind in self._dm_sim._dmsim_kraus_half:
            target = sites[0] if kind == "single" else list(sites)
            circuit.kraus_channel(target, kraus_ops)

        # 3. Reverse palindromic pass (Lindblad channels as Kraus).
        for sites, kraus_ops, kind in reversed(self._dm_sim._dmsim_kraus_half):
            target = sites[0] if kind == "single" else list(sites)
            circuit.kraus_channel(target, kraus_ops)

        # 4. Closing Hamiltonian half-step.
        circuit.cu_multi(
            list(range(n)),
            self._dm_sim._U_H_half.astype(np.complex128),
        )

        # 4b. Pair-depolarisation noise on every neighbour pair (closing).
        if self._pair_depol_kraus is not None:
            for i, j in self.params.neighbors:
                circuit.kraus_channel([i, j], self._pair_depol_kraus)

        job = self._dm_sim._dmsim_backend.run(circuit, initial_density_matrix=rho)
        return job.result().get_density_matrix()

    # ------------------------------------------------------------------
    # Main simulate() — runs ρ(t) on dmsim then samples shots from ρ_final
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        n_shots: int = 1000,
        seed: int | None = None,
    ) -> dict:
        """Run dmsim density-matrix evolution then Born-sample n_shots."""
        from gksl_math_utils import (
            compute_populations_from_density_matrix,
            compute_purity,
            compute_von_neumann_entropy,
        )
        from mqt.qudits.simulation.qudit_provider import MQTQuditProvider

        start = time_module.time()
        dt = t_max / n_steps

        # Pre-compute Hamiltonian half-step and Lindblad Kraus.
        self._dm_sim._precompute_unitaries(dt)
        # Acquire the dmsim backend handle (matches QuditGKSLSimulator.simulate).
        self._dm_sim._dmsim_backend = MQTQuditProvider().get_backend("dmsim")

        rho = self._dm_sim.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho, self.params)]
        entropies = [compute_von_neumann_entropy(rho)]
        purities = [compute_purity(rho)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step_dmsim_with_noise(rho)
            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        # Born sampling at the end.
        counts = sample_counts_from_density_matrix(rho, n_shots, seed=seed)

        elapsed = time_module.time() - start

        # High-level gate count (palindromic 2nd-order Trotter).
        n_lindblad = len(self._dm_sim.lindblad_ops)
        n_pairs = len(self.params.neighbors)
        n_pair_kraus = (2 * n_pairs) if self.p_depol > 0.0 else 0
        gates_per_step = (
            2 * (self._dm_sim.n_system_qudits + n_pairs)
            + n_lindblad * 2
            + n_pair_kraus
        )

        result: dict = {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": (
                "qudit_gksl_dmsim_shot_noisy"
                if self.p_depol > 0.0
                else "qudit_gksl_dmsim_shot"
            ),
            "backend": "mqt_qudits:dmsim",
            "execution": "density_matrix_trajectory + final_state_born_sampling",
            "params": self.params.to_dict(),
            "n_system_qudits": self._dm_sim.n_system_qudits,
            "n_ancilla_qudits": self._dm_sim.n_ancilla_qudits,
            "n_total_qudits": self._dm_sim.n_total_qudits,
            "d_anc": self._dm_sim.d_anc,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "n_shots": n_shots,
            "counts": counts,
            "seed": seed,
        }
        if self.p_depol > 0.0:
            result["noise_params"] = {
                "p_depol": self.p_depol,
                "p_dephasing": self.p_dephasing,
                "depol_pair_only": self.depol_pair_only,
                "noise_model": "pair_depolarisation_kraus",
            }
        return result


# ---------------------------------------------------------------------------
# Qiskit-Aer qubit shot simulator (replaces 3b / 3c)
# ---------------------------------------------------------------------------


class QiskitQubitShotSimulator:
    """Backend-executed "shot-based" qubit GKSL simulator on Qiskit Aer.

    Runs the full Trotter trajectory on Qiskit Aer's
    ``density_matrix`` method, then performs exact Born sampling of
    ``diag(ρ_final)`` to produce ``n_shots`` computational-basis
    outcomes in qubit-pair encoding.  Optionally injects 4-qubit
    pair-depolarisation Kraus channels after each Hamiltonian
    half-step (one per nearest-neighbour molecule pair), matching the
    placement of :class:`QubitGKSLNoisySimulator` for
    ``depol_pair_only=True, p_dephasing=0.0``.

    The forbidden state ``|11⟩`` (per molecule) is reported via
    :func:`compute_forbidden_state_population` and forbidden outcomes
    are counted from the final-state Born sampling.

    Parameters
    ----------
    params:
        Physical parameters (must be ``with_boson=False``).
    p_depol:
        Depolarisation probability per pair gate.  ``0.0`` for cell 3b;
        ``0.001`` for cell 3c.
    p_dephasing, depol_pair_only:
        Same restriction as :class:`QuditDMSimShotSimulator` (must be
        ``0.0`` and ``True`` respectively); enforced explicitly to
        prevent silent semantic drift.
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.0,
        p_dephasing: float = 0.0,
        depol_pair_only: bool = True,
    ) -> None:
        if params.with_boson:
            raise ValueError(
                "QiskitQubitShotSimulator targets the non-boson model only"
            )
        if p_dephasing != 0.0:
            raise NotImplementedError(
                "p_dephasing != 0 is not implemented in the Qiskit Aer "
                "qubit shot simulator. Refusing to silently drop it."
            )
        if not depol_pair_only:
            raise NotImplementedError(
                "depol_pair_only=False is not implemented in the Qiskit "
                "Aer qubit shot simulator."
            )
        if not 0.0 <= p_depol <= 1.0:
            raise ValueError(f"p_depol must be in [0,1], got {p_depol}")

        from qiskit_qubit_gksl_simulator import QiskitQubitGKSLSimulator

        self.params = params
        self.p_depol = float(p_depol)
        self.p_dephasing = 0.0
        self.depol_pair_only = True

        # Underlying density-matrix backend simulator.
        self._dm_sim = QiskitQubitGKSLSimulator(params)

        # Pre-build 4-qubit (16-dim) pair-depolarisation Kraus operators.
        # Each "molecule pair" spans 4 qubits = 16 dims; we reuse the
        # generic pair_depolarisation_kraus() with d=4 so that the
        # local Hilbert space matches the qubit-pair-encoded molecule pair.
        if self.p_depol > 0.0:
            self._pair_depol_kraus = pair_depolarisation_kraus(
                d=4, p=self.p_depol
            )
        else:
            self._pair_depol_kraus = None

    # ------------------------------------------------------------------
    # Override the per-step circuit to inject pair-depol Kraus.
    # ------------------------------------------------------------------

    def _build_trotter_circuit_with_noise(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import Kraus

        sim = self._dm_sim
        qc = QuantumCircuit(sim.n_qubits)

        all_qubits = list(range(sim.n_qubits))
        qc.append(UnitaryGate(sim._U_H_half), all_qubits)

        # Pair-depolarisation noise on every neighbour molecule pair.
        if self._pair_depol_kraus is not None:
            kraus_inst = Kraus(self._pair_depol_kraus).to_instruction()
            for i, j in self.params.neighbors:
                qbits = list(sim._mol_qubits[i]) + list(sim._mol_qubits[j])
                qc.append(kraus_inst, qbits)

        # Forward palindromic Lindblad pass.
        for sites, kraus_ops, kind in sim._kraus_per_channel:
            sim._append_kraus(qc, Kraus, sites, kraus_ops, kind)
        # Reverse palindromic Lindblad pass.
        for sites, kraus_ops, kind in reversed(sim._kraus_per_channel):
            sim._append_kraus(qc, Kraus, sites, kraus_ops, kind)

        qc.append(UnitaryGate(sim._U_H_half), all_qubits)
        # Closing pair-depolarisation pass.
        if self._pair_depol_kraus is not None:
            kraus_inst = Kraus(self._pair_depol_kraus).to_instruction()
            for i, j in self.params.neighbors:
                qbits = list(sim._mol_qubits[i]) + list(sim._mol_qubits[j])
                qc.append(kraus_inst, qbits)

        return qc

    # ------------------------------------------------------------------
    # Main simulate() — runs ρ(t) on Aer then samples shots from ρ_final
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        n_shots: int = 1000,
        seed: int | None = None,
    ) -> dict:
        from gksl_math_utils import (
            compute_populations_from_density_matrix,
            compute_purity,
            compute_von_neumann_entropy,
        )
        from qiskit import QuantumCircuit, transpile
        from qiskit.quantum_info import DensityMatrix
        from qubit_gksl_simulator import extract_density_matrix_from_qubit_space

        start = time_module.time()
        dt = t_max / n_steps

        sim = self._dm_sim
        sim._precompute(dt)

        rho = sim.prepare_initial_state(initial_state)
        rho_qt = extract_density_matrix_from_qubit_space(
            rho, sim._mapping, sim.dim_qutrit
        )

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_qt, self.params)]
        entropies = [compute_von_neumann_entropy(rho_qt)]
        purities = [compute_purity(rho_qt)]
        traces = [float(np.real(np.trace(rho_qt)))]
        forbidden_pops: list[float] = [
            compute_forbidden_state_population(rho, sim._mapping)
        ]

        step_circ = self._build_trotter_circuit_with_noise()

        for step in range(n_steps):
            wrapped = QuantumCircuit(sim.n_qubits)
            wrapped.set_density_matrix(DensityMatrix(rho))
            wrapped.compose(step_circ, inplace=True)
            wrapped.save_density_matrix()
            result = sim._aer.run(transpile(wrapped, sim._aer)).result()
            rho = np.asarray(result.data(0)["density_matrix"], dtype=np.complex128)

            rho_qt = extract_density_matrix_from_qubit_space(
                rho, sim._mapping, sim.dim_qutrit
            )
            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_qt))))
            populations.append(
                compute_populations_from_density_matrix(rho_qt, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_qt))
            purities.append(compute_purity(rho_qt))
            forbidden_pops.append(
                compute_forbidden_state_population(rho, sim._mapping)
            )

        # Born sampling on the qubit-space density matrix.  We sample in
        # qubit space (so the "forbidden" indices are reachable and
        # countable), report counts per qubit-basis index, and also count
        # how many shots landed on a forbidden molecular state.
        counts_qb = sample_counts_from_density_matrix(rho, n_shots, seed=seed)
        forbidden_count = _count_forbidden_shots(counts_qb, sim._mapping, sim.dim_qubit)

        elapsed = time_module.time() - start

        result_dict: dict = {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_qt,
            "rho_final_qubit": rho,
            "elapsed_time": elapsed,
            "method": (
                "qiskit_qubit_gksl_shot_noisy"
                if self.p_depol > 0.0
                else "qiskit_qubit_gksl_shot"
            ),
            "backend": "qiskit_aer:density_matrix",
            "execution": "density_matrix_trajectory + final_state_born_sampling",
            "params": self.params.to_dict(),
            "n_sys_qubits": sim.n_qubits,
            "n_ancilla": 0,
            "n_total_qubits": sim.n_qubits,
            "dim_qubit_space": sim.dim_qubit,
            "n_shots": n_shots,
            "counts": counts_qb,
            "forbidden_count": forbidden_count,
            "forbidden_state_population": forbidden_pops,
            "seed": seed,
            "estimated_gates_per_step": _qubit_gates_per_step(sim, p_depol=self.p_depol),
            "total_estimated_gates": _qubit_gates_per_step(sim, p_depol=self.p_depol) * n_steps,
        }
        if self.p_depol > 0.0:
            result_dict["noise_params"] = {
                "p_depol": self.p_depol,
                "p_dephasing": self.p_dephasing,
                "depol_pair_only": self.depol_pair_only,
                "noise_model": "pair_depolarisation_kraus_4qubit",
            }
        return result_dict


def _qubit_gates_per_step(sim, *, p_depol: float) -> int:
    """High-level instruction count per Trotter step in the qubit circuit.

    * 2 ``UnitaryGate`` (Hamiltonian half-steps)
    * 2 × |kraus_per_channel| Kraus instructions (palindromic Lindblad)
    * + 2 × |neighbours| Kraus instructions (pair depolarisation) when
      noise is enabled.
    """
    n_lindblad = len(sim._kraus_per_channel)
    n_pairs = len(sim.params.neighbors)
    n_noise = (2 * n_pairs) if p_depol > 0.0 else 0
    return 2 + 2 * n_lindblad + n_noise


def _count_forbidden_shots(
    counts: dict[int, int], mapping: dict[int, int], dim_qubit: int
) -> int:
    """Count Born samples whose qubit-basis index lies outside the
    embedded physical (molecular) subspace.

    ``mapping`` maps every physical qutrit-product basis index to the
    corresponding qubit-product basis index.  Any qubit basis index
    *not* in ``mapping.values()`` corresponds to at least one molecule
    in the forbidden ``|11⟩`` state.
    """
    physical_qubit_indices = set(mapping.values())
    return int(sum(c for i, c in counts.items() if i not in physical_qubit_indices))
