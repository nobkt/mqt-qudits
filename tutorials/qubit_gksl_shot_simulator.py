"""Qubit GKSL shot-based simulator using quantum trajectories.

Implements shot-based (quantum trajectory) simulation for the qubit-encoded
GKSL-Lindblad dynamics. Each molecule is encoded in 2 qubits
(|00>=S0, |01>=T1, |10>=S1, |11>=forbidden). The computation runs in the
4^N = 256 dimensional qubit Hilbert space, faithfully representing what a
qubit quantum computer would compute.

Two classes are provided:
  - QubitGKSLShotSimulator: ideal shot-based (no hardware noise)
  - QubitGKSLNoisyShotSimulator: shot-based with per-gate stochastic noise
    using d=2 Pauli operators per qubit (with forbidden-state leakage)
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
from gksl_physical_parameters import GKSLPhysicalParameters
from gksl_validation import PhysicsViolationError
from qubit_gksl_simulator import (
    build_qubit_qutrit_mapping,
    compute_forbidden_state_population,
    embed_operator_in_qubit_space,
    embed_statevector_in_qubit_space,
    extract_density_matrix_from_qubit_space,
    extract_statevector_from_qubit_space,
)
from stinespring_utils import stinespring_unitary_from_lindblad

# ---------------------------------------------------------------------------
# 2-qubit Pauli matrices for qubit noise model
# ---------------------------------------------------------------------------

_PAULI_2X2 = [
    np.eye(2, dtype=np.complex128),  # I
    np.array([[0, 1], [1, 0]], dtype=np.complex128),  # X
    np.array([[0, -1j], [1j, 0]], dtype=np.complex128),  # Y
    np.array([[1, 0], [0, -1]], dtype=np.complex128),  # Z
]

# Pre-build all 16 two-qubit Pauli matrices (4x4) for per-molecule noise
_PAULI_4X4: list[np.ndarray] = []
for _pa in _PAULI_2X2:
    for _pb in _PAULI_2X2:
        _PAULI_4X4.append(np.kron(_pa, _pb))


def _apply_pauli_on_molecule(
    psi: np.ndarray, mol: int, pauli_idx: int, N: int
) -> np.ndarray:
    """Apply a 4x4 Pauli operator on molecule mol in the 4^N qubit space.

    pauli_idx: 0..15 indexing into _PAULI_4X4 (0 = identity).
    """
    d_local = 4
    P = _PAULI_4X4[pauli_idx]
    psi_tensor = psi.reshape([d_local] * N)
    result = np.zeros_like(psi_tensor)

    for i in range(d_local):
        for j in range(d_local):
            if abs(P[i, j]) > 1e-15:
                idx_src = [slice(None)] * N
                idx_src[mol] = j
                idx_dst = [slice(None)] * N
                idx_dst[mol] = i
                result[tuple(idx_dst)] += P[i, j] * psi_tensor[tuple(idx_src)]

    return result.reshape(d_local**N)


def _apply_stochastic_qubit_depolarization_single(
    psi: np.ndarray, mol: int, N: int, p: float, rng: np.random.Generator
) -> np.ndarray:
    """Stochastically apply 2-qubit Pauli depolarization on a single molecule.

    Depolarization channel on d=4 local space:
      E(rho) = (1-p) rho + (p/d) I Tr(rho)

    Stochastic (unraveled) form:
      - With prob (1 - p + p/d^2): identity (no error)
      - With prob p/d^2 each: apply Pauli P_k for k=1..15

    NOTE: Pauli errors can map physical states into the forbidden |11> state
    (leakage). This is a fundamental property of qubit encoding.
    """
    if p <= 0.0:
        return psi
    d_local = 4
    p_identity = 1.0 - p + p / (d_local * d_local)  # 1-p+p/16
    if rng.random() < p_identity:
        return psi
    # Select random non-identity Pauli (1..15)
    pauli_idx = rng.integers(1, d_local * d_local)
    return _apply_pauli_on_molecule(psi, mol, pauli_idx, N)


def _apply_stochastic_qubit_depolarization_pair(
    psi: np.ndarray, mol_a: int, mol_b: int, N: int, p: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stochastically apply 4-qubit Pauli depolarization on a molecule pair.

    Depolarization channel on d=16 joint space (4 qubits):
      E(rho) = (1-p) rho + (p/d^2) I Tr(rho)

    Stochastic form:
      - With prob (1 - p + p/d^4): identity
      - With prob p/d^4 each: apply P_a x P_b for non-identity pairs

    where d=16 for the 4-qubit space, giving 256 Pauli operators.
    """
    if p <= 0.0:
        return psi
    d_pair = 16  # 4 qubits
    p_identity = 1.0 - p + p / (d_pair * d_pair)  # 1-p+p/256
    if rng.random() < p_identity:
        return psi
    # Select random non-identity 4-qubit Pauli (1..255)
    error_idx = rng.integers(1, d_pair * d_pair)
    # Decompose: error_idx = pauli_a * 16 + pauli_b
    pauli_a = error_idx // 16  # 0..15 for molecule a
    pauli_b = error_idx % 16  # 0..15 for molecule b
    psi = _apply_pauli_on_molecule(psi, mol_a, pauli_a, N)
    psi = _apply_pauli_on_molecule(psi, mol_b, pauli_b, N)
    return psi


def _apply_stochastic_qubit_thermal_relaxation(
    psi: np.ndarray, mol: int, N: int, p_reset: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stochastically apply thermal relaxation on a molecule's 2 qubits.

    Models T1 relaxation on each physical qubit: each qubit independently
    decays to |0> with probability p_reset.

    For the 2-qubit encoding |b1 b0>:
      K_00 = diag(1, sqrt(1-p), sqrt(1-p), (1-p))  (no decay)
      K_01 = sqrt(p) |b1,0><b1,1|   (second qubit decays)
      K_10 = sqrt(p) |0,b0><1,b0|   (first qubit decays)
      K_11 = p |00><11|              (both qubits decay)

    Stochastic: sample Kraus operator by Born rule.
    """
    if p_reset <= 0.0:
        return psi

    d_local = 4
    dim = d_local**N
    psi_tensor = psi.reshape([d_local] * N)
    sq = np.sqrt(1.0 - p_reset)
    sqp = np.sqrt(p_reset)

    # Compute probabilities for each outcome
    # K_00: no decay (both qubits survive)
    k00_diag = np.array([1.0, sq, sq, sq * sq])
    p_k00 = 0.0
    for k in range(d_local):
        idx = [slice(None)] * N
        idx[mol] = k
        block = psi_tensor[tuple(idx)]
        p_k00 += k00_diag[k] ** 2 * float(np.real(np.vdot(block, block)))

    # K_01: second qubit (bit 0) decays: |x1> -> |x0>, applies sqrt(p) factor
    # Maps: |01>->|00>, |11>->|10>
    idx_01 = [slice(None)] * N
    idx_01[mol] = 1  # |01>
    block_01 = psi_tensor[tuple(idx_01)]
    idx_11 = [slice(None)] * N
    idx_11[mol] = 3  # |11>
    block_11 = psi_tensor[tuple(idx_11)]
    p_k01 = p_reset * (
        float(np.real(np.vdot(block_01, block_01)))
        + (1.0 - p_reset) * float(np.real(np.vdot(block_11, block_11)))
    )

    # K_10: first qubit (bit 1) decays: |1x> -> |0x>, applies sqrt(p) factor
    # Maps: |10>->|00>, |11>->|01>
    idx_10 = [slice(None)] * N
    idx_10[mol] = 2  # |10>
    block_10 = psi_tensor[tuple(idx_10)]
    p_k10 = p_reset * (
        float(np.real(np.vdot(block_10, block_10)))
        + (1.0 - p_reset) * float(np.real(np.vdot(block_11, block_11)))
    )

    # K_11: both qubits decay: |11> -> |00>
    p_k11 = p_reset * p_reset * float(np.real(np.vdot(block_11, block_11)))

    total = p_k00 + p_k01 + p_k10 + p_k11
    if total < 1e-15:
        return psi

    r = rng.random() * total
    if r < p_k00:
        # Apply K_00: multiply by survival amplitudes
        result_tensor = psi_tensor.copy()
        for k in range(d_local):
            idx = [slice(None)] * N
            idx[mol] = k
            result_tensor[tuple(idx)] *= k00_diag[k]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k00)
    elif r < p_k00 + p_k01:
        # K_01: second qubit decays (bit 0: 1->0)
        result_tensor = np.zeros_like(psi_tensor)
        # |01> -> |00> with sqrt(p)
        idx_dst_00 = [slice(None)] * N
        idx_dst_00[mol] = 0
        result_tensor[tuple(idx_dst_00)] += sqp * psi_tensor[tuple(idx_01)]
        # |11> -> |10> with sqrt(p)*sqrt(1-p) (first qubit survives)
        idx_dst_10 = [slice(None)] * N
        idx_dst_10[mol] = 2
        result_tensor[tuple(idx_dst_10)] += sqp * sq * psi_tensor[tuple(idx_11)]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k01)
    elif r < p_k00 + p_k01 + p_k10:
        # K_10: first qubit decays (bit 1: 1->0)
        result_tensor = np.zeros_like(psi_tensor)
        # |10> -> |00> with sqrt(p)
        idx_dst_00 = [slice(None)] * N
        idx_dst_00[mol] = 0
        result_tensor[tuple(idx_dst_00)] += sqp * psi_tensor[tuple(idx_10)]
        # |11> -> |01> with sqrt(p)*sqrt(1-p) (second qubit survives)
        idx_dst_01 = [slice(None)] * N
        idx_dst_01[mol] = 1
        result_tensor[tuple(idx_dst_01)] += sqp * sq * psi_tensor[tuple(idx_11)]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k10)
    else:
        # K_11: both qubits decay
        result_tensor = np.zeros_like(psi_tensor)
        idx_dst_00 = [slice(None)] * N
        idx_dst_00[mol] = 0
        result_tensor[tuple(idx_dst_00)] += p_reset * psi_tensor[tuple(idx_11)]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k11)

    if norm < 1e-15:
        msg = f"Thermal relaxation produced zero state at molecule {mol}"
        raise PhysicsViolationError(msg)
    return result / norm


def _apply_stochastic_qubit_dephasing_single(
    psi: np.ndarray, mol: int, N: int, p: float, rng: np.random.Generator
) -> np.ndarray:
    """Stochastically apply dephasing on a molecule's 2-qubit encoding.

    Dephasing projects onto a random 2-qubit computational basis state
    |q1 q0> with Born rule probabilities:
      E_deph(rho) = (1-p) rho + p sum_k P_k rho P_k
    where P_k = |k><k| on the 4-dim local space (k = 0..3).

    With prob (1-p): no dephasing.
    With prob p: project onto a random local basis state.
    """
    if p <= 0.0:
        return psi

    if rng.random() >= p:
        return psi

    d_local = 4
    dim = d_local**N
    psi_tensor = psi.reshape([d_local] * N)

    # Compute local probabilities
    local_probs = np.zeros(d_local)
    for k in range(d_local):
        idx = [slice(None)] * N
        idx[mol] = k
        local_probs[k] = np.real(np.vdot(
            psi_tensor[tuple(idx)], psi_tensor[tuple(idx)]
        ))

    total = local_probs.sum()
    if total < 1e-15:
        return psi
    local_probs /= total

    # Sample which local state to project onto
    k_chosen = rng.choice(d_local, p=local_probs)

    # Apply projection: zero out all components where mol != k_chosen
    result_tensor = np.zeros_like(psi_tensor)
    idx = [slice(None)] * N
    idx[mol] = k_chosen
    result_tensor[tuple(idx)] = psi_tensor[tuple(idx)]

    result = result_tensor.reshape(dim)
    norm = np.sqrt(np.real(np.vdot(result, result)))
    if norm < 1e-15:
        msg = f"Dephasing projection produced zero state at molecule {mol}"
        raise PhysicsViolationError(msg)
    return result / norm


class QubitGKSLShotSimulator:
    """Shot-based qubit GKSL simulator using quantum trajectories.

    Each shot evolves a pure state through the Trotter decomposition in
    the 4^N = 256 dimensional qubit space, stochastically measuring
    Stinespring ancillas. At the end, the system is measured in the
    computational basis.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QubitGKSLShotSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.N = params.N_molecules
        self.n_sys_qubits = 2 * self.N
        self.dim_qutrit = params.d ** params.N_molecules  # 81
        self.dim_qubit = (2**2) ** params.N_molecules  # 256

        # Build operators in qutrit space first
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total_qutrit = self.H_0 + self.H_transfer
        self.lindblad_ops_qutrit = build_lindblad_operators(params)

        # Qubit-qutrit mapping
        self._mapping = build_qubit_qutrit_mapping(self.N, params.d)

        # Embed operators in qubit space
        self.H_total = embed_operator_in_qubit_space(
            self.H_total_qutrit, self._mapping, self.dim_qubit
        )
        self.lindblad_ops = [
            (
                embed_operator_in_qubit_space(L, self._mapping, self.dim_qubit),
                gamma,
            )
            for L, gamma in self.lindblad_ops_qutrit
        ]

        self.n_ancilla = len(self.lindblad_ops)
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries in qubit space.

        Stores half-dt Stinespring unitaries for the palindromic ordering in
        _trotter_step_trajectory.
        """
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines_half = [
            stinespring_unitary_from_lindblad(L_op, dt / 2)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _prepare_initial_statevector(self, state_type: str) -> np.ndarray:
        """Prepare initial pure state in qubit space."""
        d = self.params.d
        N = self.params.N_molecules
        dim_qt = d**N
        psi_qt = np.zeros(dim_qt, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi_qt[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d**i) for i in range(N))
            psi_qt[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d**i) for i in range(N))
            psi_qt[index] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

        return embed_statevector_in_qubit_space(psi_qt, self._mapping, self.dim_qubit)

    def _apply_stinespring_with_measurement(
        self, psi: np.ndarray, U_stine: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """Apply Stinespring unitary and stochastically measure ancilla."""
        d_sys = len(psi)
        psi_ext = np.zeros(2 * d_sys, dtype=np.complex128)
        psi_ext[:d_sys] = psi
        psi_ext = U_stine @ psi_ext

        psi_0 = psi_ext[:d_sys]
        psi_1 = psi_ext[d_sys:]
        p_0 = np.real(np.vdot(psi_0, psi_0))
        p_1 = np.real(np.vdot(psi_1, psi_1))

        if rng.random() < p_0 / (p_0 + p_1):
            return psi_0 / np.sqrt(p_0)
        else:
            return psi_1 / np.sqrt(p_1)

    def _trotter_step_trajectory(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad ordering for trajectory.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        Matches QubitGKSLSimulator._trotter_step structure so that the
        infinite-shot limit equals the density matrix result exactly.
        """
        psi = self._U_H_half @ psi
        # Forward half-step for all Lindblad channels
        for U_stine_half in self._U_stines_half:
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
        # Reverse half-step for all Lindblad channels (palindromic)
        for U_stine_half in reversed(self._U_stines_half):
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
        psi = self._U_H_half @ psi
        return psi

    def _measure_system(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> int:
        """Measure system in computational basis (qubit space)."""
        probs = np.abs(psi) ** 2
        total = float(probs.sum())

        if abs(total - 1.0) > 1e-6:
            msg = f"State norm deviation before measurement: |sum(|psi|^2) - 1| = {abs(total - 1.0):.2e}"
            raise PhysicsViolationError(msg)

        probs = probs / total
        return int(rng.choice(len(probs), p=probs))

    def _map_outcome_to_qutrit(self, outcome_qubit: int) -> int | None:
        """Map qubit-space measurement outcome to qutrit index.

        Returns None if the outcome is in the forbidden subspace.
        """
        inv_mapping = {v: k for k, v in self._mapping.items()}
        return inv_mapping.get(outcome_qubit, None)

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        n_shots: int = 1000,
        seed: int | None = None,
    ) -> dict:
        """Run shot-based GKSL simulation in 256-dim qubit space.

        Parameters:
            t_max: total simulation time
            n_steps: number of Trotter steps
            initial_state: initial state type
            n_shots: number of quantum trajectories (shots)
            seed: random seed for reproducibility

        Returns:
            dict compatible with existing visualization functions.
            rho_final is in the 81-dim qutrit space for cross-comparison.
        """
        start = time_module.time()
        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        psi_init = self._prepare_initial_statevector(initial_state)
        rng = np.random.default_rng(seed)

        n_time_points = n_steps + 1

        # Accumulate in qutrit space for observables
        diag_accum = np.zeros((n_time_points, self.dim_qutrit))
        rho_accum = np.zeros(
            (n_time_points, self.dim_qutrit, self.dim_qutrit), dtype=np.complex128
        )
        final_counts: dict[int, int] = {}
        forbidden_count = 0

        for _shot in range(n_shots):
            psi = psi_init.copy()

            # Extract qutrit-space state for accumulation
            psi_qt = extract_statevector_from_qubit_space(
                psi, self._mapping, self.dim_qutrit
            )
            diag_accum[0] += np.abs(psi_qt) ** 2
            rho_accum[0] += np.outer(psi_qt, psi_qt.conj())

            for step in range(n_steps):
                psi = self._trotter_step_trajectory(psi, rng)
                psi_qt = extract_statevector_from_qubit_space(
                    psi, self._mapping, self.dim_qutrit
                )
                diag_accum[step + 1] += np.abs(psi_qt) ** 2
                rho_accum[step + 1] += np.outer(psi_qt, psi_qt.conj())

            outcome_qb = self._measure_system(psi, rng)
            outcome_qt = self._map_outcome_to_qutrit(outcome_qb)
            if outcome_qt is not None:
                final_counts[outcome_qt] = final_counts.get(outcome_qt, 0) + 1
            else:
                forbidden_count += 1

        diag_accum /= n_shots
        rho_accum /= n_shots

        times = [s * dt for s in range(n_time_points)]
        populations = [
            _populations_from_diagonal_qubit(diag_accum[s], self.params)
            for s in range(n_time_points)
        ]
        entropies = [
            compute_von_neumann_entropy(rho_accum[s]) for s in range(n_time_points)
        ]
        purities = [compute_purity(rho_accum[s]) for s in range(n_time_points)]
        traces = [float(np.real(np.trace(rho_accum[s]))) for s in range(n_time_points)]

        elapsed = time_module.time() - start

        # Gate count estimate: matches QubitGKSLSimulator
        # ~n_sys_qubits Rz gates (H0) + (N-1)*~10 gates (H_transfer pairs)
        # + n_ancilla * 2 * ~6 gates per Stinespring channel (palindromic)
        gates_per_step = (
            self.n_sys_qubits + len(self.params.neighbors) * 10 + self.n_ancilla * 6 * 2
        )

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_accum[-1],
            "elapsed_time": elapsed,
            "method": "qubit_gksl_shot",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_sys_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "dim_qubit_space": self.dim_qubit,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "n_shots": n_shots,
            "counts": final_counts,
            "forbidden_count": forbidden_count,
            "seed": seed,
        }


def _populations_from_diagonal_qubit(
    diag: np.ndarray, params: GKSLPhysicalParameters
) -> dict:
    """Compute populations from qutrit-space probability vector."""
    N = params.N_molecules
    d = params.d
    dim = d**N

    N_S0 = 0.0
    N_T1 = 0.0
    N_S1 = 0.0
    per_molecule: list[dict[int, float]] = [
        {0: 0.0, 1: 0.0, 2: 0.0} for _ in range(N)
    ]

    for idx in range(dim):
        p = diag[idx]
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


class QubitGKSLNoisyShotSimulator(QubitGKSLShotSimulator):
    """Shot-based qubit GKSL simulator with d=2 Pauli stochastic noise.

    Applies per-gate stochastic noise using d=4 two-qubit Pauli operators
    (tensor products of d=2 single-qubit Paulis). This correctly models
    qubit hardware noise, including forbidden-state leakage.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
        p_depol: depolarization probability per 2-qubit gate (default 0.01)
        p_dephasing: dephasing probability per gate (default 0.0)
        T1: energy relaxation time (default None = no relaxation)
        t_gate: 2-qubit gate time (default 300.0)
        depol_pair_only: if True, apply noise only after pair (2+ qubit)
            interactions, skipping single-site Lindblad channels (default False)
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.01,
        p_dephasing: float = 0.0,
        T1: float | None = None,
        t_gate: float = 300.0,
        depol_pair_only: bool = False,
    ) -> None:
        super().__init__(params)
        if p_depol < 0.0 or p_depol > 1.0:
            msg = f"p_depol must be in [0, 1], got {p_depol}"
            raise ValueError(msg)
        if p_dephasing < 0.0 or p_dephasing > 1.0:
            msg = f"p_dephasing must be in [0, 1], got {p_dephasing}"
            raise ValueError(msg)
        self.p_depol = p_depol
        self.p_dephasing = p_dephasing
        self.T1 = T1
        self.t_gate = t_gate
        self.depol_pair_only = depol_pair_only

        if T1 is not None and T1 > 0:
            self.p_reset = 1.0 - np.exp(-t_gate / T1)
        else:
            self.p_reset = 0.0

        self._lindblad_sites = self._compute_lindblad_sites()

    def _compute_lindblad_sites(self) -> list[list[int]]:
        """Determine molecule indices for each Lindblad operator."""
        sites: list[list[int]] = []
        for i, j in self.params.neighbors:
            sites.append([i, j])
            sites.append([i, j])
        for _channel in range(5):
            for mol in range(self.params.N_molecules):
                sites.append([mol])
        assert len(sites) == len(self.lindblad_ops)
        return sites

    def _trotter_step_trajectory(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad ordering and d=2 Pauli qubit noise."""
        N = self.params.N_molecules

        # --- Half Hamiltonian + transfer gate noise ---
        psi = self._U_H_half @ psi
        for i, j in self.params.neighbors:
            psi = _apply_stochastic_qubit_depolarization_pair(
                psi, i, j, N, self.p_depol, rng
            )
            if self.p_dephasing > 0.0:
                psi = _apply_stochastic_qubit_dephasing_single(
                    psi, i, N, self.p_dephasing, rng
                )
                psi = _apply_stochastic_qubit_dephasing_single(
                    psi, j, N, self.p_dephasing, rng
                )

        # --- Forward half-step Lindblad channels + per-channel noise ---
        for k, U_stine_half in enumerate(self._U_stines_half):
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
            sites = self._lindblad_sites[k]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    psi = _apply_stochastic_qubit_depolarization_single(
                        psi, sites[0], N, self.p_depol, rng
                    )
                    if self.p_dephasing > 0.0:
                        psi = _apply_stochastic_qubit_dephasing_single(
                            psi, sites[0], N, self.p_dephasing, rng
                        )
            else:
                psi = _apply_stochastic_qubit_depolarization_pair(
                    psi, sites[0], sites[1], N, self.p_depol, rng
                )
                if self.p_dephasing > 0.0:
                    psi = _apply_stochastic_qubit_dephasing_single(
                        psi, sites[0], N, self.p_dephasing, rng
                    )
                    psi = _apply_stochastic_qubit_dephasing_single(
                        psi, sites[1], N, self.p_dephasing, rng
                    )

        # --- Reverse half-step Lindblad channels + per-channel noise (palindromic) ---
        n_channels = len(self._U_stines_half)
        for k_rev in range(n_channels - 1, -1, -1):
            U_stine_half = self._U_stines_half[k_rev]
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
            sites = self._lindblad_sites[k_rev]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    psi = _apply_stochastic_qubit_depolarization_single(
                        psi, sites[0], N, self.p_depol, rng
                    )
                    if self.p_dephasing > 0.0:
                        psi = _apply_stochastic_qubit_dephasing_single(
                            psi, sites[0], N, self.p_dephasing, rng
                        )
            else:
                psi = _apply_stochastic_qubit_depolarization_pair(
                    psi, sites[0], sites[1], N, self.p_depol, rng
                )
                if self.p_dephasing > 0.0:
                    psi = _apply_stochastic_qubit_dephasing_single(
                        psi, sites[0], N, self.p_dephasing, rng
                    )
                    psi = _apply_stochastic_qubit_dephasing_single(
                        psi, sites[1], N, self.p_dephasing, rng
                    )

        # --- Half Hamiltonian + transfer gate noise ---
        psi = self._U_H_half @ psi
        for i, j in self.params.neighbors:
            psi = _apply_stochastic_qubit_depolarization_pair(
                psi, i, j, N, self.p_depol, rng
            )
            if self.p_dephasing > 0.0:
                psi = _apply_stochastic_qubit_dephasing_single(
                    psi, i, N, self.p_dephasing, rng
                )
                psi = _apply_stochastic_qubit_dephasing_single(
                    psi, j, N, self.p_dephasing, rng
                )

        # --- Thermal relaxation on all molecules ---
        if self.p_reset > 0.0:
            for mol in range(N):
                psi = _apply_stochastic_qubit_thermal_relaxation(
                    psi, mol, N, self.p_reset, rng
                )

        return psi

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        n_shots: int = 1000,
        seed: int | None = None,
    ) -> dict:
        """Run noisy shot-based simulation with qubit Pauli noise."""
        result = super().simulate(
            t_max=t_max,
            n_steps=n_steps,
            initial_state=initial_state,
            n_shots=n_shots,
            seed=seed,
        )
        result["method"] = "qubit_gksl_noisy_shot"
        result["noise_params"] = {
            "p_depol": self.p_depol,
            "p_dephasing": self.p_dephasing,
            "T1": self.T1,
            "t_gate": self.t_gate,
            "p_reset": self.p_reset,
            "depol_pair_only": self.depol_pair_only,
        }
        return result
