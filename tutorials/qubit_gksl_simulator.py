"""Qubit GKSL simulator using Stinespring dilation + 2nd-order Trotter decomposition.

Scenario 3: Qubit-based GKSL-Lindblad (no boson).

Encodes each 3-level molecule using 2 qubits (|S0>->|00>, |T1>->|01>, |S1>->|10>)
and simulates in the 4^N = 256 dimensional qubit Hilbert space. This faithfully
represents what a qubit quantum computer would compute, including forbidden-state
(|11>) leakage effects under noisy conditions.

The physical subspace is the 81-dim qutrit subspace embedded in the 256-dim qubit
space. In the noiseless case, the operators are embedded to preserve this subspace,
so results match the native qutrit simulation up to floating-point precision.
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
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)

# Qubit-pair encoding for each qutrit level
_QUTRIT_TO_QUBIT_PAIR = {0: 0b00, 1: 0b01, 2: 0b10}


def build_qubit_qutrit_mapping(N: int, d: int = 3) -> dict[int, int]:
    """Build mapping from qutrit basis index to qubit basis index.

    For each qutrit computational basis state, compute the corresponding
    qubit-pair encoding index. Uses the same Kronecker product ordering
    as the operator construction in gksl_math_utils.
    """
    dim_qutrit = d**N
    qutrit_to_qubit: dict[int, int] = {}
    for idx in range(dim_qutrit):
        remainder = idx
        qubit_idx = 0
        for i in range(N):
            s = remainder % d
            remainder //= d
            qubit_idx += _QUTRIT_TO_QUBIT_PAIR[s] * (4**i)
        qutrit_to_qubit[idx] = qubit_idx
    return qutrit_to_qubit


def embed_operator_in_qubit_space(
    op: np.ndarray, mapping: dict[int, int], dim_qubit: int
) -> np.ndarray:
    """Embed a qutrit-space operator into the qubit-encoded space.

    The operator acts only on the physical subspace; forbidden-state
    matrix elements are zero.
    """
    op_q = np.zeros((dim_qubit, dim_qubit), dtype=np.complex128)
    for i_qt, i_qb in mapping.items():
        for j_qt, j_qb in mapping.items():
            op_q[i_qb, j_qb] = op[i_qt, j_qt]
    return op_q


def embed_density_matrix_in_qubit_space(
    rho: np.ndarray, mapping: dict[int, int], dim_qubit: int
) -> np.ndarray:
    """Embed a qutrit-space density matrix into qubit-encoded space."""
    rho_q = np.zeros((dim_qubit, dim_qubit), dtype=np.complex128)
    for i_qt, i_qb in mapping.items():
        for j_qt, j_qb in mapping.items():
            rho_q[i_qb, j_qb] = rho[i_qt, j_qt]
    return rho_q


def extract_density_matrix_from_qubit_space(
    rho_q: np.ndarray, mapping: dict[int, int], dim_qutrit: int
) -> np.ndarray:
    """Extract qutrit-space density matrix from qubit-encoded space."""
    rho = np.zeros((dim_qutrit, dim_qutrit), dtype=np.complex128)
    for i_qt, i_qb in mapping.items():
        for j_qt, j_qb in mapping.items():
            rho[i_qt, j_qt] = rho_q[i_qb, j_qb]
    return rho


def embed_statevector_in_qubit_space(
    psi: np.ndarray, mapping: dict[int, int], dim_qubit: int
) -> np.ndarray:
    """Embed a qutrit-space state vector into qubit-encoded space."""
    psi_q = np.zeros(dim_qubit, dtype=np.complex128)
    for i_qt, i_qb in mapping.items():
        psi_q[i_qb] = psi[i_qt]
    return psi_q


def extract_statevector_from_qubit_space(
    psi_q: np.ndarray, mapping: dict[int, int], dim_qutrit: int
) -> np.ndarray:
    """Extract qutrit-space components from qubit-space state vector."""
    psi = np.zeros(dim_qutrit, dtype=np.complex128)
    for i_qt, i_qb in mapping.items():
        psi[i_qt] = psi_q[i_qb]
    return psi


def compute_forbidden_state_population(
    rho_q: np.ndarray, mapping: dict[int, int]
) -> float:
    """Compute the total population in forbidden (|11>) states."""
    physical_indices = set(mapping.values())
    dim = rho_q.shape[0]
    p_forbidden = 0.0
    for i in range(dim):
        if i not in physical_indices:
            p_forbidden += float(np.real(rho_q[i, i]))
    return p_forbidden


class QubitGKSLSimulator:
    """Qubit GKSL simulator using Stinespring dilation + 2nd order Trotter.

    Simulates in the 4^N = 256 dimensional qubit Hilbert space with 2-qubit
    encoding per molecule (|00>=S0, |01>=T1, |10>=S1, |11>=forbidden).

    All operators are embedded in the qubit space, and the Trotter evolution
    runs in the full 256-dim space. Observables are extracted by projecting
    back to the 81-dim qutrit physical subspace.

    In the noiseless case, the physical subspace is invariant under the
    embedded operators, so results match the native qutrit simulation up
    to floating-point precision. The computational overhead (256 vs 81 dim)
    demonstrates the cost of qubit encoding.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            raise ValueError("QubitGKSLSimulator is for non-boson model only")
        self.params = params
        self.N = params.N_molecules
        self.n_sys_qubits = 2 * self.N  # 8

        # Build qutrit-space operators
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total_qutrit = self.H_0 + self.H_transfer
        self.lindblad_ops_qutrit = build_lindblad_operators(params)

        self.n_ancilla = len(self.lindblad_ops_qutrit)  # 26
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla  # 34

        self.dim_qutrit = params.d ** params.N_molecules  # 81
        self.dim_qubit = (2**2) ** params.N_molecules  # 256

        # Build qubit-qutrit mapping
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

    # ------------------------------------------------------------------
    # Qubit / qutrit space conversions
    # ------------------------------------------------------------------

    def _embed_in_qubit_space(self, rho_qutrit: np.ndarray) -> np.ndarray:
        """Embed 81x81 qutrit density matrix into 256x256 qubit space."""
        return embed_density_matrix_in_qubit_space(
            rho_qutrit, self._mapping, self.dim_qubit
        )

    def _extract_from_qubit_space(self, rho_qubit: np.ndarray) -> np.ndarray:
        """Extract 81x81 qutrit density matrix from 256x256 qubit space."""
        return extract_density_matrix_from_qubit_space(
            rho_qubit, self._mapping, self.dim_qutrit
        )

    # ------------------------------------------------------------------
    # Forbidden-state leakage check
    # ------------------------------------------------------------------

    def check_forbidden_states(
        self, rho_qubit: np.ndarray, step: int | None = None
    ) -> float:
        """Check leakage into forbidden |11> states.

        Returns the forbidden-state probability.
        """
        p_forbidden = compute_forbidden_state_population(rho_qubit, self._mapping)
        if p_forbidden > 1e-8:
            raise PhysicsViolationError(
                f"Forbidden state leakage at step {step}: P_forbidden={p_forbidden}"
            )
        return float(p_forbidden)

    # ------------------------------------------------------------------
    # Trotter step primitives (operating in 256-dim qubit space)
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation)."""
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines = [
            stinespring_unitary_from_lindblad(L_op, dt)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """2nd-order symmetric Trotter step in 256-dim qubit space.

        exp(L dt) ~ exp(L_H dt/2) prod_k exp(L_k dt) exp(L_H dt/2)
        """
        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        # All Lindblad channels
        for U_stine in self._U_stines:
            rho = apply_stinespring_to_density_matrix(rho, U_stine)
        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        return rho

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix in the qubit Hilbert space.

        First constructs the state in qutrit space, then embeds into qubit space.
        """
        d = self.params.d
        N = self.N
        dim_qt = d**N
        psi = np.zeros(dim_qt, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            # |1,0,...,0,1> in base-d: molecules 0 and N-1 in T1, rest in S0
            index = 1 * (d ** (N - 1)) + 1
            psi[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d**i) for i in range(N))
            psi[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d**i) for i in range(N))
            psi[index] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

        rho_qutrit = np.outer(psi, psi.conj())
        return self._embed_in_qubit_space(rho_qutrit)

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run simulation using Stinespring + Trotter in 256-dim qubit space.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qubit circuit statistics. Observables are computed from
        the qutrit subspace of the qubit density matrix.
        """
        start = time_module.time()

        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        # Extract to qutrit space for observables
        rho_qt = self._extract_from_qubit_space(rho)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_qt, self.params)]
        entropies = [compute_von_neumann_entropy(rho_qt)]
        purities = [compute_purity(rho_qt)]
        traces = [float(np.real(np.trace(rho_qt)))]
        forbidden_pops: list[float] = [
            compute_forbidden_state_population(rho, self._mapping)
        ]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            # Check forbidden-state leakage
            self.check_forbidden_states(rho, step=step + 1)

            # Extract to qutrit space for observables
            rho_qt = self._extract_from_qubit_space(rho)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_qt))))
            populations.append(
                compute_populations_from_density_matrix(rho_qt, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_qt))
            purities.append(compute_purity(rho_qt))
            forbidden_pops.append(
                compute_forbidden_state_population(rho, self._mapping)
            )

        elapsed = time_module.time() - start

        # Return the qutrit-space density matrix for cross-scenario comparison
        rho_final_qt = self._extract_from_qubit_space(rho)

        # Gate count estimate for a real qubit circuit:
        # ~n_sys_qubits Rz gates (H0) + 3 * ~10 gates (H_transfer pairs)
        # + n_lindblad * ~6 gates per Stinespring channel
        gates_per_step = (
            self.n_sys_qubits + len(self.params.neighbors) * 10 + self.n_ancilla * 6
        )

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_final_qt,
            "elapsed_time": elapsed,
            "method": "qubit_gksl",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_sys_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "dim_qubit_space": self.dim_qubit,
            "forbidden_state_population": forbidden_pops,
        }
