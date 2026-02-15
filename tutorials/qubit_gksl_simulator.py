"""Qubit GKSL simulator using Stinespring dilation + 2nd-order Trotter decomposition.

Scenario 3: Qubit-based GKSL-Lindblad (no boson).

Encodes each 3-level molecule using 2 qubits (|S0>->|00>, |T1>->|01>, |S1>->|10>)
and simulates the quantum circuit at the matrix level. The actual computation happens
in the 81-dim qutrit Hilbert space, but corresponds to a qubit circuit implementation.
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


class QubitGKSLSimulator:
    """Qubit GKSL simulator using Stinespring dilation + 2nd order Trotter.

    Simulates what a qubit quantum circuit would compute, using matrix-level
    Stinespring operations.  The actual computation happens in the 81-dim qutrit
    Hilbert space, but the approach corresponds to a qubit circuit implementation.
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
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.n_ancilla = len(self.lindblad_ops)  # 26
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla  # 34

        self.dim_qutrit = params.d ** params.N_molecules  # 81
        self.dim_qubit = (2 ** 2) ** params.N_molecules  # 256

        # Cache the qubit<->qutrit mapping
        self._mapping = self._build_qubit_qutrit_mapping()

    # ------------------------------------------------------------------
    # Qubit / qutrit mapping
    # ------------------------------------------------------------------

    def _build_qubit_qutrit_mapping(self) -> dict[int, int]:
        """Build mapping from qutrit basis index to qubit basis index."""
        d = self.params.d
        N = self.N
        qutrit_to_qubit: dict[int, int] = {}
        for idx in range(self.dim_qutrit):
            remainder = idx
            qubit_idx = 0
            for i in range(N):
                s = remainder % d
                remainder //= d
                qubit_idx += _QUTRIT_TO_QUBIT_PAIR[s] * (4 ** i)
            qutrit_to_qubit[idx] = qubit_idx
        return qutrit_to_qubit

    def _embed_in_qubit_space(self, rho_qutrit: np.ndarray) -> np.ndarray:
        """Embed 81x81 qutrit density matrix into 256x256 qubit space."""
        rho_qubit = np.zeros((self.dim_qubit, self.dim_qubit), dtype=np.complex128)
        for i_qt, i_qb in self._mapping.items():
            for j_qt, j_qb in self._mapping.items():
                rho_qubit[i_qb, j_qb] = rho_qutrit[i_qt, j_qt]
        return rho_qubit

    def _extract_from_qubit_space(self, rho_qubit: np.ndarray) -> np.ndarray:
        """Extract 81x81 qutrit density matrix from 256x256 qubit space."""
        rho_qutrit = np.zeros((self.dim_qutrit, self.dim_qutrit), dtype=np.complex128)
        for i_qt, i_qb in self._mapping.items():
            for j_qt, j_qb in self._mapping.items():
                rho_qutrit[i_qt, j_qt] = rho_qubit[i_qb, j_qb]
        return rho_qutrit

    # ------------------------------------------------------------------
    # Forbidden-state leakage check
    # ------------------------------------------------------------------

    def check_forbidden_states(
        self, rho_qubit: np.ndarray, step: int | None = None
    ) -> float:
        """Check leakage into forbidden |11> states.

        Returns the forbidden-state probability.
        """
        physical_indices = list(self._mapping.values())
        p_phys = sum(np.real(rho_qubit[i, i]) for i in physical_indices)
        p_forbidden = 1.0 - p_phys
        if p_forbidden > 1e-8:
            raise PhysicsViolationError(
                f"Forbidden state leakage at step {step}: P_forbidden={p_forbidden}"
            )
        return float(p_forbidden)

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _apply_hamiltonian_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """Apply unitary Hamiltonian evolution: rho -> e^{-iHdt} rho e^{iHdt}."""
        U = expm(-1j * self.H_total * dt)
        return U @ rho @ U.conj().T

    @staticmethod
    def _apply_lindblad_stinespring(
        rho: np.ndarray, L_op: np.ndarray, dt: float
    ) -> np.ndarray:
        """Apply a single Lindblad channel via Stinespring dilation."""
        U = stinespring_unitary_from_lindblad(L_op, dt)
        return apply_stinespring_to_density_matrix(rho, U)

    def _trotter_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """2nd-order symmetric Trotter step.

        exp(L dt) ~ exp(L_H dt/2) prod_k exp(L_k dt) exp(L_H dt/2)
        """
        # Half Hamiltonian
        rho = self._apply_hamiltonian_step(rho, dt / 2)
        # All Lindblad channels
        for L_op, _gamma in self.lindblad_ops:
            rho = self._apply_lindblad_stinespring(rho, L_op, dt)
        # Half Hamiltonian
        rho = self._apply_hamiltonian_step(rho, dt / 2)
        return rho

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix in the qutrit Hilbert space."""
        d = self.params.d
        N = self.N
        dim = d ** N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            # |1,0,...,0,1> in base-d: molecules 0 and N-1 in T1, rest in S0
            if N < 2:
                raise ValueError("edge_triplet requires N_molecules >= 2")
            index = (d ** (N - 1)) + 1
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
        """Run simulation using Stinespring + Trotter approach.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qubit circuit statistics.
        """
        start = time_module.time()

        dt = t_max / n_steps
        rho = self.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho, self.params)]
        entropies = [compute_von_neumann_entropy(rho)]
        purities = [compute_purity(rho)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho, dt)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        elapsed = time_module.time() - start

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
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": "qubit_gksl",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_sys_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
