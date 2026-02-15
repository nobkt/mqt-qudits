"""Qudit GKSL simulator using Stinespring dilation + 2nd-order Trotter decomposition.

Scenario 5: Qudit-based GKSL-Lindblad (no boson).

Uses native qutrit (d=3) encoding. Each molecule is naturally 1 qutrit,
so there are NO forbidden states (unlike the qubit encoding which wastes |11>).
The simulation uses the same 81-dim qutrit Hilbert space and Stinespring+Trotter
approach as the qubit simulator, but corresponds to a qudit circuit with fewer gates.
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from typing import TYPE_CHECKING

from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
)
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)

if TYPE_CHECKING:
    from gksl_physical_parameters import GKSLPhysicalParameters


class QuditGKSLSimulator:
    """Qudit GKSL simulator using Stinespring dilation + 2nd order Trotter.

    Uses native qutrit (d=3) encoding. No forbidden states exist.
    4 qutrits for system + 26 ancilla qubits for Lindblad channels.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QuditGKSLSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.n_system_qudits = params.N_molecules  # 4
        # 26 ancilla qubits: 2*3 TTA + 4*4 single-site + 4 ISC channels
        self.n_ancilla_qubits = 26

        # Build operators in native qutrit space
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        # System dimension (qutrit space)
        self.dim = params.d**params.N_molecules  # 81

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _apply_hamiltonian_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """Apply unitary Hamiltonian evolution: rho -> e^{-iHdt} rho e^{iHdt}."""
        U = expm(-1j * self.H_total * dt)
        return U @ rho @ U.conj().T

    @staticmethod
    def _apply_lindblad_stinespring(rho: np.ndarray, L_op: np.ndarray, dt: float) -> np.ndarray:
        """Apply single Lindblad channel via Stinespring dilation."""
        U = stinespring_unitary_from_lindblad(L_op, dt)
        return apply_stinespring_to_density_matrix(rho, U)

    def _trotter_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """2nd-order symmetric Trotter step.

        exp(L dt) ~ exp(L_H dt/2) prod_alpha exp(L_D_alpha dt) exp(L_H dt/2)
        """
        # Half Hamiltonian
        rho = self._apply_hamiltonian_step(rho, dt / 2)
        # All Lindblad channels via Stinespring
        for L_op, _gamma in self.lindblad_ops:
            rho = self._apply_lindblad_stinespring(rho, L_op, dt)
        # Half Hamiltonian
        return self._apply_hamiltonian_step(rho, dt / 2)

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix in qutrit space."""
        d = self.params.d
        N = self.params.N_molecules
        dim = d**N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            # Molecules 0 and N-1 in T1, rest in S0
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = (d ** (N - 1)) + 1
            psi[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d**i) for i in range(N))
            psi[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d**i) for i in range(N))
            psi[index] = 1.0
        else:
            msg = f"Unknown state type: {state_type}"
            raise ValueError(msg)

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
            populations.append(compute_populations_from_density_matrix(rho, self.params))
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        elapsed = time_module.time() - start

        # Gate count estimate for qudit circuit
        # Qutrit advantages: no encoding overhead, native 3-level operations
        # 4 VirtRz gates (H_0 diagonal for 4 molecules)
        # 3 CustomTwo gates (H_transfer for 3 nearest-neighbour pairs)
        # 26 Stinespring CustomTwo gates (one per Lindblad channel)
        gates_per_step = 4 + 3 + 26  # 33 gates (much fewer than qubit ~194)

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": "qudit_gksl",
            "params": self.params.to_dict(),
            "n_system_qudits": self.n_system_qudits,
            "n_ancilla_qubits": self.n_ancilla_qubits,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
