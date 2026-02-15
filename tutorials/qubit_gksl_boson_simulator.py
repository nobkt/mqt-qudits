"""Qubit GKSL simulator with boson (phonon) interaction using Stinespring + Trotter.

Scenario 4: Qubit-based GKSL-Lindblad with boson.

Encodes each 3-level molecule using 2 qubits and phonon modes using ceil(log2(n_max+1))
qubits per molecule. Works in the extended electronic+phonon Hilbert space
dim_total = d^N * (n_max+1)^N, using the same Stinespring + 2nd-order Trotter approach.
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
    build_H_total_boson,
    build_lindblad_operators,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    extend_lindblad_operators,
    partial_trace_phonon,
)
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)

if TYPE_CHECKING:
    from gksl_physical_parameters import GKSLPhysicalParameters


class QubitGKSLBosonSimulator:
    """Qubit GKSL simulator with boson (phonon) interaction.

    Works in the extended space dim_total = d^N * (n_max+1)^N.
    Uses Stinespring + 2nd-order Trotter in this extended space.
    After simulation, partial-traces over phonon for population calculation.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            msg = "QubitGKSLBosonSimulator requires with_boson=True"
            raise ValueError(msg)
        self.params = params
        self.N = params.N_molecules

        # Dimensions
        self.dim_el = params.d**params.N_molecules
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_total = self.dim_el * self.dim_ph

        # Qubit counts
        self.n_el_qubits = 2 * self.N
        self.n_ph_qubits = int(np.ceil(np.log2(params.n_max + 1))) * self.N
        self.n_ancilla = 2 * len(params.neighbors) + 5 * params.N_molecules  # 26 for N=4
        self.n_total_qubits = self.n_el_qubits + self.n_ph_qubits + self.n_ancilla

        # Build extended Hamiltonian
        self.H_total = build_H_total_boson(params)

        # Build extended Lindblad operators
        lindblad_ops_el = build_lindblad_operators(params)
        self.lindblad_ops = extend_lindblad_operators(lindblad_ops_el, self.dim_ph)

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _apply_hamiltonian_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """Apply unitary Hamiltonian evolution: rho -> e^{-iHdt} rho e^{iHdt}."""
        U = expm(-1j * self.H_total * dt)
        return U @ rho @ U.conj().T

    @staticmethod
    def _apply_lindblad_stinespring(rho: np.ndarray, L_op: np.ndarray, dt: float) -> np.ndarray:
        """Apply a single Lindblad channel via Stinespring dilation."""
        U = stinespring_unitary_from_lindblad(L_op, dt)
        return apply_stinespring_to_density_matrix(rho, U)

    def _trotter_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """2nd-order symmetric Trotter step in extended space.

        exp(L dt) ~ exp(L_H dt/2) prod_k exp(L_k dt) exp(L_H dt/2)
        """
        # Half Hamiltonian
        rho = self._apply_hamiltonian_step(rho, dt / 2)
        # All Lindblad channels
        for L_op, _gamma in self.lindblad_ops:
            rho = self._apply_lindblad_stinespring(rho, L_op, dt)
        # Half Hamiltonian
        return self._apply_hamiltonian_step(rho, dt / 2)

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic state ⊗ phonon vacuum |00...0⟩."""
        d = self.params.d
        N = self.N
        dim_el = d**N

        # Electronic state
        psi_el = np.zeros(dim_el, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = (d ** (N - 1)) + 1
            psi_el[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d**i) for i in range(N))
            psi_el[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d**i) for i in range(N))
            psi_el[index] = 1.0
        else:
            msg = f"Unknown state type: {state_type}"
            raise ValueError(msg)

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
        """Run simulation using Stinespring + Trotter in extended space.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qubit circuit statistics. Populations are computed by
        partial-tracing over phonon degrees of freedom.
        """
        start = time_module.time()

        dt = t_max / n_steps
        rho = self.prepare_initial_state(initial_state)

        # Partial trace for initial populations
        rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_el, self.params)]
        entropies = [compute_von_neumann_entropy(rho_el)]
        purities = [compute_purity(rho_el)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho, dt)

            rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(compute_populations_from_density_matrix(rho_el, self.params))
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        elapsed = time_module.time() - start

        # Gate count estimate for qubit circuit in extended space
        gates_per_step = self.n_el_qubits + self.n_ph_qubits + len(self.params.neighbors) * 10 + self.n_ancilla * 6

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": "qubit_gksl_boson",
            "params": self.params.to_dict(),
            "n_el_qubits": self.n_el_qubits,
            "n_ph_qubits": self.n_ph_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "dim_total": self.dim_total,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
