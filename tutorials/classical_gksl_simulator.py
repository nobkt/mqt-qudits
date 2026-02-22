"""Classical GKSL-Lindblad simulator (Scenario 1: no boson interaction).

Uses Stinespring dilation + 2nd-order Trotter decomposition to integrate
the GKSL master equation.  Each step is a CPTP (Completely Positive,
Trace Preserving) map by construction, guaranteeing that density-matrix
positivity is preserved throughout the evolution.
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
from gksl_validation import PhysicsViolationError, validate_density_matrix  # noqa: F401
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)


class ClassicalGKSLSimulator:
    """Integrate the GKSL master equation for the non-boson model.

    Uses Stinespring dilation + 2nd-order symmetric Trotter decomposition:

        exp(L dt) ≈ exp(L_H dt/2) ∏_α exp(L_D_α dt) exp(L_H dt/2)

    where each Hamiltonian step is a unitary channel and each Lindblad
    channel is implemented via Stinespring dilation.  Both are CPTP by
    construction, so density-matrix positivity is guaranteed.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            raise ValueError("ClassicalGKSLSimulator is for non-boson model only")
        self.params = params

        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _apply_hamiltonian_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """Apply unitary Hamiltonian evolution: ρ → e^{-iHdt} ρ e^{iHdt}."""
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

        exp(L dt) ≈ exp(L_H dt/2) ∏_α exp(L_D_α dt) exp(L_H dt/2)
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
    # Initial state preparation
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Return an initial density matrix ρ₀ = |ψ⟩⟨ψ|."""
        d = self.params.d
        N = self.params.N_molecules
        dim = d**N

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            # |1,0,...,0,1⟩: molecules 0 and N-1 in T1, rest in S0
            psi = np.zeros(dim, dtype=np.complex128)
            index = 1 * (d ** (N - 1)) + 1
            psi[index] = 1.0
            return np.outer(psi, psi.conj())

        if state_type == "all_triplet":
            # |1111⟩
            psi = np.zeros(dim, dtype=np.complex128)
            index = sum(1 * (d**i) for i in range(N))
            psi[index] = 1.0
            return np.outer(psi, psi.conj())

        if state_type == "all_singlet":
            # |2222⟩ – all molecules in S1
            psi = np.zeros(dim, dtype=np.complex128)
            index = sum(2 * (d**i) for i in range(N))
            psi[index] = 1.0
            return np.outer(psi, psi.conj())

        raise ValueError(f"Unknown initial state type: {state_type}")

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run the GKSL time-evolution simulation.

        Parameters
        ----------
        t_max : float
            Total simulation time.
        n_steps : int
            Number of time steps (excluding t=0).
        initial_state : str
            One of ``'edge_triplet'``, ``'all_triplet'``, ``'all_singlet'``.

        Returns
        -------
        dict
            Keys: times, populations, entropy, purity, trace, rho_final,
            elapsed_time, method, params.
        """
        start = time_module.time()

        dt = t_max / n_steps
        rho = self.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        populations: list[dict] = [
            compute_populations_from_density_matrix(rho, self.params)
        ]
        entropies: list[float] = [compute_von_neumann_entropy(rho)]
        purities: list[float] = [compute_purity(rho)]
        traces: list[float] = [float(np.real(np.trace(rho)))]

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

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": "classical_gksl",
            "params": self.params.to_dict(),
        }
