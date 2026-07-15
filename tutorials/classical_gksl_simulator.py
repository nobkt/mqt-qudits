"""Classical GKSL-Lindblad simulator (Scenario 1: no boson interaction).

Uses the exact matrix-exponential of the GKSL Liouvillian superoperator to
integrate the master equation.  The propagator exp(L·t) is a CPTP map by
construction (Lindblad–GKS theorem), so density-matrix positivity and trace
preservation are structurally guaranteed.

This provides an independent reference solution that does NOT use the
Stinespring + Trotter decomposition employed by the qubit/qudit simulators.
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    unvectorize_density_matrix,
    vectorize_density_matrix,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from gksl_validation import PhysicsViolationError, validate_density_matrix  # noqa: F401
from stinespring_utils import build_gksl_superoperator


class ClassicalGKSLSimulator:
    """Integrate the GKSL master equation for the non-boson model.

    Uses the exact Liouvillian superoperator approach:

        vec(ρ(t)) = exp(L·t) · vec(ρ(0))

    where L is the full GKSL Liouvillian.  The matrix exponential exp(L·t)
    is always a CPTP map for a valid Lindblad generator L (Lindblad–GKS
    theorem), so density-matrix positivity is guaranteed without Trotter
    splitting.  This provides an independent reference to compare against
    the Stinespring + Trotter decomposition used by the qubit/qudit simulators.

    The time evolution is computed via ``scipy.sparse.linalg.expm_multiply``
    which evaluates exp(L·t)·v efficiently without forming the full matrix
    exponential.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            raise ValueError("ClassicalGKSLSimulator is for non-boson model only")
        self.params = params

        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.dim = params.d ** params.N_molecules

        # Build the full GKSL Liouvillian superoperator (dim^2 × dim^2)
        self._L_super = build_gksl_superoperator(self.H_total, self.lindblad_ops)

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

        rho_0 = self.prepare_initial_state(initial_state)
        vec_0 = vectorize_density_matrix(rho_0)

        # Compute exp(L·t_k)·vec(ρ₀) for all time points at once.
        # This is the exact formal solution of the GKSL master equation.
        vecs = expm_multiply(
            self._L_super, vec_0,
            start=0.0, stop=t_max, num=n_steps + 1, endpoint=True,
        )

        times: list[float] = []
        populations: list[dict] = []
        entropies: list[float] = []
        purities: list[float] = []
        traces: list[float] = []

        dt = t_max / n_steps
        for k in range(n_steps + 1):
            rho = unvectorize_density_matrix(vecs[k], self.dim)

            times.append(k * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        rho_final = unvectorize_density_matrix(vecs[-1], self.dim)
        elapsed = time_module.time() - start

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_final,
            "elapsed_time": elapsed,
            "method": "classical_gksl",
            "params": self.params.to_dict(),
        }
