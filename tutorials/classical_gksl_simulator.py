"""Classical GKSL-Lindblad simulator (Scenario 1: no boson interaction).

Uses scipy.integrate.solve_ivp to integrate the GKSL master equation
directly in density-matrix form (not the full superoperator), keeping
memory usage manageable for 4-molecule, 3-level systems (dim=81).
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.integrate import solve_ivp

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


class ClassicalGKSLSimulator:
    """Integrate the GKSL master equation for the non-boson model.

    The right-hand side is evaluated with the direct matrix approach:

        dρ/dt = -i[H, ρ] + Σ_α (L_α ρ L_α† − ½{L_α†L_α, ρ})

    where each L_α already contains the √γ factor.
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
    # RHS of the master equation
    # ------------------------------------------------------------------

    def _gksl_rhs(self, rho: np.ndarray) -> np.ndarray:
        """Compute dρ/dt = -i[H, ρ] + Σ_α D[L_α](ρ)."""
        H = self.H_total
        drho = -1j * (H @ rho - rho @ H)
        for L_op, _gamma in self.lindblad_ops:
            LdL = L_op.conj().T @ L_op
            drho += L_op @ rho @ L_op.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
        return drho

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
        dim = rho_0.shape[0]

        y0 = vectorize_density_matrix(rho_0)
        t_eval = np.linspace(0, t_max, n_steps + 1)

        def ode_func(_t: float, y: np.ndarray) -> np.ndarray:
            rho = unvectorize_density_matrix(y, dim)
            drho = self._gksl_rhs(rho)
            return vectorize_density_matrix(drho)

        sol = solve_ivp(
            ode_func,
            [0, t_max],
            y0,
            t_eval=t_eval,
            method="RK45",
            rtol=1e-9,
            atol=1e-12,
        )

        if not sol.success:
            raise RuntimeError(f"ODE solver failed: {sol.message}")

        # ---- Post-processing ----
        times: list[float] = sol.t.tolist()
        populations: list[dict] = []
        entropies: list[float] = []
        purities: list[float] = []
        traces: list[float] = []

        for k in range(len(times)):
            rho = unvectorize_density_matrix(sol.y[:, k], dim)

            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        rho_final = unvectorize_density_matrix(sol.y[:, -1], dim)
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
