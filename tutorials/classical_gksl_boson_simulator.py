"""Scenario 2: Classical GKSL-Lindblad WITH boson (phonon) interaction.

Extends the electronic system with phonon degrees of freedom using
Holstein-type electron-phonon coupling.
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.integrate import solve_ivp

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


class ClassicalGKSLBosonSimulator:
    """GKSL master-equation solver on the joint electronic ⊗ phonon space."""

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            msg = "ClassicalGKSLBosonSimulator requires with_boson=True"
            raise ValueError(msg)
        self.params = params

        # Dimensions
        self.dim_el = params.d ** params.N_molecules
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_total = self.dim_el * self.dim_ph

        # Full Hamiltonian (electronic + phonon + coupling)
        self.H_total = build_H_total_boson(params)

        # Extended Lindblad operators: L_el ⊗ I_phonon
        lindblad_ops_el = build_lindblad_operators(params)
        self.lindblad_ops = extend_lindblad_operators(lindblad_ops_el, self.dim_ph)

    # ------------------------------------------------------------------
    def _gksl_rhs(self, rho: np.ndarray) -> np.ndarray:
        """Compute dρ/dt directly in matrix form."""
        H = self.H_total
        drho = -1j * (H @ rho - rho @ H)
        for L_op, _gamma in self.lindblad_ops:
            LdL = L_op.conj().T @ L_op
            drho += L_op @ rho @ L_op.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
        return drho

    # ------------------------------------------------------------------
    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic state ⊗ phonon vacuum."""
        d = self.params.d
        N = self.params.N_molecules
        dim_el = d ** N

        # Electronic initial state
        psi_el = np.zeros(dim_el, dtype=np.complex128)
        if state_type == "edge_triplet":
            # Edge molecules in T1 (1), middle ones in S0 (0): |1,0,...,0,1⟩
            index = 1 * (d ** (N - 1)) + 1
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

        # Phonon vacuum state |0…0⟩
        psi_ph = np.zeros(self.dim_ph, dtype=np.complex128)
        psi_ph[0] = 1.0

        # Full initial state
        psi_full = np.kron(psi_el, psi_ph)
        return np.outer(psi_full, psi_full.conj())

    # ------------------------------------------------------------------
    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        method: str = "BDF",
    ) -> dict:
        """Run GKSL simulation with boson coupling.

        Uses BDF (stiff solver) since phonon coupling creates stiff dynamics.
        Returns dict compatible with ClassicalGKSLSimulator output format.
        """
        start = time_module.time()
        dim = self.dim_total

        rho_0 = self.prepare_initial_state(initial_state)
        y0 = rho_0.flatten()

        t_eval = np.linspace(0, t_max, n_steps + 1)

        def ode_func(_t: float, y: np.ndarray) -> np.ndarray:
            rho = y.reshape((dim, dim))
            return self._gksl_rhs(rho).flatten()

        sol = solve_ivp(
            ode_func,
            [0, t_max],
            y0,
            t_eval=t_eval,
            method=method,
            rtol=1e-8,
            atol=1e-10,
        )

        if not sol.success:
            msg = f"ODE solver failed: {sol.message}"
            raise RuntimeError(msg)

        # Extract electronic observables via partial trace over phonon
        times: list[float] = sol.t.tolist()
        populations: list[dict] = []
        entropies: list[float] = []
        purities: list[float] = []
        traces: list[float] = []

        for k in range(len(times)):
            rho_total = sol.y[:, k].reshape((dim, dim))
            rho_el = partial_trace_phonon(rho_total, self.dim_el, self.dim_ph)

            traces.append(float(np.real(np.trace(rho_el))))
            populations.append(
                compute_populations_from_density_matrix(rho_el, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        rho_final = sol.y[:, -1].reshape((dim, dim))
        rho_final_el = partial_trace_phonon(rho_final, self.dim_el, self.dim_ph)
        elapsed = time_module.time() - start

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_final_el,
            "elapsed_time": elapsed,
            "method": "classical_gksl_boson",
            "params": self.params.to_dict(),
            "ode_solver": method,
            "n_function_evals": sol.nfev,
        }
