"""Qudit GKSL simulator with boson (phonon) interaction using Stinespring + Trotter.

Scenario 6: Qudit-based GKSL-Lindblad with boson.

Uses native qutrit (d=3) encoding for both electronic and phonon degrees of freedom.
No forbidden states for either subsystem. Works in the extended electronic+phonon
Hilbert space dim_total = d^N * (n_max+1)^N, using Stinespring + 2nd-order Trotter.
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import expm

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
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)


class QuditGKSLBosonSimulator:
    """Qudit GKSL simulator with boson (phonon) interaction.

    Uses native qutrit (d=3) encoding for electronic states and qutrit encoding
    for phonon modes (when n_max=2). No forbidden states for either subsystem.
    Works in the extended space dim_total = d^N * (n_max+1)^N.

    Since this targets qudit quantum computers, all registers — including
    Stinespring ancillas — are native d-level qudits.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            raise ValueError("QuditGKSLBosonSimulator requires with_boson=True")
        self.params = params
        self.N = params.N_molecules
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension

        # Dimensions
        self.dim_el = params.d ** params.N_molecules
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_total = self.dim_el * self.dim_ph

        # Qudit counts
        self.n_system_qudits = params.N_molecules  # electronic qutrits
        self.n_phonon_qudits = params.N_molecules   # phonon qutrits (for n_max=2)
        self.n_ancilla_qudits = 2 * len(params.neighbors) + 5 * params.N_molecules
        self.n_total_qudits = self.n_system_qudits + self.n_phonon_qudits

        # Build extended Hamiltonian
        self.H_total = build_H_total_boson(params)

        # Build extended Lindblad operators
        lindblad_ops_el = build_lindblad_operators(params)
        self.lindblad_ops = extend_lindblad_operators(lindblad_ops_el, self.dim_ph)

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation).

        Ancilla dimension is d_anc (= params.d).
        """
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines_half = [
            stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad channel ordering.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        The Hamiltonian–Dissipator Strang splitting is 2nd-order O(dt³)/step.
        The palindromic Lindblad product eliminates the Lie-Trotter commutator
        error (also 2nd-order).  The remaining dominant error is the Stinespring
        approximation: O(dt²)/step, giving **O(dt) global convergence**.
        """
        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        # Forward half-step for all Lindblad channels
        for U_stine_half in self._U_stines_half:
            rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
        # Reverse half-step for all Lindblad channels (palindromic)
        for U_stine_half in reversed(self._U_stines_half):
            rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
        # Half Hamiltonian
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        return rho

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic state ⊗ phonon vacuum |00...0⟩."""
        d = self.params.d
        N = self.N
        dim_el = d ** N

        # Electronic state
        psi_el = np.zeros(dim_el, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi_el[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d ** i) for i in range(N))
            psi_el[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d ** i) for i in range(N))
            psi_el[index] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

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
        """Run Qudit GKSL+boson simulation.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qudit circuit statistics. Populations are computed by
        partial-tracing over phonon degrees of freedom.
        """
        start = time_module.time()

        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        # Partial trace for initial populations
        rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_el, self.params)]
        entropies = [compute_von_neumann_entropy(rho_el)]
        purities = [compute_purity(rho_el)]
        traces = [float(np.real(np.trace(rho_el)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_el))))
            populations.append(
                compute_populations_from_density_matrix(rho_el, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        elapsed = time_module.time() - start

        # Return electronic-only reduced density matrix for consistency
        # with ClassicalGKSLBosonSimulator and to match time-series observables
        rho_el_final = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

        # Gate count estimate for qudit circuit in extended space
        # (palindromic 2nd-order Trotter)
        # Qudit advantage: native d-level ops, no forbidden states
        # 2 × N cu_one gates (electronic on-site, two half-steps)
        # 2 × N cu_one gates (phonon on-site, two half-steps)
        # 2 × (N-1) cu_two gates (H_transfer, two half-steps)
        # 2 × n_lindblad Stinespring gates (fwd + rev)
        gates_per_step = (
            2 * (self.n_system_qudits + self.n_phonon_qudits + len(self.params.neighbors))
            + self.n_ancilla_qudits * 2
        )

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_el_final,
            "elapsed_time": elapsed,
            "method": "qudit_gksl_boson",
            "params": self.params.to_dict(),
            "n_system_qudits": self.n_system_qudits,
            "n_phonon_qudits": self.n_phonon_qudits,
            "n_ancilla_qudits": self.n_ancilla_qudits,
            "d_anc": self.d_anc,
            "n_total_qudits": self.n_total_qudits,
            "dim_total": self.dim_total,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
