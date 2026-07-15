"""Scenario 2: Classical GKSL-Lindblad WITH boson (phonon) interaction.

Extends the electronic system with phonon degrees of freedom using
Holstein-type electron-phonon coupling.

Uses the exact matrix-exponential of the GKSL Liouvillian superoperator
in the joint electronic ⊗ phonon space.  The propagator exp(L·dt) is a
CPTP map by construction (Lindblad–GKS theorem).
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_math_utils import (
    build_H_total_boson,
    build_lindblad_operators,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    extend_lindblad_operators,
    partial_trace_phonon,
    vectorize_density_matrix,
    unvectorize_density_matrix,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from stinespring_utils import build_gksl_superoperator_sparse


class ClassicalGKSLBosonSimulator:
    """GKSL master-equation solver on the joint electronic ⊗ phonon space.

    Uses the exact Liouvillian superoperator approach in the extended space:

        ρ_total(t + dt) = unvec( exp(L·dt) · vec(ρ_total(t)) )

    where L is the full GKSL Liouvillian in the extended electronic+phonon
    Hilbert space.  Electronic observables are obtained by partial-tracing
    over phonon degrees of freedom.
    """

    # A single dense (dim_total × dim_total) complex128 operator larger than
    # this is refused up-front with a clear MemoryError instead of an opaque
    # allocation failure inside the Hamiltonian builders.
    _MAX_DENSE_OPERATOR_BYTES = 2 * 1024**3

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            msg = "ClassicalGKSLBosonSimulator requires with_boson=True"
            raise ValueError(msg)
        self.params = params
        # Dimensions
        self.dim_el = params.d ** params.N_molecules
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_total = self.dim_el * self.dim_ph

        # Full Hamiltonian (electronic + phonon + coupling) and Liouvillian.
        # When g_eph == 0 the dynamics reduce exactly to the non-boson solver
        # (see ``simulate``), so the extended-space superoperator — whose dense
        # size is (dim_total^2)^2 and reaches PiB scale already for N=4,
        # n_max=2 — must not be built at all in that case.
        if params.g_eph == 0.0:
            self.H_total = None
            self.lindblad_ops = None
            self._L_super = None
        else:
            # Dimension guard: H_total and the extended Lindblad operators are
            # built as dense (dim_total × dim_total) arrays; refuse before
            # allocation if a single such array would exceed the budget.
            dense_op_bytes = 16 * self.dim_total**2
            if dense_op_bytes > self._MAX_DENSE_OPERATOR_BYTES:
                msg = (
                    f"dense operator size {dense_op_bytes / 1024**3:.1f} GiB "
                    f"(dim_total={self.dim_total}) exceeds the "
                    f"{self._MAX_DENSE_OPERATOR_BYTES / 1024**3:.1f} GiB limit; "
                    "reduce N_molecules or n_max"
                )
                raise MemoryError(msg)

            self.H_total = build_H_total_boson(params)

            # Extended Lindblad operators: L_el ⊗ I_phonon
            lindblad_ops_el = build_lindblad_operators(params)
            self.lindblad_ops = extend_lindblad_operators(lindblad_ops_el, self.dim_ph)

            # Build the full GKSL Liouvillian superoperator
            # (dim_total^2 × dim_total^2) in sparse CSR form; its own memory
            # guard raises MemoryError with the estimated size if exceeded.
            self._L_super = build_gksl_superoperator_sparse(
                self.H_total, self.lindblad_ops
            )

    # ------------------------------------------------------------------
    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic state ⊗ phonon vacuum."""
        d = self.params.d
        N = self.params.N_molecules
        dim_el = d ** N

        # Electronic initial state
        psi_el = np.zeros(dim_el, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
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
        method: str = "expm",
    ) -> dict:
        """Run GKSL simulation with boson coupling.

        Uses exact Liouvillian expm in the extended electronic+phonon space.
        Returns dict compatible with ClassicalGKSLSimulator output format.
        Electronic observables are obtained by partial-tracing over phonon DOFs.

        When g_eph == 0 the phonon degrees of freedom decouple exactly,
        so the dynamics reduce to the non-boson ClassicalGKSLSimulator.
        """
        # Exact reduction: g_eph=0 means no electron-phonon coupling,
        # so phonon DOFs factor out and we can use the cheaper non-boson solver.
        if self.params.g_eph == 0.0:
            reduced_kwargs = self.params.to_dict()
            reduced_kwargs["with_boson"] = False
            reduced_params = GKSLPhysicalParameters(**reduced_kwargs)
            from classical_gksl_simulator import ClassicalGKSLSimulator

            result = ClassicalGKSLSimulator(reduced_params).simulate(
                t_max=t_max, n_steps=n_steps, initial_state=initial_state
            )
            result["method"] = "classical_gksl_boson (g_eph=0 exact reduction)"
            return result

        start = time_module.time()
        dim = self.dim_total

        rho = self.prepare_initial_state(initial_state)
        vec_0 = vectorize_density_matrix(rho)

        # Compute exp(L·t_k)·vec(ρ₀) for all time points at once.
        dt = t_max / n_steps
        vecs = expm_multiply(
            self._L_super, vec_0,
            start=0.0, stop=t_max, num=n_steps + 1, endpoint=True,
        )

        # Extract electronic observables via partial trace over phonon
        times: list[float] = []
        populations: list[dict] = []
        entropies: list[float] = []
        purities: list[float] = []
        traces: list[float] = []

        for k in range(n_steps + 1):
            rho_total = unvectorize_density_matrix(vecs[k], dim)
            rho_el = partial_trace_phonon(rho_total, self.dim_el, self.dim_ph)

            times.append(k * dt)
            traces.append(float(np.real(np.trace(rho_el))))
            populations.append(
                compute_populations_from_density_matrix(rho_el, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        rho_final_total = unvectorize_density_matrix(vecs[-1], dim)
        rho_final_el = partial_trace_phonon(rho_final_total, self.dim_el, self.dim_ph)
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
        }
