"""Qudit GKSL simulator using Stinespring dilation + Trotter decomposition.

Scenario 5: Qudit-based GKSL-Lindblad (no boson).

Uses native qutrit (d=3) encoding. Each molecule is naturally 1 qutrit,
so there are NO forbidden states (unlike the qubit encoding which wastes |11>).
The simulation uses the same 81-dim qutrit Hilbert space and Stinespring+Trotter
approach as the qubit simulator, but corresponds to a qudit circuit with fewer gates.

Convergence note
----------------
The Hamiltonian–Dissipator splitting uses Strang (symmetric) splitting, which
is 2nd-order for the H-D decomposition.  However, the individual Lindblad
channels are applied via Stinespring dilation, which is a 1st-order
approximation of each exact Lindblad channel exp(L_{D_α} dt).  Additionally,
the sequential product of individual channel maps introduces a Lie-Trotter
product error.  The Lindblad channels are applied in *symmetric (palindromic)
order* to eliminate the Lie-Trotter commutator error, but the Stinespring
approximation itself remains 1st-order.  As a result, the effective
convergence in trace distance is **O(dt)** (1st-order), not O(dt²).
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
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)


class QuditGKSLSimulator:
    """Qudit GKSL simulator using Stinespring dilation + Trotter splitting.

    Uses native qutrit (d=3) encoding. No forbidden states exist.
    4 qutrits for system + 26 ancilla qudits (d=3) for Lindblad channels.

    Since this simulator targets qudit-type quantum computers (e.g.
    qudit-boson ion-trap processors), all registers — including Stinespring
    ancillas — are native d-level qudits, not 2-level qubits.

    The Hamiltonian–Dissipator Strang splitting is 2nd-order, but the
    Stinespring dilation is a 1st-order approximation of each Lindblad
    channel.  Lindblad channels are applied in symmetric (palindromic)
    order to eliminate the Lie-Trotter product commutator error.
    Effective convergence in trace distance is O(dt) (1st-order).
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            raise ValueError("QuditGKSLSimulator is for non-boson model only")
        self.params = params
        self.n_system_qudits = params.N_molecules
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension

        # Build operators in native qutrit space
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.n_ancilla_qudits = len(self.lindblad_ops)

        # System dimension (qutrit space)
        self.dim = params.d ** params.N_molecules  # 81

    # ------------------------------------------------------------------
    # Trotter step primitives
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation).

        Stores half-dt Stinespring unitaries for the symmetric palindromic
        product used in _trotter_step.  Ancilla dimension is d_anc (= params.d).
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
        approximation: E_α(dt/2)(ρ) = exp(L_{D_α} dt/2)(ρ) + O(dt²), giving
        O(dt²)/step and **O(dt) global convergence** in trace distance.
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
        """Prepare initial density matrix in qutrit space."""
        d = self.params.d
        N = self.params.N_molecules
        dim = d ** N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            # Molecules 0 and N-1 in T1, rest in S0
            index = 1 * (d ** (N - 1)) + 1
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
        """Run Qudit GKSL simulation.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qudit circuit statistics.
        """
        start = time_module.time()

        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho, self.params)]
        entropies = [compute_von_neumann_entropy(rho)]
        purities = [compute_purity(rho)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        elapsed = time_module.time() - start

        # Gate count estimate for qudit circuit
        # Qutrit advantages: no encoding overhead, native 3-level operations
        # N VirtRz gates (H_0 diagonal for N molecules)
        # (N-1) CustomTwo gates (H_transfer for nearest-neighbour pairs)
        # n_lindblad×2 Stinespring CustomTwo gates (palindromic: forward + reverse)
        n_lindblad = len(self.lindblad_ops)
        gates_per_step = self.n_system_qudits + len(self.params.neighbors) + n_lindblad * 2

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
            "n_ancilla_qudits": self.n_ancilla_qudits,
            "d_anc": self.d_anc,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
