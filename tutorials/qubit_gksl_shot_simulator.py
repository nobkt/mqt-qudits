"""Qubit GKSL shot-based simulator using quantum trajectories.

Implements shot-based (quantum trajectory) simulation for the qubit-encoded
GKSL-Lindblad dynamics. Each molecule is encoded in 2 qubits
(|00>=S0, |01>=T1, |10>=S1, |11>=forbidden). The computation runs in the
81-dim qutrit Hilbert space (same as the density matrix version), but
corresponds to a qubit circuit implementation.

Two classes are provided:
  - QubitGKSLShotSimulator: ideal shot-based (no hardware noise)
  - QubitGKSLNoisyShotSimulator: shot-based with per-gate stochastic noise
    (depolarization + thermal relaxation)
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
    compute_purity,
    compute_von_neumann_entropy,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from gksl_validation import PhysicsViolationError
from qudit_gksl_shot_simulator import (
    _apply_stochastic_depolarization_pair,
    _apply_stochastic_depolarization_single,
    _populations_from_diagonal,
)
from stinespring_utils import stinespring_unitary_from_lindblad


def _apply_stochastic_thermal_relaxation_single(
    psi: np.ndarray, site: int, d: int, N: int, p_reset: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stochastically apply thermal relaxation on a single molecule.

    Kraus operators (qutrit encoding |0>=S0, |1>=T1, |2>=S1):
      K0 = diag(1, sqrt(1-p), sqrt(1-p))  (no decay)
      K1 = sqrt(p) * |0><1|               (T1 -> S0)
      K2 = sqrt(p) * |0><2|               (S1 -> S0)

    Stochastic application: sample which Kraus operator to apply with Born rule.
    """
    if p_reset <= 0.0:
        return psi

    dim = d ** N
    psi_tensor = psi.reshape([d] * N)

    # Compute probabilities for each Kraus operator
    sq = np.sqrt(1.0 - p_reset)
    sqp = np.sqrt(p_reset)

    # Prob(K0): sum over all basis states of |K0[site_state] * amplitude|^2
    # K0 multiplies |0> by 1, |1> by sq, |2> by sq
    k0_diag = np.array([1.0, sq, sq])
    p_k0 = 0.0
    for k in range(d):
        idx = [slice(None)] * N
        idx[site] = k
        block = psi_tensor[tuple(idx)]
        p_k0 += k0_diag[k] ** 2 * float(np.real(np.vdot(block, block)))

    # Prob(K1): p_reset * |<1|_site psi>|^2
    idx1 = [slice(None)] * N
    idx1[site] = 1
    block1 = psi_tensor[tuple(idx1)]
    p_k1 = p_reset * float(np.real(np.vdot(block1, block1)))

    # Prob(K2): p_reset * |<2|_site psi>|^2
    idx2 = [slice(None)] * N
    idx2[site] = 2
    block2 = psi_tensor[tuple(idx2)]
    p_k2 = p_reset * float(np.real(np.vdot(block2, block2)))

    total = p_k0 + p_k1 + p_k2
    if total < 1e-15:
        return psi

    # Sample which Kraus operator to apply
    r = rng.random() * total
    if r < p_k0:
        # Apply K0: multiply site amplitudes by k0_diag
        result_tensor = psi_tensor.copy()
        for k in range(d):
            idx = [slice(None)] * N
            idx[site] = k
            result_tensor[tuple(idx)] *= k0_diag[k]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k0)
    elif r < p_k0 + p_k1:
        # Apply K1 = sqrt(p)|0><1|: move site=1 amplitude to site=0
        result_tensor = np.zeros_like(psi_tensor)
        idx_src = [slice(None)] * N
        idx_src[site] = 1
        idx_dst = [slice(None)] * N
        idx_dst[site] = 0
        result_tensor[tuple(idx_dst)] = sqp * psi_tensor[tuple(idx_src)]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k1)
    else:
        # Apply K2 = sqrt(p)|0><2|: move site=2 amplitude to site=0
        result_tensor = np.zeros_like(psi_tensor)
        idx_src = [slice(None)] * N
        idx_src[site] = 2
        idx_dst = [slice(None)] * N
        idx_dst[site] = 0
        result_tensor[tuple(idx_dst)] = sqp * psi_tensor[tuple(idx_src)]
        result = result_tensor.reshape(dim)
        norm = np.sqrt(p_k2)

    if norm < 1e-15:
        msg = f"Thermal relaxation produced zero state at site {site}"
        raise PhysicsViolationError(msg)
    return result / norm


class QubitGKSLShotSimulator:
    """Shot-based qubit GKSL simulator using quantum trajectories.

    Each shot evolves a pure state through the Trotter decomposition,
    stochastically measuring Stinespring ancillas. At the end, the system
    is measured in the computational basis.

    The simulation runs in the 81-dim qutrit Hilbert space but corresponds
    to a qubit circuit with 2-qubit encoding per molecule.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QubitGKSLShotSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.N = params.N_molecules
        self.n_sys_qubits = 2 * self.N
        self.dim = params.d ** params.N_molecules  # 81

        # Build operators in qutrit space
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.n_ancilla = len(self.lindblad_ops)
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries."""
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines = [
            stinespring_unitary_from_lindblad(L_op, dt)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _prepare_initial_statevector(self, state_type: str) -> np.ndarray:
        """Prepare initial pure state in qutrit space."""
        d = self.params.d
        N = self.params.N_molecules
        dim = d ** N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
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

        return psi

    def _apply_stinespring_with_measurement(
        self, psi: np.ndarray, U_stine: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """Apply Stinespring unitary and stochastically measure ancilla."""
        d_sys = len(psi)
        psi_ext = np.zeros(2 * d_sys, dtype=np.complex128)
        psi_ext[:d_sys] = psi
        psi_ext = U_stine @ psi_ext

        psi_0 = psi_ext[:d_sys]
        psi_1 = psi_ext[d_sys:]
        p_0 = np.real(np.vdot(psi_0, psi_0))
        p_1 = np.real(np.vdot(psi_1, psi_1))

        if rng.random() < p_0 / (p_0 + p_1):
            return psi_0 / np.sqrt(p_0)
        else:
            return psi_1 / np.sqrt(p_1)

    def _trotter_step_trajectory(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """Apply one 2nd-order Trotter step to a pure state trajectory."""
        psi = self._U_H_half @ psi
        for U_stine in self._U_stines:
            psi = self._apply_stinespring_with_measurement(psi, U_stine, rng)
        psi = self._U_H_half @ psi
        return psi

    def _measure_system(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> int:
        """Measure system in computational basis."""
        probs = np.abs(psi) ** 2
        total = float(probs.sum())

        if abs(total - 1.0) > 1e-6:
            msg = f"State norm deviation before measurement: |sum(|psi|^2) - 1| = {abs(total - 1.0):.2e}"
            raise PhysicsViolationError(msg)

        probs = probs / total
        return int(rng.choice(len(probs), p=probs))

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        n_shots: int = 1000,
        seed: int | None = None,
    ) -> dict:
        """Run shot-based GKSL simulation using quantum trajectories.

        Parameters:
            t_max: total simulation time
            n_steps: number of Trotter steps
            initial_state: initial state type
            n_shots: number of quantum trajectories (shots)
            seed: random seed for reproducibility

        Returns:
            dict compatible with existing visualization functions.
        """
        start = time_module.time()
        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        psi_init = self._prepare_initial_statevector(initial_state)
        rng = np.random.default_rng(seed)

        n_time_points = n_steps + 1

        diag_accum = np.zeros((n_time_points, self.dim))
        rho_accum = np.zeros((n_time_points, self.dim, self.dim), dtype=np.complex128)
        final_counts: dict[int, int] = {}

        for _shot in range(n_shots):
            psi = psi_init.copy()
            diag_accum[0] += np.abs(psi) ** 2
            rho_accum[0] += np.outer(psi, psi.conj())

            for step in range(n_steps):
                psi = self._trotter_step_trajectory(psi, rng)
                diag_accum[step + 1] += np.abs(psi) ** 2
                rho_accum[step + 1] += np.outer(psi, psi.conj())

            outcome = self._measure_system(psi, rng)
            final_counts[outcome] = final_counts.get(outcome, 0) + 1

        diag_accum /= n_shots
        rho_accum /= n_shots

        times = [s * dt for s in range(n_time_points)]
        populations = [
            _populations_from_diagonal(diag_accum[s], self.params)
            for s in range(n_time_points)
        ]
        entropies = [
            compute_von_neumann_entropy(rho_accum[s]) for s in range(n_time_points)
        ]
        purities = [compute_purity(rho_accum[s]) for s in range(n_time_points)]
        traces = [float(np.real(np.trace(rho_accum[s]))) for s in range(n_time_points)]

        elapsed = time_module.time() - start

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_accum[-1],
            "elapsed_time": elapsed,
            "method": "qubit_gksl_shot",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_sys_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "n_shots": n_shots,
            "counts": final_counts,
            "seed": seed,
        }


class QubitGKSLNoisyShotSimulator(QubitGKSLShotSimulator):
    """Shot-based qubit GKSL simulator with stochastic hardware noise.

    Adds per-gate stochastic depolarization and thermal relaxation noise
    to each trajectory, matching the noise model of QubitGKSLNoisySimulator.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
        p_depol: depolarization probability per 2-qubit gate (default 0.01)
        T1: energy relaxation time (default None = no relaxation)
        t_gate: 2-qubit gate time (default 300.0)
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.01,
        T1: float | None = None,
        t_gate: float = 300.0,
    ) -> None:
        super().__init__(params)
        if p_depol < 0.0 or p_depol > 1.0:
            msg = f"p_depol must be in [0, 1], got {p_depol}"
            raise ValueError(msg)
        self.p_depol = p_depol
        self.T1 = T1
        self.t_gate = t_gate

        if T1 is not None and T1 > 0:
            self.p_reset = 1.0 - np.exp(-t_gate / T1)
        else:
            self.p_reset = 0.0

        self._lindblad_sites = self._compute_lindblad_sites()

    def _compute_lindblad_sites(self) -> list[list[int]]:
        """Determine molecule indices for each Lindblad operator."""
        sites: list[list[int]] = []
        for i, j in self.params.neighbors:
            sites.append([i, j])
            sites.append([i, j])
        for _channel in range(5):
            for mol in range(self.params.N_molecules):
                sites.append([mol])
        assert len(sites) == len(self.lindblad_ops)
        return sites

    def _trotter_step_trajectory(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """2nd-order Trotter step with stochastic noise."""
        d = self.params.d
        N = self.params.N_molecules

        # --- Half Hamiltonian + transfer gate noise ---
        psi = self._U_H_half @ psi
        for i, j in self.params.neighbors:
            psi = _apply_stochastic_depolarization_pair(
                psi, i, j, d, N, self.p_depol, rng
            )

        # --- All Lindblad channels + per-channel noise ---
        for k, U_stine in enumerate(self._U_stines):
            psi = self._apply_stinespring_with_measurement(psi, U_stine, rng)
            sites = self._lindblad_sites[k]
            if len(sites) == 1:
                psi = _apply_stochastic_depolarization_single(
                    psi, sites[0], d, N, self.p_depol, rng
                )
            else:
                psi = _apply_stochastic_depolarization_pair(
                    psi, sites[0], sites[1], d, N, self.p_depol, rng
                )

        # --- Half Hamiltonian + transfer gate noise ---
        psi = self._U_H_half @ psi
        for i, j in self.params.neighbors:
            psi = _apply_stochastic_depolarization_pair(
                psi, i, j, d, N, self.p_depol, rng
            )

        # --- Thermal relaxation on all molecules ---
        if self.p_reset > 0.0:
            for mol in range(N):
                psi = _apply_stochastic_thermal_relaxation_single(
                    psi, mol, d, N, self.p_reset, rng
                )

        return psi

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        n_shots: int = 1000,
        seed: int | None = None,
    ) -> dict:
        """Run noisy shot-based simulation."""
        result = super().simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state,
            n_shots=n_shots, seed=seed,
        )
        result["method"] = "qubit_gksl_noisy_shot"
        result["noise_params"] = {
            "p_depol": self.p_depol,
            "T1": self.T1,
            "t_gate": self.t_gate,
            "p_reset": self.p_reset,
        }
        return result
