"""Qudit GKSL shot-based simulator using quantum trajectories.

Implements shot-based (quantum trajectory) simulation for GKSL-Lindblad dynamics
using the Stinespring dilation approach. Instead of evolving the full density matrix,
each shot evolves a pure state through the Trotter steps, stochastically measuring
ancilla qudits at each Stinespring channel, then measures the system at the end.

Since this simulator targets qudit-type quantum computers (e.g. qudit-boson
ion-trap processors), all registers — including Stinespring ancillas — are
native d-level qudits, not 2-level qubits.

Two classes are provided:
  - QuditGKSLShotSimulator: ideal shot-based (no hardware noise)
  - QuditGKSLNoisyShotSimulator: shot-based with per-gate stochastic noise

The quantum trajectory method is mathematically equivalent to the density matrix
approach (QuditGKSLSimulator) in the limit of infinite shots:
  rho = lim_{N->inf} (1/N) sum_k |psi_k><psi_k|

Both use the same palindromic half-step Trotter structure for Lindblad channels.
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
from stinespring_utils import stinespring_unitary_from_lindblad


def _populations_from_diagonal(diag: np.ndarray, params: GKSLPhysicalParameters) -> dict:
    """Compute N_S0, N_T1, N_S1 populations from probability vector (diagonal of rho)."""
    N = params.N_molecules
    d = params.d
    dim = d ** N

    N_S0 = 0.0
    N_T1 = 0.0
    N_S1 = 0.0
    per_molecule: list[dict[int, float]] = [{0: 0.0, 1: 0.0, 2: 0.0} for _ in range(N)]

    for idx in range(dim):
        p = diag[idx]
        remainder = idx
        for mol in range(N - 1, -1, -1):
            local_state = remainder % d
            remainder //= d
            if local_state == 0:
                N_S0 += p
            elif local_state == 1:
                N_T1 += p
            else:
                N_S1 += p
            per_molecule[mol][local_state] += p

    return {
        "N_S0": float(N_S0),
        "N_T1": float(N_T1),
        "N_S1": float(N_S1),
        "per_molecule_populations": per_molecule,
    }


class QuditGKSLShotSimulator:
    """Shot-based qudit GKSL simulator using quantum trajectories.

    Each shot evolves a pure state |psi> through the Trotter decomposition,
    stochastically measuring Stinespring ancilla qudits via the Born rule.
    At the end, the system is measured in the computational basis.

    Since this targets qudit quantum computers, Stinespring ancillas are
    d-level qudits (same dimension as system qudits), not 2-level qubits.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QuditGKSLShotSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.dim = params.d ** params.N_molecules
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension

        # Build operators in native qutrit space
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.n_system_qudits = params.N_molecules
        self.n_ancilla_qudits = len(self.lindblad_ops)

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (shared across all shots).

        Stores half-dt Stinespring unitaries for the palindromic ordering in
        _trotter_step_trajectory.  Ancilla dimension is d_anc (= params.d).
        """
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines_half = [
            stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _prepare_initial_statevector(self, state_type: str) -> np.ndarray:
        """Prepare initial pure state |psi> in qutrit space."""
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
        """Apply Stinespring unitary and stochastically measure ancilla qudit.

        Extends |psi> with ancilla |0>, applies U_stine, then measures
        the ancilla in the {|0>, |1>, ..., |d_anc-1>} basis using the Born rule.

        Returns the post-measurement system state (normalized).
        """
        d_sys = len(psi)
        d_anc = self.d_anc
        # |Psi_ext> = |0>_anc ⊗ |psi>_sys
        psi_ext = np.zeros(d_anc * d_sys, dtype=np.complex128)
        psi_ext[:d_sys] = psi

        # Apply Stinespring unitary
        psi_ext = U_stine @ psi_ext

        # Project onto ancilla outcomes and compute Born probabilities
        probs = np.zeros(d_anc)
        psi_branches = []
        for k in range(d_anc):
            branch = psi_ext[k * d_sys : (k + 1) * d_sys]
            psi_branches.append(branch)
            probs[k] = np.real(np.vdot(branch, branch))

        # Born rule: sample measurement outcome
        total = probs.sum()
        probs_norm = probs / total
        outcome = rng.choice(d_anc, p=probs_norm)
        return psi_branches[outcome] / np.sqrt(probs[outcome])

    def _trotter_step_trajectory(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad ordering for trajectory.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        Matches QuditGKSLSimulator._trotter_step structure so that the
        infinite-shot limit equals the density matrix result exactly.
        """
        # Half Hamiltonian
        psi = self._U_H_half @ psi
        # Forward half-step for all Lindblad channels
        for U_stine_half in self._U_stines_half:
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
        # Reverse half-step for all Lindblad channels (palindromic)
        for U_stine_half in reversed(self._U_stines_half):
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
        # Half Hamiltonian
        psi = self._U_H_half @ psi
        return psi

    def _measure_system(
        self, psi: np.ndarray, rng: np.random.Generator
    ) -> int:
        """Measure system in computational basis, returning basis state index.

        Raises PhysicsViolationError if state norm deviates significantly from 1.
        """
        probs = np.abs(psi) ** 2
        total = float(probs.sum())

        if abs(total - 1.0) > 1e-6:
            msg = f"State norm deviation before measurement: |sum(|psi|^2) - 1| = {abs(total - 1.0):.2e}"
            raise PhysicsViolationError(msg)

        # Correct floating-point rounding for np.random.choice
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
            dict with time series of populations, entropy, purity, trace,
            reconstructed final density matrix, measurement counts, etc.
        """
        start = time_module.time()
        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        psi_init = self._prepare_initial_statevector(initial_state)
        rng = np.random.default_rng(seed)

        n_time_points = n_steps + 1

        # Accumulators: diagonal probabilities at each time step
        diag_accum = np.zeros((n_time_points, self.dim))
        # Density matrix accumulator at each time step (for entropy/purity)
        rho_accum = np.zeros((n_time_points, self.dim, self.dim), dtype=np.complex128)
        # Final measurement counts
        final_counts: dict[int, int] = {}

        for _shot in range(n_shots):
            psi = psi_init.copy()

            # Record initial state
            diag_accum[0] += np.abs(psi) ** 2
            rho_accum[0] += np.outer(psi, psi.conj())

            for step in range(n_steps):
                psi = self._trotter_step_trajectory(psi, rng)
                diag_accum[step + 1] += np.abs(psi) ** 2
                rho_accum[step + 1] += np.outer(psi, psi.conj())

            # Final measurement
            outcome = self._measure_system(psi, rng)
            final_counts[outcome] = final_counts.get(outcome, 0) + 1

        # Average over shots
        diag_accum /= n_shots
        rho_accum /= n_shots

        # Compute time series
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

        # N VirtRz (H_0 diagonal) + (N-1) CustomTwo (NN transfer) + n_lindblad×2 Stinespring (palindromic)
        gates_per_step = self.n_system_qudits + len(self.params.neighbors) + self.n_ancilla_qudits * 2

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_accum[-1],
            "elapsed_time": elapsed,
            "method": "qudit_gksl_shot",
            "params": self.params.to_dict(),
            "n_system_qudits": self.n_system_qudits,
            "n_ancilla_qudits": self.n_ancilla_qudits,
            "d_anc": self.d_anc,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "n_shots": n_shots,
            "counts": final_counts,
            "seed": seed,
        }


# ---------------------------------------------------------------------------
# Stochastic noise operators for quantum trajectories
# ---------------------------------------------------------------------------


def _apply_weyl_heisenberg_single_site(
    psi: np.ndarray, site: int, a: int, b: int, d: int, N: int
) -> np.ndarray:
    """Apply the Weyl-Heisenberg operator X^a Z^b on a single molecule site.

    X^a Z^b |k> = omega^(b*k) |(k+a) mod d>
    where omega = exp(2*pi*i/d).

    Operates on the multi-site state vector by reshaping to tensor form.
    """
    omega = np.exp(2j * np.pi / d)
    dim = d ** N
    psi_tensor = psi.reshape([d] * N)
    result = np.zeros_like(psi_tensor)

    for k in range(d):
        k_new = (k + a) % d
        phase = omega ** (b * k)
        # Select slice where site has value k
        idx_src = [slice(None)] * N
        idx_src[site] = k
        idx_dst = [slice(None)] * N
        idx_dst[site] = k_new
        result[tuple(idx_dst)] += phase * psi_tensor[tuple(idx_src)]

    return result.reshape(dim)


def _apply_stochastic_depolarization_single(
    psi: np.ndarray, site: int, d: int, N: int, p: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stochastically apply single-site depolarization noise.

    Kraus decomposition of E_S(rho) = (1-p)rho + (p/d) I_S x Tr_S[rho]:
      - With prob (1 - p + p/d^2): identity (no error)
      - With prob p/d^2 each: apply X^a Z^b for (a,b) != (0,0)
    """
    if p <= 0.0:
        return psi

    p_identity = 1.0 - p + p / (d * d)
    r = rng.random()

    if r < p_identity:
        return psi

    # Select random non-identity Weyl-Heisenberg operator X^a Z^b.
    # error_idx in [0, d^2-2] maps to ab in [1, d^2-1] (skipping identity at ab=0).
    # ab encodes (a, b) as ab = a*d + b, giving d^2-1 non-identity operators.
    error_idx = rng.integers(0, d * d - 1)
    ab = error_idx + 1
    a = ab // d
    b = ab % d
    return _apply_weyl_heisenberg_single_site(psi, site, a, b, d, N)


def _apply_stochastic_depolarization_pair(
    psi: np.ndarray, site_a: int, site_b: int, d: int, N: int, p: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stochastically apply pair depolarization noise.

    Kraus decomposition of E_{AB}(rho) = (1-p)rho + (p/d^2) I_{AB} x Tr_{AB}[rho]:
      - With prob (1 - p + p/d^4): identity
      - With prob p/d^4 each: apply X_A^a Z_A^b x X_B^c Z_B^d for (a,b,c,d) != (0,0,0,0)
    """
    if p <= 0.0:
        return psi

    d4 = d ** 4
    p_identity = 1.0 - p + p / d4
    r = rng.random()

    if r < p_identity:
        return psi

    # Select random non-identity pair operator
    error_idx = rng.integers(0, d4 - 1)
    abcd = error_idx + 1  # skip (0,0,0,0)
    a_A = abcd // (d * d * d) % d
    b_A = abcd // (d * d) % d
    a_B = abcd // d % d
    b_B = abcd % d

    psi = _apply_weyl_heisenberg_single_site(psi, site_a, a_A, b_A, d, N)
    psi = _apply_weyl_heisenberg_single_site(psi, site_b, a_B, b_B, d, N)
    return psi


def _apply_stochastic_dephasing_single(
    psi: np.ndarray, site: int, d: int, N: int, p: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Stochastically apply single-site dephasing noise.

    E_deph(rho) = (1-p) rho + p sum_k P_k rho P_k
    where P_k = |k><k| on the site.

    Stochastic: with prob (1-p) keep state, with prob p project onto
    random local basis state (chosen with Born rule probabilities).
    """
    if p <= 0.0:
        return psi

    if rng.random() >= p:
        return psi

    # Compute local probabilities at the site
    dim = d ** N
    psi_tensor = psi.reshape([d] * N)
    local_probs = np.zeros(d)
    for k in range(d):
        idx = [slice(None)] * N
        idx[site] = k
        local_probs[k] = np.real(np.vdot(
            psi_tensor[tuple(idx)], psi_tensor[tuple(idx)]
        ))

    total = local_probs.sum()
    if total < 1e-15:
        return psi
    local_probs /= total

    # Sample which local state to project onto
    k_chosen = rng.choice(d, p=local_probs)

    # Apply projection P_k: zero out all components where site != k
    result_tensor = np.zeros_like(psi_tensor)
    idx = [slice(None)] * N
    idx[site] = k_chosen
    result_tensor[tuple(idx)] = psi_tensor[tuple(idx)]

    # Renormalize
    result = result_tensor.reshape(dim)
    norm = np.sqrt(np.real(np.vdot(result, result)))
    if norm < 1e-15:
        msg = f"Dephasing projection produced zero state at site {site}"
        raise PhysicsViolationError(msg)
    return result / norm


class QuditGKSLNoisyShotSimulator(QuditGKSLShotSimulator):
    """Shot-based qudit GKSL simulator with stochastic hardware noise.

    Extends QuditGKSLShotSimulator by applying per-gate stochastic noise
    (depolarization + dephasing) after each gate operation in each trajectory.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
        p_depol: depolarization probability per 2-qudit gate (default 0.01)
        p_dephasing: dephasing probability per 2-qudit gate (default 0.0)
        depol_pair_only: if True, apply noise only after pair (2+ qudit)
            interactions, skipping single-site Lindblad channels (default False)
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.01,
        p_dephasing: float = 0.0,
        depol_pair_only: bool = False,
    ) -> None:
        super().__init__(params)
        if p_depol < 0.0 or p_depol > 1.0:
            msg = f"p_depol must be in [0, 1], got {p_depol}"
            raise ValueError(msg)
        if p_dephasing < 0.0 or p_dephasing > 1.0:
            msg = f"p_dephasing must be in [0, 1], got {p_dephasing}"
            raise ValueError(msg)
        self.p_depol = p_depol
        self.p_dephasing = p_dephasing
        self.depol_pair_only = depol_pair_only
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
        """Symmetric Trotter step with palindromic Lindblad ordering and noise."""
        d = self.params.d
        N = self.params.N_molecules

        # --- Half Hamiltonian + transfer gate noise ---
        psi = self._U_H_half @ psi
        for i, j in self.params.neighbors:
            psi = _apply_stochastic_depolarization_pair(
                psi, i, j, d, N, self.p_depol, rng
            )
            if self.p_dephasing > 0.0:
                psi = _apply_stochastic_dephasing_single(
                    psi, i, d, N, self.p_dephasing, rng
                )
                psi = _apply_stochastic_dephasing_single(
                    psi, j, d, N, self.p_dephasing, rng
                )

        # --- Forward half-step Lindblad channels + per-channel noise ---
        for k, U_stine_half in enumerate(self._U_stines_half):
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
            sites = self._lindblad_sites[k]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    psi = _apply_stochastic_depolarization_single(
                        psi, sites[0], d, N, self.p_depol, rng
                    )
                    if self.p_dephasing > 0.0:
                        psi = _apply_stochastic_dephasing_single(
                            psi, sites[0], d, N, self.p_dephasing, rng
                        )
            else:
                psi = _apply_stochastic_depolarization_pair(
                    psi, sites[0], sites[1], d, N, self.p_depol, rng
                )
                if self.p_dephasing > 0.0:
                    psi = _apply_stochastic_dephasing_single(
                        psi, sites[0], d, N, self.p_dephasing, rng
                    )
                    psi = _apply_stochastic_dephasing_single(
                        psi, sites[1], d, N, self.p_dephasing, rng
                    )

        # --- Reverse half-step Lindblad channels + per-channel noise (palindromic) ---
        n_channels = len(self._U_stines_half)
        for k_rev in range(n_channels - 1, -1, -1):
            U_stine_half = self._U_stines_half[k_rev]
            psi = self._apply_stinespring_with_measurement(psi, U_stine_half, rng)
            sites = self._lindblad_sites[k_rev]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    psi = _apply_stochastic_depolarization_single(
                        psi, sites[0], d, N, self.p_depol, rng
                    )
                    if self.p_dephasing > 0.0:
                        psi = _apply_stochastic_dephasing_single(
                            psi, sites[0], d, N, self.p_dephasing, rng
                        )
            else:
                psi = _apply_stochastic_depolarization_pair(
                    psi, sites[0], sites[1], d, N, self.p_depol, rng
                )
                if self.p_dephasing > 0.0:
                    psi = _apply_stochastic_dephasing_single(
                        psi, sites[0], d, N, self.p_dephasing, rng
                    )
                    psi = _apply_stochastic_dephasing_single(
                        psi, sites[1], d, N, self.p_dephasing, rng
                    )

        # --- Half Hamiltonian + transfer gate noise ---
        psi = self._U_H_half @ psi
        for i, j in self.params.neighbors:
            psi = _apply_stochastic_depolarization_pair(
                psi, i, j, d, N, self.p_depol, rng
            )
            if self.p_dephasing > 0.0:
                psi = _apply_stochastic_dephasing_single(
                    psi, i, d, N, self.p_dephasing, rng
                )
                psi = _apply_stochastic_dephasing_single(
                    psi, j, d, N, self.p_dephasing, rng
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
        """Run noisy shot-based simulation.

        Returns the same dict structure as QuditGKSLShotSimulator.simulate()
        plus noise_params.
        """
        result = super().simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state,
            n_shots=n_shots, seed=seed,
        )
        result["method"] = "qudit_gksl_noisy_shot"
        result["noise_params"] = {
            "p_depol": self.p_depol,
            "p_dephasing": self.p_dephasing,
            "depol_pair_only": self.depol_pair_only,
        }
        return result
