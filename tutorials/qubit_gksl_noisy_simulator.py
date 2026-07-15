"""Qubit GKSL noisy simulator with hardware noise model.

Scenario 8: Qubit-based GKSL-Lindblad with hardware depolarization + thermal relaxation.

Extends QubitGKSLSimulator by applying per-gate local noise channels after
each gate operation in the Trotter step. The simulation runs in the 4^N = 256
dimensional qubit-encoded space, using d=4 depolarization per molecule (equivalent
to d=2 Pauli depolarization on each physical qubit pair).

Physical distinction:
  - Lindblad dissipators (TTA, fluorescence, etc.) model real physical processes
  - Hardware noise (depolarization, thermal relaxation) models gate imperfections
  These are conceptually independent and both included in this simulator.

Noise channels use d=4 Pauli operators on the 2-qubit encoding of each molecule.
This correctly models forbidden-state (|11>) leakage from qubit hardware noise,
which is absent in the native qutrit (d=3) encoding.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_simulator import QubitGKSLSimulator
from stinespring_utils import apply_stinespring_to_density_matrix


def _apply_local_depolarization_single_qubit(
    rho: np.ndarray, site: int, N: int, p: float
) -> np.ndarray:
    """Apply d=4 depolarization on a single molecule in 4^N qubit space.

    E_site[rho] = (1-p) rho + p/d * I_site ⊗ Tr_site[rho]

    where d=4 is the local 2-qubit dimension per molecule.
    This can cause forbidden-state leakage via the I_site term.
    """
    if p <= 0.0:
        return rho
    d = 4  # local dimension per molecule (2 qubits)
    dim = d**N
    rho_tensor = rho.reshape([d] * (2 * N))

    # Partial trace over 'site': contract ket[site] with bra[site]
    rho_rest = np.trace(rho_tensor, axis1=site, axis2=N + site)

    # Construct I_site/d ⊗ rho_rest
    mixed = np.zeros_like(rho_tensor)
    for m in range(d):
        idx = [slice(None)] * (2 * N)
        idx[site] = m
        idx[N + site] = m
        mixed[tuple(idx)] = rho_rest / d

    return ((1 - p) * rho_tensor + p * mixed).reshape(dim, dim)


def _apply_local_depolarization_pair_qubit(
    rho: np.ndarray, site_a: int, site_b: int, N: int, p: float
) -> np.ndarray:
    """Apply d=16 depolarization on a molecule pair in 4^N qubit space.

    E_{a,b}[rho] = (1-p) rho + p/d^2 * I_{a,b} ⊗ Tr_{a,b}[rho]

    where d=4 and d^2=16 is the joint pair dimension.
    """
    if p <= 0.0:
        return rho
    d = 4
    dim = d**N
    d_pair = d * d  # 16
    rho_tensor = rho.reshape([d] * (2 * N))

    # Partial trace over sites a and b
    rho_rest_sum = None
    for ma in range(d):
        for mb in range(d):
            idx = [slice(None)] * (2 * N)
            idx[site_a] = ma
            idx[N + site_a] = ma
            idx[site_b] = mb
            idx[N + site_b] = mb
            block = rho_tensor[tuple(idx)]
            if rho_rest_sum is None:
                rho_rest_sum = block.copy()
            else:
                rho_rest_sum = rho_rest_sum + block

    # Construct I_{a,b}/d^2 ⊗ rho_rest
    mixed = np.zeros_like(rho_tensor)
    for ka in range(d):
        for kb in range(d):
            idx = [slice(None)] * (2 * N)
            idx[site_a] = ka
            idx[N + site_a] = ka
            idx[site_b] = kb
            idx[N + site_b] = kb
            mixed[tuple(idx)] = rho_rest_sum / d_pair

    return ((1 - p) * rho_tensor + p * mixed).reshape(dim, dim)


def _apply_thermal_relaxation_qubit(
    rho: np.ndarray, site: int, N: int, p_reset: float
) -> np.ndarray:
    """Apply 2-qubit thermal relaxation on a molecule in 4^N qubit space.

    Models independent T1 relaxation on each physical qubit of the molecule's
    2-qubit encoding. Each qubit independently decays to |0> with probability
    p_reset.

    4x4 Kraus operators for the 2-qubit encoding |b1 b0>:
      K_00 = diag(1, sqrt(1-p), sqrt(1-p), (1-p))  (no decay)
      K_01 = sqrt(p) * [[0,1,0,0],[0,0,0,0],[0,0,0,sqrt(1-p)],[0,0,0,0]]
      K_10 = sqrt(p) * [[0,0,1,0],[0,0,0,sqrt(1-p)],[0,0,0,0],[0,0,0,0]]
      K_11 = p * |00><11|

    These satisfy Σ K†K = I and model independent qubit relaxation.
    """
    if p_reset <= 0.0:
        return rho
    d = 4
    dim = d**N
    sq = np.sqrt(1.0 - p_reset)
    sqp = np.sqrt(p_reset)

    # Build 4x4 Kraus matrices
    K00 = np.diag(np.array([1.0, sq, sq, sq * sq], dtype=np.complex128))
    K01 = sqp * np.array(
        [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, sq], [0, 0, 0, 0]],
        dtype=np.complex128,
    )
    K10 = sqp * np.array(
        [[0, 0, 1, 0], [0, 0, 0, sq], [0, 0, 0, 0], [0, 0, 0, 0]],
        dtype=np.complex128,
    )
    K11 = p_reset * np.array(
        [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
        dtype=np.complex128,
    )

    rho_tensor = rho.reshape([d] * (2 * N))
    result = np.zeros_like(rho_tensor)

    for K in [K00, K01, K10, K11]:
        # Apply K on ket index (site), K† on bra index (N + site)
        # temp = K_site ρ K†_site
        temp = rho_tensor.copy()

        # Apply K on ket index
        new_temp = np.zeros_like(temp)
        for a in range(d):
            for c in range(d):
                if abs(K[a, c]) < 1e-15:
                    continue
                idx_in = [slice(None)] * (2 * N)
                idx_in[site] = c
                idx_out = [slice(None)] * (2 * N)
                idx_out[site] = a
                new_temp[tuple(idx_out)] += K[a, c] * temp[tuple(idx_in)]
        temp = new_temp

        # Apply K† on bra index
        new_temp = np.zeros_like(temp)
        for b in range(d):
            for e in range(d):
                kbe_conj = np.conj(K[b, e])
                if abs(kbe_conj) < 1e-15:
                    continue
                idx_in = [slice(None)] * (2 * N)
                idx_in[N + site] = e
                idx_out = [slice(None)] * (2 * N)
                idx_out[N + site] = b
                new_temp[tuple(idx_out)] += kbe_conj * temp[tuple(idx_in)]

        result += new_temp

    return result.reshape(dim, dim)


def _apply_local_dephasing_single_qubit(
    rho: np.ndarray, site: int, N: int, p: float
) -> np.ndarray:
    """Apply dephasing channel on a single molecule in 4^N qubit space.

    E_deph[rho] = (1-p) rho + p * sum_k (|k><k|_site ⊗ I_rest) rho (|k><k|_site ⊗ I_rest)

    This decoheres the off-diagonal elements between different local states
    at the specified site while leaving diagonal elements and other sites intact.
    """
    if p <= 0.0:
        return rho
    d = 4
    dim = d**N
    rho_tensor = rho.reshape([d] * (2 * N))

    dephased = np.zeros_like(rho_tensor)
    for k in range(d):
        idx = [slice(None)] * (2 * N)
        idx[site] = k
        idx[N + site] = k
        dephased[tuple(idx)] = rho_tensor[tuple(idx)]

    return ((1 - p) * rho_tensor + p * dephased).reshape(dim, dim)


class QubitGKSLNoisySimulator(QubitGKSLSimulator):
    """Qubit GKSL simulator with d=4 Pauli depolarization, dephasing and thermal relaxation.

    Inherits the Stinespring + 2nd-order Trotter approach from QubitGKSLSimulator
    (running in 4^N = 256 dim qubit space) and adds per-gate local noise channels
    using d=4 two-qubit Pauli operators to model hardware imperfections.

    The noise model correctly accounts for forbidden-state leakage: d=4 Pauli
    errors can drive the state into the |11> forbidden subspace, which is a
    fundamental property of qubit encoding that is absent in native qutrit encoding.

    Noise placement matches ``QubitGKSLNoisyShotSimulator``.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
        p_depol: depolarization probability per 2-qubit gate (default 0.01 = 1%)
        p_dephasing: dephasing probability per gate (default 0.0)
        T1: energy relaxation time in natural units (default None = no relaxation)
        T2: dephasing time in natural units (default None = no relaxation)
        t_gate: 2-qubit gate time in natural units (default 300.0 fs)
        depol_pair_only: if True, apply noise only after pair (2+ qubit)
            interactions, skipping single-site Lindblad channels (default False)
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.01,
        p_dephasing: float = 0.0,
        T1: float | None = None,
        T2: float | None = None,
        t_gate: float = 300.0,
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
        self.T1 = T1
        self.T2 = T2
        self.t_gate = t_gate
        self.depol_pair_only = depol_pair_only

        # Compute thermal relaxation probability per gate
        if T1 is not None and T1 > 0:
            self.p_reset = 1.0 - np.exp(-t_gate / T1)
        else:
            self.p_reset = 0.0

        # Pre-compute Lindblad operator site mapping
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

    def check_forbidden_states(
        self, rho_qubit: np.ndarray, step: int | None = None
    ) -> float:
        """Track leakage without raising.

        For the noisy simulator, forbidden-state leakage is an expected physical
        effect of qubit encoding with hardware noise, not a validation failure.
        """
        from qubit_gksl_simulator import compute_forbidden_state_population

        return float(compute_forbidden_state_population(rho_qubit, self._mapping))

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad ordering and d=4 qubit hardware noise.

        Gate noise is applied after each gate operation using d=4 Pauli
        depolarization in the 256-dim qubit space:
          1. Half Hamiltonian + NN pair depolarization + dephasing
          2. Forward half-step Lindblad channels + per-channel depolarization + dephasing
          3. Reverse half-step Lindblad channels + per-channel depolarization + dephasing (palindromic)
          4. Half Hamiltonian + NN pair depolarization + dephasing
          5. Thermal relaxation on all molecules

        Uses precomputed half-step unitaries from ``_precompute_unitaries``.
        Noise placement matches ``QubitGKSLNoisyShotSimulator``.
        """
        N = self.params.N_molecules

        # --- Half Hamiltonian ---
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        for i, j in self.params.neighbors:
            rho = _apply_local_depolarization_pair_qubit(
                rho, i, j, N, self.p_depol
            )
            if self.p_dephasing > 0.0:
                rho = _apply_local_dephasing_single_qubit(rho, i, N, self.p_dephasing)
                rho = _apply_local_dephasing_single_qubit(rho, j, N, self.p_dephasing)

        # --- Forward half-step Lindblad channels ---
        for k, U_stine_half in enumerate(self._U_stines_half):
            rho = apply_stinespring_to_density_matrix(rho, U_stine_half)
            sites = self._lindblad_sites[k]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    rho = _apply_local_depolarization_single_qubit(
                        rho, sites[0], N, self.p_depol
                    )
                    if self.p_dephasing > 0.0:
                        rho = _apply_local_dephasing_single_qubit(
                            rho, sites[0], N, self.p_dephasing
                        )
            else:
                rho = _apply_local_depolarization_pair_qubit(
                    rho, sites[0], sites[1], N, self.p_depol
                )
                if self.p_dephasing > 0.0:
                    rho = _apply_local_dephasing_single_qubit(
                        rho, sites[0], N, self.p_dephasing
                    )
                    rho = _apply_local_dephasing_single_qubit(
                        rho, sites[1], N, self.p_dephasing
                    )

        # --- Reverse half-step Lindblad channels (palindromic) ---
        n_channels = len(self._U_stines_half)
        for k_rev in range(n_channels - 1, -1, -1):
            U_stine_half = self._U_stines_half[k_rev]
            rho = apply_stinespring_to_density_matrix(rho, U_stine_half)
            sites = self._lindblad_sites[k_rev]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    rho = _apply_local_depolarization_single_qubit(
                        rho, sites[0], N, self.p_depol
                    )
                    if self.p_dephasing > 0.0:
                        rho = _apply_local_dephasing_single_qubit(
                            rho, sites[0], N, self.p_dephasing
                        )
            else:
                rho = _apply_local_depolarization_pair_qubit(
                    rho, sites[0], sites[1], N, self.p_depol
                )
                if self.p_dephasing > 0.0:
                    rho = _apply_local_dephasing_single_qubit(
                        rho, sites[0], N, self.p_dephasing
                    )
                    rho = _apply_local_dephasing_single_qubit(
                        rho, sites[1], N, self.p_dephasing
                    )

        # --- Half Hamiltonian ---
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        for i, j in self.params.neighbors:
            rho = _apply_local_depolarization_pair_qubit(
                rho, i, j, N, self.p_depol
            )
            if self.p_dephasing > 0.0:
                rho = _apply_local_dephasing_single_qubit(rho, i, N, self.p_dephasing)
                rho = _apply_local_dephasing_single_qubit(rho, j, N, self.p_dephasing)

        # --- Thermal relaxation (all molecules, accumulated gate time) ---
        if self.p_reset > 0.0:
            for mol in range(N):
                rho = _apply_thermal_relaxation_qubit(
                    rho, mol, N, self.p_reset
                )

        return rho

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run noisy Qubit GKSL simulation.

        Returns the same dict structure as QubitGKSLSimulator.simulate()
        plus noise_params.
        """
        result = super().simulate(t_max=t_max, n_steps=n_steps, initial_state=initial_state)
        result["method"] = "qubit_gksl_noisy"
        result["noise_params"] = {
            "p_depol": self.p_depol,
            "p_dephasing": self.p_dephasing,
            "T1": self.T1,
            "T2": self.T2,
            "t_gate": self.t_gate,
            "p_reset": self.p_reset,
            "depol_pair_only": self.depol_pair_only,
        }
        return result
