"""Qudit GKSL noisy simulator with hardware noise model.

Scenario 7: Qudit-based GKSL-Lindblad with hardware depolarization + dephasing.

Extends QuditGKSLSimulator by applying per-gate local noise channels after
each 2-qudit gate operation in the Trotter step. The noise acts on the
subsystem(s) involved in each gate, not globally.

Physical distinction:
  - Lindblad dissipators (TTA, fluorescence, etc.) model real physical processes
  - Hardware noise (depolarization, dephasing) models gate imperfections
  These are conceptually independent and both included in this simulator.

Noise channels (2-qudit gates only; 1-qudit gates are ideal):
  - Depolarization: E_S[rho] = (1-p) rho + p/d_S * I_S ⊗ Tr_S[rho]
  - Dephasing: E_S[rho] = (1-p) rho + p * sum_k (|k><k|_S ⊗ I_rest) rho (|k><k|_S ⊗ I_rest)

Gate-noise mapping for the Trotter step:
  - Hamiltonian transfer: one local pair depolarization per nearest-neighbor pair
  - Stinespring channel: one local depolarization on the subsystem the Lindblad
    operator acts on (single molecule for site-local operators, molecule pair for TTA)
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_simulator import QuditGKSLSimulator
from stinespring_utils import apply_stinespring_to_density_matrix


def _apply_local_depolarization_single(
    rho: np.ndarray, site: int, d: int, N: int, p: float
) -> np.ndarray:
    """Apply depolarization channel on a single molecule site.

    E_site[rho] = (1-p) rho + p/d * I_site ⊗ Tr_site[rho]

    Parameters:
        rho: d^N × d^N density matrix
        site: molecule index (0-based)
        d: local dimension per molecule (3 for qutrit)
        N: total number of molecules
        p: depolarization probability in [0, 1]
    """
    if p <= 0.0:
        return rho
    dim = d**N
    rho_tensor = rho.reshape([d] * (2 * N))

    # Partial trace over 'site': contract ket[site] with bra[site]
    # np.trace sums over the diagonal of the two specified axes
    rho_rest = np.trace(rho_tensor, axis1=site, axis2=N + site)

    # Construct I_site/d ⊗ rho_rest
    mixed = np.zeros_like(rho_tensor)
    for m in range(d):
        idx = [slice(None)] * (2 * N)
        idx[site] = m
        idx[N + site] = m
        mixed[tuple(idx)] = rho_rest / d

    return ((1 - p) * rho_tensor + p * mixed).reshape(dim, dim)


def _apply_local_depolarization_pair(
    rho: np.ndarray, site_a: int, site_b: int, d: int, N: int, p: float
) -> np.ndarray:
    """Apply depolarization channel on a molecule pair.

    E_{a,b}[rho] = (1-p) rho + p/d^2 * I_{a,b} ⊗ Tr_{a,b}[rho]

    Parameters:
        rho: d^N × d^N density matrix
        site_a, site_b: molecule indices (0-based, site_a != site_b)
        d: local dimension per molecule
        N: total number of molecules
        p: depolarization probability
    """
    if p <= 0.0:
        return rho
    dim = d**N
    d_pair = d * d
    rho_tensor = rho.reshape([d] * (2 * N))

    # Partial trace over sites a and b: sum over diagonal in both
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


def _apply_local_dephasing_single(
    rho: np.ndarray, site: int, d: int, N: int, p: float
) -> np.ndarray:
    """Apply dephasing channel on a single molecule site.

    E_deph[rho] = (1-p) rho + p * sum_k (|k><k|_site ⊗ I_rest) rho (|k><k|_site ⊗ I_rest)

    This decoheres the off-diagonal elements between different local states
    at the specified site while leaving diagonal elements and other sites intact.
    """
    if p <= 0.0:
        return rho
    dim = d**N
    rho_tensor = rho.reshape([d] * (2 * N))

    dephased = np.zeros_like(rho_tensor)
    for k in range(d):
        idx = [slice(None)] * (2 * N)
        idx[site] = k
        idx[N + site] = k
        dephased[tuple(idx)] = rho_tensor[tuple(idx)]

    return ((1 - p) * rho_tensor + p * dephased).reshape(dim, dim)


class QuditGKSLNoisySimulator(QuditGKSLSimulator):
    """Qudit GKSL simulator with hardware depolarization and dephasing noise.

    Inherits the Stinespring + 2nd-order Trotter approach from QuditGKSLSimulator
    and adds per-gate local noise channels to model hardware imperfections.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
        p_depol: depolarization probability per 2-qudit gate (default 0.01 = 1%)
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

        # Pre-compute which molecules each Lindblad operator acts on
        self._lindblad_sites = self._compute_lindblad_sites()

    def _compute_lindblad_sites(self) -> list[list[int]]:
        """Determine the molecule indices each Lindblad operator acts on.

        Structure from build_lindblad_operators:
          - 2*len(neighbors) TTA operators: each pair (i,j) contributes 2 ops
          - 5*N single-site operators: fluorescence, phosphorescence, IC, ISC_ST, ISC_TS
        """
        sites: list[list[int]] = []
        # TTA: 2 per neighbor pair
        for i, j in self.params.neighbors:
            sites.append([i, j])
            sites.append([i, j])
        # Single-site: 5 types × N molecules
        for _channel in range(5):
            for mol in range(self.params.N_molecules):
                sites.append([mol])
        assert len(sites) == len(self.lindblad_ops)
        return sites

    def _apply_gate_noise_single(self, rho: np.ndarray, site: int) -> np.ndarray:
        """Apply hardware noise after a gate involving a single molecule."""
        d = self.params.d
        N = self.params.N_molecules
        rho = _apply_local_depolarization_single(rho, site, d, N, self.p_depol)
        if self.p_dephasing > 0.0:
            rho = _apply_local_dephasing_single(rho, site, d, N, self.p_dephasing)
        return rho

    def _apply_gate_noise_pair(
        self, rho: np.ndarray, site_a: int, site_b: int
    ) -> np.ndarray:
        """Apply hardware noise after a gate involving a molecule pair."""
        d = self.params.d
        N = self.params.N_molecules
        rho = _apply_local_depolarization_pair(rho, site_a, site_b, d, N, self.p_depol)
        if self.p_dephasing > 0.0:
            rho = _apply_local_dephasing_single(rho, site_a, d, N, self.p_dephasing)
            rho = _apply_local_dephasing_single(rho, site_b, d, N, self.p_dephasing)
        return rho

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """Symmetric Trotter step with palindromic Lindblad ordering and hardware noise.

        Same palindromic structure as base class but applies local noise after
        each gate-like operation:
          1. Half Hamiltonian + transfer gate noise (depol + dephasing)
          2. Forward half-step Lindblad channels + per-channel gate noise
          3. Reverse half-step Lindblad channels + per-channel gate noise (palindromic)
          4. Half Hamiltonian + transfer gate noise (depol + dephasing)

        Uses precomputed half-step unitaries from ``_precompute_unitaries``.
        Noise placement matches ``QuditGKSLNoisyShotSimulator``.
        """
        d = self.params.d
        N = self.params.N_molecules

        # --- Half Hamiltonian ---
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        for i, j in self.params.neighbors:
            rho = _apply_local_depolarization_pair(
                rho, i, j, d, N, self.p_depol
            )
            if self.p_dephasing > 0.0:
                rho = _apply_local_dephasing_single(rho, i, d, N, self.p_dephasing)
                rho = _apply_local_dephasing_single(rho, j, d, N, self.p_dephasing)

        # --- Forward half-step Lindblad channels ---
        for k, U_stine_half in enumerate(self._U_stines_half):
            rho = apply_stinespring_to_density_matrix(rho, U_stine_half)
            sites = self._lindblad_sites[k]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    rho = self._apply_gate_noise_single(rho, sites[0])
            else:
                rho = self._apply_gate_noise_pair(rho, sites[0], sites[1])

        # --- Reverse half-step Lindblad channels (palindromic) ---
        n_channels = len(self._U_stines_half)
        for k_rev in range(n_channels - 1, -1, -1):
            U_stine_half = self._U_stines_half[k_rev]
            rho = apply_stinespring_to_density_matrix(rho, U_stine_half)
            sites = self._lindblad_sites[k_rev]
            if len(sites) == 1:
                if not self.depol_pair_only:
                    rho = self._apply_gate_noise_single(rho, sites[0])
            else:
                rho = self._apply_gate_noise_pair(rho, sites[0], sites[1])

        # --- Half Hamiltonian ---
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        for i, j in self.params.neighbors:
            rho = _apply_local_depolarization_pair(
                rho, i, j, d, N, self.p_depol
            )
            if self.p_dephasing > 0.0:
                rho = _apply_local_dephasing_single(rho, i, d, N, self.p_dephasing)
                rho = _apply_local_dephasing_single(rho, j, d, N, self.p_dephasing)

        return rho

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run noisy Qudit GKSL simulation.

        Returns the same dict structure as QuditGKSLSimulator.simulate()
        plus noise_params with p_depol and p_dephasing.
        """
        result = super().simulate(t_max=t_max, n_steps=n_steps, initial_state=initial_state)
        result["method"] = "qudit_gksl_noisy"
        result["noise_params"] = {
            "p_depol": self.p_depol,
            "p_dephasing": self.p_dephasing,
            "depol_pair_only": self.depol_pair_only,
        }
        return result
