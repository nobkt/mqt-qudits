"""Qubit GKSL noisy simulator with hardware noise model.

Scenario 8: Qubit-based GKSL-Lindblad with hardware depolarization + thermal relaxation.

Extends QubitGKSLSimulator by applying per-gate local noise channels after
each 2-qubit gate operation in the Trotter step. Since the matrix-level
simulation operates in the 81-dim qutrit Hilbert space, the depolarization
is applied as a local qutrit-space channel on the molecule subsystem(s)
involved in each gate.

Physical distinction:
  - Lindblad dissipators (TTA, fluorescence, etc.) model real physical processes
  - Hardware noise (depolarization, thermal relaxation) models gate imperfections
  These are conceptually independent and both included in this simulator.

Noise channels (2-qubit gates only; 1-qubit gates are ideal):
  - Depolarization: per-gate local depolarization on the molecule subsystem(s)
  - Thermal relaxation: amplitude damping toward |S0> (ground state) per molecule,
    modeling T1 energy decay during gate execution

Gate-noise mapping for the Trotter step:
  - Hamiltonian transfer: one local pair depolarization per nearest-neighbor pair
  - Stinespring channel: one local depolarization on the acted subsystem
    (single molecule for site-local operators, molecule pair for TTA)
  - Thermal relaxation: applied per molecule after each Trotter step
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_noisy_simulator import (
    _apply_local_depolarization_pair,
    _apply_local_depolarization_single,
)


def _apply_thermal_relaxation_single(
    rho: np.ndarray, site: int, d: int, N: int, p_reset: float
) -> np.ndarray:
    """Apply thermal relaxation (amplitude damping toward |0>) on a single molecule.

    Models T1 energy decay during gate execution. In the qutrit encoding
    (|0>=S0, |1>=T1, |2>=S1), thermal relaxation causes transitions
    |1> -> |0> and |2> -> |0> with probability p_reset.

    Kraus operators for molecule k (in the {|0>,|1>,|2>} basis):
      K0 = |0><0| + sqrt(1-p)*|1><1| + sqrt(1-p)*|2><2|  (no decay)
      K1 = sqrt(p) * |0><1|  (T1 -> S0)
      K2 = sqrt(p) * |0><2|  (S1 -> S0)

    Parameters:
        rho: d^N × d^N density matrix
        site: molecule index
        d: local dimension (3)
        N: number of molecules
        p_reset: reset probability = 1 - exp(-t_gate/T1)
    """
    if p_reset <= 0.0:
        return rho
    dim = d**N
    rho_tensor = rho.reshape([d] * (2 * N))
    result = np.zeros_like(rho_tensor)

    sq = np.sqrt(1.0 - p_reset)
    # Build 3x3 Kraus matrices for the single site
    # K0[a,b] coefficients for |a><b| contribution
    # K0 = diag(1, sqrt(1-p), sqrt(1-p))
    k0_diag = np.array([1.0, sq, sq], dtype=np.complex128)
    # K1 = sqrt(p) * |0><1|
    # K2 = sqrt(p) * |0><2|
    sqp = np.sqrt(p_reset)

    # Apply Kraus: rho' = K0 rho K0† + K1 rho K1† + K2 rho K2†
    # This is done element-wise on the site's local indices.
    # For the tensor with indices (..., ket_site, ..., bra_site, ...):
    # (K rho K†)_{a,b}^{site} = sum_{c,d} K_{a,c} * rho_{c,d}^{site} * conj(K_{b,d})

    for a in range(d):
        for b in range(d):
            idx_out = [slice(None)] * (2 * N)
            idx_out[site] = a
            idx_out[N + site] = b

            # K0 contribution: K0[a,c] * K0[b,d]* = k0_diag[a] * k0_diag[b] when c=a, d=b
            idx_in = [slice(None)] * (2 * N)
            idx_in[site] = a
            idx_in[N + site] = b
            result[tuple(idx_out)] += (
                k0_diag[a] * np.conj(k0_diag[b]) * rho_tensor[tuple(idx_in)]
            )

            # K1 contribution: K1 = sqp * |0><1|
            # K1[a,c] = sqp if a=0,c=1 else 0
            # K1[b,d]* = sqp if b=0,d=1 else 0
            if a == 0 and b == 0:
                idx_k1 = [slice(None)] * (2 * N)
                idx_k1[site] = 1  # c=1
                idx_k1[N + site] = 1  # d=1
                result[tuple(idx_out)] += p_reset * rho_tensor[tuple(idx_k1)]

            # K2 contribution: K2 = sqp * |0><2|
            if a == 0 and b == 0:
                idx_k2 = [slice(None)] * (2 * N)
                idx_k2[site] = 2  # c=2
                idx_k2[N + site] = 2  # d=2
                result[tuple(idx_out)] += p_reset * rho_tensor[tuple(idx_k2)]

    return result.reshape(dim, dim)


class QubitGKSLNoisySimulator(QubitGKSLSimulator):
    """Qubit GKSL simulator with hardware depolarization and thermal relaxation.

    Inherits the Stinespring + 2nd-order Trotter approach from QubitGKSLSimulator
    and adds per-gate local noise channels to model hardware imperfections.

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
        p_depol: depolarization probability per 2-qubit gate (default 0.01 = 1%)
        T1: energy relaxation time in natural units (default None = no relaxation)
        T2: dephasing time in natural units (default None = no relaxation)
        t_gate: 2-qubit gate time in natural units (default 300.0 fs)

    Note on thermal relaxation:
        For typical superconducting qubit parameters (T1=50μs, gate_time=300fs),
        p_reset = 1 - exp(-t_gate/T1) ≈ 6e-9, which is negligible compared to
        depolarization (p=0.01). Thermal relaxation is included for completeness.
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        p_depol: float = 0.01,
        T1: float | None = None,
        T2: float | None = None,
        t_gate: float = 300.0,
    ) -> None:
        super().__init__(params)
        if p_depol < 0.0 or p_depol > 1.0:
            msg = f"p_depol must be in [0, 1], got {p_depol}"
            raise ValueError(msg)
        self.p_depol = p_depol
        self.T1 = T1
        self.T2 = T2
        self.t_gate = t_gate

        # Compute thermal relaxation probability per gate
        if T1 is not None and T1 > 0:
            self.p_reset = 1.0 - np.exp(-t_gate / T1)
        else:
            self.p_reset = 0.0

        # Pre-compute Lindblad operator site mapping
        self._lindblad_sites = self._compute_lindblad_sites()

    def _compute_lindblad_sites(self) -> list[list[int]]:
        """Determine molecule indices for each Lindblad operator.

        Same structure as build_lindblad_operators:
          - 2*len(neighbors) TTA operators: each pair (i,j) -> 2 ops
          - 5*N single-site operators: fl, ph, IC, ISC_ST, ISC_TS
        """
        sites: list[list[int]] = []
        for i, j in self.params.neighbors:
            sites.append([i, j])
            sites.append([i, j])
        for _channel in range(5):
            for mol in range(self.params.N_molecules):
                sites.append([mol])
        assert len(sites) == len(self.lindblad_ops)
        return sites

    def _trotter_step(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """2nd-order Trotter step with hardware noise.

        Gate noise is applied after each 2-body gate operation:
          1. Half Hamiltonian + NN pair depolarization
          2. Lindblad channels + per-channel depolarization
          3. Half Hamiltonian + NN pair depolarization
          4. Thermal relaxation on all molecules (accumulated gate time)
        """
        d = self.params.d
        N = self.params.N_molecules

        # --- Half Hamiltonian ---
        rho = self._apply_hamiltonian_step(rho, dt / 2)
        for i, j in self.params.neighbors:
            rho = _apply_local_depolarization_pair(rho, i, j, d, N, self.p_depol)

        # --- All Lindblad channels ---
        for k, (L_op, _gamma) in enumerate(self.lindblad_ops):
            rho = self._apply_lindblad_stinespring(rho, L_op, dt)
            sites = self._lindblad_sites[k]
            if len(sites) == 1:
                rho = _apply_local_depolarization_single(
                    rho, sites[0], d, N, self.p_depol
                )
            else:
                rho = _apply_local_depolarization_pair(
                    rho, sites[0], sites[1], d, N, self.p_depol
                )

        # --- Half Hamiltonian ---
        rho = self._apply_hamiltonian_step(rho, dt / 2)
        for i, j in self.params.neighbors:
            rho = _apply_local_depolarization_pair(rho, i, j, d, N, self.p_depol)

        # --- Thermal relaxation (all molecules, accumulated gate time) ---
        if self.p_reset > 0.0:
            for mol in range(N):
                rho = _apply_thermal_relaxation_single(
                    rho, mol, d, N, self.p_reset
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
            "T1": self.T1,
            "T2": self.T2,
            "t_gate": self.t_gate,
            "p_reset": self.p_reset,
        }
        return result
