"""Physical parameters for GKSL-Lindblad quantum dynamics of TTA-UC phenomena.

This module defines the complete parameter set for open quantum system simulation
of 4-molecule linear chain with triplet-triplet annihilation upconversion (TTA-UC).

All parameters follow the implementation specification:
  tutorials/doc/GKSL/TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md

Units:
  - Energies: eV
  - Time: fs (femtoseconds)
  - Rate constants: eV/ℏ
  - ℏ = 0.6582119569 eV·fs
"""

from __future__ import annotations

import numpy as np


class GKSLPhysicalParameters:
    """Complete physical parameters for GKSL-Lindblad TTA-UC simulation.

    Attributes:
        N_molecules: Number of molecules in linear chain.
        E_T: Triplet state energy (eV).
        E_S: Excited singlet state energy (eV).
        V: Dexter energy transfer integral (eV).
        gamma_TTA: TTA dissipation rate (eV/ℏ).
        Gamma_fl: Fluorescence emission rate (eV/ℏ).
        Gamma_ph: Phosphorescence emission rate (eV/ℏ).
        k_IC: Internal conversion rate (eV/ℏ).
        k_ISC_ST: Intersystem crossing S₁→T₁ rate (eV/ℏ).
        k_ISC_TS: Intersystem crossing T₁→S₀ rate (eV/ℏ).
        hbar: Reduced Planck constant (eV·fs).
        neighbors: List of adjacent molecule pairs.
        T_total: Total simulation time (fs).
        N_steps: Number of Trotter steps.
        dt: Time step (fs).
        initial_state_type: Initial state configuration.
    """

    def __init__(
        self,
        *,
        N_molecules: int = 4,
        E_T: float = 1.5,
        E_S: float = 3.0,
        V: float = 0.1,
        gamma_TTA: float = 0.05,
        Gamma_fl: float = 0.01,
        Gamma_ph: float = 1e-6,
        k_IC: float = 0.005,
        k_ISC_ST: float = 0.003,
        k_ISC_TS: float = 1e-5,
        T_total: float = 100.0,
        N_steps: int = 100,
        initial_state_type: str = "edge_triplet",
    ) -> None:
        # System size
        self.N_molecules = N_molecules
        self.dim_single = 3  # |S0>, |T1>, |S1>
        self.dim_total = self.dim_single**N_molecules  # 3^N = 81 for N=4

        # Energy levels (eV)
        self.E_T = E_T  # T1 energy
        self.E_S = E_S  # S1 energy
        # E_S0 = 0 (reference)

        # Coherent interaction (eV)
        self.V = V  # Dexter energy transfer integral

        # Dissipation rates (eV/ℏ)
        self.gamma_TTA = gamma_TTA
        self.Gamma_fl = Gamma_fl
        self.Gamma_ph = Gamma_ph
        self.k_IC = k_IC
        self.k_ISC_ST = k_ISC_ST
        self.k_ISC_TS = k_ISC_TS

        # Physical constants
        self.hbar = 0.6582119569  # eV·fs

        # System geometry (linear chain)
        self.neighbors = [(i, i + 1) for i in range(N_molecules - 1)]

        # Simulation parameters
        self.T_total = T_total
        self.N_steps = N_steps
        self.dt = T_total / N_steps

        # Initial state
        self.initial_state_type = initial_state_type

    def get_initial_density_matrix(self) -> np.ndarray:
        """Construct initial density matrix ρ(0) = |ψ₀⟩⟨ψ₀|.

        For 'edge_triplet': |1,0,0,1⟩ = T₁-S₀-S₀-T₁
        """
        psi = np.zeros(self.dim_total, dtype=complex)

        if self.initial_state_type == "edge_triplet":
            # |1,0,0,1⟩ in qutrit basis
            # Index = 1*27 + 0*9 + 0*3 + 1 = 28
            idx = self._state_index([1, 0, 0, 1])
            psi[idx] = 1.0
        elif self.initial_state_type == "all_triplet":
            # |1,1,1,1⟩
            idx = self._state_index([1, 1, 1, 1])
            psi[idx] = 1.0
        else:
            msg = f"Unknown initial state type: {self.initial_state_type}"
            raise ValueError(msg)

        return np.outer(psi, psi.conj())

    def _state_index(self, states: list[int]) -> int:
        """Convert list of per-molecule states to global index.

        Args:
            states: List of states [s0, s1, s2, s3] where each si ∈ {0,1,2}.

        Returns:
            Global index in the tensor product space.
        """
        idx = 0
        for s in states:
            idx = idx * self.dim_single + s
        return idx

    def get_initial_state_vector(self) -> np.ndarray:
        """Get initial state vector |ψ₀⟩."""
        psi = np.zeros(self.dim_total, dtype=complex)
        if self.initial_state_type == "edge_triplet":
            idx = self._state_index([1, 0, 0, 1])
            psi[idx] = 1.0
        elif self.initial_state_type == "all_triplet":
            idx = self._state_index([1, 1, 1, 1])
            psi[idx] = 1.0
        else:
            msg = f"Unknown initial state type: {self.initial_state_type}"
            raise ValueError(msg)
        return psi

    def summary(self) -> str:
        """Return human-readable parameter summary."""
        lines = [
            "=== GKSL Physical Parameters ===",
            f"N_molecules: {self.N_molecules}",
            f"Hilbert space dim: {self.dim_total}",
            f"Energy levels: E_T={self.E_T} eV, E_S={self.E_S} eV",
            f"Transfer integral: V={self.V} eV",
            f"TTA rate: γ_TTA={self.gamma_TTA} eV/ℏ",
            f"Fluorescence: Γ_fl={self.Gamma_fl} eV/ℏ",
            f"Phosphorescence: Γ_ph={self.Gamma_ph} eV/ℏ",
            f"Internal conversion: k_IC={self.k_IC} eV/ℏ",
            f"ISC S→T: k_ISC_ST={self.k_ISC_ST} eV/ℏ",
            f"ISC T→S: k_ISC_TS={self.k_ISC_TS} eV/ℏ",
            f"Simulation: T={self.T_total} fs, N_steps={self.N_steps}, dt={self.dt} fs",
            f"Initial state: {self.initial_state_type}",
            f"Neighbor pairs: {self.neighbors}",
        ]
        return "\n".join(lines)
