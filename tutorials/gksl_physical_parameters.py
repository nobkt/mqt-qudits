"""GKSLPhysicalParameters class for TTA-UC GKSL-Lindblad quantum dynamics."""

from __future__ import annotations

from typing import Any

# hbar in eV·fs units (NIST CODATA value, truncated)
HBAR_EV_FS = 0.6582

# Tolerance for 2*E_T ≈ E_S energy relation check
ENERGY_RATIO_TOLERANCE = 0.1


class GKSLPhysicalParameters:
    """Physical parameters for TTA-UC GKSL-Lindblad quantum dynamics simulation."""

    def __init__(
        self,
        *,
        E_T: float = 1.5,
        E_S: float = 3.0,
        V: float = 0.1,
        gamma_TTA: float = 0.05,
        Gamma_fl: float = 0.01,
        Gamma_ph: float = 1e-6,
        k_IC: float = 0.005,
        k_ISC_ST: float = 0.003,
        k_ISC_TS: float = 1e-5,
        N_molecules: int = 4,
        d: int = 3,
        with_boson: bool = False,
        n_max: int = 2,
        omega_ph: float = 0.15,
        g_eph: float = 0.02,
    ) -> None:
        # Energy parameters (eV)
        self.E_T = E_T
        self.E_S = E_S
        self.V = V

        # Dissipation rates (eV/hbar)
        self.gamma_TTA = gamma_TTA
        self.Gamma_fl = Gamma_fl
        self.Gamma_ph = Gamma_ph
        self.k_IC = k_IC
        self.k_ISC_ST = k_ISC_ST
        self.k_ISC_TS = k_ISC_TS

        # System parameters
        self.N_molecules = N_molecules
        self.d = d

        # Boson parameters
        self.with_boson = with_boson
        self.n_max = n_max
        self.omega_ph = omega_ph
        self.g_eph = g_eph

    @property
    def hbar(self) -> float:
        """Natural units where hbar = 1."""
        return 1.0

    @property
    def neighbors(self) -> list[tuple[int, int]]:
        """Nearest-neighbor pairs."""
        return [(i, i + 1) for i in range(self.N_molecules - 1)]

    def validate(self) -> list[str]:
        """Validate physical parameters. Returns list of error messages (empty if valid)."""
        errors: list[str] = []

        # 1. Energy positivity
        if self.E_T <= 0:
            errors.append("E_T must be positive")
        if self.E_S <= 0:
            errors.append("E_S must be positive")
        if self.E_S <= self.E_T:
            errors.append("E_S must be greater than E_T")

        # 2. Energy relation 2*E_T ≈ E_S (10% tolerance)
        if self.E_T > 0 and self.E_S > 0:
            ratio = abs(2 * self.E_T - self.E_S) / self.E_S
            if ratio > ENERGY_RATIO_TOLERANCE:
                errors.append(f"2*E_T should approximate E_S (ratio deviation {ratio:.3f} > {ENERGY_RATIO_TOLERANCE})")

        # 3. Dissipation rate non-negativity
        errors.extend(
            f"{name} must be non-negative"
            for name in ("gamma_TTA", "Gamma_fl", "Gamma_ph", "k_IC", "k_ISC_ST", "k_ISC_TS")
            if getattr(self, name) < 0
        )

        # 4. Time scale hierarchy: gamma_TTA >= Gamma_fl >= k_IC >= k_ISC_ST
        if self.gamma_TTA < self.Gamma_fl:
            errors.append("Time scale hierarchy violated: gamma_TTA < Gamma_fl")
        if self.Gamma_fl < self.k_IC:
            errors.append("Time scale hierarchy violated: Gamma_fl < k_IC")
        if self.k_IC < self.k_ISC_ST:
            errors.append("Time scale hierarchy violated: k_IC < k_ISC_ST")

        # 5. System parameters
        if self.N_molecules < 2:
            errors.append("N_molecules must be >= 2")
        # d=3 required: ground (S0), triplet (T1), singlet (S1) states per molecule
        if self.d != 3:
            errors.append("d must be 3")

        # 6. Boson parameters
        if self.with_boson:
            if self.n_max < 1:
                errors.append("n_max must be >= 1 when with_boson is True")
            if self.omega_ph <= 0:
                errors.append("omega_ph must be positive when with_boson is True")
            if self.g_eph < 0:
                errors.append("g_eph must be non-negative when with_boson is True")

        # 7. Weak coupling: max_dissipation < 0.5 * min_energy
        # Skip when all dissipation rates are zero (pure unitary limit)
        # or when V=0 (no inter-molecular coupling; use E_T as energy scale)
        max_dissipation = max(self.gamma_TTA, self.Gamma_fl, self.Gamma_ph, self.k_IC, self.k_ISC_ST, self.k_ISC_TS)
        if max_dissipation > 0:
            energy_scales = [self.E_T]
            if self.V > 0:
                energy_scales.append(self.V)
            min_energy = min(energy_scales)
            if max_dissipation > 0.5 * min_energy:
                errors.append(
                    f"Weak coupling violated: max_dissipation ({max_dissipation}) >= 0.5 * min_energy ({0.5 * min_energy})"
                )

        return errors

    def get_hilbert_space_dim(self) -> int:
        """Return the Hilbert space dimension."""
        dim = int(self.d**self.N_molecules)
        if self.with_boson:
            dim *= int((self.n_max + 1) ** self.N_molecules)
        return dim

    def get_time_unit(self) -> float:
        """Return hbar in eV·fs units."""
        return HBAR_EV_FS

    def to_dict(self) -> dict[str, Any]:
        """Return all parameters as a dictionary."""
        return {
            "E_T": self.E_T,
            "E_S": self.E_S,
            "V": self.V,
            "gamma_TTA": self.gamma_TTA,
            "Gamma_fl": self.Gamma_fl,
            "Gamma_ph": self.Gamma_ph,
            "k_IC": self.k_IC,
            "k_ISC_ST": self.k_ISC_ST,
            "k_ISC_TS": self.k_ISC_TS,
            "N_molecules": self.N_molecules,
            "d": self.d,
            "with_boson": self.with_boson,
            "n_max": self.n_max,
            "omega_ph": self.omega_ph,
            "g_eph": self.g_eph,
        }

    def __repr__(self) -> str:
        return (
            f"GKSLPhysicalParameters("
            f"E_T={self.E_T}, E_S={self.E_S}, V={self.V}, "
            f"gamma_TTA={self.gamma_TTA}, Gamma_fl={self.Gamma_fl}, Gamma_ph={self.Gamma_ph}, "
            f"k_IC={self.k_IC}, k_ISC_ST={self.k_ISC_ST}, k_ISC_TS={self.k_ISC_TS}, "
            f"N_molecules={self.N_molecules}, d={self.d}, "
            f"with_boson={self.with_boson}, n_max={self.n_max}, "
            f"omega_ph={self.omega_ph}, g_eph={self.g_eph})"
        )
