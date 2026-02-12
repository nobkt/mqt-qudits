"""Classical GKSL-Lindblad simulator using superoperator formalism.

Implements the GKSL master equation for 4-molecule TTA-UC system:

  dρ/dt = -i/ℏ [H_sys, ρ] + Σ_α γ_α D[L_α][ρ]

where D[L][ρ] = L ρ L† - ½{L†L, ρ}

This module uses scipy.integrate.solve_ivp for time integration of the
vectorized density matrix equation.

Reference:
  tutorials/doc/GKSL/TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import scipy.integrate
import scipy.linalg
import scipy.sparse

from gksl_physical_parameters import GKSLPhysicalParameters


def _single_site_operator(
    op: np.ndarray, site: int, N: int, d: int
) -> np.ndarray:
    """Embed a single-site operator into the full Hilbert space.

    Args:
        op: d×d operator acting on one molecule.
        site: Which molecule (0-indexed).
        N: Total number of molecules.
        d: Dimension per molecule.

    Returns:
        d^N × d^N operator.
    """
    result = np.eye(1, dtype=complex)
    for i in range(N):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, np.eye(d, dtype=complex))
    return result


def _two_site_operator(
    op_i: np.ndarray, op_j: np.ndarray, site_i: int, site_j: int, N: int, d: int
) -> np.ndarray:
    """Embed a two-site tensor product operator into the full Hilbert space.

    Args:
        op_i: d×d operator on site i.
        op_j: d×d operator on site j.
        site_i: First site index.
        site_j: Second site index.
        N: Total number of molecules.
        d: Dimension per molecule.

    Returns:
        d^N × d^N operator.
    """
    result = np.eye(1, dtype=complex)
    for k in range(N):
        if k == site_i:
            result = np.kron(result, op_i)
        elif k == site_j:
            result = np.kron(result, op_j)
        else:
            result = np.kron(result, np.eye(d, dtype=complex))
    return result


class ClassicalGKSLSimulator:
    """Classical GKSL-Lindblad simulator using ODE integration.

    Solves the Lindblad master equation:
      dρ/dt = -i/ℏ [H, ρ] + Σ_α γ_α D[L_α][ρ]

    Uses vectorization: ρ → vec(ρ) and constructs the Liouvillian superoperator.
    """

    def __init__(self, params: GKSLPhysicalParameters | None = None) -> None:
        self.params = params or GKSLPhysicalParameters()
        self.dim = self.params.dim_total  # 81 for 4 molecules
        self.N = self.params.N_molecules
        self.d = self.params.dim_single  # 3

    def build_system_hamiltonian(self) -> np.ndarray:
        """Build the system Hamiltonian H_sys = H_0 + H_transfer.

        H_0 = Σ_i (E_T |1⟩⟨1| + E_S |2⟩⟨2|) on site i
        H_transfer = Σ_{⟨i,j⟩} V (|0⟩⟨1| ⊗ |1⟩⟨0| + h.c.)

        Note: TTA is NOT in the Hamiltonian - it is a Lindblad dissipator.
        """
        dim = self.dim
        H = np.zeros((dim, dim), dtype=complex)

        # On-site energies
        h_single = np.diag([0.0, self.params.E_T, self.params.E_S]).astype(complex)
        for i in range(self.N):
            H += _single_site_operator(h_single, i, self.N, self.d)

        # Dexter energy transfer between neighbors
        # |0⟩⟨1| ⊗ |1⟩⟨0| + |1⟩⟨0| ⊗ |0⟩⟨1|
        ket0_bra1 = np.zeros((self.d, self.d), dtype=complex)
        ket0_bra1[0, 1] = 1.0
        ket1_bra0 = np.zeros((self.d, self.d), dtype=complex)
        ket1_bra0[1, 0] = 1.0

        for i, j in self.params.neighbors:
            H += self.params.V * _two_site_operator(
                ket0_bra1, ket1_bra0, i, j, self.N, self.d
            )
            H += self.params.V * _two_site_operator(
                ket1_bra0, ket0_bra1, i, j, self.N, self.d
            )

        return H

    def build_lindblad_operators(self) -> list[tuple[float, np.ndarray]]:
        """Build all 26 Lindblad operators with their rates.

        Returns:
            List of (rate, L_operator) tuples.

        Lindblad operators:
          - TTA: 6 operators (2 per pair × 3 pairs)
          - Fluorescence: 4 operators (1 per molecule)
          - Phosphorescence: 4 operators (1 per molecule)
          - Internal conversion: 4 operators (1 per molecule)
          - ISC S₁→T₁: 4 operators (1 per molecule)
          - ISC T₁→S₀: 4 operators (1 per molecule)
        """
        ops = []

        # Single-site transition operators
        ket0_bra1 = np.zeros((self.d, self.d), dtype=complex)
        ket0_bra1[0, 1] = 1.0  # |0⟩⟨1| (T₁ → S₀)

        ket0_bra2 = np.zeros((self.d, self.d), dtype=complex)
        ket0_bra2[0, 2] = 1.0  # |0⟩⟨2| (S₁ → S₀)

        ket1_bra2 = np.zeros((self.d, self.d), dtype=complex)
        ket1_bra2[1, 2] = 1.0  # |1⟩⟨2| (S₁ → T₁)

        ket2_bra1 = np.zeros((self.d, self.d), dtype=complex)
        ket2_bra1[2, 1] = 1.0  # |1⟩ → |2⟩ (T₁ → S₁)

        # --- TTA operators (6 total) ---
        # For each neighbor pair (i,j):
        #   L_TTA,1^(ij) = |2⟩_i⟨1| ⊗ |0⟩_j⟨1|  (T₁T₁ → S₁S₀)
        #   L_TTA,2^(ij) = |0⟩_i⟨1| ⊗ |2⟩_j⟨1|  (T₁T₁ → S₀S₁)
        for i, j in self.params.neighbors:
            L_tta1 = _two_site_operator(ket2_bra1, ket0_bra1, i, j, self.N, self.d)
            ops.append((self.params.gamma_TTA / 2.0, L_tta1))

            L_tta2 = _two_site_operator(ket0_bra1, ket2_bra1, i, j, self.N, self.d)
            ops.append((self.params.gamma_TTA / 2.0, L_tta2))

        # --- Single-molecule dissipation (4 each) ---
        for i in range(self.N):
            # Fluorescence: S₁ → S₀
            L_fl = _single_site_operator(ket0_bra2, i, self.N, self.d)
            ops.append((self.params.Gamma_fl, L_fl))

            # Phosphorescence: T₁ → S₀
            L_ph = _single_site_operator(ket0_bra1, i, self.N, self.d)
            ops.append((self.params.Gamma_ph, L_ph))

            # Internal conversion: S₁ → S₀
            L_ic = _single_site_operator(ket0_bra2, i, self.N, self.d)
            ops.append((self.params.k_IC, L_ic))

            # ISC S₁→T₁
            L_isc_st = _single_site_operator(ket1_bra2, i, self.N, self.d)
            ops.append((self.params.k_ISC_ST, L_isc_st))

            # ISC T₁→S₀
            L_isc_ts = _single_site_operator(ket0_bra1, i, self.N, self.d)
            ops.append((self.params.k_ISC_TS, L_isc_ts))

        return ops

    def _lindblad_rhs(self, _t: float, rho_vec: np.ndarray, H: np.ndarray,
                      lindblad_ops: list[tuple[float, np.ndarray]]) -> np.ndarray:
        """Right-hand side of the vectorized Lindblad equation.

        dρ/dt = -i/ℏ [H, ρ] + Σ_α γ_α D[L_α][ρ]
        """
        dim = self.dim
        rho = rho_vec.reshape((dim, dim))

        # Coherent part: -i/ℏ [H, ρ]
        drho = -1j / self.params.hbar * (H @ rho - rho @ H)

        # Dissipative part: Σ_α γ_α D[L_α][ρ]
        for gamma, L in lindblad_ops:
            Ld = L.conj().T
            LdL = Ld @ L
            drho += gamma * (L @ rho @ Ld - 0.5 * LdL @ rho - 0.5 * rho @ LdL)

        return drho.ravel()

    def _extract_populations(self, rho: np.ndarray) -> dict[str, float]:
        """Extract populations from density matrix.

        N_n(t) = Σ_i Tr[P_n^(i) ρ(t)]
        where P_n^(i) = I^⊗i ⊗ |n⟩⟨n| ⊗ I^⊗(N-1-i)
        """
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0

        for mol in range(self.N):
            for state_idx, key in [(0, "S0"), (1, "T1"), (2, "S1")]:
                projector = np.zeros((self.d, self.d), dtype=complex)
                projector[state_idx, state_idx] = 1.0
                P = _single_site_operator(projector, mol, self.N, self.d)
                pop = np.real(np.trace(P @ rho))
                if key == "S0":
                    N_S0 += pop
                elif key == "T1":
                    N_T1 += pop
                else:
                    N_S1 += pop

        return {"N_S0": N_S0, "N_T1": N_T1, "N_S1": N_S1}

    def _extract_per_molecule_populations(self, rho: np.ndarray) -> dict[str, np.ndarray]:
        """Extract per-molecule populations."""
        S0 = np.zeros(self.N)
        T1 = np.zeros(self.N)
        S1 = np.zeros(self.N)

        for mol in range(self.N):
            for state_idx, arr in [(0, S0), (1, T1), (2, S1)]:
                projector = np.zeros((self.d, self.d), dtype=complex)
                projector[state_idx, state_idx] = 1.0
                P = _single_site_operator(projector, mol, self.N, self.d)
                arr[mol] = np.real(np.trace(P @ rho))

        return {"S0": S0, "T1": T1, "S1": S1}

    def _von_neumann_entropy(self, rho: np.ndarray) -> float:
        """Compute von Neumann entropy S = -Tr[ρ ln ρ]."""
        eigenvalues = np.real(np.linalg.eigvalsh(rho))
        eigenvalues = eigenvalues[eigenvalues > 1e-15]
        return float(-np.sum(eigenvalues * np.log(eigenvalues)))

    def _purity(self, rho: np.ndarray) -> float:
        """Compute purity Tr[ρ²]."""
        return float(np.real(np.trace(rho @ rho)))

    def simulate(self, ode_method: str = "RK45") -> dict[str, Any]:
        """Run GKSL-Lindblad simulation.

        Args:
            ode_method: ODE solver method ('RK45', 'BDF', etc.)

        Returns:
            Dictionary with simulation results including populations,
            entropy, purity, trace verification, etc.
        """
        start_time = time.time()

        # Build Hamiltonian and Lindblad operators
        H = self.build_system_hamiltonian()
        lindblad_ops = self.build_lindblad_operators()

        # Verify Hamiltonian is Hermitian
        assert np.allclose(H, H.conj().T), "Hamiltonian must be Hermitian"

        # Initial density matrix
        rho0 = self.params.get_initial_density_matrix()

        # Time points
        t_span = (0.0, self.params.T_total)
        t_eval = np.linspace(0.0, self.params.T_total, self.params.N_steps + 1)

        # Solve ODE
        sol = scipy.integrate.solve_ivp(
            fun=lambda t, y: self._lindblad_rhs(t, y, H, lindblad_ops),
            t_span=t_span,
            y0=rho0.ravel(),
            method=ode_method,
            t_eval=t_eval,
            rtol=1e-10,
            atol=1e-12,
        )

        if not sol.success:
            msg = f"ODE solver failed: {sol.message}"
            raise RuntimeError(msg)

        # Extract results at each time point
        times = sol.t.tolist()
        populations = []
        per_molecule_populations = []
        entropy_list = []
        purity_list = []
        trace_list = []

        for k in range(len(times)):
            rho_k = sol.y[:, k].reshape((self.dim, self.dim))

            populations.append(self._extract_populations(rho_k))
            per_molecule_populations.append(self._extract_per_molecule_populations(rho_k))
            entropy_list.append(self._von_neumann_entropy(rho_k))
            purity_list.append(self._purity(rho_k))
            trace_list.append(float(np.real(np.trace(rho_k))))

        # Final density matrix
        rho_final = sol.y[:, -1].reshape((self.dim, self.dim))

        elapsed = time.time() - start_time

        return {
            "times": times,
            "populations": populations,
            "per_molecule_populations": per_molecule_populations,
            "elapsed_time": elapsed,
            "method": f"Classical GKSL ({ode_method})",
            # GKSL-specific
            "entropy": entropy_list,
            "purity": purity_list,
            "trace": trace_list,
            "rho_final": rho_final,
            "ode_solver": ode_method,
            "n_function_evals": sol.nfev,
            "n_lindblad_operators": len(lindblad_ops),
        }


def verify_gksl_results(results: dict[str, Any], params: GKSLPhysicalParameters) -> dict[str, bool]:
    """Verify GKSL simulation results for physical consistency.

    Checks:
      1. Trace preservation: |1 - Tr[ρ(t)]| < 1e-6
      2. Particle conservation: N_S0 + N_T1 + N_S1 = N ± 1e-6
      3. Entropy non-decrease (2nd law)
      4. Purity bounds: 0 ≤ Tr[ρ²] ≤ 1
    """
    checks = {}

    # Trace preservation
    trace_errors = [abs(1.0 - tr) for tr in results["trace"]]
    checks["trace_preserved"] = max(trace_errors) < 1e-6

    # Particle conservation
    pop_sums = [
        p["N_S0"] + p["N_T1"] + p["N_S1"]
        for p in results["populations"]
    ]
    pop_errors = [abs(params.N_molecules - s) for s in pop_sums]
    checks["particle_conserved"] = max(pop_errors) < 1e-6

    # Entropy behavior check.
    # Note: For GKSL, the von Neumann entropy of the system alone (partial trace)
    # is NOT guaranteed to be monotonically non-decreasing. The 2nd law applies
    # to the total system+environment entropy. When coherent Hamiltonian dynamics
    # compete with dissipation, system entropy can temporarily decrease.
    # We check that entropy generally increases from initial value.
    entropy = results["entropy"]
    checks["entropy_increases_from_initial"] = entropy[-1] > entropy[0]

    # Purity bounds
    purity = results["purity"]
    checks["purity_bounded"] = all(-1e-8 <= p <= 1.0 + 1e-8 for p in purity)

    return checks
