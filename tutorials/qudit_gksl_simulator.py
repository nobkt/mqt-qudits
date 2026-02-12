"""Qudit GKSL-Lindblad simulator using Stinespring dilation.

Implements the GKSL master equation using qutrit (3-level qudit) quantum circuits
via Stinespring dilation. Each molecule is naturally represented as a single qutrit,
avoiding the forbidden state problem of qubit encoding.

Qutrit encoding (per molecule):
  |S₀⟩ ↔ |0⟩, |T₁⟩ ↔ |1⟩, |S₁⟩ ↔ |2⟩

System: 4 qutrits (81-dimensional Hilbert space)
Ancilla: 26 qubits (1 per Lindblad operator)

Reference:
  tutorials/doc/GKSL/GKSL-Lindblad量子ダイナミクスQudit完全実装理論.md
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import scipy.linalg

from gksl_physical_parameters import GKSLPhysicalParameters


def _single_site_operator(
    op: np.ndarray, site: int, N: int, d: int
) -> np.ndarray:
    """Embed a single-site operator into the full Hilbert space."""
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
    """Embed a two-site tensor product operator into the full Hilbert space."""
    result = np.eye(1, dtype=complex)
    for k in range(N):
        if k == site_i:
            result = np.kron(result, op_i)
        elif k == site_j:
            result = np.kron(result, op_j)
        else:
            result = np.kron(result, np.eye(d, dtype=complex))
    return result


def _build_stinespring_unitary(
    L: np.ndarray, gamma: float, dt: float, hbar: float
) -> np.ndarray:
    """Build Stinespring dilation unitary for a Lindblad operator.

    U = exp(-i θ G) where:
      G = L ⊗ |1⟩⟨0|_E + L† ⊗ |0⟩⟨1|_E
      θ = √(γ dt / ℏ)
    """
    dim = L.shape[0]
    theta = np.sqrt(gamma * dt / hbar)

    G = np.zeros((2 * dim, 2 * dim), dtype=complex)
    G[:dim, dim:] = L.conj().T
    G[dim:, :dim] = L

    return scipy.linalg.expm(-1j * theta * G)


class QuditGKSLSimulator:
    """Qudit-based GKSL-Lindblad simulator using Stinespring dilation.

    Uses natural qutrit representation (3-level per molecule):
      |S₀⟩ = |0⟩, |T₁⟩ = |1⟩, |S₁⟩ = |2⟩

    Advantages over qubit encoding:
      - No forbidden states
      - Natural partial-space rotations
      - Compact representation (4 qutrits vs 8 qubits)
    """

    def __init__(self, params: GKSLPhysicalParameters | None = None) -> None:
        self.params = params or GKSLPhysicalParameters()
        self.N = self.params.N_molecules
        self.d = self.params.dim_single  # 3
        self.dim = self.params.dim_total  # 81

    def build_system_hamiltonian(self) -> np.ndarray:
        """Build system Hamiltonian in qutrit space (81×81)."""
        dim = self.dim
        H = np.zeros((dim, dim), dtype=complex)

        # On-site energies
        h_single = np.diag([0.0, self.params.E_T, self.params.E_S]).astype(complex)
        for i in range(self.N):
            H += _single_site_operator(h_single, i, self.N, self.d)

        # Dexter energy transfer
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
        """Build all 26 Lindblad operators in qutrit space."""
        ops = []

        ket0_bra1 = np.zeros((self.d, self.d), dtype=complex)
        ket0_bra1[0, 1] = 1.0
        ket0_bra2 = np.zeros((self.d, self.d), dtype=complex)
        ket0_bra2[0, 2] = 1.0
        ket1_bra2 = np.zeros((self.d, self.d), dtype=complex)
        ket1_bra2[1, 2] = 1.0
        ket2_bra1 = np.zeros((self.d, self.d), dtype=complex)
        ket2_bra1[2, 1] = 1.0

        # TTA operators (6)
        for i, j in self.params.neighbors:
            L1 = _two_site_operator(ket2_bra1, ket0_bra1, i, j, self.N, self.d)
            ops.append((self.params.gamma_TTA / 2.0, L1))
            L2 = _two_site_operator(ket0_bra1, ket2_bra1, i, j, self.N, self.d)
            ops.append((self.params.gamma_TTA / 2.0, L2))

        # Single-molecule operators (20)
        for i in range(self.N):
            ops.append((self.params.Gamma_fl,
                        _single_site_operator(ket0_bra2, i, self.N, self.d)))
            ops.append((self.params.Gamma_ph,
                        _single_site_operator(ket0_bra1, i, self.N, self.d)))
            ops.append((self.params.k_IC,
                        _single_site_operator(ket0_bra2, i, self.N, self.d)))
            ops.append((self.params.k_ISC_ST,
                        _single_site_operator(ket1_bra2, i, self.N, self.d)))
            ops.append((self.params.k_ISC_TS,
                        _single_site_operator(ket0_bra1, i, self.N, self.d)))

        return ops

    def _precompute_kraus_operators(
        self, lindblad_ops: list[tuple[float, np.ndarray]]
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        """Precompute Kraus operators from Stinespring unitaries."""
        kraus_list = []
        for gamma, L in lindblad_ops:
            dim = L.shape[0]
            U_se = _build_stinespring_unitary(L, gamma, self.params.dt, self.params.hbar)
            K0 = U_se[:dim, :dim]
            K1 = U_se[dim:, :dim]
            kraus_list.append((K0, K1))
        return kraus_list

    def _apply_kraus_channel(
        self, rho: np.ndarray, K0: np.ndarray, K1: np.ndarray
    ) -> np.ndarray:
        """Apply Kraus channel: ρ → K₀ ρ K₀† + K₁ ρ K₁†."""
        return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T

    def _extract_populations(self, rho: np.ndarray) -> dict[str, float]:
        """Extract populations from density matrix."""
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0

        for mol in range(self.N):
            for state_idx, key in [(0, "S0"), (1, "T1"), (2, "S1")]:
                proj = np.zeros((self.d, self.d), dtype=complex)
                proj[state_idx, state_idx] = 1.0
                P = _single_site_operator(proj, mol, self.N, self.d)
                pop = np.real(np.trace(P @ rho))
                if key == "S0":
                    N_S0 += pop
                elif key == "T1":
                    N_T1 += pop
                else:
                    N_S1 += pop

        return {"N_S0": N_S0, "N_T1": N_T1, "N_S1": N_S1}

    def _extract_per_molecule_populations(
        self, rho: np.ndarray
    ) -> dict[str, np.ndarray]:
        """Extract per-molecule populations."""
        S0 = np.zeros(self.N)
        T1 = np.zeros(self.N)
        S1 = np.zeros(self.N)

        for mol in range(self.N):
            for state_idx, arr in [(0, S0), (1, T1), (2, S1)]:
                proj = np.zeros((self.d, self.d), dtype=complex)
                proj[state_idx, state_idx] = 1.0
                P = _single_site_operator(proj, mol, self.N, self.d)
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

    def simulate(self, shots: int | None = None) -> dict[str, Any]:
        """Run Qudit GKSL simulation via Stinespring dilation.

        Uses 2nd-order symmetric Trotter splitting (Strang splitting):
          exp(L_H dt/2) → exp(L_diss dt) → exp(L_H dt/2)

        Args:
            shots: If provided, add shot noise to final populations.

        Returns:
            Dictionary with simulation results.
        """
        start_time = time.time()

        H = self.build_system_hamiltonian()
        lindblad_ops = self.build_lindblad_operators()

        assert np.allclose(H, H.conj().T), "Hamiltonian must be Hermitian"

        rho = self.params.get_initial_density_matrix()

        times = [0.0]
        populations = [self._extract_populations(rho)]
        per_molecule_populations = [self._extract_per_molecule_populations(rho)]
        entropy_list = [self._von_neumann_entropy(rho)]
        purity_list = [self._purity(rho)]
        trace_list = [float(np.real(np.trace(rho)))]

        # Precompute
        U_half = scipy.linalg.expm(
            -1j * H * 0.5 * self.params.dt / self.params.hbar
        )
        kraus_ops = self._precompute_kraus_operators(lindblad_ops)

        for step in range(self.params.N_steps):
            # Strang splitting: H/2 → Dissipation → H/2

            rho = U_half @ rho @ U_half.conj().T

            for K0, K1 in kraus_ops:
                rho = self._apply_kraus_channel(rho, K0, K1)

            rho = U_half @ rho @ U_half.conj().T

            times.append((step + 1) * self.params.dt)
            populations.append(self._extract_populations(rho))
            per_molecule_populations.append(
                self._extract_per_molecule_populations(rho)
            )
            entropy_list.append(self._von_neumann_entropy(rho))
            purity_list.append(self._purity(rho))
            trace_list.append(float(np.real(np.trace(rho))))

        elapsed = time.time() - start_time

        result = {
            "times": times,
            "populations": populations,
            "per_molecule_populations": per_molecule_populations,
            "elapsed_time": elapsed,
            "method": "Qudit GKSL (Stinespring)",
            "entropy": entropy_list,
            "purity": purity_list,
            "trace": trace_list,
            "rho_final": rho,
            "n_lindblad_operators": len(lindblad_ops),
            "n_system_qutrits": self.N,
            "n_ancilla_qubits": len(lindblad_ops),
            "equivalent_qubits": self.N * 2 + len(lindblad_ops),
        }

        if shots is not None:
            result["shots"] = shots
            result["populations_shot"] = self._add_shot_noise(
                populations, shots
            )

        return result

    def _add_shot_noise(
        self, populations: list[dict[str, float]], shots: int
    ) -> list[dict[str, float]]:
        """Add shot noise to populations."""
        rng = np.random.default_rng(42)
        noisy_pops = []
        for pop in populations:
            total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
            if total < 1e-10:
                noisy_pops.append(pop.copy())
                continue

            probs = np.array([
                pop["N_S0"] / self.N,
                pop["N_T1"] / self.N,
                pop["N_S1"] / self.N,
            ])
            probs = np.clip(probs, 0, None)
            probs /= probs.sum()

            n_s0 = 0.0
            n_t1 = 0.0
            n_s1 = 0.0
            for _ in range(self.N):
                counts = rng.multinomial(shots, probs)
                n_s0 += counts[0] / shots
                n_t1 += counts[1] / shots
                n_s1 += counts[2] / shots

            noisy_pops.append({"N_S0": n_s0, "N_T1": n_t1, "N_S1": n_s1})

        return noisy_pops

    def get_circuit_info(self) -> dict[str, Any]:
        """Get quantum circuit resource estimates."""
        n_lindblad = 6 + 4 * 5
        return {
            "n_system_qutrits": self.N,
            "n_ancilla_qubits": n_lindblad,
            "equivalent_qubits": self.N * 2 + n_lindblad,
            "n_lindblad_operators": n_lindblad,
            "encoding": "|S₀⟩=|0⟩, |T₁⟩=|1⟩, |S₁⟩=|2⟩",
            "trotter_order": "2nd (Strang splitting)",
            "stinespring_error": "O(γ² dt²)",
            "advantage": "No forbidden states (natural 3-level encoding)",
        }

    def build_stinespring_unitary_info(self) -> list[dict[str, Any]]:
        """Get detailed Stinespring unitary information for each channel.

        This provides circuit-level information about the Stinespring
        construction for documentation and visualization purposes.
        """
        lindblad_ops = self.build_lindblad_operators()
        info = []

        channel_names = []
        for i, j in self.params.neighbors:
            channel_names.append(f"TTA_1({i},{j})")
            channel_names.append(f"TTA_2({i},{j})")
        for mol in range(self.N):
            channel_names.extend([
                f"Fluorescence({mol})",
                f"Phosphorescence({mol})",
                f"IC({mol})",
                f"ISC_ST({mol})",
                f"ISC_TS({mol})",
            ])

        for idx, (gamma, L) in enumerate(lindblad_ops):
            U_se = _build_stinespring_unitary(
                L, gamma, self.params.dt, self.params.hbar
            )

            # Verify unitarity
            unitarity_error = np.linalg.norm(
                U_se @ U_se.conj().T - np.eye(U_se.shape[0])
            )

            dim = L.shape[0]
            K0 = U_se[:dim, :dim]
            K1 = U_se[dim:, :dim]
            cptp_error = np.linalg.norm(
                K0.conj().T @ K0 + K1.conj().T @ K1 - np.eye(dim)
            )

            info.append({
                "name": channel_names[idx] if idx < len(channel_names) else f"Channel_{idx}",
                "rate": gamma,
                "L_dim": L.shape,
                "U_dim": U_se.shape,
                "theta": np.sqrt(gamma * self.params.dt / self.params.hbar),
                "unitarity_error": unitarity_error,
                "cptp_error": cptp_error,
                "L_nonzero": np.count_nonzero(np.abs(L) > 1e-15),
            })

        return info
