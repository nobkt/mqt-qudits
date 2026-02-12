"""Qubit GKSL-Lindblad simulator using Stinespring dilation.

Implements the GKSL master equation on qubit quantum circuits via Stinespring
dilation of each Lindblad operator. Each dissipative channel is implemented as
a unitary operation on system + ancilla qubits, followed by partial trace
(ancilla measurement and reset).

Qubit encoding (per molecule, 2 qubits):
  |S₀⟩ ↔ |00⟩, |T₁⟩ ↔ |01⟩, |S₁⟩ ↔ |10⟩, |11⟩ = forbidden

System: 4 molecules × 2 qubits = 8 system qubits
Ancilla: 26 qubits (1 per Lindblad operator)
Total: 34 qubits

Reference:
  tutorials/doc/GKSL/GKSL-Lindblad量子ダイナミクスQiskit-Qubit完全実装理論.md
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import scipy.linalg

from gksl_physical_parameters import GKSLPhysicalParameters


def _embed_in_qubit_space(op_3x3: np.ndarray) -> np.ndarray:
    """Embed a 3×3 qutrit operator into 4×4 qubit space.

    Mapping: |0⟩→|00⟩, |1⟩→|01⟩, |2⟩→|10⟩
    The |11⟩ state is forbidden and left as identity.
    """
    op_4x4 = np.eye(4, dtype=complex)
    # Map qutrit indices to qubit indices
    # qutrit |0⟩ = qubit index 0 (|00⟩)
    # qutrit |1⟩ = qubit index 1 (|01⟩)
    # qutrit |2⟩ = qubit index 2 (|10⟩)
    for i in range(3):
        for j in range(3):
            op_4x4[i, j] = op_3x3[i, j]
    return op_4x4


def _single_site_op_qubit(op_3x3: np.ndarray, site: int, N: int) -> np.ndarray:
    """Embed single-site 3×3 operator into full 4^N qubit space."""
    op_4x4 = _embed_in_qubit_space(op_3x3)
    result = np.eye(1, dtype=complex)
    for i in range(N):
        if i == site:
            result = np.kron(result, op_4x4)
        else:
            result = np.kron(result, np.eye(4, dtype=complex))
    return result


def _two_site_op_qubit(
    op_i_3x3: np.ndarray, op_j_3x3: np.ndarray,
    site_i: int, site_j: int, N: int
) -> np.ndarray:
    """Embed two-site tensor product into full 4^N qubit space."""
    op_i = _embed_in_qubit_space(op_i_3x3)
    op_j = _embed_in_qubit_space(op_j_3x3)
    result = np.eye(1, dtype=complex)
    for k in range(N):
        if k == site_i:
            result = np.kron(result, op_i)
        elif k == site_j:
            result = np.kron(result, op_j)
        else:
            result = np.kron(result, np.eye(4, dtype=complex))
    return result


def _build_stinespring_unitary(
    L: np.ndarray, gamma: float, dt: float, hbar: float
) -> np.ndarray:
    """Build Stinespring dilation unitary for a Lindblad operator.

    U = exp(-i θ G) where:
      G = (L ⊗ |1⟩⟨0|_E + L† ⊗ |0⟩⟨1|_E) / √2
      θ = √(2 γ dt / ℏ)

    This implements: E[ρ] ≈ Tr_E[U(ρ ⊗ |0⟩⟨0|_E)U†]
    which reproduces the Lindblad dissipator to O(γ²dt²).

    Args:
        L: Lindblad operator (dim × dim).
        gamma: Dissipation rate (eV/ℏ).
        dt: Time step (fs).
        hbar: ℏ in eV·fs.

    Returns:
        Unitary matrix of dimension 2*dim × 2*dim.
    """
    dim = L.shape[0]
    theta = np.sqrt(gamma * dt / hbar)

    # Build Hermitian generator G
    # G = L ⊗ |1⟩⟨0| + L† ⊗ |0⟩⟨1|
    # In the basis [system⊗|0⟩_E, system⊗|1⟩_E]:
    G = np.zeros((2 * dim, 2 * dim), dtype=complex)
    G[:dim, dim:] = L.conj().T   # L† ⊗ |0⟩⟨1|
    G[dim:, :dim] = L             # L ⊗ |1⟩⟨0|

    U = scipy.linalg.expm(-1j * theta * G)
    return U


class QubitGKSLSimulator:
    """Qubit-based GKSL-Lindblad simulator using Stinespring dilation.

    Uses density matrix simulation in the qubit representation.
    Each Lindblad operator is implemented via Stinespring dilation:
      1. Append ancilla qubit in |0⟩
      2. Apply Stinespring unitary on system + ancilla
      3. Trace out ancilla (partial trace)

    This is a statevector-level simulation that computes the density matrix
    evolution exactly (no shot noise), equivalent to what a quantum computer
    would produce with infinite shots.
    """

    def __init__(self, params: GKSLPhysicalParameters | None = None) -> None:
        self.params = params or GKSLPhysicalParameters()
        self.N = self.params.N_molecules
        self.d = self.params.dim_single  # 3
        self.dim_qubit = 4 ** self.N  # 256 for 4 molecules

    def build_system_hamiltonian_qubit(self) -> np.ndarray:
        """Build system Hamiltonian in qubit encoding (4^N × 4^N)."""
        dim = self.dim_qubit
        H = np.zeros((dim, dim), dtype=complex)

        # On-site energies
        h_single = np.diag([0.0, self.params.E_T, self.params.E_S]).astype(complex)
        for i in range(self.N):
            H += _single_site_op_qubit(h_single, i, self.N)

        # Dexter energy transfer
        ket0_bra1 = np.zeros((3, 3), dtype=complex)
        ket0_bra1[0, 1] = 1.0
        ket1_bra0 = np.zeros((3, 3), dtype=complex)
        ket1_bra0[1, 0] = 1.0

        for i, j in self.params.neighbors:
            H += self.params.V * _two_site_op_qubit(ket0_bra1, ket1_bra0, i, j, self.N)
            H += self.params.V * _two_site_op_qubit(ket1_bra0, ket0_bra1, i, j, self.N)

        return H

    def build_lindblad_operators_qubit(self) -> list[tuple[float, np.ndarray]]:
        """Build all 26 Lindblad operators in qubit encoding."""
        ops = []

        ket0_bra1 = np.zeros((3, 3), dtype=complex)
        ket0_bra1[0, 1] = 1.0
        ket0_bra2 = np.zeros((3, 3), dtype=complex)
        ket0_bra2[0, 2] = 1.0
        ket1_bra2 = np.zeros((3, 3), dtype=complex)
        ket1_bra2[1, 2] = 1.0
        ket2_bra1 = np.zeros((3, 3), dtype=complex)
        ket2_bra1[2, 1] = 1.0

        # TTA operators (6)
        for i, j in self.params.neighbors:
            L1 = _two_site_op_qubit(ket2_bra1, ket0_bra1, i, j, self.N)
            ops.append((self.params.gamma_TTA / 2.0, L1))
            L2 = _two_site_op_qubit(ket0_bra1, ket2_bra1, i, j, self.N)
            ops.append((self.params.gamma_TTA / 2.0, L2))

        # Single-molecule operators (20)
        for i in range(self.N):
            ops.append((self.params.Gamma_fl, _single_site_op_qubit(ket0_bra2, i, self.N)))
            ops.append((self.params.Gamma_ph, _single_site_op_qubit(ket0_bra1, i, self.N)))
            ops.append((self.params.k_IC, _single_site_op_qubit(ket0_bra2, i, self.N)))
            ops.append((self.params.k_ISC_ST, _single_site_op_qubit(ket1_bra2, i, self.N)))
            ops.append((self.params.k_ISC_TS, _single_site_op_qubit(ket0_bra1, i, self.N)))

        return ops

    def _precompute_kraus_operators(
        self, lindblad_ops: list[tuple[float, np.ndarray]]
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        """Precompute Kraus operators from Stinespring unitaries.

        For each Lindblad operator L with rate γ, the Stinespring unitary
        U_SE gives Kraus operators K₀, K₁ via:
          K₀ = ⟨0_E|U_SE|0_E⟩  (no-jump)
          K₁ = ⟨1_E|U_SE|0_E⟩  (jump)

        Returns:
            List of (K0, K1) tuples for each channel.
        """
        kraus_list = []
        for gamma, L in lindblad_ops:
            dim = L.shape[0]
            U_se = _build_stinespring_unitary(L, gamma, self.params.dt, self.params.hbar)
            K0 = U_se[:dim, :dim]   # ⟨0_E|U|0_E⟩
            K1 = U_se[dim:, :dim]   # ⟨1_E|U|0_E⟩
            kraus_list.append((K0, K1))
        return kraus_list

    def _apply_kraus_channel(
        self, rho: np.ndarray, K0: np.ndarray, K1: np.ndarray
    ) -> np.ndarray:
        """Apply Kraus channel: ρ → K₀ ρ K₀† + K₁ ρ K₁†."""
        return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T

    def _apply_hamiltonian_step(
        self, rho: np.ndarray, H: np.ndarray, dt_factor: float
    ) -> np.ndarray:
        """Apply unitary Hamiltonian evolution for time dt_factor * dt.

        ρ → U ρ U† where U = exp(-i H Δt / ℏ)
        """
        U = scipy.linalg.expm(
            -1j * H * dt_factor * self.params.dt / self.params.hbar
        )
        return U @ rho @ U.conj().T

    def _extract_populations_qubit(self, rho: np.ndarray) -> dict[str, float]:
        """Extract populations from qubit-encoded density matrix."""
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0

        for mol in range(self.N):
            for state_idx, key in [(0, "S0"), (1, "T1"), (2, "S1")]:
                proj_3x3 = np.zeros((3, 3), dtype=complex)
                proj_3x3[state_idx, state_idx] = 1.0
                P = _single_site_op_qubit(proj_3x3, mol, self.N)
                pop = np.real(np.trace(P @ rho))
                if key == "S0":
                    N_S0 += pop
                elif key == "T1":
                    N_T1 += pop
                else:
                    N_S1 += pop

        return {"N_S0": N_S0, "N_T1": N_T1, "N_S1": N_S1}

    def _extract_per_molecule_populations_qubit(
        self, rho: np.ndarray
    ) -> dict[str, np.ndarray]:
        """Extract per-molecule populations from qubit-encoded density matrix."""
        S0 = np.zeros(self.N)
        T1 = np.zeros(self.N)
        S1 = np.zeros(self.N)

        for mol in range(self.N):
            for state_idx, arr in [(0, S0), (1, T1), (2, S1)]:
                proj_3x3 = np.zeros((3, 3), dtype=complex)
                proj_3x3[state_idx, state_idx] = 1.0
                P = _single_site_op_qubit(proj_3x3, mol, self.N)
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
        """Run Qubit GKSL simulation via Stinespring dilation.

        Uses 2nd-order symmetric Trotter splitting (Strang splitting):
          exp(L_H dt/2) → exp(L_diss dt) → exp(L_H dt/2)

        Args:
            shots: If provided, add shot noise to final populations.

        Returns:
            Dictionary with simulation results.
        """
        start_time = time.time()

        H = self.build_system_hamiltonian_qubit()
        lindblad_ops = self.build_lindblad_operators_qubit()

        # Build initial density matrix in qubit space
        psi_qutrit = self.params.get_initial_state_vector()
        psi_qubit = np.zeros(self.dim_qubit, dtype=complex)
        for idx in range(self.params.dim_total):
            # Convert qutrit index to qubit index
            states = []
            temp = idx
            for _ in range(self.N):
                states.append(temp % 3)
                temp //= 3
            states.reverse()
            qubit_idx = 0
            for s in states:
                qubit_idx = qubit_idx * 4 + s
            psi_qubit[qubit_idx] = psi_qutrit[idx]

        rho = np.outer(psi_qubit, psi_qubit.conj())

        times = [0.0]
        populations = [self._extract_populations_qubit(rho)]
        per_molecule_populations = [self._extract_per_molecule_populations_qubit(rho)]
        entropy_list = [self._von_neumann_entropy(rho)]
        purity_list = [self._purity(rho)]
        trace_list = [float(np.real(np.trace(rho)))]

        # Precompute half-step Hamiltonian unitary
        U_half = scipy.linalg.expm(
            -1j * H * 0.5 * self.params.dt / self.params.hbar
        )

        # Precompute all Kraus operators (much faster than per-step Stinespring)
        kraus_ops = self._precompute_kraus_operators(lindblad_ops)

        for step in range(self.params.N_steps):
            # Strang splitting: H/2 → Dissipation → H/2

            # 1. Half-step Hamiltonian
            rho = U_half @ rho @ U_half.conj().T

            # 2. Apply all Kraus channels (dissipation)
            for K0, K1 in kraus_ops:
                rho = self._apply_kraus_channel(rho, K0, K1)

            # 3. Half-step Hamiltonian
            rho = U_half @ rho @ U_half.conj().T

            # Record observables
            times.append((step + 1) * self.params.dt)
            populations.append(self._extract_populations_qubit(rho))
            per_molecule_populations.append(
                self._extract_per_molecule_populations_qubit(rho)
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
            "method": "Qubit GKSL (Stinespring)",
            "entropy": entropy_list,
            "purity": purity_list,
            "trace": trace_list,
            "rho_final": rho,
            "n_lindblad_operators": len(lindblad_ops),
            "n_system_qubits": 2 * self.N,
            "n_ancilla_qubits": len(lindblad_ops),
            "n_total_qubits": 2 * self.N + len(lindblad_ops),
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
        """Add shot noise to populations by sampling from multinomial distribution."""
        rng = np.random.default_rng(42)
        noisy_pops = []
        for pop in populations:
            total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
            if total < 1e-10:
                noisy_pops.append(pop.copy())
                continue

            # Per-molecule average probabilities
            probs = np.array([
                pop["N_S0"] / self.N,
                pop["N_T1"] / self.N,
                pop["N_S1"] / self.N,
            ])
            probs = np.clip(probs, 0, None)
            probs /= probs.sum()

            # Sample for each molecule
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
        n_lindblad = 6 + 4 * 5  # 6 TTA + 20 single-molecule
        return {
            "n_system_qubits": 2 * self.N,
            "n_ancilla_qubits": n_lindblad,
            "n_total_qubits": 2 * self.N + n_lindblad,
            "n_lindblad_operators": n_lindblad,
            "encoding": "|S0⟩=|00⟩, |T1⟩=|01⟩, |S1⟩=|10⟩",
            "trotter_order": "2nd (Strang splitting)",
            "stinespring_error": "O(γ² dt²)",
        }
