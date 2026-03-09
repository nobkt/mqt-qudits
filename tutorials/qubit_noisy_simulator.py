#!/usr/bin/env python3
"""
Qubit Molecular Dynamics Simulator with Noise Support.

This module implements a density-matrix-based noisy qubit simulator for the
4-molecule TTA-UC system.  Depolarizing noise is applied **only** to
2-molecule interaction gates (H_transfer and H_TTA), exactly matching the
approach used by NoisyQuditMolecularDynamicsSimulator for qutrits.

Noise model:
- Each molecule is treated as a d=4 "super-qudit" (|S0>=0, |T1>=1, |S1>=2,
  |11>=3 forbidden).  The full system is N d=4 molecules → dimension 4^N.
- H0 (on-site, single-molecule, 4×4 diagonal): **no noise**
- H_transfer (pair 16×16 unitary): depolarizing on the molecule pair
- H_TTA    (pair 16×16 unitary): depolarizing on the molecule pair
- Depolarizing channel: ε(ρ)=(1-p)ρ + p·Tr_{ij}(ρ)⊗I_{16}/16

This gives 12 noise events per Suzuki-Trotter step (same as the qudit case),
which is physically meaningful and avoids the catastrophic noise accumulation
that occurs when Qiskit Aer applies noise to every transpiled cx gate (≈4330
cx gates per step with 1% noise each ≈ 100% noise per step).
"""

import numpy as np
import time
from typing import Dict, List

# Import exact Hamiltonian builders
from exact_qubit_hamiltonians import (
    build_H_transfer_qubit_unitary,
    build_H_TTA_qubit_unitary,
)


class QubitMolecularDynamicsSimulatorNoisy:
    """Qubit-based noisy simulator using density-matrix evolution.

    Each molecule is modelled as a d=4 "super-qudit":
        state 0 → |00⟩ = S0
        state 1 → |01⟩ = T1  (q0=1, q1=0, i.e. q0 is the LSB of each molecule)
        state 2 → |10⟩ = S1
        state 3 → |11⟩ = forbidden

    Total Hilbert-space dimension: 4^N = 256 for N=4 molecules.

    Trotter structure (2nd-order symmetric Suzuki-Trotter):
        H0(dt/2) → H_tr(dt/2) → H_TTA(dt/2) →
        H_TTA(dt/2, rev) → H_tr(dt/2, rev) → H0(dt/2)

    Noise model:
        • H0  (single-molecule 4×4 diagonal): **no noise**
        • H_transfer (molecule-pair 16×16 unitary): 2-molecule depolarizing
        • H_TTA      (molecule-pair 16×16 unitary): 2-molecule depolarizing

    Depolarizing channel on molecule pair (i, j):
        ε(ρ) = (1-p)ρ + p · Tr_{i,j}(ρ) ⊗ I_{16}/16
    """

    def __init__(self, params):
        self.params = params
        self.N = params.N_molecules
        self.d_mol = 4          # dimension per molecule (2 qubits → 4 states)
        self.D = self.d_mol ** self.N   # total dimension (256 for N=4)

        print("Qubit Noisy シミュレータを初期化しました（密度行列方式）")
        print(f"  分子数: {self.N}")
        print(f"  必要Qubit数: {2 * self.N}")
        print(f"  密度行列次元: {self.D} × {self.D}")

    # ------------------------------------------------------------------
    # Unitary builders
    # ------------------------------------------------------------------

    def _build_H0_molecule_unitary(self, dt_half: float) -> np.ndarray:
        """Build the 4×4 single-molecule H0 unitary for a half step dt/2.

        H0 = diag([0, E_T, E_S, 0]) in the molecule basis
        → U_H0 = diag([1, exp(-i E_T dt/(2ℏ)), exp(-i E_S dt/(2ℏ)), 1])

        The forbidden state (index 3) is assigned the same energy as S0 (0).
        """
        E_T = self.params.E_T
        E_S = self.params.E_S
        hbar = self.params.hbar
        return np.diag([
            1.0,
            np.exp(-1j * E_T * dt_half / hbar),
            np.exp(-1j * E_S * dt_half / hbar),
            1.0,
        ])

    def _build_pair_unitary_full(
        self,
        U_pair: np.ndarray,
        mol_i: int,
        mol_j: int,
    ) -> np.ndarray:
        """Embed a 16×16 molecule-pair unitary in the full D×D space.

        Index convention for U_pair (from build_H_*_qubit_unitary):
            index = mol_i_state + 4 * mol_j_state
            (mol_i_state = q_{2i} + 2*q_{2i+1}, same for mol_j)

        Parameters
        ----------
        U_pair : ndarray, shape (16, 16)
        mol_i, mol_j : int
            Molecule indices (mol_i < mol_j assumed by neighbors)
        """
        d = self.d_mol
        N = self.N
        D = self.D
        U_full = np.zeros((D, D), dtype=complex)

        for col in range(D):
            # Decompose col → molecule states
            tmp = col
            states_col = []
            for _ in range(N):
                states_col.append(tmp % d)
                tmp //= d

            s_i_old = states_col[mol_i]
            s_j_old = states_col[mol_j]
            col_pair = s_i_old + d * s_j_old   # index into 16×16

            for row_pair in range(d * d):
                u_elem = U_pair[row_pair, col_pair]
                if abs(u_elem) < 1e-14:
                    continue
                s_i_new = row_pair % d
                s_j_new = row_pair // d

                states_row = list(states_col)
                states_row[mol_i] = s_i_new
                states_row[mol_j] = s_j_new

                row = sum(s * d**k for k, s in enumerate(states_row))
                U_full[row, col] = u_elem

        return U_full

    def _build_single_mol_unitary_full(
        self,
        U_mol: np.ndarray,
        mol_idx: int,
    ) -> np.ndarray:
        """Embed a 4×4 single-molecule unitary in the full D×D space."""
        # U_full = I^{⊗(N-1-mol_idx)} ⊗ U_mol ⊗ I^{⊗mol_idx}
        # (mol_0 is least-significant, so it sits at the right)
        n_high = self.N - 1 - mol_idx   # molecules above mol_idx
        n_low = mol_idx                 # molecules below mol_idx
        left = np.eye(self.d_mol ** n_high, dtype=complex)
        right = np.eye(self.d_mol ** n_low, dtype=complex)
        return np.kron(np.kron(left, U_mol), right)

    # ------------------------------------------------------------------
    # Depolarizing channel
    # ------------------------------------------------------------------

    def _apply_2mol_depolarizing_dm(
        self,
        rho: np.ndarray,
        mol_i: int,
        mol_j: int,
        depol_prob: float,
    ) -> np.ndarray:
        """Apply 2-molecule depolarizing channel on molecule pair (i, j).

        ε(ρ) = (1-p)ρ + p · Tr_{mol_i, mol_j}(ρ) ⊗ I_{d²}/d²

        where d = d_mol = 4 (treating each molecule as a d=4 qudit).

        This is the exact same mathematical structure as
        NoisyQuditMolecularDynamicsSimulator._apply_2qudit_depolarizing_dm,
        evaluated with d=4.
        """
        if depol_prob <= 0.0:
            return rho

        d = self.d_mol   # 4
        N = self.N        # 4
        D = self.D        # 256

        # Reshape to tensor [d]*N row indices + [d]*N col indices
        rho_r = rho.reshape([d] * (2 * N))

        # Step 1: Partial trace over mol_i and mol_j
        row_chars = [chr(ord('a') + k) for k in range(N)]
        col_chars = [chr(ord('a') + N + k) for k in range(N)]

        trace_col_chars = list(col_chars)
        trace_col_chars[mol_i] = row_chars[mol_i]
        trace_col_chars[mol_j] = row_chars[mol_j]

        remaining = [k for k in range(N) if k != mol_i and k != mol_j]

        out_chars = [row_chars[k] for k in remaining] + [col_chars[k] for k in remaining]
        ein_trace = (
            ''.join(row_chars)
            + ''.join(trace_col_chars)
            + '->'
            + ''.join(out_chars)
        )
        rho_rest_r = np.einsum(ein_trace, rho_r)

        # Step 2: Rebuild mixed state rho_rest ⊗ I_d/d ⊗ I_d/d
        I_d = np.eye(d, dtype=complex) / d

        rest_row = [row_chars[k] for k in remaining]
        rest_col = [col_chars[k] for k in remaining]
        Ii_in = row_chars[mol_i] + col_chars[mol_i]
        Ij_in = row_chars[mol_j] + col_chars[mol_j]
        full_out = ''.join(row_chars) + ''.join(col_chars)
        ein_build = (
            ''.join(rest_row)
            + ''.join(rest_col)
            + ','
            + Ii_in
            + ','
            + Ij_in
            + '->'
            + full_out
        )
        mixed_r = np.einsum(ein_build, rho_rest_r, I_d, I_d)
        mixed = mixed_r.reshape(D, D)

        # Step 3: ε(ρ) = (1-p)ρ + p·mixed
        return (1.0 - depol_prob) * rho + depol_prob * mixed

    # ------------------------------------------------------------------
    # Population extraction from density matrix diagonal
    # ------------------------------------------------------------------

    def _populations_from_dm(self, rho: np.ndarray) -> Dict:
        """Compute N_S0, N_T1, N_S1 (and unphysical) from density matrix."""
        d = self.d_mol
        N = self.N
        probs = np.real(np.diag(rho))   # length D

        N_S0 = N_T1 = N_S1 = unphysical = 0.0
        for idx in range(len(probs)):
            p = probs[idx]
            if p <= 0.0:
                continue
            tmp = idx
            for _ in range(N):
                mol_state = tmp % d
                tmp //= d
                if mol_state == 0:
                    N_S0 += p
                elif mol_state == 1:
                    N_T1 += p
                elif mol_state == 2:
                    N_S1 += p
                else:               # state 3 = forbidden
                    unphysical += p
        return {'N_S0': N_S0, 'N_T1': N_T1, 'N_S1': N_S1, 'unphysical': unphysical}

    def _sample_populations(self, rho: np.ndarray, shots: int) -> Dict:
        """Sample from the density-matrix diagonal and return populations."""
        probs = np.real(np.diag(rho))
        probs = np.maximum(probs, 0.0)
        total = probs.sum()
        if total > 0.0:
            probs /= total

        d = self.d_mol
        N = self.N
        samples = np.random.choice(len(probs), size=shots, p=probs)

        N_S0 = N_T1 = N_S1 = unphysical = 0.0
        for s in samples:
            tmp = s
            for _ in range(N):
                mol_state = tmp % d
                tmp //= d
                if mol_state == 0:
                    N_S0 += 1
                elif mol_state == 1:
                    N_T1 += 1
                elif mol_state == 2:
                    N_S1 += 1
                else:
                    unphysical += 1

        scale = 1.0 / shots
        return {
            'N_S0': N_S0 * scale,
            'N_T1': N_T1 * scale,
            'N_S1': N_S1 * scale,
            'unphysical': unphysical * scale,
        }

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def _build_initial_density_matrix(self, state_type: str) -> np.ndarray:
        """Build the initial pure-state density matrix."""
        d = self.d_mol
        N = self.N
        D = self.D

        # T1 = state index 1, S0 = state index 0
        if state_type == 'edge_triplet':
            # mol 0 = T1, mol 1..N-2 = S0, mol N-1 = T1
            mol_states = [0] * N
            mol_states[0] = 1           # T1
            mol_states[N - 1] = 1       # T1
        elif state_type == 'all_triplet':
            mol_states = [1] * N        # all T1
        else:
            mol_states = [0] * N        # all S0

        idx = sum(s * d**k for k, s in enumerate(mol_states))
        psi = np.zeros(D, dtype=complex)
        psi[idx] = 1.0
        return np.outer(psi, psi.conj())

    # ------------------------------------------------------------------
    # Main simulation
    # ------------------------------------------------------------------

    def simulate(
        self,
        T_total: float,
        N_steps: int,
        initial_state_type: str = 'edge_triplet',
        shots: int = 10000,
        noise_params: Dict = None,
    ) -> Dict:
        """Run density-matrix simulation with per-2-molecule-gate depolarizing noise.

        Parameters
        ----------
        T_total : float
            Total simulation time (fs)
        N_steps : int
            Number of Suzuki-Trotter steps
        initial_state_type : str
            'edge_triplet' (default), 'all_triplet'
        shots : int
            Number of samples drawn from the density-matrix diagonal at each
            time point to compute populations (default 10000)
        noise_params : dict or None
            Accepted keys:
              'depol_1q'  – kept for API compatibility, not used
              'depol_2q'  – depolarizing probability per 2-molecule gate (default 0.01)
        """
        print("\n" + "=" * 70)
        print("Qubitベースシミュレーション開始（密度行列ノイズモデル付き）")
        print("=" * 70)
        print(f"ショット数: {shots}")

        if noise_params is None:
            noise_params = {}
        depol_2q = noise_params.get('depol_2q', 0.01)

        print("\nノイズモデルパラメータ:")
        print("  H0（単一分子ゲート）: 理想的（ノイズなし）")
        print(f"  H_transfer / H_TTA（2-分子ゲート）脱分極エラー: {depol_2q * 100:.3f}%")
        print("  その他のノイズ: なし（脱分極のみ）")

        dt = T_total / N_steps
        dt_half = dt / 2.0

        # ---------------------------------------------------------------
        # Pre-build unitaries
        # ---------------------------------------------------------------
        print("\n全ユニタリ行列を事前構築中...")

        # H0 (single-molecule, d=4, no noise)
        U_H0_mol = self._build_H0_molecule_unitary(dt_half)
        U_H0_full_list = [
            self._build_single_mol_unitary_full(U_H0_mol, mol)
            for mol in range(self.N)
        ]

        # H_transfer per neighbor pair
        U_transfer_half_list = []
        for mol_i, mol_j in self.params.neighbors:
            U16 = build_H_transfer_qubit_unitary(
                self.params.V, dt_half, self.params.hbar
            )
            U_full = self._build_pair_unitary_full(U16, mol_i, mol_j)
            U_transfer_half_list.append((U_full, (mol_i, mol_j)))

        # H_TTA per neighbor pair
        U_TTA_half_list = []
        for mol_i, mol_j in self.params.neighbors:
            U16 = build_H_TTA_qubit_unitary(
                self.params.J, dt_half, self.params.hbar
            )
            U_full = self._build_pair_unitary_full(U16, mol_i, mol_j)
            U_TTA_half_list.append((U_full, (mol_i, mol_j)))

        n_pair_gates_per_step = 2 * (len(U_transfer_half_list) + len(U_TTA_half_list))
        effective_noise = 1.0 - (1.0 - depol_2q) ** n_pair_gates_per_step
        print(f"  2-分子ゲート数/ステップ: {n_pair_gates_per_step}")
        print(f"  有効ステップノイズ: 1-(1-{depol_2q})^{n_pair_gates_per_step} = {effective_noise:.4f}")

        # ---------------------------------------------------------------
        # Initial density matrix
        # ---------------------------------------------------------------
        rho = self._build_initial_density_matrix(initial_state_type)
        pop_0 = self._sample_populations(rho, shots)

        print(f"\n初期状態: {initial_state_type}")
        print(f"  N_S0 = {pop_0['N_S0']:.4f}")
        print(f"  N_T1 = {pop_0['N_T1']:.4f}")
        print(f"  N_S1 = {pop_0['N_S1']:.4f}")

        times = [0.0]
        populations = [pop_0]

        start_time = time.time()

        # ---------------------------------------------------------------
        # Time evolution
        # ---------------------------------------------------------------
        print(f"\n時間発展を実行中（{N_steps}ステップ、密度行列シミュレーション）...")
        print(f"  ノイズ: 各2-分子ゲートに{depol_2q * 100:.1f}%脱分極エラーを適用")

        for step in range(1, N_steps + 1):
            # --- Forward half: H0 → H_transfer → H_TTA ---

            # H0 (no noise)
            for U_H0 in U_H0_full_list:
                rho = U_H0 @ rho @ U_H0.conj().T

            # H_transfer (2-molecule gates, with noise)
            for U_tr, (mi, mj) in U_transfer_half_list:
                rho = U_tr @ rho @ U_tr.conj().T
                rho = self._apply_2mol_depolarizing_dm(rho, mi, mj, depol_2q)

            # H_TTA (2-molecule gates, with noise)
            for U_TTA, (mi, mj) in U_TTA_half_list:
                rho = U_TTA @ rho @ U_TTA.conj().T
                rho = self._apply_2mol_depolarizing_dm(rho, mi, mj, depol_2q)

            # --- Backward half: H_TTA(rev) → H_transfer(rev) → H0 ---

            for U_TTA, (mi, mj) in reversed(U_TTA_half_list):
                rho = U_TTA @ rho @ U_TTA.conj().T
                rho = self._apply_2mol_depolarizing_dm(rho, mi, mj, depol_2q)

            for U_tr, (mi, mj) in reversed(U_transfer_half_list):
                rho = U_tr @ rho @ U_tr.conj().T
                rho = self._apply_2mol_depolarizing_dm(rho, mi, mj, depol_2q)

            # H0 (no noise)
            for U_H0 in reversed(U_H0_full_list):
                rho = U_H0 @ rho @ U_H0.conj().T

            # Enforce Hermiticity and normalization (numerical stability)
            rho = (rho + rho.conj().T) * 0.5
            tr = np.real(np.trace(rho))
            if tr > 0.0:
                rho /= tr

            pop = self._sample_populations(rho, shots)
            t = step * dt
            times.append(t)
            populations.append(pop)

            if step % max(1, N_steps // 10) == 0:
                purity = np.real(np.trace(rho @ rho))
                print(
                    f"  ステップ {step}/{N_steps}: t = {t:.2f} fs, "
                    f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}, "
                    f"非物理的: {pop['unphysical']:.4f}, 純度: {purity:.4f}"
                )

        elapsed = time.time() - start_time
        final_pop = populations[-1]
        final_purity = float(np.real(np.trace(rho @ rho)))

        print("\n" + "=" * 70)
        print("シミュレーション完了")
        print("=" * 70)
        print("最終個体数:")
        print(f"  N_S0 = {final_pop['N_S0']:.4f}")
        print(f"  N_T1 = {final_pop['N_T1']:.4f}")
        print(f"  N_S1 = {final_pop['N_S1']:.4f}")
        print(f"  非物理的状態（禁止状態）: {final_pop['unphysical']:.6f}")
        print(f"  最終純度: {final_purity:.6f}")
        print(f"\n回路統計:")
        print(f"  2-分子ゲート数/ステップ: {n_pair_gates_per_step}")
        print(f"  総2-分子ゲート数: {n_pair_gates_per_step * N_steps}")
        print(f"  実行時間: {elapsed:.2f}秒")

        return {
            'times': times,
            'populations': populations,
            'elapsed_time': elapsed,
            'method': f'Qubit (Noisy DM, depol_2q={depol_2q:.4f})',
            'n_pair_gates_per_step': n_pair_gates_per_step,
            'shots': shots,
            'noise_params': {'depol_2q': depol_2q},
        }
