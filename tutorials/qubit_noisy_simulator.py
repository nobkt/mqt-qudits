#!/usr/bin/env python3
"""
Qubit Molecular Dynamics Simulator with Noise Support.

This module extends the exact Qubit simulator to support noise models
using density matrix simulation with pair-level depolarizing noise.

The noise model matches the qudit noisy simulator (mqt_qudits_noisy_simulator.py):
- Depolarizing noise applied at the 2-molecule pair-gate level
- 12 noise events per 2nd-order symmetric Trotter step
- H0 single-molecule gates receive no noise

Depolarizing channel for a pair of molecules (each d=4, 2 qubits):
    ε(ρ) = (1-p)ρ + p · Tr_{ij}(ρ) ⊗ I_{16}/16

Note: Since d=4 includes the forbidden state |11⟩, the depolarizing channel
inherently causes forbidden-state leakage. This is a physical property of
the qubit encoding, not a bug.
"""

import numpy as np
import time
from typing import Dict

from exact_qubit_hamiltonians import (
    build_H_transfer_qubit_unitary,
    build_H_TTA_qubit_unitary,
)


class QubitMolecularDynamicsSimulatorNoisy:
    """Qubit-based simulator with density matrix pair-level noise.

    Uses the same noise granularity as the qudit noisy simulator:
    depolarizing noise applied after each 2-molecule pair gate
    (H_transfer and H_TTA), with no noise on H0 single-molecule gates.
    """

    def __init__(self, params):
        self.params = params
        self.N = params.N_molecules
        self.d = 4  # 2 qubits per molecule -> d=4 local dimension
        self.dim = self.d ** self.N  # 4^N (256 for 4 molecules)
        self.n_qubits = 2 * self.N

        print(f"Qubit Noisy シミュレータを初期化しました")
        print(f"  分子数: {self.N}")
        print(f"  必要Qubit数: {self.n_qubits}")

    def _build_per_pair_unitaries(self, dt: float):
        """Build per-molecule and per-pair unitary matrices for the Trotter step.

        Returns:
            Tuple of (U_H0_half_list, U_transfer_half_list, U_TTA_half_list)
        """
        I4 = np.eye(self.d, dtype=complex)

        def build_single_molecule_operator(mol_idx, op_4x4):
            """Build dim x dim operator from a 4x4 single-molecule operator."""
            operators = [I4] * self.N
            operators[mol_idx] = op_4x4
            result = operators[0]
            for i in range(1, self.N):
                result = np.kron(result, operators[i])
            return result

        def build_two_molecule_operator(mol_i, mol_j, op_16x16):
            """Build dim x dim operator from a 16x16 pair operator."""
            if mol_i == 0 and mol_j == 1:
                result = op_16x16
                for k in range(2, self.N):
                    result = np.kron(result, I4)
            elif mol_i == 1 and mol_j == 2:
                result = np.kron(I4, op_16x16)
                for k in range(3, self.N):
                    result = np.kron(result, I4)
            elif mol_i == 2 and mol_j == 3:
                result = np.kron(np.kron(I4, I4), op_16x16)
            else:
                msg = f"Unsupported molecule pair: ({mol_i}, {mol_j})"
                raise ValueError(msg)
            return result

        # H0 per-molecule unitaries (dt/2)
        # Qubit encoding: |00>=S0(E=0), |01>=T1(E=E_T), |10>=S1(E=E_S), |11>=forbidden(E=0)
        U_H0_half_list = []
        for mol_idx in range(self.N):
            U_H0_mol = np.diag([
                1.0,
                np.exp(-1j * self.params.E_T * dt / (2 * self.params.hbar)),
                np.exp(-1j * self.params.E_S * dt / (2 * self.params.hbar)),
                1.0,  # forbidden state: no evolution
            ])
            U_mol_full = build_single_molecule_operator(mol_idx, U_H0_mol)
            U_H0_half_list.append(U_mol_full)

        # H_transfer per-pair unitaries (dt/2) -- 16x16 from exact_qubit_hamiltonians
        U_transfer_half_list = []
        for pair_idx, (mol_i, mol_j) in enumerate(self.params.neighbors):
            V = self.params.V[pair_idx] if isinstance(self.params.V, (list, np.ndarray)) else self.params.V
            U_16x16 = build_H_transfer_qubit_unitary(V, dt / 2, self.params.hbar)
            U_full = build_two_molecule_operator(mol_i, mol_j, U_16x16)
            U_transfer_half_list.append((U_full, (mol_i, mol_j)))

        # H_TTA per-pair unitaries (dt/2) -- 16x16 from exact_qubit_hamiltonians
        U_TTA_half_list = []
        for pair_idx, (mol_i, mol_j) in enumerate(self.params.neighbors):
            J = self.params.J[pair_idx] if isinstance(self.params.J, (list, np.ndarray)) else self.params.J
            U_16x16 = build_H_TTA_qubit_unitary(J, dt / 2, self.params.hbar)
            U_full = build_two_molecule_operator(mol_i, mol_j, U_16x16)
            U_TTA_half_list.append((U_full, (mol_i, mol_j)))

        return U_H0_half_list, U_transfer_half_list, U_TTA_half_list

    def _apply_2qubit_pair_depolarizing_dm(self, rho: np.ndarray,
                                            mol_i: int, mol_j: int,
                                            depol_prob: float) -> np.ndarray:
        """Apply 2-molecule pair depolarizing noise to density matrix.

        Implements the exact quantum channel:
            E(rho) = (1-p)*rho + p * Tr_{ij}(rho) (x) I_{ij}/d^2

        where d = 4 (2-qubit molecule dimension), d^2 = 16 (pair dimension).
        """
        if depol_prob <= 0:
            return rho

        d = self.d  # 4
        N = self.N
        D = self.dim  # 4^N

        # Reshape density matrix to tensor form [d]*2N
        rho_r = rho.reshape([d] * (2 * N))

        # Step 1: Partial trace over molecules i and j using einsum
        row_chars = [chr(ord('a') + k) for k in range(N)]
        col_chars = [chr(ord('a') + N + k) for k in range(N)]

        trace_col_chars = list(col_chars)
        trace_col_chars[mol_i] = row_chars[mol_i]
        trace_col_chars[mol_j] = row_chars[mol_j]

        remaining = [k for k in range(N) if k != mol_i and k != mol_j]

        output_chars = []
        for k in remaining:
            output_chars.append(row_chars[k])
        for k in remaining:
            output_chars.append(col_chars[k])

        einsum_trace = ''.join(row_chars) + ''.join(trace_col_chars) + '->' + ''.join(output_chars)
        rho_rest_r = np.einsum(einsum_trace, rho_r)

        # Step 2: Build mixed state: rho_rest (x) I_pair/d^2
        I_d = np.eye(d, dtype=complex) / d

        rest_row_chars = [row_chars[k] for k in remaining]
        rest_col_chars = [col_chars[k] for k in remaining]
        rest_input = ''.join(rest_row_chars) + ''.join(rest_col_chars)

        Ii_input = row_chars[mol_i] + col_chars[mol_i]
        Ij_input = row_chars[mol_j] + col_chars[mol_j]

        full_output = ''.join(row_chars) + ''.join(col_chars)

        einsum_build = f'{rest_input},{Ii_input},{Ij_input}->{full_output}'
        mixed_r = np.einsum(einsum_build, rho_rest_r, I_d, I_d)
        mixed = mixed_r.reshape(D, D)

        # Step 3: Apply depolarizing channel
        rho_new = (1 - depol_prob) * rho + depol_prob * mixed

        return rho_new

    def _calculate_populations_from_dm(self, rho: np.ndarray) -> Dict[str, float]:
        """Calculate populations directly from density matrix diagonal."""
        N_S0 = N_T1 = N_S1 = 0.0
        unphysical = 0.0

        probs = np.diag(rho).real

        for idx in range(self.dim):
            prob = probs[idx]
            if prob < 1e-15:
                continue

            is_forbidden = False
            mol_S0 = mol_T1 = mol_S1 = 0
            temp = idx
            for _mol in range(self.N):
                mol_state = temp % self.d
                temp //= self.d
                if mol_state == 3:  # |11> = forbidden
                    is_forbidden = True
                    break
                if mol_state == 0:  # |00> = S0
                    mol_S0 += 1
                elif mol_state == 1:  # |01> = T1
                    mol_T1 += 1
                elif mol_state == 2:  # |10> = S1
                    mol_S1 += 1

            if is_forbidden:
                unphysical += prob
            else:
                N_S0 += prob * mol_S0
                N_T1 += prob * mol_T1
                N_S1 += prob * mol_S1

        return {
            'N_S0': N_S0,
            'N_T1': N_T1,
            'N_S1': N_S1,
            'unphysical': unphysical
        }

    def _build_initial_statevector(self, initial_state_type: str) -> np.ndarray:
        """Build initial statevector in the 4^N dimensional space."""
        psi = np.zeros(self.dim, dtype=complex)

        if initial_state_type == 'edge_triplet':
            # Molecule 0 = T1 (|01> = index 1), Molecule N-1 = T1 (|01> = index 1)
            # Others = S0 (|00> = index 0)
            mol_states = [0] * self.N  # all S0
            mol_states[0] = 1  # T1
            mol_states[self.N - 1] = 1  # T1
            idx = 0
            for i in range(self.N):
                idx += mol_states[i] * (self.d ** i)
            psi[idx] = 1.0
        elif initial_state_type == 'all_triplet':
            mol_states = [1] * self.N  # all T1
            idx = 0
            for i in range(self.N):
                idx += mol_states[i] * (self.d ** i)
            psi[idx] = 1.0
        else:
            msg = f"Unknown initial state type: {initial_state_type}"
            raise ValueError(msg)

        return psi

    def _calculate_populations_from_samples(self, samples: np.ndarray,
                                             shots: int) -> Dict[str, float]:
        """Calculate populations from Monte Carlo samples."""
        N_S0 = N_T1 = N_S1 = 0.0
        unphysical = 0.0

        for idx in samples:
            is_forbidden = False
            mol_S0 = mol_T1 = mol_S1 = 0
            temp = int(idx)
            for _mol in range(self.N):
                mol_state = temp % self.d
                temp //= self.d
                if mol_state == 3:
                    is_forbidden = True
                    break
                if mol_state == 0:
                    mol_S0 += 1
                elif mol_state == 1:
                    mol_T1 += 1
                elif mol_state == 2:
                    mol_S1 += 1

            if is_forbidden:
                unphysical += 1.0 / shots
            else:
                N_S0 += mol_S0 / shots
                N_T1 += mol_T1 / shots
                N_S1 += mol_S1 / shots

        return {
            'N_S0': N_S0,
            'N_T1': N_T1,
            'N_S1': N_S1,
            'unphysical': unphysical
        }

    def simulate(self, T_total: float, N_steps: int,
                 initial_state_type: str = 'edge_triplet',
                 shots: int = 10000,
                 noise_params: Dict = None) -> Dict:
        """Run noisy simulation using density matrix with pair-level noise.

        Uses density matrix formalism with per-pair-gate 2-molecule depolarizing
        noise. Each 2-molecule pair interaction (H_transfer, H_TTA) receives
        independent depolarizing noise, matching the per-gate noise model used
        in the qudit noisy simulator.

        Parameters
        ----------
        T_total : float
            Total simulation time in fs
        N_steps : int
            Number of Trotter steps
        initial_state_type : str
            Initial state configuration
        shots : int
            Number of measurement shots per time step (for Monte Carlo sampling)
        noise_params : dict or None
            Noise parameters. Keys: 'depol_1q' (unused), 'depol_2q'
        """
        print("\n" + "="*70)
        print("Qubitベースシミュレーション開始（ノイズモデル付き）")
        print("="*70)
        print(f"ショット数: {shots}")

        dt = T_total / N_steps

        if noise_params is None:
            noise_params = {}

        depol_1q = noise_params.get('depol_1q', 0.001)
        depol_2q = noise_params.get('depol_2q', 0.01)

        print("\nノイズモデルパラメータ:")
        print(f"  1量子ビットゲート: 理想的（ノイズなし）")
        print(f"  2量子ビットゲート脱分極エラー: {depol_2q*100:.3f}%")
        print(f"  ノイズ適用: ペアゲートレベル（H_transfer, H_TTA）")
        print(f"  熱緩和: 無効")

        start_time = time.time()

        # Build per-pair unitaries
        print("\nBuilding per-pair unitaries for density matrix simulation...")
        U_H0_half_list, U_transfer_half_list, U_TTA_half_list = \
            self._build_per_pair_unitaries(dt)

        n_2qubit_gates = 2 * (len(U_transfer_half_list) + len(U_TTA_half_list))
        print(f"  2-molecule pair gates per Trotter step: {n_2qubit_gates}")
        print(f"  Effective per-step noise: 1-(1-{depol_2q})^{n_2qubit_gates} = "
              f"{1-(1-depol_2q)**n_2qubit_gates:.4f}")

        # Initial state
        psi_0 = self._build_initial_statevector(initial_state_type)
        rho = np.outer(psi_0, psi_0.conj())

        pop_0 = self._calculate_populations_from_dm(rho)

        print(f"\n初期状態: {initial_state_type}")
        print(f"  N_S0 = {pop_0['N_S0']:.4f}")
        print(f"  N_T1 = {pop_0['N_T1']:.4f}")
        print(f"  N_S1 = {pop_0['N_S1']:.4f}")

        times = [0.0]
        populations = [pop_0]

        # Time evolution
        print(f"\n時間発展を実行中（{N_steps}ステップ、密度行列シミュレーション）...")
        print(f"  ノイズ: 各ペアゲートに{depol_2q*100:.3f}%脱分極エラーを適用")

        for step in range(1, N_steps + 1):
            # Forward half: H0 -> H_transfer -> H_TTA
            for U_H0 in U_H0_half_list:
                rho = U_H0 @ rho @ U_H0.conj().T

            for U_tr, pair in U_transfer_half_list:
                rho = U_tr @ rho @ U_tr.conj().T
                rho = self._apply_2qubit_pair_depolarizing_dm(rho, pair[0], pair[1], depol_2q)

            for U_TTA, pair in U_TTA_half_list:
                rho = U_TTA @ rho @ U_TTA.conj().T
                rho = self._apply_2qubit_pair_depolarizing_dm(rho, pair[0], pair[1], depol_2q)

            # Backward half: H_TTA -> H_transfer -> H0 (reverse)
            for U_TTA, pair in reversed(U_TTA_half_list):
                rho = U_TTA @ rho @ U_TTA.conj().T
                rho = self._apply_2qubit_pair_depolarizing_dm(rho, pair[0], pair[1], depol_2q)

            for U_tr, pair in reversed(U_transfer_half_list):
                rho = U_tr @ rho @ U_tr.conj().T
                rho = self._apply_2qubit_pair_depolarizing_dm(rho, pair[0], pair[1], depol_2q)

            for U_H0 in reversed(U_H0_half_list):
                rho = U_H0 @ rho @ U_H0.conj().T

            # Numerical stability: ensure Hermitian and normalized
            rho = (rho + rho.conj().T) / 2
            trace_val = np.trace(rho).real
            if trace_val > 0:
                rho /= trace_val

            # Calculate populations from Monte Carlo sampling
            probabilities = np.diag(rho).real
            probabilities = np.maximum(probabilities, 0)
            probabilities /= np.sum(probabilities)
            samples = np.random.choice(self.dim, size=shots, p=probabilities)

            pop = self._calculate_populations_from_samples(samples, shots)

            t = step * dt
            times.append(t)
            populations.append(pop)

            if step % max(1, N_steps // 10) == 0:
                print(f"  ステップ {step}/{N_steps}: t = {t:.2f} fs, "
                      f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}, "
                      f"非物理的: {pop['unphysical']:.4f}")

        elapsed = time.time() - start_time

        print("\n" + "="*70)
        print("シミュレーション完了")
        print("="*70)
        print(f"最終個体数:")
        print(f"  N_S0 = {populations[-1]['N_S0']:.4f}")
        print(f"  N_T1 = {populations[-1]['N_T1']:.4f}")
        print(f"  N_S1 = {populations[-1]['N_S1']:.4f}")
        print(f"  非物理的状態: {populations[-1]['unphysical']:.6f}")
        print(f"\n回路統計:")
        print(f"  ペアゲート数/ステップ: {n_2qubit_gates}")
        print(f"  総ペアゲート数: {n_2qubit_gates * N_steps}")
        print(f"  実行時間: {elapsed:.2f}秒")

        return {
            'times': times,
            'populations': populations,
            'elapsed_time': elapsed,
            'method': f'Qubit (Noisy, depol_2q={depol_2q:.4f})',
            'n_2qubit_gates_per_step': n_2qubit_gates,
            'total_gates': n_2qubit_gates * N_steps,
            'shots': shots,
            'noise_params': noise_params
        }
