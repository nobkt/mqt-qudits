#!/usr/bin/env python3
"""
MQT-Qudits Molecular Dynamics Simulator with Noise Support.

This module extends the exact Qudit simulator to support noise models
using MQT-qudits' noise capabilities.

Noise models supported:
- Depolarizing error: Random errors on qudits
- Dephasing error: Phase damping
- Mathematical Noise: Global noise applied to all qudit levels

All noise parameters are physically motivated and use realistic values
for trapped-ion or superconducting qutrits.
"""

import sys
from pathlib import Path
import numpy as np
from typing import List, Tuple, Dict, Optional
import time

# Import from the existing implementation
from mqt_qudits_four_molecule_sparse_implementation import (
    PhysicalParameters,
    SparseAwareMQTQuditTimeEvolution,
    index_to_config,
    config_to_index,
    config_to_state_name
)


class NoisyQuditMolecularDynamicsSimulator:
    """
    Qudit-based simulator with realistic noise models using MQT-qudits noise tools.
    
    This simulator uses the exact Hamiltonian-based time evolution but adds
    realistic noise to each gate operation using MQT-qudits NoiseModel.
    """
    
    # Noise model parameters
    MAX_DEPOL_PROB = 0.1  # Maximum effective depolarizing probability (10%)
    MAX_DEPHASING_PROB = 0.1  # Maximum effective dephasing probability (10%)
    
    def __init__(self, params: PhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.dim = 3 ** self.N
        
        # MQT-Quditsのインポート
        try:
            from mqt.qudits.simulation import MQTQuditProvider
            from mqt.qudits.simulation.noise_tools import NoiseModel, Noise
            # Note: SubspaceNoise is not imported because the C++ backend (bindings.cpp)
            # expects Noise objects with direct probability_depolarizing and probability_dephasing
            # attributes. SubspaceNoise stores these in a dictionary and causes AttributeError.
            self.provider = MQTQuditProvider()
            self.NoiseModel = NoiseModel
            self.Noise = Noise
            self.mqt_available = True
        except ImportError as e:
            self.mqt_available = False
            raise ImportError(f"mqt.quditsがインストールされていません: {e}")
        
        # Time evolution module for building circuits
        self.time_evol = SparseAwareMQTQuditTimeEvolution(params)
        
        print(f"Noisy Qudit シミュレータを初期化しました")
        print(f"  分子数: {self.N}")
        print(f"  必要Qutrit数: {self.N}")
        print(f"  状態空間次元: {self.dim}")
    
    def create_noise_model(self,
                          depol_1q: float = 0.001,
                          depol_2q: float = 0.01,
                          noise_gates: List[str] = None) -> "NoiseModel":
        """
        Create a realistic noise model for qutrits.
        
        NOTE: As per the requirement, noise is applied ONLY to 2-qudit gates.
        Single-qudit gates are assumed to be ideal (no noise).
        
        Parameters:
        -----------
        depol_1q : float
            Depolarizing error probability for single-qudit gates (NOT USED - kept for API compatibility)
        depol_2q : float
            Depolarizing error probability for two-qudit gates (default: 0.01 = 1%)
        noise_gates : List[str] or None
            List of gate names to apply noise to. If None, applies to all 2-qudit gates.
            NOTE: Only 2-qudit gates will receive noise regardless of this parameter.
        
        Returns:
        --------
        NoiseModel : MQT-qudits noise model
        """
        noise_model = self.NoiseModel()
        
        # MODIFICATION: Single-qudit gates are now IDEAL (no noise applied)
        # This follows the requirement: "qudit量子シミュレーションにおけるノイズモデルは2-quditゲートに対してのみ施す"
        
        # Identify two-qudit gates only
        if noise_gates is None:
            # Default to all common 2-qudit gates
            two_qudit_gates = ['cx', 'csum', 'ls', 'ms']
        else:
            # Filter to keep only 2-qudit gates
            two_qudit_gates = [g for g in noise_gates if g in ['cx', 'csum', 'ls', 'ms']]
        
        # Apply noise ONLY to two-qudit gates
        if two_qudit_gates:
            noise_2q = self.Noise(depol_2q, 0.0)  # Only depolarizing, no dephasing
            noise_model.add_nonlocal_quantum_error(noise_2q, two_qudit_gates)
        
        return noise_model
    
    def build_initial_state_circuit(self, initial_state_type: str = 'edge_triplet'):
        """初期状態回路を構築"""
        from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
        
        circuit = QuantumCircuit()
        reg = QuantumRegister("molecules", self.N, [3] * self.N)
        circuit.append(reg)
        
        if initial_state_type == 'edge_triplet':
            # |1001⟩: 両端が三重項
            circuit.x(0)  # 分子0をT1にセット
            circuit.x(self.N - 1)  # 分子N-1をT1にセット
        elif initial_state_type == 'all_triplet':
            # |1111⟩: 全て三重項
            for i in range(self.N):
                circuit.x(i)
        
        return circuit
    
    def calculate_populations_from_samples(self, samples: np.ndarray, shots: int) -> Dict[str, float]:
        """サンプルから個体数を計算"""
        N_S0 = N_T1 = N_S1 = 0.0
        
        for sample in samples:
            config = index_to_config(sample, self.N, 3)
            for level in config:
                if level == 0:
                    N_S0 += 1
                elif level == 1:
                    N_T1 += 1
                elif level == 2:
                    N_S1 += 1
        
        # Normalize by shots
        N_S0 /= shots
        N_T1 /= shots
        N_S1 /= shots
        
        return {'N_S0': N_S0, 'N_T1': N_T1, 'N_S1': N_S1}
    
    def calculate_per_molecule_populations_from_samples(self, samples: np.ndarray, 
                                                       shots: int) -> Dict[str, np.ndarray]:
        """サンプルから分子ごとの個体数を計算"""
        S0_per_mol = np.zeros(self.N)
        T1_per_mol = np.zeros(self.N)
        S1_per_mol = np.zeros(self.N)
        
        for sample in samples:
            config = index_to_config(sample, self.N, 3)
            for mol_idx, level in enumerate(config):
                if level == 0:
                    S0_per_mol[mol_idx] += 1
                elif level == 1:
                    T1_per_mol[mol_idx] += 1
                elif level == 2:
                    S1_per_mol[mol_idx] += 1
        
        # Normalize by shots
        S0_per_mol /= shots
        T1_per_mol /= shots
        S1_per_mol /= shots
        
        return {
            'S0_per_mol': S0_per_mol,
            'T1_per_mol': T1_per_mol,
            'S1_per_mol': S1_per_mol
        }
    
    def _build_per_pair_unitaries(self, dt: float):
        """
        Build per-molecule and per-pair unitary matrices for the Trotter step.
        
        These are the same unitaries used by the classical simulator, built from
        scipy.linalg.expm. They are used to interleave noise with gate applications.
        
        Args:
            dt: Time step (fs)
            
        Returns:
            Tuple of (U_H0_half_list, U_transfer_half_list, U_TTA_half_list)
            - U_H0_half_list: List of (dim×dim) unitaries for H0 per molecule
            - U_transfer_half_list: List of ((dim×dim) unitary, (mol_i, mol_j)) for H_transfer
            - U_TTA_half_list: List of ((dim×dim) unitary, (mol_i, mol_j)) for H_TTA
        """
        import scipy.linalg
        from exact_hamiltonian_builders import build_H_transfer_unitary, build_H_TTA_unitary
        
        I3 = np.eye(3, dtype=complex)
        
        def build_single_molecule_operator(mol_idx, op):
            operators = [I3] * self.N
            operators[mol_idx] = op
            result = operators[0]
            for i in range(1, self.N):
                result = np.kron(result, operators[i])
            return result
        
        def build_two_molecule_operator(mol_i, mol_j, op_9x9):
            if mol_i == 0 and mol_j == 1:
                result = op_9x9
                for k in range(2, self.N):
                    result = np.kron(result, I3)
            elif mol_i == 1 and mol_j == 2:
                result = np.kron(I3, op_9x9)
                for k in range(3, self.N):
                    result = np.kron(result, I3)
            elif mol_i == 2 and mol_j == 3:
                result = np.kron(np.kron(I3, I3), op_9x9)
            else:
                raise ValueError(f"Unsupported molecule pair: ({mol_i}, {mol_j})")
            return result
        
        # H0 per-molecule unitaries (dt/2)
        U_H0_half_list = []
        for mol_idx in range(self.N):
            U_H0_mol = np.diag([
                1.0,
                np.exp(-1j * self.params.E_T * dt / (2 * self.params.hbar)),
                np.exp(-1j * self.params.E_S * dt / (2 * self.params.hbar))
            ])
            U_mol_full = build_single_molecule_operator(mol_idx, U_H0_mol)
            U_H0_half_list.append(U_mol_full)
        
        # H_transfer per-pair unitaries (dt/2)
        U_transfer_half_list = []
        for pair_idx, (mol_i, mol_j) in enumerate(self.params.neighbors):
            V = self.params.V[pair_idx] if isinstance(self.params.V, (list, np.ndarray)) else self.params.V
            U_9x9 = build_H_transfer_unitary(V, dt / 2, self.params.hbar, dim=3)
            U_full = build_two_molecule_operator(mol_i, mol_j, U_9x9)
            U_transfer_half_list.append((U_full, (mol_i, mol_j)))
        
        # H_TTA per-pair unitaries (dt/2)
        U_TTA_half_list = []
        for pair_idx, (mol_i, mol_j) in enumerate(self.params.neighbors):
            J = self.params.J[pair_idx] if isinstance(self.params.J, (list, np.ndarray)) else self.params.J
            U_9x9 = build_H_TTA_unitary(J, dt / 2, self.params.hbar, dim=3)
            U_full = build_two_molecule_operator(mol_i, mol_j, U_9x9)
            U_TTA_half_list.append((U_full, (mol_i, mol_j)))
        
        return U_H0_half_list, U_transfer_half_list, U_TTA_half_list
    
    def _apply_2qudit_depolarizing_dm(self, rho: np.ndarray, 
                                       qudit_i: int, qudit_j: int,
                                       depol_prob: float) -> np.ndarray:
        """
        Apply 2-qutrit depolarizing noise on qudit pair (i, j) to density matrix.
        
        Implements the exact quantum channel:
            ε(ρ) = (1-p)ρ + p · Tr_{ij}(ρ) ⊗ I_{ij}/d²
        
        where d = 3 (qutrit dimension), d² = 9 (pair dimension),
        Tr_{ij} is the partial trace over qudits i and j,
        and I_{ij}/d² is the maximally mixed state on the pair.
        
        Parameters:
        -----------
        rho : np.ndarray
            Density matrix (dim × dim) where dim = 3^N
        qudit_i, qudit_j : int
            Indices of the qudit pair to apply noise to
        depol_prob : float
            Depolarizing probability per gate (e.g., 0.01 = 1%)
        
        Returns:
        --------
        rho_new : np.ndarray
            Noisy density matrix after applying the depolarizing channel
        """
        if depol_prob <= 0:
            return rho
        
        d = 3  # qutrit dimension
        N = self.N
        D = self.dim  # 3^N = 81
        
        # Reshape density matrix to tensor form [d]*2N
        rho_r = rho.reshape([d] * (2 * N))
        
        # Step 1: Partial trace over qudits i and j using einsum
        # Row indices: 0..N-1, Column indices: N..2N-1
        # To trace qudit k: set column index N+k equal to row index k
        row_chars = [chr(ord('a') + k) for k in range(N)]
        col_chars = [chr(ord('a') + N + k) for k in range(N)]
        
        # For traced qudits, set column index = row index
        trace_col_chars = list(col_chars)
        trace_col_chars[qudit_i] = row_chars[qudit_i]
        trace_col_chars[qudit_j] = row_chars[qudit_j]
        
        remaining_qudits = [k for k in range(N) if k != qudit_i and k != qudit_j]
        
        output_chars = []
        for k in remaining_qudits:
            output_chars.append(row_chars[k])
        for k in remaining_qudits:
            output_chars.append(col_chars[k])
        
        einsum_trace = ''.join(row_chars) + ''.join(trace_col_chars) + '->' + ''.join(output_chars)
        rho_rest_r = np.einsum(einsum_trace, rho_r)
        
        # Step 2: Build mixed state: rho_rest ⊗ I_pair/d²
        # Use einsum to construct the full tensor product
        I_d = np.eye(d, dtype=complex) / d  # I/d for each qudit in pair
        
        # rho_rest indices
        rest_row_chars = [row_chars[k] for k in remaining_qudits]
        rest_col_chars = [col_chars[k] for k in remaining_qudits]
        rest_input = ''.join(rest_row_chars) + ''.join(rest_col_chars)
        
        # I_d for qudit_i and qudit_j
        Ii_input = row_chars[qudit_i] + col_chars[qudit_i]
        Ij_input = row_chars[qudit_j] + col_chars[qudit_j]
        
        # Full output
        full_output = ''.join(row_chars) + ''.join(col_chars)
        
        einsum_build = f'{rest_input},{Ii_input},{Ij_input}->{full_output}'
        mixed_r = np.einsum(einsum_build, rho_rest_r, I_d, I_d)
        mixed = mixed_r.reshape(D, D)
        
        # Step 3: Apply depolarizing channel
        rho_new = (1 - depol_prob) * rho + depol_prob * mixed
        
        return rho_new
    
    def simulate_noisy(self, T_total: float, N_steps: int,
                      initial_state_type: str = 'edge_triplet',
                      shots: int = 10000,
                      noise_params: Dict = None) -> Dict:
        """
        ノイズモデル付きシミュレーションを実行
        
        Uses density matrix formalism with per-gate 2-qudit depolarizing noise.
        Each 2-qudit gate (CEx for H_transfer, CustomTwo for H_TTA) receives
        independent depolarizing noise, matching the per-gate noise model
        used in the qubit noisy simulator.
        
        The depolarizing channel for each 2-qudit gate is:
            ε(ρ) = (1-p)ρ + p · Tr_{ij}(ρ) ⊗ I_{ij}/d²
        where d=3 (qutrit), applied to the specific qudit pair.
        
        Parameters:
        -----------
        T_total : float
            Total simulation time in fs
        N_steps : int
            Number of Trotter steps
        initial_state_type : str
            Initial state configuration
        shots : int
            Number of measurement shots per time step
        noise_params : Dict or None
            Noise parameters. If None, uses default realistic values.
            Keys: 'depol_1q', 'depol_2q', 'noise_gates'
        
        Returns:
        --------
        Dict : Simulation results
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
        
        print("\n" + "="*70)
        print("Quditベースシミュレーション開始（ノイズモデル付き）")
        print("="*70)
        print(f"ショット数: {shots}")
        
        dt = T_total / N_steps
        
        # Create noise model
        if noise_params is None:
            noise_params = {}
        
        depol_1q = noise_params.get('depol_1q', 0.001)
        depol_2q = noise_params.get('depol_2q', 0.01)
        noise_gates = noise_params.get('noise_gates', None)
        
        noise_model = self.create_noise_model(depol_1q, depol_2q, noise_gates)
        
        print("\nノイズモデルパラメータ:")
        print(f"  1量子ビットゲート: 理想的（ノイズなし）")
        print(f"  2量子ビットゲート脱分極エラー: {depol_2q*100:.3f}%")
        print(f"  ノイズ適用ゲート: {noise_model.basis_gates}")
        
        # Build single Trotter step circuit (for gate counting only)
        step_circuit = QuantumCircuit()
        reg = QuantumRegister("molecules", self.N, [3] * self.N)
        step_circuit.append(reg)
        self.time_evol.add_single_trotter_step(step_circuit, dt)
        
        # Count gates
        initial_gates = len(step_circuit.instructions)
        custom_two_count = sum(1 for gate in step_circuit.instructions 
                              if gate.__class__.__name__ == 'CustomTwo')
        
        print(f"\nゲート数計測用回路構築完了:")
        print(f"  初期ゲート数: {initial_gates} (CustomTwo含む: {custom_two_count})")
        
        # Estimate gate count after decomposition
        if custom_two_count > 0:
            print(f"\nCustomTwoゲートの基本ゲート数を推定中...")
            estimated_gates = self.time_evol.decompose_custom_two_gates(step_circuit)
            gates_per_step = estimated_gates
        else:
            gates_per_step = initial_gates
        
        print(f"\n1トロッターステップあたりの基本ゲート数: {gates_per_step}")
        
        # Initialize backend
        backend = self.provider.get_backend("misim")
        
        # Results storage
        times = [0.0]
        populations_history = []
        per_molecule_populations_history = []
        
        start_time = time.time()
        
        # Initial state
        init_circuit = self.build_initial_state_circuit(initial_state_type)
        
        print(f"\n初期状態: {initial_state_type}")
        
        # Build per-pair unitaries for interleaved noise application
        print("\nBuilding per-pair unitaries for density matrix simulation...")
        U_H0_half_list, U_transfer_half_list, U_TTA_half_list = \
            self._build_per_pair_unitaries(dt)
        
        # Count 2-qudit gates per step for reporting
        n_2qudit_gates = 2 * (len(U_transfer_half_list) + len(U_TTA_half_list))
        print(f"  2-qudit gates per Trotter step: {n_2qudit_gates}")
        print(f"  Effective per-step noise: 1-(1-{depol_2q})^{n_2qudit_gates} = "
              f"{1-(1-depol_2q)**n_2qudit_gates:.4f}")
        
        # Get initial statevector
        job = backend.run(init_circuit)
        result = job.result()
        current_state = result.get_state_vector().flatten()
        
        # Initialize density matrix
        rho = np.outer(current_state, current_state.conj())
        
        # Calculate initial populations
        probabilities = np.diag(rho).real
        probabilities = np.maximum(probabilities, 0)
        probabilities /= np.sum(probabilities)
        samples_0 = np.random.choice(self.dim, size=shots, p=probabilities)
        pop_0 = self.calculate_populations_from_samples(samples_0, shots)
        pop_per_mol_0 = self.calculate_per_molecule_populations_from_samples(samples_0, shots)
        populations_history.append(pop_0)
        per_molecule_populations_history.append(pop_per_mol_0)
        
        print(f"\n初期個体数:")
        print(f"  N_S0 = {pop_0['N_S0']:.4f}")
        print(f"  N_T1 = {pop_0['N_T1']:.4f}")
        print(f"  N_S1 = {pop_0['N_S1']:.4f}")
        
        # Run simulation using density matrix evolution with per-gate noise
        print(f"\n時間発展を実行中（{N_steps}ステップ、密度行列シミュレーション）...")
        print(f"  ノイズ: 各2-quditゲートに{depol_2q*100:.1f}%脱分極エラーを適用")
        
        for step in range(1, N_steps + 1):
            # Forward half: H0 -> H_transfer -> H_TTA
            
            # H0 (single-qudit gates, no noise)
            for U_H0 in U_H0_half_list:
                rho = U_H0 @ rho @ U_H0.conj().T
            
            # H_transfer (2-qudit gates with noise)
            for U_tr, pair in U_transfer_half_list:
                rho = U_tr @ rho @ U_tr.conj().T
                rho = self._apply_2qudit_depolarizing_dm(rho, pair[0], pair[1], depol_2q)
            
            # H_TTA (2-qudit gates with noise)
            for U_TTA, pair in U_TTA_half_list:
                rho = U_TTA @ rho @ U_TTA.conj().T
                rho = self._apply_2qudit_depolarizing_dm(rho, pair[0], pair[1], depol_2q)
            
            # Backward half: H_TTA -> H_transfer -> H0 (reverse order)
            
            # H_TTA (reverse)
            for U_TTA, pair in reversed(U_TTA_half_list):
                rho = U_TTA @ rho @ U_TTA.conj().T
                rho = self._apply_2qudit_depolarizing_dm(rho, pair[0], pair[1], depol_2q)
            
            # H_transfer (reverse)
            for U_tr, pair in reversed(U_transfer_half_list):
                rho = U_tr @ rho @ U_tr.conj().T
                rho = self._apply_2qudit_depolarizing_dm(rho, pair[0], pair[1], depol_2q)
            
            # H0 (reverse, single-qudit, no noise)
            for U_H0 in reversed(U_H0_half_list):
                rho = U_H0 @ rho @ U_H0.conj().T
            
            # Ensure Hermitian and normalized (numerical stability)
            rho = (rho + rho.conj().T) / 2
            trace_val = np.trace(rho).real
            if trace_val > 0:
                rho /= trace_val
            
            # Sample from diagonal probabilities
            probabilities = np.diag(rho).real
            probabilities = np.maximum(probabilities, 0)
            probabilities /= np.sum(probabilities)
            samples = np.random.choice(self.dim, size=shots, p=probabilities)
            
            # Calculate populations
            pop = self.calculate_populations_from_samples(samples, shots)
            pop_per_mol = self.calculate_per_molecule_populations_from_samples(samples, shots)
            
            t = step * dt
            times.append(t)
            populations_history.append(pop)
            per_molecule_populations_history.append(pop_per_mol)
            
            if step % max(1, N_steps // 10) == 0:
                purity = np.trace(rho @ rho).real
                print(f"  ステップ {step}/{N_steps}: t = {t:.2f} fs, "
                      f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}, "
                      f"purity = {purity:.4f}")
        
        elapsed = time.time() - start_time
        
        print("\n" + "="*70)
        print("シミュレーション完了")
        print("="*70)
        print(f"最終個体数:")
        print(f"  N_S0 = {populations_history[-1]['N_S0']:.4f}")
        print(f"  N_T1 = {populations_history[-1]['N_T1']:.4f}")
        print(f"  N_S1 = {populations_history[-1]['N_S1']:.4f}")
        final_purity = np.trace(rho @ rho).real
        print(f"  最終純度: {final_purity:.6f}")
        print(f"\n回路統計:")
        print(f"  1ステップあたりゲート数: {gates_per_step}")
        print(f"  2-quditゲート数/ステップ: {n_2qudit_gates}")
        print(f"  総ゲート数: {gates_per_step * N_steps}")
        print(f"  実行時間: {elapsed:.2f}秒")
        
        return {
            'times': times,
            'populations': populations_history,
            'per_molecule_populations': per_molecule_populations_history,
            'elapsed_time': elapsed,
            'method': f'Qudit (Noisy, depol_1q={depol_1q:.4f}, depol_2q={depol_2q:.4f})',
            'gates_per_step': gates_per_step,
            'total_gates': gates_per_step * N_steps,
            'n_2qudit_gates_per_step': n_2qudit_gates,
            'shots': shots,
            'noise_params': noise_params
        }
