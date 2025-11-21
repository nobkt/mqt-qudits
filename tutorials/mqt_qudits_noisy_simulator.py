#!/usr/bin/env python3
"""
MQT-Qudits Molecular Dynamics Simulator with Noise Support.

This module extends the exact Qudit simulator to support noise models
using MQT-qudits' noise capabilities.

Noise models supported:
- Depolarizing error: Random errors on qudits
- Dephasing error: Phase damping
- SubspaceNoise: Physical noise on specific level transitions

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
    SparseAwareMQTTimeEvolution,
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
    
    def __init__(self, params: PhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.dim = 3 ** self.N
        
        # MQT-Quditsのインポート
        try:
            from mqt.qudits.simulation import MQTQuditProvider
            from mqt.qudits.simulation.noise_tools import NoiseModel, SubspaceNoise, Noise
            self.provider = MQTQuditProvider()
            self.NoiseModel = NoiseModel
            self.SubspaceNoise = SubspaceNoise
            self.Noise = Noise
            self.mqt_available = True
        except ImportError as e:
            self.mqt_available = False
            raise ImportError(f"mqt.quditsがインストールされていません: {e}")
        
        # Time evolution module for building circuits
        self.time_evol = SparseAwareMQTTimeEvolution(params)
        
        print(f"Noisy Qudit シミュレータを初期化しました")
        print(f"  分子数: {self.N}")
        print(f"  必要Qutrit数: {self.N}")
        print(f"  状態空間次元: {self.dim}")
    
    def create_noise_model(self,
                          depol_prob: float = 0.001,
                          dephasing_prob: float = 0.001,
                          noise_gates: List[str] = None) -> "NoiseModel":
        """
        Create a realistic noise model for qutrits.
        
        Parameters:
        -----------
        depol_prob : float
            Depolarizing error probability (default: 0.001 = 0.1%)
        dephasing_prob : float
            Dephasing error probability (default: 0.001 = 0.1%)
        noise_gates : List[str] or None
            List of gate names to apply noise to. If None, applies to all gates.
        
        Returns:
        --------
        NoiseModel : MQT-qudits noise model
        """
        noise_model = self.NoiseModel()
        
        if noise_gates is None:
            # Apply noise to all common gates
            noise_gates = ['virtrz', 'r', 'rz', 'rh', 'cx', 'h', 'x', 'z', 's']
        
        # Create noise for all level transitions in a qutrit (0-1, 0-2, 1-2)
        # Use SubspaceNoise to specify noise for each transition
        subspace_noise = self.SubspaceNoise(depol_prob, dephasing_prob, 
                                           [(0, 1), (0, 2), (1, 2)])
        
        # Apply noise locally to all qudits
        noise_model.add_quantum_error_locally(subspace_noise, noise_gates)
        
        # For two-qudit gates, add stronger noise
        two_qudit_gates = [g for g in noise_gates if g in ['cx', 'csum', 'ls', 'ms']]
        if two_qudit_gates:
            two_qudit_noise = self.SubspaceNoise(depol_prob * 5, dephasing_prob * 3,
                                                [(0, 1), (0, 2), (1, 2)])
            noise_model.add_nonlocal_quantum_error(two_qudit_noise, two_qudit_gates)
        
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
    
    def simulate_noisy(self, T_total: float, N_steps: int,
                      initial_state_type: str = 'edge_triplet',
                      shots: int = 10000,
                      noise_params: Dict = None) -> Dict:
        """
        ノイズモデル付きシミュレーションを実行
        
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
            Keys: 'depol_prob', 'dephasing_prob', 'noise_gates'
        
        Returns:
        --------
        Dict : Simulation results
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
        from mqt.qudits.simulation.noise_tools import NoisyCircuitFactory
        
        print("\n" + "="*70)
        print("Quditベースシミュレーション開始（ノイズモデル付き）")
        print("="*70)
        print(f"ショット数: {shots}")
        
        dt = T_total / N_steps
        
        # Create noise model
        if noise_params is None:
            noise_params = {}
        
        depol_prob = noise_params.get('depol_prob', 0.001)
        dephasing_prob = noise_params.get('dephasing_prob', 0.001)
        noise_gates = noise_params.get('noise_gates', None)
        
        noise_model = self.create_noise_model(depol_prob, dephasing_prob, noise_gates)
        
        print("\nノイズモデルパラメータ:")
        print(f"  脱分極エラー確率: {depol_prob*100:.3f}%")
        print(f"  位相緩和エラー確率: {dephasing_prob*100:.3f}%")
        print(f"  ノイズ適用ゲート: {noise_model.basis_gates}")
        
        # Build single Trotter step circuit
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
        
        # Initialize backend with noise
        backend = self.provider.get_backend("misim")
        
        # Results storage
        times = [0.0]
        populations_history = []
        per_molecule_populations_history = []
        
        start_time = time.time()
        
        # Initial state
        init_circuit = self.build_initial_state_circuit(initial_state_type)
        
        print(f"\n初期状態: {initial_state_type}")
        
        # Run simulation for each time step
        print(f"\n時間発展を実行中（{N_steps}ステップ、各ステップ{shots}ショット）...")
        
        for step in range(N_steps + 1):
            # Build circuit up to current step
            circuit = self.build_initial_state_circuit(initial_state_type)
            
            for _ in range(step):
                circuit.compose(step_circuit, inplace=True)
            
            # Add noise to the circuit
            noisy_factory = NoisyCircuitFactory(noise_model, circuit)
            noisy_circuit = noisy_factory.generate_circuit()
            
            # Run with noise model on backend
            job = backend.run(noisy_circuit, noise_model=noise_model, shots=shots)
            result = job.result()
            
            # Get measurement counts
            counts = result.get_counts()
            
            # Convert counts to samples
            samples = []
            for measurement, count in counts:
                # measurement is the state index
                samples.extend([measurement] * count)
            
            samples = np.array(samples)
            
            # Calculate populations from samples
            pop = self.calculate_populations_from_samples(samples, shots)
            pop_per_mol = self.calculate_per_molecule_populations_from_samples(samples, shots)
            
            t = step * dt
            times.append(t)
            populations_history.append(pop)
            per_molecule_populations_history.append(pop_per_mol)
            
            if step % max(1, N_steps // 10) == 0:
                print(f"  ステップ {step}/{N_steps}: t = {t:.2f} fs, "
                      f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}")
        
        elapsed = time.time() - start_time
        
        print("\n" + "="*70)
        print("シミュレーション完了")
        print("="*70)
        print(f"最終個体数:")
        print(f"  N_S0 = {populations_history[-1]['N_S0']:.4f}")
        print(f"  N_T1 = {populations_history[-1]['N_T1']:.4f}")
        print(f"  N_S1 = {populations_history[-1]['N_S1']:.4f}")
        print(f"\n回路統計:")
        print(f"  1ステップあたりゲート数: {gates_per_step}")
        print(f"  総ゲート数: {gates_per_step * N_steps}")
        print(f"  実行時間: {elapsed:.2f}秒")
        
        return {
            'times': times,
            'populations': populations_history,
            'per_molecule_populations': per_molecule_populations_history,
            'elapsed_time': elapsed,
            'method': f'Qudit (Noisy, depol={depol_prob:.4f}, dephasing={dephasing_prob:.4f})',
            'gates_per_step': gates_per_step,
            'total_gates': gates_per_step * N_steps,
            'shots': shots,
            'noise_params': noise_params
        }
