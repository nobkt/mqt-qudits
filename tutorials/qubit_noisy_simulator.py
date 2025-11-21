#!/usr/bin/env python3
"""
Qubit Molecular Dynamics Simulator with Noise Support.

This module extends the exact Qubit simulator to support noise models
using Qiskit Aer's noise capabilities.

Noise models supported:
- Depolarizing error: Random Pauli errors on qubits
- Amplitude damping: Energy relaxation (T1 decay)
- Thermal relaxation: Combined T1 and T2 decay

All noise parameters are physically motivated and use realistic values
for superconducting qubits.
"""

import numpy as np
import time
from typing import Dict, List
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit.primitives import StatevectorSampler
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, thermal_relaxation_error, amplitude_damping_error

# Import exact Hamiltonian builders (UnitaryGate version)
from exact_qubit_hamiltonians import (
    apply_exact_H_transfer_qubit,
    apply_exact_H_TTA_qubit
)


class QubitMolecularDynamicsSimulatorNoisy:
    """Qubit-based simulator with realistic noise models"""
    
    def __init__(self, params):
        self.params = params
        self.N = params.N_molecules
        self.n_qubits = 2 * self.N
        
        print(f"Qubit Noisy シミュレータを初期化しました")
        print(f"  分子数: {self.N}")
        print(f"  必要Qubit数: {self.n_qubits}")
    
    def create_noise_model(self, 
                          depol_1q: float = 0.001,
                          depol_2q: float = 0.01,
                          t1: float = 50000.0,  # 50 μs in fs
                          t2: float = 70000.0,  # 70 μs in fs
                          gate_time_1q: float = 50.0,  # 50 fs
                          gate_time_2q: float = 300.0) -> NoiseModel:
        """
        Create a realistic noise model for superconducting qubits.
        
        Parameters:
        -----------
        depol_1q : float
            Depolarizing error probability for single-qubit gates (default: 0.001 = 0.1%)
        depol_2q : float
            Depolarizing error probability for two-qubit gates (default: 0.01 = 1%)
        t1 : float
            Energy relaxation time in fs (default: 50 μs = 50000 fs)
        t2 : float
            Dephasing time in fs (default: 70 μs = 70000 fs)
        gate_time_1q : float
            Single-qubit gate time in fs (default: 50 fs)
        gate_time_2q : float
            Two-qubit gate time in fs (default: 300 fs)
        
        Returns:
        --------
        NoiseModel : Qiskit Aer noise model
        """
        noise_model = NoiseModel()
        
        # Single-qubit gates
        error_1q = depolarizing_error(depol_1q, 1)
        noise_model.add_all_qubit_quantum_error(error_1q, ['u1', 'u2', 'u3', 'rz', 'ry', 'rx', 'x'])
        
        # Two-qubit gates
        error_2q = depolarizing_error(depol_2q, 2)
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
        
        # Thermal relaxation for single-qubit gates
        if t1 > 0 and t2 > 0:
            error_thermal_1q = thermal_relaxation_error(t1, t2, gate_time_1q)
            noise_model.add_all_qubit_quantum_error(error_thermal_1q, ['u1', 'u2', 'u3'])
            
            error_thermal_2q = thermal_relaxation_error(t1, t2, gate_time_2q)
            noise_model.add_all_qubit_quantum_error(error_thermal_2q, ['cx'])
        
        return noise_model
    
    def prepare_initial_state(self, circuit: QuantumCircuit, state_type: str = 'edge_triplet'):
        """初期状態を準備"""
        if state_type == 'edge_triplet':
            circuit.x(0)  # 分子0の右側qubit
            circuit.x(2 * (self.N - 1))  # 分子N-1の右側qubit
        elif state_type == 'all_triplet':
            for i in range(self.N):
                circuit.x(2 * i)
    
    def apply_H0_evolution(self, circuit: QuantumCircuit, mol_idx: int, dt: float):
        """対角ハミルトニアン H0 の時間発展"""
        q0 = 2 * mol_idx
        q1 = 2 * mol_idx + 1
        
        E_T = self.params.E_T
        E_S = self.params.E_S
        hbar = self.params.hbar
        
        # Pauli分解による実装
        beta = (E_S - E_T) / 4
        gamma = (E_T - E_S) / 4
        delta = -(E_T + E_S) / 4
        
        theta_0 = -2 * beta * dt / hbar
        theta_1 = -2 * gamma * dt / hbar
        theta_zz = -2 * delta * dt / hbar
        
        circuit.rz(theta_0, q0)
        circuit.rz(theta_1, q1)
        
        circuit.cx(q0, q1)
        circuit.rz(theta_zz, q1)
        circuit.cx(q0, q1)
    
    def build_single_trotter_step(self, dt: float) -> QuantumCircuit:
        """1トロッターステップの回路を構築（UnitaryGate使用）"""
        circuit = QuantumCircuit(self.n_qubits)
        
        # 前半: H0, H_transfer, H_TTA
        for i in range(self.N):
            self.apply_H0_evolution(circuit, i, dt/2)
        
        for i, j in self.params.neighbors:
            # UnitaryGate版のH_transfer
            apply_exact_H_transfer_qubit(circuit, i, j, 
                                        self.params.V, dt/2, self.params.hbar)
        
        for i, j in self.params.neighbors:
            # UnitaryGate版のH_TTA
            apply_exact_H_TTA_qubit(circuit, i, j,
                                   self.params.J, dt/2, self.params.hbar)
        
        # 後半: 逆順
        for i, j in reversed(self.params.neighbors):
            apply_exact_H_TTA_qubit(circuit, i, j,
                                   self.params.J, dt/2, self.params.hbar)
        
        for i, j in reversed(self.params.neighbors):
            apply_exact_H_transfer_qubit(circuit, i, j,
                                        self.params.V, dt/2, self.params.hbar)
        
        for i in reversed(range(self.N)):
            self.apply_H0_evolution(circuit, i, dt/2)
        
        return circuit
    
    def calculate_populations_from_counts(self, counts: dict, shots: int) -> Dict[str, float]:
        """測定カウントから個体数を計算"""
        N_S0 = N_T1 = N_S1 = 0.0
        unphysical = 0.0
        
        for bitstring, count in counts.items():
            if count == 0:
                continue
            
            prob = count / shots
            bits = bitstring
            
            is_unphysical = False
            mol_count_S0 = mol_count_T1 = mol_count_S1 = 0
            
            for mol in range(self.N):
                q0_bit = int(bits[-(2*mol+1)])
                q1_bit = int(bits[-(2*mol+2)])
                
                if q0_bit == 1 and q1_bit == 1:
                    is_unphysical = True
                    break
                
                if q1_bit == 0 and q0_bit == 0:
                    mol_count_S0 += 1
                elif q1_bit == 0 and q0_bit == 1:
                    mol_count_T1 += 1
                elif q1_bit == 1 and q0_bit == 0:
                    mol_count_S1 += 1
            
            if is_unphysical:
                unphysical += prob
            else:
                N_S0 += prob * mol_count_S0
                N_T1 += prob * mol_count_T1
                N_S1 += prob * mol_count_S1
        
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
        """
        完全なシミュレーションを実行（ノイズモデル付き）
        
        Parameters:
        -----------
        T_total : float
            Total simulation time
        N_steps : int
            Number of Trotter steps
        initial_state_type : str
            Initial state configuration
        shots : int
            Number of measurement shots
        noise_params : Dict or None
            Noise parameters. If None, uses default realistic values.
            Keys: 'depol_1q', 'depol_2q', 't1', 't2', 'gate_time_1q', 'gate_time_2q'
        """
        print("\n" + "="*70)
        print("Qubitベースシミュレーション開始（ノイズモデル付き）")
        print("="*70)
        print(f"ショット数: {shots}")
        
        # Create noise model
        if noise_params is None:
            noise_params = {}
        noise_model = self.create_noise_model(**noise_params)
        
        print("\nノイズモデルパラメータ:")
        depol_1q = noise_params.get('depol_1q', 0.001)
        depol_2q = noise_params.get('depol_2q', 0.01)
        t1 = noise_params.get('t1', 50000.0)
        t2 = noise_params.get('t2', 70000.0)
        print(f"  1量子ビットゲート脱分極エラー: {depol_1q*100:.3f}%")
        print(f"  2量子ビットゲート脱分極エラー: {depol_2q*100:.3f}%")
        print(f"  T1 (エネルギー緩和時間): {t1:.1f} fs")
        print(f"  T2 (位相緩和時間): {t2:.1f} fs")
        
        start_time = time.time()
        dt = T_total / N_steps
        
        # Use Aer simulator with noise
        simulator = AerSimulator(noise_model=noise_model)
        
        step_circuit = self.build_single_trotter_step(dt)
        
        # Transpile for the noisy simulator
        transpiled_step = transpile(step_circuit, simulator)
        
        print(f"\n1トロッターステップあたりのゲート数: {len(transpiled_step.data)}")
        print(f"回路深さ: {transpiled_step.depth()}\n")
        
        # 初期状態の確認（理想的な状態）
        init_circuit = QuantumCircuit(self.n_qubits)
        self.prepare_initial_state(init_circuit, initial_state_type)
        state_0 = Statevector(init_circuit)
        
        probabilities = state_0.probabilities_dict()
        N_S0 = N_T1 = N_S1 = 0.0
        for bitstring, prob in probabilities.items():
            if prob < 1e-15:
                continue
            bits = bitstring[::-1]
            for mol in range(self.N):
                q0_bit = int(bits[2*mol])
                q1_bit = int(bits[2*mol+1])
                if q1_bit == 0 and q0_bit == 0:
                    N_S0 += prob
                elif q1_bit == 0 and q0_bit == 1:
                    N_T1 += prob
                elif q1_bit == 1 and q0_bit == 0:
                    N_S1 += prob
        
        pop_0 = {'N_S0': N_S0, 'N_T1': N_T1, 'N_S1': N_S1, 'unphysical': 0.0}
        
        print(f"初期状態: {initial_state_type}")
        print(f"  N_S0 = {pop_0['N_S0']:.4f}")
        print(f"  N_T1 = {pop_0['N_T1']:.4f}")
        print(f"  N_S1 = {pop_0['N_S1']:.4f}")
        
        times = [0.0]
        populations = [pop_0]
        
        # 時間発展
        print(f"\n時間発展を実行中（{N_steps}ステップ、各ステップ{shots}ショット）...")
        for step in range(1, N_steps + 1):
            circuit = QuantumCircuit(self.n_qubits, self.n_qubits)
            self.prepare_initial_state(circuit, initial_state_type)
            
            for _ in range(step):
                circuit.compose(transpiled_step, inplace=True)
            
            circuit.measure(range(self.n_qubits), range(self.n_qubits))
            
            # Run on noisy simulator
            transpiled_circuit = transpile(circuit, simulator)
            job = simulator.run(transpiled_circuit, shots=shots)
            result = job.result()
            counts = result.get_counts()
            
            pop = self.calculate_populations_from_counts(counts, shots)
            
            t = step * dt
            times.append(t)
            populations.append(pop)
            
            if step % max(1, N_steps // 10) == 0:
                print(f"  ステップ {step}/{N_steps}: t = {t:.2f} fs, "
                      f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}, "
                      f"非物理的: {pop['unphysical']:.4f}")
        
        elapsed = time.time() - start_time
        
        # 回路統計
        circuit_no_measure = QuantumCircuit(self.n_qubits)
        self.prepare_initial_state(circuit_no_measure, initial_state_type)
        for _ in range(N_steps):
            circuit_no_measure.compose(transpiled_step, inplace=True)
        
        total_gates = len(circuit_no_measure.data)
        total_depth = circuit_no_measure.depth()
        
        print("\n" + "="*70)
        print("シミュレーション完了")
        print("="*70)
        print(f"最終個体数:")
        print(f"  N_S0 = {populations[-1]['N_S0']:.4f}")
        print(f"  N_T1 = {populations[-1]['N_T1']:.4f}")
        print(f"  N_S1 = {populations[-1]['N_S1']:.4f}")
        print(f"  非物理的状態: {populations[-1]['unphysical']:.6f}")
        print(f"\n回路統計:")
        print(f"  総ゲート数: {total_gates}")
        print(f"  総回路深さ: {total_depth}")
        print(f"  実行時間: {elapsed:.2f}秒")
        
        return {
            'times': times,
            'populations': populations,
            'elapsed_time': elapsed,
            'method': f'Qubit (Noisy, depol_1q={depol_1q:.4f}, depol_2q={depol_2q:.4f})',
            'total_gates': total_gates,
            'total_depth': total_depth,
            'shots': shots,
            'noise_params': noise_params
        }
