#!/usr/bin/env python3
"""
Qubit Molecular Dynamics Simulator using UnitaryGate.

This is an alternative to the basic gate implementation that uses
Qiskit's UnitaryGate for H_transfer and H_TTA evolution.

The physics is identical, but the gate representation differs:
- Basic gate version: Decomposes unitaries into CNOT, Rz, Ry, Rx
- UnitaryGate version: Uses high-level unitary gates directly

Both are mathematically exact (no approximations).
"""

import numpy as np
import time
from typing import Dict, List
from qiskit import QuantumCircuit, ClassicalRegister
from qiskit.primitives import StatevectorSampler
from qiskit.quantum_info import Statevector

# Import exact Hamiltonian builders (UnitaryGate version)
from exact_qubit_hamiltonians import (
    apply_exact_H_transfer_qubit,
    apply_exact_H_TTA_qubit
)


class QubitMolecularDynamicsSimulatorUnitary:
    """Qubit-based simulator using UnitaryGate for Hamiltonian evolution"""
    
    def __init__(self, params):
        self.params = params
        self.N = params.N_molecules
        self.n_qubits = 2 * self.N
        
        print(f"Qubit UnitaryGate シミュレータを初期化しました")
        print(f"  分子数: {self.N}")
        print(f"  必要Qubit数: {self.n_qubits}")
    
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
        
        # Pauli分解による実装（UnitaryGate版でも同じ）
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
                 shots: int = 10000) -> Dict:
        """完全なシミュレーションを実行（UnitaryGate版）"""
        print("\n" + "="*70)
        print("Qubitベースシミュレーション開始（UnitaryGate版）")
        print("="*70)
        print(f"ショット数: {shots}")
        
        start_time = time.time()
        dt = T_total / N_steps
        
        sampler = StatevectorSampler()
        step_circuit = self.build_single_trotter_step(dt)
        
        print(f"\n1トロッターステップあたりのゲート数: {len(step_circuit.data)}")
        print(f"回路深さ: {step_circuit.depth()}\n")
        
        # 初期状態の確認
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
                circuit.compose(step_circuit, inplace=True)
            
            circuit.measure(range(self.n_qubits), range(self.n_qubits))
            
            job = sampler.run([circuit], shots=shots)
            result = job.result()
            counts = result[0].data.c.get_counts()
            
            pop = self.calculate_populations_from_counts(counts, shots)
            
            t = step * dt
            times.append(t)
            populations.append(pop)
            
            if step % max(1, N_steps // 10) == 0:
                print(f"  ステップ {step}/{N_steps}: t = {t:.2f} fs, "
                      f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}")
        
        elapsed = time.time() - start_time
        
        # 回路統計
        circuit_no_measure = QuantumCircuit(self.n_qubits)
        self.prepare_initial_state(circuit_no_measure, initial_state_type)
        for _ in range(N_steps):
            circuit_no_measure.compose(step_circuit, inplace=True)
        
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
            'circuit_final': circuit_no_measure,
            'step_circuit': step_circuit,
            'elapsed_time': elapsed,
            'total_gates': total_gates,
            'total_depth': total_depth,
            'gates_per_step': len(step_circuit.data),
            'depth_per_step': step_circuit.depth(),
            'method': 'Qubit (Qiskit - UnitaryGate)',
            'shots': shots
        }
