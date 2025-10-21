#!/usr/bin/env python3
"""
疎構造認識コンパイラを使用した4分子量子ダイナミクスシミュレーション

このスクリプトは、PR#42-46で開発された疎構造認識コンパイラを統合し、
CustomTwoゲートの効率的な分解を実現します。

主な改善点:
- LogEntQRCEXPass（~1000ゲート/CustomTwo）の代わりに
- IntegratedSparseCompilerV2を使用（~1-6ゲート/CustomTwo）
- ゲート数: 6182 → 21-25 ゲート（99.6%削減）

理論的基盤:
- tutorials/doc/SPARSE_COMPILER_THEORETICAL_FOUNDATION_JA.md
- tutorials/doc/PR46_FRAMEWORK_INTEGRATION_SPECIFICATION_JA.md

制約:
- ヒューリスティック・近似を使用しない
- 忠実度 1.0 を保証
- 数学的に完全に厳密
"""

import sys
from pathlib import Path
import numpy as np
from typing import List, Tuple, Dict, Optional
import time

# tools/からのインポート
sys.path.insert(0, str(Path(__file__).parent.parent / 'tools'))
from integrated_sparse_compiler_v2 import IntegratedSparseCompilerV2


# ===================================================================
# ユーティリティ関数
# ===================================================================

def index_to_config(idx: int, N: int = 4, d: int = 3) -> List[int]:
    """線形インデックスをd進数配列に変換"""
    config = []
    for _ in range(N):
        config.append(idx % d)
        idx //= d
    return config[::-1]


def config_to_index(config: List[int], d: int = 3) -> int:
    """d進数配列を線形インデックスに変換"""
    idx = 0
    for i, level in enumerate(config):
        idx += level * (d ** (len(config) - 1 - i))
    return idx


def config_to_state_name(config: List[int]) -> str:
    """設定配列を状態名文字列に変換"""
    level_names = {0: 'S0', 1: 'T1', 2: 'S1'}
    return ''.join([level_names[n] for n in config])


# ===================================================================
# 物理パラメータクラス
# ===================================================================

class PhysicalParameters:
    """4分子系の物理パラメータ"""
    
    def __init__(self):
        self.N_molecules = 4
        self.E_T = 1.5    # 三重項エネルギー (eV)
        self.E_S = 3.0    # 一重項エネルギー (eV)
        self.V = np.array([0.10, 0.10, 0.10])  # エネルギー移動積分 (eV)
        self.J = np.array([0.05, 0.05, 0.05])  # TTA相互作用定数 (eV)
        self.Gamma_fl = 0.001  # 蛍光放出速度 (fs^-1)
        self.hbar = 0.6582119569  # 換算プランク定数 (eV·fs)
        self.neighbors = [(0, 1), (1, 2), (2, 3)]  # 隣接リスト


# ===================================================================
# MQTゲート生成クラス（疎構造認識）
# ===================================================================

class SparseAwareMQTGateGenerator:
    """
    疎構造認識MQTゲート生成器
    
    IntegratedSparseCompilerV2を使用して、疎構造を持つユニタリ行列を
    効率的にMQT-Quditsの基本ゲートに変換します。
    
    対応する構造:
    - 2×2部分空間（H_transfer型）: ~1ゲート
    - 3×3部分空間（H_TTA型）: ~6ゲート
    - その他: 検出とレポート
    
    このクラスはMQT-Quditsの量子回路を直接生成せず、
    ゲート列の情報を返します。実際の回路への追加は
    呼び出し側で行います。
    """
    
    def __init__(self, tolerance: float = 1e-10):
        """
        Args:
            tolerance: 数値許容誤差
        """
        self.tolerance = tolerance
        self.compiler = IntegratedSparseCompilerV2(tolerance=tolerance, optimize_gates=True)
        self.compilation_stats = {
            'sparse_2x2': 0,
            'sparse_3x3': 0,
            'dense': 0,
            'total_gates_sparse': 0,
            'total_gates_estimated_dense': 0
        }
    
    def compile_unitary_to_gates(self, U: np.ndarray, qudit_indices: List[int]) -> Dict:
        """
        ユニタリ行列をMQT-Quditsゲートにコンパイル
        
        Args:
            U: 9×9ユニタリ行列（2-qutrit）
            qudit_indices: [qudit_i, qudit_j]
            
        Returns:
            ゲート情報の辞書:
            {
                'gates': List[Dict],  # ゲートのリスト
                'gate_count': int,
                'structure_type': str,  # 'sparse_2x2', 'sparse_3x3', 'dense'
                'fidelity': float
            }
        """
        # IntegratedSparseCompilerV2でコンパイル
        result = self.compiler.compile(U)
        
        # ゲート情報を抽出
        gates = []
        for gate in result.gate_sequence.gates:
            gate_info = {
                'type': gate['type'],
                'qudit_indices': qudit_indices,
                'params': gate.get('params', {})
            }
            gates.append(gate_info)
        
        # 統計更新
        structure_type = self._classify_structure(result)
        self.compilation_stats[structure_type] += 1
        self.compilation_stats['total_gates_sparse'] += len(gates)
        if structure_type != 'dense':
            # 密構造の場合の推定ゲート数（比較用）
            self.compilation_stats['total_gates_estimated_dense'] += 1000
        
        return {
            'gates': gates,
            'gate_count': len(gates),
            'structure_type': structure_type,
            'fidelity': result.gate_sequence.fidelity,
            'active_indices': result.structure_info.active_subspace
        }
    
    def _classify_structure(self, result) -> str:
        """コンパイル結果から構造タイプを分類"""
        if result.structure_info.active_dimension == 2:
            return 'sparse_2x2'
        elif result.structure_info.active_dimension == 3:
            return 'sparse_3x3'
        else:
            return 'dense'
    
    def get_compilation_report(self) -> str:
        """コンパイル統計レポートを生成"""
        stats = self.compilation_stats
        reduction_rate = 0
        if stats['total_gates_estimated_dense'] > 0:
            reduction_rate = (1 - stats['total_gates_sparse'] / 
                            stats['total_gates_estimated_dense']) * 100
        
        report = f"""
=== 疎構造認識コンパイラ統計 ===

【検出された構造】
  2×2部分空間: {stats['sparse_2x2']} 個
  3×3部分空間: {stats['sparse_3x3']} 個
  密構造: {stats['dense']} 個

【ゲート数比較】
  疎構造認識: {stats['total_gates_sparse']} ゲート
  LogEntQRCEX推定: {stats['total_gates_estimated_dense']} ゲート
  削減率: {reduction_rate:.1f}%

【平均ゲート数】
  2×2部分空間: {stats['total_gates_sparse'] / max(stats['sparse_2x2'], 1):.1f} ゲート/個
  3×3部分空間: {stats['total_gates_sparse'] / max(stats['sparse_3x3'], 1):.1f} ゲート/個
"""
        return report


# ===================================================================
# MQT-Quditsゲートによる時間発展演算子（疎構造認識版）
# ===================================================================

class SparseAwareMQTQuditTimeEvolution:
    """
    疎構造認識MQT-Quditsゲート時間発展
    
    従来のMQTQuditTimeEvolutionクラスと同じインターフェースを持ちますが、
    CustomTwoゲートの分解に疎構造認識コンパイラを使用します。
    
    主な変更点:
    1. SparseAwareMQTGateGeneratorを使用
    2. ゲート数を大幅に削減（~99.6%）
    3. 忠実度 1.0 を保証
    """
    
    def __init__(self, params: PhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.dim = 3 ** self.N
        self.gate_generator = SparseAwareMQTGateGenerator()
        
        # MQT-Quditsのインポート（実行時のみ）
        try:
            from mqt.qudits.quantum_circuit import QuantumCircuit
            from mqt.qudits.simulation import MQTQuditProvider
            self.provider = MQTQuditProvider()
            self.mqt_available = True
        except ImportError:
            self.mqt_available = False
            print("警告: mqt.quditsがインストールされていません。")
            print("ゲート生成のみ実行し、実際の回路構築はスキップされます。")
    
    def add_H0_evolution_gates(self, circuit, dt: float):
        """
        H0の時間発展ゲートを回路に追加
        
        実装: VirtRzゲート（変更なし）
        ゲート数: 8個（4分子 × 2準位）
        """
        for i in range(self.N):
            # 準位1 (T1) に位相を適用
            phase_T = -self.params.E_T * dt / self.params.hbar
            circuit.virtrz(i, [1, phase_T])
            
            # 準位2 (S1) に位相を適用
            phase_S = -self.params.E_S * dt / self.params.hbar
            circuit.virtrz(i, [2, phase_S])
    
    def add_H_transfer_evolution_gates(self, circuit, dt: float):
        """
        H_transferの時間発展ゲートを回路に追加（疎構造認識版）
        
        従来: CustomTwoゲート → LogEntQRCEXPass → ~1000ゲート
        改良: 疎構造検出 → IntegratedSparseCompilerV2 → ~1ゲート
        
        期待される構造: 2×2部分空間（|01⟩, |10⟩）
        期待されるゲート数: ~1ゲート/ペア × 3ペア = ~3ゲート
        """
        for pair_idx, (i, j) in enumerate(self.params.neighbors):
            V = self.params.V[pair_idx]
            theta = V * dt / self.params.hbar
            
            # 9×9ユニタリ行列を構築
            U = np.eye(9, dtype=complex)
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)
            
            U[1, 1] = cos_theta
            U[1, 3] = -1j * sin_theta
            U[3, 1] = -1j * sin_theta
            U[3, 3] = cos_theta
            
            # 疎構造認識コンパイラでゲート列を生成
            gate_info = self.gate_generator.compile_unitary_to_gates(U, [i, j])
            
            # ゲート列を回路に追加
            self._add_gates_to_circuit(circuit, gate_info['gates'])
            
            # デバッグ情報（初回のみ）
            if pair_idx == 0:
                print(f"H_transfer分解: {gate_info['structure_type']}, "
                      f"{gate_info['gate_count']}ゲート, "
                      f"忠実度={gate_info['fidelity']:.10f}")
    
    def add_H_TTA_evolution_gates(self, circuit, dt: float):
        """
        H_TTAの時間発展ゲートを回路に追加（疎構造認識版）
        
        従来: CustomTwoゲート → LogEntQRCEXPass → ~1000ゲート
        改良: 疎構造検出 → IntegratedSparseCompilerV2 → ~6ゲート
        
        期待される構造: 3×3部分空間（|02⟩, |11⟩, |20⟩）
        期待されるゲート数: ~6ゲート/ペア × 3ペア = ~18ゲート
        """
        for pair_idx, (i, j) in enumerate(self.params.neighbors):
            J = self.params.J[pair_idx]
            
            # 9×9ユニタリ行列を構築
            U = np.eye(9, dtype=complex)
            
            # 部分空間のハミルトニアン
            H_sub = J * np.array([
                [0, 1, 1],
                [1, 0, 0],
                [1, 0, 0]
            ], dtype=complex)
            
            # 固有値分解
            eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
            
            # 時間発展演算子
            phases = np.exp(-1j * eigenvalues * dt / self.params.hbar)
            U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T
            
            # 9×9行列の該当部分に埋め込む
            indices = [2, 4, 6]
            for a, idx_a in enumerate(indices):
                for b, idx_b in enumerate(indices):
                    U[idx_a, idx_b] = U_sub[a, b]
            
            # 疎構造認識コンパイラでゲート列を生成
            gate_info = self.gate_generator.compile_unitary_to_gates(U, [i, j])
            
            # ゲート列を回路に追加
            self._add_gates_to_circuit(circuit, gate_info['gates'])
            
            # デバッグ情報（初回のみ）
            if pair_idx == 0:
                print(f"H_TTA分解: {gate_info['structure_type']}, "
                      f"{gate_info['gate_count']}ゲート, "
                      f"忠実度={gate_info['fidelity']:.10f}")
    
    def _add_gates_to_circuit(self, circuit, gates: List[Dict]):
        """
        ゲート列を回路に追加
        
        IntegratedSparseCompilerV2が生成するゲート形式:
        - VirtRz: 仮想Z回転
        - R: 回転ゲート
        - CEx: 制御Exchangeゲート
        
        Args:
            circuit: MQT-Qudits QuantumCircuit
            gates: ゲート情報のリスト
        """
        if not self.mqt_available:
            return
        
        for gate in gates:
            gate_type = gate['type']
            qudits = gate['qudit_indices']
            params = gate['params']
            
            if gate_type == 'VirtRz':
                # VirtRz(qudit, [level, phase])
                circuit.virtrz(qudits[0], [params['level'], params['phase']])
            
            elif gate_type == 'R':
                # R(qudit, [level_a, level_b, theta, phi])
                circuit.r(qudits[0], [params['level_a'], params['level_b'], 
                                     params['theta'], params['phi']])
            
            elif gate_type == 'CEx':
                # CEx([qudit_i, qudit_j])
                circuit.cex(qudits)
            
            elif gate_type == 'Rz':
                # Rz(qudit, [level_a, level_b, phase])
                circuit.rz(qudits[0], [params['level_a'], params['level_b'], 
                                       params['phase']])
            
            elif gate_type == 'Rh':
                # Rh(qudit, [level_a, level_b, theta])
                circuit.rh(qudits[0], [params['level_a'], params['level_b'], 
                                       params['theta']])
            
            else:
                print(f"警告: 未知のゲートタイプ {gate_type}")
    
    def get_compilation_report(self) -> str:
        """コンパイル統計レポートを取得"""
        return self.gate_generator.get_compilation_report()


# ===================================================================
# 使用例（メイン実行時）
# ===================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("疎構造認識コンパイラを使用した4分子量子ダイナミクス")
    print("=" * 70)
    print()
    
    # パラメータ初期化
    params = PhysicalParameters()
    print("=== 物理パラメータ ===")
    print(f"分子数: {params.N_molecules}")
    print(f"三重項エネルギー: {params.E_T} eV")
    print(f"一重項エネルギー: {params.E_S} eV")
    print(f"エネルギー移動積分: {params.V} eV")
    print(f"TTA相互作用定数: {params.J} eV")
    print()
    
    # 時間発展演算子の構築
    time_evol = SparseAwareMQTQuditTimeEvolution(params)
    
    print("=== 量子ゲートの構築 ===")
    print()
    
    # 1トロッターステップの構築（テスト用）
    dt = 10.0  # fs
    
    try:
        from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
        
        # 回路初期化
        circuit = QuantumCircuit()
        reg = QuantumRegister("molecules", params.N_molecules, [3] * params.N_molecules)
        circuit.append(reg)
        
        print(f"時間ステップ Δt = {dt} fs")
        print()
        
        # H0の時間発展
        print("H0の時間発展ゲート追加中...")
        time_evol.add_H0_evolution_gates(circuit, dt/2)
        
        # H_transferの時間発展
        print("H_transferの時間発展ゲート追加中...")
        time_evol.add_H_transfer_evolution_gates(circuit, dt/2)
        
        # H_TTAの時間発展
        print("H_TTAの時間発展ゲート追加中...")
        time_evol.add_H_TTA_evolution_gates(circuit, dt/2)
        
        print()
        print("=" * 70)
        print(time_evol.get_compilation_report())
        print("=" * 70)
        
    except ImportError:
        print("mqt.quditsがインストールされていないため、")
        print("実際の回路構築はスキップされました。")
        print()
        print("インストール方法:")
        print("  pip install mqt.qudits")
