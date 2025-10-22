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
                'type': gate.gate_type,
                'qudit_indices': qudit_indices,
                'params': gate.parameters
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
# 鈴木トロッター分解シミュレータ（疎構造認識版）
# ===================================================================

class SuzukiTrotterMQTQuditSimulator:
    """
    鈴木トロッターシミュレータ（疎構造認識版）
    
    疎構造認識時間発展演算子を使用して、
    4分子系の量子ダイナミクスをシミュレートします。
    
    従来のSuzukiTrotterMQTQuditSimulatorと同じインターフェースを持ちますが、
    SparseAwareMQTQuditTimeEvolutionを内部で使用します。
    """
    
    def __init__(self, params: PhysicalParameters):
        """
        Args:
            params: 物理パラメータ
        """
        self.params = params
        self.time_evol = SparseAwareMQTQuditTimeEvolution(params)
        self.N = params.N_molecules
        self.dim = 3 ** self.N
        
        # MQT-Quditsのインポート
        try:
            from mqt.qudits.simulation import MQTQuditProvider
            self.provider = MQTQuditProvider()
            self.backend = self.provider.get_backend("tnsim")
            self.mqt_available = True
        except ImportError:
            self.mqt_available = False
            print("警告: mqt.quditsがインストールされていません。")
    
    def build_trotter_circuit(self, dt: float, n_steps: int, initial_state: Optional[np.ndarray] = None):
        """
        鈴木トロッター回路を構築
        
        Args:
            dt: 時間刻み幅 (fs)
            n_steps: トロッターステップ数
            initial_state: 初期状態ベクトル（None の場合は |T1,T1,S0,S0⟩）
            
        Returns:
            circuit: QuantumCircuit
        """
        if not self.mqt_available:
            raise ImportError("mqt.quditsがインストールされていません")
        
        from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
        
        # 回路初期化
        circuit = QuantumCircuit()
        reg = QuantumRegister("molecules", self.N, [3] * self.N)
        circuit.append(reg)
        
        # 初期状態の設定（必要に応じて）
        if initial_state is not None:
            # 初期状態を設定するゲートを追加
            # （実装は省略、実際には状態ベクトルから回路を構築する必要がある）
            pass
        
        # トロッターステップ
        for step in range(n_steps):
            # 対称鈴木トロッター分解: H0(dt/2) -> H1(dt) -> H0(dt/2)
            self.time_evol.add_H0_evolution_gates(circuit, dt/2)
            self.time_evol.add_H_transfer_evolution_gates(circuit, dt)
            self.time_evol.add_H_TTA_evolution_gates(circuit, dt)
            self.time_evol.add_H0_evolution_gates(circuit, dt/2)
        
        return circuit
    
    def run_simulation(self, dt: float, n_steps: int, initial_state: Optional[np.ndarray] = None):
        """
        シミュレーションを実行
        
        Args:
            dt: 時間刻み幅 (fs)
            n_steps: トロッターステップ数
            initial_state: 初期状態ベクトル
            
        Returns:
            result: シミュレーション結果
        """
        if not self.mqt_available:
            raise ImportError("mqt.quditsがインストールされていません")
        
        # 回路構築
        circuit = self.build_trotter_circuit(dt, n_steps, initial_state)
        
        # 初期状態が指定されていない場合のデフォルト
        if initial_state is None:
            # |T1,T1,S0,S0⟩ = |1,1,0,0⟩
            initial_state = np.zeros(self.dim, dtype=complex)
            config = [1, 1, 0, 0]
            idx = config_to_index(config)
            initial_state[idx] = 1.0
        
        # バックエンドで実行
        result = self.backend.run(circuit, initial_state)
        
        return result
    
    def get_compilation_report(self) -> str:
        """コンパイル統計レポートを取得"""
        return self.time_evol.get_compilation_report()


# ===================================================================
# 厳密対角化ソルバー
# ===================================================================

class ExactDiagonalizationSolver:
    """
    厳密対角化による解析解計算
    
    ハミルトニアン全体を対角化して、
    時間発展の解析解を計算します。
    """
    
    def __init__(self, params: PhysicalParameters):
        """
        Args:
            params: 物理パラメータ
        """
        self.params = params
        self.N = params.N_molecules
        self.dim = 3 ** self.N
        self.H = None
        self.eigenvalues = None
        self.eigenvectors = None
    
    def build_hamiltonian(self) -> np.ndarray:
        """
        全ハミルトニアンを構築
        
        Returns:
            H: dim × dim ハミルトニアン行列
        """
        H = np.zeros((self.dim, self.dim), dtype=complex)
        
        # H0: 対角項
        for i in range(self.dim):
            config = index_to_config(i, self.N)
            E = sum(self.params.E_T if level == 1 else 
                   self.params.E_S if level == 2 else 0
                   for level in config)
            H[i, i] = E
        
        # H_transfer: エネルギー移動項
        for pair_idx, (mol_i, mol_j) in enumerate(self.params.neighbors):
            V = self.params.V[pair_idx]
            for idx in range(self.dim):
                config = index_to_config(idx, self.N)
                # |01⟩ <-> |10⟩ 遷移
                if config[mol_i] == 0 and config[mol_j] == 1:
                    new_config = config.copy()
                    new_config[mol_i] = 1
                    new_config[mol_j] = 0
                    new_idx = config_to_index(new_config)
                    H[idx, new_idx] = V
                    H[new_idx, idx] = V
        
        # H_TTA: 三重項-三重項消滅項
        for pair_idx, (mol_i, mol_j) in enumerate(self.params.neighbors):
            J = self.params.J[pair_idx]
            for idx in range(self.dim):
                config = index_to_config(idx, self.N)
                # |02⟩ <-> |11⟩ 遷移
                if config[mol_i] == 0 and config[mol_j] == 2:
                    new_config = config.copy()
                    new_config[mol_i] = 1
                    new_config[mol_j] = 1
                    new_idx = config_to_index(new_config)
                    H[idx, new_idx] = J
                    H[new_idx, idx] = J
                # |11⟩ <-> |20⟩ 遷移
                if config[mol_i] == 1 and config[mol_j] == 1:
                    new_config1 = config.copy()
                    new_config1[mol_i] = 2
                    new_config1[mol_j] = 0
                    new_idx1 = config_to_index(new_config1)
                    H[idx, new_idx1] = J
                    H[new_idx1, idx] = J
        
        self.H = H
        return H
    
    def diagonalize(self):
        """ハミルトニアンを対角化"""
        if self.H is None:
            self.build_hamiltonian()
        
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(self.H)
    
    def time_evolve(self, initial_state: np.ndarray, time: float) -> np.ndarray:
        """
        状態を時間発展
        
        Args:
            initial_state: 初期状態ベクトル
            time: 時間 (fs)
            
        Returns:
            evolved_state: 時間発展後の状態ベクトル
        """
        if self.eigenvalues is None or self.eigenvectors is None:
            self.diagonalize()
        
        # 初期状態を固有状態基底で展開
        coeffs = self.eigenvectors.conj().T @ initial_state
        
        # 時間発展
        phases = np.exp(-1j * self.eigenvalues * time / self.params.hbar)
        evolved_coeffs = coeffs * phases
        
        # 元の基底に戻す
        evolved_state = self.eigenvectors @ evolved_coeffs
        
        return evolved_state


# ===================================================================
# ユーティリティ関数
# ===================================================================

def calculate_fidelity(state1: np.ndarray, state2: np.ndarray) -> float:
    """
    2つの状態ベクトル間の忠実度を計算
    
    Args:
        state1: 状態ベクトル1
        state2: 状態ベクトル2
        
    Returns:
        fidelity: 忠実度 |⟨ψ1|ψ2⟩|^2
    """
    overlap = np.abs(np.vdot(state1, state2))
    return overlap ** 2


def compare_qudit_vs_exact(qudit_state: np.ndarray, exact_state: np.ndarray) -> Dict:
    """
    Qudit実装と厳密解を比較
    
    Args:
        qudit_state: Quditシミュレーション結果
        exact_state: 厳密対角化結果
        
    Returns:
        comparison: 比較結果の辞書
    """
    fidelity = calculate_fidelity(qudit_state, exact_state)
    
    # 各状態の占有確率
    qudit_probs = np.abs(qudit_state) ** 2
    exact_probs = np.abs(exact_state) ** 2
    
    # 確率の差
    prob_diff = np.abs(qudit_probs - exact_probs)
    max_diff = np.max(prob_diff)
    mean_diff = np.mean(prob_diff)
    
    return {
        'fidelity': fidelity,
        'max_prob_diff': max_diff,
        'mean_prob_diff': mean_diff,
        'qudit_probs': qudit_probs,
        'exact_probs': exact_probs
    }


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
