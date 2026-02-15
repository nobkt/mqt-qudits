# PR#39 Phase 3詳細仕様書：MQT-Qudits統合

## 文書の目的

本仕様書は、PR#39のPhase 3（MQT-Quditsフレームワークへの統合）の完全な詳細仕様を提供します。
Phase 2で実装したゲート変換を、MQT-Quditsの量子回路フレームワークに統合し、
エンドツーエンドの最適化パイプラインを構築します。

## 前提条件

### Phase 1の成果（完了）

✅ **integrated_sparse_compiler.py**:

- 疎構造解析
- PR#37分解器の統合
- 忠実度 1.0 の達成

### Phase 2の成果（予定）

⏳ **gate_converter.py**:

- 2×2変換: ZYZ → MQT-Quditsゲート
- 3×3変換: Givens → MQT-Quditsゲート
- 忠実度 1.0 の保持

## Phase 3の目的

### 主要目標

1. **MQT-Qudits量子回路への統合**

   - QuantumCircuitクラスとの統合
   - ゲートの追加機能
   - 回路の構築と実行

2. **CompilerPassの実装**

   - CustomTwoゲートの最適化パス
   - 疎構造認識コンパイラパス
   - 既存のコンパイラパイプラインへの統合

3. **エンドツーエンドテスト**

   - 4分子鎖の完全な回路
   - 100ステップのトロッター分解
   - ゲート数削減の実測

4. **パフォーマンス最適化**
   - キャッシング
   - 並列処理
   - プロファイリングと最適化

### 期待される最終成果

```
現状（Phase 0）:
  - ゲート数: 約6,000/トロッターステップ
  - 忠実度: 0.24/0.63（不合格）
  - 計算時間: 長い

最終状態（Phase 3完了後）:
  - ゲート数: 約150/トロッターステップ（97.5%削減）
  - 忠実度: 1.0（完璧）
  - 計算時間: 大幅削減
  - エンドツーエンドで動作
```

## MQT-Quditsフレームワークの理解

### 基本構造

```python
from mqt.qudits import QuantumCircuit

# 量子回路の作成
circuit = QuantumCircuit(num_qudits=4, dimensions=[3, 3, 3, 3])

# ゲートの追加
circuit.virtrz(qudit_index=0, level=1, phase=np.pi / 4)
circuit.r(qudit_index=0, level1=0, level2=1, theta=np.pi / 2, phi=0.0)
circuit.cex(control_qudit=0, target_qudit=1, level1=0, level2=1)

# コンパイラパスの適用
from mqt.qudits.compiler import LogEntQRCEXPass

optimized_circuit = LogEntQRCEXPass().run(circuit)

# 実行
from mqt.qudits.simulation import MQTQuditProvider

provider = MQTQuditProvider()
backend = provider.get_backend("qasm_simulator")
result = backend.run(optimized_circuit)
```

### 既存のCompilerPass

```python
class LogEntQRCEXPass(CompilerPass):
    """
    既存のコンパイラパス

    問題点:
    - CustomTwoゲートを一般的に分解（約810ゲート）
    - 疎構造を無視
    - 忠実度が低い場合がある
    """

    def run(self, circuit):
        # CustomTwoゲートを基本ゲートに分解
        pass
```

## Phase 3の実装設計

### アーキテクチャ

```
入力: CustomTwoゲートを含む量子回路
  ↓
1. SparseStructureOptimizationPass（新規）
  ├─ CustomTwoゲートを検出
  ├─ 疎構造を解析
  ├─ Phase 1: 部分空間抽出と分解
  ├─ Phase 2: MQT-Quditsゲートに変換
  └─ 最適化されたゲートシーケンスを挿入
  ↓
出力: 最適化された量子回路（ゲート数97.5%削減）
```

### 1. MQT-Quditsゲートラッパー

```python
# tools/mqt_gate_builder.py

import numpy as np
from typing import List, Dict, Any


class MQTGateSequence:
    """MQT-Quditsゲートシーケンスの抽象表現"""

    def __init__(self):
        self.gates: List[Dict[str, Any]] = []

    def add_virtrz(self, qudit_index: int, level: int, phase: float):
        """VirtRzゲートを追加"""
        self.gates.append(
            {
                "type": "VirtRz",
                "qudit_index": qudit_index,
                "level": level,
                "phase": phase,
            }
        )

    def add_r(
        self, qudit_index: int, level1: int, level2: int, theta: float, phi: float = 0.0
    ):
        """Rゲートを追加"""
        self.gates.append(
            {
                "type": "R",
                "qudit_index": qudit_index,
                "level1": level1,
                "level2": level2,
                "theta": theta,
                "phi": phi,
            }
        )

    def get_gate_count(self) -> int:
        """ゲート数を取得"""
        return len(self.gates)

    def get_physical_gate_count(self) -> int:
        """物理ゲート数を取得（VirtRz以外）"""
        return sum(1 for g in self.gates if g["type"] != "VirtRz")

    def apply_to_circuit(self, circuit):
        """MQT-Qudits回路にゲートを追加"""
        for gate in self.gates:
            if gate["type"] == "VirtRz":
                circuit.virtrz(
                    qudit_index=gate["qudit_index"],
                    level=gate["level"],
                    phase=gate["phase"],
                )
            elif gate["type"] == "R":
                circuit.r(
                    qudit_index=gate["qudit_index"],
                    level1=gate["level1"],
                    level2=gate["level2"],
                    theta=gate["theta"],
                    phi=gate["phi"],
                )

        return circuit


class OptimizedGateBuilder:
    """
    最適化されたゲートシーケンスを構築

    Phase 1とPhase 2の結果を統合してMQT-Quditsゲートシーケンスを生成
    """

    @staticmethod
    def build_from_2x2_decomposition(
        params: Dict, qudit_index: int, active_indices: List[int]
    ) -> MQTGateSequence:
        """
        2×2分解からゲートシーケンスを構築

        Args:
            params: Phase 1の分解結果（ZYZ）
            qudit_index: quditのインデックス
            active_indices: 部分空間のインデックス

        Returns:
            MQTGateSequence
        """
        sequence = MQTGateSequence()

        theta = params["theta"]
        phi = params["phi"]
        lam = params["lambda"]
        alpha = params["global_phase"]

        level1, level2 = active_indices

        # VirtRz(α+φ, level1)
        phase1 = alpha + phi
        if abs(phase1) > 1e-10:
            sequence.add_virtrz(qudit_index, level1, phase1)

        # R(θ, 0, level1, level2)
        if abs(theta) > 1e-10:
            sequence.add_r(qudit_index, level1, level2, theta, 0.0)

        # VirtRz(λ, level2)
        if abs(lam) > 1e-10:
            sequence.add_virtrz(qudit_index, level2, lam)

        return sequence

    @staticmethod
    def build_from_3x3_decomposition(
        params: Dict, qudit_index: int, active_indices: List[int]
    ) -> MQTGateSequence:
        """
        3×3分解からゲートシーケンスを構築

        Args:
            params: Phase 1の分解結果（Givens）
            qudit_index: quditのインデックス
            active_indices: 部分空間のインデックス

        Returns:
            MQTGateSequence
        """
        sequence = MQTGateSequence()

        # Givens回転をゲートに変換
        for local_level1, local_level2, theta, phi in params["rotations"]:
            global_level1 = active_indices[local_level1]
            global_level2 = active_indices[local_level2]

            # Givens回転: Rz(φ/2) Ry(θ) Rz(-φ/2)
            alpha = phi / 2.0
            beta = -phi / 2.0

            if abs(alpha) > 1e-10:
                sequence.add_virtrz(qudit_index, global_level1, alpha)

            if abs(theta) > 1e-10:
                sequence.add_r(qudit_index, global_level1, global_level2, theta, 0.0)

            if abs(beta) > 1e-10:
                sequence.add_virtrz(qudit_index, global_level2, beta)

        # 対角位相
        for local_level, phase in enumerate(params["diagonal_phases"]):
            global_level = active_indices[local_level]
            if abs(phase) > 1e-10:
                sequence.add_virtrz(qudit_index, global_level, phase)

        return sequence
```

### 2. SparseStructureOptimizationPass

```python
# tools/sparse_optimization_pass.py

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from integrated_sparse_compiler import IntegratedSparseCompiler
from mqt_gate_builder import OptimizedGateBuilder


class SparseStructureOptimizationPass:
    """
    疎構造認識型最適化コンパイラパス

    CustomTwoゲートを疎構造を利用して効率的に分解します。
    MQT-QuditsのCompilerPassとして実装。
    """

    def __init__(self, tolerance: float = 1e-10):
        """
        Args:
            tolerance: 数値誤差の許容範囲
        """
        self.compiler = IntegratedSparseCompiler(tolerance)
        self.tolerance = tolerance
        self.statistics = {
            "original_gate_count": 0,
            "optimized_gate_count": 0,
            "customtwo_count": 0,
            "optimization_ratio": 0.0,
        }

    def run(self, circuit):
        """
        量子回路を最適化

        Args:
            circuit: MQT-Quditsの量子回路

        Returns:
            最適化された量子回路
        """
        from mqt.qudits import QuantumCircuit

        # 新しい回路を作成
        optimized_circuit = QuantumCircuit(
            num_qudits=circuit.num_qudits, dimensions=circuit.dimensions
        )

        # 統計情報の初期化
        self.statistics["original_gate_count"] = len(circuit.gates)
        optimized_gate_count = 0
        customtwo_count = 0

        # 各ゲートを処理
        for gate in circuit.gates:
            if self._is_customtwo_gate(gate):
                # CustomTwoゲート: 最適化
                customtwo_count += 1
                optimized_gates = self._optimize_customtwo_gate(gate)
                optimized_gate_count += optimized_gates.get_gate_count()
                optimized_gates.apply_to_circuit(optimized_circuit)
            else:
                # その他のゲート: そのまま追加
                optimized_circuit.append(gate)
                optimized_gate_count += 1

        # 統計情報の更新
        self.statistics["optimized_gate_count"] = optimized_gate_count
        self.statistics["customtwo_count"] = customtwo_count
        if self.statistics["original_gate_count"] > 0:
            self.statistics["optimization_ratio"] = (
                1.0 - optimized_gate_count / self.statistics["original_gate_count"]
            ) * 100.0

        return optimized_circuit

    def _is_customtwo_gate(self, gate) -> bool:
        """CustomTwoゲートかどうかを判定"""
        # MQT-Quditsの実装に応じて調整
        return hasattr(gate, "unitary") and gate.unitary.shape == (9, 9)

    def _optimize_customtwo_gate(self, gate):
        """
        CustomTwoゲートを最適化

        Args:
            gate: CustomTwoゲート

        Returns:
            MQTGateSequence
        """
        # ユニタリ行列を取得
        U = gate.unitary

        # Phase 1: 疎構造コンパイル
        result = self.compiler.compile(U)

        # Phase 2: MQT-Quditsゲートに変換
        qudit_index = gate.qudit_index  # ゲートが作用するquditのインデックス

        if result.structure_info.active_dimension == 2:
            # 2×2部分空間
            sequence = OptimizedGateBuilder.build_from_2x2_decomposition(
                result.decomposition_params,
                qudit_index,
                result.structure_info.active_subspace,
            )
        elif result.structure_info.active_dimension == 3:
            # 3×3部分空間
            sequence = OptimizedGateBuilder.build_from_3x3_decomposition(
                result.decomposition_params,
                qudit_index,
                result.structure_info.active_subspace,
            )
        else:
            # その他: 未実装
            raise NotImplementedError(
                f"部分空間次元 {result.structure_info.active_dimension} は未実装"
            )

        return sequence

    def get_statistics(self) -> Dict[str, Any]:
        """最適化統計情報を取得"""
        return self.statistics.copy()
```

### 3. エンドツーエンドテスト

```python
# tools/end_to_end_test.py


def test_four_molecule_chain():
    """4分子鎖の完全なエンドツーエンドテスト"""
    from mqt.qudits import QuantumCircuit
    from sparse_optimization_pass import SparseStructureOptimizationPass

    print("=" * 70)
    print("4分子鎖エンドツーエンドテスト")
    print("=" * 70)

    # パラメータ
    num_trotter_steps = 10
    dt = 0.1

    # 量子回路の構築
    circuit = construct_four_molecule_circuit(num_trotter_steps, dt)

    print(f"\n元の回路:")
    print(f"  ゲート数: {len(circuit.gates)}")
    print(f"  トロッターステップ数: {num_trotter_steps}")

    # 最適化パスの適用
    optimizer = SparseStructureOptimizationPass()
    optimized_circuit = optimizer.run(circuit)

    print(f"\n最適化された回路:")
    print(f"  ゲート数: {len(optimized_circuit.gates)}")

    # 統計情報
    stats = optimizer.get_statistics()
    print(f"\n統計情報:")
    print(f"  元のゲート数: {stats['original_gate_count']}")
    print(f"  最適化後: {stats['optimized_gate_count']}")
    print(f"  CustomTwoゲート数: {stats['customtwo_count']}")
    print(f"  削減率: {stats['optimization_ratio']:.1f}%")

    # 検証: ゲート数が期待値以下
    expected_max_gates = 200 * num_trotter_steps  # 約150-200/ステップ
    assert len(optimized_circuit.gates) <= expected_max_gates

    # 検証: 削減率が95%以上
    assert stats["optimization_ratio"] >= 95.0

    print(f"\n✓ エンドツーエンドテスト合格")
    print(f"  ゲート数削減: {stats['optimization_ratio']:.1f}%")
    print(f"  期待値（97.5%）に近い")

    return True


def construct_four_molecule_circuit(num_trotter_steps: int, dt: float):
    """
    4分子鎖の量子回路を構築

    H_transfer（3個）とH_TTA（3個）のCustomTwoゲートを含む
    """
    from mqt.qudits import QuantumCircuit
    import numpy as np

    # 4つのqutrit
    circuit = QuantumCircuit(num_qudits=4, dimensions=[3, 3, 3, 3])

    # トロッターステップ
    for step in range(num_trotter_steps):
        # H_transfer（3個）
        for i in range(3):
            U_transfer = construct_h_transfer_unitary(dt)
            circuit.custom_two(qudit_index=i, unitary=U_transfer)

        # H_TTA（3個）
        for i in range(3):
            U_tta = construct_h_tta_unitary(dt)
            circuit.custom_two(qudit_index=i, unitary=U_tta)

    return circuit


def construct_h_transfer_unitary(dt: float) -> np.ndarray:
    """H_transferユニタリを構築（2×2部分空間）"""
    theta = 0.1 * dt
    U = np.eye(9, dtype=complex)
    U[1, 1] = np.cos(theta)
    U[1, 3] = -1j * np.sin(theta)
    U[3, 1] = -1j * np.sin(theta)
    U[3, 3] = np.cos(theta)
    return U


def construct_h_tta_unitary(dt: float) -> np.ndarray:
    """H_TTAユニタリを構築（3×3部分空間）"""
    J = 0.05
    hbar = 0.6582119569

    H_sub = J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=complex)

    eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
    phases = np.exp(-1j * eigenvalues * dt / hbar)
    U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T

    U = np.eye(9, dtype=complex)
    active_indices = [0, 1, 2]
    for i, idx_i in enumerate(active_indices):
        for j, idx_j in enumerate(active_indices):
            U[idx_i, idx_j] = U_sub[i, j]

    return U
```

## 実装スケジュール

### Phase 3全体（4-8週間）

#### Week 1-2: 基盤実装

- MQTGateSequenceクラス
- OptimizedGateBuilderクラス
- 基本的なテスト

#### Week 3-4: CompilerPass実装

- SparseStructureOptimizationPass
- CustomTwoゲートの検出と最適化
- 統計情報の収集

#### Week 5-6: エンドツーエンドテスト

- 4分子鎖の回路構築
- 最適化パスの適用
- ゲート数削減の実測

#### Week 7: パフォーマンス最適化

- プロファイリング
- キャッシング
- 並列処理

#### Week 8: ドキュメントと最終化

- 完全なドキュメント
- チュートリアル
- 論文執筆の準備

## 成功基準

### 必須基準

1. ✅ **ゲート数削減**: 95%以上（目標97.5%）
2. ✅ **忠実度**: 1.0（完璧）
3. ✅ **エンドツーエンド動作**: 4分子鎖で完璧に動作
4. ✅ **MQT-Qudits互換**: 完全に統合

### 望ましい基準

5. ✅ **パフォーマンス**: 高速な最適化
6. ✅ **スケーラビリティ**: 大規模な回路に対応
7. ✅ **ドキュメント**: 完全で正確
8. ✅ **コミュニティ**: MQT-Quditsへの貢献

## デリバラブル

### コード（tools/下）

1. **mqt_gate_builder.py**

   - MQTGateSequence
   - OptimizedGateBuilder

2. **sparse_optimization_pass.py**

   - SparseStructureOptimizationPass
   - 統計情報

3. **end_to_end_test.py**
   - 包括的なテスト
   - 4分子鎖の実装

### ドキュメント（tutorials/doc/下）

1. **pr39_phase3_implementation_report_ja.md**（実装後）

   - 実装の詳細
   - テスト結果
   - 最終的な成果

2. **mqt_qudits_optimization_guide_ja.md**
   - ユーザーガイド
   - APIリファレンス
   - チュートリアル

## リスクと緩和策

### 技術的リスク

1. **MQT-Qudits互換性**

   - 緩和策: 早期の統合テスト
   - 代替案: スタンドアロン実装

2. **パフォーマンス**

   - 緩和策: プロファイリングと最適化
   - 代替案: 並列処理の導入

3. **複雑性**
   - 緩和策: 段階的な実装とレビュー
   - 代替案: 範囲の調整

### プロジェクトリスク

1. **工数超過**

   - 緩和策: 週次レビューと調整
   - 代替案: Phase分割

2. **品質問題**
   - 緩和策: 包括的なテストスイート
   - 代替案: より多くのレビュー

## 最終目標

Phase 3完了後の最終状態:

```
量子回路の性能:
  - ゲート数: 6,000 → 150（97.5%削減）
  - 忠実度: 0.24/0.63 → 1.0（完璧）
  - 計算時間: 大幅削減

実用性:
  - Qubitと競争力のある性能
  - 実問題で実証済み
  - MQT-Quditsへの貢献

コミュニティへの影響:
  - 研究論文の執筆
  - オープンソースへの貢献
  - Qudit量子計算の実用化
```

---

**文書作成日**: 2025年10月21日
**作成者**: GitHub Copilot AI分析システム
**バージョン**: 1.0
**ステータス**: Phase 3詳細仕様
