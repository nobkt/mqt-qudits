# PR#38 実装ロードマップ

## 文書の目的

本文書は、PR#37の成果（完璧なユニタリ分解器）を基に、
完全なQuditゲート最適化パイプラインを構築するための
包括的な実装ロードマップを提供します。

## エグゼクティブサマリー

### PR#37の成果

- ✅ **2×2分解**: 忠実度 1.0（100/100テスト合格）
- ✅ **3×3分解**: 忠実度 1.0（100/100テスト合格）
- ✅ **数学的厳密性**: ヒューリスティックゼロ、近似ゼロ

### PR#38の目標

**短期目標**（1-2週間）:
- sparse_structure_compiler.pyとの統合
- 忠実度 0.24/0.63 → 1.0 の達成
- ゲート数見積もりの精度向上

**中期目標**（1-2ヶ月）:
- MQT-Quditsゲートへの変換実装
- 完全な最適化パイプライン構築
- エンドツーエンドのテスト

**長期目標**（3-6ヶ月）:
- MQT-Quditsフレームワークへの貢献
- 一般的な問題への適用拡大
- 論文発表とコミュニティへの公開

### 期待される最終成果

```
現状:
  - ゲート数: 約6,000/トロッターステップ
  - 忠実度: 0.24-0.63（不合格）
  
最終状態:
  - ゲート数: 約150/トロッターステップ（97.5%削減）
  - 忠実度: 1.0（完璧）
  - 計算時間: 大幅削減
```

## フェーズ1: 統合実装（1-2週間、35-56時間）

### 目的

PR#37の完璧な分解器をsparse_structure_compiler.pyに統合し、
忠実度を1.0に改善する。

### 詳細仕様

詳細は`pr38_integration_specification_ja.md`を参照。

### タスク一覧

#### タスク1.1: 環境準備（2-3時間）

```bash
# 依存関係のインストール
pip install numpy scipy matplotlib

# PR#37分解器の動作確認
python tools/improved_unitary_decomposition.py
python tools/perfect_3x3_decomposition.py

# 期待される出力: すべてのテストで忠実度=1.0
```

**成功基準**:
- ✓ すべてのテストが合格
- ✓ 忠実度 > 0.9999（実際は1.0）

#### タスク1.2: TwoLevelRotationDecomposerの統合（8-12時間）

**実装方法**（推奨）:

```python
# sparse_structure_compiler.py

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from improved_unitary_decomposition import ImprovedTwoQubitDecomposer

class TwoLevelRotationDecomposer:
    """PR#37統合版2×2分解器"""
    
    @staticmethod
    def decompose_2x2_unitary(U: np.ndarray) -> Dict:
        """忠実度1.0の分解"""
        decomposer = ImprovedTwoQubitDecomposer()
        result = decomposer.decompose_zyz(U)
        
        return {
            'theta': result.theta,
            'phi': result.phi,
            'lambda': result.lam,
            'global_phase': result.global_phase,
            'fidelity': result.fidelity
        }
```

**テスト**:

```python
def test_2x2_integration():
    # ランダムユニタリでテスト
    for i in range(100):
        U = generate_random_2x2_unitary()
        result = TwoLevelRotationDecomposer.decompose_2x2_unitary(U)
        assert result['fidelity'] > 0.9999
```

**成功基準**:
- ✓ 100/100テストで忠実度 > 0.9999
- ✓ H_transferで正しく動作

#### タスク1.3: ThreeLevelRotationDecomposerの統合（10-15時間）

**実装方法**:

```python
from perfect_3x3_decomposition import Perfect3x3Decomposer

class ThreeLevelRotationDecomposer:
    """PR#37統合版3×3分解器"""
    
    @staticmethod
    def decompose_3x3_unitary(U: np.ndarray) -> Dict:
        """忠実度1.0の分解"""
        decomposer = Perfect3x3Decomposer()
        result = decomposer.decompose(U)
        
        # QからGivens回転を抽出
        rotations = extract_givens_from_q(result.Q)
        phases = decomposer.extract_diagonal_phases(result.R)
        
        return {
            'rotations': rotations,
            'diagonal_phases': phases,
            'fidelity': result.fidelity,
            'Q': result.Q,
            'R': result.R
        }
```

**Givens抽出の実装**:

```python
@staticmethod
def extract_givens_from_q(Q: np.ndarray) -> List[Tuple[int, int, float, float]]:
    """
    QをGivens回転のリストに分解
    
    Q = G(0,1) @ G(0,2) @ G(1,2)
    """
    rotations = []
    Q_work = Q.conj().T.copy()
    
    # G(0,1): Q†[1,0]をゼロにする
    a, b = Q_work[0, 0], Q_work[1, 0]
    if abs(b) > 1e-10:
        theta, phi = compute_givens_params(a, b)
        rotations.append((0, 1, theta, phi))
        G_dag = construct_givens(3, 0, 1, theta, phi).conj().T
        Q_work = G_dag @ Q_work
    
    # G(0,2): Q†[2,0]をゼロにする
    a, b = Q_work[0, 0], Q_work[2, 0]
    if abs(b) > 1e-10:
        theta, phi = compute_givens_params(a, b)
        rotations.append((0, 2, theta, phi))
        G_dag = construct_givens(3, 0, 2, theta, phi).conj().T
        Q_work = G_dag @ Q_work
    
    # G(1,2): Q†[2,1]をゼロにする
    a, b = Q_work[1, 1], Q_work[2, 1]
    if abs(b) > 1e-10:
        theta, phi = compute_givens_params(a, b)
        rotations.append((1, 2, theta, phi))
    
    return rotations
```

**テスト**:

```python
def test_3x3_integration():
    # ランダムユニタリでテスト
    for i in range(100):
        U = generate_random_3x3_unitary()
        result = ThreeLevelRotationDecomposer.decompose_3x3_unitary(U)
        assert result['fidelity'] > 0.9999
    
    # H_TTA実問題でテスト
    U_tta = construct_h_tta_unitary()
    result = ThreeLevelRotationDecomposer.decompose_3x3_unitary(U_tta)
    assert result['fidelity'] > 0.9999
```

**成功基準**:
- ✓ 100/100テストで忠実度 > 0.9999
- ✓ H_TTAで正しく動作

#### タスク1.4: SubspaceRotationOptimizerの更新（3-5時間）

```python
class SubspaceRotationOptimizer:
    """改良版ゲート数見積もり"""
    
    def estimate_gate_count(self, structure: SparseStructureInfo) -> int:
        """正確なゲート数見積もり"""
        if structure.active_dimension == 2:
            return 3  # VirtRz + R + VirtRz
        elif structure.active_dimension == 3:
            return 12  # 3×(VirtRz+R+VirtRz) + 3×VirtRz
        else:
            return self._fallback_estimate(structure)
```

#### タスク1.5: 統合テストと検証（10-15時間）

**テストスイート**:

```python
def test_integration_h_transfer():
    """H_transfer完全統合テスト"""
    # 9×9ユニタリを構築
    U_9x9 = construct_h_transfer_unitary()
    
    # 疎構造解析
    analyzer = SparseStructureAnalyzer()
    structure = analyzer.analyze(U_9x9)
    
    # 部分空間抽出
    U_sub = analyzer.extract_subspace_unitary(U_9x9, structure.active_subspace)
    
    # 分解
    decomposer = TwoLevelRotationDecomposer()
    result = decomposer.decompose_2x2_unitary(U_sub)
    
    # 検証
    assert result['fidelity'] > 0.9999
    
    # ゲート数確認
    optimizer = SubspaceRotationOptimizer()
    gate_count = optimizer.estimate_gate_count(structure)
    assert gate_count == 3

def test_integration_h_tta():
    """H_TTA完全統合テスト"""
    # 同様のテスト
    pass

def test_integration_regression():
    """回帰テスト"""
    # sparse_structure_compiler.pyの全機能が動作することを確認
    pass
```

**成功基準**:
- ✓ すべてのユニットテストが合格
- ✓ 回帰テストが合格
- ✓ H_transfer/H_TTAで忠実度=1.0

### デリバラブル

1. ✅ 統合されたsparse_structure_compiler.py
2. ✅ 包括的なテストスイート
3. ✅ 統合レポート（忠実度、ゲート数、性能）
4. ✅ 更新されたREADME.md

## フェーズ2: ゲート変換実装（2-4週間、60-80時間）

### 目的

QR分解の結果をMQT-Quditsの基本ゲート（CEx, R, Rz, VirtRz）に
変換する完全な実装を提供する。

### 詳細仕様

詳細は`pr38_gate_conversion_theory_ja.md`を参照。

### タスク一覧

#### タスク2.1: 2×2変換の実装（15-20時間）

**ZYZ → MQT-Quditsゲート変換**:

```python
def zyz_to_mqt_gates(theta: float, phi: float, lam: float,
                     global_phase: float,
                     level1: int, level2: int) -> List[MQTGate]:
    """
    ZYZ分解パラメータをMQT-Quditsゲートに変換
    
    Returns:
        [VirtRz(α+φ, level1), R(θ, 0, level1, level2), VirtRz(λ, level2)]
    """
    gates = []
    
    # VirtRz(α+φ)
    phase1 = global_phase + phi
    if abs(phase1) > 1e-10:
        gates.append(VirtRzGate(phase=phase1, level=level1))
    
    # R(θ, 0)
    if abs(theta) > 1e-10:
        gates.append(RGate(theta=theta, phi=0.0, 
                          level1=level1, level2=level2))
    
    # VirtRz(λ)
    if abs(lam) > 1e-10:
        gates.append(VirtRzGate(phase=lam, level=level2))
    
    return gates
```

**テスト**:

```python
def test_2x2_gate_conversion():
    """2×2ゲート変換のテスト"""
    # ランダムユニタリ
    U = generate_random_2x2_unitary()
    
    # ZYZ分解
    decomposer = ImprovedTwoQubitDecomposer()
    result = decomposer.decompose_zyz(U)
    
    # ゲート変換
    gates = zyz_to_mqt_gates(
        result.theta, result.phi, result.lam, result.global_phase,
        level1=0, level2=1
    )
    
    # ゲートから行列を再構築
    U_reconstructed = reconstruct_from_gates(gates)
    
    # 忠実度検証
    fidelity = compute_fidelity(U, U_reconstructed)
    assert fidelity > 0.9999
```

**成功基準**:
- ✓ 100/100テストで忠実度 > 0.9999
- ✓ ゲート数 = 3（期待値）
- ✓ 物理ゲート数 = 1（Rゲートのみ）

#### タスク2.2: 3×3変換の実装（25-35時間）

**QR → Givens → MQT-Quditsゲート変換**:

```python
def convert_3x3_to_mqt_gates(U: np.ndarray) -> List[MQTGate]:
    """
    3×3ユニタリをMQT-Quditsゲートに変換
    
    U = QR → Q = G1 G2 G3 → 各Gi → MQT-Quditsゲート
    """
    # Step 1: QR分解
    decomposer = Perfect3x3Decomposer()
    result = decomposer.decompose(U)
    
    # Step 2: QからGivens回転を抽出
    givens_rotations = extract_givens_from_q(result.Q)
    
    # Step 3: 各Givens回転をMQT-Quditsゲートに変換
    gates = []
    for level1, level2, theta, phi in givens_rotations:
        givens_gates = givens_to_mqt_gates(level1, level2, theta, phi)
        gates.extend(givens_gates)
    
    # Step 4: 対角位相をVirtRzゲートに変換
    phases = decomposer.extract_diagonal_phases(result.R)
    for level, phase in enumerate(phases):
        if abs(phase) > 1e-10:
            gates.append(VirtRzGate(phase=phase, level=level))
    
    return gates


def givens_to_mqt_gates(level1: int, level2: int,
                        theta: float, phi: float) -> List[MQTGate]:
    """
    Givens回転をMQT-Quditsゲートに変換
    
    G(i,j; θ, φ) → Rz(α) Ry(θ) Rz(β) → VirtRz + R + VirtRz
    """
    # ZYZパラメータに変換
    alpha = phi / 2.0
    beta = -phi / 2.0
    
    gates = []
    
    if abs(alpha) > 1e-10:
        gates.append(VirtRzGate(phase=alpha, level=level1))
    
    if abs(theta) > 1e-10:
        gates.append(RGate(theta=theta, phi=0.0,
                          level1=level1, level2=level2))
    
    if abs(beta) > 1e-10:
        gates.append(VirtRzGate(phase=beta, level=level2))
    
    return gates
```

**テスト**:

```python
def test_3x3_gate_conversion():
    """3×3ゲート変換のテスト"""
    # 100個のランダムユニタリでテスト
    for _ in range(100):
        U = generate_random_3x3_unitary()
        
        # ゲート変換
        gates = convert_3x3_to_mqt_gates(U)
        
        # 再構築
        U_reconstructed = reconstruct_from_gates(gates)
        
        # 忠実度検証
        fidelity = compute_fidelity(U, U_reconstructed)
        assert fidelity > 0.9999
    
    # H_TTA実問題でテスト
    U_tta = construct_h_tta_unitary()
    gates = convert_3x3_to_mqt_gates(U_tta)
    U_reconstructed = reconstruct_from_gates(gates)
    fidelity = compute_fidelity(U_tta, U_reconstructed)
    assert fidelity > 0.9999
```

**成功基準**:
- ✓ 100/100テストで忠実度 > 0.9999
- ✓ ゲート数 = 12（期待値）
- ✓ 物理ゲート数 = 3（Rゲートのみ）

#### タスク2.3: 部分空間への埋め込み（10-15時間）

```python
def embed_gates_in_9x9(gates: List[MQTGate],
                      active_indices: List[int],
                      total_dimension: int = 9) -> List[MQTGate]:
    """
    部分空間のゲートを完全な空間に埋め込み
    
    Args:
        gates: 部分空間のゲートリスト
        active_indices: [1, 3] など
        total_dimension: 完全な空間の次元（9）
        
    Returns:
        埋め込まれたゲートリスト
    """
    embedded_gates = []
    
    for gate in gates:
        if isinstance(gate, VirtRzGate):
            # 準位インデックスを変換
            local_level = gate.level
            global_level = active_indices[local_level]
            embedded_gates.append(VirtRzGate(
                phase=gate.phase,
                level=global_level
            ))
        
        elif isinstance(gate, RGate):
            # 準位インデックスを変換
            local_level1 = gate.level1
            local_level2 = gate.level2
            global_level1 = active_indices[local_level1]
            global_level2 = active_indices[local_level2]
            embedded_gates.append(RGate(
                theta=gate.theta,
                phi=gate.phi,
                level1=global_level1,
                level2=global_level2
            ))
    
    return embedded_gates
```

#### タスク2.4: 統合テスト（10-10時間）

```python
def test_full_conversion_pipeline():
    """完全な変換パイプラインのテスト"""
    
    # H_transferのテスト
    U_transfer_9x9 = construct_h_transfer_unitary()
    
    # 疎構造解析
    analyzer = SparseStructureAnalyzer()
    structure = analyzer.analyze(U_transfer_9x9)
    
    # 部分空間抽出
    U_sub = analyzer.extract_subspace_unitary(
        U_transfer_9x9, structure.active_subspace
    )
    
    # ZYZ分解
    decomposer = ImprovedTwoQubitDecomposer()
    result = decomposer.decompose_zyz(U_sub)
    
    # MQT-Quditsゲート変換
    gates = zyz_to_mqt_gates(
        result.theta, result.phi, result.lam, result.global_phase,
        level1=0, level2=1
    )
    
    # 9×9空間に埋め込み
    embedded_gates = embed_gates_in_9x9(
        gates, structure.active_subspace, total_dimension=9
    )
    
    # 検証
    print(f"H_transfer完全変換:")
    print(f"  元のゲート数: 810")
    print(f"  最適化後: {len(embedded_gates)}ゲート")
    print(f"  物理ゲート: {sum(g.cost for g in embedded_gates)}")
    print(f"  削減率: {(1 - len(embedded_gates)/810)*100:.1f}%")
    
    assert len(embedded_gates) <= 3
    assert sum(g.cost for g in embedded_gates) == 1
```

### デリバラブル

1. ✅ 完全なゲート変換実装
2. ✅ 包括的なテストスイート（100+テストケース）
3. ✅ 実問題での検証（H_transfer, H_TTA）
4. ✅ ゲート変換ドキュメント

## フェーズ3: MQT-Qudits統合（4-8週間、120-160時間）

### 目的

最適化されたゲートシーケンスをMQT-Quditsフレームワークに統合し、
エンドツーエンドの最適化パイプラインを構築する。

### タスク一覧

#### タスク3.1: MQT-Quditsゲートクラスの実装（20-30時間）

```python
# MQT-Quditsとの互換性を持つゲートクラス

from mqt.qudits import QuantumCircuit

class OptimizedGateBuilder:
    """最適化されたゲートシーケンスを構築"""
    
    @staticmethod
    def add_virtrz(circuit: QuantumCircuit, 
                   qudit_index: int, 
                   level: int, 
                   phase: float):
        """VirtRzゲートを追加"""
        circuit.virtrz(qudit_index, level, phase)
    
    @staticmethod
    def add_r_gate(circuit: QuantumCircuit,
                   qudit_index: int,
                   level1: int,
                   level2: int,
                   theta: float,
                   phi: float):
        """Rゲートを追加"""
        circuit.r(qudit_index, level1, level2, theta, phi)
```

#### タスク3.2: CompilerPassの実装（30-40時間）

```python
from mqt.qudits.compiler import CompilerPass

class SparseStructureOptimizationPass(CompilerPass):
    """疎構造を認識した最適化パス"""
    
    def run(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """
        CustomTwoゲートを最適化されたゲートシーケンスに置き換え
        """
        optimized_circuit = QuantumCircuit(circuit.num_qudits)
        
        for gate in circuit.gates:
            if isinstance(gate, CustomTwoGate):
                # 疎構造解析
                structure = self._analyze_structure(gate.unitary)
                
                if structure.structure_type == 'sparse_subspace':
                    # 最適化されたゲートシーケンスを使用
                    optimized_gates = self._optimize_gate(gate, structure)
                    for opt_gate in optimized_gates:
                        optimized_circuit.append(opt_gate)
                else:
                    # 標準的な分解を使用
                    optimized_circuit.append(gate)
            else:
                optimized_circuit.append(gate)
        
        return optimized_circuit
    
    def _optimize_gate(self, gate: CustomTwoGate, 
                      structure: SparseStructureInfo) -> List[Gate]:
        """ゲートの最適化"""
        # 部分空間抽出
        U_sub = extract_subspace(gate.unitary, structure.active_subspace)
        
        # 適切な分解器を選択
        if structure.active_dimension == 2:
            gates = self._optimize_2x2(U_sub, structure)
        elif structure.active_dimension == 3:
            gates = self._optimize_3x3(U_sub, structure)
        else:
            gates = self._fallback_decomposition(gate)
        
        return gates
```

#### タスク3.3: エンドツーエンドテスト（30-40時間）

```python
def test_end_to_end_optimization():
    """4分子鎖の完全な最適化テスト"""
    # 回路を構築
    circuit = construct_four_molecule_circuit(num_trotter_steps=10)
    
    # 最適化前のゲート数
    gate_count_before = circuit.count_gates()
    print(f"最適化前: {gate_count_before}ゲート")
    
    # 最適化パスを適用
    optimizer = SparseStructureOptimizationPass()
    optimized_circuit = optimizer.run(circuit)
    
    # 最適化後のゲート数
    gate_count_after = optimized_circuit.count_gates()
    print(f"最適化後: {gate_count_after}ゲート")
    
    # 削減率
    reduction = (1 - gate_count_after / gate_count_before) * 100
    print(f"削減率: {reduction:.1f}%")
    
    # 検証
    assert reduction > 95  # 95%以上の削減
    
    # 忠実度検証
    # （実際の量子状態のシミュレーションが必要）
```

#### タスク3.4: パフォーマンスチューニング（20-30時間）

- ゲートシーケンスの最適化
- キャッシュの実装
- 並列処理の導入
- プロファイリングと最適化

#### タスク3.5: ドキュメントとチュートリアル（20-20時間）

- ユーザーガイドの作成
- APIリファレンスの作成
- チュートリアルノートブックの作成
- 論文執筆の準備

### デリバラブル

1. ✅ MQT-Qudits互換のCompilerPass
2. ✅ エンドツーエンドのテストスイート
3. ✅ パフォーマンスベンチマーク
4. ✅ 完全なドキュメント
5. ✅ チュートリアル

## 成功基準

### 必須基準

1. ✅ **忠実度**: すべての分解で > 0.9999（実際は1.0）
2. ✅ **ゲート数削減**: > 95%（目標97.5%）
3. ✅ **テスト合格率**: 100%
4. ✅ **実問題検証**: 4分子鎖で完璧に動作

### 望ましい基準

5. ✅ **計算時間**: 大幅削減
6. ✅ **コード品質**: 高い可読性と保守性
7. ✅ **ドキュメント**: 完全で正確
8. ✅ **コミュニティ**: MQT-Quditsへの貢献

## 工数とスケジュール

### 総工数見積もり

| フェーズ | 期間 | 工数（時間） |
|---------|------|-------------|
| 1. 統合実装 | 1-2週間 | 35-56 |
| 2. ゲート変換 | 2-4週間 | 60-80 |
| 3. MQT-Qudits統合 | 4-8週間 | 120-160 |
| **合計** | **3-6ヶ月** | **215-296** |

### 推奨スケジュール（フルタイム）

**月1: フェーズ1完了**
- 週1: 環境準備とTwoLevel統合
- 週2: ThreeLevel統合
- 週3: Optimizer更新とテスト
- 週4: ドキュメントとレビュー

**月2-3: フェーズ2完了**
- 月2週1-2: 2×2変換実装
- 月2週3-4: 3×3変換実装（前半）
- 月3週1-2: 3×3変換実装（後半）
- 月3週3-4: 埋め込みとテスト

**月4-6: フェーズ3完了**
- 月4: ゲートクラスとCompilerPass実装
- 月5: エンドツーエンドテストと最適化
- 月6: パフォーマンスチューニングとドキュメント

## リスク管理

### 技術的リスク

1. **Givens抽出の数値安定性**
   - 緩和策: 厳密なテストと検証
   - 代替案: 直接QRを使用

2. **MQT-Qudits互換性**
   - 緩和策: 早期の統合テスト
   - 代替案: スタンドアロン実装

3. **パフォーマンス**
   - 緩和策: プロファイリングと最適化
   - 代替案: 並列処理の導入

### プロジェクトリスク

1. **工数超過**
   - 緩和策: 段階的な実装とレビュー
   - 代替案: 範囲の調整

2. **品質問題**
   - 緩和策: 包括的なテストスイート
   - 代替案: より多くのレビュー

## まとめ

本ロードマップは、PR#37の成果を基に、完全なQuditゲート最適化
パイプラインを構築するための包括的な計画を提供しました。

### 重要なマイルストーン

1. **フェーズ1完了**: 忠実度 1.0 達成
2. **フェーズ2完了**: MQT-Quditsゲートへの変換完了
3. **フェーズ3完了**: 97.5%のゲート数削減実現

### 期待される最終成果

```
性能改善:
  - ゲート数: 6,000 → 150（97.5%削減）
  - 忠実度: 0.24/0.63 → 1.0（完璧）
  - 計算時間: 大幅削減

実用性:
  - Qubitと競争力のある性能
  - 実問題で実証済み
  - MQT-Quditsへの貢献
```

### 次のアクション

1. **即時**: フェーズ1の開始（統合実装）
2. **1-2週間後**: フェーズ2の開始（ゲート変換）
3. **2-3ヶ月後**: フェーズ3の開始（MQT-Qudits統合）

---

**文書作成日**: 2025年10月20日  
**作成者**: GitHub Copilot AI分析システム  
**バージョン**: 1.0  
**ステータス**: 実装ロードマップ完成
