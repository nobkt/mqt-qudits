# PR#40 継続作業詳細仕様書

## 文書の目的

本仕様書は、PR#40で未完了の作業について、完全な技術仕様と実装ロードマップを提供します。
主に以下の2つのタスクに焦点を当てます:

1. **ゲート数最適化**: VirtRzゲートの自動結合
2. **3×3変換の修正**: Givens回転のMQT-Quditsゲートへの正確な変換

## 前提条件

### PR#40で達成したこと

✅ **gate_converter.py の基本実装**:

- TwoLevelGateConverter: 2×2変換（忠実度 1.0）
- ThreeLevelGateConverter: 3×3変換（要改善、忠実度 0.68）
- MQT-Qudits Rゲートの正しい定義を使用
- 行列積の順序を正確に実装

✅ **重要な発見**:

- MQT-Qudits R(θ, φ) ≠ 標準 Ry(θ)
- Ry(θ) = R(-θ, 0) の関係
- ZYZ分解のRz(φ)は半角位相を使用

### 継続課題

1. ⚠️ ゲート数最適化（2×2）: 5ゲート → 3ゲート
2. ⚠️ 3×3変換の忠実度: 0.68 → 1.0
3. ⚠️ ランダムユニタリテスト: 合格率 0% → 100%

## タスク 1: ゲート数最適化

### 1.1 目的

連続するVirtRzゲートを自動的に結合し、ゲート数を最小化する。

### 1.2 現状の問題

**H_transferの例**:

```python
gates = [
    VirtRz(level=1, phase=0.7854),
    VirtRz(level=3, phase=-0.7854),
    R(level1=1, level2=3, theta=-0.2000),
    VirtRz(level=1, phase=-0.7854),
    VirtRz(level=3, phase=0.7854),
]
# 総ゲート数: 5
# 物理ゲート数: 1
```

**最適化後**:

```python
# level=1: 0.7854 + (-0.7854) = 0 → 削除
# level=3: -0.7854 + 0.7854 = 0 → 削除
gates = [R(level1=1, level2=3, theta=-0.2000)]
# 総ゲート数: 1
# 物理ゲート数: 1
```

### 1.3 実装設計

#### クラス設計

```python
class GateSequenceOptimizer:
    """ゲートシーケンスの最適化器"""

    def __init__(self, tolerance: float = 1e-10):
        """
        Args:
            tolerance: 位相が実質的にゼロとみなす閾値
        """
        self.tolerance = tolerance

    def optimize(self, gates: List[MQTGate]) -> List[MQTGate]:
        """
        ゲートシーケンスを最適化

        最適化戦略:
        1. VirtRzゲートの結合: 同じレベルの連続VirtRzを1つに
        2. ゼロ位相の削除: 位相が実質的にゼロのVirtRzを削除
        3. 恒等変換の削除: R(θ≈0)などを削除

        Args:
            gates: 最適化前のゲートシーケンス

        Returns:
            最適化後のゲートシーケンス
        """
        pass

    def _combine_consecutive_virtrz(self, gates: List[MQTGate]) -> List[MQTGate]:
        """連続するVirtRzゲートを結合"""
        pass

    def _remove_identity_gates(self, gates: List[MQTGate]) -> List[MQTGate]:
        """恒等変換のゲートを削除"""
        pass
```

#### アルゴリズム

```python
def optimize(self, gates: List[MQTGate]) -> List[MQTGate]:
    """最適化のメインアルゴリズム"""

    optimized = []
    virtrz_buffer = {}  # {level: accumulated_phase}

    for gate in gates:
        if gate.gate_type == "VirtRz":
            # VirtRzゲート: バッファに累積
            level = gate.parameters["level"]
            phase = gate.parameters["phase"]
            virtrz_buffer[level] = virtrz_buffer.get(level, 0.0) + phase

        else:
            # 非VirtRzゲート: バッファをフラッシュ
            self._flush_virtrz_buffer(optimized, virtrz_buffer)

            # 恒等変換でない場合のみ追加
            if not self._is_identity_gate(gate):
                optimized.append(gate)

    # 最後のバッファをフラッシュ
    self._flush_virtrz_buffer(optimized, virtrz_buffer)

    return optimized


def _flush_virtrz_buffer(
    self, optimized: List[MQTGate], virtrz_buffer: Dict[int, float]
):
    """VirtRzバッファを出力してクリア"""
    for level, phase in virtrz_buffer.items():
        # 位相を[-π, π]に正規化
        phase = np.angle(np.exp(1j * phase))

        # ゼロでない位相のみ出力
        if abs(phase) > self.tolerance:
            optimized.append(
                MQTGate(
                    gate_type="VirtRz",
                    parameters={"level": level, "phase": phase},
                    cost=0,
                )
            )

    virtrz_buffer.clear()


def _is_identity_gate(self, gate: MQTGate) -> bool:
    """ゲートが恒等変換かどうかを判定"""
    if gate.gate_type == "R":
        theta = gate.parameters["theta"]
        return abs(theta) < self.tolerance

    return False
```

### 1.4 テスト設計

```python
def test_virtrz_combination():
    """VirtRz結合のテスト"""

    # テストケース1: 完全キャンセル
    gates = [
        MQTGate("VirtRz", {"level": 0, "phase": 0.5}, 0),
        MQTGate("VirtRz", {"level": 0, "phase": -0.5}, 0),
    ]
    optimized = optimizer.optimize(gates)
    assert len(optimized) == 0  # すべてキャンセル

    # テストケース2: 部分結合
    gates = [
        MQTGate("VirtRz", {"level": 0, "phase": 0.5}, 0),
        MQTGate("VirtRz", {"level": 0, "phase": 0.3}, 0),
        MQTGate("R", {"level1": 0, "level2": 1, "theta": 0.2, "phi": 0.0}, 1),
        MQTGate("VirtRz", {"level": 0, "phase": 0.1}, 0),
    ]
    optimized = optimizer.optimize(gates)
    assert len(optimized) == 3  # VirtRz(0.8) + R + VirtRz(0.1)

    # テストケース3: H_transfer
    gates = [
        MQTGate("VirtRz", {"level": 1, "phase": 0.7854}, 0),
        MQTGate("VirtRz", {"level": 3, "phase": -0.7854}, 0),
        MQTGate("R", {"level1": 1, "level2": 3, "theta": -0.2, "phi": 0.0}, 1),
        MQTGate("VirtRz", {"level": 1, "phase": -0.7854}, 0),
        MQTGate("VirtRz", {"level": 3, "phase": 0.7854}, 0),
    ]
    optimized = optimizer.optimize(gates)
    assert len(optimized) == 1  # Rゲートのみ

    print("✓ VirtRz結合テスト合格")


def test_h_transfer_optimization():
    """H_transferの完全な最適化テスト"""
    from improved_unitary_decomposition import ImprovedTwoQubitDecomposer

    # H_transfer 2×2ユニタリ
    theta = 0.1
    U = np.array(
        [[np.cos(theta), -1j * np.sin(theta)], [-1j * np.sin(theta), np.cos(theta)]],
        dtype=complex,
    )

    # Phase 1: ZYZ分解
    decomposer = ImprovedTwoQubitDecomposer()
    result = decomposer.decompose_zyz(U)

    # Phase 2: ゲート変換
    converter = TwoLevelGateConverter()
    params = {
        "theta": result.theta,
        "phi": result.phi,
        "lambda": result.lam,
        "global_phase": result.global_phase,
    }
    gates = converter.convert(params, [1, 3])

    # Phase 2.5: 最適化
    optimizer = GateSequenceOptimizer()
    optimized_gates = optimizer.optimize(gates.gates)

    # 検証
    assert len(optimized_gates) <= 3  # 最大3ゲート
    assert sum(1 for g in optimized_gates if g.cost > 0) == 1  # 物理ゲート1個

    # 忠実度検証
    optimized_sequence = MQTGateSequence(optimized_gates)
    fidelity = converter.verify_conversion(U, optimized_sequence, [1, 3])
    assert fidelity > 0.9999

    print(f"✓ H_transfer最適化テスト合格")
    print(f"  元のゲート数: {len(gates.gates)}")
    print(f"  最適化後: {len(optimized_gates)}")
    print(f"  忠実度: {fidelity:.10f}")
```

### 1.5 統合

```python
class TwoLevelGateConverter:
    """更新版: 最適化機能を統合"""

    def __init__(self, tolerance: float = 1e-10, optimize: bool = True):
        """
        Args:
            tolerance: 数値誤差の許容範囲
            optimize: ゲートシーケンスを自動最適化するか
        """
        self.tolerance = tolerance
        self.optimize_flag = optimize
        self.optimizer = GateSequenceOptimizer(tolerance) if optimize else None

    def convert(self, params: Dict, active_indices: List[int]) -> MQTGateSequence:
        """ZYZ分解結果をMQT-Quditsゲートに変換（最適化付き）"""

        # 基本的な変換
        gates = self._convert_basic(params, active_indices)

        # 最適化
        if self.optimize_flag and self.optimizer:
            gates = self.optimizer.optimize(gates)

        return MQTGateSequence(gates=gates, fidelity=1.0)
```

### 1.6 期待される結果

**H_transfer**:

- 元: 5ゲート（1物理 + 4仮想）
- 最適化後: 1ゲート（1物理）
- 忠実度: 1.0（維持）

**一般的な2×2ユニタリ**:

- 元: 5ゲート
- 最適化後: 1-3ゲート（状況による）
- 忠実度: 1.0（維持）

## タスク 2: 3×3変換の修正

### 2.1 目的

Givens回転のMQT-Quditsゲートへの正確な変換を実現し、忠実度 1.0 を達成する。

### 2.2 問題の診断

**現状**:

- H_TTA: 忠実度 0.6795
- ランダムユニタリ: 平均忠実度 0.3612

**可能性のある原因**:

1. **Givens回転の定義の誤解**

   - givens_rotation_theory_ja.md の定義
   - integrated_sparse_compiler.py の抽出方法
   - 両者の一致を確認する必要

2. **MQT-Qudits Rゲートでの実装誤り**

   - 2×2の場合と同様の符号問題
   - 位相の扱い
   - 行列積の順序

3. **対角位相の扱い**
   - QR分解のR行列からの抽出
   - MQT-QuditsゲートのGivens回転 + 同一レベルのVirtRzの順序

### 2.3 診断手順

#### Step 1: Givens回転の定義を確認

```python
def test_givens_definition():
    """Givens回転の定義を確認"""

    # givens_rotation_theory_ja.md の定義
    def construct_givens_theory(i: int, j: int, theta: float, phi: float) -> np.ndarray:
        """理論的なGivens行列"""
        c = np.cos(theta / 2) * np.exp(1j * phi / 2)
        s = np.sin(theta / 2) * np.exp(-1j * phi / 2)

        G = np.eye(3, dtype=complex)
        G[i, i] = c
        G[i, j] = s
        G[j, i] = -np.conj(s)
        G[j, j] = np.conj(c)

        return G

    # integrated_sparse_compiler.py の抽出結果と比較
    from integrated_sparse_compiler import IntegratedThreeLevelDecomposer

    # テスト用3×3ユニタリ
    U = np.array([[0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0, 0.0, 1.0]], dtype=complex)

    # QR分解
    Q, R = np.linalg.qr(U)

    # Givens抽出
    decomposer = IntegratedThreeLevelDecomposer()
    rotations = decomposer._extract_givens_from_q(Q)

    # 再構築して比較
    Q_reconstructed = np.eye(3, dtype=complex)
    for i, j, theta, phi in rotations:
        G = construct_givens_theory(i, j, theta, phi)
        Q_reconstructed = G @ Q_reconstructed

    fidelity = abs(np.trace(Q.conj().T @ Q_reconstructed)) / 3.0
    print(f"Givens定義検証: 忠実度 = {fidelity:.10f}")
    assert fidelity > 0.9999
```

#### Step 2: ZYZ分解の対応を確認

```python
def test_givens_zyz_decomposition():
    """Givens回転のZYZ分解を確認"""

    # Givens回転: G(i,j; θ, φ)
    # ZYZ分解: G = Rz(φ/2)_i Ry(θ)_{i,j} Rz(-φ/2)_j

    theta = 0.5
    phi = 1.0

    # 理論的なGivens
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    G_theory = np.array(
        [[c, s, 0], [-np.conj(s), np.conj(c), 0], [0, 0, 1]], dtype=complex
    )

    # ZYZ分解
    Rz_phi_2 = np.diag([np.exp(1j * phi / 2), 1.0, 1.0])

    # Ry(θ) on levels (0,1)
    cos_t2 = np.cos(theta / 2)
    sin_t2 = np.sin(theta / 2)
    Ry = np.array([[cos_t2, -sin_t2, 0], [sin_t2, cos_t2, 0], [0, 0, 1]], dtype=complex)

    Rz_minus_phi_2 = np.diag([1.0, np.exp(-1j * phi / 2), 1.0])

    G_zyz = Rz_phi_2 @ Ry @ Rz_minus_phi_2

    fidelity = abs(np.trace(G_theory.conj().T @ G_zyz)) / 3.0
    print(f"Givens ZYZ分解検証: 忠実度 = {fidelity:.10f}")

    if fidelity < 0.9999:
        print("ERROR: ZYZ分解が正しくありません")
        print(f"G_theory:\n{G_theory}")
        print(f"G_zyz:\n{G_zyz}")
        return False

    return True
```

#### Step 3: MQT-Quditsゲートへの変換を確認

```python
def test_givens_to_mqt_gates():
    """Givens回転のMQT-Quditsゲート変換を確認"""

    theta = 0.5
    phi = 1.0

    # 理論的なGivens
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    G_theory = np.array(
        [[c, s, 0], [-np.conj(s), np.conj(c), 0], [0, 0, 1]], dtype=complex
    )

    # MQT-Quditsゲート変換
    # 重要: 2×2の場合と同様に符号反転が必要かも
    gates = [
        MQTGate("VirtRz", {"level": 0, "phase": phi / 2}, 0),
        MQTGate(
            "R", {"level1": 0, "level2": 1, "theta": -theta, "phi": 0.0}, 1
        ),  # 符号反転?
        MQTGate("VirtRz", {"level": 1, "phase": -phi / 2}, 0),
    ]

    # 再構築（逆順で左から掛ける）
    G_mqt = np.eye(3, dtype=complex)
    for gate in reversed(gates):
        if gate.gate_type == "VirtRz":
            level = gate.parameters["level"]
            phase = gate.parameters["phase"]
            D = np.eye(3, dtype=complex)
            D[level, level] = np.exp(1j * phase)
            G_mqt = D @ G_mqt

        elif gate.gate_type == "R":
            level1 = gate.parameters["level1"]
            level2 = gate.parameters["level2"]
            theta_r = gate.parameters["theta"]

            # MQT-Qudits R gate
            c_r = np.cos(theta_r / 2)
            s_r = np.sin(theta_r / 2)

            R = np.eye(3, dtype=complex)
            R[level1, level1] = c_r
            R[level1, level2] = s_r
            R[level2, level1] = -s_r
            R[level2, level2] = c_r

            G_mqt = R @ G_mqt

    fidelity = abs(np.trace(G_theory.conj().T @ G_mqt)) / 3.0
    print(f"Givens → MQT-Quditsゲート検証: 忠実度 = {fidelity:.10f}")

    if fidelity < 0.9999:
        print("ERROR: MQT-Quditsゲート変換が正しくありません")
        print(f"G_theory:\n{G_theory}")
        print(f"G_mqt:\n{G_mqt}")

        # 符号反転を試す
        print("\n符号反転を試行...")
        gates_alt = [
            MQTGate("VirtRz", {"level": 0, "phase": phi / 2}, 0),
            MQTGate(
                "R", {"level1": 0, "level2": 1, "theta": theta, "phi": 0.0}, 1
            ),  # 符号そのまま
            MQTGate("VirtRz", {"level": 1, "phase": -phi / 2}, 0),
        ]

        G_mqt_alt = reconstruct_from_gates_3x3(gates_alt, [0, 1, 2])
        fidelity_alt = abs(np.trace(G_theory.conj().T @ G_mqt_alt)) / 3.0
        print(f"符号そのままの忠実度: {fidelity_alt:.10f}")

        return False

    return True
```

### 2.4 修正実装

診断結果に基づいて、`ThreeLevelGateConverter.convert()` を修正:

```python
def convert(self, params: Dict, active_indices: List[int]) -> MQTGateSequence:
    """
    3×3 Givens分解結果をMQT-Quditsゲートに変換（修正版）

    修正内容:
    1. Givens回転の定義を再確認
    2. MQT-Qudits Rゲートの符号を確認
    3. 対角位相の適用順序を確認
    """
    gates = []

    # Givens回転をゲートに変換
    rotations = params.get("rotations", [])
    for local_level1, local_level2, theta, phi in rotations:
        global_level1 = active_indices[local_level1]
        global_level2 = active_indices[local_level2]

        # Givens回転のZYZ分解:
        # G(i,j; θ, φ) = Rz(φ/2)_i Ry(θ)_{i,j} Rz(-φ/2)_j
        #
        # 重要: 診断結果に基づいて符号を決定
        # もし test_givens_to_mqt_gates() で符号反転が必要と判明したら、
        # theta_mqt = -theta を使用

        alpha = phi / 2.0
        beta = -phi / 2.0
        theta_mqt = -theta  # または theta（診断結果による）

        if abs(alpha) > self.tolerance:
            gates.append(
                MQTGate(
                    gate_type="VirtRz",
                    parameters={"level": global_level1, "phase": alpha},
                    cost=0,
                )
            )

        if abs(theta) > self.tolerance:
            gates.append(
                MQTGate(
                    gate_type="R",
                    parameters={
                        "level1": global_level1,
                        "level2": global_level2,
                        "theta": theta_mqt,
                        "phi": 0.0,
                    },
                    cost=1,
                )
            )

        if abs(beta) > self.tolerance:
            gates.append(
                MQTGate(
                    gate_type="VirtRz",
                    parameters={"level": global_level2, "phase": beta},
                    cost=0,
                )
            )

    # 対角位相をVirtRzゲートに変換
    diagonal_phases = params.get("diagonal_phases", [])
    for local_level, phase in enumerate(diagonal_phases):
        global_level = active_indices[local_level]

        if abs(phase) > self.tolerance:
            gates.append(
                MQTGate(
                    gate_type="VirtRz",
                    parameters={"level": global_level, "phase": phase},
                    cost=0,
                )
            )

    return MQTGateSequence(gates=gates, fidelity=1.0)
```

### 2.5 包括的なテスト

```python
def test_3x3_comprehensive():
    """3×3変換の包括的テスト"""

    print("=" * 70)
    print("3×3変換包括的テスト")
    print("=" * 70)

    # テスト1: 単一のGivens回転
    print("\nテスト1: 単一のGivens回転")
    theta, phi = 0.5, 1.0
    G_single = construct_givens_3x3(0, 1, theta, phi)
    gates_single = convert_givens_to_gates(0, 1, theta, phi)
    fidelity_single = verify_3x3_conversion(G_single, gates_single, [0, 1, 2])
    print(f"  忠実度: {fidelity_single:.10f}")
    assert fidelity_single > 0.9999

    # テスト2: H_TTA
    print("\nテスト2: H_TTA")
    U_tta = construct_h_tta()
    result_tta = convert_3x3_full(U_tta)
    fidelity_tta = verify_3x3_conversion(U_tta, result_tta.gates, [0, 1, 2])
    print(f"  忠実度: {fidelity_tta:.10f}")
    assert fidelity_tta > 0.9999

    # テスト3: ランダムユニタリ（100個）
    print("\nテスト3: ランダムユニタリ")
    fidelities = []
    for i in range(100):
        U_random = generate_random_3x3_unitary()
        result_random = convert_3x3_full(U_random)
        fidelity = verify_3x3_conversion(U_random, result_random.gates, [0, 1, 2])
        fidelities.append(fidelity)

    min_fid = min(fidelities)
    avg_fid = np.mean(fidelities)
    pass_count = sum(1 for f in fidelities if f > 0.9999)

    print(f"  最小忠実度: {min_fid:.10f}")
    print(f"  平均忠実度: {avg_fid:.10f}")
    print(f"  合格率: {pass_count}/100 ({pass_count}%)")

    assert pass_count == 100

    print("\n✓✓✓ すべてのテストに合格")
```

### 2.6 期待される結果

**単一のGivens回転**:

- 忠実度: 1.0

**H_TTA**:

- 元: 忠実度 0.68
- 修正後: 忠実度 1.0
- ゲート数: 12（変わらず）

**ランダムユニタリテスト**:

- 元: 合格率 0%
- 修正後: 合格率 100%
- 平均忠実度: 1.0

## 実装スケジュール

### Week 1: ゲート数最適化（2-3日）

**Day 1**:

- GateSequenceOptimizerクラスの実装
- VirtRz結合アルゴリズム
- 基本的なテスト

**Day 2**:

- TwoLevelGateConverterへの統合
- H_transferでの検証
- ランダムユニタリでの検証

**Day 3**:

- ThreeLevelGateConverterへの統合（3×3修正後）
- 包括的なテスト
- ドキュメント更新

### Week 2: 3×3変換の修正（3-5日）

**Day 1-2**:

- 診断テストの実装と実行
  - test_givens_definition()
  - test_givens_zyz_decomposition()
  - test_givens_to_mqt_gates()
- 問題の特定

**Day 3-4**:

- ThreeLevelGateConverter.convert()の修正
- 単一Givens回転での検証
- H_TTAでの検証

**Day 5**:

- ランダムユニタリテスト（100個）
- 包括的なテスト
- ドキュメント更新

### Week 3: 統合とドキュメント（2-3日）

**Day 1**:

- すべての変更の統合
- エンドツーエンドテスト
- パフォーマンス測定

**Day 2**:

- ドキュメントの完成
  - gate_converter.py のdocstring
  - チュートリアルの追加
  - README.mdの更新

**Day 3**:

- PR#40の最終化
- レビュー対応
- Phase 3への準備

## 成功基準

### 必須基準

1. ✅ **ゲート数最適化**:

   - H_transfer: 5ゲート → 1-3ゲート
   - 忠実度: 1.0 を維持

2. ✅ **3×3変換の忠実度**:

   - H_TTA: 1.0
   - ランダムユニタリ: すべて > 0.9999
   - 合格率: 100%

3. ✅ **数学的厳密性**:
   - ヒューリスティックゼロ
   - 近似ゼロ
   - すべて厳密な線形代数

### 望ましい基準

4. ✅ **パフォーマンス**:

   - 高速な変換
   - 効率的な最適化

5. ✅ **コードの品質**:

   - 高い可読性
   - 完全なドキュメント
   - 包括的なテスト

6. ✅ **実用性**:
   - 即座に使用可能
   - 他のツールとの統合が容易

## 次のステップ（Phase 3）

Phase 2完了後、Phase 3（MQT-Quditsフレームワーク統合）に移行:

1. **MQTGateSequenceクラスの拡張**

   - 実際のMQT-Qudits回路への適用
   - ゲートの追加機能

2. **SparseStructureOptimizationPassの実装**

   - CustomTwoゲートの検出
   - 疎構造解析とゲート変換の統合
   - CompilerPassとしての実装

3. **エンドツーエンドテスト**

   - 4分子鎖の完全な回路
   - 100ステップのトロッター分解
   - ゲート数削減の実測

4. **最終目標**:
   - ゲート数: 6,000 → 150 (97.5%削減)
   - 忠実度: 1.0 (完璧)
   - Qubitと競争力のある性能

---

**文書作成日**: 2025年10月21日
**作成者**: GitHub Copilot AI分析システム
**バージョン**: 1.0
**ステータス**: PR#40継続作業詳細仕様
