# Quditゲート最適化：継続実装仕様書

## 文書の目的

本仕様書は、PR#36で開始されたquditゲート最適化の継続実装に必要な詳細仕様を提供します。
既存の分析と初期実装を基に、完全な実装に必要なすべての技術的詳細を記述します。

## 実装状況のサマリー

### 完了した作業

1. ✅ **問題分析**

   - 根本原因の特定: CustomTwoゲートの一般的分解が疎構造を無視
   - ゲート数: 現在6,182 → 理論値158（97.4%削減可能）

2. ✅ **疎構造解析器** (`tools/sparse_structure_compiler.py`)

   - H_transfer: 2×2部分空間を正確に検出
   - H_TTA: 3×3部分空間を正確に検出
   - 恒等変換の検出機能

3. ✅ **詳細ドキュメント**
   - 問題分析レポート4件
   - 実装計画書
   - 理論的基礎文書

### 未完了の作業（本PR#37で実施）

1. ⏳ **厳密なユニタリ分解器の改良**

   - 現在の実装には数値精度の問題あり
   - ZYZ/ZXZ分解の正確な実装が必要
   - Givens分解の改良が必要

2. ⏳ **MQT-Quditsフレームワークとの統合**

   - 基本ゲート（CEx, R, Rz, VirtRz）への変換
   - CompilerPassとしての実装
   - 既存回路への適用

3. ⏳ **H_transfer/H_TTA専用ゲートシーケンス**
   - 物理的に意味のある操作の直接実装
   - CustomTwoゲートを完全に回避
   - 15-35ゲートでの実装

## 技術的課題と解決方針

### 課題1: 2×2ユニタリ分解の数値精度

**問題**:

- 現在の`TwoLevelRotationDecomposer`は忠実度0.24（不合格）
- ZYZ分解の実装に数学的誤りがある可能性
- グローバル位相の処理が不適切

**解決方針**:

#### オプションA: Qiskitの実装を参考にする（推奨）

Qiskit の `TwoQubitBasisDecomposer` は高精度な2×2分解を提供しています。
同じアルゴリズムを使用することで、忠実度 > 0.9999 を達成できます。

```python
def decompose_2x2_qiskit_style(U: np.ndarray) -> Tuple[float, float, float]:
    """
    Qiskitスタイルの2×2ユニタリ分解

    参考: qiskit.quantum_info.synthesis.two_qubit_decompose

    U = e^(iα) Rz(φ) Ry(θ) Rz(λ)
    """
    # Step 1: det(U) = e^(2iα) を使ってグローバル位相を抽出
    det_U = np.linalg.det(U)
    alpha = np.angle(det_U) / 2.0

    # Step 2: SU(2)に正規化
    U_su2 = U / np.exp(1j * alpha)

    # Step 3: パラメータ抽出（Qiskitの方法）
    # U_su2 = [[a, -b*], [b, a*]] の形式（SU(2)の一般形）

    # θ: |a|から計算
    theta = 2 * np.arccos(np.clip(abs(U_su2[0, 0]), 0, 1))

    # φとλ: 複素位相から計算
    if abs(np.sin(theta / 2)) > 1e-10:
        # 一般的なケース
        exp_i_phi = U_su2[1, 0] / np.sin(theta / 2)
        exp_i_lambda = -U_su2[0, 1] / np.sin(theta / 2)
        phi = np.angle(exp_i_phi)
        lam = np.angle(exp_i_lambda)
    else:
        # θ ≈ 0 または π の場合
        phi = 0.0
        lam = np.angle(U_su2[0, 0] / U_su2[1, 1])

    return theta, phi, lam, alpha
```

**実装要件**:

- 忠実度 > 0.9999 を保証
- すべての特異点（θ=0, π）で正しく動作
- グローバル位相を正確に処理

#### オプションB: 固有値分解を使用

より数学的に厳密な方法として、固有値分解を使用できます。

```python
def decompose_2x2_via_eigenvectors(U: np.ndarray) -> Dict:
    """
    固有値分解による2×2ユニタリの分解

    理論:
    U = V @ diag(e^(iλ₁), e^(iλ₂)) @ V†

    ここで V はユニタリ（SU(2)まで）
    """
    # Step 1: 固有値分解
    eigenvalues, eigenvectors = np.linalg.eig(U)

    # Step 2: 固有値から回転角を抽出
    lambda1, lambda2 = np.angle(eigenvalues[0]), np.angle(eigenvalues[1])

    # Step 3: 固有ベクトルからパラメータを抽出
    # V をさらに分解

    # 実装の詳細は複雑だが、数学的に厳密
```

**実装要件**:

- np.linalg.eigのみ使用（厳密な線形代数）
- 特異値分解は使用しない（ヒューリスティックの可能性）

### 課題2: 3×3ユニタリ分解の数値精度

**問題**:

- 現在のGivens分解は忠実度0.63（不合格）
- 対角位相の処理が不完全
- 回転の順序が最適でない可能性

**解決方針**:

#### オプションA: 完全なCosine-Sine分解（推奨）

```python
def decompose_3x3_complete(U: np.ndarray) -> Dict:
    """
    完全な3×3ユニタリ分解

    手法:
    1. QR分解で上三角化
    2. 対角要素をユニタリ化
    3. 各Givens回転を正確に記録
    """
    rotations = []
    U_work = U.copy()

    # Phase 1: 下三角の要素を順次ゼロ化
    # G(0,1): U[1,0] → 0
    # G(0,2): U[2,0] → 0
    # G(1,2): U[2,1] → 0

    for i, j in [(1, 0), (2, 0), (2, 1)]:
        if abs(U_work[i, j]) > 1e-10:
            # Givens回転パラメータを計算
            theta, phi = compute_givens_to_zero(U_work, i, j)

            # 回転を記録
            rotations.append((i, j, theta, phi))

            # 適用
            G = construct_givens(3, i, j, theta, phi)
            U_work = G.conj().T @ U_work

    # Phase 2: 対角位相の抽出と調整
    diagonal = np.diag(U_work)
    phases = np.angle(diagonal)

    # Phase 3: 検証
    U_reconstructed = reconstruct_from_decomposition(rotations, phases)
    assert np.allclose(U, U_reconstructed, atol=1e-8)

    return {
        "rotations": rotations,
        "diagonal_phases": phases,
        "fidelity": compute_fidelity(U, U_reconstructed),
    }
```

**重要な注意点**:

- Givens回転の正確な定義を使用
- 各ステップで数値誤差をチェック
- 再構築時に同じ行列が得られることを検証

#### オプションB: Schur分解の利用

より安定した数値アルゴリズムとして、Schur分解を使用できます。

```python
def decompose_3x3_via_schur(U: np.ndarray) -> Dict:
    """
    Schur分解による3×3ユニタリの分解

    理論:
    U = Q @ T @ Q†

    ここで T は上三角、Q はユニタリ
    """
    # scipy.linalg.schur を使用
    # ただし、これはヒューリスティックではなく、
    # 厳密な線形代数アルゴリズム

    from scipy.linalg import schur

    T, Q = schur(U, output="complex")

    # Q をさらに2準位回転に分解
    # ...
```

**注意**:

- `scipy.linalg.schur` は厳密なアルゴリズム（QR法）を使用
- ヒューリスティックではない
- ただし、問題文の意図によっては避けるべき

### 課題3: MQT-Quditsフレームワークへの統合

**問題**:

- MQT-Quditsの内部APIが複雑
- 基本ゲート（CEx, R, Rz）への変換方法が不明確
- CompilerPassとしての実装方法が未定義

**解決方針**:

#### ステップ1: 既存のCRotGenを参考にする

MQT-Quditsには既に`CRotGen`というゲートシーケンス生成器があります。
これと同様の構造を使用します。

```python
# 参考: src/mqt/qudits/compiler/twodit/blocks/crot.py


class SparseStructureGateGen:
    """
    疎構造を活用したゲートシーケンス生成器

    CRotGenと同様の構造だが、疎構造に特化
    """

    def __init__(self, circuit, name, qudits, dimensions):
        self.circuit = circuit
        self.name = name
        self.qudits = qudits
        self.dimensions = dimensions

    def generate_2level_rotation(self, subspace_indices, theta, phi, lam):
        """
        2準位回転のゲートシーケンスを生成

        Args:
            subspace_indices: 作用する部分空間 [idx1, idx2]
            theta, phi, lam: ZYZ分解パラメータ

        Returns:
            ゲートのリスト
        """
        gates = []

        # Step 1: 部分空間を標準位置に移動
        gates.extend(self._prepare_subspace(subspace_indices))

        # Step 2: Rz(φ) を VirtRz で実装
        gates.append(self._create_virtrz_gate(self.qudits[1], phi))

        # Step 3: Ry(θ) を CEx, R, Rz の組み合わせで実装
        gates.extend(self._create_controlled_y_rotation(theta))

        # Step 4: Rz(λ) を VirtRz で実装
        gates.append(self._create_virtrz_gate(self.qudits[1], lam))

        # Step 5: 部分空間を元に戻す
        gates.extend(self._restore_subspace(subspace_indices))

        return gates

    def _create_controlled_y_rotation(self, theta):
        """
        制御Y回転を基本ゲートで実装

        理論:
        CRy(θ) = CEx @ Rz(-θ/2) @ CEx @ Rz(θ/2)

        ただし、実際にはより複雑な準位操作が必要
        """
        gates = []

        # フレーム調整: Ry → Rz 変換
        # Ry = Rz(-π/2) Rx Rz(π/2)
        gates.append(self._R_gate([0, 1, np.pi / 2, -np.pi / 2]))

        # 制御Exchange
        gates.append(self._CEx_gate())

        # Z回転
        gates.append(self._Rz_gate([0, 1, -theta / 2]))

        # 再びCEx
        gates.append(self._CEx_gate())

        # Z回転
        gates.append(self._Rz_gate([0, 1, theta / 2]))

        # フレーム復元
        gates.append(self._R_gate([0, 1, -np.pi / 2, np.pi / 2]))

        return gates
```

**実装要件**:

- `mqt.qudits.quantum_circuit.gates` モジュールを使用
- 各ゲートは `Gate` オブジェクトとして生成
- `circuit.append()` で回路に追加

#### ステップ2: CompilerPassの実装

```python
from mqt.qudits.compiler import CompilerPass


class SparseStructureAwareCompilerPass(CompilerPass):
    """
    疎構造認識型コンパイラパス

    CustomTwoゲートを検出し、疎構造を活用して効率的に分解
    """

    def __init__(self, backend):
        super().__init__(backend)
        self.analyzer = SparseStructureAnalyzer()

    def transpile(self, circuit):
        """
        回路を変換

        Args:
            circuit: 入力QuantumCircuit

        Returns:
            変換されたQuantumCircuit
        """
        new_instructions = []

        for gate in circuit.instructions:
            if gate.gate_type == GateTypes.TWO:
                # CustomTwoゲートを分解
                decomposed = self._decompose_custom_two(gate)
                new_instructions.extend(decomposed)
            else:
                # その他はそのまま
                new_instructions.append(gate)

        # 新しい回路を構築
        transpiled = circuit.copy()
        transpiled.instructions = new_instructions
        return transpiled

    def _decompose_custom_two(self, gate):
        """CustomTwoゲートを疎構造を活用して分解"""
        U = gate.to_matrix(identities=0)
        structure = self.analyzer.analyze(U)

        if structure.structure_type == "sparse_subspace":
            # 疎構造を活用した分解
            return self._decompose_sparse(gate, structure)
        else:
            # 一般的分解にフォールバック
            from mqt.qudits.compiler.twodit.entanglement_qr import LogEntQRCEXPass

            fallback_compiler = LogEntQRCEXPass(self.backend)
            return fallback_compiler.transpile_gate(gate)
```

**実装要件**:

- `CompilerPass` を正しく継承
- `transpile()` メソッドを実装
- 既存のコンパイラとの互換性を保持

### 課題4: H_transfer/H_TTA専用ゲートシーケンス

**問題**:

- 物理的に意味のある操作を直接実装する必要
- CustomTwoゲートを完全に回避
- MQT-Quditsの基本ゲートのみを使用

**解決方針**:

#### H_transfer専用実装

```python
class HTransferOptimizedSequence:
    """
    H_transfer時間発展の最適化実装

    物理:
    H_transfer = V (|01⟩⟨10| + |10⟩⟨01|)

    部分空間: {|01⟩, |10⟩} ⊂ ℋ^9

    ゲート数: 約15ゲート（vs 現行1000ゲート）
    """

    def __init__(self, circuit, qudits, V, dt, hbar):
        self.circuit = circuit
        self.qudits = qudits
        self.theta = V * dt / hbar

    def generate(self):
        """最適化されたゲートシーケンスを生成"""
        i, j = self.qudits
        gates = []

        # 戦略:
        # |01⟩ = |i:0⟩ ⊗ |j:1⟩
        # |10⟩ = |i:1⟩ ⊗ |j:0⟩
        #
        # この2つの状態間で回転を実装

        # Step 1: 準位マッピング
        # qudit i の準位 0,1 を使用
        # qudit j の準位 0,1 を使用

        # Step 2: 制御回転の実装
        # qudit i が |0⟩ の時、qudit j を回転
        # qudit i が |1⟩ の時、qudit j を逆回転

        # 実装（簡略版）:

        # 2a. フレーム設定（qudit jの準位0-1）
        gates.append(create_R_gate(j, [0, 1, np.pi / 2, -np.pi / 2]))

        # 2b. 制御Exchange（i-jペア）
        gates.append(create_CEx_gate([i, j]))

        # 2c. Z回転（qudit j、準位0-1、角度-θ/2）
        gates.append(create_Rz_gate(j, [0, 1, -self.theta / 2]))

        # 2d. 再びCEx
        gates.append(create_CEx_gate([i, j]))

        # 2e. Z回転（角度+θ/2）
        gates.append(create_Rz_gate(j, [0, 1, self.theta / 2]))

        # 2f. フレーム復元
        gates.append(create_R_gate(j, [0, 1, -np.pi / 2, np.pi / 2]))

        # 注: 実際には、準位2を含む場合の処理も必要
        # ここでは簡略化して記述

        return gates
```

**実装要件**:

- ゲート数 ≤ 20
- 数学的厳密性を保持
- すべてMQT-Quditsの基本ゲートのみ

#### H_TTA専用実装

```python
class HTTAOptimizedSequence:
    """
    H_TTA時間発展の最適化実装

    物理:
    H_TTA = J (|02⟩⟨11| + |11⟩⟨02| + |20⟩⟨11| + |11⟩⟨20|)

    部分空間: {|02⟩, |11⟩, |20⟩} ⊂ ℋ^9

    ゲート数: 約35ゲート（vs 現行1000ゲート）
    """

    def __init__(self, circuit, qudits, J, dt, hbar):
        self.circuit = circuit
        self.qudits = qudits
        self.J = J
        self.dt = dt
        self.hbar = hbar

    def generate(self):
        """最適化されたゲートシーケンスを生成"""
        # Step 1: 部分空間のハミルトニアンを対角化
        H_sub = self.J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])

        eigenvalues, eigenvectors = np.linalg.eigh(H_sub)

        # Step 2: 時間発展演算子を計算
        phases = np.exp(-1j * eigenvalues * self.dt / self.hbar)
        U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T

        # Step 3: U_subをGivens分解
        decomposer = RigorousThreeQuditDecomposer()
        decomp = decomposer.decompose_givens(U_sub)

        # Step 4: 各Givens回転をゲートシーケンスに変換
        gates = []

        subspace_levels = [2, 4, 6]  # |02⟩, |11⟩, |20⟩の位置

        for level1, level2, theta, phi in decomp.rotations:
            # 物理的準位のマッピング:
            # level 0 → |02⟩: (i=0, j=2)
            # level 1 → |11⟩: (i=1, j=1)
            # level 2 → |20⟩: (i=2, j=0)

            physical_level1 = subspace_levels[level1]
            physical_level2 = subspace_levels[level2]

            # この2準位間の回転を実装
            gates.extend(
                self._implement_two_level_rotation(
                    physical_level1, physical_level2, theta, phi
                )
            )

        return gates

    def _implement_two_level_rotation(self, level1, level2, theta, phi):
        """
        2つの物理準位間の回転を実装

        レベルのエンコーディング:
        level = i * 3 + j  (i, j ∈ {0, 1, 2})

        例:
        |02⟩ → level 2 → i=0, j=2
        |11⟩ → level 4 → i=1, j=1
        |20⟩ → level 6 → i=2, j=0
        """
        gates = []

        # level1, level2 から (i1, j1), (i2, j2) を計算
        i1, j1 = level1 // 3, level1 % 3
        i2, j2 = level2 // 3, level2 % 3

        # 2-qudit回転を実装
        # これは複雑だが、約12ゲートで実装可能

        # 簡略化された実装（実際はより複雑）:
        # ...

        return gates
```

**実装要件**:

- ゲート数 ≤ 40
- 固有値分解は厳密（np.linalg.eighのみ）
- scipy.linalg.expmは使用しない

## 実装の優先順位

### フェーズ1: 基礎の改良（最優先）

**目標**: 数学的厳密性を確保

1. **厳密な2×2分解の実装**

   - Qiskitの方法を参考
   - 忠実度 > 0.9999 を達成
   - 工数: 30-40時間

2. **厳密な3×3分解の実装**
   - Givens分解の完全実装
   - 対角位相の正確な処理
   - 工数: 40-50時間

**成果物**:

- `tools/unitary_decomposition_rigorous.py` の完成版
- すべてのテストケースで忠実度 > 0.9999

### フェーズ2: MQT-Qudits統合（高優先）

**目標**: 実用的なツールとして機能

3. **基本ゲートへの変換実装**

   - CEx, R, Rz, VirtRzへのマッピング
   - ゲートシーケンス生成器
   - 工数: 50-70時間

4. **CompilerPassの実装**
   - SparseStructureAwareCompilerPass
   - 既存のLogEntQRCEXPassとの統合
   - 工数: 40-60時間

**成果物**:

- `tools/mqt_qudits_integration.py`
- `tools/sparse_compiler_pass.py`

### フェーズ3: 専用シーケンスの実装（中優先）

**目標**: 最大のゲート数削減

5. **H_transfer専用実装**

   - 15ゲート以内での実装
   - 完全な検証
   - 工数: 60-80時間

6. **H_TTA専用実装**
   - 35ゲート以内での実装
   - 完全な検証
   - 工数: 80-100時間

**成果物**:

- `tools/h_transfer_optimized.py`
- `tools/h_tta_optimized.py`

### フェーズ4: 統合テストと文書化（必須）

**目標**: 品質保証と知識共有

7. **完全な統合テスト**

   - 4分子系での実測
   - ゲート数の確認（目標: 6,182 → 158）
   - 工数: 40-60時間

8. **ドキュメント完成**
   - ユーザーガイド
   - APIリファレンス
   - 実装例とチュートリアル
   - 工数: 30-40時間

**成果物**:

- `tests/test_sparse_optimization.py`
- `tutorials/doc/sparse_optimization_user_guide_ja.md`

## 総工数見積もり

| フェーズ | タスク           | 工数（時間） | 状態      |
| -------- | ---------------- | ------------ | --------- |
| 1        | 厳密な2×2分解    | 30-40        | ⏳ 未着手 |
| 1        | 厳密な3×3分解    | 40-50        | ⏳ 未着手 |
| 2        | 基本ゲート変換   | 50-70        | ⏳ 未着手 |
| 2        | CompilerPass実装 | 40-60        | ⏳ 未着手 |
| 3        | H_transfer専用   | 60-80        | ⏳ 未着手 |
| 3        | H_TTA専用        | 80-100       | ⏳ 未着手 |
| 4        | 統合テスト       | 40-60        | ⏳ 未着手 |
| 4        | ドキュメント     | 30-40        | ⏳ 未着手 |
| **合計** |                  | **370-500**  |           |

**注記**:

- 各フェーズは並行して作業可能な部分あり
- 経験豊富な開発者なら下限値が達成可能
- 完璧な品質を求めるなら上限値を想定

## 数学的厳密性の保証方法

### 禁止事項の再確認

❌ **絶対に使用してはいけないもの**:

1. `scipy.linalg.expm` - 行列指数関数（ヒューリスティック）
2. トロッター次数の削減 - 精度低下
3. 小さな行列要素の無視 - 数学的に不正確
4. 近似的な時間発展 - 物理的に不正確
5. 任意のfallback処理 - ごまかしの可能性

✅ **使用可能なもの**:

1. `np.linalg.eigh` - エルミート行列の固有値分解（厳密）
2. `np.linalg.eig` - 一般行列の固有値分解（厳密）
3. `np.linalg.qr` - QR分解（厳密）
4. `scipy.linalg.schur` - Schur分解（厳密なQR法ベース）
5. 三角関数（arccos, arctan2など）- 厳密な数学関数
6. 量子ゲートの組み合わせ - 厳密な量子操作

### 検証方法

すべての実装において、以下の検証を実施する必要があります：

```python
def verify_mathematical_rigor(original_U, decomposed_gates):
    """
    数学的厳密性を検証

    Args:
        original_U: 元のユニタリ行列
        decomposed_gates: 分解されたゲートのリスト

    Raises:
        AssertionError: 厳密性が保たれていない場合
    """
    # Step 1: ゲートから実際のユニタリを再構築
    U_reconstructed = reconstruct_unitary_from_gates(decomposed_gates)

    # Step 2: ユニタリ性の検証
    assert is_unitary(
        U_reconstructed, tolerance=1e-10
    ), "再構築されたユニタリがユニタリ性を失っています"

    # Step 3: 忠実度の計算
    fidelity = compute_fidelity(original_U, U_reconstructed)
    assert fidelity > 0.9999, f"忠実度が不十分です: {fidelity:.6f} < 0.9999"

    # Step 4: 固有値の保存確認
    eigenvals_orig = np.linalg.eigvals(original_U)
    eigenvals_recon = np.linalg.eigvals(U_reconstructed)
    eigenvals_orig_sorted = np.sort(np.angle(eigenvals_orig))
    eigenvals_recon_sorted = np.sort(np.angle(eigenvals_recon))
    assert np.allclose(
        eigenvals_orig_sorted, eigenvals_recon_sorted, atol=1e-8
    ), "固有値が保存されていません"

    print("✓ 数学的厳密性検証: 合格")
```

## 成功基準

### 最小成功基準（必須）

1. ✅ 疎構造解析が正確に動作（完了）
2. ⏳ 2×2ユニタリ分解の忠実度 > 0.9999
3. ⏳ 3×3ユニタリ分解の忠実度 > 0.9999
4. ⏳ 数学的厳密性を完全に保持

### 望ましい成功基準

5. ⏳ MQT-Quditsフレームワークへの統合完了
6. ⏳ H_transfer/H_TTAの専用実装完了
7. ⏳ ゲート数削減: 6,182 → 158 (97.4%削減)
8. ⏳ 4分子系での実測値で検証

### 理想的な成功基準

9. ⏳ 一般的な疎構造ユニタリに対応
10. ⏳ MQT-Quditsのメインブランチへの貢献
11. ⏳ 学術論文としての発表可能性
12. ⏳ 他の物理系への応用可能性

## 次のステップ

### 即時実施（この PR）

1. 本継続仕様書の作成（完了）
2. 理論的基礎の詳細化
3. 実装設計書の完成
4. テストケースの定義

### 短期（1-2週間）

5. 厳密な2×2分解の実装
6. 厳密な3×3分解の実装
7. 初期テストと検証

### 中期（1-2ヶ月）

8. MQT-Qudits統合
9. H_transfer専用実装
10. 部分的な性能測定

### 長期（3-6ヶ月）

11. H_TTA専用実装
12. 完全な統合テスト
13. ドキュメント完成
14. コミュニティへの貢献

## 参考資料

### 既存ドキュメント

1. `tutorials/doc/qudit_gate_cost_analysis.md` - 問題分析
2. `tutorials/doc/qudit_gate_optimization_implementation_plan.md` - 実装計画
3. `tutorials/doc/qudit_gate_optimization_detailed_specification.md` - 詳細仕様
4. `tutorials/doc/qudit_gate_optimization_theoretical_foundation.md` - 理論基礎

### 外部参考文献

1. **Qiskit Quantum Information Module**

   - `qiskit.quantum_info.synthesis.two_qubit_decompose`
   - 高精度な2-qubitゲート分解の実装

2. **Nielsen & Chuang "Quantum Computation and Quantum Information"**

   - Chapter 4.2: Single qubit operations
   - Chapter 4.3: Controlled operations

3. **MQT-Qudits Documentation**

   - https://mqt-qudits.readthedocs.io/
   - CompilerPassの実装方法

4. **Givens Rotations**
   - Golub & Van Loan "Matrix Computations" Chapter 5
   - 数値的に安定したユニタリ分解

## まとめ

本継続仕様書は、PR#36で開始されたquditゲート最適化の完全な実装に必要な
すべての技術的詳細を提供します。

**重要なポイント**:

- 数学的厳密性を最優先
- ヒューリスティックや近似は一切使用しない
- 段階的な実装で品質を確保
- 完全なテストと検証を実施

**推定工数**: 370-500時間
**推定期間**: 3-6ヶ月（フルタイム換算）

次のPRでは、フェーズ1（基礎の改良）の実装を開始することを推奨します。

---

**文書作成日**: 2025年10月20日
**作成者**: GitHub Copilot AI分析システム
**バージョン**: 1.0
**ステータス**: 継続仕様完成
