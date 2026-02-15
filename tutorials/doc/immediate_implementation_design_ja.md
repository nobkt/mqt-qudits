# Quditゲート最適化：即時実装設計書

## 文書の目的

本文書は、PR#37で実施可能な即時実装タスクの詳細設計を提供します。
既存の問題（2×2/3×3ユニタリ分解の数値精度）を解決し、
次のフェーズへの基盤を構築することを目的とします。

## 現状の問題分析

### 既存実装の動作確認

```bash
$ python tools/sparse_structure_compiler.py

H_transfer 疎構造解析テスト
======================================================================
✓ 構造タイプ: sparse_subspace
✓ 作用する部分空間の次元: 2
✓ 作用する基底インデックス: [1, 3]
✓ ゲート数見積もり: 現在810 → 最適化後15（削減率98.1%）

H_TTA 疎構造解析テスト
======================================================================
✓ 構造タイプ: sparse_subspace
✓ 作用する部分空間の次元: 3
✓ 作用する基底インデックス: [2, 4, 6]
✓ ゲート数見積もり: 現在810 → 最適化後35（削減率95.7%）

数学的厳密性検証テスト
======================================================================
✗ 2×2ユニタリの分解と再構築:
  - 忠実度: 0.237508 （要求: > 0.9999）
  - 最大誤差: 1.76e+00
  - 厳密性: ✗ 不合格

✗ 3×3ユニタリの分解と検証:
  - 忠実度: 0.628500 （要求: > 0.9999）
  - 上三角化誤差: 1.45e+00
  - 厳密性: ✗ 不合格
```

### 問題の特定

1. **sparse_structure_compiler.py**:
   - ✅ 疎構造解析: 正常に動作
   - ✅ ゲート数見積もり: 正確
   - ❌ 2×2分解: 忠実度0.24 << 0.9999
   - ❌ 3×3分解: 忠実度0.63 << 0.9999

2. **unitary_decomposition_rigorous.py**:
   - ❌ ZYZ分解: 実装に誤りあり
   - ❌ ZXZ分解: 実装に誤りあり
   - ❌ Givens分解: 数値精度不足

### 根本原因

#### 2×2分解の問題

現在の実装（`TwoLevelRotationDecomposer.decompose_2x2_unitary`）:

```python
# 問題のある実装
theta = 2 * np.arccos(min(1.0, abs(U[0, 0])))

if abs(np.sin(theta/2)) > 1e-10:
    alpha = np.angle(-U[1, 0] / np.sin(theta/2))
    beta = np.angle(U[0, 1] / np.sin(theta/2))
```

**問題点**:
1. グローバル位相の処理が不適切
2. ZYZ分解の公式が数学的に正しくない
3. 再構築時の行列が元の行列と一致しない

#### 3×3分解の問題

現在の実装（`ThreeLevelRotationDecomposer.decompose_3x3_unitary`）:

```python
# Givens回転の構築に問題
G[level1, level1] = c * np.exp(1j * phi)
G[level1, level2] = -s
G[level2, level1] = s
G[level2, level2] = c
```

**問題点**:
1. Givens回転行列の定義が不正確
2. 対角位相の処理が不完全
3. 累積誤差が大きい

## 即時実装タスク

### タスク1: 参考実装の調査と理解

**目的**: 高精度な2×2/3×3分解の正しい実装方法を理解

**手順**:

1. **Qiskitの2-qubitゲート分解を調査**:
   ```python
   # qiskit.quantum_info.synthesis.two_qubit_decompose
   # 特に TwoQubitBasisDecomposer クラス
   
   from qiskit.quantum_info.synthesis import TwoQubitBasisDecomposer
   from qiskit.quantum_info import Operator
   
   # 任意の2×2ユニタリUに対して:
   decomposer = TwoQubitBasisDecomposer(...)
   result = decomposer(U)
   # result には theta, phi, lambda などのパラメータが含まれる
   ```

2. **NumPyの公式ドキュメントを参照**:
   - `np.linalg.eigh` の使用方法
   - `np.linalg.qr` の数値的安定性

3. **線形代数の教科書を参照**:
   - Golub & Van Loan "Matrix Computations" Chapter 5 (Givens Rotations)
   - Horn & Johnson "Matrix Analysis" Chapter 2 (Unitary Matrices)

**成果物**:
- 参考実装のサマリー文書（Markdown）
- 正しいアルゴリズムの疑似コード

**工数**: 10-15時間

### タスク2: 厳密な2×2ユニタリ分解の実装

**目的**: 忠実度 > 0.9999 を達成する2×2分解の実装

**実装方針**:

#### 方法A: Qiskitスタイルの実装（推奨）

```python
class ImprovedTwoQubitDecomposer:
    """
    Qiskitの実装を参考にした高精度2×2分解
    
    参考:
    - qiskit.quantum_info.synthesis.two_qubit_decompose
    - Nielsen & Chuang "Quantum Computation and Quantum Information" Section 4.2
    """
    
    @staticmethod
    def decompose_zyz_accurate(U: np.ndarray) -> Tuple[float, float, float, float]:
        """
        高精度ZYZ分解
        
        U = e^(iα) Rz(φ) Ry(θ) Rz(λ)
        
        Returns:
            (theta, phi, lambda, alpha)
        """
        # Step 1: グローバル位相の抽出
        det_U = np.linalg.det(U)
        alpha = np.angle(det_U) / 2.0
        
        # Step 2: SU(2)に正規化
        U_su2 = U * np.exp(-1j * alpha)
        
        # Step 3: パラメータ抽出（Qiskitの方法）
        # U_su2[0,0] = e^(i(φ+λ)/2) cos(θ/2)
        # U_su2[0,1] = -e^(i(φ-λ)/2) sin(θ/2)
        # U_su2[1,0] = e^(-i(φ-λ)/2) sin(θ/2)
        # U_su2[1,1] = e^(-i(φ+λ)/2) cos(θ/2)
        
        a = U_su2[0, 0]
        b = U_su2[0, 1]
        c = U_su2[1, 0]
        d = U_su2[1, 1]
        
        # θ: Y回転角
        # |a|² + |b|² = 1 より |a| = cos(θ/2)
        cos_theta_2 = abs(a)
        cos_theta_2 = np.clip(cos_theta_2, 0.0, 1.0)  # 数値誤差対策
        theta = 2.0 * np.arccos(cos_theta_2)
        
        # φとλの計算
        sin_theta_2 = np.sin(theta / 2.0)
        
        if sin_theta_2 > 1e-10:
            # 一般的なケース
            # a / d = e^(i(φ+λ))
            # -b / c = e^(i(φ+λ))
            # c / a = e^(-i(φ-λ)) tan(θ/2)
            # -b / d = e^(i(φ-λ)) tan(θ/2)
            
            # より安定した方法:
            # (φ+λ) / 2 = arg(a) = arg(d*)
            # (φ-λ) / 2 = arg(c) = arg(-b*)
            
            phi_plus_lambda = np.angle(a) - np.angle(d)
            phi_minus_lambda = np.angle(c) - np.angle(-b)
            
            phi = (phi_plus_lambda + phi_minus_lambda) / 2.0
            lam = (phi_plus_lambda - phi_minus_lambda) / 2.0
            
        else:
            # 特異点: θ ≈ 0 または π
            if abs(cos_theta_2 - 1.0) < 1e-10:
                # θ ≈ 0: U ≈ e^(i(φ+λ)/2) I
                phi = 0.0
                lam = np.angle(a) - np.angle(d)
            else:
                # θ ≈ π: U ≈ e^(i(φ+λ)/2) σ_y
                phi = np.angle(b) - np.angle(c)
                lam = 0.0
        
        # 正規化（-π ~ π）
        phi = np.angle(np.exp(1j * phi))
        lam = np.angle(np.exp(1j * lam))
        
        return theta, phi, lam, alpha
    
    @staticmethod
    def reconstruct_zyz(theta: float, phi: float, lam: float, 
                        alpha: float = 0.0) -> np.ndarray:
        """
        ZYZ分解からユニタリを再構築
        
        検証用: 元のUと一致することを確認
        """
        # Rz(φ)
        Rz_phi = np.array([
            [np.exp(1j * phi / 2), 0],
            [0, np.exp(-1j * phi / 2)]
        ], dtype=complex)
        
        # Ry(θ)
        Ry_theta = np.array([
            [np.cos(theta / 2), -np.sin(theta / 2)],
            [np.sin(theta / 2), np.cos(theta / 2)]
        ], dtype=complex)
        
        # Rz(λ)
        Rz_lam = np.array([
            [np.exp(1j * lam / 2), 0],
            [0, np.exp(-1j * lam / 2)]
        ], dtype=complex)
        
        # 組み合わせ
        U = Rz_phi @ Ry_theta @ Rz_lam
        U = U * np.exp(1j * alpha)
        
        return U
    
    @staticmethod
    def verify_decomposition(U_original: np.ndarray, 
                            theta: float, phi: float, lam: float, alpha: float) -> float:
        """
        分解結果を検証
        
        Returns:
            fidelity: 忠実度（> 0.9999 を期待）
        """
        U_reconstructed = ImprovedTwoQubitDecomposer.reconstruct_zyz(
            theta, phi, lam, alpha
        )
        
        # 忠実度計算
        fidelity = abs(np.trace(U_original.conj().T @ U_reconstructed)) / 2.0
        
        return fidelity
```

**テストケース**:

```python
def test_improved_2x2_decomposition():
    """改良された2×2分解のテスト"""
    from scipy.stats import unitary_group
    
    print("改良された2×2ユニタリ分解テスト")
    print("="*70)
    
    # ランダムなユニタリでテスト
    num_tests = 100
    fidelities = []
    
    for i in range(num_tests):
        U = unitary_group.rvs(2)
        
        decomposer = ImprovedTwoQubitDecomposer()
        theta, phi, lam, alpha = decomposer.decompose_zyz_accurate(U)
        fidelity = decomposer.verify_decomposition(U, theta, phi, lam, alpha)
        
        fidelities.append(fidelity)
        
        if fidelity < 0.9999:
            print(f"⚠ テスト{i}: 忠実度 = {fidelity:.6f} < 0.9999")
    
    min_fidelity = min(fidelities)
    avg_fidelity = sum(fidelities) / len(fidelities)
    
    print(f"\n結果:")
    print(f"  - 最小忠実度: {min_fidelity:.10f}")
    print(f"  - 平均忠実度: {avg_fidelity:.10f}")
    print(f"  - 合格率: {sum(1 for f in fidelities if f > 0.9999) / num_tests * 100:.1f}%")
    
    assert min_fidelity > 0.9999, f"最小忠実度が不十分: {min_fidelity}"
    print("\n✓ すべてのテストに合格")
```

**成果物**:
- `tools/unitary_decomposition_rigorous.py` の改良版
- 100個のランダムユニタリでテスト済み
- すべてのテストケースで忠実度 > 0.9999

**工数**: 30-40時間

### タスク3: 厳密な3×3ユニタリ分解の実装

**目的**: 忠実度 > 0.9999 を達成する3×3分解の実装

**実装方針**:

```python
class ImprovedThreeQuditDecomposer:
    """
    改良された3×3ユニタリ分解
    
    手法: 修正Givens分解 + 対角位相調整
    """
    
    @staticmethod
    def decompose_accurate(U: np.ndarray) -> Dict:
        """
        高精度Givens分解
        
        Args:
            U: 3×3ユニタリ行列
        
        Returns:
            {
                'rotations': [(level1, level2, theta, phi), ...],
                'diagonal_phases': [φ₀, φ₁, φ₂],
                'fidelity': float
            }
        """
        # Step 1: ユニタリ性の検証
        assert ImprovedThreeQuditDecomposer._is_unitary(U), \
            "入力行列がユニタリではありません"
        
        rotations = []
        U_work = U.copy()
        
        # Step 2: QR分解で下三角要素をゼロ化
        # G(0,1): U[1,0] → 0
        theta_01, phi_01 = ImprovedThreeQuditDecomposer._compute_givens_accurate(
            U_work[0, 0], U_work[1, 0]
        )
        if abs(theta_01) > 1e-10:
            rotations.append((0, 1, theta_01, phi_01))
            G_01 = ImprovedThreeQuditDecomposer._construct_givens_accurate(
                3, 0, 1, theta_01, phi_01
            )
            U_work = G_01.conj().T @ U_work
        
        # G(0,2): U[2,0] → 0
        theta_02, phi_02 = ImprovedThreeQuditDecomposer._compute_givens_accurate(
            U_work[0, 0], U_work[2, 0]
        )
        if abs(theta_02) > 1e-10:
            rotations.append((0, 2, theta_02, phi_02))
            G_02 = ImprovedThreeQuditDecomposer._construct_givens_accurate(
                3, 0, 2, theta_02, phi_02
            )
            U_work = G_02.conj().T @ U_work
        
        # G(1,2): U[2,1] → 0
        theta_12, phi_12 = ImprovedThreeQuditDecomposer._compute_givens_accurate(
            U_work[1, 1], U_work[2, 1]
        )
        if abs(theta_12) > 1e-10:
            rotations.append((1, 2, theta_12, phi_12))
            G_12 = ImprovedThreeQuditDecomposer._construct_givens_accurate(
                3, 1, 2, theta_12, phi_12
            )
            U_work = G_12.conj().T @ U_work
        
        # Step 3: 対角位相の抽出
        diagonal_phases = np.angle(np.diag(U_work))
        
        # Step 4: 再構築して検証
        U_reconstructed = ImprovedThreeQuditDecomposer._reconstruct(
            rotations, diagonal_phases
        )
        
        fidelity = ImprovedThreeQuditDecomposer._compute_fidelity(
            U, U_reconstructed
        )
        
        return {
            'rotations': rotations,
            'diagonal_phases': diagonal_phases,
            'fidelity': fidelity
        }
    
    @staticmethod
    def _compute_givens_accurate(a: complex, b: complex) -> Tuple[float, float]:
        """
        高精度Givensパラメータ計算
        
        目標: G(θ,φ) @ [a, b]ᵀ = [r, 0]ᵀ
        
        G = [[c,  -s*],
             [s,   c ]]
        
        c = cos(θ/2) e^(iφ/2)
        s = sin(θ/2) e^(-iφ/2)
        """
        r = np.sqrt(abs(a)**2 + abs(b)**2)
        
        if r < 1e-15:
            return 0.0, 0.0
        
        # 正規化
        a_norm = a / r
        b_norm = b / r
        
        # θの計算
        # c a_norm - s* b_norm = 1
        # |c| = cos(θ/2), |s| = sin(θ/2)
        # tan(θ/2) = |b_norm| / |a_norm|
        
        theta = 2.0 * np.arctan2(abs(b_norm), abs(a_norm))
        
        # φの計算
        if abs(b_norm) > 1e-10:
            # c / s* = (a_norm / b_norm*) * (1 / tan(θ/2))
            # arg(c) - arg(s*) = arg(a_norm) - arg(-b_norm)
            phi = np.angle(a_norm) - np.angle(-b_norm)
            phi = np.angle(np.exp(1j * phi))  # 正規化
        else:
            phi = 0.0
        
        return theta, phi
    
    @staticmethod
    def _construct_givens_accurate(d: int, level1: int, level2: int,
                                   theta: float, phi: float) -> np.ndarray:
        """
        高精度Givens行列構築
        """
        G = np.eye(d, dtype=complex)
        
        c = np.cos(theta / 2) * np.exp(1j * phi / 2)
        s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
        
        G[level1, level1] = c
        G[level1, level2] = -s.conj()
        G[level2, level1] = s
        G[level2, level2] = c.conj()
        
        # 検証: Gがユニタリであることを確認
        assert ImprovedThreeQuditDecomposer._is_unitary(G), \
            "構築されたGivens行列がユニタリではありません"
        
        return G
    
    # 他のメソッド: _reconstruct, _is_unitary, _compute_fidelity...
```

**成果物**:
- `tools/unitary_decomposition_rigorous.py` のさらなる改良
- 100個のランダム3×3ユニタリでテスト済み
- すべてのテストケースで忠実度 > 0.9999

**工数**: 40-50時間

### タスク4: 統合テストと検証

**目的**: 改良された分解器が実際の問題（H_transfer, H_TTA）で機能することを検証

**テストケース**:

```python
def test_h_transfer_with_improved_decomposer():
    """H_transferに対する改良された分解のテスト"""
    print("H_transfer 改良版分解テスト")
    print("="*70)
    
    # H_transferのユニタリを構築
    theta = 0.1
    U_9x9 = np.eye(9, dtype=complex)
    U_9x9[1, 1] = np.cos(theta)
    U_9x9[1, 3] = -1j * np.sin(theta)
    U_9x9[3, 1] = -1j * np.sin(theta)
    U_9x9[3, 3] = np.cos(theta)
    
    # 疎構造解析
    analyzer = SparseStructureAnalyzer()
    structure = analyzer.analyze(U_9x9)
    
    print(f"1. 構造解析:")
    print(f"   - タイプ: {structure.structure_type}")
    print(f"   - 部分空間次元: {structure.active_dimension}")
    
    # 部分空間抽出
    U_sub = analyzer.extract_subspace_unitary(U_9x9, structure.active_subspace)
    
    print(f"\n2. 部分空間ユニタリ:")
    print(f"   - 形状: {U_sub.shape}")
    
    # 改良された分解
    decomposer = ImprovedTwoQuditDecomposer()
    theta_zyz, phi, lam, alpha = decomposer.decompose_zyz_accurate(U_sub)
    
    print(f"\n3. ZYZ分解:")
    print(f"   - θ = {theta_zyz:.6f} (期待値: {2*theta:.6f})")
    print(f"   - φ = {phi:.6f}")
    print(f"   - λ = {lam:.6f}")
    
    # 検証
    fidelity = decomposer.verify_decomposition(U_sub, theta_zyz, phi, lam, alpha)
    
    print(f"\n4. 検証:")
    print(f"   - 忠実度 = {fidelity:.10f}")
    print(f"   - 合格: {'✓' if fidelity > 0.9999 else '✗'}")
    
    assert fidelity > 0.9999, f"忠実度が不十分: {fidelity}"
    
    print("\n✓ H_transfer改良版分解テスト完了")
    print("="*70)


def test_h_tta_with_improved_decomposer():
    """H_TTAに対する改良された分解のテスト"""
    print("H_TTA 改良版分解テスト")
    print("="*70)
    
    # H_TTAのユニタリを構築
    J = 0.05
    dt = 1.0
    hbar = 0.6582119569
    
    H_sub = J * np.array([
        [0, 1, 1],
        [1, 0, 0],
        [1, 0, 0]
    ], dtype=complex)
    
    eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
    phases = np.exp(-1j * eigenvalues * dt / hbar)
    U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T
    
    print(f"1. H_TTA部分空間ユニタリ:")
    print(f"   - 形状: {U_sub.shape}")
    
    # 改良された分解
    decomposer = ImprovedThreeQuditDecomposer()
    result = decomposer.decompose_accurate(U_sub)
    
    print(f"\n2. Givens分解:")
    print(f"   - 回転数: {len(result['rotations'])}")
    for i, (l1, l2, theta, phi) in enumerate(result['rotations']):
        print(f"     [{i}] 準位{l1}-{l2}: θ={theta:.6f}, φ={phi:.6f}")
    
    print(f"\n3. 検証:")
    print(f"   - 忠実度 = {result['fidelity']:.10f}")
    print(f"   - 合格: {'✓' if result['fidelity'] > 0.9999 else '✗'}")
    
    assert result['fidelity'] > 0.9999, f"忠実度が不十分: {result['fidelity']}"
    
    print("\n✓ H_TTA改良版分解テスト完了")
    print("="*70)
```

**成果物**:
- 実際の物理問題での検証完了
- テストスクリプト `tests/test_improved_decomposition.py`

**工数**: 20-30時間

## 実装スケジュール

### 即時実施（1週間以内）

| タスク | 工数 | 担当 | 状態 |
|--------|------|------|------|
| タスク1: 参考実装調査 | 10-15h | 実装者 | 📋 計画中 |
| タスク2: 2×2分解改良 | 30-40h | 実装者 | 📋 計画中 |

### 短期実施（2-3週間以内）

| タスク | 工数 | 担当 | 状態 |
|--------|------|------|------|
| タスク3: 3×3分解改良 | 40-50h | 実装者 | 📋 計画中 |
| タスク4: 統合テスト | 20-30h | 実装者 | 📋 計画中 |

**総工数**: 100-135時間
**推定期間**: 2-3週間（フルタイム）

## 成功基準

### 必須基準

1. ✅ 2×2ユニタリ分解の忠実度 > 0.9999（すべてのテストケース）
2. ✅ 3×3ユニタリ分解の忠実度 > 0.9999（すべてのテストケース）
3. ✅ H_transferでの動作検証完了
4. ✅ H_TTAでの動作検証完了

### 望ましい基準

5. ✅ 100個のランダムユニタリでテスト済み
6. ✅ 数値的安定性の確保（特異点での動作確認）
7. ✅ ドキュメント完備（使用方法、アルゴリズム説明）

## 次フェーズへの準備

本タスク完了後、次のフェーズ（MQT-Qudits統合）に移行できます：

1. ✅ 高精度な2×2/3×3分解が利用可能
2. ⏳ MQT-Quditsの基本ゲート（CEx, R, Rz）への変換実装
3. ⏳ CompilerPassの実装
4. ⏳ 完全な最適化パイプライン

## まとめ

本即時実装設計書は、PR#37で実施可能な具体的なタスクを定義しました。

**重点事項**:
- 数学的厳密性を最優先
- 参考実装（Qiskit等）を活用
- 段階的な検証で品質を確保
- 忠実度 > 0.9999 を厳守

**推定工数**: 100-135時間
**推定期間**: 2-3週間（フルタイム）

これらのタスクが完了すれば、次のフェーズ（MQT-Qudits統合）への
強固な基盤が確立されます。

---

**文書作成日**: 2025年10月20日  
**作成者**: GitHub Copilot AI分析システム  
**バージョン**: 1.0  
**ステータス**: 即時実装設計完成
