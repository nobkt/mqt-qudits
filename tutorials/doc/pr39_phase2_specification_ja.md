# PR#39 Phase 2詳細仕様書：ゲート変換実装

## 文書の目的

本仕様書は、PR#39のPhase 2（ゲート変換実装）の完全な詳細仕様を提供します。
Phase 1で統合した疎構造コンパイラの分解結果を、MQT-Quditsの基本ゲートに変換する
実装を定義します。

## Phase 1の成果（前提条件）

### 完了した実装

✅ **integrated_sparse_compiler.py**:
- 疎構造解析
- PR#37分解器の統合
- 2×2分解（ZYZ、忠実度 1.0）
- 3×3分解（QR→Givens、忠実度 1.0）
- ゲート数見積もり

### 現在の出力形式

**2×2分解の結果**:
```python
{
    'theta': float,      # Y軸回転角
    'phi': float,        # 第1のZ軸回転角
    'lambda': float,     # 第2のZ軸回転角
    'global_phase': float,
    'fidelity': 1.0
}
```

**3×3分解の結果**:
```python
{
    'Q': np.ndarray,     # ユニタリ行列
    'R': np.ndarray,     # 上三角行列
    'rotations': [       # Givens回転のリスト
        (level1, level2, theta, phi),
        ...
    ],
    'diagonal_phases': np.ndarray,  # 対角位相
    'fidelity': 1.0
}
```

## Phase 2の目的

### 主要目標

1. **MQT-Quditsゲートへの変換**
   - ZYZパラメータ → MQT-Quditsゲート
   - Givens回転 → MQT-Quditsゲート
   - 対角位相 → VirtRzゲート

2. **部分空間への埋め込み**
   - 2×2 → 9×9空間
   - 3×3 → 9×9空間
   - インデックスマッピング

3. **忠実度の保持**
   - 変換前後で忠実度 1.0 を維持
   - 数学的厳密性の保証

4. **ゲート数の最小化**
   - 仮想ゲート（VirtRz）のコストはゼロ
   - 物理ゲート（R, Rz, CEx）を最小化

### 期待される成果

```
入力: Phase 1の分解結果（ZYZ or Givens）
出力: MQT-Quditsゲートシーケンス

ゲート数:
  - 2×2: 3ゲート（VirtRz×2 + R×1）
  - 3×3: 12ゲート（Givens×3の3ゲート + VirtRz×3）

忠実度: 1.0（完璧に保持）
```

## MQT-Qudits基本ゲートセット

### ゲート定義

PR#38のpr38_gate_conversion_theory_ja.mdより:

1. **VirtRz（仮想Z回転）**
   ```
   VirtRz(θ, level) = |level⟩⟨level| e^(iθ) + Σ_{k≠level} |k⟩⟨k|
   コスト: 0（仮想ゲート）
   ```

2. **R（Y回転）**
   ```
   R(θ, φ, level1, level2) = I + (cos(θ/2)-1)(|level1⟩⟨level1| + |level2⟩⟨level2|)
                             + sin(θ/2)e^(iφ)(|level1⟩⟨level2| - |level2⟩⟨level1|)
   コスト: 1（物理ゲート）
   ```

3. **Rz（実Z回転）**
   ```
   Rz(θ, level) = VirtRzと同じだが物理的に実装
   コスト: 1（物理ゲート）
   注: 通常はVirtRzを使用するのでRzは不要
   ```

4. **CEx（制御交換）**
   ```
   CEx(control, target1, target2) = 制御付きスワップゲート
   コスト: 1（物理ゲート）
   注: 本実装では使用しない（2準位回転のみ）
   ```

### ゲートの数学的表現

**MQT-Quditsの2準位ユニタリ**:
```
U(level1, level2) = exp(i(α + φ/2)) × 
                    [VirtRz(α+φ, level1) × R(θ, 0, level1, level2) × VirtRz(λ, level2)]
```

これは以下と等価:
```
U = e^(i(α+φ)) Rz(α+φ)_{level1} Ry(θ)_{level1,level2} Rz(λ)_{level2}
  = VirtRz(α+φ, level1) R(θ, 0, level1, level2) VirtRz(λ, level2)
```

## 2×2ゲート変換の詳細設計

### 入力形式

```python
{
    'theta': float,      # Y軸回転角
    'phi': float,        # 第1のZ軸回転角
    'lambda': float,     # 第2のZ軸回転角
    'global_phase': float,
    'active_indices': [i, j],  # 部分空間のインデックス
    'fidelity': 1.0
}
```

### 変換アルゴリズム

```python
def convert_2x2_to_mqt_gates(params: Dict, 
                             active_indices: List[int]) -> List[MQTGate]:
    """
    ZYZパラメータをMQT-Quditsゲートに変換
    
    U = e^(iα) Rz(φ) Ry(θ) Rz(λ)
      → VirtRz(α+φ, level1) + R(θ, 0, level1, level2) + VirtRz(λ, level2)
    
    Args:
        params: Phase 1の分解結果
        active_indices: 部分空間のインデックス [i, j]
        
    Returns:
        MQT-Quditsゲートのリスト
    """
    gates = []
    
    theta = params['theta']
    phi = params['phi']
    lam = params['lambda']
    alpha = params['global_phase']
    
    level1, level2 = active_indices
    
    # VirtRz(α+φ, level1)
    phase1 = alpha + phi
    if abs(phase1) > 1e-10:
        gates.append({
            'type': 'VirtRz',
            'phase': phase1,
            'level': level1,
            'cost': 0
        })
    
    # R(θ, 0, level1, level2)
    if abs(theta) > 1e-10:
        gates.append({
            'type': 'R',
            'theta': theta,
            'phi': 0.0,  # MQT-QuditsのRゲートはφ=0を使用
            'level1': level1,
            'level2': level2,
            'cost': 1
        })
    
    # VirtRz(λ, level2)
    if abs(lam) > 1e-10:
        gates.append({
            'type': 'VirtRz',
            'phase': lam,
            'level': level2,
            'cost': 0
        })
    
    return gates
```

### 検証手順

```python
def verify_2x2_conversion(original_U: np.ndarray,
                         gates: List[MQTGate],
                         active_indices: List[int]) -> float:
    """
    変換の正確性を検証
    
    Returns:
        忠実度（1.0であるべき）
    """
    # ゲートから行列を再構築
    U_reconstructed = np.eye(2, dtype=complex)
    
    for gate in gates:
        if gate['type'] == 'VirtRz':
            # VirtRzは位相回転
            if gate['level'] == active_indices[0]:
                U_reconstructed = np.diag([np.exp(1j * gate['phase']), 1.0]) @ U_reconstructed
            else:
                U_reconstructed = np.diag([1.0, np.exp(1j * gate['phase'])]) @ U_reconstructed
        
        elif gate['type'] == 'R':
            # Rは2準位回転
            theta = gate['theta']
            c = np.cos(theta / 2)
            s = np.sin(theta / 2)
            R = np.array([[c, -s], [s, c]], dtype=complex)
            U_reconstructed = R @ U_reconstructed
    
    # 忠実度計算
    fidelity = abs(np.trace(original_U.conj().T @ U_reconstructed)) / 2.0
    return float(fidelity)
```

### テスト設計

```python
def test_2x2_gate_conversion():
    """2×2ゲート変換のテスト"""
    # ランダムな2×2ユニタリを生成
    for i in range(100):
        U = generate_random_2x2_unitary()
        
        # Phase 1: ZYZ分解
        decomposer = ImprovedTwoQubitDecomposer()
        result = decomposer.decompose_zyz(U)
        
        # Phase 2: MQT-Quditsゲート変換
        gates = convert_2x2_to_mqt_gates(result, active_indices=[0, 1])
        
        # 検証
        fidelity = verify_2x2_conversion(U, gates, [0, 1])
        assert fidelity > 0.9999, f"忠実度 {fidelity} が要求値を下回っています"
        
        # ゲート数確認
        assert len(gates) <= 3, f"ゲート数 {len(gates)} が期待値を上回っています"
    
    print("✓ 2×2ゲート変換テスト合格")
```

## 3×3ゲート変換の詳細設計

### 入力形式

```python
{
    'Q': np.ndarray,     # ユニタリ行列（3×3）
    'R': np.ndarray,     # 上三角行列（3×3）
    'rotations': [       # Givens回転のリスト
        (0, 1, theta_01, phi_01),  # G(0,1)
        (0, 2, theta_02, phi_02),  # G(0,2)
        (1, 2, theta_12, phi_12),  # G(1,2)
    ],
    'diagonal_phases': [φ0, φ1, φ2],  # 対角位相
    'active_indices': [i, j, k],       # 部分空間のインデックス
    'fidelity': 1.0
}
```

### 変換アルゴリズム

```python
def convert_3x3_to_mqt_gates(params: Dict,
                             active_indices: List[int]) -> List[MQTGate]:
    """
    Givens回転をMQT-Quditsゲートに変換
    
    U = Q R = [G(0,1) G(0,2) G(1,2)] [diag(e^(iφ0), e^(iφ1), e^(iφ2))]
    
    各Givens回転:
      G(i,j; θ, φ) → VirtRz(φ/2, i) + R(θ, 0, i, j) + VirtRz(-φ/2, j)
    
    Args:
        params: Phase 1の分解結果
        active_indices: 部分空間のインデックス [i, j, k]
        
    Returns:
        MQT-Quditsゲートのリスト
    """
    gates = []
    
    # Givens回転をゲートに変換
    for level1, level2, theta, phi in params['rotations']:
        # グローバルインデックスに変換
        global_level1 = active_indices[level1]
        global_level2 = active_indices[level2]
        
        # Givens回転のZYZ分解:
        # G(i,j; θ, φ) = Rz(φ/2)_i Ry(θ)_{i,j} Rz(-φ/2)_j
        
        alpha = phi / 2.0
        beta = -phi / 2.0
        
        # VirtRz(α, level1)
        if abs(alpha) > 1e-10:
            gates.append({
                'type': 'VirtRz',
                'phase': alpha,
                'level': global_level1,
                'cost': 0
            })
        
        # R(θ, 0, level1, level2)
        if abs(theta) > 1e-10:
            gates.append({
                'type': 'R',
                'theta': theta,
                'phi': 0.0,
                'level1': global_level1,
                'level2': global_level2,
                'cost': 1
            })
        
        # VirtRz(β, level2)
        if abs(beta) > 1e-10:
            gates.append({
                'type': 'VirtRz',
                'phase': beta,
                'level': global_level2,
                'cost': 0
            })
    
    # 対角位相をVirtRzゲートに変換
    for local_level, phase in enumerate(params['diagonal_phases']):
        global_level = active_indices[local_level]
        
        if abs(phase) > 1e-10:
            gates.append({
                'type': 'VirtRz',
                'phase': phase,
                'level': global_level,
                'cost': 0
            })
    
    return gates
```

### Givens回転の理論的根拠

**Givens行列の定義**:
```
G(i,j; θ, φ) = I + (c-1)(|i⟩⟨i| + |j⟩⟨j|) - s*|i⟩⟨j| + s|j⟩⟨i|
```
ここで:
```
c = cos(θ/2) e^(iφ/2)
s = sin(θ/2) e^(-iφ/2)
```

**ZYZ分解**:
```
G(i,j; θ, φ) = Rz(φ/2)_i Ry(θ)_{i,j} Rz(-φ/2)_j
```

**証明**（簡略版）:
```
Rz(φ/2)_i = diag(e^(iφ/2), 1) （i,j部分空間）
Ry(θ)_{i,j} = [[cos(θ/2), -sin(θ/2)],
               [sin(θ/2),  cos(θ/2)]]
Rz(-φ/2)_j = diag(1, e^(-iφ/2))

Rz(φ/2) Ry(θ) Rz(-φ/2) = [[cos(θ/2)e^(iφ/2), -sin(θ/2)e^(-iφ/2)],
                           [sin(θ/2)e^(iφ/2),   cos(θ/2)e^(-iφ/2)]]
                        = G(i,j; θ, φ)  ✓
```

### 検証手順

```python
def verify_3x3_conversion(original_U: np.ndarray,
                         gates: List[MQTGate],
                         active_indices: List[int]) -> float:
    """
    3×3変換の正確性を検証
    
    Returns:
        忠実度（1.0であるべき）
    """
    # ゲートから行列を再構築
    U_reconstructed = np.eye(3, dtype=complex)
    
    for gate in gates:
        if gate['type'] == 'VirtRz':
            # 対角行列
            D = np.eye(3, dtype=complex)
            local_level = active_indices.index(gate['level'])
            D[local_level, local_level] = np.exp(1j * gate['phase'])
            U_reconstructed = D @ U_reconstructed
        
        elif gate['type'] == 'R':
            # 2準位回転
            local_level1 = active_indices.index(gate['level1'])
            local_level2 = active_indices.index(gate['level2'])
            
            theta = gate['theta']
            c = np.cos(theta / 2)
            s = np.sin(theta / 2)
            
            R = np.eye(3, dtype=complex)
            R[local_level1, local_level1] = c
            R[local_level1, local_level2] = -s
            R[local_level2, local_level1] = s
            R[local_level2, local_level2] = c
            
            U_reconstructed = R @ U_reconstructed
    
    # 忠実度計算
    fidelity = abs(np.trace(original_U.conj().T @ U_reconstructed)) / 3.0
    return float(fidelity)
```

### テスト設計

```python
def test_3x3_gate_conversion():
    """3×3ゲート変換のテスト"""
    # ランダムな3×3ユニタリでテスト
    for i in range(100):
        U = generate_random_3x3_unitary()
        
        # Phase 1: QR分解 + Givens抽出
        decomposer = Perfect3x3Decomposer()
        result = decomposer.decompose(U)
        rotations = extract_givens_from_q(result.Q)
        phases = decomposer.extract_diagonal_phases(result.R)
        
        params = {
            'Q': result.Q,
            'R': result.R,
            'rotations': rotations,
            'diagonal_phases': phases
        }
        
        # Phase 2: MQT-Quditsゲート変換
        gates = convert_3x3_to_mqt_gates(params, active_indices=[0, 1, 2])
        
        # 検証
        fidelity = verify_3x3_conversion(U, gates, [0, 1, 2])
        assert fidelity > 0.9999, f"忠実度 {fidelity} が要求値を下回っています"
        
        # ゲート数確認
        assert len(gates) <= 12, f"ゲート数 {len(gates)} が期待値を上回っています"
    
    print("✓ 3×3ゲート変換テスト合格")
```

## 部分空間への埋め込み

### 9×9空間への埋め込み

```python
def embed_gates_in_9x9(gates: List[MQTGate],
                      active_indices: List[int],
                      total_dimension: int = 9) -> List[MQTGate]:
    """
    部分空間のゲートを完全な空間に埋め込み
    
    注意: インデックスの変換のみ。ゲートの定義は変わらない。
    
    Args:
        gates: 部分空間のゲートリスト
        active_indices: 部分空間のインデックス（例: [1, 3]）
        total_dimension: 完全な空間の次元（9）
        
    Returns:
        埋め込まれたゲートリスト
    """
    # すでにグローバルインデックスを使用しているため、
    # このステップは実際には不要。
    # ただし、検証のために実装する。
    
    for gate in gates:
        # インデックスが正しい範囲にあることを確認
        if gate['type'] == 'VirtRz':
            assert 0 <= gate['level'] < total_dimension
        elif gate['type'] == 'R':
            assert 0 <= gate['level1'] < total_dimension
            assert 0 <= gate['level2'] < total_dimension
    
    return gates  # すでに正しいインデックス
```

## 実装スケジュール

### Phase 2全体（2-4週間）

#### Week 1: 2×2変換実装
- Day 1-2: 基本クラスとデータ構造
- Day 3-4: convert_2x2_to_mqt_gates実装
- Day 5: 検証関数とテスト

#### Week 2: 3×3変換実装
- Day 1-2: Givens→MQT-Quditsゲート変換
- Day 3-4: 対角位相の処理
- Day 5: 検証関数とテスト

#### Week 3: 統合とテスト
- Day 1-2: Phase 1との統合
- Day 3: H_transfer/H_TTAでのテスト
- Day 4-5: ドキュメント作成

#### Week 4: バッファ（必要に応じて）

## 成功基準

### 必須基準

1. ✅ **忠実度**: すべての変換で > 0.9999（実際は1.0）
2. ✅ **ゲート数**: 
   - 2×2: ≤ 3ゲート
   - 3×3: ≤ 12ゲート
3. ✅ **テスト合格率**: 100%
4. ✅ **実問題検証**: H_transfer/H_TTAで正しく動作

### 望ましい基準

5. ✅ **数値安定性**: 特異点やエッジケースでも正しく動作
6. ✅ **コードの品質**: 高い可読性と保守性
7. ✅ **ドキュメント**: 完全で正確
8. ✅ **パフォーマンス**: 効率的な実装

## デリバラブル

### コード（tools/下）

1. **gate_converter.py**（予定）
   - `convert_2x2_to_mqt_gates()`
   - `convert_3x3_to_mqt_gates()`
   - `verify_conversion()`
   - 包括的なテスト

### ドキュメント（tutorials/doc/下）

1. **pr39_phase2_implementation_report_ja.md**（実装後）
   - 実装の詳細
   - テスト結果
   - 次のステップ

## 次のステップ（Phase 3）

Phase 2完了後、Phase 3に進みます:
- MQT-Quditsフレームワークへの統合
- CompilerPassの実装
- エンドツーエンドテスト

詳細は `pr39_phase3_specification_ja.md` を参照。

---

**文書作成日**: 2025年10月21日  
**作成者**: GitHub Copilot AI分析システム  
**バージョン**: 1.0  
**ステータス**: Phase 2詳細仕様
