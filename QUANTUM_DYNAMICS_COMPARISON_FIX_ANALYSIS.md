# 量子ダイナミクス完全比較ノートブック修正分析

## 問題の特定

### 1. 理論文書の分析結果

`tutorials/doc/theory_quantum_dynamics_complete_comparison.md`を再分析した結果、以下の重要な点を確認:

#### 1.1 H_TTA ユニタリ行列の正しい構造

理論文書 (Section 7.8.3, lines 3104-3107) より:

```
U_TTA^{subspace}(t) = [
    [ 0.5*(1+cos(ω)),  -i*sin(ω)/√2,  -0.5*(1-cos(ω))],
    [-i*sin(ω)/√2,     cos(ω),        -i*sin(ω)/√2   ],
    [-0.5*(1-cos(ω)),  -i*sin(ω)/√2,   0.5*(1+cos(ω))]
]
```

ここで `ω = √2·J·t/ℏ`

**重要な性質**:
- **対角要素は実数**: `0.5*(1+cos(ω))` と `cos(ω)`
- **非対角要素は純虚数**: `-i*sin(ω)/√2`
- **U[0,2]とU[2,0]は負**: `-0.5*(1-cos(ω))`

理論文書 lines 3038-3041 にも明記:

```
**注意**: 以前のバージョンの文書では、オフ対角要素を実数（$\frac{\sqrt{2}\sin\omega}{2}$）と誤って
記載していたが、これは**非ユニタリ行列**となり誤りである。正しくは純虚数でなければならない。
```

#### 1.2 scipy.linalg.expm による厳密計算の必要性

理論文書 lines 3043-3062 の検証コードより:

```python
from scipy.linalg import expm
import numpy as np

# パラメータ
J = 0.05  # eV
dt = 10.0  # fs
hbar = 0.6582  # eV·fs

# ハミルトニアン
H_TTA = J * np.array([[0, 1, 0],
                      [1, 0, 1],
                      [0, 1, 0]])

# 厳密なユニタリ行列
U_exact = expm(-1j * H_TTA * dt / hbar)

# ユニタリ性検証: ||U†U - I|| < 10^-15
error = np.linalg.norm(U_exact @ U_exact.conj().T - np.eye(3))
print(f"ユニタリ性誤差: {error:.2e}")  # ~10^-16
```

**この方法が唯一の正しいアプローチ**:
- 解析的な公式ではなく、scipy.linalg.expmを使用する
- これにより数値誤差を機械精度（~10^-15）に抑える
- ヒューリスティックや近似は一切含まない

#### 1.3 per-pair Trotter分解の統一性

理論文書 Section 3.8 (lines 597-618) より:

**重要**: 古典シミュレーション、Qubitシミュレーション、Quditシミュレーションの全てで、
同じper-pair分解を使用することで、公平な比較を実現している。

すなわち、`H_transfer`と`H_TTA`については:

```
exp(-i H_transfer t) ≈ ∏_{pairs} exp(-i H_transfer^{pair} t)
exp(-i H_TTA t) ≈ ∏_{pairs} exp(-i H_TTA^{pair} t)
```

という形で、ペアごとに演算子を逐次適用する。

**Trotter分解の順序** (理論文書 lines 456-459):

```
順序: H0(forward) -> H_transfer(forward) -> H_TTA(forward) ->
      H_TTA(backward) -> H_transfer(backward) -> H0(backward)
```

### 2. 実装ファイルの検証

#### 2.1 exact_qudit_basic_gates.py の問題

`tutorials/exact_qudit_basic_gates.py` lines 68-167:

```python
def apply_H_TTA_basic_gates(circuit, qudit_i: int, qudit_j: int,
                             J: float, dt: float, hbar: float):
    """
    Apply TTA Hamiltonian evolution using EXACT unitary from scipy.linalg.expm.
    ...
    """
    # Calculate fundamental parameter
    omega = np.sqrt(2) * J * dt / hbar
    
    # Define the Hamiltonian in the 3D subspace
    H_TTA = J * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    
    # Compute EXACT 3×3 unitary using matrix exponential
    U_3x3 = expm(-1j * H_TTA * dt / hbar)
    
    # CRITICAL: Verify unitarity - this is required by PR#86 fix
    unitarity_error = np.linalg.norm(U_3x3 @ U_3x3.conj().T - np.eye(3))
    if unitarity_error > 1e-10:
        raise ValueError(f"H_TTA unitary is not unitary! Error: {unitarity_error:.2e}")
    
    # Verify against analytical formula
    ...
```

**この実装は正しい**: scipy.linalg.expmを使用し、厳密なユニタリ行列を計算している。

#### 2.2 exact_qubit_hamiltonians.py の問題

`tutorials/exact_qubit_hamiltonians.py` lines 70-108:

```python
def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """
    Build exact H_TTA unitary for 4-qubit system (2 molecules).
    ...
    """
    H = np.zeros((16, 16), dtype=complex)
    
    # |T1⟩_i|T1⟩_j = |01⟩_i|01⟩_j = |0101⟩
    # Little-endian: q3=0, q2=1, q1=0, q0=1
    idx_01_01 = 0b0101  # = 5
    
    # |S0⟩_i|S1⟩_j = |00⟩_i|10⟩_j = |0010⟩
    # Little-endian: q3=0, q2=0, q1=1, q0=0
    idx_00_10 = 0b0010  # = 2
    
    # H_TTA = J(|0010⟩⟨0101| + |0101⟩⟨0010|)
    # This implements: J(|S0⟩_i|S1⟩_j⟨T1|_i⟨T1|_j + h.c.)
    H[idx_00_10, idx_01_01] = J
    H[idx_01_01, idx_00_10] = J
    
    # Time evolution
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    
    return U
```

**問題**: この実装は**不完全**である。

H_TTAは2つの独立な2D部分空間で作用する:
1. `{|0101⟩, |0010⟩}`: `|T1⟩_i|T1⟩_j ↔ |S0⟩_i|S1⟩_j`
2. `{|0101⟩, |1000⟩}`: `|T1⟩_i|T1⟩_j ↔ |S1⟩_i|S0⟩_j`

現在の実装は最初の遷移のみを実装しており、2番目の遷移（`|S1⟩_i|S0⟩_j`への遷移）が欠けている。

**正しい実装**:

```python
def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    H = np.zeros((16, 16), dtype=complex)
    
    # State 1: |T1⟩_i|T1⟩_j = |0101⟩
    idx_01_01 = 0b0101  # = 5
    
    # State 2a: |S0⟩_i|S1⟩_j = |0010⟩  
    idx_00_10 = 0b0010  # = 2
    
    # State 2b: |S1⟩_i|S0⟩_j = |1000⟩
    idx_10_00 = 0b1000  # = 8
    
    # H_TTA couples:
    # |T1 T1⟩ ↔ |S0 S1⟩
    H[idx_00_10, idx_01_01] = J
    H[idx_01_01, idx_00_10] = J
    
    # |T1 T1⟩ ↔ |S1 S0⟩  
    H[idx_10_00, idx_01_01] = J
    H[idx_01_01, idx_10_00] = J
    
    # Time evolution
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    
    return U
```

これにより、H_TTAは3次元部分空間 `{|0010⟩, |0101⟩, |1000⟩}` で作用するようになる。

### 3. ノートブック内の古典シミュレータの問題

`quantum_dynamics_complete_comparison.ipynb` Cell 3 (lines 252-525) の `ClassicalSuzukiTrotterSimulator`:

#### 3.1 per-pair分解の実装

現在の実装 (lines 414-441) は正しい:

```python
# H0: per-molecule unitaries
U_H0_half_list = []
for mol_idx in range(self.N):
    H0_mol = self.build_H0_single_molecule(mol_idx)
    U = scipy.linalg.expm(-1j * H0_mol * dt / (2 * self.params.hbar))
    U_H0_half_list.append(U)

# H_transfer: per-pair unitaries
U_transfer_half_list = []
for mol_i, mol_j in self.params.neighbors:
    H_tr_pair = self.build_H_transfer_pair(mol_i, mol_j)
    U = scipy.linalg.expm(-1j * H_tr_pair * dt / (2 * self.params.hbar))
    U_transfer_half_list.append(U)

# H_TTA: per-pair unitaries
U_TTA_half_list = []
for mol_i, mol_j in self.params.neighbors:
    H_TTA_pair = self.build_H_TTA_pair(mol_i, mol_j)
    U = scipy.linalg.expm(-1j * H_TTA_pair * dt / (2 * self.params.hbar))
    U_TTA_half_list.append(U)
```

#### 3.2 適用順序の実装

現在の実装 (lines 459-485) は正しい:

```python
# Forward: H0, H_transfer, H_TTA (same order as quantum)
# H0: forward
for mol_idx in range(self.N):
    state = U_H0_half_list[mol_idx] @ state

# H_transfer: forward
for pair_idx in range(len(self.params.neighbors)):
    state = U_transfer_half_list[pair_idx] @ state

# H_TTA: forward
for pair_idx in range(len(self.params.neighbors)):
    state = U_TTA_half_list[pair_idx] @ state

# Backward: reverse order (same as quantum)
# H_TTA: backward
for pair_idx in reversed(range(len(self.params.neighbors))):
    state = U_TTA_half_list[pair_idx] @ state

# H_transfer: backward
for pair_idx in reversed(range(len(self.params.neighbors))):
    state = U_transfer_half_list[pair_idx] @ state

# H0: backward
for mol_idx in reversed(range(self.N)):
    state = U_H0_half_list[mol_idx] @ state
```

### 4. 根本的な問題

上記の分析から、以下の問題が特定された:

1. **Qubit実装のH_TTA**:
   - `exact_qubit_hamiltonians.py`の`build_H_TTA_qubit_unitary`が不完全
   - `|S1⟩_i|S0⟩_j`への遷移が欠けている
   - これにより、Qubitシミュレーションが古典シミュレーションと一致しない

2. **exact_qubit_basic_gates.pyの依存**:
   - `apply_H_TTA_basic_gates`は`build_H_TTA_qubit_unitary`を使用している
   - したがって、この関数も不完全なH_TTAを実装している

3. **ノートブック内のQubitシミュレータ**:
   - `exact_qubit_basic_gates.py`を使用しているため、上記の問題を引き継いでいる

### 5. 修正方針

#### 5.1 優先度1: exact_qubit_hamiltonians.pyの修正

`build_H_TTA_qubit_unitary`関数を修正し、完全なH_TTA実装を提供する:

```python
def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """
    Build exact H_TTA unitary for 4-qubit system (2 molecules).
    
    H_TTA couples THREE states in the full TTA process:
    1. |T1⟩_i|T1⟩_j ↔ |S0⟩_i|S1⟩_j  (idx 5 ↔ 2)
    2. |T1⟩_i|T1⟩_j ↔ |S1⟩_i|S0⟩_j  (idx 5 ↔ 8)
    
    This creates a 3D subspace {|0010⟩, |0101⟩, |1000⟩} = {2, 5, 8}
    analogous to the qutrit case {|02⟩, |11⟩, |20⟩}
    
    Args:
        J: TTA coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    
    Returns:
        U: 16×16 unitary matrix
    """
    H = np.zeros((16, 16), dtype=complex)
    
    # Define the three states involved in TTA
    # |T1⟩_i|T1⟩_j = |01⟩_i|01⟩_j = |0101⟩
    idx_T1_T1 = 0b0101  # = 5
    
    # |S0⟩_i|S1⟩_j = |00⟩_i|10⟩_j = |0010⟩
    idx_S0_S1 = 0b0010  # = 2
    
    # |S1⟩_i|S0⟩_j = |10⟩_i|00⟩_j = |1000⟩
    idx_S1_S0 = 0b1000  # = 8
    
    # Build H_TTA in the 3D subspace
    # Following the theory document section 2.4.3, we have:
    # H_TTA = J(|S0 S1⟩⟨T1 T1| + |T1 T1⟩⟨S0 S1| + |S1 S0⟩⟨T1 T1| + |T1 T1⟩⟨S1 S0|)
    
    # Coupling: |S0 S1⟩ ↔ |T1 T1⟩
    H[idx_S0_S1, idx_T1_T1] = J
    H[idx_T1_T1, idx_S0_S1] = J
    
    # Coupling: |S1 S0⟩ ↔ |T1 T1⟩
    H[idx_S1_S0, idx_T1_T1] = J
    H[idx_T1_T1, idx_S1_S0] = J
    
    # Compute exact time evolution using scipy.linalg.expm
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    
    return U
```

#### 5.2 優先度2: 検証テストの追加

修正後、以下を検証する:

1. **構造の検証**: H_TTAが3D部分空間で作用することを確認
2. **ユニタリ性の検証**: `||U†U - I|| < 1e-10`
3. **固有値の検証**: 期待される固有値 `{-√2·J, 0, +√2·J}` と一致
4. **古典シミュレータとの比較**: per-pair分解後の結果が一致

#### 5.3 優先度3: ノートブックの検証

修正後のノートブックで以下を確認:

1. **個体数保存**: `N_S0 + N_T1 + N_S1 = 4` が常に成立
2. **3手法の一致**: 古典・Qubit・Quditの結果が許容誤差内で一致
3. **誤差の評価**: 各手法の誤差が `< 1e-3` 程度（Trotter誤差）

## 結論

**根本原因**: Qubit実装のH_TTAが不完全であり、完全なTTAプロセスを実装していない。

**解決策**: `exact_qubit_hamiltonians.py`の`build_H_TTA_qubit_unitary`を修正し、
`|S1⟩_i|S0⟩_j`への遷移を追加する。

**検証**: 修正後、ノートブックを実行し、3手法の結果が一致することを確認する。
