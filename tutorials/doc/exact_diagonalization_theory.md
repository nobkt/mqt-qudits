# 厳密対角化法による量子ダイナミクス解析理論

## Exact Diagonalization Theory for Quantum Dynamics Analysis

**作成日 / Date**: 2025-10-17
**バージョン / Version**: 1.0.0
**目的 / Purpose**: Qudit量子アルゴリズムの検証のための厳密解計算理論

---

## 1. 概要 / Overview

### 1.1 厳密対角化法とは

厳密対角化法（Exact Diagonalization, ED）は、量子多体系のハミルトニアンを完全に対角化することで、系の厳密な固有状態と固有エネルギーを求める手法である。

本文書では、4分子直線配置系の量子ダイナミクスにおいて：

- **Qudit量子アルゴリズム** による時間発展（鈴木トロッター分解）
- **厳密対角化法** による時間発展（解析解）

の2つを比較し、量子アルゴリズムの精度を検証する方法を詳述する。

### 1.2 なぜ厳密対角化が必要か

Qudit量子アルゴリズムの実装検証において、厳密対角化は以下の役割を果たす：

1. **精度検証**: トロッター分解の誤差評価
2. **ベンチマーク**: アルゴリズムの正しさの確認
3. **収束性テスト**: 時間刻み $\Delta t$ の依存性評価
4. **物理的妥当性**: 結果の物理的解釈の確認

### 1.3 適用範囲

厳密対角化は状態空間の次元 $d$ に対して $O(d^3)$ の計算量を要するため、小規模系に限定される：

| 系のサイズ | 状態空間次元    | 厳密対角化の実行可能性 |
| ---------- | --------------- | ---------------------- |
| N=2 分子   | $3^2 = 9$       | ✅ 容易                |
| N=3 分子   | $3^3 = 27$      | ✅ 容易                |
| N=4 分子   | $3^4 = 81$      | ✅ 可能（本実装）      |
| N=5 分子   | $3^5 = 243$     | ⚠️ やや困難            |
| N=6 分子   | $3^6 = 729$     | ❌ 困難                |
| N≥7 分子   | $3^N \geq 2187$ | ❌ 実用的でない        |

本実装では **N=4分子系**（81次元）を対象とする。

---

## 2. 数学的定式化 / Mathematical Formulation

### 2.1 ハミルトニアン行列の構築

#### 2.1.1 全ハミルトニアン

4分子系の全ハミルトニアンは以下の3項から成る：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}}
$$

計算基底 $\{|n_0, n_1, n_2, n_3\rangle\}_{n_i \in \{0,1,2\}}$ において、$81 \times 81$ 行列 $\mathbf{H}_{\text{total}}$ を構築する。

#### 2.1.2 対角項 $\hat{H}_0$

$$
\hat{H}_0 = \sum_{i=0}^{3} \left( E_T |1\rangle_i\langle 1| + E_S |2\rangle_i\langle 2| \right)
$$

行列要素：

$$
\langle n_0, n_1, n_2, n_3 | \hat{H}_0 | n_0', n_1', n_2', n_3' \rangle = \delta_{nn'} \sum_{i=0}^{3} E(n_i)
$$

ここで：

$$
E(n) = \begin{cases}
0 & (n = 0, \ \text{基底状態}) \\
E_T & (n = 1, \ \text{三重項}) \\
E_S & (n = 2, \ \text{励起一重項})
\end{cases}
$$

#### 2.1.3 エネルギー移動項 $\hat{H}_{\text{transfer}}$

隣接ペア $\langle i, j \rangle \in \{(0,1), (1,2), (2,3)\}$ について：

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + \text{h.c.} \right)
$$

行列要素（$i < j$）：

$$
\begin{align}
&\langle \ldots, 0, \ldots, 1, \ldots | \hat{H}_{\text{transfer}} | \ldots, 1, \ldots, 0, \ldots \rangle = V_{ij} \\
&\langle \ldots, 1, \ldots, 0, \ldots | \hat{H}_{\text{transfer}} | \ldots, 0, \ldots, 1, \ldots \rangle = V_{ij}
\end{align}
$$

その他の要素は0。

#### 2.1.4 三重項-三重項消滅項 $\hat{H}_{\text{TTA}}$

隣接ペア $\langle i, j \rangle$ について：

$$
\hat{H}_{\text{TTA}} = \sum_{\langle i,j \rangle} J_{ij} \left( |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1| + |0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1| + \text{h.c.} \right)
$$

行列要素：

$$
\begin{align}
&\langle \ldots, 2, \ldots, 0, \ldots | \hat{H}_{\text{TTA}} | \ldots, 1, \ldots, 1, \ldots \rangle = J_{ij} \\
&\langle \ldots, 0, \ldots, 2, \ldots | \hat{H}_{\text{TTA}} | \ldots, 1, \ldots, 1, \ldots \rangle = J_{ij} \\
&(\text{およびエルミート共役})
\end{align}
$$

### 2.2 固有値分解

構築した $81 \times 81$ エルミート行列 $\mathbf{H}_{\text{total}}$ を対角化：

$$
\mathbf{H}_{\text{total}} = \mathbf{V} \mathbf{\Lambda} \mathbf{V}^\dagger
$$

ここで：

- $\mathbf{\Lambda} = \text{diag}(\lambda_0, \lambda_1, \ldots, \lambda_{80})$：固有値（エネルギー準位）
- $\mathbf{V} = [|\phi_0\rangle, |\phi_1\rangle, \ldots, |\phi_{80}\rangle]$：固有ベクトル（エネルギー固有状態）

NumPyでの実装：

```python
eigenvalues, eigenvectors = np.linalg.eigh(H_total)
```

`np.linalg.eigh` はエルミート行列に対して **厳密な** 固有値分解を実行（ヒューリスティックではない）。

### 2.3 厳密な時間発展

#### 2.3.1 時間発展演算子

時間発展演算子：

$$
\hat{U}(t) = e^{-i \hat{H}_{\text{total}} t / \hbar}
$$

固有基底での表現：

$$
\hat{U}(t) = \sum_{k=0}^{80} e^{-i \lambda_k t / \hbar} |\phi_k\rangle \langle \phi_k|
$$

行列形式：

$$
\mathbf{U}(t) = \mathbf{V} \cdot \text{diag}(e^{-i \lambda_0 t / \hbar}, \ldots, e^{-i \lambda_{80} t / \hbar}) \cdot \mathbf{V}^\dagger
$$

#### 2.3.2 初期状態の固有基底展開

初期状態 $|\Psi_0\rangle$ を固有基底に展開：

$$
|\Psi_0\rangle = \sum_{k=0}^{80} c_k |\phi_k\rangle
$$

係数：

$$
c_k = \langle \phi_k | \Psi_0 \rangle
$$

行列形式：

$$
\mathbf{c} = \mathbf{V}^\dagger \cdot |\Psi_0\rangle
$$

#### 2.3.3 時間発展後の状態

時間 $t$ での状態：

$$
|\Psi(t)\rangle = \hat{U}(t) |\Psi_0\rangle = \sum_{k=0}^{80} c_k e^{-i \lambda_k t / \hbar} |\phi_k\rangle
$$

行列形式：

$$
|\Psi(t)\rangle = \mathbf{V} \cdot \text{diag}(e^{-i \lambda_0 t / \hbar}, \ldots) \cdot \mathbf{c}
$$

これは **厳密解** であり、近似を含まない。

### 2.4 個体数の計算

状態 $|\Psi(t)\rangle$ から各準位の個体数を計算：

$$
\begin{align}
N_{S_0}(t) &= \sum_{i=0}^{3} \sum_{n: n_i=0} |\langle n | \Psi(t) \rangle|^2 \\
N_{T_1}(t) &= \sum_{i=0}^{3} \sum_{n: n_i=1} |\langle n | \Psi(t) \rangle|^2 \\
N_{S_1}(t) &= \sum_{i=0}^{3} \sum_{n: n_i=2} |\langle n | \Psi(t) \rangle|^2
\end{align}
$$

---

## 3. 実装アルゴリズム / Implementation Algorithm

### 3.1 厳密対角化シミュレータのクラス設計

```python
class ExactDiagonalizationSolver:
    """
    厳密対角化による量子ダイナミクスの解析解計算

    ヒューリスティックな手法を使用せず、数学的に厳密な計算のみを実行：
    - np.linalg.eigh: エルミート行列の固有値分解（厳密）
    - 行列-ベクトル積: 厳密な線形代数演算
    - 指数関数: 数学的に定義された演算
    """

    def __init__(self, params: PhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.dim = 3**self.N
        self.H_total = None
        self.eigenvalues = None
        self.eigenvectors = None

    def build_total_hamiltonian(self) -> np.ndarray:
        """全ハミルトニアン行列を構築（厳密）"""
        pass

    def diagonalize(self):
        """ハミルトニアンを対角化（厳密）"""
        pass

    def time_evolution(self, t: float, initial_state: np.ndarray) -> np.ndarray:
        """時間発展を計算（厳密解）"""
        pass

    def calculate_populations(self, state: np.ndarray) -> Dict[str, float]:
        """個体数を計算"""
        pass

    def simulate(self, T_total: float, N_points: int, initial_state_type: str) -> Dict:
        """完全なシミュレーション（厳密解）"""
        pass
```

### 3.2 ハミルトニアン構築アルゴリズム

```python
def build_total_hamiltonian(self) -> np.ndarray:
    """
    全ハミルトニアン行列の構築

    81×81 エルミート行列を構築
    """
    dim = self.dim
    H = np.zeros((dim, dim), dtype=complex)

    # 1. H₀項（対角）
    for idx in range(dim):
        config = index_to_config(idx, self.N, 3)
        energy = 0.0
        for level in config:
            if level == 1:
                energy += self.params.E_T
            elif level == 2:
                energy += self.params.E_S
        H[idx, idx] += energy

    # 2. H_transfer項（非対角）
    for pair_idx, (i, j) in enumerate(self.params.neighbors):
        V = self.params.V[pair_idx]

        for idx1 in range(dim):
            config1 = index_to_config(idx1, self.N, 3)

            # |...0...1...⟩ ↔ |...1...0...⟩
            if config1[i] == 0 and config1[j] == 1:
                config2 = config1.copy()
                config2[i] = 1
                config2[j] = 0
                idx2 = config_to_index(config2, 3)
                H[idx1, idx2] += V
                H[idx2, idx1] += V

    # 3. H_TTA項（非対角）
    for pair_idx, (i, j) in enumerate(self.params.neighbors):
        J = self.params.J[pair_idx]

        for idx1 in range(dim):
            config1 = index_to_config(idx1, self.N, 3)

            # |...1...1...⟩ → |...2...0...⟩
            if config1[i] == 1 and config1[j] == 1:
                config2 = config1.copy()
                config2[i] = 2
                config2[j] = 0
                idx2 = config_to_index(config2, 3)
                H[idx1, idx2] += J
                H[idx2, idx1] += J

                # |...1...1...⟩ → |...0...2...⟩
                config3 = config1.copy()
                config3[i] = 0
                config3[j] = 2
                idx3 = config_to_index(config3, 3)
                H[idx1, idx3] += J
                H[idx3, idx1] += J

    return H
```

### 3.3 時間発展アルゴリズム

```python
def time_evolution(self, t: float, initial_state: np.ndarray) -> np.ndarray:
    """
    厳密な時間発展

    Args:
        t: 時間
        initial_state: 初期状態ベクトル（81次元）

    Returns:
        時間発展後の状態ベクトル（厳密解）
    """
    if self.eigenvalues is None or self.eigenvectors is None:
        raise RuntimeError("先にdiagonalize()を実行してください")

    # 初期状態を固有基底に展開
    coeffs = self.eigenvectors.conj().T @ initial_state.flatten()

    # 時間発展
    time_evolved_coeffs = coeffs * np.exp(-1j * self.eigenvalues * t / self.params.hbar)

    # 元の基底に戻す
    state_final = self.eigenvectors @ time_evolved_coeffs

    return state_final
```

---

## 4. Qudit量子アルゴリズムとの比較手法

### 4.1 フィデリティによる精度評価

2つの状態ベクトル間のフィデリティ：

$$
F = \left| \langle \Psi_{\text{exact}}(t) | \Psi_{\text{Trotter}}(t) \rangle \right|^2
$$

- $F = 1$：完全一致（トロッター誤差なし）
- $F < 1$：トロッター分解による誤差あり

### 4.2 個体数の比較

各時刻での個体数の差：

$$
\Delta N_{\alpha}(t) = N_{\alpha}^{\text{Trotter}}(t) - N_{\alpha}^{\text{exact}}(t), \quad \alpha \in \{S_0, T_1, S_1\}
$$

### 4.3 収束性テスト

時間刻み $\Delta t$ を変化させて、フィデリティの変化を評価：

$$
F(\Delta t) \approx 1 - C \cdot (\Delta t)^p
$$

2次トロッター分解では $p \approx 3$ が期待される。

---

## 5. 実装上の注意事項

### 5.1 厳密計算とヒューリスティック計算の区別

#### ✅ 厳密計算（使用可能）

以下は数学的に定義された厳密な演算：

1. **`np.linalg.eigh`**: エルミート行列の固有値分解

   - ハウスホルダー変換とQRアルゴリズムによる厳密解
   - 数値誤差は機械精度のみ（$\sim 10^{-16}$）

2. **行列積**: $\mathbf{A} \cdot \mathbf{B}$

   - 線形代数の基本演算

3. **複素指数関数**: $e^{-i \lambda t / \hbar}$

   - 数学的に定義された関数

4. **内積**: $\langle \psi | \phi \rangle$
   - 線形代数の基本演算

#### ❌ ヒューリスティック計算（使用禁止）

以下は近似計算であり、本実装では使用しない：

1. **`scipy.linalg.expm`**: 行列指数関数

   - パデ近似やスケーリング＆二乗法による近似
   - 本実装では固有値分解による厳密計算を使用

2. **時間刻み $\Delta t$ による近似**
   - 厳密対角化では時間刻みは不要（連続時間での厳密解）

### 5.2 計算量とメモリ使用量

#### 計算量

- ハミルトニアン構築: $O(d^2)$ （$d = 81$）
- 固有値分解: $O(d^3)$
- 時間発展1ステップ: $O(d^2)$

#### メモリ使用量

- ハミルトニアン行列: $81 \times 81 \times 16$ bytes = 104 KB (complex128)
- 固有ベクトル: $81 \times 81 \times 16$ bytes = 104 KB
- 状態ベクトル: $81 \times 16$ bytes = 1.3 KB

合計: 約 210 KB （十分小さい）

### 5.3 数値安定性

固有値分解の数値安定性を保つため：

1. **エルミート性の保証**:

   ```python
   H = (H + H.conj().T) / 2  # 数値誤差によるエルミート性の破れを修正
   ```

2. **規格化の確認**:
   ```python
   norm = np.linalg.norm(state)
   assert abs(norm - 1.0) < 1e-10
   ```

---

## 6. 検証と妥当性

### 6.1 保存則の確認

#### エネルギー保存

時間発展中、エネルギー期待値は保存される：

$$
\langle \Psi(t) | \hat{H}_{\text{total}} | \Psi(t) \rangle = \text{const.}
$$

実装での確認：

```python
E_initial = np.vdot(state_initial, H_total @ state_initial).real
E_final = np.vdot(state_final, H_total @ state_final).real
assert abs(E_final - E_initial) < 1e-10
```

#### 規格化保存

状態ベクトルのノルムは常に1：

$$
\langle \Psi(t) | \Psi(t) \rangle = 1
$$

### 6.2 極限ケースのテスト

#### テスト1: 初期状態が固有状態の場合

初期状態を固有状態 $|\phi_k\rangle$ にした場合：

$$
|\Psi(t)\rangle = e^{-i \lambda_k t / \hbar} |\phi_k\rangle
$$

位相因子を除いて状態は変化しない。

#### テスト2: 相互作用がない場合（$V=J=0$）

ハミルトニアンが対角：

$$
\hat{H}_{\text{total}} = \hat{H}_0 = \text{diag}(\ldots)
$$

時間発展は各基底状態に独立な位相：

$$
|\Psi(t)\rangle = \sum_n e^{-i E_n t / \hbar} \langle n | \Psi_0 \rangle |n\rangle
$$

---

## 7. 結論 / Conclusion

### 7.1 厳密対角化法の役割

本理論により、以下が可能となる：

1. **Qudit量子アルゴリズムの精度検証**

   - トロッター分解の誤差評価
   - フィデリティによる定量的評価

2. **収束性の確認**

   - 時間刻み $\Delta t$ 依存性
   - トロッター次数の効果

3. **物理的妥当性の保証**
   - エネルギー保存則
   - 個体数の時間発展

### 7.2 理論的厳密性

本実装は以下を保証する：

✅ **厳密な数学的手法のみ使用**

- 固有値分解による厳密解
- ヒューリスティック不使用
- 近似を含まない時間発展

✅ **数値的安定性**

- エルミート行列の対称性保持
- 規格化の厳密な保存
- 保存則の検証

✅ **実装可能性**

- N=4分子系で実用的な計算量
- メモリ使用量: ~210 KB
- 計算時間: 秒オーダー

---

**実装完了後の検証事項 / Verification Checklist**

- [ ] ハミルトニアン行列のエルミート性確認
- [ ] 固有値が実数であることの確認
- [ ] エネルギー保存則の確認
- [ ] 状態ベクトルの規格化保存の確認
- [ ] Quditアルゴリズムとのフィデリティ計算
- [ ] 時間刻み依存性の評価
- [ ] 個体数動態の比較プロット

---

**参考文献 / References**

1. `quantum_dynamics_molecular_triplet_states.md` - 基礎理論
2. `suzuki_trotter_decomposition_theory.md` - トロッター分解理論
3. `qudit_quantum_algorithm_for_molecular_triplet_dynamics.md` - Quditアルゴリズム
4. NumPy Documentation - `numpy.linalg.eigh`

---

**END OF THEORY DOCUMENT**
