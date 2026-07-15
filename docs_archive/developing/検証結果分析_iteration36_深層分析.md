# Iteration 36: 深層分析結果

- 分析日時: 2026-03-02
- 分析対象: Iteration 35検証結果、ノートブック全セル、コード全体、数学的理論

## 1. 分析スコープ

以下を網羅的に分析した：

1. `run_tta_uc_gksl_verification_iteration35.py` の実行結果（9チェック全PASS）
2. `quantum_dynamics_gksl_comparison.ipynb` の全40セルの実行結果
3. コアモジュール: `stinespring_utils.py`, `gksl_math_utils.py`, `classical_gksl_simulator.py`, `qudit_gksl_simulator.py`, `gksl_physical_parameters.py`
4. Stinespring dilation の数学的正当性
5. GKSL superoperator の構成
6. Lindblad演算子の物理的正確性
7. 収束次数の理論的根拠と数値的検証

## 2. 追加数値検証結果

既存の検証スクリプト（iteration 35）では検証されていない項目について、独自に数値検証を実施した。

### 2.1 粒子数保存（N_S0 + N_T1 + N_S1 = N_molecules）

| シミュレータ | 最大粒子数誤差 |
|:---|:---|
| ClassicalGKSLSimulator | 2.22e-14 |
| QuditGKSLSimulator | 1.24e-14 |

**結論**: 両シミュレータとも機械精度で粒子数を保存。問題なし。

### 2.2 トレース保存（全時刻ステップ）

| シミュレータ | 最大トレース誤差 |
|:---|:---|
| ClassicalGKSLSimulator | 5.33e-15 |
| QuditGKSLSimulator | 2.89e-15 |

**結論**: 機械精度でトレースを保存。問題なし。

### 2.3 中間時刻での人口動態比較（4分子系、n_steps=100）

| 時刻 t | |ΔN_S0| | |ΔN_T1| | |ΔN_S1| |
|-------:|:------:|:------:|:------:|
| 0.0 | 0.0e+00 | 0.0e+00 | 0.0e+00 |
| 10.0 | 1.21e-04 | 2.95e-04 | 1.74e-04 |
| 25.0 | 7.93e-04 | 1.27e-03 | 4.75e-04 |
| 50.0 | 1.21e-03 | 1.60e-03 | 3.87e-04 |
| 75.0 | 1.33e-03 | 1.47e-03 | 1.37e-04 |
| 100.0 | 1.27e-03 | 1.16e-03 | 1.04e-04 |

**結論**: 中間時刻での人口差は最大 ~1.6e-3 であり、n_steps=100 (dt=1.0) での O(dt) 収束と整合。全時刻で一貫した精度。問題なし。

### 2.4 Stinespring チャネルの CPTP 性

全26チャネルについて dt=1.0 で検証（最大混合状態 ρ=I/81 に対して）：

| 項目 | 結果 |
|:---|:---|
| トレース保存 | Tr(E(ρ)) = 1.00000000000000 (全チャネル) |
| 最小固有値 | > 1.2e-02 (全チャネル) |
| エルミート誤差 | 0.0e+00 (全チャネル) |

**結論**: 全Stinespringチャネルが CPTP 写像として正しく機能。問題なし。

### 2.5 GKSL Liouvillian superoperator のトレース保存性

$$\text{vec}(I)^\top \cdot \mathcal{L} = 0$$

数値検証結果: $\|\text{vec}(I)^\top \mathcal{L}\|_\infty = 3.47 \times 10^{-17}$

**結論**: superoperator が数学的にトレース保存生成子であることを数値確認。問題なし。

### 2.6 単一チャネル Stinespring 近似誤差（2分子系、dim=9）

最初のTTAチャネルについて、厳密チャネル exp(L_D dt) と Stinespring 近似 E(dt) の差を計測：

| dt | T(exact, stine) | T/dt² |
|-----:|:---:|:---:|
| 10.00 | 8.650e-03 | 8.650e-05 |
| 5.00 | 2.375e-03 | 9.498e-05 |
| 2.00 | 4.016e-04 | 1.004e-04 |
| 1.00 | 1.023e-04 | 1.023e-04 |
| 0.50 | 2.580e-05 | 1.032e-04 |
| 0.10 | 1.040e-06 | 1.040e-04 |

**結論**: T/dt² が定数 (~1.04e-4) に収束しており、Stinespring 近似がステップあたり O(dt²) であることを数値的に確認。これは理論予測と完全に一致。

### 2.7 全 Trotter ステップの誤差（2分子系）

| dt | T(exact, trotter) | T/dt² |
|-----:|:---:|:---:|
| 10.00 | 6.608e-03 | 6.608e-05 |
| 5.00 | 2.073e-03 | 8.292e-05 |
| 2.00 | 3.803e-04 | 9.507e-05 |
| 1.00 | 9.952e-05 | 9.952e-05 |
| 0.50 | 2.545e-05 | 1.018e-04 |
| 0.10 | 1.037e-06 | 1.037e-04 |

**結論**: 全 Trotter ステップの誤差も O(dt²)/step で、Stinespring 近似が支配的誤差源であることを確認。Strang 分割および回文順序積の誤差はこれより小さい。

### 2.8 ユニタリ（散逸なし）のみの収束（2分子系）

全散逸率を 0 に設定した場合のハミルトニアン発展：

| n_steps | dt | T(exact, trotter) |
|--------:|-----:|:---:|
| 10 | 10.0 | 1.60e-14 |
| 50 | 2.0 | 1.48e-14 |
| 100 | 1.0 | 2.22e-14 |
| 200 | 0.5 | 4.15e-14 |
| 500 | 0.2 | 0.00e+00 |

**結論**: ハミルトニアン発展は機械精度で厳密。これは `expm(-iH dt/2)` が厳密な行列指数関数であるため。Trotter 分割誤差はハミルトニアン部分には存在しない（H_0 と H_transfer の分割はしていない）。**近似誤差は100% Stinespring dilation に由来する。**

### 2.9 2分子系の収束次数

| n_steps 区間 | Rate(T) |
|:---|:---|
| 10→20 | 1.0653 |
| 20→50 | 1.0217 |
| 50→100 | 1.0095 |
| 100→200 | 1.0048 |
| 200→500 | 1.0022 |

**結論**: 2分子系でも Rate(T) は漸近的に 1.0 に収束するが、厳密に 1.0000 には至らない。iteration 34 の作業ログにある「Rate(T) = 1.0000 for 2-molecule systems」は微小な高次項を無視した近似的表現であり、厳密には「Rate(T) → 1.0（dt → 0 で）」が正確。ただし、これは実質的に O(dt) 収束を確認するものであり、本質的な問題ではない。

### 2.10 誤差スケーリングの精密フィット（4分子系）

$T(\Delta t) = c_1 \Delta t + c_2 \Delta t^2$ のフィッティング（dt ≤ 2 のデータ使用）：

- $c_1 = 7.19 \times 10^{-4}$
- $c_2 = 1.3 \times 10^{-5}$
- $c_2 / c_1 = 0.018$

**結論**: 2次補正項は1次項の約1.8%であり、dt ≤ 2 では支配的でない。大きな dt（dt=10）では3次以上の項も寄与するため、Rate(T) が 1.0 から外れる。これは標準的な漸近展開挙動であり、問題ない。

## 3. コード分析結果

### 3.1 `stinespring_utils.py`

#### `stinespring_unitary_from_lindblad(L, dt)`

- **生成子**: $G = \begin{pmatrix} 0 & L^\dagger \\ L & 0 \end{pmatrix}$, $\theta = \sqrt{\Delta t}$, $U = \exp(-i\theta G)$ ✓
- **ユニタリ性検証**: $\|U^\dagger U - I\|_F < 10^{-10}$ ✓
- **数学的正当性**: 摂動展開で $\mathcal{E}(\rho) = \rho + \Delta t(L\rho L^\dagger - \frac{1}{2}\{L^\dagger L, \rho\}) + O(\Delta t^2)$ を確認 ✓

#### `apply_stinespring_to_density_matrix(rho, U)`

- **Kronecker 順序**: `np.kron(env0, rho)` → 環境 ⊗ 系 の順序 ✓
- **G 行列のブロック構造と整合**: `G[0:d, d:] = L†`, `G[d:, 0:d] = L` → `|0⟩⟨1| ⊗ L† + |1⟩⟨0| ⊗ L` ✓
- **部分トレース**: `rho_out = rho_prime[:d, :d] + rho_prime[d:, d:]` ✓

#### `build_gksl_superoperator(H, lindblad_ops)`

- **ハミルトニアン部分**: $\mathcal{L}_H = -\frac{i}{\hbar}(I \otimes H - H^T \otimes I)$ ✓
- **散逸部分**: $\mathcal{L}_D = \sum_\alpha [\overline{L_\alpha} \otimes L_\alpha - \frac{1}{2}(I \otimes L_\alpha^\dagger L_\alpha + (L_\alpha^\dagger L_\alpha)^T \otimes I)]$ ✓
- **列優先ベクトル化**: `vec(AXB) = (B^T ⊗ A) vec(X)` と整合 ✓

### 3.2 `gksl_math_utils.py`

#### `build_lindblad_operators(params)`

- **TTA チャネル**: 3対 × 2チャネル = 6 ✓
- **蛍光**: 4分子 × 1 = 4 ✓
- **燐光**: 4分子 × 1 = 4 ✓
- **内部転換**: 4分子 × 1 = 4 ✓
- **ISC S→T**: 4分子 × 1 = 4 ✓
- **ISC T→S**: 4分子 × 1 = 4 ✓
- **合計**: 26 ✓
- **√γ の包含**: 各 L_alpha に sqrt(gamma) が含まれている ✓

#### `build_onsite_hamiltonian(params)` / `build_transfer_hamiltonian(params)`

- **H_0**: $\sum_i (E_T|1\rangle\langle 1| + E_S|2\rangle\langle 2|)_i$ ✓
- **H_transfer**: $\sum_{\langle i,j\rangle} V(|0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + \text{h.c.})$ ✓
- **エルミート性検証**: コード内で `np.allclose(H, H.conj().T)` をチェック ✓

### 3.3 `qudit_gksl_simulator.py`

#### `_precompute_unitaries(dt)`

- **ハミルトニアン**: `expm(-1j * H_total * dt/2)` — H_0 + H_transfer の厳密行列指数関数 ✓
- **Stinespring**: 半ステップ `dt/2` と全ステップ `dt` の両方を事前計算 ✓

#### `_trotter_step(rho)`

- **Strang 分割**: `exp(L_H dt/2) · [prod E_α(dt/2) · prod_rev E_α(dt/2)] · exp(L_H dt/2)` ✓
- **回文順序**: forward + reversed で Lie-Trotter 交換子誤差を2次消去 ✓
- **各チャネルの半ステップ**: dt/2 で2回適用 → 実効的に dt の散逸 ✓

### 3.4 `classical_gksl_simulator.py`

- **厳密解**: `expm_multiply(L_super, vec_0, ...)` で $\exp(\mathcal{L}t)$ を計算 ✓
- **CPTP 保証**: Lindblad-GKS 定理により構造的に保証 ✓

### 3.5 `gksl_physical_parameters.py`

- **パラメータ検証**: エネルギー正値性、散逸率非負性、時間スケール階層、弱結合条件 ✓
- **デフォルト値**: 物理的に妥当 ✓

## 4. ノートブック全セル分析

### 全40セルの実行状態

全セルがエラーなしで実行完了。出力値はコードおよび理論と整合。

### 主要セルの出力確認

| セル | 内容 | 出力 | 判定 |
|:---:|:---|:---|:---:|
| 4 | 古典GKSL (NB) | N_S0=3.2762, Tr誤差=5.33e-15 | ✓ |
| 6 | Qudit GKSL (NB) | Tr誤差=2.02e-14 | ✓ |
| 8 | Qubit GKSL (NB) | Tr誤差=2.04e-14 | ✓ |
| 22 | トレース・粒子数保存 | 全シナリオ < 1e-13 | ✓ |
| 24 | 忠実度比較 | F(Cl,Qd)=0.999999, F(Cl,Qb)=0.999999 | ✓ |
| 26 | 収束解析 | Rate(T) → 1.0 確認 | ✓ |

### Cell 25, 27, 39 の文書修正確認

| セル | 修正内容 | 確認 |
|:---:|:---|:---:|
| 25 | 「34回のイテレーション検証」 | ✓ |
| 27 | n_steps=200 → 「~3.6e-04」 | ✓ |
| 39 | 「34回の検証イテレーション」 | ✓ |

## 5. 発見事項

### 5.1 コードバグ

**発見されなかった。** 数学的実装は正確であり、理論と一致する。

### 5.2 軽微な文書記載

iteration 34 の作業ログにおける「Rate(T) = 1.0000 for 2-molecule systems」は、漸近極限での表現としては正しいが、厳密には有限 dt での Rate は 1.0 よりわずかに大きい（例: dt=0.5 で Rate ≈ 1.005）。ただし、これは O(dt) 収束の結論に影響しない。

### 5.3 検証スクリプトの改善余地

既存のスクリプト（iteration 35）には以下のチェック項目が含まれていない：

1. **粒子数保存の検証**
2. **中間時刻での人口動態比較**
3. **2分子系での収束検証**
4. **ユニタリ限界（散逸なし）での精度検証**
5. **固有値スペクトルの比較**

これらを iteration 36 の検証スクリプトに追加する。

## 6. 結論

**コードの数学的実装は正確であり、問題は発見されなかった。**

- Stinespring dilation: 正しく GKSL 散逸子を O(dt²)/step の精度で近似
- Strang 分割: 正しく H-D 分割を O(dt³)/step で実装
- 回文順序積: Lie-Trotter 交換子誤差を 2次消去
- 全体収束次数: O(dt)（1次）— Stinespring 近似が支配的
- 保存則: トレース保存・粒子数保存とも機械精度

追加検証項目を含む iteration 36 検証スクリプトを作成し、検証フレームワークの網羅性を向上させる。
