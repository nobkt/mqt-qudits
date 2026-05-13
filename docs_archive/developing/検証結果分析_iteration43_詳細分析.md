# Iteration 43 詳細分析: Iteration 42 検証結果の網羅的レビュー

## 分析日時
2026-03-04

## 分析対象
- `tutorials/run_tta_uc_gksl_verification_iteration42.py` 実行結果
- `developing/verification_results/iteration42_verification_20260303T110735Z.json`
- `developing/verification_results/iteration42_verification_20260303T110735Z.md`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb` の全40セル出力
- GKSL シミュレーターコード全体（16ファイル）
- 詳細理論・設計・仕様
- Stinespring dilation の数学的正確性の独立検証
- GKSL 超演算子の構成の独立検証
- Lindblad 演算子の物理的正確性の独立検証
- T/dt 誤差係数の漸近挙動の独立分析

## 1. Iteration 42 検証結果の確認

全33チェックが PASS ✓

### チェック結果一覧

| # | チェック項目 | 結果 | 備考 |
|---|---|---|---|
| 1 | トレース距離単調減少 | ✓ | |
| 2 | 収束次数 ≈ 1.0 (平均: 1.2010) | ✓ | O(dt) 収束を確認 |
| 3 | 忠実度単調増加 | ✓ | |
| 4 | 密度行列品質 | ✓ | |
| 5 | n_steps=100 で T < 1e-3 | ✓ | T=7.319e-04 |
| 6 | Cell 27 テーブル値の整合性 | ✓ | |
| 7 | Cell 25 「41回」表記 | ✓ | |
| 8 | Cell 39 「41回」表記 | ✓ | |
| 9 | Cell 27 n_steps=200 値 (~3.6e-04) | ✓ | |
| 10 | 粒子数保存 | ✓ | Cl: 2.22e-14, Qd: 8.17e-14 |
| 11 | 中間時刻人口動態 | ✓ | max diff: 1.598e-03 |
| 12 | 2分子系収束次数 | ✓ | Rate: 1.0022 |
| 13 | ユニタリ限界精度 | ✓ | T: 1.34e-14 |
| 14 | 固有値スペクトル一致 | ✓ | max diff: 6.42e-04 |
| 15 | Qudit Shot-DM一致 | ✓ | T: 2.81e-03, F: 0.9999 |
| 16 | Qubit Shot-DM一致 | ✓ | T: 2.81e-03, F: 0.9999 |
| 17 | Qudit Noisy(p=0)=DM | ✓ | T: 0.00e+00 |
| 18 | Qubit Noisy(p=0)=DM | ✓ | T: 0.00e+00 |
| 19 | Qudit DM アンシラ数 (4分子) | ✓ | 26=26 |
| 20 | Qudit DM アンシラ数 (2分子) | ✓ | 12=12 |
| 21 | Qudit DM ゲート数 (4分子) | ✓ | 59=59 |
| 22 | Qudit DM ゲート数 (2分子) | ✓ | 27=27 |
| 23 | Qudit Shot ゲート数 (4分子) | ✓ | 59=59 |
| 24 | Qudit Shot ゲート数 (2分子) | ✓ | 27=27 |
| 25 | Qubit Shot ゲート数推定あり | ✓ | 158=158 |
| 26 | Qudit-Qubit Shot 同一結果 | ✓ | T: 3.34e-19 |
| 27 | Boson Qudit 収束次数 | ✓ | Rate: 1.0009 |
| 28 | Boson Qubit 収束次数 | ✓ | Rate: 1.0009 |
| 29 | Boson Qudit-Qubit一致 | ✓ | T: 0.00e+00 |
| 30 | Boson 密度行列品質 | ✓ | tr_err: 4.44e-16 |
| 31 | Cell 25/39 「41回」確認 | ✓ | 40回なし確認 |
| 32 | 非Bosonに_U_stinesなし | ✓ | 4シミュレータ全て確認 |
| 33 | 誤差係数 T/dt の安定性 | ✓ | 偏差: 9.02e-03 |

## 2. 独立数学検証

### 2.1 Stinespring dilation の数学的正確性（再確認）

**検証**: `stinespring_unitary_from_lindblad(L, dt)` のコードを独立にトレースした。

生成子構成:
```python
G[:d_sys, d_sys:] = L.conj().T  # = L†
G[d_sys:, :d_sys] = L
```

$G = \begin{pmatrix} 0 & L^\dagger \\ L & 0 \end{pmatrix}$ はエルミート（$(L^\dagger)^\dagger = L$）であるため、$U = \exp(-i\sqrt{\Delta t} \cdot G)$ は常にユニタリ。✓

**テイラー展開**:

$G^2 = \begin{pmatrix} L^\dagger L & 0 \\ 0 & LL^\dagger \end{pmatrix}$

$U \approx I - i\sqrt{\Delta t} G - \frac{\Delta t}{2} G^2 + ...$

$K_0 = U_{00} = I - \frac{\Delta t}{2} L^\dagger L + O(\Delta t^2)$, $K_1 = U_{10} = -i\sqrt{\Delta t} \cdot L + O(\Delta t^{3/2})$

$\mathcal{E}(\rho) = K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger = \rho + \Delta t (L\rho L^\dagger - \frac{1}{2}\{L^\dagger L, \rho\}) + O(\Delta t^2)$ ✓

### 2.2 `apply_stinespring_to_density_matrix` の正確性

**テンソル積順序の確認**:
```python
rho_ext = np.kron(env0, rho)  # (env ⊗ sys) 順序
rho_prime = U @ rho_ext @ U.conj().T
rho_out = rho_prime[:d_sys, :d_sys] + rho_prime[d_sys:, d_sys:]
```

`kron(env0, rho)` により、インデックス `[0, d_sys)` が env=|0⟩、`[d_sys, 2*d_sys)` が env=|1⟩ に対応。
部分トレースは:
$\text{Tr}_{\text{env}}(\rho') = \langle 0|\rho'|0\rangle + \langle 1|\rho'|1\rangle = \rho'_{00} + \rho'_{11}$
これは `rho_prime[:d_sys, :d_sys] + rho_prime[d_sys:, d_sys:]` と一致。✓

### 2.3 GKSL 超演算子の正確性（再確認）

列優先ベクトル化 $\text{vec}(AXB) = (B^T \otimes A) \text{vec}(X)$ を使用:

- ハミルトニアン部: $\mathcal{L}_H = -\frac{i}{\hbar}(I \otimes H - H^T \otimes I)$ ✓
- 散逸部: $\mathcal{L}_D = \bar{L} \otimes L - \frac{1}{2}(I \otimes L^\dagger L) - \frac{1}{2}((L^\dagger L)^T \otimes I)$ ✓
  - `L_op.conj()` は $\bar{L}$（要素ごとの共役 = $(L^\dagger)^T$）であることを確認 ✓

### 2.4 Lindblad 演算子の物理的正確性（再確認）

各 Lindblad 演算子に $\sqrt{\gamma}$ が含まれていることを `build_lindblad_operators` のコードで確認:
```python
ops.append((sq * reduce(np.kron, op_list), gamma))  # sq = sqrt(gamma)
```

TTA 2チャネル × 3ペア = 6、蛍光 × 4、リン光 × 4、内部変換 × 4、ISC S→T × 4、ISC T→S × 4 = 合計26。✓

### 2.5 回文順序積 Trotter ステップの正確性

全シミュレーター（QuditGKSLSimulator, QubitGKSLSimulator, QuditGKSLBosonSimulator, QubitGKSLBosonSimulator, QuditGKSLShotSimulator, QubitGKSLShotSimulator）の `_trotter_step` / `_trotter_step_trajectory` を確認:

```
exp(L dt/2) · prod_{α=1..n} E_α(dt/2) · prod_{α=n..1} E_α(dt/2) · exp(L dt/2)
```

全て同一の回文構造を使用。✓

### 2.6 T/dt 誤差係数の漸近挙動分析

Check 33 の T/dt 値を詳細に分析した:

| n_steps | dt | T/dt |
|--------:|------:|--------:|
| 10 | 10.0 | 0.001294 |
| 20 | 5.0 | 0.000816 |
| 50 | 2.0 | 0.000746 |
| 100 | 1.0 | 0.000732 |
| 200 | 0.5 | 0.000725 |

$T = C_1 \cdot dt + C_2 \cdot dt^2 + O(dt^3)$ であるため、$T/dt = C_1 + C_2 \cdot dt + O(dt^2)$。

最後の2点から線形外挿:
- 傾き $C_2 \approx (0.000732 - 0.000725) / (1.0 - 0.5) = 0.000014$
- 切片 $C_1 \approx 0.000725 - 0.000014 \times 0.5 = 0.000718$

検証:
- dt=2.0: 予測 0.000718 + 0.000028 = 0.000746、実測 0.000746 ✓（完全一致）
- dt=5.0: 予測 0.000788、実測 0.000816（3.5% 偏差 → 2次以上の寄与）
- dt=10.0: 予測 0.000858、実測 0.001294（50.8% 偏差 → 漸近領域外）

**結論**: dt ≤ 2 (n_steps ≥ 50) で漸近的 O(dt) 収束が高精度で成立。dt > 2 では高次補正項が支配的。これは数学的に期待される挙動であり、バグではない。✓

### 2.7 ノートブック出力の整合性

Cell 26 の収束解析出力が検証スクリプトの結果と完全一致:

| n_steps | ノートブック T | 検証スクリプト T | 一致 |
|--------:|--------:|--------:|----:|
| 10 | 1.294305e-02 | 1.294305e-02 | ✓ |
| 20 | 4.081287e-03 | 4.081287e-03 | ✓ |
| 50 | 1.492513e-03 | 1.492513e-03 | ✓ |
| 100 | 7.319238e-04 | 7.319238e-04 | ✓ |
| 200 | 3.626913e-04 | 3.626913e-04 | ✓ |

Cell 22 の粒子数保存もCheck 10と整合。✓

Cell 24 の忠実度:
- Classical vs Qudit: F = 0.999999 ✓
- Classical vs Qubit: F = 0.999999 ✓
- Qubit vs Qudit: F = 1.000000 ✓

Cell 38 のショット比較:
- Qudit DM vs Qudit Shot (no noise): F = 0.998056 ✓
- Qubit DM vs Qubit Shot (no noise): F = 0.998056 ✓（Qudit と同一 → Check 26 と整合）

### 2.8 Qubit ノイズシミュレータの禁止状態リーケージ

Cell 36 出力:
- Forbidden state count (final measurement): 702/1000
- Max trace deficit (leakage): 0.7217 (72.2%)

p_depol=0.01 で100ステップ（各ステップ ~60ゲート = 合計 ~6000ゲート適用）において、ほぼ全トラジェクトリで複数回のパウリエラーが発生する。パウリエラーは物理状態を禁止状態 |11⟩ にマッピングできるため、高いリーケージは期待される挙動。

Qudit ノイズシミュレータ (Cell 32): Tr(ρ) = 1.000000（リーケージなし）- Weyl-Heisenberg ノイズは qutrit 空間内で完結するため。✓

## 3. 新たに発見された問題点

### 問題1（低）: ノートブックのイテレーション回数表記

**影響ファイル**: `tutorials/quantum_dynamics_gksl_comparison.ipynb`

**詳細**:
- Cell 25: 「41回のイテレーション検証」→「42回」に更新すべき
- Cell 39: 「41回の検証イテレーション」→「42回」に更新すべき

Iteration 42 が全33チェック PASS で完了したため。

## 4. 重大なバグ・数学的誤り

**なし。**

Stinespring dilation、GKSL 超演算子、Lindblad 演算子構造、Trotter 分割、回文順序積、ボソン拡張空間、テンソル積順序、部分トレース、射影測定（ショットシミュレータ）の全てにおいて、数学的正確性が独立検証により確認された。

## 5. コード全体の数学的正確性（総括）

### 5.1 Stinespring dilation の実装 ✓
- 生成子 $G = [[0, L^\dagger], [L, 0]]$ はエルミート → $U$ は常にユニタリ
- チャネルは常に CPTP（$\sum_k K_k^\dagger K_k = I$）
- 1次近似誤差は $O(\Delta t^2)$/チャネル適用
- `apply_stinespring_to_density_matrix` のテンソル積順序が正しい

### 5.2 GKSL 超演算子 ✓
- 列優先ベクトル化と整合する超演算子構造
- トレース保存（$6.94 \times 10^{-18}$）
- `L_op.conj()` は $\bar{L}$（要素ごとの共役）

### 5.3 Lindblad 演算子構造 ✓
- TTA 2チャネル × 3ペア + 5種 × 4分子 = 26 演算子
- 各演算子に $\sqrt{\gamma}$ が含まれている
- 減衰レートが物理的期待値と完全一致

### 5.4 Qubit-Qutrit マッピング ✓
- Check 26 で T ≈ 3.34e-19 により Qudit-Qubit 一致を確認

### 5.5 収束特性 ✓
- Rate(T) → 1.0 で $O(\Delta t)$ 収束を確認
- $T/\Delta t \to C_1$ で線形漸近を確認

### 5.6 ボソンシミュレーター ✓
- 回文半ステップ Trotter 構造に統一済み
- Check 27-30 で収束・品質を確認
- Qudit-Qubit 一致（T=0.0）

### 5.7 ショットシミュレーター ✓
- Born 則に基づく確率的アンシラ測定
- 回文順序の _trotter_step_trajectory が密度行列版と同一構造
- 無限ショット極限で密度行列結果に一致

### 5.8 ノイズシミュレーター ✓
- p=0 でノイズなし結果に一致（Check 17, 18）
- Qudit: Weyl-Heisenberg 脱分極（リーケージなし）
- Qubit: Pauli 脱分極（禁止状態リーケージあり → 期待される挙動）

### 5.9 ヒューリスティック・フォールバック ✓
- 密度行列の強制正規化、固有値クリッピング等は一切不使用
- `quantum_fidelity` の `np.maximum(evals, 0.0)` は浮動小数点精度対応であり、ヒューリスティックではない

## 6. 修正計画

| # | 修正内容 | 優先度 | ファイル |
|---|----------|--------|----------|
| 1 | ノートブックのイテレーション回数更新 (41→42) | 低 | `quantum_dynamics_gksl_comparison.ipynb` |
| 2 | Iteration 43 検証スクリプト作成 | 中 | `run_tta_uc_gksl_verification_iteration43.py` |

## 7. 結論

Iteration 42 の全33チェックは PASS であり、コードベース全体の数学的正確性が独立した数値検証およびコード精読によって確認された。

唯一の発見事項は、ノートブックのイテレーション回数表記の更新（41回→42回）であり、これは機能的な問題ではなく文書の保守作業である。

T/dt 誤差係数の漸近解析により、dt ≤ 2 (4分子系で n_steps ≥ 50) で高精度な O(dt) 収束が確認された。大きな dt では高次補正項が支配的になるが、これは数学的に期待される挙動である。
