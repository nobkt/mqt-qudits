# Iteration 38 検証結果詳細分析

- 分析日時: 2026-03-03
- 対象: Iteration 37 検証結果（全18チェック）、ノートブック全40セル、シミュレーターコード全体

## 1. Iteration 37 検証結果の確認

Iteration 37 で実施された全18チェックは全て PASS しています。

| # | チェック項目 | 結果 | 判定 |
|---|---|---|---|
| 1 | トレース距離単調減少 | ✓ | 正常 |
| 2 | 収束次数 ≈ 1.0 (平均: 1.2010) | ✓ | 理論と一致 |
| 3 | 忠実度単調増加 | ✓ | 正常 |
| 4 | 密度行列品質 | ✓ | 機械精度 |
| 5 | n_steps=100 で T < 1e-3 | ✓ | 正常 |
| 6 | Cell 27 テーブル値の整合性 | ✓ | 正常 |
| 7 | Cell 25 「34回」表記 | ✓ | 正常 |
| 8 | Cell 39 「34回」表記 | ✓ | 正常 |
| 9 | Cell 27 n_steps=200 値 (~3.6e-04) | ✓ | 正常 |
| 10 | 粒子数保存 | ✓ | 機械精度 |
| 11 | 中間時刻人口動態 | ✓ | 正常 |
| 12 | 2分子系収束次数 | ✓ | 正常 |
| 13 | ユニタリ限界精度 | ✓ | 機械精度 |
| 14 | 固有値スペクトル一致 | ✓ | 正常 |
| 15 | Qudit Shot-DM一致 | ✓ | 統計的範囲内 |
| 16 | Qubit Shot-DM一致 | ✓ | 統計的範囲内 |
| 17 | Qudit Noisy(p=0)=DM | ✓ | 完全一致 |
| 18 | Qubit Noisy(p=0)=DM | ✓ | 完全一致 |

## 2. 数学的整合性の検証

### 2.1 Stinespring dilation の実装

Stinespring ユニタリの構成:

```
G = [[0, L†], [L, 0]]
θ = √dt
U = exp(-iθG)
```

展開すると:
- K₀ = U₀₀ ≈ I - (dt/2)L†L + O(dt²)
- K₁ = U₁₀ ≈ -i√dt L + O(dt^{3/2})

チャネル E(ρ) = K₀ρK₀† + K₁ρK₁†:
```
E(ρ) ≈ ρ + dt(LρL† - ½{L†L, ρ}) + O(dt²)
```

これは GKSL 散逸項の 1 次近似であり、**正しく実装されている**。

### 2.2 Trotter 分割の構造

```
Φ_step(ρ) = U_H(dt/2) · [∏_{α=1..n} E_α(dt/2)] · [∏_{α=n..1} E_α(dt/2)] · U_H(dt/2)
```

- **Strang H-D 分割**: O(dt³)/ステップ → 全体 O(dt²)
- **回文順序 Lindblad 積**: Lie-Trotter 交換子誤差の 2 次消去 → O(dt³)/ステップ
- **Stinespring 近似**: O(dt²/4)/チャネル × 52 チャネル = O(dt²)/ステップ → **全体 O(dt) が支配的**

### 2.3 収束次数の理論的整合性

実測データ:
| n_steps | dt | Rate(T) |
|---------|------|---------|
| 10→20 | 10→5 | 1.665 |
| 20→50 | 5→2 | 1.098 |
| 50→100 | 2→1 | 1.028 |
| 100→200 | 1→0.5 | 1.013 |

Rate(T) → 1.0 への収束が確認された。大きな dt での Rate > 1.0 は高次項（O(dt²), O(dt³)）の寄与によるもので、これは理論的に予想される正常な挙動。

### 2.4 リウビリアン超演算子の構成

`build_gksl_superoperator` は列優先ベクトル化（Fortran 順序）を使用:
```
L_H = (-i/ℏ)(I⊗H - H^T⊗I)
L_D = Σ_α [L_α*⊗L_α - ½(I⊗L†_αL_α + (L†_αL_α)^T⊗I)]
```

`vectorize_density_matrix`（列優先 `flatten(order='F')`）と整合しており、**正しい**。

### 2.5 Shot-DM 等価性

量子軌道法の数学的等価性: CPTP マップの線形性により、無限ショット極限で密度行列結果に収束する。

```
⟨|ψ_final⟩⟨ψ_final|⟩_shots = E_n(E_{n-1}(...E_1(|ψ⟩⟨ψ|)...))
```

**チェック 15 と 16 が完全に同一の値**（T: 2.8069e-03, F: 0.999949）を示す理由:
- Qubit シミュレーターは Qudit 演算子を大きな空間に埋め込むが、物理的部分空間の進化は同一
- 同一の RNG シード(seed=42)により、同一の測定確率から同一の測定結果が生成される
- 抽出された qutrit 空間の密度行列は完全に一致する

これは**正常な挙動であり、バグではない**。ただし、ドキュメントに説明がないため、混乱を招く可能性がある。

## 3. 発見された問題点

### 3.1 【バグ】QuditGKSLSimulator のアンシラ数がハードコード（重要度: 低）

**ファイル**: `tutorials/qudit_gksl_simulator.py`

```python
# 現在のコード（行67-68）
self.n_ancilla_qubits = 26
```

`n_ancilla_qubits` が 26 にハードコードされている。4 分子系（デフォルト）では正しいが、2 分子系では実際のアンシラ数は 12 であるべき。

**正しい値**: `len(self.lindblad_ops)` を使用すべき。

**影響**: シミュレーション結果自体には影響なし。メタデータ（レポートされるアンシラ数）が誤る。

### 3.2 【バグ】QuditGKSLSimulator のゲート数推定がハードコード（重要度: 低）

**ファイル**: `tutorials/qudit_gksl_simulator.py`

```python
# 現在のコード（行200）
gates_per_step = 4 + 3 + 26 * 2  # 59 gates
```

ゲート数が N=4 分子系の値 59 にハードコードされている。

**正しい計算**: 
```python
n_lindblad = len(self.lindblad_ops)
gates_per_step = self.n_system_qudits + len(self.params.neighbors) + n_lindblad * 2
```

| N_molecules | 正しいゲート数 | ハードコード値 | 差異 |
|-------------|--------------|-------------|------|
| 4 | 4+3+52=59 | 59 | なし |
| 2 | 2+1+24=27 | 59 | **+32 (誤)** |
| 3 | 3+2+36=41 | 59 | **+18 (誤)** |

### 3.3 【バグ】QuditGKSLShotSimulator のゲート数推定がハードコード（重要度: 低）

**ファイル**: `tutorials/qudit_gksl_shot_simulator.py`

```python
# 現在のコード（行287）
gates_per_step = 4 + 3 + 26 * 2
```

上記 3.2 と同一の問題。`n_ancilla_qubits` は正しく `len(self.lindblad_ops)` を使用しているにもかかわらず、ゲート数は依然としてハードコード。

### 3.4 【バグ】QubitGKSLShotSimulator にゲート数推定がない（重要度: 低）

**ファイル**: `tutorials/qubit_gksl_shot_simulator.py`

`QubitGKSLShotSimulator.simulate()` の返却値に `estimated_gates_per_step` と `total_estimated_gates` が含まれていない。

- `QuditGKSLShotSimulator` にはゲート数推定がある（ハードコードだが）
- `QubitGKSLSimulator`（DM版）にはゲート数推定がある（動的計算）
- `QubitGKSLShotSimulator` にだけゲート数推定がない → **不整合**

## 4. 設計上の観察事項（バグではない）

### 4.1 冗長な Lindblad 演算子

以下の遷移演算子が重複している:
- **蛍光 (Γ_fl)** と **内部変換 (k_IC)**: どちらも `|0⟩⟨2|` 遷移
- **燐光 (Γ_ph)** と **ISC T→S (k_ISC_TS)**: どちらも `|0⟩⟨1|` 遷移

これらは数学的に `√(γ₁+γ₂) |a⟩⟨b|` に統合可能で、Lindblad チャネル数を 26 から 18 に削減できる。

**影響**:
- 量子回路のゲート数が約 30% 削減可能（52 → 36 Stinespring ゲート/ステップ）
- シミュレーション結果の O(dt) 収束次数は変わらないが、誤差の係数が改善される可能性がある
- 現在の実装は数学的に正しく、バグではない

### 4.2 ノートブック Cell 25, 27, 39 の文書整合性

Iteration 37 で修正された文書修正（「34回のイテレーション検証」表記、n_steps=200 の値更新）は正しく適用されている。

## 5. 修正方針

以下の修正を実施する:

1. **QuditGKSLSimulator**: `n_ancilla_qubits` を動的に計算
2. **QuditGKSLSimulator**: `gates_per_step` を動的に計算
3. **QuditGKSLShotSimulator**: `gates_per_step` を動的に計算
4. **QubitGKSLShotSimulator**: `estimated_gates_per_step` と `total_estimated_gates` を返却値に追加

これらはメタデータの修正であり、シミュレーションの物理的正確性には影響しない。

## 6. 次のイテレーションの検証スクリプト

Iteration 38 の検証スクリプトで以下を追加検証:
- 各シミュレーターのアンシラ数が Lindblad 演算子数と一致することの確認
- ゲート数推定が動的に計算されることの確認（2 分子系と 4 分子系の比較）
- QubitGKSLShotSimulator にゲート数推定が含まれることの確認
