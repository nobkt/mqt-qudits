# 検証結果分析: Quantum Dynamics Comparison - Iteration 17

## 1. 背景

### 1.1 Iteration 16の結果

Iteration 16ではBug L-O（ショットベースおよびボソンシミュレータの
`estimated_gates_per_step`における回文Trotter 2×因子の欠落）が修正され、
18/18チェックが全てPASSした。

### 1.2 Iteration 16の残課題

Iteration 16の分析ドキュメント（Section 6, item 3）で、以下の問題が
明示的に指摘されていた：

> ボソンシミュレータの`estimated_gates_per_step`に電子-フォノン結合ゲートが
> 含まれていない件の検討（回路ボソンシミュレータでは`n_eph_gates`が含まれる）

## 2. 発見されたバグの詳細

### 2.1 Bug P: quditボソンシミュレータのH_ephゲート欠落

**ファイル**: `qudit_gksl_boson_simulator.py`

**問題**: ボソンモデルのハミルトニアンは3つの項から成る：

```
H_total = H_el_ext + H_phonon + H_eph
```

ここで：
- `H_el_ext` = H_onsite ⊗ I_ph + H_transfer ⊗ I_ph（電子的ハミルトニアン）
- `H_phonon` = I_el ⊗ Σ_i ω_ph n̂_i（フォノン周波数項）
- `H_eph` = g_eph × Σ_i |1⟩_i⟨1| ⊗ (a_i + a†_i)（ホルシュタイン電子-フォノン結合）

しかし`estimated_gates_per_step`の計算では、H_eph結合ゲートが
含まれていなかった。H_ephはN個の分子について、各分子の電子自由度と
フォノン自由度を結合する2体ゲート（cu_two）を必要とする。

```python
# 修正前（74）— H_ephゲート欠落
gates_per_step = (
    2 * (self.n_system_qudits + self.n_phonon_qudits + len(self.params.neighbors))
    + self.n_ancilla_qudits * 2
)

# 修正後（82）— H_ephゲート追加
gates_per_step = (
    2 * (self.n_system_qudits + self.n_phonon_qudits
         + len(self.params.neighbors) + self.params.N_molecules)
    + self.n_ancilla_qudits * 2
)
```

**内訳変更**:
| 項目 | 半ステップあたり | 2半ステップ |
|------|-----------------|------------|
| 電子オンサイト (cu_one) | N = 4 | 8 |
| フォノンオンサイト (cu_one) | N = 4 | 8 |
| 転送 (cu_two) | n_pairs = 3 | 6 |
| **電子-フォノン結合 (cu_two)** | **N = 4** | **8** ← **新規追加** |
| Stinespring (fwd+rev) | n_lindblad = 26 | 52 |
| **合計** | | **82** |

### 2.2 Bug Q: qubitボソンシミュレータのH_ephゲート欠落

**ファイル**: `qubit_gksl_boson_simulator.py`

**問題**: 同様にH_eph結合ゲートが欠落。qubitエンコーディングでは、
各H_eph結合は2 el qubits + 2 ph qubits = 4 qubit unitary として
分解され、転送ゲートと同様の複雑さ（~10基本ゲート）を持つ。

```python
# 修正前（404）— H_ephゲート欠落
gates_per_step = (
    2 * (self.n_el_qubits + self.n_ph_qubits + len(self.params.neighbors) * 10)
    + self.n_ancilla * 6 * 2
)

# 修正後（484）— H_ephゲート追加
gates_per_step = (
    2 * (self.n_el_qubits + self.n_ph_qubits
         + len(self.params.neighbors) * 10 + self.params.N_molecules * 10)
    + self.n_ancilla * 6 * 2
)
```

**内訳変更**:
| 項目 | 半ステップあたり | 2半ステップ |
|------|-----------------|------------|
| 電子Rz (1 per qubit) | n_el = 8 | 16 |
| フォノンRz (1 per qubit) | n_ph = 8 | 16 |
| 転送 (~10 per pair) | 3×10 = 30 | 60 |
| **電子-フォノン結合 (~10 per mol)** | **4×10 = 40** | **80** ← **新規追加** |
| Stinespring (~6 per channel, fwd+rev) | 26×6 = 156 | 312 |
| **合計** | | **484** |

## 3. バグの根本原因

### 3.1 ハミルトニアン項の不完全な列挙

ボソンモデルのハミルトニアン `H_total = H_el_ext + H_phonon + H_eph` は
3つの独立した物理的相互作用を含む。しかし`estimated_gates_per_step`の
コメントと計算式では、最初の2つ（電子ハミルトニアンとフォノン周波数）のみが
列挙され、3番目の電子-フォノン結合（H_eph）が漏れていた。

### 3.2 密度行列シミュレータでは隠蔽される

密度行列シミュレータでは、`_U_H_half = expm(-i * H_total * dt/2)` として
ハミルトニアン全体を1つの行列指数関数で適用する。H_ephは物理的に正しく
シミュレーションに含まれるが、ゲートカウント推定では明示的に数えられない。

### 3.3 H_ephの物理的意味

ホルシュタイン電子-フォノン結合は、電子状態（特にT1三重項状態）と
フォノンモード（格子振動）を結合する。TTA-UC過程において、この結合は
エネルギー散逸と振動緩和に重要な役割を果たす。

量子回路実装では、各分子について1つのcu_two（qudit）または
~10基本ゲート（qubit）が追加で必要となる。

## 4. 修正内容

### 4.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `qudit_gksl_boson_simulator.py` | `estimated_gates_per_step`: 74 → 82; H_ephコメント追加 |
| `qubit_gksl_boson_simulator.py` | `estimated_gates_per_step`: 404 → 484; H_ephコメント追加 |

### 4.2 変更不要

- 非ボソンシミュレータ（`*_gksl_simulator.py`、`*_shot_simulator.py`等）:
  H_ephを含まないため変更不要
- 回路シミュレータ（`*_circuit_simulator.py`）:
  ボソンモデル非対応（`with_boson=True`で`ValueError`）
- `_trotter_step()`: 全シミュレータで正しく実装済み
  （H_ephは`H_total`に含まれ、`expm(-i*H_total*dt/2)`で正しく適用される）

### 4.3 影響範囲

全てのバグは報告値（`estimated_gates_per_step`）のみの問題であり、
物理的なシミュレーション結果（人口動態、トレース保存等）には影響しない。

## 5. 検証結果

`tutorials/run_iteration17_verification.py` で16項目を検証：

| # | チェック項目 | 結果 |
|---|-------------|------|
| 1 | quditボソン ソースにH_eph項あり | PASS |
| 2 | quditボソン GPS ≠ 74（旧値） | PASS |
| 3 | quditボソン GPS = 82 | PASS |
| 4 | qubitボソン ソースにH_eph項あり | PASS |
| 5 | qubitボソン GPS ≠ 404（旧値） | PASS |
| 6 | qubitボソン GPS = 484 | PASS |
| 7 | H_eph ≠ 0 (ノルム > 0) | PASS |
| 8 | H_eph エルミート性 | PASS |
| 9 | H_eph N結合項の確認 | PASS |
| 10 | quditボソン 2×因子保存（回帰） | PASS |
| 11 | qubitボソン 2×因子保存（回帰） | PASS |
| 12 | 回帰: qudit密度行列 GPS = 66 | PASS |
| 13 | 回帰: qubit密度行列 GPS = 388 | PASS |
| 14 | 回帰: quditショット GPS = 66 | PASS |
| 15 | 回帰: qubitショット GPS = 388 | PASS |
| 16 | ボソン qubit/qudit比率 = 5.90× | PASS |

Iteration 16の全チェック(18/18)も引き続きPASSすることを確認済み。

## 6. ゲートカウント一覧（全シミュレータ）

### 6.1 非ボソンモデル

| シミュレータ | GPS | 内訳 |
|-------------|-----|------|
| qudit DM | 66 | 2×(4+3) + 26×2 |
| qubit DM | 388 | 2×(8+30) + 26×12 |
| qudit shot | 66 | 同上 |
| qubit shot | 388 | 同上 |
| qudit circuit | 66 | 2×7 + 2×20 + 2×6 |
| qubit circuit | 66 | 同上（複合ゲート単位） |

### 6.2 ボソンモデル

| シミュレータ | GPS | 内訳 |
|-------------|-----|------|
| qudit boson | **82** | 2×(4+4+3+4) + 26×2 |
| qubit boson | **484** | 2×(8+8+30+40) + 26×12 |

### 6.3 注意事項

- qubit回路シミュレータの`gates_per_step = 66`は複合ゲート（UnitaryGate）
  単位のカウント。qubit DM/shotの388は基本ゲート（CX+1qゲート）推定値。
  これらは異なる抽象化レベルでの測定であり、直接比較には注意が必要。
- ボソン回路シミュレータは現在非対応（`with_boson=True`で`ValueError`）。

## 7. 次のステップ

1. ユーザーが検証スクリプトをローカルで実行して結果を確認
2. ノートブック（`quantum_dynamics_gksl_comparison.ipynb`）でボソン
   シミュレーション結果のゲートカウント表示を確認
3. qubit回路シミュレータの`gates_per_step`（複合ゲート単位66）と
   qubit DM/shotの`estimated_gates_per_step`（基本ゲート単位388）の
   不整合について検討が必要（直接比較のためには統一が望ましい）
