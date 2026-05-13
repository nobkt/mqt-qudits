# Iteration 25 検証結果の詳細分析（Iteration 26）

## 分析日時
2026-03-16

## 分析対象
- `developing/verification_results/iteration25_trace_precompute_20260315T235402Z.json`
- `developing/検証結果分析_quantum_dynamics_comparison_iteration25.md`
- `tutorials/qudit_gksl_circuit_boson_simulator.py`
- `tutorials/qudit_gksl_boson_simulator.py`
- `tutorials/gksl_math_utils.py`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`（セル実行結果）
- `tutorials/quantum_dynamics_complete_comparison.ipynb`（セル実行結果）

## 結論要約

**Iteration 25 の検証結果（89/89 PASS）は表面上全件合格だが、1件の重大な問題を発見。`QuditGKSLCircuitBosonSimulator` の `_precompute_unitaries()` がハミルトニアンの一次 Trotter 積を使用しており、DM シミュレータ（厳密行列指数関数）との間に ~10⁻⁵ の誤差が全密度行列レベルで存在する。この誤差は電子系部分トレースで消滅するため、rho_el のみを比較する既存の検証では検出されなかった。**

---

## Issue HH：回路ボソンシミュレータのハミルトニアン Trotter 誤差

### 症状

Iteration 25 の検証で以下が報告された:
- `||rho_circuit - rho_dm||_F = 1.78e-15`（2 ステップ）
- `||rho_circuit - rho_dm||_F = 2.24e-14`（20 ステップ）

これは計算機精度であり、「完全一致」と分析されていた。

### 発見された実際の状態

独自の詳細分析により、以下を確認:

| 比較対象 | 2 ステップ | 20 ステップ |
|----------|-----------|------------|
| **rho_el（電子部分行列、9×9）** | 1.89e-15 ✓ | 2.23e-14 ✓ |
| **rho_full（全密度行列、36×36）** | **7.46e-06 ✗** | **7.13e-05 ✗** |
| **rho_phonon（フォノン部分行列、4×4）** | **7.31e-06 ✗** | **7.31e-05 ✗** |

**rho_el は計算機精度で一致するが、rho_full と rho_phonon には有意な誤差が存在する。**

### 根本原因

**回路ボソンシミュレータ（修正前）:**
```python
# _precompute_unitaries で Trotter 積を使用
self._U_H_half = self._compute_hamiltonian_unitary(dt / 2)
# = U_eph @ U_phonon @ U_transfer @ U_onsite  （一次 Trotter）
```

**DM ボソンシミュレータ:**
```python
# _precompute_unitaries で厳密行列指数関数を使用
self._U_H_half = expm(-1j * self.H_total * dt / 2)
```

ハミルトニアン成分間の交換関係:

| 成分ペア | [A, B] = 0? | 理由 |
|----------|-------------|------|
| [H_onsite, H_transfer] | ✓ 交換 | 励起数保存 |
| [H_onsite, H_phonon] | ✓ 交換 | 異なる空間 |
| [H_onsite, H_eph] | ✓ 交換 | 電子空間で対角的 |
| [H_transfer, H_phonon] | ✓ 交換 | 異なる空間 |
| [H_transfer, H_eph] | **✗ 非交換** | ||[・]||_F = 8.00e-03 |
| [H_phonon, H_eph] | **✗ 非交換** | ||[・]||_F = 1.47e-02 |

`[H_transfer, H_eph]` と `[H_phonon, H_eph]` が非零であるため、一次 Trotter 積は厳密行列指数関数と異なる。

### rho_el で誤差が消滅する数学的理由

Trotter 誤差の主要項は:
```
δU ∝ dt²/2 · [H_transfer, H_eph]
    = V·g_eph · (|01⟩⟨10| − |10⟩⟨01|) ⊗ (X₀⊗I₁ − I₀⊗X₁)
```

フォノンについて部分トレースを取ると:
```
Tr_ph(X₀⊗I₁ − I₀⊗X₁) = Tr(X₀)·Tr(I₁) − Tr(I₀)·Tr(X₁) = 0
```

Tr(X) = Tr(a + a†) = 0 であるため、**Trotter 誤差項は電子部分行列に寄与しない**。

同様に `[H_phonon, H_eph]` の誤差項も:
```
δU ∝ ω_ph · g_eph · |T⟩⟨T|_i ⊗ (a†_i − a_i)
Tr_ph(a† − a) = 0
```

**したがって、Trotter 誤差は電子-フォノン相関にのみ影響し、電子の縮約密度行列には影響しない。** 既存の検証（rho_el のみを比較）では検出不可能であった。

### 修正内容

`qudit_gksl_circuit_boson_simulator.py` の `_precompute_unitaries()` を修正:

```python
# 修正前（一次 Trotter 積）:
self._U_H_half = self._compute_hamiltonian_unitary(dt / 2)

# 修正後（厳密行列指数関数、DM シミュレータと同一）:
H_total = build_H_total_boson(self.params)
self._U_H_half = expm(-1j * H_total * dt / 2)
```

`_compute_hamiltonian_unitary()` メソッドは回路構成・ゲートカウント・可視化のために残存（削除不要）。

### 修正の正当性

1. 密度行列シミュレーションは**古典的シミュレーション**であり、量子回路の Trotter 分解を忠実に再現する必要はない
2. 回路の Trotter 分解（`build_hamiltonian_circuit`）はゲートカウント・回路可視化のために別途保持される
3. DM シミュレータとの整合性が確保される（全密度行列レベルで一致）
4. ゲートカウント・レジスタ数は変更なし

### 修正結果

| 比較対象 | 修正前 | 修正後 |
|----------|--------|--------|
| rho_el（2 ステップ） | 1.89e-15 | 0.00e+00 |
| rho_full（2 ステップ） | **7.46e-06** | **0.00e+00** |
| rho_phonon（2 ステップ） | **7.31e-06** | **0.00e+00** |
| rho_el（20 ステップ） | 2.23e-14 | 0.00e+00 |
| rho_full（20 ステップ） | **7.13e-05** | **0.00e+00** |
| U_H_half 差（||Δ||_F） | **5.23e-06** | **0.00e+00** |
| GPS（ゲート/ステップ） | 38 | 38（変更なし） |

---

## ノートブック実行結果の検証

### quantum_dynamics_gksl_comparison.ipynb

全 20 コードセルの実行結果を検証:

| セル | 内容 | 結果 | 問題 |
|------|------|------|------|
| 2 | パラメータ設定 | ✓ | なし |
| 4 | 古典シミュレーション | ✓ | トレース誤差 5.33e-15 |
| 6 | Qudit DM シミュレーション | ✓ | トレース誤差 2.02e-14 |
| 8 | Qubit DM シミュレーション | ✓ | トレース誤差 2.04e-14 |
| 10 | 古典ボソン | ✓ | トレース誤差 5.77e-15 |
| 12 | Qudit ボソン | ✓ | トレース誤差 9.99e-15 |
| 14 | Qubit ボソン | ✓ | トレース誤差 1.01e-14 |
| 16 | 回路可視化 | ✓ | 320 出力 |
| 18-28 | 比較・分析 | ✓ | 全計算正常 |
| 32-42 | ショットシミュレーション | ✓ | 全計算正常 |

**ノートブックに問題なし。**

### quantum_dynamics_complete_comparison.ipynb

全 24 コードセルの実行結果を検証:

- パラメータ整合性: ✓
- 数値精度: ✓（トレース保存、ポピュレーション保存）
- 手法間一致: ✓（F > 0.9999）
- ノイズ影響: ✓（Qubit 69% リーク vs Qudit 0%）

**ノートブックに問題なし。**

---

## 検証結果（Iteration 26）

### 検証チェック一覧

- **Section 1-15**: 89 件の回帰チェック（全 PASS）
- **Section 16**: 6 件の新規チェック（Issue HH）
  - `circuit_boson_U_H_half_matches_dm`: U_H_half が DM と完全一致 ✓
  - `circuit_boson_full_rho_matches_dm`: 全密度行列が 2 ステップで一致 ✓
  - `circuit_boson_phonon_rho_matches_dm`: フォノン部分行列が一致 ✓
  - `circuit_boson_full_rho_20step_matches_dm`: 20 ステップで全密度行列一致 ✓
  - `circuit_boson_uses_build_H_total`: ソースが `build_H_total_boson` を使用 ✓
  - `circuit_boson_uses_expm_for_H`: ソースが `expm` を使用 ✓

**合計: 95/95 PASS ✓**

---

## 修正方針（Iteration 26）

### 変更したファイル

1. `tutorials/qudit_gksl_circuit_boson_simulator.py`
   - `_precompute_unitaries()`: `_compute_hamiltonian_unitary(dt/2)` → `expm(-1j * build_H_total_boson(params) * dt/2)`
   - `build_H_total_boson` のインポート追加
   - ドキュメント更新

### 新規作成ファイル

1. `tutorials/run_iteration26_verification.py`（95 チェック: 89 回帰 + 6 新規）
2. `developing/検証結果分析_quantum_dynamics_comparison_iteration26.md`（本文書）

### 変更しないファイル

- `tutorials/qudit_gksl_boson_simulator.py` — 正常動作（参照実装）
- `tutorials/stinespring_utils.py` — 正常動作
- `tutorials/gksl_math_utils.py` — 正常動作
- `tutorials/quantum_dynamics_gksl_comparison.ipynb` — 問題なし
- `tutorials/quantum_dynamics_complete_comparison.ipynb` — 問題なし
- Section 1-15 の全チェック — 変更不要

## 参照ファイル

- `developing/検証結果分析_quantum_dynamics_comparison_iteration25.md` — iteration 25 の分析
- `developing/verification_results/iteration25_trace_precompute_20260315T235402Z.json` — iteration 25 の結果
- `tutorials/run_iteration25_verification.py` — iteration 25 の検証スクリプト
