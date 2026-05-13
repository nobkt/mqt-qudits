# Iteration 26 検証結果の詳細分析（Iteration 27）

## 分析日時
2026-03-16

## 分析対象
- `developing/verification_results/iteration26_hamiltonian_fix_20260316T035105Z.json`
- `developing/検証結果分析_quantum_dynamics_comparison_iteration26.md`
- `tutorials/qudit_gksl_circuit_boson_simulator.py`
- `tutorials/qudit_gksl_boson_simulator.py`
- `tutorials/gksl_math_utils.py`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`（セル実行結果）
- `tutorials/quantum_dynamics_complete_comparison.ipynb`（セル実行結果）
- PR#251 履歴

## 結論要約

**Iteration 26 の修正（Issue HH）は正確かつ完全であり、95/95 の検証チェックが PASS している。ノートブック実行結果にもエラー・NaN・Inf は存在しない。**

**追加分析により、1件の軽微な API 不整合（Issue II）を発見：`QuditGKSLCircuitBosonSimulator.simulate()` の `g_eph=0` exact reduction パスで、返却辞書のキーが `g_eph≠0` パスと不一致であった。修正適用済み。**

---

## Issue HH の検証確認（Iteration 26 修正の再確認）

### 修正内容
`qudit_gksl_circuit_boson_simulator.py` の `_precompute_unitaries()` が `expm(-1j * build_H_total_boson(params) * dt/2)` を使用（一次 Trotter 積から変更）。

### 独自検証結果

| 検証項目 | 結果 |
|----------|------|
| 検証スクリプト実行（95 チェック） | 95/95 PASS ✓ |
| U_H_half 差 (DM vs Circuit) | 0.00e+00 ✓ |
| rho_full 差 (2ステップ) | 0.00e+00 ✓ |
| rho_full 差 (20ステップ) | 0.00e+00 ✓ |
| rho_phonon 差 (2ステップ) | 0.00e+00 ✓ |
| rho_el 差 (100ステップ) | 0.00e+00 ✓ |
| トレース保存 (100ステップ) | max\|Tr-1\| = 1.02e-14 ✓ |
| n_max=2 クロス一致 | 0.00e+00 ✓ |
| Trotter 誤差スケーリング | O(dt²) 確認済み（比率 9.98e+05 ≈ 10⁶） |

### ハミルトニアン構成の検証

| 検証項目 | 結果 |
|----------|------|
| H_total 分解チェック | 0.00e+00 ✓ |
| H_total エルミート性 | 0.00e+00 ✓ |
| H_phonon エルミート性 | 0.00e+00 ✓ |
| H_eph エルミート性 | 0.00e+00 ✓ |
| \|\|[H_transfer, H_eph]\|\|_F | 8.00e-03（分析文書と一致） |
| \|\|[H_phonon, H_eph]\|\|_F | 1.47e-02（分析文書と一致） |

**Issue HH の修正は数学的に正確であり、物理的に妥当である。**

---

## ノートブック実行結果の検証

### quantum_dynamics_gksl_comparison.ipynb

全 44 セルの実行結果を独自に検証：

| セル | 内容 | エラー | NaN/Inf | 備考 |
|------|------|--------|---------|------|
| 4 | 古典 DM | なし | なし | Trace: 5.33e-15 |
| 6 | Qudit DM | なし | なし | n_total=30, GPS=66 |
| 8 | Qubit DM | なし | なし | n_total=34, GPS=388 |
| 10 | 古典ボソン | なし | なし | Trace: 5.77e-15 |
| 12 | Qudit ボソン | なし | なし | n_total=16, GPS=38 |
| 14 | Qubit ボソン | なし | なし | n_total=18, GPS=216 |
| 20-28 | 比較・収束分析 | なし | なし | F > 0.999999 |
| 32 | Qudit Shot | なし | なし | Shot noise のみ |
| 36 | Qubit Shot | なし | なし | forbidden_count=0 |
| 38 | Qubit Noisy | なし | なし | 70.4% リーク（物理的に正常） |

**全セルにエラー・NaN・Inf なし。数値結果は全て物理的に妥当。**

### quantum_dynamics_complete_comparison.ipynb

全 39 セルの実行結果を独自に検証：

| セル | 内容 | エラー | NaN/Inf |
|------|------|--------|---------|
| 3 | パラメータ | なし | なし |
| 5 | 古典 | なし | なし |
| 8 | Qubit Shot | なし | なし |
| 11 | Qubit Unitary | なし | なし |
| 17 | Qudit Shot | なし | なし |
| 23 | Qudit Noisy | なし | なし |
| 30-35 | 比較・表 | なし | なし |

**全セルにエラー・NaN・Inf なし。**

---

## Issue II: g_eph=0 exact reduction パスの API 不整合

### 症状

`QuditGKSLCircuitBosonSimulator.simulate()` で `g_eph=0` のとき、`QuditGKSLCircuitSimulator` に委譲する exact reduction パスが存在する。この際、返却辞書のキーが `g_eph≠0` のときと不一致であった。

### 具体的な差異

| キー | g_eph≠0 | g_eph=0（修正前） |
|------|---------|------------------|
| `n_total_qudits` | 16 | **欠落** |
| `n_phonon_qudits` | 2 | **欠落** |
| `dim_total` | 36 | **欠落** |
| `estimated_gates_per_step` | 38 | **欠落**（`gates_per_step=30` あり） |
| `total_estimated_gates` | 76 | **欠落**（`total_gates=60` あり） |

### 根本原因

`QuditGKSLCircuitSimulator.simulate()` の返却辞書に `n_total_qudits`, `n_phonon_qudits`, `dim_total` キーが含まれず、ゲートカウントのキー名が異なる（`gates_per_step` vs `estimated_gates_per_step`）。

### 修正内容

`qudit_gksl_circuit_boson_simulator.py` の `simulate()` 内 `g_eph=0` パスで、委譲結果に不足キーを追加：

```python
result.setdefault("n_phonon_qudits", 0)
result.setdefault("n_total_qudits", n_sys + n_anc)
result.setdefault("dim_total", self.dim_el)
if "gates_per_step" in result and "estimated_gates_per_step" not in result:
    result["estimated_gates_per_step"] = result["gates_per_step"]
if "total_gates" in result and "total_estimated_gates" not in result:
    result["total_estimated_gates"] = result["total_gates"]
```

### 修正の正当性

1. `n_phonon_qudits=0` は正確：g_eph=0 では完全にフォノンが分離するため、フォノンレジスタ不要
2. `dim_total=dim_el` は正確：フォノン空間が完全に分離
3. キー名正規化はダウンストリーム互換性のため
4. 元の `gates_per_step` / `total_gates` キーも保持（既存コードとの互換性）

### 修正結果（修正後）

| キー | g_eph≠0 | g_eph=0 |
|------|---------|---------|
| `n_total_qudits` | 16 | 14 ✓ |
| `n_phonon_qudits` | 2 | 0 ✓ |
| `dim_total` | 36 | 9 ✓ |
| `estimated_gates_per_step` | 38 | 30 ✓ |
| `total_estimated_gates` | 76 | 60 ✓ |

---

## _apply_hamiltonian_step デッドコードの確認

`qudit_gksl_circuit_boson_simulator.py` の `_apply_hamiltonian_step()` メソッド（行 414-419）は `_compute_hamiltonian_unitary()` を使用するが、`simulate()` および `_trotter_step()` からは呼び出されていない（デッドコード）。

- `_trotter_step()` はプリコンピュートされた `_U_H_half` を使用
- `_apply_hamiltonian_step()` はどのファイルからも呼び出されていない
- docstring は "Apply Hamiltonian evolution using circuit-decomposed unitaries" で正確

**対応不要**：削除ではなくそのまま残置。将来のデバッグ・テスト用途に使用可能。

---

## 検証結果（Iteration 27）

### 検証チェック一覧

- **Section 1-16**: 95 件の回帰チェック（全 PASS）
- **Section 17**: 9 件の新規チェック（Issue II: g_eph=0 API 整合性）
  - `geph0_has_n_total_qudits`: n_total_qudits キーが存在 ✓
  - `geph0_n_phonon_zero`: n_phonon_qudits=0 ✓
  - `geph0_has_dim_total`: dim_total キーが存在 ✓
  - `geph0_has_estimated_gps`: estimated_gates_per_step キーが存在 ✓
  - `geph0_has_total_estimated_gates`: total_estimated_gates キーが存在 ✓
  - `geph0_n_total_sum_correct`: n_total = n_sys + n_ph + n_anc ✓
  - `geph0_rho_final_shape`: rho_final は (9,9) ✓
  - `geph0_method_indicates_reduction`: method 文字列に exact reduction を含む ✓
  - `geph0_trace_conservation`: Tr(rho_final) = 1.0 ✓

**合計: 104/104 PASS ✓**

---

## 修正方針（Iteration 27）

### 変更したファイル

1. `tutorials/qudit_gksl_circuit_boson_simulator.py`
   - `simulate()` の g_eph=0 パスに API 整合性のためのキー追加

### 新規作成ファイル

1. `tutorials/run_iteration27_verification.py`（104 チェック: 95 回帰 + 9 新規）
2. `developing/検証結果分析_quantum_dynamics_comparison_iteration27.md`（本文書）

### 変更しないファイル

- `tutorials/qudit_gksl_boson_simulator.py` — g_eph=0 パスなし、変更不要
- `tutorials/qudit_gksl_circuit_simulator.py` — 根本原因だが本 PR スコープ外
- `tutorials/stinespring_utils.py` — 正常動作
- `tutorials/gksl_math_utils.py` — 正常動作
- ノートブック — g_eph=0 は使用されないため変更不要

## 参照ファイル

- `developing/検証結果分析_quantum_dynamics_comparison_iteration26.md` — iteration 26 の分析
- `developing/verification_results/iteration26_hamiltonian_fix_20260316T035105Z.json` — iteration 26 の結果
- `tutorials/run_iteration26_verification.py` — iteration 26 の検証スクリプト
- PR#251 — Issue HH 修正の PR（closed, not merged）
