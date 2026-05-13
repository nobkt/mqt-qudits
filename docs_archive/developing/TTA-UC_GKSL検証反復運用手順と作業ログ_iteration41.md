# TTA-UC GKSL 検証反復運用手順と作業ログ - Iteration 41

## 実施日
2026-03-03

## 前回（Iteration 40）の結果
- 全26チェック PASS ✓
- Cell 25/39 のイテレーション回数を「38回」→「39回」に更新
- `build_trotter_step_classical` デッドコード削除

## Iteration 41 の修正内容

### 修正1（重大）: ボソンシミュレーターの Trotter 分割を回文半ステップに統一

**問題**: `qudit_gksl_boson_simulator.py` と `qubit_gksl_boson_simulator.py` が非回文・フルステップの Lindblad チャネル適用を使用しており、非ボソンシミュレーターの回文半ステップと不整合。

**修正**:
- `_precompute_unitaries`: `_U_stines`（dt）→ `_U_stines_half`（dt/2）を計算
- `_trotter_step`: 回文半ステップ構造に変更
  ```
  exp(L_H dt/2) · Π_{fwd} E_α(dt/2) · Π_{rev} E_α(dt/2) · exp(L_H dt/2)
  ```
- docstring を正確な記述に更新

**影響**:
- 同じ dt でより高精度な結果（Lie-Trotter 交換子誤差が2次消去されるため）
- 全シミュレーターで統一された Trotter 分割構造

### 修正2（中）: 非ボソンシミュレーターの `_U_stines` デッドコード削除

**問題**: 4つの非ボソンシミュレーターが `_precompute_unitaries` でフルステップ Stinespring ユニタリ `_U_stines` を計算しているが、`_trotter_step` / `_trotter_step_trajectory` では使用されていない。

**修正**: 4ファイルから `_U_stines` 計算を削除
- `qudit_gksl_simulator.py`
- `qubit_gksl_simulator.py`
- `qudit_gksl_shot_simulator.py`
- `qubit_gksl_shot_simulator.py`

### 修正3（低）: ノートブックのイテレーション回数更新

**修正**: Cell 25 と Cell 39 の「39回」→「40回」に更新

## 検証項目（Iteration 41）

Iteration 40 の全26チェック + 新規チェック:

| # | チェック項目 |
|---|---|
| 1-26 | Iteration 40 全チェック（回帰テスト） |
| 27 | ボソン Qudit シミュレーター収束テスト |
| 28 | ボソン Qubit シミュレーター収束テスト |
| 29 | ボソン Qudit-Qubit 一致テスト |
| 30 | ボソン Qudit zero-noise = ideal テスト |
| 31 | Cell 25/39 「40回」表記確認 |
| 32 | 非ボソンシミュレーターに `_U_stines` 属性が無いことの確認 |

## ファイル変更一覧

| ファイル | 変更内容 |
|----------|----------|
| `tutorials/qudit_gksl_boson_simulator.py` | 回文半ステップ Trotter に変更 |
| `tutorials/qubit_gksl_boson_simulator.py` | 回文半ステップ Trotter に変更 |
| `tutorials/qudit_gksl_simulator.py` | `_U_stines` デッドコード削除 |
| `tutorials/qubit_gksl_simulator.py` | `_U_stines` デッドコード削除 |
| `tutorials/qudit_gksl_shot_simulator.py` | `_U_stines` デッドコード削除 |
| `tutorials/qubit_gksl_shot_simulator.py` | `_U_stines` デッドコード削除 |
| `tutorials/quantum_dynamics_gksl_comparison.ipynb` | Cell 25/39 「40回」更新 |
| `developing/検証結果分析_iteration41_詳細分析.md` | 詳細分析ドキュメント |
| `tutorials/run_tta_uc_gksl_verification_iteration41.py` | 検証スクリプト |

## 次のステップ

1. ユーザーが `run_tta_uc_gksl_verification_iteration41.py` をローカル実行
2. 結果を push
3. 結果に基づいて次のイテレーションへ
