# Iteration 37: TTA-UC GKSL 検証反復運用手順と作業ログ

## 作業概要

- 作業日: 2026-03-03
- 前回: Iteration 36（全14チェックPASS）
- 今回: Iteration 37（ショットベース回文順序一貫性修正 + 検証チェック追加）

## 参照ドキュメント

- `developing/検証結果分析_iteration36_深層分析.md`: Iteration 36 の分析結果
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration36.md`: 前回の作業ログ

## ① 検証結果の分析と問題点の特定

### 分析スコープ

1. Iteration 36 の全14チェック検証結果（全PASS確認）
2. `quantum_dynamics_gksl_comparison.ipynb` の全40セルの実行結果
3. 全シミュレータの数学的構造比較

### 発見された問題点

**問題1: Trotterステップ構造の不整合（重大度: 中）**

理想DM版シミュレータと理想ショットベースシミュレータが異なるTrotterステップ構造を使用していた。

- DM版: 半ステップ回文順序 `E_α(dt/2) forward + E_α(dt/2) reverse`
- ショットベース: 全ステップ前方のみ `E_α(dt) forward only`

docstringでは「無限ショット極限でDM版と等価」と記述されていたが、上記構造の差異のため不正確だった。

**問題2: ノイズDMシミュレータの不整合（重大度: 中）**

ノイズDMシミュレータも `_U_stines`（全ステップ）を使用しており、ノイズ=0の極限でも理想DM版と一致しなかった。

**問題3: ゲート数見積もりの不正確性（重大度: 低）**

回文順序による2回適用を考慮していなかった（26 → 52 Stinespring gates/step）。

詳細は `developing/検証結果分析_iteration37_深層分析.md` を参照。

## ② コード修正

### 修正対象ファイル

1. `tutorials/qudit_gksl_shot_simulator.py`
   - `_precompute_unitaries`: `_U_stines_half` の追加
   - `_trotter_step_trajectory`: 回文順序化（半ステップ forward + reverse）
   - `QuditGKSLNoisyShotSimulator._trotter_step_trajectory`: 同上

2. `tutorials/qubit_gksl_shot_simulator.py`
   - `_precompute_unitaries`: `_U_stines_half` の追加
   - `_trotter_step_trajectory`: 回文順序化
   - `QubitGKSLNoisyShotSimulator._trotter_step_trajectory`: 同上

3. `tutorials/qudit_gksl_noisy_simulator.py`
   - `_trotter_step`: `_U_stines` → `_U_stines_half` + 回文順序化

4. `tutorials/qubit_gksl_noisy_simulator.py`
   - `_trotter_step`: `_U_stines` → `_U_stines_half` + 回文順序化

5. `tutorials/qudit_gksl_simulator.py`
   - ゲート数見積もり修正: `4 + 3 + 26` → `4 + 3 + 26 * 2`

6. `tutorials/qubit_gksl_simulator.py`
   - ゲート数見積もり修正: `n_ancilla * 6` → `n_ancilla * 6 * 2`

### テスト確認結果

- `test_gksl_simulators.py`: 全GKSL関連テスト PASS
  - `TestQuditGKSLNoisySimulator::test_zero_noise_matches_ideal` PASSED ✓
  - `TestQubitGKSLNoisySimulator::test_zero_noise_matches_ideal` PASSED ✓
  - 全ノイズチャネルプロパティテスト PASSED ✓
  - 全収束テスト PASSED ✓
- 失敗テストは全て既存の環境依存問題（mqt/qiskit未インストール、メモリ不足）

## ③ 検証スクリプト

`tutorials/run_tta_uc_gksl_verification_iteration37.py` を作成。

### 新規チェック項目

| # | チェック項目 | 目的 |
|---|---|---|
| 15 | Qudit Shot-DM一致性 | 無限ショット極限でDM版と一致するか |
| 16 | Qubit Shot-DM一致性 | 同上（qubit版） |
| 17 | Qudit Noisy(p=0)=DM | ノイズゼロ極限で理想DM版と一致するか |
| 18 | Qubit Noisy(p=0)=DM | 同上（qubit版） |

### 実行方法

```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration37.py
```

結果は `developing/verification_results/` に JSON + Markdown で出力される。

## ④ 次回への引き継ぎ事項

1. ユーザーにてローカルで検証スクリプトを実行し、結果をリポジトリにpush
2. 全18チェックの結果を確認
3. 特にチェック15-18（回文順序一貫性の新規チェック）の結果を確認
4. ノートブック `quantum_dynamics_gksl_comparison.ipynb` の再実行は任意（コアのDM版は変更なし）
