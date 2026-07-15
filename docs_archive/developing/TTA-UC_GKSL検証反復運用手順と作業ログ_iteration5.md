# TTA-UC GKSL検証反復運用手順と作業ログ（iteration5）

## 反復の概要

- **反復番号**: 5
- **日時**: 2026-02-22
- **前回反復**: iteration4（全11シナリオPASS、qubit 256次元空間対応完了）
- **目的**: iteration4の検証結果とコード全体を詳細分析し、残存する問題を修正

## 手順1: 検証結果分析

### 入力
- `developing/verification_results/tta_uc_gksl_verification_20260222T135417Z.json`
- `developing/verification_results/tta_uc_gksl_verification_20260222T135417Z.md`

### 分析結果

iteration4の最終検証は全11シナリオPASS。qubitシミュレータが256次元空間で正しく
動作していることを、エントロピー値の浮動小数点精度差（~1e-15レベル）で確認した。

以下の4つの問題を発見:

1. **問題1（重要）**: 密度行列シミュレータ4つでStinespringユニタリとハミルトニアンexpmが
   毎Trotterステップで冗長に再計算されている
2. **問題2（中程度）**: QubitGKSLNoisyShotSimulatorにディフェージングチャネルがない
3. **問題3（軽微）**: リーケージシナリオのトレース誤差報告がmax_trace_error=0と誤解を招く
4. **問題4（軽微）**: qubit/quditノイズシミュレータで異なるノイズ種別を使用

詳細は `developing/検証結果分析_20260222_iteration5.md` に記載。

## 手順2: コード修正

### 修正1: Stinespring/Hamiltonian expm事前計算の追加

**対象ファイル（4ファイル）:**
- `tutorials/qudit_gksl_simulator.py`
- `tutorials/qubit_gksl_simulator.py`
- `tutorials/qudit_gksl_boson_simulator.py`
- `tutorials/qubit_gksl_boson_simulator.py`

**変更内容:**
- `_apply_hamiltonian_step(rho, dt)` と `_apply_lindblad_stinespring(rho, L_op, dt)` を削除
- `_precompute_unitaries(dt)` メソッドを追加（shot系シミュレータと同様のインターフェース）
- `_trotter_step(rho, dt)` を `_trotter_step(rho)` に変更（キャッシュされたユニタリを使用）
- `simulate()` メソッド内で `_precompute_unitaries(dt)` を1回だけ呼び出し

**効果:**
- expm呼び出し回数: 140回 → 27回（約80%削減）
- 計算時間の大幅な短縮が期待される

### 修正2: QubitGKSLNoisyShotSimulatorへのディフェージング追加

**対象ファイル:** `tutorials/qubit_gksl_shot_simulator.py`

**変更内容:**
- `_apply_stochastic_qubit_dephasing_single()` 関数を追加
  - 4次元局所空間（2量子ビット符号化）での計算基底への確率的射影
  - QuditGKSLNoisyShotSimulatorの `_apply_stochastic_dephasing_single()` と対称的な設計
- `QubitGKSLNoisyShotSimulator.__init__()` に `p_dephasing` パラメータを追加
- `_trotter_step_trajectory()` にディフェージング処理を追加（quditと同構造）
- `noise_params` に `p_dephasing` を含むように更新

### 修正3: 検証スクリプトのトレース報告改善

**対象ファイル:** `tutorials/run_tta_uc_gksl_verification.py`

**変更内容:**
- `_validate_result()` に `max_trace_deficit` メトリクスを追加
  - リーケージシナリオで `1 - trace` の最大値を報告
- `_build_markdown_report()` に `max_trace_deficit` の表示を追加
- qubit_noisy_shotのノイズパラメータを更新: `p_dephasing=0.005` を追加

## 手順3: 検証

### ローカル検証結果

全シミュレータの基本動作を確認:
- QuditGKSLSimulator: ✅ (precomputed unitaries)
- QubitGKSLSimulator: ✅ (precomputed unitaries, 256-dim confirmed)
- QuditGKSLBosonSimulator: ✅ (precomputed unitaries)
- QubitGKSLBosonSimulator: ✅ (precomputed unitaries)
- QubitGKSLNoisyShotSimulator (p_dephasing=0): ✅ (backward compatible)
- QubitGKSLNoisyShotSimulator (p_dephasing=0.005): ✅
- max_trace_deficit metric: ✅ (qubit_noisy_shot: 0.259, qudit_noisy_shot: N/A)

### 検証スクリプト実行
- `--scenarios classical qudit qubit`: 全PASS
- `--scenarios qubit_noisy_shot qudit_noisy_shot`: 全PASS
- qubit_noisy_shot: max_trace_deficit正常に報告

## 次のステップ

ユーザーによるフル検証スクリプトの実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification.py
```

全11シナリオの結果を `developing/verification_results/` にpushし、
iteration6で結果を確認する。
