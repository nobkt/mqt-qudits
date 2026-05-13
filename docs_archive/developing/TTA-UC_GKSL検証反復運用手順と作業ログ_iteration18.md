# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 18

## 基本情報

- 日時: 2026-02-24
- 前回: Iteration 17（修正版ノイズ診断）
- 目的: Iteration 17 結果の詳細分析で発見された検証方法論の欠陥を修正

## 前回（Iteration 17）からの引き継ぎ

### 確認済み事項
- ✅ シミュレータコード（8ファイル）にバグなし
- ✅ Test B（ノイズバジェット分解）正常
- ✅ Test C（ε(n) 追跡）正常 — 収穫逓減効果を正しく実証
- ✅ Test D（Qubit忠実度追跡）正常

### 発見された問題
1. **重大**: Stinespring近似誤差が未検証（Test A は Trotter 分割誤差のみ測定）
2. **軽微**: Test C の dt 依存性報告がハードコードで不正確（~0.1% vs 実際 ~0.9%）
3. **中程度**: 密度行列の妥当性検証（トレース・正値性等）が欠如

## 作業手順

### ステップ 1: 分析ドキュメント作成 ✅
- `developing/検証結果分析_20260224_iteration18.md` を作成
- 3 件の問題を詳細に記述

### ステップ 2: 検証スクリプト作成
- `tutorials/run_tta_uc_gksl_verification_iteration18.py` を作成
- 変更内容:
  - **Test A〜D**: Iteration 17 と同一（変更なし）
  - **Test C**: dt 変動量をデータから計算するように修正（ハードコード削除）
  - **Test E（新規）**: 正確な GKSL Liouvillian 解との比較
    - `ClassicalGKSLSimulator` を参照解として使用
    - 複数 dt 値で Stinespring+Trotter の全体誤差を測定
    - 密度行列診断（トレース、エルミート性、最小固有値）を記録
    - 収束次数を数値的に算出

### ステップ 3: ユーザーによるスクリプト実行
- ユーザーがローカルで `run_tta_uc_gksl_verification_iteration18.py` を実行
- 結果は `developing/verification_results/` に自動保存

### ステップ 4: 結果確認と修正
- Iteration 18 の結果を分析
- Stinespring 近似誤差の大きさに応じて追加修正が必要か判断

## コード変更方針

### 変更するファイル
- `tutorials/run_tta_uc_gksl_verification_iteration18.py`（新規作成）

### 変更しないファイル（シミュレータコア）
- `tutorials/qudit_gksl_simulator.py` — 変更なし
- `tutorials/qudit_gksl_noisy_simulator.py` — 変更なし
- `tutorials/qubit_gksl_simulator.py` — 変更なし
- `tutorials/qubit_gksl_noisy_simulator.py` — 変更なし
- `tutorials/stinespring_utils.py` — 変更なし
- `tutorials/gksl_math_utils.py` — 変更なし
- `tutorials/gksl_physical_parameters.py` — 変更なし
- `tutorials/classical_gksl_simulator.py` — 変更なし
