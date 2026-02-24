# TTA-UC GKSL 検証反復運用手順と作業ログ Iteration 15

## 概要

- 反復番号: 15
- 日時: 2026-02-24
- 前回: iteration 14 (developing/検証結果分析_20260224_iteration14.md)
- 目的: Trotter収束性検証、厳密GKSL解との比較、DMノイジーvsショットノイジー整合性検証

## 背景

Iteration 14の検証結果を詳細に分析した結果、コード実装自体にはバグは
検出されなかったが、検証方法に以下の問題点が発見された:

1. Trotter近似（dt=1.0）の収束性が未検証
2. DM ノイジーとショットノイジーの整合性が未検証（特に多ステップ）
3. 厳密GKSL解との比較が未実施
4. F_norm メトリクスがleakageを隠蔽する問題
5. forbidden state population の追跡不足

## 作業手順

### ステップ1: 分析ドキュメント作成

`developing/検証結果分析_20260224_iteration15.md` を作成:
- iteration 14 結果の詳細分析
- 5つの問題点の特定と説明
- コード実装の正確性確認

### ステップ2: 検証スクリプト作成

`tutorials/run_tta_uc_gksl_verification_iteration15.py` を作成:

- **Test 1: Trotter収束テスト** — n_steps=10,20,50,100 でノイズレス qudit を
  実行し、populations の収束を確認
- **Test 2: 厳密GKSL解との比較** — `build_gksl_superoperator` + `expm(L·t)` で
  厳密解を計算し、Trotter解と比較
- **Test 3: DMノイジー10ステップ忠実度** — DMノイジーシミュレータで10ステップを
  実行し、各ステップで忠実度・forbidden state population を計算
- **Test 4: DMノイジー vs ショットノイジー整合性** — 同じパラメータでDMノイジー
  とショットベースを比較し、整合性を確認

### ステップ3: 作業ログ作成

本ファイル `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration15.md`

### ステップ4: ユーザーによるローカル検証（次回）

ユーザーがローカルで以下を実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration15.py
```

結果が `developing/verification_results/` に出力される。

## 修正ファイル一覧

| ファイル | 変更種別 | 内容 |
|---------|---------|------|
| tutorials/run_tta_uc_gksl_verification_iteration15.py | 新規 | iteration15 検証スクリプト |
| developing/検証結果分析_20260224_iteration15.md | 新規 | iteration14 結果分析 |
| developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration15.md | 新規 | 本作業ログ |

## コード変更について

シミュレータ本体のコード（qudit_gksl_noisy_simulator.py 等）には変更なし。
iteration 14 のコード修正は正確であり、バグは検出されなかった。
今回は検証方法の改善（新しい検証スクリプト）のみ。

## 次回の作業

1. ユーザーがローカルで iteration15 検証スクリプトを実行
2. 結果をリポジトリにpush
3. 検証結果を確認し、Trotter誤差が問題であれば:
   - dt を小さくした（n_steps を増やした）ノイジーシミュレーションを追加
   - 必要に応じてデフォルトの n_steps 推奨値をドキュメントに追記
4. DM ノイジー vs ショットノイジーの整合性に問題があれば:
   - 該当シミュレータのコードを修正
