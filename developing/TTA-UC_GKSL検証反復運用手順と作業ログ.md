# TTA-UC現象 GKSL-Lindblad 検証反復運用手順と作業ログ

## 目的

TTA-UC現象のGKSL-Lindblad量子ダイナミクスについて、
ヒューリスティックやfallbackを使わず、実行結果ファイルに基づいて検証・修正を反復する。

## 反復手順（①〜④）

1. **検証スクリプト実行（①）**

   ```bash
   python /home/runner/work/mqt-qudits/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py \
     --output-dir /home/runner/work/mqt-qudits/mqt-qudits/developing/verification_results
   ```

   - 実行結果を以下2ファイルとして必ず出力する。
     - `tta_uc_gksl_verification_<timestamp>.json`
     - `tta_uc_gksl_verification_<timestamp>.md`

2. **ユーザーが結果をpush（②）**

   - ①で生成された結果ファイルをコミットし、対象ブランチにpushする。

3. **結果確認PR作成（③）**

   - push済み結果ファイルを根拠に、修正依頼PRをagent経由で作成する。
   - PR本文には、対象結果ファイル名とFAIL項目を明記する。

4. **結果に基づく修正実施（④）**

   - ③の結果ファイルを確認し、必要な修正を実施する。
   - 修正後は①に戻る。

## 作業ログ運用ルール

- 本ファイル末尾に毎回ログを追記する（日本語Markdown）。
- 新しい作業開始時に、直前ログを参照してから着手する。
- 各ログに最低限以下を記録する。
  - 実行日時
  - 実行コマンド
  - 生成ファイル名
  - 判定（PASS/FAIL）
  - 次アクション

---

## 作業ログ

### 2026-02-16 初回整備

- 実施内容
  - 検証スクリプト `tutorials/run_tta_uc_gksl_verification.py` を追加。
  - 実行時にJSON/Markdownの結果ファイルを出力する仕様を追加。
  - 本手順書を作成し、①〜④反復とログ運用ルールを明記。
- 次アクション
  - ローカルで検証スクリプトを実行し、`developing/verification_results`配下の結果をpushして次の修正サイクルへ進む。
