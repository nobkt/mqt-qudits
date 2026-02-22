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

### 2026-02-22 検証スクリプト大幅拡張（①対応）

- 実施内容
  - `tutorials/run_tta_uc_gksl_verification.py` を `quantum_dynamics_gksl_comparison.ipynb` の全シナリオに対応するように拡張。
  - 対応シナリオ一覧:
    - Scenario 1: classical（非ボソン、N=4）
    - Scenario 2: classical_boson（ボソン、N=2、n_max=1）
    - Scenario 3: qubit（非ボソン、N=4）
    - Scenario 4: qubit_boson（ボソン、N=2、n_max=1）
    - Scenario 5: qudit（非ボソン、N=4）
    - Scenario 6: qudit_boson（ボソン、N=2、n_max=1）
    - Scenario 3b: qubit_shot（ショットベース、ノイズなし）
    - Scenario 3c: qubit_noisy_shot（ショットベース、ノイズあり）
    - Scenario 5b: qudit_shot（ショットベース、ノイズなし）
    - Scenario 5c: qudit_noisy_shot（ショットベース、ノイズあり）
    - unitary（ユニタリ極限、散逸なし）
  - 追加機能:
    - シナリオ間の量子忠実度（fidelity）比較を自動計算
    - Unitary vs GKSL比較分析
    - 詳細なMarkdownテーブル形式のレポート出力
    - トレースバック付きエラーレポート
    - CLI引数: `--n-shots`, `--seed` を追加
  - ヒューリスティックやfallbackは一切使用していない。物理ベースの検証のみ。
- 生成ファイル: なし（スクリプト改修のみ、ユーザーがローカルで実行する）
- 次アクション
  - ユーザーがローカルで `python tutorials/run_tta_uc_gksl_verification.py` を実行し、結果を `developing/verification_results` 配下にpush（②）。
  - push後にagent側で結果を確認し修正PR作成（③④）。
