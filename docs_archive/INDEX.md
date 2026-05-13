# docs_archive/ — 退避された文書・スクリプトの索引

このディレクトリは「過去の作業経緯を保全しつつ、リポジトリ表面の文書堆積を整理する」目的で
作られた退避先です。**削除はしていません**（過去の検討経緯と照合できる状態を維持しています）。

E-1（文書の堆積）と E-2（文書間の相互矛盾）への本 PR での対応は以下のとおりです。

## 何を退避したか

| 退避先サブディレクトリ | 内容 | 退避前の場所 | 件数 |
|---|---|---|---|
| `root_md/` | リポジトリ直下に堆積していた "*_FIX_REPORT*.md", "*_COMPLETION_*.md", "*_SUMMARY*.md", "*_OLD.md", "H_TTA_*", "PR<N>_*完了報告.md" 等の解決済み報告系 .md（55+ 本） | リポジトリ直下 | 69 |
| `developing/` | iteration ごとの作業ログ `TTA-UC_GKSL検証反復運用手順と作業ログ_iterationN.md` および古い設計書改訂履歴・進捗報告書 | `developing/` | 130 |
| `iteration_scripts/` | `run_iteration*_verification.py`、`run_tta_uc_gksl_verification_iterationN.py`、および iteration29 の issue 別検証スクリプト | `tutorials/` | 62 |
| `notebook_backups/` | ノートブックの `.backup*` 系（4 本） | `tutorials/` | 4 |
| `tutorials_md/` | `tutorials/` 直下の "*_SUMMARY*.md", "*_REPORT*.md", "PR41_IMPLEMENTATION_SUMMARY.md" 等 | `tutorials/` | 15 |
| `tutorials_doc/` | `tutorials/doc/` 配下の "PR<N>_COMPLETION_REPORT*.md", "PR<N>_IMPLEMENTATION_SUMMARY*.md" など PR 単位の完了報告群 | `tutorials/doc/` | 64 |

## 何を退避していないか（リポジトリに残しているもの）

- `README.md`、`QUICK_START.md`：ユーザー向けの導入資料。
- `LICENSE`、`pyproject.toml`、ソースコード（`src/`, `include/`, `cmake/`, `tools/`, `test/`）。
- `developing/TTA-UC現象GKSL-Lindblad量子ダイナミクス詳細実装計画書.md`：能動的な実装計画書として残置。
- `developing/verification_results/`：検証結果ファイル（成果物）として残置。
- `tutorials/README.md`：チュートリアル一覧として残置。
- `tutorials/doc/GKSL/`：理論・仕様の参照ドキュメント群として残置。
- `tutorials/run_quantum_dynamics_verification.py`、`run_qudit_ancilla_verification.py`、
  `run_tta_uc_gksl_verification.py`、`run_tta_uc_gksl_regression_test.py`：iteration 番号を持たない
  能動的な検証エントリポイントとして残置。

## なぜ削除ではなく退避なのか

- 退避された文書群の中には、過去の修正で**「解決済み」と宣言された問題が後の検討で覆っているケース**が
  少なからずあり、それらは「リポジトリ全体としての主張」を再評価する際の重要な一次資料になります（E-2）。
- `git rm` で消すと検索性は上がりますが、過去の検討経緯と現状の食い違いを後から検証する手段が
  失われます。そのコストは、削除で得られる利点を上回ると判断しました。

## 新規 "完了報告" .md を増やさない方針

本 PR ではこのインデックスと、`STATUS_HONEST_2026-05.md`（本 PR で「解決した／部分的に対処した／
未解決のまま明示する」事項を表で記載）の **2 ファイルだけ** を新規追加しています。
従来のような "PR<N>_COMPLETION_REPORT_JA.md" 形式の文書は本 PR では新規作成しません（E-1 への対応）。
