# STATUS_HONEST_2026-05.md — 本 PR で解決した／していない事項の正直な記録

このリポジトリは過去に多くの「完了報告」.md を生成してきましたが、それらの主張は
互いに矛盾しているものが少なくありません（E-1/E-2）。本ファイルでは、ご指摘いただいた
**A-1 / A-2 / A-3 / A-4 / B-1 / D-1 / E-1 / E-2** の各項目について、本 PR で何を
解決し、何を解決していないかを「嘘・誤魔化し・ヒューリスティックを使わずに」記載します。

## 本 PR の対応状況一覧

| ID | 内容 | 状況 | 詳細 |
|---|---|---|---|
| **A-1** | MQT-Qudits を実際には使っていない（NumPy `expm` のみで計算） | ❌ **本 PR では対応しない** | 6 シナリオの本体は依然 `scipy.linalg.expm` 直接計算であり、MQT-Qudits backend (`tnsim`/`misim`) 上での実行は本 PR でも行っていない。`*KrausSimulator.build_combined_trotter_step_circuit` で組まれる回路は mid-circuit ancilla reset を含むため backend で実行不可（その旨を当該 docstring に明記）。次 PR で N=2 の最小構成について `tnsim` 経由の検証テストを追加予定。 |
| **A-2** | Qudit/Qubit 比較に独立性が無い | ✅ **正直化のみ実施** | Qubit/Qudit 両シミュレータの docstring と NB Cell 7/13/25/35 に「両者は同一の Stinespring/Trotter コアを共有し、qutrit→qubit 写像で埋め込んだもの。両者一致は写像の正しさの確認であり、独立な相互検証ではない」と明記した。**真に独立な 2 実装を作ることは本 PR では行っていない**。 |
| **A-3** | ゲート数が `2*(N+pairs)+...` のハードコード式 | ❌ **本 PR では対応しない** | 当該ハードコード式はそのまま残っている。`*KrausSimulator` の docstring に「報告するゲート数は回路に追加した高レベルゲート数であり、コンパイラ計測値ではない」と明記。次 PR で MQT-Qudits compiler を実際に呼んだ実測値への置換を試みる予定。 |
| **A-4** | ノイズ比較の根拠不足・非対称（qubit 側だけ `cx_per_pair_gate=46`） | ✅ **解決** | NB Cell 38 から CX 倍率を撤去し、qubit と qudit を**同一規約**（ペアゲート 1 つあたり `p_depol` を 1 回適用）で比較するよう変更。Cell 37 markdown に「片側だけの増幅は qudit 有利の結論を恣意的に作り出すバイアス」と明記。`cx_per_pair_gate` パラメータ自体は API として残し、明示的にオプトインで利用可能。 |
| **B-1** | Stinespring が 1 次精度のため Strang 2 次が潰れて O(dt) | ❌ **本 PR では対応しない** | 既存の docstring と NB Cell 27/28 における「O(dt)（1次）が支配的」「これはバグではない」の記述は数学的事実なのでそのまま残す。次 PR で 2 次化（exp(L_α dt) を Choi-Kraus 経由で正確適用 等）の実装を試行し、収束次数が実測 2 になった場合のみ Cell 28 を差し替える。**実測で 2 にならなかった場合は NB を変更しない**。 |
| **D-1** | "circuit simulator" の名称詐称 | ✅ **解決** | `QuditGKSLCircuitSimulator` → `QuditGKSLKrausSimulator`、`QubitGKSLCircuitSimulator` → `QubitGKSLKrausSimulator`、`QuditGKSLCircuitBosonSimulator` → `QuditGKSLKrausBosonSimulator` に改名。各クラス・モジュールの docstring 冒頭に「これは MQT-Qudits / Qiskit 回路を構築するが、回路実行ではなく Kraus 演算子を抽出して NumPy で密度行列に適用するシミュレータである。`build_combined_trotter_step_circuit` で組まれる統合回路は MQT-Qudits backend では実行不可（mid-circuit ancilla reset 未対応）であり、可視化目的に限定される」を明記。後方互換のため旧名のエイリアスを残置。 |
| **E-1** | 文書の堆積 | ✅ **解決** | リポジトリ直下 .md（69 本）、`developing/` の iteration 系（130 本）、`tutorials/` の `run_iteration*` / `run_tta_uc_gksl_verification_iteration*` スクリプト（62 本）、`tutorials/*.md` の summary 系（15 本）、`tutorials/doc/` の PR 完了報告群（64 本）、ノートブック `.backup*`（4 本）の合計 **344 ファイル** を `docs_archive/` 配下へ `git mv` で退避。詳細は [`docs_archive/INDEX.md`](docs_archive/INDEX.md)。**削除はしていない**（過去の検討経緯の保全）。 |
| **E-2** | 文書間の相互矛盾 | ⚠️ **完全には解決していない** | `docs_archive/` への退避により、**現役を主張するドキュメント（リポジトリ表面）は `README.md` / `QUICK_START.md` / 本ファイル / `docs_archive/INDEX.md` および `developing/`・`tutorials/doc/GKSL/` の参照系のみ**になり、相互矛盾の表面積は大幅に縮小した。ただし `docs_archive/` 内の旧 .md の主張同士は依然として矛盾しており、これは保全のため意図的にそのままにしている。新規の "PR<N>_COMPLETION_REPORT" 系 .md は本 PR では作成しない方針を明示する（本ファイルは唯一の例外）。 |

## 用語

- **解決**：本 PR の変更で当該問題が事実として解消された、または検証可能な形で正直化された。
- **正直化のみ実施**：根本実装は変えていないが、誤った主張・誇大な主張をコード/文書から削除し、
  事実に即した記述に置き換えた。
- **本 PR では対応しない**：当該問題には触れていない。次 PR で扱うことを明示。
- **完全には解決していない**：部分的に対処したが残課題があることを明示。

## 嘘をつかないために — 残課題の明示

1. 本ノートブック (`tutorials/quantum_dynamics_gksl_comparison.ipynb`) が「MQT-Qudits を使った
   qudit 量子計算の実例」であるかのように見える可能性は依然残っています。実態は「MQT-Qudits の
   `QuantumCircuit` API を使って回路オブジェクトを*構成*はするが、時間発展の数値計算は NumPy
   と `scipy.linalg.expm` で行っている」です。これは A-1 の解決を待つべき問題です。
2. ゲート数の表示（`Estimated gates/step: 66` 等）は実際にコンパイラを通した値ではなく、
   コードがハードコードした式の結果です。これは A-3 の解決を待つべき問題です。
3. 収束テーブル（NB Cell 28）の rate≈1.0 は事実であり、「これはバグではない」という記述も
   事実ですが、それは「現在の実装では 2 次精度が原理的に出ない」ことを意味しており、
   「実用上の精度は十分」を意味しません。これは B-1 の解決を待つべき問題です。
