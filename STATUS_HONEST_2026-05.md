# STATUS_HONEST_2026-05.md — 本 PR で解決した／していない事項の正直な記録

このリポジトリは過去に多くの「完了報告」.md を生成してきましたが、それらの主張は
互いに矛盾しているものが少なくありません（E-1/E-2）。本ファイルでは、ご指摘いただいた
**A-1 / A-2 / A-3 / A-4 / B-1 / D-1 / E-1 / E-2** の各項目について、本 PR で何を
解決し、何を解決していないかを「嘘・誤魔化し・ヒューリスティックを使わずに」記載します。

## 本 PR の対応状況一覧

| ID | 内容 | 状況 | 詳細 |
|---|---|---|---|
| **A-1** | MQT-Qudits を実際には使っていない（NumPy `expm` のみで計算） | ⚠️ **部分解決（次 PR で継続）** | 6 シナリオの本体は依然 `scipy.linalg.expm` 直接計算。**ただし**、`tutorials/run_n2_tnsim_verification.py` で N=2 最小構成における Trotter ステップの **全部品回路**（Hamiltonian 半ステップ + 単一サイト Stinespring 10 個 + TTA-pair Stinespring 2 個 = 計 13 サブ回路）を MQT-Qudits の `tnsim` バックエンドで実行し、NumPy `expm` 期待値と完全一致 (||Δsv|| = 0.0, tol 1e-10) を確認・自動テスト化 (`test_gksl_simulators.py::TestN2TnsimVerification`)。組合せステップ回路は依然 mid-circuit ancilla reset 不在のため `tnsim`/`misim` 上で実行不可 — これが残課題として明示。 |
| **A-2** | Qudit/Qubit 比較に独立性が無い | ✅ **正直化のみ実施** | Qubit/Qudit 両シミュレータの docstring と NB Cell 7/13/25/35 に「両者は同一の Stinespring/Trotter コアを共有し、qutrit→qubit 写像で埋め込んだもの。両者一致は写像の正しさの確認であり、独立な相互検証ではない」と明記した。**真に独立な 2 実装を作ることは本 PR では行っていない**。 |
| **A-3** | ゲート数が `2*(N+pairs)+...` のハードコード式 | ✅ **解決（コンパイラ実測値を提供）** | `qudit_gksl_simulator.py` のハードコード式は backward-compat で残るが、`simulate()` 戻り値に `n_high_level_gates_per_step` および `gate_count_method="high_level_count"` を追加し「これは高レベルゲートのオブジェクト数であってコンパイラ計測値ではない」ことを honestly に明示。新規メソッド `QuditGKSLSimulator.compute_compiler_measured_gate_counts(dt, optimization_level)` を追加し、MQT-Qudits の `compileO0/compileO1` を実際に呼んでネイティブゲート（VirtRz/R/Rh/Rz/CEx）の内訳を返す。N=2, dt=0.5, opt=0 では native gates per step = 10,836（ハードコード値 30 と一桁以上違うことが実測で示される）。`cu_multi` (TTA-pair Stinespring, 27x27) は MQT-Qudits compiler 現状で decompose 不可という制約を `per_step_summary['uncompiled_cu_multi']` として正直に記録（ヒューリスティック推定では埋めない）。NB に専用セルを追加し全シナリオの内訳を表示。テスト: `test_gksl_simulators.py::TestCompilerMeasuredGateCounts`。 |
| **A-4** | ノイズ比較の根拠不足・非対称（qubit 側だけ `cx_per_pair_gate=46`） | ✅ **解決** | NB Cell 38 から CX 倍率を撤去し、qubit と qudit を**同一規約**（ペアゲート 1 つあたり `p_depol` を 1 回適用）で比較するよう変更。Cell 37 markdown に「片側だけの増幅は qudit 有利の結論を恣意的に作り出すバイアス」と明記。`cx_per_pair_gate` パラメータ自体は API として残し、明示的にオプトインで利用可能。 |
| **B-1** | Stinespring が 1 次精度のため Strang 2 次が潰れて O(dt) | ✅ **解決（exact_local_channels アルゴリズムで 2 次収束を実測）** | 新規モジュール `tutorials/exact_local_channels.py` を追加し、各 Lindblad チャネルを Stinespring 近似ではなく**局所超演算子の厳密 exponentiation** `expm(L_D_α^local · dt)`（単一サイト 9×9, ペア 81×81 の行列）として直接適用する scheme を実装。`QuditGKSLSimulator(algorithm="exact_local_channels")` の opt-in パラメータで切替可能（default は `"stinespring"` のまま、後方互換完全保持）。**実測** (N=2, t_max=10, n_steps=20→200): rate = 2.000 を**全範囲で確認**（既存 `"stinespring"` は同条件で rate ≈ 1.000）。N=4, t_max=100, n_steps=100 では Stinespring T ≈ 7e-4 に対し exact T ≈ 4.7e-5（**約 15 倍精度改善**）。テスト: `test_gksl_simulators.py::TestExactLocalChannelsConvergence` (rate ≥ 1.9 を assert)。NB Cell 28 を両アルゴリズム並列出力に書き換え、Cell 27/29 を更新。**注意**: `"exact_local_channels"` は古典シミュレータ上で正確な dissipator 適用を行うものであり、対応する量子回路の構築（mid-circuit ancilla reset 等）は A-1 の残課題のまま。 |
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

1. 本ノートブック (`tutorials/quantum_dynamics_gksl_comparison.ipynb`) の主要な
   時間発展（Cell 6 の `QuditGKSLSimulator.simulate`、Cell 8 の Qubit GKSL、ボソン
   付き Cell 12/14 等）は依然 NumPy/`scipy.linalg.expm` ベースです。本 PR で
   `tutorials/run_n2_tnsim_verification.py` により N=2 の **部品回路レベル**では
   MQT-Qudits backend 上で正しく実行できることを**実測で**示しましたが、
   **組合せ per-step 回路の MQT-Qudits backend 実行**は mid-circuit ancilla
   reset の不在により依然不可能です。これが A-1 の残課題です。
2. ゲート数の表示には現在 2 系統あります：
   - `result["estimated_gates_per_step"]` / `result["n_high_level_gates_per_step"]`
     は依然「ハードコードした高レベルゲート数式」の値で、A-3 の旧表示と同じ値です。
     後方互換のため残してあります。`result["gate_count_method"] = "high_level_count"`
     によりこれが**コンパイラ計測値ではない**ことを明示しています。
   - `QuditGKSLSimulator.compute_compiler_measured_gate_counts(dt)` を呼べば
     MQT-Qudits の `compileO0`/`compileO1` を実際に通したネイティブゲート
     (`VirtRz`/`R`/`Rh`/`Rz`/`CEx`) の内訳と、`cu_multi` (TTA-pair Stinespring) の
     未分解個数（MQT-Qudits compiler が現状 `cu_multi` を 2-qudit gate に
     decompose できないため）を**そのまま**返します。NB に内訳を表示する
     セルを追加済みです。これが A-3 の対応です。
3. 収束テーブル（NB Cell 28〜30）には `algorithm="stinespring"` (1 次, default)
   と `algorithm="exact_local_channels"` (2 次, 本 PR で追加) の両方の rate が
   並列出力されます。Stinespring 側の `rate ≈ 1.0` は依然として事実
   （Stinespring dilation の数学的性質）であり、Cell 28 のその記述は維持。
   Exact local channel 側の `rate ≈ 2.0` も**実測**であり、嘘・ヒューリスティック
   なしで Strang+palindromic 構造が本来期待する 2 次収束を取り戻します。
   ただしこれは古典シミュレータ上での厳密 dissipator 適用であり、対応する
   量子回路化（Choi-Kraus からの Kraus operator 抽出 + ancilla reset 不要な
   実装等）は依然として未完了で、A-1 の残課題に含まれます。
