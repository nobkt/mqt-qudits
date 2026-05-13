# STATUS_HONEST_2026-05.md — 累積残課題の正直な記録

このリポジトリは過去に多くの「完了報告」.md を生成してきましたが、それらの主張は
互いに矛盾しているものが少なくありません（E-1/E-2）。本ファイルは、ご指摘いただいた
**A-1 / A-2 / A-3 / A-4 / B-1 / D-1 / E-1 / E-2** の各項目について、PR#257 〜
本セッションまでに何を解決し、何を解決していないかを「嘘・誤魔化し・ヒューリスティックを
使わずに」記載した累積記録です。

## 本ファイルの履歴（嘘をつかないため）

| 改訂 | 日付 | 反映範囲 |
|---|---|---|
| 初版 | PR#257 | D-1 / A-2 (正直化) / A-4 / E-1・E-2 を解決、A-1 / A-3 / B-1 を未解決として記録 |
| 第 2 版 | PR#258 | A-1 を「N=2 部品回路を tnsim で実測検証 (||Δsv||=0)」、A-3 を「コンパイラ実測値 API 追加」、B-1 を「`exact_local_channels` で rate=2 実測」に更新 |
| 第 3 版 | 本セッション (2026-05-13) | PR#259 (QuTiP 独立リファレンス) と PR#260 (fresh-ancilla per-step + `interleaved_n2` で N=2 full per-step を tnsim 実行成功) を反映、本セッションで再実測した全数値を追記、N>2 の interleaved 拡張は state-vector 自体が infeasible で実装しても無効である事実を追記 |

## 累積対応状況一覧

| ID | 内容 | 状況 | 詳細 |
|---|---|---|---|
| **A-1** | MQT-Qudits を実際には使っていない（NumPy `expm` のみで計算） | ⚠️ **大幅前進・但し根本未解決** | (i) 6 シナリオの本体時間発展は依然 `scipy.linalg.expm` 直接計算。(ii) **PR#258**: `tutorials/run_n2_tnsim_verification.py` で N=2 の Trotter ステップ全 13 サブ回路（Hamiltonian 半ステップ + 単一 Stinespring×10 + TTA-pair Stinespring×2）を `tnsim` で実行し ‖Δsv‖=0.0 (tol 1e-10) を確認 (`test_gksl_simulators.py::TestN2TnsimVerification` パス)。(iii) **PR#259**: `QuditGKSLKrausSimulator.build_executable_per_step_circuit_fresh_ancillas` を新設し、各 Stinespring チャネルに新規 ancilla を割当てて mid-circuit reset を不要化。Part A (math equivalence) で ‖Δρ_sys‖_F = 0.0 を確認 (`run_n2_full_step_fresh_ancilla_verification.py`)。(iv) **PR#260**: `ancilla_layout="interleaved_n2"` を追加（**qudit 添字の純粋な permutation のみ**でヒューリスティック圧縮ではない）、worst gate range を 13→6 に削減。これにより **N=2 forward-only full per-step (14 qudits, 12 channels) を tnsim 上で実行成功**。本セッションで end-to-end 再実測：Part A 1.2s ✓ (‖Δρ‖_F=0.0)、Part C 1.3s ✓ (‖ρ_sys(tnsim) − ρ_sys(Kraus)‖_F = 0.0、Tr=1.0)、Part B 11 分・5 中 3 成功 + 残り 2 は MemoryError (sequential layout の 12-/14-qudit 中間行列 468 GiB)。**根本的な未解決**: `palindromic=True` fresh-ancilla は state-vector が 3^26≈2.5e12 amplitudes (76 TB) で全 state-vector backend で infeasible、組合せ Strang 半ステップ済 backend 実行は依然 mid-circuit ancilla reset 不在で不可能。N>2 の interleaved 拡張は forward-only でも state-vector が N=3:3^22≈500GB / N=4:3^30≈3PB で **layout を変えても tnsim 実行不可**（本セッションで定量確認）であり、A-1 の根本解決には backend 側 mid-circuit ancilla reset の実装か、新たな量子回路化アルゴリズムが必要。 |
| **A-2** | Qudit/Qubit 比較に独立性が無い | ✅ **解決（独立リファレンスを追加）** | (i) **PR#257**: docstring と NB Cell 7/13/25/35 に「Qubit/Qudit 両シミュレータは同一 `stinespring_utils`/Trotter コア共有、独立検証ではない」を明記。(ii) **PR#259**: `tutorials/qutip_gksl_reference.py` を新設し、in-tree GKSL スタックとコード共有なしの `qutip.mesolve` (適応 Adams/BDF ODE) を独立リファレンスとして導入。同じ物理 Hamiltonian と 5N+2|neighbors| Lindblad collapse operator を QuTiP プリミティブから独立に再構築。(iii) `tutorials/run_qutip_cross_validation.py` で実測。**本セッションで再実測**: N=2/n_steps=200/t_max=10 で `Tr|ρ_qutip − ρ_exact_local|/2 = 1.786e-08`、N=4 同条件で `7.683e-07`（PR#260 主張値 1.79e-08, 7.68e-07 と一致）。`exact_local_channels` の妥当性を独立スタックで確認。 |
| **A-3** | ゲート数が `2*(N+pairs)+...` のハードコード式 | ✅ **解決（コンパイラ実測値を提供）** | `qudit_gksl_simulator.py` のハードコード式は backward-compat で残るが、`simulate()` 戻り値に `n_high_level_gates_per_step` と `gate_count_method="high_level_count"` を追加し「これは高レベルゲートのオブジェクト数であってコンパイラ計測値ではない」ことを honestly に明示。新規メソッド `QuditGKSLSimulator.compute_compiler_measured_gate_counts(dt, optimization_level)` が MQT-Qudits の `compileO0/compileO1` を実際に呼んでネイティブゲート (VirtRz/R/Rh/Rz/CEx) の内訳を返す。**本セッションで再実測** (N=2, dt=0.5, opt=0): `per_step_summary = {'native_gates': 10836, 'uncompiled_cu_multi': 4, 'optimization_level': 0, 'backend': 'faketraps2trits'}`（ハードコード値 30 と 360 倍以上の乖離が実測で示される）。`cu_multi` (TTA-pair Stinespring, 27×27) は MQT-Qudits compiler が現状 decompose 不可という制約を `per_step_summary['uncompiled_cu_multi']` として正直に記録（ヒューリスティック推定では埋めない）。テスト: `TestCompilerMeasuredGateCounts` パス。 |
| **A-4** | ノイズ比較の根拠不足・非対称（qubit 側だけ `cx_per_pair_gate=46`） | ✅ **解決** | NB Cell 38 から CX 倍率を撤去し、qubit と qudit を**同一規約**（ペアゲート 1 つあたり `p_depol` を 1 回適用）で比較するよう変更。Cell 37 markdown に「片側だけの増幅は qudit 有利の結論を恣意的に作り出すバイアス」と明記。`cx_per_pair_gate` パラメータ自体は API として残し、明示的にオプトインで利用可能。 |
| **B-1** | Stinespring が 1 次精度のため Strang 2 次が潰れて O(dt) | ✅ **解決（`exact_local_channels` で rate=2 実測）** | 新規モジュール `tutorials/exact_local_channels.py` で各 Lindblad チャネルを Stinespring 近似ではなく**局所超演算子の厳密 exponentiate** `expm(L_D_α^local · dt)`（単一サイト 9×9, ペア 81×81）として直接適用。`QuditGKSLSimulator(algorithm="exact_local_channels")` の opt-in パラメータで切替（default は `"stinespring"` のまま、後方互換完全保持）。**実測** (N=2, t_max=10, n_steps=20→200): rate = 2.000 を**全範囲で確認**（既存 `"stinespring"` は同条件で rate ≈ 1.000、これは Stinespring dilation の数学的性質として事実）。N=4, t_max=100, n_steps=100 では Stinespring T ≈ 7e-4 に対し exact T ≈ 4.7e-5（**約 15 倍精度改善**）。本セッションでテスト `TestExactLocalChannelsConvergence` (6 件) パスを再確認。**注意（残課題）**: `"exact_local_channels"` は古典シミュレータ上の厳密 dissipator 適用であり、対応する量子回路化（Choi-Kraus からの Kraus operator 抽出 + ancilla reset 不要な実装等）は A-1 の残課題のまま。 |
| **D-1** | "circuit simulator" の名称詐称 | ✅ **解決** | `QuditGKSLCircuitSimulator` → `QuditGKSLKrausSimulator`、`QubitGKSLCircuitSimulator` → `QubitGKSLKrausSimulator`、`QuditGKSLCircuitBosonSimulator` → `QuditGKSLKrausBosonSimulator` に改名。各クラス・モジュールの docstring 冒頭に「これは MQT-Qudits / Qiskit 回路を構築するが、回路実行ではなく Kraus 演算子を抽出して NumPy で密度行列に適用するシミュレータである。`build_combined_trotter_step_circuit` で組まれる統合回路は MQT-Qudits backend では実行不可（mid-circuit ancilla reset 未対応）であり、可視化目的に限定される」を明記。後方互換のため旧名のエイリアスを残置。 |
| **E-1** | 文書の堆積 | ✅ **解決** | リポジトリ直下 .md（69 本）、`developing/` の iteration 系（130 本）、`tutorials/` の `run_iteration*` / `run_tta_uc_gksl_verification_iteration*` スクリプト（62 本）、`tutorials/*.md` の summary 系（15 本）、`tutorials/doc/` の PR 完了報告群（64 本）、ノートブック `.backup*`（4 本）の合計 **344 ファイル** を `docs_archive/` 配下へ `git mv` で退避。詳細は [`docs_archive/INDEX.md`](docs_archive/INDEX.md)。**削除はしていない**（過去の検討経緯の保全）。 |
| **E-2** | 文書間の相互矛盾 | ⚠️ **完全には解決していない** | `docs_archive/` への退避により、**現役を主張するドキュメント（リポジトリ表面）は `README.md` / `QUICK_START.md` / 本ファイル / `docs_archive/INDEX.md` および `developing/`・`tutorials/doc/GKSL/` の参照系のみ**になり、相互矛盾の表面積は大幅に縮小した。ただし `docs_archive/` 内の旧 .md の主張同士は依然として矛盾しており、これは保全のため意図的にそのままにしている。新規の "PR<N>_COMPLETION_REPORT" 系 .md は作成しない方針を維持（本ファイルは唯一の例外。PR#258 以降の新規追加はゼロ）。 |

## 用語

- **解決**：本 PR の変更で当該問題が事実として解消された、または検証可能な形で正直化された。
- **正直化のみ実施**：根本実装は変えていないが、誤った主張・誇大な主張をコード/文書から削除し、
  事実に即した記述に置き換えた。
- **本 PR では対応しない**：当該問題には触れていない。次 PR で扱うことを明示。
- **完全には解決していない**：部分的に対処したが残課題があることを明示。

## 嘘をつかないために — 残課題の明示（累積版）

### 解決していない根本問題（量子回路としての backend 実行の欠如）

1. **本ノートブック (`tutorials/quantum_dynamics_gksl_comparison.ipynb`) の主要な
   時間発展（Cell 6 の `QuditGKSLSimulator.simulate`、Cell 8 の Qubit GKSL、ボソン
   付き Cell 12/14 等）は依然 NumPy/`scipy.linalg.expm` ベース**です。
   PR#258 で N=2 の **部品回路レベル**で MQT-Qudits backend (tnsim) 上で正しく実行
   できることを実測、PR#260 で N=2 の **forward-only full per-step (12 channels,
   14 qudits) 回路**を `interleaved_n2` 純 permutation レイアウト下で tnsim 実行成功
   (本セッションで再実測、‖Δρ_sys‖_F = 0.0)。しかし以下は依然未解決：
    - **組合せ Strang per-step (palindromic=True) を ancilla 再利用付き形態で
      backend 実行することは依然不可能** — mid-circuit ancilla reset が
      `tnsim`/`misim` で未実装のため。
    - **fresh-ancilla 版の palindromic=True**（24 ancillas + 2 system = 26 qudits、
      state vector 3^26 ≈ 2.5e12 amplitudes ≈ 76 TB）はあらゆる state-vector
      backend で root cause として infeasible。
    - **N>2 の fresh-ancilla forward-only 回路**は state vector が N=3:3^22≈500GB /
      N=4:3^30≈3PB であり、interleaved layout を実装しても state vector 自体が
      backend で持てない（本セッションで定量確認）。よって "N>2 用 interleaved
      layout が未実装" は **layout を実装しても backend 実行可能にはならない**ため、
      実質的な意味はない。
    - **実シナリオの NumPy `expm` 経路を量子回路に置換すること**には、上記の通り
      backend 側 mid-circuit ancilla reset の実装か、新たな量子回路化アルゴリズム
      （例：Kraus operator を直接 backend に渡す API、Lindblad-aware compiler）の
      導入が必要。本フォークの tutorials スコープ外。

2. **ゲート数の表示には現在 2 系統あります**：
    - `result["estimated_gates_per_step"]` / `result["n_high_level_gates_per_step"]`
      は依然「ハードコードした高レベルゲート数式」の値で、A-3 の旧表示と同じです。
      後方互換のため残してあります。`result["gate_count_method"] = "high_level_count"`
      によりこれが**コンパイラ計測値ではない**ことを明示しています。
    - `QuditGKSLSimulator.compute_compiler_measured_gate_counts(dt)` を呼べば
      MQT-Qudits の `compileO0`/`compileO1` を実際に通したネイティブゲート
      (`VirtRz`/`R`/`Rh`/`Rz`/`CEx`) の内訳と、`cu_multi` (TTA-pair Stinespring) の
      未分解個数（MQT-Qudits compiler が現状 `cu_multi` を 2-qudit gate に
      decompose できないため）を**そのまま**返します。本セッション再実測で
      `per_step_summary = {'native_gates': 10836, 'uncompiled_cu_multi': 4,
      'optimization_level': 0, 'backend': 'faketraps2trits'}`（N=2, dt=0.5）。
      `cu_multi` の compiler decompose 制約は MQT-Qudits 本体側の課題で、本フォーク
      tutorials のスコープ外。

3. **収束テーブル（NB Cell 28〜30）には `algorithm="stinespring"` (1 次, default)
   と `algorithm="exact_local_channels"` (2 次, PR#258 で追加) の両方の rate が
   並列出力**されます。Stinespring 側の `rate ≈ 1.0` は依然として事実
   （Stinespring dilation の数学的性質）であり、Cell 28 のその記述は維持。
   Exact local channel 側の `rate ≈ 2.0` も**実測**であり、嘘・ヒューリスティック
   なしで Strang+palindromic 構造が本来期待する 2 次収束を取り戻します。
   ただしこれは古典シミュレータ上での厳密 dissipator 適用であり、対応する
   量子回路化（Choi-Kraus からの Kraus operator 抽出 + ancilla reset 不要な
   実装等）は依然として未完了で、A-1 の残課題に含まれます。

4. **A-2 の独立性**：PR#259/#260 で `tutorials/qutip_gksl_reference.py`
   （`qutip.mesolve` ベース、in-tree GKSL スタックとコード共有なし）と
   `tutorials/run_qutip_cross_validation.py` を追加。これは **`exact_local_channels`
   の妥当性**を独立検証するもので、Stinespring パスや組合せ回路の正しさを直接
   検証するものではない。Qubit/Qudit シミュレータ間の比較は、依然として同一
   Stinespring/Trotter コアを共有している事実は変わらない（その旨は docstring と
   NB に明記済）。

### 本セッション (2026-05-13) で実施した再実測

PR 説明にある数値が現コードベース上で再現可能であることを以下の手順で個別に確認した：

| 検証項目 | 実行物 | 実測結果 | PR 主張 | 一致 |
|---|---|---|---|---|
| A-1 Part A (math equivalence) | `run_n2_full_step_fresh_ancilla_verification.py::part_a_math_equivalence` | 1.2s, ‖Δρ_sys‖_F = 0.0 | PR#259 | ✓ |
| A-1 Part B (subset probes) | 同 `part_b_backend_probe` | 11 分, 5/5 中 3 成功 + 2 MemoryError (468 GiB) | PR#259/#260 | ✓ |
| A-1 Part C (interleaved_n2 full per-step) | 同 `part_c_full_step_interleaved` | 1.3s, tnsim_ok=True, ‖Δρ_sys‖_F = 0.0, Tr=1.0, worst gate range = 6, intermediate ≈ 8 MB | PR#260 | ✓ |
| A-2 QuTiP cross-validation N=2 | `run_qutip_cross_validation.py` | `Tr|ρ_qutip − ρ_exact_local|/2 = 1.786e-08` | PR#260 (1.79e-08) | ✓ |
| A-2 QuTiP cross-validation N=4 | 同 | `7.683e-07` | PR#260 (7.68e-07) | ✓ |
| A-3 compiler-measured gates | `QuditGKSLSimulator.compute_compiler_measured_gate_counts(0.5, 0)` | `native_gates=10836, uncompiled_cu_multi=4` | PR#258 (10836, +cu_multi) | ✓ |
| B-1 `exact_local_channels` rate=2 | `pytest TestExactLocalChannelsConvergence` | 6/6 パス | PR#258 | ✓ |
| その他テスト | `pytest TestCompilerMeasuredGateCounts`, `TestN2TnsimVerification` | 各 1/1 パス | PR#258 | ✓ |

PR 説明と実環境の挙動はすべて一致した。PR#260 が「本ツリー統合後の verification
完走は時間切れで未確認」と記録した Part C は、本セッションで Part C 単体 1.3s で
完走、‖Δρ_sys‖_F = 0.0 を実測確認した。フルスクリプトの 11 分の所要時間は、Part B
の MemoryError ケースで NumPy が 468 GiB の確保を試みる時間が支配的であって、
Part C 自体の問題ではない。

### honestly に残る制約（実装で取り除くことのできない事項）

- **`palindromic=True` fresh-ancilla N=2**：state vector 3^26 ≈ 76 TB。state-vector
  backend の根本的容量制約。Density-matrix backend や Kraus-channel backend の
  導入が必要。
- **N≥3 の fresh-ancilla forward-only**：state vector が N=3 で 500 GB を越える。
  layout 工夫では救えない。
- **`cu_multi` (3-qudit Stinespring) の compiler decompose**：MQT-Qudits 本体の
  compiler が現状 `cu_multi` を 2-qudit gate に分解できない。本フォーク tutorials
  からは制御できない（compiler 側の改修が必要）。
- **6 シナリオ本体の量子回路化**：上記 A-1 セクションの通り、backend 側に
  mid-circuit ancilla reset または Kraus-aware API が無いと根本的に不可能。
  ヒューリスティックで「動いているふり」をすることは本ファイルの方針上禁止。
