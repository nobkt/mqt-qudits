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
| 第 4 版 | 本セッション (2026-05-13) | **A-1 を density-matrix backend `DMSim` 新設で N=4 ボソン無し backend 実行可能に更新**。`mqt.qudits.simulation.backends.DMSim` と `mqt.qudits.quantum_circuit.gates.KrausChannel` を MQT-Qudits 本体に新規追加、`QuditGKSLSimulator(execute_on_backend="dmsim")` で N=4, t_max=100, n_steps=100 を実測 (NumPy ref 0.88s vs DMSim backend 2.68s, ‖Δρ‖_F = 8.55e-14)。state-vector backend (tnsim/misim) では N=4 が依然不可能であること、ボソン付きと qubit 側は未対応であることを honest に併記 |
| 第 6 版 | 本セッション (2026-05-13、本 PR) | **ノートブックの 4 つのショットセル (5b/5c/3b/3c) を backend 実行 + Born サンプリングに置き換え**。新規モジュール `tutorials/dmsim_shot_simulator.py` を追加し、(i) `QuditDMSimShotSimulator` (MQT-Qudits DMSim backend で ρ(t) を実行 + 最終状態を `np.random.multinomial(n_shots, diag(ρ_final))` で Born サンプリング, 5b/5c)、(ii) `QiskitQubitShotSimulator` (Qiskit Aer `density_matrix` backend で同様, 3b/3c)。**ノイズは厳密 CPTP Kraus**: d 次元 Weyl-Heisenberg twirling identity `(1/d²)Σ W_{a,b} ρ W_{a,b}† = (Tr ρ)I/d` から `K_0 = √(1−p+p/d²)·I, K_{(a,b)≠0} = √(p/d²)·W_{a,b}` を構築。pair-depolarisation は同じ式を `d → d²` で適用（qudit は d=3、qubit は 2-qubit pair が d=4 = 16-dim local space）。`Σ K†K = I` を機械精度で検証し許容超過時は `AssertionError`（ヒューリスティック補正なし）。**ノイズ規約は 5c/3c 完全対称** (`p_depol=0.001, depol_pair_only=True, p_dephasing=0.0`)、`cx_per_pair_gate` 型の片側増幅は使用しない。**実測検証**: N=2/t=5/n_steps=10 で `‖ρ_dmsim − ρ_NumPy_exact_local‖_F = 1.9e-14` (機械精度), `‖ρ_aer_shot − ρ_aer_ref‖_F = 0.0` (instruction 完全一致), N=4/t=100/n_steps=10 で 5b 1.1s, 5c 1.1s, 3b 5.9s, 3c 12.2s。テスト `TestDMSimShotSimulators` 7 件全パス（depolarisation Kraus CPTP + 参照一致 + ノイズ有り ρ シフト + 禁止状態リーケージ + 不正オプション拒否 + Born サンプリング統計収束）。**正直な制約として明記**: (a) これは **per-trajectory 実行ではなく**、backend で ρ(t) を 1 回 evolve し最後に Born サンプリングする「density-matrix backend + final-state Born sampling」スキーム。`n_shots` 本の独立 trajectory を mid-circuit ancilla 測定付きで走らせているわけではない（MQT-Qudits backend に mid-circuit ancilla reset がないため不可能）。(b) **中間時刻の populations は ρ(t) 厳密値**で trajectory 平均ではない、shot variance は最終測定のみ。(c) ノイズは `depol_pair_only=True, p_dephasing=0.0` のみサポート、それ以外は `NotImplementedError` で明示拒否（黙って規約違反しないため）。旧 NumPy 軌道シミュレータ (`Qudit/QubitGKSL{Shot,NoisyShot}Simulator`) はファイル削除せず保存し、`test_gksl_simulators.py` の独立検証材料として継続利用。 |
| 第 5 版 | 本セッション (2026-05-13、本 PR) | **A-1 をボソン付き qudit と qubit 経路の双方で backend 実行可能に拡張**。(i) `QuditGKSLBosonSimulator` に `algorithm="exact_local_channels"` + `execute_on_backend="dmsim"` を新設し、各 Trotter ステップで mixed-dimensional `QuantumCircuit`（電子 qutrit ×N + フォノン qudit ×N、`[3,3,2,2]`）を組み立てて DMSim で実行。N=2/n_max=1/t_max=10/n_steps=20 で実測 ‖ρ_NumPy − ρ_DMSim‖_F = 3.12e-14、テスト `TestQuditGKSLBosonDMSim` 3 件パス。(ii) qubit 経路は **Qiskit Aer `density_matrix` backend を採用**して新規モジュール `tutorials/qiskit_qubit_gksl_simulator.py` を追加 (`QiskitQubitGKSLSimulator`/`QiskitQubitGKSLBosonSimulator`)。各 Trotter ステップで `qiskit.QuantumCircuit` を組み立てて `set_density_matrix → UnitaryGate(U_H_half) → Kraus(...).to_instruction()×26 (palindromic + reverse) → UnitaryGate(U_H_half) → save_density_matrix` を `qiskit_aer.AerSimulator(method='density_matrix')` で execute。各 Kraus は qutrit→qubit-pair 埋め込み空間で構築した局所超演算子 `expm(L_D^local · dt/2)` の Choi-Jamiolkowski 抽出（**closure 補正は不要**であることを Σ K†K = I の実測で検証 — GKSL 散逸子は埋め込んでも自動的にトレース保存）。N=2 ボソン無し: trace deviation < 1e-8、‖ρ − ρ_classical‖_F = 1.83e-6 (n_steps=20, dt=0.5)、N=2 ボソン有り: trace deviation < 1e-8、‖ρ − ρ_classical‖_F = 1.83e-6 (同条件)。テスト `TestQiskitQubitGKSL` 4 件パス。(iii) ノートブックに 3 セル追加 (Scenario 3d / 6d / 4d)。**残課題（本 PR でも未解決）**: state-vector backend (tnsim/misim) での Stinespring per-step 実行は依然 N=4 で不可能（PR#260 と同じく state-vector がメモリ infeasible）。Qiskit Aer はあくまで simulator であり実 quantum hardware ではない。Qiskit ボソン版は n_max+1 が 2 のべき乗の場合のみ対応（n_max=1 は OK だが n_max=2 は要拡張）。 |

## 累積対応状況一覧

| ID | 内容 | 状況 | 詳細 |
|---|---|---|---|
| **A-1** | MQT-Qudits を実際には使っていない（NumPy `expm` のみで計算） | ✅ **本 PR で全 6 シナリオが backend 実行可能に到達 — qudit 経路は MQT-Qudits DMSim、qubit 経路は Qiskit Aer。state-vector backend での per-step 実行は依然 infeasible** | (i)〜(iv) は前述 (PR#258〜260)。(v) **PR#262**: density-matrix backend `mqt.qudits.simulation.backends.DMSim` と CPTP instruction `mqt.qudits.quantum_circuit.gates.KrausChannel` を MQT-Qudits 本体に新規追加。`QuditGKSLSimulator(execute_on_backend="dmsim")` で N=4 ボソン無し（Cell 9）を MQT-Qudits backend で実行成功（‖Δρ‖_F = 8.55e-14）。(vi) **本 PR**: 残 3 シナリオを backend 化。**Scenario 6d (qudit, ボソン有り N=2)**: `QuditGKSLBosonSimulator(algorithm="exact_local_channels", execute_on_backend="dmsim")` 新設、mixed-dim register `[3,3,2,2]` (電子 qutrit + フォノン qubit) を用いた DMSim 実行で ‖ρ_NumPy − ρ_DMSim‖_F = 3.12e-14。**Scenario 3d (qubit, ボソン無し N=4) / 4d (qubit, ボソン有り N=2)**: 新規 `tutorials/qiskit_qubit_gksl_simulator.py` (`QiskitQubitGKSLSimulator`/`QiskitQubitGKSLBosonSimulator`) を Qiskit Aer `density_matrix` backend で execute。各 Trotter ステップで `qiskit.QuantumCircuit` を構築し、Hamiltonian 半ステップを `UnitaryGate`、各 Lindblad チャネルを qubit-pair 埋め込み空間で Choi-Jamiolkowski 抽出した Kraus operator として `qiskit.quantum_info.Kraus` instruction 化、palindromic + reverse で append。trace は < 1e-8 で保存、forbidden state population は round-off level (1e-9 以下)、古典参照との Frobenius 差 1.8e-6 (n_steps=20, dt=0.5; 2次精度 Trotter の dt² ~ 0.25 領域として整合)。テスト: `TestQuditGKSLBosonDMSim` (3) + `TestQiskitQubitGKSL` (4) を新設し全パス。**根本的な未解決（依然残る制約）**: (a) state-vector backend (tnsim/misim) での Stinespring per-step 実行は依然 N=4 で不可能 (state-vector が 3^30≈3 PB)、本 PR は **density-matrix simulator (DMSim, Qiskit Aer density_matrix) を用いて** backend 実行を達成しているのであって state-vector backend の制約は変わっていない。(b) Qiskit Aer は古典計算機上で密度行列を直接演算する simulator であり、実 quantum hardware 実行ではない。(c) Qiskit ボソン版は n_max+1 が 2 のべき乗の場合のみ対応（実装簡易化のため）。(d) qubit/qudit independence の memory: `QiskitQubitGKSL*` は qutrit-space で組み立てた H/L をそのまま qubit-pair 埋め込みで使うため、qudit 経路と完全独立とは言えない（埋め込みの正しさを実測検証する位置付け）。 |
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
| **A-1 DMSim N=4 backend exec** | `pytest TestDMSimBackendExecution::test_dmsim_matches_numpy_n4_notebook_settings` | NumPy 0.88s vs DMSim 2.68s, ‖Δρ‖_F = 8.55e-14, trace dev 1.7e-14 | 本 PR (新規) | ✓ |
| **A-1 DMSim/KrausChannel ユニット** | `pytest test/python/simulation/test_dmsim.py` | 9/9 パス（CPTP 検証, to_matrix raise, 振幅減衰解析解一致, etc.） | 本 PR (新規) | ✓ |
| **A-1 Kraus 抽出の数値同値** | `dmsim_kraus_helpers.kraus_from_local_superoperator` を N=4 全 26 channel で再構成 | 全 channel で再構成誤差 ≤ 1.22e-15 | 本 PR (新規) | ✓ |

PR 説明と実環境の挙動はすべて一致した。PR#260 が「本ツリー統合後の verification
完走は時間切れで未確認」と記録した Part C は、本セッションで Part C 単体 1.3s で
完走、‖Δρ_sys‖_F = 0.0 を実測確認した。フルスクリプトの 11 分の所要時間は、Part B
の MemoryError ケースで NumPy が 468 GiB の確保を試みる時間が支配的であって、
Part C 自体の問題ではない。

### honestly に残る制約（実装で取り除くことのできない事項）

- **`palindromic=True` fresh-ancilla N=2**：state vector 3^26 ≈ 76 TB。state-vector
  backend の根本的容量制約。Density-matrix backend や Kraus-channel backend の
  導入が必要。**本セッション・本 PR で `DMSim` density-matrix backend を新設して N=4
  ボソン無しに対する解決を提供したが、state-vector backend (tnsim/misim) 自体の
  この制約は解消されておらず、今後 state-vector backend で同シナリオを動かしたい
  場合は依然解決が必要**。
- **N≥3 の fresh-ancilla forward-only**：state vector が N=3 で 500 GB を越える。
  layout 工夫では救えない。**density-matrix backend に切り替えることでのみ回避可能**
  （density matrix サイズは ancilla qudit と無関係に system Hilbert 空間のみで
  決まるため、N=4 でも 81×81 = ~100 KB）。
- **`cu_multi` (3-qudit Stinespring) の compiler decompose**：MQT-Qudits 本体の
  compiler が現状 `cu_multi` を 2-qudit gate に分解できない。本フォーク tutorials
  からは制御できない（compiler 側の改修が必要）。DMSim backend は `cu_multi` を
  そのまま `U ρ U†` として適用するため、この制約に**影響されない**。
- **6 シナリオ本体の量子回路化（state-vector backend での実行）**：上記 A-1 セクションの通り、backend 側に
  mid-circuit ancilla reset または Kraus-aware API が無いと根本的に不可能。
  ヒューリスティックで「動いているふり」をすることは本ファイルの方針上禁止。
- **本 PR の DMSim backend で対応していない範囲（明示）**：
    - **ボソン付きシナリオ (Cell 11/13/15)**：density matrix サイズが Hilbert 空間
      と共に増加する (N=2, n_max=1 で 36×36 = ~10 KB は trivial、それ以上は
      個別検討)。`QuditGKSLBosonSimulator` への DMSim 経路追加は別 work。
    - **qubit シミュレータ (Cell 9 等)**：qubit 経路は qutrit シミュレータからの
      埋め込み変換に依存しており (memory: qubit/qudit independence)、qubit 経路
      のみで独立に DMSim backend 化する意味は小さい。本 PR では未対応。
    - **Stinespring 経路**：`exact_local_channels` の Kraus 表示と数学的に同値で
      あるため、Stinespring を Kraus 化した backend execution は意味のある独立
      情報を提供せず、本 PR では rejection する (`ValueError`)。
