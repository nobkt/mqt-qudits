# `tutorials/quantum_dynamics_gksl_comparison.ipynb` 詳細分析・改修報告（2026-07-13）

本書は、`tutorials/quantum_dynamics_gksl_comparison.ipynb`（全 53 セル、約 9.2 MB）の
詳細分析結果、実施した改修内容、および次ステップの計画をまとめたものである。
**事実ベースで記載しており、実施していないことを実施したとは書かない。**
各項目には検証方法（実測コマンド・実測値）を併記する。

---

## 1. 分析方法（実際に行ったこと）

1. ノートブックの全 53 セル（コード 27 / Markdown 26）の JSON を展開し、
   各セルのソース・出力・`execution_count` を静的に精査した。
2. 依存モジュール（`qudit_gksl_simulator.py`, `qubit_gksl_simulator.py`,
   `qiskit_qubit_gksl_simulator.py`, `dmsim_shot_simulator.py`,
   `qudit_gksl_circuit_simulator.py` 等）の API（クラス名・引数・戻り値のキー）と
   ノートブック側の呼び出しの整合性を grep ベースで照合した。
3. サンドボックス環境に `pip install -e .`（本フォークの MQT-Qudits）、
   `qiskit`, `qiskit-aer`, `ipython` をインストールし、
   **全 27 コードセルを縮小パラメータ（`n_steps=10`, `n_steps_quantum=10`、
   セル 11b は `n_steps_list=[5,10,20]`）で順次実行するスモークテスト**を行った。
4. スモークテストで検出した唯一の実行時エラー（後述の `pylatexenc` 欠如）を解消後、
   該当セル（追補: 回路可視化セル）の再実行で正常終了することを確認した。
5. Markdown 修正後、`n_steps=100` 等の**本来のパラメータでノートブック全体を
   `jupyter nbconvert --execute --inplace` によりフル実行**した
   （結果は本書 §4 に記載）。

---

## 2. 検出した問題点

### 2.1 実行状態の不整合（未実行セルと旧出力の混在）【重大・修正済】

改修前のノートブックは、以下の 10 個のコードセルが **`execution_count: null`・出力なし**
（＝追加後一度も実行されていない）のまま保存されていた:

- セル 7（A-3: コンパイラ実測ゲート数）
- セル 9（シナリオ 5d: DMSim backend）
- セル 13（シナリオ 3d: Qiskit Aer）
- セル 19（シナリオ 6d: DMSim mixed-dim）
- セル 23（シナリオ 4d: Qiskit Aer boson）
- セル 37（11b: 収束特性分析）
- セル 41/43/45/47（シナリオ 5b/5c/3b/3c: ショットベース）

一方、残りのセルには旧実行時の出力が残っており、`execution_count` が
1〜13, 19, 20 と**欠番だらけの不整合状態**であった。特にセル 49・51
（ショットベース比較プロット）は `result5b` / `result3c` 等の**未実行セルの変数に
依存する出力を保持**しており、「掲載されている出力がどのコードから生成されたのか
検証不能」という、本リポジトリの正直性ポリシー（`STATUS_HONEST_2026-05.md`）に
反する状態だった。

**対処**: ノートブック全体を本来のパラメータでフル実行し、全セルの出力と
`execution_count` を一貫した状態に更新した（§4）。

### 2.2 まとめセル（セル 52）の記述が旧実装のままで事実と不一致【重大・修正済】

セル 12b（Markdown、セル 39）には「ショットベース 4 セルは
**density-matrix backend 実行 + 最終状態 Born サンプリング**に置き換えた。
per-trajectory ではない」と正直に明記されているにもかかわらず、
まとめセル（セル 52）は**置き換え前の旧実装の記述のまま**であった:

| 箇所 | 旧記述（誤り） | 事実（修正後） |
|---|---|---|
| 3b/5b の方式 | 「量子軌道法、アンシラ測定」 | backend ρ 実行 + 最終 Born サンプリング（per-trajectory ではない） |
| ショットベースの説明 | 「純粋状態 \|ψ⟩ を発展、アンシラを Born 則で確率的に測定、N_shots 回繰り返し」 | ρ(t) を backend で 1 回発展、最終時刻のみ multinomial Born サンプリング。中間時刻は ρ(t) 厳密値 |
| ノイズモデル | 「確率 p_depol で Weyl-Heisenberg 演算子を**確率的適用**」 | **厳密 CPTP Kraus チャネル**として決定論的に ρ に作用（Σ K†K = I を機械精度検証） |
| シナリオ 4 の量子資源 | 「42 qubits」 | **18 qubits**（4 electronic + 2 phonon + 12 ancilla; セル 7b の実行出力と一致） |
| シナリオ 6 の量子資源 | 「8 qutrits + 26 ancilla」 | **16 qudits**（2 qutrits + 2 phonon qudits + 12 ancilla; セル 6 の実行出力と一致） |
| 3b/3c/5b/5c の量子資源 | 「34 qubits / 4 qutrits + 26 ancilla」 | backend 回路はアンシラなし（Kraus instruction / KrausChannel 使用; `dmsim_shot_simulator.py` の `n_ancilla: 0` と回路構成で確認） |

このほか「密度行列（Statevector）」という自己矛盾した用語
（密度行列シミュレーションを Statevector と呼んでいた）を
「密度行列（NumPy 直接演算）」に統一した。

**対処**: セル 52 を上表の事実に基づいて全面的に書き換えた。あわせて
収束特性の節に `exact_local_channels` の $O(\Delta t^2)$ 収束
（セル 11b で実測される内容）を追記した。

### 2.3 シナリオ一覧（セル 0）の欠落【中・修正済】

冒頭のシナリオ一覧表に、ノートブック本体に存在するシナリオ
**3d / 4d / 5d / 6d（backend 実行系 4 セル）が記載されていなかった**。
4 行を追加し、全 14 シナリオを網羅する表に更新した。

### 2.4 実行要件の未記載（`pylatexenc`）【小・修正済】

「追補: 1トロッターステップ量子回路の可視化」セル（セル 25）は
`qiskit.visualization.circuit_drawer(output="mpl")` を使用しており、
`pylatexenc` がないと `MissingOptionalLibraryError` で失敗する
（スモークテストで実測）。セル 0 に「実行要件」節を追加し、
`pylatexenc` / `qiskit` / `qiskit-aer` の必要性を明記した。

### 2.5 問題なしと確認した点（誤検出の防止のため記録）

- **コードセルの API 整合性**: 旧クラス名 `QubitGKSLCircuitSimulator` /
  `QuditGKSLCircuitSimulator` / `QuditGKSLCircuitBosonSimulator`（セル 25 で使用）は
  D-1 改名後も後方互換エイリアスとして各モジュール末尾に残置されており
  （`qudit_gksl_circuit_simulator.py:1636` 等）、動作する。
- **全 27 コードセルが縮小パラメータで例外なく完走**（`pylatexenc` 導入後）。
  各セル内の `assert`（DMSim vs NumPy の ‖Δρ‖_F < 1e-10 等）もすべて通過。
- セル 2 の `os.path.dirname(os.path.abspath('__file__'))` は文字列リテラル
  `'__file__'` を用いた慣用句であり、カーネルの CWD（tutorials/）を返すため動作上の
  問題はない（美観の問題のみ。挙動を変えないため今回は修正しない）。
- セル 9b（ユニタリ TTA）の `H_TTA` はエルミート性を `assert` 済みで、
  GKSL 側の Lindblad TTA 演算子と同一の遷移構造を持つ。
  「物理的に等価な比較ではない」ことはセル内コメントに明記済み。
- セル 12b の正直な記述（per-trajectory ではない等）は
  `dmsim_shot_simulator.py` の実装・`STATUS_HONEST_2026-05.md` 第 6 版と一致。

---

## 3. 実施した改修の一覧

1. **セル 0（Markdown）**: シナリオ一覧に 3d/4d/5d/6d を追加、
   「密度行列（Statevector）」→「密度行列（NumPy 直接演算）」に用語修正、
   「実行要件」節（`pylatexenc` 等）を新設。
2. **セル 52（Markdown, まとめ）**: §2.2 の表のとおり、旧実装の記述・誤った量子資源数を
   現実装の事実に合わせて全面修正。backend 実行系（3d/4d/5d/6d）の行を追加。
   ノイズモデル・収束特性・設計原則の各節を現状に合わせて更新。
3. **ノートブック全体のフル実行**（本来のパラメータ: `t_max=100`, `n_steps=100`,
   `n_steps_quantum=100`, `n_shots=1000`）により、全セルの出力を一貫した状態に再生成。

コード（.py モジュール側）は今回の分析で不具合が見つからなかったため変更していない。

---

## 4. フル実行の結果

フル実行は `jupyter nbconvert --to notebook --execute --inplace` で実施した。
実行結果の要点（トレース保存・参照一致等の数値）は、実行後のノートブック本体の
各セル出力を一次情報とすること。本書には実行の成否のみを記録する:

- 実行環境: Python 3.12 / numpy 2.5.1 / scipy 1.18.0 / qiskit 2.5.0 /
  qiskit-aer 0.17.2 / 本フォークの mqt.qudits（editable install）
- 実行結果: 全 53 セルがエラーなく完走（`execution_count` 1〜27 の連番、
  セル内 `assert`（‖ρ_ref − ρ_dmsim‖_F < 1e-10 等）もすべて通過）

---

## 5. 次ステップの計画

### 5.1 短期（次の PR で対応可能）

1. **`QuditDMSimShotSimulator` の結果 dict の `n_ancilla_qudits` 表示の見直し**:
   現状、内部で保持する `QuditGKSLSimulator` の属性（Stinespring 換算の 26）を
   そのまま返しているが、DMSim 実行回路にはアンシラは存在しない
   （KrausChannel を使用）。誤解を招くため、`n_ancilla_qudits_backend_circuit: 0`
   のような backend 回路実態のキーを併記するか、キー名を明確化する。
2. **CI へのノートブックスモークテスト追加**: 今回の「未実行セルと旧出力の混在」は、
   セル追加時にフル実行せずコミットしたことが原因。縮小パラメータ
   （`n_steps=10` 程度）で全コードセルを実行する pytest
   （`test_pr89_notebook_fix.py` と同系統）を追加し、セル追加時の実行漏れを機械的に
   検出する。
3. **`pylatexenc` を含む依存の宣言**: `pyproject.toml` の optional-dependencies
   （例: `[project.optional-dependencies] tutorials = [...]`）として
   `qiskit`, `qiskit-aer`, `pylatexenc`, `ipython` を宣言する。

### 5.2 中期（既知の残課題、`STATUS_HONEST_2026-05.md` と共通）

4. **A-1 の残課題**: state-vector backend（tnsim/misim）での Stinespring per-step
   実行は N=4 で依然不可能（3^30 amplitudes ≈ 3 PB）。mid-circuit ancilla reset の
   backend 実装、または Lindblad-aware compiler の導入が必要。
5. **`cu_multi`（TTA-pair Stinespring, 27×27）のコンパイラ分解**: 現状
   MQT-Qudits compiler が decompose できず、セル A-3 は未分解個数を正直に報告する
   のみ。分解アルゴリズムの実装は MQT-Qudits 本体側の課題。
6. **Qiskit ボソン版の `n_max` 制約**: `n_max+1` が 2 のべき乗の場合のみ対応
   （`n_max=1` は可、`n_max=2` は要拡張）。
7. **収束テーブル（セル 11b 直後の Markdown）の自動生成化**: 現状は PR#258 時点の
   実測値を手書きした静的表であり、コード変更時に古くなるリスクがある。
   セル 11b の実測結果から表を生成する形に変更する。

### 5.3 実施しないと判断した事項（理由付き）

- **セル 2 の `abspath('__file__')` イディオムの書き換え**: 挙動に問題がなく、
  変更は差分を増やすだけのため今回は見送り。
- **旧 NumPy 軌道シミュレータ（`qudit/qubit_gksl_shot_simulator.py` 等）の削除**:
  `test_gksl_simulators.py` の独立検証材料として意図的に残置されている
  （セル 12b に明記）ため、削除しない。

---

## 6. 検証コマンド（再現手順）

```bash
# 依存インストール
pip install -e .
pip install qiskit qiskit-aer ipython jupyter pylatexenc

# フル実行（tutorials/ ディレクトリで）
cd tutorials
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=3600 quantum_dynamics_gksl_comparison.ipynb
```
