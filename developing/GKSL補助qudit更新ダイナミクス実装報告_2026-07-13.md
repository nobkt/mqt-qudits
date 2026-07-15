# TTA-UC GKSL-Lindblad 量子ダイナミクス: 補助qudit（アンシラ）を更新しながらの実行 —— 実装報告と次ステップ計画（2026-07-13）

本書は、`tutorials/quantum_dynamics_gksl_comparison.ipynb` の GKSL-Lindblad 量子ダイナミクスを
「**Lindblad 項の処理に補助量子ビット（アンシラ）を更新しながら量子ダイナミクスを実行する本来の処理**」
となるように改修した作業の報告である。
**事実ベースで記載し、実施していないことを実施したとは書かない。**
各項目には検証方法（実測コマンド・実測値）を併記する。

---

## 1. 改修前の状態（問題点の確認）

指摘のとおり、改修前のノートブックの backend 実行系 GKSL ダイナミクスは
**アンシラを一切使わずに** Lindblad 項を処理していた:

| 改修前のパス | Lindblad 項の処理 | アンシラ |
|---|---|---|
| シナリオ 5（NumPy 参照） | Stinespring 膨張 → NumPy で部分トレース | 数学的操作としてのみ存在（回路・backend 実行なし） |
| シナリオ 5d（DMSim backend） | 厳密 Kraus（`KrausChannel`、Choi 分解） | **なし**（回路はシステム 4 qutrit のみ） |
| シナリオ 5b/5c/3b/3c（ショット系） | backend ρ 実行（Kraus） + 最終 Born サンプリング | **なし**（per-trajectory でもない） |

さらに `QuditGKSLSimulator` は `algorithm='stinespring'` と `execute_on_backend='dmsim'` の
組合せを **`ValueError` で明示拒否** していた（改修前 `qudit_gksl_simulator.py` の
「Stinespring dilation has no Kraus representation on the system alone」）。
つまり「アンシラを回路に持ち、更新しながら backend で時間発展する」パスは存在しなかった。

確認方法: 改修前コードの `grep -n "exact_local_channels" tutorials/qudit_gksl_simulator.py`、
ノートブックのセル 9（シナリオ 5d）のソース精査。

## 2. 参照した実装

### 2.1 tensorcircuit-ng（https://github.com/tensorcircuit/tensorcircuit-ng）

tensorcircuit-ng における一般量子チャネルの扱いを参照した:

- `Circuit.general_kraus` / `unitary_kraus`（`tensorcircuit/channels.py`, `basecircuit.py`）:
  純粋状態シミュレータ上で Kraus チャネル $\{K_k\}$ を適用する際、
  **各 Kraus 分岐の Born 確率 $p_k = \langle\psi|K_k^\dagger K_k|\psi\rangle$ で結果 $k$ を確率的に選択し、
  $|\psi\rangle \to K_k|\psi\rangle/\sqrt{p_k}$ と射影・規格化する**。
  これは「膨張アンシラを付加 → 膨張ユニタリ → アンシラを測定して破棄（=更新）」の
  Stinespring 描像と数学的に同一（モンテカルロ波動関数法 / 量子軌道法）。
- `DMCircuit`（密度行列シミュレータ）: 同じチャネルを
  $\rho \to \sum_k K_k \rho K_k^\dagger$ として決定論的に適用する。
  これは「アンシラを付加 → 膨張ユニタリ → アンシラを部分トレース（=更新して再初期化）」と同一。

### 2.2 MQT-Qudits 側（本フォーク）

- `src/mqt/qudits/simulation/backends/dmsim.py`: 密度行列 backend。
  `CustomTwo`/`CustomMulti` ユニタリと `KrausChannel`（CPTP 検証付き）を実行できる。
- `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`: CPTP 条件
  $\sum_k K_k^\dagger K_k = I$ を構築時に機械精度で検証する Kraus 命令。
- MQT-Qudits には **専用の reset / mid-circuit 測定命令は存在しない**（`grep -rn "reset" src/mqt/qudits` で確認）。
  ただし「アンシラを $|0\rangle$ に初期化する」操作は厳密 CPTP チャネル
  $\{K_k = |0\rangle\langle k|\}_{k=0..d-1}$（$\sum_k K_k^\dagger K_k = I$）として表現でき、
  既存の `KrausChannel` で **mid-circuit のアンシラ更新を厳密に実行できる**。
  本改修はこの事実を利用した（MQT-Qudits 本体のコード変更は不要だった）。

## 3. 実施した改修（実際に行ったこと）

### 3.1 `tutorials/qudit_gksl_simulator.py`: アンシラ更新付き Stinespring backend 実行モード

`QuditGKSLSimulator(algorithm='stinespring', execute_on_backend='dmsim')` を新規に有効化した。

- **`_prepare_stinespring_ancilla_circuit(dt)`（新設）**: 1 トロッターステップの回路を
  **N=4 システム qutrit + 1 物理アンシラ qudit（d=3, 位置 N）** の
  MQT-Qudits `QuantumCircuit` として一度だけ構築する:
  1. `cu_multi`: $\exp(-iH_{total}\,dt/2)$（システムのみ）
  2. 各 Lindblad チャネル（palindromic 順、26 チャネル × 2 パス = 52 回）ごとに
     - **局所 Stinespring 膨張ユニタリ** $U_\alpha = \exp(-i\sqrt{dt/2}\,G_\alpha)$
       （単一サイト: 9×9 を `cu_two([anc, site])`、TTA ペア: 27×27 を
       `cu_multi([anc, site_i, site_j])`。アンシラ脚が先頭 = env⊗sys ブロック構造）
     - **mid-circuit アンシラ更新（reset）**: `KrausChannel(anc, {K_k=|0><k|})`。
       アンシラを測定して結果を破棄し $|0\rangle$ に初期化する操作の密度行列表現であり、
       tensorcircuit-ng の `general_kraus` の measure-and-discard と同じ意味論。
  3. `cu_multi`: 閉じの $\exp(-iH_{total}\,dt/2)$
- **`_trotter_step_dmsim_stinespring_ancilla(rho)`（新設）**: 各ステップで
  $\rho_{tot} = \rho_{sys} \otimes |0\rangle\langle 0|_{anc}$ を `initial_density_matrix` として
  DMSim backend で上記回路を実行し、アンシラの部分トレースで $\rho_{sys}$ を回収する。
- 局所 Lindblad 演算子（3×3 / 9×9）とチャネル順序は
  `qudit_gksl_circuit_simulator.QuditGKSLKrausSimulator.lindblad_local_info` を再利用
  （`gksl_math_utils.build_lindblad_operators` と同一順序であることをコード精査で確認、
  チャネル数一致は実行時 `ValueError` チェックで防御）。
- 結果 dict に `execute_on_backend` と `n_ancilla_qudits_backend_circuit`
  （stinespring+dmsim → 1、exact_local_channels+dmsim → 0、backend なし → None）を追加し、
  backend 回路の実態を正直に報告するようにした
  （既報告書 §5.1-1 の「`n_ancilla_qudits` 表示の見直し」にも部分対応）。

**数学的正当性**: 「アンシラ $|0\rangle$ 付加 → 膨張ユニタリ → 部分トレース」（NumPy 参照実装
`stinespring_utils.apply_stinespring_to_density_matrix`）と
「アンシラ $|0\rangle$ → 膨張ユニタリ → reset チャネル $\{|0\rangle\langle k|\}$（トレースアウト + 再初期化）」は
チャネルとして厳密に同一。局所ユニタリの埋め込みは、生成子 $G_\alpha$ が
（アンシラ, 対象サイト）にのみ台を持つため大域 `expm` と厳密一致する。

### 3.2 `tutorials/test_gksl_simulators.py`: 検証テスト

`TestDMSimBackendExecution` を更新:

- 旧テスト `test_dmsim_requires_exact_local_channels_algorithm`
  （stinespring+dmsim が拒否されることを確認するテスト）を削除し、以下に置換:
- `test_dmsim_stinespring_ancilla_matches_numpy_reference`:
  N=2, t=10, 20 ステップで NumPy Stinespring 参照と `‖Δρ‖_F < 1e-12` の一致、
  トレース保存、`n_ancilla_qudits_backend_circuit == 1` を検証。
- `test_dmsim_stinespring_circuit_contains_ancilla_and_resets`:
  backend 回路が本当に N+1 qudit で構成され、チャネル適用回数（52）と同数の
  reset `KrausChannel` がアンシラのみに作用し、その Kraus 集合が厳密に
  $\{|0\rangle\langle k|\}$ であることを instruction レベルで検証。

実測: `python -m pytest test_gksl_simulators.py -k "DMSimBackendExecution" -q`
→ **5 passed**（既存 2 テストも回帰なし）。
`-k "QuditGKSLSimulator or Stinespring"` → **14 passed**（回帰なし）。

### 3.3 ノートブック `quantum_dynamics_gksl_comparison.ipynb`: 新規 4 セル + Markdown 更新

- **シナリオ 5e（新規、セル 10-11）**: 「補助qudit を回路内で更新しながらの
  Stinespring backend 実行」。セル内で以下を実測検証:
  - 回路実態: qudit 構成 `[3,3,3,3,3]`、Stinespring ユニタリ 52 個/ステップ、
    アンシラ reset 52 個/ステップ、reset がアンシラのみに作用（`assert`）
  - NumPy Stinespring 参照（シナリオ 5 の `result5`）との一致
    `‖Δρ‖_F < 1e-10`（`assert`; 実測は round-off レベル ~1e-16）
  - 厳密解（`res5d_ref`）との差 = Stinespring 1 次近似の $O(\Delta t)$ 誤差を**隠さず表示**
- **シナリオ 5f（新規、セル 12-13）**: tensorcircuit-ng `general_kraus` と同じ確率的意味論の
  量子軌道法（既存 `QuditGKSLShotSimulator` を使用）。チャネル毎に
  アンシラ付加 → Stinespring ユニタリ → **アンシラ Born 測定 → 射影・規格化 → 破棄（更新）**。
  軌道平均と密度行列解の差が統計誤差 $O(1/\sqrt{n_{shots}})$ の範囲内であることを `assert`。
  **NumPy 純粋状態演算であり backend 実行ではない**ことをセル・Markdown 両方に明記
  （状態ベクトル backend は mid-circuit 測定未実装のため; §5 残課題）。
- セル 0 のシナリオ一覧に 5e / 5f の行を追加。
- まとめセル（旧セル 52）の表に 5e / 5f を追加し、
  「アンシラ更新付き Stinespring 実行」の節を新設。
- ノートブック全体を本来のパラメータ（`t_max=100`, `n_steps=100`,
  `n_steps_quantum=100`, 5f は `n_shots=200`）で
  `jupyter nbconvert --to notebook --execute --inplace` によりフル実行し、
  全セルの出力・`execution_count` を一貫した状態にした（実行結果は §4）。

### 3.4 変更しなかったもの（理由付き）

- **MQT-Qudits 本体（`src/`）**: 変更なし。アンシラ reset は既存 `KrausChannel` で
  厳密に表現できるため、新規命令の追加は不要だった（§2.2）。
  専用 `reset` 命令や状態ベクトル backend の mid-circuit 測定は次ステップ（§5）。
- **シナリオ 5d（KrausChannel 方式）**: 削除せず併存。アンシラ無しの厳密 Kraus 実行は
  2 次精度の参照値として引き続き有用であり、5e との差（Stinespring 1 次近似誤差）を
  実測表示する比較対象として使用。
- **tensorcircuit-ng の依存追加**: 行わない。参照したのはアルゴリズムの意味論
  （`general_kraus` の測定・破棄、`DMCircuit` のチャネル適用）であり、
  qudit（d=3）ネイティブな本問題には既存の MQT-Qudits/NumPy 実装で
  同じ意味論を厳密に実装できるため。

## 4. 実測検証結果

実行環境: Python 3.12 / numpy / scipy / 本フォークの mqt.qudits（editable install）。

1. **プロトタイプ検証**（`/tmp/proto/proto_ancilla.py`、N=4, 5 ステップ）:
   アンシラ更新回路（DMSim）と NumPy Stinespring 参照の差
   `‖Δρ‖_F = 7.5e-17`（round-off）。
2. **単体テスト**: `pytest test_gksl_simulators.py -k DMSimBackendExecution` → 5 passed。
3. **N=4, 10 ステップのセル・スモークテスト**: `‖ρ_5e − ρ_5‖_F = 1.8e-16`、
   trace deviation `1.2e-14`、厳密解との差 `8.9e-3`（$O(\Delta t)$、$dt=10$）。
4. **軌道法（5f）**: n_shots=30, 10 ステップで
   `‖ρ̄_軌道 − ρ_密度行列‖_F = 1.6e-1 < 5/√30`（統計誤差の範囲内）。
5. **ノートブックフル実行**: `jupyter nbconvert --execute --inplace`（本来のパラメータ）で
   全 57 セル（コード 27 セル、`execution_count` 1〜27 連番）がエラーなく完走
   （セル内 `assert` すべて通過）。主要実測値:
   - シナリオ 5e: `‖ρ_5e − ρ_5(NumPy Stinespring 参照)‖_F = 2.63e-16`（round-off 一致）、
     trace deviation `0.0`、厳密解との差 `8.78e-4`（$O(\Delta t)$、$dt=1$）、
     backend 実行 12.1 s / 100 ステップ
   - シナリオ 5f: n_shots=200 で `‖ρ̄_軌道 − ρ_密度行列‖_F = 5.34e-2 < 1/√200 = 7.07e-2`
     （統計誤差の範囲内）、60.0 s
   その他の数値の一次情報は実行後のノートブック本体の各セル出力を参照のこと。

## 5. 次ステップの計画

### 5.1 短期

1. **状態ベクトル backend（tnsim/misim）への mid-circuit 測定・reset 命令の実装**:
   現状、per-trajectory のアンシラ測定・更新（シナリオ 5f）は NumPy 実装であり、
   MQT-Qudits backend 上では実行できない。tensorcircuit-ng の
   `cond_measure`/`general_kraus` に相当する「計算基底測定 → 分岐選択 → 規格化 →
   （必要なら）$|0\rangle$ 再準備」命令を `src/mqt/qudits` の命令セットと
   tnsim/misim に追加すれば、5f を本物の backend 実行にできる。
   これは MQT-Qudits 本体への API 追加を伴うため独立 PR とすべき。
2. **専用 `Reset` 命令の追加**: 現在 reset は `KrausChannel({|0><k|})` として
   発行しており動作は厳密だが、回路の可読性・コンパイラ対応の観点では
   専用命令（QASM の `reset` 相当）が望ましい。DMSim では Kraus と同一実装で済む。
3. **qubit 側（シナリオ 3 系）への同等機能の展開**: 本改修は qudit 側のみ。
   Qiskit Aer は `reset` 命令をネイティブにサポートするため、
   `qiskit_qubit_gksl_simulator.py` に「アンシラ qubit + mid-circuit reset」版を
   追加するのは比較的容易（Aer の `density_matrix` method は reset を扱える）。
4. **5f の n_shots 収束スタディ**: n_shots ∈ {50, 200, 800, 3200} で
   $O(1/\sqrt{n_{shots}})$ 収束を系統的に実測するセルの追加。

### 5.2 中期（既知の残課題、`STATUS_HONEST_2026-05.md` と共通）

5. **A-1 の残課題**: N=4 の per-step Stinespring を状態ベクトル backend で
   実行するには、上記 1 の mid-circuit reset 実装が前提
  （fresh-ancilla 方式は 3^30 amplitudes で不可能。reset 再利用方式なら 3^5 で済む —
   本改修の DMSim 版と同じ回路構造が使える）。
6. **アンシラ更新回路のコンパイラ分解**: 5e の回路には 27×27 の `cu_multi`
   （TTA ペア Stinespring）が含まれ、MQT-Qudits compiler は現状これを分解できない
   （既知課題、A-3）。ネイティブゲート実行への道はコンパイラ側の拡張が必要。

## 6. 検証コマンド（再現手順）

```bash
pip install -e .
pip install qiskit qiskit-aer ipython jupyter pylatexenc pytest

# 単体テスト
cd tutorials
python -m pytest test_gksl_simulators.py -k "DMSimBackendExecution" -q

# ノートブックフル実行
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=3600 quantum_dynamics_gksl_comparison.ipynb
```
