# GKSL補助qudit更新ダイナミクス 次ステップ実装報告（2026-07-13）

本報告は `developing/GKSL補助qudit更新ダイナミクス実装報告_2026-07-13.md` の「5. 次ステップの計画」（§5.1 短期項目）の実装結果を、事実ベースでまとめたものである。誇張・虚偽は一切含まない。未完了項目は未完了と明記する。

## 1. 実装完了項目

### 1.1 §5.1-2: 専用 `Reset` 命令（完了）

- 新規ファイル `src/mqt/qudits/quantum_circuit/gates/reset.py`
  - `Reset(KrausChannel)` サブクラス。Krausセット {K_k = |0><k|, k=0..d-1}、`qasm_tag="reset"`。
- `src/mqt/qudits/quantum_circuit/gates/__init__.py` に登録。
- `QuantumCircuit.reset_qudit(qudit)` ビルダーを `circuit.py` に追加（既存の `reset()` は回路全体クリアのため別名にした）。
- `Reset` は `KrausChannel` のサブクラスであるため、DMSim は無変更で動作する（`isinstance` 判定を利用）。
- テスト: `test/python/simulation/test_reset.py` を新規作成、**7件全て合格**（実測）。

### 1.2 §5.1-1: TNSim での mid-circuit Reset/Kraus チャネル対応（完了）

- `src/mqt/qudits/simulation/backends/tnsim.py` を改修:
  - 回路を `KrausChannel` 命令で区切り、ユニタリ区間は tensornetwork で縮約、チャネルは Born 確率サンプリング＋射影＋再正規化で確率的に適用（1回の `execute` = 1トラジェクトリ）。
  - `initial_state` オプション（任意の初期状態ベクトルから開始）と `seed` オプション（`backendv2.py` の `DefaultOptions` に追加）を実装。
  - Kraus確率の合計が 1±1e-8 を外れた場合はごまかさず `ValueError` を送出する。
- 既存の tnsim/dmsim/circuit テスト **20件合格**（回帰なし、実測）。
- **misim は非対応のまま**（C++カーネル `state_vector_simulation` の改修が必要なため。`QuditGKSLSimulator` は misim 指定時に明示的に `ValueError` を出す）。

### 1.3 §5.1-1 継続: `QuditGKSLSimulator` の tnsim トラジェクトリモード（完了）

- `tutorials/qudit_gksl_simulator.py` に `execute_on_backend='tnsim'`（`algorithm='stinespring'` 限定）、`n_trajectories`、`seed` を追加。
- `_simulate_tnsim_trajectories`: トラジェクトリごとに純粋状態を TNSim 上でステップ実行し、補助quditのリーク（>1e-10 で `RuntimeError`）を検証したうえで |ψ⟩⟨ψ| を平均。
- 検証（実測）: N=2、n_trajectories=100 で ‖ρ_traj − ρ_dmsim‖_F = **0.353** < 5/√100 = 0.5（統計誤差の理論オーダー内）。
- テスト: `tutorials/test_gksl_simulators.py` に `TestTNSimTrajectoryBackendExecution`（4件）を追加。既存 DMSim テストと合わせ **9件合格**（実測）。

### 1.4 §5.1-3: Qiskit Aer 量子ビット側 ancilla+reset シミュレータ（完了）

- `tutorials/qiskit_qubit_gksl_simulator.py` に `QiskitQubitGKSLStinespringSimulator` を追加。
  - 2N 系qubit + 補助qubit 1個。チャネルごとに局所 Stinespring 拡張ユニタリ（8×8/32×32、Qiskit little-endian 規約で ancilla を最後に配置）＋ネイティブ `qc.reset(anc)`。Aer density_matrix 法で実行し、各ステップで ancilla を部分トレース。
- 検証（実測）: N=2 で NumPy 参照実装 `QubitGKSLSimulator` と最大差 **3.7e-15** で一致。トレース≈1、禁制ポピュレーション 0。
- テスト: `TestQiskitQubitGKSL` に3件追加、**10件合格**（実測）。

## 2. 未完了項目（正直な記載）

### 2.1 §5.1-4: n_shots 収束スタディ（スクリプト作成済み・実行未完了）

- `tutorials/run_5f_shots_convergence_study.py` を作成した（N=4, t_max=100, n_steps=100, n_shots ∈ {50, 200, 800, 3200}、NumPy Stinespring 決定論参照との ‖Δρ‖_F を log-log 傾き（理論値 −0.5）でフィット、tnsim トラジェクトリ（50本）とのクロスチェック付き）。
- **実行はセッション時間内に完了しなかった**。約30分以上走らせたが出力ゼロの段階で時間切れとなり停止した。**収束表・傾きの実測値は得られていない**。この点を偽らずに明記する。
- スクリプト自体は成果物としてコミットしており、後述の手順で再実行可能。

### 2.2 その他の未実施事項

- misim（C++）への reset 対応（§5.1-1 の一部）: 未実施。C++ 側の改修が必要。
- ノートブック `tutorials/quantum_dynamics_gksl_comparison.ipynb` のシナリオ 5e/5f セルの新機能対応・再実行: 未実施。
- 元報告 §5.2 の中期項目（A-1, A-3 等）: 未着手。

## 3. 検証コマンド（実際に実行し合格を確認したもの）

```bash
pip install -e .        # 加えて qiskit, qiskit-aer, pytest が必要
python -m pytest test/python/simulation/test_reset.py -v                     # 7 passed
python -m pytest test/python/simulation/test_tnsim.py test/python/simulation/test_dmsim.py test/python/simulation/test_circuit.py  # 20 passed
python -m pytest tutorials/test_gksl_simulators.py -k "DMSimBackendExecution or TNSimTrajectory" -v  # 9 passed
python -m pytest tutorials/test_gksl_simulators.py -k "QiskitQubitGKSL" -v   # 10 passed
```

収束スタディの再実行（長時間: 目安30分〜1時間以上）:

```bash
python tutorials/run_5f_shots_convergence_study.py
```

## 4. 変更ファイル一覧

- 新規: `src/mqt/qudits/quantum_circuit/gates/reset.py`, `test/python/simulation/test_reset.py`, `tutorials/run_5f_shots_convergence_study.py`
- 変更: `src/mqt/qudits/quantum_circuit/gates/__init__.py`, `src/mqt/qudits/quantum_circuit/circuit.py`, `src/mqt/qudits/simulation/backends/backendv2.py`, `src/mqt/qudits/simulation/backends/tnsim.py`, `tutorials/qudit_gksl_simulator.py`, `tutorials/qiskit_qubit_gksl_simulator.py`, `tutorials/test_gksl_simulators.py`

## 5. 次ステップの計画

### 5.1 短期

1. **収束スタディの実行完了**: `run_5f_shots_convergence_study.py` を計算資源のある環境で完走させ、n_shots 収束表と log-log 傾き（理論 −0.5 との比較）を本報告に追記する。必要なら n_steps を減らした縮小版を先に走らせ、パラメータ変更は正直に記録する。
2. **ノートブック更新**: `quantum_dynamics_gksl_comparison.ipynb` のシナリオ 5e/5f を `execute_on_backend='tnsim'` および `Reset` 命令ベースに書き換え、`jupyter nbconvert --execute` で再実行（約13分＋新セル分）。
3. **リント整備**: 新規コードは ruff の新規指摘を解消済み。リポジトリ全体の既存 ruff エラー（例: docstring 内の ρ 等）は既存スタイルであり本タスクでは触れていない。

### 5.2 中期

1. misim（C++カーネル）への mid-circuit reset 対応の設計検討。
2. tnsim トラジェクトリ法の並列化（現状は逐次実行で N=4 では非常に遅い）。
3. 元報告 §5.2 の A-1/A-3（物理モデル拡張）への着手。
