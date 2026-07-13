# GKSL補助qudit更新ダイナミクス 収束スタディ実行報告（2026-07-13）

本報告は `developing/GKSL補助qudit更新ダイナミクス次ステップ実装報告_2026-07-13.md` の「5. 次ステップの計画」§5.1 の実施結果を事実ベースでまとめたものである。誇張・虚偽は含まない。未完了項目は未完了と明記する。

## 1. 実施内容と実測結果

### 1.1 事前回帰テスト（完了・実測）

- `pip install -e .` + qiskit, qiskit-aer, pylatexenc, pytest, jupyter を導入。
- `python -m pytest test/python/simulation/test_reset.py test/python/simulation/test_tnsim.py test/python/simulation/test_dmsim.py -q` → **23 passed**（55.5 秒、実測）。
- `tutorials/test_gksl_simulators.py` フルスイート: **148 passed, 3 failed, 2 skipped**（49分50秒、実測）。失敗3件は本タスクと無関係の既存問題（`TestBosonGephZeroReduction::test_boson_g_eph_zero_traces` の 26.3 PiB メモリ確保エラー、`TestQuditGKSLCircuitBosonSimulator::test_gate_breakdown` の KeyError 'n_el_qutrits' 等）であり、reset/tnsim トラジェクトリ関連テストの失敗はない。

### 1.2 §5.1-1 収束スタディ: 性能問題の発見と対処（実測）

- 4プロセス並列実行時、BLAS のマルチスレッドと ProcessPoolExecutor が競合し、1ショットあたり約15秒（単独実行の約15倍遅い）となる問題を実測で確認した。
- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` を設定することで解消（1ショットあたり約1.2秒に改善、実測）。
- **推奨実行コマンド**（今後の再実行用）:

```bash
cd tutorials
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python run_5f_shots_convergence_study.py --processes 4
```

### 1.3 §5.1-1 収束スタディ: 縮小版（n_steps=20、完了・実測）

設定: t_max=100, n_steps=20, n_shots ∈ {50, 200, 800, 3200}, seed=12345, processes=4。

| n_shots | ‖Δρ‖_F | 1/√n | ratio | time [s] |
|---|---|---|---|---|
| 50 | 1.6424e-01 | 1.4142e-01 | 1.161 | 4.3 |
| 200 | 5.1200e-02 | 7.0711e-02 | 0.724 | 13.1 |
| 800 | 2.4506e-02 | 3.5355e-02 | 0.693 | 50.3 |
| 3200 | 1.1718e-02 | 1.7678e-02 | 0.663 | 191.8 |

- log-log フィット傾き **p = −0.625**（理論 −0.5）。O(1/√n) 収束と整合的（有限サンプルゆえのばらつきあり）。
- tnsim backend クロスチェック（n_trajectories=50）: ‖Δρ‖_F = **3.3947e-01** < 上界 5/√50 = 7.0711e-01、367.1 秒 → **PASSED**（実測）。

### 1.4 §5.1-1 収束スタディ: フル本番版（n_steps=100、ショット部完了・実測）

設定: t_max=100, n_steps=100, n_shots ∈ {50, 200, 800, 3200}, seed=12345, processes=4（デフォルトパラメータ）。

| n_shots | ‖Δρ‖_F | 1/√n | ratio | time [s] |
|---|---|---|---|---|
| 50 | 1.4517e-01 | 1.4142e-01 | 1.027 | 26.9 |
| 200 | 1.1793e-01 | 7.0711e-02 | 1.668 | 96.0 |
| 800 | 3.9350e-02 | 3.5355e-02 | 1.113 | 304.5 |
| 3200 | 2.1447e-02 | 1.7678e-02 | 1.213 | 1173.4 |

- log-log フィット傾き **p = −0.493**（理論 −0.5）。**理論値とよく一致**（実測）。
- 決定論参照（NumPy Stinespring、n_steps=100）: 53.2 秒。
- **tnsim backend クロスチェック（n_trajectories=50, n_steps=100）はセッション終了時点で実行中のまま完了していない**。数値は得られていない。縮小版（§1.3）では同クロスチェックが合格していることを付記する。

## 2. 未完了項目（正直な記載）

1. フル版（n_steps=100）の tnsim クロスチェック: 実行中のままセッション時間切れ。上記推奨コマンドで再実行可能（所要目安: ショット部含め全体で約60〜90分。ショット部の結果は再現される）。
2. ノートブック `tutorials/quantum_dynamics_gksl_comparison.ipynb`（前回コミットで 5e/5f/5f-2 セル改修済み）の `jupyter nbconvert --execute` による再実行: **未実施**。収束スタディを優先したため着手できなかった。
3. `tutorials/test_gksl_simulators.py`: 148 passed, 3 failed（**既存の boson 系問題で本タスクとは無関係**）, 2 skipped。reset/tnsim 関連の失敗はゼロ。
4. misim（C++）reset 対応、元報告 §5.2 中期項目: 未着手。

## 3. 次ステップの計画

### 3.1 短期

1. フル版収束スタディの tnsim クロスチェック完走（上記コマンド、`--n-shots 50` などでショット部を省略short-circuitはできないため、`--n-trajectories 50` のみ確認したい場合は `--n-shots 50` と組み合わせて時間短縮可能）。
2. ノートブックの `jupyter nbconvert --to notebook --execute --inplace tutorials/quantum_dynamics_gksl_comparison.ipynb` 実行（約13分＋5f-2セル約10分の見込み。`OMP_NUM_THREADS` 制限は不要）。
3. `test_gksl_simulators.py` の既存失敗3件（boson 系: 26.3 PiB メモリエラー、KeyError 'n_el_qutrits'）の原因調査・修正（本タスクの reset/tnsim 実装とは無関係の既存問題）。

### 3.2 中期

1. `run_5f_shots_convergence_study.py` に BLAS スレッド数の自動制限（`--processes > 1` 時に環境変数を子プロセスに設定）を組み込む（本セッションで実測した15倍の速度低下を恒久的に防ぐため）。
2. misim（C++カーネル）への mid-circuit reset 対応の設計検討。
3. 元報告 §5.2 の A-1/A-3（物理モデル拡張）への着手。
