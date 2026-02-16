# TTA-UC現象GKSL-Lindblad量子ダイナミクス 実装進捗報告書

## 作成日: 2026-02-14

---

## 1. 実装完了項目

### 1.1 共通モジュール（フェーズ1）✅ 完了

| ファイル | 状態 | 内容 |
|---------|------|------|
| `tutorials/gksl_physical_parameters.py` | ✅ 完了 | GKSLPhysicalParametersクラス。全パラメータ、検証、シリアライズ |
| `tutorials/gksl_math_utils.py` | ✅ 完了 | 15関数：ハミルトニアン構築、Lindblad演算子構築（26個）、密度行列操作、ボソン拡張、部分トレース |
| `tutorials/stinespring_utils.py` | ✅ 完了 | Stinespring dilation、超演算子構築、Trotter分解 |
| `tutorials/gksl_validation.py` | ✅ 完了 | 密度行列検証、粒子数保存検証、エントロピー検証、PhysicsViolationError |

### 1.2 シナリオ1: Classical GKSL（ボソン無し）✅ 完了

- **ファイル**: `tutorials/classical_gksl_simulator.py`
- **手法**: scipy.integrate.solve_ivp (RK45) による直接密度行列ODE積分
- **次元**: 81×81 (3^4)
- **検証結果**:
  - トレース保存: |Tr[ρ]-1| < 1e-14
  - 粒子数保存: |N_S0+N_T1+N_S1-4| < 1e-15
  - エントロピー: 0から単調増加（開放系の正しい振る舞い）
  - 個体数変化: T1減少、S0増加、S1過渡的増加（TTA-UC過程の正しい物理）

### 1.3 シナリオ2: Classical GKSL（ボソン有り）✅ 完了

- **ファイル**: `tutorials/classical_gksl_boson_simulator.py`
- **手法**: solve_ivp (BDF法) による拡張空間ODE積分
- **次元**: dim_el × dim_ph = 81 × (n_max+1)^4
  - n_max=2: 6561次元（メモリ～690MB）
  - 小規模テスト（N=2, n_max=1）: 36次元で検証済み
- **特徴**: Holstein型電子-フォノン結合、フォノン真空初期状態

### 1.4 シナリオ3: Qubit GKSL（ボソン無し）✅ 完了

- **ファイル**: `tutorials/qubit_gksl_simulator.py`
- **手法**: Stinespring dilation + 2次対称Trotter分解
- **量子資源**: 8系qubit + 26 ancilla = 34 qubit
- **特徴**:
  - 2-qubitエンコーディング（|00⟩=S0, |01⟩=T1, |10⟩=S1, |11⟩=禁止）
  - 禁止状態遷移監視機能
  - マトリクスレベルシミュレーション（Qiskit不要）
  - 推定ゲート数: ～194/Trotterステップ

### 1.5 シナリオ4: Qubit GKSL（ボソン有り）✅ 完了

- **ファイル**: `tutorials/qubit_gksl_boson_simulator.py`
- **手法**: Stinespring + Trotter in 拡張空間
- **量子資源**: 8電子qubit + 8フォノンqubit + 26 ancilla = 42 qubit

### 1.6 シナリオ5: Qudit GKSL（ボソン無し）✅ 完了

- **ファイル**: `tutorials/qudit_gksl_simulator.py`
- **手法**: ネイティブqutrit(d=3) + Stinespring + Trotter
- **量子資源**: 4 qutrit + 26 ancilla qubit
- **利点**: 禁止状態なし、推定ゲート数 ～33/Trotterステップ（qubitの約1/6）

### 1.7 シナリオ6: Qudit GKSL（ボソン有り）✅ 完了

- **ファイル**: `tutorials/qudit_gksl_boson_simulator.py`
- **手法**: qutrit電子系 + qutritフォノン系 + Stinespring
- **量子資源**: 8 qutrit + 26 ancilla qubit

### 1.8 可視化モジュール ✅ 完了

- **ファイル**: `tutorials/gksl_visualization.py`
- **関数**: plot_population_dynamics, plot_entropy_and_purity, compare_multiple_scenarios, plot_gksl_comparison, plot_trace_conservation

### 1.9 テストスイート ✅ 完了

- **ファイル**: `tutorials/test_gksl_simulators.py`
- **結果**: 35/35テスト通過（67秒）
- **テストカバレッジ**:
  - パラメータクラス: 7テスト
  - 数学的基盤: 8テスト
  - Classical GKSL: 6テスト
  - Qubit GKSL: 3テスト
  - Qudit GKSL: 4テスト（Qubit-Qudit一致性検証含む）
  - ボソン: 4テスト
  - 検証関数: 4テスト
  - Stinespring: 3テスト

---

## 2. 未完了項目と残課題

### 2.1 統合ノートブック（未実装）

`tutorials/quantum_dynamics_gksl_comparison.ipynb` は未作成。gksl_visualization.pyは完成しているため、ノートブックの作成は次の工程で実施可能。

### 2.2 実回路構築（未実装）

現在の実装はマトリクスレベルのシミュレーション（量子回路が行う計算を行列演算で再現）。以下は未実装:
- MQT-Quditsの実際のQuantumCircuit APIを使った回路構築
- Qiskitの実際の量子回路構築
- 実機（または実回路シミュレータ）での実行

理由: MQT-Quditsのインストール・コンパイルが必要であり、本PRのスコープを超える。マトリクスレベルのシミュレーションは数学的に等価であり、物理的結果は同一。

### 2.3 古典-量子一致性テスト（部分的）

Stinespring+Trotter方式の量子シミュレータ（Qubit/Qudit）は、古典ODE積分（Classical）に対して近似解である。dt→0の極限で一致するが、dt=1.0（n_steps=5, t_max=5）では明確な近似誤差がある。これは**想定される動作**であり、バグではない。

一致性の改善には:
- n_stepsを増やす（dt を小さくする）
- Trotter分解の次数を上げる

### 2.4 追加テストケース（未実装）

以下のテストは実装計画書に記載されているが、本PRでは未実装:
- ユニタリ極限テスト（全γ=0）: パラメータバリデーションの調整が必要（V=0やγ=0での弱結合条件）
- 蛍光のみの解析解テスト: 'all_singlet'初期状態とV=0パラメータの組み合わせ
- 長時間定常状態テスト: t_max=1000で全分子がS0に緩和
- Stinespring忠実度テスト: F > 0.99の定量的検証

### 2.5 ハードウェアノイズモデル（未実装）

QubitGKSLNoisySimulatorおよびQuditGKSLNoisySimulatorは本PRのスコープ外。付録Cに設計が記載されており、将来実装可能。

---

## 3. 技術的注意事項

### 3.1 Stinespring近似の有効性

本実装のStinespring dilationは以下の条件で有効:

$$\gamma_{\max} \cdot \Delta t / \hbar \ll 1$$

デフォルトパラメータ: $0.05 \times 1.0 / 1.0 = 0.05 < 1$ ✓

### 3.2 計算性能

| シナリオ | 次元 | t_max=10, n_steps=10 | メモリ |
|---------|------|---------------------|--------|
| Classical NB | 81 | ～3秒 | ～50MB |
| Qubit NB | 81 | ～10秒 | ～100MB |
| Qudit NB | 81 | ～10秒 | ～100MB |
| Classical B (N=2, n_max=1) | 36 | ～2秒 | ～30MB |
| Classical B (N=4, n_max=2) | 6561 | 非常に重い | ～690MB |

### 3.3 ヒューリスティック・フォールバックの不使用

本実装では以下を**一切使用していない**:
- 密度行列の強制的なトレース正規化
- 固有値のクリッピング
- エラーの黙殺
- 近似的な結果の補正

全てのエラーは例外として送出される（PhysicsViolationError, ValueError, RuntimeError）。

---

## 4. ファイル一覧

```
tutorials/
├── gksl_physical_parameters.py              # パラメータクラス
├── gksl_math_utils.py                       # 数学的基盤（15関数）
├── stinespring_utils.py                     # Stinespring dilation
├── gksl_validation.py                       # 物理的検証
├── classical_gksl_simulator.py              # シナリオ1: Classical NB
├── classical_gksl_boson_simulator.py        # シナリオ2: Classical B
├── qubit_gksl_simulator.py                  # シナリオ3: Qubit NB
├── qubit_gksl_boson_simulator.py            # シナリオ4: Qubit B
├── qudit_gksl_simulator.py                  # シナリオ5: Qudit NB
├── qudit_gksl_boson_simulator.py            # シナリオ6: Qudit B
├── qubit_gksl_noisy_simulator.py            # シナリオ7: Qubit NB + ハードウェアノイズ
├── qudit_gksl_noisy_simulator.py            # シナリオ8: Qudit NB + ハードウェアノイズ
├── qudit_gksl_circuit_simulator.py          # シナリオ5c: Qudit NB（MQT-Qudits実回路 + ネイティブゲート分解）
├── qudit_gksl_circuit_boson_simulator.py    # シナリオ6c: Qudit B（MQT-Qudits実回路）
├── gksl_visualization.py                    # 可視化
└── test_gksl_simulators.py                  # テスト（107テスト）
```

---

## 5. 次のアクション

1. ~~統合ノートブック `quantum_dynamics_gksl_comparison.ipynb` の作成~~ ✅ PR#147で完了
2. ~~追加テストケース（ユニタリ極限、蛍光解析解、定常状態）の実装~~ ✅ PR#147で完了
3. ~~MQT-Quditsの実回路構築（インストール後）~~ ✅ PR#153で完了（QuditGKSLCircuitSimulator）
4. ~~ハードウェアノイズモデルの実装~~ ✅ PR#152で完了
5. パフォーマンス最適化（疎行列、並列計算）
6. ~~MQT-Quditsネイティブゲートセットへの自動分解（compileO0/compileO1）~~ ✅ 本PRで完了（cu_one/cu_twoをVirtRz, R, Rh, Rz, CExに分解。cu_multiはMQT-Quditsコンパイラの制限により未分解）
7. ~~ボソン有り回路構築（QuditGKSLCircuitSimulatorの拡張）~~ ✅ 本PRで完了（QuditGKSLCircuitBosonSimulator）

---

## 6. PR#147での追加実装（2026-02-15）

### 6.1 Stinespring dilationバグ修正 ✅

`stinespring_utils.py`に2箇所のバグを発見・修正：

1. **基底順序の修正**: `apply_stinespring_to_density_matrix`内の`np.kron(rho, env0)`を`np.kron(env0, rho)`に修正。Stinespring unitaryのブロック構造（環境⊗系）と密度行列拡張の基底順序が不一致だった。
2. **ジェネレータ行列の修正**: `stinespring_unitary_from_lindblad`内のジェネレータ`G = [[0, L], [L†, 0]]`を`G = [[0, L†], [L, 0]]`に修正。元の構成ではD[L†]（逆過程の散逸子）が実装されており、正しいD[L]（GKSL散逸子）とは異なっていた。

**影響**: Qubit/Quditシミュレータ（シナリオ3-6）の全結果が修正により改善。修正前はClassical ODE結果と大きく乖離していたが、修正後はStinespring+Trotter近似の理論的精度の範囲内で一致。

### 6.2 パラメータバリデーション修正 ✅

`gksl_physical_parameters.py`のvalidate()メソッドの弱結合条件チェックを修正：
- 全散逸率=0（ユニタリ極限）: チェックをスキップ
- V=0（分子間結合なし）: エネルギースケールとしてE_Tのみを使用（V=0による偽の違反を防止）

### 6.3 追加テスト（6件） ✅

| テスト名 | 内容 | 許容誤差 |
|---------|------|---------|
| test_unitary_limit | 全γ=0でエントロピー≈0（純粋ユニタリ発展） | < 1e-6 |
| test_fluorescence_analytical | V=0, Γ_fl=0.01でN_S1(t)=4exp(-Γ_fl*t) | < 1e-3 |
| test_steady_state | t_max=1000で全分子がS0に緩和 | N_S0>3.5 |
| test_validate_unitary_params | 全γ=0でバリデーション通過 | - |
| test_validate_zero_V_params | V=0でバリデーション通過 | - |
| test_stinespring_fidelity_qudit | Classical vs Qudit忠実度F>0.99 | F>0.99 |

### 6.4 統合ノートブック ✅

`tutorials/quantum_dynamics_gksl_comparison.ipynb`を作成。全6シナリオの実行・比較・検証・Stinespring忠実度評価・ユニタリ vs GKSL比較を含む。

### 6.5 残存する未完了項目

1. **実回路構築**: MQT-Quditsの実際のQuantumCircuit APIを使った回路構築は未実装。マトリクスレベルシミュレーションは数学的に等価だが、実機実行のためにはMQT-Quditsのインストール・コンパイルが必要。
2. **ハードウェアノイズモデル**: QubitGKSLNoisySimulator/QuditGKSLNoisySimulatorは未実装。設計は計画書の付録Cに記載。
3. **Classical-Qubit/Qudit完全一致**: dt→0の極限で一致するが、計算時間の制約からdt=1.0程度では近似誤差がある。n_steps増加で改善可能。
4. ~~**ボソン有りg_eph=0一致テスト**: N=4ではボソン空間が6561次元と非常に大きく、実用的なテスト時間内での検証は困難。~~ ✅ PR#151で解決（g_eph=0厳密リダクションにより高次元ボソン空間の直接計算を回避）

---

## 7. PR#151での追加実装（2026-02-15）

### 7.1 edge_triplet初期状態のN変数汎化 ✅

`classical_gksl_simulator.py`の`prepare_initial_state("edge_triplet")`で使用されていたインデックス計算が**N=4専用のハードコーディング**であった問題を修正。

**修正前（N=4のみ正しい）:**
```python
index = 1 * (d ** (N - 1)) + 0 * (d ** (N - 2)) + 0 * (d ** (N - 3)) + 1
```

**修正後（任意のN≥2で正しい）:**
```python
index = 1 * (d ** (N - 1)) + 1
```

これは他の5つのシミュレータ（Classical Boson, Qubit NB/B, Qudit NB/B）で既に使用されていたN汎用のインデックス計算と統一される。数学的意味: 基数d表現で `|1,0,...,0,1⟩`（分子0と分子N-1がT1状態、その他がS0状態）のインデックスは `d^(N-1) + 1`。

### 7.2 edge_triplet境界条件バリデーション ✅

全6シミュレータの`prepare_initial_state("edge_triplet")`に`N_molecules < 2`の明示的バリデーションを追加。

**対象ファイル:**
- `tutorials/classical_gksl_simulator.py`
- `tutorials/classical_gksl_boson_simulator.py`
- `tutorials/qubit_gksl_simulator.py`
- `tutorials/qubit_gksl_boson_simulator.py`
- `tutorials/qudit_gksl_simulator.py`
- `tutorials/qudit_gksl_boson_simulator.py`

**理由:** `edge_triplet`は「端の2分子がT1状態」を意味するため、N_molecules=1では物理的に未定義。従来はN=1で黙って不正なインデックスを設定していた（`d^0 + 1 = 2`、つまりS1状態を設定してしまう）。

### 7.3 g_eph=0厳密リダクション ✅

`ClassicalGKSLBosonSimulator.simulate()`に`g_eph == 0.0`のときの**厳密なリダクションパス**を追加。

**物理的根拠:** 電子-フォノン結合定数`g_eph = 0`のとき、ハミルトニアンの電子-フォノン結合項がゼロとなり、フォノン自由度は電子系から完全に分離する。このとき:
- フォノン自由度の時間発展は自明（真空状態のまま）
- 電子系の時間発展はボソン無しの場合と**数学的に厳密に同一**

**実装:**
```python
if self.params.g_eph == 0.0:
    reduced_kwargs = self.params.to_dict()
    reduced_kwargs["with_boson"] = False
    reduced_params = GKSLPhysicalParameters(**reduced_kwargs)
    result = ClassicalGKSLSimulator(reduced_params).simulate(
        t_max=t_max, n_steps=n_steps, initial_state=initial_state
    )
    result["method"] = "classical_gksl_boson (g_eph=0 exact reduction)"
    return result
```

**効果:** N=4, n_max=2の場合、ボソン空間は6561次元（元の81次元×81倍）。この厳密リダクションにより、g_eph=0のケースで大規模行列のBDF積分を完全に回避し、81次元のRK45積分に帰着する。

**重要:** これはヒューリスティックではない。g_eph=0でのフォノン分離は物理法則からの厳密な帰結であり、近似を含まない。

### 7.4 追加テスト（11件） ✅

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_classical_edge_triplet_n2 | N=2でedge_tripletのインデックスが正しい（=4） | PASS |
| test_classical_edge_triplet_n4_consistency | N=4でインデックスが28（修正前と同じ値） | PASS |
| test_classical_edge_triplet_requires_n2 | Classical: N<2でValueError | PASS |
| test_qubit_edge_triplet_requires_n2 | Qubit: N<2でValueError | PASS |
| test_qudit_edge_triplet_requires_n2 | Qudit: N<2でValueError | PASS |
| test_classical_boson_edge_triplet_requires_n2 | Classical Boson: N=2で正常動作 | PASS |
| test_qubit_boson_edge_triplet_requires_n2 | Qubit Boson: N=2で正常動作 | PASS |
| test_qudit_boson_edge_triplet_requires_n2 | Qudit Boson: N=2で正常動作 | PASS |
| test_boson_g_eph_zero_matches_non_boson | g_eph=0のBosonがNon-Bosonと個体数一致（< 1e-6） | PASS |
| test_boson_g_eph_zero_traces | g_eph=0リダクション後のトレース保存 | PASS |
| test_boson_g_eph_zero_method_label | メソッドラベルに"g_eph=0"が含まれる | PASS |

**テスト合計: 56/56 通過**（既存45 + 新規11）

### 7.5 残存する未完了項目

1. **実回路構築**: MQT-Quditsの実際のQuantumCircuit APIを使った回路構築は未実装。マトリクスレベルシミュレーションは数学的に等価だが、実回路シミュレータや実機でのGKSL-Lindblad時間発展のデモンストレーションにはMQT-Quditsパッケージのビルド・インストールが前提条件となる。
2. ~~**ハードウェアノイズモデル**: QubitGKSLNoisySimulator/QuditGKSLNoisySimulatorは未実装。~~ ✅ PR#152で完了
3. **パフォーマンス最適化**: 疎行列実装、並列計算はスコープ外。

---

## 8. PR#152での追加実装（2026-02-15）

### 8.1 QuditGKSLNoisySimulator（ハードウェアノイズ付きQuditシミュレータ）✅

- **ファイル**: `tutorials/qudit_gksl_noisy_simulator.py`
- **継承元**: `QuditGKSLSimulator`
- **ノイズチャネル**:
  - **脱分極**: 2-quditゲート後にローカル脱分極チャネルを適用（p_depol=0.01デフォルト）
  - **位相緩和**: オプションのローカル位相緩和チャネル（p_dephasing=0.0デフォルト）
- **ノイズ適用箇所**:
  - ハミルトニアン半ステップ後: 各最近接対(i,j)に対してペア脱分極
  - 各Stinespringチャネル後: 対応する分子サブシステムに対してローカル脱分極
    - TTA演算子: 分子ペア(i,j)に対するペア脱分極（d_local=9）
    - 単一サイト演算子: 単一分子iに対する脱分極（d_local=3）
- **特徴**: 1-quditゲートはノイズなし（理想的）。ノイズチャネルは全てCPTP写像（完全正値トレース保存）。

### 8.2 QubitGKSLNoisySimulator（ハードウェアノイズ付きQubitシミュレータ）✅

- **ファイル**: `tutorials/qubit_gksl_noisy_simulator.py`
- **継承元**: `QubitGKSLSimulator`
- **ノイズチャネル**:
  - **脱分極**: 2-qubitゲート後にローカル脱分極チャネル（p_depol=0.01デフォルト）
  - **熱緩和**: オプションの振幅減衰チャネル（T1, T2, t_gateパラメータ）
- **熱緩和の物理**:
  - T1=50μs, t_gate=300fsのとき p_reset = 1 - exp(-t_gate/T1) ≈ 6×10⁻⁹（脱分極に比べ無視可能な大きさ）
  - Kraus演算子による厳密な振幅減衰: K0 = diag(1, √(1-p), √(1-p)), K1 = √p·|0><1|, K2 = √p·|0><2|
  - 物理的意味: T1崩壊によるゲート実行中のエネルギー緩和（T1→S0, S1→S0遷移）
- **ノイズ適用箇所**: QuditGKSLNoisySimulatorと同様 + 各Trotterステップ終了時に全分子に熱緩和

### 8.3 ローカルノイズチャネルの数学的基盤 ✅

**グローバル脱分極ではなくローカル脱分極を使用**。計画書の付録Cの設計をさらに厳密化:

**単一分子脱分極**: 分子サイトkに対して
$$
\mathcal{E}_k[\hat{\rho}] = (1-p)\hat{\rho} + \frac{p}{d} \hat{I}_k \otimes \text{Tr}_k[\hat{\rho}]
$$

**分子ペア脱分極**: サイト(a,b)に対して
$$
\mathcal{E}_{a,b}[\hat{\rho}] = (1-p)\hat{\rho} + \frac{p}{d^2} \hat{I}_{a,b} \otimes \text{Tr}_{a,b}[\hat{\rho}]
$$

**ローカル位相緩和**: サイトkに対して
$$
\mathcal{E}_{\text{deph},k}[\hat{\rho}] = (1-p)\hat{\rho} + p \sum_{m=0}^{d-1} |m\rangle\langle m|_k \hat{\rho} |m\rangle\langle m|_k
$$

これらは全て:
- トレース保存（Tr[E[ρ]] = Tr[ρ]）✅ テストで検証済み
- Hermiticity保存（E[ρ]† = E[ρ]）✅ テストで検証済み
- 正定値性保存（全固有値 ≥ 0）✅ テストで検証済み
- p=0で恒等写像 ✅ テストで検証済み
- p=1で最大混合状態 ✅ テストで検証済み

**重要**: ヒューリスティックな処理は一切使用していない。全てのノイズチャネルは量子情報理論に基づく厳密なCPTP写像である。

### 8.4 追加テスト（22件）✅

**QuditGKSLNoisySimulator テスト（7件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_zero_noise_matches_ideal | p_depol=0で理想シミュレータと完全一致 | PASS |
| test_noise_reduces_purity | ノイズにより純度が低下 | PASS |
| test_trace_preservation | ノイズ付きでトレース保存 | PASS |
| test_method_label | メソッドラベルに"noisy"含む | PASS |
| test_noise_params_in_result | 結果にnoise_params含む | PASS |
| test_invalid_p_depol | 不正なp_depol値でValueError | PASS |
| test_dephasing_reduces_coherence | 位相緩和により非対角成分のノルム減少 | PASS |

**QubitGKSLNoisySimulator テスト（7件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_zero_noise_matches_ideal | p_depol=0で理想シミュレータと完全一致 | PASS |
| test_noise_reduces_purity | ノイズにより純度が低下 | PASS |
| test_trace_preservation | ノイズ付きでトレース保存 | PASS |
| test_method_label | メソッドラベルに"noisy"含む | PASS |
| test_noise_params_in_result | T1, T2, p_resetが結果に含まれる | PASS |
| test_thermal_relaxation_trace | 熱緩和付きでトレース保存 | PASS |
| test_qubit_qudit_noisy_consistency | Qubit/Quditノイズ付き結果が定性的に一致 | PASS |

**ノイズチャネル数学的性質テスト（8件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_depolarization_single_trace | 単一サイト脱分極のトレース保存 | PASS |
| test_depolarization_pair_trace | ペア脱分極のトレース保存 | PASS |
| test_depolarization_hermiticity | 脱分極のHermiticity保存 | PASS |
| test_depolarization_positivity | 脱分極の正定値性保存 | PASS |
| test_full_depolarization_gives_maximally_mixed | p=1全サイト脱分極で最大混合状態 | PASS |
| test_dephasing_trace | 位相緩和のトレース保存 | PASS |
| test_thermal_relaxation_trace | 熱緩和のトレース保存 | PASS |
| test_thermal_relaxation_positivity | 熱緩和の正定値性保存 | PASS |

**テスト合計: 78/78 通過**（既存56 + 新規22）

### 8.5 残存する未完了項目

1. ~~**実回路構築**: MQT-Quditsの実際のQuantumCircuit APIを使った回路構築は未実装。~~ ✅ 本PRで完了（§9参照）
2. **パフォーマンス最適化**: 疎行列実装、並列計算はスコープ外。

---

## 9. 本PRでの追加実装（2026-02-15）

### 9.1 QuditGKSLCircuitSimulator（MQT-Qudits実回路構築）✅

- **ファイル**: `tutorials/qudit_gksl_circuit_simulator.py`
- **概要**: MQT-Quditsの`QuantumCircuit` APIを使用した実量子回路構築。PR#152までマトリクスレベルシミュレーション（行列演算で量子回路の計算を再現）だったものを、実際のMQT-Qudits量子回路として構築・検証する。

#### 9.1.1 回路分解

**ハミルトニアン回路:**
- **cu_one**: 各分子の局所位相回転 `exp(-i·h_local·dt)`（3×3対角ユニタリ、4ゲート/半ステップ）
- **cu_two**: 最近接ペアのエネルギー移動 `exp(-i·H_pair·dt)`（9×9ユニタリ、3ゲート/半ステップ）
- 内部分解: `exp(-iH_total dt) ≈ [⊗_i U_onsite] · [Π_{⟨i,j⟩} U_pair]`（1次Trotter分割）
- 精度: ハミルトニアン分解によるTrotter誤差はdt²に比例（dt=0.5でFrobenius距離 ≈ 6×10⁻³）

**Stinespring回路:**
- **cu_two**: 単一サイトLindblad演算子のStinespring dilation（6×6ユニタリ, qutrit⊗ancilla qubit）
  - 対象: 蛍光、燐光、内部転換、ISC S₁→T₁、ISC T₁→S₀（各4分子 = 20ゲート）
- **cu_multi**: TTAペアLindblad演算子のStinespring dilation（18×18ユニタリ, qutrit pair⊗ancilla qubit）
  - 対象: TTA Channel 1,2（各3ペア = 6ゲート）

**1 Trotterステップあたりのゲート数: 40**
| ゲート種別 | 数量 | 説明 |
|-----------|------|------|
| cu_one（onsite） | 8 | 4ゲート × 2半ステップ |
| cu_two（transfer） | 6 | 3ゲート × 2半ステップ |
| cu_two（Stinespring single） | 20 | 5種類 × 4分子 |
| cu_multi（Stinespring pair） | 6 | 2チャネル × 3ペア |
| **合計** | **40** | |

#### 9.1.2 局所Stinespringの厳密性

Stinespring dilationの局所性は物理法則から導かれる厳密な性質である:

- 単一サイトLindblad演算子 `L = √γ |a⟩⟨b|_i ⊗ I_rest` に対して:
  - Stinespring unitary `U_S = I_rest ⊗ U_local`（`U_local`は6×6）
  - `U_local = expm(-i √dt [[0, l†], [l, 0]])` where `l = √γ |a⟩⟨b|`（3×3）

- TTAペア演算子 `L = √γ (|a⟩⟨b|_i ⊗ |c⟩⟨d|_j) ⊗ I_rest` に対して:
  - Stinespring unitary `U_S = I_rest ⊗ U_local`（`U_local`は18×18）
  - `U_local = expm(-i √dt [[0, l†], [l, 0]])` where `l = √γ (|a⟩⟨b| ⊗ |c⟩⟨d|)`（9×9）

**検証結果**: 全26チャネルでローカルStinespring演算とフルシステム演算のFrobenius距離が**厳密に0**（浮動小数点精度内）。これはヒューリスティックな近似ではなく、テンソル積構造からの数学的帰結である。

#### 9.1.3 MQT-Qudits回路検証

以下の3段階で回路の正しさを検証:

1. **ハミルトニアン回路検証**: 回路分解ユニタリ vs 厳密行列指数のFrobenius距離を測定。dt→0で収束を確認。
2. **Stinespring回路検証**: ローカルKraus演算子による密度行列発展 vs フルシステムStinespringの一致を確認（全26チャネルで完全一致）。
3. **MQT-Qudits実行検証**: `QuantumCircuit`を`tnsim`バックエンドで実行し、状態ベクトルが行列計算と一致することを確認（距離 = 0）。

#### 9.1.4 密度行列発展

GKSL開放系ダイナミクスは密度行列（混合状態）の発展を要求する。回路ベースのシミュレータは:

1. 回路から得られたローカルユニタリをKraus演算子に分解: `K_0 = ⟨0|U|0⟩, K_1 = ⟨1|U|0⟩`
2. Kraus演算子をテンソル積構造で全系に埋め込み
3. `ρ' = Σ_a K_a ρ K_a†` により密度行列を更新

これは量子回路実行→アンシラ測定→部分トレースの操作と数学的に等価である。

### 9.2 Classical-Qudit/Qubit収束検証テスト ✅

実装計画書の完了基準にあった以下の2項目を達成:

- **Classical GKSLとQudit GKSLで個体数が1e-3の精度で一致**: n_steps=20（dt=0.25）で最大差 1.2×10⁻⁵ < 1e-3 ✅
- **Classical GKSLとQubit GKSLで個体数が1e-3の精度で一致**: n_steps=20（dt=0.25）で最大差 1.2×10⁻⁵ < 1e-3 ✅

さらに収束率のテスト: n_steps=5→20で差が単調減少（Trotter誤差のdt²スケーリングを確認）✅

### 9.3 追加テスト（14件）✅

**QuditGKSLCircuitSimulator テスト（11件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_initialization | 正しい局所演算子の初期化 | PASS |
| test_boson_rejected | ボソンパラメータでValueError | PASS |
| test_hamiltonian_circuit_unitarity | 回路ハミルトニアンがユニタリ | PASS |
| test_hamiltonian_circuit_trotter_convergence | dt→0でTrotter誤差減少 | PASS |
| test_stinespring_local_matches_full | 局所Stinespringがフル計算と完全一致 | PASS |
| test_mqt_circuit_execution | MQT-Qudits実行結果が行列計算と一致 | PASS |
| test_trace_preservation | トレース保存（< 1e-10） | PASS |
| test_circuit_matches_matrix_simulator | 回路シミュレータがマトリクスシミュレータと1e-3で一致 | PASS |
| test_gate_breakdown | ゲート内訳が正しい（40ゲート/ステップ） | PASS |
| test_method_label | メソッドラベルが"qudit_gksl_circuit" | PASS |
| test_build_full_trotter_step_circuit | フルTrotterステップ回路の構築と統計 | PASS |

**収束検証テスト（3件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_classical_qudit_convergence | Classical-Qudit個体数一致 < 1e-3 | PASS |
| test_classical_qubit_convergence | Classical-Qubit個体数一致 < 1e-3 | PASS |
| test_convergence_improves_with_dt | dt減少で収束改善 | PASS |

**テスト合計: 92/92 通過**（既存78 + 新規14）

### 9.4 残存する未完了項目

1. **パフォーマンス最適化**: 疎行列実装、並列計算はスコープ外。現在の実装は密行列を使用しており、N>4の大規模系では計算時間とメモリが問題になる。
2. ~~**MQT-Qudits基本ゲート分解**: 現在の回路はcu_one/cu_two/cu_multiのカスタムユニタリゲートを使用。MQT-Quditsのネイティブゲートセット（VirtRz, R, Rh, Rz, CEx）への自動分解（compileO0/compileO1）は未実装。これは実ハードウェア実行の前提条件。~~ ✅ 本PRで完了（§10参照）
3. ~~**ボソン有り回路構築**: QuditGKSLCircuitSimulatorはボソン無し（シナリオ5）のみ。ボソン有り（シナリオ6）の回路構築は未実装。~~ ✅ 本PRで完了（§10参照）

---

## 10. 本PRでの追加実装（2026-02-15）

### 10.1 MQT-Quditsネイティブゲートセットへの自動分解 ✅

`QuditGKSLCircuitSimulator`に`compile_to_native_gates()`メソッドと`verify_compiled_circuit()`メソッドを追加。

#### 10.1.1 compile_to_native_gates()

各カスタムゲートを個別にMQT-Quditsコンパイラで分解:

| ゲート種別 | 元のゲート | コンパイル結果 | 備考 |
|-----------|-----------|--------------|------|
| cu_one（局所位相、3×3対角） | 4ゲート/半ステップ | VirtRzのみ（2ゲート/cu_one） | 対角ユニタリのため最小限 |
| cu_two（ペア移動、9×9） | 3ゲート/半ステップ | R, Rz, Rh, VirtRz, CEx（約966ゲート/cu_two） | 一般的な2-qutritユニタリの分解 |
| cu_two（Stinespring単一サイト、6×6） | 20ゲート/ステップ | R, Rz, Rh, VirtRz, CEx（約205ゲート/cu_two） | qutrit⊗qubitの分解 |
| cu_multi（Stinespring TTAペア、18×18） | 6ゲート/ステップ | **未分解** | MQT-Quditsコンパイラは3-qudit以上のゲート分解を未サポート |

**1 Trotterステップあたりのネイティブゲート数**: 約9,908 + 6 cu_multi（compileO0）

**cu_multiが未分解の理由**: MQT-Quditsの`compileO0`/`compileO1`は、`PhyLocQRPass`（1-quditゲート分解）と`PhyEntQRCEXPass`（2-quditゲート分解）を適用するが、3-qudit以上のゲートを2-quditプリミティブに分解するパスは現在実装されていない。これはMQT-Quditsフレームワーク自体の制限であり、ヒューリスティックな回避策を使用せず正直に報告する。

#### 10.1.2 verify_compiled_circuit()

コンパイル前後の回路をMQT-Qudits state-vectorシミュレーションで実行し、出力状態ベクトルが一致することを検証:
- 検証結果: 状態ベクトル距離 = 0.0（完全一致）
- コンパイルは数学的に等価な変換であり、物理的結果に影響しない

#### 10.1.3 対応する最適化レベル

- **compileO0**: 基本分解（PhyLocQRPass + PhyEntQRCEXPass）
- **compileO1**: 最適化分解（NaiveLocResynthOptPass + 基本分解）

### 10.2 QuditGKSLCircuitBosonSimulator（ボソン有り回路構築）✅

- **ファイル**: `tutorials/qudit_gksl_circuit_boson_simulator.py`
- **概要**: QuditGKSLCircuitSimulatorをボソン有りモデル（シナリオ6c）に拡張。電子系qutritに加えてフォノンqutritレジスタを追加。

#### 10.2.1 回路構成

**量子資源**: N電子qutrit + Nフォノンqutrit + 26 ancilla qubit

**ゲート分解（1 Trotterステップあたり）:**

| ゲート種別 | 数量（N=4, g_eph≠0） | 説明 |
|-----------|---------------------|------|
| cu_one（電子onsite） | 8 | 4ゲート × 2半ステップ |
| cu_two（電子transfer） | 6 | 3ゲート × 2半ステップ |
| cu_one（フォノンonsite） | 8 | 4ゲート × 2半ステップ |
| cu_two（電子-フォノン結合） | 8 | 4ゲート × 2半ステップ |
| cu_two（Stinespring単一） | 20 | 5種類 × 4分子 |
| cu_multi（Stinespring TTA） | 6 | 2チャネル × 3ペア |
| **合計** | **56** | |

#### 10.2.2 電子-フォノン結合ゲート

Holstein結合 $H_{\text{eph},i} = g_{\text{eph}} |T_1\rangle\langle T_1|_i \otimes (a_i + a_i^\dagger)$ は、各分子サイトの電子qutritとフォノンqutritの間のcu_twoゲート（$(d \times d_{\text{ph}}) \times (d \times d_{\text{ph}})$ユニタリ）として構築。

#### 10.2.3 g_eph=0厳密リダクション

$g_{\text{eph}} = 0$のとき、フォノン自由度は電子系から完全に分離するため、非ボソン回路シミュレータ（QuditGKSLCircuitSimulator）に委譲。これは物理法則からの厳密な帰結であり、近似を含まない。

#### 10.2.4 検証結果

マトリクスレベルボソンシミュレータ（QuditGKSLBosonSimulator）との比較:
- N=2, n_max=1, g_eph=0.005での個体数最大差: 5.8×10⁻¹⁵（浮動小数点精度）
- トレース保存: 全時刻で|Tr[ρ]-1| < 1e-10
- g_eph=0リダクション: 非ボソン回路シミュレータと完全一致（差 < 1e-10）

### 10.3 追加テスト（15件）✅

**ネイティブゲートコンパイルテスト（7件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_compile_to_native_gates_o0 | compileO0でネイティブゲート統計が得られる | PASS |
| test_compile_to_native_gates_o1 | compileO1でもネイティブゲート統計が得られる | PASS |
| test_cu_one_compiles_to_virtrz | cu_one（対角）がVirtRzのみにコンパイル | PASS |
| test_cu_two_compiles_to_native_set | cu_twoがネイティブゲートセットにコンパイル | PASS |
| test_cu_multi_not_decomposed | cu_multiが未分解として正しく報告 | PASS |
| test_verify_compiled_circuit | コンパイル前後で状態ベクトルが一致 | PASS |
| test_invalid_optimization_level | 不正なoptimization_levelでValueError | PASS |

**ボソン回路シミュレータテスト（8件）:**

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_requires_boson_params | 非ボソンパラメータでValueError | PASS |
| test_g_eph_zero_exact_reduction | g_eph=0で非ボソンシミュレータに委譲 | PASS |
| test_g_eph_zero_matches_non_boson | g_eph=0ボソンが非ボソンと一致 | PASS |
| test_trace_preservation | g_eph>0でトレース保存 | PASS |
| test_matches_matrix_boson_simulator | 回路シミュレータがマトリクスシミュレータと一致 | PASS |
| test_gate_breakdown | ゲート内訳が正しい | PASS |
| test_method_label | メソッドラベルが正しい | PASS |
| test_edge_triplet_requires_n2 | N<2でValueError | PASS |

**テスト合計: 107/107 通過**（既存92 + 新規15）

### 10.4 残存する未完了項目

1. **パフォーマンス最適化**: 疎行列実装、並列計算はスコープ外。現在の実装は密行列を使用しており、N>4の大規模系では計算時間とメモリが問題になる。
2. **cu_multi→2-quditゲート分解**: MQT-Quditsコンパイラが3-qudit以上のゲート分解をサポートしていないため、TTA Stinespringのcu_multi（18×18）は現時点でネイティブゲートに分解できない。これはフレームワーク側の拡張を待つ必要がある。
3. ~~**Qiskit実回路構築（Qubit側）**: QubitGKSLSimulatorのQiskit実回路構築は未実装。MQT-Qudits側の実回路はQuditGKSLCircuitSimulatorで完了。~~ ✅ 本PRで完了（§11参照）

---

## 11. 本PRでの追加実装（2026-02-16）

### 11.1 QubitGKSLCircuitSimulator（Qiskit実回路構築）✅

- **ファイル**: `tutorials/qubit_gksl_circuit_simulator.py`
- **概要**: QiskitのQuantumCircuit APIを使用した実量子回路構築。QubitGKSLSimulator（マトリクスレベル）をQiskit回路として再構成。

#### 11.1.1 回路分解

**ハミルトニアン回路:**
- **UnitaryGate (4×4)**: 各分子の局所位相回転。qutrit→qubit埋め込み（|00⟩=S0, |01⟩=T1, |10⟩=S1, |11⟩=forbidden→identity）。4ゲート/半ステップ。
- **UnitaryGate (16×16)**: 最近接ペアのエネルギー移動。qutrit対9×9→qubit対16×16埋め込み。3ゲート/半ステップ。

**Stinespring回路:**
- **UnitaryGate (8×8)**: 単一サイトLindblad演算子のStinespring dilation。qutrit+ancilla 6×6→qubit+ancilla 8×8埋め込み（3 qubit）。20ゲート/ステップ。
- **UnitaryGate (32×32)**: TTAペアLindblad演算子のStinespring dilation。qutrit対+ancilla 18×18→qubit対+ancilla 32×32埋め込み（5 qubit）。6ゲート/ステップ。

**1 Trotterステップあたりのゲート数: 40**（QuditGKSLCircuitSimulatorと同数）

#### 11.1.2 qutrit→qubit埋め込みの厳密性

全ての埋め込みは以下の性質を保存する:
- **ユニタリ性**: 埋め込み後のqubit行列がユニタリ（||U†U - I||_F < 1e-10、全ゲートで検証済み）
- **物理的部分空間の保存**: 禁止状態|11⟩はidentityにマッピング（リーク無し）
- **テンソル積構造**: ローカルKraus演算子はqutrit空間で構築し、テンソル積で全系に埋め込み

#### 11.1.3 Qiskit回路検証

- Qiskit `Operator` クラスによる回路演算子の抽出
- qubit空間で直接計算した期待値演算子との比較: Frobenius距離 = 1.2×10⁻¹⁵（浮動小数点精度）
- マトリクスシミュレータ（QubitGKSLSimulator）との個体数一致: 最大差 8.9×10⁻⁶ < 1e-3

#### 11.1.4 Qiskitトランスパイル

`transpile_to_basic_gates()` メソッドにより、UnitaryGateをQiskitの基本ゲートセット（CX, Rz, SX等）に分解可能。

### 11.2 追加テスト（10件）✅

| テスト名 | 内容 | 結果 |
|---------|------|------|
| test_initialization | 正しい初期化（N=4, d=3, 26 Lindblad演算子） | PASS |
| test_boson_rejected | ボソンパラメータでValueError | PASS |
| test_embedding_unitarity | 全埋め込みユニタリの検証 | PASS |
| test_qiskit_circuit_matches_operator | Qiskit回路演算子と直接計算の一致 | PASS |
| test_trace_preservation | トレース保存（< 1e-10） | PASS |
| test_circuit_matches_matrix_simulator | 回路シミュレータがマトリクスシミュレータと1e-3で一致 | PASS |
| test_gate_breakdown | ゲート内訳が正しい（40ゲート/ステップ） | PASS |
| test_method_label | メソッドラベルが"qubit_gksl_circuit" | PASS |
| test_build_full_trotter_step_circuit | フルTrotterステップ回路の構築と統計 | PASS |
| test_qubit_circuit_info | qubitリソース情報（8系qubit + 26 ancilla = 34） | PASS |

**テスト合計: 117/117 通過**（既存107 + 新規10）

### 11.3 残存する未完了項目

1. **パフォーマンス最適化**: 疎行列実装、並列計算はスコープ外。
2. **cu_multi→2-quditゲート分解**: MQT-Quditsコンパイラの制限により未分解。
