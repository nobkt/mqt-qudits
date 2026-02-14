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
├── gksl_physical_parameters.py     # パラメータクラス
├── gksl_math_utils.py              # 数学的基盤（15関数）
├── stinespring_utils.py            # Stinespring dilation
├── gksl_validation.py              # 物理的検証
├── classical_gksl_simulator.py     # シナリオ1: Classical NB
├── classical_gksl_boson_simulator.py  # シナリオ2: Classical B
├── qubit_gksl_simulator.py         # シナリオ3: Qubit NB
├── qubit_gksl_boson_simulator.py   # シナリオ4: Qubit B
├── qudit_gksl_simulator.py         # シナリオ5: Qudit NB
├── qudit_gksl_boson_simulator.py   # シナリオ6: Qudit B
├── gksl_visualization.py           # 可視化
└── test_gksl_simulators.py         # テスト（35テスト）
```

---

## 5. 次のアクション

1. 統合ノートブック `quantum_dynamics_gksl_comparison.ipynb` の作成
2. 追加テストケース（ユニタリ極限、蛍光解析解、定常状態）の実装
3. MQT-Quditsの実回路構築（インストール後）
4. ハードウェアノイズモデルの実装
5. パフォーマンス最適化（疎行列、並列計算）
