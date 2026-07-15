# Iteration 48: 統合トロッターステップ量子回路の可視化

## 日付
2026-03-23

## 目標
`quantum_dynamics_gksl_comparison.ipynb` の量子回路可視化に、1トロッターステップ分の量子回路をハミルトニアンとStinespringで分けずに**1つの統合量子回路図**として表示する機能を追加する。

## 背景
- 現状: Cell 16 で各シナリオの1トロッターステップを複数のサブ回路（ハミルトニアン半ステップ、個別Stinespringチャネル）として分割表示
- 要望: 全ゲートを1つの回路図に統合して表示できるようにする（既存の分割表示も維持）

## 技術的課題分析

### レジスタサイズの不一致
各サブ回路は異なるqubit/qudit数を持つ:
- ハミルトニアン回路: システムqubit/quditのみ（補助qubitなし）
- Stinespring単一サイト: 2-3 qubit/qudit（システム＋補助）
- Stinespringペア: 3-5 qubit/qudit（システム＋補助）

### 解決策: 統合レジスタ
システムqubit/qudit + 1つの補助qubit/quditを持つ単一回路を構築:

| シナリオ | 方式 | システム | 補助 | 合計 |
|----------|------|----------|------|------|
| 3 (Qubit非ボソン) | Qiskit | 2N=8 qubit | 1 qubit | 9 qubit |
| 4 (Qubitボソン) | Qiskit | n_el+n_ph qubit | 1 qubit | n_sys+1 qubit |
| 5 (Qudit非ボソン) | MQT-Qudits | N=4 qutrit | 1 qutrit | 5 qutrit |
| 6 (Quditボソン) | MQT-Qudits | 2N=8 qudit | 1 qutrit | 9 qudit |

### 補助qubitリセット
- **Qiskit**: `circuit.reset(ancilla_idx)` で各Stinespringチャネル前に明示的リセット
- **MQT-Qudits**: ミッドサーキットリセット未対応。回路図はゲート構造を表示し、補助quditは各チャネル間で概念的に|0⟩にリセットされることを前提とする

## 実施内容

### 1. シミュレータクラスへのメソッド追加

#### `QubitGKSLCircuitSimulator.build_combined_trotter_step_circuit(dt)`
- ファイル: `tutorials/qubit_gksl_circuit_simulator.py`
- 9 qubit (8システム + 1補助) の単一 `QuantumCircuit` を構築
- ハミルトニアンゲートをシステムqubitに配置
- Stinespringゲートをシステム+補助qubitに配置（リセット付き）

#### `QuditGKSLCircuitSimulator.build_combined_trotter_step_circuit(dt)`
- ファイル: `tutorials/qudit_gksl_circuit_simulator.py`
- 5 qutrit (4システム + 1補助) の単一 MQT-Qudits 回路を構築
- `cu_one`/`cu_two`/`cu_multi` でゲート配置

#### `QuditGKSLCircuitBosonSimulator.build_combined_trotter_step_circuit(dt)`
- ファイル: `tutorials/qudit_gksl_circuit_boson_simulator.py`
- 9 qudit (4電子 + 4フォノン + 1補助) の単一回路を構築

### 2. ノートブック Cell 16 の更新
- 既存の分割表示の後に統合回路表示セクションを追加
- 各シナリオの統合回路図をPNG保存＋表示

### 3. 検証スクリプト作成
- `tutorials/run_iteration48_verification.py`
- 全シナリオの統合回路構築・ゲート数・レジスタサイズ検証
- 結果: `developing/verification_results/iteration48_results.txt`

## ゲート配置の正当性

### Qiskitの場合
- `circuit.append(gate, qubit_list)` はゲートのローカルqubitをグローバルqubitにマッピング
- ローカルqubit 0 → グローバル ancilla_idx
- ローカルqubit 1,2,... → グローバル システムqubit

### MQT-Quditsの場合
- `cu_two([system_k, ancilla_N], U_local)` は `cu_two([0, 1], U_local)` と同じユニタリ行列を使用
- ローカルqudit 0 → グローバル system_k
- ローカルqudit 1 → グローバル ancilla_N
- Kronecker積の順序が保存されるため、同一のユニタリ行列が正しく適用される

## 検証項目 (48チェック予定)

### Section 1: Scenario 3 (21チェック)
- 辞書キー確認 (8)
- 分割回路とのゲート数一致 (3)
- レジスタサイズ (3)
- 回路オブジェクト (4)
- 回帰値 (3)

### Section 2: Scenario 5 (17チェック)
- 辞書・ゲート数・レジスタ確認 (8)
- 回路構造 (4)
- 回帰値 (3)
- ゲートターゲット確認 (2)

### Section 3: Scenario 6 (11チェック)
- 辞書・レジスタ・ゲート数 (10)
- 補助quditターゲット確認 (1)

### Section 4: Scenario 4 (7チェック)
- 回路構築成功・サイズ・ゲート数 (7)

### Section 5: クロスシナリオ一貫性 (4チェック)
- S3とS5の物理的等価性 (3)
- qudit対qubitレジスタ効率 (1)

### Section 6: 可視化可能性 (6チェック)
- 各回路の描画可能性確認 (6)

## ステータス
- [x] コード実装完了
- [ ] ユーザーによるローカル検証実行待ち
- [ ] 検証結果に基づく修正（必要に応じて）
