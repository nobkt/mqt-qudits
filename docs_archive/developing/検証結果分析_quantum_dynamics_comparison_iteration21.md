# 検証結果分析: Quantum Dynamics Comparison - Iteration 21

## 1. 背景

### 1.1 Iteration 20の結果

Iteration 20ではCell 32/36にレジスタ情報表示を追加し（Issue U/V）、
GPS回帰テストの誤検証を修正し（Issue W）、検証項目を拡充した（Issue X）。
34/34チェックがPASSした。

### 1.2 Iteration 20の検証結果分析で発見された問題

Iteration 20の検証結果`iteration20_shot_display_20260313T061734Z.json`と
ノートブックセル出力を詳細に分析した結果、以下の問題を発見した。

## 2. 発見された問題の詳細

### 2.1 Issue Y: Cell 36（Qubit Shot）の測定結果表示欠落

**ファイル**: `quantum_dynamics_gksl_comparison.ipynb` Cell 36

**問題**: Cell 32（QuditGKSLShotSimulator）は測定カウント上位5件を表示しているが、
Cell 36（QubitGKSLShotSimulator）では表示されていなかった。
また、QubitGKSLShotSimulatorは`forbidden_count`を返すが、Cell 36では
表示されておらず、ノイズなしケースで禁止状態遷移がゼロであることを
読者が確認できなかった。

| 項目 | Cell 32 (Qudit Shot) | Cell 36 (Qubit Shot) 修正前 | Cell 36 修正後 |
|------|---------------------|---------------------------|---------------|
| System qudits/qubits | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Ancilla | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Total | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Gates/step | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Final populations | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Trace conservation | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Measurement counts | ✓ 表示 | ✗ 未表示 | ✓ 表示 |
| Forbidden count | N/A（quditは禁止状態なし） | ✗ 未表示 | ✓ 表示 |

**影響**: qubit shot simulatorの測定結果分布が読者に伝わらず、
ノイズなしケースで禁止状態がゼロであることを確認できなかった。

### 2.2 Issue Z: Section 5検証項目の不足

**ファイル**: `run_iteration20_verification.py` Section 5

**問題**: Section 5の結果辞書完全性チェックで以下が不足していた：

| 欠落していた検証 | 期待値 |
|---------------|--------|
| QubitGKSL result n_total_qubits | 34 |
| QubitGKSLBoson result n_total_qubits | 20 |
| QubitGKSLShot forbidden_count（ノイズなし） | 0 |
| QubitGKSLShot result dict 'counts' key | 存在すること |
| QubitGKSLShot result dict 'forbidden_count' key | 存在すること |

## 3. 修正内容

### 3.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `quantum_dynamics_gksl_comparison.ipynb` Cell 36 | 測定カウント・forbidden_count表示追加（Issue Y修正） |
| `run_iteration21_verification.py` | 42項目の検証スクリプト（新規、Issue Z修正） |

### 3.2 変更不要

- シミュレータのPythonファイル: コードに問題なし（コードバグなし）
- `quantum_dynamics_complete_comparison.ipynb`: GKSLシミュレータを直接使用しない
- Cell 32: 変更不要（既に正しい表示）
- `run_iteration20_verification.py`: 旧バージョンとして保持

## 4. 検証結果

`tutorials/run_iteration21_verification.py` で42項目を検証：

### Section 1: n_total正確性（8項目） - 回帰

| # | チェック項目 | 結果 |
|---|-------------|------|
| 1 | qudit DM n_total = 30 | PASS |
| 2 | qubit DM n_total = 34 | PASS |
| 3 | qudit boson N=2 n_total = 16 | PASS |
| 4 | qubit boson N=2 n_total = 20 | PASS |
| 5 | qudit shot n_total = 30 | PASS |
| 6 | qubit shot n_total = 34 | PASS |
| 7 | qudit boson N=4 n_total = 34 | PASS |
| 8 | qubit boson N=4 n_total = 42（数式計算） | PASS |

### Section 2: クロスチェック（4項目） - 回帰

| # | チェック項目 | 結果 |
|---|-------------|------|
| 9 | 非ボソンアンシラ一致 (26=26) | PASS |
| 10 | ボソンアンシラ一致 (12=12) | PASS |
| 11 | qudit DM/Shot n_total一致 (30=30) | PASS |
| 12 | qubit DM/Shot n_total一致 (34=34) | PASS |

### Section 3: GPSレジスタ回帰テスト（8項目） - 回帰

| # | チェック項目 | 結果 |
|---|-------------|------|
| 13 | qudit DM GPS = 66 | PASS |
| 14 | qubit DM GPS = 388 | PASS |
| 15 | qudit boson N=2 GPS = 38 | PASS |
| 16 | qubit boson N=2 GPS = 220 | PASS |
| 17 | qudit shot GPS = 66 | PASS |
| 18 | qubit shot GPS = 388 | PASS |
| 19 | qudit boson N=4 GPS = 82 | PASS |
| 20 | qubit boson N=4 GPS = 484（数式計算） | PASS |

### Section 4: ノートブックセル検査（13項目） - Issue Y修正確認含む

| # | チェック項目 | 結果 |
|---|-------------|------|
| 21 | Cell 6 n_total_qudits表示 | PASS |
| 22 | Cell 8 n_total_qubits表示 | PASS |
| 23 | Cell 12 n_total_qudits表示 | PASS |
| 24 | Cell 14 n_total_qubits表示 | PASS |
| 25 | Cell 32 n_total_qudits表示 | PASS |
| 26 | Cell 32 n_system_qudits表示 | PASS |
| 27 | Cell 32 n_ancilla_qudits表示 | PASS |
| 28 | Cell 32 counts表示 | PASS |
| 29 | Cell 36 n_total_qubits表示 | PASS |
| 30 | Cell 36 n_ancilla表示 | PASS |
| 31 | Cell 36 Trace conservation表示 | PASS |
| 32 | Cell 36 counts表示 **[NEW]** | PASS |
| 33 | Cell 36 forbidden_count表示 **[NEW]** | PASS |

### Section 5: 結果辞書完全性（6項目） - Issue Z修正含む

| # | チェック項目 | 結果 |
|---|-------------|------|
| 34 | qudit DM result n_total = 30 | PASS |
| 35 | qubit DM result n_total = 34 **[NEW]** | PASS |
| 36 | qudit boson result n_total = 16 | PASS |
| 37 | qubit boson result n_total = 20 **[NEW]** | PASS |
| 38 | qudit shot result n_total = 30 | PASS |
| 39 | qubit shot result n_total = 34 | PASS |

### Section 6: Qubit shot forbidden_count検証（3項目） - 全て新規

| # | チェック項目 | 結果 |
|---|-------------|------|
| 40 | qubit shot forbidden_count = 0（ノイズなし） **[NEW]** | PASS |
| 41 | qubit shot result dict has 'counts' **[NEW]** | PASS |
| 42 | qubit shot result dict has 'forbidden_count' **[NEW]** | PASS |

## 5. 詳細分析

### 5.1 コードバグの不在

全シミュレータコードを詳細に分析した結果、以下を確認した：

- **GPS計算**: DM版とShot版で完全に一致（qudit: 66, qubit: 388, boson N=2: 38/220, boson N=4: 82/484）
- **n_total計算**: 全8シミュレータで正しいアンシラ込みの値
- **d_anc設定**: qudit=3, qubit=2（各シミュレータ内で一貫性あり）
- **per_molecule_populations順序**: Kronecker積の慣例（big-endian）に沿っており正しい
- **Stinespring測定**: 物理サブスペースを保存（ノイズなし時forbidden_count=0で確認）
- **Trace保存**: ノイズなし shot sim で max|Tr-1| = 1.15e-14（DM sim: 2.0e-14より良好）

### 5.2 ノートブック出力結果の分析

Cell 32とCell 36のノイズなし結果を比較：

| 項目 | Cell 32 (Qudit Shot) | Cell 36 (Qubit Shot) |
|------|---------------------|---------------------|
| Final N_S0 | 3.2400 | 3.2400 |
| Final N_T1 | 0.4600 | 0.4600 |
| Final N_S1 | 0.3000 | 0.3000 |
| max|Tr-1| | 1.15e-14 | 1.15e-14 |
| GPS | 66 | 388 |
| Total gates | 6,600 | 38,800 |
| Elapsed time | 194.56s | 260.52s |

母集団動態が一致する（同じ物理を記述）ことは正しく、
GPS差（66 vs 388 = 5.9倍）はqubit encodingのオーバーヘッドを正確に反映。

## 6. 残課題

### 6.1 qubit回路シミュレータの抽象化レベル不整合（継続）

Iteration 17/18/19/20から継続の課題。qubit回路シミュレータの`gates_per_step`（複合
ゲート単位66）とqubit DM/shotの`estimated_gates_per_step`（基本ゲート推定388）
は異なる抽象化レベルの測定値。

### 6.2 ユーザーによるノートブック再実行

修正後のノートブックをユーザーがローカルで実行し、Cell 36で
`Forbidden count: 0 / 1000`と
`Measurement counts (top 5): {...}`が正しく表示されることを確認する必要がある。

## 7. 次のステップ

1. ユーザーがiteration 21の検証スクリプトをローカルで実行し結果を確認
2. ノートブックを再実行してCell 36の測定カウント・forbidden_count表示を目視確認
3. qubit回路シミュレータの抽象化レベル統一を検討（次イテレーション候補）
