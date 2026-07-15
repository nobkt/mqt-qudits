# 検証結果分析: Quantum Dynamics Comparison - Iteration 20

## 1. 背景

### 1.1 Iteration 19の結果

Iteration 19では`n_total_qudits`がアンシラを含むよう修正し（Bug R）、
非ボソン・ショットシミュレータに`n_total_qudits`を追加し（Issue S）、
ノートブックのDM/ボソンセル（Cell 6, 12）に「Total qudits」表示を追加した（Issue T）。
30/30チェックがPASSした。

### 1.2 Iteration 19の検証結果分析で発見された問題

Iteration 19の検証結果`iteration19_n_total_qudits_20260313T045845Z.json`と
ノートブックセルを注意深く分析した結果、以下の4つの問題を発見した。

## 2. 発見された問題の詳細

### 2.1 Issue U: Cell 32（Qudit Shot）のレジスタ情報表示欠落

**ファイル**: `quantum_dynamics_gksl_comparison.ipynb` Cell 32

**問題**: Cell 32（QuditGKSLShotSimulator）で`n_system_qudits`、`n_ancilla_qudits`、
`n_total_qudits`が表示されていなかった。

| 項目 | Cell 6 (Qudit DM) | Cell 32 (Qudit Shot) 修正前 | Cell 32 修正後 |
|------|------------------|---------------------------|---------------|
| System qutrits | ✓ 表示 | ✗ 未表示 | ✓ 表示 |
| Ancilla qudits | ✓ 表示 | ✗ 未表示 | ✓ 表示 |
| Total qudits | ✓ 表示 | ✗ 未表示 | ✓ 表示 |

**影響**: ショットベースシミュレーションのレジスタ数が読者に伝わらず、
密度行列版との比較が困難だった。

### 2.2 Issue V: Cell 36（Qubit Shot）のレジスタ情報表示欠落

**ファイル**: `quantum_dynamics_gksl_comparison.ipynb` Cell 36

**問題**: Cell 36（QubitGKSLShotSimulator）で`n_ancilla`と`n_total_qubits`が
表示されておらず、`Trace conservation`も欠落していた。

| 項目 | Cell 8 (Qubit DM) | Cell 36 (Qubit Shot) 修正前 | Cell 36 修正後 |
|------|------------------|---------------------------|---------------|
| System qubits | ✓ 表示 | ✓ 表示 | ✓ 表示 |
| Ancilla qubits | ✓ 表示 | ✗ 未表示 | ✓ 表示 |
| Total qubits | ✓ 表示 | ✗ 未表示 | ✓ 表示 |
| Trace conservation | ✓ 表示 | ✗ 未表示 | ✓ 表示 |

### 2.3 Issue W: Iteration 19検証スクリプトのチェック#30が誤った項目を検証

**ファイル**: `run_iteration19_verification.py` 行288-292

**問題**: チェック名が`regression_qudit_boson_n4_gps_82`（GPS=82の回帰テスト）
だが、実際には`n_total_qudits == 34`を検証していた。

```python
# 修正前（iteration 19）
add(
    "regression_qudit_boson_n4_gps_82",      # ← GPSのテストを名乗る
    sim_qudit_boson_n4.n_total_qudits == 34,  # ← 実際はn_totalをテスト
    ...
)
```

**影響**: QuditGKSLBoson N=4のGPS=82の回帰テストが実質的に未実行だった。
30/30 PASSの結果が、GPS回帰については偽陽性を含んでいた。

### 2.4 Issue X: 検証項目の不足

Iteration 19の検証スクリプトに以下の項目が不足していた：

| 欠落していた検証 | 期待値 |
|---------------|--------|
| Cell 32のn_total_qudits表示 | n_total_quditsを含むこと |
| Cell 36のn_total_qubits表示 | n_total_qubitsを含むこと |
| Cell 36のTrace conservation表示 | trace情報を含むこと |
| QubitGKSLShot n_total_qubits | 34 |
| QubitGKSLShot GPS | 388 |
| QubitGKSLBoson N=4 n_total_qubits | 42 |
| QubitGKSLBoson N=4 GPS | 484 |
| Qubit DM/Shot n_total一致 | DM == Shot |

## 3. 修正内容

### 3.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `quantum_dynamics_gksl_comparison.ipynb` Cell 32 | レジスタ情報表示追加（Issue U修正） |
| `quantum_dynamics_gksl_comparison.ipynb` Cell 36 | n_ancilla, n_total_qubits, Trace conservation追加（Issue V修正） |
| `run_iteration20_verification.py` | 34項目の検証スクリプト（新規、Issue W/X修正） |

### 3.2 変更不要

- シミュレータのPythonファイル: コードに問題なし
- `quantum_dynamics_complete_comparison.ipynb`: GKSLシミュレータを直接使用しない
- `run_iteration19_verification.py`: 旧バージョンとして保持

## 4. 検証結果

`tutorials/run_iteration20_verification.py` で34項目を検証：

### Section 1: n_total正確性（8項目）

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

### Section 2: クロスチェック（4項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 9 | 非ボソンアンシラ一致 (26=26) | PASS |
| 10 | ボソンアンシラ一致 (12=12) | PASS |
| 11 | qudit DM/Shot n_total一致 (30=30) | PASS |
| 12 | qubit DM/Shot n_total一致 (34=34) | PASS |

### Section 3: GPSレジスタ回帰テスト（8項目）

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

### Section 4: ノートブックセル検査（10項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 21 | Cell 6 n_total_qudits表示 | PASS |
| 22 | Cell 8 n_total_qubits表示 | PASS |
| 23 | Cell 12 n_total_qudits表示 | PASS |
| 24 | Cell 14 n_total_qubits表示 | PASS |
| 25 | Cell 32 n_total_qudits表示 | PASS |
| 26 | Cell 32 n_system_qudits表示 | PASS |
| 27 | Cell 32 n_ancilla_qudits表示 | PASS |
| 28 | Cell 36 n_total_qubits表示 | PASS |
| 29 | Cell 36 n_ancilla表示 | PASS |
| 30 | Cell 36 Trace conservation表示 | PASS |

### Section 5: 結果辞書完全性（4項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 31 | qudit DM result n_total = 30 | PASS |
| 32 | qudit boson result n_total = 16 | PASS |
| 33 | qudit shot result n_total = 30 | PASS |
| 34 | qubit shot result n_total = 34 | PASS |

## 5. 残課題

### 5.1 qubit回路シミュレータの抽象化レベル不整合（継続）

Iteration 17/18/19から継続の課題。qubit回路シミュレータの`gates_per_step`（複合
ゲート単位66）とqubit DM/shotの`estimated_gates_per_step`（基本ゲート推定388）
は異なる抽象化レベルの測定値。

### 5.2 ユーザーによるノートブック再実行

修正後のノートブックをユーザーがローカルで実行し、Cell 32で
`System qutrits: 4 / Ancilla qudits: 26 / Total qudits: 30`、
Cell 36で`Ancilla qubits: 26 / Total qubits: 34 / Trace conservation`
が正しく表示されることを確認する必要がある。

## 6. 次のステップ

1. ユーザーがiteration 20の検証スクリプトをローカルで実行し結果を確認
2. ノートブックを再実行してCell 32, 36の表示を目視確認
3. qubit回路シミュレータの抽象化レベル統一を検討（次イテレーション候補）
