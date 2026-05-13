# 検証結果分析: Quantum Dynamics Comparison - Iteration 19

## 1. 背景

### 1.1 Iteration 18の結果

Iteration 18では6つのシミュレータセルすべてに`estimated_gates_per_step`の
表示を追加し、ボソンセルにはレジスタ情報も追加した。
31/31チェックがPASSした。

### 1.2 Iteration 18の残課題

Iteration 18の分析ドキュメント（Section 6）で以下が挙げられていた：

1. qubit回路シミュレータの抽象化レベル不整合 → **次イテレーション以降**
2. ユーザーによるノートブック再実行確認 → **完了**

## 2. 発見された問題の詳細

### 2.1 Bug R: `n_total_qudits`がアンシラを含んでいない

**ファイル**: `qudit_gksl_boson_simulator.py`（行63）

**問題**: `n_total_qudits`が`n_system_qudits + n_phonon_qudits`と定義されており、
**Stinespringアンシラキュディットを含んでいなかった**。

一方、対応するqubitシミュレータ（`qubit_gksl_boson_simulator.py`）では
`n_total_qubits = n_el_qubits + n_ph_qubits + n_ancilla`と定義されており、
**アンシラを含んでいた**。

| シミュレータ | 修正前 | 修正後 |
|-------------|--------|--------|
| qudit_gksl_boson_simulator | n_total = sys + ph = **4** (N=2) | n_total = sys + ph + anc = **16** (N=2) |
| qubit_gksl_boson_simulator | n_total = el + ph + anc = **20** (N=2) | 変更なし |

**影響**: `result["n_total_qudits"]`を使用して比較を行うコードがあった場合、
quditの総レジスタ数が過小評価され、不公平な比較となる。

### 2.2 Issue S: qudit非ボソン・ショットシミュレータに`n_total_qudits`が未定義

| シミュレータ | 修正前 | 修正後 |
|-------------|--------|--------|
| qudit_gksl_simulator | `n_total_qudits`未定義 | 追加: sys + anc = **30** (N=4) |
| qudit_gksl_shot_simulator | `n_total_qudits`未定義 | 追加: sys + anc = **30** (N=4) |
| qubit_gksl_simulator | `n_total_qubits`=**34** (既存) | 変更なし |
| qubit_gksl_shot_simulator | `n_total_qubits`定義済み (既存) | 変更なし |

### 2.3 Issue T: ノートブック表示の非対称性

qubitシミュレータのセル（Cell 8, 14）では「Total qubits」を表示していたが、
対応するquditシミュレータのセル（Cell 6, 12）では「Total qudits」を
表示していなかった。

| セル | シミュレータ | 修正前 | 修正後 |
|------|-------------|--------|--------|
| Cell 6 | QuditGKSLSimulator | Total未表示 | Total qudits: 30 を追加 |
| Cell 8 | QubitGKSLSimulator | Total qubits: 34 | 変更なし |
| Cell 12 | QuditGKSLBosonSimulator | Total未表示 | Total qudits: 16 を追加 |
| Cell 14 | QubitGKSLBosonSimulator | Total qubits: 20 | 変更なし |

## 3. 修正内容

### 3.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `qudit_gksl_boson_simulator.py` | `n_total_qudits`にアンシラを追加（Bug R修正） |
| `qudit_gksl_simulator.py` | `n_total_qudits`を新規追加（Issue S修正） |
| `qudit_gksl_shot_simulator.py` | `n_total_qudits`を新規追加（Issue S修正） |
| `quantum_dynamics_gksl_comparison.ipynb` | Cell 6, 12に「Total qudits」表示追加（Issue T修正） |
| `run_iteration19_verification.py` | 30項目の検証スクリプト（新規） |

### 3.2 変更不要

- qubit系シミュレータ: 既に正しく`n_total_qubits`にアンシラを含んでいた
- ゲートカウント計算式: 変更なし（`estimated_gates_per_step`は正しい）
- `quantum_dynamics_complete_comparison.ipynb`: GKSLシミュレータを直接使用しない

## 4. n_total一覧（修正後の確認済み値）

### 4.1 非ボソンモデル（N=4）

| シミュレータ | sys | anc | n_total |
|-------------|-----|-----|---------|
| qudit DM | 4 qutrits | 26 qudits | **30 qudits** |
| qubit DM | 8 qubits | 26 qubits | **34 qubits** |
| qudit shot | 4 qutrits | 26 qudits | **30 qudits** |
| qubit shot | 8 qubits | 26 qubits | **34 qubits** |

### 4.2 ボソンモデル（N=2）

| シミュレータ | sys | ph | anc | n_total |
|-------------|-----|-----|-----|---------|
| qudit boson | 2 qutrits | 2 qutrits | 12 qudits | **16 qudits** |
| qubit boson | 4 qubits | 4 qubits | 12 qubits | **20 qubits** |

### 4.3 ボソンモデル（N=4）

| シミュレータ | sys | ph | anc | n_total |
|-------------|-----|-----|-----|---------|
| qudit boson | 4 qutrits | 4 qutrits | 26 qudits | **34 qudits** |
| qubit boson | 8 qubits | 8 qubits | 26 qubits | **42 qubits** |

### 4.4 qubit/qudit比率（n_total）

| モデル | N | qudit | qubit | 比率 |
|--------|---|-------|-------|------|
| 非ボソン | 4 | 30 | 34 | 1.13× |
| ボソン | 2 | 16 | 20 | 1.25× |
| ボソン | 4 | 34 | 42 | 1.24× |

**注**: アンシラ数は同一（Lindblad演算子の数に依存）。差はシステム+フォノン
レジスタの次元差（qudit: d=3で1レジスタ vs qubit: d=4で2レジスタ/分子）に起因。

## 5. 検証結果

`tutorials/run_iteration19_verification.py` で30項目を検証：

### Section 1: シミュレータクラスのn_total検査（11項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 1 | QuditGKSLSimulatorがn_total_quditsを持つ | PASS |
| 2 | n_total_quditsにアンシラ含む | PASS |
| 3 | qudit DM n_total = 30 | PASS |
| 4 | qubit DM n_total = 34 | PASS |
| 5 | qudit boson n_totalにアンシラ含む | PASS |
| 6 | qudit boson n_total = 16 | PASS |
| 7 | qubit boson n_totalにアンシラ含む | PASS |
| 8 | qubit boson n_total = 20 | PASS |
| 9 | QuditGKSLShotSimulatorがn_total_quditsを持つ | PASS |
| 10 | qudit shot n_totalにアンシラ含む | PASS |
| 11 | qudit shot n_total = 30 | PASS |

### Section 2: 結果辞書のn_total検査（6項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 12 | qudit DM結果にn_total_quditsあり | PASS |
| 13 | qudit DM結果のn_total = 30 | PASS |
| 14 | qudit boson結果にn_total_quditsあり | PASS |
| 15 | qudit boson結果のn_total = 16 | PASS |
| 16 | qudit shot結果にn_total_quditsあり | PASS |
| 17 | qudit shot結果のn_total = 30 | PASS |

### Section 3: クロスチェック（3項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 18 | 非ボソン: qudit/qubitアンシラ数一致 (26=26) | PASS |
| 19 | ボソン: qudit/qubitアンシラ数一致 (12=12) | PASS |
| 20 | qudit DM/shotのn_total一致 (30=30) | PASS |

### Section 4: ノートブックセル検査（4項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 21 | Cell 6にn_total_qudits表示あり | PASS |
| 22 | Cell 8にn_total_qubits表示あり | PASS |
| 23 | Cell 12にn_total_qudits表示あり | PASS |
| 24 | Cell 14にn_total_qubits表示あり | PASS |

### Section 5: ゲートカウント回帰テスト（6項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 25 | qudit DM GPS = 66 | PASS |
| 26 | qubit DM GPS = 388 | PASS |
| 27 | qudit boson GPS = 38 | PASS |
| 28 | qubit boson GPS = 220 | PASS |
| 29 | qudit shot GPS = 66 | PASS |
| 30 | qudit boson N=4 n_total = 34 | PASS |

## 6. 残課題

### 6.1 qubit回路シミュレータの抽象化レベル不整合（継続）

Iteration 17/18から継続の課題。qubit回路シミュレータの`gates_per_step`（複合
ゲート単位66）とqubit DM/shotの`estimated_gates_per_step`（基本ゲート推定
388）は異なる抽象化レベルの測定値。

### 6.2 ユーザーによるノートブック再実行

修正後のノートブックをユーザーがローカルで実行し、各セルで
n_total_quditsが正しく表示されることを確認する必要がある。

## 7. 次のステップ

1. ユーザーがiteration 19の検証スクリプトをローカルで実行し結果を確認
2. ノートブックを再実行してn_total_qudits表示を目視確認
3. qubit回路シミュレータの抽象化レベル統一を検討（次イテレーション候補）
