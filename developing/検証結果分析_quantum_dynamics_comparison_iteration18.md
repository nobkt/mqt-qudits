# 検証結果分析: Quantum Dynamics Comparison - Iteration 18

## 1. 背景

### 1.1 Iteration 17の結果

Iteration 17ではBug P-Q（ボソンシミュレータの`estimated_gates_per_step`に
H_eph電子-フォノン結合ゲートが含まれていない問題）が修正され、
16/16チェックが全てPASSした。

### 1.2 Iteration 17の残課題

Iteration 17の分析ドキュメント（Section 7）で、以下の3項目が
次のステップとして挙げられていた：

1. ユーザーが検証スクリプトをローカルで実行して結果を確認 → **完了**
2. ノートブックでボソンシミュレーション結果のゲートカウント表示を確認 → **問題発見**
3. qubit回路シミュレータの`gates_per_step`（複合ゲート単位66）と
   qubit DM/shotの`estimated_gates_per_step`（基本ゲート単位388）の
   不整合について検討 → **次イテレーションで対応**

## 2. 発見された問題の詳細

### 2.1 ノートブックのゲートカウント表示欠落

**ファイル**: `quantum_dynamics_gksl_comparison.ipynb`

**問題**: 6つの量子シミュレータセルのうち、`estimated_gates_per_step`を
表示していたのはCell 6（QuditGKSLSimulator）のみ。残り5つのセルでは、
シミュレータの結果辞書にゲートカウントが含まれているにもかかわらず、
表示されていなかった。

| セル | シミュレータ | 修正前 | 修正後 |
|------|-------------|--------|--------|
| Cell 6 | QuditGKSLSimulator | ✓ 表示あり (66) | ✓ 変更なし |
| Cell 8 | QubitGKSLSimulator | ✗ 表示なし | ✓ 表示追加 (388) |
| Cell 12 | QuditGKSLBosonSimulator | ✗ 表示なし | ✓ 表示追加 (38) |
| Cell 14 | QubitGKSLBosonSimulator | ✗ 表示なし | ✓ 表示追加 (220) |
| Cell 32 | QuditGKSLShotSimulator | ✗ 表示なし | ✓ 表示追加 (66) |
| Cell 36 | QubitGKSLShotSimulator | ✗ 表示なし | ✓ 表示追加 (388) |

### 2.2 ボソンセルのレジスタ情報欠落

ボソンシミュレータのセル（Cell 12, 14）では、経過時間とトレース保存の
2行のみが表示され、レジスタ情報（qudit/qubit数）が欠落していた。
非ボソンセル（Cell 6, 8）では表示されていたため、情報の非対称性があった。

**Cell 12に追加した情報**:
- System qutrits（電子系キュートリット数）
- Phonon qutrits（フォノン系キュートリット数）
- Ancilla qudits（アンシラキュディット数）

**Cell 14に追加した情報**:
- Electronic qubits（電子系キュービット数）
- Phonon qubits（フォノン系キュービット数）
- Ancilla qubits（アンシラキュービット数）
- Total qubits（総キュービット数）

## 3. 修正内容

### 3.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `quantum_dynamics_gksl_comparison.ipynb` | Cell 8, 12, 14, 32, 36にゲートカウント表示を追加 |
| `run_iteration18_verification.py` | 31項目の検証スクリプト（新規） |

### 3.2 変更不要

- シミュレータのPythonコード: 変更なし（`estimated_gates_per_step`は
  全てのシミュレータで正しく計算・返却されていた）
- `quantum_dynamics_complete_comparison.ipynb`: GKSLシミュレータを
  直接使用しないため変更不要

## 4. ゲートカウント一覧（確認済み値）

### 4.1 非ボソンモデル（N=4）

| シミュレータ | GPS | 内訳 |
|-------------|-----|------|
| qudit DM | 66 | 2×(4+3) + 26×2 |
| qubit DM | 388 | 2×(8+30) + 26×12 |
| qudit shot | 66 | 同上 |
| qubit shot | 388 | 同上 |

### 4.2 ボソンモデル（N=2, ノートブックのデフォルト）

| シミュレータ | GPS | 内訳 |
|-------------|-----|------|
| qudit boson | **38** | 2×(2+2+1+2) + 12×2 |
| qubit boson | **220** | 2×(4+4+10+20) + 12×12 |

### 4.3 ボソンモデル（N=4, 回帰テスト用）

| シミュレータ | GPS | 内訳 |
|-------------|-----|------|
| qudit boson | **82** | 2×(4+4+3+4) + 26×2 |
| qubit boson | **484** | 2×(8+8+30+40) + 26×12 |

### 4.4 qubit/qudit比率

| モデル | N | 比率 |
|--------|---|------|
| 非ボソン | 4 | 5.88× (388/66) |
| ボソン | 2 | 5.79× (220/38) |
| ボソン | 4 | 5.90× (484/82) |

## 5. 検証結果

`tutorials/run_iteration18_verification.py` で31項目を検証：

### Section 1: ノートブックセル検査（13項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 1 | Cell 6 gates_per_step表示 | PASS |
| 2 | Cell 8 gates_per_step表示 | PASS |
| 3 | Cell 8 total_gates表示 | PASS |
| 4 | Cell 12 gates_per_step表示 | PASS |
| 5 | Cell 12 system_qudits表示 | PASS |
| 6 | Cell 12 phonon_qudits表示 | PASS |
| 7 | Cell 12 ancilla表示 | PASS |
| 8 | Cell 14 gates_per_step表示 | PASS |
| 9 | Cell 14 el_qubits表示 | PASS |
| 10 | Cell 14 ph_qubits表示 | PASS |
| 11 | Cell 14 total_qubits表示 | PASS |
| 12 | Cell 32 gates_per_step表示 | PASS |
| 13 | Cell 36 gates_per_step表示 | PASS |

### Section 2: ゲートカウント値検証（6項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 14 | qudit DM GPS = 66 | PASS |
| 15 | qubit DM GPS = 388 | PASS |
| 16 | qudit boson GPS = 38 (N=2) | PASS |
| 17 | qubit boson GPS = 220 (N=2) | PASS |
| 18 | qudit shot GPS = 66 | PASS |
| 19 | qubit shot GPS = 388 | PASS |

### Section 3: クロスチェック（4項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 20 | qudit DM == qudit shot | PASS |
| 21 | qubit DM == qubit shot | PASS |
| 22 | ボソンN=2 比率 = 5.79× | PASS |
| 23 | 非ボソンN=4 比率 = 5.88× | PASS |

### Section 4: ボソンレジスタ数検証（6項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 24 | qudit boson system_qudits = 2 | PASS |
| 25 | qudit boson phonon_qudits = 2 | PASS |
| 26 | qudit boson ancilla = 12 | PASS |
| 27 | qubit boson el_qubits = 4 | PASS |
| 28 | qubit boson ph_qubits = 4 | PASS |
| 29 | qubit boson total = 20 | PASS |

### Section 5: 回帰テスト（2項目）

| # | チェック項目 | 結果 |
|---|-------------|------|
| 30 | qudit boson N=4 GPS = 82 | PASS |
| 31 | qubit boson N=4 GPS = 484 | PASS |

## 6. 残課題

### 6.1 qubit回路シミュレータの抽象化レベル不整合（継続）

Iteration 17で指摘されたqubit回路シミュレータの`gates_per_step`（複合ゲート
単位66）とqubit DM/shotの`estimated_gates_per_step`（基本ゲート推定388）の
不整合は、今回のスコープ外とした。

**問題の本質**:
- qudit回路シミュレータ: 66ゲート（ネイティブqudit操作 = 複合ゲート）
- qubit回路シミュレータ: 66ゲート（UnitaryGate複合ゲート）
- qubit DM/shotシミュレータ: 388ゲート（CX+1qゲート推定）

quditでは複合ゲート = ネイティブゲートのため一致するが、qubitでは
UnitaryGateを基本ゲートに分解する必要があり、66と388は異なる抽象化レベルの
測定値である。公平な比較のためには統一が望ましい。

### 6.2 ユーザーによるノートブック再実行

修正後のノートブックをユーザーがローカルで実行し、各セルで
ゲートカウントが正しく表示されることを確認する必要がある。

## 7. 次のステップ

1. ユーザーがiteration 18の検証スクリプトをローカルで実行し結果を確認
2. ノートブックを再実行してゲートカウント表示を目視確認
3. qubit回路シミュレータの抽象化レベル統一を検討（次イテレーション候補）
