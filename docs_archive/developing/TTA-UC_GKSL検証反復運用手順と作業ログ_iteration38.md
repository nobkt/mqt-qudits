# TTA-UC GKSL 検証反復運用手順と作業ログ: Iteration 38

## 作業概要

- 日時: 2026-03-03
- 目的: Iteration 37 の検証結果・ノートブック・コード全体の詳細分析、メタデータバグの修正
- 前提: Iteration 37 で回文順序一貫性が確認され、全18チェック PASS

## 問題の特定

Iteration 37 の検証結果、ノートブック全40セル、シミュレーターコード全体、数学的理論を精査した結果、以下の問題を特定:

### 問題1: QuditGKSLSimulator の n_ancilla_qubits ハードコード（重要度: 低）
- **場所**: `tutorials/qudit_gksl_simulator.py` 行67-68
- **内容**: `self.n_ancilla_qubits = 26` がハードコード。2分子系では12であるべき
- **影響**: メタデータのみ。シミュレーション計算には影響なし

### 問題2: QuditGKSLSimulator のゲート数推定ハードコード（重要度: 低）
- **場所**: `tutorials/qudit_gksl_simulator.py` 行200
- **内容**: `gates_per_step = 4 + 3 + 26 * 2` がN=4固定。N≠4で誤った値
- **影響**: メタデータのみ

### 問題3: QuditGKSLShotSimulator のゲート数推定ハードコード（重要度: 低）
- **場所**: `tutorials/qudit_gksl_shot_simulator.py` 行287
- **内容**: 問題2と同一。`n_ancilla_qubits` は動的だがゲート数はハードコード
- **影響**: メタデータのみ

### 問題4: QubitGKSLShotSimulator にゲート数推定が欠如（重要度: 低）
- **場所**: `tutorials/qubit_gksl_shot_simulator.py`
- **内容**: simulate()の返却値にゲート数推定が含まれていない
- **影響**: 他シミュレーターとの不整合

## 修正内容

### ファイル1: `tutorials/qudit_gksl_simulator.py`
1. `n_ancilla_qubits` を `len(self.lindblad_ops)` に変更
2. `gates_per_step` を動的に計算:
   ```python
   n_lindblad = len(self.lindblad_ops)
   gates_per_step = self.n_system_qudits + len(self.params.neighbors) + n_lindblad * 2
   ```

### ファイル2: `tutorials/qudit_gksl_shot_simulator.py`
1. `gates_per_step` を動的に計算（`n_ancilla_qubits` は既に正しい）

### ファイル3: `tutorials/qubit_gksl_shot_simulator.py`
1. simulate() にゲート数推定を追加

## テスト

既存テスト `tutorials/test_gksl_simulators.py` の以下のテストが影響:
- `TestQuditGKSLSimulator.test_initialization`: `assert sim.n_ancilla_qubits == 26`
  → デフォルト4分子系では依然として26なので変更不要

新規検証: Iteration 38 検証スクリプトで2分子系と4分子系のアンシラ数・ゲート数を比較

## 次のステップ

1. ユーザーがローカルで iteration 38 検証スクリプトを実行
2. 結果を確認し、全チェック PASS を確認
