# quditアンシラ修正 作業ログ iteration1

## 作業日時
2026-03-11T07:20:00Z

## 課題
`quantum_dynamics_gksl_comparison.ipynb`のqudit表現は、qudit量子演算型の量子コンピュータ（例：qudit-boson型イオントラップ量子コンピュータ）で実行することを想定しているため、Stinespring拡張の補助ビット（アンシラ）はqubit(d=2)ではなくqudit(d=3)を使用する必要がある。

## 分析結果

### 修正前の状態
全quditシミュレータ（QuditGKSLSimulator, QuditGKSLShotSimulator, QuditGKSLNoisySimulator, QuditGKSLBosonSimulator, QuditGKSLCircuitSimulator, QuditGKSLCircuitBosonSimulator）において、Stinespring拡張のアンシラが**d=2（qubit）**として実装されていた。

具体的には：
- `stinespring_utils.py`: `stinespring_unitary_from_lindblad(L, dt)` が常に `2*d_sys × 2*d_sys` のユニタリを生成
- `apply_stinespring_to_density_matrix(rho, U)` が常に2レベルアンシラ `|0><0|` (2×2) を使用
- 各シミュレータに `n_ancilla_qubits` 属性があり、qubitとして表現
- 回路シミュレータで `QuantumCircuit(2, [d, 2], 0)` とアンシラがd=2で登録

### qubit vs qudit アンシラの物理的違い
- **qubitアンシラ (d=2)**: Stinespring unitaryは `2*d_sys × 2*d_sys` 次元。Kraus演算子K₀, K₁の2つ。
- **quditアンシラ (d=3)**: Stinespring unitaryは `3*d_sys × 3*d_sys` 次元。Kraus演算子K₀, K₁, K₂の3つ。ただしK₂=0（生成子Gの第3ブロックが0であるため）。
- **物理結果は完全に同一**: K₀, K₁はd_anc=2とd_anc=3で完全に一致し、K₂=0なので量子チャネルは数学的に同一。

### 修正の正当性
qudit量子コンピュータでは全てのレジスタがd次元のquditである。Stinespring拡張のアンシラも物理的にはd次元のquditレジスタを使用する（第3レベル|2⟩は未使用のまま）。これはqudit量子コンピュータのハードウェア構成を正確に反映している。

## 修正内容

### 1. `stinespring_utils.py`
- `stinespring_unitary_from_lindblad(L, dt, d_anc=2)`: `d_anc`パラメータ追加
- `apply_stinespring_to_density_matrix(rho, U, d_anc=2)`: `d_anc`パラメータ追加
- デフォルトd_anc=2で後方互換性維持（qubitシミュレータはそのまま動作）

### 2. `qudit_gksl_simulator.py`
- `self.d_anc = params.d` 追加
- `n_ancilla_qubits` → `n_ancilla_qudits` 名前変更
- `_precompute_unitaries`: `d_anc=self.d_anc` で呼び出し
- `_trotter_step`: `d_anc=self.d_anc` で呼び出し
- 結果辞書に `d_anc` キー追加

### 3. `qudit_gksl_shot_simulator.py`
- `self.d_anc = params.d` 追加
- `n_ancilla_qubits` → `n_ancilla_qudits`
- `_apply_stinespring_with_measurement`: d_ancレベルのアンシラ測定に拡張
  - `d_anc`個の分岐でBorn確率を計算し、ランダムに射影
- 結果辞書に `d_anc`, `n_ancilla_qudits` キー追加

### 4. `qudit_gksl_noisy_simulator.py`
- 親クラスからd_anc継承
- `_trotter_step`: `d_anc=self.d_anc` で呼び出し

### 5. `qudit_gksl_boson_simulator.py`
- `self.d_anc = params.d` 追加
- `n_ancilla_qubits` → `n_ancilla_qudits`
- Stinespring呼び出しに `d_anc` パラメータ追加

### 6. `qudit_gksl_circuit_simulator.py`
- `self.d_anc = params.d` 追加
- `_build_local_stinespring_unitary`: `d_anc*d_local` 次元の生成子を構築
- 回路構築: `QuantumCircuit(2, [d, d_anc], 0)` でアンシラをd_ancレベルに
- `_extract_kraus_from_local_stinespring`: d_anc個のKraus演算子をリストで返す
- `_apply_single_site_channel`, `_apply_pair_channel`: Kraus演算子リストを受け取る形式に変更

### 7. `qudit_gksl_circuit_boson_simulator.py`
- 同上の変更をボソン版にも適用

### 8. ノートブック・テスト
- `quantum_dynamics_gksl_comparison.ipynb`: 「ancilla qubit」→「ancilla qudit(d=3)」
- `test_gksl_simulators.py`: `n_ancilla_qubits` → `n_ancilla_qudits` + d_ancチェック追加

## 検証結果
32/32チェック全てパス。

主要検証項目：
1. ✓ Stinespring unitary次元: d_anc=3で9×9（3×3系の場合）、d_anc=2で6×6
2. ✓ K₀, K₁はd_anc=2とd_anc=3で完全一致（||差||_F = 0.0）
3. ✓ K₂ = 0（||K₂||_F = 0.0）
4. ✓ トレース保存: ||ΣK_k†K_k - I||_F = 0.0
5. ✓ 密度行列チャネル出力: d_anc=2とd_anc=3で完全一致
6. ✓ 全シミュレータに`n_ancilla_qudits`属性あり（`n_ancilla_qubits`なし）
7. ✓ 全シミュレータの`d_anc == 3`
8. ✓ qubitシミュレータは変更なし（デフォルトd_anc=2を使用）
9. ✓ 81次元系のフルStinespring: 243×243ユニタリ、ユニタリ性残差4.75e-16
10. ✓ 5ステップシミュレーション: トレース保存誤差6.66e-16

## 変更されないもの
- qubitシミュレータ（QubitGKSLSimulator等）: qubit量子コンピュータ向けなのでd_anc=2のまま
- 物理結果: d_anc=2とd_anc=3で量子チャネルは数学的に完全同一
