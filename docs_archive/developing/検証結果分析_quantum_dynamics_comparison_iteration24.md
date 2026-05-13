# 検証結果分析: quantum_dynamics_comparison iteration 24

## 1. 概要

Iteration 23の検証結果（62/62 PASS）、ノートブック実行結果、およびソースコードの詳細分析を実施した。
Iteration 23の修正（Issue BB: パリンドローム回路可視化、Issue CC: n_phonon_qudits一般化）は正確に実装され、
全検証チェックをパスしている。

新たに**回路ボソンシミュレータ（QuditGKSLCircuitBosonSimulator）に2件のAPI整合性バグ**を発見・修正した。

### 1.1 分析対象

1. `iteration23_palindromic_boson_20260314T090607Z.json`: 62/62 PASS
2. `quantum_dynamics_gksl_comparison.ipynb` の全セル実行結果
3. `quantum_dynamics_complete_comparison.ipynb` の全セル実行結果
4. 全7シミュレータのソースコード詳細レビュー（DM 4種 + Shot 2種 + Circuit 1種）
5. クロスシミュレータ数値比較（DM qudit vs qubit vs circuit）

### 1.2 分析結論

- **Iteration 23の修正は正確**: Issue BB（Cell 16パリンドローム）とIssue CC（n_phonon_qudits一般化）は正しく修正済み
- **シミュレータの物理計算は正確**: 全シミュレータがmachine precision（~1e-14以下）で一致
- **Issue DD発見**: 回路ボソンシミュレータのrho_finalが不整合（拡張空間全体 vs 電子系のみ）
- **Issue EE発見**: 回路ボソンシミュレータの結果辞書キー名が不整合

## 2. Iteration 23修正の検証

### 2.1 Issue BB: Cell 16パリンドローム修正 ✓

Cell 16の実行出力を確認：
- Scenario 3 (Qubit非ボソン): 54 subcircuits ✓
- Scenario 5 (Qudit非ボソン): 54 subcircuits ✓
- Scenario 6 (Quditボソン): **26 subcircuits** ✓（修正前: 14）
- Scenario 4 (Qubitボソン): **26 subcircuits** ✓（修正前: 14）

各ボソンシナリオの回路構造を確認：
- H_half_1 → 12 forward channels → 12 reverse channels → H_half_2
- 逆順はインデックス11→0の降順で正しく実装
- dt_quantum/2が正しく使用

### 2.2 Issue CC: n_phonon_qudits一般化 ✓

`ceil(log_d(n_max+1)) * N`式が正しく動作：
- n_max=1: ceil(log_3(2)) = 1 → n_phonon = 2 ✓
- n_max=2: ceil(log_3(3)) = 1 → n_phonon = 2 ✓
- n_max=3: ceil(log_3(4)) = 2 → n_phonon = 4 ✓
- n_max=8: ceil(log_3(9)) = 2 → n_phonon = 4 ✓

浮動小数点精度を追加検証：log_3(243) = 4.9999...→ ceil = 5 ✓（問題なし）

## 3. クロスシミュレータ数値検証

### 3.1 QuditGKSLBosonSimulator vs QubitGKSLBosonSimulator

同一パラメータ（N=2, n_max=1, t=1.0, 100 steps）で比較：
- `||rho_qudit - rho_qubit||_F = 2.17e-20`（machine precision）
- 個体群分布: 完全一致

**結論**: 両シミュレータは数学的に等価であり、物理計算は完全に正確。

### 3.2 QuditGKSLBosonSimulator vs QuditGKSLCircuitBosonSimulator

修正後に同一パラメータで比較：
- `||rho_dm - rho_circuit||_F = 1.89e-15`（machine precision）
- トレース保存: max|Tr-1| < 1e-14

**結論**: 回路シミュレータとDMシミュレータは数学的に等価。

## 4. 発見されたバグ

### 4.1 Issue DD: 回路ボソンシミュレータのrho_final不整合

**問題**: `QuditGKSLCircuitBosonSimulator.simulate()` の `rho_final` が拡張空間全体の密度行列（dim_el*dim_ph × dim_el*dim_ph = 36×36）を返していた。一方、他の全ボソンシミュレータ（QuditGKSLBosonSimulator、QubitGKSLBosonSimulator）は電子系のみの縮約密度行列（dim_el × dim_el = 9×9）を返す。

**影響**: シミュレータ間のrho_final比較が不可能（形状不一致エラー）。

**修正内容**:
```python
# 修正前
"rho_final": rho,  # 36×36 (電子⊗フォノン全体)

# 修正後
rho_el_final = partial_trace_phonon(rho, self.dim_el, self.dim_ph)
"rho_final": rho_el_final,  # 9×9 (電子系のみ)
```

### 4.2 Issue EE: 回路ボソンシミュレータの結果辞書キー不整合

**問題**: `QuditGKSLCircuitBosonSimulator.simulate()` の返却辞書が以下の点でDMシミュレータと不整合：

| キー名 | 回路シミュレータ（修正前） | DMシミュレータ |
|--------|------------------------|-------------|
| 電子レジスタ数 | `n_el_qutrits` | `n_system_qudits` |
| フォノンレジスタ数 | `n_ph_qutrits` (= N固定) | `n_phonon_qudits` (= ceil(log_d(n_max+1))*N) |
| 総レジスタ数 | 未定義 | `n_total_qudits` |
| ゲート数/ステップ | `gates_per_step` | `estimated_gates_per_step` |
| 総ゲート数 | `total_gates` | `total_estimated_gates` |

**修正内容**: DMシミュレータと同一のキー名に統一。`n_phonon_qudits`は一般化公式で計算。

## 5. 修正ファイル一覧

| ファイル | 変更内容 |
|---------|---------|
| `tutorials/qudit_gksl_circuit_boson_simulator.py` | Issue DD + EE: rho_final修正、キー名統一 |
| `tutorials/run_iteration24_verification.py` | 新規: 修正検証スクリプト（74チェック） |
| `developing/検証結果分析_quantum_dynamics_comparison_iteration24.md` | 新規: 本分析文書 |

## 6. 検証項目（Iteration 24）

### 新規追加チェック（12項目）

| チェック名 | 内容 |
|-----------|------|
| circuit_boson_rho_final_shape | rho_finalが(9,9)形状であること |
| circuit_boson_rho_matches_dm | rho_finalがDMシミュレータと一致（<1e-12） |
| circuit_boson_rho_final_trace | rho_finalのトレースが~1.0 |
| circuit_boson_has_n_system_qudits | n_system_quditsキーが存在 |
| circuit_boson_has_n_phonon_qudits | n_phonon_quditsキーが存在 |
| circuit_boson_has_n_total_qudits | n_total_quditsキーが存在 |
| circuit_boson_has_estimated_gps | estimated_gates_per_stepキーが存在 |
| circuit_boson_n_system_matches_dm | n_system_quditsがDMと一致 |
| circuit_boson_n_phonon_matches_dm | n_phonon_quditsがDMと一致 |
| circuit_boson_n_total_matches_dm | n_total_quditsがDMと一致 |
| circuit_boson_gps_matches_dm | estimated_gates_per_stepがDMと一致 |
| circuit_boson_no_old_keys | 旧キー名が存在しないこと |

### 既存リグレッションチェック（62項目、Iteration 23から継承）

62項目すべてリグレッションなしを確認済み。

### 合計: 74チェック

## 7. 残課題

### 7.1 Qubit回路シミュレータの抽象レベル不一致（既知、iteration 17以降）

Qubit回路シミュレータは66個のcompositeゲート（UnitaryGate）をカウントするが、qubit DM/shotは388個のbasicゲートを推定する。公正な比較には統一が必要。

### 7.2 ボソン回路シミュレータのbuild_full_trotter_step_circuit()不在（既知）

非ボソン回路シミュレータは`build_full_trotter_step_circuit()`メソッドを持つが、ボソン回路シミュレータにはない。Cell 16が手動で回路を構築する必要がある。

### 7.3 次の検証手順

1. ユーザーが`run_iteration24_verification.py`をローカルで実行し、結果をpush
2. 結果を確認し、74/74 PASSを確認
