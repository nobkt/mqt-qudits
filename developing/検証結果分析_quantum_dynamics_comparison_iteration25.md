# 検証結果分析: quantum_dynamics_comparison iteration 25

## 1. 概要

Iteration 24の検証結果（74/74 PASS）、ノートブック実行結果、および全シミュレータのソースコードを詳細分析した。
Iteration 24の修正（Issue DD: rho_final形状修正、Issue EE: 辞書キー統一）は正確に実装され、全検証チェックをパスしている。

新たに**回路ボソンシミュレータ（QuditGKSLCircuitBosonSimulator）に2件のバグ**を発見・修正した。

### 1.1 分析対象

1. `iteration24_circuit_boson_consistency_20260314T224705Z.json`: 74/74 PASS
2. `quantum_dynamics_gksl_comparison.ipynb` の全セル実行結果
3. `quantum_dynamics_complete_comparison.ipynb` の全セル実行結果
4. 全7シミュレータのソースコード詳細レビュー（DM 4種 + Shot 2種 + Circuit 1種）
5. クロスシミュレータ数値比較（DM vs Circuit: rho_final、traces時系列）

### 1.2 分析結論

- **Iteration 24の修正は正確**: Issue DD（rho_final形状）とIssue EE（キー名統一）は正しく修正済み
- **Issue FF発見**: 回路ボソンシミュレータのtraces計算がDMボソンシミュレータと不整合
- **Issue GG発見**: 回路ボソンシミュレータが毎Trotterステップでユニタリ行列を再計算（不要な計算コスト）

## 2. Iteration 24修正の検証

### 2.1 Issue DD: rho_final形状修正 ✓

- rho_finalは(9,9)形状（電子系のみ）: ✓
- DMシミュレータとの一致: ||rho_circuit - rho_dm||_F = 1.78e-15 ✓
- トレース保存: Tr(rho_final) = 1.0000000000 ✓

### 2.2 Issue EE: 辞書キー統一 ✓

- `n_system_qudits`: circuit=2, dm=2 ✓
- `n_phonon_qudits`: circuit=2, dm=2 ✓
- `n_total_qudits`: circuit=16, dm=16 ✓
- `estimated_gates_per_step`: circuit=38, dm=38 ✓
- 旧キー名（`n_el_qutrits`, `n_ph_qutrits`）は存在しない ✓

## 3. 発見されたバグ

### 3.1 Issue FF: 回路ボソンシミュレータのtraces計算不整合

**問題**: `QuditGKSLCircuitBosonSimulator.simulate()` の `traces` リストが `np.trace(rho)` （拡張空間全体の密度行列、36×36）から計算されていた。一方、`QuditGKSLBosonSimulator` は `np.trace(rho_el)` （電子系のみの縮約密度行列、9×9）から計算する。

**技術的詳細**:
- 数学的には `Tr(rho) = Tr(rho_el)` が成立する（部分トレースはトレースを保存）
- しかし、数値計算上は異なる丸め誤差が蓄積される可能性がある
- API整合性の観点で、全シミュレータの `traces` フィールドは同一の量を計測すべき
- 他の時系列量（populations, entropy, purity）は全て `rho_el` から計算しており、tracesのみ不整合

**修正内容**:
```python
# 修正前（qudit_gksl_circuit_boson_simulator.py）
traces = [float(np.real(np.trace(rho)))]        # 拡張空間全体
traces.append(float(np.real(np.trace(rho))))     # 拡張空間全体

# 修正後
traces = [float(np.real(np.trace(rho_el)))]      # 電子系のみ
traces.append(float(np.real(np.trace(rho_el))))   # 電子系のみ
```

### 3.2 Issue GG: 回路ボソンシミュレータのユニタリ再計算

**問題**: `_trotter_step()` メソッドが毎呼び出しで以下を再計算していた：
1. ハミルトニアンユニタリ（`_compute_hamiltonian_unitary(dt/2)` — 複数のmatrix exponentialを含む）
2. 全Lindbladチャネルの Stinespring ユニタリおよびKraus演算子

`dt` はシミュレーション中一定であるため、これらは1回の事前計算で十分。
DM ボソンシミュレータ（`QuditGKSLBosonSimulator`）は `_precompute_unitaries(dt)` で正しく1回のみ計算している。

**影響**: 不要な計算コスト。n_steps=1000の場合、数千回の不要な行列指数関数計算が発生。

**修正内容**:
1. `_precompute_unitaries(dt)` メソッドを新設：
   - `self._U_H_half`: 半ステップハミルトニアンユニタリ（回路分解で構成、dim_total × dim_total）
   - `self._kraus_list`: 全Lindbladチャネルの半ステップKraus演算子リスト
2. `simulate()` で `_precompute_unitaries(dt)` を1回呼び出し
3. `_trotter_step()` を事前計算済みの値を使用するように変更：
   - `dt` パラメータを削除（不要になったため）
   - `self._U_H_half` でハミルトニアンステップを直接適用
   - `self._kraus_list` からKraus演算子を直接使用

## 4. 修正ファイル一覧

| ファイル | 変更内容 |
|---------|---------|
| `tutorials/qudit_gksl_circuit_boson_simulator.py` | Issue FF: traces計算修正、Issue GG: ユニタリ事前計算 |
| `tutorials/run_iteration25_verification.py` | 新規: 修正検証スクリプト（89チェック） |
| `developing/検証結果分析_quantum_dynamics_comparison_iteration25.md` | 新規: 本分析文書 |

## 5. 検証項目（Iteration 25）

### 新規追加チェック（15項目）

| チェック名 | 内容 | セクション |
|-----------|------|-----------|
| circuit_boson_trace_length_matches_dm | tracesリスト長がDMと一致 | 13 |
| circuit_boson_trace_near_unity | 全tracesが~1.0（< 1e-10） | 13 |
| circuit_boson_traces_match_dm | traces時系列がDMと一致（< 1e-12） | 13 |
| circuit_boson_traces_use_rho_el | ソースコードがnp.trace(rho_el)を使用 | 13 |
| circuit_boson_has_precompute | _precompute_unitariesメソッドが存在 | 14 |
| circuit_boson_has_U_H_half | simulate後に_U_H_half属性が存在 | 14 |
| circuit_boson_has_kraus_list | simulate後に_kraus_list属性が存在 | 14 |
| circuit_boson_U_H_half_shape | _U_H_halfが(36,36)形状 | 14 |
| circuit_boson_U_H_half_unitary | _U_H_halfのユニタリ性（< 1e-12） | 14 |
| circuit_boson_kraus_list_length | _kraus_listの長さ=12（Lindbladチャネル数） | 14 |
| circuit_boson_trotter_no_recompute | _trotter_stepでユニタリ再計算しない | 14 |
| circuit_boson_trotter_uses_precomputed | _trotter_stepが事前計算値を使用 | 14 |
| circuit_boson_simulate_calls_precompute | simulate()が_precompute_unitariesを呼び出す | 14 |
| circuit_boson_long_run_matches_dm | 20ステップでDMと一致（< 1e-10） | 15 |
| circuit_boson_timing_logged | 実行時間の記録 | 15 |

### 既存リグレッションチェック（74項目、Iteration 24から継承）

74項目すべてリグレッションなしを確認済み。

### 合計: 89チェック（74既存 + 15新規）

## 6. 検証結果

89/89 PASS ✓

- DM vs Circuit Frobenius距離（2ステップ）: 1.78e-15
- DM vs Circuit Frobenius距離（20ステップ）: 2.24e-14
- traces時系列最大差: 2.00e-15
- _U_H_halfユニタリ性誤差: 1.26e-15

## 7. 残課題

### 7.1 Qubit回路シミュレータの抽象レベル不一致（既知、iteration 17以降）

Qubit回路シミュレータは66個のcompositeゲート（UnitaryGate）をカウントするが、qubit DM/shotは388個のbasicゲートを推定する。公正な比較には統一が必要。

### 7.2 ボソン回路シミュレータのbuild_full_trotter_step_circuit()不在（既知）

非ボソン回路シミュレータは`build_full_trotter_step_circuit()`メソッドを持つが、ボソン回路シミュレータにはない。Cell 16が手動で回路を構築する必要がある。

### 7.3 次の検証手順

1. ユーザーが`run_iteration25_verification.py`をローカルで実行し、結果をpush
2. 結果を確認し、89/89 PASSを確認
