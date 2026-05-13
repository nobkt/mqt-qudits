# 検証結果分析: quantum_dynamics_comparison iteration 23

## 1. 概要

Iteration 22の検証結果（49/49 PASS）およびノートブック実行結果を真実ベースで詳細に分析した。
シミュレータの物理的計算結果は正確であるが、**回路可視化コードに2件のバグ**と**レジスタ計算式に1件の一般性欠如**を発見・修正した。

### 1.1 分析対象

1. `iteration22_boson_params_20260313T124031Z.json`: 49/49 PASS
2. `quantum_dynamics_gksl_comparison.ipynb` の全セル実行結果（特にCell 16）
3. `quantum_dynamics_complete_comparison.ipynb` の全セル実行結果
4. 全6シミュレータおよび回路可視化コードのソースコード詳細レビュー

### 1.2 分析結論

- **シミュレータコードの物理計算は正確**: 全シミュレータの密度行列演算・個体群計算・トレース保存は正しい
- **Issue BB発見**: ノートブックCell 16のボソン回路可視化がパリンドローム構造を欠如
- **Issue CC発見**: qudit bosonのn_phonon_qudits計算式がn_max > d-1で不正確

## 2. 発見されたバグ

### 2.1 Issue BB: Cell 16ボソン回路可視化のパリンドローム構造欠如（重大）

**問題**: `quantum_dynamics_gksl_comparison.ipynb` のCell 16は、ボソンシナリオ（Scenario 4, 6）の回路を手動構築しているが、以下の2つの問題がある：

1. **パリンドローム逆順が欠如**: Forward Lindblad channelsのみ表示し、Reverse（palindromic）パスが無い
2. **Stinespring dtが不正**: 全dt（`dt_quantum`）を使用しているが、正しくは半dt（`dt_quantum/2`）

**影響範囲**:

| シナリオ | 修正前subcircuits | 正しいsubcircuits | dt使用 |
|----------|------------------|------------------|--------|
| Scenario 6 (Qudit Boson) | 14 (= 2 + 12) | 26 (= 2 + 12 + 12) | dt → dt/2 |
| Scenario 4 (Qubit Boson) | 14 (= 2 + 12) | 26 (= 2 + 12 + 12) | dt → dt/2 |
| Scenario 3 (Qubit 非ボソン) | 54 ✓ | 54 ✓ | dt/2 ✓ |
| Scenario 5 (Qudit 非ボソン) | 54 ✓ | 54 ✓ | dt/2 ✓ |

**根本原因**: 非ボソンのScenario 3/5は`build_full_trotter_step_circuit()`メソッドを使用し、正しいパリンドローム構造を生成する。しかしボソン回路シミュレータには`build_full_trotter_step_circuit()`メソッドが存在せず、Cell 16のコードが回路を手動構築した際にパリンドローム構造を省略し、dtも誤った値を使用していた。

**重要**: これは**回路可視化のみのバグ**であり、実際のシミュレーション結果（Cell 12, 14の密度行列計算）には影響しない。実際のシミュレータ（`QuditGKSLBosonSimulator._trotter_step()`、`QubitGKSLBosonSimulator._trotter_step()`）はパリンドローム構造を正しく実装している。

**修正内容**:

Scenario 6（Qudit Boson）:
```python
# 修正前: 単一パス、全dt
for idx, ... in enumerate(lindblad_local_info):
    circ = build_stinespring_circuit_single(L_local, dt_quantum, ...)  # 不正
    qudit_boson_step_circuits.append((f"stinespring_single_{idx}", circ))

# 修正後: Forward + Reverse、半dt
# Forward half-step
for idx, ... in enumerate(lindblad_local_info):
    circ = build_stinespring_circuit_single(L_local, dt_quantum / 2, ...)  # 修正
    qudit_boson_step_circuits.append((f"stinespring_fwd_single_{idx}", circ))
# Reverse half-step (palindromic)
for idx, ... in reversed(list(enumerate(lindblad_local_info))):
    circ = build_stinespring_circuit_single(L_local, dt_quantum / 2, ...)  # 修正
    qudit_boson_step_circuits.append((f"stinespring_rev_single_{idx}", circ))
```

Scenario 4（Qubit Boson）:
```python
# 修正前: 単一パス、全dt
U_stinespring = stinespring_unitary_from_lindblad(L_op, dt_quantum)  # 不正

# 修正後: Forward + Reverse、半dt
U_stinespring = stinespring_unitary_from_lindblad(L_op, dt_quantum / 2)  # 修正
```

### 2.2 Issue CC: n_phonon_qudits計算式の一般性欠如

**問題**: `qudit_gksl_boson_simulator.py` の `n_phonon_qudits = params.N_molecules` はn_max ≤ d-1（d=3でn_max ≤ 2）の場合のみ正確。

**修正前**:
```python
self.n_phonon_qudits = params.N_molecules   # phonon qutrits (for n_max=2)
```

**修正後**:
```python
# Number of d-level qudits needed per phonon mode: ceil(log_d(n_max+1))
n_ph_per_mol = int(np.ceil(np.log(params.n_max + 1) / np.log(params.d)))
self.n_phonon_qudits = n_ph_per_mol * params.N_molecules
```

**影響**: 現在のテストパラメータ（n_max=1,2, d=3）では結果が変わらない（ceil(log_3(2))=1, ceil(log_3(3))=1）。n_max=3以上で正しい値を返すようになる（ceil(log_3(4))=2）。

## 3. 修正ファイル一覧

| ファイル | 変更内容 |
|---------|---------|
| `tutorials/quantum_dynamics_gksl_comparison.ipynb` Cell 16 | Issue BB: ボソン回路のパリンドローム修正 |
| `tutorials/qudit_gksl_boson_simulator.py` | Issue CC: n_phonon_qudits一般化 |
| `tutorials/run_iteration23_verification.py` | 新規: 修正検証スクリプト |
| `developing/検証結果分析_quantum_dynamics_comparison_iteration23.md` | 新規: 本分析文書 |

## 4. 検証項目（Iteration 23）

### 新規追加チェック

| チェック名 | 内容 |
|-----------|------|
| cell16_scenario6_subcircuits_26 | Scenario 6が26 subcircuitsを含むか |
| cell16_scenario4_subcircuits_26 | Scenario 4が26 subcircuitsを含むか |
| cell16_scenario6_has_fwd_rev | Scenario 6がfwd/revラベルを含むか |
| cell16_scenario4_has_fwd_rev | Scenario 4がfwd/revラベルを含むか |
| cell16_scenario6_uses_half_dt | Scenario 6がdt/2を使用しているか |
| cell16_scenario4_uses_half_dt | Scenario 4がdt/2を使用しているか |
| n_phonon_qudits_nmax1 | n_max=1でn_phonon_qudits = N ✓ |
| n_phonon_qudits_nmax2 | n_max=2でn_phonon_qudits = N ✓ |
| n_phonon_qudits_nmax3 | n_max=3でn_phonon_qudits = 2*N（一般化修正） |

### 既存リグレッションチェック（Iteration 22から継承）

49項目すべてリグレッションなしを確認する。

## 5. 残課題

### 5.1 Qubit回路シミュレータの抽象レベル不一致

既知問題（iteration 17以降）: qubit回路シミュレータは66個のcompositeゲート（UnitaryGate）をカウントするが、qubit DM/shotは388個のbasicゲートを推定する。これは公正な比較のための統一が依然必要。

### 5.2 ボソン回路シミュレータの`build_full_trotter_step_circuit()`不在

非ボソン回路シミュレータ（`QubitGKSLCircuitSimulator`、`QuditGKSLCircuitSimulator`）は`build_full_trotter_step_circuit()`メソッドを持ち、パリンドローム構造を正しく生成する。しかしボソン回路シミュレータ（`QuditGKSLCircuitBosonSimulator`）にはこのメソッドが存在せず、ノートブックCell 16が回路を手動構築する必要がある。将来的にはこのメソッドを追加することが望ましい。

### 5.3 次の検証手順

1. ユーザーが`run_iteration23_verification.py`をローカルで実行し、結果をpush
2. 結果を確認し、ノートブックを再実行してCell 16の出力がsubcircuits=26になることを確認
