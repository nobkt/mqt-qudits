# 検証結果分析: Quantum Dynamics Comparison - Iteration 14

## 1. 背景

### 1.1 Iteration 13の結果

Iteration 13では6つのバグ（A〜F）が修正され、33/33チェックが全てPASSした。
しかし、回文Trotter修正に伴うゲートカウント報告の整合性が一部更新されていなかった。

### 1.2 Iteration 14で発見された問題

Iteration 13の検証結果（`iteration13_palindromic_trotter_20260311T112739Z.json`）と
ソースコードを詳細に分析した結果、以下の3つの不整合が発見された。

## 2. 発見されたバグの詳細

### 2.1 Bug G: クラスdocstringのゲートカウント不整合

**問題**: 3つの回路シミュレータのクラスdocstringが、回文Trotter修正前の
ゲートカウント（40ゲート/ステップ）を記載していた。

**影響ファイル**:
- `qudit_gksl_circuit_simulator.py`: "Total: 40 gates per Trotter step"
- `qubit_gksl_circuit_simulator.py`: "Total: 40 gates per Trotter step"
- `qudit_gksl_circuit_boson_simulator.py`: Stinespring 20/6（1回分）

**修正後**: 回文構造を反映し、Stinespringゲートを2倍（forward+reverse）、
合計66ゲート/ステップに更新。

### 2.2 Bug H: `compile_to_native_gates()` per_step_summary の不整合

**ファイル**: `qudit_gksl_circuit_simulator.py`

**問題**: `per_step_summary`のネイティブゲート数計算で、Stinespringチャネルの
ネイティブゲート数が回文の2倍を反映していなかった。

```python
# 修正前（非回文カウント）
total_native_per_step = 2 * ham_native_half + single_native_total
total_uncompiled_per_step = pair_uncompiled_total

# 修正後（回文カウント）
total_native_per_step = 2 * ham_native_half + 2 * single_native_total
total_uncompiled_per_step = 2 * pair_uncompiled_total
```

### 2.3 Bug I: `transpile_to_basic_gates()` per_step_summary の不整合

**ファイル**: `qubit_gksl_circuit_simulator.py`

**問題**: `per_step_summary`の基本ゲート数計算で、Stinespringチャネルの
基本ゲート数が回文の2倍を反映していなかった。

```python
# 修正前（非回文カウント）
total_single = sum(r["basic_gates"] for r in single_results)
total_pair = sum(r["basic_gates"] for r in pair_results)

# 修正後（回文カウント）
total_single = sum(r["basic_gates"] for r in single_results) * 2  # fwd+rev
total_pair = sum(r["basic_gates"] for r in pair_results) * 2  # fwd+rev
```

## 3. バグ発見の経緯

### 3.1 なぜIteration 13で見落とされたか

1. Iteration 13の検証スクリプトは`simulate()`の`gates_per_step`と`gate_breakdown`を
   検証したが、`compile_to_native_gates()`と`transpile_to_basic_gates()`の
   `per_step_summary`は検証対象外だった
2. クラスdocstringの整合性は自動検証されていなかった
3. 既存テスト（`test_compile_to_native_gates_o0`）は`per_step_summary`の値を
   `> 0`でしかチェックしておらず、正確な値を検証していなかった

## 4. 修正内容

### 4.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `qudit_gksl_circuit_simulator.py` | docstring: 66ゲート/ステップ; `compile_to_native_gates`: 2× |
| `qubit_gksl_circuit_simulator.py` | docstring: 66ゲート/ステップ; `transpile_to_basic_gates`: 2× |
| `qudit_gksl_circuit_boson_simulator.py` | docstring: Stinespring 40/12 |

### 4.2 変更なし

- `_trotter_step_circuit`: Iteration 13で正しく修正済み
- `simulate()`: `gates_per_step`と`gate_breakdown`は正しい（2×済み）
- `build_full_trotter_step_circuit`: Iteration 13で正しく修正済み

## 5. 次のステップ

1. ユーザーが検証スクリプトをローカルで実行して結果を確認
2. `compile_to_native_gates()`と`transpile_to_basic_gates()`の`per_step_summary`値を
   検証するテストを追加（Iteration 15で検討）
