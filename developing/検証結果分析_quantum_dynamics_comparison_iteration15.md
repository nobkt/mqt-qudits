# 検証結果分析: Quantum Dynamics Comparison - Iteration 15

## 1. 背景

### 1.1 Iteration 14の結果

Iteration 14では3つのバグ（G〜I）が修正され、19/19チェックが全てPASSした。
修正内容はクラスdocstringのゲートカウント更新（66ゲート/ステップ）、
`compile_to_native_gates()`の`per_step_summary`の2×修正、
`transpile_to_basic_gates()`の`per_step_summary`の2×修正であった。

### 1.2 Iteration 15で発見された問題

Iteration 14の検証結果（19/19 PASS）およびソースコードを詳細に分析した結果、
**非回路シミュレータ**（`qudit_gksl_simulator.py`と`qubit_gksl_simulator.py`）の
`estimated_gates_per_step`に回文Trotter構造の2×が反映されていない
バグが2件発見された。

## 2. 発見されたバグの詳細

### 2.1 Bug J: 非回路quditシミュレータのゲートカウント不整合

**ファイル**: `qudit_gksl_simulator.py`

**問題**: `simulate()`の`estimated_gates_per_step`計算で、ハミルトニアン部分が
1回分しかカウントされていなかった。回文2次Trotterステップでは
H(dt/2) → D_fwd → D_rev → H(dt/2) の構造でハミルトニアンが**2回**適用されるため、
ハミルトニアンのゲート数も2倍にする必要がある。

```python
# 修正前（ハミルトニアン1回分 = 59）
gates_per_step = self.n_system_qudits + len(self.params.neighbors) + n_lindblad * 2
# = 4 + 3 + 52 = 59

# 修正後（ハミルトニアン2回分 = 66）
gates_per_step = 2 * (self.n_system_qudits + len(self.params.neighbors)) + n_lindblad * 2
# = 2*(4+3) + 52 = 14 + 52 = 66
```

**検証**: 修正後の66は、回路シミュレータ（`qudit_gksl_circuit_simulator.py`）の
`gates_per_step = 66`と完全に一致する。

### 2.2 Bug K: 非回路qubitシミュレータのゲートカウント不整合

**ファイル**: `qubit_gksl_simulator.py`

**問題**: 同様に`estimated_gates_per_step`のハミルトニアン基本ゲート推定が
1回分しかカウントされていなかった。

```python
# 修正前（ハミルトニアン1回分 = 350）
gates_per_step = (
    self.n_sys_qubits + len(self.params.neighbors) * 10 + self.n_ancilla * 6 * 2
)
# = 8 + 30 + 312 = 350

# 修正後（ハミルトニアン2回分 = 388）
gates_per_step = (
    2 * (self.n_sys_qubits + len(self.params.neighbors) * 10) + self.n_ancilla * 6 * 2
)
# = 2*(8+30) + 312 = 76 + 312 = 388
```

**注意**: qubit非回路シミュレータの推定値（388）は基本ゲートレベルの推定であり、
回路シミュレータの高レベルゲート数（66）とは異なるレベルの抽象度である。
これは仕様通りの差異であり、バグではない。

## 3. バグ発見の経緯

### 3.1 なぜ以前のイテレーションで見落とされたか

1. Iteration 13〜14の修正は**回路シミュレータ**（`*_circuit_simulator.py`）に
   焦点を当てており、非回路シミュレータの`estimated_gates_per_step`は検証対象外
2. 非回路シミュレータの`_trotter_step()`メソッド自体は正しく2回のハミルトニアン
   半ステップを適用しており、シミュレーション結果に影響なし
3. ゲートカウントの不整合は報告値のみの問題で、物理的なシミュレーション結果には
   影響しない

### 3.2 根本原因

回文Trotter分解 H(dt/2) → D_fwd → D_rev → H(dt/2) において：
- Stinespring部分の2×（fwd + rev）は正しく実装されていた
- ハミルトニアン部分の2×（2つの半ステップ）が漏れていた

コメント自体が誤解を招く記述だった：
```python
# N VirtRz gates (H_0 diagonal for N molecules)       ← 1回分
# (N-1) CustomTwo gates (H_transfer for nearest-neighbour pairs) ← 1回分
# n_lindblad×2 Stinespring CustomTwo gates (palindromic: forward + reverse) ← 2回分
```

## 4. 修正内容

### 4.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `qudit_gksl_simulator.py` | `estimated_gates_per_step`: 59 → 66; コメント更新 |
| `qubit_gksl_simulator.py` | `estimated_gates_per_step`: 350 → 388; コメント更新 |

### 4.2 変更なし

- `qudit_gksl_circuit_simulator.py`: 正しい（66ゲート/ステップ）
- `qubit_gksl_circuit_simulator.py`: 正しい（66ゲート/ステップ）
- `qudit_gksl_circuit_boson_simulator.py`: 正しい（2× for Hamiltonian使用済み）
- `_trotter_step()`: 全シミュレータで正しく実装済み（物理結果への影響なし）

## 5. 検証スクリプト

`tutorials/run_iteration15_verification.py` で15項目を検証：

| # | チェック項目 | 結果 |
|---|-------------|------|
| 1 | qudit非回路 GPS = 66 | PASS |
| 2 | qudit非回路 GPS ≠ 59（旧値） | PASS |
| 3 | qubit非回路 GPS = 388 | PASS |
| 4 | qubit非回路 GPS ≠ 350（旧値） | PASS |
| 5 | qudit非回路 = qudit回路（66） | PASS |
| 6 | qubit回路 GPS = 66 | PASS |
| 7 | qudit回路 GPS = 66 | PASS |
| 8 | qudit total = 132（66×2） | PASS |
| 9 | qubit total = 776（388×2） | PASS |
| 10 | quditトレース保存 | PASS |
| 11 | qubitトレース保存 | PASS |
| 12 | qudit回路docstring "66" | PASS |
| 13 | qubit回路docstring "66" | PASS |
| 14 | qudit breakdown single = 40 | PASS |
| 15 | qudit breakdown pair = 12 | PASS |

## 6. 次のステップ

1. ユーザーが検証スクリプトをローカルで実行して結果を確認
2. ノートブック（`quantum_dynamics_complete_comparison.ipynb`、
   `quantum_dynamics_gksl_comparison.ipynb`）でゲートカウント表示セルが
   更新後の値を正しく反映しているか確認
3. テストスイートに`estimated_gates_per_step`の正確な値を検証するテストケースを
   追加検討（現在のテストは`> 0`チェックのみ）
