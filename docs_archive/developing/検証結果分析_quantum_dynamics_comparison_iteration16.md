# 検証結果分析: Quantum Dynamics Comparison - Iteration 16

## 1. 背景

### 1.1 Iteration 15の結果

Iteration 15ではBug J（qudit_gksl_simulator.pyのestimated_gates_per_step: 59→66）と
Bug K（qubit_gksl_simulator.pyのestimated_gates_per_step: 350→388）が修正され、
15/15チェックが全てPASSした。

### 1.2 Iteration 16で発見された問題

Iteration 15の修正は**密度行列シミュレータ**（`qudit_gksl_simulator.py`と
`qubit_gksl_simulator.py`）のみに適用されていた。同一の回文Trotter構造を持つ
他の4つのシミュレータでは、同じ2×因子の欠落バグが修正されていなかった。

## 2. 発見されたバグの詳細

### 2.1 Bug L: ショットベースquditシミュレータのゲートカウント不整合

**ファイル**: `qudit_gksl_shot_simulator.py`

**問題**: `simulate()`の`estimated_gates_per_step`でハミルトニアン部分が
1回分しかカウントされていなかった。`_trotter_step_trajectory()`メソッドでは
ハミルトニアンが2回適用されている（191行目と200行目）にもかかわらず。

```python
# 修正前（59）
gates_per_step = self.n_system_qudits + len(self.params.neighbors) + self.n_ancilla_qudits * 2

# 修正後（66）
gates_per_step = 2 * (self.n_system_qudits + len(self.params.neighbors)) + self.n_ancilla_qudits * 2
```

### 2.2 Bug M: ショットベースqubitシミュレータのゲートカウント不整合

**ファイル**: `qubit_gksl_shot_simulator.py`

**問題**: 同様に`estimated_gates_per_step`のハミルトニアン基本ゲート推定が
1回分しかカウントされていなかった。

```python
# 修正前（350）
gates_per_step = self.n_sys_qubits + len(self.params.neighbors) * 10 + self.n_ancilla * 6 * 2

# 修正後（388）
gates_per_step = 2 * (self.n_sys_qubits + len(self.params.neighbors) * 10) + self.n_ancilla * 6 * 2
```

### 2.3 Bug N: ボソンquditシミュレータのゲートカウント不整合

**ファイル**: `qudit_gksl_boson_simulator.py`

**問題**: `estimated_gates_per_step`でハミルトニアン部分もStinespring部分も
1回分しかカウントされていなかった。`_trotter_step()`メソッド（87-110行目）では
ハミルトニアンが2回（101行目と109行目）、Stinespringが2回（forward 103-104行目、
reverse 106-107行目）適用されている。

```python
# 修正前
gates_per_step = (
    self.n_system_qudits + self.n_phonon_qudits
    + len(self.params.neighbors) + self.n_ancilla_qudits
)

# 修正後
gates_per_step = (
    2 * (self.n_system_qudits + self.n_phonon_qudits + len(self.params.neighbors))
    + self.n_ancilla_qudits * 2
)
```

### 2.4 Bug O: ボソンqubitシミュレータのゲートカウント不整合

**ファイル**: `qubit_gksl_boson_simulator.py`

**問題**: 同様にハミルトニアン部分もStinespring部分も1回分しかカウントされていなかった。

```python
# 修正前
gates_per_step = (
    self.n_el_qubits + self.n_ph_qubits
    + len(self.params.neighbors) * 10 + self.n_ancilla * 6
)

# 修正後
gates_per_step = (
    2 * (self.n_el_qubits + self.n_ph_qubits + len(self.params.neighbors) * 10)
    + self.n_ancilla * 6 * 2
)
```

## 3. バグの根本原因

### 3.1 コピー&ペーストの不完全な伝播

Iteration 15でBug J/Kが修正された際、修正は密度行列シミュレータ（`qudit_gksl_simulator.py`、
`qubit_gksl_simulator.py`）のみに適用された。同じ回文Trotter構造を持つ
4つの派生シミュレータ（ショットベース2種、ボソン2種）には修正が伝播されなかった。

### 3.2 共通パターン

全6つの非回路シミュレータが同じ回文Trotter構造を使用:
```
H(dt/2) → D₁...Dₙ(dt/2) → Dₙ...D₁(dt/2) → H(dt/2)
```

したがって全てのゲートカウントで:
- **ハミルトニアン**: 2× (2つの半ステップ)
- **Stinespring**: 2× (forward + reverse)

### 3.3 影響の範囲

| シミュレータ | Bug | 影響 |
|-------------|-----|------|
| `qudit_gksl_simulator.py` | Bug J (IT15で修正済) | 報告値のみ |
| `qubit_gksl_simulator.py` | Bug K (IT15で修正済) | 報告値のみ |
| `qudit_gksl_shot_simulator.py` | Bug L (本IT16で修正) | 報告値のみ |
| `qubit_gksl_shot_simulator.py` | Bug M (本IT16で修正) | 報告値のみ |
| `qudit_gksl_boson_simulator.py` | Bug N (本IT16で修正) | 報告値のみ |
| `qubit_gksl_boson_simulator.py` | Bug O (本IT16で修正) | 報告値のみ |
| `qudit_gksl_noisy_simulator.py` | 継承で自動修正済 (IT15) | - |
| `qubit_gksl_noisy_simulator.py` | 継承で自動修正済 (IT15) | - |
| `qudit_gksl_noisy_shot_simulator.py` | 継承で自動修正済 (本IT16) | - |
| `qubit_gksl_noisy_shot_simulator.py` | 継承で自動修正済 (本IT16) | - |

全てのバグは報告値（`estimated_gates_per_step`）のみの問題であり、
物理的なシミュレーション結果（人口動態、トレース保存等）には影響しない。

## 4. 修正内容

### 4.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `qudit_gksl_shot_simulator.py` | `estimated_gates_per_step`: 59 → 66; コメント更新 |
| `qubit_gksl_shot_simulator.py` | `estimated_gates_per_step`: 350 → 388; コメント更新 |
| `qudit_gksl_boson_simulator.py` | `estimated_gates_per_step`: 2×因子追加; コメント更新 |
| `qubit_gksl_boson_simulator.py` | `estimated_gates_per_step`: 2×因子追加; コメント更新 |

### 4.2 変更不要

- 回路シミュレータ（`*_circuit_simulator.py`）: 正しい（Iteration 14で検証済み）
- Noisy密度行列シミュレータ: 親クラスから正しい値を継承
- Noisyショットシミュレータ: 親クラスから正しい値を継承（本修正後）
- `_trotter_step()` / `_trotter_step_trajectory()`: 全シミュレータで正しく実装済み

## 5. 検証結果

`tutorials/run_iteration16_verification.py` で18項目を検証:

| # | チェック項目 | 結果 |
|---|-------------|------|
| 1 | quditショット GPS = 66 | PASS |
| 2 | quditショット GPS ≠ 59（旧値） | PASS |
| 3 | qubitショット GPS = 388 | PASS |
| 4 | qubitショット GPS ≠ 350（旧値） | PASS |
| 5 | quditショット = qudit密度行列（66） | PASS |
| 6 | qubitショット = qubit密度行列（388） | PASS |
| 7 | quditノイジーショット GPS = 66 | PASS |
| 8 | qubitノイジーショット GPS = 388 | PASS |
| 9 | quditボソン ソースに2×因子あり | PASS |
| 10 | qubitボソン ソースに2×因子あり | PASS |
| 11 | quditショット total = 132 | PASS |
| 12 | qubitショット total = 776 | PASS |
| 13 | 回帰: qudit密度行列 GPS = 66 | PASS |
| 14 | 回帰: qubit密度行列 GPS = 388 | PASS |
| 15 | 回帰: qudit回路 GPS = 66 | PASS |
| 16 | 回帰: qubit回路 GPS = 66 | PASS |
| 17 | 回帰: quditショット トレース保存 | PASS |
| 18 | 回帰: qubitショット トレース保存 | PASS |

## 6. 次のステップ

1. ユーザーが検証スクリプトをローカルで実行して結果を確認
2. ノートブック（`quantum_dynamics_gksl_comparison.ipynb`）でショットベース
   シミュレーション結果のゲートカウント表示を確認
3. ボソンシミュレータの`estimated_gates_per_step`に電子-フォノン結合ゲートが
   含まれていない件の検討（回路ボソンシミュレータでは`n_eph_gates`が含まれる）
