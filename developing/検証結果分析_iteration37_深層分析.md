# Iteration 37: 深層分析結果

- 分析日時: 2026-03-03
- 分析対象: Iteration 36検証結果、ノートブック全セル、コード全体（全シミュレータ）

## 1. 分析スコープ

以下を網羅的に分析した：

1. Iteration 36 の全14チェック検証結果（全PASS）
2. `quantum_dynamics_gksl_comparison.ipynb` の全40セルの実行結果
3. コアモジュール: `stinespring_utils.py`, `gksl_math_utils.py`, `classical_gksl_simulator.py`
4. DM版シミュレータ: `qudit_gksl_simulator.py`, `qubit_gksl_simulator.py`
5. ショットベースシミュレータ: `qudit_gksl_shot_simulator.py`, `qubit_gksl_shot_simulator.py`
6. ノイズDMシミュレータ: `qudit_gksl_noisy_simulator.py`, `qubit_gksl_noisy_simulator.py`
7. 数学的理論: Stinespring dilation, GKSL superoperator, Trotter分割, 回文順序積

## 2. Iteration 36 結果の確認

全14チェックがPASSしており、コア実装（ClassicalGKSLSimulator, QuditGKSLSimulator）の数学的正確性は確認済み。

## 3. 発見された問題点

### 問題1: ショットベースとDM版のTrotterステップ構造の不整合

**重大度: 中**

#### 概要

理想DM版シミュレータ（`QuditGKSLSimulator`, `QubitGKSLSimulator`）と理想ショットベースシミュレータ（`QuditGKSLShotSimulator`, `QubitGKSLShotSimulator`）が**異なるTrotterステップ構造**を使用している。

#### DM版（理想、ノイズなし）のTrotterステップ

```
exp(L_H dt/2) · [∏_{α=1..n} E_α(dt/2)] · [∏_{α=n..1} E_α(dt/2)] · exp(L_H dt/2)
```

- `_U_stines_half`（半ステップ dt/2）を使用
- **回文順序**（forward + reverse）で適用
- 各チャネルを2回（各 dt/2）適用 → 実効的に dt の散逸

#### ショットベース（理想、ノイズなし）のTrotterステップ

```
exp(L_H dt/2) · [∏_{α=1..n} E_α(dt)] · exp(L_H dt/2)
```

- `_U_stines`（全ステップ dt）を使用
- **前方のみ**（回文順序なし）で適用
- 各チャネルを1回（全 dt）適用

#### 理論的影響

1. **Stinespring近似誤差**:
   - DM版: 各チャネルの誤差 ≈ C₂·(dt/2)² × 2回 = C₂·dt²/2 per step
   - ショット版: 各チャネルの誤差 ≈ C₂·dt² per step
   - → ショット版はDM版の約**2倍**のStinespring近似誤差

2. **Lie-Trotter積誤差**:
   - DM版: 回文順序により O(dt³)/step（2次消去）
   - ショット版: 前方のみで O(dt²)/step（1次）

3. **全体収束次数**: 両者とも O(dt) だが、ショット版は誤差定数が大きい

#### コードの証拠

**DM版** (`qudit_gksl_simulator.py` L100-124):
```python
def _trotter_step(self, rho):
    rho = self._U_H_half @ rho @ self._U_H_half.conj().T
    for U_stine_half in self._U_stines_half:        # dt/2
        rho = apply_stinespring_to_density_matrix(rho, U_stine_half)
    for U_stine_half in reversed(self._U_stines_half):  # 回文
        rho = apply_stinespring_to_density_matrix(rho, U_stine_half)
    rho = self._U_H_half @ rho @ self._U_H_half.conj().T
    return rho
```

**ショットベース** (`qudit_gksl_shot_simulator.py` L162-176):
```python
def _trotter_step_trajectory(self, psi, rng):
    psi = self._U_H_half @ psi
    for U_stine in self._U_stines:              # dt（全ステップ）
        psi = self._apply_stinespring_with_measurement(psi, U_stine, rng)
    psi = self._U_H_half @ psi                  # 回文なし
    return psi
```

#### 影響

ショットベースシミュレータのdocstringは「The quantum trajectory method is mathematically equivalent to the density matrix approach in the limit of infinite shots」と記載しているが、これは**不正確**。無限ショット極限でも異なるTrotterステップ構造に起因する系統的差異が残る。

### 問題2: ノイズDMシミュレータもDM版と異なるTrotterステップ構造

**重大度: 中**

ノイズDM版シミュレータ（`QuditGKSLNoisySimulator`, `QubitGKSLNoisySimulator`）も `_U_stines`（全ステップ、前方のみ）を使用しており、理想DM版と異なるTrotterステップ構造を持つ。

**コードの証拠** (`qudit_gksl_noisy_simulator.py` L246):
```python
for k, U_stine in enumerate(self._U_stines):  # 全ステップ、前方のみ
    rho = apply_stinespring_to_density_matrix(rho, U_stine)
```

これは「ノイズなし極限（p_depol=0, p_dephasing=0）でもDM版と一致しない」ことを意味する。

### 問題3: ゲート数見積もりの不正確性

**重大度: 低**

`QuditGKSLSimulator` のゲート数見積もり:
```python
gates_per_step = 4 + 3 + 26  # 33 gates
```

しかし実際のTrotterステップでは26チャネルを2回（前方＋逆方向）適用しているため、正確なゲート数は:
```
4 (VirtRz) + 3 (CustomTwo) + 26×2 (Stinespring) = 59 gates
```

### 問題4: docstringの不正確な記述

**重大度: 低**

1. `qudit_gksl_shot_simulator.py` L15: 「mathematically equivalent to the density matrix approach in the limit of infinite shots」→ 現状のTrotterステップ構造では不正確
2. `qudit_gksl_shot_simulator.py` L167: 「2nd-order Trotter step」→ H-D分割は2次だが、チャネル積は1次（前方のみ）

## 4. 修正方針

### 方針

すべてのシミュレータ（理想、ショットベース、ノイズDM、ノイズショットベース）で**同一のTrotterステップ構造**（半ステップ回文順序）を使用する。

#### 理由

1. **数学的一貫性**: 理想DM版が最も正確な近似であり、他のシミュレータはそれを基盤とすべき
2. **等価性の保証**: ショットベース → DM版（無限ショット極限）、ノイズDM版 → DM版（ノイズ=0極限）
3. **物理的正確性**: 量子回路が回文順序を使用する場合、ノイズも各ゲートに対して適用されるべき

#### 変更点

1. ショットベースシミュレータ: `_precompute_unitaries` に `_U_stines_half` を追加、`_trotter_step_trajectory` を回文順序に変更
2. ノイズDMシミュレータ: `_trotter_step` を `_U_stines_half` + 回文順序に変更
3. ノイズショットベースシミュレータ: `_trotter_step_trajectory` を回文順序に変更
4. ゲート数見積もりの修正
5. docstringの修正

#### ノイズシミュレータへの影響

回文順序化により各Lindbladチャネルゲートが2回/ステップ適用されるため、ゲートノイズも2倍/ステップになる。これは物理的に正しい（より多くのゲート → より多くのノイズ）が、既存のノイズシミュレーション結果が変化する。

## 5. 結論

コア実装（ClassicalGKSLSimulator, QuditGKSLSimulator）は数学的に正確であり、Iteration 36の全検証チェックはPASS。しかし、ショットベースシミュレータおよびノイズシミュレータには**Trotterステップ構造の不整合**があり、docstringの記述と矛盾する。修正はヒューリスティックやfallbackではなく、数学的一貫性を回復するための構造的修正である。
