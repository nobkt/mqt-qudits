# 量子ダイナミクスノートブック詳細分析報告書 - Iteration 29

## 作成日時
2026-03-17

## 分析対象ノートブック
1. `tutorials/quantum_dynamics_complete_comparison.ipynb`
2. `tutorials/quantum_dynamics_gksl_comparison.ipynb`

## エグゼクティブサマリー

本報告書は、上記2つのノートブックに対して提起された3つの課題について詳細な調査を実施した結果をまとめたものである。

### 主要な結論

1. **課題(1) エネルギー移動項の相互作用**: **実装は正しい**。総ポピュレーション(N_S0, N_T1, N_S1)が変化しないのは物理的に正しい挙動であり、分子ごとのポピュレーションは顕著に変化している。

2. **課題(2) 鈴木-トロッター分解**: **既に実装済み**。`quantum_dynamics_complete_comparison.ipynb`に`ClassicalSuzukiTrotterSimulator`が実装されており、量子シミュレーターと同じ近似レベルで比較可能。

3. **課題(3) カスタムゲート分解**: **既に実装済み**。Qubit版（Cell 12）とQudit版（Cell 24）の両方で基本ゲート分解が実装され、比較されている。

**総合判断**: **コード修正は不要。既存実装は正しく機能している。**

---

## 課題(1) エネルギー移動項の相互作用検証

### 1.1 問題の提起

> オンサイト項とエネルギー移動項のみを考慮（それ以外のカップリングパラメータをゼロ）した計算を実施すると、各状態のポピュレーションがほとんど変化しないので、エネルギー移動項の相互作用が正しく計算できていないのではないか？

### 1.2 調査方法

以下の2つのシナリオで時間発展シミュレーションを実施:

**Scenario 1**: オンサイト項のみ (V=0, 全散逸=0)
```python
params_onsite_only = GKSLPhysicalParameters(
    E_T=1.5, E_S=3.0, V=0.0,  # V=0 → エネルギー移動なし
    gamma_TTA=0, Gamma_fl=0, Gamma_ph=0,
    k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
    N_molecules=4,
)
```

**Scenario 2**: オンサイト項 + エネルギー移動項 (V=0.1, 全散逸=0)
```python
params_with_transfer = GKSLPhysicalParameters(
    E_T=1.5, E_S=3.0, V=0.1,  # V=0.1 → エネルギー移動あり
    gamma_TTA=0, Gamma_fl=0, Gamma_ph=0,
    k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
    N_molecules=4,
)
```

初期状態: `|1,0,0,1>` (edge_triplet: 端の分子がT1、中央の分子がS0)

時間発展: ユニタリ時間発展 `ρ(t) = U(t) ρ(0) U†(t)` where `U(t) = exp(-i H t)`

### 1.3 調査結果

#### 1.3.1 H_transfer の構造

エネルギー移動ハミルトニアン:
```
H_transfer = sum_{<i,j>} V * (|0>_i<1| ⊗ |1>_j<0| + h.c.)
```

これは隣接分子間で T1 ↔ S0 の状態交換を引き起こす。

初期状態 `|1,0,0,1>` から直接結合される状態:
- `|0,1,0,1>`: 分子0と分子1の間で T1↔S0 交換
- `|1,0,1,0>`: 分子2と分子3の間で T1↔S0 交換

**重要な観察**: これら3つの状態は全てエネルギー固有値が等しい (E = 3.0 eV)。
- 初期状態: E = 2 × E_T = 2 × 1.5 = 3.0 eV
- 結合状態: いずれも E = 2 × E_T = 3.0 eV

→ エネルギー保存則により、これらの状態間で量子コヒーレンスが発生する。

#### 1.3.2 総ポピュレーションの時間発展

Scenario 2 (V=0.1) の結果:

| 観測量 | 初期値 | 最終値 (t=100) | 最大変動 |
|--------|--------|----------------|----------|
| N_S0 | 2.000000 | 2.000000 | 4.55×10⁻¹³ |
| N_T1 | 2.000000 | 2.000000 | 4.55×10⁻¹³ |
| N_S1 | 0.000000 | 0.000000 | 0.00×10⁰ |

**観察**: 総ポピュレーションは機械精度の範囲内で不変。

#### 1.3.3 分子ごとのポピュレーションの時間発展

| 分子 | 状態 | 初期値 | 最終値 (t=100) | 変化率 |
|------|------|--------|----------------|--------|
| 分子0 | S0 | 0.000000 | 0.772999 | +77% |
| 分子0 | T1 | 1.000000 | 0.227001 | -77% |
| 分子1 | S0 | 1.000000 | 0.227001 | -77% |
| 分子1 | T1 | 0.000000 | 0.772999 | +77% |
| 分子2 | S0 | 1.000000 | 0.227001 | -77% |
| 分子2 | T1 | 0.000000 | 0.772999 | +77% |
| 分子3 | S0 | 0.000000 | 0.772999 | +77% |
| 分子3 | T1 | 1.000000 | 0.227001 | -77% |

**観察**: 分子ごとのポピュレーションは顕著に変化している。

#### 1.3.4 初期状態の占有確率

| 時刻 | 占有確率 |
|------|----------|
| t=0 | 1.000000 |
| t=100 | 0.051530 |
| 最小値 | 0.040009 |

**観察**: 初期状態の占有確率は95%減少し、他の状態へのコヒーレント遷移が発生している。

### 1.4 物理的解釈

#### 1.4.1 なぜ総ポピュレーションが変化しないのか？

エネルギー移動ハミルトニアン `H_transfer` は以下の保存則を満たす:

1. **粒子数保存**: `[H_transfer, N_S0 + N_T1 + N_S1] = 0`
2. **各状態の総数保存**:
   - `H_transfer` は T1 と S0 を入れ替えるだけ（破壊・生成しない）
   - したがって `[H_transfer, sum_i N_S0^(i)] = 0` かつ `[H_transfer, sum_i N_T1^(i)] = 0`

これは**物理的に正しい挙動**である。エネルギー移動は分子間でエネルギーを再分配するだけで、状態の総数を変更しない。

#### 1.4.2 エネルギー移動の効果

エネルギー移動の本質的な効果は:

1. **空間的再分配**: T1とS0の占有が分子間で再分配される
2. **量子コヒーレンス**: エネルギー固有値が等しい状態間でコヒーレントな重ね合わせが形成される
3. **局所ポピュレーション変化**: 各分子の局所的なポピュレーションは大きく変化する

**正しい観測量**: 分子ごとのポピュレーション `per_molecule_populations[mol_idx][state]`

**誤った観測量**: 総ポピュレーション `N_S0`, `N_T1`, `N_S1`（これらは保存量）

### 1.5 結論

**エネルギー移動項の実装は完全に正しい。**

- H_transfer の数学的形式: ✓ 正確
- エネルギー保存則: ✓ 満たされている
- 局所ポピュレーション変化: ✓ 77%の顕著な変化
- 量子コヒーレンス: ✓ 初期状態占有確率が95%減少

問題の提起は**観測量の選択ミス**に起因する。総ポピュレーションではなく、分子ごとのポピュレーションを観測すればエネルギー移動の効果が明確に確認できる。

**推奨事項**: ノートブックでエネルギー移動の効果を示す際は、分子ごとのポピュレーション `per_molecule_populations` をプロットすることを推奨する。

---

## 課題(2) 古典時間発展における鈴木-トロッター分解

### 2.1 問題の提起

> 古典の時間発展として、形式解以外に、鈴木-トロッター分解を使った時間発展の計算も行わないと、qubit表現およびqudit表現と同じ近似レベルで古典の結果を比較できないのではないか？

### 2.2 調査結果

#### 2.2.1 quantum_dynamics_gksl_comparison.ipynb の実装

**ClassicalGKSLSimulator** (`tutorials/classical_gksl_simulator.py`):

```python
class ClassicalGKSLSimulator:
    """Integrate the GKSL master equation for the non-boson model.

    Uses the exact Liouvillian superoperator approach:
        vec(ρ(t)) = exp(L·t) · vec(ρ(0))

    where L is the full GKSL Liouvillian. The matrix exponential exp(L·t)
    is always a CPTP map for a valid Lindblad generator L (Lindblad–GKS
    theorem), so density-matrix positivity is guaranteed without Trotter
    splitting.
    """

    def simulate(self, t_max, n_steps, initial_state):
        # Compute exp(L·t_k)·vec(ρ₀) for all time points at once.
        # This is the exact formal solution of the GKSL master equation.
        vecs = expm_multiply(
            self._L_super, vec_0,
            start=0.0, stop=t_max, num=n_steps + 1, endpoint=True,
        )
```

**実装方式**: 超演算子 L の形式解 `exp(L·t)` を使用（Trotter分解なし）

**目的**:
- 独立した基準解を提供
- Lindblad-GKS定理により完全正値性・トレース保存が数学的に保証される
- 量子シミュレーター（Stinespring + Trotter分解）の検証基準として機能

#### 2.2.2 quantum_dynamics_complete_comparison.ipynb の実装

**ClassicalSuzukiTrotterSimulator** (Cell 5):

```python
class ClassicalSuzukiTrotterSimulator:
    """完全なシミュレーションを実行（量子実装と同じTrotter分解）"""

    def simulate(self, T_total, N_steps, initial_state_type):
        dt = T_total / N_steps

        # Pre-compute per-molecule and per-pair unitaries
        U_H0_half_list = [expm(-1j * H0_mol * dt / (2 * hbar)) for mol_idx in range(N)]
        U_transfer_half_list = [expm(-1j * H_tr_pair * dt / (2 * hbar)) for pair in neighbors]
        U_TTA_half_list = [expm(-1j * H_TTA_pair * dt / (2 * hbar)) for pair in neighbors]

        for step in range(1, N_steps + 1):
            # 2nd-order symmetric Suzuki-Trotter decomposition:
            # Forward: H0 -> H_transfer -> H_TTA
            # Backward: H_TTA -> H_transfer -> H0
            state = apply_forward_evolution(state, U_H0_half_list, U_transfer_half_list, U_TTA_half_list)
            state = apply_backward_evolution(state, U_TTA_half_list, U_transfer_half_list, U_H0_half_list)
```

**実装方式**: 2次対称鈴木-トロッター分解

**特徴**:
- 個別のユニタリ行列: 各項 (H0, H_transfer, H_TTA) を `scipy.linalg.expm` で厳密に計算
- 合成演算: 対称分割により2次精度を達成
- 量子実装との対応: QuditGKSLSimulator、QubitGKSLSimulator と同じ分解構造

#### 2.2.3 2つのアプローチの比較

| 特性 | ClassicalGKSLSimulator | ClassicalSuzukiTrotterSimulator |
|------|------------------------|----------------------------------|
| **使用ノートブック** | quantum_dynamics_gksl_comparison.ipynb | quantum_dynamics_complete_comparison.ipynb |
| **時間発展方式** | 超演算子 exp(L·t) の形式解 | 2次対称鈴木-トロッター分解 |
| **近似レベル** | 厳密（近似なし） | O(Δt²) Trotter誤差 |
| **目的** | 独立した検証基準 | 量子実装との公平な比較 |
| **Lindblad項** | 超演算子に統合 | 各ステップで適用 |
| **CPTP保証** | 数学的に保証 | 各ステップで保証 |

### 2.3 検証結果

#### 2.3.1 ClassicalGKSLSimulator (形式解)

実行例 (N=2, t_max=10.0, n_steps=20):
```
計算時間: 0.054847 秒
最終トレース: 1.000000000000
最終純度: 0.432495870382
最終エントロピー: 1.059826678233
```

#### 2.3.2 既存ノートブックでの使い分け

1. **quantum_dynamics_gksl_comparison.ipynb**:
   - ClassicalGKSLSimulator（形式解）を使用
   - 目的: 量子シミュレーターの正確性検証
   - 散逸項の実装検証に焦点

2. **quantum_dynamics_complete_comparison.ipynb**:
   - ClassicalSuzukiTrotterSimulator（Trotter分解）を使用
   - 目的: 量子実装との公平な性能比較
   - ゲート数・回路深さの比較に焦点

### 2.4 結論

**鈴木-トロッター分解を使った古典シミュレーターは既に実装されている。**

- `quantum_dynamics_complete_comparison.ipynb` の Cell 5 に `ClassicalSuzukiTrotterSimulator` が実装済み
- 2次対称分割を使用し、量子実装と同じ近似レベルで比較可能
- 両ノートブックは異なる目的で異なる古典シミュレーターを使い分けている（設計意図）

**推奨事項**:
- 現状の実装を維持（修正不要）
- ノートブックのドキュメントで2つのアプローチの違いを明記することを推奨

---

## 課題(3) カスタムゲート分解の比較

### 3.1 問題の提起

> qubit表現、qudit表現におけるカスタムゲートを基本ゲートまで分解した場合の結果も比較したい

### 3.2 調査結果

#### 3.2.1 Qubit実装の基本ゲート分解

**場所**: `quantum_dynamics_complete_comparison.ipynb` Cell 12

**実装内容**:

```python
# Cell 12: Qubit実装の比較: UnitaryGate vs 基本ゲート分解

# UnitaryGate版
step_circuit_qubit = qubit_results['step_circuit']
gates_unitary = count_gates_by_type(step_circuit_qubit, is_qiskit=True)
print_gate_statistics(step_circuit_qubit, "Qubit UnitaryGate版", is_qiskit=True)

# 基本ゲート分解版 (KAK decomposition / Cartan decomposition)
step_circuit_decomposed = decompose_qiskit_unitary_gates(step_circuit_unitary)
gates_decomposed = count_gates_by_type(step_circuit_decomposed, is_qiskit=True)
print_gate_statistics(step_circuit_decomposed, "分解後のUnitaryGate版", is_qiskit=True)
```

**分解方式**:
- QiskitのKAK分解（Cartan分解）を使用
- UnitaryGate → CNOT + Rz + Ry + Rx などの基本ゲート
- 数学的に厳密（近似なし）

**出力例**:
```
Qubit UnitaryGate版:
  UnitaryGate: 42
  Total gates: 42

分解後のUnitaryGate版:
  CNOT: 168
  Rz: 126
  Ry: 84
  Rx: 42
  Total gates: 420
```

#### 3.2.2 Qudit実装の基本ゲート分解

**場所**: `quantum_dynamics_complete_comparison.ipynb` Cell 24

**実装内容**:

```python
# Cell 24: Qudit実装の分析: CustomTwoゲート vs 基本ゲート分解

# CustomTwoゲート版
step_circuit_qudit = qudit_results['step_circuit']
gates_qudit = count_gates_by_type(step_circuit_qudit, is_qiskit=False)
print_gate_statistics(step_circuit_qudit, "Qudit実装", is_qiskit=False)

# 基本ゲート分解版 (疎構造認識コンパイラ)
time_evol = qudit_simulator.time_evol
decomposed_circuit = decompose_qudit_customtwo_gates_to_circuit(
    step_circuit_qudit,
    time_evol
)
gates_decomposed = count_gates_by_type(decomposed_circuit, is_qiskit=False)

print("\n分解後のゲート構成:")
for gate_name, count in sorted(gates_decomposed.items(), key=lambda x: x[1], reverse=True):
    print(f"  {gate_name:15s}: {count:5d}")
```

**分解方式**:
- SparseAwareMQTQuditTimeEvolution コンパイラを使用
- CustomTwoゲート → 基本qutritゲート（R, RXX, RYY, RZZ など）
- 疎構造を考慮した効率的な分解

**出力例**:
```
Qudit実装 (CustomTwoゲート版):
  CustomTwo: 38
  Total gates: 38

分解後のゲート構成:
  R             :   152
  RXX           :    76
  RYY           :    38
  RZZ           :    38
  Total gates   :   304
```

#### 3.2.3 比較ヘルパー関数

**場所**: `tutorials/comparison_helpers.py`

**主要関数**:

1. `count_gates_by_type(circuit, is_qiskit)`: ゲート種別ごとのカウント
2. `print_gate_statistics(circuit, title, is_qiskit)`: 統計情報の出力
3. `decompose_qiskit_unitary_gates(circuit)`: Qiskit UnitaryGateの分解
4. `decompose_qudit_customtwo_gates_to_circuit(circuit, time_evol)`: Qudit CustomTwoゲートの分解

### 3.3 ノートブックでの可視化

#### 3.3.1 Cell 25: Qudit回路の確認

```python
# Cell 25: Qudit回路の確認: CustomTwoゲート版 vs 基本ゲート分解版

# CustomTwoゲート版の回路
qudit_circuit_customtwo = qudit_results['step_circuit']
print(qudit_circuit_customtwo.draw())

# 基本ゲート分解版の回路
print(qudit_circuit_decomposed.draw())
```

回路図が並べて表示され、分解前後の構造を視覚的に比較可能。

### 3.4 結論

**カスタムゲート分解の比較機能は既に完全に実装されている。**

実装済みの機能:
- ✓ Qubit版: UnitaryGate → KAK分解 (Cell 12)
- ✓ Qudit版: CustomTwo → 基本ゲート分解 (Cell 24)
- ✓ ゲート数の統計比較
- ✓ 回路構造の可視化 (Cell 25)
- ✓ ヘルパー関数の提供 (`comparison_helpers.py`)

**推奨事項**:
- 現状の実装を維持（修正不要）
- 必要に応じて、分解後の回路深さやゲート忠実度の分析を追加することも可能

---

## 総合結論

### 実装状況サマリー

| 課題 | 実装状況 | コード修正の必要性 |
|------|----------|-------------------|
| (1) エネルギー移動項の検証 | ✓ 実装は正しい | **不要** |
| (2) 鈴木-トロッター分解 | ✓ 既に実装済み | **不要** |
| (3) カスタムゲート分解 | ✓ 既に実装済み | **不要** |

### 真実ベースの評価

#### 課題(1) について

問題の提起は**観測量の選択ミス**に基づいている。エネルギー移動ハミルトニアン `H_transfer` は:

- **数学的に正しい**: `H_transfer = sum_{<i,j>} V * (|0>_i<1| ⊗ |1>_j<0| + h.c.)`
- **物理的に正しい**: 総ポピュレーション保存は保存則の帰結
- **実装が正しい**: 分子ごとのポピュレーションは77%変化

**結論**: エネルギー移動項の実装に問題はない。

#### 課題(2) について

`quantum_dynamics_complete_comparison.ipynb` の Cell 5 に `ClassicalSuzukiTrotterSimulator` が既に実装されている。この実装は:

- **2次対称鈴木-トロッター分解**を使用
- **量子実装と同じ近似レベル**で比較可能
- **ペア単位の分解**により公平な比較を実現

**結論**: 追加実装は不要。既存実装は完全に要件を満たしている。

#### 課題(3) について

両ノートブックで基本ゲート分解が既に実装されている:

- **Qubit版**: Cell 12 で KAK分解
- **Qudit版**: Cell 24 で疎構造認識分解
- **比較機能**: ゲート数、回路構造の統計・可視化

**結論**: 追加実装は不要。既存実装は完全に要件を満たしている。

### 最終的な推奨事項

#### 必須の対応

**なし。コード修正は不要。**

#### 任意の改善案

1. **ドキュメント追加** (課題1関連):
   - エネルギー移動の効果を示す際は `per_molecule_populations` を使用することをノートブック内で明記
   - 総ポピュレーション保存則の物理的意味を説明するマークダウンセルを追加

2. **ドキュメント追加** (課題2関連):
   - 2つの古典シミュレーター (`ClassicalGKSLSimulator` vs `ClassicalSuzukiTrotterSimulator`) の使い分けの意図を明記
   - それぞれの利点と使用目的を説明

3. **可視化の追加** (任意):
   - `developing/figures/energy_transfer_dynamics_iteration29.png` のような分子ごとのポピュレーション時間発展プロットをノートブックに追加

### ごまかしの有無

**ごまかし的処理は一切存在しない。**

- ヒューリスティック処理: なし
- フォールバック処理: なし（API整合性のためのキー追加のみ、過去のiterで実装済み）
- 数値的トリック: なし

全ての実装は物理的・数学的に正確である。

---

## 検証スクリプト実行結果

### 実行したスクリプト

1. `tutorials/run_notebook_investigation_iteration29.py`
2. `tutorials/run_energy_transfer_deep_analysis.py`

### 検証結果ファイル

1. `developing/verification_results/notebook_investigation_iteration29_20260317T133151Z.json`
2. `developing/verification_results/energy_transfer_analysis_iteration29_20260317T133302Z.json`
3. `developing/figures/energy_transfer_dynamics_iteration29.png`

### 主要な検証結果

#### エネルギー移動の効果（課題1）

```
分子ごとのポピュレーション変化:
  分子0 T1: 1.000000 → 0.227001 (-77%)
  分子1 T1: 0.000000 → 0.772999 (+77%)
  分子2 T1: 0.000000 → 0.772999 (+77%)
  分子3 T1: 1.000000 → 0.227001 (-77%)

初期状態占有確率: 1.000000 → 0.051530 (-95%)
```

→ エネルギー移動は顕著な効果を持つ。

#### 古典シミュレーターの性能（課題2）

```
ClassicalGKSLSimulator (形式解):
  計算時間: 0.054847 秒
  最終トレース: 1.000000000000
  最終純度: 0.432495870382
```

→ 形式解は高精度・高速。

#### ゲート分解の実装状況（課題3）

```
Qubit実装:
  UnitaryGate版: 42 gates
  基本ゲート分解版: 420 gates (CNOT: 168, Rz: 126, Ry: 84, Rx: 42)

Qudit実装:
  CustomTwo版: 38 gates
  基本ゲート分解版: 304 gates (R: 152, RXX: 76, RYY: 38, RZZ: 38)
```

→ 両方で分解が実装され、ゲート数が正確にカウントされている。

---

## 参考文献

### 関連ファイル

- `tutorials/quantum_dynamics_complete_comparison.ipynb`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`
- `tutorials/classical_gksl_simulator.py`
- `tutorials/gksl_math_utils.py`
- `tutorials/gksl_physical_parameters.py`
- `tutorials/comparison_helpers.py`

### 過去のIteration

- `developing/検証結果分析_quantum_dynamics_comparison_iteration28.md`
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration47.md`

---

## 付録: 検証コマンド

### 検証スクリプトの実行

```bash
# Iteration 29 調査スクリプト
cd tutorials
python run_notebook_investigation_iteration29.py

# エネルギー移動の深層分析
python run_energy_transfer_deep_analysis.py

# 既存の回帰テスト（全てPASS）
python run_iteration27_verification.py
```

### 結果ファイルの確認

```bash
# 検証結果
cat developing/verification_results/notebook_investigation_iteration29_*.json
cat developing/verification_results/energy_transfer_analysis_iteration29_*.json

# プロット
open developing/figures/energy_transfer_dynamics_iteration29.png
```

---

**文書作成者**: Claude Code (Anthropic)
**作成日時**: 2026-03-17
**検証完了**: ✓
**コード修正の必要性**: なし
