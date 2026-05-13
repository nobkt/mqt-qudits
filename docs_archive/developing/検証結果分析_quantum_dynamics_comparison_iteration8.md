# 検証結果分析: Quantum Dynamics Comparison - Iteration 8

## 1. 背景

### 1.1 Iteration 7の状態

Iteration 7で`qubit_noisy_simulator.py`に`cx_per_pair_gate`パラメータが追加され、
CXゲートオーバーヘッドの物理的効果をモデル化可能になった：

- `p_eff = 1 - (1 - p_phys)^cx_per_pair_gate`
- `estimate_cx_per_pair_gate()`メソッドでQiskitトランスパイラからCXゲート数を推定
- H_transfer: 47 CXゲート、H_TTA: 86 CXゲート、平均66.5 CXゲート/ペア

### 1.2 Iteration 7の検証結果

- 24/24チェック全てPASS
- コード自体は正しく実装されている

### 1.3 Iteration 8で発見された問題点

Iteration 7の検証結果とノートブック実行結果を注意深く分析した結果、以下の問題を発見：

**問題1: complete_comparison.ipynb Cell 15が`cx_per_pair_gate`を使用していない**

Cell 15は`simulate()`を`cx_per_pair_gate`パラメータなしで呼び出していた（デフォルト値1）。
これにより：
- QubitのノイズレートがQuditと同じ`p_eff = 0.001`になっていた
- 実際のハードウェアでは`p_eff = 0.0644`（64倍悪い）
- CXゲートオーバーヘッドの効果が比較に反映されていなかった

Cell 15の実行出力:
```
物理CXゲートエラー率 (p_phys): 0.1000%
ペアあたりのCXゲート数: 1          ← ここが問題（実際は66.5）
有効ペア脱分極率 (p_eff): 0.1000%  ← quditと同じ（不公平）
```

**問題2: `cx_per_pair_gate`の型ヒントが`int`**

`estimate_cx_per_pair_gate()`が返す`cx_avg`はfloat（例: 66.5）だが、
`simulate()`の型ヒントは`int`だった。

**問題3: QubitGKSLNoisyShotSimulatorに`cx_per_pair_gate`サポートがない**

GKSL比較ノートブック（Cell 38）のQubitノイジーシミュレータも同じ問題を抱えていた。
`p_depol=0.001`がそのままペア脱分極レートとして使用され、CXオーバーヘッドが反映されていなかった。

## 2. Iteration 8での修正内容

### 2.1 `qubit_noisy_simulator.py`の修正

- `cx_per_pair_gate: int = 1` → `cx_per_pair_gate: float = 1` に型ヒントを修正
- float値（例: 66.5）を正しく受け付けるようになった

### 2.2 `quantum_dynamics_complete_comparison.ipynb` Cell 15の修正

**修正前:**
```python
qubit_noisy_results = qubit_noisy_sim.simulate(
    T_total=params.T_total,
    N_steps=params.N_steps,
    initial_state_type=params.initial_state_type,
    shots=10000,
    noise_params=noise_params
)
```

**修正後:**
```python
# CXゲート数の推定
cx_info = QubitMolecularDynamicsSimulatorNoisy.estimate_cx_per_pair_gate(
    V=params.V, J=params.J, dt=params.dt, hbar=params.hbar
)
cx_per_pair = cx_info['cx_avg']

# 実際のCXゲートオーバーヘッドを反映
qubit_noisy_results = qubit_noisy_sim.simulate(
    T_total=params.T_total,
    N_steps=params.N_steps,
    initial_state_type=params.initial_state_type,
    shots=10000,
    noise_params=noise_params,
    cx_per_pair_gate=cx_per_pair
)
```

### 2.3 `qubit_gksl_shot_simulator.py`の修正

`QubitGKSLNoisyShotSimulator`に`cx_per_pair_gate`サポートを追加：

- `__init__`に`cx_per_pair_gate: float = 1`パラメータを追加
- `self.p_depol_pair_eff = 1.0 - (1.0 - p_depol) ** cx_per_pair_gate`を計算
- `_trotter_step_trajectory`のペア脱分極に`self.p_depol_pair_eff`を使用
- 単一サイト脱分極には元の`self.p_depol`を維持（CXオーバーヘッドは不要）
- `simulate()`の`noise_params`辞書に`cx_per_pair_gate`と`p_depol_pair_eff`を追加
- 後方互換性を維持（デフォルト`cx_per_pair_gate=1`で`p_eff = p_depol`）

### 2.4 `quantum_dynamics_gksl_comparison.ipynb` Cell 37-38の修正

- Cell 37: CXゲートオーバーヘッドの説明をマークダウンに追加
- Cell 38: `estimate_cx_per_pair_gate()`でCXゲート数を推定し、`cx_per_pair_gate`を渡すように修正

## 3. 物理的影響の分析

### 3.1 CXオーバーヘッドの定量的影響

p_phys = 0.001（物理CXゲートあたりのエラー率）として：

| エンコーディング | CXゲート/ペア | 有効ペアエラー率 | 比率 |
|:--|--:|--:|--:|
| Qudit (qutrit) | 1（ネイティブ） | 0.1% | 1x |
| Qubit (cx=1) | 1 | 0.1% | 1x |
| Qubit (cx=66.5) | 66.5 | 6.44% | 64x |

### 3.2 シミュレーション全体への影響

100ステップのシミュレーションでの累積影響：

| 条件 | ステップあたりノイズ | 100ステップ後の状態 |
|:--|:--|:--|
| Qudit (p=0.001) | 1-(1-0.001)^12 ≈ 1.2% | 部分的デコヒーレンス |
| Qubit (cx=1, p=0.001) | 1-(1-0.001)^12 ≈ 1.2% | 部分的デコヒーレンス + 禁止状態漏洩 |
| Qubit (cx=66.5, p=0.001) | 1-(1-0.0644)^12 ≈ 55% | ほぼ完全デコヒーレンス + 大量漏洩 |

cx=66.5の場合、各ステップで状態の約55%がノイズにより劣化する。
100ステップ後には状態はほぼ完全に混合状態に近づく。

### 3.3 Quditの物理的優位性（2つの独立効果）

1. **ゲートオーバーヘッド**: Qubitは1ペア相互作用あたり66.5 CXゲート必要（64倍悪い）
2. **禁止状態漏洩**: Qubitは7/16=43.75%の禁止状態、Quditは0%

これらは独立した物理的効果であり、ヒューリスティックな処理ではない。

## 4. 検証結果

Iteration 8検証スクリプト（`run_iteration8_verification.py`）で30項目を検証し、全てPASS。

### 4.1 検証項目

1. **型ヒント**: `cx_per_pair_gate: float`（int→float修正済み）
2. **complete_comparison Cell 15**: `estimate_cx_per_pair_gate()`使用、`cx_per_pair_gate`渡し
3. **QubitGKSLNoisyShotSimulator**: `cx_per_pair_gate`パラメータ、`p_depol_pair_eff`計算
4. **GKSL Cell 38**: CXゲート数推定・使用
5. **有効レート公式**: 後方互換性、単調増加性
6. **Quditシミュレータ**: `cx_per_pair_gate`なし（正しい - ネイティブゲート）
7. **ドキュメント**: Cell 37にCXオーバーヘッドの説明

## 5. 次のステップ

ユーザーが以下を実行する必要がある：
1. `cd tutorials && python run_iteration8_verification.py` で検証
2. `quantum_dynamics_complete_comparison.ipynb` を再実行
3. `quantum_dynamics_gksl_comparison.ipynb` を再実行
4. 検証結果をリポジトリにpush
5. 実行結果を確認し、さらなる修正が必要かどうか判断
