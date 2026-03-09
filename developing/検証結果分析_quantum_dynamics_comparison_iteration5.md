# 量子ダイナミクス比較ノートブック 検証結果分析レポート（Iteration 5）

## 作成日: 2026-03-09
## 対象ファイル:
- `tutorials/quantum_dynamics_complete_comparison.ipynb`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`
## 前回参照: `developing/検証結果分析_quantum_dynamics_comparison_iteration4.md`

---

## 1. セル実行結果の詳細分析

### 1.1 `quantum_dynamics_complete_comparison.ipynb`

#### Cell 15: Qubitノイズモデル付きシミュレーション（実行済み出力）

```
1トロッターステップあたりのゲート数: 4330
回路深さ: 2374

ステップ 10/100: t = 10.00 fs, N_T1 = 0.4190, N_S1 = 0.4313, 非物理的: 0.6779
ステップ 20/100: t = 20.00 fs, N_T1 = 0.4272, N_S1 = 0.4176, 非物理的: 0.6793
...
ステップ 100/100: t = 100.00 fs, N_T1 = 0.4255, N_S1 = 0.4268, 非物理的: 0.6778

最終個体数:
  N_S0 = 0.4365  N_T1 = 0.4255  N_S1 = 0.4268
  非物理的状態: 0.677800
```

**問題**: 非物理的状態が67.8%! 物理的状態の合計は0.4365+0.4255+0.4268 = 1.29 (4分子で正規化すると1.29/4 = 32.2%のみが物理的)

#### Cell 23: Quditノイズモデル付きシミュレーション（実行済み出力）

```
2-qudit gates per Trotter step: 12
Effective per-step noise: 1-(1-0.01)^12 = 0.1136

最終個体数:
  N_S0 = 1.3404  N_T1 = 1.3309  N_S1 = 1.3287
  最終純度: 0.012346
```

**状態**: N_S0+N_T1+N_S1 = 4.0 → トレース保存 ✓。純度=0.012は強いデコヒーレンスだが物理的に妥当。

### 1.2 `quantum_dynamics_gksl_comparison.ipynb`

#### Cell 34: Quditショットベース（ノイズ有り, 実行済み出力）

```
Elapsed time: 105.97s
Noise: p_depol=0.01, p_dephasing=0.005
Final populations: N_S0=1.5107, N_T1=1.2143, N_S1=1.2751
```

合計: 1.5107+1.2143+1.2751 = 4.0001 ≈ 4.0 ✓（トレース保存）

**問題**: `p_dephasing=0.005`（位相緩和ノイズ）が適用されている。要件は脱分極ノイズのみ。  
**問題**: `depol_pair_only=False`（デフォルト）により、単一サイトのLindblad チャネルにもノイズが適用されている。

#### Cell 38: Qubitショットベース（ノイズ有り, 実行済み出力）

```
Elapsed time: 262.94s
Noise: p_depol=0.0100, p_dephasing=0.005
Final populations: N_S0=0.4417, N_T1=0.3834, N_S1=0.3773
Tr(rho_final): 0.300600 (1.0 = no leakage)
Forbidden state count (final measurement): 702/1000
Max trace deficit (leakage): 0.7217 (72.2%)
```

**重大問題**: トレースが0.3006 → 70%の確率が禁止状態にリークしている！
禁止状態への72.2%のリークは完全に非物理的な結果をもたらしている。

---

## 2. 根本原因分析

### 問題1: `qubit_noisy_simulator.py` のQiskit Aerノイズモデルが過剰（最重要・重大度：致命的）

**原因**: `QubitMolecularDynamicsSimulatorNoisy.simulate()`メソッドはQiskit Aerのノイズモデルを使用し、`'cx'`ゲートに脱分極ノイズを適用する。

#### ノイズ蓄積の問題:

| 項目 | 値 |
|------|-----|
| 1トロッターステップあたりの基本ゲート数（変換後） | 4,330 |
| `cx`ゲートに適用する脱分極確率 | 1% |
| 有効1ステップノイズ | 1-(0.99)^4330 ≈ 100% |
| 100ステップ後の非物理的状態 | 67.8% |

`UnitaryGate`（H_transferおよびH_TTA用の4-qubit行列）がQiskit AerによってCX+RZゲートに変換される際、1ステップあたり約4330個のcxゲートが生成される。各cxゲートが1%の脱分極ノイズを受けると、1ステップ後にすでに100%近いノイズが蓄積し、状態が完全にランダム化される。

#### quditノイズモデルとの比較:

| 項目 | Qudit（正常） | Qubit（問題あり） |
|------|------------|------------|
| 2-分子ゲート数/ステップ | 12 | 4,330 (cx変換後) |
| 有効ステップノイズ | 11.4% | ≈100% |
| 最終純度 | 0.012 | N/A（非物理的） |

**正しいアプローチ**: Quditシミュレータと同様に、密度行列形式で2-分子相互作用（H_transfer, H_TTA）のみに脱分極ノイズを適用する。

#### Qubitエンコーディングにおける2-分子ゲートの定義:
- 1分子 = 2物理qubit（|S0⟩=|00⟩, |T1⟩=|01⟩, |S1⟩=|10⟩, 禁止=|11⟩）
- 2-分子インタラクション = 4物理qubitゲート（16×16ユニタリ）
- H_transfer, H_TTA は各分子ペアに1回の脱分極ノイズ（quditの場合と同様）

### 問題2: gksl_comparison.ipynb の位相緩和ノイズと単一サイトノイズ（重大度：高）

**原因**: `QuditGKSLNoisyShotSimulator`と`QubitGKSLNoisyShotSimulator`が:
1. `p_dephasing=0.005`（位相緩和ノイズ）を適用 → 要件は「脱分極のみ」に反する
2. `depol_pair_only=False`（デフォルト）→ 単一サイトLindblad チャネルにもノイズが適用される

**Qubit GKSLでの影響**:
- 単一サイトLindblad チャネルにも`_apply_stochastic_qubit_depolarization_single`が適用
- Weyl-Heisenberg演算子が2-qubitエンコーディングの禁止状態|11⟩を生成
- 結果: 72.2%の禁止状態リーク → 完全に非物理的な結果

---

## 3. 修正内容

### 修正1: `tutorials/qubit_noisy_simulator.py` の全面的な書き直し

**現在の実装**:
- Qiskit Aerのショットシミュレータ（回路レベルノイズモデル）
- `'cx'`, `'cz'`ゲートに脱分極ノイズ → 変換後4,330 cxゲート/ステップ

**新実装**:
- Quditノイズシミュレータと同様の密度行列シミュレーション
- 各分子を「d=4次元量子ビット」として扱う（|S0⟩=0, |T1⟩=1, |S1⟩=2, 禁止=3）
- 全体次元: 4^4 = 256（= 2^8 ✓）
- H0: 単一分子（4×4対角行列）→ ノイズなし
- H_transfer: 2-分子ペア（16×16ユニタリ）→ 脱分極ノイズ
- H_TTA: 2-分子ペア（16×16ユニタリ）→ 脱分極ノイズ
- 脱分極チャネル: ε(ρ) = (1-p)ρ + p·Tr_{mol_i,mol_j}(ρ) ⊗ I_{16}/16

#### 期待される改善:

| 項目 | 修正前 | 修正後 |
|------|--------|--------|
| 2-分子ゲート数/ステップ | 4,330 (cx変換後) | 12 |
| 有効ステップノイズ | ≈100% | 11.4% |
| 非物理的状態 | 67.8% | <1%（禁止状態への微小リーク） |

### 修正2: `quantum_dynamics_gksl_comparison.ipynb` Cell 34

```python
# 修正前
sim5c = QuditGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.005)

# 修正後
sim5c = QuditGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.0, depol_pair_only=True)
```

**理由**: 「2-quditゲートのみ脱分極ノイズ」の要件に対応。
- `p_dephasing=0.0`: 位相緩和ノイズを無効化（脱分極のみ）
- `depol_pair_only=True`: 単一サイトLindblad チャネルへのノイズを無効化

### 修正3: `quantum_dynamics_gksl_comparison.ipynb` Cell 38

```python
# 修正前
sim3c = QubitGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.005)

# 修正後
sim3c = QubitGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.0, depol_pair_only=True)
```

**理由**: 同上。`depol_pair_only=True`により禁止状態リークが大幅に減少する見込み。

### 修正4: マークダウンセルの更新

- `complete_comparison.ipynb` Cell 14: 密度行列ベースの新ノイズモデルを説明するよう更新
- `gksl_comparison.ipynb` Cells 33, 37: 「脱分極のみ、ペアのみ」に更新

---

## 4. 技術的詳細

### 4.1 新しいQubitノイズシミュレータの密度行列実装

#### 状態空間の再解釈:
```
1分子 = d=4次元量子ビット (|S0⟩=0, |T1⟩=1, |S1⟩=2, 禁止=3)
4分子 = 4^4 = 256次元
密度行列 ρ: 256×256
```

#### H0ユニタリ（分子i, 半ステップ dt/2）:
```python
U_H0_mol = diag([1, exp(-iE_T dt/(2ℏ)), exp(-iE_S dt/(2ℏ)), 1])
```
（禁止状態はH0=0として扱う）

#### H_transfer/H_TTA ユニタリ（分子ペア(i,j), 半ステップ）:
- `build_H_transfer_qubit_unitary(V, dt/2, hbar)` → 16×16 U_transfer
- `build_H_TTA_qubit_unitary(J, dt/2, hbar)` → 16×16 U_TTA

#### 2-分子脱分極チャネル（分子ペア(i,j), d=4）:
```
ε(ρ) = (1-p)ρ + p·Tr_{mol_i,mol_j}(ρ) ⊗ I_{4×4}/4 ⊗ I_{4×4}/4
```
quditノイズシミュレータの`_apply_2qudit_depolarizing_dm`と同じ実装（d=4で使用）

#### 1トロッターステップあたりのノイズ適用:
- H_transfer × 3ペア × 2（前進+後退）= 6回
- H_TTA × 3ペア × 2 = 6回
- 計12回（quditと同じ）

### 4.2 GKSLノイズシミュレータの修正

#### `depol_pair_only=True`の効果（qudit版）:
- `_trotter_step_trajectory`内:
  - Hamiltonianペアノイズ（transfer gates）: 適用 ✓
  - ペアLindbladチャネル（TTA）: 適用 ✓
  - 単一サイトLindbladチャネル（緩和）: スキップ ✓

#### `p_dephasing=0.0`の効果:
- すべての位相緩和ノイズを無効化
- 脱分極ノイズのみが残る ✓

---

## 5. 期待される結果

### `quantum_dynamics_complete_comparison.ipynb`:

| シミュレーション | 修正前 | 修正後（期待） |
|----------------|--------|----------------|
| Qubit noisy 非物理状態 | 67.8% | <1% |
| Qubit noisy N_T1 最終値 | ~0.43 (非物理) | 物理的な値 |
| Qudit noisy N_T1 最終値 | ~1.33 (変化なし) | 同じ（既に正しい） |

### `quantum_dynamics_gksl_comparison.ipynb`:

| シミュレーション | 修正前 | 修正後（期待） |
|----------------|--------|----------------|
| Qubit noisy 禁止状態リーク | 72.2% | <5% |
| Qudit noisy N_S0 最終値 | 1.5107 (dephasing含む) | より小さい変化 |

---

## 6. 作業ログ

| ステップ | 内容 |
|--------|------|
| ① 検証スクリプト作成（Iteration 1-3） | ノートブック実行・検証 |
| ② ユーザーによるローカル実行・push | Iteration 4の結果が確認 |
| ③ PR作成 (PR#218) | Iteration 4の問題を確認 |
| ④ Iteration 5修正（本ドキュメント） | ノイズモデルの根本的な修正 |

### Iteration 5での修正方針:
- **ごまかしゼロ**: ヒューリスティックなfallbackや近似は使用しない
- **物理的正確性**: qubitとquditで対称なノイズモデル（2-分子ゲートのみ）
- **密度行列**: Qiskit Aer回路レベルノイズから密度行列方式に変更（quditと同様）
- **真実ベース**: 現状の67-72%非物理状態を正直に報告し、根本から修正

---

## 7. 次回の検証で確認すべき事項

1. `complete_comparison.ipynb` Cell 15: 非物理的状態が<5%に改善されること
2. `complete_comparison.ipynb` Cell 29: QubitノイジーとQuditノイジーの比較プロットが正常であること
3. `gksl_comparison.ipynb` Cell 38: 禁止状態リークが<10%に改善されること
4. `gksl_comparison.ipynb` Cell 42 (9c): ショットベース比較プロットが更新されること
5. 修正後のノイズレベルが物理的に妥当であること（Qubit Noisy ≒ Qudit Noisy in noise strength）
