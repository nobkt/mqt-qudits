# TTA-UC現象 GKSL-Lindblad量子ダイナミクス 詳細設計書

## 文書情報

**作成日**: 2026年2月12日  
**バージョン**: 1.0.0  
**上位文書**:  
- `TTA-UC現象のGKSL-Lindblad量子ダイナミクス完全理論書.md`（理論的基礎）  
- `TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md`（実装仕様）  

**目的**: 上記2文書の要件を完全に満たし、誰でも仕様書通りにTTA-UC現象のGKSL-Lindblad量子ダイナミクスを実装できる解像度の詳細設計を提供する  
**対象ノートブック**: `tutorials/quantum_dynamics_complete_comparison.ipynb`  
**対象フレームワーク**: MQT-Qudits / Qiskit  

---

## 厳密性の原則

本設計書は以下の原則に厳格に従う：

**✅ 許容される手法**：
- GKSL定理に基づく数学的に厳密な導出
- Stinespring dilationによる非ユニタリ演算の量子回路表現
- 実験的に測定可能なパラメータのみの使用
- 制御可能な近似誤差（Trotter分解）の明示的評価

**❌ 禁止される手法**：
- ヒューリスティックな近似や経験的フィッティング
- Fallback処理（計算失敗時の「適当な値」への置き換え）
- 物理的根拠のない簡略化
- ごまかしや真実を隠蔽する記述

---

## 目次

1. [設計概要](#1-設計概要)
2. [データ構造設計](#2-データ構造設計)
3. [数学的基盤モジュールの設計](#3-数学的基盤モジュールの設計)
4. [シナリオ1: 古典GKSL・ボソン無し 詳細設計](#4-シナリオ1-古典gkslボソン無し-詳細設計)
5. [シナリオ2: 古典GKSL・ボソン有り 詳細設計](#5-シナリオ2-古典gkslボソン有り-詳細設計)
6. [シナリオ3: Qubit GKSL・ボソン無し 詳細設計](#6-シナリオ3-qubit-gkslボソン無し-詳細設計)
7. [シナリオ4: Qubit GKSL・ボソン有り 詳細設計](#7-シナリオ4-qubit-gkslボソン有り-詳細設計)
8. [シナリオ5: Qudit GKSL・ボソン無し 詳細設計](#8-シナリオ5-qudit-gkslボソン無し-詳細設計)
9. [シナリオ6: Qudit GKSL・ボソン有り 詳細設計](#9-シナリオ6-qudit-gkslボソン有り-詳細設計)
10. [検証設計](#10-検証設計)
11. [可視化・比較フレームワーク設計](#11-可視化比較フレームワーク設計)
12. [全体統合フローチャート](#12-全体統合フローチャート)
13. [付録A: 記号一覧](#付録a-記号一覧)
14. [付録B: 単位系変換](#付録b-単位系変換)
15. [付録C: エラーハンドリング設計](#付録c-エラーハンドリング設計)
16. [付録D: Lindblad超演算子の要素展開](#付録d-lindblad超演算子の要素展開)

---

## 1. 設計概要

### 1.1 本設計書のスコープ

本設計書は、4分子直線配置TTA-UC系のGKSL-Lindblad量子ダイナミクスを以下の6シナリオで実装するための詳細設計を定義する：

| # | 計算方式 | ボソン相互作用 | 略称 | 実装優先度 |
|---|---------|--------------|------|-----------|
| 1 | 古典計算 | 無し | Classical-NB | 1（最優先） |
| 2 | 古典計算 | 有り | Classical-B | 4 |
| 3 | Qubit量子計算 | 無し | Qubit-NB | 3 |
| 4 | Qubit量子計算 | 有り | Qubit-B | 6 |
| 5 | Qudit量子計算 | 無し | Qudit-NB | 2 |
| 6 | Qudit量子計算 | 有り | Qudit-B | 5 |

### 1.2 前提条件

#### 1.2.1 物理系の定義

- **分子数**: $N = 4$（直線配置）
- **各分子の電子状態数**: $d = 3$（$|S_0\rangle$, $|T_1\rangle$, $|S_1\rangle$）
- **状態空間次元**: $d^N = 3^4 = 81$
- **密度行列サイズ**: $81 \times 81 = 6561$ 要素
- **超演算子サイズ**: $6561 \times 6561 = 43046721$ 要素
- **隣接ペア**: $\{(0,1), (1,2), (2,3)\}$（3ペア）

#### 1.2.2 Qutrit基底の対応

$$
|S_0\rangle_i \longleftrightarrow |0\rangle_i, \quad |T_1\rangle_i \longleftrightarrow |1\rangle_i, \quad |S_1\rangle_i \longleftrightarrow |2\rangle_i
$$

#### 1.2.3 GKSL方程式（本設計の中核）

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}_{\text{sys}}, \hat{\rho}] + \sum_{\alpha=1}^{26} \gamma_\alpha \mathcal{D}[\hat{L}_\alpha][\hat{\rho}]
$$

ここで Lindblad超演算子は：

$$
\mathcal{D}[\hat{L}][\hat{\rho}] = \hat{L}\hat{\rho}\hat{L}^\dagger - \frac{1}{2}\hat{L}^\dagger\hat{L}\hat{\rho} - \frac{1}{2}\hat{\rho}\hat{L}^\dagger\hat{L}
$$

ハミルトニアン（ユニタリ部分）：

$$
\hat{H}_{\text{sys}} = \hat{H}_0 + \hat{H}_{\text{transfer}}
$$

**重要**: TTA過程はハミルトニアンから**除外**し、Lindblad散逸項として実装する。

#### 1.2.4 ソフトウェア前提

| ライブラリ | 用途 |
|-----------|------|
| NumPy | 行列演算 |
| SciPy | ODE積分、行列指数関数、疎行列 |
| Qiskit | Qubit量子回路 |
| MQT-Qudits | Qudit量子回路 |
| Matplotlib | 可視化 |

### 1.3 設計方針

1. **現行ノートブックとの互換性**: 出力データ構造を上位互換とし、可視化関数を再利用可能にする
2. **段階的実装**: 優先度順にPhase 1〜5で実装（§1.1参照）
3. **モジュール構成**: 各シナリオを独立したクラスとして実装し、共通基盤を共有する

### 1.4 ファイル構成設計

```
tutorials/
├── gksl_physical_parameters.py           # パラメータクラス（全シナリオ共通）
├── stinespring_utils.py                  # Stinespring共通ユーティリティ
├── classical_gksl_simulator.py           # シナリオ1
├── classical_gksl_boson_simulator.py     # シナリオ2
├── qubit_gksl_simulator.py              # シナリオ3
├── qubit_gksl_boson_simulator.py        # シナリオ4
├── qudit_gksl_simulator.py              # シナリオ5
├── qudit_gksl_boson_simulator.py        # シナリオ6
├── test_gksl_simulators.py              # テスト
├── quantum_dynamics_gksl_comparison.ipynb # 統合ノートブック
└── doc/GKSL/
    ├── TTA-UC現象のGKSL-Lindblad量子ダイナミクス完全理論書.md
    ├── TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md
    └── TTA-UC現象GKSL-Lindblad量子ダイナミクス詳細設計書.md  ← 本文書
```

---

## 2. データ構造設計

### 2.1 GKSLPhysicalParametersクラス

#### 2.1.1 クラス図

```
┌──────────────────────────────────────────┐
│         GKSLPhysicalParameters           │
├──────────────────────────────────────────┤
│ + N_molecules: int = 4                   │
│ + E_T: float = 1.5           # eV       │
│ + E_S: float = 3.0           # eV       │
│ + V: float = 0.1             # eV       │
│ + gamma_TTA: float = 0.05    # eV/ℏ     │
│ + Gamma_fl: float = 0.01     # eV/ℏ     │
│ + Gamma_ph: float = 1e-6     # eV/ℏ     │
│ + k_IC: float = 0.005        # eV/ℏ     │
│ + k_ISC_ST: float = 0.003    # eV/ℏ     │
│ + k_ISC_TS: float = 1e-5     # eV/ℏ     │
│ + hbar: float = 0.6582119569 # eV·fs    │
│ + neighbors: List[Tuple]                 │
│ + T_total: float = 100.0     # fs       │
│ + N_steps: int = 100                     │
│ + dt: float                              │
│ + initial_state_type: str                │
├──────────────────────────────────────────┤
│ + validate() → bool                      │
│ + get_dissipation_rates() → Dict         │
│ + check_time_scale_hierarchy() → bool    │
└──────────────────────────────────────────┘
         ▲
         │ 継承
┌──────────────────────────────────────────┐
│    GKSLPhysicalParametersWithBoson       │
├──────────────────────────────────────────┤
│ + n_max: int = 2                         │
│ + omega_ph: float = 0.15     # eV/ℏ     │
│ + g_eph: float = 0.02        # eV       │
│ + huang_rhys: float (computed)           │
├──────────────────────────────────────────┤
│ + get_phonon_dim() → int                 │
│ + get_total_dim() → int                  │
└──────────────────────────────────────────┘
```

#### 2.1.2 validate()メソッドの設計

以下のチェックを実施し、全て通過した場合のみTrueを返す：

```
validate():
  1. 2 * E_T >= E_S であることを確認（TTA条件）
  2. 全速度定数 > 0 であることを確認
  3. N_molecules >= 2 であることを確認
  4. N_steps >= 1 であることを確認
  5. T_total > 0 であることを確認
  6. hbar > 0 であることを確認
  7. dt = T_total / N_steps が正であることを確認
  違反時: ValueError を raise（Fallbackは禁止）
```

#### 2.1.3 パラメータの物理的根拠と階層構造

各パラメータ値の設計根拠を明示する：

| パラメータ | 値 | 設計根拠 | 実験的範囲 |
|-----------|-----|---------|-----------|
| $\gamma_{\text{TTA}} = 0.05$ | eV/ℏ | 現行 $J = 0.05$ eV からの対応。特性時間 $\tau_{\text{TTA}} \approx 13.2$ fs | $10^{-3}$–$10^{-1}$ |
| $\Gamma_{\text{fl}} = 0.01$ | eV/ℏ | 現行ノートブックと同値。100fs窓で効果を観測するための教育的選択 | $\sim 10^{-7}$ |
| $\Gamma_{\text{ph}} = 10^{-6}$ | eV/ℏ | 蛍光より4–6桁遅い（スピン禁制） | $\sim 10^{-10}$ |
| $k_{\text{IC}} = 0.005$ | eV/ℏ | 蛍光と同程度のオーダー | $10^{-8}$–$10^{-4}$ |
| $k_{\text{ISC}}^{S\to T} = 0.003$ | eV/ℏ | スピン-軌道結合による | $10^{-7}$–$10^{-4}$ |
| $k_{\text{ISC}}^{T\to S} = 10^{-5}$ | eV/ℏ | 大エネルギーギャップのため遅い | $\sim 10^{-10}$ |

時間スケール階層構造（保存すべき不等式）：

$$
\gamma_{\text{TTA}} > \Gamma_{\text{fl}} > k_{\text{IC}} > k_{\text{ISC}}^{S \to T} \gg \Gamma_{\text{ph}},\ k_{\text{ISC}}^{T \to S}
$$

$$
0.05 > 0.01 > 0.005 > 0.003 \gg 10^{-6},\ 10^{-5}
$$

### 2.2 初期状態の設計

#### 2.2.1 初期密度行列の構築手順

`edge_triplet` 初期状態の構築フロー：

```
1. 初期状態ベクトルの定義:
   |ψ(0)⟩ = |T₁⟩₀ ⊗ |S₀⟩₁ ⊗ |S₀⟩₂ ⊗ |T₁⟩₃ = |1,0,0,1⟩

2. Qutrit基底でのインデックス計算:
   idx = 1×3³ + 0×3² + 0×3¹ + 1×3⁰ = 27 + 0 + 0 + 1 = 28

3. 状態ベクトルの構築:
   psi = np.zeros(81, dtype=complex)
   psi[28] = 1.0

4. 密度行列の構築:
   rho_0 = np.outer(psi, psi.conj())
   # 結果: 81×81行列で (28,28)要素のみが1、他は全て0
```

#### 2.2.2 密度行列の検証

構築直後に以下を検証する：

$$
\text{Tr}[\hat{\rho}_0] = 1, \quad \hat{\rho}_0 = \hat{\rho}_0^\dagger, \quad \hat{\rho}_0^2 = \hat{\rho}_0 \quad (\text{純粋状態})
$$

### 2.3 出力データ構造設計

#### 2.3.1 統一出力辞書

全シナリオで返す辞書の構造（現行ノートブックと上位互換）：

```python
GKSLResult = {
    # ── 現行ノートブック互換フィールド ──
    'times': List[float],                      # [0, dt, 2*dt, ..., T_total]
    'populations': List[Dict[str, float]],     # [{'N_S0':f, 'N_T1':f, 'N_S1':f}, ...]
    'per_molecule_populations': List[Dict[str, np.ndarray]],
    # [{'S0_per_mol':ndarray(4,), 'T1_per_mol':ndarray(4,), 'S1_per_mol':ndarray(4,)}, ...]
    'elapsed_time': float,                     # 実行時間（秒）
    'method': str,                             # 手法名

    # ── GKSL版追加フィールド ──
    'entropy': List[float],                    # von Neumannエントロピー S(t)
    'purity': List[float],                     # 純度 Tr[ρ²]
    'trace': List[float],                      # トレース Tr[ρ]（検証用、常に≈1）
    'rho_final': np.ndarray,                   # 最終密度行列 (dim × dim)
    'coherences': List[Dict[str, float]],      # 主要コヒーレンスの時間発展

    # ── シナリオ固有フィールド（存在する場合のみ） ──
    'ode_solver': Optional[str],               # 古典: 'RK45', 'BDF' 等
    'n_function_evals': Optional[int],         # 古典: 関数評価回数
    'step_circuit': Optional[object],          # 量子: 1トロッターステップの回路
    'total_gates': Optional[int],              # 量子: 総ゲート数
    'total_depth': Optional[int],              # 量子: 総回路深さ
    'gates_per_step': Optional[int],           # 量子: ステップ当たりゲート数
    'depth_per_step': Optional[int],           # 量子: ステップ当たり深さ
    'n_ancilla': Optional[int],                # 量子: ancilla数
    'shots': Optional[int],                    # 量子: ショット数
}
```

#### 2.3.2 populationsの計算設計

密度行列 $\hat{\rho}(t)$ から個体数を計算する手順：

$$
N_{X}(t) = \sum_{i=0}^{3} \text{Tr}\left[\hat{P}_{X}^{(i)} \hat{\rho}(t)\right]
$$

射影演算子 $\hat{P}_{X}^{(i)}$ の構築：

$$
\hat{P}_{S_0}^{(i)} = \hat{I}_3^{\otimes i} \otimes |0\rangle\langle 0| \otimes \hat{I}_3^{\otimes (3-i)}
$$

$$
\hat{P}_{T_1}^{(i)} = \hat{I}_3^{\otimes i} \otimes |1\rangle\langle 1| \otimes \hat{I}_3^{\otimes (3-i)}
$$

$$
\hat{P}_{S_1}^{(i)} = \hat{I}_3^{\otimes i} \otimes |2\rangle\langle 2| \otimes \hat{I}_3^{\otimes (3-i)}
$$

**効率的な計算**: $\text{Tr}[\hat{P}_{X}^{(i)} \hat{\rho}]$ は密度行列の対角要素のみから計算可能。

```
calculate_populations(rho, N=4, d=3):
  dim = 81
  N_S0 = N_T1 = N_S1 = 0.0
  S0_per_mol = zeros(4)
  T1_per_mol = zeros(4)
  S1_per_mol = zeros(4)

  FOR idx IN 0..80:
    prob = real(rho[idx, idx])  # 対角要素 = 確率
    config = index_to_config(idx, N=4, d=3)
    FOR mol_idx, level IN enumerate(config):
      IF level == 0: S0_per_mol[mol_idx] += prob; N_S0 += prob
      IF level == 1: T1_per_mol[mol_idx] += prob; N_T1 += prob
      IF level == 2: S1_per_mol[mol_idx] += prob; N_S1 += prob

  RETURN {N_S0, N_T1, N_S1, S0_per_mol, T1_per_mol, S1_per_mol}
```

`index_to_config` の設計：

```
index_to_config(idx, N=4, d=3) → List[int]:
  config = []
  FOR i IN 0..N-1:
    config.append(idx % d)
    idx = idx // d
  RETURN reversed(config)
  # 例: idx=28 → [1, 0, 0, 1]（MSB first）

  注: 実際にはMSB-first表記に合わせるため:
  config = []
  remaining = idx
  FOR i IN range(N):
    config.insert(0, remaining % d)
    remaining = remaining // d
  RETURN config
```

**重要**: `index_to_config` の基底順序は現行ノートブックと一致させる。すなわち、インデックス $\text{idx}$ は以下で計算される：

$$
\text{idx} = \sum_{i=0}^{N-1} c_i \times d^{N-1-i}
$$

ここで $c_i$ は分子 $i$ の状態（0, 1, 2）。

#### 2.3.3 追加観測量の計算設計

**von Neumannエントロピー**:

$$
S(t) = -\text{Tr}[\hat{\rho}(t) \ln \hat{\rho}(t)] = -\sum_k \lambda_k \ln \lambda_k
$$

```
compute_entropy(rho):
  eigenvalues = eigvalsh(rho)  # エルミート行列の実固有値
  eigenvalues = eigenvalues[eigenvalues > 1e-15]  # 数値的ゼロを除外
  RETURN -sum(eigenvalues * log(eigenvalues))
```

**純度**:

$$
P(t) = \text{Tr}[\hat{\rho}(t)^2]
$$

```
compute_purity(rho):
  RETURN real(trace(rho @ rho))
```

**コヒーレンス**（分子 $i$ の $T_1$-$S_1$ コヒーレンス）:

$$
C_{\text{T1-S1}}^{(i)}(t) = |\langle 1|_i \hat{\rho}_{\text{reduced}}^{(i)}(t) |2\rangle_i|
$$

```
compute_coherences(rho, N=4, d=3):
  coherences = {}
  FOR mol_idx IN 0..3:
    rho_reduced = partial_trace_except(rho, mol_idx, N, d)  # 3×3 縮約密度行列
    coherences[f'C_T1S1_mol{mol_idx}'] = abs(rho_reduced[1, 2])
  RETURN coherences
```

---

## 3. 数学的基盤モジュールの設計

### 3.1 ハミルトニアン構築

#### 3.1.1 オンサイトエネルギー $\hat{H}_0$ の構築手順

単一分子のオンサイトハミルトニアン（$3 \times 3$行列）：

$$
\hat{H}_0^{(i)} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1.5 & 0 \\ 0 & 0 & 3.0 \end{pmatrix} \text{ eV}
$$

全空間（$81 \times 81$）への拡張：

$$
\hat{H}_0 = \sum_{i=0}^{3} \hat{I}_3^{\otimes i} \otimes \hat{H}_0^{(i)} \otimes \hat{I}_3^{\otimes (3-i)}
$$

```
build_H0(params):
  N = params.N_molecules  # 4
  d = 3
  dim = d^N  # 81
  H0 = zeros(dim, dim, dtype=complex)

  h_single = diag([0, params.E_T, params.E_S])  # 3×3

  FOR mol_idx IN 0..N-1:
    # テンソル積の構築
    H_full = identity(1)
    FOR i IN 0..N-1:
      IF i == mol_idx:
        H_full = kron(H_full, h_single)
      ELSE:
        H_full = kron(H_full, eye(d))
    H0 += H_full

  RETURN H0  # 81×81 エルミート対角行列
```

**構築後の検証**: $\hat{H}_0 = \hat{H}_0^\dagger$（エルミート性）

#### 3.1.2 エネルギー移動 $\hat{H}_{\text{transfer}}$ の構築手順

2分子間の三重項エネルギー移動演算子（$9 \times 9$ 部分空間）：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V \left( |01\rangle\langle 10| + |10\rangle\langle 01| \right)_{ij}
$$

$9 \times 9$ 行列表現（基底: $|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle$）：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V \begin{pmatrix} 0&0&0&0&0&0&0&0&0 \\ 0&0&0&1&0&0&0&0&0 \\ 0&0&0&0&0&0&0&0&0 \\ 0&1&0&0&0&0&0&0&0 \\ 0&0&0&0&0&0&0&0&0 \\ 0&0&0&0&0&0&0&0&0 \\ 0&0&0&0&0&0&0&0&0 \\ 0&0&0&0&0&0&0&0&0 \\ 0&0&0&0&0&0&0&0&0 \end{pmatrix}
$$

```
build_H_transfer(params):
  N = params.N_molecules  # 4
  d = 3
  dim = d^N  # 81
  H_transfer = zeros(dim, dim, dtype=complex)

  # 局所遷移演算子の定義
  # |0⟩⟨1| : S0←T1 transition (3×3)
  op_01 = zeros(3,3); op_01[0,1] = 1.0
  # |1⟩⟨0| : T1←S0 transition (3×3)
  op_10 = zeros(3,3); op_10[1,0] = 1.0

  FOR (i, j) IN params.neighbors:
    # |0⟩_i⟨1| ⊗ |1⟩_j⟨0| + h.c.
    h_pair_fwd = build_two_site_operator(op_01, op_10, i, j, N, d)
    h_pair_bwd = build_two_site_operator(op_10, op_01, i, j, N, d)
    H_transfer += params.V * (h_pair_fwd + h_pair_bwd)

  RETURN H_transfer  # 81×81 エルミート行列
```

```
build_two_site_operator(op_i, op_j, site_i, site_j, N, d):
  """サイト i と j に op_i, op_j を適用する全空間演算子を構築"""
  result = identity(1)
  FOR k IN 0..N-1:
    IF k == site_i:
      result = kron(result, op_i)
    ELIF k == site_j:
      result = kron(result, op_j)
    ELSE:
      result = kron(result, eye(d))
  RETURN result
```

**構築後の検証**: $\hat{H}_{\text{transfer}} = \hat{H}_{\text{transfer}}^\dagger$（エルミート性）

#### 3.1.3 系ハミルトニアンの組立

$$
\hat{H}_{\text{sys}} = \hat{H}_0 + \hat{H}_{\text{transfer}}
$$

```
build_H_sys(params):
  H0 = build_H0(params)
  H_transfer = build_H_transfer(params)
  H_sys = H0 + H_transfer
  # 検証: H_sys がエルミートであること
  assert norm(H_sys - H_sys.conj().T) < 1e-12
  RETURN H_sys
```

### 3.2 Lindblad演算子の構築

#### 3.2.1 演算子の総数と分類

| 過程 | 演算子形式 | 速度定数 | 1分子/ペアあたり | 対象数 | 小計 |
|------|----------|---------|----------------|--------|------|
| TTA | $\|2\rangle_i\langle 1\| \otimes \|0\rangle_j\langle 1\|$ | $\gamma_{\text{TTA}}/2$ | 2 | 3ペア | 6 |
| 蛍光 | $\|0\rangle_i\langle 2\|$ | $\Gamma_{\text{fl}}$ | 1 | 4分子 | 4 |
| 燐光 | $\|0\rangle_i\langle 1\|$ | $\Gamma_{\text{ph}}$ | 1 | 4分子 | 4 |
| 内部転換 | $\|0\rangle_i\langle 2\|$ | $k_{\text{IC}}$ | 1 | 4分子 | 4 |
| ISC S→T | $\|1\rangle_i\langle 2\|$ | $k_{\text{ISC}}^{S\to T}$ | 1 | 4分子 | 4 |
| ISC T→S | $\|0\rangle_i\langle 1\|$ | $k_{\text{ISC}}^{T\to S}$ | 1 | 4分子 | 4 |
| **合計** | | | | | **26** |

#### 3.2.2 各Lindblad演算子の行列定義と構築手順

**TTA Lindblad演算子 $\hat{L}_{\text{TTA},1}^{(ij)}$**:

物理過程: $|T_1\rangle_i|T_1\rangle_j \to |S_1\rangle_i|S_0\rangle_j$

局所表現: $|2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|$

$$
\hat{L}_{\text{TTA},1}^{(ij)} = |20\rangle_{ij}\langle 11|
$$

$9 \times 9$ 行列表現:

$$
\hat{L}_{\text{TTA},1}^{(ij)} = \begin{pmatrix} 0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&1&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0 \end{pmatrix}
$$

（行6列4の要素のみ1、0-indexed）

**TTA Lindblad演算子 $\hat{L}_{\text{TTA},2}^{(ij)}$**:

物理過程: $|T_1\rangle_i|T_1\rangle_j \to |S_0\rangle_i|S_1\rangle_j$

局所表現: $|0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1|$

$$
\hat{L}_{\text{TTA},2}^{(ij)} = |02\rangle_{ij}\langle 11|
$$

$9 \times 9$ 行列表現:

$$
\hat{L}_{\text{TTA},2}^{(ij)} = \begin{pmatrix} 0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&1&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0\\0&0&0&0&0&0&0&0&0 \end{pmatrix}
$$

（行2列4の要素のみ1、0-indexed）

**蛍光 Lindblad演算子 $\hat{L}_{\text{fl}}^{(i)}$**:

物理過程: $|S_1\rangle_i \to |S_0\rangle_i + h\nu_{\text{fl}}$

単一分子表現 ($3 \times 3$):

$$
\hat{L}_{\text{fl}}^{(i)} = |0\rangle\langle 2| = \begin{pmatrix} 0&0&1\\0&0&0\\0&0&0 \end{pmatrix}
$$

**燐光 Lindblad演算子 $\hat{L}_{\text{ph}}^{(i)}$**:

物理過程: $|T_1\rangle_i \to |S_0\rangle_i + h\nu_{\text{ph}}$

単一分子表現 ($3 \times 3$):

$$
\hat{L}_{\text{ph}}^{(i)} = |0\rangle\langle 1| = \begin{pmatrix} 0&1&0\\0&0&0\\0&0&0 \end{pmatrix}
$$

**内部転換 Lindblad演算子 $\hat{L}_{\text{IC}}^{(i)}$**:

物理過程: $|S_1\rangle_i \to |S_0\rangle_i + \text{phonons}$

単一分子表現 ($3 \times 3$): 蛍光と同形式

$$
\hat{L}_{\text{IC}}^{(i)} = |0\rangle\langle 2| = \begin{pmatrix} 0&0&1\\0&0&0\\0&0&0 \end{pmatrix}
$$

**ISC S₁→T₁ Lindblad演算子 $\hat{L}_{\text{ISC}}^{S\to T,(i)}$**:

物理過程: $|S_1\rangle_i \to |T_1\rangle_i$

単一分子表現 ($3 \times 3$):

$$
\hat{L}_{\text{ISC}}^{S\to T,(i)} = |1\rangle\langle 2| = \begin{pmatrix} 0&0&0\\0&0&1\\0&0&0 \end{pmatrix}
$$

**ISC T₁→S₀ Lindblad演算子 $\hat{L}_{\text{ISC}}^{T\to S,(i)}$**:

物理過程: $|T_1\rangle_i \to |S_0\rangle_i + \text{phonons}$

単一分子表現 ($3 \times 3$): 燐光と同形式

$$
\hat{L}_{\text{ISC}}^{T\to S,(i)} = |0\rangle\langle 1| = \begin{pmatrix} 0&1&0\\0&0&0\\0&0&0 \end{pmatrix}
$$

#### 3.2.3 全Lindblad演算子リストの構築フロー

```
build_lindblad_operators(params):
  N = 4; d = 3; dim = 81
  operators = []  # List of (gamma, L_full)

  # 局所遷移演算子の定義 (3×3)
  L_02 = |0⟩⟨2| = zeros(3,3); L_02[0,2] = 1  # S1→S0
  L_01 = |0⟩⟨1| = zeros(3,3); L_01[0,1] = 1  # T1→S0
  L_12 = |1⟩⟨2| = zeros(3,3); L_12[1,2] = 1  # S1→T1
  L_21 = |2⟩⟨1| = zeros(3,3); L_21[2,1] = 1  # T1→S1
  L_10 = |1⟩⟨0| = zeros(3,3); L_10[1,0] = 1  # S0→T1（L_01の転置）

  # ──── TTA（6個） ────
  FOR (i, j) IN params.neighbors:  # 3ペア
    # チャネル1: |T1⟩_i|T1⟩_j → |S1⟩_i|S0⟩_j
    L_TTA1 = build_two_site_operator(L_21, L_01, i, j, N, d)
    operators.append( (params.gamma_TTA / 2, L_TTA1) )

    # チャネル2: |T1⟩_i|T1⟩_j → |S0⟩_i|S1⟩_j
    L_TTA2 = build_two_site_operator(L_01, L_21, i, j, N, d)
    operators.append( (params.gamma_TTA / 2, L_TTA2) )

  # ──── 単一分子散逸（20個） ────
  FOR mol_idx IN 0..3:
    # 蛍光: S1→S0
    L_fl = build_single_site_operator(L_02, mol_idx, N, d)
    operators.append( (params.Gamma_fl, L_fl) )

    # 燐光: T1→S0
    L_ph = build_single_site_operator(L_01, mol_idx, N, d)
    operators.append( (params.Gamma_ph, L_ph) )

    # 内部転換: S1→S0
    L_IC = build_single_site_operator(L_02, mol_idx, N, d)
    operators.append( (params.k_IC, L_IC) )

    # ISC S1→T1
    L_ISC_ST = build_single_site_operator(L_12, mol_idx, N, d)
    operators.append( (params.k_ISC_ST, L_ISC_ST) )

    # ISC T1→S0
    L_ISC_TS = build_single_site_operator(L_01, mol_idx, N, d)
    operators.append( (params.k_ISC_TS, L_ISC_TS) )

  ASSERT len(operators) == 26
  RETURN operators
```

```
build_single_site_operator(op_local, site_idx, N, d):
  """単一サイトの局所演算子を全空間に拡張"""
  result = identity(1)
  FOR k IN 0..N-1:
    IF k == site_idx:
      result = kron(result, op_local)
    ELSE:
      result = kron(result, eye(d))
  RETURN result  # dim×dim 行列
```

**注意**: $\hat{L}_{\text{fl}}^{(i)}$ と $\hat{L}_{\text{IC}}^{(i)}$ は同一の行列形式（$|0\rangle\langle 2|$）を持つが、異なる速度定数（$\Gamma_{\text{fl}}$ vs $k_{\text{IC}}$）を持つ別個の物理過程として扱う。同様に $\hat{L}_{\text{ph}}^{(i)}$ と $\hat{L}_{\text{ISC}}^{T\to S,(i)}$ も同一行列形式（$|0\rangle\langle 1|$）だが別速度定数。

### 3.3 完全なGKSL方程式の構築

系のGKSL方程式の全項展開：

$$
\frac{d\hat{\rho}}{dt} = \underbrace{-\frac{i}{\hbar}[\hat{H}_0 + \hat{H}_{\text{transfer}}, \hat{\rho}]}_{\text{コヒーレント項}}
$$

$$
+ \underbrace{\sum_{(i,j) \in \text{neighbors}} \frac{\gamma_{\text{TTA}}}{2}\left(\mathcal{D}[\hat{L}_{\text{TTA},1}^{(ij)}][\hat{\rho}] + \mathcal{D}[\hat{L}_{\text{TTA},2}^{(ij)}][\hat{\rho}]\right)}_{\text{TTA散逸（6演算子）}}
$$

$$
+ \underbrace{\sum_{i=0}^{3} \Gamma_{\text{fl}} \cdot \mathcal{D}[\hat{L}_{\text{fl}}^{(i)}][\hat{\rho}]}_{\text{蛍光散逸（4演算子）}}
+ \underbrace{\sum_{i=0}^{3} \Gamma_{\text{ph}} \cdot \mathcal{D}[\hat{L}_{\text{ph}}^{(i)}][\hat{\rho}]}_{\text{燐光散逸（4演算子）}}
$$

$$
+ \underbrace{\sum_{i=0}^{3} k_{\text{IC}} \cdot \mathcal{D}[\hat{L}_{\text{IC}}^{(i)}][\hat{\rho}]}_{\text{内部転換散逸（4演算子）}}
+ \underbrace{\sum_{i=0}^{3} k_{\text{ISC}}^{S\to T} \cdot \mathcal{D}[\hat{L}_{\text{ISC}}^{S\to T,(i)}][\hat{\rho}]}_{\text{ISC S→T散逸（4演算子）}}
$$

$$
+ \underbrace{\sum_{i=0}^{3} k_{\text{ISC}}^{T\to S} \cdot \mathcal{D}[\hat{L}_{\text{ISC}}^{T\to S,(i)}][\hat{\rho}]}_{\text{ISC T→S散逸（4演算子）}}
$$

### 3.4 Stinespring Dilationの共通設計

#### 3.4.1 理論的基礎

**Stinespring表現定理**: 任意のCPTP写像 $\mathcal{E}[\hat{\rho}_S]$ は、補助系（ancilla）$E$ を追加したユニタリ演算と部分トレースで実現できる：

$$
\mathcal{E}[\hat{\rho}_S] = \text{Tr}_E\left[\hat{U}_{SE}(\hat{\rho}_S \otimes |0\rangle_E\langle 0|)\hat{U}_{SE}^\dagger\right]
$$

#### 3.4.2 Lindblad演算子からのStinespringユニタリ構築

単一Lindblad演算子 $\hat{L}$ に対する生成子：

$$
\hat{G} = \hat{L} \otimes |1\rangle_E\langle 0| + \hat{L}^\dagger \otimes |0\rangle_E\langle 1|
$$

Stinespringユニタリ：

$$
\hat{U}(\theta) = e^{-i\theta \hat{G}}
$$

パラメータ：

$$
\theta = \sqrt{\gamma \Delta t}
$$

近似精度：

$$
\text{Tr}_E[\hat{U}(\theta)(\hat{\rho}\otimes|0\rangle\langle 0|)\hat{U}^\dagger(\theta)] = \hat{\rho} + \gamma\Delta t \cdot \mathcal{D}[\hat{L}][\hat{\rho}] + O(\gamma^2 \Delta t^2)
$$

#### 3.4.3 Stinespringユニタリの具体的構築手順

```
build_stinespring_unitary(L_operator, gamma, dt):
  """
  入力:
    L_operator: dim_sys × dim_sys のLindblad演算子
    gamma: 散逸速度定数
    dt: 時間刻み
  出力:
    U: (dim_sys * 2) × (dim_sys * 2) のユニタリ行列
  """
  dim_sys = L_operator.shape[0]
  theta = sqrt(gamma * dt)

  # 補助系（qubit）の演算子
  sigma_minus = array([[0, 0], [1, 0]])  # |1⟩⟨0|
  sigma_plus  = array([[0, 1], [0, 0]])  # |0⟩⟨1|

  # 生成子 G の構築
  G = kron(L_operator, sigma_minus) + kron(L_operator.conj().T, sigma_plus)
  # G は (dim_sys * 2) × (dim_sys * 2) のエルミート行列

  # 検証: G がエルミートであること
  assert norm(G - G.conj().T) < 1e-12

  # ユニタリの計算
  U = expm(-1j * theta * G)

  # 検証: U がユニタリであること
  assert norm(U @ U.conj().T - eye(dim_sys * 2)) < 1e-10

  RETURN U
```

#### 3.4.4 Stinespring近似の有効条件

Stinespring近似が有効であるための条件：

$$
\gamma_{\max} \cdot \Delta t \ll 1
$$

本仕様のパラメータでの確認：

$$
\gamma_{\max} = \gamma_{\text{TTA}} = 0.05 \text{ eV/ℏ}, \quad \Delta t = 1 \text{ fs}
$$

$$
\gamma_{\max} \cdot \Delta t / \hbar = 0.05 \times 1 / 0.658 \approx 0.076 \ll 1 \quad \checkmark
$$

累積誤差（$N_{\text{steps}}$ ステップ後）：

$$
\epsilon_{\text{total}} \sim N_{\text{steps}} \cdot \gamma^2 (\Delta t)^2 = \gamma^2 T \Delta t
$$

本仕様では：

$$
\epsilon_{\text{total}} \sim (0.05)^2 \times (100/0.658) \times (1/0.658) \approx 0.058
$$


---

## 4. シナリオ1: 古典GKSL・ボソン無し 詳細設計

### 4.1 設計概要

密度行列のGKSL-Lindblad時間発展を超演算子形式で直接計算する古典的シミュレータ。他の全シナリオの検証基準となる。

### 4.2 超演算子の構築設計

#### 4.2.1 ベクトル化の設計

$81 \times 81$ の密度行列 $\hat{\rho}$ を $6561$ 次元の列ベクトル $|\hat{\rho}\rangle\rangle$ に変換する。

変換規則（列優先、column-major / Fortranスタイル）：

$$
\text{vec}(\hat{\rho})_{i + j \cdot d} = \hat{\rho}_{ij}, \quad d = 81
$$

```
vectorize(rho):
  RETURN rho.flatten(order='F')  # column-major

unvectorize(rho_vec, d=81):
  RETURN rho_vec.reshape((d, d), order='F')
```

#### 4.2.2 ハミルトニアン超演算子の構築

$$
\mathcal{L}_H = -\frac{i}{\hbar}\left(\hat{H} \otimes \hat{I}_d - \hat{I}_d \otimes \hat{H}^T\right)
$$

ここで $\otimes$ は Kronecker積、$d = 81$。

```
build_hamiltonian_superoperator(H_sys, hbar):
  d = H_sys.shape[0]  # 81
  I_d = eye(d)
  L_H = -1j / hbar * (kron(H_sys, I_d) - kron(I_d, H_sys.T))
  # L_H は 6561×6561 の複素行列
  RETURN L_H
```

**検証**: $\text{Tr}[\mathcal{L}_H] = 0$（トレース保存を暗示）

#### 4.2.3 Lindblad超演算子の構築

各Lindblad演算子 $\hat{L}_\alpha$（$81 \times 81$）と速度定数 $\gamma_\alpha$ に対して：

$$
\mathcal{L}_\alpha = \gamma_\alpha \left(\hat{L}_\alpha \otimes \bar{\hat{L}}_\alpha - \frac{1}{2}\hat{L}_\alpha^\dagger\hat{L}_\alpha \otimes \hat{I}_d - \frac{1}{2}\hat{I}_d \otimes (\hat{L}_\alpha^\dagger\hat{L}_\alpha)^T\right)
$$

ここで $\bar{\hat{L}}_\alpha$ は $\hat{L}_\alpha$ の要素ごとの複素共役（転置なし）。

```
build_lindblad_superoperator_single(gamma, L, d):
  I_d = eye(d)
  LdL = L.conj().T @ L  # L†L
  L_super = gamma * (
      kron(L, L.conj())           # L ⊗ L̄
    - 0.5 * kron(LdL, I_d)        # -½ L†L ⊗ I
    - 0.5 * kron(I_d, LdL.T)      # -½ I ⊗ (L†L)ᵀ
  )
  RETURN L_super  # d²×d² 行列
```

#### 4.2.4 全超演算子の構築

$$
\mathcal{L}_{\text{super}} = \mathcal{L}_H + \sum_{\alpha=1}^{26} \mathcal{L}_\alpha
$$

```
build_full_superoperator(params):
  H_sys = build_H_sys(params)
  lindblad_ops = build_lindblad_operators(params)
  d = 81

  # ハミルトニアン超演算子
  L_super = build_hamiltonian_superoperator(H_sys, params.hbar)

  # Lindblad超演算子の加算
  FOR (gamma, L) IN lindblad_ops:
    L_super += build_lindblad_superoperator_single(gamma, L, d)

  RETURN L_super  # 6561×6561 複素行列
```

**メモリ見積もり**: $6561^2 \times 16$ bytes（complex128）$\approx 689$ MB

**疎行列最適化**: 超演算子は高度に疎であるため、`scipy.sparse.csr_matrix` を使用してメモリと計算時間を大幅に削減可能。

### 4.3 時間発展の設計

#### 4.3.1 方法A: ODEソルバー（推奨）

```
simulate_ode(L_super, rho_0, T_total, N_steps, method='RK45'):
  d = 81
  rho_0_vec = vectorize(rho_0)
  t_eval = linspace(0, T_total, N_steps + 1)

  def rhs(t, rho_vec):
    RETURN L_super @ rho_vec

  sol = solve_ivp(
    rhs,
    [0, T_total],
    rho_0_vec,
    method=method,     # 'RK45'（非スティッフ）or 'BDF'（スティッフ）
    t_eval=t_eval,
    rtol=1e-10,
    atol=1e-12
  )

  IF NOT sol.success:
    raise RuntimeError(f"ODE solver failed: {sol.message}")
    # Fallbackは禁止

  # 各時刻の密度行列を復元
  results = []
  FOR k IN 0..N_steps:
    rho_k = unvectorize(sol.y[:, k], d)
    results.append(rho_k)

  RETURN results, sol.nfev
```

**方法の選択基準**:
- 速度定数の比が $10^4$ 以上（例: $\gamma_{\text{TTA}} / k_{\text{ISC}}^{T\to S} = 5000$）→ スティッフ系 → `method='BDF'`
- 速度定数の比が $10^4$ 未満 → `method='RK45'`

#### 4.3.2 方法B: 行列指数関数法

```
simulate_expm(L_super, rho_0, T_total, N_steps):
  d = 81
  dt = T_total / N_steps
  rho_vec = vectorize(rho_0)

  # ステップ行列を1回だけ計算
  step_matrix = expm(L_super * dt)  # 6561×6561

  results = [rho_0.copy()]
  FOR step IN 1..N_steps:
    rho_vec = step_matrix @ rho_vec
    rho = unvectorize(rho_vec, d)
    results.append(rho)

  RETURN results
```

計算量: $O(d^6)$（expmの計算）+ $O(N_{\text{steps}} \cdot d^4)$（行列-ベクトル積）

#### 4.3.3 方法C: 疎行列＋Krylov法

```
simulate_sparse(L_super_sparse, rho_0, T_total, N_steps):
  d = 81
  dt = T_total / N_steps
  rho_vec = vectorize(rho_0)

  results = [rho_0.copy()]
  FOR step IN 1..N_steps:
    rho_vec = expm_multiply(L_super_sparse * dt, rho_vec)
    rho = unvectorize(rho_vec, d)
    results.append(rho)

  RETURN results
```

### 4.4 ClassicalGKSLSimulator 詳細設計

#### 4.4.1 フローチャート

```
┌─────────────────────────┐
│   simulate() 開始        │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 1. パラメータ検証         │
│    params.validate()     │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 2. ハミルトニアン構築     │
│    H_sys = build_H_sys() │
│    検証: H = H†           │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 3. Lindblad演算子構築     │
│    ops = build_lindblad_ │
│          operators()     │
│    検証: len(ops) == 26   │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 4. 超演算子構築           │
│    L_super = build_full_ │
│              superop()   │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 5. 初期密度行列構築       │
│    rho_0 = prepare_      │
│    initial_density_      │
│    matrix()              │
│    検証: Tr[ρ₀]=1,       │
│          ρ₀=ρ₀†, ρ₀²=ρ₀ │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 6. 時間発展の実行         │
│    FOR step IN 1..N:     │
│      ρ(t+dt) = evolve(ρ) │
│      ├─ 検証: Tr[ρ]=1    │
│      ├─ 検証: λₖ≥-ε     │
│      ├─ 計算: populations│
│      ├─ 計算: entropy    │
│      ├─ 計算: purity     │
│      └─ 計算: coherences │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│ 7. 結果の集約と出力       │
│    RETURN GKSLResult     │
└─────────────────────────┘
```

#### 4.4.2 密度行列の物理性検証（各ステップ）

```
validate_density_matrix(rho, step, tolerance=1e-8):
  # 1. トレース保存
  tr = real(trace(rho))
  IF abs(1 - tr) > tolerance:
    raise PhysicsViolationError(
      f"Step {step}: Trace = {tr}, deviation = {abs(1-tr)}")

  # 2. エルミート性
  herm_err = norm(rho - rho.conj().T, 'fro')
  IF herm_err > tolerance:
    raise PhysicsViolationError(
      f"Step {step}: Hermiticity violation = {herm_err}")

  # 3. 正定値性
  eigenvalues = eigvalsh(rho)
  min_eig = min(eigenvalues)
  IF min_eig < -tolerance:
    raise PhysicsViolationError(
      f"Step {step}: Negative eigenvalue = {min_eig}")

  # 4. 粒子数保存
  pops = calculate_populations(rho)
  N_total = pops['N_S0'] + pops['N_T1'] + pops['N_S1']
  IF abs(N_total - N_molecules) > tolerance:
    raise PhysicsViolationError(
      f"Step {step}: Particle number = {N_total}")
```

**重要**: 検証に失敗した場合は `PhysicsViolationError` を raise する。Fallback処理（「適当な値」への置き換えなど）は厳禁。

---

## 5. シナリオ2: 古典GKSL・ボソン有り 詳細設計

### 5.1 設計概要

電子系にフォノン（格子振動）モードを明示的に結合させたモデル。Holstein型電子-フォノン結合を採用。

### 5.2 拡張ヒルベルト空間の設計

#### 5.2.1 空間構造

$$
\mathcal{H}_{\text{total}} = \mathcal{H}_{\text{el}} \otimes \mathcal{H}_{\text{phonon}}
$$

$$
\mathcal{H}_{\text{el}} = \bigotimes_{i=0}^{3} \mathbb{C}^3 \quad (\dim = 81)
$$

$$
\mathcal{H}_{\text{phonon}} = \bigotimes_{i=0}^{3} \mathbb{C}^{n_{\max}+1} \quad (\dim = (n_{\max}+1)^4)
$$

$n_{\max} = 2$ の場合の次元：

| 部分空間 | 次元 |
|---------|------|
| 電子系 | $81$ |
| フォノン系 | $3^4 = 81$ |
| 全系 | $81 \times 81 = 6561$ |
| 密度行列 | $6561 \times 6561 = 43046721$ |

#### 5.2.2 基底の順序

全空間の基底状態は以下の順序で列挙する：

$$
|e_0, e_1, e_2, e_3, n_0, n_1, n_2, n_3\rangle
$$

ここで $e_i \in \{0, 1, 2\}$ は分子 $i$ の電子状態、$n_i \in \{0, 1, \ldots, n_{\max}\}$ はフォノン数。

インデックスの計算：

$$
\text{idx} = \left(\sum_{i=0}^{3} e_i \cdot 3^{3-i}\right) \cdot (n_{\max}+1)^4 + \sum_{i=0}^{3} n_i \cdot (n_{\max}+1)^{3-i}
$$

### 5.3 拡張ハミルトニアンの構築設計

#### 5.3.1 電子系ハミルトニアン（フォノン空間への拡張）

$$
\hat{H}_{\text{el}}^{\text{ext}} = \hat{H}_{\text{sys}} \otimes \hat{I}_{\text{phonon}}
$$

```
build_H_el_extended(H_sys, dim_phonon):
  RETURN kron(H_sys, eye(dim_phonon))
```

#### 5.3.2 フォノンハミルトニアン

$$
\hat{H}_{\text{phonon}} = \hat{I}_{\text{el}} \otimes \sum_{i=0}^{3} \hbar\omega_{\text{ph}} \hat{n}_i
$$

フォノン数演算子の行列表現（$(n_{\max}+1) \times (n_{\max}+1)$）：

$$
\hat{n} = \hat{a}^\dagger\hat{a} = \text{diag}(0, 1, 2, \ldots, n_{\max})
$$

フォノン消滅演算子（$(n_{\max}+1) \times (n_{\max}+1)$）：

$$
\hat{a} = \begin{pmatrix} 0 & \sqrt{1} & 0 & \cdots \\ 0 & 0 & \sqrt{2} & \cdots \\ \vdots & & & \ddots \\ 0 & \cdots & 0 & \sqrt{n_{\max}} \\ 0 & \cdots & & 0 \end{pmatrix}
$$

```
build_phonon_operators(n_max):
  d_ph = n_max + 1
  a = zeros(d_ph, d_ph)
  FOR n IN 0..n_max-1:
    a[n, n+1] = sqrt(n + 1)
  a_dag = a.conj().T
  n_op = a_dag @ a
  RETURN a, a_dag, n_op
```

```
build_H_phonon(params, dim_el):
  d_ph = params.n_max + 1
  dim_ph = d_ph^4
  a, a_dag, n_op = build_phonon_operators(params.n_max)

  H_phonon_local = zeros(dim_ph, dim_ph)
  FOR i IN 0..3:
    n_i = build_single_site_operator(n_op, i, N=4, d=d_ph)
    H_phonon_local += params.hbar * params.omega_ph * n_i

  RETURN kron(eye(dim_el), H_phonon_local)
```

#### 5.3.3 電子-フォノン結合（Holstein型）

$$
\hat{H}_{e\text{-ph}} = g \sum_{i=0}^{3} |1\rangle_i\langle 1| \otimes (\hat{a}_i + \hat{a}_i^\dagger)
$$

```
build_H_eph(params):
  d_el = 3; d_ph = params.n_max + 1
  dim_el = d_el^4; dim_ph = d_ph^4

  # 電子系の|1⟩⟨1|射影演算子 (3×3)
  proj_T1 = diag([0, 1, 0])

  # フォノンの (a + a†) 演算子
  a, a_dag, _ = build_phonon_operators(params.n_max)
  x_op = a + a_dag  # 変位演算子

  H_eph = zeros(dim_el * dim_ph, dim_el * dim_ph, dtype=complex)

  FOR i IN 0..3:
    # 電子部分: I⊗...⊗|1⟩⟨1|⊗...⊗I
    el_part = build_single_site_operator(proj_T1, i, N=4, d=d_el)

    # フォノン部分: I⊗...⊗(a+a†)⊗...⊗I
    ph_part = build_single_site_operator(x_op, i, N=4, d=d_ph)

    # テンソル積
    H_eph += params.g_eph * kron(el_part, ph_part)

  RETURN H_eph
```

#### 5.3.4 全ハミルトニアンの組立

$$
\hat{H}_{\text{total}} = \hat{H}_{\text{el}}^{\text{ext}} + \hat{H}_{\text{phonon}} + \hat{H}_{e\text{-ph}}
$$

```
build_H_total_boson(params):
  H_sys = build_H_sys(params)
  dim_el = 81
  d_ph = params.n_max + 1
  dim_ph = d_ph^4

  H_el_ext = build_H_el_extended(H_sys, dim_ph)
  H_phonon = build_H_phonon(params, dim_el)
  H_eph = build_H_eph(params)

  H_total = H_el_ext + H_phonon + H_eph
  # 検証: H_total がエルミート
  assert norm(H_total - H_total.conj().T) < 1e-10
  RETURN H_total
```

### 5.4 Lindblad演算子の拡張

Lindblad演算子は電子系にのみ作用するが、フォノン空間への拡張が必要：

$$
\hat{L}_\alpha^{\text{ext}} = \hat{L}_\alpha^{\text{el}} \otimes \hat{I}_{\text{phonon}}
$$

```
extend_lindblad_operators(lindblad_ops_el, dim_phonon):
  extended_ops = []
  FOR (gamma, L_el) IN lindblad_ops_el:
    L_ext = kron(L_el, eye(dim_phonon))
    extended_ops.append((gamma, L_ext))
  RETURN extended_ops
```

### 5.5 個体数計算（部分トレース）

電子系の縮約密度行列を計算し、そこから個体数を得る：

$$
\hat{\rho}_{\text{el}}(t) = \text{Tr}_{\text{phonon}}[\hat{\rho}_{\text{total}}(t)]
$$

```
partial_trace_phonon(rho_total, dim_el=81, dim_ph):
  """全密度行列からフォノン部分をトレースアウト"""
  rho_el = zeros(dim_el, dim_el, dtype=complex)
  FOR i IN 0..dim_el-1:
    FOR j IN 0..dim_el-1:
      FOR k IN 0..dim_ph-1:
        rho_el[i, j] += rho_total[i*dim_ph + k, j*dim_ph + k]
  RETURN rho_el
```

### 5.6 数値手法の選択

$n_{\max} = 2$ の場合、$d_{\text{total}} = 6561$。超演算子サイズは $6561^2 \approx 4.3 \times 10^7$ であるため、明示的な超演算子行列の構築は困難。

**推奨手法**: 暗黙的ODE積分（超演算子の行列-ベクトル積を関数として定義）

```
simulate_boson(params):
  H_total = build_H_total_boson(params)
  ext_ops = extend_lindblad_operators(...)
  d_total = H_total.shape[0]

  def lindblad_rhs(t, rho_vec):
    rho = rho_vec.reshape((d_total, d_total), order='F')
    drho = -1j / params.hbar * (H_total @ rho - rho @ H_total)
    FOR (gamma, L) IN ext_ops:
      LdL = L.conj().T @ L
      drho += gamma * (L @ rho @ L.conj().T
                       - 0.5 * LdL @ rho
                       - 0.5 * rho @ LdL)
    RETURN drho.flatten(order='F')

  # BDF法（スティッフ系向け）を使用
  sol = solve_ivp(lindblad_rhs, [0, T_total], rho_0_vec,
                  method='BDF', t_eval=t_eval, rtol=1e-8, atol=1e-10)
```

### 5.7 ClassicalGKSLBosonSimulator フローチャート

```
┌──────────────────────────────┐
│    simulate() 開始            │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│ 1. パラメータ検証             │
│    + ボソンパラメータの検証    │
│    (n_max >= 1, omega_ph > 0, │
│     g_eph >= 0)               │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│ 2. 拡張ハミルトニアン構築     │
│    H_total = H_el⊗I + I⊗H_ph │
│              + H_eph          │
│    検証: H_total = H_total†   │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│ 3. 拡張Lindblad演算子構築     │
│    L_ext = L_el ⊗ I_phonon   │
│    26個の演算子を拡張          │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│ 4. 初期状態構築               │
│    ρ_0 = |ψ_el⟩⟨ψ_el| ⊗     │
│          |0000⟩⟨0000|_phonon │
│    （フォノン真空状態、T=0）   │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│ 5. ODE積分（暗黙的rhs関数）   │
│    method='BDF'               │
│    各ステップで:               │
│    ├─ 部分トレース→ ρ_el      │
│    ├─ 個体数計算               │
│    ├─ エントロピー計算         │
│    └─ 物理性検証               │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│ 6. 結果の集約と出力           │
│    + フォノン占有数も記録      │
└──────────────────────────────┘
```

---

## 6. シナリオ3: Qubit GKSL・ボソン無し 詳細設計

### 6.1 設計概要

現行ノートブックのQubit実装を拡張し、Stinespring dilationを用いてLindblad散逸項を量子回路に実装する。

### 6.2 Qubitエンコーディング設計

各分子 $i$ を2個のqubit $(q_{2i}, q_{2i+1})$ で表現：

$$
|S_0\rangle_i \leftrightarrow |00\rangle_{2i,2i+1}, \quad |T_1\rangle_i \leftrightarrow |01\rangle_{2i,2i+1}, \quad |S_1\rangle_i \leftrightarrow |10\rangle_{2i,2i+1}
$$

禁止状態: $|11\rangle_{2i,2i+1}$（物理的意味なし）

| リソース | 数 | 説明 |
|---------|-----|------|
| 系 qubit | 8 | 4分子 × 2 qubit |
| Ancilla qubit | 26 | 26個のLindblad演算子 |
| **合計** | **34** | |

### 6.3 量子回路の構成設計

#### 6.3.1 1 Trotterステップの回路構造

2次対称Trotter分解：

$$
e^{\mathcal{L}_{\text{GKSL}} \Delta t} \approx e^{\mathcal{L}_H \Delta t/2} \cdot e^{\mathcal{L}_{\text{diss}} \Delta t} \cdot e^{\mathcal{L}_H \Delta t/2}
$$

```
1トロッターステップの回路構成:

┌────────────────────────────────────────────────────┐
│ ユニタリ前半 (Δt/2)                                 │
│  ├─ H0: 各分子の対角位相 (Rz ゲート × 4分子)         │
│  └─ H_transfer: 2分子エネルギー移動 (ペア (0,1),     │
│                  (1,2), (2,3))                       │
├────────────────────────────────────────────────────┤
│ 散逸ステップ (Δt)                                    │
│  ├─ TTA Stinespring: ペア(0,1)×2ch, (1,2)×2ch,     │
│  │                    (2,3)×2ch → 6 ancilla          │
│  ├─ 蛍光 Stinespring: 分子 0,1,2,3 → 4 ancilla      │
│  ├─ 燐光 Stinespring: 分子 0,1,2,3 → 4 ancilla      │
│  ├─ IC Stinespring: 分子 0,1,2,3 → 4 ancilla        │
│  ├─ ISC S→T Stinespring: 分子 0,1,2,3 → 4 ancilla   │
│  └─ ISC T→S Stinespring: 分子 0,1,2,3 → 4 ancilla   │
├────────────────────────────────────────────────────┤
│ ユニタリ後半 (Δt/2)                                  │
│  ├─ H_transfer: ペア (2,3), (1,2), (0,1) [逆順]     │
│  └─ H0: 分子 3,2,1,0 [逆順]                         │
└────────────────────────────────────────────────────┘
```

#### 6.3.2 各Lindblad演算子のQubit回路設計

**蛍光 $\hat{L}_{\text{fl}}^{(i)} = |0\rangle_i\langle 2|$** (Qubit: $|00\rangle\langle 10|$)

Stinespring生成子：

$$
\hat{G}_{\text{fl}} = |00\rangle\langle 10|_{(q_0,q_1)} \otimes |1\rangle\langle 0|_E + |10\rangle\langle 00|_{(q_0,q_1)} \otimes |0\rangle\langle 1|_E
$$

回路手順：
1. $q_0 = 1, q_1 = 0$（$|10\rangle$ 状態 = $S_1$）を検出
2. 条件付きで ancilla $q_E$ に $R_Y(2\arcsin(\sqrt{\Gamma_{\text{fl}} \Delta t}))$ 適用
3. ancilla が $|1\rangle$ なら系を $|10\rangle \to |00\rangle$ に遷移

**燐光 $\hat{L}_{\text{ph}}^{(i)} = |0\rangle_i\langle 1|$** (Qubit: $|00\rangle\langle 01|$)

Stinespring生成子：

$$
\hat{G}_{\text{ph}} = |00\rangle\langle 01|_{(q_0,q_1)} \otimes |1\rangle\langle 0|_E + |01\rangle\langle 00|_{(q_0,q_1)} \otimes |0\rangle\langle 1|_E
$$

回路手順：
1. $q_1 = 1, q_0 = 0$（$|01\rangle$ 状態 = $T_1$）を検出
2. 条件付き ancilla 回転
3. ancilla が $|1\rangle$ なら $|01\rangle \to |00\rangle$

**ISC S₁→T₁ $\hat{L}_{\text{ISC}}^{S\to T,(i)} = |1\rangle_i\langle 2|$** (Qubit: $|01\rangle\langle 10|$)

Stinespring生成子：

$$
\hat{G}_{\text{ISC}} = |01\rangle\langle 10|_{(q_0,q_1)} \otimes |1\rangle\langle 0|_E + |10\rangle\langle 01|_{(q_0,q_1)} \otimes |0\rangle\langle 1|_E
$$

回路手順：
1. $|10\rangle$（$S_1$）を検出
2. 条件付き ancilla 回転
3. ancilla が $|1\rangle$ なら $|10\rangle \to |01\rangle$（両qubit反転）

**TTA $\hat{L}_{\text{TTA},1}^{(ij)} = |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|$** (Qubit: $|10,00\rangle\langle 01,01|$)

4-qubit系 + 1 ancilla の操作。

Stinespring生成子：

$$
\hat{G}_{\text{TTA},1} = |10,00\rangle\langle 01,01| \otimes |1\rangle\langle 0|_E + |01,01\rangle\langle 10,00| \otimes |0\rangle\langle 1|_E
$$

回路手順：
1. 4-qubit状態 $|01,01\rangle$（$T_1, T_1$）を検出
2. 条件付き ancilla 回転（$\theta = \sqrt{(\gamma_{\text{TTA}}/2) \Delta t}$）
3. ancilla が $|1\rangle$ なら $|01,01\rangle \to |10,00\rangle$

#### 6.3.3 ゲート分解戦略

2段階アプローチ：

**レベル1 (UnitaryGate表現)**:
- Stinespringユニタリ行列を直接 `UnitaryGate` としてQiskit回路に挿入
- ゲート数: 少ない（検証用）

**レベル2 (基本ゲート分解)**:
- `UnitaryGate` → KAK分解 → CNOT + Rz + Ry + Rx
- ゲート数: 多い（実機実行用）

### 6.4 密度行列の再構成

Statevectorシミュレータ使用時、ancillaを含む全系の状態ベクトルから系の密度行列を部分トレースで取得：

$$
\hat{\rho}_{\text{sys}} = \text{Tr}_{\text{ancilla}}[|\Psi_{\text{total}}\rangle\langle\Psi_{\text{total}}|]
$$

```
reconstruct_density_matrix(statevector, n_sys_qubits=8, n_ancilla=26):
  # Qiskitのpartial_traceを使用
  from qiskit.quantum_info import Statevector, partial_trace
  sv = Statevector(statevector)
  ancilla_indices = list(range(n_sys_qubits, n_sys_qubits + n_ancilla))
  rho_sys = partial_trace(sv, ancilla_indices)
  RETURN rho_sys.data  # 2^8 × 2^8 = 256×256 行列

  # 注: 物理的部分空間は3^4=81次元。
  # 256×256行列から81×81の物理的密度行列を抽出する必要がある。
```

### 6.5 禁止状態への遷移監視

各ステップで禁止状態（$|11\rangle$）への遷移確率を検証：

$$
P_{\text{forbidden}}(t) = 1 - \text{Tr}[\hat{P}_{\text{phys}} \hat{\rho}_{\text{sys}}(t)] < 10^{-8}
$$

```
check_forbidden_states(rho_qubit, N=4):
  P_phys = construct_physical_projector(N)
  P_forbidden = 1 - real(trace(P_phys @ rho_qubit))
  IF P_forbidden > 1e-8:
    raise PhysicsViolationError(
      f"Forbidden state leakage: {P_forbidden}")
```

### 6.6 QubitGKSLSimulator フローチャート

```
┌──────────────────────────────────┐
│     simulate() 開始               │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ 1. パラメータ検証                 │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ 2. 量子回路の構築                 │
│  ├─ qubit レジスタ割当て          │
│  │   (sys: 8, ancilla: 26)       │
│  ├─ 初期状態準備                  │
│  │   |ψ₀⟩ = |01,00,00,01⟩_sys   │
│  │         ⊗ |0...0⟩_ancilla     │
│  └─ N_stepsループ:               │
│      ├─ build_unitary_step(dt/2) │
│      ├─ build_lindblad_step(dt)  │
│      ├─ ancilla リセット（|0⟩）   │
│      └─ build_unitary_step(dt/2) │
│          [逆順]                   │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ 3. シミュレーション実行           │
│  ├─ Statevector: 密度行列再構成   │
│  └─ Shot-based: 測定統計から推定  │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ 4. 各ステップの後処理             │
│  ├─ 部分トレース（ancilla除去）   │
│  ├─ 禁止状態遷移チェック          │
│  ├─ 個体数計算                    │
│  ├─ エントロピー計算              │
│  └─ 純度計算                      │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ 5. 結果出力                       │
│  + 回路情報（ゲート数、深さ）      │
└──────────────────────────────────┘
```


---

## 7. シナリオ4: Qubit GKSL・ボソン有り 詳細設計

### 7.1 設計概要

シナリオ3にフォノンモードのqubit表現を追加。電子-フォノン結合を量子回路で実現する。

### 7.2 フォノンモードのQubitエンコーディング

バイナリエンコーディング: $n_{\max} = 2$ の場合

$$
|0\rangle_{\text{Fock}} \leftrightarrow |00\rangle, \quad |1\rangle_{\text{Fock}} \leftrightarrow |01\rangle, \quad |2\rangle_{\text{Fock}} \leftrightarrow |10\rangle
$$

必要qubit数: $\lceil\log_2(n_{\max}+1)\rceil = 2$ qubit/フォノンモード

禁止状態: $|11\rangle$（Fock状態 $|3\rangle$ に対応するが $n_{\max} = 2$ で切断）

### 7.3 量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 電子系 qubit | 8 | 4分子 × 2 qubit |
| フォノン qubit | 8 | 4分子 × 2 qubit ($n_{\max} = 2$) |
| Ancilla qubit (Lindblad) | 26 | 26個のLindblad演算子 |
| **合計** | **42** | |

### 7.4 電子-フォノン結合の量子回路設計

Holstein型結合 $g(\hat{a}+\hat{a}^\dagger)|1\rangle\langle 1|$ の回路実装：

```
電子-フォノン結合の1ステップ回路:

1. 電子系の |T₁⟩ = |01⟩ を制御条件として検出
2. フォノン qubit に対して条件付きインクリメント/デクリメント回路を適用

回路の概略:
  IF electron_state == |01⟩ THEN:
    フォノンqubitに exp(-ig·Δt·(a+a†)/ℏ) を適用
```

フォノン昇降演算子のqubit表現（$n_{\max} = 2$の場合の $3 \times 3$ → $4 \times 4$ 拡張）：

$$
\hat{a}_{\text{qubit}} = \begin{pmatrix} 0&1&0&0\\0&0&\sqrt{2}&0\\0&0&0&0\\0&0&0&0 \end{pmatrix}
$$

これを2-qubitゲートに分解する必要がある。

### 7.5 QubitGKSLBosonSimulator フローチャート

```
┌────────────────────────────────────┐
│     simulate() 開始                 │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 1. qubit割当て                      │
│   電子: q0-q7 (8 qubit)            │
│   フォノン: q8-q15 (8 qubit)       │
│   Ancilla: q16-q41 (26 qubit)      │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 2. 1 Trotterステップの構築          │
│  ├─ ユニタリ前半:                   │
│  │   H0(el) + H_transfer + H_phonon│
│  │   + H_eph                        │
│  ├─ Lindblad散逸ステップ            │
│  │   (電子系qubit + ancilla qubitのみ│
│  │    フォノンqubitには作用しない)    │
│  └─ ユニタリ後半                    │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 3. N_stepsループ実行               │
│    各ステップ後に部分トレース       │
│    （ancilla + フォノンを除去）      │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 4. 結果出力                         │
└────────────────────────────────────┘
```

---

## 8. シナリオ5: Qudit GKSL・ボソン無し 詳細設計

### 8.1 設計概要

MQT-Quditsフレームワークを使用し、各分子を1 qutrit（$d = 3$）で直接表現。Stinespring dilationで散逸項を実装する。Qubitエンコーディングと比較して禁止状態が存在しないことが最大の利点。

### 8.2 Qutritエンコーディング

$$
|S_0\rangle_i \leftrightarrow |0\rangle_i, \quad |T_1\rangle_i \leftrightarrow |1\rangle_i, \quad |S_1\rangle_i \leftrightarrow |2\rangle_i
$$

状態空間: $3^4 = 81$（全て物理状態、禁止状態なし）

### 8.3 Stinespring演算子の具体的行列設計

#### 8.3.1 蛍光のStinespringユニタリ

$\hat{L}_{\text{fl}}^{(i)} = |0\rangle_i\langle 2|$、系qutrit（$d=3$）+ ancilla qubit（$d=2$）

基底: $\{|0,0\rangle, |0,1\rangle, |1,0\rangle, |1,1\rangle, |2,0\rangle, |2,1\rangle\}$

生成子：

$$
\hat{G}_{\text{fl}} = \begin{pmatrix} 0&0&0&0&0&1\\0&0&0&0&0&0\\0&0&0&0&0&0\\0&0&0&0&0&0\\0&0&0&0&0&0\\1&0&0&0&0&0 \end{pmatrix}
$$

ユニタリ（$\theta = \sqrt{\Gamma_{\text{fl}} \Delta t}$）：

$$
\hat{U}_{\text{fl}}(\theta) = \begin{pmatrix} \cos\theta&0&0&0&0&-i\sin\theta\\0&1&0&0&0&0\\0&0&1&0&0&0\\0&0&0&1&0&0\\0&0&0&0&1&0\\-i\sin\theta&0&0&0&0&\cos\theta \end{pmatrix}
$$

**検証**: $\hat{U}_{\text{fl}}^\dagger \hat{U}_{\text{fl}} = \hat{I}_6$

#### 8.3.2 燐光のStinespringユニタリ

$\hat{L}_{\text{ph}}^{(i)} = |0\rangle_i\langle 1|$

基底: $\{|0,0\rangle, |0,1\rangle, |1,0\rangle, |1,1\rangle, |2,0\rangle, |2,1\rangle\}$

$$
\hat{U}_{\text{ph}}(\theta) = \begin{pmatrix} \cos\theta&0&0&-i\sin\theta&0&0\\0&1&0&0&0&0\\0&0&1&0&0&0\\-i\sin\theta&0&0&\cos\theta&0&0\\0&0&0&0&1&0\\0&0&0&0&0&1 \end{pmatrix}
$$

$\theta = \sqrt{\Gamma_{\text{ph}} \Delta t}$

#### 8.3.3 ISC S₁→T₁ のStinespringユニタリ

$\hat{L}_{\text{ISC}}^{S\to T,(i)} = |1\rangle_i\langle 2|$

$$
\hat{U}_{\text{ISC}}(\theta) = \begin{pmatrix} 1&0&0&0&0&0\\0&1&0&0&0&0\\0&0&\cos\theta&0&0&-i\sin\theta\\0&0&0&1&0&0\\0&0&0&0&1&0\\0&0&-i\sin\theta&0&0&\cos\theta \end{pmatrix}
$$

$\theta = \sqrt{k_{\text{ISC}}^{S\to T} \Delta t}$

#### 8.3.4 TTA チャネル1のStinespringユニタリ

$\hat{L}_{\text{TTA},1}^{(ij)} = |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|$

2 qutrit（$d=3$）+ 1 ancilla qubit（$d=2$） → $3 \times 3 \times 2 = 18$ 次元

基底: $\{|mn,e\rangle\}$ where $m,n \in \{0,1,2\}$, $e \in \{0,1\}$

ユニタリは18×18行列で、$|20,0\rangle \leftrightarrow |11,1\rangle$ の部分空間でのみ非自明な回転を行う：

$$
\hat{U}_{\text{TTA},1}(\theta) = \hat{I}_{18} + (\cos\theta - 1)(|20,0\rangle\langle 20,0| + |11,1\rangle\langle 11,1|) - i\sin\theta(|20,0\rangle\langle 11,1| + |11,1\rangle\langle 20,0|)
$$

$\theta = \sqrt{(\gamma_{\text{TTA}}/2) \Delta t}$

**$|20,0\rangle$ のインデックス**: $2 \times 3 \times 2 + 0 \times 2 + 0 = 12$
**$|11,1\rangle$ のインデックス**: $1 \times 3 \times 2 + 1 \times 2 + 1 = 9$

（注: インデックスは基底の具体的な列挙順序に依存するため、実装時に確認が必要）

### 8.4 MQT-Qudits回路での実装設計

#### 8.4.1 ゲート分解

各Stinespringユニタリは、MQT-Quditsの基本ゲートセットに分解する：

**単一分子Lindblad演算子**（蛍光、燐光、IC、ISC）の場合：
- 6×6 ユニタリ行列（系qutrit + ancilla qubit）
- `CustomTwo` ゲートとして実装
- `IntegratedSparseCompilerV2` により基本ゲートに分解
- 典型的には ~6 基本ゲートに分解

**TTA Lindblad演算子**の場合：
- 18×18 ユニタリ行列（2 qutrit + 1 ancilla qubit）
- より多くの基本ゲートが必要

#### 8.4.2 基本ゲートセット

| ゲート名 | 記号 | 数学的定義 |
|---------|------|----------|
| 仮想Z回転 | `VirtRz` | $\text{diag}(e^{i\phi_0}, e^{i\phi_1}, \ldots)$ |
| 部分空間回転 | `R` | $\hat{R}_{mn}(\theta, \phi) = e^{-i\theta(\cos\phi\,\hat{\sigma}_{mn}^x + \sin\phi\,\hat{\sigma}_{mn}^y)/2}$ |
| 制御交換 | `CEx` | 制御付き部分空間交換（2-quditゲート） |

Givens回転：

$$
G_{mn}(\theta, \phi) = I + (\cos\theta - 1)(|m\rangle\langle m| + |n\rangle\langle n|) + \sin\theta(e^{i\phi}|m\rangle\langle n| - e^{-i\phi}|n\rangle\langle m|)
$$

### 8.5 1 Trotterステップの回路フロー

```
1 Trotterステップ（Qudit GKSL）:

┌──────────────────────────────────────────────┐
│ ユニタリ前半 (Δt/2)                           │
│  ├─ H0: VirtRz ゲート（各 qutrit）            │
│  │   diag(1, e^{-iE_T·Δt/(2ℏ)},              │
│  │        e^{-iE_S·Δt/(2ℏ)})                  │
│  └─ H_transfer: CEx 系ゲート（各ペア）         │
│     exp(-iV·Δt·(|01⟩⟨10|+|10⟩⟨01|)/(2ℏ))     │
├──────────────────────────────────────────────┤
│ Lindblad散逸ステップ (Δt)                      │
│  ├─ TTA Stinespring × 6                       │
│  │   各: CustomTwo(qutrit_i, qutrit_j, anc)   │
│  ├─ 蛍光 Stinespring × 4                      │
│  │   各: CustomTwo(qutrit_i, anc)              │
│  ├─ 燐光 Stinespring × 4                      │
│  ├─ IC Stinespring × 4                         │
│  ├─ ISC S→T Stinespring × 4                   │
│  └─ ISC T→S Stinespring × 4                   │
│  各Stinespring後: ancilla → |0⟩ リセット       │
├──────────────────────────────────────────────┤
│ ユニタリ後半 (Δt/2) [逆順]                     │
│  ├─ H_transfer [逆順ペア]                      │
│  └─ H0 [逆順分子]                              │
└──────────────────────────────────────────────┘
```

### 8.6 量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 系 qutrit | 4 | 4分子 × 1 qutrit |
| Ancilla qubit | 26 | 26個のLindblad演算子 |
| **合計** | **4 qutrit + 26 qubit** | |
| 等価 qubit 数 | $4 \times 2 + 26 = 34$ | 1 qutrit ≈ 2 qubit |

### 8.7 QuditGKSLSimulator フローチャート

```
┌────────────────────────────────────┐
│     simulate_shot_based() 開始      │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 1. MQT-Qudits回路初期化            │
│   系: 4 qutrit (d=3)              │
│   ancilla: 26 qubit (d=2)         │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 2. 初期状態準備                     │
│   qutrit 0: |1⟩                   │
│   qutrit 1: |0⟩                   │
│   qutrit 2: |0⟩                   │
│   qutrit 3: |1⟩                   │
│   ancilla: all |0⟩                │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 3. Stinespringユニタリの構築        │
│   FOR 各Lindblad演算子:             │
│     U = build_stinespring_unitary  │
│         (L, gamma, dt)             │
│     検証: U†U = I                  │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 4. N_stepsループ                   │
│   FOR step IN 1..N_steps:          │
│     ├─ build_unitary_step(dt/2)    │
│     ├─ build_lindblad_step(dt)     │
│     ├─ ancilla リセット            │
│     └─ build_unitary_step(dt/2)    │
│         [逆順]                     │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 5. 測定とポスト処理                 │
│   shots回の測定を実行               │
│   統計から個体数を推定              │
└──────────┬─────────────────────────┘
           ▼
┌────────────────────────────────────┐
│ 6. 結果出力                         │
│   + MQT-Qudits回路情報             │
└────────────────────────────────────┘
```

### 8.8 Qudit実装の利点（GKSL文脈）

1. **禁止状態なし**: Lindblad演算子が禁止状態に遷移するリスクがゼロ
2. **自然な部分空間回転**: $|0\rangle\langle 2|$ 等の遷移演算子が部分空間回転 $R_{02}(\theta)$ で直接実装可能
3. **少ないゲート数**: Stinespringユニタリの疎構造を認識し、効率的に分解可能
4. **Qubit比較**: $N=4$ で Qubit-NB ~34 qubit vs Qudit-NB ~34 等価qubitだが、ゲート数でQudit版が優位

---

## 9. シナリオ6: Qudit GKSL・ボソン有り 詳細設計

### 9.1 設計概要

シナリオ5にフォノンモードの qudit 表現を追加。フォノンモードを高次元 qudit で自然に表現する。

### 9.2 フォノンの Qudit エンコーディング

$$
|n\rangle_{\text{Fock}} \leftrightarrow |n\rangle_{d_{\text{ph}}}
$$

$n_{\max} = 2$ の場合: $d_{\text{ph}} = 3$（qutrit）

Qubitのバイナリエンコーディングと異なり禁止状態が存在しない。

### 9.3 量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 電子系 qutrit | 4 | $d = 3$ |
| フォノン qutrit | 4 | $d = 3$ ($n_{\max} = 2$) |
| Ancilla qubit | 26 | Lindblad演算子数 |
| **合計** | **8 qutrit + 26 qubit** | |
| 等価 qubit 数 | $8 \times 2 + 26 = 42$ | |

### 9.4 電子-フォノン結合の Qudit 回路設計

電子qutrit（$d=3$）とフォノンqutrit（$d=3$）間の2-quditゲート：

$$
\hat{U}_{e\text{-ph}} = |0\rangle\langle 0|_{\text{el}} \otimes \hat{I}_{\text{ph}} + |1\rangle\langle 1|_{\text{el}} \otimes e^{-ig\Delta t(\hat{a}+\hat{a}^\dagger)/\hbar} + |2\rangle\langle 2|_{\text{el}} \otimes \hat{I}_{\text{ph}}
$$

$9 \times 9$ のユニタリ行列を `CustomTwo` ゲートとして MQT-Qudits に実装。

フォノンqutrit上の $e^{-ig\Delta t(\hat{a}+\hat{a}^\dagger)/\hbar}$ の行列表現：

$$
\hat{a} + \hat{a}^\dagger = \begin{pmatrix} 0&1&0\\1&0&\sqrt{2}\\0&\sqrt{2}&0 \end{pmatrix}
$$

$$
e^{-i\alpha(\hat{a}+\hat{a}^\dagger)} \text{ を行列指数関数で計算} \quad (\alpha = g\Delta t / \hbar)
$$

### 9.5 QuditGKSLBosonSimulator フローチャート

```
┌────────────────────────────────────────┐
│     simulate_shot_based() 開始          │
└──────────┬─────────────────────────────┘
           ▼
┌────────────────────────────────────────┐
│ 1. qudit/qubit割当て                    │
│   電子 qutrit: q0-q3 (d=3)            │
│   フォノン qutrit: q4-q7 (d=3)        │
│   ancilla qubit: q8-q33 (d=2)         │
└──────────┬─────────────────────────────┘
           ▼
┌────────────────────────────────────────┐
│ 2. 初期状態準備                         │
│   電子: |1,0,0,1⟩                      │
│   フォノン: |0,0,0,0⟩（真空）           │
│   ancilla: |0...0⟩                     │
└──────────┬─────────────────────────────┘
           ▼
┌────────────────────────────────────────┐
│ 3. 1 Trotterステップの構築              │
│  ├─ ユニタリ前半:                       │
│  │   H0(el) + H_transfer               │
│  │   + H_phonon + H_eph                │
│  ├─ Lindblad散逸ステップ                │
│  │   (電子qutrit + ancillaのみ)         │
│  └─ ユニタリ後半                        │
└──────────┬─────────────────────────────┘
           ▼
┌────────────────────────────────────────┐
│ 4. N_stepsループ→測定→結果出力          │
└────────────────────────────────────────┘
```

---

## 10. 検証設計

### 10.1 数学的整合性の検証体系

全シナリオで以下の検証を実施する。

#### 10.1.1 トレース保存

$$
|1 - \text{Tr}[\hat{\rho}(t)]| < 10^{-8} \quad \forall t
$$

実装: 各時間ステップで `abs(1 - trace(rho))` を計算。違反時は `PhysicsViolationError` を raise。

#### 10.1.2 エルミート性

$$
\|\hat{\rho} - \hat{\rho}^\dagger\|_F < 10^{-10}
$$

実装: `norm(rho - rho.conj().T, 'fro')` を計算。

#### 10.1.3 正定値性

$$
\lambda_k \geq -10^{-10} \quad \forall k
$$

実装: `eigvalsh(rho)` の最小値を確認。

#### 10.1.4 エントロピー非減少

$$
S(t + \Delta t) \geq S(t) - 10^{-8}
$$

散逸のみのGKSL方程式ではエントロピーは単調増加する。微小な数値的減少（$\sim 10^{-8}$）は許容。

#### 10.1.5 粒子数保存

$$
N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = N = 4 \pm 10^{-8}
$$

GKSL方程式のトレース保存より自動的に成立するが、数値的に検証。

### 10.2 物理法則の検証

#### 10.2.1 エネルギー変化の方向

散逸により系のエネルギーは一般に減少する（基底状態 $|0000\rangle$ への緩和）。

$$
\langle \hat{H}_0 \rangle(t) = \text{Tr}[\hat{H}_0 \hat{\rho}(t)]
$$

長時間極限: $\langle \hat{H}_0 \rangle \to 0$

#### 10.2.2 定常状態の確認

十分長い時間（$t \gg 1/\gamma_{\min}$）で：

$$
\hat{\rho}(\infty) \to |0000\rangle\langle 0000|
$$

（全分子が基底状態に緩和）

### 10.3 テストケース設計

| テストID | テスト名 | 条件 | 期待結果 | 許容誤差 |
|---------|---------|------|---------|---------|
| T1 | ユニタリ極限 | $\gamma_\alpha = 0$ 全て | 現行ユニタリ版と個体数一致 | $10^{-6}$ |
| T2 | 蛍光のみ | $V=0$, $\Gamma_{\text{fl}}>0$, 初期$S_1$ | $N_{S_1} \sim e^{-\Gamma_{\text{fl}} t / \hbar}$ | $10^{-3}$ |
| T3 | 純TTA | $V=0$, $\gamma_{\text{TTA}}>0$, 初期$T_1T_1$ | $N_{T_1}$ 単調減少 | 定性的 |
| T4 | トレース保存 | 全条件 | $\text{Tr}[\hat{\rho}] = 1$ | $10^{-8}$ |
| T5 | 正定値性 | 全条件 | $\lambda_{\min} \geq -10^{-10}$ | $10^{-10}$ |
| T6 | 古典-Qubit一致 | 同パラメータ | 個体数一致（Statevector） | $10^{-3}$ |
| T7 | 古典-Qudit一致 | 同パラメータ | 個体数一致（Statevector） | $10^{-3}$ |
| T8 | 長時間緩和 | $t \to \infty$ | $\hat{\rho} \to |0000\rangle\langle 0000|$ | $10^{-2}$ |

#### 10.3.1 テストT1の詳細設計（ユニタリ極限）

```
test_unitary_limit():
  params = GKSLPhysicalParameters()
  # 全散逸速度を0に設定
  params.gamma_TTA = 0
  params.Gamma_fl = 0
  params.Gamma_ph = 0
  params.k_IC = 0
  params.k_ISC_ST = 0
  params.k_ISC_TS = 0
  # 注: ハミルトニアンにTTA項が無いため、
  # 現行のH_transferのみの時間発展と比較する

  gksl_result = ClassicalGKSLSimulator(params).simulate(...)
  unitary_result = ClassicalSuzukiTrotterSimulator(params_unitary).simulate(...)
  # 注: 現行版はH_TTAを含むため、GKSLのユニタリ極限とは直接比較できない
  # → H_transfer のみの現行版と比較するか、
  #   新たにH_TTA無しの現行版を構築して比較

  FOR t IN times:
    assert abs(gksl_N_T1[t] - unitary_N_T1[t]) < 1e-6
```

#### 10.3.2 テストT2の詳細設計（蛍光のみの解析解との比較）

$V = 0$（エネルギー移動なし）、蛍光のみの場合、分子 $i$ が $|S_1\rangle$ にある確率の解析解：

$$
P_{S_1}^{(i)}(t) = P_{S_1}^{(i)}(0) \cdot e^{-\Gamma_{\text{fl}} t / \hbar}
$$

```
test_fluorescence_only():
  params = GKSLPhysicalParameters()
  params.V = 0  # エネルギー移動なし
  params.gamma_TTA = 0
  params.Gamma_ph = 0
  params.k_IC = 0
  params.k_ISC_ST = 0
  params.k_ISC_TS = 0
  # 初期状態: 全分子が S1
  params.initial_state_type = 'all_S1'
  # → |2,2,2,2⟩

  result = ClassicalGKSLSimulator(params).simulate(...)

  FOR t IN times:
    expected_N_S1 = 4 * exp(-params.Gamma_fl * t / params.hbar)
    assert abs(result.N_S1[t] - expected_N_S1) < 1e-3
```

### 10.4 Trotter誤差の評価設計

#### 10.4.1 理論的誤差限界

2次Trotter誤差（本仕様で使用）：

$$
\|e^{(\hat{A}+\hat{B})t} - (e^{\hat{A}t/2N}e^{\hat{B}t/N}e^{\hat{A}t/2N})^N\| \leq \frac{t^3}{12N^2}\|[[\hat{A},\hat{B}],\hat{A}+\hat{B}]\|
$$

#### 10.4.2 本仕様のパラメータでの見積もり

$\hat{A} = -i\hat{H}_0/\hbar$, $\hat{B} = -i\hat{H}_{\text{transfer}}/\hbar$

$$
\|[\hat{H}_0, \hat{H}_{\text{transfer}}]\| \sim V \cdot (E_T - E_{S_0}) = 0.1 \times 1.5 = 0.15 \text{ eV}^2
$$

交換子ノルムの上界に基づく粗い上界。実際の個体数の系統誤差は $N_{\text{steps}} = 100$ で $\sim 10^{-3}$ 程度。

### 10.5 Stinespring忠実度検証設計

理想的Lindblad時間発展（古典GKSLソルバー）とStinespring量子回路実装の比較：

$$
F = \left(\text{Tr}\sqrt{\sqrt{\hat{\rho}_{\text{ideal}}}\hat{\rho}_{\text{Stinespring}}\sqrt{\hat{\rho}_{\text{ideal}}}}\right)^2
$$

目標: $F > 0.99$（$N_{\text{steps}} = 100$）

```
verify_stinespring_fidelity(rho_classical, rho_quantum):
  # scipy.linalg.sqrtm を使用
  sqrt_rho_c = sqrtm(rho_classical)
  inner = sqrtm(sqrt_rho_c @ rho_quantum @ sqrt_rho_c)
  F = real(trace(inner))^2
  IF F < 0.99:
    log_warning(f"Fidelity {F} < 0.99 target")
  RETURN F
```

### 10.6 検証フローチャート

```
┌─────────────────────────────────────┐
│          検証パイプライン             │
└──────────┬──────────────────────────┘
           ▼
┌─────────────────────────────────────┐
│ Phase 1: 古典GKSL単体テスト          │
│  ├─ T1: ユニタリ極限                 │
│  ├─ T2: 蛍光解析解                   │
│  ├─ T3: TTA定性的挙動               │
│  ├─ T4: トレース保存                 │
│  └─ T5: 正定値性                     │
└──────────┬──────────────────────────┘
           ▼
┌─────────────────────────────────────┐
│ Phase 2: 量子-古典一致検証            │
│  ├─ T6: 古典 vs Qubit GKSL          │
│  └─ T7: 古典 vs Qudit GKSL          │
└──────────┬──────────────────────────┘
           ▼
┌─────────────────────────────────────┐
│ Phase 3: 物理的整合性                 │
│  ├─ T8: 長時間緩和                   │
│  ├─ エントロピー非減少               │
│  └─ エネルギー変化の方向              │
└──────────┬──────────────────────────┘
           ▼
┌─────────────────────────────────────┐
│ Phase 4: Stinespring忠実度           │
│  ├─ 各Lindblad演算子の忠実度         │
│  └─ 全体の累積忠実度                 │
└─────────────────────────────────────┘
```


---

## 11. 可視化・比較フレームワーク設計

### 11.1 比較対象の全体像

現行ノートブック（3手法）に加え、GKSL版（6シナリオ）を含む計9手法を比較する。

| # | 手法 | 方式 | ボソン | 状態記述 | TTA実装 |
|---|------|------|--------|---------|---------|
| 1 | 古典 Suzuki-Trotter | ユニタリ | 無し | 状態ベクトル | ハミルトニアン |
| 2 | 古典 GKSL (NB) | GKSL | 無し | 密度行列 | Lindblad |
| 3 | 古典 GKSL (B) | GKSL | 有り | 密度行列 | Lindblad |
| 4 | Qubit ユニタリ | ユニタリ | 無し | 状態ベクトル | ハミルトニアン |
| 5 | Qubit GKSL (NB) | GKSL | 無し | 密度行列 | Stinespring |
| 6 | Qubit GKSL (B) | GKSL | 有り | 密度行列 | Stinespring |
| 7 | Qudit ユニタリ | ユニタリ | 無し | 状態ベクトル | ハミルトニアン |
| 8 | Qudit GKSL (NB) | GKSL | 無し | 密度行列 | Stinespring |
| 9 | Qudit GKSL (B) | GKSL | 有り | 密度行列 | Stinespring |

### 11.2 可視化関数の設計

#### 11.2.1 plot_population_dynamics（現行互換）

```
plot_population_dynamics(results, title):
  入力: results (GKSLResult辞書), title (str)
  出力: matplotlib Figure

  fig, ax = subplots(1, 1, figsize=(10, 6))
  times = results['times']
  N_S0 = [p['N_S0'] for p in results['populations']]
  N_T1 = [p['N_T1'] for p in results['populations']]
  N_S1 = [p['N_S1'] for p in results['populations']]

  ax.plot(times, N_S0, 'b-o', label='N_S0', markersize=3)
  ax.plot(times, N_T1, 'r-s', label='N_T1', markersize=3)
  ax.plot(times, N_S1, 'g-^', label='N_S1', markersize=3)
  ax.set_xlabel('Time (fs)')
  ax.set_ylabel('Population')
  ax.set_ylim(0, 4.5)
  ax.legend()
  ax.set_title(title)
```

#### 11.2.2 plot_per_molecule_populations（現行互換）

```
plot_per_molecule_populations(results, title):
  fig, axes = subplots(2, 2, figsize=(14, 10))
  FOR mol_idx IN 0..3:
    ax = axes[mol_idx // 2][mol_idx % 2]
    S0 = [p['S0_per_mol'][mol_idx] for p in results['per_molecule_populations']]
    T1 = [p['T1_per_mol'][mol_idx] for p in results['per_molecule_populations']]
    S1 = [p['S1_per_mol'][mol_idx] for p in results['per_molecule_populations']]
    ax.plot(times, S0, 'b-', label='S0')
    ax.plot(times, T1, 'r-', label='T1')
    ax.plot(times, S1, 'g-', label='S1')
    ax.set_ylim(0, 1.1)
    ax.set_title(f'Molecule {mol_idx}')
```

#### 11.2.3 plot_gksl_comparison（新規）

```
plot_gksl_comparison(unitary_results, gksl_results, title):
  fig, axes = subplots(1, 2, figsize=(16, 6))

  # 左パネル: 個体数の比較
  ax = axes[0]
  ax.plot(times, N_S0_unitary, 'b-', label='N_S0 (unitary)')
  ax.plot(times, N_T1_unitary, 'r-', label='N_T1 (unitary)')
  ax.plot(times, N_S1_unitary, 'g-', label='N_S1 (unitary)')
  ax.plot(times, N_S0_gksl, 'b--', label='N_S0 (GKSL)')
  ax.plot(times, N_T1_gksl, 'r--', label='N_T1 (GKSL)')
  ax.plot(times, N_S1_gksl, 'g--', label='N_S1 (GKSL)')
  ax.set_title('Population: Unitary vs GKSL')

  # 右パネル: エントロピーと純度
  ax = axes[1]
  ax.plot(times, gksl_results['entropy'], 'k-', label='Entropy')
  ax2 = ax.twinx()
  ax2.plot(times, gksl_results['purity'], 'm-', label='Purity')
  ax.set_title('Entropy and Purity')
```

#### 11.2.4 plot_entropy_dynamics（新規）

```
plot_entropy_dynamics(results_list, labels, title):
  fig, axes = subplots(1, 2, figsize=(14, 5))

  # 左: von Neumannエントロピー
  FOR result, label IN zip(results_list, labels):
    axes[0].plot(result['times'], result['entropy'], label=label)
  axes[0].set_ylabel('von Neumann Entropy')
  axes[0].set_xlabel('Time (fs)')

  # 右: 純度
  FOR result, label IN zip(results_list, labels):
    axes[1].plot(result['times'], result['purity'], label=label)
  axes[1].set_ylabel('Purity Tr[ρ²]')
```

#### 11.2.5 plot_6scenario_comparison（新規）

```
plot_6scenario_comparison(results_dict, title):
  fig, axes = subplots(2, 3, figsize=(20, 12))
  scenarios = [
    ('Classical NB', 'Classical B'),
    ('Qubit NB', 'Qubit B'),
    ('Qudit NB', 'Qudit B'),
  ]
  FOR col, (nb_key, b_key) IN enumerate(scenarios):
    # 行1: ボソン無し
    plot_pops(axes[0][col], results_dict[nb_key])
    # 行2: ボソン有り
    plot_pops(axes[1][col], results_dict[b_key])
```

### 11.3 包括的比較表の設計

```python
comparison_data = {
    'シナリオ': ['古典ユニタリ', '古典GKSL(NB)', '古典GKSL(B)',
                'Qubitユニタリ', 'QubitGKSL(NB)', 'QubitGKSL(B)',
                'Quditユニタリ', 'QuditGKSL(NB)', 'QuditGKSL(B)'],
    '状態記述': ['状態ベクトル', '密度行列', '密度行列',
                '状態ベクトル', '密度行列', '密度行列',
                '状態ベクトル', '密度行列', '密度行列'],
    'TTA実装': ['ハミルトニアン', 'Lindblad', 'Lindblad',
               'ハミルトニアン', 'Stinespring', 'Stinespring',
               'ハミルトニアン', 'Stinespring', 'Stinespring'],
    '散逸過程数': [0, 26, 26, 0, 26, 26, 0, 26, 26],
    'ヒルベルト空間次元': [81, 81, 6561, 256, 256, '16384+',
                          81, 81, 6561],
    '量子リソース': ['N/A', 'N/A', 'N/A',
                    '8 qubit', '34 qubit', '42 qubit',
                    '4 qutrit', '4 qutrit+26 qubit', '8 qutrit+26 qubit'],
    'N_T1 (最終)': [],   # 実行時に記録
    'N_S1 (最終)': [],   # 実行時に記録
    'エントロピー (最終)': [],
    '純度 (最終)': [],
    '実行時間 (秒)': [],
}
```

### 11.4 量子回路可視化設計

#### 11.4.1 Qubit GKSL回路

```python
from qiskit.visualization import circuit_drawer
circuit_drawer(step_circuit, output='mpl', fold=80)
```

対象:
- 1トロッターステップ全体の回路
- UnitaryGate版 vs 基本ゲート分解版の並列比較

#### 11.4.2 Qudit GKSL回路

```python
from mqt.qudits.visualisation import visualize_circuit
visualize_circuit(step_circuit)
```

対象:
- CustomTwo版 vs 基本ゲート分解版の並列比較
- 系qutrit と ancilla qubit の区別

### 11.5 スケーラビリティ比較設計

| シナリオ | ヒルベルト空間次元 | 密度行列要素数 | $N$依存性 | 実用限界 |
|---------|-----------------|-------------|-----------|---------|
| 古典GKSL (NB) | $3^N$ | $9^N$ | $O(9^N)$ | $N \sim 6$ |
| 古典GKSL (B) | $3^N \cdot d_{\text{ph}}^N$ | $(3 d_{\text{ph}})^{2N}$ | $O((3d_{\text{ph}})^{2N})$ | $N \sim 4$ |
| Qubit GKSL (NB) | $2^{2N}$ | — | $O(N^2)$ qubit | $N \sim 100+$ |
| Qubit GKSL (B) | $2^{2N+N\lceil\log_2 d_{\text{ph}}\rceil}$ | — | $O(N^2 \log n_{\max})$ | $N \sim 50+$ |
| Qudit GKSL (NB) | $3^N$ | — | $O(N)$ qutrit | $N \sim 100+$ |
| Qudit GKSL (B) | $3^N \cdot d_{\text{ph}}^N$ | — | $O(N)$ qudit | $N \sim 50+$ |

---

## 12. 全体統合フローチャート

### 12.1 メインノートブックの実行フロー

```
┌──────────────────────────────────────────────────────────────┐
│   quantum_dynamics_gksl_comparison.ipynb 全体フロー           │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル1: パラメータ設定                                          │
│  ├─ GKSLPhysicalParameters()                                 │
│  ├─ GKSLPhysicalParametersWithBoson()                        │
│  └─ validate() 実行                                          │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル2: 現行ユニタリシミュレーション（比較基準）                   │
│  ├─ ClassicalSuzukiTrotterSimulator.simulate()                │
│  ├─ QubitMolecularDynamicsSimulator.simulate()                │
│  └─ SuzukiTrotterMQTQuditSimulator.simulate_shot_based()      │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル3: 古典GKSL（ボソン無し）                                  │
│  ├─ ClassicalGKSLSimulator(params).simulate()                 │
│  ├─ 検証: トレース保存、正定値性、粒子数保存                     │
│  └─ plot_gksl_comparison(unitary, gksl)                       │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル4: 古典GKSL（ボソン有り）                                  │
│  ├─ ClassicalGKSLBosonSimulator(params_boson).simulate()       │
│  └─ ボソン無しとの比較                                         │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル5: Qubit GKSL（ボソン無し）                                │
│  ├─ QubitGKSLSimulator(params).simulate()                     │
│  ├─ 古典GKSLとの一致検証                                      │
│  └─ 回路情報の表示                                             │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル6: Qubit GKSL（ボソン有り）                                │
│  ├─ QubitGKSLBosonSimulator(params_boson).simulate()           │
│  └─ 古典GKSL(B)との一致検証                                   │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル7: Qudit GKSL（ボソン無し）                                │
│  ├─ QuditGKSLSimulator(params).simulate_shot_based()           │
│  ├─ 古典GKSLとの一致検証                                       │
│  └─ 回路情報の表示                                              │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル8: Qudit GKSL（ボソン有り）                                │
│  ├─ QuditGKSLBosonSimulator(params_boson).simulate_shot_based()│
│  └─ 古典GKSL(B)との一致検証                                    │
└──────────┬───────────────────────────────────────────────────┘
           ▼
┌──────────────────────────────────────────────────────────────┐
│ セル9: 包括的比較                                               │
│  ├─ plot_6scenario_comparison()                                │
│  ├─ plot_entropy_dynamics()                                    │
│  ├─ 比較表の表示                                                │
│  └─ スケーラビリティ分析                                        │
└──────────────────────────────────────────────────────────────┘
```

### 12.2 シナリオ1（古典GKSL・NB）の詳細実行フロー

```
┌──────────────────────┐
│ ClassicalGKSLSimulator│
│     .simulate()       │
└──────────┬───────────┘
           ▼
┌──────────────────────────────────┐
│ ① H_sys構築 (81×81)              │
│  H0 = Σ I⊗...⊗h_single⊗...⊗I   │
│  H_tr = Σ V(|01⟩⟨10|+h.c.)     │
│  H_sys = H0 + H_tr              │
│  検証: H_sys = H_sys†            │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ ② Lindblad演算子構築 (26個)       │
│  FOR each process:               │
│    L_full = I⊗...⊗L_local⊗...⊗I │
│  operators = [(γ, L)]×26         │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ ③ 超演算子構築                    │
│    (6561×6561)                   │
│  L_H = -i/ℏ (H⊗I - I⊗H^T)      │
│  FOR (γ,L) IN operators:         │
│    LdL = L†L                     │
│    L_α = γ(L⊗L̄ - ½LdL⊗I        │
│           - ½I⊗(LdL)^T)          │
│    L_super += L_α                 │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ ④ 初期密度行列構築                │
│  |ψ₀⟩ = |1,0,0,1⟩               │
│  ρ₀ = |ψ₀⟩⟨ψ₀|                  │
│  ρ₀_vec = vec(ρ₀) (列優先)       │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ ⑤ ODE積分                        │
│  dρ_vec/dt = L_super · ρ_vec     │
│  solve_ivp(method='RK45' or 'BDF│
│    ', rtol=1e-10, atol=1e-12)    │
│                                  │
│  FOR each t_eval point:          │
│    ρ(t) = reshape(sol.y[:, k])   │
│    ├─ validate: Tr=1, λ≥0, H†   │
│    ├─ populations                │
│    ├─ entropy                    │
│    ├─ purity                     │
│    └─ coherences                 │
└──────────┬───────────────────────┘
           ▼
┌──────────────────────────────────┐
│ ⑥ 結果辞書の構築                  │
│  RETURN {times, populations,     │
│    per_molecule_populations,      │
│    entropy, purity, trace,       │
│    rho_final, coherences,        │
│    elapsed_time, method,         │
│    ode_solver, n_function_evals} │
└──────────────────────────────────┘
```

### 12.3 シナリオ5（Qudit GKSL・NB）の詳細回路構築フロー

```
┌───────────────────────────────────┐
│ QuditGKSLSimulator                │
│   .simulate_shot_based()          │
└──────────┬────────────────────────┘
           ▼
┌───────────────────────────────────┐
│ ① 回路初期化                       │
│  circuit = QuantumCircuit(         │
│    dimensions = [3,3,3,3,  ← 系    │
│                  2,2,...,2] ← anc  │
│    num_qudits = 30)                │
└──────────┬────────────────────────┘
           ▼
┌───────────────────────────────────┐
│ ② 初期状態準備                     │
│  X gate: qutrit 0 → |1⟩           │
│  X gate: qutrit 3 → |1⟩           │
│  (他は|0⟩のまま)                   │
└──────────┬────────────────────────┘
           ▼
┌───────────────────────────────────┐
│ ③ Stinespringユニタリ事前構築      │
│  FOR 各 (γ, L) IN 26 operators:  │
│    θ = √(γ·dt)                    │
│    G = L⊗σ⁻ + L†⊗σ⁺              │
│    U = expm(-iθG)                  │
│    検証: U†U = I                   │
│  → U_stinespring_list (26個)       │
└──────────┬────────────────────────┘
           ▼
┌───────────────────────────────────────────────────┐
│ ④ N_stepsループ                                    │
│  FOR step IN 1..N_steps:                           │
│                                                    │
│  ┌─────────────────────────────────────────────┐   │
│  │ ④-a ユニタリ前半 (dt/2)                      │   │
│  │  FOR mol IN 0..3:                            │   │
│  │    VirtRz(mol, phases=[0, -E_T·dt/(2ℏ),     │   │
│  │                        -E_S·dt/(2ℏ)])        │   │
│  │  FOR (i,j) IN [(0,1),(1,2),(2,3)]:           │   │
│  │    U_tr = expm(-i·V·dt/(2ℏ)·H_tr_pair)      │   │
│  │    CustomTwo(qutrit_i, qutrit_j, U_tr)       │   │
│  └─────────────────────────────────────────────┘   │
│                                                    │
│  ┌─────────────────────────────────────────────┐   │
│  │ ④-b Lindblad散逸ステップ (dt)                │   │
│  │  ancilla_idx = 4  (ancillaの開始インデックス) │   │
│  │  FOR k, U_k IN enumerate(U_stinespring_list):│   │
│  │    IF TTA演算子:                              │   │
│  │      CustomMulti(qutrit_i, qutrit_j,         │   │
│  │                  ancilla[ancilla_idx], U_k)   │   │
│  │    ELSE:                                      │   │
│  │      CustomTwo(qutrit_i,                      │   │
│  │               ancilla[ancilla_idx], U_k)      │   │
│  │    ancilla_idx += 1                           │   │
│  │  Ancilla全体をリセット → |0⟩                  │   │
│  └─────────────────────────────────────────────┘   │
│                                                    │
│  ┌─────────────────────────────────────────────┐   │
│  │ ④-c ユニタリ後半 (dt/2) [逆順]               │   │
│  │  FOR (i,j) IN [(2,3),(1,2),(0,1)]:           │   │
│  │    CustomTwo(qutrit_i, qutrit_j, U_tr†)      │   │
│  │  FOR mol IN 3..0:                            │   │
│  │    VirtRz(mol, 逆位相)                        │   │
│  └─────────────────────────────────────────────┘   │
│                                                    │
│  ┌─────────────────────────────────────────────┐   │
│  │ ④-d 測定（track_dynamics=True時）             │   │
│  │  circuit.measure(系qutrit)                   │   │
│  │  shots回の測定結果から個体数を推定             │   │
│  └─────────────────────────────────────────────┘   │
└──────────┬────────────────────────────────────────┘
           ▼
┌───────────────────────────────────┐
│ ⑤ 結果集約                        │
│  各ステップの個体数を辞書に格納    │
│  回路情報（ゲート数、深さ）を記録   │
│  RETURN GKSLResult                 │
└───────────────────────────────────┘
```

### 12.4 全6シナリオの相互検証フロー

```
┌────────────────────────────────────────────────────────────┐
│                  相互検証パイプライン                         │
└──────────┬─────────────────────────────────────────────────┘
           ▼
┌─────────────────┐
│ 古典GKSL (NB)    │ ←── 基準実装（最も信頼性が高い）
│ ClassicalGKSL    │
└──────┬──────────┘
       │
       ├──── 一致検証 ────► Qubit GKSL (NB)
       │     |ΔN_X| < ε       QubitGKSL
       │
       ├──── 一致検証 ────► Qudit GKSL (NB)
       │     |ΔN_X| < ε       QuditGKSL
       │
       └──── ユニタリ極限 ──► 古典ユニタリ
             γ_α = 0             ClassicalSuzukiTrotter

┌─────────────────┐
│ 古典GKSL (B)     │ ←── ボソン有りの基準
│ ClassicalGKSLBos │
└──────┬──────────┘
       │
       ├──── 一致検証 ────► Qubit GKSL (B)
       │                      QubitGKSLBoson
       │
       ├──── 一致検証 ────► Qudit GKSL (B)
       │                      QuditGKSLBoson
       │
       └──── g_eph=0極限 ──► 古典GKSL (NB)
             ボソン無しと一致     ClassicalGKSL
```

---

## 付録A: 記号一覧

### A.1 物理量

| 記号 | 意味 | 単位 | デフォルト値 |
|------|------|------|------------|
| $N$ | 分子数 | — | 4 |
| $d$ | 各分子の電子状態数 | — | 3 |
| $E_T$ | 三重項エネルギー | eV | 1.5 |
| $E_S$ | 一重項エネルギー | eV | 3.0 |
| $V$ | エネルギー移動積分 | eV | 0.1 |
| $\gamma_{\text{TTA}}$ | TTA散逸速度定数 | eV/ℏ | 0.05 |
| $\Gamma_{\text{fl}}$ | 蛍光放出速度 | eV/ℏ | 0.01 |
| $\Gamma_{\text{ph}}$ | 燐光放出速度 | eV/ℏ | $10^{-6}$ |
| $k_{\text{IC}}$ | 内部転換速度 | eV/ℏ | 0.005 |
| $k_{\text{ISC}}^{S\to T}$ | ISC S₁→T₁ 速度 | eV/ℏ | 0.003 |
| $k_{\text{ISC}}^{T\to S}$ | ISC T₁→S₀ 速度 | eV/ℏ | $10^{-5}$ |
| $\hbar$ | 換算プランク定数 | eV·fs | 0.6582119569 |
| $T_{\text{total}}$ | 総シミュレーション時間 | fs | 100 |
| $N_{\text{steps}}$ | 時間ステップ数 | — | 100 |
| $\Delta t$ | 時間刻み | fs | 1.0 |
| $\omega_{\text{ph}}$ | フォノン振動周波数 | eV/ℏ | 0.15 |
| $g$ | 電子-フォノン結合定数 | eV | 0.02 |
| $n_{\max}$ | フォノンFock空間の切断 | — | 2 |

### A.2 数学的記号

| 記号 | 意味 |
|------|------|
| $\hat{\rho}$ | 密度演算子（密度行列） |
| $\hat{H}$ | ハミルトニアン演算子 |
| $\hat{L}_\alpha$ | Lindblad演算子（ジャンプ演算子） |
| $\mathcal{L}$ | リンドブラディアン（超演算子） |
| $\mathcal{D}[\hat{L}]$ | Lindblad超演算子 $\hat{L}\hat{\rho}\hat{L}^\dagger - \frac{1}{2}\{\hat{L}^\dagger\hat{L}, \hat{\rho}\}$ |
| $\text{Tr}[\cdot]$ | トレース |
| $\text{Tr}_E[\cdot]$ | 環境系に対する部分トレース |
| $[\hat{A}, \hat{B}]$ | 交換子 $\hat{A}\hat{B} - \hat{B}\hat{A}$ |
| $\{\hat{A}, \hat{B}\}$ | 反交換子 $\hat{A}\hat{B} + \hat{B}\hat{A}$ |
| $\hat{A}^\dagger$ | エルミート共役 |
| $\bar{\hat{A}}$ | 要素ごとの複素共役（転置なし） |
| $\hat{A}^T$ | 転置 |
| $\otimes$ | テンソル積 / Kronecker積 |
| $\|\cdot\|_F$ | フロベニウスノルム |
| $\|\cdot\|_1$ | トレースノルム |
| $S(\hat{\rho})$ | von Neumannエントロピー $-\text{Tr}[\hat{\rho}\ln\hat{\rho}]$ |
| $P(\hat{\rho})$ | 純度 $\text{Tr}[\hat{\rho}^2]$ |
| $F(\hat{\rho}, \hat{\sigma})$ | 忠実度（フィデリティ） |

### A.3 状態表記

| 表記 | 意味 |
|------|------|
| $\|S_0\rangle = \|0\rangle$ | 一重項基底状態 |
| $\|T_1\rangle = \|1\rangle$ | 三重項第一励起状態 |
| $\|S_1\rangle = \|2\rangle$ | 一重項第一励起状態 |
| $\|1,0,0,1\rangle$ | 分子0がT₁、分子1,2がS₀、分子3がT₁ |

---

## 付録B: 単位系変換

### B.1 eV/ℏ → fs⁻¹ 変換

$$
\gamma [\text{fs}^{-1}] = \frac{\gamma [\text{eV/ℏ}]}{\hbar [\text{eV·fs}]} = \frac{\gamma}{0.6582119569}
$$

### B.2 eV/ℏ → s⁻¹ 変換

$$
\gamma [\text{s}^{-1}] = \gamma [\text{eV/ℏ}] \times 1.519 \times 10^{15}
$$

### B.3 主要パラメータの単位変換表

| パラメータ | eV/ℏ | fs⁻¹ | s⁻¹ | 特性時間 |
|-----------|------|-------|-----|---------|
| $\gamma_{\text{TTA}}$ | 0.05 | 0.0760 | $7.60 \times 10^{13}$ | 13.2 fs |
| $\Gamma_{\text{fl}}$ | 0.01 | 0.0152 | $1.52 \times 10^{13}$ | 65.8 fs |
| $\Gamma_{\text{ph}}$ | $10^{-6}$ | $1.52 \times 10^{-6}$ | $1.52 \times 10^{9}$ | $6.58 \times 10^{5}$ fs |
| $k_{\text{IC}}$ | 0.005 | $7.60 \times 10^{-3}$ | $7.60 \times 10^{12}$ | 131.6 fs |
| $k_{\text{ISC}}^{S\to T}$ | 0.003 | $4.56 \times 10^{-3}$ | $4.56 \times 10^{12}$ | 219.4 fs |
| $k_{\text{ISC}}^{T\to S}$ | $10^{-5}$ | $1.52 \times 10^{-5}$ | $1.52 \times 10^{10}$ | $6.58 \times 10^{4}$ fs |

---

## 付録C: エラーハンドリング設計

### C.1 カスタム例外クラス

```python
class PhysicsViolationError(Exception):
    """物理法則の違反を検出したときに送出する例外"""
    pass

class NumericalInstabilityError(Exception):
    """数値的不安定性を検出したときに送出する例外"""
    pass

class ParameterValidationError(ValueError):
    """パラメータの検証に失敗したときに送出する例外"""
    pass
```

### C.2 エラーハンドリングの原則

| 状況 | 対応 | Fallbackの有無 |
|------|------|---------------|
| パラメータ値が物理的範囲外 | `ParameterValidationError` を raise | **なし** |
| $\text{Tr}[\hat{\rho}] \neq 1$ | `PhysicsViolationError` を raise | **なし** |
| 負の固有値（$< -10^{-10}$） | `PhysicsViolationError` を raise | **なし** |
| ODE積分の失敗 | `RuntimeError` を raise | **なし** |
| ユニタリ性の破れ | `NumericalInstabilityError` を raise | **なし** |
| ゲート分解の失敗 | `RuntimeError` を raise | **なし** |

**重要**: いかなる場合もFallback処理（「適当な値への置き換え」、「結果の補正」、「エラーの黙殺」等）は行わない。全てのエラーは明示的に例外として送出する。

### C.3 ログ出力の設計

```python
import logging

logger = logging.getLogger('gksl_simulator')

# 各ステップで出力する情報
logger.info(f"Step {step}/{N_steps}: "
            f"Tr[ρ]={tr:.12f}, "
            f"min(λ)={min_eig:.2e}, "
            f"S={entropy:.6f}, "
            f"P={purity:.6f}")

# 警告（エラーではないが注意が必要な場合）
IF entropy_decrease > 1e-10:
    logger.warning(f"Step {step}: Entropy decreased by {entropy_decrease:.2e}")
```

---

## 付録D: Lindblad超演算子の要素展開

### D.1 蛍光 $\hat{L} = |0\rangle\langle 2|$（単一分子 $3 \times 3$）

$$
\hat{L}^\dagger\hat{L} = |2\rangle\langle 2|
$$

$$
\hat{L}\hat{\rho}\hat{L}^\dagger = \rho_{22} |0\rangle\langle 0|
$$

対角要素の変化率:

$$
\dot{\rho}_{00}^{(\text{fl})} = +\Gamma_{\text{fl}} \rho_{22}
$$

$$
\dot{\rho}_{11}^{(\text{fl})} = 0
$$

$$
\dot{\rho}_{22}^{(\text{fl})} = -\Gamma_{\text{fl}} \rho_{22}
$$

オフ対角要素の変化率:

$$
\dot{\rho}_{01}^{(\text{fl})} = 0, \quad \dot{\rho}_{02}^{(\text{fl})} = -\frac{\Gamma_{\text{fl}}}{2} \rho_{02}, \quad \dot{\rho}_{12}^{(\text{fl})} = -\frac{\Gamma_{\text{fl}}}{2} \rho_{12}
$$

物理的意味: $S_1$ 準位の確率が $S_0$ に移行し、$S_1$ に関連するコヒーレンスが減衰する。

### D.2 燐光 $\hat{L} = |0\rangle\langle 1|$（単一分子 $3 \times 3$）

$$
\hat{L}^\dagger\hat{L} = |1\rangle\langle 1|
$$

対角要素の変化率:

$$
\dot{\rho}_{00}^{(\text{ph})} = +\Gamma_{\text{ph}} \rho_{11}
$$

$$
\dot{\rho}_{11}^{(\text{ph})} = -\Gamma_{\text{ph}} \rho_{11}
$$

$$
\dot{\rho}_{22}^{(\text{ph})} = 0
$$

オフ対角要素の変化率:

$$
\dot{\rho}_{01}^{(\text{ph})} = -\frac{\Gamma_{\text{ph}}}{2} \rho_{01}, \quad \dot{\rho}_{02}^{(\text{ph})} = 0, \quad \dot{\rho}_{12}^{(\text{ph})} = -\frac{\Gamma_{\text{ph}}}{2} \rho_{12}
$$

### D.3 ISC S₁→T₁ $\hat{L} = |1\rangle\langle 2|$（単一分子 $3 \times 3$）

$$
\hat{L}^\dagger\hat{L} = |2\rangle\langle 2|
$$

対角要素の変化率:

$$
\dot{\rho}_{00}^{(\text{ISC})} = 0
$$

$$
\dot{\rho}_{11}^{(\text{ISC})} = +k_{\text{ISC}}^{S\to T} \rho_{22}
$$

$$
\dot{\rho}_{22}^{(\text{ISC})} = -k_{\text{ISC}}^{S\to T} \rho_{22}
$$

オフ対角要素の変化率:

$$
\dot{\rho}_{01}^{(\text{ISC})} = -\frac{k_{\text{ISC}}^{S\to T}}{2} \rho_{01} \quad (\text{注: }\dot{\rho}_{01} = 0 \text{ ではない。正しくは...})
$$

**修正**: ISC $|1\rangle\langle 2|$ の場合:

$$
\hat{L}^\dagger\hat{L} = |2\rangle\langle 2|, \quad \hat{L}\hat{L}^\dagger = |1\rangle\langle 1|
$$

$$
\dot{\rho}_{01} = -\frac{k_{\text{ISC}}^{S\to T}}{2}(0 \cdot \rho_{01}) = 0
$$

いいえ、正確に展開する:

$$
\mathcal{D}[|1\rangle\langle 2|][\hat{\rho}]_{01} = (|1\rangle\langle 2|\hat{\rho}|2\rangle\langle 1|)_{01} - \frac{1}{2}(|2\rangle\langle 2|\hat{\rho})_{01} - \frac{1}{2}(\hat{\rho}|2\rangle\langle 2|)_{01}
$$

$$
= 0 - \frac{1}{2} \cdot 0 - \frac{1}{2} \rho_{02}\delta_{21} = 0
$$

結果:

$$
\dot{\rho}_{01}^{(\text{ISC})} = 0, \quad \dot{\rho}_{02}^{(\text{ISC})} = -\frac{k_{\text{ISC}}^{S\to T}}{2}\rho_{02}, \quad \dot{\rho}_{12}^{(\text{ISC})} = -\frac{k_{\text{ISC}}^{S\to T}}{2}\rho_{12}
$$

### D.4 TTA $\hat{L} = |20\rangle\langle 11|$（2分子 $9 \times 9$ 部分空間）

$$
\hat{L}^\dagger\hat{L} = |11\rangle\langle 11|
$$

対角要素の変化率:

$$
\dot{\rho}_{20,20}^{(\text{TTA})} = +\frac{\gamma_{\text{TTA}}}{2} \rho_{11,11}
$$

$$
\dot{\rho}_{11,11}^{(\text{TTA})} = -\frac{\gamma_{\text{TTA}}}{2} \rho_{11,11}
$$

$$
\dot{\rho}_{mn,mn}^{(\text{TTA})} = 0 \quad (mn \neq 11, 20)
$$

オフ対角要素: $|11\rangle$ に関連するコヒーレンス（例: $\rho_{11,10}, \rho_{11,01}$ 等）が $\gamma_{\text{TTA}}/4$ の速度で減衰する。

---

**文書終了**

作成日: 2026年2月12日  
バージョン: 1.0.0  
対象リポジトリ: nobkt/mqt-qudits  
ライセンス: MIT License

