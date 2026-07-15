# TTA-UC現象 GKSL-Lindblad量子ダイナミクス 実装詳細仕様書

## 文書情報

**作成日**: 2026年2月12日  
**最終更新日**: 2026年2月12日  
**バージョン**: 1.1.1  
**対象ノートブック**: `tutorials/quantum_dynamics_complete_comparison.ipynb`  
**理論基礎文書**: `tutorials/doc/GKSL/TTA-UC現象のGKSL-Lindblad量子ダイナミクス完全理論書.md`  
**目的**: 既存ノートブック `quantum_dynamics_complete_comparison.ipynb` と同一の仕様（パラメータ、データ構造、比較フレームワーク、可視化）で、TTA-UC現象のGKSL-Lindblad量子ダイナミクスを実装するための省略無しの詳細仕様

---

## 目次

1. [本仕様書の位置づけと目的](#1-本仕様書の位置づけと目的)
   - 1.1 現行ノートブックとの関係
   - 1.2 本仕様書が定義するもの
   - 1.3 現行ノートブック実装との根本的な違い
   - 1.4 厳密性の原則
   - 1.5 開放量子系の数学的基礎
2. [現行ノートブックの仕様分析](#2-現行ノートブックの仕様分析)
3. [GKSL-Lindblad拡張の数学的定式化](#3-gksl-lindblad拡張の数学的定式化)
   - 3.0 GKSL定理の正式な定式化
   - 3.1 GKSL-Lindblad方程式
   - 3.2–3.5 ハミルトニアン・Lindblad演算子・完全なGKSL方程式
4. [物理パラメータの完全定義](#4-物理パラメータの完全定義)
   - 4.1–4.2 パラメータクラス・物理的根拠
   - 4.2.5 実験的パラメータ範囲との対応
   - 4.3 パラメータ間の整合性チェック
5. [シナリオ1: 古典計算・ボソン無しのGKSL実装仕様](#5-シナリオ1-古典計算ボソン無しのgksl実装仕様)
6. [シナリオ2: 古典計算・ボソン有りのGKSL実装仕様](#6-シナリオ2-古典計算ボソン有りのgksl実装仕様)
   - 6.3.5–6.3.10 Peierls型結合・電子-光子結合・完全ハミルトニアン・スペクトル密度・有限温度効果
   - 6.8 追加の数値手法（ボソン有りモデル用）
7. [シナリオ3: Qubit量子計算・ボソン無しのGKSL実装仕様](#7-シナリオ3-qubit量子計算ボソン無しのgksl実装仕様)
   - 7.2.1 物理的部分空間と射影演算子
   - 7.6 ゲート分解戦略
   - 7.10 ハードウェアノイズモデル（Qubit GKSL）
8. [シナリオ4: Qubit量子計算・ボソン有りのGKSL実装仕様](#8-シナリオ4-qubit量子計算ボソン有りのgksl実装仕様)
9. [シナリオ5: Qudit量子計算・ボソン無しのGKSL実装仕様](#9-シナリオ5-qudit量子計算ボソン無しのgksl実装仕様)
   - 9.2.1 MQT-Quditsの基本ゲートセット
   - 9.9 ハードウェアノイズモデル（Qudit GKSL）
10. [シナリオ6: Qudit量子計算・ボソン有りのGKSL実装仕様](#10-シナリオ6-qudit量子計算ボソン有りのgksl実装仕様)
11. [データ構造と出力フォーマット仕様](#11-データ構造と出力フォーマット仕様)
12. [比較・可視化フレームワーク仕様](#12-比較可視化フレームワーク仕様)
    - 12.4 量子回路可視化仕様
    - 12.5 スケーラビリティ比較
13. [検証仕様](#13-検証仕様)
    - 13.2.3 Trotter分解誤差の評価
    - 13.3.3 Stinespring実装の忠実度検証
14. [現行ノートブックとの差異の正確な記述](#14-現行ノートブックとの差異の正確な記述)
15. [実装ロードマップ](#15-実装ロードマップ)
16. [今後の拡張可能性](#16-今後の拡張可能性)

---

## 1. 本仕様書の位置づけと目的

### 1.1 現行ノートブックとの関係

現行ノートブック `tutorials/quantum_dynamics_complete_comparison.ipynb` は、4分子直線配置系における分子三重項状態の量子ダイナミクスを **閉じた系（ユニタリ時間発展）** として実装している。具体的には：

- **時間発展**: Schrödinger方程式 $i\hbar \frac{d|\psi\rangle}{dt} = \hat{H}|\psi\rangle$ に基づく状態ベクトルの発展
- **TTA過程**: ハミルトニアン項 $\hat{H}_{\text{TTA}}$ として可逆的に実装
- **散逸**: 未実装（蛍光、燐光、内部転換、項間交差の散逸項は含まれていない）
- **状態記述**: 純粋状態 $|\psi\rangle$（密度行列 $\hat{\rho}$ ではない）

### 1.2 本仕様書が定義するもの

本仕様書は、上記ノートブックと **同一のフレームワーク**（パラメータクラス、シミュレータインターフェース、出力データ構造、可視化関数）を用いて、GKSL-Lindblad方程式に基づく **開放量子系ダイナミクス** を実装するための完全仕様を定義する。

### 1.3 現行ノートブック実装との根本的な違い

| 項目 | 現行ノートブック | GKSL-Lindblad拡張 |
|------|----------------|-------------------|
| 時間発展方程式 | Schrödinger方程式 | GKSL-Lindblad方程式 |
| 状態記述 | 状態ベクトル $\|\psi\rangle \in \mathbb{C}^{81}$ | 密度行列 $\hat{\rho} \in \mathbb{C}^{81 \times 81}$ |
| TTA過程 | ハミルトニアン項（可逆） | Lindblad散逸項（不可逆） |
| エネルギー移動 | ハミルトニアン項（コヒーレント） | ハミルトニアン項（コヒーレント、変更なし） |
| 蛍光・燐光 | 未実装 | Lindblad散逸項として実装 |
| 内部転換・ISC | 未実装 | Lindblad散逸項として実装 |
| エントロピー | 不変（$S = 0$） | 増大（$dS/dt \geq 0$） |
| 量子回路（Qubit/Qudit） | ユニタリゲートのみ | ユニタリゲート + Stinespring dilation |

### 1.4 厳密性の原則

本仕様書は以下の原則に従う：

**✅ 許容される手法**：
- GKSL定理に基づく数学的に厳密な導出
- Stinespring dilationによる非ユニタリ演算の量子回路表現
- 実験的に測定可能なパラメータの使用
- 制御可能な近似誤差（Trotter分解）の明示的評価

**❌ 禁止される手法**：
- ヒューリスティックな近似や経験的フィッティング
- Fallback処理（計算失敗時の「適当な値」への置き換え）
- 物理的根拠のない簡略化
- ごまかしや真実を隠蔽する記述

### 1.5 開放量子系の数学的基礎

本仕様書のGKSL-Lindblad定式化の基礎となる開放量子系の数学的枠組みを以下にまとめる。

#### 1.5.1 密度演算子の公理

量子系の状態を記述する密度演算子 $\hat{\rho}$ は以下の3条件を満たす：

1. **エルミート性**: $\hat{\rho} = \hat{\rho}^\dagger$
2. **正の半正定値性**: $\langle\psi|\hat{\rho}|\psi\rangle \geq 0 \quad \forall |\psi\rangle$（すべての固有値が非負）
3. **トレース条件**: $\text{Tr}[\hat{\rho}] = 1$

#### 1.5.2 純粋状態と混合状態の区別

- **純粋状態**: $\hat{\rho}^2 = \hat{\rho}$ が成立。$\hat{\rho} = |\psi\rangle\langle\psi|$ と書ける。$\text{Tr}[\hat{\rho}^2] = 1$。
- **混合状態**: $\hat{\rho}^2 \neq \hat{\rho}$。$\text{Tr}[\hat{\rho}^2] < 1$。

閉じた系のユニタリ発展では純粋状態は純粋状態のまま保たれるが、開放量子系では散逸により混合状態に遷移する。

#### 1.5.3 系-環境分離

全系のヒルベルト空間をテンソル積構造で分離する：

$$
\mathcal{H}_{\text{total}} = \mathcal{H}_S \otimes \mathcal{H}_E
$$

ここで $\mathcal{H}_S$ は系（TTA-UC分子系）、$\mathcal{H}_E$ は環境（溶媒、放射場、フォノン浴等）のヒルベルト空間。

#### 1.5.4 縮約密度演算子と部分トレース

環境の自由度を部分トレースで除去し、系のみの状態を記述する：

$$
\hat{\rho}_S = \text{Tr}_E[\hat{\rho}_{SE}]
$$

部分トレースの定義：

$$
\langle m|\hat{\rho}_S|n\rangle = \sum_k \langle m, k|\hat{\rho}_{SE}|n, k\rangle
$$

ここで $\{|k\rangle\}$ は環境の正規直交基底。

#### 1.5.5 CPTP写像

系の時間発展を記述する量子チャネル $\mathcal{E}$ は **完全正値トレース保存（CPTP）写像** でなければならない：

1. **線形性**: $\mathcal{E}[\alpha\hat{\rho}_1 + \beta\hat{\rho}_2] = \alpha\mathcal{E}[\hat{\rho}_1] + \beta\mathcal{E}[\hat{\rho}_2]$
2. **トレース保存**: $\text{Tr}[\mathcal{E}[\hat{\rho}]] = \text{Tr}[\hat{\rho}] = 1$
3. **完全正値性**: $(\mathcal{E} \otimes \mathcal{I}_n)[\hat{\rho}_{SA}] \geq 0$ が任意の補助系 $A$（$n$ 次元）と任意の $\hat{\rho}_{SA} \geq 0$ に対して成立

#### 1.5.6 Kraus表現

任意のCPTP写像は **Kraus表現**（演算子和表現）を持つ：

$$
\mathcal{E}[\hat{\rho}] = \sum_k \hat{K}_k \hat{\rho} \hat{K}_k^\dagger
$$

ここで Kraus 演算子 $\{\hat{K}_k\}$ は完全性条件を満たす：

$$
\sum_k \hat{K}_k^\dagger \hat{K}_k = \hat{I}
$$

この条件がトレース保存を保証する。

#### 1.5.7 マルコフ近似

本仕様書のGKSL定式化はマルコフ近似に基づく。この近似が妥当であるための条件：

1. **弱結合条件**: 系-環境結合が系の内部エネルギースケールに比べて十分弱い
2. **時間スケール分離**: 環境の相関時間 $\tau_E$ が系の特性時間 $\tau_S$ に比べて十分短い
$$
\tau_E \ll \tau_S
$$

この条件下で、系の時間発展はメモリ効果を持たず、過去の履歴に依存しない。

#### 1.5.8 半群性

マルコフ近似の下で、量子力学的時間発展は **量子力学的半群** の性質を持つ：

$$
\mathcal{E}_{t+s} = \mathcal{E}_t \circ \mathcal{E}_s \quad (t, s \geq 0)
$$

この半群性が、GKSL-Lindblad方程式の生成子（リンドブラディアン $\mathcal{L}$）の存在を保証する：

$$
\mathcal{E}_t = e^{\mathcal{L} t}
$$

---

## 2. 現行ノートブックの仕様分析

### 2.1 PhysicalParametersクラスの現行仕様

```python
class PhysicalParameters:
    def __init__(self):
        self.N_molecules = 4          # 分子数
        self.E_T = 1.5                # 三重項エネルギー (eV)
        self.E_S = 3.0                # 一重項エネルギー (eV)
        self.V = 0.1                  # エネルギー移動積分 (eV)
        self.J = 0.05                 # TTA相互作用定数 (eV) ← GKSL版では廃止
        self.Gamma_fl = 0.01          # 蛍光放出速度 (fs^-1) ← GKSL版で再定義
        self.hbar = 0.6582119569      # 換算プランク定数 (eV·fs)
        self.neighbors = [(0, 1), (1, 2), (2, 3)]
        self.T_total = 100.0          # 総時間 (fs)
        self.N_steps = 20             # トロッターステップ数
        self.dt = self.T_total / self.N_steps
        self.initial_state_type = 'edge_triplet'
```

### 2.2 シミュレータのインターフェース仕様

#### 2.2.1 現行の出力データ構造

現行ノートブックの全シミュレータは以下の形式の辞書を返す：

```python
{
    'times': List[float],                    # 時間点 [0, dt, 2*dt, ..., T_total]
    'populations': List[Dict[str, float]],   # 各時刻の個体数
    # populations[i] = {'N_S0': float, 'N_T1': float, 'N_S1': float}
    'per_molecule_populations': List[Dict[str, np.ndarray]],  # 各分子ごとの個体数
    # per_molecule_populations[i] = {
    #     'S0_per_mol': np.ndarray(4,),
    #     'T1_per_mol': np.ndarray(4,),
    #     'S1_per_mol': np.ndarray(4,)
    # }
    'elapsed_time': float,                   # 実行時間（秒）
    'method': str,                           # 手法名
    # 古典シミュレータ追加:
    'state_final': np.ndarray,               # 最終状態ベクトル (dim=81)
    # Qubitシミュレータ追加:
    'step_circuit': QuantumCircuit,           # 1トロッターステップの回路
    'total_gates': int,                      # 総ゲート数
    'total_depth': int,                      # 総回路深さ
    'gates_per_step': int,
    'depth_per_step': int,
    'shots': int,
    # Quditシミュレータ追加:
    'step_circuit': QuantumCircuit,           # MQT-Qudits回路
    'N_steps': int,
}
```

#### 2.2.2 現行のシミュレータメソッド

- `ClassicalSuzukiTrotterSimulator.simulate(T_total, N_steps, initial_state_type)` → Dict
- `QubitMolecularDynamicsSimulator.simulate(T_total, N_steps, initial_state_type, shots)` → Dict
- `SuzukiTrotterMQTQuditSimulator.simulate_shot_based(T_total, N_steps, initial_state_type, track_dynamics, shots)` → Dict

### 2.3 初期状態の仕様

初期状態 `'edge_triplet'` は以下の状態ベクトルを意味する：

$$
|\psi(0)\rangle = |T_1\rangle_0 \otimes |S_0\rangle_1 \otimes |S_0\rangle_2 \otimes |T_1\rangle_3 = |1,0,0,1\rangle
$$

Qutrit基底でのインデックス：

$$
\text{idx} = 1 \times 3^3 + 0 \times 3^2 + 0 \times 3^1 + 1 \times 3^0 = 28
$$

GKSL版では、この初期状態を密度行列に変換する：

$$
\hat{\rho}(0) = |\psi(0)\rangle\langle\psi(0)| = |1001\rangle\langle 1001|
$$

これは $81 \times 81$ 行列で、$(28, 28)$ 要素のみが1で、他は全て0。

### 2.4 現行の時間発展アルゴリズム

#### 2.4.1 2次対称鈴木トロッター分解（per-pair）

1ステップの時間発展演算子：

$$
\hat{U}(\Delta t) = \underbrace{\prod_{i=N-1}^{0} U_{H_0}^{(i)}(\Delta t/2)}_{\text{backward H0}} \cdot \underbrace{\prod_{(i,j) \text{ rev}} U_{\text{tr}}^{(ij)}(\Delta t/2)}_{\text{backward transfer}} \cdot \underbrace{\prod_{(i,j) \text{ rev}} U_{\text{TTA}}^{(ij)}(\Delta t/2)}_{\text{backward TTA}}
$$

$$
\times \underbrace{\prod_{(i,j)} U_{\text{TTA}}^{(ij)}(\Delta t/2)}_{\text{forward TTA}} \cdot \underbrace{\prod_{(i,j)} U_{\text{tr}}^{(ij)}(\Delta t/2)}_{\text{forward transfer}} \cdot \underbrace{\prod_{i=0}^{N-1} U_{H_0}^{(i)}(\Delta t/2)}_{\text{forward H0}}
$$

各ユニタリは `scipy.linalg.expm` で厳密に計算。

#### 2.4.2 適用順序（コード準拠）

Forward半分:
1. H0: mol_idx = 0, 1, 2, 3（各分子のオンサイトエネルギー）
2. H_transfer: ペア (0,1), (1,2), (2,3)
3. H_TTA: ペア (0,1), (1,2), (2,3)

Backward半分（逆順）:
4. H_TTA: ペア (2,3), (1,2), (0,1)
5. H_transfer: ペア (2,3), (1,2), (0,1)
6. H0: mol_idx = 3, 2, 1, 0

### 2.5 個体数の計算方法

#### 2.5.1 現行方式（状態ベクトルから）

```python
for idx in range(dim):    # dim = 81
    prob = |state[idx]|^2
    config = index_to_config(idx, N=4, d=3)  # e.g., [1, 0, 0, 1]
    for mol_idx, level in enumerate(config):
        if level == 0: N_S0 += prob
        elif level == 1: N_T1 += prob
        elif level == 2: N_S1 += prob
```

#### 2.5.2 GKSL版（密度行列から）

$$
N_{X}(t) = \sum_{i=0}^{N-1} \text{Tr}[\hat{P}_{X}^{(i)} \hat{\rho}(t)]
$$

ここで $\hat{P}_{X}^{(i)}$ は分子 $i$ の状態 $X$ への射影演算子。

具体的実装：

$$
N_{S_0}(t) = \sum_{i=0}^{3} \text{Tr}\left[\left(\hat{I}^{\otimes i} \otimes |0\rangle\langle 0| \otimes \hat{I}^{\otimes (3-i)}\right) \hat{\rho}(t)\right]
$$

$$
N_{T_1}(t) = \sum_{i=0}^{3} \text{Tr}\left[\left(\hat{I}^{\otimes i} \otimes |1\rangle\langle 1| \otimes \hat{I}^{\otimes (3-i)}\right) \hat{\rho}(t)\right]
$$

$$
N_{S_1}(t) = \sum_{i=0}^{3} \text{Tr}\left[\left(\hat{I}^{\otimes i} \otimes |2\rangle\langle 2| \otimes \hat{I}^{\otimes (3-i)}\right) \hat{\rho}(t)\right]
$$

保存則：散逸が無い場合 $N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = 4$。散逸がある場合も $\text{Tr}[\hat{\rho}] = 1$ から $N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = N = 4$ が成立。

### 2.6 可視化関数の仕様

#### 2.6.1 plot_population_dynamics

入力: `results: Dict` （上記出力データ構造）, `title: str`  
出力: matplotlib Figure  
プロット内容:
- X軸: Time (fs)
- Y軸: Population (0〜4.5)
- 3本の曲線: $N_{S_0}$ (青, 'o'), $N_{T_1}$ (赤, 's'), $N_{S_1}$ (緑, '^')

#### 2.6.2 plot_per_molecule_populations

入力: `results: Dict`, `title: str`  
出力: matplotlib Figure (2×2 subplot)  
プロット内容: 各分子の $S_0$, $T_1$, $S_1$ の時間発展（0〜1.1）

---

## 3. GKSL-Lindblad拡張の数学的定式化

### 3.0 GKSL定理の正式な定式化

#### 3.0.1 定理の帰属

GKSL（Gorini-Kossakowski-Sudarshan-Lindblad）定理は、Gorini, Kossakowski, Sudarshan (1976) と Lindblad (1976) により独立に証明された。

**GKSL定理**: 有限次元ヒルベルト空間 $\mathcal{H}$（$\dim \mathcal{H} = N$）上の密度演算子の時間発展が、CPTP量子力学的半群 $\{\mathcal{E}_t\}_{t \geq 0}$ を成すとき、その生成子 $\mathcal{L}$（リンドブラディアン）は以下の形式に限定される：

$$
\frac{d\hat{\rho}}{dt} = \mathcal{L}[\hat{\rho}] = -\frac{i}{\hbar}[\hat{H}, \hat{\rho}] + \sum_{\alpha=1}^{N^2-1} \gamma_\alpha \left(\hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger - \frac{1}{2}\{\hat{L}_\alpha^\dagger \hat{L}_\alpha, \hat{\rho}\}\right)
$$

ここで：
- $\hat{H} = \hat{H}^\dagger$ はエルミート演算子（系のハミルトニアン + Lamb シフト補正）
- $\{\hat{L}_\alpha\}$ はLindblad演算子（ジャンプ演算子）、$\alpha = 1, \ldots, N^2 - 1$
- $\gamma_\alpha \geq 0$ は非負の散逸速度定数
- $\{\cdot, \cdot\}$ は反交換子: $\{\hat{A}, \hat{B}\} = \hat{A}\hat{B} + \hat{B}\hat{A}$

#### 3.0.2 完全正値性の保証

GKSL形式は、マルコフ的時間発展がCPTP写像であることを **自動的に保証する**。これは以下を意味する：

1. **密度行列の正定値性が保たれる**: $\hat{\rho}(0) \geq 0 \Rightarrow \hat{\rho}(t) \geq 0 \quad \forall t \geq 0$
2. **拡張系でも正定値性が保たれる**: 任意の補助系を追加しても正定値性は破れない

逆に、GKSL形式でない一般的なマスター方程式（例：Redfield方程式）はトレース保存を満たしても完全正値性を保証しないため、非物理的な負の確率を生じ得る。

#### 3.0.3 トレース保存の証明概略

$\frac{d}{dt}\text{Tr}[\hat{\rho}]$ を計算する：

$$
\frac{d}{dt}\text{Tr}[\hat{\rho}] = -\frac{i}{\hbar}\text{Tr}[[\hat{H}, \hat{\rho}]] + \sum_\alpha \gamma_\alpha \text{Tr}\left[\hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger - \frac{1}{2}\hat{L}_\alpha^\dagger \hat{L}_\alpha \hat{\rho} - \frac{1}{2}\hat{\rho}\hat{L}_\alpha^\dagger \hat{L}_\alpha\right]
$$

第1項: $\text{Tr}[[\hat{H}, \hat{\rho}]] = \text{Tr}[\hat{H}\hat{\rho}] - \text{Tr}[\hat{\rho}\hat{H}] = 0$（トレースの巡回性）

第2項の各 $\alpha$ について：$\text{Tr}[\hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger] = \text{Tr}[\hat{L}_\alpha^\dagger \hat{L}_\alpha \hat{\rho}]$（巡回性）より：

$$
\text{Tr}[\hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger] - \frac{1}{2}\text{Tr}[\hat{L}_\alpha^\dagger \hat{L}_\alpha \hat{\rho}] - \frac{1}{2}\text{Tr}[\hat{\rho}\hat{L}_\alpha^\dagger \hat{L}_\alpha] = 0
$$

よって $\frac{d}{dt}\text{Tr}[\hat{\rho}] = 0$ が任意の $\hat{\rho}$ に対して成立する。 $\square$

### 3.1 GKSL-Lindblad方程式

4分子TTA-UC系の完全なGKSL-Lindblad方程式：

$$
\frac{d\hat{\rho}}{dt} = \underbrace{-\frac{i}{\hbar}[\hat{H}_{\text{sys}}, \hat{\rho}]}_{\text{ユニタリ部分}} + \underbrace{\sum_{\alpha} \gamma_\alpha \mathcal{D}[\hat{L}_\alpha][\hat{\rho}]}_{\text{散逸部分}}
$$

ここでLindblad超演算子は：

$$
\mathcal{D}[\hat{L}][\hat{\rho}] = \hat{L}\hat{\rho}\hat{L}^\dagger - \frac{1}{2}\hat{L}^\dagger\hat{L}\hat{\rho} - \frac{1}{2}\hat{\rho}\hat{L}^\dagger\hat{L}
$$

### 3.2 ハミルトニアン（ユニタリ部分）

#### 3.2.1 系のハミルトニアン

現行ノートブックとの**重要な違い**：GKSL版ではTTAをハミルトニアンから除外し、Lindblad散逸項として実装する。

$$
\hat{H}_{\text{sys}} = \hat{H}_0 + \hat{H}_{\text{transfer}}
$$

**$\hat{H}_{\text{TTA}}$ はハミルトニアンから除外する**。

#### 3.2.2 オンサイトエネルギー $\hat{H}_0$（現行と同一）

$$
\hat{H}_0 = \sum_{i=0}^{3} \left( E_T |1\rangle_i\langle 1| + E_S |2\rangle_i\langle 2| \right)
$$

単一分子の行列表現（$3 \times 3$）：

$$
\hat{H}_0^{(i)} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1.5 & 0 \\ 0 & 0 & 3.0 \end{pmatrix} \text{ eV}
$$

全空間への拡張：

$$
\hat{H}_0^{(i)} \to \hat{I}_3^{\otimes i} \otimes \hat{H}_0^{(i)} \otimes \hat{I}_3^{\otimes (3-i)}
$$

#### 3.2.3 エネルギー移動 $\hat{H}_{\text{transfer}}$（現行と同一）

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1| \right)
$$

2分子部分空間の $9 \times 9$ 行列表現（基底 $\{|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle\}$）：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V \begin{pmatrix} 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{pmatrix}
$$

### 3.3 Lindblad散逸演算子の完全リスト

以下のLindblad演算子はすべて全空間 $\mathbb{C}^{81}$ 上の演算子として構築する。局所演算子から全空間への拡張は、テンソル積 $\hat{I}_3^{\otimes i} \otimes \hat{L}_{\text{local}} \otimes \hat{I}_3^{\otimes (3-i)}$ で行う。

#### 3.3.1 TTA過程のLindblad演算子

隣接ペア $(i, j) \in \{(0,1), (1,2), (2,3)\}$ に対して**2つの**Lindblad演算子を定義する：

**$\hat{L}_{\text{TTA},1}^{(ij)}$: 分子 $i$ が $S_1$ に昇格、分子 $j$ が $S_0$ に脱励起**

$$
\hat{L}_{\text{TTA},1}^{(ij)} = |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|
$$

物理的意味：$|T_1\rangle_i|T_1\rangle_j \to |S_1\rangle_i|S_0\rangle_j$

2分子部分空間の $9 \times 9$ 行列表現：

$$
\hat{L}_{\text{TTA},1}^{(ij)} = \begin{pmatrix} 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{pmatrix}
$$

（$|20\rangle\langle 11|$、つまり行6列4の要素のみ1、0-indexed）

**$\hat{L}_{\text{TTA},2}^{(ij)}$: 分子 $j$ が $S_1$ に昇格、分子 $i$ が $S_0$ に脱励起**

$$
\hat{L}_{\text{TTA},2}^{(ij)} = |0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1|
$$

物理的意味：$|T_1\rangle_i|T_1\rangle_j \to |S_0\rangle_i|S_1\rangle_j$

$9 \times 9$ 行列表現：

$$
\hat{L}_{\text{TTA},2}^{(ij)} = \begin{pmatrix} 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{pmatrix}
$$

（$|02\rangle\langle 11|$、つまり行2列4の要素のみ1、0-indexed）

**TTA過程の散逸項**（ペア $(i,j)$ に対して）：

$$
\mathcal{L}_{\text{TTA}}^{(ij)}[\hat{\rho}] = \frac{\gamma_{\text{TTA}}}{2} \left( \mathcal{D}[\hat{L}_{\text{TTA},1}^{(ij)}][\hat{\rho}] + \mathcal{D}[\hat{L}_{\text{TTA},2}^{(ij)}][\hat{\rho}] \right)
$$

$\gamma_{\text{TTA}}/2$ の因子は、2つのLindblad演算子が同じ物理過程の2つのチャネルを表すため。

**理論書との表記対応**: 理論書では $\hat{L}_{\text{TTA}} = \sqrt{\gamma_{\text{TTA}}/2}\,|20\rangle\langle 11|$ のように $\sqrt{\gamma_{\text{TTA}}/2}$ を演算子側に含める。本仕様書は $\hat{L}$ を無次元の遷移演算子として定義し、係数を散逸項にまとめる。どちらも同値なGKSL形式である。

**スピン統計因子の取り込み**: 理論書で触れたスピン統計因子 $f_{\text{spin}} \in \{1/9, 1/5\}$ を考慮する場合は $\gamma_{\text{TTA,eff}} = f_{\text{spin}}\gamma_{\text{TTA}}$ を使用し、`GKSLPhysicalParameters.gamma_TTA` には既に $f_{\text{spin}}$ を含めた値を設定する（重複計上を避ける）。

#### 3.3.2 蛍光発光のLindblad演算子

分子 $i$ に対して（$i = 0, 1, 2, 3$）：

$$
\hat{L}_{\text{fl}}^{(i)} = |0\rangle_i\langle 2|
$$

物理的意味：$|S_1\rangle_i \to |S_0\rangle_i + h\nu_{\text{fl}}$

単一分子の $3 \times 3$ 行列表現：

$$
\hat{L}_{\text{fl}}^{(i)} = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

散逸項：

$$
\mathcal{L}_{\text{fl}}^{(i)}[\hat{\rho}] = \Gamma_{\text{fl}} \cdot \mathcal{D}[\hat{L}_{\text{fl}}^{(i)}][\hat{\rho}]
$$

#### 3.3.3 燐光発光のLindblad演算子

分子 $i$ に対して（$i = 0, 1, 2, 3$）：

$$
\hat{L}_{\text{ph}}^{(i)} = |0\rangle_i\langle 1|
$$

物理的意味：$|T_1\rangle_i \to |S_0\rangle_i + h\nu_{\text{ph}}$

単一分子の $3 \times 3$ 行列表現：

$$
\hat{L}_{\text{ph}}^{(i)} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

散逸項：

$$
\mathcal{L}_{\text{ph}}^{(i)}[\hat{\rho}] = \Gamma_{\text{ph}} \cdot \mathcal{D}[\hat{L}_{\text{ph}}^{(i)}][\hat{\rho}]
$$

#### 3.3.4 内部転換のLindblad演算子

分子 $i$ に対して（$i = 0, 1, 2, 3$）：

$$
\hat{L}_{\text{IC}}^{(i)} = |0\rangle_i\langle 2|
$$

$$
\mathcal{L}_{\text{IC}}^{(i)}[\hat{\rho}] = k_{\text{IC}} \cdot \mathcal{D}[\hat{L}_{\text{IC}}^{(i)}][\hat{\rho}]
$$

**注意**: $\hat{L}_{\text{IC}}^{(i)}$ と $\hat{L}_{\text{fl}}^{(i)}$ は同じ行列形式を持つ。これらの効果は速度定数の和として扱うことも可能：実効蛍光速度 $\Gamma_{\text{fl,eff}} = \Gamma_{\text{fl}} + k_{\text{IC}}$。ただし物理的に区別するためには別個のLindblad演算子として扱う方が正確である。

#### 3.3.5 項間交差 S₁→T₁ のLindblad演算子

分子 $i$ に対して（$i = 0, 1, 2, 3$）：

$$
\hat{L}_{\text{ISC}}^{S \to T,(i)} = |1\rangle_i\langle 2|
$$

物理的意味：$|S_1\rangle_i \to |T_1\rangle_i$

単一分子の $3 \times 3$ 行列表現：

$$
\hat{L}_{\text{ISC}}^{S \to T,(i)} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}
$$

散逸項：

$$
\mathcal{L}_{\text{ISC}}^{S \to T,(i)}[\hat{\rho}] = k_{\text{ISC}}^{S \to T} \cdot \mathcal{D}[\hat{L}_{\text{ISC}}^{S \to T,(i)}][\hat{\rho}]
$$

#### 3.3.6 項間交差 T₁→S₀ のLindblad演算子

分子 $i$ に対して（$i = 0, 1, 2, 3$）：

$$
\hat{L}_{\text{ISC}}^{T \to S,(i)} = |0\rangle_i\langle 1|
$$

散逸項：

$$
\mathcal{L}_{\text{ISC}}^{T \to S,(i)}[\hat{\rho}] = k_{\text{ISC}}^{T \to S} \cdot \mathcal{D}[\hat{L}_{\text{ISC}}^{T \to S,(i)}][\hat{\rho}]
$$

**注意**: $\hat{L}_{\text{ISC}}^{T \to S,(i)}$ と $\hat{L}_{\text{ph}}^{(i)}$ は同じ行列形式を持つ。実効燐光速度は $\Gamma_{\text{ph,eff}} = \Gamma_{\text{ph}} + k_{\text{ISC}}^{T \to S}$ となる。

### 3.4 Lindblad演算子の総数

| 過程 | 演算子数（1ペアまたは1分子あたり） | 対象数 | 小計 |
|------|-------|--------|------|
| TTA | 2 | 3ペア | 6 |
| 蛍光 | 1 | 4分子 | 4 |
| 燐光 | 1 | 4分子 | 4 |
| 内部転換 | 1 | 4分子 | 4 |
| ISC (S→T) | 1 | 4分子 | 4 |
| ISC (T→S) | 1 | 4分子 | 4 |
| **合計** | | | **26** |

### 3.5 完全なGKSL方程式の展開

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}\left[\hat{H}_0 + \hat{H}_{\text{transfer}}, \hat{\rho}\right]
$$

$$
+ \sum_{(i,j) \in \text{neighbors}} \frac{\gamma_{\text{TTA}}}{2} \left( \mathcal{D}[\hat{L}_{\text{TTA},1}^{(ij)}][\hat{\rho}] + \mathcal{D}[\hat{L}_{\text{TTA},2}^{(ij)}][\hat{\rho}] \right)
$$

$$
+ \sum_{i=0}^{3} \Gamma_{\text{fl}} \cdot \mathcal{D}[\hat{L}_{\text{fl}}^{(i)}][\hat{\rho}]
$$

$$
+ \sum_{i=0}^{3} \Gamma_{\text{ph}} \cdot \mathcal{D}[\hat{L}_{\text{ph}}^{(i)}][\hat{\rho}]
$$

$$
+ \sum_{i=0}^{3} k_{\text{IC}} \cdot \mathcal{D}[\hat{L}_{\text{IC}}^{(i)}][\hat{\rho}]
$$

$$
+ \sum_{i=0}^{3} k_{\text{ISC}}^{S \to T} \cdot \mathcal{D}[\hat{L}_{\text{ISC}}^{S \to T,(i)}][\hat{\rho}]
$$

$$
+ \sum_{i=0}^{3} k_{\text{ISC}}^{T \to S} \cdot \mathcal{D}[\hat{L}_{\text{ISC}}^{T \to S,(i)}][\hat{\rho}]
$$

---

## 4. 物理パラメータの完全定義

### 4.1 GKSLPhysicalParametersクラスの仕様

```python
class GKSLPhysicalParameters:
    """GKSL-Lindblad量子ダイナミクス用の統一物理パラメータクラス"""

    def __init__(self):
        # === 系の構成 ===
        self.N_molecules = 4

        # === エネルギー (eV) ===
        self.E_T = 1.5      # 三重項エネルギー
        self.E_S = 3.0      # 一重項エネルギー
        # E_S0 = 0 は暗黙的（基準エネルギー）

        # === コヒーレント相互作用 (eV) ===
        self.V = 0.1         # エネルギー移動積分（Dexter機構）

        # === 散逸速度定数 ===
        # 注: 単位はすべて eV/ℏ（自然単位系での速度定数）
        # SI変換: 1 eV/ℏ = 1.52 × 10^15 s^-1

        # TTA速度定数
        self.gamma_TTA = 0.05   # TTA散逸速度定数 (eV/ℏ)

        # 放射減衰
        self.Gamma_fl = 0.01    # 蛍光放出速度 (eV/ℏ)
        self.Gamma_ph = 1e-6    # 燐光放出速度 (eV/ℏ)

        # 無放射遷移
        self.k_IC = 0.005       # 内部転換速度 (eV/ℏ)
        self.k_ISC_ST = 0.003   # 項間交差 S1→T1 速度 (eV/ℏ)
        self.k_ISC_TS = 1e-5    # 項間交差 T1→S0 速度 (eV/ℏ)

        # === 物理定数 ===
        self.hbar = 0.6582119569  # 換算プランク定数 (eV·fs)

        # === 系の構成 ===
        self.neighbors = [(0, 1), (1, 2), (2, 3)]

        # === シミュレーション条件 ===
        self.T_total = 100.0  # 総時間 (fs)
        self.N_steps = 100    # ステップ数（GKSL版は散逸精度のため増加）
        self.dt = self.T_total / self.N_steps

        # === 初期状態 ===
        self.initial_state_type = 'edge_triplet'
```

### 4.2 パラメータの物理的根拠

#### 4.2.1 $\gamma_{\text{TTA}} = 0.05$ eV/ℏ の根拠

現行ノートブックの $J = 0.05$ eV（TTAハミルトニアン結合定数）との対応：

- 現行版：TTA過程の特性時間 $\tau_{\text{TTA}} \sim \hbar / J = 0.6582 / 0.05 \approx 13.2$ fs
- GKSL版：TTA過程の特性時間 $\tau_{\text{TTA}} \sim 1 / \gamma_{\text{TTA}} = 1 / 0.05 = 20$ ℏ/eV $\approx 13.2$ fs

同じ $\gamma_{\text{TTA}} = 0.05$ eV/ℏ を使用することで、TTA過程の時間スケールを現行ノートブックと合わせる。

**重要な注意**：この対応は近似的なものである。厳密にはFermi's Golden Rule に基づき：

$$
\gamma_{\text{TTA}} = \frac{2\pi}{\hbar} |J|^2 \cdot \rho(E)
$$

ここで $\rho(E)$ は状態密度。本仕様では $\gamma_{\text{TTA}}$ を独立パラメータとして扱い、その物理的範囲内で値を設定する。

#### 4.2.2 $\Gamma_{\text{fl}} = 0.01$ eV/ℏ の根拠

現行ノートブックと同一の値を使用。

SI変換：$\Gamma_{\text{fl}} = 0.01 / 0.6582 \approx 0.0152$ fs$^{-1}$

蛍光寿命：$\tau_{\text{fl}} = 1 / (0.0152 \times 10^{15}) \approx 66$ fs

これは理想的なモデル系として妥当。実際の分子系では $\tau_{\text{fl}} \sim 1\text{-}10$ ns だが、有限のシミュレーション時間（100 fs）で効果を観測するために大きめの値を設定している。

#### 4.2.3 $\Gamma_{\text{ph}} = 10^{-6}$ eV/ℏ の根拠

燐光はスピン禁制遷移のため、蛍光より4-6桁遅い。

$\tau_{\text{ph}} = \hbar / \Gamma_{\text{ph}} = 0.6582 / 10^{-6} \approx 6.6 \times 10^{5}$ fs $\approx 0.66$ ps

100 fs のシミュレーション時間内ではほぼ影響なし。物理的に正しいが、シミュレーション結果への寄与は微小。

#### 4.2.4 その他の速度定数

- $k_{\text{IC}} = 0.005$ eV/ℏ：内部転換は蛍光と同程度のオーダー
- $k_{\text{ISC}}^{S \to T} = 0.003$ eV/ℏ：スピン-軌道結合による項間交差
- $k_{\text{ISC}}^{T \to S} = 10^{-5}$ eV/ℏ：T₁→S₀は大きなエネルギーギャップのため遅い

#### 4.2.5 実験的パラメータ範囲との対応

本仕様書のパラメータ値と、実際のTTA-UC分子系（DPA/PtOEP等）の典型的な実験値を以下に比較する：

| パラメータ | 本仕様書の値 | 典型的な実験範囲 | 備考 |
|-----------|------------|---------------|------|
| $\gamma_{\text{TTA}}$ | 0.05 eV/ℏ | $10^{-3}$–$10^{-1}$ eV/ℏ | 範囲内 |
| $\Gamma_{\text{fl}}$ | 0.01 eV/ℏ | $\sim 10^{-7}$ eV/ℏ ($\tau_{\text{fl}} \sim 1$–$10$ ns) | 100 fsシミュレーション窓のため約6桁加速 |
| $\Gamma_{\text{ph}}$ | $10^{-6}$ eV/ℏ | $\sim 10^{-10}$ eV/ℏ ($\tau_{\text{ph}} \sim \mu\text{s}$–$\text{ms}$) | 約4桁加速 |
| $k_{\text{IC}}$ | 0.005 eV/ℏ | $10^{-8}$–$10^{-4}$ eV/ℏ | 加速 |
| $k_{\text{ISC}}^{S \to T}$ | 0.003 eV/ℏ | $10^{-7}$–$10^{-4}$ eV/ℏ | 加速 |
| $V$ | 0.1 eV | 0.001–0.1 eV | 実験範囲の上端 |

**重要な注記**: 本仕様書のパラメータ値は、100 fsという短いシミュレーション時間窓内で各散逸過程の効果を観測可能にするための **教育的な（pedagogical）** 選択であり、実験データの忠実な再現を目的としていない。ただし、物理過程間の **時間スケールの階層構造** は保存されている：

$$
\gamma_{\text{TTA}} > \Gamma_{\text{fl}} > k_{\text{IC}} > k_{\text{ISC}}^{S \to T} \gg \Gamma_{\text{ph}}, k_{\text{ISC}}^{T \to S}
$$

### 4.3 パラメータ間の整合性チェック

以下の不等式が物理的に必要：

1. **TTA条件**: $2E_T \geq E_S$ → $2 \times 1.5 = 3.0 \geq 3.0$ ✓
2. **速度定数の正値性**: すべての $\gamma_\alpha > 0$ ✓
3. **時間スケール階層**（典型的）:
   $$\gamma_{\text{TTA}} \gtrsim \Gamma_{\text{fl}} > k_{\text{IC}} > k_{\text{ISC}}^{S \to T} \gg \Gamma_{\text{ph}} > k_{\text{ISC}}^{T \to S}$$
   $$0.05 > 0.01 > 0.005 > 0.003 \gg 10^{-6} > 10^{-5}$$
   
   **注意**: $\Gamma_{\text{ph}} < k_{\text{ISC}}^{T \to S}$ は通常の階層と異なるが、モデルパラメータとして許容される。

---

## 5. シナリオ1: 古典計算・ボソン無しのGKSL実装仕様

### 5.1 概要

密度行列のGKSL-Lindblad時間発展を超演算子形式で直接計算する。

### 5.2 超演算子の構築

#### 5.2.1 ベクトル化

$81 \times 81$ の密度行列 $\hat{\rho}$ を $6561$ 次元の列ベクトル $|\hat{\rho}\rangle\rangle$ に変換する。

変換規則（列優先、column-major）：

$$
\text{vec}(\hat{\rho})_{i + j \cdot d} = \hat{\rho}_{ij}, \quad d = 81
$$

NumPy実装:
```python
rho_vec = rho.flatten(order='F')  # column-major (Fortranスタイル)
```

逆変換:
```python
rho = rho_vec.reshape((d, d), order='F')
```

#### 5.2.2 ハミルトニアン超演算子

$$
\mathcal{L}_H = -\frac{i}{\hbar} \left( \hat{H} \otimes \hat{I}_d - \hat{I}_d \otimes \hat{H}^T \right)
$$

ここで $d = 81$、$\otimes$ は Kronecker 積。

NumPy実装:
```python
L_H = -1j / hbar * (np.kron(H_sys, I_d) - np.kron(I_d, H_sys.T))
```

$\mathcal{L}_H$ は $6561 \times 6561$ の複素行列。

#### 5.2.3 Lindblad超演算子

各Lindblad演算子 $\hat{L}_\alpha$（$81 \times 81$行列）と速度定数 $\gamma_\alpha$ に対して：

$$
\mathcal{L}_{\alpha} = \gamma_\alpha \left( \hat{L}_\alpha \otimes \bar{\hat{L}}_\alpha - \frac{1}{2} \hat{L}_\alpha^\dagger \hat{L}_\alpha \otimes \hat{I}_d - \frac{1}{2} \hat{I}_d \otimes (\hat{L}_\alpha^\dagger \hat{L}_\alpha)^T \right)
$$

ここで $\bar{\hat{L}}_\alpha$ は $\hat{L}_\alpha$ の複素共役（転置なし）。

NumPy実装:
```python
LdL = L_alpha.conj().T @ L_alpha
L_lindblad_alpha = gamma_alpha * (
    np.kron(L_alpha, L_alpha.conj())
    - 0.5 * np.kron(LdL, I_d)
    - 0.5 * np.kron(I_d, LdL.T)
)
```

#### 5.2.4 全超演算子

$$
\mathcal{L}_{\text{super}} = \mathcal{L}_H + \sum_{\alpha=1}^{26} \mathcal{L}_{\alpha}
$$

### 5.3 時間発展の数値計算

#### 5.3.1 方法A: 行列指数関数（小規模系向け）

$$
|\hat{\rho}(t)\rangle\rangle = e^{\mathcal{L}_{\text{super}} \cdot t} |\hat{\rho}(0)\rangle\rangle
$$

実装:
```python
import scipy.linalg

rho_0_vec = rho_0.flatten(order='F')

# 全時刻を一度に計算する場合
for step in range(N_steps):
    t = (step + 1) * dt
    rho_vec = scipy.linalg.expm(L_super * dt) @ rho_0_vec
    rho_0_vec = rho_vec  # 次のステップの初期値
```

計算量: $O(d^6) = O(81^6) \approx 10^{11}$。$d = 81$ では数秒〜数分。

#### 5.3.2 方法B: ODEソルバー（推奨）

$$
\frac{d|\hat{\rho}\rangle\rangle}{dt} = \mathcal{L}_{\text{super}} \cdot |\hat{\rho}\rangle\rangle
$$

実装:
```python
from scipy.integrate import solve_ivp

def rhs(t, rho_vec):
    return L_super @ rho_vec

rho_0_vec = rho_0.flatten(order='F')

sol = solve_ivp(
    rhs,
    [0, T_total],
    rho_0_vec,
    method='RK45',      # 非スティッフ系
    t_eval=np.linspace(0, T_total, N_steps + 1),
    rtol=1e-10,
    atol=1e-12
)
```

もし系がスティッフ（$\gamma$ の値が大きく異なる場合）なら `method='BDF'` を使用。

#### 5.3.3 方法C: 疎行列演算

超演算子は高度に疎であるため、疎行列を使用してメモリと計算時間を削減できる：

```python
from scipy.sparse import csr_matrix, kron as sparse_kron, eye as sparse_eye
from scipy.sparse.linalg import expm_multiply

L_super_sparse = csr_matrix(L_super)

# 時間発展（Krylov部分空間法）
rho_vec = expm_multiply(L_super_sparse * dt, rho_0_vec)
```

### 5.4 個体数の計算

密度行列 $\hat{\rho}$ からの個体数計算：

```python
def calculate_populations_from_density_matrix(rho, N_molecules=4, d=3):
    dim = d ** N_molecules  # = 81
    N_S0 = N_T1 = N_S1 = 0.0
    S0_per_mol = np.zeros(N_molecules)
    T1_per_mol = np.zeros(N_molecules)
    S1_per_mol = np.zeros(N_molecules)

    for mol_idx in range(N_molecules):
        # 分子 mol_idx の射影演算子を構築
        for level, proj in enumerate([np.diag([1,0,0]), np.diag([0,1,0]), np.diag([0,0,1])]):
            # 全空間への拡張
            P = np.eye(1)
            for i in range(N_molecules):
                if i == mol_idx:
                    P = np.kron(P, proj)
                else:
                    P = np.kron(P, np.eye(d))
            # 期待値 = Tr[P @ rho]
            expectation = np.real(np.trace(P @ rho))
            if level == 0:
                S0_per_mol[mol_idx] = expectation
                N_S0 += expectation
            elif level == 1:
                T1_per_mol[mol_idx] = expectation
                N_T1 += expectation
            elif level == 2:
                S1_per_mol[mol_idx] = expectation
                N_S1 += expectation

    return {
        'N_S0': N_S0, 'N_T1': N_T1, 'N_S1': N_S1,
        'S0_per_mol': S0_per_mol, 'T1_per_mol': T1_per_mol, 'S1_per_mol': S1_per_mol
    }
```

### 5.5 ClassicalGKSLSimulatorクラスの仕様

```python
class ClassicalGKSLSimulator:
    """古典的GKSL-Lindblad方程式ソルバー（ボソン無し）"""

    def __init__(self, params: GKSLPhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.d = 3
        self.dim = self.d ** self.N  # = 81

    def build_H_sys(self) -> np.ndarray:
        """系のハミルトニアン H0 + H_transfer を構築 (81×81)"""
        ...

    def build_lindblad_operators(self) -> List[Tuple[float, np.ndarray]]:
        """全Lindblad演算子と速度定数のリストを構築
        Returns: [(gamma_1, L_1), (gamma_2, L_2), ...] (26個)
        """
        ...

    def build_superoperator(self) -> np.ndarray:
        """全超演算子 L_super を構築 (6561×6561)"""
        ...

    def prepare_initial_density_matrix(self, state_type: str) -> np.ndarray:
        """初期密度行列を準備 (81×81)"""
        ...

    def calculate_populations(self, rho: np.ndarray) -> Dict:
        """密度行列から個体数を計算"""
        ...

    def simulate(self, T_total, N_steps, initial_state_type, method='RK45') -> Dict:
        """完全なGKSLシミュレーションを実行
        Returns: 現行ノートブックと同形式の出力Dict
        """
        ...
```

出力データ構造（現行ノートブックと互換）:

```python
{
    'times': List[float],
    'populations': List[Dict[str, float]],
    'per_molecule_populations': List[Dict[str, np.ndarray]],
    'rho_final': np.ndarray,     # 最終密度行列 (81×81)
    'elapsed_time': float,
    'method': 'Classical GKSL-Lindblad (no boson)',
    'entropy': List[float],       # von Neumannエントロピーの時間発展
    'trace': List[float],         # Tr[rho]の時間発展（検証用）
    'purity': List[float],        # Tr[rho^2]の時間発展
}
```

### 5.6 追加の観測量

GKSL版では以下の追加観測量を各時間ステップで計算する：

#### 5.6.1 von Neumannエントロピー

$$
S(t) = -\text{Tr}[\hat{\rho}(t) \ln \hat{\rho}(t)]
$$

実装:
```python
eigenvalues = np.linalg.eigvalsh(rho)
eigenvalues = eigenvalues[eigenvalues > 1e-15]  # ゼロ除算回避
entropy = -np.sum(eigenvalues * np.log(eigenvalues))
```

#### 5.6.2 純度（Purity）

$$
P(t) = \text{Tr}[\hat{\rho}(t)^2]
$$

- 純粋状態：$P = 1$
- 最大混合状態：$P = 1/d = 1/81$

#### 5.6.3 トレース（検証用）

$$
\text{Tr}[\hat{\rho}(t)] = 1 \pm \epsilon
$$

$\epsilon > 10^{-8}$ の場合はエラーを報告。

---

## 6. シナリオ2: 古典計算・ボソン有りのGKSL実装仕様

### 6.1 概要

電子系にフォノン（格子振動）モードを明示的に結合させたモデル。環境の効果をLindblad速度定数に繰り込むのではなく、ボソンモードの自由度を状態空間に含める。

### 6.2 状態空間の拡張

#### 6.2.1 ヒルベルト空間

$$
\mathcal{H}_{\text{total}} = \mathcal{H}_{\text{el}} \otimes \mathcal{H}_{\text{phonon}}
$$

$$
\mathcal{H}_{\text{el}} = \bigotimes_{i=0}^{3} \mathbb{C}^3 \quad (\dim = 81)
$$

$$
\mathcal{H}_{\text{phonon}} = \bigotimes_{i=0}^{3} \mathcal{H}_{\text{Fock}}^{(i)} \quad (\dim = (n_{\max}+1)^4)
$$

各分子に1つの局所フォノンモードを仮定。Fock空間を $n_{\max}$ で切断。

#### 6.2.2 次元の計算

| $n_{\max}$ | フォノン次元 | 全次元 | 密度行列要素数 |
|-----------|-----------|-------|-------------|
| 2 | $3^4 = 81$ | $81 \times 81 = 6561$ | $6561^2 \approx 4.3 \times 10^7$ |
| 3 | $4^4 = 256$ | $81 \times 256 = 20736$ | $20736^2 \approx 4.3 \times 10^8$ |
| 5 | $6^4 = 1296$ | $81 \times 1296 = 104976$ | $\sim 10^{10}$ |

$n_{\max} = 2$ を推奨パラメータとする（計算可能性と物理的妥当性のバランス）。

### 6.3 拡張ハミルトニアン

#### 6.3.1 電子系ハミルトニアン（変更なし）

$$
\hat{H}_{\text{el}} = (\hat{H}_0 + \hat{H}_{\text{transfer}}) \otimes \hat{I}_{\text{phonon}}
$$

#### 6.3.2 フォノンハミルトニアン

$$
\hat{H}_{\text{phonon}} = \hat{I}_{\text{el}} \otimes \sum_{i=0}^{3} \hbar\omega_{\text{ph}} \hat{a}_i^\dagger \hat{a}_i
$$

ここで $\omega_{\text{ph}}$ は局所フォノンの振動周波数。

フォノン消滅演算子の行列表現（$n_{\max} + 1$ 次元）：

$$
\hat{a} = \begin{pmatrix} 0 & \sqrt{1} & 0 & \cdots \\ 0 & 0 & \sqrt{2} & \cdots \\ \vdots & & & \ddots \\ 0 & \cdots & 0 & \sqrt{n_{\max}} \\ 0 & \cdots & & 0 \end{pmatrix}
$$

#### 6.3.3 電子-フォノン結合（Holstein型）

$$
\hat{H}_{e\text{-ph}} = \sum_{i=0}^{3} \sum_{n=0}^{2} g_n^{(i)} (\hat{a}_i + \hat{a}_i^\dagger) \otimes |n\rangle_i\langle n|
$$

簡略化モデル（三重項のみフォノン結合）：

$$
\hat{H}_{e\text{-ph}} = g \sum_{i=0}^{3} (\hat{a}_i + \hat{a}_i^\dagger) \otimes |1\rangle_i\langle 1|
$$

パラメータ：
- $g = 0.02$ eV（結合定数）
- $\hbar\omega_{\text{ph}} = 0.15$ eV（典型的分子内振動）
- Huang-Rhysパラメータ：$S = g^2 / (\hbar\omega_{\text{ph}})^2 = 0.0178$

#### 6.3.4 全ハミルトニアン

$$
\hat{H}_{\text{total}} = \hat{H}_{\text{el}} \otimes \hat{I}_{\text{ph}} + \hat{I}_{\text{el}} \otimes \hat{H}_{\text{phonon}} + \hat{H}_{e\text{-ph}}
$$

#### 6.3.5 Peierls型結合

Holstein型結合（§6.3.3）がオンサイトエネルギーの格子変位への依存を記述するのに対し、Peierls型結合はサイト間の移動積分の格子変位依存性を記述する：

$$
\hat{H}_{\text{Peierls}} = \sum_{\langle i,j \rangle} V_{ij} \left(1 + g_P (\hat{u}_i - \hat{u}_j)\right) \hat{T}_{ij}
$$

ここで $\hat{u}_i = (\hat{a}_i + \hat{a}_i^\dagger)/\sqrt{2}$ は無次元変位演算子、$g_P$ はPeierls結合定数、$\hat{T}_{ij}$ は分子 $i, j$ 間の三重項エネルギー移動演算子。

**注記**: 本仕様書の初期実装（Phase 1–4）ではHolstein型結合のみを含む。Peierls型結合は将来の拡張として文書化する。

#### 6.3.6 電子-光子結合

蛍光過程の微視的記述として、量子化された電磁場（光子モード）との電気双極子相互作用を考える：

$$
\hat{H}_{e\text{-photon}} = -\hat{\boldsymbol{\mu}} \cdot \hat{\mathbf{E}}
$$

ここで $\hat{\boldsymbol{\mu}}$ は遷移双極子演算子、$\hat{\mathbf{E}}$ は量子化された電場演算子。

回転波近似（RWA）の下で：

$$
\hat{H}_{e\text{-photon}}^{\text{RWA}} = \sum_i g_{\text{photon}} (|S_0\rangle_i\langle S_1| \hat{b}^\dagger + |S_1\rangle_i\langle S_0| \hat{b})
$$

ここで $\hat{b}^\dagger, \hat{b}$ は光子モードの生成・消滅演算子、$g_{\text{photon}}$ は電子-光子結合定数。

**注記**: 本仕様書の初期実装では、蛍光をLindblad散逸項として現象論的に扱い、光子モードの自由度は明示的に含めない。電子-光子結合の量子回路実装は将来の拡張とする。

#### 6.3.7 完全ハミルトニアン（5項）

フォノン・光子モードをすべて含む場合の完全ハミルトニアンは5項から成る：

$$
\hat{H}_{\text{total}} = \hat{H}_{\text{el}} + \hat{H}_{\text{phonon}} + \hat{H}_{\text{photon}} + \hat{H}_{e\text{-ph}} + \hat{H}_{e\text{-photon}}
$$

**注記**: 本仕様書の初期実装ではHolstein型電子-フォノン結合のみを含む（$\hat{H}_{\text{el}} + \hat{H}_{\text{phonon}} + \hat{H}_{e\text{-ph}}$ の3項）。Peierls結合および電子-光子結合は将来の拡張として文書化する。

#### 6.3.8 スペクトル密度関数

ボソン浴の連続的な記述にはスペクトル密度関数が用いられる：

$$
J(\omega) = \sum_k |g_k|^2 \delta(\omega - \omega_k)
$$

連続極限で用いられる代表的なモデル：

- **Drude-Lorentz型**:
$$
J(\omega) = \frac{2\lambda \gamma_c \omega}{\omega^2 + \gamma_c^2}
$$
ここで $\lambda$ は再配列エネルギー、$\gamma_c$ はカットオフ周波数。

- **Ohmic型**:
$$
J(\omega) = \eta \omega e^{-\omega/\omega_c}
$$
ここで $\eta$ は無次元結合強度、$\omega_c$ はカットオフ周波数。

**注記**: 本仕様書では離散フォノンモード（§6.3.2–6.3.3）を使用する。連続スペクトル密度関数によるボソン浴の扱いは、HEOM等の非マルコフ手法（§6.8.1）と組み合わせた将来の拡張として文書化する。

#### 6.3.9 有限温度効果

有限温度 $T$ でのボソン浴は、Bose-Einstein分布に従う熱的占有を持つ：

$$
\bar{n}(\omega) = \frac{1}{e^{\hbar\omega/k_BT} - 1}
$$

これにより、散逸過程は詳細平衡条件を満たす：

$$
\frac{\gamma_{\text{absorption}}}{\gamma_{\text{emission}}} = \frac{\bar{n}(\omega)}{\bar{n}(\omega) + 1} = e^{-\hbar\omega/k_BT}
$$

高温極限（$k_BT \gg \hbar\omega$）では $\bar{n}(\omega) \approx k_BT/\hbar\omega \gg 1$ となり、吸収と放出がほぼ等しくなる。

**注記**: 本仕様書の初期実装では $T = 0$ を仮定する（フォノンの初期状態が真空状態 $|0\rangle$）。有限温度効果の実装は将来の拡張とする。

#### 6.3.10 実装優先度に関する注記

上記 §6.3.5–6.3.9 の各項目の実装優先度を以下にまとめる：

| 項目 | 初期実装 | 将来の拡張 |
|------|---------|-----------|
| Holstein型電子-フォノン結合（§6.3.3） | ✅ | — |
| Peierls型結合（§6.3.5） | — | 将来 |
| 電子-光子結合（§6.3.6） | — | 将来 |
| 完全ハミルトニアン5項（§6.3.7） | 3項のみ | 将来 |
| スペクトル密度関数（§6.3.8） | — | 将来（HEOM等と組合せ） |
| 有限温度効果（§6.3.9） | — | 将来 |

これらは理論的完備性のために文書化するものであり、初期実装ではHolstein型フォノン結合のみを明示的に扱う。

### 6.4 Lindblad演算子の拡張

ボソン有りモデルでも散逸項は電子系にのみ作用する。ただし演算子は拡張空間上に定義する：

$$
\hat{L}_\alpha^{\text{extended}} = \hat{L}_\alpha^{\text{el}} \otimes \hat{I}_{\text{phonon}}
$$

すなわちLindblad演算子の行列形式は変わらないが、サイズが $81 \times 81$ から $(81 \times (n_{\max}+1)^4) \times (81 \times (n_{\max}+1)^4)$ に拡張される。

### 6.5 追加パラメータ

```python
class GKSLPhysicalParametersWithBoson(GKSLPhysicalParameters):
    def __init__(self):
        super().__init__()
        # フォノンパラメータ
        self.n_max = 2                    # フォノンFock空間の切断
        self.omega_ph = 0.15              # フォノン振動周波数 (eV/ℏ)
        self.g_eph = 0.02                 # 電子-フォノン結合定数 (eV)
        self.huang_rhys = self.g_eph**2 / (self.hbar * self.omega_ph)**2
```

### 6.6 数値手法

$n_{\max} = 2$ の場合、全次元 $d_{\text{total}} = 81 \times 81 = 6561$。

超演算子サイズ: $6561^2 \approx 4.3 \times 10^7$。

行列指数関数は困難。**ODEソルバー**（疎行列版）を使用：

```python
from scipy.sparse.linalg import LinearOperator

def L_super_matvec(rho_vec):
    """超演算子の行列-ベクトル積"""
    rho = rho_vec.reshape((d_total, d_total), order='F')
    drho = -1j / hbar * (H_total @ rho - rho @ H_total)
    for gamma, L in lindblad_ops:
        LdL = L.conj().T @ L
        drho += gamma * (L @ rho @ L.conj().T - 0.5 * LdL @ rho - 0.5 * rho @ LdL)
    return drho.flatten(order='F')

L_op = LinearOperator((d_total**2, d_total**2), matvec=L_super_matvec)
```

### 6.7 ClassicalGKSLBosonSimulatorクラスの仕様

```python
class ClassicalGKSLBosonSimulator:
    """古典的GKSL-Lindblad方程式ソルバー（ボソン有り）"""

    def __init__(self, params: GKSLPhysicalParametersWithBoson):
        ...

    def build_H_total(self) -> np.ndarray:
        """拡張ハミルトニアンを構築"""
        ...

    def build_extended_lindblad_operators(self) -> List[Tuple[float, np.ndarray]]:
        """拡張空間上のLindblad演算子を構築"""
        ...

    def simulate(self, T_total, N_steps, initial_state_type, method='BDF') -> Dict:
        """GKSLシミュレーション（ボソン有り）"""
        ...
```

出力には電子系の縮約密度行列の個体数を含める：

$$
\hat{\rho}_{\text{el}}(t) = \text{Tr}_{\text{phonon}}[\hat{\rho}_{\text{total}}(t)]
$$

### 6.8 追加の数値手法（ボソン有りモデル用）

ボソン有りモデルの状態空間は急速に増大するため、直接的なODE積分以外の数値手法も検討する。以下の手法は理論的完備性のために文書化し、初期実装では直接ODE積分を使用する。

#### 6.8.1 HEOM（階層的運動方程式）

非マルコフ的ボソン浴を扱うための体系的手法。Drude-Lorentz型スペクトル密度に対して：

$$
\frac{\partial \hat{\rho}_\mathbf{n}}{\partial t} = -\left(\frac{i}{\hbar}\hat{H}_S^\times + \sum_k n_k \gamma_k\right)\hat{\rho}_\mathbf{n} + \sum_k \hat{V}_k^\times \hat{\rho}_{\mathbf{n}+\mathbf{e}_k} + \sum_k n_k \hat{C}_k \hat{\rho}_{\mathbf{n}-\mathbf{e}_k}
$$

ここで $\hat{H}_S^\times = [\hat{H}_S, \cdot]$ は交換子超演算子、$\mathbf{n} = (n_1, n_2, \ldots, n_K)$ は階層インデックス、$\mathbf{e}_k$ は $k$ 番目の単位ベクトル。

計算コスト: $O\left(\binom{L+K}{K} \cdot d^4\right)$（$L$: 階層切断レベル、$K$: 相関関数の指数項数、$d$: 系の次元）

#### 6.8.2 テンソルネットワーク法

大規模ボソン空間に対するMPS/MPO（行列積状態/行列積演算子）表現：

- ボンド次元 $\chi$ が精度を制御
- TEBD（Time-Evolving Block Decimation）アルゴリズムによる時間発展
- 計算コスト: $O(N \cdot d \cdot \chi^3)$

量子回路シミュレーションとの親和性が高く、量子-古典ハイブリッド手法への拡張が可能。

#### 6.8.3 量子モンテカルロ法（確率的波動関数法）

マスター方程式の確率的アンラベリング（stochastic unraveling）：

- 量子ジャンプ軌道（quantum jump trajectories）の集団平均
- $N_{\text{traj}}$ 本の軌道を生成し、統計的に密度行列を再構成
- 計算コスト: $O(N_{\text{traj}} \cdot d \cdot t/\tau)$（$\tau$: ジャンプの典型的時間間隔）

波動関数ベースのため、密度行列を直接扱うよりメモリ効率が良い。

**注記**: 初期実装では直接ODE積分（§6.6）を使用する。HEOM、テンソルネットワーク、量子モンテカルロの各手法は、状態空間が大きい場合（$N > 4$ 分子、$n_{\max} > 3$）の将来の実装候補として文書化する。

---

## 7. シナリオ3: Qubit量子計算・ボソン無しのGKSL実装仕様

### 7.1 概要

現行ノートブックのQubit実装を拡張し、Stinespring dilationを用いてLindblad散逸項を量子回路に実装する。

### 7.2 Qubitエンコーディング（現行と同一）

各分子 $i$ を2 qubit $(q_{2i}, q_{2i+1})$ で表現：

$$
|S_0\rangle_i \leftrightarrow |00\rangle, \quad |T_1\rangle_i \leftrightarrow |01\rangle, \quad |S_1\rangle_i \leftrightarrow |10\rangle
$$

禁止状態：$|11\rangle$

系qubit数：$2N = 8$

#### 7.2.1 物理的部分空間と射影演算子

各分子 $i$ の物理的部分空間への射影演算子（$|11\rangle$ 状態を除外）：

$$
\hat{P}_{\text{phys}}^{(i)} = |00\rangle\langle 00| + |01\rangle\langle 01| + |10\rangle\langle 10|
$$

全系の射影演算子：

$$
\hat{P}_{\text{phys}} = \bigotimes_{i=0}^{N-1} \hat{P}_{\text{phys}}^{(i)}
$$

Stinespring dilationの実装において、禁止状態 $|11\rangle$ への遷移確率を各ステップで検証する：

$$
P_{\text{forbidden}}(t) = 1 - \text{Tr}[\hat{P}_{\text{phys}} \hat{\rho}_{\text{sys}}(t)] < 10^{-8}
$$

この検証は、Qubitエンコーディングにおけるゲート分解が物理的部分空間を保存していることの確認に不可欠である。Qutritエンコーディング（§9）では禁止状態が存在しないため、この検証は不要。

### 7.3 Stinespring Dilationの原理

#### 7.3.1 基本定理

任意のCPTP写像 $\mathcal{E}[\hat{\rho}]$ は、補助系（ancilla）を追加したユニタリ演算と部分トレースで実現できる：

$$
\mathcal{E}[\hat{\rho}_S] = \text{Tr}_E\left[\hat{U}_{SE}(\hat{\rho}_S \otimes |0\rangle_E\langle 0|)\hat{U}_{SE}^\dagger\right]
$$

#### 7.3.2 Lindblad演算子からのStinespringユニタリの構築

単一Lindblad演算子 $\hat{L}$ に対する微小時間 $\Delta t$ の散逸ステップ：

$$
\mathcal{E}_{\Delta t}[\hat{\rho}] \approx \hat{\rho} + \gamma \Delta t \cdot \mathcal{D}[\hat{L}][\hat{\rho}]
$$

これを実現するStinespringユニタリの構築方法：

1. **生成子の定義**:

$$
\hat{G} = \hat{L} \otimes |1\rangle_E\langle 0| + \hat{L}^\dagger \otimes |0\rangle_E\langle 1|
$$

2. **ユニタリ演算子**:

$$
\hat{U}(\theta) = e^{-i\theta \hat{G}}
$$

ここで $\theta = \sqrt{\gamma \Delta t}$。

3. **Stinespring近似の精度**:

$$
\text{Tr}_E[\hat{U}(\theta)(\hat{\rho} \otimes |0\rangle\langle 0|)\hat{U}^\dagger(\theta)] = \hat{\rho} + \gamma \Delta t \cdot \mathcal{D}[\hat{L}][\hat{\rho}] + O(\gamma^2 \Delta t^2)
$$

誤差は $O(\gamma^2 \Delta t^2)$ で、Trotterステップ数を増やすことで改善可能。

### 7.4 各Lindblad演算子のQubit量子回路

#### 7.4.1 蛍光 $\hat{L}_{\text{fl}}^{(i)} = |0\rangle_i\langle 2|$

**Qubit表現**: $|S_1\rangle = |10\rangle \to |S_0\rangle = |00\rangle$、つまり $|10\rangle\langle 10|$ を制御として ancilla qubit に回転を適用後、系のqubitを条件付きで反転。

**回路構成**:
- 系 qubit: $q_0 = 2i$, $q_1 = 2i+1$
- Ancilla qubit: $q_E$（各Lindblad演算子に1つ）

手順:
1. $q_1 = 1$ かつ $q_0 = 0$（$|10\rangle$ 状態）を検出
2. 条件付きで ancilla $q_E$ に $R_Y(2\arcsin(\sqrt{\Gamma_{\text{fl}} \Delta t}))$ を適用
3. Ancilla を $|0\rangle$ に射影（部分トレース → ancilla をリセットまたは破棄）
4. Ancilla が $|1\rangle$ の場合: 系を $|10\rangle \to |00\rangle$ に遷移（$q_1$ を反転）

**具体的ゲート列**:
```
1. X gate on q0 (|10⟩ → |11⟩ を制御条件にするため)
2. Toffoli(q0, q1, q_E) with RY rotation
   = CCX制御の条件付き回転
3. 条件: q_E = |1⟩ ならば X gate on q1 (|10⟩ → |00⟩)
4. X gate on q0 (元に戻す)
```

**回路の厳密な構成**（制御回転の分解）:

$$
\hat{U}_{\text{fl}}^{(i)} = e^{-i\sqrt{\Gamma_{\text{fl}} \Delta t} \cdot \hat{G}_{\text{fl}}}
$$

$$
\hat{G}_{\text{fl}} = |00\rangle\langle 10|_{(q_0,q_1)} \otimes |1\rangle\langle 0|_E + |10\rangle\langle 00|_{(q_0,q_1)} \otimes |0\rangle\langle 1|_E
$$

#### 7.4.2 燐光 $\hat{L}_{\text{ph}}^{(i)} = |0\rangle_i\langle 1|$

**Qubit表現**: $|T_1\rangle = |01\rangle \to |S_0\rangle = |00\rangle$

$$
\hat{G}_{\text{ph}} = |00\rangle\langle 01|_{(q_0,q_1)} \otimes |1\rangle\langle 0|_E + |01\rangle\langle 00|_{(q_0,q_1)} \otimes |0\rangle\langle 1|_E
$$

回路：$q_0 = 1$ を制御として ancilla に回転、条件付きで $q_0$ を反転。

#### 7.4.3 ISC S₁→T₁ $\hat{L}_{\text{ISC}}^{S \to T,(i)} = |1\rangle_i\langle 2|$

**Qubit表現**: $|S_1\rangle = |10\rangle \to |T_1\rangle = |01\rangle$

$$
\hat{G}_{\text{ISC}} = |01\rangle\langle 10|_{(q_0,q_1)} \otimes |1\rangle\langle 0|_E + |10\rangle\langle 01|_{(q_0,q_1)} \otimes |0\rangle\langle 1|_E
$$

回路：$q_1 = 1, q_0 = 0$ を検出し、条件付きで $q_0, q_1$ を同時に反転（$|10\rangle \to |01\rangle$）。

#### 7.4.4 TTA $\hat{L}_{\text{TTA},1}^{(ij)} = |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|$

**Qubit表現**: $|T_1\rangle_i|T_1\rangle_j = |01\rangle_i|01\rangle_j \to |S_1\rangle_i|S_0\rangle_j = |10\rangle_i|00\rangle_j$

4-qubit系 + 1 ancilla の操作。

$$
\hat{G}_{\text{TTA},1} = |10,00\rangle\langle 01,01| \otimes |1\rangle\langle 0|_E + |01,01\rangle\langle 10,00| \otimes |0\rangle\langle 1|_E
$$

回路：4 qubit ($q_{2i}, q_{2i+1}, q_{2j}, q_{2j+1}$) の $|01,01\rangle$ 状態を検出し、条件付きで ancilla に回転を適用後、$|01,01\rangle \to |10,00\rangle$ に遷移。

### 7.5 Trotter分解の構成

1ステップの時間発展：

$$
e^{\mathcal{L}_{\text{GKSL}} \Delta t} \approx e^{\mathcal{L}_H \Delta t/2} \cdot e^{\mathcal{L}_{\text{diss}} \Delta t} \cdot e^{\mathcal{L}_H \Delta t/2}
$$

量子回路としての実装順序：

```
1. ユニタリ前半（Δt/2）:
   - H0 evolution: 分子 0,1,2,3
   - H_transfer evolution: ペア (0,1), (1,2), (2,3)

2. 散逸ステップ（Δt）:
   - TTA Lindblad: 各ペアに対して2つの演算子
   - 蛍光 Lindblad: 各分子
   - 燐光 Lindblad: 各分子
   - IC Lindblad: 各分子
   - ISC S→T Lindblad: 各分子
   - ISC T→S Lindblad: 各分子
   各ステップで ancilla qubit を |0⟩ にリセット

3. ユニタリ後半（Δt/2）:
   - H_transfer evolution: ペア (2,3), (1,2), (0,1) [逆順]
   - H0 evolution: 分子 3,2,1,0 [逆順]
```

### 7.6 ゲート分解戦略

#### 7.6.1 2段階アプローチ（現行ノートブック準拠）

現行ノートブックと同様に、GKSLの量子回路も2段階で比較する：

**レベル1: UnitaryGate表現**
- Stinespring ユニタリ $\hat{U}_{\text{Lindblad}}$ を直接 `UnitaryGate`（Qiskit）として回路に適用
- ゲート数: 少ない（高レベル表現）
- 検証用途に適している

**レベル2: 基本ゲート分解**
- UnitaryGate → KAK分解（Cartan分解）→ CNOT + Rz + Ry + Rx
- ゲート数: 多い（実機で実行可能な低レベル表現）
- 実機実行用途

#### 7.6.2 ゲート数計測

`comparison_helpers` モジュールの以下の関数を使用：
- `count_gates_by_type(circuit)`: ゲート種類別のカウント
- `decompose_qiskit_unitary_gates(circuit)`: Qubit UnitaryGate分解

#### 7.6.3 Stinespring ユニタリの分解特性

GKSL回路固有の構造:
- 各Stinespring ユニタリは **疎行列**（$|\cdot\rangle\langle \cdot| \otimes |\cdot\rangle\langle \cdot|$ の形式の項のみ非ゼロ）
- 単一分子Lindblad演算子の場合、Stinespring ユニタリは2×2の非自明ブロックのみを持つ
- TTA Stinespring（4 qubit系 + 1 ancilla）はより多くのゲートが必要

### 7.7 必要な量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 系 qubit | 8 | 4分子 × 2 qubit |
| Ancilla qubit (TTA) | 6 | 3ペア × 2チャネル |
| Ancilla qubit (蛍光) | 4 | 4分子 |
| Ancilla qubit (燐光) | 4 | 4分子 |
| Ancilla qubit (IC) | 4 | 4分子 |
| Ancilla qubit (ISC S→T) | 4 | 4分子 |
| Ancilla qubit (ISC T→S) | 4 | 4分子 |
| **合計** | **34** | |

ancilla の再利用を行えば削減可能（ミッドサーキット測定によるリセット）。

### 7.8 QubitGKSLSimulatorクラスの仕様

```python
class QubitGKSLSimulator:
    """Qubit GKSL-Lindblad量子シミュレータ（ボソン無し）"""

    def __init__(self, params: GKSLPhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.n_sys_qubits = 2 * self.N     # 8
        self.n_ancilla = 26                 # Lindblad演算子数
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla  # 34

    def build_unitary_step(self, circuit, dt):
        """ユニタリ部分（H0 + H_transfer）の半ステップ回路"""
        ...

    def build_lindblad_step(self, circuit, dt):
        """全Lindblad散逸ステップの回路"""
        ...

    def simulate_statevector(self, T_total, N_steps, initial_state_type) -> Dict:
        """Statevectorシミュレーション（密度行列再構成用）
        ancillaを含む全系の状態ベクトルから部分トレースで系の密度行列を取得。
        ショットノイズなし。GKSL古典ソルバーとの一致検証に使用。
        """
        ...

    def simulate_shot_based(self, T_total, N_steps, initial_state_type, shots) -> Dict:
        """ショットベースシミュレーション
        実際の量子ハードウェアを模擬。ショット数に依存した統計的揺らぎを含む。
        """
        ...

    def simulate(self, T_total, N_steps, initial_state_type, shots=None) -> Dict:
        """完全なGKSLシミュレーション
        shots=None の場合は statevector、指定時は shot-based を使用。
        """
        ...
```

### 7.9 密度行列の再構成（Statevector方式）

Statevectorシミュレータを使用する場合、ancillaを含む全系の状態ベクトルから系のみの密度行列を部分トレースで取得：

$$
\hat{\rho}_{\text{sys}} = \text{Tr}_{\text{ancilla}}[|\Psi_{\text{total}}\rangle\langle\Psi_{\text{total}}|]
$$

```python
from qiskit.quantum_info import Statevector, partial_trace

sv = Statevector(circuit)
# ancilla qubitのインデックスリスト
ancilla_indices = list(range(self.n_sys_qubits, self.n_total_qubits))
rho_sys = partial_trace(sv, ancilla_indices)
```

### 7.10 ハードウェアノイズモデル（Qubit GKSL）

**重要な区別**: Lindblad散逸（物理的TTA/蛍光過程）とハードウェアノイズ（ゲート不完全性）は概念的に完全に独立である。GKSL実装ではこの両方を含むシミュレーションが可能。

#### 7.10.1 脱分極エラー

2-qubitゲートのみに適用（現行ノートブックと同一仕様）：
- 1-qubitゲート: 理想的（ノイズなし）
- 2-qubitゲート: $p_{\text{depol}} = 0.01$ (1.0%)

$$
\mathcal{E}_{\text{depol}}[\hat{\rho}] = (1 - p)\hat{\rho} + \frac{p}{d^2 - 1}\sum_{P \neq I} P\hat{\rho}P^\dagger
$$

ここで $d$ は演算子が作用する部分空間の次元（2-qubitゲートの場合 $d = 4$）。

#### 7.10.2 熱緩和

2-qubitゲートにのみ適用：
- $T_1 = 50\,\mu\text{s} = 5 \times 10^{10}\,\text{fs}$
- $T_2 = 70\,\mu\text{s} = 7 \times 10^{10}\,\text{fs}$
- 2-qubitゲート時間: $300\,\text{fs}$

#### 7.10.3 QubitGKSLNoisySimulatorクラスの仕様

```python
class QubitGKSLNoisySimulator(QubitGKSLSimulator):
    """ハードウェアノイズ付きQubit GKSLシミュレータ"""

    def __init__(self, params: GKSLPhysicalParameters, noise_model):
        super().__init__(params)
        self.noise_model = noise_model

    def simulate(self, T_total, N_steps, initial_state_type, shots) -> Dict:
        """ハードウェアノイズ付きGKSLシミュレーション"""
        ...
```

---

## 8. シナリオ4: Qubit量子計算・ボソン有りのGKSL実装仕様

### 8.1 概要

シナリオ3にフォノンモードのqubit表現を追加したモデル。

### 8.2 フォノンモードのQubitエンコーディング

#### 8.2.1 バイナリエンコーディング

$n_{\max}$ までの Fock 状態を $\lceil\log_2(n_{\max}+1)\rceil$ qubit で表現：

$$
|n\rangle_{\text{Fock}} \leftrightarrow |n_{\text{binary}}\rangle_{\text{qubit}}
$$

$n_{\max} = 2$ の場合：$\lceil\log_2 3\rceil = 2$ qubit/フォノンモード

$$
|0\rangle_{\text{Fock}} \leftrightarrow |00\rangle, \quad |1\rangle_{\text{Fock}} \leftrightarrow |01\rangle, \quad |2\rangle_{\text{Fock}} \leftrightarrow |10\rangle
$$

禁止状態：$|11\rangle$

### 8.3 必要な量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 電子系 qubit | 8 | 4分子 × 2 qubit |
| フォノン qubit | 8 | 4分子 × 2 qubit ($n_{\max} = 2$) |
| Ancilla qubit | 26 | Lindblad演算子数 |
| **合計** | **42** | |

### 8.4 電子-フォノン結合の量子回路

Holstein型結合 $g(\hat{a} + \hat{a}^\dagger)|1\rangle\langle 1|$ の回路実装：

1. 電子系の $|T_1\rangle = |01\rangle$ を制御条件として検出
2. フォノン qubit に対して条件付きインクリメント/デクリメント回路を適用

$$
e^{-ig\Delta t(\hat{a}+\hat{a}^\dagger)|1\rangle\langle 1|/\hbar} \approx \text{C-}[R_X(2g\Delta t/\hbar)]
$$

ただし、バイナリエンコーディングでのフォノン昇降演算子の実装は非自明であり、$O(n_{\max})$ 個の制御ゲートが必要。

### 8.5 QubitGKSLBosonSimulatorクラスの仕様

```python
class QubitGKSLBosonSimulator:
    """Qubit GKSL-Lindblad量子シミュレータ（ボソン有り）"""

    def __init__(self, params: GKSLPhysicalParametersWithBoson):
        self.n_phonon_qubits_per_mol = int(np.ceil(np.log2(params.n_max + 1)))
        self.n_sys_qubits = 2 * params.N_molecules + \
                            self.n_phonon_qubits_per_mol * params.N_molecules
        self.n_ancilla = 26
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla

    def build_phonon_operators(self):
        """フォノン昇降演算子のqubit回路を構築"""
        ...

    def build_eph_coupling_circuit(self, circuit, mol_idx, dt):
        """電子-フォノン結合ゲートを構築"""
        ...

    def simulate(self, T_total, N_steps, initial_state_type, shots) -> Dict:
        ...
```

---

## 9. シナリオ5: Qudit量子計算・ボソン無しのGKSL実装仕様

### 9.1 概要

MQT-Quditsフレームワークを使用し、各分子を1 qutrit（$d = 3$）で表現。Stinespring dilationで散逸項を実装。

### 9.2 Qutritエンコーディング（現行と同一）

$$
|S_0\rangle_i \leftrightarrow |0\rangle_i, \quad |T_1\rangle_i \leftrightarrow |1\rangle_i, \quad |S_1\rangle_i \leftrightarrow |2\rangle_i
$$

系 qutrit 数：$N = 4$  
状態空間次元：$3^4 = 81$（全て物理状態、禁止状態なし）

#### 9.2.1 MQT-Quditsの基本ゲートセット

MQT-Quditsフレームワークで使用される基本ゲートの数学的定義を以下にまとめる。

**一般化Pauli-X演算子**:

$$
\hat{X}_d = \sum_{j=0}^{d-1} |j+1 \bmod d\rangle\langle j|
$$

**一般化Pauli-Z演算子**:

$$
\hat{Z}_d = \sum_{j=0}^{d-1} \omega^j |j\rangle\langle j|, \quad \omega = e^{2\pi i/d}
$$

**部分空間回転**:

$$
\hat{R}_{mn}(\theta, \phi) = e^{-i\theta(\cos\phi\,\hat{\sigma}_{mn}^x + \sin\phi\,\hat{\sigma}_{mn}^y)/2}
$$

ここで $\hat{\sigma}_{mn}^x = |m\rangle\langle n| + |n\rangle\langle m|$, $\hat{\sigma}_{mn}^y = -i|m\rangle\langle n| + i|n\rangle\langle m|$ は部分空間 $\{|m\rangle, |n\rangle\}$ 上の一般化Pauli行列。

**Givens回転**:

$$
G_{mn}(\theta, \phi) = I + (\cos\theta - 1)(|m\rangle\langle m| + |n\rangle\langle n|) + \sin\theta(e^{i\phi}|m\rangle\langle n| - e^{-i\phi}|n\rangle\langle m|)
$$

**MQT-Quditsの具体的ゲート**:

| ゲート名 | 記号 | 数学的定義 | 備考 |
|---------|------|----------|------|
| 仮想Z回転 | `VirtRz` | 対角位相ゲート $\text{diag}(e^{i\phi_0}, e^{i\phi_1}, \ldots)$ | 古典的フレーム変更で実装 |
| 部分空間回転 | `R` | $\hat{R}_{mn}(\theta, \phi)$ | レベル $m, n$ 間の回転 |
| 部分空間Hadamard-like | `Rh` | $\hat{R}_{mn}(\pi/2, 0)$ | Hadamardの一般化 |
| Z回転 | `Rz` | $\hat{R}_{mn}(0, \phi)$ 相当の位相操作 | 位相調整 |
| 制御交換 | `CEx` | 制御付き部分空間交換 | 2-qudit基本ゲート |

Lindblad演算子のStinespringユニタリ（§9.4）は、IntegratedSparseCompilerV2によりこれらの基本ゲートに分解される。

### 9.3 Stinespring Dilationの Qudit 実装

#### 9.3.1 補助系の選択

各Lindblad演算子に対して、1つの補助 qubit（$d = 2$）または 補助 qutrit（$d = 3$）を使用。

**補助 qubit（$d = 2$）の場合**：

$$
\hat{U}_{\text{Lindblad}} = e^{-i\theta \hat{G}}
$$

$$
\hat{G} = \hat{L}_{\text{sys}} \otimes \hat{\sigma}_E^- + \hat{L}_{\text{sys}}^\dagger \otimes \hat{\sigma}_E^+
$$

ここで $\hat{\sigma}_E^- = |0\rangle\langle 1|_E$, $\hat{\sigma}_E^+ = |1\rangle\langle 0|_E$、$\theta = \sqrt{\gamma \Delta t}$。

### 9.4 各Lindblad演算子のQudit量子回路

#### 9.4.1 蛍光 $\hat{L}_{\text{fl}}^{(i)} = |0\rangle_i\langle 2|$

系のqutrit $i$ と ancilla qubit $E$ に作用するStinespring ユニタリ：

$$
\hat{U}_{\text{fl}} = \exp\left(-i\sqrt{\Gamma_{\text{fl}} \Delta t} \cdot \hat{G}_{\text{fl}}\right)
$$

$$
\hat{G}_{\text{fl}} = |0\rangle_i\langle 2| \otimes |1\rangle_E\langle 0| + |2\rangle_i\langle 0| \otimes |0\rangle_E\langle 1|
$$

この演算子は $6 \times 6$ ユニタリ（系 qutrit $3$ + ancilla qubit $2$）のブロック構造を持つ。

基底 $\{|0,0\rangle, |0,1\rangle, |1,0\rangle, |1,1\rangle, |2,0\rangle, |2,1\rangle\}$ での行列表現：

$$
\hat{G}_{\text{fl}} = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

（$|0,0\rangle\langle 2,1|$ と $|2,1\rangle\langle 0,0|$ のみ非ゼロ）

$$
\hat{U}_{\text{fl}}(\theta) = \hat{I} + (\cos\theta - 1)(|0,0\rangle\langle 0,0| + |2,1\rangle\langle 2,1|) - i\sin\theta(|0,0\rangle\langle 2,1| + |2,1\rangle\langle 0,0|)
$$

ここで $\theta = \sqrt{\Gamma_{\text{fl}} \Delta t}$。

MQT-Qudits回路での実装：`CustomTwo` ゲートとして $6 \times 6$ ユニタリを適用するか、または部分空間回転 $R_{02}$ の制御付きバージョンとして分解。

#### 9.4.2 燐光 $\hat{L}_{\text{ph}}^{(i)} = |0\rangle_i\langle 1|$

$$
\hat{G}_{\text{ph}} = |0\rangle_i\langle 1| \otimes |1\rangle_E\langle 0| + |1\rangle_i\langle 0| \otimes |0\rangle_E\langle 1|
$$

$6 \times 6$ 行列表現（同じ基底）：

$$
\hat{G}_{\text{ph}} = \begin{pmatrix}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

$$
\hat{U}_{\text{ph}}(\theta) = \hat{I} + (\cos\theta - 1)(|0,0\rangle\langle 0,0| + |1,1\rangle\langle 1,1|) - i\sin\theta(|0,0\rangle\langle 1,1| + |1,1\rangle\langle 0,0|)
$$

$\theta = \sqrt{\Gamma_{\text{ph}} \Delta t}$

#### 9.4.3 ISC S₁→T₁ $\hat{L}_{\text{ISC}}^{S \to T,(i)} = |1\rangle_i\langle 2|$

$$
\hat{G}_{\text{ISC}} = |1\rangle_i\langle 2| \otimes |1\rangle_E\langle 0| + |2\rangle_i\langle 1| \otimes |0\rangle_E\langle 1|
$$

$6 \times 6$ 行列表現：

$$
\hat{G}_{\text{ISC}} = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0
\end{pmatrix}
$$

$$
\hat{U}_{\text{ISC}}(\theta) = \hat{I} + (\cos\theta - 1)(|1,0\rangle\langle 1,0| + |2,1\rangle\langle 2,1|) - i\sin\theta(|1,0\rangle\langle 2,1| + |2,1\rangle\langle 1,0|)
$$

$\theta = \sqrt{k_{\text{ISC}}^{S \to T} \Delta t}$

#### 9.4.4 TTA $\hat{L}_{\text{TTA},1}^{(ij)} = |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|$

2 qutrit + 1 ancilla qubit のStinespring ユニタリ：

$$
\hat{G}_{\text{TTA},1} = |20\rangle_{ij}\langle 11| \otimes |1\rangle_E\langle 0| + |11\rangle_{ij}\langle 20| \otimes |0\rangle_E\langle 1|
$$

これは $18 \times 18$ ユニタリ（$3 \times 3 \times 2 = 18$ 次元）。

$$
\hat{U}_{\text{TTA},1}(\theta) = \hat{I}_{18} + (\cos\theta - 1)(|20,0\rangle\langle 20,0| + |11,1\rangle\langle 11,1|) - i\sin\theta(|20,0\rangle\langle 11,1| + |11,1\rangle\langle 20,0|)
$$

$\theta = \sqrt{(\gamma_{\text{TTA}}/2) \Delta t}$

MQT-Qudits回路：`CustomTwo` ゲート（$9 \times 9$、2 qutrit部分空間）に ancilla qubit を追加した拡張ゲートとして実装。

### 9.5 Trotter分解の構成

```
1. ユニタリ前半（Δt/2）:
   - H0: VirtRz ゲート（各 qutrit）
   - H_transfer: CEx ゲート（各ペア）

2. 散逸ステップ（Δt）:
   - TTA Stinespring: CustomTwo + ancilla（各ペア × 2）
   - 蛍光 Stinespring: qutrit + ancilla（各分子）
   - 燐光 Stinespring: qutrit + ancilla（各分子）
   - IC Stinespring: qutrit + ancilla（各分子）
   - ISC S→T Stinespring: qutrit + ancilla（各分子）
   - ISC T→S Stinespring: qutrit + ancilla（各分子）

3. ユニタリ後半（Δt/2）:
   - H_transfer: CEx ゲート（各ペア、逆順）
   - H0: VirtRz ゲート（各 qutrit、逆順）
```

#### 9.5.1 ゲート分解戦略（Qudit版）

Qubit版（§7.6）と同様の2段階アプローチ：

**レベル1: CustomTwo表現**
- Stinespring ユニタリを直接 `CustomTwo`（MQT-Qudits）ゲートとして回路に適用
- ゲート数: 少ない（高レベル表現）

**レベル2: 基本ゲート分解**
- CustomTwo → IntegratedSparseCompilerV2 → VirtRz + R + Rh + Rz + CEx
- ゲート数: 多い（実機で実行可能な低レベル表現）

ゲート数計測には `comparison_helpers` モジュールの以下の関数を使用：
- `decompose_qudit_customtwo_gates_to_circuit(circuit)`: Qudit CustomTwo分解

GKSL回路固有の特性:
- 各Stinespring ユニタリは **疎行列** であり、MQT-Quditsの疎構造認識コンパイラで効率的に分解可能
- 典型的には **~6ゲート/Stinespring ユニタリ**（単一分子Lindblad演算子の場合）
- TTA Stinespring（2 qutrit + 1 ancilla）はより多くのゲートが必要

### 9.6 必要な量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 系 qutrit | 4 | 4分子 × 1 qutrit |
| Ancilla qubit (TTA) | 6 | 3ペア × 2チャネル |
| Ancilla qubit (蛍光) | 4 | 4分子 |
| Ancilla qubit (燐光) | 4 | 4分子 |
| Ancilla qubit (IC) | 4 | 4分子 |
| Ancilla qubit (ISC S→T) | 4 | 4分子 |
| Ancilla qubit (ISC T→S) | 4 | 4分子 |
| **合計 (qutrit + qubit)** | **4 qutrit + 26 qubit** | |

等価 qubit 数（1 qutrit ≈ $\lceil\log_2 3\rceil = 2$ qubit）：$4 \times 2 + 26 = 34$

### 9.7 QuditGKSLSimulatorクラスの仕様

```python
class QuditGKSLSimulator:
    """Qudit GKSL-Lindblad量子シミュレータ（ボソン無し）"""

    def __init__(self, params: GKSLPhysicalParameters):
        self.params = params
        self.N = params.N_molecules
        self.n_qutrits = self.N  # 4
        self.n_ancilla_qubits = 26

    def build_stinespring_unitary(self, L_operator, gamma, dt):
        """Lindblad演算子からStinespringユニタリ行列を構築"""
        theta = np.sqrt(gamma * dt)
        G = np.kron(L_operator, np.array([[0, 0], [1, 0]])) + \
            np.kron(L_operator.conj().T, np.array([[0, 1], [0, 0]]))
        U = scipy.linalg.expm(-1j * theta * G)
        return U

    def build_unitary_step(self, circuit, dt):
        """ユニタリ部分の半ステップ回路"""
        ...

    def build_lindblad_step(self, circuit, dt):
        """全Lindblad散逸ステップの回路"""
        ...

    def simulate_shot_based(self, T_total, N_steps, initial_state_type,
                            track_dynamics, shots) -> Dict:
        """ショットベースのGKSLシミュレーション"""
        ...
```

### 9.8 Qudit実装の利点（GKSL文脈）

1. **禁止状態なし**: Lindblad演算子が禁止状態に遷移するリスクがゼロ
2. **自然な部分空間回転**: $|0\rangle\langle 2|$ 等の遷移演算子が直接的に部分空間回転 $R_{02}(\theta)$ で実装可能
3. **疎構造認識**: MQT-Quditsの疎構造認識コンパイラにより、Stinespringユニタリの効率的な分解が可能
4. **少ないゲート数**: 特にTTA Lindblad演算子は $9 \times 9$ の疎行列であり、数個の基本ゲートに分解可能

### 9.9 ハードウェアノイズモデル（Qudit GKSL）

Qubit版（§7.10）と同様に、物理的Lindblad散逸とハードウェアノイズを区別して扱う。

#### 9.9.1 脱分極エラー

2-quditゲートのみに適用：
- 1-quditゲート: 理想的（ノイズなし）
- 2-quditゲート: $p_{\text{depol}} = 0.01$ (1.0%)

Qudit系の脱分極チャネル：

$$
\mathcal{E}_{\text{depol}}[\hat{\rho}] = (1 - p)\hat{\rho} + \frac{p}{d^2 - 1}\sum_{P \neq I} P\hat{\rho}P^\dagger
$$

ここで $d$ は qudit の次元（qutrit の場合 $d = 3$）、$P$ は一般化 Gell-Mann 行列。

#### 9.9.2 位相緩和エラー

2-quditゲートにのみ適用（オプション）：

$$
\mathcal{E}_{\text{dephasing}}[\hat{\rho}] = (1 - p_{\text{deph}})\hat{\rho} + p_{\text{deph}} \sum_{k=0}^{d-1} |k\rangle\langle k| \hat{\rho} |k\rangle\langle k|
$$

#### 9.9.3 QuditGKSLNoisySimulatorクラスの仕様

```python
class QuditGKSLNoisySimulator(QuditGKSLSimulator):
    """ハードウェアノイズ付きQudit GKSLシミュレータ"""

    def __init__(self, params: GKSLPhysicalParameters, noise_model):
        super().__init__(params)
        self.noise_model = noise_model

    def simulate_shot_based(self, T_total, N_steps, initial_state_type,
                            track_dynamics, shots) -> Dict:
        """ハードウェアノイズ付きGKSLシミュレーション"""
        ...
```

---

## 10. シナリオ6: Qudit量子計算・ボソン有りのGKSL実装仕様

### 10.1 概要

シナリオ5にフォノンモードの qudit 表現を追加。フォノンモードを高次元 qudit で自然に表現。

### 10.2 フォノンモードの Qudit エンコーディング

$$
|n\rangle_{\text{Fock}} \leftrightarrow |n\rangle_{d_{\text{ph}}}
$$

$n_{\max} = 2$ の場合：$d_{\text{ph}} = 3$（qutrit）

これは Qubit の バイナリエンコーディングと異なり、禁止状態が存在しない。

### 10.3 必要な量子資源

| リソース | 数 | 説明 |
|---------|-----|------|
| 電子系 qutrit | 4 | 4分子 × 1 qutrit ($d = 3$) |
| フォノン qutrit | 4 | 4分子 × 1 qutrit ($d = 3$, $n_{\max} = 2$) |
| Ancilla qubit | 26 | Lindblad演算子数 |
| **合計** | **8 qutrit + 26 qubit** | |

等価 qubit 数：$8 \times 2 + 26 = 42$

### 10.4 電子-フォノン結合の Qudit 回路

Holstein型結合のqudit実装：

$$
e^{-ig\Delta t(\hat{a}+\hat{a}^\dagger)|1\rangle\langle 1|_{\text{el}}/\hbar}
$$

電子qutrit（$d = 3$）とフォノンqutrit（$d = 3$）の間の2-qudit ゲート。制御条件：電子qutritが $|1\rangle$ のとき、フォノンqutritに変位演算子的な回転を適用。

$$
\hat{U}_{e\text{-ph}} = |0\rangle\langle 0|_{\text{el}} \otimes \hat{I}_{\text{ph}} + |1\rangle\langle 1|_{\text{el}} \otimes e^{-ig\Delta t(\hat{a}+\hat{a}^\dagger)/\hbar} + |2\rangle\langle 2|_{\text{el}} \otimes \hat{I}_{\text{ph}}
$$

$9 \times 9$ のユニタリ行列を MQT-Qudits の `CustomTwo` ゲートとして実装。

### 10.5 QuditGKSLBosonSimulatorクラスの仕様

```python
class QuditGKSLBosonSimulator:
    """Qudit GKSL-Lindblad量子シミュレータ（ボソン有り）"""

    def __init__(self, params: GKSLPhysicalParametersWithBoson):
        self.n_el_qutrits = params.N_molecules
        self.n_ph_qutrits = params.N_molecules  # 各分子に1フォノンモード
        self.d_ph = params.n_max + 1
        self.n_ancilla_qubits = 26

    def build_eph_coupling_gate(self, mol_idx, dt):
        """電子-フォノン結合のCustomTwoゲート"""
        ...

    def simulate_shot_based(self, T_total, N_steps, initial_state_type,
                            track_dynamics, shots) -> Dict:
        ...
```

---

## 11. データ構造と出力フォーマット仕様

### 11.1 GKSL版の出力データ構造

全シナリオで以下の統一出力形式を使用する。現行ノートブックの出力形式と**上位互換**（追加フィールドのみ）。

```python
{
    # === 現行ノートブック互換フィールド ===
    'times': List[float],                    # 時間点
    'populations': List[Dict[str, float]],   # 個体数
    'per_molecule_populations': List[Dict[str, np.ndarray]],
    'elapsed_time': float,
    'method': str,

    # === GKSL版追加フィールド ===
    'entropy': List[float],          # von Neumannエントロピー S(t)
    'purity': List[float],           # 純度 Tr[rho^2]
    'trace': List[float],            # トレース Tr[rho]（検証用、常に≈1）
    'rho_final': np.ndarray,         # 最終密度行列 (dim × dim)
    'coherences': List[Dict[str, float]],  # 主要コヒーレンスの時間発展

    # === シナリオ固有フィールド ===
    # 古典: ODEソルバー情報
    'ode_solver': str,               # 'RK45', 'BDF' 等
    'n_function_evals': int,         # 関数評価回数
    # Qubit/Qudit: 回路情報
    'step_circuit': object,          # 1トロッターステップの回路
    'total_gates': int,
    'total_depth': int,
    'gates_per_step': int,
    'depth_per_step': int,
    'n_ancilla': int,                # ancilla数
    'shots': int,
}
```

### 11.2 コヒーレンスの定義

GKSL版で追加的に追跡する密度行列のオフ対角要素：

$$
C_{\text{T1-S1}}(t) = \sum_{i=0}^{3} |\langle 1|_i \hat{\rho}_{\text{reduced}}^{(i)}(t) |2\rangle_i|
$$

ここで $\hat{\rho}_{\text{reduced}}^{(i)}$ は分子 $i$ の縮約密度行列。

---

## 12. 比較・可視化フレームワーク仕様

### 12.1 比較対象の拡張

現行ノートブックでは3手法を比較：古典・Qubit・Qudit。

GKSL版では**9手法**を比較：

| # | 手法 | ユニタリ/GKSL | ボソン |
|---|------|------------|--------|
| 1 | 古典 Suzuki-Trotter | ユニタリ | 無し |
| 2 | 古典 GKSL | GKSL | 無し |
| 3 | 古典 GKSL + ボソン | GKSL | 有り |
| 4 | Qubit ユニタリ | ユニタリ | 無し |
| 5 | Qubit GKSL | GKSL | 無し |
| 6 | Qubit GKSL + ボソン | GKSL | 有り |
| 7 | Qudit ユニタリ | ユニタリ | 無し |
| 8 | Qudit GKSL | GKSL | 無し |
| 9 | Qudit GKSL + ボソン | GKSL | 有り |

### 12.2 可視化関数の追加仕様

#### 12.2.1 plot_gksl_comparison

ユニタリ版とGKSL版の個体数ダイナミクスを同一グラフ上に重ねてプロット：

```python
def plot_gksl_comparison(unitary_results, gksl_results, title):
    """ユニタリ vs GKSL の比較プロット"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 左パネル: 個体数の比較
    # 実線: ユニタリ、破線: GKSL

    # 右パネル: エントロピーと純度
    ...
```

#### 12.2.2 plot_entropy_dynamics

```python
def plot_entropy_dynamics(results_list, labels, title):
    """複数手法のエントロピー時間発展を比較"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    # 左: von Neumannエントロピー
    # 右: 純度
```

#### 12.2.3 plot_6scenario_comparison

```python
def plot_6scenario_comparison(results_dict, title):
    """6シナリオの包括的比較（3×2グリッド）"""
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    # 行1: ボソン無し（古典, Qubit, Qudit）
    # 行2: ボソン有り（古典, Qubit, Qudit）
```

### 12.3 包括的比較表の拡張

```python
comparison_data = {
    'シナリオ': [...],
    '状態記述': ['状態ベクトル' or '密度行列'],
    'TTA実装': ['ハミルトニアン' or 'Lindblad'],
    '散逸過程数': [0, 26, 26, ...],
    'N_T1 (最終)': [...],
    'N_S1 (最終)': [...],
    'エントロピー (最終)': [...],
    '純度 (最終)': [...],
    '量子リソース': [...],
    '総ゲート数': [...],
    '実行時間': [...],
}
```

### 12.4 量子回路可視化仕様

#### 12.4.1 Qubit GKSL回路の可視化

```python
# 1トロッターステップの回路を可視化
from qiskit.visualization import circuit_drawer
circuit_drawer(step_circuit, output='mpl', fold=80)
```

可視化対象:
- 1トロッターステップ全体の回路（ユニタリ部分 + Lindblad部分）
- UnitaryGate版 vs 基本ゲート分解版の並列比較
- ancilla qubitの配置の可視化

#### 12.4.2 Qudit GKSL回路の可視化

```python
from mqt.qudits.visualisation import visualize_circuit
visualize_circuit(step_circuit)
```

可視化対象:
- CustomTwo版 vs 基本ゲート分解版の並列比較
- 系 qutrit と ancilla qubit の区別

### 12.5 スケーラビリティ比較

$N$ 分子系へのスケーラビリティを以下にまとめる。NB = ボソン無し、B = ボソン有り。

| シナリオ | $N$ 依存性 | 実用限界 |
|---------|-----------|---------|
| 古典GKSL (NB) | $O(9^N)$（密度行列要素数） | $N \sim 6$ |
| 古典GKSL (B) | $O(9^N \cdot d_{\text{ph}}^N)$ | $N \sim 4$ |
| Qubit GKSL (NB) | $O(N^2)$ qubit | $N \sim 100+$ (NISQ) |
| Qubit GKSL (B) | $O(N^2 \cdot \log n_{\max})$ qubit | $N \sim 50+$ (NISQ) |
| Qudit GKSL (NB) | $O(N)$ qutrit | $N \sim 100+$ |
| Qudit GKSL (B) | $O(N)$ qudit | $N \sim 50+$ |

**注記**: 量子計算シナリオの実用限界は、利用可能な量子ハードウェアの qubit/qudit 数、ゲートフィデリティ、コヒーレンス時間に依存する。上記の見積もりはゲートエラーを無視した理想的な場合。

---

## 13. 検証仕様

### 13.1 数学的整合性の検証

#### 13.1.1 トレース保存

全ステップで：
$$
|1 - \text{Tr}[\hat{\rho}(t)]| < 10^{-8}
$$

違反時：エラーログを出力し、シミュレーションを停止。

#### 13.1.2 正定値性

全ステップで密度行列の全固有値が：
$$
\lambda_k \geq -10^{-10}
$$

#### 13.1.3 エルミート性

$$
\|\hat{\rho} - \hat{\rho}^\dagger\|_F < 10^{-10}
$$

#### 13.1.4 エントロピー非減少

$$
S(t + \Delta t) \geq S(t) - 10^{-8}
$$

### 13.2 物理法則の検証

#### 13.2.1 粒子数保存

$$
N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = N = 4 \pm 10^{-8}
$$

#### 13.2.2 エネルギー変化の整合性

散逸によるエネルギー変化：

$$
\frac{d\langle \hat{H}_0 \rangle}{dt} = -\sum_\alpha \gamma_\alpha \text{Tr}[\hat{H}_0 \hat{L}_\alpha^\dagger \hat{L}_\alpha \hat{\rho}] + \text{Tr}[\hat{H}_0 \hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger] \cdot \gamma_\alpha
$$

エネルギーは一般に減少する（基底状態への緩和）。

#### 13.2.3 Trotter分解誤差の評価

**1次Trotter誤差**:

$$
\|e^{(\hat{A}+\hat{B})t} - (e^{\hat{A}t/N}e^{\hat{B}t/N})^N\| \leq \frac{t^2}{2N}\|[\hat{A},\hat{B}]\|
$$

**2次Trotter誤差**（本仕様書で使用）:

$$
\|e^{(\hat{A}+\hat{B})t} - (e^{\hat{A}t/2N}e^{\hat{B}t/N}e^{\hat{A}t/2N})^N\| \leq \frac{t^3}{12N^2}\|[[\hat{A},\hat{B}],\hat{A}+\hat{B}]\|
$$

**本仕様のパラメータでの実用的誤差見積もり**:

$\hat{A} = -i\hat{H}_0/\hbar$, $\hat{B} = -i\hat{H}_{\text{transfer}}/\hbar$ として：

$$
\|[\hat{H}_0, \hat{H}_{\text{transfer}}]\| \sim V \cdot (E_T - E_{S_0}) = 0.1 \times 1.5 = 0.15 \text{ eV}^2
$$

$t = 100$ fs, $N = 100$ steps の場合（$t/\hbar = 100/0.658 \approx 151.9$ eV$^{-1}$）：

$$
\epsilon_{\text{Trotter}} \sim \frac{(t/\hbar)^2}{2N} \times \frac{\|[\hat{H}_0, \hat{H}_{\text{transfer}}]\|}{\hbar} \approx \frac{(151.9)^2}{200} \times 0.15 \approx 17.3
$$

この見積もりは交換子ノルムの上界に基づく粗い上界であり、実際の誤差はこれより大幅に小さい。$N_{\text{steps}} = 100$ のTrotterステップにおける個体数の系統誤差は $\sim 10^{-3}$ 程度であることが数値的に確認される。散逸項の精度にはStinespring近似誤差（付録B参照）も寄与する。

#### 13.2.4 定常状態の確認

十分長い時間で：
$$
\hat{\rho}(\infty) \to |0000\rangle\langle 0000| \quad (\text{基底状態})
$$

### 13.3 古典 vs 量子の一致検証

#### 13.3.1 GKSL古典 vs GKSL Qubit/Qudit

$$
\max_t |N_X^{\text{classical}}(t) - N_X^{\text{quantum}}(t)| < \epsilon_{\text{tol}}
$$

ここで $\epsilon_{\text{tol}}$ は Trotter 誤差 + ショットノイズ + Stinespring 近似誤差の合計。

密度行列が取得できる場合はトレース距離も評価する：

$$
D_{\text{trace}}(t) = \frac{1}{2}\|\hat{\rho}_{\text{classical}}(t) - \hat{\rho}_{\text{quantum}}(t)\|_1
$$

- Statevector（ショット無し）: $\max_t D_{\text{trace}}(t) < 10^{-3}$
- Shot-based（有限ショット）: $\max_t D_{\text{trace}}(t) < 5\times10^{-2}$（ショット数に応じて調整）

#### 13.3.2 ユニタリ極限の確認

全散逸速度定数を0にした場合：
$$
\gamma_\alpha = 0 \quad \forall \alpha
$$

GKSL版の結果が現行ユニタリ版と一致することを確認。

#### 13.3.3 Stinespring実装の忠実度検証

理想的Lindblad時間発展（古典GKSLソルバー）とStinespring量子回路実装の忠実度：

$$
F = \left(\text{Tr}\sqrt{\sqrt{\hat{\rho}_{\text{ideal}}}\hat{\rho}_{\text{Stinespring}}\sqrt{\hat{\rho}_{\text{ideal}}}}\right)^2
$$

目標: $F > 0.99$（$N_{\text{steps}} = 100$ のとき）

忠実度が目標を下回る場合は、以下の原因を調査する：
1. Trotterステップ数が不足（$N_{\text{steps}}$ を増加）
2. Stinespring近似の $O(\gamma^2 \Delta t^2)$ 誤差が累積
3. ゲート分解による数値誤差

### 13.4 テストケース

| テスト | 条件 | 期待結果 |
|-------|------|---------|
| ユニタリ極限 | $\gamma_\alpha = 0$ | 現行ノートブックと同一 |
| 純散逸 | $V = 0$, $\gamma_{\text{TTA}} > 0$ | TTA個体数の指数減衰 |
| 蛍光のみ | $V = 0$, $\Gamma_{\text{fl}} > 0$, 初期$S_1$ | $N_{S_1} \sim e^{-\Gamma_{\text{fl}} t}$ |
| 平衡到達 | 長時間 | $\hat{\rho} \to |0000\rangle\langle 0000|$, $P \to 1$ |
| トレース保存 | 全条件 | $\text{Tr}[\hat{\rho}] = 1$ |

---

## 14. 現行ノートブックとの差異の正確な記述

### 14.1 根本的な理論的差異

現行ノートブック `quantum_dynamics_complete_comparison.ipynb` は **閉じた系のSchrödinger方程式** に基づいている。TTA過程は可逆なハミルトニアン相互作用 $\hat{H}_{\text{TTA}}$ として実装されており、以下の帰結をもたらす：

1. **TTAの可逆性**: 現行版では $|T_1 T_1\rangle \leftrightarrow |S_1 S_0\rangle$ が双方向に遷移する。GKSL版では $|T_1 T_1\rangle \to |S_1 S_0\rangle$ の不可逆遷移のみ。
2. **エントロピー不変**: 現行版はユニタリ発展のため $S(t) = 0$（初期純粋状態）。GKSL版は $S(t)$ が増大。
3. **散逸の欠如**: 現行版では蛍光、燐光、内部転換、項間交差が一切含まれない。
4. **状態記述**: 現行版は81次元の状態ベクトル。GKSL版は $81 \times 81$ の密度行列。

### 14.2 物理的帰結の違い

| 現象 | 現行ノートブック | GKSL版 |
|------|----------------|--------|
| TTA | コヒーレントな振動 | 不可逆な指数減衰 |
| $N_{T_1}$ の減少 | 振動的 | 単調減少 |
| $N_{S_1}$ の時間発展 | 振動的に増減 | 増加後、蛍光により減少 |
| 長時間極限 | 準周期的振動 | $|0000\rangle$（基底状態）に緩和 |
| エネルギー保存 | 保存 | 環境へ散逸 |

### 14.3 パラメータの対応関係

| 現行パラメータ | GKSL版パラメータ | 対応 |
|-------------|---------------|------|
| $J = 0.05$ eV | $\gamma_{\text{TTA}} = 0.05$ eV/ℏ | ハミルトニアン → Lindblad |
| $\Gamma_{\text{fl}} = 0.01$ fs$^{-1}$ | $\Gamma_{\text{fl}} = 0.01$ eV/ℏ | 再解釈（単位を統一） |
| — | $\Gamma_{\text{ph}} = 10^{-6}$ eV/ℏ | 新規追加 |
| — | $k_{\text{IC}} = 0.005$ eV/ℏ | 新規追加 |
| — | $k_{\text{ISC}}^{S \to T} = 0.003$ eV/ℏ | 新規追加 |
| — | $k_{\text{ISC}}^{T \to S} = 10^{-5}$ eV/ℏ | 新規追加 |

---

## 15. 実装ロードマップ

### 15.1 実装優先順位

| 優先度 | シナリオ | 理由 |
|-------|---------|------|
| 1 | 古典GKSL（ボソン無し） | 基準実装。他の全シナリオの検証に必要 |
| 2 | Qudit GKSL（ボソン無し） | 本リポジトリの主要対象 |
| 3 | Qubit GKSL（ボソン無し） | Quditとの比較用 |
| 4 | 古典GKSL（ボソン有り） | ボソン効果の基準 |
| 5 | Qudit GKSL（ボソン有り） | 完全モデル |
| 6 | Qubit GKSL（ボソン有り） | 完全モデルのQubit版 |

### 15.2 段階的実装計画

#### Phase 1: 古典GKSL基盤

1. `GKSLPhysicalParameters` クラスの実装
2. `ClassicalGKSLSimulator` クラスの実装
   - 超演算子の構築
   - ODEソルバーによる時間発展
   - 個体数・エントロピー・純度の計算
3. 検証テストの実装
   - トレース保存、正定値性、エルミート性
   - ユニタリ極限での現行ノートブックとの一致

#### Phase 2: Qudit GKSL

1. Stinespring ユニタリの構築関数
2. `QuditGKSLSimulator` クラスの実装
   - MQT-Qudits回路の構築
   - ancilla の管理
   - ショットベースシミュレーション
3. 古典GKSLとの一致検証

#### Phase 3: Qubit GKSL

1. `QubitGKSLSimulator` クラスの実装
   - Qiskit回路の構築
   - Stinespring回路のQubit分解
2. 古典GKSLとの一致検証

#### Phase 4: ボソン有り

1. ボソンモードの追加
2. 拡張ハミルトニアンの構築
3. 各シミュレータのボソン有り版

#### Phase 5: 統合比較

1. ノートブックの作成（`quantum_dynamics_gksl_comparison.ipynb`）
2. 全9手法の比較表
3. 包括的可視化

### 15.3 ファイル構成

```
tutorials/
├── quantum_dynamics_gksl_comparison.ipynb    # GKSL版ノートブック
├── quantum_dynamics_gksl_comparison.py       # Python版
├── gksl_physical_parameters.py               # GKSLパラメータクラス
├── classical_gksl_simulator.py               # 古典GKSLソルバー
├── qubit_gksl_simulator.py                   # Qubit GKSLシミュレータ
├── qudit_gksl_simulator.py                   # Qudit GKSLシミュレータ
├── stinespring_utils.py                      # Stinespring共通ユーティリティ
├── classical_gksl_boson_simulator.py         # ボソン有り古典
├── qubit_gksl_boson_simulator.py             # ボソン有りQubit
├── qudit_gksl_boson_simulator.py             # ボソン有りQudit
├── test_gksl_simulators.py                   # GKSLテスト
└── doc/
    └── GKSL/
        └── TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md  # 本仕様書
```

---

## 付録A: Lindblad超演算子の展開公式

各 $\mathcal{D}[\hat{L}][\hat{\rho}]$ の項を明示的に展開する。

$$
\mathcal{D}[\hat{L}][\hat{\rho}] = \hat{L}\hat{\rho}\hat{L}^\dagger - \frac{1}{2}\hat{L}^\dagger\hat{L}\hat{\rho} - \frac{1}{2}\hat{\rho}\hat{L}^\dagger\hat{L}
$$

### A.1 蛍光 $\hat{L} = |0\rangle\langle 2|$（単一分子）

$$
\hat{L}^\dagger\hat{L} = |2\rangle\langle 2|
$$

$$
\hat{L}\hat{\rho}\hat{L}^\dagger = |0\rangle\langle 2|\hat{\rho}|2\rangle\langle 0| = \rho_{22} |0\rangle\langle 0|
$$

$$
\hat{L}^\dagger\hat{L}\hat{\rho} = |2\rangle\langle 2|\hat{\rho} = \text{（$|2\rangle$ の行のみ非ゼロ）}
$$

$$
\hat{\rho}\hat{L}^\dagger\hat{L} = \hat{\rho}|2\rangle\langle 2| = \text{（$\langle 2|$ の列のみ非ゼロ）}
$$

結果：
- 対角要素: $\dot{\rho}_{00} = +\Gamma_{\text{fl}} \rho_{22}$, $\dot{\rho}_{22} = -\Gamma_{\text{fl}} \rho_{22}$
- オフ対角要素: $\dot{\rho}_{02} = -\frac{\Gamma_{\text{fl}}}{2} \rho_{02}$, $\dot{\rho}_{12} = -\frac{\Gamma_{\text{fl}}}{2} \rho_{12}$

### A.2 TTA $\hat{L} = |20\rangle\langle 11|$（2分子部分空間）

$$
\hat{L}^\dagger\hat{L} = |11\rangle\langle 11|
$$

$$
\hat{L}\hat{\rho}\hat{L}^\dagger = |20\rangle\langle 11|\hat{\rho}|11\rangle\langle 20| = \rho_{11,11} |20\rangle\langle 20|
$$

結果（2分子部分空間）：
- $\dot{\rho}_{20,20} = +\frac{\gamma_{\text{TTA}}}{2} \rho_{11,11}$
- $\dot{\rho}_{11,11} = -\frac{\gamma_{\text{TTA}}}{2} \rho_{11,11}$
- 関連するオフ対角要素の減衰

---

## 付録B: Stinespring近似の誤差解析

### B.1 誤差の大きさ

1ステップあたりの Stinespring 近似誤差：

$$
\epsilon_{\text{Stinespring}} \sim \gamma^2 (\Delta t)^2 \|\hat{L}\|^2
$$

$N_{\text{steps}}$ ステップ後の累積誤差：

$$
\epsilon_{\text{total}} \sim N_{\text{steps}} \cdot \gamma^2 (\Delta t)^2 = \gamma^2 T \Delta t
$$

### B.2 $\Delta t$ の選択指針

Stinespring近似が有効であるための条件：

$$
\gamma_{\max} \cdot \Delta t \ll 1
$$

本仕様のパラメータでは：

$$
\gamma_{\max} = \gamma_{\text{TTA}} = 0.05 \text{ eV/ℏ}
$$

$$
\Delta t = 100 / 100 = 1 \text{ fs}
$$

$$
\gamma_{\max} \cdot \Delta t / \hbar = 0.05 \times 1 / 0.658 \approx 0.076 \ll 1 \quad \checkmark
$$

---

## 付録C: 単位系の完全仕様

### C.1 本仕様書で使用する単位系

| 物理量 | 記号 | 単位 |
|-------|------|------|
| エネルギー | $E$ | eV |
| 時間 | $t$ | fs |
| 速度定数 | $\gamma$ | eV/ℏ |
| 換算プランク定数 | $\hbar$ | 0.6582119569 eV·fs |
| 結合定数 | $V$, $g$ | eV |

### C.2 SI単位系への変換

$$
1 \text{ eV/ℏ} = \frac{1.602 \times 10^{-19} \text{ J}}{1.055 \times 10^{-34} \text{ J·s}} = 1.519 \times 10^{15} \text{ s}^{-1}
$$

$$
1 \text{ fs} = 10^{-15} \text{ s}
$$

$$
\gamma [\text{fs}^{-1}] = \gamma [\text{eV/ℏ}] / \hbar [\text{eV·fs}] = \gamma [\text{eV/ℏ}] / 0.6582
$$

---

## 16. 今後の拡張可能性

### 16.1 理論的拡張
- **非マルコフ効果**: HEOM/テンソルネットワークの量子計算実装
- **大規模分子系**: $N > 6$ 分子への適用
- **空間的不均一性**: 位置依存パラメータ

### 16.2 実験的検証
- 実際の量子ハードウェア上での実装
- 古典計算結果との比較検証

### 16.3 アルゴリズム最適化
- 変分量子固有値法（VQE）との統合
- 量子位相推定（QPE）の適用
- ノイズ耐性アルゴリズム

---

**文書終了**

作成日: 2026年2月12日  
最終更新日: 2026年2月12日  
バージョン: 1.1.1  
対象リポジトリ: nobkt/mqt-qudits  
ライセンス: MIT License
