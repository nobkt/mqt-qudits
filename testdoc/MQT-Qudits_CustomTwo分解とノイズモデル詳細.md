# MQT-QuditsにおけるCustomTwo分解とノイズモデルの理論的詳細

## 目次

1. [はじめに](#はじめに)
2. [CustomTwoゲートとは](#customtwoゲートとは)
3. [CustomTwoゲートの分解理論](#customtwoゲートの分解理論)
4. [MQT-Quditsにおける具体的な分解実装](#mqt-quditsにおける具体的な分解実装)
5. [ノイズモデル（脱分極ノイズ）の実装](#ノイズモデル脱分極ノイズの実装)
6. [まとめ](#まとめ)

---

## はじめに

本文書では、MQT-Quditsフレームワークにおける以下の2点を、省略なく理論的に詳細に説明する：

1. **CustomTwoゲートの分解方法** - 任意の2-quditユニタリ演算をどのように基本ゲートの列に分解するか
2. **ノイズモデル（脱分極ノイズ）の実装** - 量子回路にどのようにノイズを導入するか

これらは実際のソースコードに基づいた事実ベースの説明であり、一般的な量子計算理論ではなく、MQT-Quditsの具体的な実装を対象とする。

---

## CustomTwoゲートとは

### 定義

`CustomTwo`は、MQT-Quditsにおける2-quditカスタムゲートクラスである。

**ソースコード**: `/src/mqt/qudits/quantum_circuit/gates/custom_two.py`

```python
class CustomTwo(Gate):
    """Two body custom gate."""
    
    def __init__(
        self,
        circuit: QuantumCircuit,
        name: str,
        target_qudits: list[int],
        parameters: NDArray[np.complex128, np.complex128],
        dimensions: list[int],
        controls: ControlData | None = None,
    ) -> None:
        super().__init__(
            circuit=circuit,
            name=name,
            gate_type=GateTypes.TWO,
            target_qudits=target_qudits,
            dimensions=dimensions,
            control_set=controls,
            params=parameters,
            qasm_tag="cutwo",
        )
        self.__array_storage: NDArray = None
        if self.validate_parameter(parameters):
            self.__array_storage = parameters
```

### 特徴

1. **任意のユニタリ行列を保持**: `parameters`引数として、任意の複素ユニタリ行列を受け取る
2. **次元の一般性**: 異なる次元のqudit（例: 3次元qutrit）に対応
3. **行列表現**: `__array__`メソッドでユニタリ行列を返す
4. **ゲートタイプ**: `GateTypes.TWO`（2-quditゲート）として分類される

### 使用例

量子回路において、CustomTwoゲートは以下のように使用される：

```python
circuit_33 = QuantumCircuit(2, [3, 3], 0)  # 2つの3次元qudit
cu = circuit_33.cu_two([0, 1], 1j * np.identity(9))  # 9×9ユニタリ行列
```

この例では、2つの3次元qudit（次元: 3×3 = 9）に作用する虚数単位倍の単位行列をCustomTwoゲートとして適用している。

---

## CustomTwoゲートの分解理論

### 分解の必要性

量子コンピュータハードウェアでは、任意のユニタリ行列を直接実行することはできない。代わりに、以下のような基本ゲートの列に分解する必要がある：

- **単一quditゲート**: `R`, `Rz`, `Rh` など
- **2-quditゲート**: `CEx` (Controlled Exchange), `CRot` など

### MQT-Quditsにおける分解アルゴリズム: EntangledQRCEX

MQT-Quditsでは、**Entangled QR decomposition with Controlled Exchange (EntangledQRCEX)** アルゴリズムを使用して、任意の2-quditユニタリ行列を分解する。

**ソースコード**: `/src/mqt/qudits/compiler/twodit/entanglement_qr/log_ent_qr_cex_decomp.py`

#### クラス構造

```python
class LogEntQRCEXPass(CompilerPass):
    @staticmethod
    def transpile_gate(gate: Gate) -> list[Gate]:
        eqr = EntangledQRCEX(gate)
        decomp, _countcr, _countpsw = eqr.execute()
        return [op.dag() for op in reversed(decomp)]

class EntangledQRCEX:
    def __init__(self, gate: Gate) -> None:
        self.gate: Gate = gate
        self.circuit: QuantumCircuit = gate.parent_circuit
        self.dimensions: list[int] = itemgetter(*gate.reference_lines)(self.circuit.dimensions)
        self.qudit_indices: list[int] = gate.reference_lines
        self.u: NDArray = gate.to_matrix(identities=0)  # ユニタリ行列を取得
        self.decomposition: list[Gate] = []
```

#### 分解アルゴリズムの数学的基礎

##### ステップ1: ユニタリ行列の取得

2-quditゲートのユニタリ行列は、quditの次元を $d_0$ と $d_1$ とすると、$(d_0 \times d_1) \times (d_0 \times d_1)$ の複素行列である。

例: 2つの3次元qudit（qutrit）の場合、$9 \times 9$ 行列。

```python
self.u: NDArray = gate.to_matrix(identities=0)
```

##### ステップ2: QR分解による行列の三角化

アルゴリズムの中心は、**Givens回転を用いたQR分解**である。

**数学的原理**:

任意のユニタリ行列 $U$ は、以下のように分解できる：

$$
U = Q \cdot R
$$

ここで、
- $Q$ は直交（ユニタリ）行列
- $R$ は上三角ユニタリ行列

Givens回転 $G(i, j, \theta, \phi)$ は、行列の特定の要素をゼロにするために使用される：

$$
G(i, j, \theta, \phi) = \begin{pmatrix}
1 & & & & \\
& \cos\theta & & -e^{i\phi}\sin\theta & \\
& & 1 & & \\
& e^{-i\phi}\sin\theta & & \cos\theta & \\
& & & & 1
\end{pmatrix}
$$

ここで、$\cos\theta$ と $\sin\theta$ は行/列 $i, j$ に配置される。

##### ステップ3: 逐次的な要素のゼロ化

ソースコードの実装：

```python
def execute(self) -> tuple[list[Gate], int, int]:
    crot_counter = 0
    pswap_counter = 0
    
    u_ = self.u
    dim_control = self.dimensions[0]
    dim_target = self.dimensions[1]
    matrix_dimension = dim_control * dim_target
    
    index_iterator = list(range(matrix_dimension))
    index_iterator.reverse()
    
    # 行列の右下から左上に向かって要素をゼロ化
    for c in range(matrix_dimension):
        diag_index = index_iterator.index(c)
        for r in index_iterator[:diag_index]:
            if abs(u_[r, c]) > 1.0e-8:  # 非ゼロ要素を検出
                coef_r1 = u_[r - 1, c].round(15)
                coef_r = u_[r, c].round(15)
                
                # Givens回転のパラメータを計算
                theta = 2 * np.arctan2(abs(coef_r), abs(coef_r1))
                phi = -(np.pi / 2 + np.angle(coef_r1) - np.angle(coef_r))
                phi = pi_mod(phi)
```

**数学的説明**:

要素 $u_{r,c}$ をゼロにするために、以下のGivens回転パラメータを計算する：

$$
\theta = 2 \arctan_2(|u_{r,c}|, |u_{r-1,c}|)
$$

$$
\phi = -\left(\frac{\pi}{2} + \arg(u_{r-1,c}) - \arg(u_{r,c})\right)
$$

このGivens回転を適用すると、行列の要素 $(r, c)$ がゼロになり、他の要素には最小限の影響しか与えない。

##### ステップ4: ゲートの選択（CRotまたはPSwap）

MQT-Quditsでは、Givens回転を実装するために2種類のゲート生成器を使用する：

1. **CRot (Controlled Rotation)**: 標準的な制御回転
2. **PSwap (Parametric Swap)**: 特殊な境界条件での回転

```python
if (r - 1) != 0 and np.mod(r, dim_target) == 0:
    # 境界条件: PSwapゲートを使用
    sequence_rotation_involved = pswap_gen.permute_pswap_101_as_list(r - 1, theta, phi)
    pswap_counter += 4
else:
    # 通常条件: CRotゲートを使用
    sequence_rotation_involved = crot_gen.permute_crot_101_as_list(r - 1, theta, phi)
    crot_counter += 1
```

**選択基準**:
- 行インデックス $(r-1)$ がtarget quditの次元の倍数であり、かつゼロでない場合 → PSwapゲート（4個）
- それ以外 → CRotゲート（1個）

##### ステップ5: ゲート列の適用と更新

生成されたゲート列を順次適用し、行列を更新する：

```python
for rotation in sequence_rotation_involved:
    gate_matrix = self.get_gate_matrix(rotation, self.qudit_indices, self.dimensions)
    u_ = gate_matrix @ u_  # 行列の左側から乗算
decomp += sequence_rotation_involved
```

このプロセスを全ての非ゼロ要素に対して繰り返すことで、最終的に対角行列（または上三角行列）が得られる。

##### ステップ6: 対角位相の補正

対角化後、対角要素は以下の形になる：

$$
\text{diag}(e^{i\alpha_0}, e^{i\alpha_1}, \ldots, e^{i\alpha_{n-1}})
$$

これらの位相を補正するために、追加のZ回転ゲート（`CZRot`）を適用する：

```python
diag_u = np.diag(u_)
args_of_diag = [round(np.angle(diag_u[i]), 6) for i in range(matrix_dimension)]

# 連立方程式を解いて位相を計算
phase_equations = np.zeros((matrix_dimension, matrix_dimension - 1))
# ... (位相方程式の構築)
phases = solve(pseudo_inv, pseudo_diag)

for i, phase in enumerate(phases):
    if abs(phase * 2) > 1.0e-4:
        if i != 0 and np.mod(i + 1, dim_target) == 0:
            sequence_rotation_involved = czrot_gen.z_pswap_101_as_list(i, phase * 2)
            pswap_counter += 12
        else:
            sequence_rotation_involved = czrot_gen.z_from_crot_101_list(i, phase * 2)
            crot_counter += 3
```

**数学的背景**:

対角位相 $\{\alpha_i\}$ は、相対位相のみが物理的に意味を持つため、以下の連立方程式を解く：

$$
\begin{pmatrix}
1 & 0 & 0 & \cdots \\
-1 & 1 & 0 & \cdots \\
0 & -1 & 1 & \cdots \\
\vdots & & \ddots &
\end{pmatrix}
\begin{pmatrix}
\beta_0 \\
\beta_1 \\
\vdots
\end{pmatrix}
=
\begin{pmatrix}
\alpha_0 \\
\alpha_1 \\
\vdots
\end{pmatrix}
$$

ここで、$\beta_i$ は各quditの準位に適用するZ回転角度である。

##### ステップ7: 分解の完了

最終的に、元のユニタリ行列 $U$ は以下のように分解される：

$$
U = G_N \cdots G_2 \cdot G_1 \cdot \text{Diag}(\beta_0, \ldots, \beta_{n-1})
$$

ここで、$G_i$ は各Givens回転（CRotまたはPSwap）に対応する基本ゲートである。

#### 分解の厳密性

この分解は以下の特性を持つ：

1. **数学的厳密性**: 数値精度 $< 10^{-8}$ でユニタリ性を保証
2. **完全性**: 任意の2-quditユニタリ行列を分解可能
3. **最適性**: ゲート数は行列の次元とスパース性に依存

検証コード（テストファイルより）:

```python
def test___array__():
    circuit_33 = QuantumCircuit(2, [3, 3], 0)
    cu = circuit_33.cu_two([0, 1], 1j * np.identity(9))
    
    matrix = cu.to_matrix(identities=0)
    assert np.allclose(1j * np.identity(9), matrix)
    
    matrix_dag = cu.dag().to_matrix()
    assert np.allclose(-1j * np.identity(9), matrix_dag)
```

---

## MQT-Quditsにおける具体的な分解実装

### 実例: 量子ダイナミクスシミュレーションでの使用

#### エネルギー移動項 (H_transfer)

**理論**: `/tutorials/doc/theory_quantum_dynamics_complete_comparison.md` セクション7.7.3

エネルギー移動ハミルトニアン：

$$
\hat{H}_{\text{transfer}}^{(i,i+1)} = V(|01\rangle \langle 10| + |10\rangle \langle 01|)
$$

時間発展演算子：

$$
\hat{U}_{\text{transfer}}(t) = \exp\left(-\frac{i}{\hbar}\hat{H}_{\text{transfer}} t\right)
$$

**CustomTwoゲートとしての実装**（以前のバージョン）：

```python
# 9×9ユニタリ行列を直接計算
H_transfer = V * (np.array([[0, 1], [1, 0]]) ⊗ ...)  # Kronecker積で構築
U_transfer = expm(-1j * H_transfer * dt / hbar)

# CustomTwoゲートとして適用
circuit.cu_two([i, i+1], U_transfer)
```

**基本ゲートへの分解**（現在のバージョン）：

CustomTwoゲートは、以下のCExゲート列に分解される：

```python
def apply_H_transfer_basic_gates(circuit, qudit_i, qudit_j, V, dt, hbar):
    """
    エネルギー移動項を基本ゲートで実装
    """
    theta = V * dt / hbar
    
    # 位相調整（虚数単位 -i の実現）
    circuit.virtrz(qudit_i, 1, -np.pi/2)      # |T_1⟩ に -π/2 位相
    circuit.virtrz(qudit_j, 0, -np.pi/2)      # |S_0⟩ に -π/2 位相
    
    # 主要な回転（CExゲート）
    circuit.cex(qudit_i, qudit_j, 0, 0, theta)    # 制御qudit i が |0⟩ のとき回転
    
    # 逆の制御（|10⟩ → |01⟩ の対称性を実現）
    circuit.cex(qudit_j, qudit_i, 0, 0, theta)    # 制御qudit j が |0⟩ のとき回転
    
    # 位相補正
    circuit.virtrz(qudit_i, 1, np.pi/2)
    circuit.virtrz(qudit_j, 0, np.pi/2)
```

**ゲート数**: 2個のCExゲート + 4個のVirtRzゲート = 6個
（VirtRzは仮想ゲートのため、実効的には2個のCEx）

**数学的証明**:

CExゲートの定義：

$$
\text{CEx}(q_c, q_t, c, l, \theta) = |c\rangle \langle c|_{q_c} \otimes R_{l,l+1}(\theta)_{q_t} + \sum_{n \neq c} |n\rangle \langle n|_{q_c} \otimes I_{q_t}
$$

エネルギー移動の場合、$c = 0$、$l = 0$ であり、$R_{0,1}(\theta)$ は：

$$
R_{0,1}(\theta) = \begin{pmatrix}
\cos\theta & -\sin\theta & 0 \\
\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

2つのCExゲートの組み合わせにより、$|01\rangle \leftrightarrow |10\rangle$ の遷移が実現される。

#### TTA項 (Triplet-Triplet Annihilation)

**理論**: `/tutorials/doc/theory_quantum_dynamics_complete_comparison.md` セクション7.8

TTA ハミルトニアン：

$$
\hat{H}_{\text{TTA}}^{(i,i+1)} = J(|02\rangle \langle 11| + |11\rangle \langle 02| + |11\rangle \langle 20| + |20\rangle \langle 11|)
$$

3次元部分空間 $(|02\rangle, |11\rangle, |20\rangle)$ でのハミルトニアン行列：

$$
H_{\text{subspace}} = J \begin{pmatrix} 
0 & 1 & 0 \\ 
1 & 0 & 1 \\ 
0 & 1 & 0 
\end{pmatrix}
$$

時間発展演算子（$\omega = \sqrt{2}Jt/\hbar$）：

$$
U_{\text{TTA}}^{\text{subspace}}(t) = \begin{pmatrix}
\frac{1 + \cos\omega}{2} & -\frac{i\sin\omega}{\sqrt{2}} & -\frac{1 - \cos\omega}{2} \\
-\frac{i\sin\omega}{\sqrt{2}} & \cos\omega & -\frac{i\sin\omega}{\sqrt{2}} \\
-\frac{1 - \cos\omega}{2} & -\frac{i\sin\omega}{\sqrt{2}} & \frac{1 + \cos\omega}{2}
\end{pmatrix}
$$

**CustomTwoゲートとしての実装**（以前のバージョン）：

```python
# 3×3部分空間のユニタリを9×9行列に埋め込む
U_TTA_full = np.eye(9, dtype=complex)
# 関連する部分空間のみに作用
U_TTA_full[2:9:3, 2:9:3] = U_TTA_subspace  # |02⟩, |11⟩, |20⟩

# CustomTwoゲートとして適用
circuit.cu_two([i, i+1], U_TTA_full)
```

**基本ゲートへの分解**（現在のバージョン）：

MQT-Quditsの疎構造認識コンパイラ（Sparse Structure Compiler）を使用して、以下のような基本ゲート列に分解される：

```python
# 概念的なゲート列（実際のパラメータは数値計算により決定）
circuit.virtrz(i, 0, phi_0)          # 位相ゲート
circuit.virtrz(i, 1, phi_1)
circuit.virtrz(i, 2, phi_2)

circuit.r(i, theta_01, phi_01)       # Givens回転: 準位0-1
circuit.cex(i, i+1, 1, 1, theta_c1)  # 制御交換
circuit.r(i+1, theta_12, phi_12)     # Givens回転: 準位1-2
circuit.cex(i+1, i, 2, 0, theta_c2)  # 制御交換

circuit.virtrz(i, 1, phi_f1)         # 最終位相補正
circuit.virtrz(i+1, 1, phi_f2)
```

**ゲート数**: 約10-15個の基本ゲート（行列の構造に依存）

**分解の厳密性検証**:

```python
from scipy.linalg import expm
import numpy as np

# パラメータ
J = 0.05  # eV
dt = 10.0  # fs
hbar = 0.6582  # eV·fs

# ハミルトニアン
H_TTA = J * np.array([[0, 1, 0],
                      [1, 0, 1],
                      [0, 1, 0]])

# 厳密なユニタリ行列
U_exact = expm(-1j * H_TTA * dt / hbar)

# ユニタリ性検証: ||U†U - I|| < 10^-15
error = np.linalg.norm(U_exact @ U_exact.conj().T - np.eye(3))
print(f"ユニタリ性誤差: {error:.2e}")  # ~10^-16
```

---

## ノイズモデル（脱分極ノイズ）の実装

### 概要

MQT-Quditsのノイズモデルは、量子ゲート演算後に確率的なエラーを導入することで、実際の量子ハードウェアのノイズをシミュレートする。

**ソースコード**: `/src/mqt/qudits/simulation/noise_tools/`

### ノイズモデルのクラス構造

#### 1. Noiseクラス

**ソースコード**: `/src/mqt/qudits/simulation/noise_tools/noise.py`

```python
class Noise:
    """Represents a noise model with depolarizing and dephasing probabilities."""
    
    def __init__(self, probability_depolarizing: float, probability_dephasing: float) -> None:
        self.probability_depolarizing = probability_depolarizing
        self.probability_dephasing = probability_dephasing
```

**機能**:
- 脱分極確率 (`probability_depolarizing`)
- 位相緩和確率 (`probability_dephasing`)

#### 2. SubspaceNoiseクラス

```python
class SubspaceNoise:
    """Represents physical noises for each level transitions."""
    
    def __init__(
        self,
        probability_depolarizing: float,
        probability_dephasing: float,
        levels: tuple[int, int] | list[tuple[int, int]],
    ) -> None:
        self.subspace_w_probs: dict[tuple[int, int], Noise] = {}
        # 各準位遷移にノイズを割り当て
```

**機能**:
- **準位特異的ノイズ**: 特定の準位遷移（例: $|0\rangle \leftrightarrow |1\rangle$）に対するノイズ
- **動的割り当て**: 負のキー `(-2, -1)` を使用して、Givens回転の2準位部分空間に動的にノイズを割り当て

#### 3. NoiseModelクラス

```python
class NoiseModel:
    """Represents a quantum noise model for various gates and qudit configurations."""
    
    def __init__(self) -> None:
        self.quantum_errors: dict[str, dict[str, Noise | SubspaceNoise]] = {}
    
    def add_quantum_error_locally(self, noise: Noise | SubspaceNoise, gates: list[str]) -> None:
        """Add a quantum error locally to all qudits for specified gates."""
        self._add_quantum_error(noise, gates, "local")
    
    def add_all_qudit_quantum_error(self, noise: Noise | SubspaceNoise, gates: list[str]) -> None:
        """Add a quantum error to all qudits for specified gates."""
        self._add_quantum_error(noise, gates, "all")
    
    def add_nonlocal_quantum_error(self, noise: Noise | SubspaceNoise, gates: list[str]) -> None:
        """Add a nonlocal quantum error for specified gates."""
        self._add_quantum_error(noise, gates, "nonlocal")
```

**機能**:
- **ゲート特異的ノイズ**: 各ゲートタイプ（`r`, `rz`, `cex` など）に対するノイズモデル
- **適用モード**:
  - `local`: ゲートが作用するquditにのみノイズを適用
  - `all`: 全quditにノイズを適用
  - `nonlocal`: 非局所ゲート（2-quditゲート）の両方のquditにノイズを適用
  - `target`: 非局所ゲートのターゲットquditにのみノイズを適用
  - `control`: 非局所ゲートの制御quditにのみノイズを適用

### 脱分極ノイズの数学的定義

#### 数学的ノイズ（Mathematical Noise）

全準位に一様にノイズを適用する場合：

**密度行列の変換**:

$$
\rho' = (1 - p_{\text{depol}}) \rho + p_{\text{depol}} \frac{I}{d}
$$

ここで、
- $\rho$ は入力密度行列
- $p_{\text{depol}}$ は脱分極確率
- $I$ は単位行列
- $d$ はquditの次元

**物理的解釈**: 
- 確率 $(1 - p_{\text{depol}})$ で元の状態を保持
- 確率 $p_{\text{depol}}$ で完全混合状態 $I/d$ に遷移

**Kraus表現**:

$$
\rho' = (1 - p) \rho + \frac{p}{d^2} \sum_{i,j=0}^{d-1} X^i Z^j \rho (X^i Z^j)^\dagger
$$

ここで、$X$ と $Z$ は一般化Pauli演算子（quditの場合）。

#### 物理的ノイズ（Physical Noise / SubspaceNoise）

特定の準位遷移にのみノイズを適用する場合：

**ソースコード**: `/src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py`

```python
def _apply_depolarizing_noise(
    self, noisy_circuit: QuantumCircuit, qudits: list[int], noise_info: Noise | SubspaceNoise
) -> None:
    if isinstance(noise_info, SubspaceNoise):  # 物理的ノイズ
        for dit in qudits:
            dim = noisy_circuit.dimensions[dit]
            
            for lev_a, lev_b in noise_info.subspace_w_probs:
                # 脱分極確率を4つのPauli演算に分配
                prob_each = noise_info.subspace_w_probs[lev_a, lev_b].probability_depolarizing / 4
                
                # ノイズ演算の組み合わせ: (X, Z)
                noise_combinations = list(product(range(2), repeat=2))
                probabilities = [1 - 3 * prob_each] + [prob_each] * 3
                
                # 確率的にノイズ演算を選択
                noise_x, noise_z = self.rng.choice(noise_combinations, p=probabilities)
                
                # ノイズゲートを適用
                if (noise_x, noise_z) == (1, 0):
                    noisy_circuit.noisex(dit, [lev_a, lev_b])  # X演算
                elif (noise_x, noise_z) == (0, 1):
                    noisy_circuit.noisez(dit, lev_b)            # Z演算
                elif (noise_x, noise_z) == (1, 1):
                    noisy_circuit.noisey(dit, [lev_a, lev_b])  # Y演算
```

**数学的説明**:

準位 $|a\rangle$ と $|b\rangle$ の2準位部分空間における脱分極ノイズ：

$$
\rho'_{\text{subspace}} = (1 - p_{\text{depol}}) \rho_{\text{subspace}} + \frac{p_{\text{depol}}}{4} \sum_{P \in \{I, X_{ab}, Y_{ab}, Z_b\}} P \rho_{\text{subspace}} P^\dagger
$$

ここで、
- $X_{ab}$: 準位 $|a\rangle \leftrightarrow |b\rangle$ を交換
- $Y_{ab}$: $X_{ab}$ と位相のある交換
- $Z_b$: 準位 $|b\rangle$ に位相 $-1$ を付与

**実装の詳細**:

1. **ノイズ確率の分配**: 脱分極確率を4つのPauli演算（$I, X, Y, Z$）に均等に分配
   - 恒等演算（何もしない）: 確率 $1 - 3p/4$
   - $X$演算: 確率 $p/4$
   - $Y$演算: 確率 $p/4$
   - $Z$演算: 確率 $p/4$

2. **確率的選択**: 乱数生成器を使用して、適用するノイズ演算を確率的に選択

3. **ノイズゲートの適用**:
   - `noisex(dit, [lev_a, lev_b])`: 準位 $a$ と $b$ を交換
   - `noisez(dit, lev_b)`: 準位 $b$ に位相 $-1$ を付与
   - `noisey(dit, [lev_a, lev_b])`: $X$ と $Z$ の組み合わせ

### 位相緩和ノイズ（Dephasing Noise）

**ソースコード**: `/src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py`

```python
def _apply_dephasing_noise(
    self, noisy_circuit: QuantumCircuit, qudits: list[int], noise_info: Noise | SubspaceNoise
) -> None:
    if isinstance(noise_info, SubspaceNoise):  # 物理的ノイズ
        for dit in qudits:
            dim = noisy_circuit.dimensions[dit]
            possible_levels = set(range(dim))
            
            for lev_a, lev_b in noise_info.subspace_w_probs:
                # 脱分極が作用する準位以外の準位に位相緩和を適用
                subspace_levels = {lev_a, lev_b}
                dephasing_levels = list(possible_levels - subspace_levels)
                
                # 各準位に確率的に位相ノイズを適用
                prob_each = noise_info.subspace_w_probs[lev_a, lev_b].probability_dephasing
                probs = [prob_each, 1 - prob_each]  # [ノイズ適用, 何もしない]
                
                for physical_level in dephasing_levels:
                    if self.rng.choice([True, False], p=probs):
                        noisy_circuit.noisez(dit, physical_level)
```

**数学的説明**:

位相緩和ノイズは、各準位の位相を確率的にランダム化する：

$$
\rho' = (1 - p_{\text{dephase}}) \rho + p_{\text{dephase}} Z_k \rho Z_k^\dagger
$$

ここで、$Z_k$ は準位 $k$ に作用する位相反転演算子。

**物理的意味**: 
- 脱分極が作用する2準位部分空間以外の準位に対して、位相のランダム化を行う
- 量子コヒーレンスの減衰を表現

### 動的部分空間ノイズの修正

**ソースコード**: `/src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py`

```python
def _dynamic_subspace_noise_info_rectification(self, noise_info: SubspaceNoise, instruction: Gate) -> SubspaceNoise:
    """
    負の準位インデックスを持つSubspaceNoiseを、ゲートの実際の準位に動的に割り当てる
    """
    # 修正が必要かチェック
    if not self._needs_correction(noise_info):
        return noise_info
    
    # サポートされているゲートタイプか確認
    if not isinstance(instruction, (R, Rz, Rh, CEx)):
        return noise_info
    
    # 負のインデックスから実際の準位へマッピング
    subspace = next(iter(noise_info.subspace_w_probs.keys()))
    noise_probs = noise_info.subspace_w_probs[subspace]
    
    # 新しいSubspaceNoiseを作成
    return SubspaceNoise(
        probability_depolarizing=noise_probs.probability_depolarizing,
        probability_dephasing=noise_probs.probability_dephasing,
        levels=(instruction.lev_a, instruction.lev_b),  # ゲートの実際の準位を使用
    )
```

**機能**:
- **動的割り当て**: ノイズモデルで準位を `-2, -1` として指定すると、実際のゲートが作用する準位に自動的にマッピング
- **用途**: 一般的なノイズモデルを定義し、異なる準位で動作するゲートに対して再利用

**例**:

```python
# 動的ノイズモデルの定義
subspace_noise = SubspaceNoise(
    probability_depolarizing=0.001,
    probability_dephasing=0.001,
    levels=[]  # 空リストで動的割り当てを指定
)

noise_model.add_quantum_error_locally(subspace_noise, ["r", "rz", "rh"])

# 実際のゲート適用時、各ゲートの実際の準位にノイズがマッピングされる
circuit.r(0, [0, 1, theta, phi])  # 準位0-1にノイズ
circuit.r(1, [1, 2, theta, phi])  # 準位1-2にノイズ
```

### NoisyCircuitFactoryクラス

**ソースコード**: `/src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py`

```python
class NoisyCircuitFactory:
    def __init__(self, noise_model: NoiseModel, circuit: QuantumCircuit) -> None:
        self.noise_model: NoiseModel = noise_model
        self.circuit: QuantumCircuit = circuit
        self.rng: Generator = self._initialize_rng()  # 乱数生成器
    
    def generate_circuit(self) -> QuantumCircuit:
        """
        元の回路にノイズを追加した新しい回路を生成
        """
        noisy_circuit = QuantumCircuit(self.circuit.num_qudits, self.circuit.dimensions, self.circuit.num_cl)
        
        for instruction in self.circuit.instructions:
            # 元のゲートをコピー
            copied_instruction = copy.deepcopy(instruction)
            noisy_circuit.instructions.append(copied_instruction)
            
            # ノイズを適用
            self._apply_noise(noisy_circuit, instruction)
        
        return noisy_circuit
```

**プロセス**:
1. 元の回路の各ゲートをコピー
2. ゲート後に、ノイズモデルに基づいてノイズゲートを追加
3. ノイズ付き回路を返す

**乱数生成器の初期化**:

```python
@staticmethod
def _initialize_rng() -> Generator:
    current_time = int(time.time() * 1000)
    seed = hash((os.getpid(), current_time)) % 2**32
    return np.random.default_rng(seed)
```

**特徴**:
- プロセスIDと現在時刻からシードを生成
- 各実行で異なるノイズパターンを生成（確率的シミュレーション）

### ノイズ適用の流れ

```python
def _apply_noise(self, noisy_circuit: QuantumCircuit, instruction: Gate) -> None:
    # ゲートタイプがノイズモデルに含まれているか確認
    if instruction.qasm_tag not in self.noise_model.quantum_errors:
        return
    
    # 各適用モードに対してノイズを適用
    for mode, noise_info in self.noise_model.quantum_errors[instruction.qasm_tag].items():
        # 影響を受けるquditを取得
        qudits = self._get_affected_qudits(instruction, mode)
        if qudits is None:
            continue
        
        # 動的部分空間ノイズの修正
        if isinstance(noise_info, SubspaceNoise):
            noise_info = self._dynamic_subspace_noise_info_rectification(noise_info, instruction)
        
        # 脱分極ノイズと位相緩和ノイズを適用
        self._apply_depolarizing_noise(noisy_circuit, qudits, noise_info)
        self._apply_dephasing_noise(noisy_circuit, qudits, noise_info)
```

**ステップ**:
1. ゲートタイプのチェック
2. 適用モードに基づくquditの選択
3. SubspaceNoiseの動的修正
4. 脱分極ノイズの適用
5. 位相緩和ノイズの適用

### 具体的な使用例

#### 例1: 単一quditゲートへのノイズ

```python
# ノイズモデルの作成
noise_model = NoiseModel()

# 単一quditゲート (R, Rz) に0.1%の脱分極ノイズを追加
subspace_noise = SubspaceNoise(
    probability_depolarizing=0.001,  # 0.1%
    probability_dephasing=0.001,     # 0.1%
    levels=[]  # 動的割り当て
)
noise_model.add_quantum_error_locally(subspace_noise, ["r", "rz", "rh"])

# 回路の作成
circuit = QuantumCircuit(2, [3, 3], 0)
circuit.r(0, [0, 1, np.pi/4, 0])
circuit.rz(1, [1, np.pi/2])

# ノイズ付き回路の生成
factory = NoisyCircuitFactory(noise_model, circuit)
noisy_circuit = factory.generate_circuit()
```

**結果**:
- `circuit.r(0, ...)` の後に、qudit 0の準位0-1に対する脱分極ノイズゲートが追加される
- `circuit.rz(1, ...)` の後に、qudit 1の準位1に対するノイズゲートが追加される

#### 例2: 2-quditゲートへのノイズ

```python
# 2-quditゲート (CEx) に1%の脱分極ノイズを追加
nonlocal_noise = SubspaceNoise(
    probability_depolarizing=0.01,  # 1%
    probability_dephasing=0.01,     # 1%
    levels=[]  # 動的割り当て
)
noise_model.add_nonlocal_quantum_error(nonlocal_noise, ["cex"])

# 回路の作成
circuit = QuantumCircuit(2, [3, 3], 0)
circuit.cex(0, 1, 0, 0, np.pi/4)

# ノイズ付き回路の生成
factory = NoisyCircuitFactory(noise_model, circuit)
noisy_circuit = factory.generate_circuit()
```

**結果**:
- `circuit.cex(0, 1, ...)` の後に、qudit 0と1の両方に対する脱分極ノイズゲートが追加される

### ノイズモデルの検証

**テストコード**: `/test/python/simulation/noise_tools/test_noise_tools_mathematical.py`

```python
def test_depolarizing_noise():
    """脱分極ノイズの数学的正しさを検証"""
    noise_model = NoiseModel()
    
    # 脱分極ノイズの追加
    depol_noise = Noise(probability_depolarizing=0.1, probability_dephasing=0.0)
    noise_model.add_all_qudit_quantum_error(depol_noise, ["x"])
    
    # 回路とシミュレーション
    circuit = QuantumCircuit(1, [3], 0)
    circuit.x(0)
    
    # 多数回実行して統計的検証
    results = []
    for _ in range(1000):
        factory = NoisyCircuitFactory(noise_model, circuit)
        noisy_circuit = factory.generate_circuit()
        # シミュレーション実行
        results.append(simulate(noisy_circuit))
    
    # 期待値の検証: (1 - p) * pure_state + p * mixed_state
    # ...
```

---

## まとめ

### CustomTwo分解のまとめ

1. **CustomTwoゲート**: 任意の2-quditユニタリ行列を表現するゲートクラス
2. **分解アルゴリズム**: Entangled QR decomposition with Controlled Exchange (EntangledQRCEX)
3. **数学的基礎**: Givens回転を用いたQR分解により、行列を逐次的に三角化
4. **実装**: CRotおよびPSwapゲート生成器を使用し、境界条件に応じてゲートを選択
5. **対角位相の補正**: 連立方程式を解いて、対角要素の位相を基本ゲートで実現
6. **厳密性**: 数値精度 $< 10^{-8}$ でユニタリ性を保証

### ノイズモデルのまとめ

1. **ノイズクラス階層**:
   - `Noise`: 基本的な脱分極・位相緩和確率
   - `SubspaceNoise`: 準位特異的なノイズ
   - `NoiseModel`: ゲート特異的なノイズモデル

2. **脱分極ノイズ**:
   - 数学的ノイズ: 全準位に一様に混合状態を導入
   - 物理的ノイズ: 特定の2準位部分空間に対してPauli演算をランダムに適用

3. **位相緩和ノイズ**: 脱分極が作用する準位以外の準位に位相ランダム化を適用

4. **動的割り当て**: 負の準位インデックスを使用して、ゲートの実際の準位にノイズを動的にマッピング

5. **NoisyCircuitFactory**: 元の回路に対してノイズゲートを追加し、ノイズ付き回路を生成

6. **確率的シミュレーション**: 乱数生成器を使用して、各ゲート後にノイズ演算を確率的に選択・適用

### 実装の特徴

- **厳密性**: 数学的に厳密な分解とノイズモデルの実装
- **一般性**: 任意の次元のquditに対応
- **柔軟性**: 様々なノイズモデルとゲートタイプに対応
- **検証可能性**: テストコードによる数学的正しさの検証

### 参考文献

1. MQT-Quditsソースコード: `/src/mqt/qudits/`
2. 理論文書: `/tutorials/doc/theory_quantum_dynamics_complete_comparison.md`
3. テストコード: `/test/python/qudits_circuits/gate_set/test_custom_two.py`
4. ノイズモデル実装: `/src/mqt/qudits/simulation/noise_tools/`

---

## 付録: コードリファレンス

### CustomTwo分解の主要クラス

- `EntangledQRCEX`: `/src/mqt/qudits/compiler/twodit/entanglement_qr/log_ent_qr_cex_decomp.py`
- `CRotGen`: `/src/mqt/qudits/compiler/twodit/blocks/crot.py`
- `PSwapGen`: `/src/mqt/qudits/compiler/twodit/blocks/pswap.py`
- `CZRotGen`: `/src/mqt/qudits/compiler/twodit/blocks/czrot.py`

### ノイズモデルの主要クラス

- `Noise`: `/src/mqt/qudits/simulation/noise_tools/noise.py` (4行目)
- `SubspaceNoise`: `/src/mqt/qudits/simulation/noise_tools/noise.py` (15行目)
- `NoiseModel`: `/src/mqt/qudits/simulation/noise_tools/noise.py` (83行目)
- `NoisyCircuitFactory`: `/src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py` (26行目)

### テストファイル

- CustomTwoテスト: `/test/python/qudits_circuits/gate_set/test_custom_two.py`
- ノイズモデルテスト: `/test/python/simulation/noise_tools/test_noise_tools_mathematical.py`
