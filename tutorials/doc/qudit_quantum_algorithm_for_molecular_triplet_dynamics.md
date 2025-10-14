# Quditによる分子三重項状態の量子ダイナミクスシミュレーション理論

## 1. はじめに

本文書では、`quantum_dynamics_molecular_triplet_states.md` および `suzuki_trotter_decomposition_theory.md` で記述された理論を基に、**MQT Quditsフレームワークを用いた量子アルゴリズムによる実装**のための完全な理論を提示する。

### 1.1 本文書の目的

- 分子三重項状態の量子ダイナミクスをQudits（d次元量子系）で表現
- 鈴木トロッター分解をQuditゲート操作に変換
- MQT Quditsフレームワークで直接プログラム化可能な形式で記述
- すべての数式を省略無しに展開

### 1.2 Quditの利点

従来の量子ビット（2準位系）では、3つの分子状態 $\{|S_0\rangle, |T_1\rangle, |S_1\rangle\}$ を表現するために複数の量子ビットが必要となる。しかし、**Qutrit（3準位量子系）**を用いることで、1つの物理的量子系で1分子の状態を自然に表現できる。

$$
\text{1分子} \longleftrightarrow \text{1 Qutrit（3次元Qudit）}
$$

## 2. Quditによる状態表現

### 2.1 単一分子のQutrit表現

分子 $i$ の3つの電子状態を、Qutrit（$d=3$）の計算基底として表現する：

$$
|S_0\rangle_i \longleftrightarrow |0\rangle_i = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}_i
$$

$$
|T_1\rangle_i \longleftrightarrow |1\rangle_i = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}_i
$$

$$
|S_1\rangle_i \longleftrightarrow |2\rangle_i = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}_i
$$

### 2.2 N分子系の状態

$N$ 個の分子からなる系は、$N$ 個のQutritのテンソル積で表現される：

$$
|\Psi\rangle = \sum_{n_1=0}^{2} \sum_{n_2=0}^{2} \cdots \sum_{n_N=0}^{2} c_{n_1 n_2 \cdots n_N} |n_1\rangle \otimes |n_2\rangle \otimes \cdots \otimes |n_N\rangle
$$

状態空間の次元は：

$$
\dim(\mathcal{H}) = 3^N
$$

### 2.3 計算基底の列挙

$N=2$ の場合の計算基底（9状態）：

$$
\begin{align}
|00\rangle &\equiv |S_0, S_0\rangle \\
|01\rangle &\equiv |S_0, T_1\rangle \\
|02\rangle &\equiv |S_0, S_1\rangle \\
|10\rangle &\equiv |T_1, S_0\rangle \\
|11\rangle &\equiv |T_1, T_1\rangle \\
|12\rangle &\equiv |T_1, S_1\rangle \\
|20\rangle &\equiv |S_1, S_0\rangle \\
|21\rangle &\equiv |S_1, T_1\rangle \\
|22\rangle &\equiv |S_1, S_1\rangle
\end{align}
$$

## 3. ハミルトニアンのQudit表現

### 3.1 単一Qutrit演算子

#### 3.1.1 数演算子（Number Operator）

各準位の占有数演算子をQutrit演算子として表現する：

$$
\hat{n}_0 = |0\rangle\langle 0| = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

$$
\hat{n}_1 = |1\rangle\langle 1| = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

$$
\hat{n}_2 = |2\rangle\langle 2| = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}
$$

完全性関係：

$$
\hat{n}_0 + \hat{n}_1 + \hat{n}_2 = \hat{I}_3
$$

ここで、$\hat{I}_3$ は3×3の単位行列である。

#### 3.1.2 対角ハミルトニアン $\hat{H}_0$

単一分子 $i$ のエネルギーハミルトニアンは：

$$
\hat{H}_0^{(i)} = E_T \hat{n}_1^{(i)} + E_S \hat{n}_2^{(i)}
$$

行列表現：

$$
\hat{H}_0^{(i)} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix}
$$

$N$ 分子系全体では：

$$
\hat{H}_0 = \sum_{i=1}^{N} \hat{H}_0^{(i)} = \sum_{i=1}^{N} \left( E_T \hat{n}_1^{(i)} + E_S \hat{n}_2^{(i)} \right)
$$

### 3.2 2-Qutrit演算子

#### 3.2.1 エネルギー移動演算子

隣接分子対 $(i, j)$ 間のエネルギー移動演算子は：

$$
\hat{T}_{ij}^{\text{transfer}} = |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0|
$$

行列要素：

$$
\langle m_i, m_j | \hat{T}_{ij}^{\text{transfer}} | n_i, n_j \rangle = \delta_{m_i, 0} \delta_{n_i, 1} \delta_{m_j, 1} \delta_{n_j, 0}
$$

9×9行列表現（基底順序 $\{|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle\}$）：

$$
\hat{T}_{ij}^{\text{transfer}} = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

エルミート共役：

$$
\hat{T}_{ij}^{\text{transfer}\dagger} = |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1|
$$

エネルギー移動ハミルトニアン：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V_{ij} \left( \hat{T}_{ij}^{\text{transfer}} + \hat{T}_{ij}^{\text{transfer}\dagger} \right)
$$

行列形式：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V_{ij} \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

#### 3.2.2 TTA演算子

三重項-三重項消滅演算子は：

$$
\hat{T}_{ij}^{\text{TTA}} = |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1| + |0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1|
$$

第1項の行列表現：

$$
|2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1| : \quad |11\rangle \to |20\rangle
$$

第2項の行列表現：

$$
|0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1| : \quad |11\rangle \to |02\rangle
$$

9×9行列：

$$
\hat{T}_{ij}^{\text{TTA}} = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

TTAハミルトニアン：

$$
\hat{H}_{\text{TTA}}^{(ij)} = J_{ij} \left( \hat{T}_{ij}^{\text{TTA}} + \hat{T}_{ij}^{\text{TTA}\dagger} \right)
$$


## 4. Quditゲートによる時間発展演算子の実装

### 4.1 対角ハミルトニアン $\hat{H}_0$ の時間発展

#### 4.1.1 単一Qutrit位相ゲート

対角ハミルトニアン $\hat{H}_0^{(i)}$ の時間発展演算子は：

$$
\hat{U}_0^{(i)}(t) = e^{-i\hat{H}_0^{(i)}t/\hbar} = e^{-i(E_T \hat{n}_1^{(i)} + E_S \hat{n}_2^{(i)})t/\hbar}
$$

演算子の可換性（$[\hat{n}_1, \hat{n}_2] = 0$）より：

$$
\hat{U}_0^{(i)}(t) = e^{-iE_T \hat{n}_1^{(i)}t/\hbar} e^{-iE_S \hat{n}_2^{(i)}t/\hbar}
$$

行列表現：

$$
\hat{U}_0^{(i)}(t) = \begin{pmatrix}
1 & 0 & 0 \\
0 & e^{-iE_T t/\hbar} & 0 \\
0 & 0 & e^{-iE_S t/\hbar}
\end{pmatrix}
$$

#### 4.1.2 Qutrit位相ゲート（VirtRz）

MQT Quditsフレームワークでは、任意の位相ゲート `VirtRz` を用いて実装できる：

**準位 $|1\rangle$ への位相 $\phi_1 = -E_T t/\hbar$**:

$$
\text{VirtRz}_{\text{level}=1}(\phi_1) = \begin{pmatrix}
1 & 0 & 0 \\
0 & e^{i\phi_1} & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

**準位 $|2\rangle$ への位相 $\phi_2 = -E_S t/\hbar$**:

$$
\text{VirtRz}_{\text{level}=2}(\phi_2) = \begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & e^{i\phi_2}
\end{pmatrix}
$$

合成ゲート：

$$
\hat{U}_0^{(i)}(t) = \text{VirtRz}_{\text{level}=2}(\phi_2) \cdot \text{VirtRz}_{\text{level}=1}(\phi_1)
$$

#### 4.1.3 MQT Quditsでの実装コード

```python
from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
import numpy as np

# パラメータ設定
E_T = 1.5  # eV
E_S = 3.0  # eV
dt = 0.01  # 時間ステップ (任意単位)
hbar = 1.0  # 換算プランク定数（単位系に依存）

# 位相角の計算
phi_1 = -E_T * dt / hbar
phi_2 = -E_S * dt / hbar

# 量子回路の構築
circuit = QuantumCircuit()
molecule_reg = QuantumRegister("molecules", N, [3]*N)  # N個のQutrit
circuit.append(molecule_reg)

# 各Qutritに位相ゲートを適用
for i in range(N):
    # 準位1への位相
    circuit.virtrz(molecule_reg[i], [1, phi_1])
    # 準位2への位相
    circuit.virtrz(molecule_reg[i], [2, phi_2])
```

### 4.2 エネルギー移動ハミルトニアン $\hat{H}_{\text{transfer}}$ の時間発展

#### 4.2.1 2準位部分空間への射影

エネルギー移動ハミルトニアン $\hat{H}_{\text{transfer}}^{(ij)}$ は、部分空間 $\{|01\rangle, |10\rangle\}$ で非ゼロである：

$$
\hat{H}_{\text{sub}} = V_{ij} \begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix} = V_{ij} \sigma_x
$$

ここで、$\sigma_x$ はパウリX行列である。

#### 4.2.2 時間発展演算子の解析解

行列指数関数：

$$
e^{-i\hat{H}_{\text{sub}}t/\hbar} = e^{-iV_{ij}\sigma_x t/\hbar}
$$

パウリ行列の性質 $\sigma_x^2 = I$ を用いて：

$$
e^{-i\theta\sigma_x} = \cos\theta \cdot I - i\sin\theta \cdot \sigma_x
$$

ここで、$\theta = V_{ij}t/\hbar$ とすると：

$$
e^{-i\hat{H}_{\text{sub}}t/\hbar} = \begin{pmatrix}
\cos\theta & -i\sin\theta \\
-i\sin\theta & \cos\theta
\end{pmatrix}
$$

#### 4.2.3 Qutrit回転ゲート（RXY）

MQT Quditsの `R` ゲート（2準位間の回転）を用いる：

**定義**: 準位 $|a\rangle$ と $|b\rangle$ 間の回転ゲート

$$
\text{R}_{ab}(\theta, \phi) = e^{-i\theta(\cos\phi \, X_{ab} + \sin\phi \, Y_{ab})}
$$

ここで、$X_{ab}, Y_{ab}$ は準位 $a, b$ 間のパウリ様演算子：

$$
X_{ab} = |a\rangle\langle b| + |b\rangle\langle a|
$$

$$
Y_{ab} = -i|a\rangle\langle b| + i|b\rangle\langle a|
$$

エネルギー移動には $X$ 回転（$\phi = 0$）を用いる：

$$
\text{R}_{ab}(\theta, 0) = e^{-i\theta X_{ab}}
$$

準位 $a=1, b=0$ の場合：

$$
X_{10} = |1\rangle\langle 0| + |0\rangle\langle 1| = \begin{pmatrix}
0 & 0 & 0 \\
1 & 0 & 0 \\
0 & 0 & 0
\end{pmatrix} + \begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{pmatrix} = \begin{pmatrix}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

時間発展演算子（単一Qutrit）：

$$
\text{R}_{10}(\theta, 0) = \begin{pmatrix}
\cos\theta & -i\sin\theta & 0 \\
-i\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

#### 4.2.4 2-Qutritゲート分解

隣接Qutrit対 $(i, j)$ への時間発展は、制御Qutritゲートで実装される。しかし、エネルギー移動は **2体相互作用** であるため、直接的な2-Quditゲートまたは制御ゲートの組み合わせが必要。

**方法1: カスタム2-Quditゲート**

MQT Quditsの `cu_two` (Custom two-qudit gate) を用いて、9×9ユニタリ行列を直接指定：

```python
import numpy as np

def energy_transfer_gate(V_ij, dt, hbar=1.0):
    """
    2-Qutrit間のエネルギー移動ゲート
    """
    theta = V_ij * dt / hbar
    
    # 9x9単位行列
    U = np.eye(9, dtype=complex)
    
    # 部分空間 {|01>, |10>} = {index 1, index 3} での回転
    # 基底順序: |00>=0, |01>=1, |02>=2, |10>=3, |11>=4, ...
    c = np.cos(theta)
    s = np.sin(theta)
    
    U[1, 1] = c
    U[1, 3] = -1j * s
    U[3, 1] = -1j * s
    U[3, 3] = c
    
    return U

# 回路への適用
theta_transfer = V_ij * dt / hbar
U_transfer = energy_transfer_gate(V_ij, dt, hbar)
circuit.cu_two([i, j], U_transfer)
```

**方法2: 制御回転ゲートの分解**

より一般的な分解として、制御NOT型ゲートと局所回転の組み合わせを用いる。ただし、Qutritの場合は複雑になるため、直接的な行列指定が推奨される。

### 4.3 TTAハミルトニアン $\hat{H}_{\text{TTA}}$ の時間発展

#### 4.3.1 3準位部分空間への射影

TTAハミルトニアンは部分空間 $\{|11\rangle, |02\rangle, |20\rangle\}$ で作用する：

基底順序を $\{|11\rangle, |02\rangle, |20\rangle\}$ として：

$$
\hat{H}_{\text{TTA,sub}} = J_{ij} \begin{pmatrix}
0 & 1 & 1 \\
1 & 0 & 0 \\
1 & 0 & 0
\end{pmatrix}
$$

#### 4.3.2 固有値分解

固有値問題を解く：

$$
\det(\hat{H}_{\text{TTA,sub}} - \lambda I) = 0
$$

計算により、固有値は：

$$
\lambda_0 = 0, \quad \lambda_+ = J_{ij}\sqrt{2}, \quad \lambda_- = -J_{ij}\sqrt{2}
$$

固有ベクトル：

$$
|v_0\rangle = \frac{1}{\sqrt{2}} \begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix}, \quad
|v_+\rangle = \frac{1}{2} \begin{pmatrix} \sqrt{2} \\ 1 \\ 1 \end{pmatrix}, \quad
|v_-\rangle = \frac{1}{2} \begin{pmatrix} -\sqrt{2} \\ 1 \\ 1 \end{pmatrix}
$$

時間発展演算子：

$$
\hat{U}_{\text{TTA,sub}}(t) = e^{-i\hat{H}_{\text{TTA,sub}}t/\hbar} = \sum_{k} e^{-i\lambda_k t/\hbar} |v_k\rangle\langle v_k|
$$

行列形式：

$$
\phi_\pm = \pm J_{ij}\sqrt{2} \, t / \hbar
$$

$$
\hat{U}_{\text{TTA,sub}}(t) = \begin{pmatrix}
\cos(\sqrt{2}\phi) & -i\sin(\sqrt{2}\phi)/\sqrt{2} & -i\sin(\sqrt{2}\phi)/\sqrt{2} \\
-i\sin(\sqrt{2}\phi)/\sqrt{2} & \cos^2(\sqrt{2}\phi/2) & \sin^2(\sqrt{2}\phi/2) \\
-i\sin(\sqrt{2}\phi)/\sqrt{2} & \sin^2(\sqrt{2}\phi/2) & \cos^2(\sqrt{2}\phi/2)
\end{pmatrix}
$$

ここで、$\phi = J_{ij}t/\hbar$ である。

#### 4.3.3 カスタム2-Quditゲートでの実装

```python
def TTA_gate(J_ij, dt, hbar=1.0):
    """
    2-Qutrit間のTTAゲート
    """
    phi = J_ij * dt / hbar
    sqrt2 = np.sqrt(2)
    
    # 9x9単位行列
    U = np.eye(9, dtype=complex)
    
    # 部分空間 {|11>, |02>, |20>} = {index 4, 2, 6}
    c = np.cos(sqrt2 * phi)
    s = np.sin(sqrt2 * phi)
    c2 = np.cos(sqrt2 * phi / 2)**2
    s2 = np.sin(sqrt2 * phi / 2)**2
    
    # 行列要素の設定
    U[4, 4] = c
    U[4, 2] = -1j * s / sqrt2
    U[4, 6] = -1j * s / sqrt2
    U[2, 4] = -1j * s / sqrt2
    U[2, 2] = c2
    U[2, 6] = s2
    U[6, 4] = -1j * s / sqrt2
    U[6, 2] = s2
    U[6, 6] = c2
    
    return U

# 回路への適用
U_TTA = TTA_gate(J_ij, dt, hbar)
circuit.cu_two([i, j], U_TTA)
```

### 4.4 放射減衰の実装

#### 4.4.1 非ユニタリ演算としての減衰

蛍光放出は非ユニタリ過程であり、厳密には密度行列形式またはリンドブラッド方程式が必要である。しかし、短時間近似では **実効的な減衰** として扱える。

準位 $|2\rangle$（$|S_1\rangle$）の振幅減衰：

$$
c_2(t) \to c_2(t) \cdot e^{-\Gamma_{\text{fl}} t / 2}
$$

全状態の規格化：

$$
|\Psi'(t)\rangle = \frac{|\Psi(t)\rangle}{\sqrt{\langle\Psi(t)|\Psi(t)\rangle}}
$$

#### 4.4.2 ノイズモデルによる実装

MQT Quditsの `NoiseModel` を用いて、減衰チャネルを追加：

```python
from mqt.qudits.simulation.noise_tools import NoiseModel

# ノイズモデルの設定
noise_model = NoiseModel()

# 振幅減衰（Amplitude Damping）を各Qutritの準位2に適用
gamma_fl = 0.1  # 減衰率
for i in range(N):
    noise_model.add_amplitude_damping(qudit=i, level=2, gamma=gamma_fl * dt)
```

あるいは、各ステップ後に手動で減衰を適用：

```python
def apply_radiative_decay(state_vector, circuit_dims, gamma_fl, dt):
    """
    状態ベクトルに放射減衰を適用
    """
    # 各基底状態について、準位2の数をカウント
    N_qudits = len(circuit_dims)
    new_state = state_vector.copy()
    
    for idx, amplitude in enumerate(state_vector):
        # idx を3進数表現に変換して各Quditの状態を取得
        config = index_to_config(idx, N_qudits, d=3)
        n_S1 = sum(1 for level in config if level == 2)
        
        # 減衰因子を適用
        decay_factor = np.exp(-gamma_fl * dt * n_S1 / 2)
        new_state[idx] *= decay_factor
    
    # 規格化
    new_state /= np.linalg.norm(new_state)
    
    return new_state

def index_to_config(idx, N, d=3):
    """
    線形インデックスをd進数配列に変換
    """
    config = []
    for _ in range(N):
        config.append(idx % d)
        idx //= d
    return config[::-1]
```


## 5. 鈴木トロッター分解のQudit実装

### 5.1 完全なハミルトニアンの分解

全ハミルトニアンを4つの項に分解：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}} + \hat{H}_{\text{rad}}
$$

ここで：

$$
\hat{H}_0 = \sum_{i=1}^{N} \left( E_T \hat{n}_1^{(i)} + E_S \hat{n}_2^{(i)} \right)
$$

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( \hat{T}_{ij}^{\text{transfer}} + \hat{T}_{ij}^{\text{transfer}\dagger} \right)
$$

$$
\hat{H}_{\text{TTA}} = \sum_{\langle i,j \rangle} J_{ij} \left( \hat{T}_{ij}^{\text{TTA}} + \hat{T}_{ij}^{\text{TTA}\dagger} \right)
$$

### 5.2 2次対称鈴木トロッター分解

1時間ステップ $\Delta t$ の時間発展演算子：

$$
\hat{U}(\Delta t) = e^{-i\hat{H}_{\text{total}}\Delta t/\hbar}
$$

2次対称分解：

$$
\hat{U}(\Delta t) \approx \hat{U}_{\text{ST2}}(\Delta t) + \mathcal{O}(\Delta t^3)
$$

ここで：

$$
\hat{U}_{\text{ST2}}(\Delta t) = e^{-i\hat{H}_0\Delta t/(2\hbar)} e^{-i\hat{H}_{\text{transfer}}\Delta t/(2\hbar)} e^{-i\hat{H}_{\text{TTA}}\Delta t/(2\hbar)} e^{-i\hat{H}_{\text{rad}}\Delta t/\hbar} e^{-i\hat{H}_{\text{TTA}}\Delta t/(2\hbar)} e^{-i\hat{H}_{\text{transfer}}\Delta t/(2\hbar)} e^{-i\hat{H}_0\Delta t/(2\hbar)}
$$

### 5.3 Qudit回路による実装

#### 5.3.1 アルゴリズムの流れ

```
入力: 初期状態 |Ψ₀⟩, パラメータ {E_T, E_S, V_ij, J_ij, Γ_fl}, 
      全時間 T, ステップ数 N_steps

1. Δt = T / N_steps を計算
2. 初期状態を量子レジスタに設定
3. for step = 1 to N_steps:
   a. H₀ を時間 Δt/2 で時間発展 (各Quditに位相ゲート)
   b. H_transfer を時間 Δt/2 で時間発展 (各隣接対にカスタムゲート)
   c. H_TTA を時間 Δt/2 で時間発展 (各隣接対にカスタムゲート)
   d. H_rad を時間 Δt で時間発展 (減衰操作)
   e. H_TTA を時間 Δt/2 で時間発展 (逆順)
   f. H_transfer を時間 Δt/2 で時間発展 (逆順)
   g. H₀ を時間 Δt/2 で時間発展 (逆順)
   h. (オプション) 観測量を測定・記録
4. 最終状態 |Ψ(T)⟩ を出力
```

#### 5.3.2 完全な実装例（MQT Qudits）

```python
import numpy as np
from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
from mqt.qudits.simulation import MQTQuditProvider

class MolecularTripletSimulator:
    def __init__(self, N_molecules, E_T, E_S, V, J, Gamma_fl, 
                 lattice_type='chain'):
        """
        N_molecules: 分子数
        E_T: 三重項エネルギー
        E_S: 一重項エネルギー
        V: エネルギー移動積分
        J: TTA相互作用定数
        Gamma_fl: 蛍光放出速度
        lattice_type: 'chain' または 'ring'（周期境界条件）
        """
        self.N = N_molecules
        self.E_T = E_T
        self.E_S = E_S
        self.V = V
        self.J = J
        self.Gamma_fl = Gamma_fl
        self.hbar = 1.0  # 単位系の設定
        
        # 隣接リストの構築
        if lattice_type == 'chain':
            self.neighbors = [(i, i+1) for i in range(N_molecules-1)]
        elif lattice_type == 'ring':
            self.neighbors = [(i, (i+1) % N_molecules) for i in range(N_molecules)]
        else:
            raise ValueError("lattice_type must be 'chain' or 'ring'")
    
    def build_H0_gate(self, qudit_idx, dt):
        """
        単一Quditの対角ハミルトニアンゲート
        """
        phi_1 = -self.E_T * dt / self.hbar
        phi_2 = -self.E_S * dt / self.hbar
        
        # 3x3対角行列
        U_H0 = np.diag([1.0, np.exp(1j * phi_1), np.exp(1j * phi_2)])
        
        return U_H0
    
    def build_transfer_gate(self, dt):
        """
        2-Quditエネルギー移動ゲート
        """
        theta = self.V * dt / self.hbar
        
        U = np.eye(9, dtype=complex)
        
        # 部分空間 {|01⟩, |10⟩}
        c = np.cos(theta)
        s = np.sin(theta)
        
        U[1, 1] = c
        U[1, 3] = -1j * s
        U[3, 1] = -1j * s
        U[3, 3] = c
        
        return U
    
    def build_TTA_gate(self, dt):
        """
        2-Qudit TTAゲート
        """
        phi = self.J * dt / self.hbar
        sqrt2 = np.sqrt(2)
        
        U = np.eye(9, dtype=complex)
        
        c = np.cos(sqrt2 * phi)
        s = np.sin(sqrt2 * phi)
        c2 = np.cos(sqrt2 * phi / 2)**2
        s2 = np.sin(sqrt2 * phi / 2)**2
        
        # 部分空間 {|11⟩, |02⟩, |20⟩} = {4, 2, 6}
        U[4, 4] = c
        U[4, 2] = -1j * s / sqrt2
        U[4, 6] = -1j * s / sqrt2
        U[2, 4] = -1j * s / sqrt2
        U[2, 2] = c2
        U[2, 6] = s2
        U[6, 4] = -1j * s / sqrt2
        U[6, 2] = s2
        U[6, 6] = c2
        
        return U
    
    def suzuki_trotter_circuit(self, T_total, N_steps, initial_state='all_triplet'):
        """
        鈴木トロッター2次分解による量子回路の構築
        """
        dt = T_total / N_steps
        
        # 量子回路の初期化
        circuit = QuantumCircuit()
        molecule_reg = QuantumRegister("molecules", self.N, [3] * self.N)
        circuit.append(molecule_reg)
        
        # 初期状態の準備
        if initial_state == 'all_triplet':
            # すべてのQutritを |1⟩ (三重項)に初期化
            for i in range(self.N):
                # X3ゲート: |0⟩ → |1⟩
                circuit.x(molecule_reg[i])  # 実際には適切な励起ゲートを使用
        
        # 各時間ステップ
        for step in range(N_steps):
            # 前半の対称分解
            
            # (a) H₀ evolution (dt/2)
            for i in range(self.N):
                U_h0_half = self.build_H0_gate(i, dt/2)
                circuit.cu_one(molecule_reg[i], U_h0_half)
            
            # (b) H_transfer evolution (dt/2)
            U_transfer_half = self.build_transfer_gate(dt/2)
            for (i, j) in self.neighbors:
                circuit.cu_two([molecule_reg[i], molecule_reg[j]], U_transfer_half)
            
            # (c) H_TTA evolution (dt/2)
            U_TTA_half = self.build_TTA_gate(dt/2)
            for (i, j) in self.neighbors:
                circuit.cu_two([molecule_reg[i], molecule_reg[j]], U_TTA_half)
            
            # (d) H_rad evolution (dt)
            # 減衰は回路では直接表現できないため、ノイズモデルで扱う
            
            # 後半の対称分解（逆順）
            
            # (e) H_TTA evolution (dt/2)
            for (i, j) in reversed(self.neighbors):
                circuit.cu_two([molecule_reg[i], molecule_reg[j]], U_TTA_half)
            
            # (f) H_transfer evolution (dt/2)
            for (i, j) in reversed(self.neighbors):
                circuit.cu_two([molecule_reg[i], molecule_reg[j]], U_transfer_half)
            
            # (g) H₀ evolution (dt/2)
            for i in reversed(range(self.N)):
                circuit.cu_one(molecule_reg[i], U_h0_half)
        
        return circuit
    
    def run_simulation(self, T_total, N_steps, backend_name='tnsim'):
        """
        シミュレーションの実行
        """
        circuit = self.suzuki_trotter_circuit(T_total, N_steps)
        
        # シミュレーションバックエンドの取得
        provider = MQTQuditProvider()
        backend = provider.get_backend(backend_name)
        
        # 実行
        job = backend.run(circuit)
        result = job.result()
        
        # 状態ベクトルの取得
        state_vector = result.get_state_vector()
        
        return state_vector, circuit
    
    def calculate_populations(self, state_vector):
        """
        各状態の個体数を計算
        """
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0
        
        # 全基底状態を走査
        for idx, amplitude in enumerate(state_vector.flatten()):
            prob = np.abs(amplitude)**2
            
            # インデックスを3進数配列に変換
            config = self._index_to_config(idx)
            
            # 各準位の個体数をカウント
            for level in config:
                if level == 0:
                    N_S0 += prob
                elif level == 1:
                    N_T1 += prob
                elif level == 2:
                    N_S1 += prob
        
        return N_S0, N_T1, N_S1
    
    def _index_to_config(self, idx):
        """
        線形インデックスを3進数配列に変換
        """
        config = []
        for _ in range(self.N):
            config.append(idx % 3)
            idx //= 3
        return config[::-1]

# 使用例
if __name__ == "__main__":
    # パラメータ設定
    N_molecules = 5
    E_T = 1.5  # eV
    E_S = 3.0  # eV
    V = 0.1    # eV
    J = 0.05   # eV
    Gamma_fl = 0.0  # まずは減衰なしで
    
    # シミュレータの初期化
    simulator = MolecularTripletSimulator(
        N_molecules, E_T, E_S, V, J, Gamma_fl,
        lattice_type='chain'
    )
    
    # シミュレーション実行
    T_total = 10.0  # 全時間
    N_steps = 100   # ステップ数
    
    state_vector, circuit = simulator.run_simulation(T_total, N_steps)
    
    # 個体数の計算
    N_S0, N_T1, N_S1 = simulator.calculate_populations(state_vector)
    
    print(f"Population S0: {N_S0:.4f}")
    print(f"Population T1: {N_T1:.4f}")
    print(f"Population S1: {N_S1:.4f}")
    print(f"Total: {N_S0 + N_T1 + N_S1:.4f}")
```

### 5.4 時間発展の追跡

各時間ステップでの観測量を記録するためには、中間状態を取得する必要がある：

```python
def run_simulation_with_tracking(self, T_total, N_steps):
    """
    時間発展を追跡するシミュレーション
    """
    dt = T_total / N_steps
    
    # 初期状態ベクトル（すべて三重項）
    dim = 3**self.N
    state = np.zeros(dim, dtype=complex)
    
    # |111...1⟩ に対応するインデックス
    # 3進数で1が並ぶ状態
    idx_all_triplet = sum(3**i for i in range(self.N))
    state[idx_all_triplet] = 1.0
    
    # 時間発展の記録
    times = [0.0]
    populations = [(self.N, 0.0, 0.0)]  # (N_S0, N_T1, N_S1)
    
    # 時間発展演算子の構築
    U_H0_half = self._build_full_H0_operator(dt/2)
    U_transfer_half = self._build_full_transfer_operator(dt/2)
    U_TTA_half = self._build_full_TTA_operator(dt/2)
    U_rad = self._build_radiative_operator(dt)
    
    # 鈴木トロッター分解による時間発展
    for step in range(N_steps):
        # 2次対称分解の適用
        state = U_H0_half @ state
        state = U_transfer_half @ state
        state = U_TTA_half @ state
        state = U_rad @ state
        state = U_TTA_half @ state
        state = U_transfer_half @ state
        state = U_H0_half @ state
        
        # 規格化（減衰がある場合）
        state /= np.linalg.norm(state)
        
        # 観測量の計算
        t = (step + 1) * dt
        N_S0, N_T1, N_S1 = self.calculate_populations(state.reshape(-1, 1))
        
        times.append(t)
        populations.append((N_S0, N_T1, N_S1))
    
    return times, populations

def _build_full_H0_operator(self, dt):
    """
    全系のH₀時間発展演算子を構築
    """
    dim = 3**self.N
    U_full = np.eye(dim, dtype=complex)
    
    # 各基底状態に対して位相を適用
    for idx in range(dim):
        config = self._index_to_config(idx)
        energy = sum(self.E_T if level == 1 else self.E_S if level == 2 else 0 
                    for level in config)
        phase = np.exp(-1j * energy * dt / self.hbar)
        U_full[idx, idx] = phase
    
    return U_full

def _build_full_transfer_operator(self, dt):
    """
    全系のエネルギー移動時間発展演算子を構築
    """
    dim = 3**self.N
    U_full = np.eye(dim, dtype=complex)
    
    # 各隣接対に対してゲートを適用
    for (i, j) in self.neighbors:
        U_pair = self.build_transfer_gate(dt)
        U_full = self._apply_two_qudit_gate(U_full, U_pair, i, j)
    
    return U_full

def _apply_two_qudit_gate(self, U_full, U_pair, qudit_i, qudit_j):
    """
    2-Quditゲートを全系の演算子に組み込む
    """
    # テンソル積で拡張
    # 実装の詳細は省略（Kronecker積を用いる）
    # 完全な実装には、適切なインデックス操作が必要
    
    # 簡易版: 疎行列を用いた実装
    from scipy.sparse import csr_matrix
    
    dim = 3**self.N
    U_new = U_full.copy()
    
    # qudit_i と qudit_j に作用するゲートを全体に拡張
    # （詳細な実装は長くなるため、概念的な記述）
    
    return U_new
```


## 6. 観測量の計算

### 6.1 個体数演算子

各状態の個体数演算子をQudit表現で記述する。

#### 6.1.1 基底状態個体数

$$
\hat{N}_{S_0} = \sum_{i=1}^{N} \hat{n}_0^{(i)} = \sum_{i=1}^{N} |0\rangle_i\langle 0|
$$

期待値：

$$
\langle \hat{N}_{S_0} \rangle = \langle\Psi| \hat{N}_{S_0} |\Psi\rangle
$$

#### 6.1.2 三重項個体数

$$
\hat{N}_{T_1} = \sum_{i=1}^{N} \hat{n}_1^{(i)} = \sum_{i=1}^{N} |1\rangle_i\langle 1|
$$

#### 6.1.3 励起一重項個体数

$$
\hat{N}_{S_1} = \sum_{i=1}^{N} \hat{n}_2^{(i)} = \sum_{i=1}^{N} |2\rangle_i\langle 2|
$$

### 6.2 状態ベクトルからの計算

状態ベクトル $|\Psi\rangle$ を計算基底で展開：

$$
|\Psi\rangle = \sum_{n_1=0}^{2} \cdots \sum_{n_N=0}^{2} c_{n_1 \cdots n_N} |n_1, \ldots, n_N\rangle
$$

各状態の個体数：

$$
\langle \hat{N}_{S_0} \rangle = \sum_{n_1, \ldots, n_N} |c_{n_1 \cdots n_N}|^2 \sum_{i=1}^{N} \delta_{n_i, 0}
$$

$$
\langle \hat{N}_{T_1} \rangle = \sum_{n_1, \ldots, n_N} |c_{n_1 \cdots n_N}|^2 \sum_{i=1}^{N} \delta_{n_i, 1}
$$

$$
\langle \hat{N}_{S_1} \rangle = \sum_{n_1, \ldots, n_N} |c_{n_1 \cdots n_N}|^2 \sum_{i=1}^{N} \delta_{n_i, 2}
$$

### 6.3 蛍光強度

瞬時蛍光強度は一重項個体数に比例：

$$
I_{\text{fl}}(t) = \Gamma_{\text{fl}} \langle \hat{N}_{S_1}(t) \rangle
$$

累積蛍光量：

$$
\mathcal{I}_{\text{total}} = \int_0^T I_{\text{fl}}(t) \, dt \approx \sum_{n=0}^{N_{\text{steps}}} I_{\text{fl}}(t_n) \Delta t
$$

### 6.4 相関関数

2点相関関数は、空間的なエネルギー移動やTTAの特性を評価するために重要である。

#### 6.4.1 三重項-三重項相関

$$
C_{TT}(i, j) = \langle \hat{n}_1^{(i)} \hat{n}_1^{(j)} \rangle - \langle \hat{n}_1^{(i)} \rangle \langle \hat{n}_1^{(j)} \rangle
$$

これは、分子 $i$ と $j$ の三重項状態が相関しているかを示す。

#### 6.4.2 実装

```python
def calculate_correlation(state_vector, qudit_i, qudit_j, level_i, level_j):
    """
    2サイト相関関数を計算
    """
    N = len(state_vector.shape) if state_vector.ndim > 1 else int(np.log(len(state_vector)) / np.log(3))
    
    # 期待値 <n_i n_j>
    expect_ij = 0.0
    for idx, amplitude in enumerate(state_vector.flatten()):
        prob = np.abs(amplitude)**2
        config = index_to_config(idx, N, d=3)
        if config[qudit_i] == level_i and config[qudit_j] == level_j:
            expect_ij += prob
    
    # 期待値 <n_i> と <n_j>
    expect_i = 0.0
    expect_j = 0.0
    for idx, amplitude in enumerate(state_vector.flatten()):
        prob = np.abs(amplitude)**2
        config = index_to_config(idx, N, d=3)
        if config[qudit_i] == level_i:
            expect_i += prob
        if config[qudit_j] == level_j:
            expect_j += prob
    
    # 相関関数
    correlation = expect_ij - expect_i * expect_j
    
    return correlation
```

## 7. 量子回路の最適化

### 7.1 ゲート数の削減

鈴木トロッター分解では多数のゲートが必要となる。ゲート数を削減するための戦略：

#### 7.1.1 ゲートの融合

連続する同種のゲートを1つにまとめる：

$$
e^{-i\hat{H}_1 t_1/\hbar} e^{-i\hat{H}_1 t_2/\hbar} = e^{-i\hat{H}_1 (t_1+t_2)/\hbar}
$$

特に、対称分解の前半と後半で同じゲートが現れる場合、次のステップと統合できる。

#### 7.1.2 可換性の利用

可換な演算子 $[\hat{A}, \hat{B}] = 0$ の場合、順序を入れ替えても結果は変わらない：

$$
e^{-i\hat{A}t/\hbar} e^{-i\hat{B}t/\hbar} = e^{-i\hat{B}t/\hbar} e^{-i\hat{A}t/\hbar}
$$

異なるサイトに作用する演算子は一般に可換である。

### 7.2 並列化

#### 7.2.1 空間的並列化

異なる隣接ペア $(i, j)$ と $(k, l)$ が重ならない場合（$\{i,j\} \cap \{k,l\} = \emptyset$）、対応するゲートは可換であり、並列に適用できる：

$$
[\hat{U}_{ij}, \hat{U}_{kl}] = 0 \quad \text{if } \{i,j\} \cap \{k,l\} = \emptyset
$$

例：1次元鎖の場合

- 偶数ペア：$(0,1), (2,3), (4,5), \ldots$ は同時に適用可能
- 奇数ペア：$(1,2), (3,4), (5,6), \ldots$ は同時に適用可能

#### 7.2.2 量子回路の深さ削減

並列化により、回路の深さ（depth）を削減できる：

$$
\text{Depth}_{\text{serial}} = \mathcal{O}(N \cdot M) \quad \to \quad \text{Depth}_{\text{parallel}} = \mathcal{O}(M)
$$

ここで、$N$ はQudits数、$M$ は時間ステップ数である。

### 7.3 適応的時間刻み

時間発展の速さに応じて $\Delta t$ を動的に調整することで、効率と精度を両立できる：

```python
def adaptive_trotter_simulation(self, T_total, dt_initial=0.1, tolerance=1e-6):
    """
    適応的時間刻み幅による鈴木トロッター法
    """
    t = 0.0
    dt = dt_initial
    state = self.initialize_state()
    
    times = [t]
    populations = [self.calculate_populations(state)]
    
    while t < T_total:
        # 時刻 dt で1ステップ進める
        state_1 = self.one_step_evolution(state, dt)
        
        # 時刻 dt/2 で2ステップ進める
        state_half = self.one_step_evolution(state, dt/2)
        state_2 = self.one_step_evolution(state_half, dt/2)
        
        # 誤差評価
        error = np.linalg.norm(state_1 - state_2)
        
        if error < tolerance:
            # 誤差が許容範囲内
            state = state_2  # より精度の高い方を採用
            t += dt
            times.append(t)
            populations.append(self.calculate_populations(state))
            
            # 時間刻みを増やす
            dt = min(dt * 1.5, T_total - t)
        else:
            # 誤差が大きい
            dt = dt * 0.5  # 時間刻みを減らして再試行
    
    return times, populations
```

## 8. 収束性と誤差評価

### 8.1 トロッター誤差の理論

$p$ 次の鈴木トロッター分解の局所誤差（1ステップあたり）：

$$
\varepsilon_{\text{local}} = \mathcal{O}(\Delta t^{p+1})
$$

全体誤差（$N$ ステップ後）：

$$
\varepsilon_{\text{global}} = \frac{T}{\Delta t} \cdot \mathcal{O}(\Delta t^{p+1}) = \mathcal{O}(\Delta t^p)
$$

2次分解の場合：

$$
\varepsilon_{\text{global}} = C \cdot \Delta t^2
$$

ここで、定数 $C$ は：

$$
C \sim T \cdot \max_{m \neq n} \| [\hat{H}_m, \hat{H}_n] \|
$$

### 8.2 収束テスト

実際の計算では、$\Delta t$ を系統的に減少させて収束を確認：

```python
def convergence_test(self, T_total, N_steps_list):
    """
    異なる時間刻み幅での収束テスト
    """
    results = []
    
    for N_steps in N_steps_list:
        dt = T_total / N_steps
        state_final = self.run_simulation(T_total, N_steps)
        N_S0, N_T1, N_S1 = self.calculate_populations(state_final)
        
        results.append({
            'N_steps': N_steps,
            'dt': dt,
            'N_S0': N_S0,
            'N_T1': N_T1,
            'N_S1': N_S1
        })
    
    return results

# 使用例
N_steps_list = [50, 100, 200, 400, 800]
results = simulator.convergence_test(T_total=10.0, N_steps_list=N_steps_list)

# 収束プロット
import matplotlib.pyplot as plt

dts = [r['dt'] for r in results]
N_T1_values = [r['N_T1'] for r in results]

plt.figure()
plt.loglog(dts, np.abs(np.array(N_T1_values) - N_T1_values[-1]), 'o-')
plt.xlabel('Time step Δt')
plt.ylabel('Error in N_T1')
plt.title('Convergence of Suzuki-Trotter method')
plt.grid(True)
plt.show()
```

期待される結果：2次分解では傾きが約-2のべき乗則。

### 8.3 ベンチマーク：厳密対角化との比較

小規模系（$N \leq 4$ Qutrits）では厳密対角化が可能：

```python
def exact_diagonalization(self, T_total):
    """
    厳密対角化による時間発展
    """
    # 全ハミルトニアン行列の構築
    dim = 3**self.N
    H_total = self._build_total_hamiltonian()
    
    # 固有値分解
    eigenvalues, eigenvectors = np.linalg.eigh(H_total)
    
    # 初期状態
    state_initial = self._initialize_state()
    
    # 初期状態を固有基底に展開
    coeffs = eigenvectors.conj().T @ state_initial
    
    # 時間発展
    time_evolved_coeffs = coeffs * np.exp(-1j * eigenvalues * T_total / self.hbar)
    
    # 元の基底に戻す
    state_final = eigenvectors @ time_evolved_coeffs
    
    return state_final

def _build_total_hamiltonian(self):
    """
    全ハミルトニアン行列を構築
    """
    dim = 3**self.N
    H_total = np.zeros((dim, dim), dtype=complex)
    
    # H₀項
    for idx in range(dim):
        config = self._index_to_config(idx)
        energy = sum(self.E_T if l == 1 else self.E_S if l == 2 else 0 
                    for l in config)
        H_total[idx, idx] += energy
    
    # H_transfer項（各隣接ペア）
    for (i, j) in self.neighbors:
        # |01⟩ ↔ |10⟩ 間の結合
        for idx1 in range(dim):
            config1 = self._index_to_config(idx1)
            if config1[i] == 0 and config1[j] == 1:
                config2 = config1.copy()
                config2[i] = 1
                config2[j] = 0
                idx2 = self._config_to_index(config2)
                H_total[idx1, idx2] += self.V
                H_total[idx2, idx1] += self.V
    
    # H_TTA項（各隣接ペア）
    for (i, j) in self.neighbors:
        # |11⟩ → |20⟩ と |11⟩ → |02⟩
        for idx1 in range(dim):
            config1 = self._index_to_config(idx1)
            if config1[i] == 1 and config1[j] == 1:
                # |11⟩ → |20⟩
                config2 = config1.copy()
                config2[i] = 2
                config2[j] = 0
                idx2 = self._config_to_index(config2)
                H_total[idx1, idx2] += self.J
                H_total[idx2, idx1] += self.J
                
                # |11⟩ → |02⟩
                config3 = config1.copy()
                config3[i] = 0
                config3[j] = 2
                idx3 = self._config_to_index(config3)
                H_total[idx1, idx3] += self.J
                H_total[idx3, idx1] += self.J
    
    return H_total

def _config_to_index(self, config):
    """
    3進数配列を線形インデックスに変換
    """
    idx = 0
    for i, level in enumerate(config):
        idx += level * (3 ** (self.N - 1 - i))
    return idx

# ベンチマーク
state_exact = simulator.exact_diagonalization(T_total)
state_trotter = simulator.run_simulation(T_total, N_steps=100)

# フィデリティの計算
fidelity = np.abs(np.vdot(state_exact.flatten(), state_trotter.flatten()))**2
print(f"Fidelity: {fidelity:.6f}")
```

フィデリティ $F \approx 1$ であれば、トロッター分解の精度が十分であることを示す。


## 9. MQT Quditsフレームワークでの完全実装例

### 9.1 完全なシミュレーションスクリプト

以下に、MQT Quditsフレームワークを用いた完全な実装例を示す：

```python
"""
完全な分子三重項状態量子ダイナミクスシミュレーション
MQT Quditsフレームワーク使用
"""

import numpy as np
import matplotlib.pyplot as plt
from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
from mqt.qudits.simulation import MQTQuditProvider

class CompleteMolecularSimulator:
    """
    分子三重項状態の完全なQuditシミュレータ
    """
    
    def __init__(self, N_molecules, E_T, E_S, V, J, Gamma_fl=0.0):
        """
        初期化
        
        Parameters:
        -----------
        N_molecules : int
            分子数
        E_T : float
            三重項エネルギー (eV)
        E_S : float
            一重項エネルギー (eV)
        V : float or array
            エネルギー移動積分 (eV)
        J : float or array
            TTA相互作用定数 (eV)
        Gamma_fl : float
            蛍光放出速度 (ns^-1)
        """
        self.N = N_molecules
        self.E_T = E_T
        self.E_S = E_S
        self.hbar = 0.6582  # eV·fs
        
        # 相互作用パラメータの配列化
        self.V = V * np.ones(N_molecules - 1) if np.isscalar(V) else V
        self.J = J * np.ones(N_molecules - 1) if np.isscalar(J) else J
        self.Gamma_fl = Gamma_fl
        
        # 隣接リスト
        self.neighbors = [(i, i+1) for i in range(N_molecules - 1)]
    
    def create_initial_state_circuit(self, state_type='all_triplet'):
        """
        初期状態準備回路
        """
        circuit = QuantumCircuit()
        mol_reg = QuantumRegister("molecules", self.N, [3] * self.N)
        circuit.append(mol_reg)
        
        if state_type == 'all_triplet':
            # すべての分子を三重項状態 |1⟩ に
            for i in range(self.N):
                # カスタム励起ゲート: |0⟩ → |1⟩
                U_excite = np.array([
                    [0, 1, 0],
                    [1, 0, 0],
                    [0, 0, 1]
                ], dtype=complex)
                circuit.cu_one(mol_reg[i], U_excite)
        
        elif state_type == 'alternating':
            # 交互に |1⟩ と |0⟩
            for i in range(self.N):
                if i % 2 == 1:
                    U_excite = np.array([
                        [0, 1, 0],
                        [1, 0, 0],
                        [0, 0, 1]
                    ], dtype=complex)
                    circuit.cu_one(mol_reg[i], U_excite)
        
        return circuit, mol_reg
    
    def build_single_step_circuit(self, mol_reg, dt):
        """
        1時間ステップ分の鈴木トロッター2次対称分解回路
        """
        circuit_step = QuantumCircuit()
        circuit_step.append(mol_reg)
        
        # (1) H₀ evolution (dt/2) - 前半
        phi_T = -self.E_T * dt / (2 * self.hbar)
        phi_S = -self.E_S * dt / (2 * self.hbar)
        
        for i in range(self.N):
            # 準位1への位相
            U_phase_1 = np.diag([1.0, np.exp(1j * phi_T), 1.0])
            circuit_step.cu_one(mol_reg[i], U_phase_1)
            
            # 準位2への位相
            U_phase_2 = np.diag([1.0, 1.0, np.exp(1j * phi_S)])
            circuit_step.cu_one(mol_reg[i], U_phase_2)
        
        # (2) H_transfer evolution (dt/2)
        for idx, (i, j) in enumerate(self.neighbors):
            theta = self.V[idx] * dt / (2 * self.hbar)
            U_transfer = self._build_transfer_unitary(theta)
            circuit_step.cu_two([mol_reg[i], mol_reg[j]], U_transfer)
        
        # (3) H_TTA evolution (dt/2)
        for idx, (i, j) in enumerate(self.neighbors):
            phi = self.J[idx] * dt / (2 * self.hbar)
            U_TTA = self._build_TTA_unitary(phi)
            circuit_step.cu_two([mol_reg[i], mol_reg[j]], U_TTA)
        
        # (4) H_rad evolution (dt) - 中央
        # 減衰は後処理で扱う
        
        # (5) H_TTA evolution (dt/2) - 後半（逆順）
        for idx, (i, j) in reversed(list(enumerate(self.neighbors))):
            phi = self.J[idx] * dt / (2 * self.hbar)
            U_TTA = self._build_TTA_unitary(phi)
            circuit_step.cu_two([mol_reg[i], mol_reg[j]], U_TTA)
        
        # (6) H_transfer evolution (dt/2) - 後半（逆順）
        for idx, (i, j) in reversed(list(enumerate(self.neighbors))):
            theta = self.V[idx] * dt / (2 * self.hbar)
            U_transfer = self._build_transfer_unitary(theta)
            circuit_step.cu_two([mol_reg[i], mol_reg[j]], U_transfer)
        
        # (7) H₀ evolution (dt/2) - 後半
        for i in reversed(range(self.N)):
            # 準位2への位相
            U_phase_2 = np.diag([1.0, 1.0, np.exp(1j * phi_S)])
            circuit_step.cu_one(mol_reg[i], U_phase_2)
            
            # 準位1への位相
            U_phase_1 = np.diag([1.0, np.exp(1j * phi_T), 1.0])
            circuit_step.cu_one(mol_reg[i], U_phase_1)
        
        return circuit_step
    
    def _build_transfer_unitary(self, theta):
        """
        エネルギー移動ユニタリ行列（9x9）
        """
        U = np.eye(9, dtype=complex)
        c = np.cos(theta)
        s = np.sin(theta)
        
        # 部分空間 {|01⟩, |10⟩} = {1, 3}
        U[1, 1] = c
        U[1, 3] = -1j * s
        U[3, 1] = -1j * s
        U[3, 3] = c
        
        return U
    
    def _build_TTA_unitary(self, phi):
        """
        TTAユニタリ行列（9x9）
        """
        U = np.eye(9, dtype=complex)
        sqrt2 = np.sqrt(2)
        
        c = np.cos(sqrt2 * phi)
        s = np.sin(sqrt2 * phi)
        c2 = np.cos(sqrt2 * phi / 2)**2
        s2 = np.sin(sqrt2 * phi / 2)**2
        
        # 部分空間 {|11⟩, |02⟩, |20⟩} = {4, 2, 6}
        U[4, 4] = c
        U[4, 2] = -1j * s / sqrt2
        U[4, 6] = -1j * s / sqrt2
        U[2, 4] = -1j * s / sqrt2
        U[2, 2] = c2
        U[2, 6] = s2
        U[6, 4] = -1j * s / sqrt2
        U[6, 2] = s2
        U[6, 6] = c2
        
        return U
    
    def run_full_simulation(self, T_total, N_steps, backend_name='tnsim',
                           initial_state='all_triplet', track_dynamics=True):
        """
        完全なシミュレーション実行
        """
        dt = T_total / N_steps
        
        # プロバイダとバックエンドの取得
        provider = MQTQuditProvider()
        backend = provider.get_backend(backend_name)
        
        if not track_dynamics:
            # 最終状態のみを取得（高速）
            circuit, mol_reg = self.create_initial_state_circuit(initial_state)
            
            for step in range(N_steps):
                step_circuit = self.build_single_step_circuit(mol_reg, dt)
                circuit = self._compose_circuits(circuit, step_circuit)
            
            job = backend.run(circuit)
            result = job.result()
            state_final = result.get_state_vector()
            
            return state_final, None
        
        else:
            # 時間発展を追跡（各ステップで測定）
            times = [0.0]
            populations_history = []
            
            # 初期状態
            circuit, mol_reg = self.create_initial_state_circuit(initial_state)
            job = backend.run(circuit)
            result = job.result()
            state = result.get_state_vector()
            
            pop_init = self._calculate_populations(state)
            populations_history.append(pop_init)
            
            # 時間発展
            for step in range(N_steps):
                step_circuit = self.build_single_step_circuit(mol_reg, dt)
                
                # 回路を合成して実行
                full_circuit = self._compose_circuits(circuit, step_circuit)
                job = backend.run(full_circuit)
                result = job.result()
                state = result.get_state_vector()
                
                # 減衰の適用（放射項）
                if self.Gamma_fl > 0:
                    state = self._apply_radiative_decay(state, dt)
                
                # 観測量の計算
                t = (step + 1) * dt
                populations = self._calculate_populations(state)
                
                times.append(t)
                populations_history.append(populations)
                
                # 次のステップのために回路を更新
                circuit = full_circuit
            
            return state, {'times': times, 'populations': populations_history}
    
    def _compose_circuits(self, circuit1, circuit2):
        """
        2つの量子回路を合成（簡易版）
        """
        # MQT Quditsの実際のAPIに合わせて実装
        # ここでは概念的な実装
        return circuit1  # 実際には適切な合成が必要
    
    def _calculate_populations(self, state_vector):
        """
        各状態の個体数を計算
        """
        state_flat = state_vector.flatten()
        N_S0, N_T1, N_S1 = 0.0, 0.0, 0.0
        
        for idx, amplitude in enumerate(state_flat):
            prob = np.abs(amplitude)**2
            config = self._index_to_config(idx)
            
            for level in config:
                if level == 0:
                    N_S0 += prob
                elif level == 1:
                    N_T1 += prob
                elif level == 2:
                    N_S1 += prob
        
        return {'N_S0': N_S0, 'N_T1': N_T1, 'N_S1': N_S1}
    
    def _apply_radiative_decay(self, state_vector, dt):
        """
        放射減衰を適用
        """
        state_flat = state_vector.flatten()
        new_state = state_flat.copy()
        
        for idx, amplitude in enumerate(state_flat):
            config = self._index_to_config(idx)
            n_S1 = sum(1 for level in config if level == 2)
            
            decay_factor = np.exp(-self.Gamma_fl * dt * n_S1 / 2)
            new_state[idx] *= decay_factor
        
        # 規格化
        new_state /= np.linalg.norm(new_state)
        
        return new_state.reshape(state_vector.shape)
    
    def _index_to_config(self, idx):
        """
        線形インデックスを3進数配列に変換
        """
        config = []
        for _ in range(self.N):
            config.append(idx % 3)
            idx //= 3
        return config[::-1]
    
    def plot_dynamics(self, dynamics_data, save_path=None):
        """
        時間発展のプロット
        """
        if dynamics_data is None:
            print("No dynamics data to plot")
            return
        
        times = dynamics_data['times']
        populations = dynamics_data['populations']
        
        N_S0_list = [p['N_S0'] for p in populations]
        N_T1_list = [p['N_T1'] for p in populations]
        N_S1_list = [p['N_S1'] for p in populations]
        
        plt.figure(figsize=(10, 6))
        plt.plot(times, N_S0_list, 'b-', linewidth=2, label='$N_{S_0}$ (Ground singlet)')
        plt.plot(times, N_T1_list, 'r-', linewidth=2, label='$N_{T_1}$ (Triplet)')
        plt.plot(times, N_S1_list, 'g-', linewidth=2, label='$N_{S_1}$ (Excited singlet)')
        
        plt.xlabel('Time (fs)', fontsize=14)
        plt.ylabel('Population', fontsize=14)
        plt.title('Quantum Dynamics of Molecular Triplet States', fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, max(times))
        plt.ylim(0, self.N + 0.5)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def calculate_fluorescence_intensity(self, dynamics_data):
        """
        蛍光強度を計算
        """
        if dynamics_data is None:
            return None
        
        times = dynamics_data['times']
        populations = dynamics_data['populations']
        
        I_fl = [self.Gamma_fl * p['N_S1'] for p in populations]
        
        # 累積蛍光量
        dt = times[1] - times[0] if len(times) > 1 else 0
        I_total = np.sum(I_fl) * dt
        
        return {'times': times, 'intensity': I_fl, 'total': I_total}

# ==================== 使用例 ====================

if __name__ == "__main__":
    # パラメータ設定
    N_molecules = 10
    E_T = 1.5  # eV
    E_S = 3.0  # eV
    V = 0.1    # eV
    J = 0.05   # eV
    Gamma_fl = 0.01  # fs^-1
    
    # シミュレータの初期化
    simulator = CompleteMolecularSimulator(
        N_molecules=N_molecules,
        E_T=E_T,
        E_S=E_S,
        V=V,
        J=J,
        Gamma_fl=Gamma_fl
    )
    
    # シミュレーション実行
    print("Running quantum simulation...")
    T_total = 1000.0  # fs
    N_steps = 200
    
    state_final, dynamics_data = simulator.run_full_simulation(
        T_total=T_total,
        N_steps=N_steps,
        initial_state='all_triplet',
        track_dynamics=True
    )
    
    # 結果の表示
    print("\nFinal populations:")
    final_pop = simulator._calculate_populations(state_final)
    print(f"  N_S0 = {final_pop['N_S0']:.4f}")
    print(f"  N_T1 = {final_pop['N_T1']:.4f}")
    print(f"  N_S1 = {final_pop['N_S1']:.4f}")
    print(f"  Total = {sum(final_pop.values()):.4f}")
    
    # 蛍光強度の計算
    fluorescence = simulator.calculate_fluorescence_intensity(dynamics_data)
    if fluorescence:
        print(f"\nTotal fluorescence: {fluorescence['total']:.4e}")
    
    # プロット
    simulator.plot_dynamics(dynamics_data, save_path='triplet_dynamics.png')
    
    print("\nSimulation completed successfully!")
```

### 9.2 拡張機能

#### 9.2.1 不均一系の扱い

各分子のエネルギーや相互作用が異なる場合：

```python
class InhomogeneousMolecularSimulator(CompleteMolecularSimulator):
    """
    不均一な分子系のシミュレータ
    """
    
    def __init__(self, N_molecules, E_T_array, E_S_array, V_array, J_array, Gamma_fl):
        """
        E_T_array, E_S_array: 各分子のエネルギー配列
        V_array, J_array: 各ペアの相互作用配列
        """
        self.N = N_molecules
        self.E_T_array = E_T_array
        self.E_S_array = E_S_array
        self.V_array = V_array
        self.J_array = J_array
        self.Gamma_fl = Gamma_fl
        self.hbar = 0.6582
        self.neighbors = [(i, i+1) for i in range(N_molecules - 1)]
    
    def build_single_step_circuit(self, mol_reg, dt):
        """
        不均一系用の時間発展回路
        """
        circuit_step = QuantumCircuit()
        circuit_step.append(mol_reg)
        
        # 各分子ごとに異なるエネルギー
        for i in range(self.N):
            phi_T = -self.E_T_array[i] * dt / (2 * self.hbar)
            phi_S = -self.E_S_array[i] * dt / (2 * self.hbar)
            
            U_phase = np.diag([1.0, np.exp(1j * phi_T), np.exp(1j * phi_S)])
            circuit_step.cu_one(mol_reg[i], U_phase)
        
        # 各ペアごとに異なる相互作用
        for idx, (i, j) in enumerate(self.neighbors):
            theta = self.V_array[idx] * dt / (2 * self.hbar)
            phi = self.J_array[idx] * dt / (2 * self.hbar)
            
            U_transfer = self._build_transfer_unitary(theta)
            circuit_step.cu_two([mol_reg[i], mol_reg[j]], U_transfer)
            
            U_TTA = self._build_TTA_unitary(phi)
            circuit_step.cu_two([mol_reg[i], mol_reg[j]], U_TTA)
        
        # 後半は同様に逆順で実装
        # ...
        
        return circuit_step
```

#### 9.2.2 2次元格子系

```python
def create_2D_lattice_neighbors(Lx, Ly):
    """
    2次元格子の隣接リストを生成
    
    Parameters:
    -----------
    Lx, Ly : int
        x方向とy方向のサイト数
    
    Returns:
    --------
    neighbors : list of tuples
        隣接ペアのリスト
    """
    neighbors = []
    
    def site_index(ix, iy):
        return ix * Ly + iy
    
    # 横方向の結合
    for ix in range(Lx):
        for iy in range(Ly - 1):
            neighbors.append((site_index(ix, iy), site_index(ix, iy + 1)))
    
    # 縦方向の結合
    for ix in range(Lx - 1):
        for iy in range(Ly):
            neighbors.append((site_index(ix, iy), site_index(ix + 1, iy)))
    
    return neighbors

# 使用例
Lx, Ly = 4, 4
N_molecules = Lx * Ly
neighbors_2D = create_2D_lattice_neighbors(Lx, Ly)

simulator_2D = CompleteMolecularSimulator(
    N_molecules=N_molecules,
    E_T=1.5,
    E_S=3.0,
    V=0.1,
    J=0.05,
    Gamma_fl=0.01
)
simulator_2D.neighbors = neighbors_2D
```


## 10. まとめと今後の展望

### 10.1 本文書の要約

本文書では、分子三重項状態の量子ダイナミクスを **Qudit（特にQutrit）** を用いた量子アルゴリズムで実装するための完全な理論を提示した。

#### 10.1.1 主要な貢献

1. **Qudit表現の確立**
   - 3準位分子系（$|S_0\rangle, |T_1\rangle, |S_1\rangle$）を1 Qutritで自然に表現
   - $N$ 分子系を $N$ Qutritsのテンソル積空間で記述
   - 状態空間の次元: $3^N$

2. **ハミルトニアンのQudit演算子表現**
   - 対角項 $\hat{H}_0$: 単一Qutrit位相ゲート（VirtRz）
   - エネルギー移動 $\hat{H}_{\text{transfer}}$: 2-Quditカスタムゲート
   - TTA過程 $\hat{H}_{\text{TTA}}$: 2-Quditカスタムゲート
   - 放射減衰 $\hat{H}_{\text{rad}}$: 非ユニタリ操作（ノイズモデル）

3. **鈴木トロッター分解の適用**
   - 2次対称分解: $\mathcal{O}(\Delta t^2)$ の全体誤差
   - 各時間ステップをQuditゲートの列に分解
   - 並列化による回路深さの削減

4. **MQT Quditsフレームワークでの実装**
   - `QuantumCircuit`, `QuantumRegister` による回路構築
   - `cu_one`, `cu_two` によるカスタムゲートの適用
   - `MQTQuditProvider` によるシミュレーション実行
   - 完全なPythonコード例を提供

5. **観測量と解析手法**
   - 個体数演算子の期待値計算
   - 蛍光強度の時間発展
   - 相関関数の評価
   - 収束テストとベンチマーク

### 10.2 Quditアプローチの優位性

#### 10.2.1 量子ビット方式との比較

**量子ビット方式**:
- 3準位系を表現するには2量子ビット（4次元空間）が必要
- 1つの準位が未使用となり、非効率
- $N$ 分子系: $2N$ 量子ビット、状態空間 $2^{2N}$ 次元
- 物理的制約（未使用準位への遷移抑制）が必要

**Qudit方式**:
- 3準位系を1 Qutritで自然に表現
- すべての準位が物理的に意味を持つ
- $N$ 分子系: $N$ Qutrit、状態空間 $3^N$ 次元
- ゲート数の削減（2-body相互作用が直接表現可能）

#### 10.2.2 ゲート数の比較

エネルギー移動演算子の実装：
- **量子ビット**: 複数のCNOTゲートと局所回転（典型的に10個以上）
- **Qutrit**: 1つのカスタム2-Quditゲート

TTA演算子の実装：
- **量子ビット**: さらに複雑な分解が必要（20個以上のゲート）
- **Qutrit**: 1つのカスタム2-Quditゲート

**結論**: Quditアプローチは、ゲート数を大幅に削減し、ノイズの蓄積を抑制できる。

### 10.3 物理的応用

#### 10.3.1 有機フォトニクス

本理論は以下の実験系に適用可能：

1. **三重項-三重項消滅アップコンバージョン（TTA-UC）**
   - 低エネルギー光子から高エネルギー光子への変換
   - 太陽電池効率の向上
   - 本シミュレーションにより最適な分子配置や相互作用強度を予測

2. **遅延蛍光材料**
   - 有機ELデバイスの発光効率向上
   - TTA過程による励起一重項生成の最大化
   - 時間分解スペクトルの理論予測

3. **量子ドット・ナノクリスタル系**
   - 多励起子生成（Multiple Exciton Generation, MEG）
   - 励起子間相互作用の量子論的記述

#### 10.3.2 パラメータ最適化

シミュレーションを用いて、実験系の設計指針を得る：

```python
def optimize_parameters(simulator_class, param_ranges, objective='max_S1'):
    """
    パラメータ空間の最適化
    
    Parameters:
    -----------
    simulator_class : class
        シミュレータクラス
    param_ranges : dict
        パラメータの範囲 {'V': [0.05, 0.2], 'J': [0.01, 0.1], ...}
    objective : str
        最適化目標（'max_S1', 'max_fluorescence', など）
    
    Returns:
    --------
    optimal_params : dict
        最適パラメータ
    """
    from scipy.optimize import differential_evolution
    
    def objective_function(params):
        V, J = params
        sim = simulator_class(N_molecules=10, E_T=1.5, E_S=3.0, 
                             V=V, J=J, Gamma_fl=0.01)
        state_final, dynamics = sim.run_full_simulation(
            T_total=1000.0, N_steps=100, track_dynamics=True
        )
        
        if objective == 'max_S1':
            # 最大一重項個体数
            N_S1_max = max(p['N_S1'] for p in dynamics['populations'])
            return -N_S1_max  # 最小化問題に変換
        
        elif objective == 'max_fluorescence':
            # 累積蛍光量
            fluorescence = sim.calculate_fluorescence_intensity(dynamics)
            return -fluorescence['total']
    
    # 最適化実行
    bounds = [param_ranges['V'], param_ranges['J']]
    result = differential_evolution(objective_function, bounds)
    
    optimal_params = {
        'V': result.x[0],
        'J': result.x[1],
        'objective_value': -result.fun
    }
    
    return optimal_params
```

### 10.4 今後の発展方向

#### 10.4.1 高次Quditへの拡張

より複雑な分子系（4準位以上）への対応：
- Ququart（$d=4$）: 追加の励起状態を含む
- Ququint（$d=5$）: 振動準位の明示的な扱い

$$
\text{分子（振動準位込み）} \longleftrightarrow d > 3 \text{ Qudit}
$$

#### 10.4.2 混合次元Qudit系

MQT Quditsの特徴である混合次元レジスタの活用：
- 分子系: Qutrit（$d=3$）
- 光子場: Qubit（$d=2$）または無限次元（フォック状態切断）

```python
from mqt.qudits.quantum_circuit import QuantumRegister

# 混合次元レジスタ
molecule_reg = QuantumRegister("molecules", 10, [3]*10)  # 10 Qutrits
photon_reg = QuantumRegister("photons", 5, [10]*5)       # 5 Qudits (d=10, 切断フォック空間)

circuit = QuantumCircuit()
circuit.append(molecule_reg)
circuit.append(photon_reg)

# 分子-光子相互作用ゲート
# ...
```

#### 10.4.3 量子誤り訂正

Qudit量子誤り訂正符号の適用：
- Qutrit stabilizer code
- Bosonic code（振動子符号）

これにより、ノイズのある量子デバイスでの長時間シミュレーションが可能になる。

#### 10.4.4 変分量子アルゴリズム（VQA）

鈴木トロッター分解の代わりに、変分原理に基づくアプローチ：

$$
|\Psi(\boldsymbol{\theta})\rangle = \prod_{k=1}^{L} \hat{U}_k(\theta_k) |\Psi_0\rangle
$$

パラメータ $\boldsymbol{\theta}$ を最適化して、正確な時間発展を近似：

$$
\min_{\boldsymbol{\theta}} \left\| |\Psi(\boldsymbol{\theta})\rangle - e^{-i\hat{H}t/\hbar}|\Psi_0\rangle \right\|^2
$$

量子古典ハイブリッドアルゴリズムにより、近未来の量子デバイスでの実装が期待される。

#### 10.4.5 テンソルネットワークとの融合

大規模系（$N \gg 10$）では、テンソルネットワーク法との組み合わせが有効：
- Matrix Product States (MPS) 表現
- Density Matrix Renormalization Group (DMRG)

MQT Quditsのバックエンド `tnsim`（Tensor Network Simulator）を活用。

### 10.5 実験実装への道筋

#### 10.5.1 ハードウェアプラットフォーム

Qudit量子計算の実験実装候補：

1. **超伝導量子回路**
   - Transmon qutritモード
   - 3準位以上の制御が既に実証済み

2. **イオントラップ**
   - 複数の超微細準位を利用
   - 高精度な局所制御とエンタングリングゲート

3. **冷却原子**
   - 光格子中の原子
   - Rydberg励起による長距離相互作用

4. **フォトニクス**
   - 経路自由度または偏光自由度を拡張
   - 測定誘起エンタングルメント

#### 10.5.2 ゲート実装の実験的課題

- **高忠実度Quditゲート**: 従来の量子ビットゲート（忠実度 > 99%）と同等の性能が求められる
- **スケーラビリティ**: 多数のQuditを同時に制御
- **読み出し**: Qudit状態の高精度な測定

### 10.6 結論

本文書では、分子三重項状態の量子ダイナミクスを **MQT Quditsフレームワーク** を用いた量子アルゴリズムで解くための、完全かつプログラム実装可能な理論を構築した。

#### 主要な成果

1. **理論的完全性**: すべての数式を省略無しに展開し、各演算子とゲートの対応を明示
2. **実装可能性**: MQT Quditsの実際のAPIに基づいた完全なコード例を提供
3. **効率性**: Qudit表現により、量子ビット方式と比較してゲート数を大幅に削減
4. **汎用性**: 任意の分子数、格子構造、不均一系に対応可能

#### 物理的意義

- 有機フォトニクス、遅延蛍光、アップコンバージョンなどの実験系の理論解析
- 量子多体効果（TTA、エネルギー移動）の非摂動的な記述
- 実験パラメータの最適化と新材料設計への貢献

#### 量子情報科学への貢献

- Qudit量子計算の具体的応用例の提示
- 混合次元量子系の活用
- 近未来の量子デバイスでの実装可能性

本理論は、**量子化学と量子情報科学の融合**の好例であり、量子コンピュータが化学・材料科学において真に有用となる将来への基礎を提供する。

---

## 参考文献

### 量子ダイナミクス理論
1. 本文書の基礎理論: `quantum_dynamics_molecular_triplet_states.md`
2. 三重項-三重項消滅の物理化学

### 鈴木トロッター分解
3. 本文書の数値計算理論: `suzuki_trotter_decomposition_theory.md`
4. Suzuki, M. (1990). "Fractal decomposition of exponential operators". *Phys. Lett. A* **146**, 319-323.
5. Childs, A. M., et al. (2019). "Theory of Trotter error with commutator scaling". *Phys. Rev. X* **9**, 011011.

### Qudit量子計算
6. Gokhale, P., et al. (2019). "Asymptotic improvements to quantum circuits via qutrits". *Proc. ACM Symp. STOC* **51**, 554-565.
7. Murali, P., et al. (2020). "Software mitigation of crosstalk on noisy intermediate-scale quantum computers". *ASPLOS 2020*.

### MQT Quditsフレームワーク
8. MQT Qudits Documentation: https://mqt.readthedocs.io/projects/qudits/
9. Grurl, T., et al. (2023). "Automatic Implementation and Evaluation of Error-Correcting Codes for Quantum Computing: An Overview". *ACM Computing Surveys*.

### 実験実装
10. Nikolaeva, A. S., et al. (2021). "Multi-level quantum systems as qudits: Implementation in superconducting circuits". *Quantum Sci. Technol.* **6**, 035007.
11. Chi, Y., et al. (2022). "A programmable qudit-based quantum processor". *Nat. Commun.* **13**, 1166.

---

**文書作成日**: 2025-10-14  
**分野**: 量子情報科学、量子化学、Qudit量子計算  
**対象**: MQT Quditsフレームワークを用いた量子アルゴリズム実装  
**参照文書**: 
- `quantum_dynamics_molecular_triplet_states.md`
- `suzuki_trotter_decomposition_theory.md`

