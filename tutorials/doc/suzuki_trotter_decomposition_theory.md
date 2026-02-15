# 鈴木トロッター分解による分子三重項状態の量子ダイナミクス数値計算理論

## 1. はじめに

本文書では、`quantum_dynamics_molecular_triplet_states.md` で記述された分子系の量子ダイナミクスを、鈴木トロッター分解（Suzuki-Trotter decomposition）を用いて数値的に解くための詳細な理論を、プログラム実装可能なレベルで記述する。

### 1.1 目的

時間依存シュレーディンガー方程式：

$$
i\hbar \frac{\partial |\Psi(t)\rangle}{\partial t} = \hat{H}_{\text{total}} |\Psi(t)\rangle
$$

の形式解は時間発展演算子（time evolution operator）を用いて：

$$
|\Psi(t)\rangle = \hat{U}(t) |\Psi(0)\rangle
$$

と表される。ここで：

$$
\hat{U}(t) = e^{-i\hat{H}_{\text{total}}t/\hbar}
$$

しかし、ハミルトニアンが複数の非可換な項の和である場合：

$$
\hat{H}_{\text{total}} = \hat{H}_1 + \hat{H}_2 + \cdots + \hat{H}_M
$$

一般に以下は成立しない：

$$
e^{-i(\hat{H}_1 + \hat{H}_2 + \cdots + \hat{H}_M)t/\hbar} \neq e^{-i\hat{H}_1 t/\hbar} e^{-i\hat{H}_2 t/\hbar} \cdots e^{-i\hat{H}_M t/\hbar}
$$

鈴木トロッター分解は、この問題を系統的に解決し、時間発展演算子を近似的に分解する手法である。

## 2. Baker-Campbell-Hausdorffの公式と非可換性の影響

### 2.1 BCH公式

2つの演算子 $\hat{A}$, $\hat{B}$ に対して、Baker-Campbell-Hausdorff（BCH）公式は：

$$
e^{\hat{A}} e^{\hat{B}} = e^{\hat{A} + \hat{B} + \frac{1}{2}[\hat{A}, \hat{B}] + \frac{1}{12}[\hat{A}, [\hat{A}, \hat{B}]] - \frac{1}{12}[\hat{B}, [\hat{A}, \hat{B}]] + \cdots}
$$

ここで、$[\hat{A}, \hat{B}] = \hat{A}\hat{B} - \hat{B}\hat{A}$ は交換子である。

### 2.2 素朴な積分解の誤差

時間発展を単純に分解すると：

$$
e^{-i(\hat{H}_1 + \hat{H}_2)t/\hbar} \approx e^{-i\hat{H}_1 t/\hbar} e^{-i\hat{H}_2 t/\hbar}
$$

BCH公式より、誤差項は：

$$
\Delta = -\frac{t^2}{2\hbar^2}[\hat{H}_1, \hat{H}_2] + \mathcal{O}(t^3)
$$

この誤差は $\mathcal{O}(t^2)$ であり、長時間発展では蓄積して大きくなる。

## 3. 鈴木トロッター分解の基礎理論

### 3.1 時間刻み分割

時間区間 $[0, T]$ を $N$ 個の小区間に分割する：

$$
\Delta t = \frac{T}{N}
$$

各時刻を：

$$
t_n = n \Delta t, \quad n = 0, 1, 2, \ldots, N
$$

と定義する。

### 3.2 時間発展演算子の分解

全時間発展演算子は：

$$
\hat{U}(T) = \hat{U}(\Delta t)^N + \mathcal{O}(\Delta t^{p+1})
$$

ここで、$p$ は分解の次数（order）である。

## 4. 1次の鈴木トロッター分解（Lie-Trotter分解）

### 4.1 定義

1次の鈴木トロッター分解（Lie-Trotter splitting）は：

$$
e^{-i(\hat{H}_1 + \hat{H}_2)t/\hbar} = \lim_{N \to \infty} \left( e^{-i\hat{H}_1 \Delta t/\hbar} e^{-i\hat{H}_2 \Delta t/\hbar} \right)^N
$$

ここで、$\Delta t = t/N$ である。

### 4.2 誤差評価

1回の時間ステップ $\Delta t$ における誤差は：

$$
\left\| e^{-i(\hat{H}_1 + \hat{H}_2)\Delta t/\hbar} - e^{-i\hat{H}_1 \Delta t/\hbar} e^{-i\hat{H}_2 \Delta t/\hbar} \right\| = \mathcal{O}(\Delta t^2)
$$

$N$ ステップの累積誤差は：

$$
\varepsilon_{\text{total}} = N \cdot \mathcal{O}(\Delta t^2) = \frac{T}{\Delta t} \cdot \mathcal{O}(\Delta t^2) = \mathcal{O}(\Delta t)
$$

したがって、1次分解の全体誤差は $\mathcal{O}(\Delta t)$ である。

### 4.3 具体的なアルゴリズム（1次）

初期状態 $|\Psi(0)\rangle$ から時刻 $T$ の状態 $|\Psi(T)\rangle$ を求める：

**入力**:

- 初期状態: $|\Psi_0\rangle = |\Psi(0)\rangle$
- ハミルトニアン: $\hat{H}_1, \hat{H}_2, \ldots, \hat{H}_M$
- 全時間: $T$
- 時間ステップ数: $N$

**計算**:

1. 時間刻み幅を計算: $\Delta t = T / N$
2. $|\Psi\rangle \leftarrow |\Psi_0\rangle$ と初期化
3. $n = 1$ から $N$ まで繰り返し:
   - $m = 1$ から $M$ まで繰り返し:
     - $|\Psi\rangle \leftarrow e^{-i\hat{H}_m \Delta t/\hbar} |\Psi\rangle$ を計算
4. 出力: $|\Psi(T)\rangle = |\Psi\rangle$

### 4.4 実装の要点

各 $e^{-i\hat{H}_m \Delta t/\hbar}$ の計算方法：

**a) ハミルトニアンが対角化可能な場合**

$\hat{H}_m$ が固有値 $E_k^{(m)}$ と固有ベクトル $|k^{(m)}\rangle$ を持つ場合：

$$
\hat{H}_m = \sum_k E_k^{(m)} |k^{(m)}\rangle\langle k^{(m)}|
$$

時間発展演算子は：

$$
e^{-i\hat{H}_m \Delta t/\hbar} = \sum_k e^{-iE_k^{(m)} \Delta t/\hbar} |k^{(m)}\rangle\langle k^{(m)}|
$$

状態への作用は：

$$
|\Psi'\rangle = e^{-i\hat{H}_m \Delta t/\hbar} |\Psi\rangle = \sum_k e^{-iE_k^{(m)} \Delta t/\hbar} \langle k^{(m)}|\Psi\rangle |k^{(m)}\rangle
$$

**実装ステップ**:

1. $|\Psi\rangle$ を固有基底 $\{|k^{(m)}\rangle\}$ に展開: $c_k = \langle k^{(m)}|\Psi\rangle$
2. 各係数に位相を掛ける: $c_k' = e^{-iE_k^{(m)} \Delta t/\hbar} c_k$
3. 固有基底から元の基底に戻す: $|\Psi'\rangle = \sum_k c_k' |k^{(m)}\rangle$

**b) ハミルトニアンが行列表現される場合**

ハミルトニアン行列 $H_m$ に対して：

$$
U_m(\Delta t) = \exp\left(-i\frac{\Delta t}{\hbar} H_m\right)
$$

行列指数関数の計算方法：

1. **固有値分解法**: $H_m = V D V^\dagger$ と分解し、$U_m = V e^{-iD\Delta t/\hbar} V^\dagger$
2. **Padé近似法**: 有理関数近似を用いる
3. **テイラー展開法**: $e^X = \sum_{n=0}^\infty \frac{X^n}{n!}$ を有限項で打ち切る

## 5. 2次の鈴木トロッター分解（対称分解）

### 5.1 定義

2次の対称鈴木トロッター分解（Strange splitting）は：

$$
e^{-i(\hat{H}_1 + \hat{H}_2)t/\hbar} = \lim_{N \to \infty} \left( e^{-i\hat{H}_1 \Delta t/(2\hbar)} e^{-i\hat{H}_2 \Delta t/\hbar} e^{-i\hat{H}_1 \Delta t/(2\hbar)} \right)^N + \mathcal{O}(\Delta t^2)
$$

### 5.2 対称性の利点

対称形式により、奇数次の誤差項がキャンセルされる：

$$
e^{-i(\hat{H}_1 + \hat{H}_2)\Delta t/\hbar} = e^{-i\hat{H}_1 \Delta t/(2\hbar)} e^{-i\hat{H}_2 \Delta t/\hbar} e^{-i\hat{H}_1 \Delta t/(2\hbar)} + \mathcal{O}(\Delta t^3)
$$

### 5.3 誤差評価

1ステップの誤差: $\mathcal{O}(\Delta t^3)$

$N$ ステップの累積誤差: $\mathcal{O}(\Delta t^2)$

### 5.4 具体的なアルゴリズム（2次、2項の場合）

**入力**:

- 初期状態: $|\Psi_0\rangle$
- ハミルトニアン: $\hat{H}_1, \hat{H}_2$
- 全時間: $T$
- 時間ステップ数: $N$

**計算**:

1. $\Delta t = T / N$
2. $|\Psi\rangle \leftarrow |\Psi_0\rangle$
3. $n = 1$ から $N$ まで繰り返し:
   - $|\Psi\rangle \leftarrow e^{-i\hat{H}_1 \Delta t/(2\hbar)} |\Psi\rangle$
   - $|\Psi\rangle \leftarrow e^{-i\hat{H}_2 \Delta t/\hbar} |\Psi\rangle$
   - $|\Psi\rangle \leftarrow e^{-i\hat{H}_1 \Delta t/(2\hbar)} |\Psi\rangle$
4. 出力: $|\Psi(T)\rangle = |\Psi\rangle$

### 5.5 一般の多項ハミルトニアンへの拡張

$M$ 個のハミルトニアン項 $\hat{H}_1, \hat{H}_2, \ldots, \hat{H}_M$ に対する2次対称分解：

$$
\hat{S}_2(\Delta t) = e^{-i\hat{H}_1 \Delta t/(2\hbar)} e^{-i\hat{H}_2 \Delta t/(2\hbar)} \cdots e^{-i\hat{H}_M \Delta t/\hbar} \cdots e^{-i\hat{H}_2 \Delta t/(2\hbar)} e^{-i\hat{H}_1 \Delta t/(2\hbar)}
$$

具体的な順序：

$$
\hat{S}_2(\Delta t) = \prod_{m=1}^{M-1} e^{-i\hat{H}_m \Delta t/(2\hbar)} \cdot e^{-i\hat{H}_M \Delta t/\hbar} \cdot \prod_{m=M-1}^{1} e^{-i\hat{H}_m \Delta t/(2\hbar)}
$$

**アルゴリズム**:

1. $\Delta t = T / N$
2. $|\Psi\rangle \leftarrow |\Psi_0\rangle$
3. $n = 1$ から $N$ まで繰り返し:
   - 前半: $m = 1$ から $M-1$ まで $|\Psi\rangle \leftarrow e^{-i\hat{H}_m \Delta t/(2\hbar)} |\Psi\rangle$
   - 中央: $|\Psi\rangle \leftarrow e^{-i\hat{H}_M \Delta t/\hbar} |\Psi\rangle$
   - 後半: $m = M-1$ から $1$ まで $|\Psi\rangle \leftarrow e^{-i\hat{H}_m \Delta t/(2\hbar)} |\Psi\rangle$
4. 出力: $|\Psi(T)\rangle$

## 6. 4次の鈴木トロッター分解

### 6.1 4次分解の構成

4次の鈴木分解は、2次分解を再帰的に組み合わせて構成される。鈴木の公式により：

$$
\hat{S}_4(\Delta t) = \hat{S}_2(p\Delta t) \hat{S}_2(p\Delta t) \hat{S}_2((1-4p)\Delta t) \hat{S}_2(p\Delta t) \hat{S}_2(p\Delta t)
$$

ここで：

$$
p = \frac{1}{4 - 4^{1/3}} \approx 0.414907
$$

### 6.2 具体的な係数

数値的には：

$$
p_1 = p_2 = p_4 = p_5 = \frac{1}{4 - 4^{1/3}}
$$

$$
p_3 = 1 - 4p = \frac{-4^{1/3}}{4 - 4^{1/3}}
$$

注意: $p_3 < 0$ であるため、時間を「逆行」させる演算が含まれる。

### 6.3 アルゴリズム（4次）

**入力**: $|\Psi_0\rangle$, $\hat{H}_1, \hat{H}_2$, $T$, $N$

**計算**:

1. $\Delta t = T / N$
2. $p = 1/(4 - 4^{1/3})$
3. 時間刻み: $\tau_1 = \tau_2 = \tau_4 = \tau_5 = p\Delta t$, $\tau_3 = (1-4p)\Delta t$
4. $|\Psi\rangle \leftarrow |\Psi_0\rangle$
5. $n = 1$ から $N$ まで繰り返し:
   - $j = 1$ から $5$ まで繰り返し:
     - 2次対称分解 $\hat{S}_2(\tau_j)$ を適用:
       - $|\Psi\rangle \leftarrow e^{-i\hat{H}_1 \tau_j/(2\hbar)} |\Psi\rangle$
       - $|\Psi\rangle \leftarrow e^{-i\hat{H}_2 \tau_j/\hbar} |\Psi\rangle$
       - $|\Psi\rangle \leftarrow e^{-i\hat{H}_1 \tau_j/(2\hbar)} |\Psi\rangle$
6. 出力: $|\Psi(T)\rangle$

### 6.4 誤差評価

4次分解の全体誤差は $\mathcal{O}(\Delta t^4)$ である。

## 7. 分子三重項系への適用

### 7.1 ハミルトニアンの分解

`quantum_dynamics_molecular_triplet_states.md` の全ハミルトニアン：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}} + \hat{H}_{\text{rad}}
$$

ここで：

$$
\hat{H}_0 = \sum_{i=1}^{N} \left( E_T \hat{n}_{T_1}^{(i)} + E_S \hat{n}_{S_1}^{(i)} \right)
$$

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( \hat{T}_{ij}^{\text{transfer}} + \hat{T}_{ij}^{\text{transfer}\dagger} \right)
$$

$$
\hat{H}_{\text{TTA}} = \sum_{\langle i,j \rangle} J_{ij} \left( \hat{T}_{ij}^{\text{TTA}} + \hat{T}_{ij}^{\text{TTA}\dagger} \right)
$$

$$
\hat{H}_{\text{rad}} = \sum_i \sum_{\mathbf{k},\lambda} g_{\mathbf{k}\lambda} \left( \hat{A}_i \hat{a}_{\mathbf{k}\lambda}^\dagger + \hat{A}_i^\dagger \hat{a}_{\mathbf{k}\lambda} \right)
$$

### 7.2 分解戦略

以下のようにハミルトニアンを4つの部分に分解する：

- $\hat{H}_1 = \hat{H}_0$ （対角項、簡単に計算可能）
- $\hat{H}_2 = \hat{H}_{\text{transfer}}$ （エネルギー移動）
- $\hat{H}_3 = \hat{H}_{\text{TTA}}$ （TTA過程）
- $\hat{H}_4 = \hat{H}_{\text{rad}}$ （放射過程）

### 7.3 各項の時間発展演算子

#### 7.3.1 $\hat{H}_0$ の時間発展

$\hat{H}_0$ は対角的なので：

$$
e^{-i\hat{H}_0 t/\hbar} = \exp\left(-i\frac{t}{\hbar} \sum_{i=1}^{N} \left( E_T \hat{n}_{T_1}^{(i)} + E_S \hat{n}_{S_1}^{(i)} \right)\right)
$$

$$
= \prod_{i=1}^{N} e^{-i E_T \hat{n}_{T_1}^{(i)} t/\hbar} e^{-i E_S \hat{n}_{S_1}^{(i)} t/\hbar}
$$

各分子 $i$ の基底 $\{|S_0\rangle_i, |T_1\rangle_i, |S_1\rangle_i\}$ において：

$$
e^{-i\hat{H}_0 t/\hbar} |S_0\rangle_i = |S_0\rangle_i
$$

$$
e^{-i\hat{H}_0 t/\hbar} |T_1\rangle_i = e^{-iE_T t/\hbar} |T_1\rangle_i
$$

$$
e^{-i\hat{H}_0 t/\hbar} |S_1\rangle_i = e^{-iE_S t/\hbar} |S_1\rangle_i
$$

**実装**:

状態 $|\Psi\rangle$ を基底 $|s_1, s_2, \ldots, s_N\rangle$ で展開（$s_i \in \{S_0, T_1, S_1\}$）：

$$
|\Psi\rangle = \sum_{s_1, \ldots, s_N} c_{s_1 \ldots s_N} |s_1, \ldots, s_N\rangle
$$

時間発展後：

$$
|\Psi'\rangle = \sum_{s_1, \ldots, s_N} c_{s_1 \ldots s_N} \exp\left(-i\frac{t}{\hbar} \sum_{i=1}^N E(s_i)\right) |s_1, \ldots, s_N\rangle
$$

ここで：

$$
E(s_i) = \begin{cases}
0 & \text{if } s_i = S_0 \\
E_T & \text{if } s_i = T_1 \\
E_S & \text{if } s_i = S_1
\end{cases}
$$

**疑似コード**:

```
function apply_H0_evolution(state, dt):
    for each basis_state in state.basis:
        energy = 0
        for i = 1 to N:
            if basis_state[i] == T1:
                energy += E_T
            elif basis_state[i] == S1:
                energy += E_S
        phase = exp(-1j * energy * dt / hbar)
        state.coefficients[basis_state] *= phase
    return state
```

#### 7.3.2 $\hat{H}_{\text{transfer}}$ の時間発展

エネルギー移動項は：

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( |S_0\rangle_i\langle T_1| \otimes |T_1\rangle_j\langle S_0| + |T_1\rangle_i\langle S_0| \otimes |S_0\rangle_j\langle T_1| \right)
$$

各隣接ペア $(i,j)$ に対する局所ハミルトニアンは：

$$
\hat{h}_{ij}^{\text{transfer}} = V_{ij} \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & V_{ij} & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & V_{ij} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

ここで、基底の順序は $\{|S_0,S_0\rangle, |S_0,T_1\rangle, |S_0,S_1\rangle, |T_1,S_0\rangle, |T_1,T_1\rangle, |T_1,S_1\rangle, |S_1,S_0\rangle, |S_1,T_1\rangle, |S_1,S_1\rangle\}$ である。

この2×2ブロック対角形式は固有値解析が容易である：

関連する部分空間 $\{|S_0,T_1\rangle, |T_1,S_0\rangle\}$ において：

$$
\hat{h}_{ij}^{\text{transfer}} = V_{ij} \begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix} = V_{ij} \sigma_x
$$

固有値: $\pm V_{ij}$

固有ベクトル: $|\pm\rangle = \frac{1}{\sqrt{2}}(|S_0,T_1\rangle \pm |T_1,S_0\rangle)$

時間発展演算子:

$$
e^{-i\hat{h}_{ij}^{\text{transfer}} t/\hbar} = \cos(V_{ij}t/\hbar) \hat{I} - i\sin(V_{ij}t/\hbar) \sigma_x
$$

基底表現では：

$$
e^{-i\hat{h}_{ij}^{\text{transfer}} t/\hbar} = \begin{pmatrix}
\cos\theta & -i\sin\theta \\
-i\sin\theta & \cos\theta
\end{pmatrix}
$$

ここで、$\theta = V_{ij}t/\hbar$ である。

**実装ステップ**:

1. 隣接ペア $(i,j)$ を列挙
2. 各ペアに対して、関連する4状態 $\{|S_0,T_1\rangle, |T_1,S_0\rangle\}$ の係数を抽出
3. 2×2行列の時間発展を適用
4. 係数を更新

**疑似コード**:

```
function apply_transfer_evolution(state, dt):
    for each pair (i, j) in adjacent_pairs:
        V = V_ij
        theta = V * dt / hbar

        # 抽出
        c_01 = state.coefficient(i=S0, j=T1)
        c_10 = state.coefficient(i=T1, j=S0)

        # 時間発展
        c_01_new = cos(theta) * c_01 - 1j * sin(theta) * c_10
        c_10_new = -1j * sin(theta) * c_01 + cos(theta) * c_10

        # 更新
        state.set_coefficient(i=S0, j=T1, c_01_new)
        state.set_coefficient(i=T1, j=S0, c_10_new)

    return state
```

#### 7.3.3 $\hat{H}_{\text{TTA}}$ の時間発展

TTA項は：

$$
\hat{H}_{\text{TTA}} = \sum_{\langle i,j \rangle} J_{ij} \left( |S_1\rangle_i\langle T_1| \otimes |S_0\rangle_j\langle T_1| + |S_0\rangle_i\langle T_1| \otimes |S_1\rangle_j\langle T_1| + \text{h.c.} \right)
$$

関連する部分空間は $\{|T_1,T_1\rangle, |S_1,S_0\rangle, |S_0,S_1\rangle\}$ である。

9×9行列の該当ブロック：

$$
\hat{h}_{ij}^{\text{TTA}} = J_{ij} \begin{pmatrix}
0 & 0 & 0 \\
0 & 0 & J_{ij} \\
0 & J_{ij} & 0
\end{pmatrix}
$$

ただし、基底順序は $\{|T_1,T_1\rangle, |S_1,S_0\rangle, |S_0,S_1\rangle\}$ である。

この行列の固有値問題を解く：

$$
\begin{pmatrix}
0 & 0 & 0 \\
0 & 0 & 1 \\
0 & 1 & 0
\end{pmatrix}
$$

固有値: $0, \pm 1$

固有ベクトル:

- $\lambda = 0$: $|0\rangle = (1, 0, 0)^T$ すなわち $|T_1,T_1\rangle$
- $\lambda = +1$: $|+\rangle = (0, 1, 1)^T/\sqrt{2}$ すなわち $(|S_1,S_0\rangle + |S_0,S_1\rangle)/\sqrt{2}$
- $\lambda = -1$: $|-\rangle = (0, 1, -1)^T/\sqrt{2}$ すなわち $(|S_1,S_0\rangle - |S_0,S_1\rangle)/\sqrt{2}$

時間発展演算子:

$$
e^{-i\hat{h}_{ij}^{\text{TTA}} t/\hbar} = |0\rangle\langle 0| + e^{-iJ_{ij}t/\hbar} |+\rangle\langle +| + e^{+iJ_{ij}t/\hbar} |-\rangle\langle -|
$$

行列表現に戻すと：

$$
e^{-i\hat{h}_{ij}^{\text{TTA}} t/\hbar} = \begin{pmatrix}
1 & 0 & 0 \\
0 & \cos\phi & -i\sin\phi \\
0 & -i\sin\phi & \cos\phi
\end{pmatrix}
$$

ここで、$\phi = J_{ij}t/\hbar$ である。

**疑似コード**:

```
function apply_TTA_evolution(state, dt):
    for each pair (i, j) in adjacent_pairs:
        J = J_ij
        phi = J * dt / hbar

        # 抽出
        c_TT = state.coefficient(i=T1, j=T1)
        c_S0 = state.coefficient(i=S1, j=S0)
        c_0S = state.coefficient(i=S0, j=S1)

        # 時間発展（|T1,T1>は変化しない）
        c_TT_new = c_TT
        c_S0_new = cos(phi) * c_S0 - 1j * sin(phi) * c_0S
        c_0S_new = -1j * sin(phi) * c_S0 + cos(phi) * c_0S

        # 更新
        state.set_coefficient(i=T1, j=T1, c_TT_new)
        state.set_coefficient(i=S1, j=S0, c_S0_new)
        state.set_coefficient(i=S0, j=S1, c_0S_new)

    return state
```

#### 7.3.4 $\hat{H}_{\text{rad}}$ の時間発展（簡略化）

完全な光子場を含む場合は非常に複雑になるため、実用的には **散逸過程として扱う** か、**Weisskopf-Wigner近似** を用いる。

**マルコフ近似による実効ハミルトニアン**:

$$
\hat{H}_{\text{eff}} = \hat{H}_{\text{rad}} - i\frac{\hbar \Gamma_{\text{fl}}}{2} \sum_i \hat{n}_{S_1}^{(i)}
$$

時間発展演算子:

$$
e^{-i\hat{H}_{\text{eff}} t/\hbar} \approx e^{-\Gamma_{\text{fl}} t / 2} \text{ on } |S_1\rangle \text{ states}
$$

**実装**:

```
function apply_radiative_evolution(state, dt):
    gamma_fl = radiation_decay_rate
    decay_factor = exp(-gamma_fl * dt / 2)

    for each basis_state in state.basis:
        count_S1 = number of S1 states in basis_state
        state.coefficients[basis_state] *= decay_factor^count_S1

    # 規格化
    norm = sqrt(sum(|c|^2 for c in state.coefficients))
    state.coefficients /= norm

    return state
```

注意: 厳密には量子ジャンプ法やリンドブラッド方程式を用いるべきだが、短時間では上記の近似が有効である。

### 7.4 完全な鈴木トロッター法の実装（2次対称）

**入力**:

- 初期状態: $|\Psi_0\rangle$
- パラメータ: $E_T, E_S, V_{ij}, J_{ij}, \Gamma_{\text{fl}}$
- 分子数: $N$
- 全時間: $T$
- 時間ステップ数: $N_{\text{steps}}$

**アルゴリズム**:

```
function suzuki_trotter_simulation_2and_order(Psi0, params, T, N_steps):
    dt = T / N_steps
    Psi = Psi0

    for step = 1 to N_steps:
        # 前半の対称分解
        Psi = apply_H0_evolution(Psi, dt/2)
        Psi = apply_transfer_evolution(Psi, dt/2)
        Psi = apply_TTA_evolution(Psi, dt/2)

        # 中央（放射項は全時間）
        Psi = apply_radiative_evolution(Psi, dt)

        # 後半の対称分解（逆順）
        Psi = apply_TTA_evolution(Psi, dt/2)
        Psi = apply_transfer_evolution(Psi, dt/2)
        Psi = apply_H0_evolution(Psi, dt/2)

    return Psi
```

### 7.5 観測量の計算

#### 7.5.1 個体数の計算

各時刻での各状態の個体数：

$$
N_{S_0}(t) = \sum_{i=1}^N \langle\Psi(t)| \hat{n}_{S_0}^{(i)} |\Psi(t)\rangle
$$

$$
N_{T_1}(t) = \sum_{i=1}^N \langle\Psi(t)| \hat{n}_{T_1}^{(i)} |\Psi(t)\rangle
$$

$$
N_{S_1}(t) = \sum_{i=1}^N \langle\Psi(t)| \hat{n}_{S_1}^{(i)} |\Psi(t)\rangle
$$

**実装**:

```
function calculate_populations(Psi):
    N_S0 = 0
    N_T1 = 0
    N_S1 = 0

    for each basis_state, coeff in Psi:
        prob = |coeff|^2
        for i = 1 to N:
            if basis_state[i] == S0:
                N_S0 += prob
            elif basis_state[i] == T1:
                N_T1 += prob
            elif basis_state[i] == S1:
                N_S1 += prob

    return N_S0, N_T1, N_S1
```

#### 7.5.2 蛍光強度の計算

瞬時蛍光強度：

$$
I_{\text{fl}}(t) = \Gamma_{\text{fl}} N_{S_1}(t)
$$

累積蛍光量：

$$
\mathcal{I}_{\text{total}} = \sum_{n=0}^{N_{\text{steps}}} I_{\text{fl}}(t_n) \Delta t
$$

## 8. 実装上の最適化

### 8.1 疎行列表現

多体系では状態空間の次元は $3^N$ と指数関数的に増大する。しかし、各時間発展演算子は疎（sparse）であるため、疎行列を用いることで計算を効率化できる。

### 8.2 対称性の利用

- **全粒子数保存**: $N_{\text{total}} = N_{S_0} + N_{T_1} + N_{S_1} = \text{const.}$
- **空間並進対称性**: 周期境界条件の場合、運動量が保存量
- **ブロック対角化**: 保存量ごとに状態空間を分割

### 8.3 並列化

各隣接ペアへの時間発展演算は独立なので、並列化が可能：

```
function apply_transfer_evolution_parallel(state, dt):
    Parallel for each pair (i, j) in adjacent_pairs:
        apply_local_evolution(state, i, j, dt)

    return state
```

### 8.4 適応的時間刻み

誤差評価に基づいて時間刻み幅 $\Delta t$ を動的に調整する：

```
function adaptive_time_step(Psi, H, dt_initial, tolerance):
    dt = dt_initial

    while time < T:
        # 時間刻み dt で1ステップ進める
        Psi1 = one_step(Psi, dt)

        # 時間刻み dt/2 で2ステップ進める
        Psi_half = one_step(Psi, dt/2)
        Psi2 = one_step(Psi_half, dt/2)

        # 誤差評価
        error = norm(Psi1 - Psi2)

        if error < tolerance:
            Psi = Psi2  # より精度の高い方を採用
            dt = dt * 1.5  # 時間刻みを増やす
        else:
            dt = dt * 0.5  # 時間刻みを減らす

    return Psi
```

## 9. 誤差解析と収束性

### 9.1 全体誤差の見積もり

$p$ 次の鈴木トロッター分解の全体誤差は：

$$
\varepsilon_{\text{total}} = C T \Delta t^p
$$

ここで、$C$ は系依存の定数であり、以下に依存する：

$$
C \sim \max_{m \neq n} \| [\hat{H}_m, \hat{H}_n] \|
$$

### 9.2 収束テスト

実際の計算では、$\Delta t$ を減少させながら結果が収束するかを確認する：

```
for N_steps in [100, 200, 400, 800, 1600]:
    dt = T / N_steps
    Psi_final = suzuki_trotter_simulation(Psi0, T, N_steps)
    observable = calculate_observable(Psi_final)
    print(N_steps, observable)
```

観測量が $\Delta t$ の冪乗則に従って収束することを確認：

$$
|\mathcal{O}(\Delta t) - \mathcal{O}_{\text{exact}}| \propto \Delta t^p
$$

### 9.3 ベンチマーク

小規模系（$N = 2$ 分子など）では、厳密対角化と比較：

$$
|\Psi_{\text{exact}}(T)\rangle = e^{-i\hat{H}_{\text{total}}T/\hbar} |\Psi_0\rangle
$$

フィデリティ:

$$
F = |\langle\Psi_{\text{exact}}(T) | \Psi_{\text{Trotter}}(T)\rangle|^2
$$

が1に近いことを確認する。

## 10. 高度な話題

### 10.1 Magnus展開との比較

Magnus展開は別の時間発展近似手法である：

$$
|\Psi(t)\rangle = \exp\left(-\frac{i}{\hbar} \Omega(t)\right) |\Psi(0)\rangle
$$

ここで、$\Omega(t)$ はMagnus級数で与えられる。鈴木トロッター分解と比較して、特定の系では優れた性能を示すことがある。

### 10.2 Krylov部分空間法

大規模系では、行列指数関数 $e^{-iHt}$ の計算にKrylov部分空間法（特にLanczos法）が有効である。

### 10.3 量子回路への変換

量子コンピュータで実装する場合、各 $e^{-i\hat{H}_m \Delta t/\hbar}$ を量子ゲートの列に分解する必要がある（Trotterization）。

## 11. プログラム実装例（疑似コード全体像）

```python
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply


class MolecularTripletSystem:
    def __init__(self, N, E_T, E_S, V, J, Gamma_fl):
        self.N = N  # 分子数
        self.E_T = E_T
        self.E_S = E_S
        self.V = V  # エネルギー移動積分
        self.J = J  # TTA相互作用
        self.Gamma_fl = Gamma_fl  # 蛍光放出速度

        # 状態空間の次元
        self.dim = 3**N

        # ハミルトニアン行列の構築（疎行列）
        self.H0 = self.build_H0()
        self.H_transfer = self.build_H_transfer()
        self.H_TTA = self.build_H_TTA()

    def build_H0(self):
        # 対角的なハミルトニアン H0 を構築
        H0 = csr_matrix((self.dim, self.dim), dtype=complex)
        # ... 実装 ...
        return H0

    def build_H_transfer(self):
        # エネルギー移動ハミルトニアンを構築
        H_transfer = csr_matrix((self.dim, self.dim), dtype=complex)
        # ... 実装 ...
        return H_transfer

    def build_H_TTA(self):
        # TTAハミルトニアンを構築
        H_TTA = csr_matrix((self.dim, self.dim), dtype=complex)
        # ... 実装 ...
        return H_TTA

    def apply_H0_evolution(self, psi, dt):
        # 対角項の時間発展（効率的）
        psi_new = psi.copy()
        # ... 対角位相の掛け算 ...
        return psi_new

    def apply_H_transfer_evolution(self, psi, dt):
        # 行列指数関数を疎行列で計算
        psi_new = expm_multiply(-1j * self.H_transfer * dt, psi)
        return psi_new

    def apply_H_TTA_evolution(self, psi, dt):
        psi_new = expm_multiply(-1j * self.H_TTA * dt, psi)
        return psi_new

    def apply_radiative_evolution(self, psi, dt):
        # 減衰項の適用
        decay_factor = np.exp(-self.Gamma_fl * dt / 2)
        psi_new = psi.copy()
        # S1状態に対して減衰を適用
        # ... 実装 ...
        # 規格化
        psi_new /= np.linalg.norm(psi_new)
        return psi_new

    def suzuki_trotter_2and_order(self, psi0, T, N_steps):
        dt = T / N_steps
        psi = psi0.copy()

        # 時間発展の記録
        times = []
        populations = []

        for step in range(N_steps):
            # 2次対称分解
            psi = self.apply_H0_evolution(psi, dt / 2)
            psi = self.apply_H_transfer_evolution(psi, dt / 2)
            psi = self.apply_H_TTA_evolution(psi, dt / 2)

            psi = self.apply_radiative_evolution(psi, dt)

            psi = self.apply_H_TTA_evolution(psi, dt / 2)
            psi = self.apply_H_transfer_evolution(psi, dt / 2)
            psi = self.apply_H0_evolution(psi, dt / 2)

            # 観測量の計算
            t = (step + 1) * dt
            N_S0, N_T1, N_S1 = self.calculate_populations(psi)
            times.append(t)
            populations.append((N_S0, N_T1, N_S1))

        return times, populations

    def calculate_populations(self, psi):
        # 各状態の個体数を計算
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0
        # ... 実装 ...
        return N_S0, N_T1, N_S1


# 使用例
if __name__ == "__main__":
    # パラメータ設定
    N_molecules = 10
    E_T = 1.5  # eV
    E_S = 3.0  # eV
    V = 0.1  # eV
    J = 0.05  # eV
    Gamma_fl = 1.0  # ns^-1

    # システム初期化
    system = MolecularTripletSystem(N_molecules, E_T, E_S, V, J, Gamma_fl)

    # 初期状態（すべての分子が三重項状態）
    psi0 = np.zeros(system.dim, dtype=complex)
    # ... 初期状態の設定 ...

    # シミュレーション実行
    T_final = 100.0  # ns
    N_steps = 1000
    times, populations = system.suzuki_trotter_2and_order(psi0, T_final, N_steps)

    # 結果の可視化
    import matplotlib.pyplot as plt

    N_S0_list = [p[0] for p in populations]
    N_T1_list = [p[1] for p in populations]
    N_S1_list = [p[2] for p in populations]

    plt.figure(figsize=(10, 6))
    plt.plot(times, N_S0_list, label="$N_{S_0}$")
    plt.plot(times, N_T1_list, label="$N_{T_1}$")
    plt.plot(times, N_S1_list, label="$N_{S_1}$")
    plt.xlabel("Time (ns)")
    plt.ylabel("Population")
    plt.legend()
    plt.title("Quantum Dynamics of Molecular Triplet States")
    plt.savefig("triplet_dynamics.png")
    plt.show()
```

## 12. まとめ

本文書では、分子三重項状態の量子ダイナミクスを鈴木トロッター分解を用いて数値的に解くための完全な理論とアルゴリズムを提示した。

### 12.1 主要な結果

1. **鈴木トロッター分解の理論的基礎**

   - 1次分解: $\mathcal{O}(\Delta t)$ の全体誤差
   - 2次対称分解: $\mathcal{O}(\Delta t^2)$ の全体誤差
   - 4次分解: $\mathcal{O}(\Delta t^4)$ の全体誤差

2. **ハミルトニアンの具体的分解**

   - $\hat{H}_0$: 対角的、位相の掛け算で計算
   - $\hat{H}_{\text{transfer}}$: 2サイト間の回転、解析的に計算可能
   - $\hat{H}_{\text{TTA}}$: 3状態間の混合、解析的に計算可能
   - $\hat{H}_{\text{rad}}$: 散逸過程として実効的に扱う

3. **実装可能なアルゴリズム**
   - 疑似コードレベルで各ステップを明示
   - 疎行列を用いた効率化
   - 並列化の可能性

### 12.2 今後の発展

- より高次の鈴木分解（6次、8次）
- 適応的時間刻み制御
- 量子-古典ハイブリッド法
- テンソルネットワーク法との組み合わせ

### 12.3 物理的意義

鈴木トロッター分解は、複雑な量子多体系のダイナミクスを、単純な演算子の積に分解することで、実用的な数値計算を可能にする。本手法は：

- 有機フォトニクスデバイスの設計
- 光アップコンバージョン効率の最適化
- 遅延蛍光材料の理論解析

などに応用できる。

## 参考文献

鈴木トロッター分解の理論と応用に関する文献：

1. **Trotter, H. F.** (1959). "On the product of semi-groups of operators". _Proc. Amer. Math. Soc._ **10**, 545-551.

2. **Suzuki, M.** (1990). "Fractal decomposition of exponential operators with applications to many-body theories and Monte Carlo simulations". _Phys. Lett. A_ **146**, 319-323.

3. **Suzuki, M.** (1991). "General theory of fractal path integrals with applications to many-body theories and statistical physics". _J. Math. Phys._ **32**, 400-407.

4. **Hatano, N. & Suzuki, M.** (2005). "Finding exponential product formulas of higher orders". _Lecture Notes in Physics_ **679**, 37-68.

5. **Lloyd, S.** (1996). "Universal Quantum Simulators". _Science_ **273**, 1073-1078.

6. **Children, A. M. & Su, Y.** (2019). "Nearly optimal lattice simulation by product formulas". _Phys. Rev. Lett._ **123**, 050503.

---

**文書作成日**: 2025-10-14
**分野**: 計算量子力学、数値シミュレーション、鈴木トロッター分解
**対象**: MQT Qudits フレームワークにおける量子ダイナミクスシミュレーション実装
