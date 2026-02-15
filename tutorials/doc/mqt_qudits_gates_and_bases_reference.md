# MQT-Qudits 量子ゲートと基底の完全リファレンス

## 目次

1. [概要](#概要)
2. [計算基底](#計算基底)
3. [単一Quditゲート](#単一quditゲート)
4. [2Quditゲート](#2quditゲート)
5. [多Quditゲート](#多quditゲート)
6. [カスタムゲート](#カスタムゲート)
7. [補助ゲート](#補助ゲート)

---

## 概要

MQT-Quditsは、量子ビット(qubit, d=2)を超えた高次元量子システム(qudit, d>2)を扱うためのフレームワークです。このドキュメントでは、フレームワークで利用可能なすべての量子ゲートと基底を、詳細な数式を用いて説明します。

### 表記規則

- $d$: Quditの次元
- $|i\rangle$: 計算基底状態 ($i = 0, 1, \ldots, d-1$)
- $\omega_d = e^{2\pi i/d}$: $d$次の単位根
- $\theta, \phi$: 回転角度パラメータ
- $a, b$: エネルギーレベル ($0 \leq a < b < d$)

---

## 計算基底

### 標準計算基底

$d$次元quditの計算基底は、以下の$d$個の正規直交状態で構成されます：

$$
\{|0\rangle, |1\rangle, |2\rangle, \ldots, |d-1\rangle\}
$$

各基底状態はベクトル表現で以下のように表されます：

$$
|0\rangle = \begin{pmatrix} 1 \\ 0 \\ \vdots \\ 0 \end{pmatrix}, \quad
|1\rangle = \begin{pmatrix} 0 \\ 1 \\ \vdots \\ 0 \end{pmatrix}, \quad \ldots, \quad
|d-1\rangle = \begin{pmatrix} 0 \\ 0 \\ \vdots \\ 1 \end{pmatrix}
$$

### 一般的な状態

一般的な$d$次元qudit状態は、計算基底の線形結合として表されます：

$$
|\psi\rangle = \sum_{i=0}^{d-1} \alpha_i |i\rangle
$$

ここで、$\alpha_i \in \mathbb{C}$は確率振幅であり、規格化条件$\sum_{i=0}^{d-1} |\alpha_i|^2 = 1$を満たします。

---

## 単一Quditゲート

### 1. 一般化Hadamardゲート (H)

**説明**: qubitのHadamardゲートを$d$次元に一般化したゲートです。すべての計算基底状態を均等な重ね合わせ状態に変換します。

**行列表現**:

$$
H_d = \frac{1}{\sqrt{d}} \sum_{j=0}^{d-1} \sum_{k=0}^{d-1} \omega_d^{jk} |j\rangle\langle k|
$$

展開形式：

$$
H_d = \frac{1}{\sqrt{d}} \begin{pmatrix}
1 & 1 & 1 & \cdots & 1 \\
1 & \omega_d & \omega_d^2 & \cdots & \omega_d^{d-1} \\
1 & \omega_d^2 & \omega_d^4 & \cdots & \omega_d^{2(d-1)} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & \omega_d^{d-1} & \omega_d^{2(d-1)} & \cdots & \omega_d^{(d-1)^2}
\end{pmatrix}
$$

ここで、$\omega_d = e^{2\pi i/d}$は$d$次の単位根です。

**性質**:

- ユニタリ性: $H_d^\dagger H_d = I_d$
- 自己随伴: $H_d^\dagger = H_d$
- 周期性: $H_d^2 = I_d$ (次元が2の場合)

**実装コード**:

```python
def H(dimension):
    matrix = np.zeros((dimension, dimension), dtype=complex)
    for e0 in range(dimension):
        for e1 in range(dimension):
            omega = np.exp(2j * np.pi * e0 * e1 / dimension)
            matrix[e0, e1] = omega
    return matrix / np.sqrt(dimension)
```

**例 (d=3)**:

$$
H_3 = \frac{1}{\sqrt{3}} \begin{pmatrix}
1 & 1 & 1 \\
1 & e^{2\pi i/3} & e^{4\pi i/3} \\
1 & e^{4\pi i/3} & e^{8\pi i/3}
\end{pmatrix}
$$

---

### 2. 一般化Pauli-Xゲート (X)

**説明**: 計算基底状態を巡回的にシフトさせるゲートです。qubitの場合は$|0\rangle \leftrightarrow |1\rangle$の反転ですが、quditでは$|i\rangle \rightarrow |i+1 \bmod d\rangle$の巡回シフトになります。

**行列表現**:

$$
X_d = \sum_{i=0}^{d-1} |i+1 \bmod d\rangle\langle i|
$$

展開形式：

$$
X_d = \begin{pmatrix}
0 & 0 & 0 & \cdots & 0 & 1 \\
1 & 0 & 0 & \cdots & 0 & 0 \\
0 & 1 & 0 & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & 1 & 0
\end{pmatrix}
$$

**性質**:

- 巡回性: $X_d^d = I_d$
- ユニタリ性: $X_d^\dagger X_d = I_d$
- 固有値: $\{e^{2\pi i k/d} \mid k=0,1,\ldots,d-1\}$

**作用例**:

$$
\begin{align}
X_d |0\rangle &= |1\rangle \\
X_d |1\rangle &= |2\rangle \\
&\vdots \\
X_d |d-1\rangle &= |0\rangle
\end{align}
$$

**例 (d=3)**:

$$
X_3 = \begin{pmatrix}
0 & 0 & 1 \\
1 & 0 & 0 \\
0 & 1 & 0
\end{pmatrix}
$$

---

### 3. 一般化Pauli-Zゲート (Z)

**説明**: 各計算基底状態に位相をかけるゲートです。状態$|i\rangle$に対して位相$e^{2\pi i \cdot i/d}$を乗じます。

**行列表現**:

$$
Z_d = \sum_{i=0}^{d-1} \omega_d^i |i\rangle\langle i|
$$

展開形式：

$$
Z_d = \begin{pmatrix}
1 & 0 & 0 & \cdots & 0 \\
0 & \omega_d & 0 & \cdots & 0 \\
0 & 0 & \omega_d^2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \omega_d^{d-1}
\end{pmatrix}
$$

ここで、$\omega_d = e^{2\pi i/d}$です。

**性質**:

- 対角: すべてが対角行列
- 巡回性: $Z_d^d = I_d$
- $X$ゲートとの関係: $Z_d X_d Z_d^\dagger = \omega_d X_d$

**作用例**:

$$
Z_d |i\rangle = e^{2\pi i \cdot i/d} |i\rangle
$$

**例 (d=3)**:

$$
Z_3 = \begin{pmatrix}
1 & 0 & 0 \\
0 & e^{2\pi i/3} & 0 \\
0 & 0 & e^{4\pi i/3}
\end{pmatrix}
$$

---

### 4. Sゲート (S)

**説明**: 素数次元のquditに対して定義される位相ゲートです。$d=2$の場合は標準的なSゲート（$\pi/2$位相ゲート）となります。

**行列表現**:

$$
S_d = \sum_{i=0}^{d-1} \omega_d^{i(i+1)/2} |i\rangle\langle i|
$$

**d=2の場合**:

$$
S_2 = \begin{pmatrix}
1 & 0 \\
0 & i
\end{pmatrix}
$$

**一般のd（素数）の場合**:

$$
S_d = \begin{pmatrix}
1 & 0 & 0 & \cdots & 0 \\
0 & \omega_d & 0 & \cdots & 0 \\
0 & 0 & \omega_d^3 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \omega_d^{d(d-1)/2}
\end{pmatrix}
$$

**条件**: 次元$d$は素数でなければなりません。

**例 (d=3)**:

$$
S_3 = \begin{pmatrix}
1 & 0 & 0 \\
0 & e^{2\pi i/3} & 0 \\
0 & 0 & e^{4\pi i/3}
\end{pmatrix} = \begin{pmatrix}
1 & 0 & 0 \\
0 & \omega_3 & 0 \\
0 & 0 & \omega_3^3
\end{pmatrix}
$$

---

### 5. 一般化回転ゲート (R)

**説明**: 指定された2つのエネルギーレベル$a$と$b$の間で回転を行うゲートです。Gell-Mann行列を生成子として使用します。

**パラメータ**:

- $a, b$: エネルギーレベル ($0 \leq a < b < d$)
- $\theta$: 回転角度
- $\phi$: 回転軸の方位角

**行列表現**:

$$
R_{a,b}(\theta, \phi) = \cos\left(\frac{\theta}{2}\right) I - i\sin\left(\frac{\theta}{2}\right) \left[\sin(\phi) \lambda_{a,b}^{(a)} + \cos(\phi) \lambda_{a,b}^{(s)}\right]
$$

ここで、$\lambda_{a,b}^{(s)}$と$\lambda_{a,b}^{(a)}$はGell-Mann行列です：

**対称Gell-Mann行列** $\lambda_{a,b}^{(s)}$:

$$
\lambda_{a,b}^{(s)} = |a\rangle\langle b| + |b\rangle\langle a|
$$

**反対称Gell-Mann行列** $\lambda_{a,b}^{(a)}$:

$$
\lambda_{a,b}^{(a)} = -i|a\rangle\langle b| + i|b\rangle\langle a|
$$

**具体的な行列要素**:

$$
R_{a,b}(\theta, \phi)_{ij} = \begin{cases}
\cos(\theta/2) & \text{if } i = j \notin \{a, b\} \\
\cos(\theta/2) & \text{if } i = j \in \{a, b\} \\
-i\sin(\theta/2)[\sin(\phi)(-i) + \cos(\phi)] = -i\sin(\theta/2)e^{-i\phi} & \text{if } (i,j) = (a, b) \\
-i\sin(\theta/2)[\sin(\phi)(i) + \cos(\phi)] = -i\sin(\theta/2)e^{i\phi} & \text{if } (i,j) = (b, a) \\
0 & \text{otherwise}
\end{cases}
$$

**性質**:

- $R_{a,b}(0, \phi) = I$
- $R_{a,b}(2\pi, \phi) = -I$
- $R_{a,b}(\theta, 0)$は$X$タイプの回転
- $R_{a,b}(\theta, \pi/2)$は$Y$タイプの回転

**例 (d=3, a=0, b=1)**:

$$
R_{0,1}(\theta, \phi) = \begin{pmatrix}
\cos(\theta/2) & -i\sin(\theta/2)e^{-i\phi} & 0 \\
-i\sin(\theta/2)e^{i\phi} & \cos(\theta/2) & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

---

### 6. Z軸回転ゲート (Rz)

**説明**: 2つのエネルギーレベル間でZ軸周りの回転を実行するゲートです。3つのRゲートの合成で実装されます。

**パラメータ**:

- $a, b$: エネルギーレベル ($0 \leq a < b < d$)
- $\phi$: 回転角度

**行列表現**:

Rzゲートは以下の3つのRゲートの積として定義されます：

$$
\text{Rz}_{a,b}(\phi) = R_{a,b}(\pi/2, 0) \cdot R_{a,b}(\phi, \pi/2) \cdot R_{a,b}(-\pi/2, 0)
$$

これは実質的に以下の位相回転を実行します：

$$
\text{Rz}_{a,b}(\phi) = \exp\left(-i\frac{\phi}{2}(|a\rangle\langle a| - |b\rangle\langle b|)\right)
$$

展開すると：

$$
\text{Rz}_{a,b}(\phi) = I + (e^{-i\phi/2} - 1)|a\rangle\langle a| + (e^{i\phi/2} - 1)|b\rangle\langle b|
$$

**具体的な行列要素**:

$$
\text{Rz}_{a,b}(\phi)_{ij} = \begin{cases}
1 & \text{if } i = j \notin \{a, b\} \\
e^{-i\phi/2} & \text{if } i = j = a \\
e^{i\phi/2} & \text{if } i = j = b \\
0 & \text{if } i \neq j
\end{cases}
$$

**性質**:

- 対角行列
- $\text{Rz}_{a,b}(0) = I$
- $\text{Rz}_{a,b}(2\pi) = -I$

**例 (d=3, a=0, b=1)**:

$$
\text{Rz}_{0,1}(\phi) = \begin{pmatrix}
e^{-i\phi/2} & 0 & 0 \\
0 & e^{i\phi/2} & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

---

### 7. Hadamard型回転ゲート (Rh)

**説明**: 2つのエネルギーレベル間でHadamard変換を実行するゲートです。SU(2)部分空間でのHadamardゲートに相当します。

**パラメータ**:

- $a, b$: エネルギーレベル ($0 \leq a < b < d$)

**行列表現**:

Rhゲートは2つのRゲートの積として定義されます：

$$
\text{Rh}_{a,b} = R_{a,b}(-\pi, 0) \cdot R_{a,b}(\pi/2, \pi/2)
$$

これは実質的に以下の変換を実行します：

$$
\text{Rh}_{a,b} = \frac{1}{\sqrt{2}} \begin{pmatrix}
1 & 1 \\
1 & -1
\end{pmatrix} \text{ on subspace } \{|a\rangle, |b\rangle\}
$$

**具体的な行列要素** (d次元空間で):

$$
\text{Rh}_{a,b|ij} = \begin{cases}
1 & \text{if } i = j \notin \{a, b\} \\
1/\sqrt{2} & \text{if } i = a, j = a \\
1/\sqrt{2} & \text{if } i = a, j = b \\
1/\sqrt{2} & \text{if } i = b, j = a \\
-1/\sqrt{2} & \text{if } i = b, j = b \\
0 & \text{otherwise}
\end{cases}
$$

**作用例**:

$$
\begin{align}
\text{Rh}_{a,b} |a\rangle &= \frac{1}{\sqrt{2}}(|a\rangle + |b\rangle) \\
\text{Rh}_{a,b} |b\rangle &= \frac{1}{\sqrt{2}}(|a\rangle - |b\rangle) \\
\text{Rh}_{a,b} |k\rangle &= |k\rangle \quad \text{for } k \notin \{a, b\}
\end{align}
$$

**例 (d=3, a=0, b=1)**:

$$
\text{Rh}_{0,1} = \begin{pmatrix}
1/\sqrt{2} & 1/\sqrt{2} & 0 \\
1/\sqrt{2} & -1/\sqrt{2} & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

---

### 8. 仮想Z回転ゲート (VirtRz)

**説明**: 単一のエネルギーレベルに対して位相を適用するゲートです。計算には含まれない「仮想」の位相回転として扱われることがあります。

**パラメータ**:

- $a$: エネルギーレベル ($0 \leq a < d$)
- $\phi$: 位相角度

**行列表現**:

$$
\text{VirtRz}_a(\phi) = \sum_{i=0}^{d-1} \delta_{ia} e^{-i\phi} |i\rangle\langle i| + \sum_{i \neq a} |i\rangle\langle i|
$$

展開形式：

$$
\text{VirtRz}_a(\phi) = I + (e^{-i\phi} - 1)|a\rangle\langle a|
$$

**具体的な行列要素**:

$$
\text{VirtRz}_a(\phi)_{ij} = \begin{cases}
e^{-i\phi} & \text{if } i = j = a \\
1 & \text{if } i = j \neq a \\
0 & \text{if } i \neq j
\end{cases}
$$

**性質**:

- 対角行列
- グローバル位相を除いて効果的
- $\text{VirtRz}_a(0) = I$

**例 (d=3, a=1)**:

$$
\text{VirtRz}_1(\phi) = \begin{pmatrix}
1 & 0 & 0 \\
0 & e^{-i\phi} & 0 \\
0 & 0 & 1
\end{pmatrix}
$$

---

### 9. 置換ゲート (Perm)

**説明**: 計算基底状態を任意の順序で入れ替えるゲートです。指定された置換に従って基底状態を並び替えます。

**パラメータ**:

- $\sigma$: 置換 (長さ$d$のリスト)

**行列表現**:

置換$\sigma = [\sigma(0), \sigma(1), \ldots, \sigma(d-1)]$に対して：

$$
\text{Perm}_\sigma = \sum_{i=0}^{d-1} |\sigma(i)\rangle\langle i|
$$

行列としては、単位行列の列を置換に従って並び替えたものになります：

$$
\text{Perm}_\sigma = [e_{\sigma(0)}, e_{\sigma(1)}, \ldots, e_{\sigma(d-1)}]
$$

ここで、$e_i$は$i$番目の標準基底ベクトルです。

**性質**:

- 置換行列
- $\text{Perm}_\sigma^{-1} = \text{Perm}_{\sigma^{-1}}$
- $\text{Perm}_\sigma^\dagger = \text{Perm}_{\sigma^{-1}}$

**作用例**:

$$
\text{Perm}_\sigma |i\rangle = |\sigma(i)\rangle
$$

**例 (d=3, σ=[2,0,1])**:

$$
\text{Perm}_{[2,0,1]} = \begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{pmatrix}
$$

この置換は $|0\rangle \to |2\rangle$, $|1\rangle \to |0\rangle$, $|2\rangle \to |1\rangle$ を実行します。

---

### 10. ノイズXゲート (NoiseX)

**説明**: 2つの指定されたエネルギーレベル間でビット反転（X回転）を行うゲートです。qubitのPauli-Xゲートの一般化です。

**パラメータ**:

- $a, b$: エネルギーレベル ($0 \leq a < b < d$)

**行列表現**:

$$
\text{NoiseX}_{a,b} = \sum_{i \notin \{a,b\}} |i\rangle\langle i| + |a\rangle\langle b| + |b\rangle\langle a|
$$

**具体的な行列要素**:

$$
\text{NoiseX}_{a,b|ij} = \begin{cases}
1 & \text{if } i = j \notin \{a, b\} \\
0 & \text{if } i = j \in \{a, b\} \\
1 & \text{if } (i, j) = (a, b) \text{ or } (i, j) = (b, a) \\
0 & \text{otherwise}
\end{cases}
$$

**作用例**:

$$
\begin{align}
\text{NoiseX}_{a,b} |a\rangle &= |b\rangle \\
\text{NoiseX}_{a,b} |b\rangle &= |a\rangle \\
\text{NoiseX}_{a,b} |k\rangle &= |k\rangle \quad \text{for } k \notin \{a, b\}
\end{align}
$$

**性質**:

- $\text{NoiseX}_{a,b}^2 = I$
- エルミート: $\text{NoiseX}_{a,b}^\dagger = \text{NoiseX}_{a,b}$
- ユニタリ: $\text{NoiseX}_{a,b}^\dagger \text{NoiseX}_{a,b} = I$

**例 (d=3, a=0, b=2)**:

$$
\text{NoiseX}_{0,2} = \begin{pmatrix}
0 & 0 & 1 \\
0 & 1 & 0 \\
1 & 0 & 0
\end{pmatrix}
$$

---

### 11. ノイズYゲート (NoiseY)

**説明**: 2つの指定されたエネルギーレベル間でY回転を行うゲートです。qubitのPauli-Yゲートの一般化です。

**パラメータ**:

- $a, b$: エネルギーレベル ($0 \leq a < b < d$)

**行列表現**:

$$
\text{NoiseY}_{a,b} = \sum_{i \notin \{a,b\}} |i\rangle\langle i| - i|a\rangle\langle b| + i|b\rangle\langle a|
$$

**具体的な行列要素**:

$$
\text{NoiseY}_{a,b|ij} = \begin{cases}
1 & \text{if } i = j \notin \{a, b\} \\
0 & \text{if } i = j \in \{a, b\} \\
-i & \text{if } (i, j) = (a, b) \\
i & \text{if } (i, j) = (b, a) \\
0 & \text{otherwise}
\end{cases}
$$

**作用例**:

$$
\begin{align}
\text{NoiseY}_{a,b} |a\rangle &= -i|b\rangle \\
\text{NoiseY}_{a,b} |b\rangle &= i|a\rangle \\
\text{NoiseY}_{a,b} |k\rangle &= |k\rangle \quad \text{for } k \notin \{a, b\}
\end{align}
$$

**性質**:

- $\text{NoiseY}_{a,b}^2 = I$
- エルミート: $\text{NoiseY}_{a,b}^\dagger = \text{NoiseY}_{a,b}$
- ユニタリ: $\text{NoiseY}_{a,b}^\dagger \text{NoiseY}_{a,b} = I$
- XおよびZ型ゲートとの関係: $\text{NoiseY}_{a,b} = -i \text{NoiseX}_{a,b} \text{Rz}_{a,b}(\pi)$

**例 (d=3, a=0, b=2)**:

$$
\text{NoiseY}_{0,2} = \begin{pmatrix}
0 & 0 & -i \\
0 & 1 & 0 \\
i & 0 & 0
\end{pmatrix}
$$

---

## 2Quditゲート

### 12. 制御Exchangeゲート (CEx)

**説明**: 制御quditが特定の状態にあるとき、ターゲットquditの2つのレベル間でスワップ（または位相付きスワップ）を実行するゲートです。

**パラメータ**:

- $a, b$: ターゲットquditのスワップするレベル ($0 \leq a < b < d_{\text{target}}$)
- $c$: 制御レベル ($0 \leq c < d_{\text{ctrl}}$)
- $\phi$: 位相角度

**行列表現**:

制御quditが状態$|c\rangle$にあるとき、ターゲットquditに以下の変換を適用します：

$$
U_{a,b}(\phi) = I + (\text{swap\_op} - I) \cdot |c\rangle\langle c|_{\text{ctrl}}
$$

ここで、スワップ演算子は：

$$
\text{swap\_op} = \begin{pmatrix}
\ddots & & & & \\
& 0 & -ie^{i\phi} & & \\
& -ie^{-i\phi} & 0 & & \\
& & & \ddots &
\end{pmatrix}
$$

位置$(a,a), (a,b), (b,a), (b,b)$に以下の要素を持ちます：

- $(a,a)$: $0$
- $(a,b)$: $-i\cos(\phi) - \sin(\phi)$
- $(b,a)$: $-i\cos(\phi) + \sin(\phi)$
- $(b,b)$: $0$

**完全な行列** (制御次元$d_c$、ターゲット次元$d_t$):

$$
\text{CEx}_{a,b,c}(\phi) = \sum_{i \neq c} |i\rangle\langle i|_{\text{ctrl}} \otimes I_{d_t} + |c\rangle\langle c|_{\text{ctrl}} \otimes U_{a,b}(\phi)
$$

**例 (制御:d=2, ターゲット:d=2, a=0, b=1, c=1, φ=0)**:

$$
\text{CEx}_{0,1,1}(0) = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0 & -i \\
0 & 0 & -i & 0
\end{pmatrix}
$$

これは標準的な制御NOTゲート（CNOT）の変形です。

---

### 13. 制御Sumゲート (CSum)

**説明**: 制御quditの値に応じて、ターゲットquditにX演算を複数回適用するゲートです。加算modulo $d$を実行します。

**パラメータ**:

- 制御quditとターゲットqudit（それぞれ次元$d_{\text{ctrl}}$と$d_{\text{target}}$）

**行列表現**:

$$
\text{CSum} = \sum_{i=0}^{d_{\text{ctrl}}-1} |i\rangle\langle i|_{\text{ctrl}} \otimes X_{\text{target}}^i
$$

ここで、$X_{\text{target}}$はターゲットquditの一般化Pauli-Xゲートです。

**作用**:

$$
\text{CSum} |i\rangle_{\text{ctrl}} |j\rangle_{\text{target}} = |i\rangle_{\text{ctrl}} |(j+i) \bmod d_{\text{target}}\rangle_{\text{target}}
$$

**性質**:

- 量子加算を実行
- 可逆演算: $\text{CSum}^{d_{\text{target}}} = I$

**例 (制御:d=3, ターゲット:d=3)**:

$$
\text{CSum}_{3,3} = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0
\end{pmatrix}
$$

基底順序は$|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle$です。

---

### 14. LSゲート (Local Spin)

**説明**: 2つのqudit間でローカルスピン相互作用を実装するゲートです。対角要素が等しい状態（$|ii\rangle$）に対して作用します。

**パラメータ**:

- $\theta$: 相互作用の強さ

**行列表現**:

射影演算子を定義します：

$$
P = \sum_{i=0}^{\min(d_0, d_1)-1} |ii\rangle\langle ii|
$$

LSゲートは以下のように定義されます：

$$
\text{LS}(\theta) = e^{-i\theta P}
$$

展開すると：

$$
\text{LS}(\theta) = I + (e^{-i\theta} - 1) P
$$

**作用**:

$$
\text{LS}(\theta) |ij\rangle = \begin{cases}
e^{-i\theta} |ij\rangle & \text{if } i = j < \min(d_0, d_1) \\
|ij\rangle & \text{otherwise}
\end{cases}
$$

**例 (d₀=3, d₁=3)**:

$$
\text{LS}(\theta) = \begin{pmatrix}
e^{-i\theta} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & e^{-i\theta} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & e^{-i\theta}
\end{pmatrix}
$$

---

### 15. MSゲート (Mølmer-Sørensen)

**説明**: イオントラップ量子コンピュータで使用されるMølmer-Sørensenゲートの一般化です。2つのqudit間でもつれを生成します。

**パラメータ**:

- $\theta$: 相互作用の強さ

**行列表現**:

Gell-Mann行列$\lambda_{0,1}^{(s)}$を用いて、ハミルトニアンを定義します：

$$
H_{\text{MS}} = (\sigma_0 + \sigma_1)^2
$$

ここで、$\sigma_i = \lambda_{0,1}^{(s)}$は各quditに作用します。

MSゲートは以下のように定義されます：

$$
\text{MS}(\theta) = \exp\left(-i\frac{\theta}{4}(\sigma_0 + \sigma_1)^2\right)
$$

展開すると：

$$
\text{MS}(\theta) = \exp\left(-i\frac{\theta}{4}(
\sigma_0 \otimes \sigma_0 +
\sigma_0 \otimes I +
I \otimes \sigma_0 +
I \otimes I
)\right)
$$

より具体的には：

$$
(\sigma_0 \otimes I + I \otimes \sigma_1)^2 = \sigma_0 \otimes \sigma_0 + \sigma_0 \otimes I + I \otimes \sigma_1 + \text{cross terms}
$$

**性質**:

- もつれゲート
- $\theta = \pi/2$で最大もつれ状態を生成

**例** ($d_0=d_1=2$の場合、$\theta=\pi/2$):

MSゲートは以下のようなベル状態を生成します：

$$
\text{MS}(\pi/2) |00\rangle \approx \frac{1}{\sqrt{2}}(|00\rangle + i|11\rangle)
$$

---

## 多Quditゲート

### 16. ランダムユニタリゲート (RandU)

**説明**: 複数のquditに作用するランダムなユニタリ行列を生成するゲートです。Haar測度に従って一様にサンプリングされます。

**パラメータ**:

- 対象quditのリスト（次元のリスト）

**行列表現**:

$n$個のquditに対して、全体の次元は$D = \prod_{i=1}^n d_i$となります。

$$
\text{RandU} = U \in U(D)
$$

ここで、$U$はHaar測度に従ってランダムにサンプリングされた$D \times D$ユニタリ行列です。

**性質**:

- ユニタリ性: $U^\dagger U = I$
- Haar測度に従う一様分布
- 完全にランダムな量子演算

**サンプリング方法**:

Haar測度からのサンプリングは以下の手順で行われます：

1. $D \times D$の複素ガウス行列$G$を生成（各要素が独立に標準正規分布に従う）
2. QR分解: $G = QR$
3. $Q$を対角行列で正規化して$U$を得る

**使用例**:

- ランダム量子回路のベンチマーク
- 量子カオスの研究
- ランダム化コンパイル技術

---

## カスタムゲート

### 17. 単一Quditカスタムゲート (CustomOne)

**説明**: ユーザーが指定した任意のユニタリ行列を単一のquditに適用するゲートです。

**パラメータ**:

- $U$: $d \times d$のユニタリ行列

**行列表現**:

$$
\text{CustomOne}(U) = U
$$

ここで、$U \in U(d)$はユーザーが指定した$d \times d$のユニタリ行列です。

**制約**:

- $U$はユニタリ行列でなければならない: $U^\dagger U = I_d$

**使用例**:

```python
import numpy as np
from mqt.qudits import QuantumCircuit

# 3次元quditの回路を作成
qc = QuantumCircuit(1, [3])

# カスタムユニタリ行列を定義
U = np.array(
    [
        [1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)],
        [
            1 / np.sqrt(3),
            np.exp(2j * np.pi / 3) / np.sqrt(3),
            np.exp(4j * np.pi / 3) / np.sqrt(3),
        ],
        [
            1 / np.sqrt(3),
            np.exp(4j * np.pi / 3) / np.sqrt(3),
            np.exp(8j * np.pi / 3) / np.sqrt(3),
        ],
    ]
)

# カスタムゲートを適用
qc.cu_one(0, U)
```

---

### 18. 2Quditカスタムゲート (CustomTwo)

**説明**: ユーザーが指定した任意のユニタリ行列を2つのquditに適用するゲートです。

**パラメータ**:

- $U$: $(d_0 \cdot d_1) \times (d_0 \cdot d_1)$のユニタリ行列

**行列表現**:

$$
\text{CustomTwo}(U) = U
$$

ここで、$U \in U(d_0 \cdot d_1)$はユーザーが指定したユニタリ行列です。

**制約**:

- $U$はユニタリ行列でなければならない: $U^\dagger U = I_{d_0 \cdot d_1}$
- 行列のサイズは2つのquditの次元の積と一致する必要がある

**使用例**:

```python
import numpy as np
from mqt.qudits import QuantumCircuit

# 2つの3次元quditの回路を作成
qc = QuantumCircuit(2, [3, 3])

# カスタム2quditユニタリ行列を定義（9x9行列）
# 例：制御スワップゲート
U = np.eye(9, dtype=complex)
# ... 行列要素を設定 ...

# カスタムゲートを適用
qc.cu_two([0, 1], U)
```

---

### 19. 多Quditカスタムゲート (CustomMulti)

**説明**: ユーザーが指定した任意のユニタリ行列を複数のquditに適用するゲートです。

**パラメータ**:

- $U$: $D \times D$のユニタリ行列（$D = \prod_i d_i$）

**行列表現**:

$$
\text{CustomMulti}(U) = U
$$

ここで、$U \in U(D)$はユーザーが指定したユニタリ行列で、$D$は全quditの次元の積です。

**制約**:

- $U$はユニタリ行列でなければならない: $U^\dagger U = I_D$
- 行列のサイズはすべての対象quditの次元の積と一致する必要がある

**使用例**:

```python
import numpy as np
from mqt.qudits import QuantumCircuit

# 3つのquditの回路を作成
qc = QuantumCircuit(3, [2, 3, 2])

# カスタム多quditユニタリ行列を定義（12x12行列、2*3*2=12）
U = np.eye(12, dtype=complex)
# ... 行列要素を設定 ...

# カスタムゲートを適用
qc.cu_multi([0, 1, 2], U)
```

---

## 補助ゲート

### 20. Gell-Mann行列 (GellMann)

**説明**: SU(d)群の生成子であるGell-Mann行列です。回転ゲートの構成要素として使用されます。

**種類**:

1. **対称型** ($\lambda_{a,b}^{(s)}$):

   $$
   \lambda_{a,b}^{(s)} = |a\rangle\langle b| + |b\rangle\langle a|
   $$

2. **反対称型** ($\lambda_{a,b}^{(a)}$):

   $$
   \lambda_{a,b}^{(a)} = -i|a\rangle\langle b| + i|b\rangle\langle a|
   $$

3. **対角型** ($\lambda_k^{(d)}$):
   $$
   \lambda_k^{(d)} = \sqrt{\frac{2}{k(k+1)}} \left(\sum_{j=0}^{k-1} |j\rangle\langle j| - k|k\rangle\langle k|\right)
   $$

**性質**:

- エルミート: $\lambda^\dagger = \lambda$
- トレースレス: $\text{Tr}(\lambda) = 0$
- 正規化: $\text{Tr}(\lambda_i \lambda_j) = 2\delta_{ij}$

**例** (d=3の場合):

$$
\lambda_{0,1}^{(s)} = \begin{pmatrix}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}, \quad
\lambda_{0,1}^{(a)} = \begin{pmatrix}
0 & -i & 0 \\
i & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}
$$

$$
\lambda_1^{(d)} = \sqrt{\frac{2}{2}} \begin{pmatrix}
1 & 0 & 0 \\
0 & -1 & 0 \\
0 & 0 & 0
\end{pmatrix}, \quad
\lambda_2^{(d)} = \sqrt{\frac{2}{6}} \begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & -2
\end{pmatrix}
$$

**用途**:

- 任意の回転ゲートの生成子
- 量子状態トモグラフィー
- quditシステムの完全な記述

---

## まとめ

このドキュメントでは、MQT-Quditsフレームワークで利用可能なすべての量子ゲートと基底を詳細に説明しました。

### ゲートの分類

**単一Quditゲート** (11種類):

1. 一般化Hadamardゲート (H)
2. 一般化Pauli-Xゲート (X)
3. 一般化Pauli-Zゲート (Z)
4. Sゲート (S)
5. 一般化回転ゲート (R)
6. Z軸回転ゲート (Rz)
7. Hadamard型回転ゲート (Rh)
8. 仮想Z回転ゲート (VirtRz)
9. 置換ゲート (Perm)
10. ノイズXゲート (NoiseX)
11. ノイズYゲート (NoiseY)

**2Quditゲート** (4種類):

1. 制御Exchangeゲート (CEx)
2. 制御Sumゲート (CSum)
3. LSゲート (Local Spin)
4. MSゲート (Mølmer-Sørensen)

**多Quditゲート** (1種類):

1. ランダムユニタリゲート (RandU)

**カスタムゲート** (3種類):

1. 単一Quditカスタムゲート (CustomOne)
2. 2Quditカスタムゲート (CustomTwo)
3. 多Quditカスタムゲート (CustomMulti)

**補助** (1種類):

1. Gell-Mann行列 (GellMann)

### 使用例

```python
from mqt.qudits import QuantumCircuit

# 3次元quditの回路を作成
qc = QuantumCircuit(3, [3, 3, 3])

# 各種ゲートを適用
qc.h(0)  # Hadamardゲート
qc.x(1)  # Pauli-Xゲート
qc.z(2)  # Pauli-Zゲート
qc.r(0, [0, 1, np.pi / 4, 0])  # 回転ゲート
qc.rz(1, [0, 2, np.pi / 2])  # Rz回転
qc.cx([0, 1], [0, 1, 1, 0.0])  # 制御Exchangeゲート
qc.csum([1, 2])  # 制御Sumゲート
qc.ls([0, 1], [np.pi / 4])  # LSゲート
qc.ms([1, 2], [np.pi / 2])  # MSゲート
```

### 参考文献

- Wang, Y., Hu, Z., Sanders, B. C., & Kais, S. (2020). Qudits and high-dimensional quantum computing. Frontiers in Physics, 8, 589504.
- Gokhale, P., Baker, J. M., Duckering, C., Brown, N. C., Brown, K. R., & Chong, F. T. (2019). Asymptotic improvements to quantum circuits via qutrits. In Proceedings of the 46th International Symposium on Computer Architecture (pp. 554-566).
- MQT Qudits Documentation: https://github.com/cda-tum/mqt-qudits

---

**ドキュメント作成日**: 2025年10月15日
**バージョン**: 1.0
**フレームワーク**: MQT-Qudits
