# シナリオ5: Qudit GKSL（ボソン無し）の省略無し詳細理論 — Part 2

## Stinespring拡張の理論：Lindbladチャネルの量子回路表現

本文書は Part 1 に続き、シナリオ5で用いられる Stinespring 拡張（Stinespring dilation）の理論を省略無しに定式化する。各Lindbladチャネルをユニタリ演算子として実装するための数学的構成法を、行列の具体的な形まで記述する。

---

## 1. Stinespring拡張定理の概要

### 1.1 背景

GKSL方程式の散逸部分は、個別のLindblad散逸子 $\mathcal{D}[\tilde{L}_\alpha]$ の和として構成される。各散逸子は完全正・トレース保存（CPTP）写像の生成子であるが、そのまま量子回路として実装することはできない。量子コンピュータ上で実装するためには、各CPTP写像をユニタリ演算子に「持ち上げる」必要がある。これが Stinespring 拡張定理の役割である。

### 1.2 定理の内容（本コードで用いる形式）

系 $\mathcal{H}_S$（$d_{\mathrm{sys}}$ 次元）に対するCPTP写像 $\mathcal{E}$ は、適当な補助系（アンシラ）$\mathcal{H}_E$（$d_{\mathrm{anc}}$ 次元）を導入し、拡大空間 $\mathcal{H}_E \otimes \mathcal{H}_S$ 上のユニタリ演算子 $U$ を用いて以下のように表現できる：

$$
\mathcal{E}(\rho) = \mathrm{Tr}_E\!\left[U\left(|0\rangle\langle 0|_E \otimes \rho\right)U^\dagger\right]
$$

ここで：
- $|0\rangle_E$ はアンシラの初期状態（標準基底の第0ベクトル）
- $\mathrm{Tr}_E$ はアンシラに対する部分トレース
- $U$ は $(d_{\mathrm{anc}} \cdot d_{\mathrm{sys}}) \times (d_{\mathrm{anc}} \cdot d_{\mathrm{sys}})$ のユニタリ行列

### 1.3 Kraus表現との関係

Stinespring拡張から自然にKraus分解が導かれる。$U$ をアンシラの基底に関してブロック分割すると：

$$
U = \begin{pmatrix} U_{00} & U_{01} & \cdots & U_{0,d_{\mathrm{anc}}-1} \\ U_{10} & U_{11} & \cdots & U_{1,d_{\mathrm{anc}}-1} \\ \vdots & \vdots & \ddots & \vdots \\ U_{d_{\mathrm{anc}}-1,0} & U_{d_{\mathrm{anc}}-1,1} & \cdots & U_{d_{\mathrm{anc}}-1,d_{\mathrm{anc}}-1} \end{pmatrix}
$$

ここで各 $U_{kl}$ は $d_{\mathrm{sys}} \times d_{\mathrm{sys}}$ のブロックであり、

$$
U_{kl} = \langle k|_E\, U\, |l\rangle_E
$$

と定義される。Kraus演算子は：

$$
K_k = U_{k0} = \langle k|_E\, U\, |0\rangle_E, \qquad k = 0, 1, \ldots, d_{\mathrm{anc}} - 1
$$

であり、チャネルは：

$$
\mathcal{E}(\rho) = \sum_{k=0}^{d_{\mathrm{anc}}-1} K_k \rho K_k^\dagger
$$

$U$ のユニタリ性 $U^\dagger U = I$ から、$l = l' = 0$ の列について：

$$
\sum_{k=0}^{d_{\mathrm{anc}}-1} U_{k0}^\dagger U_{k0} = (U^\dagger U)_{00} = I_{d_{\mathrm{sys}}}
$$

すなわち：

$$
\sum_{k=0}^{d_{\mathrm{anc}}-1} K_k^\dagger K_k = I_{d_{\mathrm{sys}}}
$$

これはトレース保存条件である。

---

## 2. Stinespring ユニタリの構成法

### 2.1 生成子の構成

Lindblad演算子 $\tilde{L}_\alpha$（$d_{\mathrm{sys}} \times d_{\mathrm{sys}}$、既に $\sqrt{\gamma_\alpha}$ を含む）から、拡大空間上の生成子 $G_\alpha$ を以下のように構成する。

#### $d_{\mathrm{anc}} = 2$（qubit アンシラ）の場合

$$
G_\alpha = \begin{pmatrix} 0 & \tilde{L}_\alpha^\dagger \\ \tilde{L}_\alpha & 0 \end{pmatrix}
$$

これは $(2 d_{\mathrm{sys}}) \times (2 d_{\mathrm{sys}})$ の行列であり、各ブロックは $d_{\mathrm{sys}} \times d_{\mathrm{sys}}$ である。

$G_\alpha$ はエルミートである：

$$
G_\alpha^\dagger = \begin{pmatrix} 0 & \tilde{L}_\alpha^\dagger \\ \tilde{L}_\alpha & 0 \end{pmatrix} = G_\alpha
$$

#### $d_{\mathrm{anc}} = 3$（qutrit アンシラ、シナリオ5で使用）の場合

$$
G_\alpha = \begin{pmatrix} 0 & \tilde{L}_\alpha^\dagger & 0 \\ \tilde{L}_\alpha & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

これは $(3 d_{\mathrm{sys}}) \times (3 d_{\mathrm{sys}})$ の行列であり、各ブロックは $d_{\mathrm{sys}} \times d_{\mathrm{sys}}$ である。上左 $2 \times 2$ ブロック部分は $d_{\mathrm{anc}} = 2$ の場合と同一であり、第3ブロック行と第3ブロック列は全てゼロである。

$G_\alpha$ はエルミートである。

**コード対応** (`stinespring_utils.py`):

```python
d_sys = L.shape[0]          # 81
dim_total = d_anc * d_sys   # 3 * 81 = 243
G = np.zeros((dim_total, dim_total), dtype=np.complex128)
G[:d_sys, d_sys:2*d_sys] = L.conj().T      # (0,1)ブロック = L†
G[d_sys:2*d_sys, :d_sys] = L               # (1,0)ブロック = L
```

### 2.2 シナリオ5における $d_{\mathrm{anc}} = 3$ の理由

シナリオ5は qudit 量子コンピュータ（各レジスタが $d = 3$ のqutrit）を対象としている。そのため、Stinespring拡張で導入するアンシラも $d_{\mathrm{anc}} = 3$ のqutritとする。$d_{\mathrm{anc}} = 3$ では第3ブロック列・行が全てゼロであるため、$U$ の $(2,0)$ ブロック $K_2$ は $O(\theta^3)$ 以上の高次項にのみ寄与し、主要な物理（$K_0$ と $K_1$）は $d_{\mathrm{anc}} = 2$ の場合と同一である。

**コード対応**: `QuditGKSLSimulator.__init__` において `self.d_anc = params.d` すなわち $d_{\mathrm{anc}} = 3$。

### 2.3 ユニタリ行列の計算

生成子 $G_\alpha$ とタイムステップ $\Delta t$ から、Stinespring ユニタリを以下のように計算する：

$$
U_\alpha(\Delta t) = \exp\!\left(-i\sqrt{\Delta t}\; G_\alpha\right)
$$

ここで $\theta = \sqrt{\Delta t}$ がパラメータである。$G_\alpha$ がエルミートであるため、$U_\alpha$ はユニタリである：$U_\alpha^\dagger U_\alpha = I$。

**コード対応**:

```python
theta = np.sqrt(dt)
U = expm(-1j * theta * G)
```

`expm` は `scipy.linalg.expm` であり、行列指数関数を計算する。結果は $(d_{\mathrm{anc}} \cdot d_{\mathrm{sys}}) \times (d_{\mathrm{anc}} \cdot d_{\mathrm{sys}})$ のユニタリ行列であり、ユニタリ性がアサーション `||U†U - I||_F < 1e-10` で検証される。

#### シナリオ5における具体的な行列サイズ

- $d_{\mathrm{sys}} = 81$（$3^4$）
- $d_{\mathrm{anc}} = 3$
- $G_\alpha$: $243 \times 243$
- $U_\alpha$: $243 \times 243$

26個のLindblad演算子に対して26個のStinespring ユニタリが構成される。

---

## 3. Stinespring ユニタリの密度行列への適用

### 3.1 適用アルゴリズム

密度行列 $\rho$（$d_{\mathrm{sys}} \times d_{\mathrm{sys}}$）に Stinespring ユニタリ $U_\alpha$（$d_{\mathrm{anc}} d_{\mathrm{sys}} \times d_{\mathrm{anc}} d_{\mathrm{sys}}$）を適用する手順は以下の通り。

**ステップ1**: アンシラの初期状態 $|0\rangle\langle 0|_E$ を構成する。

$$
|0\rangle\langle 0|_E = \begin{pmatrix}1 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix} \in \mathbb{C}^{3 \times 3}
$$

**ステップ2**: 拡大密度行列を構成する。

$$
\rho_{\mathrm{ext}} = |0\rangle\langle 0|_E \otimes \rho
$$

$\rho_{\mathrm{ext}}$ は $(d_{\mathrm{anc}} d_{\mathrm{sys}}) \times (d_{\mathrm{anc}} d_{\mathrm{sys}}) = 243 \times 243$ の行列である。Kronecker積のブロック構造は：

$$
\rho_{\mathrm{ext}} = \begin{pmatrix} \rho & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

各ブロックは $81 \times 81$。

**ステップ3**: ユニタリ変換を適用する。

$$
\rho' = U_\alpha\, \rho_{\mathrm{ext}}\, U_\alpha^\dagger
$$

$\rho'$ は $243 \times 243$ の行列。

**ステップ4**: アンシラに対する部分トレースを取る。

$$
\rho_{\mathrm{out}} = \mathrm{Tr}_E[\rho'] = \sum_{k=0}^{d_{\mathrm{anc}}-1} \langle k|_E\, \rho'\, |k\rangle_E = \sum_{k=0}^{2} \rho'_{kk}^{(\mathrm{block})}
$$

ここで $\rho'_{kk}^{(\mathrm{block})}$ は $\rho'$ の $(k, k)$ 番目の $81 \times 81$ ブロックであり、行 $[k \cdot 81, (k+1) \cdot 81)$、列 $[k \cdot 81, (k+1) \cdot 81)$ のスライスに対応する。

**コード対応** (`stinespring_utils.py`):

```python
def apply_stinespring_to_density_matrix(rho, U, d_anc=2):
    d_sys = rho.shape[0]                                # 81
    env0 = np.zeros((d_anc, d_anc), dtype=np.complex128)
    env0[0, 0] = 1.0                                    # |0><0|
    rho_ext = np.kron(env0, rho)                         # 243×243
    rho_prime = U @ rho_ext @ U.conj().T                 # 243×243
    rho_out = np.zeros((d_sys, d_sys), dtype=np.complex128)
    for k in range(d_anc):
        rho_out += rho_prime[k*d_sys:(k+1)*d_sys, k*d_sys:(k+1)*d_sys]
    return rho_out
```

### 3.2 Kraus演算子による等価な表現

上記の手順は以下のKraus表現と等価である：

$$
\mathcal{E}_\alpha(\rho) = \sum_{k=0}^{d_{\mathrm{anc}}-1} K_k^{(\alpha)} \rho\, K_k^{(\alpha)\dagger}
$$

ここで $K_k^{(\alpha)} = U_\alpha[k \cdot d_{\mathrm{sys}} : (k+1) \cdot d_{\mathrm{sys}},\; 0 : d_{\mathrm{sys}}]$ は $U_\alpha$ の第 $k$ ブロック行、第0ブロック列のブロックである：

$$
K_k^{(\alpha)} = (U_\alpha)_{k,0}^{(\mathrm{block})}
$$

**等価性の証明**:

$\rho_{\mathrm{ext}}$ のブロック構造から $(\rho_{\mathrm{ext}})_{l,l'} = \delta_{l,0}\delta_{l',0}\rho$ なので：

$$
\rho'_{k,k'} = \sum_{l,l'} (U_\alpha)_{k,l}^{(\mathrm{block})} (\rho_{\mathrm{ext}})_{l,l'} (U_\alpha^\dagger)_{l',k'}^{(\mathrm{block})} = (U_\alpha)_{k,0}^{(\mathrm{block})}\, \rho\, (U_\alpha)_{k',0}^{(\mathrm{block})\dagger}
$$

部分トレースは $k = k'$ の対角ブロックの和であるから：

$$
\rho_{\mathrm{out}} = \sum_{k} \rho'_{k,k} = \sum_{k} K_k^{(\alpha)} \rho\, K_k^{(\alpha)\dagger}
$$

---

## 4. 1次近似の証明

### 4.1 目標

Stinespring拡張によるチャネル $\mathcal{E}_\alpha(\Delta t)$ が、GKSL散逸子 $\mathcal{D}[\tilde{L}_\alpha]$ の1ステップを以下の精度で近似することを示す：

$$
\mathcal{E}_\alpha(\Delta t)(\rho) = \rho + \Delta t\, \mathcal{D}[\tilde{L}_\alpha](\rho) + O(\Delta t^2)
$$

### 4.2 $G_\alpha^2$ の計算

$d_{\mathrm{anc}} = 3$ の場合：

$$
G_\alpha^2 = \begin{pmatrix} 0 & \tilde{L}_\alpha^\dagger & 0 \\ \tilde{L}_\alpha & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}^2 = \begin{pmatrix} \tilde{L}_\alpha^\dagger \tilde{L}_\alpha & 0 & 0 \\ 0 & \tilde{L}_\alpha \tilde{L}_\alpha^\dagger & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

**検証**（ブロック$(0,0)$）: $(G^2)_{00} = G_{01}G_{10} = \tilde{L}_\alpha^\dagger \cdot \tilde{L}_\alpha = \tilde{L}_\alpha^\dagger \tilde{L}_\alpha$ ✓

**検証**（ブロック$(1,1)$）: $(G^2)_{11} = G_{10}G_{01} = \tilde{L}_\alpha \cdot \tilde{L}_\alpha^\dagger = \tilde{L}_\alpha \tilde{L}_\alpha^\dagger$ ✓

### 4.3 $U_\alpha$ の $\sqrt{\Delta t}$ による展開

$\theta = \sqrt{\Delta t}$ として $U_\alpha = e^{-i\theta G_\alpha}$ を冪級数展開する：

$$
U_\alpha = I - i\theta G_\alpha + \frac{(-i\theta)^2}{2!} G_\alpha^2 + \frac{(-i\theta)^3}{3!} G_\alpha^3 + \cdots
$$

$$
= I - i\theta G_\alpha - \frac{\theta^2}{2} G_\alpha^2 + O(\theta^3)
$$

$\theta^2 = \Delta t$ であるから、$O(\theta^3) = O(\Delta t^{3/2})$ に注意する。

ブロック表示すると：

$$
U_\alpha = \begin{pmatrix} I - \frac{\Delta t}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha & -i\sqrt{\Delta t}\,\tilde{L}_\alpha^\dagger & 0 \\ -i\sqrt{\Delta t}\,\tilde{L}_\alpha & I - \frac{\Delta t}{2}\tilde{L}_\alpha \tilde{L}_\alpha^\dagger & 0 \\ 0 & 0 & I \end{pmatrix} + O(\Delta t^{3/2})
$$

### 4.4 Kraus演算子の導出

$K_k = (U_\alpha)_{k,0}^{(\mathrm{block})}$ であるから：

$$
K_0 = I - \frac{\Delta t}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha + O(\Delta t^{3/2})
$$

$$
K_1 = -i\sqrt{\Delta t}\,\tilde{L}_\alpha + O(\Delta t^{3/2})
$$

$$
K_2 = O(\Delta t^{3/2})
$$

$K_2$ は $d_{\mathrm{anc}} = 3$ の場合に追加される Kraus 演算子であるが、$G_\alpha$ の第3ブロック行が全てゼロであるため、主要な寄与は高次項にのみ現れる。

### 4.5 チャネルの計算

$$
\mathcal{E}_\alpha(\rho) = K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger + K_2 \rho K_2^\dagger
$$

**第1項**:

$$
K_0 \rho K_0^\dagger = \left(I - \frac{\Delta t}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha\right) \rho \left(I - \frac{\Delta t}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha\right) + O(\Delta t^2)
$$

$$
= \rho - \frac{\Delta t}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha\, \rho - \frac{\Delta t}{2}\rho\, \tilde{L}_\alpha^\dagger \tilde{L}_\alpha + O(\Delta t^2)
$$

**第2項**:

$$
K_1 \rho K_1^\dagger = \Delta t\, \tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger + O(\Delta t^2)
$$

**第3項**:

$$
K_2 \rho K_2^\dagger = O(\Delta t^3)
$$

**和**:

$$
\mathcal{E}_\alpha(\rho) = \rho + \Delta t\left(\tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger - \frac{1}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha\, \rho - \frac{1}{2}\rho\, \tilde{L}_\alpha^\dagger \tilde{L}_\alpha\right) + O(\Delta t^2)
$$

$$
= \rho + \Delta t\, \mathcal{D}[\tilde{L}_\alpha](\rho) + O(\Delta t^2)
$$

これはGKSL散逸子の1ステップ・1次近似である。$\square$

### 4.6 近似の意味

上記の結果は：

$$
\mathcal{E}_\alpha(\Delta t) = e^{\mathcal{L}_{D_\alpha} \Delta t} + O(\Delta t^2)
$$

すなわち、Stinespring拡張によるチャネルは真のLindblad時間発展 $e^{\mathcal{L}_{D_\alpha} \Delta t}$ を**1次精度**で近似する。ステップあたりの誤差は $O(\Delta t^2)$ であり、$n_{\mathrm{steps}} = t_{\mathrm{max}}/\Delta t$ ステップの積算で全体の誤差は $O(\Delta t)$ となる。

これがシナリオ5の有効収束次数が $O(\Delta t)$（1次）である主要な原因の1つである（Part 3で詳述）。

### 4.7 トレース保存の厳密な保証

$\mathcal{E}_\alpha$ のKraus表現 $\sum_k K_k \rho K_k^\dagger$ において、$\sum_k K_k^\dagger K_k = I$ が $U_\alpha$ のユニタリ性から厳密に成り立つ（第1.3節）。したがって、Stinespring近似は1次近似であるにもかかわらず、**トレース保存は各ステップで厳密に保証される**。これは行列指数関数の1次Taylor近似（$I + \Delta t\, \mathcal{L}$）ではトレース保存が厳密に成り立たないのと対照的である。

同様に、**完全正値性も各ステップで厳密に保証される**（Kraus形式は構造的にCP）。

---

## 5. $G_\alpha^2$ の具体的なブロック構造

### 5.1 一般的な $G_\alpha^p$ の構造

$d_{\mathrm{anc}} = 3$ の場合、$G_\alpha$ の第3ブロック行と第3ブロック列が全てゼロであるため、$G_\alpha^p$ ($p \geq 1$) の第3ブロック行と第3ブロック列も全てゼロである。したがって、$U_\alpha$ のブロック構造は：

$$
U_\alpha = \begin{pmatrix} (U_\alpha)_{00} & (U_\alpha)_{01} & 0 \\ (U_\alpha)_{10} & (U_\alpha)_{11} & 0 \\ 0 & 0 & I_{d_{\mathrm{sys}}} \end{pmatrix}
$$

ここで上左 $2 \times 2$ ブロック部分は $d_{\mathrm{anc}} = 2$ の場合のStinespring ユニタリ $U_\alpha^{(2)}$ と一致する。

**証明**: $G_\alpha$ を $G_\alpha = \begin{pmatrix} \hat{G} & 0 \\ 0 & 0 \end{pmatrix}$ とブロック分割すると（$\hat{G}$ は $2d_{\mathrm{sys}} \times 2d_{\mathrm{sys}}$）：

$$
e^{-i\theta G_\alpha} = \begin{pmatrix} e^{-i\theta \hat{G}} & 0 \\ 0 & I_{d_{\mathrm{sys}}} \end{pmatrix}
$$

これは行列指数関数のブロック対角性から直ちに従う。したがって $K_0$, $K_1$ は $d_{\mathrm{anc}} = 2$ の場合と完全に一致し、$K_2 = 0$（厳密にゼロ、近似ではない）。

### 5.2 結論

$d_{\mathrm{anc}} = 3$ の Stinespring 拡張は、$d_{\mathrm{anc}} = 2$ の場合と**物理的に完全に同一のチャネル**を実装する（$K_2 = 0$ なので）。$d_{\mathrm{anc}} = 3$ を使用する理由は純粋にハードウェア的な制約（qudit量子コンピュータのネイティブなレジスタ次元に合わせるため）である。

---

## 6. 各Lindblad演算子に対するStinespring ユニタリの具体的な行列

### 6.1 構成の手順

26個のLindblad演算子 $\tilde{L}_\alpha$（各 $81 \times 81$）に対して、以下の手順で Stinespring ユニタリを構成する。

1. $G_\alpha \in \mathbb{C}^{243 \times 243}$ を構成：
   - $G_\alpha[0:81, 81:162] = \tilde{L}_\alpha^\dagger$
   - $G_\alpha[81:162, 0:81] = \tilde{L}_\alpha$
   - その他は全てゼロ

2. $\theta = \sqrt{\Delta t/2}$（半ステップ分。コードではdt/2が使われる）

3. $U_\alpha = \mathrm{expm}(-i\theta\, G_\alpha) \in \mathbb{C}^{243 \times 243}$

4. ユニタリ性検証：$\|U_\alpha^\dagger U_\alpha - I_{243}\|_F < 10^{-10}$

**コード対応** (`QuditGKSLSimulator._precompute_unitaries`):

```python
self._U_stines_half = [
    stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
    for L_op, _gamma in self.lindblad_ops
]
```

ここで `dt/2` が渡されるのは、Trotterステップにおいて各Lindbladチャネルが「半ステップ」として適用されるためである（Part 3で詳述）。

### 6.2 計算量

各 $U_\alpha$ の計算には $243 \times 243$ の行列指数関数が必要であり、`scipy.linalg.expm` の計算量は $O(243^3) = O(d_{\mathrm{anc}}^3 d_{\mathrm{sys}}^3)$ である。26個の Lindblad 演算子に対してこれを事前計算する。事前計算は1回のみであり、各Trotterステップでは事前計算された $U_\alpha$ を再利用する。

---

## 7. 部分トレースのブロック抽出法

### 7.1 $\mathrm{kron}(E, S)$ 順序における部分トレース

本コードでは $\rho_{\mathrm{ext}} = |0\rangle\langle 0|_E \otimes \rho$ のKronecker積の順序が `kron(env, system)` であり、拡大密度行列はアンシラが「外側」、系が「内側」のブロック構造を持つ。

$d_{\mathrm{anc}} = 3$, $d_{\mathrm{sys}} = 81$ の場合、$\rho' = U \rho_{\mathrm{ext}} U^\dagger$ は $243 \times 243$ の行列であり、$3 \times 3$ のブロック構造（各ブロック $81 \times 81$）を持つ。

アンシラに対する部分トレースは、対角ブロックの和である：

$$
\mathrm{Tr}_E[\rho'] = \rho'[0:81, 0:81] + \rho'[81:162, 81:162] + \rho'[162:243, 162:243]
$$

### 7.2 一般の順序との関係

もし `kron(system, env)` の順序を使った場合、ブロック構造が異なり、部分トレースの取り方も変わる。本コードでは一貫して `kron(env, system)` を使用しており、Stinespring ユニタリの生成子 $G_\alpha$ のブロック構造もこの順序に合わせて定義されている。

---

## 8. Stinespring近似の誤差の定量的評価

### 8.1 1ステップあたりの誤差

Stinespring近似の1ステップあたりの誤差は：

$$
\left\|\mathcal{E}_\alpha(\Delta t) - e^{\mathcal{L}_{D_\alpha} \Delta t}\right\| = O(\Delta t^2)
$$

ここでノルムはダイアモンドノルム（完全有界ノルム）を想定している。

具体的には、$e^{\mathcal{L}_{D_\alpha} \Delta t}(\rho) = \rho + \Delta t\,\mathcal{D}[\tilde{L}_\alpha](\rho) + \frac{\Delta t^2}{2}\mathcal{D}[\tilde{L}_\alpha]^2(\rho) + \cdots$ であり、Stinespring近似は1次項まで一致するが、2次項以降は一般に異なる。

### 8.2 トレース距離での評価

2つの密度行列 $\rho$ と $\sigma$ のトレース距離は：

$$
D(\rho, \sigma) = \frac{1}{2}\|\rho - \sigma\|_1 = \frac{1}{2}\mathrm{Tr}\!\left[\sqrt{(\rho - \sigma)^\dagger(\rho - \sigma)}\right]
$$

Stinespring近似から生じるトレース距離誤差は、$n_{\mathrm{steps}}$ ステップ後に蓄積する。個別のLindblad チャネルのStinespring近似誤差は1ステップあたり $O(\Delta t^2)$ であるが、26個のチャネルの積とTrotter分割誤差が加わり、全体の誤差は Part 3 で詳述する。

---

## 付録A. エルミート行列の行列指数関数のユニタリ性

エルミート行列 $H$ に対して $U = e^{-iH}$ はユニタリである：

$$
U^\dagger U = e^{iH^\dagger} e^{-iH} = e^{iH} e^{-iH} = I
$$

これは $H = H^\dagger$ より $H^\dagger = H$ であるため、$(e^{-iH})^\dagger = e^{iH^\dagger} = e^{iH}$ から従う。

$G_\alpha$ はエルミートであるため、$U_\alpha = e^{-i\sqrt{\Delta t}\,G_\alpha}$ のユニタリ性が保証される。

## 付録B. $G_\alpha$ のスペクトル構造

$d_{\mathrm{anc}} = 3$ の場合、$G_\alpha$ のブロック構造から：

$$
G_\alpha = \begin{pmatrix} \hat{G}_\alpha & 0 \\ 0 & 0_{d_{\mathrm{sys}}} \end{pmatrix}, \qquad
\hat{G}_\alpha = \begin{pmatrix} 0 & \tilde{L}_\alpha^\dagger \\ \tilde{L}_\alpha & 0 \end{pmatrix}
$$

$G_\alpha$ の固有値は $\hat{G}_\alpha$ の固有値と $d_{\mathrm{sys}}$ 個のゼロ固有値からなる。

$\hat{G}_\alpha$ の固有値は $\tilde{L}_\alpha$ の特異値 $\sigma_k$ と関係する。$\tilde{L}_\alpha = W \Sigma V^\dagger$（特異値分解）とすると、$\hat{G}_\alpha$ の固有値は $\pm \sigma_k$ である。したがって $G_\alpha$ の固有値は $\{+\sigma_1, -\sigma_1, +\sigma_2, -\sigma_2, \ldots, 0, 0, \ldots\}$ である。

本シナリオの Lindblad 演算子は全てランク1（単一の遷移演算子に比例）またはランク低い演算子であるため、$G_\alpha$ の非ゼロ固有値は少数（多くの場合1対の $\pm\sigma$）である。

---

## 参考文献

- W. F. Stinespring, "Positive functions on C*-algebras," Proc. Amer. Math. Soc. **6**, 211 (1955).
- M.-D. Choi, "Completely positive linear maps on complex matrices," Linear Algebra Appl. **10**, 285 (1975).
- K. Kraus, *States, Effects, and Operations: Fundamental Notions of Quantum Theory*, Lecture Notes in Physics **190**, Springer (1983).

---

*Part 3 では、Trotter-Suzuki 分解（Strang 分割）、回文積（palindromic product）順序、完全なシミュレーションアルゴリズム、誤差解析と収束次数、ゲート数見積を扱う。*
