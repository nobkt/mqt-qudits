# TTA-UC現象の古典GKSL-Lindbladダイナミクス：シナリオ1の省略無し詳細理論

本文書は、`tutorials/quantum_dynamics_gksl_comparison.ipynb` の「2. シナリオ1: 古典GKSL（ボソン無し）」で実行される計算の背後にある理論を、プログラムで計算可能な解像度（行列の具体的な構成法を含む）まで省略無しに定式化する。

---

## 1. 物理系の定義

### 1.1 単一分子の状態空間

各分子は3準位系（qutrit, $d=3$）としてモデル化される。基底は以下の通りである：

$$
|0\rangle = \mathrm{S}_0 \quad (\text{基底一重項状態}), \qquad
|1\rangle = \mathrm{T}_1 \quad (\text{第一励起三重項状態}), \qquad
|2\rangle = \mathrm{S}_1 \quad (\text{第一励起一重項状態})
$$

単一分子のヒルベルト空間は $\mathcal{H}_\mathrm{mol} = \mathbb{C}^3$ であり、標準基底ベクトルの列ベクトル表現は：

$$
|0\rangle = \begin{pmatrix}1\\0\\0\end{pmatrix}, \quad
|1\rangle = \begin{pmatrix}0\\1\\0\end{pmatrix}, \quad
|2\rangle = \begin{pmatrix}0\\0\\1\end{pmatrix}
$$

### 1.2 多分子系のヒルベルト空間

$N=4$ 個の分子からなる系の全ヒルベルト空間は、テンソル積空間である：

$$
\mathcal{H} = \mathcal{H}_0 \otimes \mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \mathcal{H}_3 = (\mathbb{C}^3)^{\otimes 4}
$$

次元は $\dim \mathcal{H} = 3^4 = 81$ である。

### 1.3 計算基底のラベリング

81個の計算基底状態 $|s_0, s_1, s_2, s_3\rangle$（$s_i \in \{0,1,2\}$）は、混合基数（radix-3）分解により単一の整数インデックス $n \in \{0,1,\ldots,80\}$ に対応付けられる。本コードではKronecker積の標準的な順序（最左のサイトが最上位桁）を用いる：

$$
n = s_0 \cdot 3^3 + s_1 \cdot 3^2 + s_2 \cdot 3^1 + s_3 \cdot 3^0 = 27s_0 + 9s_1 + 3s_2 + s_3
$$

逆に、インデックス $n$ から各分子の局所状態を復元するには：

$$
s_3 = n \bmod 3, \quad
s_2 = \lfloor n/3 \rfloor \bmod 3, \quad
s_1 = \lfloor n/9 \rfloor \bmod 3, \quad
s_0 = \lfloor n/27 \rfloor \bmod 3
$$

**注意**: コード内の `compute_populations_from_density_matrix` では `mol=N-1` から `mol=0` へ逆順にループし、`remainder = idx`, `local_state = remainder % d`, `remainder //= d` で分解している。これは上記の最左上位桁の定義と整合する（最初に取り出されるのは $s_3$、最後が $s_0$）。

### 1.4 物理パラメータ

| パラメータ | 記号 | 値 | 単位 | 物理的意味 |
|---|---|---|---|---|
| 三重項エネルギー | $E_\mathrm{T}$ | 1.5 | eV | T₁状態のエネルギー |
| 一重項エネルギー | $E_\mathrm{S}$ | 3.0 | eV | S₁状態のエネルギー |
| 移動結合 | $V$ | 0.1 | eV | 隣接分子間の三重項エネルギー移動結合 |
| TTA速度 | $\gamma_\mathrm{TTA}$ | 0.05 | eV/$\hbar$ | 三重項-三重項消滅速度 |
| 蛍光速度 | $\Gamma_\mathrm{fl}$ | 0.01 | eV/$\hbar$ | S₁→S₀放射遷移速度 |
| 燐光速度 | $\Gamma_\mathrm{ph}$ | $10^{-6}$ | eV/$\hbar$ | T₁→S₀放射遷移速度 |
| 内部転換速度 | $k_\mathrm{IC}$ | 0.005 | eV/$\hbar$ | S₁→S₀内部転換速度 |
| ISC (S→T) 速度 | $k_\mathrm{ISC,ST}$ | 0.003 | eV/$\hbar$ | S₁→T₁項間交差速度 |
| ISC (T→S) 速度 | $k_\mathrm{ISC,TS}$ | $10^{-5}$ | eV/$\hbar$ | T₁→S₀項間交差速度 |

本モデルでは $\hbar = 1$（自然単位系）を採用している。時間の物理単位への変換には $\hbar = 0.6582$ eV·fs を用いるが、シミュレーション内部では無次元時間で計算する。

### 1.5 隣接分子ペア

$N=4$ の鎖状配置における最近接ペアは：

$$
\text{neighbors} = \{(0,1),\;(1,2),\;(2,3)\}
$$

ペア数は $|\text{neighbors}| = N-1 = 3$ である。

---

## 2. ハミルトニアン

全ハミルトニアンは以下の2つの成分の和である：

$$
H = H_0 + H_\mathrm{transfer}
$$

### 2.1 オンサイト・ハミルトニアン $H_0$

各分子のエネルギー準位を表す。単一分子のオンサイト・ハミルトニアンは：

$$
h_\mathrm{local} = E_\mathrm{T}|1\rangle\langle 1| + E_\mathrm{S}|2\rangle\langle 2| = \begin{pmatrix}0 & 0 & 0\\0 & E_\mathrm{T} & 0\\0 & 0 & E_\mathrm{S}\end{pmatrix} = \begin{pmatrix}0 & 0 & 0\\0 & 1.5 & 0\\0 & 0 & 3.0\end{pmatrix}
$$

全系のオンサイト・ハミルトニアンは、各分子サイトへの埋め込みの総和である：

$$
H_0 = \sum_{i=0}^{N-1} h_\mathrm{local}^{(i)}
$$

ここで $h_\mathrm{local}^{(i)}$ はサイト $i$ にのみ作用し、他のサイトには恒等演算子が置かれる：

$$
h_\mathrm{local}^{(i)} = \underbrace{I_3 \otimes \cdots \otimes I_3}_{i \text{ 個}} \otimes\; h_\mathrm{local} \;\otimes \underbrace{I_3 \otimes \cdots \otimes I_3}_{N-1-i \text{ 個}}
$$

$I_3$ は $3\times 3$ 単位行列である。

**具体的な構成法（コード対応）**: `build_single_site_operator(op, site, N, d)` は、`N` 個の $d \times d$ 単位行列のリスト `[I_d, I_d, ..., I_d]` の `site` 番目を `op` に置換し、それらの順次 Kronecker 積を取る。すなわち：

```
op_list = [I_d] * N
op_list[site] = op
result = op_list[0] ⊗ op_list[1] ⊗ ... ⊗ op_list[N-1]
```

$H_0$ は $81 \times 81$ のエルミート対角行列である。対角要素 $(n, n)$ は、インデックス $n$ に対応する状態 $|s_0, s_1, s_2, s_3\rangle$ のエネルギーの総和：

$$
(H_0)_{nn} = \sum_{i=0}^{3} \varepsilon(s_i), \qquad \varepsilon(s) = \begin{cases}0 & s=0\\E_\mathrm{T} & s=1\\E_\mathrm{S} & s=2\end{cases}
$$

非対角要素は全てゼロである。

### 2.2 三重項エネルギー移動ハミルトニアン $H_\mathrm{transfer}$

隣接分子間の三重項励起子の移動（ホッピング）を記述する：

$$
H_\mathrm{transfer} = \sum_{\langle i,j \rangle} V \left( |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1| \right)
$$

ここで和は最近接ペア $\langle i,j \rangle \in \{(0,1),(1,2),(2,3)\}$ に渡る。

各項の物理的意味は：
- $|0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0|$：分子 $i$ の三重項が消滅し、分子 $j$ に三重項が生成される（$i \to j$ ホッピング）
- エルミート共役が逆方向ホッピング（$j \to i$）を与える

**単一サイト演算子の行列表現**:

$$
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}, \qquad
|1\rangle\langle 0| = \begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}
$$

**ペア $(i,j)$ に対する全系演算子の構成法**: サイト $i$ に $|0\rangle\langle 1|$、サイト $j$ に $|1\rangle\langle 0|$ を置き、その他のサイトに $I_3$ を置いた Kronecker 積を計算する：

$$
F_{ij} = \underbrace{I_3 \otimes \cdots}_{0 \text{ to } i-1} \otimes\; |0\rangle\langle 1| \;\otimes \underbrace{\cdots \otimes I_3 \otimes \cdots}_{i+1 \text{ to } j-1} \otimes\; |1\rangle\langle 0| \;\otimes \underbrace{I_3 \otimes \cdots}_{j+1 \text{ to } N-1}
$$

各ペアの寄与は $V(F_{ij} + F_{ij}^\dagger)$ である。$F_{ij}^\dagger = |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1|$ なので、合わせてエルミート演算子となる。

$H_\mathrm{transfer}$ は $81 \times 81$ のエルミート行列であり、$H_0$ と異なり非対角成分を持つ。

### 2.3 全ハミルトニアン $H$

$$
H = H_0 + H_\mathrm{transfer}
$$

$H$ は $81 \times 81$ のエルミート行列である。

---

## 3. Lindblad演算子

GKSL方程式における散逸過程は、Lindblad演算子 $L_\alpha$ で記述される。コードでは各 $L_\alpha$ に $\sqrt{\gamma_\alpha}$ が含まれた形で構成される。すなわち、返される演算子を $\tilde{L}_\alpha$ とすると：

$$
\tilde{L}_\alpha = \sqrt{\gamma_\alpha}\; L_\alpha^{(\mathrm{bare})}
$$

以下、$\tilde{L}_\alpha$ をそのまま「Lindblad演算子」と呼ぶ。全部で $2 \times |\mathrm{neighbors}| + 5 \times N = 2 \times 3 + 5 \times 4 = 26$ 個の演算子がある。

### 3.1 三重項-三重項消滅 (TTA) — 6個

TTA過程 $\mathrm{T}_1 + \mathrm{T}_1 \to \mathrm{S}_1 + \mathrm{S}_0$ を記述する。各隣接ペア $(i,j)$ に対して2つのチャネルがある。速度定数は $\gamma = \gamma_\mathrm{TTA}/2 = 0.025$（対称的に2チャネルで分配）。

**チャネル1**: 分子 $i$ が $\mathrm{T}_1 \to \mathrm{S}_1$ に励起され、分子 $j$ が $\mathrm{T}_1 \to \mathrm{S}_0$ に脱励起される。

$$
\tilde{L}_\mathrm{TTA,1}^{(i,j)} = \sqrt{\gamma}\; \left(|2\rangle_i\langle 1|\right) \otimes \left(|0\rangle_j\langle 1|\right) \otimes \bigotimes_{k \neq i,j} I_3^{(k)}
$$

**チャネル2**: 分子 $i$ が $\mathrm{T}_1 \to \mathrm{S}_0$ に脱励起され、分子 $j$ が $\mathrm{T}_1 \to \mathrm{S}_1$ に励起される。

$$
\tilde{L}_\mathrm{TTA,2}^{(i,j)} = \sqrt{\gamma}\; \left(|0\rangle_i\langle 1|\right) \otimes \left(|2\rangle_j\langle 1|\right) \otimes \bigotimes_{k \neq i,j} I_3^{(k)}
$$

**単一サイト演算子の行列表現**:

$$
|2\rangle\langle 1| = \begin{pmatrix}0&0&0\\0&0&0\\0&1&0\end{pmatrix}, \qquad
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

3ペア $(0,1),(1,2),(2,3)$ に対してそれぞれ2チャネルで、合計6個の $81 \times 81$ 行列。

### 3.2 蛍光 (Fluorescence) — 4個

$\mathrm{S}_1 \to \mathrm{S}_0$ の放射遷移（蛍光発光）。各分子 $i \in \{0,1,2,3\}$ に対して：

$$
\tilde{L}_\mathrm{fl}^{(i)} = \sqrt{\Gamma_\mathrm{fl}}\; |0\rangle_i\langle 2| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：

$$
|0\rangle\langle 2| = \begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}
$$

$\gamma_\alpha = \Gamma_\mathrm{fl} = 0.01$、$\sqrt{\gamma_\alpha} = 0.1$。

### 3.3 燐光 (Phosphorescence) — 4個

$\mathrm{T}_1 \to \mathrm{S}_0$ の放射遷移（燐光発光）。各分子 $i$ に対して：

$$
\tilde{L}_\mathrm{ph}^{(i)} = \sqrt{\Gamma_\mathrm{ph}}\; |0\rangle_i\langle 1| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：

$$
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

$\gamma_\alpha = \Gamma_\mathrm{ph} = 10^{-6}$、$\sqrt{\gamma_\alpha} = 10^{-3}$。

### 3.4 内部転換 (Internal Conversion, IC) — 4個

$\mathrm{S}_1 \to \mathrm{S}_0$ の非放射遷移。各分子 $i$ に対して：

$$
\tilde{L}_\mathrm{IC}^{(i)} = \sqrt{k_\mathrm{IC}}\; |0\rangle_i\langle 2| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：$|0\rangle\langle 2|$（蛍光と同じ遷移演算子、異なる速度定数）。

$\gamma_\alpha = k_\mathrm{IC} = 0.005$、$\sqrt{\gamma_\alpha} \approx 0.07071$。

### 3.5 項間交差 S→T (ISC S→T) — 4個

$\mathrm{S}_1 \to \mathrm{T}_1$ のスピン禁制遷移。各分子 $i$ に対して：

$$
\tilde{L}_\mathrm{ISC,ST}^{(i)} = \sqrt{k_\mathrm{ISC,ST}}\; |1\rangle_i\langle 2| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：

$$
|1\rangle\langle 2| = \begin{pmatrix}0&0&0\\0&0&1\\0&0&0\end{pmatrix}
$$

$\gamma_\alpha = k_\mathrm{ISC,ST} = 0.003$、$\sqrt{\gamma_\alpha} \approx 0.05477$。

### 3.6 項間交差 T→S (ISC T→S) — 4個

$\mathrm{T}_1 \to \mathrm{S}_0$ のスピン禁制遷移。各分子 $i$ に対して：

$$
\tilde{L}_\mathrm{ISC,TS}^{(i)} = \sqrt{k_\mathrm{ISC,TS}}\; |0\rangle_i\langle 1| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：$|0\rangle\langle 1|$（燐光と同じ遷移演算子、異なる速度定数）。

$\gamma_\alpha = k_\mathrm{ISC,TS} = 10^{-5}$、$\sqrt{\gamma_\alpha} \approx 0.003162$。

### 3.7 Lindblad演算子の一覧表

| 番号 | 種類 | ペア/サイト | $\gamma_\alpha$ | 単一サイト $i$ | 単一サイト $j$ | 個数 |
|---|---|---|---|---|---|---|
| 1–6 | TTA ch.1,2 | $(0,1),(1,2),(2,3)$ | $\gamma_\mathrm{TTA}/2 = 0.025$ | $\|2\rangle\langle 1\|$ or $\|0\rangle\langle 1\|$ | $\|0\rangle\langle 1\|$ or $\|2\rangle\langle 1\|$ | 6 |
| 7–10 | 蛍光 | $i=0,1,2,3$ | $\Gamma_\mathrm{fl} = 0.01$ | $\|0\rangle\langle 2\|$ | — | 4 |
| 11–14 | 燐光 | $i=0,1,2,3$ | $\Gamma_\mathrm{ph} = 10^{-6}$ | $\|0\rangle\langle 1\|$ | — | 4 |
| 15–18 | IC | $i=0,1,2,3$ | $k_\mathrm{IC} = 0.005$ | $\|0\rangle\langle 2\|$ | — | 4 |
| 19–22 | ISC S→T | $i=0,1,2,3$ | $k_\mathrm{ISC,ST} = 0.003$ | $\|1\rangle\langle 2\|$ | — | 4 |
| 23–26 | ISC T→S | $i=0,1,2,3$ | $k_\mathrm{ISC,TS} = 10^{-5}$ | $\|0\rangle\langle 1\|$ | — | 4 |

合計：**26個**の $81 \times 81$ Lindblad演算子。

---

## 4. GKSL-Lindblad マスター方程式

密度行列 $\rho(t)$ の時間発展は GKSL (Gorini-Kossakowski-Sudarshan-Lindblad) マスター方程式に従う：

$$
\frac{d\rho}{dt} = \mathcal{L}[\rho] = -\frac{i}{\hbar}[H, \rho] + \sum_{\alpha=1}^{26} \mathcal{D}[\tilde{L}_\alpha](\rho)
$$

ここで $\hbar = 1$（自然単位系）、各Lindblad散逸子は：

$$
\mathcal{D}[\tilde{L}_\alpha](\rho) = \tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger - \frac{1}{2}\left(\tilde{L}_\alpha^\dagger \tilde{L}_\alpha \rho + \rho \tilde{L}_\alpha^\dagger \tilde{L}_\alpha\right)
$$

$\tilde{L}_\alpha$ には既に $\sqrt{\gamma_\alpha}$ が含まれているため、散逸子を展開すると：

$$
\mathcal{D}[\tilde{L}_\alpha](\rho) = \gamma_\alpha \left( L_\alpha^{(\mathrm{bare})} \rho L_\alpha^{(\mathrm{bare})\dagger} - \frac{1}{2}\left\{L_\alpha^{(\mathrm{bare})\dagger} L_\alpha^{(\mathrm{bare})},\; \rho\right\}\right)
$$

### 4.1 ハミルトニアン交換子の成分表示

$[H, \rho] = H\rho - \rho H$ であるから、密度行列要素について：

$$
\left([H, \rho]\right)_{mn} = \sum_k H_{mk}\rho_{kn} - \sum_k \rho_{mk}H_{kn}
$$

### 4.2 散逸子の成分表示

$\tilde{L}_\alpha^\dagger \tilde{L}_\alpha$ を $M_\alpha$ と略記する。

$$
\left(\mathcal{D}[\tilde{L}_\alpha](\rho)\right)_{mn} = \sum_{k,l} (\tilde{L}_\alpha)_{mk}(\tilde{L}_\alpha^*)_{nl}\, \rho_{kl} - \frac{1}{2}\sum_k (M_\alpha)_{mk}\rho_{kn} - \frac{1}{2}\sum_k \rho_{mk}(M_\alpha)_{kn}
$$

---

## 5. リウヴィリアン超演算子のベクトル化表現

GKSL方程式を線形常微分方程式として解くため、密度行列をベクトル化する。

### 5.1 列主序（Column-major）ベクトル化

密度行列 $\rho$ ($81 \times 81$) を列主序（Fortran順序）で $81^2 = 6561$ 次元のベクトル $\mathrm{vec}(\rho)$ に変換する：

$$
\mathrm{vec}(\rho) = \begin{pmatrix}\rho_{0,0}\\\rho_{1,0}\\\vdots\\\rho_{80,0}\\\rho_{0,1}\\\vdots\\\rho_{80,80}\end{pmatrix}
$$

すなわち、$\mathrm{vec}(\rho)$ の第 $(n + 81m)$ 成分は $\rho_{n,m}$ である（$0$-indexed）。

**コード対応**: `vectorize_density_matrix(rho)` は `rho.flatten(order="F")` を実行する。`unvectorize_density_matrix(vec, dim)` は `vec.reshape((dim, dim), order="F")` を実行する。

### 5.2 ベクトル化の基本恒等式

列主序ベクトル化に対して、以下の恒等式が成り立つ：

$$
\mathrm{vec}(AXB) = (B^\mathrm{T} \otimes A)\,\mathrm{vec}(X)
$$

ここで $\otimes$ は Kronecker 積、$B^\mathrm{T}$ は転置である。

### 5.3 ハミルトニアン部分の超演算子

ハミルトニアン部分 $-i[H, \rho] = -i(H\rho - \rho H)$ をベクトル化する：

$$
\mathrm{vec}(-i H\rho) = -i(I \otimes H)\,\mathrm{vec}(\rho)
$$

$$
\mathrm{vec}(i \rho H) = i(H^\mathrm{T} \otimes I)\,\mathrm{vec}(\rho)
$$

よって：

$$
\mathcal{L}_H = -i(I_{81} \otimes H - H^\mathrm{T} \otimes I_{81})
$$

$H$ はエルミートなので $H^\mathrm{T} = H^*$ であるが、コードでは明示的に $H^\mathrm{T}$ を使用している。$\mathcal{L}_H$ は $6561 \times 6561$ の行列である。

### 5.4 散逸部分の超演算子

各Lindblad演算子 $\tilde{L}_\alpha$ に対する散逸子 $\mathcal{D}[\tilde{L}_\alpha](\rho)$ をベクトル化する。$M_\alpha = \tilde{L}_\alpha^\dagger \tilde{L}_\alpha$ とする。

$$
\mathrm{vec}(\tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger) = (\tilde{L}_\alpha^{\dagger\mathrm{T}} \otimes \tilde{L}_\alpha)\,\mathrm{vec}(\rho) = (\bar{L}_\alpha \otimes \tilde{L}_\alpha)\,\mathrm{vec}(\rho)
$$

ここで $\bar{L}_\alpha = \tilde{L}_\alpha^*$（要素ごとの複素共役）。

$$
\mathrm{vec}\left(-\frac{1}{2}M_\alpha \rho\right) = -\frac{1}{2}(I \otimes M_\alpha)\,\mathrm{vec}(\rho)
$$

$$
\mathrm{vec}\left(-\frac{1}{2}\rho M_\alpha\right) = -\frac{1}{2}(M_\alpha^\mathrm{T} \otimes I)\,\mathrm{vec}(\rho)
$$

よって、1つの散逸子の超演算子表現は：

$$
\mathcal{L}_{\alpha} = \bar{L}_\alpha \otimes \tilde{L}_\alpha - \frac{1}{2}(I \otimes M_\alpha) - \frac{1}{2}(M_\alpha^\mathrm{T} \otimes I)
$$

全散逸部分は：

$$
\mathcal{L}_D = \sum_{\alpha=1}^{26} \mathcal{L}_{\alpha}
$$

### 5.5 全リウヴィリアン超演算子

$$
\mathcal{L} = \mathcal{L}_H + \mathcal{L}_D = -i(I \otimes H - H^\mathrm{T} \otimes I) + \sum_{\alpha=1}^{26}\left[\bar{L}_\alpha \otimes \tilde{L}_\alpha - \frac{1}{2}(I \otimes M_\alpha + M_\alpha^\mathrm{T} \otimes I)\right]
$$

$\mathcal{L}$ は $6561 \times 6561$ の（一般に非エルミートな）複素行列である。

**コード対応**: `build_gksl_superoperator(H_total, lindblad_ops)` で構成される。コード内の `np.kron(L_op.conj(), L_op)` は $\bar{L}_\alpha \otimes \tilde{L}_\alpha$ に対応する。

---

## 6. 時間発展の形式解

### 6.1 線形ODE系

ベクトル化されたGKSL方程式は線形ODE系：

$$
\frac{d}{dt}\mathrm{vec}(\rho(t)) = \mathcal{L}\,\mathrm{vec}(\rho(t))
$$

### 6.2 形式解

形式解は行列指数関数で与えられる：

$$
\mathrm{vec}(\rho(t)) = e^{\mathcal{L} t}\,\mathrm{vec}(\rho(0))
$$

Lindblad-GKS定理により、正当なGKSL生成子 $\mathcal{L}$ に対して、$e^{\mathcal{L} t}$ は完全正・トレース保存（CPTP）写像を与える。したがって、密度行列の正値性とトレース保存は構造的に保証される。

### 6.3 数値計算法

コードでは `scipy.sparse.linalg.expm_multiply` を用いて、$e^{\mathcal{L} t}\mathbf{v}$ を全行列指数関数を陽に構成せずに効率的に計算する：

```python
vecs = expm_multiply(L_super, vec_0, start=0.0, stop=t_max, num=n_steps+1, endpoint=True)
```

これは $n_\mathrm{steps} + 1 = 101$ 個の等間隔時刻 $t_k = k \cdot \Delta t$（$\Delta t = t_\mathrm{max}/n_\mathrm{steps} = 1.0$、$k = 0, 1, \ldots, 100$）に対して：

$$
\mathrm{vec}(\rho(t_k)) = e^{\mathcal{L} t_k}\,\mathrm{vec}(\rho(0))
$$

を計算する。`expm_multiply` は Krylov 部分空間法（Al-Mohy & Higham, 2011）に基づき、疎行列に対して特に効率的である。

---

## 7. 初期状態

### 7.1 エッジ三重項状態

初期状態 `'edge_triplet'` は、鎖の両端（分子0と分子3）が三重項 $\mathrm{T}_1$ 状態にあり、中間の分子（1と2）が基底状態 $\mathrm{S}_0$ にある純粋状態である：

$$
|\psi_0\rangle = |1, 0, 0, 1\rangle
$$

### 7.2 基底インデックス

$$
n = 1 \cdot 27 + 0 \cdot 9 + 0 \cdot 3 + 1 \cdot 1 = 28
$$

**コード対応**: `index = 1 * (d ** (N - 1)) + 1 = 1 * 27 + 1 = 28`

ここで `1 * (d ** (N - 1))` は $s_0 = 1$ に対応する項、末尾の `+ 1` は $s_3 = 1$ に対応する項であり、$s_1 = s_2 = 0$ なので中間項は寄与しない。

### 7.3 初期密度行列

$$
\rho(0) = |\psi_0\rangle\langle\psi_0|
$$

これは $81 \times 81$ の行列で、$(28, 28)$ 成分のみが1、他は全てゼロである：

$$
\rho(0)_{mn} = \delta_{m,28}\,\delta_{n,28}
$$

### 7.4 初期ベクトル

$$
\mathrm{vec}(\rho(0)) \in \mathbb{C}^{6561}
$$

列主序ベクトル化により、$\rho(0)_{28,28} = 1$ に対応するベクトルの非ゼロ成分の位置は：

$$
\text{位置} = 28 + 81 \times 28 = 28 + 2268 = 2296
$$

すなわち、$\mathrm{vec}(\rho(0))$ の第2296成分のみが1で、他は全てゼロである。

---

## 8. 物理量の抽出

各時刻 $t_k$ において、ベクトル $\mathrm{vec}(\rho(t_k))$ を $81 \times 81$ 密度行列 $\rho(t_k)$ に逆変換した後、以下の物理量を計算する。

### 8.1 トレース（規格化の検証）

$$
\mathrm{Tr}[\rho(t_k)] = \sum_{n=0}^{80} \rho_{nn}(t_k)
$$

CPTP性により $\mathrm{Tr}[\rho(t_k)] = 1$ が保証される。数値的には $|\mathrm{Tr}[\rho] - 1| < \epsilon$（$\epsilon$ はマシン精度程度）を検証する。

### 8.2 占有数（ポピュレーション）

各電子状態（$\mathrm{S}_0$, $\mathrm{T}_1$, $\mathrm{S}_1$）の平均占有数を、密度行列の対角要素から計算する：

$$
N_{\mathrm{S}_0}(t) = \sum_{n=0}^{80} p_n(t) \cdot c_{\mathrm{S}_0}(n), \qquad
N_{\mathrm{T}_1}(t) = \sum_{n=0}^{80} p_n(t) \cdot c_{\mathrm{T}_1}(n), \qquad
N_{\mathrm{S}_1}(t) = \sum_{n=0}^{80} p_n(t) \cdot c_{\mathrm{S}_1}(n)
$$

ここで $p_n(t) = \mathrm{Re}(\rho_{nn}(t))$ は基底状態 $|n\rangle$ の占有確率であり、$c_X(n)$ は基底状態 $|n\rangle = |s_0, s_1, s_2, s_3\rangle$ において電子状態 $X$ にある分子の数：

$$
c_{\mathrm{S}_0}(n) = \#\{i : s_i = 0\}, \qquad
c_{\mathrm{T}_1}(n) = \#\{i : s_i = 1\}, \qquad
c_{\mathrm{S}_1}(n) = \#\{i : s_i = 2\}
$$

ここで $s_i$ はインデックス $n$ の混合基数分解における分子 $i$ の局所状態である（第1.3節参照）。

規格化の整合性から：

$$
N_{\mathrm{S}_0}(t) + N_{\mathrm{T}_1}(t) + N_{\mathrm{S}_1}(t) = N \cdot \mathrm{Tr}[\rho(t)] = 4
$$

**コード対応**: `compute_populations_from_density_matrix(rho, params)` では、各対角要素 `diag[idx]` に対して混合基数分解を行い、各分子の局所状態に応じて `N_S0`, `N_T1`, `N_S1` に加算する。

**初期値**: $|\psi_0\rangle = |1,0,0,1\rangle$ なので、$N_{\mathrm{T}_1}(0) = 2$、$N_{\mathrm{S}_0}(0) = 2$、$N_{\mathrm{S}_1}(0) = 0$ である。

### 8.3 フォン・ノイマンエントロピー

$$
S(t) = -\mathrm{Tr}[\rho(t)\ln\rho(t)] = -\sum_{i=0}^{80} \lambda_i(t) \ln\lambda_i(t)
$$

ここで $\{\lambda_i(t)\}$ は $\rho(t)$ の固有値であり、$\lambda_i > 0$ の項のみを含める（$\lambda_i \leq 0$ の項は数値誤差として除外）。

**コード対応**: `compute_von_neumann_entropy(rho, tol=1e-12)` では `np.linalg.eigvalsh(rho)` で固有値を計算し、`tol` より大きい固有値のみで $-\sum \lambda_i \ln \lambda_i$ を計算する。

純粋状態のとき $S = 0$、最大混合状態のとき $S = \ln 81 \approx 4.394$ である。初期状態は純粋状態なので $S(0) = 0$ であり、散逸過程によりエントロピーは増大する。

### 8.4 純度（Purity）

$$
P(t) = \mathrm{Tr}[\rho(t)^2] = \sum_{m,n} |\rho_{mn}(t)|^2
$$

**コード対応**: `compute_purity(rho)` は `np.real(np.trace(rho @ rho))` を計算する。

純粋状態のとき $P = 1$、最大混合状態のとき $P = 1/81 \approx 0.01235$ である。初期状態は純粋状態なので $P(0) = 1$ であり、散逸過程により純度は減少する。

---

## 9. GKSL方程式の構造的性質

### 9.1 トレース保存

GKSL生成子 $\mathcal{L}$ はトレースを保存する：

$$
\frac{d}{dt}\mathrm{Tr}[\rho] = \mathrm{Tr}\!\left[\mathcal{L}[\rho]\right] = 0
$$

これは以下から導かれる：
- ハミルトニアン部分：$\mathrm{Tr}[H\rho - \rho H] = 0$（トレースの巡回性）
- 散逸部分：$\mathrm{Tr}[\tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger] = \mathrm{Tr}[\tilde{L}_\alpha^\dagger \tilde{L}_\alpha \rho]$ より、$\mathrm{Tr}[\mathcal{D}[\tilde{L}_\alpha](\rho)] = 0$

### 9.2 正値性の保存

Lindblad-GKS定理により、GKSL形式の生成子 $\mathcal{L}$ から生成される動力学半群 $e^{\mathcal{L}t}$ は、任意の $t \geq 0$ に対して完全正写像（completely positive map）である。したがって $\rho(0) \geq 0$ ならば $\rho(t) \geq 0$ が保証される。

### 9.3 定常状態

長時間極限 $t \to \infty$ では、系は定常状態 $\rho_\mathrm{ss}$ に収束する：

$$
\mathcal{L}[\rho_\mathrm{ss}] = 0
$$

本モデルでは散逸過程が全ての励起状態から基底状態への遷移を含むため、$\rho_\mathrm{ss}$ は $|0,0,0,0\rangle\langle 0,0,0,0|$（全分子が $\mathrm{S}_0$）に近いと期待される（ただしハミルトニアンの非対角項やTTA過程の競合により、厳密にはこの限りではない）。

---

## 10. 各散逸過程の物理的描像

### 10.1 TTA (三重項-三重項消滅)

2つの隣接分子がともに三重項 $\mathrm{T}_1$ にあるとき、一方が励起一重項 $\mathrm{S}_1$ に遷移し、他方が基底状態 $\mathrm{S}_0$ に脱励起する：

$$
\mathrm{T}_1 + \mathrm{T}_1 \xrightarrow{\gamma_\mathrm{TTA}} \mathrm{S}_1 + \mathrm{S}_0
$$

エネルギー保存則は $2E_\mathrm{T} \approx E_\mathrm{S}$ で近似的に満たされる。本パラメータでは $2 \times 1.5 = 3.0 = E_\mathrm{S}$ と厳密に一致する。

この過程がTTA-UC（アップコンバージョン）の鍵であり、低エネルギー三重項から高エネルギー一重項を生成する。

### 10.2 蛍光

$\mathrm{S}_1 \to \mathrm{S}_0$ のスピン許容放射遷移であり、TTAにより生成された $\mathrm{S}_1$ 状態からのアップコンバージョン蛍光の放出に対応する。

### 10.3 燐光

$\mathrm{T}_1 \to \mathrm{S}_0$ のスピン禁制放射遷移であり、三重項の直接的な放射失活チャネルである。速度定数が非常に小さい（$10^{-6}$）ため、TTAと競合する時間スケールでは無視できる程度である。

### 10.4 内部転換 (IC)

$\mathrm{S}_1 \to \mathrm{S}_0$ の非放射遷移であり、TTAにより生成された $\mathrm{S}_1$ のうち蛍光として放出されない分がこの経路で失活する。蛍光と競合する。

### 10.5 項間交差 (ISC)

- $\mathrm{S}_1 \to \mathrm{T}_1$: スピン-軌道相互作用によるスピン反転遷移。$\mathrm{S}_1$ の失活経路の1つであり、三重項の再生産にも寄与する。
- $\mathrm{T}_1 \to \mathrm{S}_0$: 燐光と同じ初期・終状態だが、非放射過程としてモデル化されている。速度定数は燐光と同程度に小さい。

---

## 11. 計算の全体フロー（まとめ）

1. **パラメータ設定**: $E_\mathrm{T}, E_\mathrm{S}, V, \gamma_\mathrm{TTA}, \Gamma_\mathrm{fl}, \Gamma_\mathrm{ph}, k_\mathrm{IC}, k_\mathrm{ISC,ST}, k_\mathrm{ISC,TS}$、$N=4$、$d=3$

2. **ハミルトニアン構成**:
   - $h_\mathrm{local} = \mathrm{diag}(0, E_\mathrm{T}, E_\mathrm{S})$、$3 \times 3$
   - $H_0 = \sum_{i=0}^{3} I_3^{\otimes i} \otimes h_\mathrm{local} \otimes I_3^{\otimes(3-i)}$、$81 \times 81$
   - $H_\mathrm{transfer}$：3ペア分のホッピング項、$81 \times 81$
   - $H = H_0 + H_\mathrm{transfer}$、$81 \times 81$

3. **Lindblad演算子構成**: 26個の $\tilde{L}_\alpha$（$81 \times 81$）

4. **リウヴィリアン超演算子構成**:
   $$\mathcal{L} = -i(I_{81} \otimes H - H^\mathrm{T} \otimes I_{81}) + \sum_{\alpha=1}^{26}\left[\bar{L}_\alpha \otimes \tilde{L}_\alpha - \frac{1}{2}(I_{81} \otimes \tilde{L}_\alpha^\dagger \tilde{L}_\alpha + (\tilde{L}_\alpha^\dagger \tilde{L}_\alpha)^\mathrm{T} \otimes I_{81})\right]$$
   $\mathcal{L}$ は $6561 \times 6561$ の複素行列

5. **初期状態**: $\rho(0) = |1,0,0,1\rangle\langle 1,0,0,1|$、$\mathrm{vec}(\rho(0)) \in \mathbb{C}^{6561}$

6. **時間発展**: `expm_multiply` により $\mathrm{vec}(\rho(t_k)) = e^{\mathcal{L}t_k}\mathrm{vec}(\rho(0))$ を $t_k = 0, 1, \ldots, 100$ で計算

7. **逆ベクトル化**: $\mathrm{vec}(\rho(t_k))$ を $81 \times 81$ 密度行列 $\rho(t_k)$ に逆変換

8. **物理量計算**:
   - $\mathrm{Tr}[\rho(t_k)]$（トレース保存の検証）
   - $N_{\mathrm{S}_0}(t_k), N_{\mathrm{T}_1}(t_k), N_{\mathrm{S}_1}(t_k)$（占有数）
   - $S(t_k) = -\sum_i \lambda_i \ln \lambda_i$（フォン・ノイマンエントロピー）
   - $P(t_k) = \mathrm{Tr}[\rho(t_k)^2]$（純度）

---

## 付録A. Kronecker積の定義と性質

$m \times n$ 行列 $A$ と $p \times q$ 行列 $B$ の Kronecker 積 $A \otimes B$ は $mp \times nq$ 行列であり：

$$
A \otimes B = \begin{pmatrix}
A_{00}B & A_{01}B & \cdots & A_{0,n-1}B\\
A_{10}B & A_{11}B & \cdots & A_{1,n-1}B\\
\vdots & \vdots & \ddots & \vdots\\
A_{m-1,0}B & A_{m-1,1}B & \cdots & A_{m-1,n-1}B
\end{pmatrix}
$$

本文書で用いる主な性質：
- $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$（混合積性質）
- $(A \otimes B)^\mathrm{T} = A^\mathrm{T} \otimes B^\mathrm{T}$
- $(A \otimes B)^\dagger = A^\dagger \otimes B^\dagger$
- $\mathrm{vec}(AXB) = (B^\mathrm{T} \otimes A)\,\mathrm{vec}(X)$（列主序ベクトル化）

## 付録B. 行列指数関数

$n \times n$ 行列 $A$ の行列指数関数は冪級数で定義される：

$$
e^A = \sum_{k=0}^{\infty} \frac{A^k}{k!} = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \cdots
$$

この級数は任意の $A$ に対して絶対収束する。

`scipy.sparse.linalg.expm_multiply(A, v)` は $e^A v$ を全行列指数関数を陽に構成せずに計算する。Al-Mohy & Higham (2011) のアルゴリズムに基づき、Krylov部分空間法を用いて効率的に近似する。$6561 \times 6561$ の行列の全体の指数関数を記憶する必要がないため、メモリ効率が良い。

## 付録C. 具体的な数値例

### C.1 オンサイト・ハミルトニアン $H_0$ の対角要素の例

| インデックス $n$ | 状態 $\|s_0, s_1, s_2, s_3\rangle$ | $(H_0)_{nn}$ (eV) |
|---|---|---|
| 0 | $\|0,0,0,0\rangle$ | $0+0+0+0 = 0$ |
| 1 | $\|0,0,0,1\rangle$ | $0+0+0+1.5 = 1.5$ |
| 2 | $\|0,0,0,2\rangle$ | $0+0+0+3.0 = 3.0$ |
| 28 | $\|1,0,0,1\rangle$ | $1.5+0+0+1.5 = 3.0$ |
| 40 | $\|1,1,1,1\rangle$ | $1.5 \times 4 = 6.0$ |
| 80 | $\|2,2,2,2\rangle$ | $3.0 \times 4 = 12.0$ |

### C.2 初期状態の物理量

| 物理量 | 値 |
|---|---|
| $N_{\mathrm{S}_0}(0)$ | 2.0 |
| $N_{\mathrm{T}_1}(0)$ | 2.0 |
| $N_{\mathrm{S}_1}(0)$ | 0.0 |
| $S(0)$ | 0.0 |
| $P(0)$ | 1.0 |
| $\mathrm{Tr}[\rho(0)]$ | 1.0 |

---

## 参考文献

- V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, "Completely positive dynamical semigroups of N-level systems," J. Math. Phys. **17**, 821 (1976).
- G. Lindblad, "On the generators of quantum dynamical semigroups," Commun. Math. Phys. **48**, 119 (1976).
- A. H. Al-Mohy and N. J. Higham, "Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators," SIAM J. Sci. Comput. **33**, 488 (2011).
