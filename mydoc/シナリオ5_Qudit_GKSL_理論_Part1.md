# シナリオ5: Qudit GKSL（ボソン無し）の省略無し詳細理論 — Part 1

## 物理系の定義・ハミルトニアン・Lindblad演算子・GKSL方程式

本文書は、`tutorials/quantum_dynamics_gksl_comparison.ipynb` の「3. シナリオ5: Qudit GKSL（ボソン無し）」セルの結果を再現するために必要な理論の第1部である。物理系の定義からGKSL方程式の定式化までを、プログラムで計算可能な解像度（行列の具体的構成法を含む）で省略無しに記述する。

シナリオ5はネイティブqutrit（$d=3$）エンコーディングを用いるStinespring拡張+Trotter分解シミュレータである。Part 1では物理系とGKSL方程式の構造（シナリオ1の古典ODEと共通する部分）を定式化し、Part 2以降でStinespring拡張とTrotter分解（シナリオ5に固有の部分）を扱う。

---

## 1. 単一分子の状態空間

各分子は3準位系（qutrit, $d=3$）としてモデル化される。基底は以下の通りである：

$$
|0\rangle = \mathrm{S}_0 \quad (\text{基底一重項状態}), \qquad
|1\rangle = \mathrm{T}_1 \quad (\text{第一励起三重項状態}), \qquad
|2\rangle = \mathrm{S}_1 \quad (\text{第一励起一重項状態})
$$

単一分子のヒルベルト空間は $\mathcal{H}_{\mathrm{mol}} = \mathbb{C}^3$ であり、標準基底ベクトルの列ベクトル表現は：

$$
|0\rangle = \begin{pmatrix}1\\0\\0\end{pmatrix}, \quad
|1\rangle = \begin{pmatrix}0\\1\\0\end{pmatrix}, \quad
|2\rangle = \begin{pmatrix}0\\0\\1\end{pmatrix}
$$

**コード対応**: `gksl_physical_parameters.py` において `d = 3` として定義される。

---

## 2. 多分子系のヒルベルト空間

### 2.1 テンソル積構造

$N = 4$ 個の分子からなる系の全ヒルベルト空間は、テンソル積空間である：

$$
\mathcal{H} = \mathcal{H}_0 \otimes \mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \mathcal{H}_3 = (\mathbb{C}^3)^{\otimes 4}
$$

次元は $\dim \mathcal{H} = 3^4 = 81$ である。密度行列 $\rho$ は $81 \times 81$ のエルミート半正定値行列で、$\mathrm{Tr}[\rho] = 1$ を満たす。

### 2.2 計算基底のラベリング

81個の計算基底状態 $|s_0, s_1, s_2, s_3\rangle$（$s_i \in \{0,1,2\}$）は、Kronecker積の標準的な順序（最左のサイトが最上位桁）に対応する混合基数（radix-3）インデックス $n \in \{0,1,\ldots,80\}$ で一意にラベル付けされる：

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

**コード対応**: `compute_populations_from_density_matrix` では `mol=N-1` から `mol=0` へ逆順にループし、`remainder = idx`, `local_state = remainder % d`, `remainder //= d` で分解する。最初に取り出されるのは $s_3$（最下位桁）、最後が $s_0$（最上位桁）である。

### 2.3 quditエンコーディングにおける禁止状態

シナリオ5ではネイティブqutrit（$d=3$）エンコーディングを使用する。$d=3$ の各qutritは物理状態 $\mathrm{S}_0, \mathrm{T}_1, \mathrm{S}_1$ に自然に対応するため、**禁止状態は存在しない**。81個の計算基底状態は全て物理的に意味のある状態である。

これはqubitエンコーディング（$d=4$、シナリオ3）との本質的な違いである。qubitエンコーディングでは各分子を2 qubit（$d=4$）でエンコードするため、$4^4 = 256$ 次元の空間のうち $3^4 = 81$ 次元のみが物理的で、$256 - 81 = 175$ 個の禁止状態が生じる。

---

## 3. 物理パラメータ

### 3.1 パラメータ一覧

| パラメータ | 記号 | 値 | 単位 | 物理的意味 |
|---|---|---|---|---|
| 三重項エネルギー | $E_{\mathrm{T}}$ | 1.5 | eV | T₁状態のエネルギー |
| 一重項エネルギー | $E_{\mathrm{S}}$ | 3.0 | eV | S₁状態のエネルギー |
| 移動結合 | $V$ | 0.1 | eV | 隣接分子間の三重項エネルギー移動結合 |
| TTA速度 | $\gamma_{\mathrm{TTA}}$ | 0.05 | eV/$\hbar$ | 三重項-三重項消滅速度 |
| 蛍光速度 | $\Gamma_{\mathrm{fl}}$ | 0.01 | eV/$\hbar$ | S₁→S₀放射遷移速度 |
| 燐光速度 | $\Gamma_{\mathrm{ph}}$ | $10^{-6}$ | eV/$\hbar$ | T₁→S₀放射遷移速度 |
| 内部転換速度 | $k_{\mathrm{IC}}$ | 0.005 | eV/$\hbar$ | S₁→S₀内部転換速度 |
| ISC (S→T) 速度 | $k_{\mathrm{ISC,ST}}$ | 0.003 | eV/$\hbar$ | S₁→T₁項間交差速度 |
| ISC (T→S) 速度 | $k_{\mathrm{ISC,TS}}$ | $10^{-5}$ | eV/$\hbar$ | T₁→S₀項間交差速度 |
| 分子数 | $N$ | 4 | — | 鎖状配置の分子数 |
| 局所次元 | $d$ | 3 | — | 各分子の準位数（qutrit） |

### 3.2 自然単位系

本モデルでは $\hbar = 1$（自然単位系）を採用する。速度定数の単位は eV/$\hbar$ であり、自然単位系では eV と等価である。時間の物理単位への変換には $\hbar = 0.6582$ eV·fs を用いるが、シミュレーション内部では無次元時間で計算する。

**コード対応**: `GKSLPhysicalParameters.hbar` プロパティは `1.0` を返す。

### 3.3 隣接分子ペア

$N=4$ の鎖状配置における最近接ペアは：

$$
\mathrm{neighbors} = \{(0,1),\;(1,2),\;(2,3)\}
$$

ペア数は $|\mathrm{neighbors}| = N - 1 = 3$ である。

**コード対応**: `GKSLPhysicalParameters.neighbors` プロパティは `[(i, i+1) for i in range(N-1)]` を返す。

---

## 4. ハミルトニアン

全ハミルトニアンは以下の2つの成分の和である：

$$
H = H_0 + H_{\mathrm{transfer}}
$$

### 4.1 オンサイト・ハミルトニアン $H_0$

#### 4.1.1 定義

各分子のエネルギー準位を記述する。単一分子のオンサイト・ハミルトニアンは：

$$
h_{\mathrm{local}} = E_{\mathrm{T}}|1\rangle\langle 1| + E_{\mathrm{S}}|2\rangle\langle 2| = \begin{pmatrix}0 & 0 & 0\\0 & E_{\mathrm{T}} & 0\\0 & 0 & E_{\mathrm{S}}\end{pmatrix} = \begin{pmatrix}0 & 0 & 0\\0 & 1.5 & 0\\0 & 0 & 3.0\end{pmatrix}
$$

全系のオンサイト・ハミルトニアンは各サイトへの埋め込みの総和である：

$$
H_0 = \sum_{i=0}^{N-1} h_{\mathrm{local}}^{(i)}
$$

ここで $h_{\mathrm{local}}^{(i)}$ はサイト $i$ にのみ作用し、他のサイトには恒等演算子が置かれる：

$$
h_{\mathrm{local}}^{(i)} = \underbrace{I_3 \otimes \cdots \otimes I_3}_{i \text{ 個}} \otimes\; h_{\mathrm{local}} \;\otimes \underbrace{I_3 \otimes \cdots \otimes I_3}_{N-1-i \text{ 個}}
$$

$I_3$ は $3\times 3$ 単位行列である。

#### 4.1.2 構成アルゴリズム

**コード対応**: `build_single_site_operator(op, site, N, d)` は以下のアルゴリズムで全系演算子を構成する：

```
入力: op (d×d行列), site (整数), N (整数), d (整数)
1. op_list = [I_d, I_d, ..., I_d]  (N個)
2. op_list[site] = op
3. result = op_list[0] ⊗ op_list[1] ⊗ ... ⊗ op_list[N-1]
4. 出力: result (d^N × d^N 行列)
```

`build_onsite_hamiltonian(params)` は：

```
1. h_local = diag(0, E_T, E_S)  (3×3)
2. H_0 = 0  (81×81 ゼロ行列)
3. for i = 0 to N-1:
     H_0 += build_single_site_operator(h_local, i, N, d)
4. エルミート性の検証: H_0 == H_0†
5. 出力: H_0
```

#### 4.1.3 行列構造

$H_0$ は $81 \times 81$ の**実対角行列**（エルミートかつ対称）である。対角要素 $(n,n)$ は、インデックス $n$ に対応する状態 $|s_0, s_1, s_2, s_3\rangle$ のエネルギーの総和：

$$
(H_0)_{nn} = \sum_{i=0}^{3} \varepsilon(s_i), \qquad \varepsilon(s) = \begin{cases}0 & s=0 \\ E_{\mathrm{T}} = 1.5 & s=1 \\ E_{\mathrm{S}} = 3.0 & s=2\end{cases}
$$

非対角要素は全てゼロである。

#### 4.1.4 具体的な対角要素の例

| インデックス $n$ | 状態 $\|s_0, s_1, s_2, s_3\rangle$ | $(H_0)_{nn}$ (eV) |
|---|---|---|
| 0 | $\|0,0,0,0\rangle$ | $0+0+0+0 = 0$ |
| 1 | $\|0,0,0,1\rangle$ | $0+0+0+1.5 = 1.5$ |
| 2 | $\|0,0,0,2\rangle$ | $0+0+0+3.0 = 3.0$ |
| 28 | $\|1,0,0,1\rangle$ | $1.5+0+0+1.5 = 3.0$ |
| 40 | $\|1,1,1,1\rangle$ | $1.5 \times 4 = 6.0$ |
| 80 | $\|2,2,2,2\rangle$ | $3.0 \times 4 = 12.0$ |

### 4.2 三重項エネルギー移動ハミルトニアン $H_{\mathrm{transfer}}$

#### 4.2.1 定義

隣接分子間の三重項励起子のホッピングを記述する：

$$
H_{\mathrm{transfer}} = \sum_{\langle i,j \rangle} V \left( |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1| \right)
$$

ここで和は最近接ペア $\langle i,j \rangle \in \{(0,1),(1,2),(2,3)\}$ に渡る。

各項の物理的意味は：
- $|0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0|$：分子 $i$ の三重項が消滅（$\mathrm{T}_1 \to \mathrm{S}_0$）し、分子 $j$ に三重項が生成（$\mathrm{S}_0 \to \mathrm{T}_1$）される（$i \to j$ ホッピング）
- エルミート共役が逆方向ホッピング（$j \to i$）を与える

#### 4.2.2 単一サイト演算子の行列表現

$$
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}, \qquad
|1\rangle\langle 0| = \begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}
$$

#### 4.2.3 全系演算子の構成法

ペア $(i,j)$ に対する順方向ホッピング演算子は、サイト $i$ に $|0\rangle\langle 1|$、サイト $j$ に $|1\rangle\langle 0|$ を置き、その他のサイトに $I_3$ を置いたKronecker積である：

$$
F_{ij} = \underbrace{I_3 \otimes \cdots}_{0 \text{ to } i-1} \otimes\; |0\rangle\langle 1| \;\otimes \underbrace{\cdots \otimes I_3 \otimes \cdots}_{i+1 \text{ to } j-1} \otimes\; |1\rangle\langle 0| \;\otimes \underbrace{I_3 \otimes \cdots}_{j+1 \text{ to } N-1}
$$

各ペアの寄与は $V(F_{ij} + F_{ij}^\dagger)$ であり、$F_{ij}^\dagger = |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1| \otimes \cdots$ なので合わせてエルミート演算子となる。

**コード対応**: `build_transfer_hamiltonian(params)` は以下のアルゴリズムを実行する：

```
1. H_t = 0  (81×81 ゼロ行列)
2. ket0_bra1 = |0><1|  (3×3)
3. ket1_bra0 = |1><0|  (3×3)
4. for (i, j) in neighbors:
     op_list = [I_3] * N
     op_list[i] = ket0_bra1
     op_list[j] = ket1_bra0
     fwd = op_list[0] ⊗ op_list[1] ⊗ ... ⊗ op_list[N-1]
     H_t += V * (fwd + fwd†)
5. エルミート性の検証: H_t == H_t†
6. 出力: H_t
```

$H_{\mathrm{transfer}}$ は $81 \times 81$ のエルミート行列であり、$H_0$ と異なり非対角成分を持つ。

### 4.3 全ハミルトニアン $H$

$$
H = H_0 + H_{\mathrm{transfer}}
$$

$H$ は $81 \times 81$ のエルミート行列である。

**コード対応**: `QuditGKSLSimulator.__init__` において `self.H_total = self.H_0 + self.H_transfer` として構成される。

---

## 5. Lindblad演算子

GKSL方程式における散逸過程は Lindblad 演算子 $\tilde{L}_\alpha$ で記述される。本コードでは各 $\tilde{L}_\alpha$ に $\sqrt{\gamma_\alpha}$ が含まれた形で構成される：

$$
\tilde{L}_\alpha = \sqrt{\gamma_\alpha}\; L_\alpha^{(\mathrm{bare})}
$$

以下では $\tilde{L}_\alpha$ をそのまま「Lindblad演算子」と呼ぶ。

### 5.1 Lindblad演算子の総数

$$
n_{\mathrm{Lindblad}} = 2 \times |\mathrm{neighbors}| + 5 \times N = 2 \times 3 + 5 \times 4 = 26
$$

26個のLindblad演算子はそれぞれ $81 \times 81$ の行列である。

### 5.2 三重項-三重項消滅 (TTA) — 6個

TTA過程 $\mathrm{T}_1 + \mathrm{T}_1 \to \mathrm{S}_1 + \mathrm{S}_0$ を記述する。各隣接ペア $(i,j)$ に対して2つのチャネルがある。速度定数は $\gamma = \gamma_{\mathrm{TTA}}/2 = 0.025$（対称的に2チャネルで分配）。

**チャネル1**: 分子 $i$ が $\mathrm{T}_1 \to \mathrm{S}_1$ に励起され、分子 $j$ が $\mathrm{T}_1 \to \mathrm{S}_0$ に脱励起される。

$$
\tilde{L}_{\mathrm{TTA},1}^{(i,j)} = \sqrt{\gamma}\; \left(|2\rangle_i\langle 1|\right) \otimes \left(|0\rangle_j\langle 1|\right) \otimes \bigotimes_{k \neq i,j} I_3^{(k)}
$$

**チャネル2**: 分子 $i$ が $\mathrm{T}_1 \to \mathrm{S}_0$ に脱励起され、分子 $j$ が $\mathrm{T}_1 \to \mathrm{S}_1$ に励起される。

$$
\tilde{L}_{\mathrm{TTA},2}^{(i,j)} = \sqrt{\gamma}\; \left(|0\rangle_i\langle 1|\right) \otimes \left(|2\rangle_j\langle 1|\right) \otimes \bigotimes_{k \neq i,j} I_3^{(k)}
$$

**単一サイト演算子の行列表現**:

$$
|2\rangle\langle 1| = \begin{pmatrix}0&0&0\\0&0&0\\0&1&0\end{pmatrix}, \qquad
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

3ペア $(0,1),(1,2),(2,3)$ に対してそれぞれ2チャネルで、合計6個の $81 \times 81$ 行列。

**コード対応**: `build_lindblad_operators` の TTA 部分：

```python
for i, j in params.neighbors:
    gamma = params.gamma_TTA / 2.0     # 0.025
    sq = np.sqrt(gamma)                 # sqrt(0.025) ≈ 0.15811
    # Channel 1: |2>_i<1| ⊗ |0>_j<1|
    op_list = [eye] * N
    op_list[i] = _ket_bra(2, 1)        # 3×3
    op_list[j] = _ket_bra(0, 1)        # 3×3
    ops.append((sq * reduce(np.kron, op_list), gamma))
    # Channel 2: |0>_i<1| ⊗ |2>_j<1|
    op_list = [eye] * N
    op_list[i] = _ket_bra(0, 1)
    op_list[j] = _ket_bra(2, 1)
    ops.append((sq * reduce(np.kron, op_list), gamma))
```

### 5.3 蛍光 (Fluorescence) — 4個

$\mathrm{S}_1 \to \mathrm{S}_0$ の放射遷移。各分子 $i \in \{0,1,2,3\}$ に対して：

$$
\tilde{L}_{\mathrm{fl}}^{(i)} = \sqrt{\Gamma_{\mathrm{fl}}}\; |0\rangle_i\langle 2| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：

$$
|0\rangle\langle 2| = \begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}
$$

$\gamma_\alpha = \Gamma_{\mathrm{fl}} = 0.01$、$\sqrt{\gamma_\alpha} = 0.1$。

### 5.4 燐光 (Phosphorescence) — 4個

$\mathrm{T}_1 \to \mathrm{S}_0$ の放射遷移。各分子 $i$ に対して：

$$
\tilde{L}_{\mathrm{ph}}^{(i)} = \sqrt{\Gamma_{\mathrm{ph}}}\; |0\rangle_i\langle 1| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：

$$
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

$\gamma_\alpha = \Gamma_{\mathrm{ph}} = 10^{-6}$、$\sqrt{\gamma_\alpha} = 10^{-3}$。

### 5.5 内部転換 (Internal Conversion, IC) — 4個

$\mathrm{S}_1 \to \mathrm{S}_0$ の非放射遷移。各分子 $i$ に対して：

$$
\tilde{L}_{\mathrm{IC}}^{(i)} = \sqrt{k_{\mathrm{IC}}}\; |0\rangle_i\langle 2| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：$|0\rangle\langle 2|$（蛍光と同じ遷移演算子、異なる速度定数）。

$\gamma_\alpha = k_{\mathrm{IC}} = 0.005$、$\sqrt{\gamma_\alpha} \approx 0.07071$。

### 5.6 項間交差 S→T (ISC S→T) — 4個

$\mathrm{S}_1 \to \mathrm{T}_1$ のスピン禁制遷移。各分子 $i$ に対して：

$$
\tilde{L}_{\mathrm{ISC,ST}}^{(i)} = \sqrt{k_{\mathrm{ISC,ST}}}\; |1\rangle_i\langle 2| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：

$$
|1\rangle\langle 2| = \begin{pmatrix}0&0&0\\0&0&1\\0&0&0\end{pmatrix}
$$

$\gamma_\alpha = k_{\mathrm{ISC,ST}} = 0.003$、$\sqrt{\gamma_\alpha} \approx 0.05477$。

### 5.7 項間交差 T→S (ISC T→S) — 4個

$\mathrm{T}_1 \to \mathrm{S}_0$ のスピン禁制遷移。各分子 $i$ に対して：

$$
\tilde{L}_{\mathrm{ISC,TS}}^{(i)} = \sqrt{k_{\mathrm{ISC,TS}}}\; |0\rangle_i\langle 1| \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

単一サイト演算子：$|0\rangle\langle 1|$（燐光と同じ遷移演算子、異なる速度定数）。

$\gamma_\alpha = k_{\mathrm{ISC,TS}} = 10^{-5}$、$\sqrt{\gamma_\alpha} \approx 0.003162$。

### 5.8 Lindblad演算子の一覧表

| 演算子番号 | 種類 | ペア/サイト | $\gamma_\alpha$ | サイト $i$ の演算子 | サイト $j$ の演算子 | 個数 |
|---|---|---|---|---|---|---|
| 1–6 | TTA ch.1,2 | $(0,1),(1,2),(2,3)$ | $\gamma_{\mathrm{TTA}}/2 = 0.025$ | $\|2\rangle\langle 1\|$ or $\|0\rangle\langle 1\|$ | $\|0\rangle\langle 1\|$ or $\|2\rangle\langle 1\|$ | 6 |
| 7–10 | 蛍光 | $i=0,1,2,3$ | $\Gamma_{\mathrm{fl}} = 0.01$ | $\|0\rangle\langle 2\|$ | — | 4 |
| 11–14 | 燐光 | $i=0,1,2,3$ | $\Gamma_{\mathrm{ph}} = 10^{-6}$ | $\|0\rangle\langle 1\|$ | — | 4 |
| 15–18 | IC | $i=0,1,2,3$ | $k_{\mathrm{IC}} = 0.005$ | $\|0\rangle\langle 2\|$ | — | 4 |
| 19–22 | ISC S→T | $i=0,1,2,3$ | $k_{\mathrm{ISC,ST}} = 0.003$ | $\|1\rangle\langle 2\|$ | — | 4 |
| 23–26 | ISC T→S | $i=0,1,2,3$ | $k_{\mathrm{ISC,TS}} = 10^{-5}$ | $\|0\rangle\langle 1\|$ | — | 4 |

合計：**26個**の $81 \times 81$ Lindblad演算子。

### 5.9 Lindblad演算子の構成に関する注意

`build_lindblad_operators` は `(L_alpha, gamma_alpha)` のタプルのリストを返す。ここで `L_alpha` は既に $\sqrt{\gamma_\alpha}$ を含んだ $81 \times 81$ 行列であり、`gamma_alpha` は速度定数の値（冗長な情報だが参照用に保持される）。返されるリストの長さが $2 \times |\mathrm{neighbors}| + 5 \times N = 26$ であることがアサーションで検証される。

---

## 6. GKSL-Lindblad マスター方程式

### 6.1 方程式の定義

密度行列 $\rho(t)$ の時間発展は GKSL (Gorini-Kossakowski-Sudarshan-Lindblad) マスター方程式に従う：

$$
\frac{d\rho}{dt} = \mathcal{L}[\rho] = -\frac{i}{\hbar}[H, \rho] + \sum_{\alpha=1}^{26} \mathcal{D}[\tilde{L}_\alpha](\rho)
$$

ここで $\hbar = 1$（自然単位系）、各Lindblad散逸子は：

$$
\mathcal{D}[\tilde{L}_\alpha](\rho) = \tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger - \frac{1}{2}\left(\tilde{L}_\alpha^\dagger \tilde{L}_\alpha \rho + \rho \tilde{L}_\alpha^\dagger \tilde{L}_\alpha\right)
$$

### 6.2 ハミルトニアン部分と散逸部分の分離

GKSL生成子 $\mathcal{L}$ は以下の2つの成分に分離できる：

$$
\mathcal{L} = \mathcal{L}_H + \mathcal{L}_D
$$

ここで：

- **ハミルトニアン生成子**: $\mathcal{L}_H[\rho] = -i[H, \rho]$ — ユニタリ時間発展を生成する
- **散逸生成子**: $\mathcal{L}_D[\rho] = \sum_{\alpha=1}^{26} \mathcal{D}[\tilde{L}_\alpha](\rho)$ — 非ユニタリ散逸過程を記述する

この分離はシナリオ5のTrotter分解の基礎となる（Part 3で詳述）。

### 6.3 ハミルトニアン部分の時間発展

ハミルトニアン部分のみの時間発展は、密度行列に対するユニタリ変換として厳密に解ける：

$$
e^{\mathcal{L}_H \Delta t}[\rho] = U_H(\Delta t)\, \rho\, U_H(\Delta t)^\dagger
$$

ここで $U_H(\Delta t) = e^{-iH\Delta t}$ は $81 \times 81$ のユニタリ行列である。

**コード対応**: `QuditGKSLSimulator._precompute_unitaries(dt)` において半ステップ分のユニタリ行列が事前計算される：

$$
U_H^{(1/2)} = e^{-iH\,dt/2}
$$

`self._U_H_half = expm(-1j * self.H_total * dt / 2)` で計算される $81 \times 81$ のユニタリ行列。

### 6.4 散逸部分の個別チャネル分解

散逸生成子はさらに個別のLindblad散逸子に分解できる：

$$
\mathcal{L}_D = \sum_{\alpha=1}^{26} \mathcal{L}_{D_\alpha}
$$

ここで $\mathcal{L}_{D_\alpha}[\rho] = \mathcal{D}[\tilde{L}_\alpha](\rho)$ は第 $\alpha$ Lindbladチャネルの生成子である。

各 $\mathcal{L}_{D_\alpha}$ は個別にCPTP写像の生成子であり、その時間発展 $e^{\mathcal{L}_{D_\alpha} \Delta t}$ はCPTP写像である。シナリオ5では、これらの個別チャネルを Stinespring 拡張で近似する（Part 2で詳述）。

### 6.5 シナリオ1（古典GKSL）との違い

シナリオ1（古典GKSL）では、$\mathcal{L}$ をリウヴィリアン超演算子として $6561 \times 6561$ 行列に明示的にベクトル化し、行列指数関数 `expm_multiply` で時間発展を計算する。

シナリオ5（Qudit GKSL）では、$\mathcal{L}$ を明示的にベクトル化**しない**。代わりに：

1. ハミルトニアン部分はユニタリ行列 $U_H$ の密度行列への左右からの作用で計算
2. 散逸部分は各Lindbladチャネルを個別にStinespring拡張で近似

この方式は量子コンピュータ上で実行可能な量子回路に対応する。

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

**コード対応**: `prepare_initial_state` において `index = 1 * (d ** (N - 1)) + 1 = 1 * 27 + 1 = 28`。ここで `1 * (d ** (N - 1))` は $s_0 = 1$ に対応する項、末尾の `+ 1` は $s_3 = 1$ に対応する項であり、$s_1 = s_2 = 0$ なので中間項は寄与しない。

### 7.3 初期密度行列

$$
\rho(0) = |\psi_0\rangle\langle\psi_0|
$$

これは $81 \times 81$ の行列で、$(28, 28)$ 成分のみが1、他は全てゼロである：

$$
\rho(0)_{mn} = \delta_{m,28}\,\delta_{n,28}
$$

**コード対応**: `prepare_initial_state` は `np.outer(psi, psi.conj())` で密度行列を構成する。`psi` は長さ81のベクトルで、`psi[28] = 1.0`、他は全て0。

### 7.4 初期状態の物理量

| 物理量 | 値 |
|---|---|
| $N_{\mathrm{S}_0}(0)$ | 2.0 |
| $N_{\mathrm{T}_1}(0)$ | 2.0 |
| $N_{\mathrm{S}_1}(0)$ | 0.0 |
| $S(0)$ (フォン・ノイマンエントロピー) | 0.0 |
| $P(0)$ (純度) | 1.0 |
| $\mathrm{Tr}[\rho(0)]$ | 1.0 |

---

## 8. 物理量の抽出

各時刻 $t_k$ における密度行列 $\rho(t_k)$ から以下の物理量を計算する。

### 8.1 トレース（規格化の検証）

$$
\mathrm{Tr}[\rho(t_k)] = \sum_{n=0}^{80} \rho_{nn}(t_k)
$$

CPTP性により $\mathrm{Tr}[\rho(t_k)] = 1$ が構造的に保証される（ただしStinespring近似の誤差範囲内で）。

**コード対応**: `float(np.real(np.trace(rho)))`

### 8.2 占有数（ポピュレーション）

各電子状態（$\mathrm{S}_0$, $\mathrm{T}_1$, $\mathrm{S}_1$）の平均占有数を密度行列の対角要素から計算する：

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

ここで $s_i$ はインデックス $n$ の混合基数分解における分子 $i$ の局所状態である（第2.2節参照）。

規格化の整合性から：

$$
N_{\mathrm{S}_0}(t) + N_{\mathrm{T}_1}(t) + N_{\mathrm{S}_1}(t) = N \cdot \mathrm{Tr}[\rho(t)] = 4
$$

**コード対応**: `compute_populations_from_density_matrix(rho, params)` は各対角要素 `diag[idx]` に対して混合基数分解を行い、各分子の局所状態に応じて `N_S0`, `N_T1`, `N_S1` に加算する。さらに `per_molecule_populations` として各分子 $i$ の状態別占有数も計算する。

### 8.3 フォン・ノイマンエントロピー

$$
S(t) = -\mathrm{Tr}[\rho(t)\ln\rho(t)] = -\sum_{i=0}^{80} \lambda_i(t) \ln\lambda_i(t)
$$

ここで $\{\lambda_i(t)\}$ は $\rho(t)$ の固有値であり、$\lambda_i > \mathrm{tol}$（$\mathrm{tol} = 10^{-12}$）の項のみを含める。

**コード対応**: `compute_von_neumann_entropy(rho, tol=1e-12)` は `np.linalg.eigvalsh(rho)` で固有値を計算し、$\mathrm{tol}$ より大きい固有値のみで $-\sum \lambda_i \ln \lambda_i$ を計算する。

純粋状態のとき $S = 0$、最大混合状態のとき $S = \ln 81 \approx 4.394$。初期状態は純粋状態なので $S(0) = 0$ であり、散逸過程によりエントロピーは増大する。

### 8.4 純度（Purity）

$$
P(t) = \mathrm{Tr}[\rho(t)^2] = \sum_{m,n} |\rho_{mn}(t)|^2
$$

**コード対応**: `compute_purity(rho)` は `np.real(np.trace(rho @ rho))` を計算する。

純粋状態のとき $P = 1$、最大混合状態のとき $P = 1/81 \approx 0.01235$。初期状態は純粋状態なので $P(0) = 1$ であり、散逸過程により純度は減少する。

---

## 9. GKSL方程式の構造的性質

### 9.1 トレース保存

GKSL生成子 $\mathcal{L}$ はトレースを保存する：

$$
\frac{d}{dt}\mathrm{Tr}[\rho] = \mathrm{Tr}\!\left[\mathcal{L}[\rho]\right] = 0
$$

これは以下から導かれる：
- ハミルトニアン部分：$\mathrm{Tr}[H\rho - \rho H] = 0$（トレースの巡回性）
- 散逸部分：$\mathrm{Tr}[\tilde{L}_\alpha \rho \tilde{L}_\alpha^\dagger] = \mathrm{Tr}[\tilde{L}_\alpha^\dagger \tilde{L}_\alpha \rho]$（トレースの巡回性）より $\mathrm{Tr}[\mathcal{D}[\tilde{L}_\alpha](\rho)] = 0$

### 9.2 正値性の保存

Lindblad-GKS定理により、GKSL形式の生成子 $\mathcal{L}$ から生成される動力学半群 $e^{\mathcal{L}t}$ は、任意の $t \geq 0$ に対して完全正写像（completely positive map）である。したがって $\rho(0) \geq 0$ ならば $\rho(t) \geq 0$ が保証される。

### 9.3 各散逸過程の物理的描像

#### TTA (三重項-三重項消滅)

2つの隣接分子がともに $\mathrm{T}_1$ にあるとき、一方が $\mathrm{S}_1$ に励起され他方が $\mathrm{S}_0$ に脱励起する：

$$
\mathrm{T}_1 + \mathrm{T}_1 \xrightarrow{\gamma_{\mathrm{TTA}}} \mathrm{S}_1 + \mathrm{S}_0
$$

エネルギー保存則は $2E_{\mathrm{T}} \approx E_{\mathrm{S}}$ で近似的に満たされる。本パラメータでは $2 \times 1.5 = 3.0 = E_{\mathrm{S}}$ と厳密に一致する。

#### 蛍光

$\mathrm{S}_1 \to \mathrm{S}_0$ のスピン許容放射遷移であり、TTAにより生成された $\mathrm{S}_1$ 状態からのアップコンバージョン蛍光の放出に対応する。

#### 燐光

$\mathrm{T}_1 \to \mathrm{S}_0$ のスピン禁制放射遷移であり、三重項の直接的な放射失活チャネル。$\Gamma_{\mathrm{ph}} = 10^{-6}$ と非常に小さいため、TTAと競合する時間スケールでは無視できる。

#### 内部転換 (IC)

$\mathrm{S}_1 \to \mathrm{S}_0$ の非放射遷移。TTAにより生成された $\mathrm{S}_1$ のうち蛍光として放出されない分がこの経路で失活する。

#### 項間交差 (ISC)

- $\mathrm{S}_1 \to \mathrm{T}_1$ ($k_{\mathrm{ISC,ST}} = 0.003$): スピン-軌道相互作用によるスピン反転遷移。$\mathrm{S}_1$ の失活経路の1つであり、三重項の再生産にも寄与する。
- $\mathrm{T}_1 \to \mathrm{S}_0$ ($k_{\mathrm{ISC,TS}} = 10^{-5}$): 非放射過程としてモデル化。速度定数は非常に小さい。

---

## 参考文献

- V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, "Completely positive dynamical semigroups of N-level systems," J. Math. Phys. **17**, 821 (1976).
- G. Lindblad, "On the generators of quantum dynamical semigroups," Commun. Math. Phys. **48**, 119 (1976).

---

*Part 2 では、Stinespring拡張の理論と各Lindbladチャネルの量子回路表現を扱う。*
*Part 3 では、Trotter-Suzuki分解、回文積順序、シミュレーションアルゴリズム全体、誤差解析、ゲート数見積を扱う。*
