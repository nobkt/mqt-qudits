# シナリオ5: Qudit GKSL（ボソン無し）— 理論の詳細

本文書は `tutorials/quantum_dynamics_gksl_comparison.ipynb` の「3. シナリオ5: Qudit GKSL（ボソン無し）」セルの計算で使われている理論を、ソースコード（`qudit_gksl_simulator.py`, `stinespring_utils.py`, `gksl_math_utils.py`, `gksl_physical_parameters.py`, `qudit_gksl_circuit_simulator.py`）に基づき、省略無しで記述する。

---

## 1. 物理モデル: TTA-UC現象

シナリオ5は、TTA-UC（Triplet-Triplet Annihilation Upconversion）現象のGKSL-Lindblad量子ダイナミクスを、ネイティブqutrit（$d=3$）エンコーディングによりシミュレートする。

### 1.1 分子とヒルベルト空間

$N=4$ 個の分子を考える。各分子は3準位系（qutrit, $d=3$）であり、以下の電子状態を持つ：

| 量子状態 | ラベル | 記号 | エネルギー |
|---------|--------|------|-----------|
| $\|0\rangle$ | 基底一重項 | $\mathrm{S}_0$ | $0$ |
| $\|1\rangle$ | 三重項 | $\mathrm{T}_1$ | $E_T$ |
| $\|2\rangle$ | 励起一重項 | $\mathrm{S}_1$ | $E_S$ |

系全体のヒルベルト空間は4個のqutritのテンソル積空間であり、次元は

$$
\dim \mathcal{H} = d^N = 3^4 = 81
$$

である。計算基底は

$$
|s_0, s_1, s_2, s_3\rangle = |s_0\rangle \otimes |s_1\rangle \otimes |s_2\rangle \otimes |s_3\rangle, \quad s_k \in \{0, 1, 2\}
$$

で与えられ、添字の対応は big-endian 順序

$$
\mathrm{index} = s_0 \cdot d^3 + s_1 \cdot d^2 + s_2 \cdot d + s_3 = s_0 \cdot 27 + s_1 \cdot 9 + s_2 \cdot 3 + s_3
$$

による。qutritエンコーディングでは、qubitエンコーディングと異なり**禁止状態が存在しない**。

### 1.2 物理パラメータ

コードで使用されるデフォルトパラメータ（自然単位系 $\hbar = 1$）：

| パラメータ | 記号 | 値 | 物理的意味 |
|-----------|------|-----|-----------|
| 三重項エネルギー | $E_T$ | $1.5\;\mathrm{eV}$ | T₁状態のエネルギー |
| 励起一重項エネルギー | $E_S$ | $3.0\;\mathrm{eV}$ | S₁状態のエネルギー |
| 移動結合定数 | $V$ | $0.1\;\mathrm{eV}$ | 隣接分子間の三重項励起子移動 |
| TTA率 | $\gamma_{\mathrm{TTA}}$ | $0.05$ | 三重項-三重項消滅の散逸率 |
| 蛍光率 | $\Gamma_{\mathrm{fl}}$ | $0.01$ | S₁→S₀放射遷移率 |
| 燐光率 | $\Gamma_{\mathrm{ph}}$ | $10^{-6}$ | T₁→S₀放射遷移率 |
| 内部転換率 | $k_{\mathrm{IC}}$ | $0.005$ | S₁→S₀非放射遷移率 |
| 項間交差率 (S→T) | $k_{\mathrm{ISC,ST}}$ | $0.003$ | S₁→T₁遷移率 |
| 項間交差率 (T→S) | $k_{\mathrm{ISC,TS}}$ | $10^{-5}$ | T₁→S₀遷移率 |

エネルギー保存条件 $2E_T \approx E_S$ が満たされている（TTA過程で2つのT₁が1つのS₁を生成する条件）。

---

## 2. GKSL-Lindblad マスター方程式

系の密度行列 $\rho(t)$ の時間発展は GKSL（Gorini–Kossakowski–Sudarshan–Lindblad）方程式

$$
\frac{d\rho}{dt} = \mathcal{L}(\rho) = -i[H, \rho] + \sum_{\alpha=1}^{26} \left( L_\alpha \rho L_\alpha^\dagger - \frac{1}{2}\{L_\alpha^\dagger L_\alpha, \rho\} \right)
$$

に従う。ここで $\hbar = 1$（自然単位系）である。

右辺は2つの部分に分解される：

$$
\mathcal{L}(\rho) = \mathcal{L}_H(\rho) + \mathcal{L}_D(\rho)
$$

- **ハミルトニアン部分**（ユニタリ発展）:

$$
\mathcal{L}_H(\rho) = -i[H, \rho]
$$

- **散逸部分**（Lindblad散逸子）:

$$
\mathcal{L}_D(\rho) = \sum_{\alpha=1}^{26} \mathcal{D}_\alpha(\rho) = \sum_{\alpha=1}^{26} \left( L_\alpha \rho L_\alpha^\dagger - \frac{1}{2}\{L_\alpha^\dagger L_\alpha, \rho\} \right)
$$

---

## 3. ハミルトニアン

### 3.1 全体構造

全ハミルトニアンは局所（オンサイト）部分と移動（ホッピング）部分の和である：

$$
H = H_0 + H_{\mathrm{transfer}}
$$

### 3.2 オンサイトハミルトニアン $H_0$

各分子 $i$ の局所ハミルトニアンは

$$
h_{\mathrm{local}} = E_T |1\rangle\langle 1| + E_S |2\rangle\langle 2| = \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix}
$$

であり、全系でのオンサイトハミルトニアンはこれらのテンソル積埋め込みの和：

$$
H_0 = \sum_{i=0}^{3} \underbrace{I_d \otimes \cdots \otimes I_d}_{i\text{個}} \otimes\; h_{\mathrm{local}} \;\otimes \underbrace{I_d \otimes \cdots \otimes I_d}_{(3-i)\text{個}}
$$

$H_0$ は $81 \times 81$ の対角行列であり、計算基底 $|s_0, s_1, s_2, s_3\rangle$ に対する対角要素は

$$
\langle s_0, s_1, s_2, s_3 | H_0 | s_0, s_1, s_2, s_3 \rangle = \sum_{i=0}^{3} \varepsilon(s_i)
$$

ただし

$$
\varepsilon(s) = \begin{cases} 0 & (s=0, \;\mathrm{S}_0) \\ E_T & (s=1, \;\mathrm{T}_1) \\ E_S & (s=2, \;\mathrm{S}_1) \end{cases}
$$

### 3.3 移動ハミルトニアン $H_{\mathrm{transfer}}$

隣接分子対 $\langle i,j \rangle \in \{(0,1),\, (1,2),\, (2,3)\}$ 間の三重項励起子ホッピングを記述する：

$$
H_{\mathrm{transfer}} = V \sum_{\langle i,j \rangle} \left( |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1| \right) \otimes I_{\mathrm{others}}
$$

ここで $I_{\mathrm{others}}$ は対 $(i,j)$ 以外の全分子に対する恒等演算子のテンソル積である。

各隣接対の局所移動ハミルトニアンは $d^2 \times d^2 = 9 \times 9$ 行列：

$$
H_{\mathrm{pair}}^{(i,j)} = V \left( |01\rangle\langle 10| + |10\rangle\langle 01| \right)
$$

2-qutrit基底 $\{|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle\}$（添字 $= s_i \cdot 3 + s_j$）での行列表現：

$$
H_{\mathrm{pair}}^{(i,j)} = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & V & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & V & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

非零要素は $(|01\rangle, |10\rangle)$ と $(|10\rangle, |01\rangle)$ の位置（行列の (1,3) 成分と (3,1) 成分、0始まり）のみ。$H_{\mathrm{transfer}}$ はエルミートである。

---

## 4. リンドブラッド演算子

全26個のリンドブラッド演算子は以下の6カテゴリに分類される。**各 $L_\alpha$ は $\sqrt{\gamma_\alpha}$ 因子を含む**（コードの規約）。

### 4.1 TTA（三重項-三重項消滅）: 6演算子

各隣接対 $(i,j) \in \{(0,1),\, (1,2),\, (2,3)\}$ に対して2チャネル：

$$
L_{\mathrm{TTA},1}^{(i,j)} = \sqrt{\frac{\gamma_{\mathrm{TTA}}}{2}} \left( |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1| \right) \otimes I_{\mathrm{others}}
$$

$$
L_{\mathrm{TTA},2}^{(i,j)} = \sqrt{\frac{\gamma_{\mathrm{TTA}}}{2}} \left( |0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1| \right) \otimes I_{\mathrm{others}}
$$

物理的意味：2つのT₁状態から1つのS₁と1つのS₀を生成する。チャネル1は分子 $i$ がS₁に遷移し分子 $j$ がS₀に戻る過程、チャネル2はその逆。2チャネルの合計率は $\gamma_{\mathrm{TTA}}$ であり、各チャネルの率は $\gamma_{\mathrm{TTA}}/2$ である。

局所演算子（$9 \times 9$, 対 $(i,j)$ のテンソル積空間）の行列表現：

**チャネル1**: $L_{\mathrm{TTA},1}^{\mathrm{local}} = \sqrt{\gamma_{\mathrm{TTA}}/2} \; |20\rangle\langle 11|$

$$
\left(L_{\mathrm{TTA},1}^{\mathrm{local}}\right)_{mn} = \begin{cases} \sqrt{\gamma_{\mathrm{TTA}}/2} & (m = 6,\; n = 4) \\ 0 & (\text{otherwise}) \end{cases}
$$

ここで $|20\rangle$ の添字 $= 2 \cdot 3 + 0 = 6$、$|11\rangle$ の添字 $= 1 \cdot 3 + 1 = 4$ である。

**チャネル2**: $L_{\mathrm{TTA},2}^{\mathrm{local}} = \sqrt{\gamma_{\mathrm{TTA}}/2} \; |02\rangle\langle 11|$

$$
\left(L_{\mathrm{TTA},2}^{\mathrm{local}}\right)_{mn} = \begin{cases} \sqrt{\gamma_{\mathrm{TTA}}/2} & (m = 2,\; n = 4) \\ 0 & (\text{otherwise}) \end{cases}
$$

ここで $|02\rangle$ の添字 $= 0 \cdot 3 + 2 = 2$ である。

### 4.2 蛍光（Fluorescence, S₁→S₀）: 4演算子

各分子 $i \in \{0,1,2,3\}$ に対して：

$$
L_{\mathrm{fl}}^{(i)} = \sqrt{\Gamma_{\mathrm{fl}}} \; |0\rangle_i\langle 2| \otimes I_{\mathrm{others}}
$$

局所演算子（$3 \times 3$）：

$$
L_{\mathrm{fl}}^{\mathrm{local}} = \sqrt{\Gamma_{\mathrm{fl}}} \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

### 4.3 燐光（Phosphorescence, T₁→S₀）: 4演算子

各分子 $i \in \{0,1,2,3\}$ に対して：

$$
L_{\mathrm{ph}}^{(i)} = \sqrt{\Gamma_{\mathrm{ph}}} \; |0\rangle_i\langle 1| \otimes I_{\mathrm{others}}
$$

局所演算子（$3 \times 3$）：

$$
L_{\mathrm{ph}}^{\mathrm{local}} = \sqrt{\Gamma_{\mathrm{ph}}} \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

### 4.4 内部転換（Internal Conversion, S₁→S₀）: 4演算子

各分子 $i \in \{0,1,2,3\}$ に対して：

$$
L_{\mathrm{IC}}^{(i)} = \sqrt{k_{\mathrm{IC}}} \; |0\rangle_i\langle 2| \otimes I_{\mathrm{others}}
$$

局所演算子（$3 \times 3$）：

$$
L_{\mathrm{IC}}^{\mathrm{local}} = \sqrt{k_{\mathrm{IC}}} \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

### 4.5 項間交差 S₁→T₁（ISC）: 4演算子

各分子 $i \in \{0,1,2,3\}$ に対して：

$$
L_{\mathrm{ISC,ST}}^{(i)} = \sqrt{k_{\mathrm{ISC,ST}}} \; |1\rangle_i\langle 2| \otimes I_{\mathrm{others}}
$$

局所演算子（$3 \times 3$）：

$$
L_{\mathrm{ISC,ST}}^{\mathrm{local}} = \sqrt{k_{\mathrm{ISC,ST}}} \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}
$$

### 4.6 項間交差 T₁→S₀（ISC）: 4演算子

各分子 $i \in \{0,1,2,3\}$ に対して：

$$
L_{\mathrm{ISC,TS}}^{(i)} = \sqrt{k_{\mathrm{ISC,TS}}} \; |0\rangle_i\langle 1| \otimes I_{\mathrm{others}}
$$

局所演算子（$3 \times 3$）：

$$
L_{\mathrm{ISC,TS}}^{\mathrm{local}} = \sqrt{k_{\mathrm{ISC,TS}}} \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

### 4.7 演算子数のまとめ

| カテゴリ | 演算子数 | 局所次元 |
|---------|---------|---------|
| TTA | $3 \times 2 = 6$ | $9 \times 9$（対空間） |
| 蛍光 | $4$ | $3 \times 3$（単一サイト） |
| 燐光 | $4$ | $3 \times 3$ |
| 内部転換 | $4$ | $3 \times 3$ |
| ISC S₁→T₁ | $4$ | $3 \times 3$ |
| ISC T₁→S₀ | $4$ | $3 \times 3$ |
| **合計** | **26** | — |

---

## 5. Stinespring拡張（Stinespring Dilation）

各リンドブラッド散逸チャネル $\mathcal{D}_\alpha$ を量子回路で実装するため、Stinespring拡張を用いる。

### 5.1 構成法

リンドブラッド演算子 $L_\alpha$（$\sqrt{\gamma}$ 因子を含む $d_{\mathrm{sys}} \times d_{\mathrm{sys}}$ 行列）に対し、以下の手順でStinespringユニタリを構成する。

**Step 1: 生成子の構成**

$2d_{\mathrm{sys}} \times 2d_{\mathrm{sys}}$ のエルミート行列

$$
G_\alpha = \begin{pmatrix} 0 & L_\alpha^\dagger \\ L_\alpha & 0 \end{pmatrix}
$$

を構成する。$G_\alpha$ がエルミートであることは $(L_\alpha^\dagger)^\dagger = L_\alpha$ から直ちに従う。

**Step 2: Stinespringユニタリ**

パラメータ $\theta = \sqrt{\Delta t}$ として

$$
U_\alpha = \exp(-i\theta G_\alpha) = \exp\!\left(-i\sqrt{\Delta t}\; G_\alpha\right)
$$

$U_\alpha$ は $2d_{\mathrm{sys}} \times 2d_{\mathrm{sys}}$ のユニタリ行列であり、環境（ancilla）⊗ 系の結合空間に作用する。

**Step 3: チャネルの適用**

環境を $|0\rangle_{\mathrm{env}}$ に初期化した拡張密度行列

$$
\rho_{\mathrm{ext}} = |0\rangle\langle 0|_{\mathrm{env}} \otimes \rho_{\mathrm{sys}}
$$

にユニタリを適用し、環境を部分トレースで消去する：

$$
\mathcal{E}_\alpha(\rho) = \mathrm{Tr}_{\mathrm{env}}\!\left[ U_\alpha \left(|0\rangle\langle 0|_{\mathrm{env}} \otimes \rho\right) U_\alpha^\dagger \right]
$$

環境が2準位（$|0\rangle$, $|1\rangle$）であるため、部分トレースは

$$
\mathcal{E}_\alpha(\rho) = \langle 0|_{\mathrm{env}} U_\alpha |0\rangle_{\mathrm{env}} \; \rho \; \langle 0|_{\mathrm{env}} U_\alpha^\dagger |0\rangle_{\mathrm{env}} + \langle 1|_{\mathrm{env}} U_\alpha |0\rangle_{\mathrm{env}} \; \rho \; \langle 0|_{\mathrm{env}} U_\alpha^\dagger |1\rangle_{\mathrm{env}}
$$

$$
= K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger
$$

### 5.2 Kraus演算子

$U_\alpha$ をブロック分解すると

$$
U_\alpha = \begin{pmatrix} A & B \\ C & D \end{pmatrix}
$$

ここで各ブロックは $d_{\mathrm{sys}} \times d_{\mathrm{sys}}$ 行列。Kraus演算子は

$$
K_0 = A = U_\alpha[0:d_{\mathrm{sys}},\; 0:d_{\mathrm{sys}}] \quad (\text{ancilla } |0\rangle \to |0\rangle)
$$

$$
K_1 = C = U_\alpha[d_{\mathrm{sys}}:2d_{\mathrm{sys}},\; 0:d_{\mathrm{sys}}] \quad (\text{ancilla } |0\rangle \to |1\rangle)
$$

であり、$U_\alpha$ のユニタリ性からトレース保存条件

$$
K_0^\dagger K_0 + K_1^\dagger K_1 = I
$$

が保証される。

### 5.3 特異値分解による厳密公式

$L_\alpha$ の特異値分解（SVD）を $L_\alpha = W \Sigma V^\dagger$ とする。ここで $W$, $V$ はユニタリ行列、$\Sigma = \mathrm{diag}(\sigma_1, \sigma_2, \ldots, \sigma_{d_{\mathrm{sys}}})$ は特異値の対角行列である。

このとき $G_\alpha$ は以下のように対角化される：

$$
\begin{pmatrix} V & 0 \\ 0 & W \end{pmatrix}^\dagger G_\alpha \begin{pmatrix} V & 0 \\ 0 & W \end{pmatrix} = \begin{pmatrix} 0 & \Sigma \\ \Sigma & 0 \end{pmatrix}
$$

各特異値 $\sigma_k$ に対応する $2 \times 2$ ブロックが $\sigma_k \sigma_x$ の形であるため、

$$
\exp\!\left(-i\theta \begin{pmatrix} 0 & \sigma_k \\ \sigma_k & 0 \end{pmatrix}\right) = \begin{pmatrix} \cos(\theta\sigma_k) & -i\sin(\theta\sigma_k) \\ -i\sin(\theta\sigma_k) & \cos(\theta\sigma_k) \end{pmatrix}
$$

全体として

$$
U_\alpha = \begin{pmatrix} V & 0 \\ 0 & W \end{pmatrix} \begin{pmatrix} \cos(\theta\Sigma) & -i\sin(\theta\Sigma) \\ -i\sin(\theta\Sigma) & \cos(\theta\Sigma) \end{pmatrix} \begin{pmatrix} V^\dagger & 0 \\ 0 & W^\dagger \end{pmatrix}
$$

ブロック行列として展開すると

$$
U_\alpha = \begin{pmatrix} V\cos(\theta\Sigma)V^\dagger & -iV\sin(\theta\Sigma)W^\dagger \\ -iW\sin(\theta\Sigma)V^\dagger & W\cos(\theta\Sigma)W^\dagger \end{pmatrix}
$$

ここで $\theta = \sqrt{\Delta t}$ であり、$\cos(\theta\Sigma) = \mathrm{diag}(\cos(\theta\sigma_1), \ldots, \cos(\theta\sigma_{d_{\mathrm{sys}}}))$、$\sin(\theta\Sigma) = \mathrm{diag}(\sin(\theta\sigma_1), \ldots, \sin(\theta\sigma_{d_{\mathrm{sys}}}))$ である。

Kraus演算子の厳密形は

$$
K_0 = V \cos\!\left(\sqrt{\Delta t}\;\Sigma\right) V^\dagger
$$

$$
K_1 = -iW \sin\!\left(\sqrt{\Delta t}\;\Sigma\right) V^\dagger
$$

### 5.4 ランク1ジャンプ演算子の場合の厳密形

本シミュレーションのリンドブラッド演算子はすべてランク1であり、$L_\alpha = \sqrt{\gamma} |a\rangle\langle b|$ の形を持つ。この場合、特異値は $\sigma_1 = \sqrt{\gamma}$ のみ非零であり、$W$ の第1列は $|a\rangle$、$V$ の第1列は $|b\rangle$ である。

Kraus演算子は

$$
K_0 = I - \left(1 - \cos\!\left(\sqrt{\gamma\,\Delta t}\right)\right) |b\rangle\langle b|
$$

$$
K_1 = -i\sin\!\left(\sqrt{\gamma\,\Delta t}\right) |a\rangle\langle b|
$$

物理的解釈：
- $K_0$: 系が遷移しない場合。ソース状態 $|b\rangle$ の振幅が $\cos(\sqrt{\gamma\,\Delta t})$ に縮小される。
- $K_1$: 系がジャンプする場合。ソース $|b\rangle$ からターゲット $|a\rangle$ への振幅 $-i\sin(\sqrt{\gamma\,\Delta t})$ による遷移。

1ステップあたりのジャンプ確率は

$$
p_{\mathrm{jump}} = \sin^2\!\left(\sqrt{\gamma\,\Delta t}\right) \approx \gamma\,\Delta t \quad (\Delta t \to 0)
$$

### 5.5 GKSL散逸子の1次近似としての正当性

$\Delta t \to 0$ の極限で

$$
\cos(\sqrt{\gamma\,\Delta t}) \approx 1 - \frac{\gamma\,\Delta t}{2}, \quad \sin(\sqrt{\gamma\,\Delta t}) \approx \sqrt{\gamma\,\Delta t}
$$

より

$$
K_0 \approx I - \frac{\Delta t}{2} L^\dagger L, \quad K_1 \approx -i\sqrt{\Delta t}\; L
$$

したがって

$$
\mathcal{E}_\alpha(\rho) = K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger
$$

$$
\approx \left(I - \frac{\Delta t}{2} L^\dagger L\right) \rho \left(I - \frac{\Delta t}{2} L^\dagger L\right) + \Delta t \; L \rho L^\dagger
$$

$$
= \rho - \frac{\Delta t}{2} L^\dagger L \rho - \frac{\Delta t}{2} \rho L^\dagger L + \Delta t \; L \rho L^\dagger + O(\Delta t^2)
$$

$$
= \rho + \Delta t \left( L \rho L^\dagger - \frac{1}{2}\{L^\dagger L, \rho\} \right) + O(\Delta t^2)
$$

これはGKSL散逸子 $\mathcal{D}_\alpha(\rho)$ の1次精度の近似 $\rho + \Delta t \cdot \mathcal{D}_\alpha(\rho)$ と一致する。すなわち、Stinespring拡張は各散逸チャネルの**1次近似**である。

---

## 6. 対称Trotter分解（Strang分割）

### 6.1 分割構造

GKSL Liouvillian $\mathcal{L} = \mathcal{L}_H + \sum_\alpha \mathcal{L}_{D_\alpha}$ を、ハミルトニアン部分と散逸部分に分割し、Strang（2次対称）分割を適用する。

1ステップ $\Delta t$ の発展は：

$$
e^{\mathcal{L}\Delta t} \approx e^{\mathcal{L}_H \Delta t/2} \circ e^{\mathcal{L}_D \Delta t} \circ e^{\mathcal{L}_H \Delta t/2}
$$

ここで $e^{\mathcal{L}_D \Delta t}$ は全散逸チャネルの合成を表す。

### 6.2 回文順序（Palindromic Ordering）による散逸チャネルの合成

散逸部分 $e^{\mathcal{L}_D \Delta t}$ は26個のチャネルの積として近似されるが、個別チャネル間の非可換性によるLie-Trotter積誤差を軽減するため、**回文（パリンドローム）順序**が採用される：

$$
e^{\mathcal{L}_D \Delta t} \approx \left[\prod_{\alpha=1}^{26} \mathcal{E}_\alpha(\Delta t/2)\right] \circ \left[\prod_{\alpha=26}^{1} \mathcal{E}_\alpha(\Delta t/2)\right]
$$

すなわち、全チャネルを順方向に半ステップ $\Delta t/2$ で適用した後、逆方向に半ステップ $\Delta t/2$ で適用する。

### 6.3 完全なTrotterステップ

以上をまとめると、1 Trotterステップの完全な操作は：

$$
\rho(t + \Delta t) \approx \mathcal{T}(\rho(t))
$$

$$
\mathcal{T} = \mathcal{U}_H^{(1/2)} \circ \left[\prod_{\alpha=1}^{26} \mathcal{E}_\alpha^{(1/2)}\right] \circ \left[\prod_{\alpha=26}^{1} \mathcal{E}_\alpha^{(1/2)}\right] \circ \mathcal{U}_H^{(1/2)}
$$

具体的な手順は以下の通り：

1. **ハミルトニアン半ステップ**:

$$
\rho \leftarrow U_H(\Delta t/2) \; \rho \; U_H(\Delta t/2)^\dagger, \quad U_H(\Delta t/2) = e^{-iH\,\Delta t/2}
$$

2. **散逸チャネル順方向半スイープ**（$\alpha = 1, 2, \ldots, 26$、各 $\Delta t/2$）:

$$
\rho \leftarrow \mathcal{E}_\alpha(\Delta t/2)(\rho) = \mathrm{Tr}_{\mathrm{env}}\!\left[U_{\alpha}^{(1/2)} \left(|0\rangle\langle 0| \otimes \rho\right) {U_{\alpha}^{(1/2)}}^\dagger\right]
$$

ここで $U_{\alpha}^{(1/2)} = \exp\!\left(-i\sqrt{\Delta t/2}\; G_\alpha\right)$

3. **散逸チャネル逆方向半スイープ**（$\alpha = 26, 25, \ldots, 1$、各 $\Delta t/2$）:

$$
\rho \leftarrow \mathcal{E}_\alpha(\Delta t/2)(\rho)
$$

4. **ハミルトニアン半ステップ**:

$$
\rho \leftarrow U_H(\Delta t/2) \; \rho \; U_H(\Delta t/2)^\dagger
$$

### 6.4 コードにおける実装

`QuditGKSLSimulator` では、ハミルトニアンユニタリ $U_H(\Delta t/2)$ は全系 $81 \times 81$ の行列指数関数として一括計算される：

$$
U_H(\Delta t/2) = \exp\!\left(-i(H_0 + H_{\mathrm{transfer}}) \cdot \frac{\Delta t}{2}\right) \in \mathbb{C}^{81 \times 81}
$$

Stinespringユニタリも全系レベル（$162 \times 162$）で計算される。回路シミュレータ（`QuditGKSLCircuitSimulator`）では後述の局所ゲートに分解される。

---

## 7. 量子ゲート表現

回路シミュレータ（`QuditGKSLCircuitSimulator`）は、上記のTrotterステップを局所的な量子ゲートに分解する。各ゲートは MQT-Qudits フレームワークの API で実装される。

### 7.1 オンサイトハミルトニアンゲート: `cu_one`（$3 \times 3$ ユニタリ）

各分子 $i$ に独立に適用される単一qutritゲート：

$$
U_{\mathrm{onsite}}(\Delta t) = \exp(-i\,h_{\mathrm{local}}\,\Delta t) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & e^{-iE_T\Delta t} & 0 \\ 0 & 0 & e^{-iE_S\Delta t} \end{pmatrix}
$$

$h_{\mathrm{local}}$ が対角行列であるため、$U_{\mathrm{onsite}}$ も対角行列であり、各対角要素は対応するエネルギー固有値の位相回転である。これは仮想 Z 回転（VirtRz）に相当する。

Trotterステップ1回あたり、半ステップ $\times 2$ 回で計 $4 \times 2 = 8$ 個の `cu_one` ゲートが適用される。

### 7.2 移動ハミルトニアンゲート: `cu_two`（$9 \times 9$ ユニタリ）

各隣接対 $(i,j)$ に適用される2-qutritゲート：

$$
U_{\mathrm{pair}}^{(i,j)}(\Delta t) = \exp\!\left(-i\,H_{\mathrm{pair}}^{(i,j)}\,\Delta t\right)
$$

$H_{\mathrm{pair}}^{(i,j)}$ は $\{|01\rangle, |10\rangle\}$ 部分空間上でのみ非自明であるため、この部分空間内の $2 \times 2$ ブロック

$$
\begin{pmatrix} 0 & V \\ V & 0 \end{pmatrix} = V \sigma_x
$$

の行列指数関数を計算すれば十分である：

$$
\exp(-iV\sigma_x \Delta t) = \cos(V\Delta t)\,I_2 - i\sin(V\Delta t)\,\sigma_x = \begin{pmatrix} \cos(V\Delta t) & -i\sin(V\Delta t) \\ -i\sin(V\Delta t) & \cos(V\Delta t) \end{pmatrix}
$$

全 $9 \times 9$ 行列は、基底 $\{|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle\}$ の順序で：

$$
U_{\mathrm{pair}}^{(i,j)}(\Delta t) = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & \cos(V\Delta t) & 0 & -i\sin(V\Delta t) & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & -i\sin(V\Delta t) & 0 & \cos(V\Delta t) & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}
$$

物理的には、$|01\rangle \leftrightarrow |10\rangle$ 間のRabi振動（三重項励起子の対間ホッピング）を表す。

Trotterステップ1回あたり、半ステップ $\times 2$ 回で計 $3 \times 2 = 6$ 個の `cu_two` ゲートが適用される。

### 7.3 単一サイトStinespringゲート: `cu_two`（$6 \times 6$ ユニタリ）

単一サイトリンドブラッド演算子 $L_{\mathrm{local}}$（$3 \times 3$）に対し、$6 \times 6$ のStinespringユニタリを構成する。

ゲートは qutrit（$d=3$）と ancilla qubit（$d=2$）の結合系に作用する。基底順序は環境 $\otimes$ 系：
$\{|0_{\mathrm{env}}, 0\rangle, |0_{\mathrm{env}}, 1\rangle, |0_{\mathrm{env}}, 2\rangle, |1_{\mathrm{env}}, 0\rangle, |1_{\mathrm{env}}, 1\rangle, |1_{\mathrm{env}}, 2\rangle\}$

各タイプの局所Stinespringユニタリの厳密な $6 \times 6$ 行列を以下に示す。記法として $\varphi = \sqrt{\gamma\,\Delta t}$、$c = \cos\varphi$、$s = \sin\varphi$ と略記する。

#### 7.3.1 蛍光ゲート ($L = \sqrt{\Gamma_{\mathrm{fl}}} |0\rangle\langle 2|$, $\varphi_{\mathrm{fl}} = \sqrt{\Gamma_{\mathrm{fl}}\,\Delta t}$)

ジャンプ: $|2\rangle \to |0\rangle$（S₁→S₀）。回転が作用する部分空間は $\{|0_{\mathrm{env}}, 2\rangle, |1_{\mathrm{env}}, 0\rangle\}$（添字 2 と 3）。

$$
U_{\mathrm{fl}} = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & \cos\varphi_{\mathrm{fl}} & -i\sin\varphi_{\mathrm{fl}} & 0 & 0 \\
0 & 0 & -i\sin\varphi_{\mathrm{fl}} & \cos\varphi_{\mathrm{fl}} & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}
$$

#### 7.3.2 燐光ゲート ($L = \sqrt{\Gamma_{\mathrm{ph}}} |0\rangle\langle 1|$, $\varphi_{\mathrm{ph}} = \sqrt{\Gamma_{\mathrm{ph}}\,\Delta t}$)

ジャンプ: $|1\rangle \to |0\rangle$（T₁→S₀）。回転が作用する部分空間は $\{|0_{\mathrm{env}}, 1\rangle, |1_{\mathrm{env}}, 0\rangle\}$（添字 1 と 3）。

$$
U_{\mathrm{ph}} = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & \cos\varphi_{\mathrm{ph}} & 0 & -i\sin\varphi_{\mathrm{ph}} & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & -i\sin\varphi_{\mathrm{ph}} & 0 & \cos\varphi_{\mathrm{ph}} & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}
$$

#### 7.3.3 内部転換ゲート ($L = \sqrt{k_{\mathrm{IC}}} |0\rangle\langle 2|$, $\varphi_{\mathrm{IC}} = \sqrt{k_{\mathrm{IC}}\,\Delta t}$)

ジャンプ: $|2\rangle \to |0\rangle$（S₁→S₀）。蛍光ゲートと同一構造で角度のみ異なる。回転が作用する部分空間は $\{|0_{\mathrm{env}}, 2\rangle, |1_{\mathrm{env}}, 0\rangle\}$。

$$
U_{\mathrm{IC}} = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & \cos\varphi_{\mathrm{IC}} & -i\sin\varphi_{\mathrm{IC}} & 0 & 0 \\
0 & 0 & -i\sin\varphi_{\mathrm{IC}} & \cos\varphi_{\mathrm{IC}} & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}
$$

#### 7.3.4 ISC S₁→T₁ ゲート ($L = \sqrt{k_{\mathrm{ISC,ST}}} |1\rangle\langle 2|$, $\varphi_{\mathrm{ISC,ST}} = \sqrt{k_{\mathrm{ISC,ST}}\,\Delta t}$)

ジャンプ: $|2\rangle \to |1\rangle$（S₁→T₁）。回転が作用する部分空間は $\{|0_{\mathrm{env}}, 2\rangle, |1_{\mathrm{env}}, 1\rangle\}$（添字 2 と 4）。

$$
U_{\mathrm{ISC,ST}} = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & \cos\varphi_{\mathrm{ISC,ST}} & 0 & -i\sin\varphi_{\mathrm{ISC,ST}} & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & -i\sin\varphi_{\mathrm{ISC,ST}} & 0 & \cos\varphi_{\mathrm{ISC,ST}} & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}
$$

#### 7.3.5 ISC T₁→S₀ ゲート ($L = \sqrt{k_{\mathrm{ISC,TS}}} |0\rangle\langle 1|$, $\varphi_{\mathrm{ISC,TS}} = \sqrt{k_{\mathrm{ISC,TS}}\,\Delta t}$)

ジャンプ: $|1\rangle \to |0\rangle$（T₁→S₀）。燐光ゲートと同一構造で角度のみ異なる。回転が作用する部分空間は $\{|0_{\mathrm{env}}, 1\rangle, |1_{\mathrm{env}}, 0\rangle\}$。

$$
U_{\mathrm{ISC,TS}} = \begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & \cos\varphi_{\mathrm{ISC,TS}} & 0 & -i\sin\varphi_{\mathrm{ISC,TS}} & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & -i\sin\varphi_{\mathrm{ISC,TS}} & 0 & \cos\varphi_{\mathrm{ISC,TS}} & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}
$$

#### 7.3.6 共通構造の解釈

すべての単一サイトStinespringゲートは、$6 \times 6$ 空間内の特定の $2 \times 2$ 部分空間上でのビームスプリッター（Givens回転）として機能する。回転が作用する2状態は：

- $|0_{\mathrm{env}}, b\rangle$: ancilla $|0\rangle$、系がソース状態 $|b\rangle$ にある
- $|1_{\mathrm{env}}, a\rangle$: ancilla $|1\rangle$、系がターゲット状態 $|a\rangle$ にある

振幅 $\cos\varphi$ で系は遷移せず ancilla は $|0\rangle$ のまま、振幅 $-i\sin\varphi$ で系は $|b\rangle \to |a\rangle$ に遷移し ancilla が $|1\rangle$ に反転する（量子ジャンプの記録）。

単一サイトStinespringゲートは20個（5種類 × 4分子）あり、MQT-Quditsの `cu_two` ゲートとして `[target_qutrit(d=3), ancilla(d=2)]` に適用される。

### 7.4 TTA対Stinespringゲート: `cu_multi`（$18 \times 18$ ユニタリ）

TTA対リンドブラッド演算子 $L_{\mathrm{TTA}}^{\mathrm{local}}$（$9 \times 9$）に対し、$18 \times 18$ のStinespringユニタリを構成する。

ゲートは2つのqutrit（$d=3$）とancilla qubit（$d=2$）の結合系に作用する。基底順序は $|e_{\mathrm{env}}\rangle \otimes |s_i, s_j\rangle$ で、添字は $e \cdot 9 + s_i \cdot 3 + s_j$（$e \in \{0,1\}$、$s_i, s_j \in \{0,1,2\}$）。

#### 7.4.1 TTAチャネル1ゲート ($L = \sqrt{\gamma_{\mathrm{TTA}}/2}\;|20\rangle\langle 11|$, $\varphi_{\mathrm{TTA}} = \sqrt{\gamma_{\mathrm{TTA}}\,\Delta t / 2}$)

ジャンプ: $|1,1\rangle \to |2,0\rangle$（$\mathrm{T}_1\mathrm{T}_1 \to \mathrm{S}_1\mathrm{S}_0$）。

$L$ は添字 $(6, 4)$ のみ非零のランク1行列（$|20\rangle$ の添字 $= 6$、$|11\rangle$ の添字 $= 4$）。

回転が作用する $18 \times 18$ 空間内の2状態：
- $|0_{\mathrm{env}}, 11\rangle$: 添字 $0 \cdot 9 + 4 = 4$
- $|1_{\mathrm{env}}, 20\rangle$: 添字 $1 \cdot 9 + 6 = 15$

$18 \times 18$ 行列は恒等行列に対して添字 $(4, 4)$, $(4, 15)$, $(15, 4)$, $(15, 15)$ のみ変更される：

$$
U_{\mathrm{TTA},1}[4,4] = \cos\varphi_{\mathrm{TTA}}, \quad U_{\mathrm{TTA},1}[4,15] = -i\sin\varphi_{\mathrm{TTA}}
$$

$$
U_{\mathrm{TTA},1}[15,4] = -i\sin\varphi_{\mathrm{TTA}}, \quad U_{\mathrm{TTA},1}[15,15] = \cos\varphi_{\mathrm{TTA}}
$$

$$
U_{\mathrm{TTA},1}[k,k] = 1 \quad (k \neq 4, 15), \quad U_{\mathrm{TTA},1}[k,l] = 0 \quad (k \neq l,\; (k,l) \neq (4,15),\;(15,4))
$$

#### 7.4.2 TTAチャネル2ゲート ($L = \sqrt{\gamma_{\mathrm{TTA}}/2}\;|02\rangle\langle 11|$, $\varphi_{\mathrm{TTA}} = \sqrt{\gamma_{\mathrm{TTA}}\,\Delta t / 2}$)

ジャンプ: $|1,1\rangle \to |0,2\rangle$（$\mathrm{T}_1\mathrm{T}_1 \to \mathrm{S}_0\mathrm{S}_1$）。

$L$ は添字 $(2, 4)$ のみ非零（$|02\rangle$ の添字 $= 2$）。

回転が作用する2状態：
- $|0_{\mathrm{env}}, 11\rangle$: 添字 $4$
- $|1_{\mathrm{env}}, 02\rangle$: 添字 $1 \cdot 9 + 2 = 11$

$$
U_{\mathrm{TTA},2}[4,4] = \cos\varphi_{\mathrm{TTA}}, \quad U_{\mathrm{TTA},2}[4,11] = -i\sin\varphi_{\mathrm{TTA}}
$$

$$
U_{\mathrm{TTA},2}[11,4] = -i\sin\varphi_{\mathrm{TTA}}, \quad U_{\mathrm{TTA},2}[11,11] = \cos\varphi_{\mathrm{TTA}}
$$

他の全要素は恒等行列と同じ。

TTAゲートは6個（3対 × 2チャネル）あり、MQT-Quditsの `cu_multi` ゲートとして `[qutrit_i(d=3), qutrit_j(d=3), ancilla(d=2)]` に適用される。

### 7.5 全系への埋め込み

局所ゲートはテンソル積構造を用いて全 $81$ 次元系に埋め込まれる。

**単一サイトゲート** ($K$ が $3 \times 3$、サイト $k$ に作用):

$$
K_{\mathrm{full}} = I_{d}^{\otimes k} \otimes K \otimes I_{d}^{\otimes (N-1-k)}
$$

**対ゲート** ($K$ が $9 \times 9$、サイト $(i,j)$ に作用):

対のKraus演算子 $K$ を全系に埋め込む際、対 $(i,j)$ 以外のqutritの状態が行と列で一致する場合のみ非零要素を持つ：

$$
K_{\mathrm{full}}[\mathbf{s}, \mathbf{s}'] = \begin{cases} K[s_i d + s_j,\; s'_i d + s'_j] & \text{if } s_k = s'_k \;\forall k \neq i,j \\ 0 & \text{otherwise} \end{cases}
$$

### 7.6 ゲート数のまとめ

| ゲートタイプ | MQT-Qudits API | 次元 | 個数 |
|-------------|----------------|------|------|
| オンサイト位相 | `cu_one` | $3 \times 3$ | $4 \times 2 = 8$ |
| 対移動ユニタリ | `cu_two` | $9 \times 9$ | $3 \times 2 = 6$ |
| 単一サイトStinespring | `cu_two` | $6 \times 6$ | $20$ |
| TTA対Stinespring | `cu_multi` | $18 \times 18$ | $6$ |
| **合計** | — | — | **40** |

（`QuditGKSLSimulator` のコード内コメントではゲート数を33と見積もっているが、これはハミルトニアン半ステップの重複を数えない簡略計算 $4 + 3 + 26 = 33$ による。回路シミュレータの実カウントは上記40個。）

---

## 8. 全系Stinespringと局所Stinespringの等価性

`QuditGKSLSimulator`（セル6で使用）は全系 $81 \times 81$ のリンドブラッド演算子から $162 \times 162$ のStinespringユニタリを構成する。一方、`QuditGKSLCircuitSimulator`（回路可視化で使用）は局所演算子（$3 \times 3$ または $9 \times 9$）から $6 \times 6$ または $18 \times 18$ の局所ユニタリを構成する。

これらが数学的に等価であることは以下から示される。単一サイトの場合、全系リンドブラッド演算子は

$$
L_{\mathrm{full}} = I^{\otimes k} \otimes L_{\mathrm{local}} \otimes I^{\otimes (N-1-k)}
$$

であるから、生成子は

$$
G_{\mathrm{full}} = \begin{pmatrix} 0 & L_{\mathrm{full}}^\dagger \\ L_{\mathrm{full}} & 0 \end{pmatrix}
$$

$G_{\mathrm{full}}$ の作用はサイト $k$ と環境の結合部分空間上でのみ非自明であり、他のサイトには恒等演算子として作用する。したがって

$$
U_{\mathrm{full}} = \exp(-i\theta G_{\mathrm{full}}) = I_{\mathrm{others}} \otimes U_{\mathrm{local}}
$$

ここで $U_{\mathrm{local}} = \exp(-i\theta G_{\mathrm{local}})$ は局所Stinespringユニタリ。

全系レベルでのKraus演算子は

$$
K_{0,\mathrm{full}} = I_{\mathrm{others}} \otimes K_{0,\mathrm{local}}, \quad K_{1,\mathrm{full}} = I_{\mathrm{others}} \otimes K_{1,\mathrm{local}}
$$

であり、チャネルの結果は

$$
\mathcal{E}(\rho) = K_{0,\mathrm{full}}\,\rho\,K_{0,\mathrm{full}}^\dagger + K_{1,\mathrm{full}}\,\rho\,K_{1,\mathrm{full}}^\dagger
$$

となって、局所ゲートによる計算と全系ユニタリによる計算は同一の量子チャネルを実現する。TTA対の場合も同様に成立する。

---

## 9. 収束次数解析

### 9.1 各近似の誤差次数

Trotterステップの3つの近似層と、それぞれの誤差次数：

| 近似 | 1ステップ誤差 | 大域誤差 |
|------|-------------|---------|
| Strang分割（$\mathcal{L}_H$ vs $\mathcal{L}_D$） | $O(\Delta t^3)$ | $O(\Delta t^2)$ |
| 回文Lie-Trotter積（個別 $\mathcal{D}_\alpha$ の合成） | $O(\Delta t^3)$ | $O(\Delta t^2)$ |
| Stinespring近似（各 $\mathcal{E}_\alpha$ の精度） | $O(\Delta t^2)$ | $O(\Delta t)$ |

### 9.2 詳細

1. **Strang分割**: ハミルトニアン部分 $\mathcal{L}_H$ と散逸部分 $\mathcal{L}_D$ の分割は対称（2次）Trotter分解であり、1ステップの誤差は $O(\Delta t^3)$。$T/\Delta t$ ステップの累積で大域誤差は $O(\Delta t^2)$。

2. **回文Lie-Trotter積**: 個別の散逸チャネル $\mathcal{E}_\alpha$ を回文順序（順方向 $\to$ 逆方向）で合成することにより、Lie-Trotter積の1次交換子誤差が相殺される。これも2次精度 $O(\Delta t^3)$/ステップ。

3. **Stinespring近似**: 各Stinespringチャネル $\mathcal{E}_\alpha(\Delta t)$ は正確な散逸チャネル $e^{\mathcal{L}_{D_\alpha}\Delta t}$ の1次近似であり、$\mathcal{E}_\alpha(\Delta t)(\rho) = e^{\mathcal{L}_{D_\alpha}\Delta t}(\rho) + O(\Delta t^2)$。この誤差は1ステップあたり $O(\Delta t^2)$ であり、$T/\Delta t$ ステップの累積で大域誤差は $O(\Delta t)$。

### 9.3 結論

ボトルネックはStinespring近似（3.）であり、**トレース距離での実効的な大域収束次数は $O(\Delta t)$（1次）** である。

---

## 10. 初期状態

シミュレーションで使用されるデフォルトの初期状態は `"edge_triplet"` であり、分子0と分子3（両端）がT₁状態、分子1と分子2（中央）がS₀状態にある：

$$
|\psi_0\rangle = |1, 0, 0, 1\rangle = |T_1, S_0, S_0, T_1\rangle
$$

計算基底における添字は

$$
\mathrm{index} = 1 \cdot 3^3 + 0 \cdot 3^2 + 0 \cdot 3 + 1 = 27 + 0 + 0 + 1 = 28
$$

初期密度行列は純粋状態：

$$
\rho_0 = |\psi_0\rangle\langle\psi_0|
$$

これは81次元ベクトルの第28成分のみが1、他が0であるベクトルの外積として構成される。

---

## 11. 観測量

シミュレーション中に以下の物理量が各時間ステップで計算される。

### 11.1 集団占有数（Population）

密度行列の対角要素から各電子状態の占有数を計算する：

$$
N_{\mathrm{S}_0}(t) = \sum_{\mathbf{s}} \rho_{\mathbf{s},\mathbf{s}}(t) \cdot \#\{i : s_i = 0\}
$$

$$
N_{\mathrm{T}_1}(t) = \sum_{\mathbf{s}} \rho_{\mathbf{s},\mathbf{s}}(t) \cdot \#\{i : s_i = 1\}
$$

$$
N_{\mathrm{S}_1}(t) = \sum_{\mathbf{s}} \rho_{\mathbf{s},\mathbf{s}}(t) \cdot \#\{i : s_i = 2\}
$$

ここで和は全計算基底 $\mathbf{s} = (s_0, s_1, s_2, s_3)$ にわたり、$\#\{i : s_i = k\}$ は状態 $\mathbf{s}$ で状態 $k$ にある分子の数である。

### 11.2 フォン・ノイマンエントロピー

$$
S(\rho) = -\mathrm{Tr}[\rho \ln \rho] = -\sum_i \lambda_i \ln \lambda_i
$$

ここで $\lambda_i$ は $\rho$ の固有値（$\lambda_i > 0$ のもののみ）。

### 11.3 純度（Purity）

$$
P(\rho) = \mathrm{Tr}[\rho^2]
$$

$P = 1$ は純粋状態、$P < 1$ は混合状態を示す。

### 11.4 トレース保存

$$
\mathrm{Tr}[\rho(t)] = 1
$$

の保存をトレース偏差 $|\mathrm{Tr}[\rho(t)] - 1|$ により各ステップで検証する。

---

## 12. 量子資源の見積もり

セル5のマークダウンに記載されている量子資源：

- **系 qutrit 数**: $N = 4$（各分子に1 qutrit）
- **ancilla qubit 数**: 26（各リンドブラッドチャネルに1 qubit）
- **合計量子レジスタ**: $4 + 26 = 30$

各 qutrit は $\lceil\log_2 3\rceil = 2$ qubit のハードウェア量子ビットに相当するため、ハードウェア上の等価qubit数は $4 \times 2 + 26 = 34$ qubit となる。

ただし、Stinespring拡張の ancilla qubit は各チャネルの適用後に部分トレースで消去され、実際には逐次再利用可能であるため、同時に必要な ancilla は1個のみである。

---

## 参考: コードファイルの対応

| ファイル | 役割 |
|---------|------|
| `tutorials/qudit_gksl_simulator.py` | セル6で使用される `QuditGKSLSimulator` クラス |
| `tutorials/stinespring_utils.py` | Stinespringユニタリの構成と密度行列への適用 |
| `tutorials/gksl_math_utils.py` | ハミルトニアン・リンドブラッド演算子の構成、物理量の計算 |
| `tutorials/gksl_physical_parameters.py` | 物理パラメータクラス `GKSLPhysicalParameters` |
| `tutorials/qudit_gksl_circuit_simulator.py` | 量子ゲート分解による回路シミュレータ `QuditGKSLCircuitSimulator` |
