# シナリオ5: Qudit GKSL（ボソン無し）— 理論まとめ

---

## スライド1: GKSL-Lindblad マスター方程式とハミルトニアン・散逸演算子

### 1.1 GKSL マスター方程式

系の密度行列 $\rho$ の時間発展は GKSL (Gorini–Kossakowski–Sudarshan–Lindblad) マスター方程式に従う（$\hbar = 1$）：

$$
\frac{d\rho}{dt} = -i[H,\,\rho] + \sum_{\alpha} \left( L_\alpha \rho L_\alpha^\dagger - \frac{1}{2}\{L_\alpha^\dagger L_\alpha,\,\rho\} \right)
$$

ここで $H$ は系の全ハミルトニアン、$L_\alpha$ はリンドブラッド演算子であり、各 $L_\alpha$ には散逸率 $\sqrt{\gamma_\alpha}$ が既に含まれている。

### 1.2 系の構成

- 分子数 $N = 4$、各分子のヒルベルト空間次元 $d = 3$（qutrit）
- 各分子の局所状態: $|0\rangle = S_0$（基底一重項）、$|1\rangle = T_1$（三重項）、$|2\rangle = S_1$（励起一重項）
- 全系のヒルベルト空間次元: $d^N = 3^4 = 81$
- 禁止状態は存在しない（qubit エンコーディングの $|11\rangle$ 禁止問題が無い）

### 1.3 オンサイトハミルトニアン $H_0$

$$
H_0 = \sum_{i=0}^{N-1} \left( E_T |1\rangle_i\langle 1| + E_S |2\rangle_i\langle 2| \right) \otimes \bigotimes_{k \neq i} I_3^{(k)}
$$

各分子の局所ハミルトニアンは対角行列:

$$
h_{\mathrm{local}} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix}
$$

デフォルトパラメータ: $E_T = 1.5\;\mathrm{eV}$, $E_S = 3.0\;\mathrm{eV}$（TTA-UC の条件 $2E_T \approx E_S$）。

### 1.4 移動ハミルトニアン $H_{\mathrm{transfer}}$

最近接分子ペア $\langle i, j \rangle$（$(0,1), (1,2), (2,3)$ の 3 ペア）に対して:

$$
H_{\mathrm{transfer}} = \sum_{\langle i,j \rangle} V \left( |0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| + |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1| \right) \otimes \bigotimes_{k \neq i,j} I_3^{(k)}
$$

ここで $V = 0.1\;\mathrm{eV}$ はトリプレット励起子の移動積分である。ペア $(i,j)$ の局所 $9 \times 9$ ハミルトニアンは:

$$
H_{\mathrm{pair}}^{(i,j)} = V \left( |01\rangle\langle 10| + |10\rangle\langle 01| \right)
$$

全ハミルトニアン: $H = H_0 + H_{\mathrm{transfer}}$

### 1.5 リンドブラッド演算子（全 26 チャネル）

| チャネル | 演算子 $L_\alpha$ | 散逸率 $\gamma_\alpha$ | 個数 |
|---|---|---|---|
| TTA (チャネル1) | $\sqrt{\gamma_{\mathrm{TTA}}/2}\; |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|$ | $\gamma_{\mathrm{TTA}}/2 = 0.025$ | 3 |
| TTA (チャネル2) | $\sqrt{\gamma_{\mathrm{TTA}}/2}\; |0\rangle_i\langle 1| \otimes |2\rangle_j\langle 1|$ | $\gamma_{\mathrm{TTA}}/2 = 0.025$ | 3 |
| 蛍光 (Fluorescence) | $\sqrt{\Gamma_{\mathrm{fl}}}\; |0\rangle_i\langle 2|$ | $\Gamma_{\mathrm{fl}} = 0.01$ | 4 |
| 燐光 (Phosphorescence) | $\sqrt{\Gamma_{\mathrm{ph}}}\; |0\rangle_i\langle 1|$ | $\Gamma_{\mathrm{ph}} = 10^{-6}$ | 4 |
| 内部変換 (IC) | $\sqrt{k_{\mathrm{IC}}}\; |0\rangle_i\langle 2|$ | $k_{\mathrm{IC}} = 0.005$ | 4 |
| ISC $S_1 \to T_1$ | $\sqrt{k_{\mathrm{ISC,ST}}}\; |1\rangle_i\langle 2|$ | $k_{\mathrm{ISC,ST}} = 0.003$ | 4 |
| ISC $T_1 \to S_0$ | $\sqrt{k_{\mathrm{ISC,TS}}}\; |0\rangle_i\langle 1|$ | $k_{\mathrm{ISC,TS}} = 10^{-5}$ | 4 |

合計: $2 \times 3 + 5 \times 4 = 26$ チャネル。単体演算子の場合、全系演算子は $L_\alpha = L_{\mathrm{local}} \otimes \bigotimes_{k \neq i} I_3^{(k)}$ で構成。TTA ペア演算子は 2 分子のテンソル積で構成される。

---

## スライド2: Stinespring 拡張・対称 Trotter 分解と量子ゲート表現

### 2.1 Strang (対称) Trotter 分解

1 タイムステップ $\Delta t$ での密度行列の時間発展を、ハミルトニアン部分 $\mathcal{L}_H$ と散逸部分 $\{\mathcal{E}_\alpha\}$ に分解する:

$$
e^{\mathcal{L}\Delta t}\rho \approx e^{\mathcal{L}_H \Delta t/2} \left[\prod_{\alpha=1}^{n} \mathcal{E}_\alpha(\Delta t/2)\right] \left[\prod_{\alpha=n}^{1} \mathcal{E}_\alpha(\Delta t/2)\right] e^{\mathcal{L}_H \Delta t/2} \rho
$$

- ハミルトニアン半ステップ: $\rho \mapsto U_{H}^{(\Delta t/2)} \rho \, U_{H}^{(\Delta t/2)\dagger}$, ただし $U_{H}^{(\Delta t/2)} = e^{-iH\Delta t/2}$
- リンドブラッドチャネルは**パリンドロミック（回文的）順序**で適用：前進順 $\alpha = 1, \ldots, n$ の後に逆順 $\alpha = n, \ldots, 1$ を適用することで Lie-Trotter 積の交換子誤差を消去
- Strang 分割自体は 2 次精度 $O(\Delta t^3)$/ステップだが、後述の Stinespring 近似が 1 次であるため、全体の収束次数はトレース距離で **$O(\Delta t)$（1 次）**

### 2.2 Stinespring 拡張（各リンドブラッドチャネル）

各リンドブラッド演算子 $L_\alpha$（$\sqrt{\gamma_\alpha}$ 込み）に対し、系 + 1 補助量子ビットの拡大空間上のユニタリ $U_{\mathrm{Stine}}$ を構成する。生成子はブロック行列:

$$
G = \begin{pmatrix} 0 & L_\alpha^\dagger \\ L_\alpha & 0 \end{pmatrix}
$$

ここで $G$ は $2d_{\mathrm{sys}} \times 2d_{\mathrm{sys}}$ のエルミート行列。Stinespring ユニタリは:

$$
U_{\mathrm{Stine}} = \exp\!\left(-i\sqrt{\Delta t}\; G\right)
$$

補助量子ビットを $|0\rangle_{\mathrm{anc}}$ に初期化し、テンソル積状態 $|0\rangle\langle 0|_{\mathrm{anc}} \otimes \rho$ にユニタリを適用した後、補助系を部分トレースで除去する:

$$
\mathcal{E}_\alpha(\Delta t)(\rho) = \mathrm{Tr}_{\mathrm{anc}}\!\left[ U_{\mathrm{Stine}} \left( |0\rangle\langle 0|_{\mathrm{anc}} \otimes \rho \right) U_{\mathrm{Stine}}^\dagger \right]
$$

$U_{\mathrm{Stine}}$ を $d_{\mathrm{sys}} \times d_{\mathrm{sys}}$ のブロックに分割すると:

$$
U_{\mathrm{Stine}} = \begin{pmatrix} A & B \\ C & D \end{pmatrix}
$$

このとき部分トレースの結果は:

$$
\mathcal{E}_\alpha(\Delta t)(\rho) = A\rho A^\dagger + C\rho C^\dagger
$$

主要次数で GKSL 散逸子を再現する: $\mathcal{E}_\alpha(\Delta t)(\rho) = \rho + \Delta t\left(L_\alpha \rho L_\alpha^\dagger - \frac{1}{2}\{L_\alpha^\dagger L_\alpha, \rho\}\right) + O(\Delta t^2)$

### 2.3 量子ゲート表現（MQT-Qudits 回路）

1 Trotter ステップの回路は以下のゲート列で構成される:

**ハミルトニアン半ステップ（前半・後半で各 1 回、計 2 回適用）:**

1. **cu_one ゲート** ×4（各 qutrit $i = 0, 1, 2, 3$ に適用）:
$$
U_{\mathrm{onsite}}(\Delta t/2) = \exp\!\left(-i \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix} \frac{\Delta t}{2}\right) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & e^{-iE_T \Delta t/2} & 0 \\ 0 & 0 & e^{-iE_S \Delta t/2} \end{pmatrix}
$$

2. **cu_two ゲート** ×3（ペア $(0,1), (1,2), (2,3)$ に適用、各 $9 \times 9$）:
$$
U_{\mathrm{pair}}^{(i,j)}(\Delta t/2) = \exp\!\left(-i V \left(|01\rangle\langle 10| + |10\rangle\langle 01|\right) \frac{\Delta t}{2}\right)
$$

この $9 \times 9$ 行列は、基底 $\{|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle\}$ において、$|01\rangle$–$|10\rangle$ 部分空間でのみ非自明:

$$
U_{\mathrm{pair}}^{(i,j)} = I_9 + (\cos(V\Delta t/2) - 1)(|01\rangle\langle 01| + |10\rangle\langle 10|) - i\sin(V\Delta t/2)(|01\rangle\langle 10| + |10\rangle\langle 01|)
$$

すなわち $|01\rangle$–$|10\rangle$ 2次元部分空間では回転:

$$
\begin{pmatrix} \langle 01|U|01\rangle & \langle 01|U|10\rangle \\ \langle 10|U|01\rangle & \langle 10|U|10\rangle \end{pmatrix} = \begin{pmatrix} \cos(V\Delta t/2) & -i\sin(V\Delta t/2) \\ -i\sin(V\Delta t/2) & \cos(V\Delta t/2) \end{pmatrix}
$$

他の基底成分はすべて 1（恒等）。

**リンドブラッドチャネル（パリンドロミック順、各半ステップ $\Delta t/2$）:**

3. **単体チャネル**: cu_two ゲート ×20（4 分子 × 5 種類）

   局所リンドブラッド演算子 $L_{\mathrm{local}}$（$3 \times 3$）に対し、Stinespring ユニタリは $6 \times 6$:

$$
U_{\mathrm{Stine}}^{\mathrm{single}} = \exp\!\left(-i\sqrt{\Delta t/2}\begin{pmatrix} 0_{3\times 3} & L_{\mathrm{local}}^\dagger \\ L_{\mathrm{local}} & 0_{3\times 3} \end{pmatrix}\right) \in \mathbb{C}^{6 \times 6}
$$

   これは qutrit $i$ と 1 補助量子ビット（$d=2$）に作用する cu_two ゲートとして実装。

4. **TTA ペアチャネル**: cu_multi ゲート ×6（3 ペア × 2 チャネル）

   ペアリンドブラッド演算子 $L_{\mathrm{pair}}$（$9 \times 9$）に対し、Stinespring ユニタリは $18 \times 18$:

$$
U_{\mathrm{Stine}}^{\mathrm{pair}} = \exp\!\left(-i\sqrt{\Delta t/2}\begin{pmatrix} 0_{9\times 9} & L_{\mathrm{pair}}^\dagger \\ L_{\mathrm{pair}} & 0_{9\times 9} \end{pmatrix}\right) \in \mathbb{C}^{18 \times 18}
$$

   これは qutrit ペア $(i, j)$ と 1 補助量子ビットに作用する cu_multi ゲートとして実装。

### 2.4 1 Trotter ステップの量子資源

| ゲート種類 | サイズ | 個数（半ステップ）| 個数（合計） |
|---|---|---|---|
| cu_one（オンサイト位相） | $3 \times 3$ | 4 | 8 |
| cu_two（移動ハミルトニアン） | $9 \times 9$ | 3 | 6 |
| cu_two（単体 Stinespring） | $6 \times 6$ | 20 | 20 |
| cu_multi（TTA Stinespring） | $18 \times 18$ | 6 | 6 |
| **合計** | | | **40** |

- 系 qutrit 数: 4
- 補助量子ビット数: 26（各 Stinespring チャネルに 1 ビット）
- 全量子ビット相当: 30（= $4 \times \lceil\log_2 3\rceil + 26$）
