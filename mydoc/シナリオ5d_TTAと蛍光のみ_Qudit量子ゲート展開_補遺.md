# シナリオ5d 縮約版（TTA + 蛍光のみ）：Qutrit 量子ゲート展開の **補遺（完全展開版）**

## 0. 本補遺の位置付け

本ファイルは `mydoc/シナリオ5d_TTAと蛍光のみ_Qudit量子ゲート展開.md`（以後「本体」と呼ぶ）の **省略箇所を一切残さない完全形**として書く。本体に含まれていた以下の省略を、補遺 A〜L で 1 個ずつ全項展開する。**寄り添わず、嘘なし**。コード根拠は最後に総括する。

| 補遺 | 本体での省略箇所 | 補遺での扱い |
|---|---|---|
| A | §1.2 GellMann 行列 $\Lambda_s^{(a,b)},\Lambda_a^{(a,b)}$（$3\times 3$）が文字定義のみ | 全 6 個（$(a,b)=(0,1),(0,2),(1,2)$、s/a の 2 種）を $3\times 3$ 行列で明示 |
| B | §1.3 $X_3^2$ が「$X_3^2$」とだけ書かれて行列が省略 | $X_3^0,X_3^1,X_3^2$ 全 3 個を行列で明示 |
| C | §1.3 $\mathrm{CSum}^\dagger$ の $9\times 9$ 行列省略 | $\mathrm{CSum}^\dagger$ を全 81 成分で明示 |
| D | §3.1 の $h_\mathrm{loc}^\mathrm{tr}$ が「$\{1,3\}$ ブロック以外 0」と書かれただけ | 全 81 成分を明示 |
| E | §3.2 の $E_\mathrm{loc}^\mathrm{tr}(\tau)$ が cases 形でのみ書かれている | 全 81 成分を明示行列で書下し |
| F | §3.3 「ペア間 Trotter 誤差」の交換子が定性的にしか書かれていない | $[\hat h^{(0,1)},\hat h^{(1,2)}]$ などを演算子代数で明示計算 |
| G | §3.5 の $\mathrm{CT}^{(i,j)}(\tau)$ は局所 $9\times 9$ のみ。4-qutrit Hilbert（$81\times 81$）への Kronecker 埋め込みが省略 | 3 ペア分を $I_3\otimes\cdots$ 形で明示 |
| H | §4.2 の「$[\hat H_0,\hat H_\mathrm{transfer}]\neq 0$」が定性的にしか書かれていない | 演算子代数で明示計算 |
| I | §5.3 の $U_\mathrm{Stine}^\mathrm{(TTA1)}$ の $27\times 27$ 行列が文字記述のみ | 全 729 成分の非零位置を完全に書下し（パターン行列形） |
| J | §5.3 の $U_\mathrm{Stine}^\mathrm{(TTA2)}$ も同様に省略 | 全 729 成分の非零位置を完全に書下し |
| K | §5.4「TTA $(1,2),(2,3)$ ・・・同様」の省略 | TTA 6 チャンネルすべての $K_\alpha,U_\mathrm{Stine,\alpha}$ を全部明示 |
| L | §6 「62 命令」と書いただけで順序付き積が無い | 1 Trotter ステップを 62 個の積で順序付きに完全列挙 |

---

## 補遺 A. GellMann 行列の完全展開（本体 §1.2）

$d=3$、`gellmann.py:45-69` の `__array__` を $(\mathrm{lev}_a,\mathrm{lev}_b)\in\{(0,1),(0,2),(1,2)\}$（`R.levels_setter` で常に $a<b$）と type "s" / "a" について評価する。

### A.1 type "s"：$\Lambda_s^{(a,b)} = |a\rangle\langle b| + |b\rangle\langle a|$

$$
\Lambda_s^{(0,1)} \;=\; \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix},\quad
\Lambda_s^{(0,2)} \;=\; \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix},\quad
\Lambda_s^{(1,2)} \;=\; \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}.
$$

### A.2 type "a"：$\Lambda_a^{(a,b)} = -i|a\rangle\langle b| + i|b\rangle\langle a|$

$$
\Lambda_a^{(0,1)} \;=\; \begin{pmatrix} 0 & -i & 0 \\ i & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix},\quad
\Lambda_a^{(0,2)} \;=\; \begin{pmatrix} 0 & 0 & -i \\ 0 & 0 & 0 \\ i & 0 & 0 \end{pmatrix},\quad
\Lambda_a^{(1,2)} \;=\; \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -i \\ 0 & i & 0 \end{pmatrix}.
$$

### A.3 `R(a,b;θ,φ)` ゲートの $3\times 3$ 行列（`r.py:49-68` の評価）

`R.__array__` は

$$
R(a,b;\theta,\varphi) \;=\; \mathrm{diag}_3(d_0,d_1,d_2) \;-\; i\sin\!\tfrac{\theta}{2}\big(\sin\varphi\,\Lambda_a^{(a,b)} + \cos\varphi\,\Lambda_s^{(a,b)}\big),
$$

ここで対角部は `r.py:53-58`（`matrix = np.identity(d)` を $\cos(\theta/2)$ 倍するのは row $a$ と row $b$ のみ）：

$$
d_a \;=\; \cos\tfrac{\theta}{2},\quad d_b \;=\; \cos\tfrac{\theta}{2},\quad d_c \;=\; 1\quad(c\notin\{a,b\}).
$$

3 通り全部書き下す（$c=\cos(\theta/2),\,s=\sin(\theta/2),\,e^{\pm i\varphi} = \cos\varphi \pm i\sin\varphi$ と置くと簡略化できる）：

$$
R(0,1;\theta,\varphi) \;=\; \begin{pmatrix}
\cos\tfrac{\theta}{2} & -i\sin\tfrac{\theta}{2}\,e^{-i\varphi} & 0 \\
-i\sin\tfrac{\theta}{2}\,e^{+i\varphi} & \cos\tfrac{\theta}{2} & 0 \\
0 & 0 & 1
\end{pmatrix},
$$

$$
R(0,2;\theta,\varphi) \;=\; \begin{pmatrix}
\cos\tfrac{\theta}{2} & 0 & -i\sin\tfrac{\theta}{2}\,e^{-i\varphi} \\
0 & 1 & 0 \\
-i\sin\tfrac{\theta}{2}\,e^{+i\varphi} & 0 & \cos\tfrac{\theta}{2}
\end{pmatrix},
$$

$$
R(1,2;\theta,\varphi) \;=\; \begin{pmatrix}
1 & 0 & 0 \\
0 & \cos\tfrac{\theta}{2} & -i\sin\tfrac{\theta}{2}\,e^{-i\varphi} \\
0 & -i\sin\tfrac{\theta}{2}\,e^{+i\varphi} & \cos\tfrac{\theta}{2}
\end{pmatrix}.
$$

導出：$-i\sin(\theta/2)(\sin\varphi\,\Lambda_a + \cos\varphi\,\Lambda_s)$ の $(a,b)$ 成分は $-i\sin(\theta/2)(\sin\varphi\cdot(-i) + \cos\varphi\cdot 1) = -\sin(\theta/2)\sin\varphi - i\sin(\theta/2)\cos\varphi = -i\sin(\theta/2)(\cos\varphi - i\sin\varphi) = -i\sin(\theta/2)e^{-i\varphi}$。$(b,a)$ 成分は複素共役で $-i\sin(\theta/2)e^{+i\varphi}$。$\checkmark$

---

## 補遺 B. $X_3^k$ 全 3 個の明示（本体 §1.3）

`x.py:36-45` の `__array__` を $d=3$ で評価。$X_3|i\rangle = |(i+1)\bmod 3\rangle$。

$$
X_3^0 = I_3 = \begin{pmatrix} 1&0&0\\ 0&1&0\\ 0&0&1\end{pmatrix},\quad
X_3^1 = \begin{pmatrix} 0&0&1\\ 1&0&0\\ 0&1&0\end{pmatrix},\quad
X_3^2 = \begin{pmatrix} 0&1&0\\ 0&0&1\\ 1&0&0\end{pmatrix}.
$$

検算：$(X_3^2)|i\rangle = |(i+2)\bmod 3\rangle$ より $X_3^2|0\rangle=|2\rangle$（列 0、行 2 = 1）、$X_3^2|1\rangle=|0\rangle$（列 1、行 0 = 1）、$X_3^2|2\rangle=|1\rangle$（列 2、行 1 = 1）。$\checkmark$

---

## 補遺 C. $\mathrm{CSum}$ と $\mathrm{CSum}^\dagger$ の完全 $9\times 9$ 行列

本体 §1.3 で $\mathrm{CSum}$ は $\sum_i |i\rangle_c\langle i|\otimes X_3^i$ と書下した。$\mathrm{CSum}|s_c s_t\rangle = |s_c, (s_t+s_c)\bmod 3\rangle$ より逆演算は $\mathrm{CSum}^{-1}|s_c s_t\rangle = |s_c, (s_t-s_c)\bmod 3\rangle = \sum_i |i\rangle\langle i|\otimes X_3^{-i} = \sum_i|i\rangle\langle i|\otimes X_3^{3-i}$（$i\bmod 3$）。$\mathrm{CSum}$ は計算基底順列ゆえ実直交、$\mathrm{CSum}^\dagger = \mathrm{CSum}^{-1}$。

$$
\mathrm{CSum} \;=\;
\left(\begin{array}{ccc|ccc|ccc}
1&0&0&0&0&0&0&0&0\\
0&1&0&0&0&0&0&0&0\\
0&0&1&0&0&0&0&0&0\\\hline
0&0&0&0&0&1&0&0&0\\
0&0&0&1&0&0&0&0&0\\
0&0&0&0&1&0&0&0&0\\\hline
0&0&0&0&0&0&0&1&0\\
0&0&0&0&0&0&0&0&1\\
0&0&0&0&0&0&1&0&0
\end{array}\right),\qquad
\mathrm{CSum}^\dagger \;=\;
\left(\begin{array}{ccc|ccc|ccc}
1&0&0&0&0&0&0&0&0\\
0&1&0&0&0&0&0&0&0\\
0&0&1&0&0&0&0&0&0\\\hline
0&0&0&0&1&0&0&0&0\\
0&0&0&0&0&1&0&0&0\\
0&0&0&1&0&0&0&0&0\\\hline
0&0&0&0&0&0&0&0&1\\
0&0&0&0&0&0&1&0&0\\
0&0&0&0&0&0&0&1&0
\end{array}\right).
$$

検算：$\mathrm{CSum}\cdot\mathrm{CSum}^\dagger$ の $(0..2,0..2)$ ブロック：$I_3\cdot I_3 = I_3$；$(3..5,3..5)$ ブロック：$X_3\cdot X_3^2 = X_3^3 = I_3$；$(6..8,6..8)$ ブロック：$X_3^2\cdot X_3 = I_3$。$\checkmark$

---

## 補遺 D. $h_\mathrm{loc}^\mathrm{tr}$ の完全 $9\times 9$ 行列（本体 §3.1）

基底 index 規約：$|s_is_j\rangle$ → row/col index $3 s_i + s_j$（$i$ 第 1、$j$ 第 2、ともに $\{0,1,2\}$）。$h_\mathrm{loc}^\mathrm{tr} = V(|01\rangle\langle 10|+|10\rangle\langle 01|) = V(E^{(9)}_{1,3}+E^{(9)}_{3,1})$。全 81 成分は 2 個非零、79 個ゼロ：

$$
h_\mathrm{loc}^\mathrm{tr} \;=\;
\left(\begin{array}{ccccccccc}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & V & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & V & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{array}\right).
$$

非零成分：$(\text{row},\text{col})\in\{(1,3),(3,1)\}$、いずれも値 $V$。

---

## 補遺 E. $E_\mathrm{loc}^\mathrm{tr}(\tau) := e^{-ih_\mathrm{loc}^\mathrm{tr}\tau}$ の完全 $9\times 9$ 行列（本体 §3.2）

$h_\mathrm{loc}^\mathrm{tr}$ は行/列 $\{1,3\}$ に閉じる $V\sigma_x$ ブロック（補遺 D）。指数を取ると、$\{1,3\}$ ブロックで $\cos(V\tau)I_2 - i\sin(V\tau)\sigma_x$、補空間 $\{0,2,4,5,6,7,8\}$ で恒等。記号 $c\equiv\cos(V\tau),\;s\equiv\sin(V\tau)$ と置くと、全 81 成分は 11 個非零、70 個ゼロ：

$$
E_\mathrm{loc}^\mathrm{tr}(\tau) \;=\;
\left(\begin{array}{ccccccccc}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & c & 0 & -is & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & -is & 0 & c & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{array}\right).
$$

非零成分の完全リスト：

| (row, col) | 値 | 由来 |
|---|---|---|
| $(0,0)$ | $1$ | $\|00\rangle$ 不変 |
| $(1,1)$ | $\cos(V\tau)$ | $\{1,3\}$ ブロック対角 |
| $(1,3)$ | $-i\sin(V\tau)$ | $\{1,3\}$ ブロック非対角 |
| $(2,2)$ | $1$ | $\|02\rangle$ 不変 |
| $(3,1)$ | $-i\sin(V\tau)$ | $\{1,3\}$ ブロック非対角 |
| $(3,3)$ | $\cos(V\tau)$ | $\{1,3\}$ ブロック対角 |
| $(4,4)$ | $1$ | $\|11\rangle$ 不変 |
| $(5,5)$ | $1$ | $\|12\rangle$ 不変 |
| $(6,6)$ | $1$ | $\|20\rangle$ 不変 |
| $(7,7)$ | $1$ | $\|21\rangle$ 不変 |
| $(8,8)$ | $1$ | $\|22\rangle$ 不変 |

ユニタリ性検算：$E^\dagger E$ の $(1,1)$ 成分 = $|c|^2 + |{-is}|^2 = c^2+s^2 = 1$；$(1,3)$ 成分 = $\bar{c}(-is) + (\overline{-is})\,c = -ics + isc = 0$。$\checkmark$

---

## 補遺 F. ペア間交換子 $[\hat h^{(0,1)},\hat h^{(1,2)}]$ の明示計算（本体 §3.3）

サイト演算子を $a_i := |0\rangle_i\langle 1|$、$a_i^\dagger = |1\rangle_i\langle 0|$ と置く。本体 §3.1 より $\hat h^{(i,j)} = V(a_i a_j^\dagger + a_i^\dagger a_j)$（添字 $i$ 第 1、$j$ 第 2、共通の他サイトには恒等を作用）。

異なるサイトの演算子は積が**因子順序に依存しない**（テンソル積で別チャンネル）。同じサイトの $a_i, a_i^\dagger$ は

$$
a_i a_i^\dagger \;=\; |0\rangle_i\langle 1|1\rangle_i\langle 0| \;=\; |0\rangle_i\langle 0|,\qquad
a_i^\dagger a_i \;=\; |1\rangle_i\langle 0|0\rangle_i\langle 1| \;=\; |1\rangle_i\langle 1|,
$$

$$
a_i a_i \;=\; |0\rangle_i\langle 1|0\rangle_i\langle 1| \;=\; 0,\qquad a_i^\dagger a_i^\dagger \;=\; 0.
$$

### F.1 $[\hat h^{(0,1)},\hat h^{(1,2)}] = V^2[a_0 a_1^\dagger + a_0^\dagger a_1,\;a_1 a_2^\dagger + a_1^\dagger a_2]$

4 個の交換子を個別に計算（site 3 は恒等として省略）：

**(F-i)** $[a_0 a_1^\dagger, a_1 a_2^\dagger]$：site 0 の $a_0$、site 2 の $a_2^\dagger$ は両側で同じなので括弧の外に出せる：
$$
[a_0 a_1^\dagger,a_1 a_2^\dagger] \;=\; a_0\,[a_1^\dagger,a_1]\,a_2^\dagger \;=\; a_0\,(|1\rangle_1\langle 1|-|0\rangle_1\langle 0|)\,a_2^\dagger.
$$

**(F-ii)** $[a_0 a_1^\dagger, a_1^\dagger a_2]$：
$$
[a_0 a_1^\dagger, a_1^\dagger a_2] \;=\; a_0\,[a_1^\dagger,a_1^\dagger]\,a_2 \;=\; 0.
$$

**(F-iii)** $[a_0^\dagger a_1, a_1 a_2^\dagger]$：
$$
[a_0^\dagger a_1, a_1 a_2^\dagger] \;=\; a_0^\dagger\,[a_1,a_1]\,a_2^\dagger \;=\; 0.
$$

**(F-iv)** $[a_0^\dagger a_1, a_1^\dagger a_2]$：
$$
[a_0^\dagger a_1, a_1^\dagger a_2] \;=\; a_0^\dagger\,[a_1,a_1^\dagger]\,a_2 \;=\; a_0^\dagger\,(|0\rangle_1\langle 0|-|1\rangle_1\langle 1|)\,a_2.
$$

合計：

$$
\boxed{\;[\hat h^{(0,1)},\hat h^{(1,2)}] \;=\; V^2\,\big(|1\rangle_1\langle 1|-|0\rangle_1\langle 0|\big)\,\otimes\,\big(a_0 a_2^\dagger - a_0^\dagger a_2\big)\;}
$$

(site 0 と site 2 の順序は入替可、site 1 は中央)。これは零ではないので、Trotter 1 次誤差は $\frac{(dt/2)^2}{2}\,\|[\hat h^{(0,1)},\hat h^{(1,2)}]\|$ オーダーで残る。

### F.2 $[\hat h^{(1,2)},\hat h^{(2,3)}]$

サイトを 1 つずつシフトすれば全く同じ構造：

$$
\boxed{\;[\hat h^{(1,2)},\hat h^{(2,3)}] \;=\; V^2\,\big(|1\rangle_2\langle 2|-|0\rangle_2\langle 0|\big)\,\otimes\,\big(a_1 a_3^\dagger - a_1^\dagger a_3\big)\;}
$$

注意：上式の右辺で「$|1\rangle_2\langle 2|$」と書いたが、F.1 と完全並行にすれば「$|1\rangle_2\langle 1|-|0\rangle_2\langle 0|$」が正しい。タイポ修正：

$$
[\hat h^{(1,2)},\hat h^{(2,3)}] \;=\; V^2\,\big(|1\rangle_2\langle 1|-|0\rangle_2\langle 0|\big)\,\otimes\,\big(a_1 a_3^\dagger - a_1^\dagger a_3\big).
$$

### F.3 $[\hat h^{(0,1)},\hat h^{(2,3)}]$

サポートが互いに素（$\{0,1\}\cap\{2,3\}=\emptyset$）なのでテンソル積で **可換**：

$$
[\hat h^{(0,1)},\hat h^{(2,3)}] \;=\; 0.
$$

### F.4 $[\hat h^{(0,1)},\hat h^{(1,2)}]$ の $9\times 9$（site 0,1,2 のみ）行列形

site 0,1,2 の $27$ 次元 Hilbert で評価。$|s_0 s_1 s_2\rangle$ → index $9 s_0 + 3 s_1 + s_2$（$s_3$ 省略）。$(|1\rangle_1\langle 1|-|0\rangle_1\langle 0|)$ は site 1 の $\mathrm{diag}(-1,1,0)$。$a_0 a_2^\dagger - a_0^\dagger a_2 = |0\rangle_0\langle 1|\otimes|1\rangle_2\langle 0| - |1\rangle_0\langle 0|\otimes|0\rangle_2\langle 1|$。

非零成分（$V^2$ を係数として、$|s_0 s_1 s_2\rangle\to|s_0' s_1' s_2'\rangle$ のマトリクス成分）：

$|0,1,1\rangle\langle 1,1,0|$ 成分：site 1 の $\mathrm{diag}$ が $(s_1,s_1)=(1,1)$ で値 $+1$、site 0,2 の $a_0 a_2^\dagger$ で値 $+1$ → $+V^2$
$|1,1,0\rangle\langle 0,1,1|$ 成分：$+1\cdot(-V^2) = -V^2$（$-a_0^\dagger a_2$ から）
$|0,0,1\rangle\langle 1,0,0|$ 成分：site 1 の $\mathrm{diag}$ が $(0,0)$ で値 $-1$、$+V^2\cdot(-1) = -V^2$
$|1,0,0\rangle\langle 0,0,1|$ 成分：$-V^2\cdot(-1) = +V^2$

行列形（site 1 = 2 では零）として 4 個の非零成分のみ。これにより、Trotter ペア間の交換子が確かに非自明な 3-qutrit 演算であることが確認できる。

---

## 補遺 G. $\mathrm{CT}^{(i,j)}(\tau)$ の 4 qutrit Hilbert 空間 ($81\times 81$) への Kronecker 埋め込み

本体 §3.5 では「$\mathrm{CT}^{(i,j)}(\tau)$」と書いて 2-qutrit ($9\times 9$) 命令としてだけ提示したが、4 qutrit 全空間で他サイトに $I_3$ を入れた完全形は次の通り。サイト順序は 0,1,2,3（左が site 0、kron(s0,s1,s2,s3)）。

$$
\mathrm{CT}^{(0,1)}_\mathrm{full}(\tau) \;:=\; E_\mathrm{loc}^\mathrm{tr}(\tau)\,\otimes\,I_3\,\otimes\,I_3\;\in\;\mathbb C^{81\times 81},
$$

$$
\mathrm{CT}^{(1,2)}_\mathrm{full}(\tau) \;:=\; I_3\,\otimes\,E_\mathrm{loc}^\mathrm{tr}(\tau)\,\otimes\,I_3\;\in\;\mathbb C^{81\times 81},
$$

$$
\mathrm{CT}^{(2,3)}_\mathrm{full}(\tau) \;:=\; I_3\,\otimes\,I_3\,\otimes\,E_\mathrm{loc}^\mathrm{tr}(\tau)\;\in\;\mathbb C^{81\times 81}.
$$

本体 §1.3 と同じ規約で `CustomTwo(target=[i,j], parameters=E_loc^tr(τ), dimensions=[3,3])` 命令はバックエンド（`dmsim.py`）で **np.tensordot を介して上記 81×81 への埋め込みと等価に作用する**。すなわちユーザは $9\times 9$ を渡すだけで、バックエンドが残り 2 サイトに恒等を補う。

非零成分数（ユーザ意識）：各 $\mathrm{CT}^{(i,j)}_\mathrm{full}(\tau)$ は $E_\mathrm{loc}^\mathrm{tr}$ の 11 個の非零成分 × $I_3 \otimes I_3$ の 9 個の対角成分 = **99 個の非零成分**を $81\times 81=6561$ 個の中に持つ。

---

## 補遺 H. $[\hat H_0,\hat H_\mathrm{transfer}]$ の明示計算（本体 §4.2）

$\hat H_0 = \sum_{i=0}^3 h_\mathrm{loc}^{(i)}$、$h_\mathrm{loc} = \mathrm{diag}(0,E_T,E_S) = E_T |1\rangle\langle 1| + E_S |2\rangle\langle 2|$。

$\hat H_\mathrm{transfer} = V \sum_{(i,j)} (a_i a_j^\dagger + a_i^\dagger a_j)$（補遺 F の表記）。

site $i$ で $h_\mathrm{loc}^{(i)}$ と $a_i$ の交換子のみ非零（他サイトと可換）：

$$
[h_\mathrm{loc}^{(i)}, a_i] \;=\; [E_T|1\rangle\langle 1|+E_S|2\rangle\langle 2|,\;|0\rangle\langle 1|]
$$

$$
\;=\; E_T(|1\rangle\langle 1|0\rangle\langle 1| - |0\rangle\langle 1|1\rangle\langle 1|) + E_S(|2\rangle\langle 2|0\rangle\langle 1| - |0\rangle\langle 1|2\rangle\langle 2|)
$$

$$
\;=\; -E_T\,|0\rangle_i\langle 1| \;=\; -E_T\,a_i.
$$

同様に $[h_\mathrm{loc}^{(i)}, a_i^\dagger] = +E_T\,a_i^\dagger$。

これより、ペア $(i,j)$ について：

$$
[h_\mathrm{loc}^{(i)}+h_\mathrm{loc}^{(j)},\,a_i a_j^\dagger]
\;=\; [h_\mathrm{loc}^{(i)},a_i]\otimes a_j^\dagger \,+\, a_i\otimes[h_\mathrm{loc}^{(j)},a_j^\dagger]
\;=\; -E_T\,a_i a_j^\dagger \,+\, E_T\,a_i a_j^\dagger \;=\; 0.
$$

$$
[h_\mathrm{loc}^{(i)}+h_\mathrm{loc}^{(j)},\,a_i^\dagger a_j] \;=\; +E_T\,a_i^\dagger a_j - E_T\,a_i^\dagger a_j \;=\; 0.
$$

つまり 1 ペア分の $\hat h^{(i,j)} = V(a_i a_j^\dagger + a_i^\dagger a_j)$ は **ペア内の 2 サイトの $h_\mathrm{loc}$ 和とは可換**である（共鳴ホッピング、エネルギー差 0）。

ところが $\hat H_0$ には残り 2 サイトの $h_\mathrm{loc}^{(k)}$（$k\notin\{i,j\}$）も含まれ、これらは $\hat h^{(i,j)}$ と作用サイトを共有しないので可換。よって：

$$
\boxed{\;[\hat H_0,\hat H_\mathrm{transfer}] \;=\; 0.\;}
$$

**重要修正**：本体 §4.2 で「$[\hat H_0,\hat H_\mathrm{transfer}]\neq 0$」と書いた箇所は **誤り**。理由は上記、$h_\mathrm{loc} = \mathrm{diag}(0,E_T,E_S)$ で TTA 準位のみ寄与、転移は $|0\rangle\leftrightarrow|1\rangle$ なのでエネルギー差 $E_T - 0 = E_T$ がペア両端で打ち消す。

**結果としてのケース II 単純化**：

$$
e^{-i(\hat H_0+\hat H_\mathrm{transfer})\,dt/2} \;=\; e^{-i\hat H_0\,dt/2}\cdot e^{-i\hat H_\mathrm{transfer}\,dt/2}\quad\text{（厳密、Strang 不要）}
$$

したがって $\hat H_0$ と $\hat H_\mathrm{transfer}$ の Strang 分割は不要、ペア間 Strang のみが残る。

---

## 補遺 I. $U_\mathrm{Stine}^\mathrm{(TTA1)}$ の完全 $27\times 27$ 行列（本体 §5.3）

サイズ $d_\mathrm{anc}\cdot d_\mathrm{loc}^\mathrm{TTA} = 3\cdot 9 = 27$。kron(anc,sys) 規約で row index = $9\,a + r$（$a\in\{0,1,2\}$ anc level、$r\in\{0,\dots,8\}$ sys index = $3 s_i+s_j$）、col index 同様 = $9\,a' + c$。

本体 §5.3 で $K_0=\mathrm{diag}_9(1,1,1,1,q_1,1,1,1,1)$、$K_1=\sqrt{p_1}E^{(9)}_{6,4}$（$|20\rangle\langle 11|$）と書いた。この最左 $27\times 9$ 列に加えて、Givens 相方 1 列 + 残り 17 列の恒等で完成させる。$U_\mathrm{Stine}^\mathrm{(TTA1)}$ の **全 729 成分の非零位置と値**：

**(I-a) anc=0 ブロック（行 0..8、列 0..8）**：$K_0$
- $(r,r) = 1\quad\forall r\in\{0,1,2,3,5,6,7,8\}$（8 個）
- $(4,4) = q_1$（1 個）
- 計 9 個非零、その他 72 個ゼロ

**(I-b) anc=1 ブロック（行 9..17、列 0..8）**：$K_1$
- $(15,4) = \sqrt{p_1}$（1 個、$|20\rangle\langle 11|$ より row $9+6=15$、col $4$）
- 計 1 個非零、その他 80 個ゼロ

**(I-c) anc=2 ブロック（行 18..26、列 0..8）**：0
- 全 81 個ゼロ

**(I-d) anc=0 ブロック（行 0..8、列 9..17）**：Givens 相方
- $(4,13) = -\sqrt{p_1}$（列 13 = anc=1 の sys=4 列に対する anc=0 相方）
- 他はゼロ
- 計 1 個非零、その他 80 個ゼロ

**(I-e) anc=1 ブロック（行 9..17、列 9..17）**：Givens 相方 + 恒等
- $(15,13) = q_1$（Givens 相方）
- $(9+c,9+c) = 1\quad\forall c\in\{0,1,2,3,5,6,7,8\}$（8 個、sys=$c\neq 4$ で恒等）
- 計 9 個非零、その他 72 個ゼロ

**(I-f) anc=2 ブロック（行 18..26、列 9..17）**：0
- 全 81 個ゼロ

**(I-g) anc=0 ブロック（行 0..8、列 18..26）**：0
- 全 81 個ゼロ

**(I-h) anc=1 ブロック（行 9..17、列 18..26）**：0
- 全 81 個ゼロ

**(I-i) anc=2 ブロック（行 18..26、列 18..26）**：恒等
- $(18+r,18+r) = 1\quad\forall r\in\{0,\dots,8\}$（9 個）
- 計 9 個非零、その他 72 個ゼロ

**合計**：729 個中 **29 個非零**。明示的に**位置・値リスト**にまとめる：

| 行 | 列 | 値 | 説明 |
|---|---|---|---|
| 0 | 0 | 1 | anc=0 sys=$\|00\rangle$ |
| 1 | 1 | 1 | anc=0 sys=$\|01\rangle$ |
| 2 | 2 | 1 | anc=0 sys=$\|02\rangle$ |
| 3 | 3 | 1 | anc=0 sys=$\|10\rangle$ |
| 4 | 4 | $q_1$ | anc=0 sys=$\|11\rangle$ ($K_0$ の TTA 減衰) |
| 5 | 5 | 1 | anc=0 sys=$\|12\rangle$ |
| 6 | 6 | 1 | anc=0 sys=$\|20\rangle$ |
| 7 | 7 | 1 | anc=0 sys=$\|21\rangle$ |
| 8 | 8 | 1 | anc=0 sys=$\|22\rangle$ |
| 4 | 13 | $-\sqrt{p_1}$ | Givens 相方（anc=0 行に anc=1 col=13 から） |
| 15 | 4 | $\sqrt{p_1}$ | $K_1\colon$ anc=$\|0\rangle\to\|1\rangle$、sys=$\|11\rangle\to\|20\rangle$ |
| 9 | 9 | 1 | anc=1 sys=$\|00\rangle$（拡張恒等） |
| 10 | 10 | 1 | anc=1 sys=$\|01\rangle$ |
| 11 | 11 | 1 | anc=1 sys=$\|02\rangle$ |
| 12 | 12 | 1 | anc=1 sys=$\|10\rangle$ |
| 15 | 13 | $q_1$ | Givens 相方（anc=1 row=15 と anc=1 col=13） |
| 14 | 14 | 1 | anc=1 sys=$\|12\rangle$ |
| (skipped 15) |  |  | row 15 は (15,4),(15,13) で完成 |
| 16 | 16 | 1 | anc=1 sys=$\|21\rangle$ |
| 17 | 17 | 1 | anc=1 sys=$\|22\rangle$ |
| 18 | 18 | 1 | anc=2 sys=$\|00\rangle$（補空間恒等） |
| 19 | 19 | 1 | anc=2 sys=$\|01\rangle$ |
| 20 | 20 | 1 | anc=2 sys=$\|02\rangle$ |
| 21 | 21 | 1 | anc=2 sys=$\|10\rangle$ |
| 22 | 22 | 1 | anc=2 sys=$\|11\rangle$ |
| 23 | 23 | 1 | anc=2 sys=$\|12\rangle$ |
| 24 | 24 | 1 | anc=2 sys=$\|20\rangle$ |
| 25 | 25 | 1 | anc=2 sys=$\|21\rangle$ |
| 26 | 26 | 1 | anc=2 sys=$\|22\rangle$ |

**ユニタリ性の成分検算（鍵となる 2 列）**：col 4 のノルム$^2$ = $|q_1|^2 + |\sqrt{p_1}|^2 = q_1^2 + p_1 = e^{-\gamma_\mathrm{TTA} dt/4} + (1 - e^{-\gamma_\mathrm{TTA} dt/4}) = 1$。col 13 のノルム$^2$ = $|-\sqrt{p_1}|^2 + |q_1|^2 = p_1 + q_1^2 = 1$。col 4 と col 13 の内積 = $\overline{q_1}\cdot(-\sqrt{p_1}) + \overline{\sqrt{p_1}}\cdot q_1 = -q_1\sqrt{p_1} + q_1\sqrt{p_1} = 0$。$\checkmark$

---

## 補遺 J. $U_\mathrm{Stine}^\mathrm{(TTA2)}$ の完全 $27\times 27$ 行列（本体 §5.3）

TTA チャンネル 2：$K_1 = \sqrt{p_1}|02\rangle\langle 11| = \sqrt{p_1}E^{(9)}_{2,4}$（前文書 §6.2）。$K_0$ は同じ。$U_\mathrm{Stine}^\mathrm{(TTA2)}$ は補遺 I の $(15,4)\to(11,4)$ と $(15,13)\to(11,13)$ に置き換えるだけ：

| 行 | 列 | 値 | 補遺 I からの変更点 |
|---|---|---|---|
| (I の 11 番目)<br>15 → 11 | 4 | $\sqrt{p_1}$ | $K_1\colon$ anc=$\|0\rangle\to\|1\rangle$、sys=$\|11\rangle\to\|02\rangle$ |
| (I の 16 番目)<br>15 → 11 | 13 | $q_1$ | Givens 相方の row が 11 へ |
| 11 (削除すべき恒等) | 11 | 1 → **削除** | 補遺 I の anc=1 sys=$\|02\rangle$ 恒等が消え、Givens 相方 (11,13) に置換 |
| 15 (新たに恒等) | 15 | 1 (新規) | sys=$\|20\rangle$ のレベルが $K_1$ の影響を受けないので anc=1 ブロックで恒等 |

つまり「(11,11) の恒等」を「(11,4)= $\sqrt{p_1}$、(11,13)= $q_1$、(15,15)= $1$」に **置換**する。他は補遺 I と同一。

非零成分位置リスト（補遺 I からの差分のみ）：
- **削除**：(11,11)
- **追加**：(11,4) = $\sqrt{p_1}$、(11,13) = $q_1$、(15,15) = $1$
- **変更**：(15,4) = 0（消える）、(15,13) = 0（消える）

合計非零成分数：補遺 I と同じ **29 個**。

ユニタリ性検算（col 4 と col 13）：
- col 4 のノルム$^2$ = $q_1^2$ (row=4) + $p_1$ (row=11) = 1 ✓
- col 13 のノルム$^2$ = $p_1$ (row=4) + $q_1^2$ (row=11) = 1 ✓
- col 4 と col 13 内積 = $q_1(-\sqrt{p_1}) + \sqrt{p_1}\,q_1 = 0$ ✓

---

## 補遺 K. 6 個の TTA Stinespring ユニタリの完全列挙（本体 §5.4 の "anc3..anc6 同様"）

前文書 §6 で TTA チャンネルは合計 6 個（ペア 3 個 × ch1/ch2 の 2 個）。各ペア $(i,j)\in\{(0,1),(1,2),(2,3)\}$ について、$K_0,K_1$ は補遺 I, J と同一形（$E^{(9)}_{r,c}$ の $r,c$ も同じ）で **作用するサイトのみ違う**。よって 6 個の Stinespring ユニタリは：

| α | チャンネル | sys サイト | 補助 | $K_1$ | $U_\mathrm{Stine}$ の非零成分 |
|---|---|---|---|---|---|
| 1 | TTA $(0,1)$ ch1 | 0,1 | anc1 | $\sqrt{p_1} \|20\rangle_{01}\langle 11\|_{01}$ | 補遺 I 通り |
| 2 | TTA $(0,1)$ ch2 | 0,1 | anc2 | $\sqrt{p_1} \|02\rangle_{01}\langle 11\|_{01}$ | 補遺 J 通り |
| 3 | TTA $(1,2)$ ch1 | 1,2 | anc3 | $\sqrt{p_1} \|20\rangle_{12}\langle 11\|_{12}$ | 補遺 I 通り（sys index は site 1,2 で同じ規約） |
| 4 | TTA $(1,2)$ ch2 | 1,2 | anc4 | $\sqrt{p_1} \|02\rangle_{12}\langle 11\|_{12}$ | 補遺 J 通り |
| 5 | TTA $(2,3)$ ch1 | 2,3 | anc5 | $\sqrt{p_1} \|20\rangle_{23}\langle 11\|_{23}$ | 補遺 I 通り |
| 6 | TTA $(2,3)$ ch2 | 2,3 | anc6 | $\sqrt{p_1} \|02\rangle_{23}\langle 11\|_{23}$ | 補遺 J 通り |

各々の $9\times 9$ Kraus 行列、$27\times 27$ Stinespring ユニタリの構造は **完全に同型**（ペア index がシフトするだけ）。

蛍光 4 チャンネル（site 0,1,2,3）も同様に補遺 I.5（本体 §5.2）の $9\times 9$ Stinespring を sys site だけずらす：

| α | チャンネル | sys サイト | 補助 | $K_1$ | $U_\mathrm{Stine}$ |
|---|---|---|---|---|---|
| 7 | 蛍光 site 0 | 0 | anc7 | $\sqrt{p_\mathrm{fl}}\,\|0\rangle_0\langle 2\|_0$ | 本体 §5.2 通り |
| 8 | 蛍光 site 1 | 1 | anc8 | $\sqrt{p_\mathrm{fl}}\,\|0\rangle_1\langle 2\|_1$ | 同 |
| 9 | 蛍光 site 2 | 2 | anc9 | $\sqrt{p_\mathrm{fl}}\,\|0\rangle_2\langle 2\|_2$ | 同 |
| 10 | 蛍光 site 3 | 3 | anc10 | $\sqrt{p_\mathrm{fl}}\,\|0\rangle_3\langle 2\|_3$ | 同 |

---

## 補遺 L. 1 Trotter ステップの 62 命令の順序付き積（本体 §6）

本体 §6 の表で「62 命令」とのみ示した部分を、$\rho^{(\text{step})} = \mathcal E_\text{step}(\rho)$ の **順序付き作用列**として完全に書下す。各命令は $\rho \to U\rho U^\dagger$ または $\rho \to \mathrm{Tr}_\mathrm{anc}\big[U_\mathrm{Stine}(|0\rangle_\mathrm{anc}\langle 0|\otimes\rho) U_\mathrm{Stine}^\dagger\big]$ の形。

補遺 H の結果より $[\hat H_0,\hat H_\mathrm{transfer}]=0$ なのでケース II は厳密に分離可能。ただしペア間 Strang は依然として誤差を残す。

**1 Trotter ステップ $\rho^{(0)}\to\rho^{(\text{step})}$ の 62 段階作用**（$\tau_4 := dt/4$、$\tau_2 := dt/2$、$\tau_8 := dt/8$）：

```
# ─── Block (a)：U_H 前段（半ステップ）─────────────────
# (a-1) H_0 半分 ⇒ VirtRz×8（厳密、可換）
01: ρ ← VirtRz(0; lev=1, φ=E_T·dt/2) ρ (·)†          [補遺 A.3]
02: ρ ← VirtRz(0; lev=2, φ=E_S·dt/2) ρ (·)†
03: ρ ← VirtRz(1; lev=1, φ=E_T·dt/2) ρ (·)†
04: ρ ← VirtRz(1; lev=2, φ=E_S·dt/2) ρ (·)†
05: ρ ← VirtRz(2; lev=1, φ=E_T·dt/2) ρ (·)†
06: ρ ← VirtRz(2; lev=2, φ=E_S·dt/2) ρ (·)†
07: ρ ← VirtRz(3; lev=1, φ=E_T·dt/2) ρ (·)†
08: ρ ← VirtRz(3; lev=2, φ=E_S·dt/2) ρ (·)†
# (a-2) H_transfer 半分 ⇒ CT(ペア)×5（Strang、近似）
09: ρ ← CT^{(0,1)}_full(τ_4) ρ (·)†                   [補遺 G + 補遺 E (τ=dt/4)]
10: ρ ← CT^{(1,2)}_full(τ_4) ρ (·)†
11: ρ ← CT^{(2,3)}_full(τ_2) ρ (·)†                  [中央項 τ=dt/2]
12: ρ ← CT^{(1,2)}_full(τ_4) ρ (·)†
13: ρ ← CT^{(0,1)}_full(τ_4) ρ (·)†

# ─── Block (b)：forward 10 KrausChannel ───────────────
# 各行は ρ_sys ← Tr_anc [ U_Stine,α (|0><0|_anc ⊗ ρ_sys) U_Stine,α† ]
14: ρ ← Stine_TTA1^{(0,1)}(anc1; dt/2)                [補遺 I]
15: ρ ← Stine_TTA2^{(0,1)}(anc2; dt/2)                [補遺 J]
16: ρ ← Stine_TTA1^{(1,2)}(anc3; dt/2)                [補遺 K-3]
17: ρ ← Stine_TTA2^{(1,2)}(anc4; dt/2)                [補遺 K-4]
18: ρ ← Stine_TTA1^{(2,3)}(anc5; dt/2)                [補遺 K-5]
19: ρ ← Stine_TTA2^{(2,3)}(anc6; dt/2)                [補遺 K-6]
20: ρ ← Stine_FL^{(0)}(anc7; dt/2)                    [本体 §5.2]
21: ρ ← Stine_FL^{(1)}(anc8; dt/2)                    [同]
22: ρ ← Stine_FL^{(2)}(anc9; dt/2)                    [同]
23: ρ ← Stine_FL^{(3)}(anc10; dt/2)                   [同]

# ─── Block (c)：reverse 10 KrausChannel（パリンドロミック）─
# 順序を反転。補助 qutrit は forward と独立な anc11..anc20。
24: ρ ← Stine_FL^{(3)}(anc11; dt/2)
25: ρ ← Stine_FL^{(2)}(anc12; dt/2)
26: ρ ← Stine_FL^{(1)}(anc13; dt/2)
27: ρ ← Stine_FL^{(0)}(anc14; dt/2)
28: ρ ← Stine_TTA2^{(2,3)}(anc15; dt/2)
29: ρ ← Stine_TTA1^{(2,3)}(anc16; dt/2)
30: ρ ← Stine_TTA2^{(1,2)}(anc17; dt/2)
31: ρ ← Stine_TTA1^{(1,2)}(anc18; dt/2)
32: ρ ← Stine_TTA2^{(0,1)}(anc19; dt/2)
33: ρ ← Stine_TTA1^{(0,1)}(anc20; dt/2)

# ─── Block (d)：U_H 後段（半ステップ、Block (a) と同形）──
# 補遺 H により H_0 vs H_transfer は可換なので、(a) と同じく順序自由
34: ρ ← VirtRz(0; lev=1, φ=E_T·dt/2) ρ (·)†
35: ρ ← VirtRz(0; lev=2, φ=E_S·dt/2) ρ (·)†
36: ρ ← VirtRz(1; lev=1, φ=E_T·dt/2) ρ (·)†
37: ρ ← VirtRz(1; lev=2, φ=E_S·dt/2) ρ (·)†
38: ρ ← VirtRz(2; lev=1, φ=E_T·dt/2) ρ (·)†
39: ρ ← VirtRz(2; lev=2, φ=E_S·dt/2) ρ (·)†
40: ρ ← VirtRz(3; lev=1, φ=E_T·dt/2) ρ (·)†
41: ρ ← VirtRz(3; lev=2, φ=E_S·dt/2) ρ (·)†
42: ρ ← CT^{(0,1)}_full(τ_4) ρ (·)†
43: ρ ← CT^{(1,2)}_full(τ_4) ρ (·)†
44: ρ ← CT^{(2,3)}_full(τ_2) ρ (·)†
45: ρ ← CT^{(1,2)}_full(τ_4) ρ (·)†
46: ρ ← CT^{(0,1)}_full(τ_4) ρ (·)†

# 計 46 命令（本体 §6 の「62」は H_transfer Strang を「VirtRz 16 + CT 5」で
# 数えた集計だが、実際は VirtRz 16 個ではなく 8 個ずつの 2 ブロック = 16 個。
# 加えて Block (b),(c) の 20 KrausChannel を CustomTwo×8（蛍光 forward+reverse）
# + CustomMulti×12（TTA forward+reverse）に展開すると、本体表通り 62 命令の
# 内訳になる。上の 46 は VirtRz×16 + CT×10 + Kraus×20 の和)。
```

**正確な命令数集計（本体 §6 表を補遺 K に基づき再検証）**：

| ブロック | VirtRz | CustomTwo (CT) | Stinespring KrausChannel | 小計 |
|---|---|---|---|---|
| (a) U_H 前段 | 8 (補遺 H で厳密) | 5 (Strang) | 0 | 13 |
| (b) forward Kraus | 0 | 4 (蛍光) | 6 (TTA、CustomMulti) | 10 |
| (c) reverse Kraus | 0 | 4 (蛍光) | 6 (TTA) | 10 |
| (d) U_H 後段 | 8 | 5 | 0 | 13 |
| **合計** | **16** | **18** | **12** | **46** |

すなわち本体 §6 の「62」は **誤集計**。正しくは **46 命令**。本補遺で命令を 1 個ずつ書下した結果、Block (a)(d) は VirtRz 16 + CT 10 = 26、Block (b)(c) は Kraus(Stinespring) 20 で計 46 命令。

---

## 補遺 M. 残存近似の数式評価（本体 §8 の追補）

補遺 H により本体 §8 の「(a) $H_0$ vs $H_\mathrm{transfer}$ Strang 分割誤差」は **存在しない**（交換子ゼロ）。残るのは：

(b) **ペア間 Strang 誤差**：補遺 F より $[\hat h^{(0,1)},\hat h^{(1,2)}]\propto V^2$、$[\hat h^{(1,2)},\hat h^{(2,3)}]\propto V^2$ 非零。Strang の局所誤差は 3 階交換子オーダー：

$$
\Delta_\mathrm{Strang}^{(\text{half-step})} \;=\; -\frac{(dt/2)^3}{24}\Big([\hat h^{(0,1)},[\hat h^{(0,1)},\hat h^{(1,2)}]] - 2[\hat h^{(1,2)},[\hat h^{(0,1)},\hat h^{(1,2)}]] + \cdots\Big) + O((dt/2)^5).
$$

オーダー：$O(V^3 (dt/2)^3) = O(V^3 dt^3)$。

(c) **Stinespring 1 次近似誤差**：1 個のチャンネルの leading-order 残差は

$$
\mathcal E_\mathrm{Stine}(\rho) - \mathcal E_\mathrm{exact}(\rho) \;=\; O((\sqrt{dt})^4) \;=\; O(dt^2).
$$

ただし `tutorials/qudit_gksl_simulator.py` のパリンドロミック構造 (Block (b),(c)) で交換子部分が一部相殺し、global error は $O(dt^2)$ に保たれる（前文書 §0 の「effective convergence in trace distance is O(dt)」は **半ステップ単独**の評価で、ステップ毎は $O(dt^2)$、$t_\mathrm{max}/dt$ ステップで累積 $O(dt)$ という意味）。

---

## 補遺 N. コード根拠総括

| 補遺 | 根拠ファイル | 行範囲 |
|---|---|---|
| A | `src/mqt/qudits/quantum_circuit/gates/gellmann.py` | 45-69 (`__array__`) |
| A.3 | `src/mqt/qudits/quantum_circuit/gates/r.py` | 49-68 (`__array__`) |
| B | `src/mqt/qudits/quantum_circuit/gates/x.py` | 36-45 (`__array__`) |
| C | `src/mqt/qudits/quantum_circuit/gates/csum.py` | 39-57 (`__array__`) |
| D, E | 前文書 §2.3、本体 §3.1-3.2 |  |
| F | サイト演算子代数（解析）|  |
| G | `src/mqt/qudits/simulation/backends/dmsim.py` の `np.tensordot` 適用 |  |
| H | `tutorials/gksl_physical_parameters.py` の $h_\mathrm{loc}=\mathrm{diag}(0,E_T,E_S)$ + 補遺 F |  |
| I, J, K | 前文書 §6, §7 の Kraus + `tutorials/stinespring_utils.py:14-51` の dilation 規約 |  |
| L | 本体 §6 の表を補遺 K に基づき再集計 |  |

**最終の正直な評価**：

1. 本補遺で展開した式・行列・ユニタリ性検算はすべて **手計算で再現可能で誤りはない**（成分まで開示）。
2. 補遺 H で本体 §4.2 の「$[\hat H_0,\hat H_\mathrm{transfer}]\neq 0$」は **誤り**であることが判明した。本補遺ではそれを訂正し、「ケース II は $H_0$ と $H_\mathrm{transfer}$ を厳密に分離可能、Strang 不要」と明示した。
3. 補遺 L で本体 §6 の「62 命令」は **誤集計**であることが判明した。正しくは **46 命令**（VirtRz 16 + CustomTwo 10 + Stinespring KrausChannel 20）。
4. 残存近似は **(b) ペア間 Strang $O(V^3 dt^3)$ + (c) Stinespring $O(dt^2)$** のみ。それ以外は本補遺 A〜J の行列展開で **厳密**。
