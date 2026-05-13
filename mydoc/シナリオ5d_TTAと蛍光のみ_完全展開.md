# シナリオ5d 縮約版：Lindblad チャンネルを TTA と蛍光のみに限定した完全展開

## 0. 本文書の前提（嘘禁止のための明示）

本文書は **Lindblad チャンネルを TTA（三重項–三重項消滅、ペア演算子）と 蛍光（単一サイト演算子）の 2 種類のみに限定した縮約モデル** について、すべての数式の項を一切省略せずに書き下す。

**実コード（`tutorials/gksl_math_utils.py::build_lindblad_operators`）が実装している Lindblad 演算子は次の 5 種類×サイト＋TTAペア の合計 26 個である**：

* TTA（ペア、$2$ チャンネル × $|\text{neighbors}|=3$ ＝ $6$ 個）
* 蛍光（単一サイト、$N=4$ 個）
* 燐光（単一サイト、$N=4$ 個）
* 内部転換 IC（単一サイト、$N=4$ 個）
* 項間交差 ISC$_{S\to T}$（単一サイト、$N=4$ 個）
* 項間交差 ISC$_{T\to S}$（単一サイト、$N=4$ 個）

本文書では**ユーザの要求に従い、燐光・IC・ISC$_{S\to T}$・ISC$_{T\to S}$ を取り除き、TTA と 蛍光のみ** の Lindblad 集合（合計 $2\cdot 3 + 1\cdot 4 = 10$ 個）を扱う。実コードの 26 個と異なるため、これは「シナリオ5d そのもの」ではなく**シナリオ5d の縮約版**であることを明記する。

ハミルトニアン $\hat H_\mathrm{total}$ は元のままとする（$E_T$、$E_S$、$V$ はいずれも残す；TTA は $|1\rangle$ から、蛍光は $|2\rangle$ からの遷移なので両エネルギーが必要）。

物理パラメータは `tutorials/gksl_physical_parameters.py` のデフォルト値を使う。$N=4$、$d=3$、$D=d^N=81$、$dt=t_\mathrm{max}/n_\mathrm{steps}=100/100=1$。隣接ペアは $\{(0,1),(1,2),(2,3)\}$。

---

## 1. ヒルベルト空間と基底

各分子は3準位系（qutrit）：

$$
|0\rangle = \mathrm{S}_0,\qquad |1\rangle = \mathrm{T}_1,\qquad |2\rangle = \mathrm{S}_1.
$$

多体ヒルベルト空間：

$$
\mathcal{H} \;=\; \mathcal{H}_0\otimes \mathcal{H}_1\otimes \mathcal{H}_2\otimes \mathcal{H}_3,\qquad \dim\mathcal{H} = 3^4 = 81.
$$

計算基底ラベル（コード `reduce(np.kron, [...])` の規約）：

$$
|s_0 s_1 s_2 s_3\rangle,\quad s_i\in\{0,1,2\},\qquad \mathrm{idx}(s_0,s_1,s_2,s_3) \;=\; 27\,s_0 + 9\,s_1 + 3\,s_2 + s_3 \;\in\;\{0,\dots,80\}.
$$

たとえば $|0000\rangle = $ idx 0、$|1001\rangle = 27+0+0+1 = 28$、$|2222\rangle = $ idx 80。

---

## 2. ハミルトニアン（全項を書き下し）

### 2.1 局所オンサイト行列 $h_\mathrm{loc}$（$3\times 3$）

$$
h_\mathrm{loc} \;=\; E_T\,|1\rangle\langle 1| \;+\; E_S\,|2\rangle\langle 2| \;=\; \begin{pmatrix} 0 & 0 & 0 \\ 0 & E_T & 0 \\ 0 & 0 & E_S \end{pmatrix}.
$$

### 2.2 オンサイト Hamiltonian $\hat H_0$（4 項の総和）

$\hat H_0$ は4分子の $h_\mathrm{loc}$ の単純な総和：

$$
\hat H_0 \;=\; h_\mathrm{loc}\otimes I_3\otimes I_3\otimes I_3 \;+\; I_3\otimes h_\mathrm{loc}\otimes I_3\otimes I_3 \;+\; I_3\otimes I_3\otimes h_\mathrm{loc}\otimes I_3 \;+\; I_3\otimes I_3\otimes I_3\otimes h_\mathrm{loc}.
$$

ここで $I_3$ は $3\times 3$ 単位行列。$\hat H_0$ は $81\times 81$ 対角行列で、対角成分は

$$
(\hat H_0)_{\mathrm{idx}(s_0 s_1 s_2 s_3),\mathrm{idx}(s_0 s_1 s_2 s_3)} \;=\; \sum_{i=0}^{3}\big( E_T\,\delta_{s_i,1} + E_S\,\delta_{s_i,2}\big) \;=\; E_T\,n_T(\vec s) + E_S\,n_S(\vec s),
$$

ここで $n_T(\vec s) = |\{i: s_i=1\}|$、$n_S(\vec s) = |\{i: s_i=2\}|$ は各構成において $\mathrm T_1$ にいる分子数、$\mathrm S_1$ にいる分子数。非対角成分はすべて 0。

### 2.3 移動 Hamiltonian $\hat H_\mathrm{transfer}$（3 ペア × 2 項 = 6 項の総和）

各隣接ペア $(i,j)$ について 2 項：

$$
\hat A^{(i,j)} \;:=\; |0\rangle_i\langle 1| \,\otimes\, |1\rangle_j\langle 0|,\qquad \hat A^{(i,j)\,\dagger} \;=\; |1\rangle_i\langle 0| \,\otimes\, |0\rangle_j\langle 1|.
$$

$\hat A^{(i,j)}$ の作用は基底状態に対して

$$
\hat A^{(i,j)}\,|s_0\cdots s_{N-1}\rangle \;=\; \delta_{s_i,1}\,\delta_{s_j,0}\,|s_0\cdots s_{i-1},0,s_{i+1},\dots,s_{j-1},1,s_{j+1},\dots\rangle.
$$

つまり「サイト $i$ にある $\mathrm T_1$ をサイト $j$ に移す」。

$$
\hat H_\mathrm{transfer} \;=\; V\,\sum_{(i,j)\in\{(0,1),(1,2),(2,3)\}}\Big(\hat A^{(i,j)} + \hat A^{(i,j)\,\dagger}\Big).
$$

これを完全に書き下すと

$$
\begin{aligned}
\hat H_\mathrm{transfer} \;=\; V\Big[\,
&\big(|0\rangle_0\langle 1|\otimes|1\rangle_1\langle 0|\otimes I_2\otimes I_3\big) + \big(|1\rangle_0\langle 0|\otimes|0\rangle_1\langle 1|\otimes I_2\otimes I_3\big) \\
+\;&\big(I_0\otimes|0\rangle_1\langle 1|\otimes|1\rangle_2\langle 0|\otimes I_3\big) + \big(I_0\otimes|1\rangle_1\langle 0|\otimes|0\rangle_2\langle 1|\otimes I_3\big) \\
+\;&\big(I_0\otimes I_1\otimes|0\rangle_2\langle 1|\otimes|1\rangle_3\langle 0|\big) + \big(I_0\otimes I_1\otimes|1\rangle_2\langle 0|\otimes|0\rangle_3\langle 1|\big)\,\Big].
\end{aligned}
$$

各 $|a\rangle_k\langle b|$ は site $k$ の qutrit 上の $3\times 3$ 行列 $E_{ab}$（$E_{ab}$ は $(a,b)$ 成分のみ 1、他 0）であり、$I_k$ はその site の $3\times 3$ 単位行列。

### 2.4 全 Hamiltonian

$$
\boxed{\;\hat H_\mathrm{total} \;=\; \hat H_0 \;+\; \hat H_\mathrm{transfer}\;}\qquad(81\times 81,\;\text{Hermitian})
$$

---

## 3. Lindblad 演算子（TTA と 蛍光のみ、合計 10 個を全部書き下す）

縮約版の GKSL 散逸子は

$$
\mathcal{L}_D \;=\; \sum_{\alpha=1}^{10} \mathcal{D}[\hat L_\alpha],\qquad \mathcal{D}[\hat L](\hat\rho) \;=\; \hat L\hat\rho\hat L^\dagger - \tfrac12\,\hat L^\dagger\hat L\,\hat\rho - \tfrac12\,\hat\rho\,\hat L^\dagger\hat L.
$$

10 個の演算子を以下に番号 $\alpha=1,\dots,10$ で固定する。

### 3.1 TTA（$\alpha=1,\dots,6$、3 ペア × 2 チャンネル）

略号：$\eta := \sqrt{\gamma_\mathrm{TTA}/2}$。各演算子は対応する 2 サイトの局所部分（$9\times 9$）と他 2 サイトの恒等のテンソル積：

| $\alpha$ | サポート $S_\alpha$ | 局所部分（site $i$ ⊗ site $j$、$9\times 9$） | 全体形 |
|---|---|---|---|
| 1 | $(i,j)=(0,1)$ | $\eta\,|2\rangle_0\langle 1| \otimes |0\rangle_1\langle 1|$ | $\hat L_1 = \eta\,(|2\rangle_0\langle 1|)\otimes(|0\rangle_1\langle 1|)\otimes I_2\otimes I_3$ |
| 2 | $(i,j)=(0,1)$ | $\eta\,|0\rangle_0\langle 1| \otimes |2\rangle_1\langle 1|$ | $\hat L_2 = \eta\,(|0\rangle_0\langle 1|)\otimes(|2\rangle_1\langle 1|)\otimes I_2\otimes I_3$ |
| 3 | $(i,j)=(1,2)$ | $\eta\,|2\rangle_1\langle 1| \otimes |0\rangle_2\langle 1|$ | $\hat L_3 = I_0\otimes(|2\rangle_1\langle 1|)\otimes(|0\rangle_2\langle 1|)\otimes I_3$ |
| 4 | $(i,j)=(1,2)$ | $\eta\,|0\rangle_1\langle 1| \otimes |2\rangle_2\langle 1|$ | $\hat L_4 = I_0\otimes(|0\rangle_1\langle 1|)\otimes(|2\rangle_2\langle 1|)\otimes I_3$ |
| 5 | $(i,j)=(2,3)$ | $\eta\,|2\rangle_2\langle 1| \otimes |0\rangle_3\langle 1|$ | $\hat L_5 = I_0\otimes I_1\otimes(|2\rangle_2\langle 1|)\otimes(|0\rangle_3\langle 1|)$ |
| 6 | $(i,j)=(2,3)$ | $\eta\,|0\rangle_2\langle 1| \otimes |2\rangle_3\langle 1|$ | $\hat L_6 = I_0\otimes I_1\otimes(|0\rangle_2\langle 1|)\otimes(|2\rangle_3\langle 1|)$ |

物理：両方とも 2 隣接サイトが共に $|11\rangle_{ij}$（双方とも $\mathrm T_1$）にあるとき発火し、片方を $|2\rangle = \mathrm S_1$、他方を $|0\rangle = \mathrm S_0$ に変える。チャンネル 1 は site $i$ が $\mathrm S_1$ になり site $j$ が $\mathrm S_0$ になる、チャンネル 2 は逆。

### 3.2 蛍光（$\alpha=7,8,9,10$、4 サイト）

略号：$\mu := \sqrt{\Gamma_\mathrm{fl}}$。各演算子は対応するサイトの局所部分（$3\times 3$）と他 3 サイトの恒等のテンソル積：

| $\alpha$ | サポート $S_\alpha$ | 局所部分（$3\times 3$） | 全体形 |
|---|---|---|---|
| 7  | site 0 | $\mu\,|0\rangle\langle 2|$ | $\hat L_7 = \mu\,(|0\rangle_0\langle 2|)\otimes I_1\otimes I_2\otimes I_3$ |
| 8  | site 1 | $\mu\,|0\rangle\langle 2|$ | $\hat L_8 = I_0\otimes \mu\,(|0\rangle_1\langle 2|)\otimes I_2\otimes I_3$ |
| 9  | site 2 | $\mu\,|0\rangle\langle 2|$ | $\hat L_9 = I_0\otimes I_1\otimes \mu\,(|0\rangle_2\langle 2|)\otimes I_3$ |
| 10 | site 3 | $\mu\,|0\rangle\langle 2|$ | $\hat L_{10} = I_0\otimes I_1\otimes I_2\otimes \mu\,(|0\rangle_3\langle 2|)$ |

物理：$\mathrm S_1\to\mathrm S_0$ の自発放射。

### 3.3 縮約 GKSL 方程式（10 項を完全に展開）

$$
\dot{\hat\rho}(t) \;=\; -i\,[\hat H_\mathrm{total},\hat\rho(t)] \;+\; \sum_{\alpha=1}^{10}\Big(\hat L_\alpha\,\hat\rho(t)\,\hat L_\alpha^\dagger \;-\; \tfrac12\,\hat L_\alpha^\dagger\hat L_\alpha\,\hat\rho(t) \;-\; \tfrac12\,\hat\rho(t)\,\hat L_\alpha^\dagger\hat L_\alpha\Big).
$$

これを $\alpha$ ごとに省略無く書き下すと：

$$
\begin{aligned}
\dot{\hat\rho} \;=\; -i\hat H_\mathrm{total}\hat\rho \;+\; i\hat\rho\hat H_\mathrm{total}\;
&+\;\hat L_1\hat\rho\hat L_1^\dagger - \tfrac12 \hat L_1^\dagger\hat L_1\hat\rho - \tfrac12 \hat\rho\hat L_1^\dagger\hat L_1\\
&+\;\hat L_2\hat\rho\hat L_2^\dagger - \tfrac12 \hat L_2^\dagger\hat L_2\hat\rho - \tfrac12 \hat\rho\hat L_2^\dagger\hat L_2\\
&+\;\hat L_3\hat\rho\hat L_3^\dagger - \tfrac12 \hat L_3^\dagger\hat L_3\hat\rho - \tfrac12 \hat\rho\hat L_3^\dagger\hat L_3\\
&+\;\hat L_4\hat\rho\hat L_4^\dagger - \tfrac12 \hat L_4^\dagger\hat L_4\hat\rho - \tfrac12 \hat\rho\hat L_4^\dagger\hat L_4\\
&+\;\hat L_5\hat\rho\hat L_5^\dagger - \tfrac12 \hat L_5^\dagger\hat L_5\hat\rho - \tfrac12 \hat\rho\hat L_5^\dagger\hat L_5\\
&+\;\hat L_6\hat\rho\hat L_6^\dagger - \tfrac12 \hat L_6^\dagger\hat L_6\hat\rho - \tfrac12 \hat\rho\hat L_6^\dagger\hat L_6\\
&+\;\hat L_7\hat\rho\hat L_7^\dagger - \tfrac12 \hat L_7^\dagger\hat L_7\hat\rho - \tfrac12 \hat\rho\hat L_7^\dagger\hat L_7\\
&+\;\hat L_8\hat\rho\hat L_8^\dagger - \tfrac12 \hat L_8^\dagger\hat L_8\hat\rho - \tfrac12 \hat\rho\hat L_8^\dagger\hat L_8\\
&+\;\hat L_9\hat\rho\hat L_9^\dagger - \tfrac12 \hat L_9^\dagger\hat L_9\hat\rho - \tfrac12 \hat\rho\hat L_9^\dagger\hat L_9\\
&+\;\hat L_{10}\hat\rho\hat L_{10}^\dagger - \tfrac12 \hat L_{10}^\dagger\hat L_{10}\hat\rho - \tfrac12 \hat\rho\hat L_{10}^\dagger\hat L_{10}.
\end{aligned}
$$

これが本縮約版の出発点となる方程式。形式解：

$$
\hat\rho(t) \;=\; e^{\mathcal{L}\,t}\hat\rho(0),\qquad \mathcal{L} = \mathcal{L}_H + \mathcal{L}_D,\quad \mathcal{L}_H(\hat\rho) = -i[\hat H_\mathrm{total},\hat\rho].
$$

---

## 4. 局所部分の積（$\hat L_\alpha^\dagger \hat L_\alpha$）を全部計算

ここから後の Strang 分割・Kraus 抽出に必要な $\hat L^\dagger\hat L$ をすべて陽に計算する。

### 4.1 TTA の局所 $\hat L^\dagger\hat L$（$9\times 9$）

2 サイト局所基底を $|ab\rangle$（$a,b\in\{0,1,2\}$、index = $3a+b\in\{0,\dots,8\}$）で並べる。順序 $|00\rangle,|01\rangle,|02\rangle,|10\rangle,|11\rangle,|12\rangle,|20\rangle,|21\rangle,|22\rangle$ が index $0,1,2,3,4,5,6,7,8$ に対応。

TTA チャンネル 1 の局所部分

$$
L_\mathrm{loc}^{(1)} \;=\; \eta\,|2\rangle\langle 1|\otimes |0\rangle\langle 1| \;=\; \eta\,|20\rangle\langle 11| \;=\; \eta\,E^{(9)}_{6,4}
$$

ここで $E^{(9)}_{a,b}$ は $9\times 9$ 行列で $(a,b)$ 成分が 1、他は 0。よって

$$
\big(L_\mathrm{loc}^{(1)}\big)^\dagger L_\mathrm{loc}^{(1)} \;=\; \eta^2\,|11\rangle\langle 20|\cdot|20\rangle\langle 11| \;=\; \eta^2\,|11\rangle\langle 11| \;=\; \tfrac{\gamma_\mathrm{TTA}}{2}\,E^{(9)}_{4,4}.
$$

これは $9\times 9$ 行列で $(4,4)$ 成分のみ $\gamma_\mathrm{TTA}/2$。

TTA チャンネル 2 の局所部分

$$
L_\mathrm{loc}^{(2)} \;=\; \eta\,|0\rangle\langle 1|\otimes |2\rangle\langle 1| \;=\; \eta\,|02\rangle\langle 11| \;=\; \eta\,E^{(9)}_{2,4},
$$

$$
\big(L_\mathrm{loc}^{(2)}\big)^\dagger L_\mathrm{loc}^{(2)} \;=\; \tfrac{\gamma_\mathrm{TTA}}{2}\,E^{(9)}_{4,4}.
$$

両者の和（同じ 2 サイトに作用する 2 チャンネルの合算）：

$$
\boxed{\;\big(L_\mathrm{loc}^{(1)}\big)^\dagger L_\mathrm{loc}^{(1)} + \big(L_\mathrm{loc}^{(2)}\big)^\dagger L_\mathrm{loc}^{(2)} \;=\; \gamma_\mathrm{TTA}\,E^{(9)}_{4,4} \;=\; \gamma_\mathrm{TTA}\,|11\rangle\langle 11|.\;}\qquad(\star_\mathrm{TTA})
$$

つまり TTA は局所 2 サイトが「両方 $\mathrm T_1$」の状態 $|11\rangle$ にあるときだけ発火し、その全脱励起率は $\gamma_\mathrm{TTA}$（2 チャンネルの和）。

### 4.2 蛍光の局所 $\hat L^\dagger\hat L$（$3\times 3$）

$$
L_\mathrm{loc}^{(\mathrm{fl})} \;=\; \mu\,|0\rangle\langle 2| \;=\; \mu\,E^{(3)}_{0,2},\qquad \big(L_\mathrm{loc}^{(\mathrm{fl})}\big)^\dagger L_\mathrm{loc}^{(\mathrm{fl})} \;=\; \mu^2\,|2\rangle\langle 2| \;=\; \Gamma_\mathrm{fl}\,E^{(3)}_{2,2} \;=\; \begin{pmatrix}0&0&0\\0&0&0\\0&0&\Gamma_\mathrm{fl}\end{pmatrix}.
$$

---

## 5. Strang + パリンドロミック分割（縮約版に合わせた書き直し）

時間発展超演算子 $\Lambda(dt) = e^{\mathcal{L}\,dt}$ を Strang 2 次分割：

$$
\Lambda(dt) \;\approx\; e^{\mathcal{L}_H\,dt/2}\;\Phi_\mathrm{pal}(dt)\;e^{\mathcal{L}_H\,dt/2},
$$

$$
\Phi_\mathrm{pal}(dt) \;=\; \Big(e^{\mathcal{D}_1\,dt/2}\,e^{\mathcal{D}_2\,dt/2}\,\cdots\,e^{\mathcal{D}_{10}\,dt/2}\Big)\;\Big(e^{\mathcal{D}_{10}\,dt/2}\,e^{\mathcal{D}_9\,dt/2}\,\cdots\,e^{\mathcal{D}_1\,dt/2}\Big),\qquad \mathcal{D}_\alpha := \mathcal{D}[\hat L_\alpha].
$$

合算で $\Lambda(dt) = e^{\mathcal{L}\,dt} + O(dt^3)$ ステップ毎、大域 $O(dt^2)$。

ハミルトニアン半ステップ：

$$
e^{\mathcal{L}_H\,dt/2}(\hat\rho) \;=\; U_H\,\hat\rho\,U_H^\dagger,\qquad U_H \;:=\; e^{-i\hat H_\mathrm{total}\,dt/2} \;\in\;\mathbb C^{81\times 81}.
$$

各局所散逸子半ステップ $e^{\mathcal{D}_\alpha\,dt/2}$ は局所空間（TTA は $9\times 9$、蛍光は $3\times 3$）でだけ非自明であり（局所性 $\mathcal{D}[L\otimes I] = \mathcal{D}^\mathrm{loc}[L]\otimes \mathrm{id}$）、補空間には恒等で作用する。よって**残る作業は、各局所散逸子の半ステップ Kraus 表現を出すこと**である。

---

## 6. TTA の半ステップ Kraus（解析形を完全に書き下す）

§3.1 の通り、TTA ペア $(i,j)$ には 2 つの局所 Lindblad 演算子 $L_\mathrm{loc}^{(1)} = \eta\,|20\rangle\langle 11|$、$L_\mathrm{loc}^{(2)} = \eta\,|02\rangle\langle 11|$ がある。これらは**同一 2 サイトに作用する 2 チャンネル**なので、**両者を合算したチャンネル**

$$
\mathcal{D}^\mathrm{loc}_\mathrm{TTA}(\rho_\mathrm{loc}) \;=\; \mathcal{D}\big[L_\mathrm{loc}^{(1)}\big](\rho_\mathrm{loc}) + \mathcal{D}\big[L_\mathrm{loc}^{(2)}\big](\rho_\mathrm{loc})
$$

の半ステップ指数 $e^{\mathcal{D}^\mathrm{loc}_\mathrm{TTA}\,dt/2}$ を出すのが最もきれい。ただし**実コードのパリンドロミック分割では各 $\alpha$ ごとに独立に $e^{\mathcal{D}_\alpha\,dt/2}$ を取る**ので、本節ではまず**1 チャンネルずつ**の Kraus を解析形で書き、§6.4 で 2 チャンネル合算版も書く。

### 6.1 単一 TTA チャンネル $L = \eta\,|20\rangle\langle 11|$ の半ステップを解く

$\rho_\mathrm{loc}$ は $9\times 9$ 行列で $\rho_{ab,cd}$（$a,b,c,d\in\{0,1,2\}$）と書く（$\rho_{(3a+b),(3c+d)}$）。$L^\dagger L = \eta^2 |11\rangle\langle 11|$、$L\rho L^\dagger = \eta^2\,\rho_{11,11}\,|20\rangle\langle 20|$。

$\dot\rho_\mathrm{loc} = \mathcal{D}[L](\rho_\mathrm{loc})$ を成分ごとに書くと（$\eta^2 = \gamma_\mathrm{TTA}/2$ と置いて略号 $\gamma' := \eta^2$）：

* $\dot\rho_{11,11} = -2\gamma'\,\rho_{11,11}\cdot\tfrac12 -\tfrac12\cdot 2\gamma'\,\rho_{11,11} = -\gamma'\rho_{11,11} \cdot \tfrac{2}{2}\cdot 1$ ←誤りやすいので**丁寧に**：

  $\mathcal{D}[L](\rho)_{ab,cd} = \eta^2\delta_{ab,20}\delta_{cd,20}\rho_{11,11} - \tfrac{\eta^2}{2}\delta_{ab,11}\rho_{11,cd} - \tfrac{\eta^2}{2}\rho_{ab,11}\delta_{cd,11}$

  $(ab,cd) = (11,11)$ のとき：$0 - \tfrac{\eta^2}{2}\rho_{11,11} - \tfrac{\eta^2}{2}\rho_{11,11} = -\eta^2\rho_{11,11} = -(\gamma_\mathrm{TTA}/2)\rho_{11,11}$.

* $(ab,cd) = (20,20)$：$\eta^2\rho_{11,11} - 0 - 0 = (\gamma_\mathrm{TTA}/2)\rho_{11,11}$.
* $(ab,cd) = (11,cd)$ で $cd\neq 11$：$0 - (\eta^2/2)\rho_{11,cd} - 0 = -(\gamma_\mathrm{TTA}/4)\rho_{11,cd}$.
* $(ab,cd) = (ab,11)$ で $ab\neq 11$：$0 - 0 - (\gamma_\mathrm{TTA}/4)\rho_{ab,11}$.
* それ以外：0。

時刻 $\tau = dt/2$ における解（線形 ODE で各 $\rho$ 成分が他成分から独立に決まる）：

$$
\rho_{11,11}(\tau) = e^{-\eta^2\tau}\,\rho_{11,11}(0),\qquad \eta^2 = \gamma_\mathrm{TTA}/2.
$$

$\rho_{20,20}(\tau)$ は $\dot\rho_{20,20} = \eta^2\,\rho_{11,11}(t) = \eta^2 e^{-\eta^2 t}\rho_{11,11}(0)$ より

$$
\rho_{20,20}(\tau) = \rho_{20,20}(0) + (1 - e^{-\eta^2\tau})\,\rho_{11,11}(0).
$$

オフ対角 $\rho_{11,cd}$（$cd\neq 11$）と $\rho_{ab,11}$（$ab\neq 11$）：$\dot\rho_{11,cd} = -(\eta^2/2)\rho_{11,cd}$ より

$$
\rho_{11,cd}(\tau) = e^{-\eta^2\tau/2}\,\rho_{11,cd}(0),\qquad \rho_{ab,11}(\tau) = e^{-\eta^2\tau/2}\,\rho_{ab,11}(0).
$$

その他の成分：$\rho_{ab,cd}(\tau) = \rho_{ab,cd}(0)$（$ab,cd \notin \{11\}$ かつ $(ab,cd)\neq(20,20)$）。

このチャンネルの Kraus 表現：定数 $p_1 := 1 - e^{-\eta^2\tau} = 1 - e^{-(\gamma_\mathrm{TTA}/2)\,dt/2} = 1 - e^{-\gamma_\mathrm{TTA}\,dt/4}$、$q_1 := \sqrt{1-p_1} = e^{-\gamma_\mathrm{TTA}\,dt/8}$ と置くと（$q_1^2 = 1-p_1$）、

$$
\boxed{\;K^{(1,\mathrm{half})}_0 \;=\; \mathrm{diag}_{9}(1,1,1,1,q_1,1,1,1,1),\qquad K^{(1,\mathrm{half})}_1 \;=\; \sqrt{p_1}\;|20\rangle\langle 11| \;=\; \sqrt{p_1}\,E^{(9)}_{6,4}.\;}
$$

完全性検証：

$$
\big(K^{(1,\mathrm{half})}_0\big)^\dagger K^{(1,\mathrm{half})}_0 + \big(K^{(1,\mathrm{half})}_1\big)^\dagger K^{(1,\mathrm{half})}_1 \;=\; \mathrm{diag}_9(1,1,1,1,q_1^2,1,1,1,1) + p_1\,E^{(9)}_{4,4} \;=\; \mathrm{diag}_9(1,1,1,1,q_1^2+p_1,1,1,1,1) \;=\; I_9.\;\checkmark
$$

作用検証（任意 $\rho_\mathrm{loc}$ に対して）：$K^{(1,\mathrm{half})}_0\rho K^{(1,\mathrm{half})\dagger}_0$ は $\rho_{ab,cd}$ を「row index が 4 なら $q_1$ 倍、column index が 4 なら $q_1$ 倍」する、つまり $\rho_{11,11}\to q_1^2\rho_{11,11} = (1-p_1)\rho_{11,11}$、$\rho_{11,cd}\to q_1\rho_{11,cd}$ ($cd\neq 4$)、$\rho_{ab,11}\to q_1\rho_{ab,11}$ ($ab\neq 4$)、その他は不変。$K^{(1,\mathrm{half})}_1\rho K^{(1,\mathrm{half})\dagger}_1 = p_1\rho_{11,11}|20\rangle\langle 20|$。両者の和は §6.1 冒頭の解析解と一致。$\checkmark$

### 6.2 単一 TTA チャンネル $L = \eta\,|02\rangle\langle 11|$ の半ステップ

§6.1 と完全に並行（$|20\rangle$ を $|02\rangle$ に置き換え、$E^{(9)}_{6,4}$ を $E^{(9)}_{2,4}$ に置き換え）：

$$
\boxed{\;K^{(2,\mathrm{half})}_0 \;=\; \mathrm{diag}_{9}(1,1,1,1,q_1,1,1,1,1),\qquad K^{(2,\mathrm{half})}_1 \;=\; \sqrt{p_1}\;|02\rangle\langle 11| \;=\; \sqrt{p_1}\,E^{(9)}_{2,4}.\;}
$$

ここで $p_1, q_1$ は §6.1 と同じ。

### 6.3 §6.1 と §6.2 を続けて適用すると何が起こるか

実コードのパリンドロミック分割は $e^{\mathcal{D}_1\,dt/2}\circ e^{\mathcal{D}_2\,dt/2}$ を適用する。これは数学的には**合算 2 チャンネル**の半ステップ $e^{(\mathcal{D}_1+\mathcal{D}_2)\,dt/2}$ とは **完全には一致しない**（$[\mathcal{D}_1,\mathcal{D}_2]\neq 0$ 由来の Trotter 誤差 $O(dt^2)$ がペア内に残る）。本縮約版でもこれを変更しない（実コード一致のため）。

### 6.4 参考：合算 2 チャンネルの半ステップ Kraus（Trotter 誤差なしの解析形）

参考までに、両 TTA チャンネルを**同時に**指数化した場合の Kraus を書く（縮約版で「もし」パリンドロミックではなく合算したらどうなるか）。$\mathcal D_\mathrm{TTA} = \mathcal D[L^{(1)}] + \mathcal D[L^{(2)}]$、$L^{(\alpha)\dagger}L^{(\alpha)}$ の和は $(\star_\mathrm{TTA})$ より $\gamma_\mathrm{TTA}|11\rangle\langle 11|$。$|11\rangle$ は 2 つの jump operator $|20\rangle\langle 11|$（重み $\eta^2$）と $|02\rangle\langle 11|$（重み $\eta^2$）でデポピュレートされる。

定数：$P := 1 - e^{-\gamma_\mathrm{TTA}\,dt/2}$、$Q := \sqrt{1-P} = e^{-\gamma_\mathrm{TTA}\,dt/4}$。

$$
K^{\mathrm{TTA,half}}_0 \;=\; \mathrm{diag}_9(1,1,1,1,Q,1,1,1,1),\qquad K^{\mathrm{TTA,half}}_1 \;=\; \sqrt{P/2}\,E^{(9)}_{6,4},\qquad K^{\mathrm{TTA,half}}_2 \;=\; \sqrt{P/2}\,E^{(9)}_{2,4}.
$$

完全性：$\mathrm{diag}_9(1,\ldots,1,Q^2,1,\ldots,1) + (P/2)E^{(9)}_{4,4} + (P/2)E^{(9)}_{4,4} = \mathrm{diag}_9(1,\ldots,1,Q^2+P,1,\ldots,1) = I_9$. $\checkmark$

**※実コードはこちらの合算形は使わない**（§6.3 の通りパリンドロミックで個別指数化）。本節 §6.4 は理解の補助のために掲示しただけである。

---

## 7. 蛍光の半ステップ Kraus（解析形を完全に書き下す）

§3.2 の通り、各サイト $i$ の局所 Lindblad は $L_\mathrm{loc}^{(\mathrm{fl})} = \mu\,|0\rangle\langle 2|$（$3\times 3$）。$L^\dagger L = \Gamma_\mathrm{fl}|2\rangle\langle 2|$、$L\rho L^\dagger = \Gamma_\mathrm{fl}\rho_{22}|0\rangle\langle 0|$。$\dot\rho_\mathrm{loc} = \mathcal D[L](\rho_\mathrm{loc})$ を成分で書くと：

* $\dot\rho_{22} = -\Gamma_\mathrm{fl}\rho_{22}$
* $\dot\rho_{00} = +\Gamma_\mathrm{fl}\rho_{22}$
* $\dot\rho_{11} = 0$
* $\dot\rho_{02} = -(\Gamma_\mathrm{fl}/2)\rho_{02}$、$\dot\rho_{20} = -(\Gamma_\mathrm{fl}/2)\rho_{20}$
* $\dot\rho_{12} = -(\Gamma_\mathrm{fl}/2)\rho_{12}$、$\dot\rho_{21} = -(\Gamma_\mathrm{fl}/2)\rho_{21}$
* $\dot\rho_{01} = 0$、$\dot\rho_{10} = 0$

時刻 $\tau = dt/2$ での解：定数 $p_\mathrm{fl} := 1 - e^{-\Gamma_\mathrm{fl}\,dt/2}$、$q_\mathrm{fl} := \sqrt{1-p_\mathrm{fl}} = e^{-\Gamma_\mathrm{fl}\,dt/4}$ と置くと（$q_\mathrm{fl}^2 = 1-p_\mathrm{fl}$）、

$$
\rho_{22}(\tau) = (1-p_\mathrm{fl})\rho_{22}(0),\quad \rho_{00}(\tau) = \rho_{00}(0) + p_\mathrm{fl}\rho_{22}(0),\quad \rho_{11}(\tau) = \rho_{11}(0),
$$

$$
\rho_{02}(\tau) = q_\mathrm{fl}\rho_{02}(0),\quad \rho_{20}(\tau) = q_\mathrm{fl}\rho_{20}(0),\quad \rho_{12}(\tau) = q_\mathrm{fl}\rho_{12}(0),\quad \rho_{21}(\tau) = q_\mathrm{fl}\rho_{21}(0),
$$

$$
\rho_{01}(\tau) = \rho_{01}(0),\quad \rho_{10}(\tau) = \rho_{10}(0).
$$

Kraus 表現（標準的な「3 準位 amplitude damping from $|2\rangle$ to $|0\rangle$」）：

$$
\boxed{\;K^{(\mathrm{fl},\mathrm{half})}_0 \;=\; \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & q_\mathrm{fl} \end{pmatrix},\qquad K^{(\mathrm{fl},\mathrm{half})}_1 \;=\; \sqrt{p_\mathrm{fl}}\,|0\rangle\langle 2| \;=\; \begin{pmatrix} 0 & 0 & \sqrt{p_\mathrm{fl}} \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.\;}
$$

完全性：

$$
\big(K^{(\mathrm{fl},\mathrm{half})}_0\big)^\dagger K^{(\mathrm{fl},\mathrm{half})}_0 + \big(K^{(\mathrm{fl},\mathrm{half})}_1\big)^\dagger K^{(\mathrm{fl},\mathrm{half})}_1 \;=\; \mathrm{diag}(1,1,q_\mathrm{fl}^2) + \mathrm{diag}(0,0,p_\mathrm{fl}) \;=\; \mathrm{diag}(1,1,1) \;=\; I_3.\;\checkmark
$$

作用検証：$K_0\rho K_0^\dagger$ は $\rho$ の row 2 と column 2 をそれぞれ $q_\mathrm{fl}$ 倍するので $\rho_{22}\to q_\mathrm{fl}^2\rho_{22} = (1-p_\mathrm{fl})\rho_{22}$、$\rho_{02}, \rho_{20}, \rho_{12}, \rho_{21}\to q_\mathrm{fl}$ 倍、$\rho_{00}, \rho_{11}, \rho_{01}, \rho_{10}\to$ 不変。$K_1\rho K_1^\dagger = p_\mathrm{fl}\rho_{22}\,|0\rangle\langle 0|$、これは $\rho_{00}$ に $p_\mathrm{fl}\rho_{22}$ を加える。両者の和が §7 冒頭の解析解と一致。$\checkmark$

---

## 8. 全 10 チャンネルの半ステップ Kraus（一覧）

§6, §7 から、縮約版の 10 個の半ステップ Kraus 演算子集合 $\{K_{\alpha,\beta}\}$ は次の通り：

| $\alpha$ | チャンネル | サポート $S_\alpha$ | 局所 Kraus $\{K_{\alpha,\beta}^\mathrm{loc}\}$ | rank $r_\alpha$ |
|---|---|---|---|---|
| 1 | TTA $(0,1)$ ch1 | $\{0,1\}$ | $K_0^{(\mathrm{TTA1})}$, $K_1^{(\mathrm{TTA1})}$（§6.1） | 2 |
| 2 | TTA $(0,1)$ ch2 | $\{0,1\}$ | $K_0^{(\mathrm{TTA2})}$, $K_1^{(\mathrm{TTA2})}$（§6.2） | 2 |
| 3 | TTA $(1,2)$ ch1 | $\{1,2\}$ | $K_0^{(\mathrm{TTA1})}$, $K_1^{(\mathrm{TTA1})}$ | 2 |
| 4 | TTA $(1,2)$ ch2 | $\{1,2\}$ | $K_0^{(\mathrm{TTA2})}$, $K_1^{(\mathrm{TTA2})}$ | 2 |
| 5 | TTA $(2,3)$ ch1 | $\{2,3\}$ | $K_0^{(\mathrm{TTA1})}$, $K_1^{(\mathrm{TTA1})}$ | 2 |
| 6 | TTA $(2,3)$ ch2 | $\{2,3\}$ | $K_0^{(\mathrm{TTA2})}$, $K_1^{(\mathrm{TTA2})}$ | 2 |
| 7 | 蛍光 site 0 | $\{0\}$ | $K_0^{(\mathrm{fl})}$, $K_1^{(\mathrm{fl})}$（§7） | 2 |
| 8 | 蛍光 site 1 | $\{1\}$ | $K_0^{(\mathrm{fl})}$, $K_1^{(\mathrm{fl})}$ | 2 |
| 9 | 蛍光 site 2 | $\{2\}$ | $K_0^{(\mathrm{fl})}$, $K_1^{(\mathrm{fl})}$ | 2 |
| 10 | 蛍光 site 3 | $\{3\}$ | $K_0^{(\mathrm{fl})}$, $K_1^{(\mathrm{fl})}$ | 2 |

すべて **rank 2**。各 $K_{\alpha,\beta}^\mathrm{loc}$ は局所空間（TTA は $9\times 9$、蛍光は $3\times 3$）の行列。全空間に持ち上げた版

$$
\tilde K_{\alpha,\beta} \;=\; \big(K_{\alpha,\beta}^\mathrm{loc}\big)_{S_\alpha}\,\otimes\,I_{\bar S_\alpha}\;\in\;\mathbb C^{81\times 81}
$$

を用いる。具体例（site レイアウト $(0,1,2,3)$）：

* $\alpha=1$（TTA $(0,1)$ ch1）：$\tilde K_{1,\beta} = K_{\beta}^\mathrm{loc}\otimes I_3\otimes I_3$
* $\alpha=3$（TTA $(1,2)$ ch1）：$\tilde K_{3,\beta} = I_3\otimes K_{\beta}^\mathrm{loc}\otimes I_3$
* $\alpha=5$（TTA $(2,3)$ ch1）：$\tilde K_{5,\beta} = I_3\otimes I_3\otimes K_{\beta}^\mathrm{loc}$
* $\alpha=7$（蛍光 site 0）：$\tilde K_{7,\beta} = K_{\beta}^\mathrm{loc}\otimes I_3\otimes I_3\otimes I_3$
* $\alpha=10$（蛍光 site 3）：$\tilde K_{10,\beta} = I_3\otimes I_3\otimes I_3\otimes K_{\beta}^\mathrm{loc}$

CPTP 検証（DMSim の `KrausChannel.__init__` が実測する条件）：各 $\alpha$ について $\sum_\beta \big(K_{\alpha,\beta}^\mathrm{loc}\big)^\dagger K_{\alpha,\beta}^\mathrm{loc} = I$ が §6.1, §6.2, §7 で解析的に証明済み。よって $\sum_\beta \tilde K_{\alpha,\beta}^\dagger \tilde K_{\alpha,\beta} = I_{81}$ も自動成立（$(A^\dagger A)\otimes I = I_{d_\mathrm{loc}}\otimes I_{\bar S} = I_{81}$）。

---

## 9. 1 Trotter ステップの量子回路（22 命令を全部書き下す）

縮約版の 1 Trotter ステップは：

* `cu_multi`：2 個（前段・後段の Hamiltonian 半ステップ）
* `kraus_channel`：forward 列 10 個 + reverse 列 10 個 ＝ 20 個
* **合計 22 命令**

`tutorials/qudit_gksl_simulator.py::_trotter_step_dmsim` の append 順に従って完全に書き下す（縮約版 10 チャンネル用に番号を再付与）：

```
QuantumRegister "sys": 4 qutrits, dimensions = [3, 3, 3, 3]

(a)   cu_multi(target=[0,1,2,3], parameters=U_H)                          ← Hamiltonian 半ステップ 前段
(b1)  kraus_channel(target=[0,1], kraus_operators=[K^TTA1_0, K^TTA1_1])    ← α=1: TTA(0,1) ch1
(b2)  kraus_channel(target=[0,1], kraus_operators=[K^TTA2_0, K^TTA2_1])    ← α=2: TTA(0,1) ch2
(b3)  kraus_channel(target=[1,2], kraus_operators=[K^TTA1_0, K^TTA1_1])    ← α=3: TTA(1,2) ch1
(b4)  kraus_channel(target=[1,2], kraus_operators=[K^TTA2_0, K^TTA2_1])    ← α=4: TTA(1,2) ch2
(b5)  kraus_channel(target=[2,3], kraus_operators=[K^TTA1_0, K^TTA1_1])    ← α=5: TTA(2,3) ch1
(b6)  kraus_channel(target=[2,3], kraus_operators=[K^TTA2_0, K^TTA2_1])    ← α=6: TTA(2,3) ch2
(b7)  kraus_channel(target=[0],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← α=7: 蛍光 site 0
(b8)  kraus_channel(target=[1],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← α=8: 蛍光 site 1
(b9)  kraus_channel(target=[2],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← α=9: 蛍光 site 2
(b10) kraus_channel(target=[3],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← α=10: 蛍光 site 3

(c1)  kraus_channel(target=[3],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← reverse: α=10
(c2)  kraus_channel(target=[2],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← reverse: α=9
(c3)  kraus_channel(target=[1],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← reverse: α=8
(c4)  kraus_channel(target=[0],   kraus_operators=[K^fl_0,   K^fl_1  ])    ← reverse: α=7
(c5)  kraus_channel(target=[2,3], kraus_operators=[K^TTA2_0, K^TTA2_1])    ← reverse: α=6
(c6)  kraus_channel(target=[2,3], kraus_operators=[K^TTA1_0, K^TTA1_1])    ← reverse: α=5
(c7)  kraus_channel(target=[1,2], kraus_operators=[K^TTA2_0, K^TTA2_1])    ← reverse: α=4
(c8)  kraus_channel(target=[1,2], kraus_operators=[K^TTA1_0, K^TTA1_1])    ← reverse: α=3
(c9)  kraus_channel(target=[0,1], kraus_operators=[K^TTA2_0, K^TTA2_1])    ← reverse: α=2
(c10) kraus_channel(target=[0,1], kraus_operators=[K^TTA1_0, K^TTA1_1])    ← reverse: α=1

(d)   cu_multi(target=[0,1,2,3], parameters=U_H)                          ← Hamiltonian 半ステップ 後段
```

ここで $U_H = e^{-i\hat H_\mathrm{total}\,dt/2}$、各 `K^{...}` は §6, §7 の解析形（局所空間の行列、全空間への持ち上げは DMSim の `apply_kraus_to_density` が `np.tensordot` で自動的に行う）。

---

## 10. 22 命令を $\rho$ に作用させた完全式（1 ステップ）

DMSim の実装（`src/mqt/qudits/simulation/backends/dmsim.py::_apply_instruction`）は、命令ごとに $\rho$ を更新する。$\rho^{(k)}$ を $k$ 番目の命令適用後の状態として、22 行を一切省略せずに書き下すと：

$$
\rho^{(0)} := \rho_n.
$$

**(a) cu_multi 前段**

$$
\rho^{(1)} \;=\; U_H\,\rho^{(0)}\,U_H^\dagger.
$$

**(b1)–(b6) TTA forward 6 命令**

$$
\rho^{(2)} \;=\; \tilde K_{1,0}\,\rho^{(1)}\,\tilde K_{1,0}^\dagger \;+\; \tilde K_{1,1}\,\rho^{(1)}\,\tilde K_{1,1}^\dagger,
$$

$$
\rho^{(3)} \;=\; \tilde K_{2,0}\,\rho^{(2)}\,\tilde K_{2,0}^\dagger \;+\; \tilde K_{2,1}\,\rho^{(2)}\,\tilde K_{2,1}^\dagger,
$$

$$
\rho^{(4)} \;=\; \tilde K_{3,0}\,\rho^{(3)}\,\tilde K_{3,0}^\dagger \;+\; \tilde K_{3,1}\,\rho^{(3)}\,\tilde K_{3,1}^\dagger,
$$

$$
\rho^{(5)} \;=\; \tilde K_{4,0}\,\rho^{(4)}\,\tilde K_{4,0}^\dagger \;+\; \tilde K_{4,1}\,\rho^{(4)}\,\tilde K_{4,1}^\dagger,
$$

$$
\rho^{(6)} \;=\; \tilde K_{5,0}\,\rho^{(5)}\,\tilde K_{5,0}^\dagger \;+\; \tilde K_{5,1}\,\rho^{(5)}\,\tilde K_{5,1}^\dagger,
$$

$$
\rho^{(7)} \;=\; \tilde K_{6,0}\,\rho^{(6)}\,\tilde K_{6,0}^\dagger \;+\; \tilde K_{6,1}\,\rho^{(6)}\,\tilde K_{6,1}^\dagger.
$$

**(b7)–(b10) 蛍光 forward 4 命令**

$$
\rho^{(8)} \;=\; \tilde K_{7,0}\,\rho^{(7)}\,\tilde K_{7,0}^\dagger \;+\; \tilde K_{7,1}\,\rho^{(7)}\,\tilde K_{7,1}^\dagger,
$$

$$
\rho^{(9)} \;=\; \tilde K_{8,0}\,\rho^{(8)}\,\tilde K_{8,0}^\dagger \;+\; \tilde K_{8,1}\,\rho^{(8)}\,\tilde K_{8,1}^\dagger,
$$

$$
\rho^{(10)} \;=\; \tilde K_{9,0}\,\rho^{(9)}\,\tilde K_{9,0}^\dagger \;+\; \tilde K_{9,1}\,\rho^{(9)}\,\tilde K_{9,1}^\dagger,
$$

$$
\rho^{(11)} \;=\; \tilde K_{10,0}\,\rho^{(10)}\,\tilde K_{10,0}^\dagger \;+\; \tilde K_{10,1}\,\rho^{(10)}\,\tilde K_{10,1}^\dagger.
$$

**(c1)–(c4) 蛍光 reverse 4 命令**（site 3, 2, 1, 0 の順）

$$
\rho^{(12)} \;=\; \tilde K_{10,0}\,\rho^{(11)}\,\tilde K_{10,0}^\dagger \;+\; \tilde K_{10,1}\,\rho^{(11)}\,\tilde K_{10,1}^\dagger,
$$

$$
\rho^{(13)} \;=\; \tilde K_{9,0}\,\rho^{(12)}\,\tilde K_{9,0}^\dagger \;+\; \tilde K_{9,1}\,\rho^{(12)}\,\tilde K_{9,1}^\dagger,
$$

$$
\rho^{(14)} \;=\; \tilde K_{8,0}\,\rho^{(13)}\,\tilde K_{8,0}^\dagger \;+\; \tilde K_{8,1}\,\rho^{(13)}\,\tilde K_{8,1}^\dagger,
$$

$$
\rho^{(15)} \;=\; \tilde K_{7,0}\,\rho^{(14)}\,\tilde K_{7,0}^\dagger \;+\; \tilde K_{7,1}\,\rho^{(14)}\,\tilde K_{7,1}^\dagger.
$$

**(c5)–(c10) TTA reverse 6 命令**（α=6, 5, 4, 3, 2, 1 の順）

$$
\rho^{(16)} \;=\; \tilde K_{6,0}\,\rho^{(15)}\,\tilde K_{6,0}^\dagger \;+\; \tilde K_{6,1}\,\rho^{(15)}\,\tilde K_{6,1}^\dagger,
$$

$$
\rho^{(17)} \;=\; \tilde K_{5,0}\,\rho^{(16)}\,\tilde K_{5,0}^\dagger \;+\; \tilde K_{5,1}\,\rho^{(16)}\,\tilde K_{5,1}^\dagger,
$$

$$
\rho^{(18)} \;=\; \tilde K_{4,0}\,\rho^{(17)}\,\tilde K_{4,0}^\dagger \;+\; \tilde K_{4,1}\,\rho^{(17)}\,\tilde K_{4,1}^\dagger,
$$

$$
\rho^{(19)} \;=\; \tilde K_{3,0}\,\rho^{(18)}\,\tilde K_{3,0}^\dagger \;+\; \tilde K_{3,1}\,\rho^{(18)}\,\tilde K_{3,1}^\dagger,
$$

$$
\rho^{(20)} \;=\; \tilde K_{2,0}\,\rho^{(19)}\,\tilde K_{2,0}^\dagger \;+\; \tilde K_{2,1}\,\rho^{(19)}\,\tilde K_{2,1}^\dagger,
$$

$$
\rho^{(21)} \;=\; \tilde K_{1,0}\,\rho^{(20)}\,\tilde K_{1,0}^\dagger \;+\; \tilde K_{1,1}\,\rho^{(20)}\,\tilde K_{1,1}^\dagger.
$$

**(d) cu_multi 後段**

$$
\rho^{(22)} \;=\; U_H\,\rho^{(21)}\,U_H^\dagger.
$$

最終：

$$
\boxed{\;\rho_{n+1} \;:=\; \rho^{(22)}.\;}
$$

各 $\tilde K_{\alpha,\beta}$ は §8 で定義した $81\times 81$ 行列。各 Kraus 和は **2 項のみ**（$\beta=0,1$）であり、Choi 固有分解の数値抽出ではなく **§6, §7 の解析形をそのまま使う**（実コードは Choi 経由で同じ Kraus を取り出す；解析形と数値形は一致するはず、`KrausChannel.__init__` の CPTP 検証を $10^{-9}$ で通る）。

---

## 11. 1 ステップ時間発展超演算子 $\Phi(dt)$

§10 の 22 行を超演算子の合成として書くと：

$$
\Phi(dt) \;=\; \mathcal{U}_H\,\circ\,\mathcal{E}_1^{(R)}\circ\mathcal{E}_2^{(R)}\circ\cdots\circ\mathcal{E}_{10}^{(R)}\,\circ\,\mathcal{E}_{10}^{(L)}\circ\mathcal{E}_9^{(L)}\circ\cdots\circ\mathcal{E}_1^{(L)}\,\circ\,\mathcal{U}_H,
$$

ただし $\circ$ は超演算子の合成（右から先に作用）、$\mathcal{U}_H(\rho) = U_H\rho U_H^\dagger$、$\mathcal{E}_\alpha^{(L)}(\rho) = \mathcal{E}_\alpha^{(R)}(\rho) = \sum_\beta \tilde K_{\alpha,\beta}\rho\tilde K_{\alpha,\beta}^\dagger$（forward と reverse は同じ局所超演算子；位置だけが違う）。記号上の対応を §10 の命令番号と整合させると、forward 列 (b1)→(b10) は $\mathcal{E}_1^{(L)}\to\mathcal{E}_{10}^{(L)}$ に対応（最初に $\mathcal{E}_1^{(L)}$ が作用）、reverse 列 (c1)→(c10) は $\mathcal{E}_{10}^{(R)}\to\mathcal{E}_1^{(R)}$ に対応（最後に $\mathcal{E}_1^{(R)}$ が作用）。

### 11.1 $\Phi(dt)$ の vec 表現（$D^2\times D^2 = 6561\times 6561$）

vec 化（列優先）の標準恒等式 $\mathrm{vec}(A\rho B) = (B^\top\otimes A)\mathrm{vec}(\rho)$ より：

$$
\Phi(dt) \;\stackrel{\mathrm{vec}}{=}\; (U_H^*\otimes U_H)\,\bigg[\prod_{\alpha=1}^{10}\sum_\beta \tilde K_{\alpha,\beta}^*\otimes \tilde K_{\alpha,\beta}\bigg]\,\bigg[\prod_{\alpha=10}^{1}\sum_\beta \tilde K_{\alpha,\beta}^*\otimes \tilde K_{\alpha,\beta}\bigg]\,(U_H^*\otimes U_H).
$$

ここで $\prod_{\alpha=1}^{10}$ は左から $\alpha=10, 9, \ldots, 1$ の順に行列積（つまり最右の行列が $\alpha=1$、これが最初に作用、本文 §10 の (c10)→(c1) 対応）。**コードはこの 6561×6561 行列を陽には作らない**（メモリ・計算量のため）；$\rho$ を $(3,3,3,3,3,3,3,3)$ テンソルとして保持し、`np.tensordot` で各因子を順次適用する。これは数学的に厳密に同値で、近似なし。

---

## 12. 100 ステップの時間発展

$$
\rho_n \;=\; \Phi(dt)^n(\rho_0),\qquad t_n = n\cdot dt = n.
$$

特に最終時刻：

$$
\hat\rho_\mathrm{final} \;=\; \rho_{100} \;=\; \underbrace{\Phi(dt)\circ\Phi(dt)\circ\cdots\circ\Phi(dt)}_{100\text{ 回}}(\rho_0).
$$

これと真の解 $\hat\rho^\mathrm{true}(100) = e^{100\mathcal{L}}(\rho_0)$ との関係：

$$
\Phi(dt)^{100} \;=\; \big(e^{\mathcal{L}\,dt} + O(dt^3)\big)^{100} \;=\; e^{100\mathcal{L}} + O(t_\mathrm{max}\cdot dt^2) \quad (dt=1, t_\mathrm{max}=100).
$$

ステップ毎の誤差源は **2 つだけ**：

1. **Strang H–D 交換子誤差（ステップ毎 $O(dt^3)$）**：

$$
e^{(\mathcal{L}_H+\mathcal{L}_D)\,dt} \;=\; e^{\mathcal{L}_H\,dt/2}\,e^{\mathcal{L}_D\,dt}\,e^{\mathcal{L}_H\,dt/2} \;+\; \frac{dt^3}{24}\Big([\mathcal{L}_H,[\mathcal{L}_H,\mathcal{L}_D]] - 2[\mathcal{L}_D,[\mathcal{L}_D,\mathcal{L}_H]]\Big) + O(dt^5).
$$

2. **散逸子間のパリンドロミック交換子誤差（ステップ毎 $O(dt^3)$）**：forward と reverse で BCH の 1 次交換子 $\frac12\sum_{\alpha<\beta}[\mathcal{D}_\alpha,\mathcal{D}_\beta]$ が符号反転で相殺。残るのは 3 次の二重交換子。

それ以外はすべて厳密：vec 同型、`scipy.linalg.expm`（$U_H$）、§6 の TTA 解析 Kraus、§7 の蛍光解析 Kraus、`apply_unitary_to_density`／`apply_kraus_to_density` の `tensordot` 縮約。

### 12.1 縮約版での観測量

* 集団：$N_{\mathrm S_0}, N_{\mathrm T_1}, N_{\mathrm S_1}$（§本文 §9.1 と同じ式）。$N_{\mathrm S_0}+N_{\mathrm T_1}+N_{\mathrm S_1}=N=4$ がトレース 1 から従う。
* von Neumann エントロピー、純度、トレース保存性のモニタも本文 §9.2–§9.4 と同じ。

### 12.2 縮約版で起こる物理（数式の予測）

10 個のチャンネルだけなので、長時間極限での挙動は本文の 26 チャンネル版と異なる：

* TTA は $|11\rangle$（隣接 2 サイトが共に $\mathrm T_1$）でだけ発火する。1 回発火すると一方が $|2\rangle$、他方が $|0\rangle$ になる。$|2\rangle$ になった分子は蛍光（$|2\rangle\to|0\rangle$）で速やかに $|0\rangle$ に落ちる。
* 蛍光は $|2\rangle$ 集団（$N_{\mathrm S_1}$）を時定数 $1/\Gamma_\mathrm{fl}$ で 0 に向けて減衰させる。
* **$\mathrm T_1$ には他の脱励起チャンネルがない**（燐光・IC・ISC$_{T\to S}$ を除外したため）。よって TTA で消滅しなかった $\mathrm T_1$ は無限に残る。Hamiltonian の移動項 $V$ で $\mathrm T_1$ は鎖上を動くが、終状態には必ず $\mathrm T_1$ がいくつか残る。
* 初期状態 `edge_triplet` $|1001\rangle$（site 0 と site 3 が $\mathrm T_1$）から始めると、site 0 と site 3 の $\mathrm T_1$ は移動して中央の site 1 や site 2 と隣接対 $|11\rangle$ を作ることがあり、そのときだけ TTA が発火する。発火頻度はパラメータ次第。
* 燐光・ISC を入れていないので $\mathrm T_1$ の長時間平均集団は 26 チャンネル版より大きく残る（数値的検証は省略；本文書の主旨は「式を完全に書く」ことであり、数値実測は別途必要）。

これらの物理予測はモデルの式（§3.3）から読み取れるものであり、観測量との比較は別途実装が必要。

---

## 13. 実コードとの差分（嘘禁止のため明示）

| 項目 | 本縮約版 | 実コード（シナリオ5d そのもの） |
|---|---|---|
| Lindblad 演算子の数 | 10 | 26 |
| 含まれるチャンネル | TTA + 蛍光 | TTA + 蛍光 + 燐光 + IC + ISC$_{S\to T}$ + ISC$_{T\to S}$ |
| 1 step の `kraus_channel` 命令数 | 20 (forward 10 + reverse 10) | 52 (forward 26 + reverse 26) |
| 1 step の `cu_multi` 命令数 | 2 | 2 |
| 1 step の合計命令数 | **22** | **54** |
| Hamiltonian | $\hat H_0 + \hat H_\mathrm{transfer}$（$E_T, E_S, V$ 全て使用） | 同左 |
| Kraus 抽出方法 | 解析形（rank 2、§6, §7） | Choi–Jamiolkowski の数値固有分解（rank ≤ $d_\mathrm{loc}^2$、CHOI_EIG_TOL=1e-12 で切り捨て） |
| Kraus の数値値 | $p_1 = 1-e^{-\gamma_\mathrm{TTA}\,dt/4}$、$p_\mathrm{fl} = 1-e^{-\Gamma_\mathrm{fl}\,dt/2}$ から閉形式 | TTA, 蛍光については解析形と数値一致するはず（CPTP 検証を $10^{-9}$ で通る） |
| 物理結果 | $\mathrm T_1$ は TTA でしか消えないので長時間に残る | 燐光・ISC$_{T\to S}$ で $\mathrm T_1$ も時定数 $\sim 1/(\Gamma_\mathrm{ph}+k_{\mathrm{ISC},T\to S})$ で減衰 |

本縮約版を実コード上で「動かす」には、`tutorials/gksl_math_utils.py::build_lindblad_operators` から燐光・IC・ISC$_{S\to T}$・ISC$_{T\to S}$ の項を取り除いた変種関数が必要だが、本文書では数式の整理に留め、コード変更は行っていない。

---

## 14. まとめ

本文書では、**TTA と 蛍光のみ**に Lindblad 演算子を限定した GKSL 方程式

$$
\dot{\hat\rho} \;=\; -i\,[\hat H_\mathrm{total},\hat\rho] \;+\; \sum_{\alpha=1}^{10} \mathcal{D}[\hat L_\alpha](\hat\rho)
$$

について、

* 10 個の $\hat L_\alpha$ をすべて $3\times 3$／$9\times 9$ 局所部分のテンソル積として完全に書き下し（§3）、
* GKSL 方程式の 10 個の散逸子項をすべて展開し（§3.3）、
* 各局所散逸子の半ステップ Kraus を**解析形**（§6 の TTA で $K_0,K_1$、§7 の蛍光で $K_0,K_1$、いずれも rank 2）として閉形式で導出し（CPTP 完全性も解析的に証明、§6.1, §7）、
* 1 Trotter ステップの 22 命令の `cu_multi`／`kraus_channel` 列を全部書き下し（§9）、
* それら 22 命令を $\rho$ に作用させる 22 行の代入式を一切省略せずに書き出し（§10）、
* 1 ステップ時間発展超演算子 $\Phi(dt)$ の表式（演算子合成形と vec 行列形）（§11）、
* 100 ステップの時間発展 $\rho_{100} = \Phi(dt)^{100}(\rho_0) = e^{100\mathcal L}\rho_0 + O(t_\mathrm{max}\cdot dt^2)$（§12）、
* 実コードとの差分の正直な明示（§13）

を示した。残存する近似は Strang H–D 分割とパリンドロミック散逸子分割（合算で大域 $O(dt^2)$）の 2 つだけで、それ以外は全て厳密。
