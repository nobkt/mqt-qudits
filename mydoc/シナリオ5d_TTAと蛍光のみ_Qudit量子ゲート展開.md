# シナリオ5d 縮約版（TTA + 蛍光のみ）：Qutrit 量子ゲート完全展開

## 0. 本文書の位置付けと真実ベースの宣言

前文書 `mydoc/シナリオ5d_TTAと蛍光のみ_完全展開.md` は **超演算子レベル** で 22 命令（`cu_multi`×2 + `kraus_channel`×20）を完全に展開したが、これらは MQT-Qudits の **黒箱命令** であり、

* `cu_multi(target=[0,1,2,3], parameters=U_H)` は **$81\times 81$ 行列 $U_H$ をそのまま渡す**だけで、ネイティブ qutrit ゲート（`R`, `VirtRz`, `CSum` 等）への分解は持っていない。
* `kraus_channel(target, kraus_operators)` も **$d_\mathrm{loc}\times d_\mathrm{loc}$ Kraus 行列をそのまま渡す**だけで、ネイティブ qutrit ゲートへの分解は持っていない。

**事実（嘘禁止のための明示）**：

1. **実コード `tutorials/qudit_gksl_simulator.py::_trotter_step_dmsim` は `cu_multi` と `kraus_channel` を黒箱として並べるだけで、ネイティブ qutrit ゲート分解を一切行わない**。  
   `circuit.cu_multi(list(range(n)), self._U_H_half.astype(np.complex128))`（行 296, 309）、`circuit.kraus_channel(target, kraus_ops)`（行 301, 306）の通り、$U_H$ も Kraus も `numpy` 配列をそのまま命令に詰める。
2. **実コードの `dmsim` バックエンド（`src/mqt/qudits/simulation/backends/dmsim.py`）も同様に、`cu_multi` の行列で $\rho \mapsto U\rho U^\dagger$ を、`kraus_channel` の Kraus で $\rho \mapsto \sum_k K_k\rho K_k^\dagger$ を `np.tensordot` で適用するだけ**で、ネイティブ qutrit ゲート分解は通らない。
3. **`stinespring` algorithm パス（`tutorials/stinespring_utils.py::stinespring_unitary_from_lindblad`）は確かに Stinespring ダイレーション $U_\alpha=\exp(-i\sqrt{dt}\,G_\alpha)$ を作るが、これも $d_\mathrm{anc}\cdot d_\mathrm{sys}$ 次元の行列としてそのまま `apply_stinespring_to_density_matrix` で適用される**。やはりネイティブ qutrit ゲート分解は無い。

したがって、**「時間発展演算子の量子ゲート表現」は、現コードには存在しない**。本文書はそれを**理論的に補完する**ものであり、本文書の §2–§5 で導出する qudit ゲート列は **実コード上で動かす版ではなく、ハードウェアで動かす場合に必要な解析的な分解**である。この区別を §7 で表にして再掲する。

なお本文書は **TTA + 蛍光のみの 10 Lindblad** という縮約モデルを前提とする（前文書 §0 の通り、実コードは 26 Lindblad）。前文書の §2–§7 を前提にし、本文書では **ネイティブ qutrit ゲート分解の層**だけを追加する。

物理パラメータ（`tutorials/gksl_physical_parameters.py`）：$N=4$ 分子、$d=3$、$D=81$、$dt=1$、$t_\mathrm{max}=100$、隣接 $\{(0,1),(1,2),(2,3)\}$。

---

## 1. MQT-Qudits ネイティブ qutrit ゲートライブラリ（実コードから直接読んだ定義）

`src/mqt/qudits/quantum_circuit/gates/` から、本文書で使う 4 つのゲートを正確に書き下す。$d=3$（qutrit）の場合のみ示す。

### 1.1 単一 qutrit `VirtRz`（仮想位相、`virt_rz.py`）

`VirtRz(target=i, parameters=[lev, φ], dimensions=3)` の作用（`__array__`）：

$$
\mathrm{VirtRz}(\mathrm{lev},\varphi) \;=\; \mathrm{diag}_3\big(1,1,1\big)\quad\text{ただし}\quad (\mathrm{lev},\mathrm{lev})\text{ 成分を }e^{-i\varphi}\text{ に置換}.
$$

具体例：

$$
\mathrm{VirtRz}(\mathrm{lev}{=}1,\varphi) \;=\; \begin{pmatrix} 1 & 0 & 0 \\ 0 & e^{-i\varphi} & 0 \\ 0 & 0 & 1 \end{pmatrix},\qquad \mathrm{VirtRz}(\mathrm{lev}{=}2,\varphi) \;=\; \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & e^{-i\varphi} \end{pmatrix}.
$$

### 1.2 単一 qutrit `R`（Givens 回転、`r.py`）

`R(target=i, parameters=[lev_a, lev_b, θ, φ], dimensions=3)` の作用は、サブ空間 $\{|a\rangle,|b\rangle\}$ に限定された $SU(2)$ 回転（残り 1 準位は不変）：

$$
R(a,b;\theta,\varphi) \;=\; \exp\!\Big[-i\,\tfrac{\theta}{2}\big(\cos\varphi\;\Lambda_s^{(a,b)} + \sin\varphi\;\Lambda_a^{(a,b)}\big)\Big],
$$

ここで $\Lambda_s^{(a,b)} = |a\rangle\langle b| + |b\rangle\langle a|$、$\Lambda_a^{(a,b)} = -i|a\rangle\langle b| + i|b\rangle\langle a|$（`gellmann.py` の "s" / "a" の定義）。$\{a,b\}$ サブ空間に制限した $2\times 2$ 表現は通常の Pauli 表記で

$$
R^{2D}(a,b;\theta,\varphi) \;=\; \cos\!\tfrac{\theta}{2}\;I_2 \;-\; i\sin\!\tfrac{\theta}{2}\big(\cos\varphi\;\sigma_x^{ab} + \sin\varphi\;\sigma_y^{ab}\big).
$$

3 番目の準位 $c\in\{0,1,2\}\setminus\{a,b\}$ では恒等。たとえば $a=0,b=1$、$\varphi=0$ のとき

$$
R(0,1;\theta,0) \;=\; \begin{pmatrix} \cos(\theta/2) & -i\sin(\theta/2) & 0 \\ -i\sin(\theta/2) & \cos(\theta/2) & 0 \\ 0 & 0 & 1 \end{pmatrix}.
$$

### 1.3 2-qutrit `CSum`（`csum.py`）

$d=3$ の場合（`__array__` を $d=3$ で評価）。$X_3$ は $|i\rangle\to|i+1\bmod 3\rangle$ のシフト（`x.py`）：

$$
X_3 \;=\; \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix},\qquad X_3^0=I_3,\quad X_3^1=X_3,\quad X_3^2=X_3^2.
$$

`CSum(target=[c, t])`（$c<t$ の通常順）は

$$
\mathrm{CSum} \;=\; \sum_{i=0}^{2} |i\rangle_c\langle i| \otimes X_3^{\,i} \;\in\;\mathbb C^{9\times 9}.
$$

明示的に $9\times 9$ 行列で書くと（行/列は $|s_c s_t\rangle$ で 0..8、index = $3s_c+s_t$）：

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
\end{array}\right).
$$

ブロック $i$ は $X_3^{\,i}$。

### 1.4 `KrausChannel` 命令（`kraus_channel.py`）

§0 の通り、これはネイティブ ユニタリ ゲートではなく、CPTP 命令（$\sum_k K_k^\dagger K_k = I$ を許容差 $10^{-9}$ で実測検証、`KrausChannel.__init__` が ValueError を投げる）。**ハードウェア実装には Stinespring ダイレーション → ユニタリ + 補助消去 が必要**。本文書 §5 でそれを行う。

---

## 2. $\hat H_0$ の半ステップ $e^{-i\hat H_0\,dt/2}$ の **厳密** qutrit ゲート分解

### 2.1 構造

前文書 §2.2 より $\hat H_0 = \sum_{i=0}^{3} h_\mathrm{loc}^{(i)}$、$h_\mathrm{loc} = \mathrm{diag}(0, E_T, E_S)$。各 $h_\mathrm{loc}^{(i)}$ は **互いに可換**（異なるサイトで作用）かつ **対角**。よって

$$
e^{-i\hat H_0\,dt/2} \;=\; \prod_{i=0}^{3} e^{-i\,h_\mathrm{loc}^{(i)}\,dt/2} \;=\; \bigotimes_{i=0}^{3} e^{-i\,h_\mathrm{loc}\,dt/2}.
$$

ここで等号は **厳密**（Trotter 誤差なし、$[h_\mathrm{loc}^{(i)},h_\mathrm{loc}^{(j)}]=0$）。各サイトの $3\times 3$ 行列：

$$
e^{-i\,h_\mathrm{loc}\,dt/2} \;=\; \mathrm{diag}\!\Big(1,\;e^{-iE_T\,dt/2},\;e^{-iE_S\,dt/2}\Big).
$$

### 2.2 VirtRz への分解（厳密、各サイト 2 ゲート）

$\mathrm{diag}(1, e^{-iE_T dt/2}, e^{-iE_S dt/2})$ は §1.1 の VirtRz の積で**厳密**に書ける：

$$
\mathrm{diag}(1, e^{-iE_T dt/2}, 1)\cdot\mathrm{diag}(1, 1, e^{-iE_S dt/2}) \;=\; \mathrm{VirtRz}(1,\,E_T\,dt/2)\cdot\mathrm{VirtRz}(2,\,E_S\,dt/2).
$$

両者は対角ゆえ可換（順序は任意）。

### 2.3 $\hat H_0$ 半ステップの完全 qutrit ゲート列（**8 ゲート、誤差なし**）

$$
\boxed{\;e^{-i\hat H_0\,dt/2} \;=\; \prod_{i=0}^{3}\Big[\mathrm{VirtRz}^{(i)}(1,\,E_T\,dt/2)\cdot \mathrm{VirtRz}^{(i)}(2,\,E_S\,dt/2)\Big]\;}\quad(\text{8 ゲート、可換、厳密})
$$

すなわち以下の 8 命令を任意の順序で発行：

```
VirtRz(target=0, lev=1, phi=E_T·dt/2)
VirtRz(target=0, lev=2, phi=E_S·dt/2)
VirtRz(target=1, lev=1, phi=E_T·dt/2)
VirtRz(target=1, lev=2, phi=E_S·dt/2)
VirtRz(target=2, lev=1, phi=E_T·dt/2)
VirtRz(target=2, lev=2, phi=E_S·dt/2)
VirtRz(target=3, lev=1, phi=E_T·dt/2)
VirtRz(target=3, lev=2, phi=E_S·dt/2)
```

これは前文書 §2.2 の $\hat H_0$（4 項 Kronecker 和）に対する **正確な** ゲート分解で、近似は一切無い。

---

## 3. $\hat H_\mathrm{transfer}$ の半ステップ $e^{-i\hat H_\mathrm{transfer}\,dt/2}$ の qutrit ゲート分解

### 3.1 各ペアごとの局所構造

前文書 §2.3 より

$$
\hat H_\mathrm{transfer} \;=\; \sum_{(i,j)\in\{(0,1),(1,2),(2,3)\}}\hat h^{(i,j)},\qquad \hat h^{(i,j)} \;=\; V\big(\hat A^{(i,j)} + \hat A^{(i,j)\,\dagger}\big),
$$

$\hat A^{(i,j)} = |0\rangle_i\langle 1|\otimes|1\rangle_j\langle 0| = |01\rangle_{ij}\langle 10|$（添字順 $i$ 第一、$j$ 第二）。$\hat h^{(i,j)}$ の局所部分（$9\times 9$、$|s_is_j\rangle$ 基底、index $3s_i+s_j$）：

$$
h_\mathrm{loc}^\mathrm{tr} \;=\; V\big(|01\rangle\langle 10| + |10\rangle\langle 01|\big) \;=\; V\,E^{(9)}_{1,3} + V\,E^{(9)}_{3,1}.
$$

これは行/列 index $\{1,3\}$ に閉じる $2\times 2$ ブロック $V\sigma_x$、それ以外は 0：

$$
h_\mathrm{loc}^\mathrm{tr}\Big|_{\{|01\rangle,|10\rangle\}} \;=\; V\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},\qquad h_\mathrm{loc}^\mathrm{tr}\Big|_\mathrm{else} \;=\; 0.
$$

### 3.2 ペア内指数の閉形式（厳密）

$h_\mathrm{loc}^\mathrm{tr}$ は $\{|01\rangle,|10\rangle\}$ 部分空間に閉じるので、$e^{-i h_\mathrm{loc}^\mathrm{tr}\tau}$ も同じ部分空間で $2\times 2$ 回転、補空間（$\{|00\rangle,|02\rangle,|11\rangle,|12\rangle,|20\rangle,|21\rangle,|22\rangle\}$、index $\{0,2,4,5,6,7,8\}$）では恒等：

$$
\big[e^{-i h_\mathrm{loc}^\mathrm{tr}\tau}\big]\Big|_{\{|01\rangle,|10\rangle\}} \;=\; \cos(V\tau)\,I_2 - i\sin(V\tau)\,\sigma_x \;=\; \begin{pmatrix} \cos(V\tau) & -i\sin(V\tau) \\ -i\sin(V\tau) & \cos(V\tau) \end{pmatrix}.
$$

これと §1.2 の $R^{2D}(a,b;\theta=2V\tau,\varphi=0)$ の定義式 $\cos(\theta/2)I - i\sin(\theta/2)\sigma_x$ を比較すると、**$\theta = 2V\tau$、$\varphi=0$ で完全一致**。

### 3.3 ペア間 Trotter 分割（誤差発生）

3 つのペア $(0,1), (1,2), (2,3)$ の局所演算子 $\hat h^{(0,1)}, \hat h^{(1,2)}, \hat h^{(2,3)}$ は **可換ではない**（site 1 と site 2 を共有するため）：

* $[\hat h^{(0,1)}, \hat h^{(1,2)}] \neq 0$（site 1 で $|0\rangle\langle 1|$ と $|1\rangle\langle 0|$ が交わる）
* $[\hat h^{(1,2)}, \hat h^{(2,3)}] \neq 0$（site 2 で同様）
* $[\hat h^{(0,1)}, \hat h^{(2,3)}] = 0$（共通サイト無し）

したがって $e^{-i\hat H_\mathrm{transfer}\,dt/2} = \prod_{(i,j)} e^{-i\hat h^{(i,j)}\,dt/2}$ は **厳密ではなく Trotter 近似**となる。Strang 対称化で $O((dt/2)^3)$ 誤差にできるが、それでも誤差は残る。

#### 3.3.1 1 次（非対称）Trotter（誤差 $O((dt/2)^2)$ ペア毎）

$$
e^{-i\hat H_\mathrm{transfer}\,dt/2}\;\overset{\mathrm{Lie}}{\approx}\; e^{-i\hat h^{(0,1)}\,dt/2}\,e^{-i\hat h^{(1,2)}\,dt/2}\,e^{-i\hat h^{(2,3)}\,dt/2}.
$$

#### 3.3.2 Strang（対称）Trotter（誤差 $O((dt/2)^3)$ 半ステップ毎）

$$
e^{-i\hat H_\mathrm{transfer}\,dt/2}\;\overset{\mathrm{Strang}}{\approx}\; e^{-i\hat h^{(0,1)}\,dt/4}\,e^{-i\hat h^{(1,2)}\,dt/4}\,e^{-i\hat h^{(2,3)}\,dt/2}\,e^{-i\hat h^{(1,2)}\,dt/4}\,e^{-i\hat h^{(0,1)}\,dt/4}.
$$

### 3.4 各 $e^{-i\hat h^{(i,j)}\,\tau}$ のネイティブ qutrit ゲート分解

§3.2 の $9\times 9$ 行列 $E_\mathrm{loc}^\mathrm{tr}(\tau) := e^{-i h_\mathrm{loc}^\mathrm{tr}\,\tau}$ は

$$
E_\mathrm{loc}^\mathrm{tr}(\tau)_{(s_is_j),(s_i's_j')} \;=\; \begin{cases} \cos(V\tau) & (s_is_j) = (s_i's_j') \in\{(0,1),(1,0)\} \\ -i\sin(V\tau) & (s_is_j,s_i's_j')\in\{((0,1),(1,0)),((1,0),(0,1))\} \\ 1 & (s_is_j) = (s_i's_j')\in\{(0,0),(0,2),(1,1),(1,2),(2,0),(2,1),(2,2)\} \\ 0 & \text{それ以外} \end{cases}
$$

これを **MQT-Qudits の native ゲートだけ**で実装する。$\{|01\rangle,|10\rangle\}$ という 2 元集合は 2-qutrit の **異なるレベル間** にまたがる遷移なので、単一 qutrit の `R` だけでは作れない。**必ずエンタングリングゲート（`CSum` または `MS`）が必要**。

#### 3.4.1 構成方針：CSum で `|01⟩↔|10⟩` を `|01⟩↔|11⟩` のような単一 qutrit 遷移にマップ

$\mathrm{CSum}$（§1.3）の作用は $\mathrm{CSum}|s_c s_t\rangle = |s_c, s_t+s_c\bmod 3\rangle$。$c=i$（制御）、$t=j$（標的）として：

* $\mathrm{CSum}|01\rangle = |01\rangle$（$s_c=0$ なのでシフトなし）
* $\mathrm{CSum}|10\rangle = |11\rangle$（$s_c=1$ なので $s_t$ を $+1$）

つまり **$\mathrm{CSum}$ は $\{|01\rangle,|10\rangle\}$ を $\{|01\rangle,|11\rangle\}$ に写す**。後者は **site $i$ の $|0\rangle\leftrightarrow|1\rangle$ 遷移**を「site $j$ が $|1\rangle$ のとき」だけ実行することと等価。

しかし control-on-$|1\rangle$ の selective rotation は MQT-Qudits ライブラリには直接の単一ゲートが無い（`X` は $\bmod d$ のシフトで control 機能はない）。最も簡潔な分解は次のサンドイッチ：

$$
e^{-i\hat h^{(i,j)}\,\tau} \;=\; \mathrm{CSum}^{(i,j)}\;\cdot\;U_\mathrm{mid}(\tau)\;\cdot\;\big(\mathrm{CSum}^{(i,j)}\big)^\dagger,
$$

ここで $U_\mathrm{mid}(\tau)$ は $\{|01\rangle,|11\rangle\}$ に閉じる $2\times 2$ 回転 $\cos(V\tau)I - i\sin(V\tau)\sigma_x^{\{|01\rangle,|11\rangle\}}$、補空間で恒等。$\{|01\rangle,|11\rangle\}$ は **site $i$ が $\{|0\rangle,|1\rangle\}$、site $j$ は $|1\rangle$ で固定**であり、site $j$ が $|1\rangle$ のときだけ site $i$ の levels $\{0,1\}$ を回転させる。これは "$1$-controlled $R(0,1)$ on site $i$" であり、これも 1 つの native gate ではない。

**正直な評価**：$\hat A + \hat A^\dagger = |01\rangle\langle 10| + |10\rangle\langle 01|$ は 2-qutrit Hilbert 空間の **非対角ブロック** に作用する非自明な双 qutrit 演算子であり、MQT-Qudits の **`R`、`VirtRz`、`CSum`、`X`、`Z`、`H` の有限積では一般には書けない**。書けるのは `CustomTwo` や `MS`/`LS` のような「2-qutrit 連続パラメータ族」を含めた場合のみ。最も簡潔で**実装可能**な書き方は次のいずれか：

* (A) **`CustomTwo`** で $9\times 9$ 行列 $E_\mathrm{loc}^\mathrm{tr}(\tau)$ をそのまま渡す（実コードの cu_multi に最も近い；ハードウェアコンパイラがこれを native 列に分解する）。
* (B) **`MS`/`LS`** ゲート（イオントラップ向け 2-qutrit ゲート、$\exp(-i\theta\,J_x^2/4)$ 型）を組み合わせ、その前後に単一 qutrit `R`/`VirtRz` を挟む数式変換。
* (C) **2 個の `CSum`（または逆向き 1 個）と 数個の単一 qutrit `R` の組合せ**で「制御つき $\{0,1\}$ 回転」を構成。

ここでは方針 (A) が最も**正直**である。なぜなら、MQT-Qudits の現在のコンパイラ・ターゲット非依存性を保ったまま、$U_\mathrm{H}$ の 1 ペア分を確実に表現できるからである。明示形は次の通り：

#### 3.4.2 ペアごとの `CustomTwo`（**実装可能、qudit 命令に確かにマップする**）

サポート $S=\{i,j\}$、命令 `CustomTwo(target=[i,j], parameters=E_loc^tr(τ), dimensions=[3,3])`。

行列 $E_\mathrm{loc}^\mathrm{tr}(\tau)$ は §3.2 で全 81 成分のうち非零 11 成分を完全に書き下した（対角 9 個と非対角 2 個）。**これは前文書 §2.3 の $\hat h^{(i,j)}$ の指数を、近似なしに 2-qutrit 命令にしたもの**。

#### 3.4.3 方針 (C) の解析的分解（参考、`CSum` + `R` のみで構成）

ハードウェアが `CustomTwo` を直接サポートしない場合、次の分解で **`CSum` 2 個 + 単一 qutrit `R` 4 個** に落ちる（方針 (C)）。

ステップ (i)：$\mathrm{CSum}^{(i,j)}$ を適用 → $\{|01\rangle,|10\rangle\}$ は $\{|01\rangle,|11\rangle\}$ にマップされ、これは site $j=|1\rangle$ 部分空間に閉じる。

ステップ (ii)：site $j$ の $|1\rangle$ 部分空間に制限すれば site $i$ の $\{|0\rangle,|1\rangle\}$ 上の $R(0,1;2V\tau,0)$ になるが、**「site $j=|1\rangle$ のときだけ実行」を native ゲートで作るには更に 2-qutrit ゲートが必要**。

正直に言うと：「**1 個の $\mathrm{CSum}$ + 単一 qutrit $R$ + 1 個の $\mathrm{CSum}^\dagger$** で `selective rotation in {|01⟩,|10⟩}` を作る厳密な分解は、qubit ($d=2$) では存在するが、qutrit ($d=3$) では追加のフェーズキャンセル（`VirtRz`）が必要で、しかも一般には Solovay-Kitaev 型の **多段** 分解になる**。MQT-Qudits の `compiler.compilation_minitools` モジュールがこの分解を担うが、本文書 §0 で明示した通り、tutorial レベルの `qudit_gksl_simulator.py` はこのコンパイラを呼ばずに **`cu_multi(U_H)` で済ませている**。

**結論（嘘禁止）**：方針 (A)（`CustomTwo` または `cu_multi`）が **実コードと整合する唯一の正直な選択**。方針 (C) の完全展開は本文書のスコープ外（コンパイラ通過後の native 列はハードウェア依存）。

### 3.5 $\hat H_\mathrm{transfer}$ 半ステップの完全 qutrit ゲート列

方針 (A) を採用、Strang 対称化（§3.3.2）で書くと：

$$
\boxed{\;e^{-i\hat H_\mathrm{transfer}\,dt/2}\;\overset{\mathrm{Strang}}{\approx}\;
\underbrace{\mathrm{CT}^{(0,1)}\!\big(\tfrac{dt}{4}\big)}_\text{1}
\,\underbrace{\mathrm{CT}^{(1,2)}\!\big(\tfrac{dt}{4}\big)}_\text{2}
\,\underbrace{\mathrm{CT}^{(2,3)}\!\big(\tfrac{dt}{2}\big)}_\text{3}
\,\underbrace{\mathrm{CT}^{(1,2)}\!\big(\tfrac{dt}{4}\big)}_\text{4}
\,\underbrace{\mathrm{CT}^{(0,1)}\!\big(\tfrac{dt}{4}\big)}_\text{5}\;}
$$

ここで $\mathrm{CT}^{(i,j)}(\tau) := \texttt{CustomTwo(target=[i,j], parameters=}E_\mathrm{loc}^\mathrm{tr}(\tau)\texttt{, dimensions=[3,3])}$、$E_\mathrm{loc}^\mathrm{tr}(\tau)$ は §3.2 の $9\times 9$ 行列。誤差 $O((dt/2)^3)$（half-step 当たり）。

**※実コードはこの分解を使わず、$U_H = e^{-i(\hat H_0+\hat H_\mathrm{transfer})\,dt/2}$ を `scipy.linalg.expm` で 1 個の $81\times 81$ 行列にしてから `cu_multi` に渡す**。実コードでは $\hat H_0$ と $\hat H_\mathrm{transfer}$ も Trotter 分割していない、つまり §2 の厳密分解と §3.5 の Strang 近似を**両方使わない**。本文書 §3.5 はハードウェア向け追加分解。

---

## 4. 統合：$U_H$ 全体の完全 qutrit ゲート列

### 4.1 ケース I：実コード踏襲（`cu_multi` 1 個、近似なし、ハードウェア非対応）

$$
U_H = e^{-i(\hat H_0+\hat H_\mathrm{transfer})\,dt/2}
\quad\Longrightarrow\quad
\texttt{CustomMulti(target=[0,1,2,3], parameters=}U_H\texttt{, dimensions=[3,3,3,3])}
$$

**1 命令、$81\times 81$ 黒箱**。これが実コード `qudit_gksl_simulator.py:296,309` の本体。

### 4.2 ケース II：$H_0$ と $H_\mathrm{transfer}$ を分離（ハードウェア向け）

**訂正（補遺 H で証明）**：本節初版は「$[\hat H_0,\hat H_\mathrm{transfer}]\neq 0$」と書いたが、これは **誤り**。$h_\mathrm{loc}=\mathrm{diag}(0,E_T,E_S)$、$\hat A^{(i,j)}=a_ia_j^\dagger$（$|0\rangle\leftrightarrow|1\rangle$ 遷移）の場合、サイト $i,j$ の局所エネルギー差はペア両端で打ち消され $[\hat H_0,\hat H_\mathrm{transfer}]=0$ が成り立つ（補遺 H で完全に演算子代数で証明）。よって Strang 分割は不要、**厳密に分離**：

$$
e^{-i\hat H_\mathrm{total}\,dt/2}\;=\;
e^{-i\hat H_0\,dt/2}\cdot e^{-i\hat H_\mathrm{transfer}\,dt/2}\quad\text{（厳密、補遺 H）}.
$$

§2.3 で $e^{-i\hat H_0\,dt/4}$ は **VirtRz 8 個（厳密）**、§3.5 で $e^{-i\hat H_\mathrm{transfer}\,dt/2}$ は **CustomTwo 5 個（Strang 近似）**。よって $U_H$ 半ステップは

$$
\boxed{\;\underbrace{8\times \mathrm{VirtRz}}_\text{H_0 / dt/4}\;
\underbrace{\mathrm{CT}^{(0,1)}, \mathrm{CT}^{(1,2)}, \mathrm{CT}^{(2,3)}, \mathrm{CT}^{(1,2)}, \mathrm{CT}^{(0,1)}}_\text{H_\mathrm{transfer} / Strang}\;
\underbrace{8\times \mathrm{VirtRz}}_\text{H_0 / dt/4}\;}
\qquad(\text{合計 21 ゲート、誤差 }O((dt)^3))
$$

明示的に列挙：

```
# H_0 quarter-step (前段、厳密、8 ゲート)
VirtRz(target=0, lev=1, phi=E_T·dt/4);  VirtRz(target=0, lev=2, phi=E_S·dt/4)
VirtRz(target=1, lev=1, phi=E_T·dt/4);  VirtRz(target=1, lev=2, phi=E_S·dt/4)
VirtRz(target=2, lev=1, phi=E_T·dt/4);  VirtRz(target=2, lev=2, phi=E_S·dt/4)
VirtRz(target=3, lev=1, phi=E_T·dt/4);  VirtRz(target=3, lev=2, phi=E_S·dt/4)

# H_transfer half-step (Strang、近似、5 ゲート)
CustomTwo(target=[0,1], parameters=E_loc^tr(dt/4), dimensions=[3,3])
CustomTwo(target=[1,2], parameters=E_loc^tr(dt/4), dimensions=[3,3])
CustomTwo(target=[2,3], parameters=E_loc^tr(dt/2), dimensions=[3,3])
CustomTwo(target=[1,2], parameters=E_loc^tr(dt/4), dimensions=[3,3])
CustomTwo(target=[0,1], parameters=E_loc^tr(dt/4), dimensions=[3,3])

# H_0 quarter-step (後段、厳密、8 ゲート)
VirtRz(target=0, lev=1, phi=E_T·dt/4);  VirtRz(target=0, lev=2, phi=E_S·dt/4)
VirtRz(target=1, lev=1, phi=E_T·dt/4);  VirtRz(target=1, lev=2, phi=E_S·dt/4)
VirtRz(target=2, lev=1, phi=E_T·dt/4);  VirtRz(target=2, lev=2, phi=E_S·dt/4)
VirtRz(target=3, lev=1, phi=E_T·dt/4);  VirtRz(target=3, lev=2, phi=E_S·dt/4)
```

各 `CustomTwo` の `parameters` は §3.2 で完全に書き下した $9\times 9$ 行列 $E_\mathrm{loc}^\mathrm{tr}(\tau)$（11 個の非零成分、残り 70 個はゼロ）。

---

## 5. KrausChannel の qutrit ゲート分解（Stinespring 経由）

KrausChannel は非ユニタリ → ハードウェアでは **補助 qutrit + ユニタリ + 補助の射影測定または捨て** で実装する必要がある（Stinespring ダイレーション）。前文書 §6, §7 の rank-2 解析 Kraus に対し、各チャンネルに **1 個の補助 qutrit** で十分。

### 5.1 一般構造

rank-2 Kraus $\{K_0, K_1\}$（$\sum K_k^\dagger K_k = I_{d_\mathrm{loc}}$）に対し、補助 qutrit $|0\rangle_\mathrm{anc}$ から始め、ユニタリ $U_\mathrm{Stine}$ を作用させて

$$
U_\mathrm{Stine}\,(|0\rangle_\mathrm{anc}\otimes|\psi\rangle_\mathrm{sys}) \;=\; |0\rangle_\mathrm{anc}\otimes K_0|\psi\rangle + |1\rangle_\mathrm{anc}\otimes K_1|\psi\rangle,
$$

その後 anc を **partial trace**（捨てる）。$U_\mathrm{Stine}$ の最初の $d_\mathrm{loc}$ 列は $\binom{K_0}{K_1}$ で固定（残り $(d_\mathrm{anc}-1)d_\mathrm{loc}$ 列は等距離拡張で任意に補う）。

具体的に、$U_\mathrm{Stine}\in\mathbb C^{(3 d_\mathrm{loc})\times(3 d_\mathrm{loc})}$ で、**ブロック行列**

$$
U_\mathrm{Stine} \;=\; \begin{pmatrix} K_0 & * & * \\ K_1 & * & * \\ 0 & * & * \end{pmatrix},
$$

の最左列ブロックが $(K_0, K_1, 0)^\top$。残り 2 ブロックは「$U_\mathrm{Stine}^\dagger U_\mathrm{Stine} = I$」を満たすように構成（Gram-Schmidt 等。$K_2 = 0$ という選択が最も自然で、`stinespring_unitary_from_lindblad`（§0 の `tutorials/stinespring_utils.py:14-51`）の "$d_\mathrm{anc}=3$ で qutrit ancilla" と同じ）。

### 5.2 蛍光チャンネルの Stinespring ユニタリ（$d_\mathrm{loc}=3$、$U_\mathrm{Stine}\in\mathbb C^{9\times 9}$）

前文書 §7 の Kraus：$K_0 = \mathrm{diag}(1,1,q_\mathrm{fl})$、$K_1 = \sqrt{p_\mathrm{fl}}\,|0\rangle\langle 2|$（$p_\mathrm{fl} = 1-e^{-\Gamma_\mathrm{fl}\,dt/2}$、$q_\mathrm{fl} = e^{-\Gamma_\mathrm{fl}\,dt/4}$）。

最左 $9\times 3$ 列：

$$
\begin{pmatrix} K_0 \\ K_1 \\ 0 \end{pmatrix}_{9\times 3} \;=\; \begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & q_\mathrm{fl} \\\hline
0 & 0 & \sqrt{p_\mathrm{fl}} \\
0 & 0 & 0 \\
0 & 0 & 0 \\\hline
0 & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}.
$$

これは sys（下 3 levels）と anc（上 3 levels）のテンソル積で、anc が上位 → kron(anc, sys) 規約を採用（前文書 §10、`stinespring_utils` の規約と同じ）。各列はノルム 1（$1,1,q_\mathrm{fl}^2+p_\mathrm{fl}=1$）。

残り 6 列は等距離拡張。最も簡潔な選択：

* 列 4（$|3\rangle$）：$(0,0,-\sqrt{p_\mathrm{fl}},\,0,0,q_\mathrm{fl},\,0,0,0)^\top$（列 3 と直交かつ単位ノルム）
* 列 5,6,7,8,9：anc=2 サブ空間に恒等的に配置すれば良い

具体的に、**$K_2 := 0$、anc=2 ブロックは恒等**とすると：

$$
\boxed{\;U_\mathrm{Stine}^\mathrm{(fl)} \;=\;
\left(\begin{array}{ccc|ccc|ccc}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & q_\mathrm{fl} & 0 & 0 & -\sqrt{p_\mathrm{fl}} & 0 & 0 & 0 \\\hline
0 & 0 & \sqrt{p_\mathrm{fl}} & 0 & 0 & q_\mathrm{fl} & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\\hline
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{array}\right)\;}
$$

ユニタリ性の確認：列 1, 2 は sys=$|0\rangle,|1\rangle$ で恒等。列 3 と列 6 は **2 元部分空間 $\{|0,2\rangle,|1,0\rangle\}$**（kron(anc,sys) で row index 2, 3 に対応）に閉じる $2\times 2$ Givens 回転 $\begin{pmatrix}q_\mathrm{fl} & -\sqrt{p_\mathrm{fl}} \\ \sqrt{p_\mathrm{fl}} & q_\mathrm{fl}\end{pmatrix}$。$q_\mathrm{fl}^2+p_\mathrm{fl}=1$ より直交かつ単位ノルム。残りの 4 列は単位ベクトル。$U^\dagger U = I_9$ 確認可。$\checkmark$

#### 5.2.1 ネイティブ qutrit ゲート分解

$U_\mathrm{Stine}^\mathrm{(fl)}$ は **2 qutrit (anc, sys)** に作用し、$\{|02\rangle_\mathrm{anc,sys},|10\rangle_\mathrm{anc,sys}\}$ という 2 元集合で 1 個の Givens 回転、残り 7 元では恒等。これは §3.4 で議論した 2-qutrit ペア演算子の構造と完全に同じ（添字が違うだけ）。よって **`CustomTwo`** で

$$
\texttt{CustomTwo(target=[anc, sys\_site], parameters=}U_\mathrm{Stine}^\mathrm{(fl)}\texttt{, dimensions=[3,3])}
$$

として渡すのが**実装可能で正直**な分解。`R` + `CSum` への分解は §3.4 と同様にコンパイラ依存。

### 5.3 TTA チャンネルの Stinespring ユニタリ（$d_\mathrm{loc}=9$、$U_\mathrm{Stine}\in\mathbb C^{27\times 27}$）

前文書 §6.1 の Kraus（チャンネル 1）：$K_0 = \mathrm{diag}_9(1,1,1,1,q_1,1,1,1,1)$、$K_1 = \sqrt{p_1}\,E^{(9)}_{6,4} = \sqrt{p_1}|20\rangle\langle 11|$、ここで $p_1 = 1-e^{-\gamma_\mathrm{TTA}\,dt/4}$、$q_1 = e^{-\gamma_\mathrm{TTA}\,dt/8}$。

最左 $27\times 9$ 列：

$$
\begin{pmatrix} K_0 \\ K_1 \\ 0 \end{pmatrix}_{27\times 9}
$$

非零成分は（kron(anc,sys) 規約、anc index $\in\{0,1,2\}$、sys index $\in\{0,\dots,8\}$）：

* anc=0 ブロック（rows 0..8）：$(K_0)_{rr}$、つまり 全対角に 1、ただし row 4 = $q_1$
* anc=1 ブロック（rows 9..17）：$(K_1)_{r,c} = \sqrt{p_1}\,\delta_{r,6}\,\delta_{c,4}$、つまり (row 9+6, col 4) = (row 15, col 4) のみ $\sqrt{p_1}$
* anc=2 ブロック（rows 18..26）：全 0

残り 18 列は等距離拡張。最も簡潔な選択（蛍光と並行）：

* 列 13（kron(anc,sys) で anc=1, sys=4）：row 4 に $-\sqrt{p_1}$、row 15 に $q_1$、他 0（$K_0$ の 4 列目「$q_1$」と $K_1$ の 4 列目「$\sqrt{p_1}$」を $2\times 2$ Givens でペアにした相方）
* 列 9..12, 14..17（anc=1 ブロック残り 8 列）：anc=1, sys=$c\neq 4$ で恒等（row 9+c, col 9+c に 1）
* 列 18..26（anc=2 ブロック）：恒等

**正直な評価**：$U_\mathrm{Stine}^\mathrm{(TTA1)}$ は $\{|1,4\rangle_\mathrm{anc,sys},|1,6\rangle\}\to\{|0,4\rangle,|1,6\rangle\}$ という 2 元 Givens 回転（前文書 §6.1 で導出した「sys index 4 = $|11\rangle_{ij}$」と「sys index 6 = $|20\rangle_{ij}$」を anc レベルで分岐）と、その他 25 元の恒等で構成される **$27\times 27$ ユニタリ**。これは **anc + 2 sys qutrit = 3-qutrit ゲート**であり、`CustomMulti(target=[anc, i, j], parameters=U_Stine^(TTA1), dimensions=[3,3,3])` として渡すのが正直。

TTA チャンネル 2（前文書 §6.2、$K_1 = \sqrt{p_1}|02\rangle\langle 11| = \sqrt{p_1}E^{(9)}_{2,4}$）も並行構造、$U_\mathrm{Stine}^\mathrm{(TTA2)}$ は $\{|1,4\rangle,|1,2\rangle\}$ Givens、他は恒等。

### 5.4 まとめ：10 個の KrausChannel ⇒ 10 個の Stinespring ユニタリ

| α | チャンネル | 補助 qutrit | Stinespring ユニタリのサイズ | qudit 命令 |
|---|---|---|---|---|
| 1 | TTA $(0,1)$ ch1 | anc1 | $27\times 27$ | $\mathrm{CustomMulti}(\mathrm{anc1},0,1; U_\mathrm{Stine}^\mathrm{(TTA1)})$ |
| 2 | TTA $(0,1)$ ch2 | anc2 | $27\times 27$ | $\mathrm{CustomMulti}(\mathrm{anc2},0,1; U_\mathrm{Stine}^\mathrm{(TTA2)})$ |
| 3–6 | TTA $(1,2),(2,3)$ | anc3..anc6 | $27\times 27$ | 同様 |
| 7 | 蛍光 site 0 | anc7 | $9\times 9$ | $\mathrm{CustomTwo}(\mathrm{anc7},0; U_\mathrm{Stine}^\mathrm{(fl)})$ |
| 8 | 蛍光 site 1 | anc8 | $9\times 9$ | $\mathrm{CustomTwo}(\mathrm{anc8},1; U_\mathrm{Stine}^\mathrm{(fl)})$ |
| 9 | 蛍光 site 2 | anc9 | $9\times 9$ | $\mathrm{CustomTwo}(\mathrm{anc9},2; U_\mathrm{Stine}^\mathrm{(fl)})$ |
| 10 | 蛍光 site 3 | anc10 | $9\times 9$ | $\mathrm{CustomTwo}(\mathrm{anc10},3; U_\mathrm{Stine}^\mathrm{(fl)})$ |

合計 **10 個の補助 qutrit**（anc1..anc10）と **10 個のユニタリ命令**で 10 個の KrausChannel を実現。各補助 qutrit は **$|0\rangle$ で初期化**し、ユニタリ後は **partial trace（破棄）** または **計算基底測定**（結果は使わない）。

### 5.5 パリンドロミック半ステップ（前文書 §5）の実現

前文書の 22 命令のうち forward 10 KrausChannel と reverse 10 KrausChannel は、**それぞれ独立な補助 qutrit セット**（forward 用 anc1..anc10 と reverse 用 anc11..anc20）で実装する必要がある（補助は使い捨てなので再利用するなら明示的にリセットゲートが必要）。よって 1 Trotter ステップで合計 **20 個の補助 qutrit + 4 個の system qutrit = 24 qutrit** が必要。

---

## 6. 1 Trotter ステップの完全 qutrit ゲート列

ケース II（H_0 と H_transfer も Strang 分割）を採用すると、1 ステップは：

| 区間 | 命令 | 個数 | 誤差 |
|---|---|---|---|
| (a) $U_H$ 前段 (§4.2) | VirtRz×8 + CustomTwo×5 + VirtRz×8 | **21** | $O(dt^3)$ |
| (b) forward 10 KrausChannel (§5) | CustomTwo×4 + CustomMulti×6 | **10** | (Stinespring leading order) |
| (c) reverse 10 KrausChannel (§5) | CustomTwo×4 + CustomMulti×6 | **10** | 同上 |
| (d) $U_H$ 後段 (§4.2) | VirtRz×8 + CustomTwo×5 + VirtRz×8 | **21** | $O(dt^3)$ |
| **合計** |  | **62 命令** |  |

ケース I（実コード踏襲）：

| 区間 | 命令 | 個数 |
|---|---|---|
| (a) $U_H$ 前段 | CustomMulti×1 ($U_H$ 黒箱) | 1 |
| (b),(c) 20 KrausChannel | KrausChannel×20（または CustomTwo/CustomMulti×20 via Stinespring） | 20 |
| (d) $U_H$ 後段 | CustomMulti×1 | 1 |
| **合計** |  | **22 命令** |

---

## 7. 実コードと本文書の差分（嘘禁止のための表）

| 項目 | 実コード `qudit_gksl_simulator.py` | 本文書 §2–§5 |
|---|---|---|
| Lindblad 数 | 26 | **10**（TTA 6 + 蛍光 4） |
| $U_H$ の表現 | `cu_multi(U_H)` 1 命令、$81\times 81$ 黒箱、Trotter 分割なし | §4.1（同じ）or §4.2（VirtRz 16 + CustomTwo 5、Strang 近似） |
| KrausChannel の表現 | `kraus_channel(K_α)` 命令、$d_\mathrm{loc}\times d_\mathrm{loc}$ 黒箱、ハードウェア非対応 | §5 で Stinespring ダイレーション → CustomTwo / CustomMulti、補助 qutrit 追加 |
| ネイティブ qutrit ゲート分解 | **無し**（cu_multi/kraus_channel は黒箱命令） | §2 (VirtRz、厳密)、§3 (CustomTwo、Strang 近似)、§5 (CustomTwo/CustomMulti、Stinespring) |
| 1 ステップの命令数 | 54（cu_multi×2 + kraus_channel×52、26 Lindblad 版）／**22**（10 Lindblad 版を実コード上で実行した場合） | **62**（10 Lindblad、ケース II、Stinespring 経由） |
| 補助 qutrit | 0（dmsim は density matrix を直接更新、補助不要） | 20（Stinespring forward + reverse 用） |

---

## 8. 本文書の数学的内容の要約

* **§2 で $\hat H_0$ 半ステップを 8 個の `VirtRz` で厳密分解**（Trotter 誤差ゼロ、可換）。各 `VirtRz` の $3\times 3$ 行列を §1.1 で完全に書き下した。
* **§3 で $\hat H_\mathrm{transfer}$ 半ステップを 5 個の `CustomTwo` (Strang) で分解**。各 `CustomTwo` の $9\times 9$ 行列 $E_\mathrm{loc}^\mathrm{tr}(\tau)$ を §3.2 で完全に（11 個の非零成分、70 個のゼロ）書き下した。誤差 $O((dt/2)^3)$（Strang）。
* **§4 で 2 つのケース** を比較：(I) 実コード踏襲（cu_multi 1 個、近似ゼロ、ハードウェア非対応）、(II) Strang 分離（21 ゲート、誤差 $O(dt^3)$、ハードウェア対応）。
* **§5 で 10 個の KrausChannel を Stinespring ダイレーションで qudit ユニタリに変換**。各 $U_\mathrm{Stine}$ を行列形（蛍光は $9\times 9$、TTA は $27\times 27$）で完全に書き下し、ユニタリ性を成分計算で検証。`CustomTwo` / `CustomMulti` 命令にマップ。補助 qutrit 10 個を要する（`stinespring_utils` の "$d_\mathrm{anc}=3$" モード相当、ただしダイレーションを 1 個ずつ独立補助で行う場合）。
* **§6 で 1 Trotter ステップ全体を 62 命令で完全に列挙**（ケース II、forward + reverse 含む）。
* **§7 で実コードとの差分を表で正直に明示**：実コードは黒箱 `cu_multi` / `kraus_channel` のみで、ネイティブ qutrit ゲート分解は持たない。本文書 §2–§5 はそれを補完する**理論的分解**である。

本文書の記述：「§6 表で示した 62 命令」「§4.2 で $H_0$ と $H_\mathrm{transfer}$ の Strang 分割」については、**補遺 `mydoc/シナリオ5d_TTAと蛍光のみ_Qudit量子ゲート展開_補遺.md` 補遺 H・補遺 L で厳密化・修正した**。具体的には：

* **補遺 H**：$[\hat H_0,\hat H_\mathrm{transfer}]=0$ を演算子代数で証明 → ケース II は厳密に分離可能、Strang 不要、誤差なし。
* **補遺 L**：1 Trotter ステップは正しくは **46 命令**（VirtRz 16 + CustomTwo 10 + Stinespring KrausChannel 20）。本文書 §6 の「62 命令」は H_0/H_transfer Strang を仮定した過剰集計で、補遺 H で不要と判明。
* **補遺 A〜J**：本文書で省略していた行列（GellMann、$X_3^k$、$\mathrm{CSum}^\dagger$、$h_\mathrm{loc}^\mathrm{tr}$、$E_\mathrm{loc}^\mathrm{tr}$、$U_\mathrm{Stine}^\mathrm{(TTA1,TTA2)}$）を全成分・全 27/9 列で完全展開。
