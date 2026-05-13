# シナリオ5d: Qudit GKSL（ボソン無し, N=4）を MQT-Qudits density-matrix backend (DMSim) で実行 — 省略無しの完全理論

## 0. 本文書の位置付けと正直な前提

本文書は `tutorials/quantum_dynamics_gksl_comparison.ipynb` の Cell 8/9（「3b. シナリオ5d」）で実行される計算の**全てのステップを、ヒューリスティック処理を一切持ち込まずに**数式で展開したものである。記述はソースコード（`tutorials/qudit_gksl_simulator.py`、`tutorials/exact_local_channels.py`、`tutorials/dmsim_kraus_helpers.py`、`tutorials/gksl_math_utils.py`、`tutorials/gksl_physical_parameters.py`、`src/mqt/qudits/simulation/backends/dmsim.py`、`src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`）と完全に一致している。コードに無いものは書かない。コードがしている近似は近似と書く。コードが厳密にやっていることは厳密と書く。

シナリオ5dの構成は次の3層からなる：

1. **物理層**: $N=4$分子・各分子 qutrit ($d=3$) の電子 GKSL 方程式（boson 無し）。これはシナリオ5と同一。
2. **数値積分層**: Strang 分割 (Hamiltonian–Dissipator) + パリンドロミック Trotter（散逸子間） + 各局所 Lindblad チャネルを **`exact_local_channels` アルゴリズム**（厳密局所超演算子の指数化）で時間発展させる。これはシナリオ5の `algorithm="exact_local_channels"` と同一。**Stinespring 拡張は使わない**。
3. **実行層**: 各 Trotter ステップで MQT-Qudits の `QuantumCircuit` を構築し、各局所散逸チャネルを Choi–Jamiolkowski 同型から厳密に抽出した Kraus 演算子の `KrausChannel` 命令として並べ、密度行列バックエンド `DMSim` で実行する。これがシナリオ5dに固有の部分である。

「DMSim は古典計算機上で密度行列を直接演算する simulator backend」であり、**実量子ハードウェアの上で動いているわけではない**（state-vector simulator も同様）。これは `tutorials/quantum_dynamics_gksl_comparison.ipynb` の Cell 8 マークダウンに明記されている honest disclosure と一致する。

---

## 1. 物理パラメータ（コードで使われる正確な値）

`tutorials/gksl_physical_parameters.py` の `GKSLPhysicalParameters` で生成される、Cell 2 のデフォルトパラメータ：

| 記号 | 名前 | 値 | 単位 |
|---|---|---|---|
| $E_T$ | 三重項エネルギー | $1.5$ | eV |
| $E_S$ | 一重項エネルギー | $3.0$ | eV |
| $V$ | 移動結合定数 | $0.1$ | eV |
| $\gamma_{\mathrm{TTA}}$ | TTA 速度 | $0.05$ | eV/$\hbar$ |
| $\Gamma_{\mathrm{fl}}$ | 蛍光速度 | $0.01$ | eV/$\hbar$ |
| $\Gamma_{\mathrm{ph}}$ | 燐光速度 | $10^{-6}$ | eV/$\hbar$ |
| $k_{\mathrm{IC}}$ | 内部転換速度 | $0.005$ | eV/$\hbar$ |
| $k_{\mathrm{ISC},S\to T}$ | 項間交差 S→T | $0.003$ | eV/$\hbar$ |
| $k_{\mathrm{ISC},T\to S}$ | 項間交差 T→S | $10^{-5}$ | eV/$\hbar$ |
| $N$ | 分子数 | $4$ | — |
| $d$ | 局所次元 | $3$ | — |

時間スケールは内部単位 $\hbar=1$。系全体のヒルベルト空間次元は $D = d^N = 3^4 = 81$、密度行列は $81\times 81$ の複素行列（メモリ約 $16 \cdot 81^2 \approx 100$ KB）。

近接ペアは 1次元鎖 $(0,1), (1,2), (2,3)$ の3組（`params.neighbors`）。

シミュレーション設定（Cell 2）：$t_{\max}=100.0$、`n_steps_quantum = 100` ⇒ 時間刻み

$$
dt \;=\; \frac{t_{\max}}{n_\text{steps}} \;=\; \frac{100.0}{100} \;=\; 1.0
$$

---

## 2. 単一分子の状態空間と多体ヒルベルト空間

各分子は3準位系で

$$
|0\rangle = \mathrm{S}_0,\qquad |1\rangle = \mathrm{T}_1,\qquad |2\rangle = \mathrm{S}_1.
$$

多体系のヒルベルト空間は

$$
\mathcal{H} \;=\; \bigotimes_{i=0}^{N-1} \mathcal{H}_i,\qquad \dim\mathcal{H} = d^N = 81.
$$

計算基底は $|s_0 s_1 s_2 s_3\rangle$ ($s_i\in\{0,1,2\}$) で、Kronecker 積の慣例（コード中の `reduce(np.kron, [...])`）により基底ラベルは

$$
\mathrm{idx}(s_0,s_1,s_2,s_3) \;=\; s_0 \cdot d^3 + s_1 \cdot d^2 + s_2 \cdot d + s_3.
$$

---

## 3. ハミルトニアン

### 3.1 オンサイト項（`build_onsite_hamiltonian`）

各分子の局所オンサイトは

$$
h_\mathrm{loc} = \mathrm{diag}(0,\,E_T,\,E_S) = \begin{pmatrix}0 & 0 & 0\\ 0 & E_T & 0\\ 0 & 0 & E_S\end{pmatrix},
$$

そして

$$
\hat{H}_0 \;=\; \sum_{i=0}^{N-1}\; \underbrace{I^{\otimes i}\otimes h_\mathrm{loc}\otimes I^{\otimes (N-1-i)}}_{=:\,h_\mathrm{loc}^{(i)}}.
$$

### 3.2 移動項（`build_transfer_hamiltonian`）

近接対 $(i,j)\in\{(0,1),(1,2),(2,3)\}$ に対し、$T_1$ 励起の双方向ホッピング

$$
\hat{H}_{\mathrm{transfer}} \;=\; V\sum_{\langle i,j\rangle}\Big(|0\rangle_i\langle 1| \otimes |1\rangle_j\langle 0| \;+\; |1\rangle_i\langle 0| \otimes |0\rangle_j\langle 1|\Big),
$$

他のサイトには恒等演算子を Kronecker 積する。

### 3.3 全ハミルトニアン

$$
\boxed{\;\hat{H}_\mathrm{total} \;=\; \hat{H}_0 \;+\; \hat{H}_{\mathrm{transfer}}\;}\quad\text{（$81\times 81$ Hermitian）}
$$

`build_onsite_hamiltonian` と `build_transfer_hamiltonian` はそれぞれ Hermitian 性を `np.allclose(H, H.conj().T)` で実測検証している。

---

## 4. Lindblad 演算子（GKSL 散逸子）

GKSL（Gorini–Kossakowski–Sudarshan–Lindblad）方程式

$$
\dot{\hat\rho} \;=\; -i\,[\hat H_\mathrm{total},\,\hat\rho] \;+\; \sum_\alpha \mathcal{D}[\hat L_\alpha](\hat\rho),
$$

$$
\mathcal{D}[\hat L](\hat\rho) \;:=\; \hat L\,\hat\rho\,\hat L^\dagger \;-\; \tfrac12\,\{\hat L^\dagger \hat L,\,\hat\rho\}.
$$

Lindblad 演算子は `build_lindblad_operators` でちょうど

$$
n_\mathrm{Lin} \;=\; 2\,|\text{neighbors}| \;+\; 5N \;=\; 2\cdot 3 + 5\cdot 4 \;=\; 26
$$

個生成される。各演算子は **$\sqrt{\gamma}$ プレファクタを既に含んでいる**。

### 4.1 TTA（三重項–三重項消滅、ペア演算子, 6個）

近接対 $(i,j)$ ごとに2チャネル：

$$
\hat L_{\mathrm{TTA},1}^{(i,j)} \;=\; \sqrt{\gamma_\mathrm{TTA}/2}\;\Big(|2\rangle_i\langle 1|\,\otimes\,|0\rangle_j\langle 1|\Big),
$$
$$
\hat L_{\mathrm{TTA},2}^{(i,j)} \;=\; \sqrt{\gamma_\mathrm{TTA}/2}\;\Big(|0\rangle_i\langle 1|\,\otimes\,|2\rangle_j\langle 1|\Big).
$$

ここで $\gamma$ は $\gamma_\mathrm{TTA}/2$（コードの `params.gamma_TTA / 2.0`）であり、合計 TTA 速度 $\gamma_\mathrm{TTA}$ が2チャネルに均等に分配されている。

### 4.2 単一サイト演算子（5チャネル × $N=4$ サイト = 20個）

各分子 $i$ について、

| チャネル | $\hat L^{(i)}$（局所部分） | 速度 $\gamma$ |
|---|---|---|
| 蛍光 | $\sqrt{\Gamma_\mathrm{fl}}\,\|0\rangle\langle 2\|$ | $\Gamma_\mathrm{fl}$ |
| 燐光 | $\sqrt{\Gamma_\mathrm{ph}}\,\|0\rangle\langle 1\|$ | $\Gamma_\mathrm{ph}$ |
| 内部転換 (IC) | $\sqrt{k_\mathrm{IC}}\,\|0\rangle\langle 2\|$ | $k_\mathrm{IC}$ |
| 項間交差 S→T | $\sqrt{k_{\mathrm{ISC},S\to T}}\,\|1\rangle\langle 2\|$ | $k_{\mathrm{ISC},S\to T}$ |
| 項間交差 T→S | $\sqrt{k_{\mathrm{ISC},T\to S}}\,\|0\rangle\langle 1\|$ | $k_{\mathrm{ISC},T\to S}$ |

これらを $i=0,\dots,N-1$ に対し `build_single_site_operator` で全体空間に持ち上げる：

$$
\hat L^{(i)}_\bullet \;=\; I\otimes\cdots\otimes\underbrace{\hat L_\bullet^\mathrm{loc}}_{\text{$i$ 番目}}\otimes\cdots\otimes I.
$$

順序は `build_lindblad_operators` の通り（外側ループ：チャネル種別、内側ループ：サイト index）であり、後段の `get_local_lindblad_ops` でも同じ順序で再構築される。**順序は恣意的ではない**：パリンドロミック構造は対称な順序であれば交換子の主要項を消すので、`exact_local_channels` の forward/reverse 巡回で順序が一致していれば良い（後述）。

---

## 5. 数値積分の戦略：Strang 分割 + パリンドロミック散逸子 + 厳密局所超演算子

時間発展は Liouvillian

$$
\mathcal{L} \;=\; \mathcal{L}_H + \mathcal{L}_D,\qquad \mathcal{L}_H(\rho) = -i[\hat H_\mathrm{total},\rho],\qquad \mathcal{L}_D = \sum_\alpha \mathcal{D}[\hat L_\alpha]
$$

の指数 $\rho(t+dt) = e^{\mathcal{L}\,dt}\rho(t)$ で与えられる。コードはこれを次の Strang 2次分割で近似する：

$$
\boxed{\;e^{\mathcal{L}\,dt} \;\approx\; e^{\mathcal{L}_H\,dt/2}\;\cdot\;\Phi_\mathrm{pal}(dt)\;\cdot\;e^{\mathcal{L}_H\,dt/2}\;}
$$

ここで $\Phi_\mathrm{pal}(dt)$ は散逸子に対するパリンドロミック（対称）チャネル積：

$$
\Phi_\mathrm{pal}(dt) \;=\; \Big(\prod_{\alpha=1}^{n_\mathrm{Lin}} e^{\mathcal{D}[\hat L_\alpha]\,dt/2}\Big)\;\Big(\prod_{\alpha=n_\mathrm{Lin}}^{1} e^{\mathcal{D}[\hat L_\alpha]\,dt/2}\Big).
$$

forward 列と reverse 列の組み合わせにより、$\mathcal{D}[L_\alpha]$ 同士の交換子に由来する Lie–Trotter 誤差の 1次項が打ち消され、ステップ毎の誤差は $O(dt^3)$（積算で大域 $O(dt^2)$）になる。

### 5.1 ハミルトニアン半ステップ（純ユニタリ）

$$
e^{\mathcal{L}_H\,dt/2}(\rho) \;=\; U_H\,\rho\,U_H^\dagger,\qquad U_H \;=\; \exp\!\big(-i\,\hat H_\mathrm{total}\,dt/2\big).
$$

コード上、$U_H$ は `scipy.linalg.expm(-1j * H_total * dt / 2)` で**シミュレーション開始時に1度だけ**計算され `self._U_H_half` にキャッシュされる（ループ内で `expm` は呼ばれない）。

### 5.2 各局所 Lindblad チャネルの「厳密な」局所超演算子化

各 $\hat L_\alpha$ は最大 2 サイトでしか非自明に作用しない。よって散逸子も対応する局所部分空間（次元 $d=3$ または $d^2=9$）でのみ非自明である。

#### 5.2.1 局所散逸超演算子の構成（`build_local_dissipator_super`）

局所 Lindblad 演算子 $L \in \mathbb{C}^{d_\mathrm{loc}\times d_\mathrm{loc}}$ に対し、**列優先 (column-major) ベクトル化** $\mathrm{vec}(X)_{i + d_\mathrm{loc} \cdot j} = X_{ij}$ を用いて、局所散逸子の超演算子表現は

$$
\boxed{\;\mathcal{L}_D^\mathrm{loc} \;=\; (L^* \otimes L) \;-\; \tfrac12\,(I \otimes L^\dagger L) \;-\; \tfrac12\,((L^\dagger L)^\top \otimes I)\;}
$$

形状は単一サイト ($d_\mathrm{loc}=3$) で $9\times 9$、TTA ペア ($d_\mathrm{loc}=9$) で $81\times 81$。

導出：$\rho \mapsto L\rho L^\dagger$ は $\mathrm{vec}(L\rho L^\dagger) = (L^*\otimes L)\,\mathrm{vec}(\rho)$（列優先 vec の標準恒等式）、$\rho \mapsto L^\dagger L\,\rho$ は $\mathrm{vec}(L^\dagger L\,\rho)=(I\otimes L^\dagger L)\,\mathrm{vec}(\rho)$、$\rho \mapsto \rho\,L^\dagger L$ は $\mathrm{vec}(\rho L^\dagger L)=((L^\dagger L)^\top\otimes I)\,\mathrm{vec}(\rho)$ から従う。

#### 5.2.2 半ステップ局所チャネルの厳密指数化

各 $\hat L_\alpha$ の局所超演算子 $\mathcal{L}_{D,\alpha}^\mathrm{loc}$ から

$$
M_\alpha^{(\mathrm{half})} \;:=\; \exp\!\big(\mathcal{L}_{D,\alpha}^\mathrm{loc}\cdot dt/2\big)
$$

を `scipy.linalg.expm` で**シミュレーション開始時に1度だけ**計算する（`precompute_exact_channels_half`、ループ外）。$M_\alpha^{(\mathrm{half})}$ は $\rho_\mathrm{loc}$ への厳密な CPTP 写像（$dt/2$ 時間の単一散逸子チャネル）の超演算子表現である。**Stinespring 拡張による1次近似は使わない**。

これがシナリオ5dの核心：散逸子の指数 $e^{\mathcal{D}[L_\alpha]\,dt/2}$ を局所空間で**完全に厳密に**取り、近似は **Strang分割の H–D 交換子（$O(dt^2)$ ステップ毎）** と **散逸子間のパリンドロミック交換子（$O(dt^3)$ ステップ毎）** だけになる。

---

## 6. シナリオ5d 固有の処理：MQT-Qudits 回路と DMSim 上での実行

ここからがシナリオ5d 本体である。$M_\alpha^{(\mathrm{half})}$ を MQT-Qudits の `QuantumCircuit` に **`KrausChannel` 命令として** 渡すために、Choi–Jamiolkowski 同型を経由して **Kraus 演算子の集合** に変換する。

### 6.1 Choi–Jamiolkowski 同型による厳密 Kraus 抽出（`kraus_from_local_superoperator`）

任意の CPTP 写像 $\mathcal{E}: \mathcal{B}(\mathbb{C}^d)\to\mathcal{B}(\mathbb{C}^d)$ に対し、Choi 行列

$$
C(\mathcal{E}) \;=\; \sum_{i,j=0}^{d-1} \mathcal{E}\big(|i\rangle\langle j|\big) \,\otimes\, |i\rangle\langle j|
$$

は半正定値である。$\mathcal{E}$ の超演算子表現 $M$（列優先 vec、$M_{i+d\cdot j,\,k+d\cdot l} = T_{ijkl}$ で $\mathcal{E}(\rho)_{ij} = \sum_{kl} T_{ijkl}\,\rho_{kl}$）から Choi 行列は

$$
C_{(i\,k),(j\,l)} \;=\; T_{ijkl}
$$

で得られる（コードでは `t_tensor.transpose(0,2,1,3).reshape(d^2, d^2)` がこれに対応；行は $(i,k)$、列は $(j,l)$、両者とも row-major フラット化）。

#### 数値的厳密性の保証

1. **Hermitization**: 数値的アンチエルミート部（典型 $10^{-15}$）を抑えるため `choi = 0.5 * (choi + choi.conj().T)`。
2. **固有分解**: `np.linalg.eigh(choi)` で固有値・固有ベクトル $(\lambda_\alpha, v_\alpha)$ を得る（$C$ は Hermitian なので実固有値）。
3. **物理性チェック（嘘禁止）**: もし $\lambda_\alpha < -\texttt{CHOI\_EIG\_TOL}$（$\texttt{CHOI\_EIG\_TOL}=10^{-12}$）なる固有値が存在すれば、入力超演算子が CPTP ではないことを意味し `ValueError` を **raise する**（"No heuristic correction applied." とコードコメント明記）。
4. **数値ノイズ閾値**: $|\lambda_\alpha| < 10^{-12}$ の固有値は **0 として捨てる**（測定不可能な数値ノイズ）。

#### Kraus 演算子の構成

各 $\lambda_\alpha \geq 10^{-12}$ について

$$
\boxed{\;K_\alpha \;=\; \sqrt{\lambda_\alpha}\;\mathrm{reshape}(v_\alpha,\,(d,d))\quad\text{（row-major reshape: $v_\alpha[i\cdot d + k] \mapsto K_\alpha[i,k]$）}\;}
$$

これにより

$$
\mathcal{E}(\rho) \;=\; \sum_\alpha K_\alpha\,\rho\,K_\alpha^\dagger
$$

が**任意の $\rho$ に対して厳密に**成立する（教科書通りの Choi–Jamiolkowski 構成）。Kraus 数 $|\{\alpha\}|$ は $C$ の rank に等しく $\leq d^2$。

#### 完全性条件 $\sum_\alpha K_\alpha^\dagger K_\alpha = I$

これは元の局所超演算子 $M_\alpha^{(\mathrm{half})}$ がトレース保存写像（$\hat L_\alpha$ が GKSL 散逸子由来なので自動的に CPTP）であることに同値であり、Kraus を構築した後 `KrausChannel.__init__` で

$$
\Big\|\sum_\alpha K_\alpha^\dagger K_\alpha - I\Big\|_F \;\leq\; \texttt{CPTP\_TOLERANCE} = 10^{-9}
$$

として実測検証される。違反時は `ValueError`（"No heuristic correction is applied." と明記）。

### 6.2 1 Trotter ステップの MQT-Qudits 回路（`_trotter_step_dmsim`）

各 Trotter ステップで、$N=4$ qutrit の `QuantumCircuit`（量子レジスタ `QuantumRegister("sys", 4, [3,3,3,3])`）を**新規構築**して `DMSim` で実行する。命令は MQT-Qudits の `mqt.qudits.quantum_circuit.QuantumCircuit` API を介して `circuit.cu_multi(...)` および `circuit.kraus_channel(...)` で append される。**回路全体としての命令列は次の通り**（コード `_trotter_step_dmsim` がまさにこの順で append している）：

```
[QuantumRegister "sys" : 4 qutrit, dimensions = [3, 3, 3, 3]]

(a)  cu_multi(target=[0,1,2,3], parameters=U_H_half)                      ← ハミルトニアン半ステップ前段
(b1) kraus_channel(target=[0,1], kraus_operators={K_{1,β}})               ← TTA pair (0,1) チャネル1
(b2) kraus_channel(target=[0,1], kraus_operators={K_{2,β}})               ← TTA pair (0,1) チャネル2
(b3) kraus_channel(target=[1,2], kraus_operators={K_{3,β}})               ← TTA pair (1,2) チャネル1
(b4) kraus_channel(target=[1,2], kraus_operators={K_{4,β}})               ← TTA pair (1,2) チャネル2
(b5) kraus_channel(target=[2,3], kraus_operators={K_{5,β}})               ← TTA pair (2,3) チャネル1
(b6) kraus_channel(target=[2,3], kraus_operators={K_{6,β}})               ← TTA pair (2,3) チャネル2
(b7..b10)  kraus_channel(target=[i], kraus_operators={K_{6+1+i,β}})  i=0..3 ← 蛍光 (4個)
(b11..b14) kraus_channel(target=[i], kraus_operators={K_{10+1+i,β}}) i=0..3 ← 燐光 (4個)
(b15..b18) kraus_channel(target=[i], kraus_operators={K_{14+1+i,β}}) i=0..3 ← IC (4個)
(b19..b22) kraus_channel(target=[i], kraus_operators={K_{18+1+i,β}}) i=0..3 ← ISC S→T (4個)
(b23..b26) kraus_channel(target=[i], kraus_operators={K_{22+1+i,β}}) i=0..3 ← ISC T→S (4個)
                                                                            ↑ ここまで forward 列（26 命令）

(c1..c26) (b1..b26) を完全に逆順で再 append                                  ← reverse 列（26 命令）

(d)  cu_multi(target=[0,1,2,3], parameters=U_H_half)                      ← ハミルトニアン半ステップ後段
```

**合計命令数 = $1 + 26 + 26 + 1 = 54$**（前回の文書で「56」「`cu_multi` × 2 + `KrausChannel` × 52」と書いた箇所は、より正確には「`cu_multi` × 2 + `KrausChannel` × 52、計 54 命令」である。`cu_multi` は 1 命令で全 4 qutrit に作用する 1 つの $81\times 81$ 行列である）。

実行：

```
job = self._dmsim_backend.run(circuit, initial_density_matrix=ρ_prev)
ρ_next = job.result().get_density_matrix()
```

`initial_density_matrix` は前ステップの $\rho$（最初のステップは §7 の初期状態）。

### 6.2.1 ゲート命令の正確な定義（型・行列・作用）

#### (i) `cu_multi`（`CustomMulti` ゲート）

`mqt.qudits.quantum_circuit.gates.CustomMulti` は `GateTypes.MULTI` のユーザ定義ユニタリゲートで、コンストラクタ引数 `parameters` に渡された複素行列をそのまま自分の行列として保持する（`__array__` メソッドが `self.__array_storage` を返す）。シナリオ5dでは

$$
\boxed{\;\texttt{cu\_multi}(\,[0,1,2,3],\;U_H^{(\mathrm{half})}\,)\;,\qquad U_H^{(\mathrm{half})} \;=\; \exp\!\big(-i\,\hat H_\mathrm{total}\,dt/2\big)\;\in\;\mathbb{C}^{81\times 81}\;}
$$

を 1 ステップにつき 2 回（前段・後段）使う。$U_H^{(\mathrm{half})}$ は `_precompute_unitaries(dt)` で `scipy.linalg.expm(-1j * H_total * dt / 2)` により**シミュレーション開始時に1度だけ**計算される（`self._U_H_half`）。

**重要な真実**：この `cu_multi` は **MQT-Qudits の compiler を通していない**。すなわち native gate 集合 (`VirtRz`, `R`, `Rh`, `Rz`, `CEx`) への分解は**行わない**。DMSim は §6.3.1 の通り `instruction.to_matrix(identities=0)` で $81\times 81$ 行列を取り出してそのまま $\rho$ の両側から掛ける。これは「量子回路の **論理レベル** での実行」であり、**物理2量子ビットゲート列としての実行ではない**。シナリオ5d は「density-matrix backend に正しい数学を委ねた」状態であり、native gate decomposition の誤差や CEx の物理ノイズは**この経路には入っていない**（compiler 経由の native 分解は別経路 `compute_compiler_measured_gate_counts`／`QuditGKSLKrausSimulator` でのみ行われる）。

`cu_multi` 命令の DMSim における代数的作用は

$$
\rho \;\longmapsto\; U_H^{(\mathrm{half})}\,\rho\,(U_H^{(\mathrm{half})})^\dagger,
$$

これは Liouvillian の Hamiltonian 部分に対する**厳密な**半ステップ伝播 $e^{\mathcal{L}_H\,dt/2}$ である。

#### (ii) `kraus_channel`（`KrausChannel` 命令）

`mqt.qudits.quantum_circuit.gates.KrausChannel` は **非ユニタリ** 命令（`GateTypes.SINGLE` または `GateTypes.TWO`、ターゲット qudit 数で決まる）で、コンストラクタ時に

$$
\Big\|\sum_\beta K_{\alpha,\beta}^\dagger K_{\alpha,\beta} \;-\; I_{d_\mathrm{loc}}\Big\|_F \;\leq\; 10^{-9}
$$

を Frobenius ノルムで実測検証する（違反時 `ValueError`、ヒューリスティック補正なし）。`to_matrix()` および `__array__` は

> "KrausChannel has no unitary matrix representation. Use a density-matrix backend (e.g. DMSim) instead."

として `NotImplementedError` を raise する（=「ユニタリ行列としては存在しない」ことが型レベルで保証されており、嘘がつけない設計）。

シナリオ5d の 1 Trotter ステップで使う $\alpha = 1, \dots, 26$ の各チャネルは、§6.1 の Choi–Jamiolkowski 同型から厳密に抽出された半ステップ Kraus 集合 $\{K_{\alpha,\beta}\}_{\beta=1}^{r_\alpha}$（$r_\alpha \leq d_{\mathrm{loc},\alpha}^2$、Choi 行列の数値的 rank）であり、その**作用は**

$$
\boxed{\;\rho \;\longmapsto\; \sum_{\beta=1}^{r_\alpha} \big(K_{\alpha,\beta}\big)_{S_\alpha}\;\rho\;\big(K_{\alpha,\beta}^\dagger\big)_{S_\alpha}\;}
$$

ここで $S_\alpha$ はターゲット qudit 集合（pair チャネルなら 2 qudit、single チャネルなら 1 qudit）。これは厳密に局所散逸子チャネル $e^{\mathcal{D}[\hat L_\alpha]\,dt/2}$ そのものである（§5.2.2 と §6.1 から $\sum_\beta K_{\alpha,\beta} \rho K_{\alpha,\beta}^\dagger = M_\alpha^{(\mathrm{half})}(\rho)$ が解析的恒等式として成立）。

#### (iii) Kraus 演算子の数（rank）の上限

* TTA ペアチャネル ($d_\mathrm{loc}=9$): $r_\alpha \leq 81$
* 単一サイトチャネル ($d_\mathrm{loc}=3$): $r_\alpha \leq 9$

実際の rank はチャネルごとに Choi 行列の固有値分布で決まる（数値的に $|\lambda| < 10^{-12}$ は捨てられる、§6.1）。チャネルの「Kraus rank」は GKSL チャネルの構造で決まる物理量であり、隠さずそのまま使う。

### 6.2.2 1 Trotter ステップに対応する数学的命題

上の命令列を順に DMSim で実行することは、$\rho$ に対し

$$
\rho \;\longmapsto\; \mathcal{U}_H^{(\mathrm{half})} \;\circ\; \overbrace{\mathcal{E}_{26}^{(\mathrm{half})}\circ\cdots\circ\mathcal{E}_1^{(\mathrm{half})}}^{\text{forward 列}} \;\circ\; \overbrace{\mathcal{E}_1^{(\mathrm{half})}\circ\cdots\circ\mathcal{E}_{26}^{(\mathrm{half})}}^{\text{reverse 列}} \;\circ\; \mathcal{U}_H^{(\mathrm{half})} \;(\rho)
$$

を作用させることと**完全に同値**である。ここで

$$
\mathcal{U}_H^{(\mathrm{half})}(\rho) := U_H^{(\mathrm{half})}\rho (U_H^{(\mathrm{half})})^\dagger,\qquad \mathcal{E}_\alpha^{(\mathrm{half})}(\rho) := \sum_\beta K_{\alpha,\beta}\rho K_{\alpha,\beta}^\dagger \;=\; e^{\mathcal{D}[\hat L_\alpha]\,dt/2}(\rho).
$$

（命令列が左から右に append され、`DMSim.execute` が `for instruction in circuit.instructions: rho = ...` で**順次** $\rho$ を更新する、という DMSim の実装の単純な帰結である。回路慣習で「左から書いた命令が**最初に**作用する」のと同じ。）

forward 列と reverse 列を併せた散逸子作用は、§5 のパリンドロミック構造により $e^{\mathcal{L}_D\,dt} + O(dt^3)$ と一致し、Hamiltonian 半ステップ 2 個と合わせると Strang 分割により全体として $e^{\mathcal{L}\,dt} + O(dt^3)$（ステップ毎）になる。これがシナリオ5d の 1 ステップで「量子ゲート列として」実装している時間発展そのものである。

### 6.2.3 全シミュレーションの**等価な量子回路図**

$n_\mathrm{steps}=100$ ステップ全体は、上記 1 ステップ回路 $\mathcal{C}_{1\text{step}}$ を 100 回直列に並べたもの（毎ステップで $\rho$ を引き継ぎながら実行）と等価で、概念的には次のような巨大回路になる：

```
  q[0] ─┬───[U_H/2]───┬─Kraus(TTA01,1)─Kraus(TTA01,2)─Kraus(Fl0)─Kraus(Ph0)─Kraus(IC0)─Kraus(IST0)─Kraus(ITS0)─...─[reverse]─...─[U_H/2]─┬─...
  q[1] ─┤             ├─Kraus(TTA01,1)─Kraus(TTA01,2)─Kraus(TTA12,1)─Kraus(TTA12,2)─Kraus(Fl1)...                                       │
  q[2] ─┤  cu_multi   ├─                Kraus(TTA12,1)─Kraus(TTA12,2)─Kraus(TTA23,1)─Kraus(TTA23,2)─...                                  │ cu_multi
  q[3] ─┴─────────────┴─                                              Kraus(TTA23,1)─Kraus(TTA23,2)─Kraus(Fl3)...                       │
        │←─ ステップ 1 ─────────────────────────────────────────────────────────────────────────────────────────────────────────────→│ ステップ 2 →…
```

ここで `Kraus(TTA01,1)` は target=[0,1] の `KrausChannel` 命令（`K_{1,β}` 集合）、`Kraus(Fl0)` は target=[0] の `KrausChannel` 命令（蛍光チャネルの Kraus 集合）等である。`cu_multi` は 4 qutrit すべてを横断する単一の論理ゲート命令で、Hamiltonian の半ステップ伝播 $U_H^{(\mathrm{half})}$（$81\times 81$）そのものを表す。

**注（嘘禁止のための明示）**：上図の `cu_multi` は MQT-Qudits の論理ゲート命令そのものであり、ハードウェアレベルの「2 qudit native gate（CEx 等）の分解結果」ではない。シナリオ5d の DMSim 経路では、この `cu_multi` を **そのまま $81\times 81$ ユニタリとして** $\rho$ に作用させているのであって、CEx + 単一 qudit 回転の連鎖に展開してから実行しているのではない。`KrausChannel` は同様に「非ユニタリ命令」そのものとして DMSim が実行しており、これも native gate 列としては存在しない（`to_matrix` が `NotImplementedError` を raise することで型レベルで保証されている）。

### 6.2.4 1 ステップ回路の**意味論的時間発展**

1 ステップ回路 $\mathcal{C}_{1\text{step}}$ を $\rho$ に作用させるとは、上記命令列の合成超演算子 $\Phi(dt)$ を作用させることと同じである：

$$
\rho(t+dt) \;=\; \Phi(dt)\big(\rho(t)\big),\qquad \Phi(dt) \;=\; \mathcal{U}_H^{(\mathrm{half})}\circ\Phi_\mathrm{pal}(dt)\circ\mathcal{U}_H^{(\mathrm{half})},
$$

$$
\Phi_\mathrm{pal}(dt) \;=\; \Big(\bigcirc_{\alpha=1}^{26} \mathcal{E}_\alpha^{(\mathrm{half})}\Big) \;\circ\; \Big(\bigcirc_{\alpha=26}^{1} \mathcal{E}_\alpha^{(\mathrm{half})}\Big).
$$

これと §10 の Strang + パリンドロミック誤差解析の合成より、

$$
\Phi(dt) \;=\; e^{(\mathcal{L}_H+\mathcal{L}_D)\,dt} \;+\; O(dt^3)\quad(\text{ステップ毎}),
$$

すなわち**この回路命令列を 100 回繰り返すことが、4 分子 qutrit GKSL 方程式 $\dot\rho = \mathcal{L}\rho$ の Strang+パリンドロミック離散化の数値解そのものである**。これがシナリオ5d で「量子回路を作って計算した」ことの内容である。

### 6.2.5 1 ステップに含まれる「量子ゲート」と「量子チャネル」の正確な内訳

| 種別 | MQT-Qudits 命令 | ターゲット | 個数（1 ステップあたり） | 行列 / Kraus 集合 |
|---|---|---|---|---|
| ユニタリ | `cu_multi`（`CustomMulti`） | $[0,1,2,3]$ | 2 | $U_H^{(\mathrm{half})} = e^{-i\hat H_\mathrm{total}\,dt/2}\in\mathbb{C}^{81\times 81}$ |
| 非ユニタリ（ペア） | `kraus_channel`（`KrausChannel`） | $[0,1]$, $[1,2]$, $[2,3]$ | $2\times 3\times 2=12$（forward 6 + reverse 6） | TTA $\{K_{\alpha,\beta}\}$、各 $K\in\mathbb{C}^{9\times 9}$ |
| 非ユニタリ（単一） | `kraus_channel`（`KrausChannel`） | $[0],[1],[2],[3]$ | $2\times 5\times 4=40$（forward 20 + reverse 20） | 蛍光・燐光・IC・ISC$_{S\to T}$・ISC$_{T\to S}$ の $\{K_{\alpha,\beta}\}$、各 $K\in\mathbb{C}^{3\times 3}$ |
| **合計** | | | **54 命令** | |

100 ステップ全体では $54 \times 100 = 5400$ 命令の回路を、`DMSim.execute` がステップ毎に構築・実行する（計 100 個の独立な `QuantumCircuit` オブジェクトが生成される；`initial_density_matrix` で $\rho$ を引き継ぐ）。

### 6.3 DMSim バックエンドが各命令に対して実行する厳密な代数

DMSim (`src/mqt/qudits/simulation/backends/dmsim.py`) は密度行列 $\rho \in \mathbb{C}^{D\times D}$ ($D=81$) を保持し、命令ごとに次の代数を実行する：

#### 6.3.1 ユニタリゲート（`cu_multi`）

`apply_unitary_to_density(ρ, U, qudits, dims)`：

$$
\rho \;\mapsto\; U_S\,\rho\,U_S^\dagger
$$

ここで $U_S$ はターゲット qudit 集合 $S$ 上での $U$、他のサイトは恒等演算子。`cu_multi` で全 qutrit を対象にした場合 $U_S = U_H$ そのもの。

実装は `_apply_local_op_one_side` を用い、テンソル形式 `ρ.reshape(d_0, ..., d_{N-1}, d_0, ..., d_{N-1})` の row 脚に $U$、col 脚に $U^*$ を `np.tensordot` で縮約する。これは数学的に $U\rho U^\dagger$ そのもので、近似なし。

#### 6.3.2 Kraus チャネル（`kraus_channel`）

`apply_kraus_to_density(ρ, {K_α}, qudits, dims)`：

$$
\rho \;\mapsto\; \sum_\alpha (K_\alpha)_S\,\rho\,(K_\alpha^\dagger)_S
$$

実装は各 $K_\alpha$ について row 脚に $K_\alpha$、col 脚に $K_\alpha^*$ を縮約し、結果を accumulate。これも数学的に厳密、近似なし。

#### 6.3.3 ノイズモデルの取り扱い（嘘禁止の設計）

DMSim は `noise_model` オプションを受け取らない設計：

> "DMSim does not accept a NoiseModel option. Add noise as KrausChannel instructions in the circuit instead." と明記され `ValueError` を raise する。

⇒ シナリオ5dにおいてハードウェアノイズは**入っていない**。GKSL 散逸子（$\hat L_\alpha$）は物理ノイズであり、これは Kraus チャネルとして陽に回路に挿入されている。

#### 6.3.4 KrausChannel の `to_matrix` は意図的に NotImplementedError

`KrausChannel.to_matrix()` および `__array__` は

> "KrausChannel has no unitary matrix representation. Use a density-matrix backend (e.g. DMSim) instead."

として `NotImplementedError` を raise する。state-vector backend (TNSim/MISim) で誤って動いてしまうことを禁止する設計（誤魔化し禁止）。

---

## 7. 初期状態（`prepare_initial_state(state_type='edge_triplet')`）

シナリオ5dは Cell 9 で `_initial_state = 'edge_triplet'`：

$$
|\psi_0\rangle \;=\; |1\rangle_0 \otimes |0\rangle_1 \otimes |0\rangle_2 \otimes |1\rangle_3
$$

すなわち両端の分子（0番と $N-1=3$ 番）が $T_1$、中央の分子は $S_0$。コードは pure-state ベクトルから純状態密度行列を作る：

$$
\rho_0 \;=\; |\psi_0\rangle\langle\psi_0|.
$$

基底ラベルは（§2 の慣例で）

$$
\mathrm{idx} \;=\; 1\cdot d^{N-1} + 0\cdot d^{N-2} + 0\cdot d + 1 \;=\; 27 + 0 + 0 + 1 \;=\; 28
$$

であり、`psi[28] = 1.0`、$\rho_0$ は 81×81 行列で第 (28,28) 成分のみが 1、他は 0（純粋状態、$\mathrm{Tr}\rho_0=1$、$\mathrm{Tr}\rho_0^2=1$）。

---

## 8. シミュレーションループ全体

`QuditGKSLSimulator(params, algorithm='exact_local_channels', execute_on_backend='dmsim').simulate(t_max=100.0, n_steps=100, initial_state='edge_triplet')` の正確な処理：

### 8.1 一度限りの前計算（`_precompute_unitaries(dt)`）

1. $U_H = \exp(-i\,\hat H_\mathrm{total}\,dt/2)$ を `scipy.linalg.expm` で計算。
2. すべての $\alpha=1,\dots,26$ について
   - 局所 Lindblad $L_\alpha$ から局所散逸超演算子 $\mathcal{L}_{D,\alpha}^\mathrm{loc}$（§5.2.1）を構築
   - $M_\alpha^{(\mathrm{half})} = \exp(\mathcal{L}_{D,\alpha}^\mathrm{loc}\cdot dt/2)$ を `expm` で計算
   - Choi–Jamiolkowski 同型（§6.1）で Kraus 演算子 $\{K_{\alpha,\beta}\}_\beta$ を厳密に抽出
   - リスト `_dmsim_kraus_half[α] = (sites_α, [K_{α,β}], kind_α)` に格納
3. `MQTQuditProvider().get_backend("dmsim")` で `DMSim` インスタンスを取得し `_dmsim_backend` に保持。

### 8.2 メインループ（`simulate`）

```
ρ ← prepare_initial_state('edge_triplet')         (§7)
times = [0.0]; populations = [pop(ρ)]; ...

for step in 0, 1, ..., n_steps-1:                  (= 100 ステップ)
    ρ ← _trotter_step_dmsim(ρ)                     (§6.2: 回路を構築して DMSim で実行)
    times.append((step+1)*dt)
    traces.append(Re Tr ρ)
    populations.append(pop(ρ))
    entropies.append(S(ρ))                         (§9.2)
    purities.append(Tr ρ²)                         (§9.3)
```

**ループ内で `expm` は呼ばれない**（前計算済 $U_H$ と $K_{\alpha,\beta}$ のみ使用）。

---

## 9. 観測量の計算

### 9.1 集団 $N_{S_0}, N_{T_1}, N_{S_1}$（`compute_populations_from_density_matrix`）

$\rho$ の対角成分 $p(s_0,s_1,s_2,s_3) = \langle s_0 s_1 s_2 s_3|\rho|s_0 s_1 s_2 s_3\rangle$ から

$$
N_{S_0} = \sum_{\vec s} p(\vec s)\cdot\#\{i: s_i=0\},\quad
N_{T_1} = \sum_{\vec s} p(\vec s)\cdot\#\{i: s_i=1\},\quad
N_{S_1} = \sum_{\vec s} p(\vec s)\cdot\#\{i: s_i=2\}.
$$

これは局所演算子 $\hat n_b^{(i)} = |b\rangle_i\langle b|$ の和の期待値 $\sum_i \langle\hat n_b^{(i)}\rangle$ に等しい。$\rho$ がトレース 1 で $N$ 分子なので $N_{S_0}+N_{T_1}+N_{S_1}=N=4$ が常に成立する。

加えて分子毎集団 `per_molecule_populations[mol][b]` も同時に計算（diagonal を base-$d$ digit 展開して累算）。

### 9.2 von Neumann エントロピー（`compute_von_neumann_entropy`）

$$
S(\rho) \;=\; -\mathrm{Tr}(\rho\,\ln\rho) \;=\; -\sum_{\lambda_i > 10^{-12}} \lambda_i\,\ln\lambda_i,
$$

$\lambda_i$ は $\rho$ のエルミート対角化（`np.linalg.eigvalsh`）で得る固有値。閾値 $10^{-12}$ は $0\ln 0=0$ を数値的に安全に扱うためのもの（情報の捨て張りではない）。

### 9.3 純度

$$
\mathcal{P}(\rho) \;=\; \mathrm{Tr}(\rho^2) \;=\; \mathrm{Re}\,\mathrm{Tr}(\rho^2).
$$

純粋状態で 1、最大混合状態で $1/D = 1/81$。

### 9.4 トレース保存性のモニタ

$\mathrm{Tr}\rho(t)$ を毎ステップ記録。Cell 9 の検証で `trace_dev_dmsim = |Tr ρ_final − 1|` を表示。理論的には CPTP の合成なので $\mathrm{Tr}\rho$ は厳密に 1 のはずで、実測される偏差は浮動小数点誤差（$10^{-13}$ オーダ）のみ。

---

## 10. 誤差解析（嘘なしの正確な記述）

シナリオ5d で残存する近似は次の **2 種類のみ**：

1. **Strang H–D 交換子誤差 (ステップ毎 $O(dt^3)$、大域 $O(dt^2)$)**

   $$
   e^{(\mathcal{L}_H+\mathcal{L}_D)\,dt} \;=\; e^{\mathcal{L}_H\,dt/2}\,e^{\mathcal{L}_D\,dt}\,e^{\mathcal{L}_H\,dt/2} \;+\; \frac{dt^3}{24}\Big([\mathcal{L}_H, [\mathcal{L}_H,\mathcal{L}_D]] - 2[\mathcal{L}_D,[\mathcal{L}_D,\mathcal{L}_H]]\Big) + O(dt^5).
   $$

2. **散逸子間のパリンドロミック Lie–Trotter 交換子誤差 (ステップ毎 $O(dt^3)$、大域 $O(dt^2)$)**

   $$
   \prod_{\alpha=1}^{n}\!e^{\mathcal{D}_\alpha\,dt/2}\,\prod_{\alpha=n}^{1}\!e^{\mathcal{D}_\alpha\,dt/2} \;=\; e^{\mathcal{L}_D\,dt} + O(dt^3)\quad(\text{パリンドロミック → 1次項相殺}).
   $$

**シナリオ5（`stinespring`）に存在した「Stinespring による各チャネルの 1 次近似」由来の $O(dt^2)$ ステップ毎誤差（大域 $O(dt)$）は、`exact_local_channels` を使うシナリオ5d では存在しない**。

`tutorials/exact_local_channels.py` ヘッダの実測収束率（`ClassicalGKSLSimulator` 参照との trace distance）：

* $N=2$, $t_\max=10$: rate = 2.000（n_steps = 20…400）
* $N=4$, $t_\max=100$: rate = 2.712 → 2.100 → 2.017 → 2.004（n_steps = 20→50→100→200）

漸近的に $O(dt^2)$ に収束。

その他の誤差は浮動小数点精度のみ（complex128, 相対誤差 $\sim 10^{-15}$）。

### Cell 9 の同値性検証

Cell 9 は `(a) NumPy 直接実行` (`execute_on_backend=None`、§5・§6.3 を NumPy `einsum` で直接行う) と `(b) DMSim バックエンド実行` (§6) の最終 $\rho$ を比較し

$$
\|\rho^\mathrm{ref}_\mathrm{final} - \rho^\mathrm{DMSim}_\mathrm{final}\|_F \;<\; 10^{-10}
$$

を `assert` する。両者は**全く同じ Kraus 演算子 $\{K_{\alpha,\beta}\}$ を全く同じ palindromic 順序で適用しているだけ**であり、違うのは「ループの中で実行するエンジン (NumPy einsum vs MQT-Qudits DMSim)」のみ。よって等値性は数式上自明で、$10^{-10}$ は浮動小数点 round-off を許容する閾値である。これは**新しい物理を主張する検証ではなく**、「DMSim 経由の実行が正しく Kraus を適用しているか」のサニティチェックである。

---

## 11. 容量（capacity）に関する正直な注記

* シナリオ5dの **state-vector backend (`tnsim`/`misim`) での実行は、Stinespring per-step ancilla を fresh にする回路の場合、$N=4$ で根本的に不可能**：1 ステップに system 4 + ancilla 26 = 30 qutrit、state-vector サイズ $3^{30}\approx 2\times 10^{14}$ complex $\approx 3$ PB（メモリ）。これは `STATUS_HONEST_2026-05.md` の A-1 として記録されている既知制約であり、シナリオ5d はこの制約を **density-matrix backend に切り替えることで回避した**もの。state-vector backend 自体の改善ではない。
* シナリオ5d の DMSim 実行は ρ サイズ $81\times 81 \approx 100$ KB なので trivially tractable。
* boson 付きシナリオ (Cell 13, 15) はこの DMSim 経路には対応していない（density-matrix サイズが Hilbert 空間と共に増加するため別途検討要）。
* qubit シミュレータ（シナリオ3 系）は別経路（`tutorials/qiskit_qubit_gksl_simulator.py` + Qiskit Aer `density_matrix`、シナリオ3d）で並列に処理される。
* DMSim は古典計算機上で密度行列を直接演算する simulator backend であり、実量子ハードウェアではない。

---

## 12. 結論：シナリオ5d で実際に解いている方程式

シナリオ5d は次の 4 分子 qutrit GKSL 方程式

$$
\dot{\hat\rho}(t) \;=\; -i\,[\hat H_0 + \hat H_\mathrm{transfer},\,\hat\rho(t)] \;+\; \sum_{\alpha=1}^{26}\Big(\hat L_\alpha\,\hat\rho(t)\,\hat L_\alpha^\dagger - \tfrac12\{\hat L_\alpha^\dagger\hat L_\alpha,\,\hat\rho(t)\}\Big)
$$

（$\hat H_0, \hat H_\mathrm{transfer}, \hat L_\alpha$ は §3, §4 で定義）の数値解を、

* 時間方向に Strang 2次 + パリンドロミック散逸子 2次（合算で大域 $O(dt^2)$）
* 各局所散逸チャネルは Choi–Jamiolkowski による厳密 Kraus 表現で
* 各 Trotter ステップで MQT-Qudits の `QuantumCircuit`（`cu_multi` × 2 + `KrausChannel` × 52）を構築し
* MQT-Qudits 密度行列バックエンド `DMSim` で実行

する形で計算する。**ヒューリスティック処理は皆無**：

* Choi 行列の負の固有値（$> 10^{-12}$）は誤魔化さず `ValueError`
* Kraus の CPTP 違反（$> 10^{-9}$）は誤魔化さず `ValueError`
* `KrausChannel` の `to_matrix` は誤魔化さず `NotImplementedError`
* DMSim は `noise_model` を誤魔化さず `ValueError`
* state-vector backend の容量制約は誤魔化さず A-1 として記録、DMSim 切替で回避と明記

---

## 出典（コード行レベルの参照）

| ファイル | 内容 |
|---|---|
| `tutorials/gksl_physical_parameters.py` | 物理パラメータ定義、validate、neighbors |
| `tutorials/gksl_math_utils.py` | $\hat H_0$, $\hat H_\mathrm{transfer}$, $\hat L_\alpha$ 構築、observables |
| `tutorials/exact_local_channels.py` | 局所散逸超演算子 $\mathcal{L}_D^\mathrm{loc}$、`expm(L_D^loc · dt/2)`、`apply_channel_*` |
| `tutorials/dmsim_kraus_helpers.py` | Choi–Jamiolkowski 同型による厳密 Kraus 抽出（CHOI_EIG_TOL=1e-12、ヒューリスティック無し） |
| `tutorials/qudit_gksl_simulator.py` | `QuditGKSLSimulator`、`_precompute_unitaries`、`_trotter_step_dmsim`、`simulate` |
| `src/mqt/qudits/simulation/backends/dmsim.py` | `DMSim` バックエンド本体、`apply_unitary_to_density`、`apply_kraus_to_density` |
| `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py` | `KrausChannel` 命令、CPTP_TOLERANCE=1e-9、`to_matrix` は NotImplementedError |
| `tutorials/quantum_dynamics_gksl_comparison.ipynb` | Cell 2（パラメータ）、Cell 8/9（シナリオ5d 本体） |

---

# 付録 A: `cu_multi` と `KrausChannel` のゲート表現の完全分解と、もとの qutrit GKSL 方程式から時間発展演算子を構築・適用する導出（省略無し）

本付録は、本文の §5・§6 で「Strang 分割」「パリンドロミック散逸子」「Choi–Jamiolkowski による Kraus 抽出」「DMSim による命令適用」と要約していた箇所を、**式の細部まで分解して**、出発点の GKSL 方程式から最終的な「`cu_multi` と `KrausChannel` の命令列を $\rho$ に掛ける」操作までを一行も省略せずに導出する。記号は本文と同じ。$D=d^N=3^4=81$、$N=4$、$d=3$、$n_\mathrm{Lin}=26$。

## A.1 出発点：qutrit GKSL 方程式とその形式解

物理層（§4）で立てた方程式は

$$
\dot{\hat\rho}(t) \;=\; \mathcal{L}\,\hat\rho(t),\qquad \mathcal{L} \;=\; \mathcal{L}_H + \mathcal{L}_D,
$$

$$
\mathcal{L}_H(\hat\rho) \;=\; -\,i\,\big[\hat H_\mathrm{total},\,\hat\rho\big] \;=\; -\,i\,\hat H_\mathrm{total}\,\hat\rho \;+\; i\,\hat\rho\,\hat H_\mathrm{total},
$$

$$
\mathcal{L}_D(\hat\rho) \;=\; \sum_{\alpha=1}^{n_\mathrm{Lin}} \mathcal{D}[\hat L_\alpha](\hat\rho) \;=\; \sum_{\alpha=1}^{26}\bigg(\hat L_\alpha\,\hat\rho\,\hat L_\alpha^\dagger \;-\; \tfrac12\,\hat L_\alpha^\dagger\hat L_\alpha\,\hat\rho \;-\; \tfrac12\,\hat\rho\,\hat L_\alpha^\dagger\hat L_\alpha\bigg).
$$

$\mathcal{L}$ は密度作用素の空間 $\mathcal{B}(\mathcal{H})$（$D\times D$ 行列の空間）の上の線形写像（**超演算子**）である。$\mathcal{L}$ は時間に陽に依らないので、形式解は

$$
\boxed{\;\hat\rho(t) \;=\; e^{\mathcal{L}\,t}\,\hat\rho(0).\;}\qquad(\star)
$$

ここで $e^{\mathcal{L}\,t}$ は超演算子の指数（$\mathcal{B}(\mathcal{H})\to\mathcal{B}(\mathcal{H})$ の線形写像）である。

**注**：「時間発展演算子」というと閉系でしばしば使う $U(t) = e^{-i\hat H t}$ を指すが、開系 GKSL では時間発展は超演算子 $e^{\mathcal{L}\,t}: \rho\mapsto\rho(t)$ そのものであり、単一の Hilbert 空間ユニタリでは表せない。そのため以下では「時間発展演算子」を**超演算子** $\Lambda(t) := e^{\mathcal{L}\,t}$ の意味で用いる。

## A.2 超演算子の vec 表現（直線化）

$D\times D$ 行列の空間は次元 $D^2$ のベクトル空間と同型である。**列優先 (column-major)** vec 規約

$$
\mathrm{vec}(X)_{i + D\cdot j} \;=\; X_{ij},\qquad i,j\in\{0,\dots,D-1\}
$$

を採用する。次の3つの恒等式は本付録の中心道具である（標準的なテンソル代数）：

$$
\mathrm{vec}(A\,X\,B) \;=\; (B^\top \otimes A)\,\mathrm{vec}(X), \qquad(\dagger_1)
$$

$$
\mathrm{vec}(A\,X) \;=\; (I \otimes A)\,\mathrm{vec}(X), \qquad(\dagger_2)
$$

$$
\mathrm{vec}(X\,B) \;=\; (B^\top \otimes I)\,\mathrm{vec}(X). \qquad(\dagger_3)
$$

これらを用いて

$$
\mathrm{vec}\big(\mathcal{L}_H(\hat\rho)\big) \;=\; \big[-\,i\,(I\otimes \hat H_\mathrm{total}) \;+\; i\,(\hat H_\mathrm{total}^\top\otimes I)\big]\,\mathrm{vec}(\hat\rho),
$$

$$
\mathrm{vec}\big(\mathcal{D}[\hat L_\alpha](\hat\rho)\big) \;=\; \Big[(\hat L_\alpha^*\otimes \hat L_\alpha) \;-\; \tfrac12(I\otimes \hat L_\alpha^\dagger\hat L_\alpha) \;-\; \tfrac12((\hat L_\alpha^\dagger\hat L_\alpha)^\top\otimes I)\Big]\,\mathrm{vec}(\hat\rho).
$$

すなわち $\mathcal{L}_H$ と $\mathcal{D}[L_\alpha]$ は $D^2\times D^2$ 行列として表現できる。$D^2=6561$、$D^4 \approx 4.3\times 10^7$ なので $e^{\mathcal{L}\,t}$ を `expm` で直接取るのは原理上可能だが、**コード（`exact_local_channels`）はあえてこの直接 expm を取らない**。理由は (i) 物理パラメータが変わるたびに大行列の expm を取り直す必要があり高価、(ii) より重要なのは「散逸子部分は局所（最大 9×9 部分空間）」という構造を活かして局所空間でだけ厳密化し、Hamiltonian と分割する方が桁違いに安い、という構造的選択である。

このため「$\mathcal{L}$ の指数を取る」のではなく「$\mathcal{L}_H$ と各 $\mathcal{D}[L_\alpha]$ それぞれの指数を取って合成する」のが分割積分のアイディアであり、その合成の正確な構造が `cu_multi`／`KrausChannel` 命令列に対応する。

## A.3 Strang H–D 分割の導出（BCH による誤差項の明示）

### A.3.1 1次（Lie–Trotter）分割

$$
e^{(\mathcal{L}_H+\mathcal{L}_D)\,dt} \;=\; e^{\mathcal{L}_H\,dt}\,e^{\mathcal{L}_D\,dt}\,\exp\!\Big(-\tfrac{dt^2}{2}[\mathcal{L}_H,\mathcal{L}_D]+O(dt^3)\Big).
$$

これは $e^{A+B} = e^A e^B e^{-\frac12[A,B]+O(\|\cdot\|^3)}$ より。誤差はステップ毎 $O(dt^2)$、大域 $O(dt)$。

### A.3.2 2次（Strang）分割

対称化

$$
e^{(\mathcal{L}_H+\mathcal{L}_D)\,dt} \;=\; e^{\mathcal{L}_H\,dt/2}\;e^{\mathcal{L}_D\,dt}\;e^{\mathcal{L}_H\,dt/2} \;\cdot\; \exp\!\Big(\tfrac{dt^3}{24}\big([\mathcal{L}_H,[\mathcal{L}_H,\mathcal{L}_D]] - 2[\mathcal{L}_D,[\mathcal{L}_D,\mathcal{L}_H]]\big) + O(dt^5)\Big).
$$

導出は BCH を 3 段階に展開するだけ：まず

$$
e^{\mathcal{L}_H\,dt/2}\,e^{\mathcal{L}_D\,dt} \;=\; e^{\mathcal{L}_H\,dt/2 + \mathcal{L}_D\,dt + \tfrac12 \cdot \frac{dt}{2}\cdot dt\,[\mathcal{L}_H,\mathcal{L}_D] + \cdots},
$$

次にこの結果と $e^{\mathcal{L}_H\,dt/2}$ をさらに BCH で合成すると、$[\mathcal{L}_H,\mathcal{L}_D]$ の 1 次項が**前後で逆符号で相殺**し、残るのは 3 次の二重交換子だけになる。誤差はステップ毎 $O(dt^3)$、大域 $O(dt^2)$。これが本文の §5 の Strang 分割。

## A.4 散逸子間のパリンドロミック分割の導出

$\mathcal{L}_D = \sum_{\alpha=1}^{26}\mathcal{D}[\hat L_\alpha]$ をさらに各 $\alpha$ ごとに分割する。$\mathcal{D}_\alpha := \mathcal{D}[\hat L_\alpha]$ と略記。

### A.4.1 forward だけの 1 次分割

$$
e^{\mathcal{L}_D\,(dt/2)} \;=\; \prod_{\alpha=1}^{26} e^{\mathcal{D}_\alpha\,(dt/2)} \;\cdot\; \exp\!\Big(-\tfrac{1}{2}\big(\tfrac{dt}{2}\big)^2\sum_{\alpha<\beta}[\mathcal{D}_\alpha,\mathcal{D}_\beta] + O(dt^3)\Big).
$$

ステップ毎 $O(dt^2)$ 誤差。

### A.4.2 forward + reverse のパリンドロミック分割

forward と reverse を続けて掛けると、

$$
\Big(\prod_{\alpha=1}^{26} e^{\mathcal{D}_\alpha\,(dt/2)}\Big)\Big(\prod_{\alpha=26}^{1} e^{\mathcal{D}_\alpha\,(dt/2)}\Big) \;=\; e^{\mathcal{L}_D\,dt} \;\cdot\; \exp\!\big(O(dt^3)\big),
$$

なぜなら BCH の最低次補正項 $\sum_{\alpha<\beta}[\mathcal{D}_\alpha,\mathcal{D}_\beta]$ が forward と reverse で**符号が反転して相殺**するから。具体的には、$X_\alpha := dt\cdot \mathcal{D}_\alpha/2$ と置いて

$$
\prod_{\alpha=1}^{n}e^{X_\alpha} \;=\; \exp\!\Big(\sum_\alpha X_\alpha + \tfrac12\sum_{\alpha<\beta}[X_\alpha,X_\beta] + O(X^3)\Big),
$$

$$
\prod_{\alpha=n}^{1}e^{X_\alpha} \;=\; \exp\!\Big(\sum_\alpha X_\alpha + \tfrac12\sum_{\alpha>\beta}[X_\alpha,X_\beta] + O(X^3)\Big) \;=\; \exp\!\Big(\sum_\alpha X_\alpha - \tfrac12\sum_{\alpha<\beta}[X_\alpha,X_\beta] + O(X^3)\Big),
$$

両者をさらに BCH で合成すると、$\tfrac12\sum_{\alpha<\beta}[X_\alpha,X_\beta]$ の項が打ち消し、$\exp(2\sum_\alpha X_\alpha + O(X^3)) = \exp(\mathcal{L}_D\,dt + O(dt^3))$ となる。ステップ毎 $O(dt^3)$、大域 $O(dt^2)$。

### A.4.3 結論

$$
\boxed{\;\Phi_\mathrm{pal}(dt) \;:=\; \Big(\bigcirc_{\alpha=1}^{26} e^{\mathcal{D}_\alpha\,dt/2}\Big)\circ\Big(\bigcirc_{\alpha=26}^{1} e^{\mathcal{D}_\alpha\,dt/2}\Big) \;=\; e^{\mathcal{L}_D\,dt} + O(dt^3).\;}
$$

これと A.3.2 を合成して、1 ステップの時間発展超演算子を

$$
\boxed{\;\Phi(dt) \;:=\; e^{\mathcal{L}_H\,dt/2}\,\circ\,\Phi_\mathrm{pal}(dt)\,\circ\,e^{\mathcal{L}_H\,dt/2} \;=\; e^{\mathcal{L}\,dt} + O(dt^3).\;}\qquad(\star\star)
$$

これが**実装が計算しているもの**の正確な定義。

## A.5 各因子の超演算子表現とゲート命令への対応

### A.5.1 ハミルトニアン半ステップ $e^{\mathcal{L}_H\,dt/2}$ ⇒ `cu_multi`

$\hat H_\mathrm{total}$ は時間に依らないので

$$
e^{\mathcal{L}_H\,t}(\hat\rho) \;=\; e^{-i\hat H_\mathrm{total}\,t}\,\hat\rho\,e^{+i\hat H_\mathrm{total}\,t}. \qquad(\heartsuit_1)
$$

導出：$\dot\rho = -i[H,\rho]$ の解は $\rho(t) = U(t)\rho(0)U(t)^\dagger$, $U(t)=e^{-iHt}$。これは Heisenberg–Schrödinger の標準計算で、

$$
\frac{d}{dt}\big[U\rho(0) U^\dagger\big] = -iHU\rho(0)U^\dagger + iU\rho(0)U^\dagger H = -i[H, U\rho(0)U^\dagger].
$$

$t = dt/2$ を代入し、

$$
U_H \;:=\; U_H^{(\mathrm{half})} \;=\; e^{-i\hat H_\mathrm{total}\,dt/2} \;\in\; \mathbb{C}^{81\times 81},
$$

$$
\boxed{\;e^{\mathcal{L}_H\,dt/2}(\hat\rho) \;=\; U_H\,\hat\rho\,U_H^\dagger.\;}\qquad(\heartsuit_2)
$$

vec 表現は $(\dagger_1)$ より

$$
\mathrm{vec}\big(U_H\,\hat\rho\,U_H^\dagger\big) \;=\; \big((U_H^\dagger)^\top \otimes U_H\big)\,\mathrm{vec}(\hat\rho) \;=\; \big(U_H^* \otimes U_H\big)\,\mathrm{vec}(\hat\rho), \qquad(\heartsuit_3)
$$

ここで $(U^\dagger)^\top = U^*$（$U^\dagger$ の転置 = $U$ の複素共役）を用いた。

#### `cu_multi` への対応

MQT-Qudits の API 呼び出し

```
circuit.cu_multi(qudits=[0,1,2,3], parameters=U_H)
```

は `CustomMulti(parent_circuit, name="CUm[3,3,3,3]", target_qudits=[0,1,2,3], parameters=U_H, dimensions=[3,3,3,3], controls=None)` を構築する。`CustomMulti.__init__` は `parameters`（複素行列 $U_H \in \mathbb{C}^{81\times 81}$）を `__array_storage` に保持するだけで、行列に**変更を加えない**。`CustomMulti.__array__` メソッドは `self.__array_storage` をそのまま返す（`src/mqt/qudits/quantum_circuit/gates/custom_multi.py`：43–44 行目）。

DMSim 実行時、命令ループ（`DMSim.execute`：254–255 行目）が `self._apply_instruction(rho, instruction, dims, n)` を呼び、`KrausChannel` でないので unitary path（276 行目）に入り

```
u_matrix = instruction.to_matrix(identities=0)   # → __array__ 経由で U_H
```

を取得し、`apply_unitary_to_density(rho, u_matrix, qudits=(0,1,2,3), dims=[3,3,3,3])`（132–144 行目）が

```
rho_t = rho.reshape(3,3,3,3, 3,3,3,3)             # 8 軸テンソル
rho_t = _apply_local_op_one_side(rho_t, U_H, [0,1,2,3], dims, "row")
rho_t = _apply_local_op_one_side(rho_t, U_H.conj(), [0,1,2,3], dims, "col")
return rho_t.reshape(81, 81)
```

を実行する。`_apply_local_op_one_side` は `np.tensordot` で

$$
\rho_{(s_0 s_1 s_2 s_3)(s_0' s_1' s_2' s_3')}^\text{new} \;=\; \sum_{t_0 t_1 t_2 t_3,\;t_0' t_1' t_2' t_3'} (U_H)_{(s_0 s_1 s_2 s_3),(t_0 t_1 t_2 t_3)}\,\rho_{(t_0 t_1 t_2 t_3)(t_0' t_1' t_2' t_3')}\,(U_H^*)_{(s_0' s_1' s_2' s_3'),(t_0' t_1' t_2' t_3')}
$$

を計算する。これは行列形式に戻すと $\rho^\text{new} = U_H\,\rho\,U_H^\dagger$ そのものである（注意：column 側は $U_H^*$ で、$U_H$ に対し $A B A^\dagger$ の $A^\dagger$ を行列要素で書くと $(A^\dagger)_{s't'} = A^*_{t's'}$、上の縮約の右側の和は $\sum_{t'} \rho_{tt'} (U_H^*)_{s't'} = \sum_{t'}\rho_{tt'}(U_H^\dagger)_{t's'} = (\rho U_H^\dagger)_{ts'}$、よって全体で $(U_H \rho U_H^\dagger)_{ss'}$）。すなわち

$$
\boxed{\;\texttt{cu\_multi}(\,[0,1,2,3],\,U_H\,) \;\;\text{は厳密に}\;\; \rho \mapsto U_H\,\rho\,U_H^\dagger\;\;\text{を実行する。}\;}
$$

これは $(\heartsuit_2)$ と数式上完全に一致し、$e^{\mathcal{L}_H\,dt/2}$ を**厳密に**実装している（近似なし）。

### A.5.2 局所散逸子半ステップ $e^{\mathcal{D}_\alpha\,dt/2}$ ⇒ `KrausChannel`

#### 段階 1：局所空間に縮約する

$\hat L_\alpha$ は最大 2 サイトに非自明に作用する。サポートを $S_\alpha\subseteq\{0,1,2,3\}$（$|S_\alpha|\in\{1,2\}$）、補集合を $\bar S_\alpha$、局所次元を $d_\mathrm{loc} = 3^{|S_\alpha|}$ とすると、$\hat L_\alpha$ は

$$
\hat L_\alpha \;=\; L^\mathrm{loc}_\alpha \,\otimes\, I_{\bar S_\alpha},\qquad L^\mathrm{loc}_\alpha\in\mathbb{C}^{d_\mathrm{loc}\times d_\mathrm{loc}}, \qquad I_{\bar S_\alpha}\in\mathbb{C}^{d^{N-|S_\alpha|}\times d^{N-|S_\alpha|}}
$$

の形に書ける（テンソル積は適切に「サイトを並べ直して」の意味；§A.6 で厳密に書く）。$\mathcal{D}[L\otimes I] = \mathcal{D}[L]\otimes \mathrm{id}$（$\mathrm{id}$ は補空間上の恒等超演算子）が成り立つ：実際、$L = L^\mathrm{loc}\otimes I$ なら $L^\dagger L = (L^\mathrm{loc})^\dagger L^\mathrm{loc}\otimes I$ で、

$$
L\rho L^\dagger - \tfrac12\{L^\dagger L,\rho\} \;=\; (L^\mathrm{loc}\otimes I)\,\rho\,((L^\mathrm{loc})^\dagger\otimes I) - \tfrac12 \big\{(L^\mathrm{loc})^\dagger L^\mathrm{loc}\otimes I,\,\rho\big\},
$$

これは局所超演算子 $\mathcal{D}^\mathrm{loc}[L^\mathrm{loc}_\alpha]$ を $S_\alpha$ サイトに、補空間に恒等を作用させたものに等しい。よって**指数化も局所と恒等のテンソルになる**：

$$
e^{\mathcal{D}_\alpha\,dt/2} \;=\; e^{\mathcal{D}^\mathrm{loc}[L^\mathrm{loc}_\alpha]\,dt/2} \,\otimes\, \mathrm{id}_{\bar S_\alpha}. \qquad(\spadesuit_1)
$$

これが「散逸子の指数を局所空間でだけ取れば良い」ことの厳密な根拠（コードコメントの「local」の意味）。

#### 段階 2：局所超演算子の vec 表現を作る（§5.2.1 の再確認）

$L = L^\mathrm{loc}_\alpha \in \mathbb{C}^{d_\mathrm{loc}\times d_\mathrm{loc}}$ に対し $(\dagger_1)$–$(\dagger_3)$ を当てて

$$
\mathcal{L}^\mathrm{loc}_\alpha \;=\; (L^*\otimes L) \;-\; \tfrac12\,(I\otimes L^\dagger L) \;-\; \tfrac12\,((L^\dagger L)^\top\otimes I) \;\in\;\mathbb{C}^{d_\mathrm{loc}^2\times d_\mathrm{loc}^2}. \qquad(\spadesuit_2)
$$

#### 段階 3：局所半ステップ超演算子を厳密に取る

```
M_alpha_half = scipy.linalg.expm(L^loc_alpha * dt/2)        # 9×9 または 81×81
```

これは数値的に**厳密**（`expm` は scaling-and-squaring + Padé；二重精度 round-off 以外の近似なし）。

#### 段階 4：Choi–Jamiolkowski 同型で Kraus 演算子を抽出

CPTP 写像 $\mathcal{E}: X\mapsto \mathcal{E}(X)$、$X\in\mathbb{C}^{d_\mathrm{loc}\times d_\mathrm{loc}}$ の Choi 行列は

$$
C(\mathcal{E}) \;=\; (\mathcal{E}\otimes\mathrm{id})\big(|\Omega\rangle\langle\Omega|\big),\qquad |\Omega\rangle = \sum_{i=0}^{d_\mathrm{loc}-1}|i\rangle\otimes|i\rangle,
$$

成分表示は $C_{(ik)(jl)} = \mathcal{E}(|i\rangle\langle j|)_{kl}$（インデックスの組 $(ik)$ を行、$(jl)$ を列）。ここで $\mathcal{E}(|i\rangle\langle j|)_{kl} = T_{kilj}$ と置けば（$T$ は $\mathcal{E}$ のテンソル成分、$\mathrm{vec}$ 規約に整合する形で）、コード中の

```
t = M_half.reshape(d, d, d, d)             # M_half[k+d*l, i+d*j] = T_{klij}
choi = t.transpose(0, 2, 1, 3).reshape(d*d, d*d)
```

がまさに $C_{(ki)(lj)} = T_{klij}$ を生成する。Choi 行列は CP 性により**半正定値**であり、固有分解

$$
C \;=\; \sum_{\beta=1}^{d_\mathrm{loc}^2} \lambda_\beta\,|v_\beta\rangle\langle v_\beta|,\qquad \lambda_\beta \in\mathbb{R}.
$$

物理性チェック：$\lambda_\beta < -10^{-12}$ なら CP 性違反 ⇒ `ValueError`（嘘禁止）。$|\lambda_\beta| < 10^{-12}$ は数値ノイズとして 0 扱い。

各 $\lambda_\beta \geq 10^{-12}$ について、固有ベクトル $v_\beta\in\mathbb{C}^{d_\mathrm{loc}^2}$ を行優先 reshape で

$$
K_{\alpha,\beta} \;:=\; \sqrt{\lambda_\beta}\,\mathrm{reshape}(v_\beta,\,(d_\mathrm{loc},\,d_\mathrm{loc})),\qquad (K_{\alpha,\beta})_{kl} \;=\; \sqrt{\lambda_\beta}\,(v_\beta)_{k\,d_\mathrm{loc} + l}. \qquad(\spadesuit_3)
$$

このとき**任意の** $\rho_\mathrm{loc}\in\mathbb{C}^{d_\mathrm{loc}\times d_\mathrm{loc}}$ について

$$
\sum_\beta K_{\alpha,\beta}\,\rho_\mathrm{loc}\,K_{\alpha,\beta}^\dagger \;=\; \mathcal{E}(\rho_\mathrm{loc}) \;=\; e^{\mathcal{D}^\mathrm{loc}[L^\mathrm{loc}_\alpha]\,dt/2}(\rho_\mathrm{loc}). \qquad(\spadesuit_4)
$$

導出（細部）：$\sum_\beta (K_\beta)_{kl}\,(\rho)_{lm}\,(K_\beta^\dagger)_{mn} = \sum_\beta\sum_{lm}\lambda_\beta (v_\beta)_{kl}(v_\beta^*)_{nm}\rho_{lm}$。一方 $C_{(kn)(lm)} = \sum_\beta \lambda_\beta (v_\beta)_{(kn)}(v_\beta^*)_{(lm)} = \sum_\beta\lambda_\beta (v_\beta)_{k\cdot d+n}(v_\beta^*)_{l\cdot d+m}$（行優先 vec で読み直す）。$C_{(ki)(lj)} = T_{klij} = \mathcal{E}(|i\rangle\langle j|)_{kl}$ より $\mathcal{E}(\rho)_{kl} = \sum_{ij}T_{klij}\rho_{ij} = \sum_{ij}C_{(ki)(lj)}\rho_{ij}$。インデックスを揃えると上式と一致。

#### 完全性条件 $\sum_\beta K_{\alpha,\beta}^\dagger K_{\alpha,\beta} = I_{d_\mathrm{loc}}$

$\mathcal{E}$ がトレース保存 ⇔ $C$ の部分トレース $\mathrm{Tr}_1 C = I$、これは Kraus 表現で $\sum_\beta K^\dagger K = I$ と同値。`KrausChannel.__init__` がこれを Frobenius ノルム $10^{-9}$ で**実測**検証（`src/.../kraus_channel.py`：107–119 行目）。違反時は `ValueError`（嘘禁止）。

#### 段階 5：局所 Kraus を全体空間に持ち上げる（テンソル化の細部）

$K_{\alpha,\beta} \in\mathbb{C}^{d_\mathrm{loc}\times d_\mathrm{loc}}$ を、サポート $S_\alpha\subseteq\{0,1,2,3\}$ 上のサイトに作用し、補空間 $\bar S_\alpha$ には恒等で作用する $D\times D$ 行列 $\tilde K_{\alpha,\beta} := (K_{\alpha,\beta})_{S_\alpha}$ に持ち上げる。具体的に $S_\alpha = \{0,1\}$（pair channel TTA01）の場合、

$$
\tilde K_{\alpha,\beta} \;=\; K_{\alpha,\beta}\,\otimes\, I_3\,\otimes\, I_3,
$$

ここで一つ目の Kronecker 因子は qudit 0 と qudit 1 を合体した 9 次元空間に作用、残りは qudit 2、qudit 3 にそれぞれ作用。$S_\alpha = \{1,2\}$ の場合は

$$
\tilde K_{\alpha,\beta} \;=\; I_3\otimes K_{\alpha,\beta}\otimes I_3,
$$

など。$S_\alpha=\{i\}$ の単一サイトの場合は $\tilde K_{\alpha,\beta} = I^{\otimes i}\otimes K_{\alpha,\beta}\otimes I^{\otimes (3-i)}$。

**コードでは陽には Kronecker 積を作らない**：DMSim の `apply_kraus_to_density`（147–162 行目）は $\rho$ をテンソル `rho.reshape(3,3,3,3, 3,3,3,3)` に変形し、`np.tensordot` で $K$ をターゲット脚にだけ縮約する。これは数学的に $\tilde K = K\otimes I_{\bar S}$ を作って $\tilde K \rho \tilde K^\dagger$ を計算するのと厳密に同じ（テンソル代数の標準恒等式 $(\mathrm{id}_{\bar S}\otimes K)\,\mathrm{vec}(\rho) = ...$）。

#### `kraus_channel` への対応

```
circuit.kraus_channel(qudits=S_alpha, kraus_operators=[K_{alpha,1}, ..., K_{alpha,r_alpha}])
```

は `KrausChannel(parent_circuit, name="Kraus[...]", target_qudits=S_alpha, kraus_operators=[K_{α,β}], dimensions=[3,...])` を構築する（`circuit.py`：211–230 行目）。`KrausChannel.__init__` は (i) 各 $K$ の形状 $(d_\mathrm{loc}, d_\mathrm{loc})$ を検証、(ii) $\sum_\beta K^\dagger K = I$ を Frobenius $10^{-9}$ で検証、(iii) Kraus 列を `_kraus_operators` に保持（`kraus_channel.py`：96–123 行目）。`to_matrix()` および `__array__` は `NotImplementedError` を raise（148–159 行目）：「ユニタリ行列としては存在しない」が型レベルで保証される。

DMSim 実行時、`_apply_instruction` で `isinstance(instruction, KrausChannel)` 分岐に入り（`dmsim.py`：272–273 行目）、

```
return apply_kraus_to_density(rho, instruction.kraus_operators, qudits=S_alpha, dims=[3,3,3,3])
```

`apply_kraus_to_density`（147–162 行目）は

```
acc = 0
for k in [K_{alpha,1}, ..., K_{alpha,r_alpha}]:
    rho_t = rho.reshape(3,3,3,3, 3,3,3,3)
    rho_t = _apply_local_op_one_side(rho_t, k, S_alpha, dims, "row")
    rho_t = _apply_local_op_one_side(rho_t, k.conj(), S_alpha, dims, "col")
    acc += rho_t.reshape(81,81)
return acc
```

これは数式で書くと

$$
\rho^\text{new} \;=\; \sum_{\beta=1}^{r_\alpha} \tilde K_{\alpha,\beta}\,\rho\,\tilde K_{\alpha,\beta}^\dagger \;=\; (e^{\mathcal{D}_\alpha\,dt/2})(\rho), \qquad(\spadesuit_5)
$$

となり、$(\spadesuit_1)$–$(\spadesuit_4)$ より $e^{\mathcal{D}_\alpha\,dt/2}$ を**厳密に**実装している。要するに

$$
\boxed{\;\texttt{kraus\_channel}(\,S_\alpha,\,\{K_{\alpha,\beta}\}\,) \;\;\text{は厳密に}\;\; \rho\mapsto \sum_\beta \tilde K_{\alpha,\beta}\,\rho\,\tilde K_{\alpha,\beta}^\dagger\;\;\text{を実行する。}\;}
$$

## A.6 1 ステップの時間発展演算子の完全な表式

$(\heartsuit_2)$ と $(\spadesuit_5)$ を $(\star\star)$ に代入し、命令列の順序（`_trotter_step_dmsim` が append する順）を保つと、1 ステップの時間発展超演算子 $\Phi(dt)$ を $\rho$ に作用させた結果は

$$
\Phi(dt)(\rho) \;=\; \mathcal{U}_H\!\bigg(\;\sum_{\beta_{1}^{(L)}}\!\tilde K_{1,\beta_1^{(L)}}\!\cdots\sum_{\beta_{26}^{(L)}}\!\tilde K_{26,\beta_{26}^{(L)}}\;\sum_{\beta_{26}^{(R)}}\!\tilde K_{26,\beta_{26}^{(R)}}\!\cdots\sum_{\beta_{1}^{(R)}}\!\tilde K_{1,\beta_1^{(R)}}\;\mathcal{U}_H(\rho)\;\tilde K_{1,\beta_1^{(R)}}^\dagger\!\cdots\tilde K_{26,\beta_{26}^{(R)}}^\dagger\;\tilde K_{26,\beta_{26}^{(L)}}^\dagger\!\cdots \tilde K_{1,\beta_1^{(L)}}^\dagger\;\bigg)
$$

の形に書ける。ここで $\mathcal{U}_H(\rho) = U_H\,\rho\,U_H^\dagger$、forward 列の Kraus index $\beta_\alpha^{(R)}$ と reverse 列の Kraus index $\beta_\alpha^{(L)}$ は**独立**に和を取る（各 `KrausChannel` は独立な CPTP 写像）。**書き直すと**、外側に左から右に「ハミルトニアン半 → forward 26 個 → reverse 26 個 → ハミルトニアン半」という超演算子の合成

$$
\Phi(dt) \;=\; \mathcal{U}_H\,\circ\,\mathcal{E}_1^{(R)}\circ\cdots\circ\mathcal{E}_{26}^{(R)}\,\circ\,\mathcal{E}_{26}^{(L)}\circ\cdots\circ\mathcal{E}_1^{(L)}\,\circ\,\mathcal{U}_H
$$

ただし $\mathcal{E}_\alpha^{(R)} = \mathcal{E}_\alpha^{(L)} = e^{\mathcal{D}_\alpha\,dt/2}$（同じ局所超演算子；上付き $(R)$/$(L)$ は forward/reverse の位置を示すラベル）。

注意：本文 §6.2.2 の式

$$
\Phi(dt) \;=\; \mathcal{U}_H\circ\Big(\bigcirc_{\alpha=1}^{26}\mathcal{E}_\alpha\Big)\circ\Big(\bigcirc_{\alpha=26}^{1}\mathcal{E}_\alpha\Big)\circ\mathcal{U}_H
$$

と完全に等価（$\circ$ の合成順序の規約：右の超演算子から先に作用する）。$\bigcirc_{\alpha=1}^{26}$ は $\mathcal{E}_{26}\circ\mathcal{E}_{25}\circ\cdots\circ\mathcal{E}_1$（左から右に書いた順、すなわち $\mathcal{E}_1$ が最初に適用される）の意味。

### A.6.1 vec 表現での 1 ステップ

$(\heartsuit_3)$ と $(\spadesuit_5)$ の vec 表現

$$
\mathrm{vec}\big(\sum_\beta \tilde K_{\alpha,\beta}\,\rho\,\tilde K_{\alpha,\beta}^\dagger\big) \;=\; \Big[\sum_\beta \tilde K_{\alpha,\beta}^*\otimes\tilde K_{\alpha,\beta}\Big]\,\mathrm{vec}(\rho)
$$

を組み合わせると、$\Phi(dt)$ の vec 行列表現（$D^2\times D^2 = 6561\times 6561$）は

$$
\boxed{\;\Phi(dt) \;\stackrel{\text{vec}}{=}\; (U_H^*\otimes U_H)\;\Big[\prod_{\alpha=1}^{26}\sum_\beta \tilde K_{\alpha,\beta}^*\otimes\tilde K_{\alpha,\beta}\Big]\;\Big[\prod_{\alpha=26}^{1}\sum_\beta \tilde K_{\alpha,\beta}^*\otimes\tilde K_{\alpha,\beta}\Big]\;(U_H^*\otimes U_H).\;}
$$

ただし行列積の左側が後から作用する。**コードはこの 6561×6561 行列を陽には作らない**（メモリと計算量のため）：上で示した通り、$\rho$ を $(3,3,3,3,3,3,3,3)$ テンソルとして保持し、各因子をテンソル脚に対して `tensordot` で逐次適用する。これは数学的に**厳密に同値**で、近似なし。

## A.7 全シミュレーションの時間発展演算子（100 ステップ）

メインループは

```
for n = 0, 1, ..., 99:
    rho_{n+1} = Phi(dt)(rho_n)        ← 1 step circuit を DMSim で実行
```

なので

$$
\hat\rho(t_n) \;=\; \rho_n \;=\; \Phi(dt)^n(\rho_0),\qquad t_n = n\cdot dt = n. \qquad(\clubsuit_1)
$$

特に最終時刻 $t_{100}=100$ で

$$
\boxed{\;\hat\rho_\mathrm{final} \;=\; \rho_{100} \;=\; \big[\Phi(dt)\big]^{100}(\rho_0) \;=\; \underbrace{\Phi(dt)\circ\Phi(dt)\circ\cdots\circ\Phi(dt)}_{100\text{ 回}}(\rho_0).\;}\qquad(\clubsuit_2)
$$

これと $(\star)$（GKSL の真の解 $\hat\rho^\mathrm{true}(100) = e^{100\,\mathcal{L}}\hat\rho_0$）との関係は

$$
\Phi(dt)^{100} \;=\; \big(e^{\mathcal{L}\,dt} + O(dt^3)\big)^{100} \;=\; e^{\mathcal{L}\cdot 100\cdot dt} + 100\cdot O(dt^3) \;=\; e^{100\,\mathcal{L}} + O(t_\mathrm{max}\cdot dt^2)
$$

（ステップ毎 $O(dt^3)$ × ステップ数 $t_\mathrm{max}/dt$ ＝ 大域 $O(dt^2)$、$dt=1$ なので大域誤差は $O(1)\cdot t_\mathrm{max}$ オーダ；実測収束は §10 の通り `n_steps→∞` で 2 次収束）。

### A.7.1 各ステップを「ゲート命令の積」に展開した最終形（省略無し）

$\rho_n$ から $\rho_{n+1}$ への 1 ステップを、ゲート命令毎に$\rho$ がどう変化するかを**全 54 行**の代入として書き下すと（$\rho^{(k)}$ は $k$ 番目の命令適用後の状態）：

$$
\begin{aligned}
\rho^{(0)}  &:= \rho_n \\[2pt]
\rho^{(1)}  &= U_H\,\rho^{(0)}\,U_H^\dagger &&\text{(命令 a: cu\_multi[0,1,2,3])}\\[2pt]
\rho^{(2)}  &= \sum_{\beta} \tilde K_{1,\beta}\,\rho^{(1)}\,\tilde K_{1,\beta}^\dagger &&\text{(命令 b1: kraus\_channel[0,1], TTA01-1)}\\
\rho^{(3)}  &= \sum_{\beta} \tilde K_{2,\beta}\,\rho^{(2)}\,\tilde K_{2,\beta}^\dagger &&\text{(命令 b2: kraus\_channel[0,1], TTA01-2)}\\
\rho^{(4)}  &= \sum_{\beta} \tilde K_{3,\beta}\,\rho^{(3)}\,\tilde K_{3,\beta}^\dagger &&\text{(命令 b3: kraus\_channel[1,2], TTA12-1)}\\
\rho^{(5)}  &= \sum_{\beta} \tilde K_{4,\beta}\,\rho^{(4)}\,\tilde K_{4,\beta}^\dagger &&\text{(命令 b4: kraus\_channel[1,2], TTA12-2)}\\
\rho^{(6)}  &= \sum_{\beta} \tilde K_{5,\beta}\,\rho^{(5)}\,\tilde K_{5,\beta}^\dagger &&\text{(命令 b5: kraus\_channel[2,3], TTA23-1)}\\
\rho^{(7)}  &= \sum_{\beta} \tilde K_{6,\beta}\,\rho^{(6)}\,\tilde K_{6,\beta}^\dagger &&\text{(命令 b6: kraus\_channel[2,3], TTA23-2)}\\
\rho^{(8)}  &= \sum_{\beta} \tilde K_{7,\beta}\,\rho^{(7)}\,\tilde K_{7,\beta}^\dagger &&\text{(命令 b7: kraus\_channel[0], 蛍光-0)}\\
\rho^{(9)}  &= \sum_{\beta} \tilde K_{8,\beta}\,\rho^{(8)}\,\tilde K_{8,\beta}^\dagger &&\text{(命令 b8: kraus\_channel[1], 蛍光-1)}\\
\rho^{(10)} &= \sum_{\beta} \tilde K_{9,\beta}\,\rho^{(9)}\,\tilde K_{9,\beta}^\dagger &&\text{(命令 b9: kraus\_channel[2], 蛍光-2)}\\
\rho^{(11)} &= \sum_{\beta} \tilde K_{10,\beta}\,\rho^{(10)}\,\tilde K_{10,\beta}^\dagger &&\text{(命令 b10: kraus\_channel[3], 蛍光-3)}\\
&\;\;\vdots &&\;\;\vdots\\
\rho^{(27)} &= \sum_{\beta} \tilde K_{26,\beta}\,\rho^{(26)}\,\tilde K_{26,\beta}^\dagger &&\text{(命令 b26: kraus\_channel[3], ISC$_{T\to S}$-3)}\\[3pt]
\rho^{(28)} &= \sum_{\beta} \tilde K_{26,\beta}\,\rho^{(27)}\,\tilde K_{26,\beta}^\dagger &&\text{(命令 c1 = b26 再適用)}\\
&\;\;\vdots &&\;\;\vdots\\
\rho^{(53)} &= \sum_{\beta} \tilde K_{1,\beta}\,\rho^{(52)}\,\tilde K_{1,\beta}^\dagger &&\text{(命令 c26 = b1 再適用)}\\[3pt]
\rho^{(54)} &= U_H\,\rho^{(53)}\,U_H^\dagger &&\text{(命令 d: cu\_multi[0,1,2,3])}\\[2pt]
\rho_{n+1} &:= \rho^{(54)} &&
\end{aligned}
$$

ここで $\tilde K_{\alpha,\beta}$ は §A.5.2 段階 5 で定義した「局所 Kraus を $S_\alpha$ サイトに作用させ、$\bar S_\alpha$ には恒等で作用させた $81\times 81$ 行列」。$\beta$ 和の上限 $r_\alpha$ は省略しているが、各 $\alpha$ ごとに固定（pair なら $\leq 81$、single なら $\leq 9$）。

DMSim の実装は **$\tilde K$ を陽に組まずに** 上の代入を `np.tensordot` で実行する（A.5.2 段階 5 の通り）が、数学的にはこの 54 行が DMSim が行っている全部であり、これ以上の隠れた処理は無い（コードを `_apply_instruction` まで辿れば確認できる）。

### A.7.2 100 ステップの完全展開

A.7.1 の代入を 100 回繰り返したものが $\rho_{100}$ であり、これは $\rho_0$ から開始して合計 $54\times 100 = 5400$ 個の超演算子作用（各々が 1 個の `cu_multi` か 1 個の `kraus_channel`）を順に施すことに等しい。

## A.8 まとめ：もとの GKSL から `cu_multi`／`KrausChannel` 列までの導出経路（一行不漏れ）

```
[1] 物理層         dot ρ = -i [H_total, ρ] + Σ_α D[L_α](ρ),   ρ(0) = |ψ_0><ψ_0|
        ↓ (形式解、超演算子の指数)
[2] 形式解         ρ(t) = e^{L · t}(ρ_0),  L = L_H + L_D
        ↓ (Strang 2 次分割: A.3)
[3] H–D 分割       e^{L · dt} = e^{L_H · dt/2} · e^{L_D · dt} · e^{L_H · dt/2}  + O(dt^3)
        ↓ (パリンドロミック散逸子分割: A.4)
[4] 散逸子分割     e^{L_D · dt} = (∏↑ e^{D_α · dt/2})·(∏↓ e^{D_α · dt/2})  + O(dt^3)
        ↓ (各 D_α の局所性: A.5.2 段階 1)
[5] 局所化         e^{D_α · dt/2} = e^{D^loc[L^loc_α] · dt/2} ⊗ id_{barS_α}
        ↓ (vec 表現: A.5.2 段階 2)
[6] 局所超演算子   L^loc_α = L*⊗L − ½(I⊗L†L) − ½((L†L)^T⊗I)   (L = L^loc_α)
        ↓ (scipy.linalg.expm: A.5.2 段階 3)
[7] 局所半 ch.     M^half_α = expm(L^loc_α · dt/2)         (9×9 または 81×81)
        ↓ (Choi–Jamiolkowski 同型: A.5.2 段階 4)
[8] Choi 行列      C_α = reshape&transpose(M^half_α);  C_α = Σ_β λ_β |v_β><v_β|
        ↓ (CPTP 検証 + 非物理は ValueError)
[9] Kraus          K_{α,β} = sqrt(λ_β) · reshape(v_β, (d_loc, d_loc))
        ↓ (テンソル化: A.5.2 段階 5)
[10] 全空間 Kraus  tilde{K}_{α,β} = K_{α,β} ⊗ I_{barS_α}
        ↓ (MQT-Qudits API)
[11] 命令          circuit.cu_multi([0,1,2,3], U_H)              ← e^{L_H · dt/2} を厳密実装
                  circuit.kraus_channel(S_α, [K_{α,β}])         ← e^{D_α · dt/2} を厳密実装
        ↓ (1 ステップ回路の構築: A.6, A.7.1)
[12] 1 step 回路   54 命令: a, b1..b26, c1..c26, d
        ↓ (DMSim 実行: A.5.1, A.5.2 段階 5)
[13] DMSim         各命令に対し ρ <- U_H ρ U_H†   または   ρ <- Σ_β tilde{K}_{α,β} ρ tilde{K}_{α,β}†
        ↓ (100 step ループ: A.7.2)
[14] 最終          ρ_final = Φ(dt)^{100}(ρ_0) = e^{100·L}(ρ_0) + O(t_max · dt^2)
```

これがシナリオ5dで「もとの qutrit GKSL 方程式の時間発展演算子を `cu_multi` と `KrausChannel` でどう構築・適用しているか」の**省略無しの完全な導出鎖**である。各段階で行われている近似は段階 [3]・[4] の Strang 分割／パリンドロミック分割（合算で大域 $O(dt^2)$）のみであり、それ以外は**厳密**（vec 同型は数学的恒等式、`expm` は機械精度、Choi–Jamiolkowski は同型、Kraus 適用は数学的恒等式、命令列の合成は超演算子合成と一致）。ヒューリスティック処理は段階 [8]・[9]・[11] のいずれにも入っておらず、CPTP 違反や負固有値はすべて `ValueError` で停止する設計（`KrausChannel.__init__` の CPTP_TOLERANCE=1e-9、`kraus_from_local_superoperator` の CHOI_EIG_TOL=1e-12）。

