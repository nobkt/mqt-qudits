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

各 Trotter ステップで、$N=4$ qutrit の `QuantumCircuit`（量子レジスタ次元 $[3,3,3,3]$）を**新規構築**して `DMSim` で実行する：

```
[システム量子レジスタ：4 qutrit]

  cu_multi(全4 qutrit, U_H_half)                           ← (a) ハミルトニアン半ステップ

  for (sites, K_α^(half), kind) in _dmsim_kraus_half:      ← (b) 散逸子 forward 列
      kraus_channel(target=sites, K_α^(half))                  （26個、構築順）

  for (sites, K_α^(half), kind) in reversed(_dmsim_kraus_half):← (c) 散逸子 reverse 列
      kraus_channel(target=sites, K_α^(half))                  （26個、逆順）

  cu_multi(全4 qutrit, U_H_half)                           ← (d) ハミルトニアン半ステップ
```

合計命令数 ＝ $2 + 26 + 26 + 2 = 56$。各 `KrausChannel` の Kraus 演算子は §6.1 で構築した $\{K_\alpha\}$ そのもの。`cu_multi` は `mqt.qudits.quantum_circuit.gates.CustomMulti` でユニタリ $U_H = \exp(-i\hat H_\mathrm{total}\,dt/2)$ を全 4 qutrit に作用させる命令。

実行は

```
job = self._dmsim_backend.run(circuit, initial_density_matrix=ρ_prev)
ρ_next = job.result().get_density_matrix()
```

ここで `initial_density_matrix` は前ステップの $\rho$（最初のステップは §7 の初期状態）。

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
