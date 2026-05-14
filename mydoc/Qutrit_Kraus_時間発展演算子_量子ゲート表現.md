# TTA と蛍光の場合における Qutrit-Kraus 表現の時間発展演算子：量子ゲート表現の完全定式化

> **真実ベース宣言**
> 本文書は `tutorials/` の実コード（Python）を直接読んだ事実のみを記述する。
> 省略や誇張は一切しない。数式に含まれない仮定は「注意」として明示する。

---

## 0. 本文書の対象と既存文書との関係

| 文書 | 対象 |
|------|------|
| `シナリオ5d_TTAと蛍光のみ_完全展開.md` | 超演算子レベル・Kraus 解析形（§6–§8） |
| `シナリオ5d_TTAと蛍光のみ_Qudit量子ゲート展開.md` | ネイティブ qutrit ゲート分解（H₀, H_transfer, Stinespring） |
| **本文書** | **コード上に存在する 2 種類の時間発展実装それぞれについて、Kraus 演算子の定義から量子ゲート（命令）への完全な写像を省略なしに書く** |

本文書は **TTA と蛍光の 2 種類の Lindblad チャンネル** に対象を絞る（実コードの 26 チャンネルの縮約版；ただし行列次元・量子ゲート種別は実コードと一致させる）。

コード上に存在する **2 つのチャンネル適用アルゴリズム** を並列に扱う：

| アルゴリズム | コード変数 | 精度 | アンシラ | コード箇所 |
|---|---|---|---|---|
| **Stinespring** | `algorithm="stinespring"` | O(dt)（1 次近似） | 必要（qutrit 1 個） | `stinespring_utils.py::stinespring_unitary_from_lindblad` |
| **完全局所チャンネル** | `algorithm="exact_local_channels"` | O(dt²) | 不要（Kraus 直接適用） | `exact_local_channels.py::precompute_exact_channels_half`, `dmsim_kraus_helpers.py::kraus_from_local_superoperator` |

---

## 1. ヒルベルト空間・パラメータ

### 1.1 qutrit 基底

各分子は 3 準位（qutrit、`d=3` は `params.d != 3` で `ValueError` が出る）：

```
|0⟩ = S₀（基底一重項）
|1⟩ = T₁（三重項励起）
|2⟩ = S₁（一重項励起）
```

コード（`gksl_physical_parameters.py:105-107`）：
```python
# d=3 required: ground (S0), triplet (T1), singlet (S1) states per molecule
if self.d != 3:
    errors.append("d must be 3")
```

### 1.2 デフォルト物理パラメータ（`gksl_physical_parameters.py`）

| 変数 | 値 | 意味 |
|------|-----|------|
| `E_T` | 1.5 eV | T₁ エネルギー |
| `E_S` | 3.0 eV | S₁ エネルギー |
| `V` | 0.1 eV | 近接分子間移動積分 |
| `gamma_TTA` | 0.05 eV/ℏ | TTA 速度定数 |
| `Gamma_fl` | 0.01 eV/ℏ | 蛍光速度定数 |
| `N_molecules` | 4 | 分子数 |
| `d_anc` | 3 | アンシラ qutrit 次元 |

---

## 2. Lindblad 演算子の明示的行列表現

### 2.1 記法

$|a\rangle\langle b|$ は `_ket_bra(d, a, b)` 関数の返り値（`gksl_math_utils.py:78-81`）：

$$
(|a\rangle\langle b|)_{ij} = \delta_{ia}\,\delta_{jb},
$$

すなわち行 $a$、列 $b$ のみ 1、他は 0 の $3\times 3$ 行列。

### 2.2 TTA Lindblad 演算子（ペア演算子、9×9）

コード（`gksl_math_utils.py:84-96`）：

```python
for i, j in params.neighbors:
    gamma = params.gamma_TTA / 2.0
    sq = np.sqrt(gamma)
    # Channel 1: |2>_i<1| x |0>_j<1|
    ops.append((sq * reduce(np.kron, op_list), gamma))
    # Channel 2: |0>_i<1| x |2>_j<1|
    ops.append((sq * reduce(np.kron, op_list), gamma))
```

隣接ペア $\{(0,1),(1,2),(2,3)\}$ に対し、各ペア $(i,j)$ に 2 チャンネル存在：

$$
L^{(i,j)}_1 = \sqrt{\frac{\gamma_\mathrm{TTA}}{2}}\;(|2\rangle_i\langle 1|)\otimes(|0\rangle_j\langle 1|), \qquad 9\times 9
$$

$$
L^{(i,j)}_2 = \sqrt{\frac{\gamma_\mathrm{TTA}}{2}}\;(|0\rangle_i\langle 1|)\otimes(|2\rangle_j\langle 1|), \qquad 9\times 9
$$

**成分表示（ペア $(i,j)$ の局所 Hilbert 空間内、行列要素 $(|s_i s_j\rangle, |s_i' s_j'\rangle)$）**：

$|a\rangle\langle b|$ の $3\times 3$ 成分：

$$
|2\rangle\langle 1| = \begin{pmatrix}0&0&0\\0&0&0\\0&1&0\end{pmatrix},\quad
|0\rangle\langle 1| = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

チャンネル 1 の $9\times 9$ 行列（Kronecker 積、添字順は $|s_i,s_j\rangle = |00\rangle,|01\rangle,|02\rangle,|10\rangle,|11\rangle,|12\rangle,|20\rangle,|21\rangle,|22\rangle$）：

$$
|2\rangle\langle 1|\otimes|0\rangle\langle 1| = \begin{pmatrix}
0&0&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0\\
0&1&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0\\
0&0&0&0&0&0&0&0&0
\end{pmatrix}
$$

行 $|20\rangle$、列 $|11\rangle$ のみ 1。$L^{(i,j)}_1 = \sqrt{\gamma_\mathrm{TTA}/2}\cdot(|2\rangle\langle 1|\otimes|0\rangle\langle 1|)$。

チャンネル 2：$L^{(i,j)}_2 = \sqrt{\gamma_\mathrm{TTA}/2}\cdot(|0\rangle\langle 1|\otimes|2\rangle\langle 1|)$（行 $|02\rangle$、列 $|11\rangle$ のみ）。

### 2.3 蛍光 Lindblad 演算子（単一サイト、3×3）

コード（`gksl_math_utils.py:98-102`）：

```python
for i in range(N):
    gamma = params.Gamma_fl
    L = np.sqrt(gamma) * build_single_site_operator(_ket_bra(0, 2), i, N, d)
    ops.append((L, gamma))
```

各サイト $i$ に対して：

$$
L^{(i)}_\mathrm{fl} = \sqrt{\Gamma_\mathrm{fl}}\;|0\rangle_i\langle 2|, \qquad 3\times 3\;\text{（局所）}
$$

$$
|0\rangle\langle 2| = \begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}
$$

すなわち、行 $|0\rangle$（S₀）、列 $|2\rangle$（S₁）のみ 1。S₁ → S₀ 遷移（光子放出）。

全系演算子は `build_single_site_operator` で テンソル積により $D\times D$（$D=3^4=81$）に拡張されるが、局所 Kraus 定式化ではこの局所 $3\times 3$ を直接使う。

---

## 3. GKSL 方程式と時間発展の枠組み

GKSL 方程式（$\hbar=1$）：

$$
\frac{d\rho}{dt} = -i[H,\rho] + \sum_\alpha\!\left(L_\alpha\,\rho\,L_\alpha^\dagger - \frac{1}{2}\{L_\alpha^\dagger L_\alpha,\,\rho\}\right) \equiv \mathcal{L}[\rho]
$$

形式解：$\rho(t+dt) = e^{\mathcal{L}\,dt}[\rho(t)]$

**数値実装では Strang splitting を使う（`qudit_gksl_simulator.py::_trotter_step`）：**

$$
\rho(t+dt) \approx e^{\mathcal{L}_H\,dt/2}\circ\!\left(\prod_{\alpha=1}^{n}\mathcal{E}_\alpha(dt/2)\right)\circ\!\left(\prod_{\alpha=n}^{1}\mathcal{E}_\alpha(dt/2)\right)\circ e^{\mathcal{L}_H\,dt/2}[\rho(t)]
$$

ここで $\mathcal{L}_H[\rho]=-i[H,\rho]$、$\mathcal{E}_\alpha(\tau)$ は Lindblad チャンネル $\alpha$ の $\tau$ 時間発展。

ハミルトニアン半ステップ：$e^{\mathcal{L}_H\,dt/2}[\rho] = U_H(dt/2)\,\rho\,U_H(dt/2)^\dagger$、$U_H(\tau)=e^{-iH\tau}$。

---

## 4. アルゴリズム A：Stinespring ダイレーション（近似、1 次）

### 4.1 一般構成（コード：`stinespring_utils.py:14-51`）

各 Lindblad 演算子 $L$（局所 $d_\mathrm{sys}\times d_\mathrm{sys}$）に対し、アンシラ（$d_\mathrm{anc}=3$）を追加して Stinespring ユニタリを構築。

**生成子 $G$（$(3d_\mathrm{sys})\times(3d_\mathrm{sys})$ Hermitian）：**

$$
G = \begin{pmatrix}
0_{d_\mathrm{sys}} & L^\dagger & 0 \\
L & 0_{d_\mathrm{sys}} & 0 \\
0 & 0 & 0_{d_\mathrm{sys}}
\end{pmatrix}
$$

ブロック $d_\mathrm{sys}\times d_\mathrm{sys}$、3 ブロック（アンシラ $|0\rangle,|1\rangle,|2\rangle$ に対応）。

コード（`stinespring_utils.py:34-38`）：
```python
G[:d_sys, d_sys:2*d_sys] = L.conj().T
G[d_sys:2*d_sys, :d_sys] = L
theta = np.sqrt(dt)  # ← dt は半ステップ dt/2 を渡すことに注意
U = expm(-1j * theta * G)
```

**Stinespring ユニタリ：**

$$
U_\alpha(\tau) = \exp\!\left(-i\sqrt{\tau}\;G_\alpha\right), \qquad \tau = dt/2
$$

$U_\alpha$ は $(3d_\mathrm{sys})\times(3d_\mathrm{sys})$ ユニタリ行列（`||U†U - I||_F < 1e-10` で検証済み）。

### 4.2 Kraus 演算子の抽出

アンシラを $|0\rangle_\mathrm{anc}$ に初期化した状態でチャンネルを適用し、アンシラについて部分トレースを取る（`stinespring_utils.py:54-81`）：

$$
\mathcal{E}_\alpha(\tau)[\rho_\mathrm{sys}] = \mathrm{Tr}_\mathrm{anc}\!\left[U_\alpha\,(|0\rangle\langle 0|_\mathrm{anc}\otimes\rho_\mathrm{sys})\,U_\alpha^\dagger\right] = \sum_{m=0}^{2} K_m\,\rho_\mathrm{sys}\,K_m^\dagger
$$

Kraus 演算子（$d_\mathrm{sys}\times d_\mathrm{sys}$）：

$$
K_m = \langle m|_\mathrm{anc}\;U_\alpha\;|0\rangle_\mathrm{anc}
= U_\alpha\!\left[m\cdot d_\mathrm{sys}:(m+1)d_\mathrm{sys},\;0:d_\mathrm{sys}\right]
$$

アンシラ–システム順の Kronecker 規約（コード：`np.kron(env0, rho)` は `env ⊗ sys` 順）。

**先頭次展開（$\theta=\sqrt{\tau}$）：**

$$
G^2 = \begin{pmatrix} L^\dagger L & 0 & 0 \\ 0 & LL^\dagger & 0 \\ 0 & 0 & 0 \end{pmatrix}
\implies U_\alpha = \begin{pmatrix}
\cos(\theta\sqrt{L^\dagger L}) & -iL^\dagger\,f(\theta,LL^\dagger) & 0 \\
-i\,g(\theta,L^\dagger L)\,L & \cos(\theta\sqrt{LL^\dagger}) & 0 \\
0 & 0 & I
\end{pmatrix}
$$

ここで $f,g$ は行列関数 $\mathrm{sinc}$ 型（詳細は §4.3 参照）。

$$
K_0 = \cos(\sqrt{\tau}\,\sqrt{L^\dagger L}),\quad K_1 = -i\,g(\sqrt{\tau},L^\dagger L)\cdot L,\quad K_2 = 0
$$

**1 次近似の確認：**

$$
\sum_m K_m\rho K_m^\dagger \approx \rho + \tau\!\left(L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho\}\right) + O(\tau^{3/2})
$$

つまり $\tau=dt/2$ として、これは $\mathcal{E}_\alpha(dt/2)$ の 1 次近似。

### 4.3 TTA チャンネル 1 の Stinespring（$d_\mathrm{sys}=9$、$U_\alpha\in\mathbb{C}^{27\times 27}$）

$L = \eta\;|20\rangle\langle 11|$（$\eta=\sqrt{\gamma_\mathrm{TTA}/2}$）。

$L^\dagger L$：

$$
L^\dagger L = \eta^2\;|11\rangle\langle 11|
= \frac{\gamma_\mathrm{TTA}}{2}\;\mathrm{diag}_9(0,0,0,0,1,0,0,0,0)
$$

（$9\times 9$ 対角行列。インデックス順：$|00\rangle,|01\rangle,|02\rangle,|10\rangle,|11\rangle,|12\rangle,|20\rangle,|21\rangle,|22\rangle$）

$LL^\dagger$：

$$
LL^\dagger = \eta^2\;|20\rangle\langle 20|
= \frac{\gamma_\mathrm{TTA}}{2}\;\mathrm{diag}_9(0,0,0,0,0,0,1,0,0)
$$

$\sqrt{L^\dagger L}=\eta\,|11\rangle\langle 11|$、$\sqrt{LL^\dagger}=\eta\,|20\rangle\langle 20|$（プロジェクタのべき乗なので自明）。

$$
K_0 = \cos(\sqrt{\tau}\cdot\sqrt{L^\dagger L})
= I_9 + (\cos(\sqrt{\tau}\cdot\eta)-1)\;|11\rangle\langle 11|
$$

$$
= \mathrm{diag}_9(1,1,1,1,\;\cos(\eta\sqrt{\tau}),\;1,1,1,1)
$$

ここで $\eta=\sqrt{\gamma_\mathrm{TTA}/2}$、$\tau=dt/2$。

$$
K_1 = -i\cdot\frac{\sin(\sqrt{\tau}\cdot\eta)}{\eta}\cdot L
= -i\sin(\eta\sqrt{\tau})\;\cdot\;|20\rangle\langle 11|
$$

すなわち $K_1$ は行 $|20\rangle$、列 $|11\rangle$ のみ $-i\sin(\eta\sqrt{dt/2})$、他 0 の $9\times 9$ 行列。

$K_2 = 0_{9\times 9}$。

**完全性確認：**

$$
K_0^\dagger K_0 + K_1^\dagger K_1 = \mathrm{diag}_9(1,\dots,1,\;\cos^2(\eta\sqrt{\tau}),\;1,\dots,1)
+ \sin^2(\eta\sqrt{\tau})\;|11\rangle\langle 11|
= I_9\quad\checkmark
$$

量子ゲート命令（`qudit_gksl_circuit_simulator.py:365-380`）：

```python
# TTA ペア (i,j)、アンシラ 1 個（d_anc=3）
circuit = QuantumCircuit(3, [d, d, d_anc], 0)
circuit.cu_multi([0, 1, 2], U_local)  # U_local は 27×27
```

`U_local` = $U_\alpha$ 行列そのもの（`scipy.linalg.expm` の出力）。**ネイティブ qutrit ゲートへの分解はコード上で行われない。**

### 4.4 蛍光チャンネルの Stinespring（$d_\mathrm{sys}=3$、$U_\alpha\in\mathbb{C}^{9\times 9}$）

$L = \sqrt{\Gamma_\mathrm{fl}}\;|0\rangle\langle 2|$。

$$
L^\dagger L = \Gamma_\mathrm{fl}\;|2\rangle\langle 2|
= \Gamma_\mathrm{fl}\;\mathrm{diag}_3(0,0,1)
$$

$$
LL^\dagger = \Gamma_\mathrm{fl}\;|0\rangle\langle 0|
= \Gamma_\mathrm{fl}\;\mathrm{diag}_3(1,0,0)
$$

$\sigma_L \equiv \sqrt{\Gamma_\mathrm{fl}}$（スカラー）として：

$$
K_0 = I_3 + (\cos(\sigma_L\sqrt{\tau})-1)\;|2\rangle\langle 2|
= \begin{pmatrix}1&0&0\\0&1&0\\0&0&\cos(\sigma_L\sqrt{\tau})\end{pmatrix}
$$

$$
K_1 = -i\sin(\sigma_L\sqrt{\tau})\;|0\rangle\langle 2|
= \begin{pmatrix}0&0&-i\sin(\sigma_L\sqrt{\tau})\\0&0&0\\0&0&0\end{pmatrix}
$$

$K_2 = 0_{3\times 3}$。

$\tau=dt/2$（半ステップ）。$\sigma_L=\sqrt{\Gamma_\mathrm{fl}}=\sqrt{0.01}\approx 0.1$ (eV/ℏ)^{1/2}。

**完全性：**

$$
K_0^\dagger K_0 + K_1^\dagger K_1 = \mathrm{diag}_3(1,1,\cos^2(\sigma_L\sqrt{\tau})) + \sin^2(\sigma_L\sqrt{\tau})\;\mathrm{diag}_3(1,0,0)\cdot |2\rangle\langle 2|
$$

整理すると $K_0^\dagger K_0 + K_1^\dagger K_1 = I_3\;\checkmark$（$(0,0)$成分：$\sin^2$、$(2,2)$成分：$\cos^2+0=\cos^2$、… 合計して $I_3$）。

量子ゲート命令（`qudit_gksl_circuit_simulator.py:336-355`）：

```python
# 蛍光サイト i、アンシラ 1 個（d_anc=3）
circuit = QuantumCircuit(2, [d, d_anc], 0)
circuit.cu_two([0, 1], U_local)  # U_local は 9×9
```

`U_local` の具体的 $9\times 9$ 成分（アンシラ $\otimes$ システム順、$d_\mathrm{anc}=3$、$c_s=\cos(\sigma_L\sqrt{\tau})$、$s_s=\sin(\sigma_L\sqrt{\tau})$）：

$G$ の非零ブロックは $G[0:3,3:6]=L^\dagger$（行 $|0_\mathrm{anc}\rangle$–システム）と $G[3:6,0:3]=L$（行 $|1_\mathrm{anc}\rangle$–システム）のみ。第三ブロック（$|2_\mathrm{anc}\rangle$）は恒等。

$U = e^{-i\sqrt{\tau}G}$ の構造：

$$
U_\alpha = \begin{pmatrix}
K_0 & -iK_1^\dagger/\text{（別形式）} & 0_3 \\
K_1 & \tilde{K}_0 & 0_3 \\
0_3 & 0_3 & I_3
\end{pmatrix}
$$

ただし $\tilde{K}_0 = \cos(\sigma_L\sqrt{\tau})\,|2\rangle\langle 2| + |0\rangle\langle 0| + |1\rangle\langle 1|$（$LL^\dagger$ 方向の余弦）。

明示すると（$c=\cos(\sigma_L\sqrt{\tau})$、$s=\sin(\sigma_L\sqrt{\tau})$、行列インデックス = アンシラ状態 $\times$ 3 + システム状態）：

$$
U_\alpha = \begin{pmatrix}
1&0&0 & 0&0&-is & 0&0&0 \\
0&1&0 & 0&0&0   & 0&0&0 \\
0&0&c & 0&0&0   & 0&0&0 \\
0&0&is& 1&0&0   & 0&0&0 \\
0&0&0 & 0&1&0   & 0&0&0 \\
0&0&0 & 0&0&c   & 0&0&0 \\
0&0&0 & 0&0&0   & 1&0&0 \\
0&0&0 & 0&0&0   & 0&1&0 \\
0&0&0 & 0&0&0   & 0&0&1
\end{pmatrix}
$$

行列のブロック位置：行 = アンシラ $|m\rangle$（上から $m=0,1,2$、各 3 行）、列 = アンシラ $|k\rangle$ 初期状態（左から $k=0,1,2$、各 3 列）。最初の 3 列（$k=0$ ブロック）が Kraus $K_m$ に対応する。

---

## 5. アルゴリズム B：完全局所チャンネル（正確、2 次）

### 5.1 一般構成（コード：`exact_local_channels.py:159-202`）

Lindblad 演算子 $L$（局所 $d_\mathrm{loc}\times d_\mathrm{loc}$、既に $\sqrt{\gamma}$ 込み）に対して局所散逸超演算子（列優先ベクトル化）：

$$
\mathcal{L}_D^\mathrm{local}\cdot\mathrm{vec}(\rho_\mathrm{loc})
= \mathrm{vec}\!\left(L\rho_\mathrm{loc}L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho_\mathrm{loc}\}\right)
$$

行列表現（$d_\mathrm{loc}^2\times d_\mathrm{loc}^2$）：

$$
M_\mathcal{L} = L^*\otimes L - \frac{1}{2}\left(I\otimes L^\dagger L + (L^\dagger L)^T\otimes I\right)
$$

コード（`exact_local_channels.py:195-202`）：
```python
LdL = L.conj().T @ L
return (
    np.kron(L.conj(), L)
    - 0.5 * np.kron(I, LdL)
    - 0.5 * np.kron(LdL.T, I)
)
```

**正確に指数化した半ステップチャンネル：**

$$
M_\alpha(\tau) = \exp\!\left(\mathcal{L}_D^\mathrm{local}\cdot\tau\right), \qquad \tau = dt/2
$$

$M_\alpha(\tau)$ は $d_\mathrm{loc}^2\times d_\mathrm{loc}^2$ 行列（`scipy.linalg.expm` で計算）。

これがチャンネル $\mathcal{E}_\alpha(\tau)$ の **完全な（近似なしの）** 超演算子表現。

### 5.2 Choi-Jamiolkowski 分解による Kraus 演算子抽出

コード（`dmsim_kraus_helpers.py:51-129`）：

1. **テンソル reshape**：$M[i+d\cdot j,\;k+d\cdot l]\to T[i,j,k,l]$（列優先 `order="F"`）

2. **Choi 行列**：$C[(ik),(jl)] = T[i,j,k,l]$

   ```python
   t_tensor = m_local.reshape(d, d, d, d, order="F")
   choi = t_tensor.transpose(0, 2, 1, 3).reshape(d*d, d*d)
   choi = 0.5 * (choi + choi.conj().T)  # 数値誤差対称化
   ```

3. **固有分解**：$C = \sum_\alpha\lambda_\alpha|\nu_\alpha\rangle\langle\nu_\alpha|$（`np.linalg.eigh`）

4. **Kraus 演算子**（$\lambda_\alpha > \epsilon_\mathrm{tol}=10^{-12}$ のみ）：
   $$
   K_\alpha = \sqrt{\lambda_\alpha}\;\mathrm{reshape}(\nu_\alpha,\;(d_\mathrm{loc},d_\mathrm{loc}))
   $$

5. **完全性検証**（コード上では暗黙）：$\sum_\alpha K_\alpha^\dagger K_\alpha = I_{d_\mathrm{loc}}$

**CPTP 条件**：Choi 行列の負の固有値が `CHOI_EIG_TOL=1e-12` を超えると `ValueError` を投げ、ヒューリスティック補正は一切行わない。

### 5.3 TTA チャンネル 1 の完全局所チャンネル（$d_\mathrm{loc}=9$）

$L = \eta\;|20\rangle\langle 11|$（$\eta=\sqrt{\gamma_\mathrm{TTA}/2}$）として、超演算子行列：

$$
\mathcal{L}_D^\mathrm{local} = \eta^2\!\left(|20\rangle\langle 11|^*\otimes|20\rangle\langle 11| - \frac{1}{2}I_9\otimes|11\rangle\langle 11| - \frac{1}{2}|11\rangle\langle 11|\otimes I_9\right)
$$

（$81\times 81$ 行列）

$M_\alpha(\tau)=\exp(\mathcal{L}_D^\mathrm{local}\cdot\tau)$ の非零成分（$\tau=dt/2$）は解析的に求まる。$L$ がランク 1 かつ 2 つのプロジェクタ $P_\mathrm{in}=|11\rangle\langle 11|$、$P_\mathrm{out}=|20\rangle\langle 20|$ で特徴付けられるため：

$\mathcal{L}_D^\mathrm{local}$ の固有値（非零）：
- $-\eta^2/2$ が `2` 重（$|11\rangle\langle 11|$ と $|20\rangle\langle 20|$ に関係）
- $-\eta^2$ が $1$ 重

指数 $M_\alpha(\tau)$ を作用させたチャンネルの Kraus 演算子（Choi 分解後）：

$$
K_0^{(B)} = I_9 + (e^{-\gamma_\mathrm{TTA}\tau/4}-1)\;|11\rangle\langle 11| + (e^{-\gamma_\mathrm{TTA}\tau/4}-1)\;|20\rangle\langle 20|
$$

（対角成分の変化のみ）

$$
K_1^{(B)} = \sqrt{1-e^{-\gamma_\mathrm{TTA}\tau/2}}\;|20\rangle\langle 11|
$$

> **注意**：この解析形は前文書 `シナリオ5d_TTAと蛍光のみ_完全展開.md` §6.1 で導出済み。ランク 1 Lindblad に限定した閉形式。数値実装では Choi 固有分解を使うため、解析形と数値形が一致することは固有分解の精度内で成立するが、コードは解析形を使わない。

### 5.4 蛍光チャンネルの完全局所チャンネル（$d_\mathrm{loc}=3$）

$L = \sigma_L\;|0\rangle\langle 2|$（$\sigma_L=\sqrt{\Gamma_\mathrm{fl}}$）。

超演算子（$9\times 9$）：

$$
\mathcal{L}_D^\mathrm{local} = \sigma_L^2\!\left(|0\rangle\langle 2|^*\otimes|0\rangle\langle 2| - \frac{1}{2}I_3\otimes|2\rangle\langle 2| - \frac{1}{2}|2\rangle\langle 2|\otimes I_3\right)
$$

固有値（非零）：$-\Gamma_\mathrm{fl}/2$（複数重複）、$-\Gamma_\mathrm{fl}$（1 重）。

解析的 Kraus（Choi 分解と一致、前文書 §7 参照）：

$$
K_0^{(B)} = \begin{pmatrix}1&0&0\\0&1&0\\0&0&e^{-\Gamma_\mathrm{fl}\tau/2}\end{pmatrix}, \qquad
K_1^{(B)} = \sqrt{1-e^{-\Gamma_\mathrm{fl}\tau}}\;\begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}
$$

$\tau=dt/2$。$K_0^{(B)}$ は $|2\rangle$ の振幅を減衰させ、$K_1^{(B)}$ は S₁→S₀ 遷移を記述。

### 5.5 量子ゲート命令（`execute_on_backend="dmsim"` 時）

コード（`qudit_gksl_simulator.py:180-204`）：

```python
for sites, M_half, kind in self._exact_channels_half:
    d_root = self.params.d if kind == "single" else self.params.d**2
    kraus = kraus_from_local_superoperator(M_half, d_root)
    kraus_list.append((sites, kraus, kind))
```

発行命令：
```python
circuit.kraus_channel(target_sites, kraus_ops)  # kraus_ops = [K_α]
```

`KrausChannel` 命令（`src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`）：
- `to_matrix()` は `NotImplementedError` → ユニタリ行列表現を持たない
- DMSim バックエンドのみが認識する（`dmsim.py`）
- DMSim は `np.tensordot` で $\rho\mapsto\sum_\alpha K_\alpha\rho_\mathrm{loc}K_\alpha^\dagger$ を直接計算
- **ネイティブ qutrit ゲートへの分解はコード上で行われない**

---

## 6. 2 アルゴリズムの比較表（TTA チャンネル 1 を例に）

| 項目 | A: Stinespring | B: 完全局所チャンネル |
|------|---------------|---------------------|
| チャンネル近似の精度 | O(τ)（1次、$\theta=\sqrt{\tau}$） | 完全（τについて機械精度） |
| 全体 Trotter 精度 | O(dt)（Lie 積誤差がO(dt²)になるが Stinespring が O(dt) で律速） | O(dt²)（Strang + palindromic） |
| Kraus の個数 | 2（$K_0,K_1$；$K_2=0$） | 2（ランク 1 Lindblad の場合；一般は最大 $d_\mathrm{loc}^2$） |
| $K_0$ の形 | $I_9+(\cos(\eta\sqrt{\tau})-1)|11\rangle\langle 11|$ | $I_9+(e^{-\gamma\tau/4}-1)(|11\rangle\langle 11|+|20\rangle\langle 20|)$ |
| $K_1$ の形 | $-i\sin(\eta\sqrt{\tau})|20\rangle\langle 11|$ | $\sqrt{1-e^{-\gamma\tau/2}}|20\rangle\langle 11|$ |
| 量子ゲート命令 | `cu_multi([i,j,anc], U_27x27)` | `kraus_channel([i,j], [K₀,K₁])` |
| アンシラ qutrit | 1 個必要（$|0\rangle$ に初期化、使い捨て） | 不要 |
| TNSim 対応 | ⛔ TNSim はアンシラ再利用不可 | ⛔ TNSim は `KrausChannel` 未対応 |

---

## 7. 1 Trotter ステップの完全な Kraus 表現（縮約 10 チャンネルモデル）

ステップの Kraus 合成（Strang + palindromic、`_trotter_step` の構造）：

$$
\mathcal{E}_\mathrm{step}(dt)
= \mathcal{U}_H(dt/2)\circ
\underbrace{\mathcal{E}_{\alpha_1}\cdots\mathcal{E}_{\alpha_{10}}}_{\text{forward}}(dt/2)\circ
\underbrace{\mathcal{E}_{\alpha_{10}}\cdots\mathcal{E}_{\alpha_1}}_{\text{reverse}}(dt/2)\circ
\mathcal{U}_H(dt/2)
$$

各 $\mathcal{U}_H(\tau)[\rho]=U_H(\tau)\,\rho\,U_H(\tau)^\dagger$ は単一の Kraus 演算子 $K=U_H$（ユニタリ）。

各 $\mathcal{E}_{\alpha_k}(dt/2)$ の Kraus は §4（Stinespring）または §5（完全局所）で与えられる。

**合成後の Kraus 演算子数**（最悪ケース）：

| アルゴリズム | 各チャンネルの Kraus 数 | 10 チャンネル forward × reverse | ハミルトニアン half × 2 | 合計（最悪） |
|------|---|---|---|---|
| A: Stinespring | 2 | $2^{10}\times 2^{10} = 2^{20}\approx 10^6$ | 各 1（ユニタリ） | $2^{20}$ |
| B: 完全局所 | 2（ランク 1 Lindblad） | 同上 | 各 1 | $2^{20}$ |

**注意**：コード上ではテンソル積を取らず局所チャンネルを逐次適用するため、$2^{20}$ 個の全系 Kraus を明示的に保持しない。

---

## 8. ハミルトニアン半ステップの Kraus（1 個）

$$
U_H(\tau) = e^{-iH_\mathrm{total}\tau}, \qquad \tau=dt/2
$$

コード（`qudit_gksl_simulator.py`）：
```python
self._U_H_half = expm(-1j * self.H_total * dt / 2)
rho = self._U_H_half @ rho @ self._U_H_half.conj().T
```

$H_\mathrm{total}=H_0+H_\mathrm{transfer}$（$81\times 81$ Hermitian）。

Kraus：$K_0=U_H(dt/2)$（単一のユニタリ行列、$81\times 81$）。

**量子ゲート命令（`execute_on_backend="dmsim"` 時）**：

```python
circuit.cu_multi(list(range(N)), U_H_half.astype(np.complex128))
```

$81\times 81$ 行列をそのまま `cu_multi` に渡す。ネイティブ分解なし。

---

## 9. 真実のまとめ

1. **「量子ゲート表現」はコード上では `cu_two`/`cu_multi`/`kraus_channel` への行列の直接代入であり、ネイティブ qutrit ゲート（`R`, `VirtRz`, `CSum`）への分解は実装されていない。**

2. **Stinespring アルゴリズム（A）** では Lindblad 演算子 $L$ から生成子 $G$ を構成し $U_\alpha=e^{-i\sqrt{\tau}G}$ を計算する。この $U_\alpha$ を `cu_two`（蛍光）または `cu_multi`（TTA）に代入し、アンシラ qutrit 1 個と `apply_stinespring_to_density_matrix` で部分トレースを取る。精度は 1 次（O(dt)）。

3. **完全局所チャンネルアルゴリズム（B）** では超演算子 $\mathcal{L}_D^\mathrm{local}$ を指数化して $M_\alpha=e^{\mathcal{L}_D\tau}$ を作り、Choi-Jamiolkowski 分解で Kraus 演算子 $\{K_\alpha\}$ を抽出し `kraus_channel` 命令で DMSim に渡す。アンシラ不要。精度は 2 次（O(dt²)）。

4. **TTA の Kraus** はランク 1 Lindblad なので解析形が閉じており、Stinespring では正弦・余弦、完全局所では指数関数で書ける（§4.3、§5.3）。

5. **蛍光の Kraus** も同様にランク 1 で、$|2\rangle$（S₁）成分のみ変化する対角行列 $K_0$ と遷移行列 $K_1$ の 2 演算子（§4.4、§5.4）。

6. **TNSim はどちらの経路でも利用できない**。`cu_multi` の行列次元が TNSim の想定を超え、かつ `KrausChannel` 命令を TNSim は認識しない（別文書 `TNSim_アンシラ制約_詳細.md` 参照）。

---

*根拠コード：*
- `tutorials/gksl_math_utils.py`（Lindblad 演算子構築）
- `tutorials/stinespring_utils.py`（Stinespring ユニタリ）
- `tutorials/exact_local_channels.py`（正確局所チャンネル）
- `tutorials/dmsim_kraus_helpers.py`（Choi-Kraus 分解）
- `tutorials/qudit_gksl_simulator.py`（Trotter ステップ、2 アルゴリズム）
- `tutorials/qudit_gksl_circuit_simulator.py`（回路命令構築）
- `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`（KrausChannel 命令）
- `src/mqt/qudits/simulation/backends/dmsim.py`（DMSim バックエンド）
