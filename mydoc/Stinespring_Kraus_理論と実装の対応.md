# Stinespring 処方・Kraus 演算子と図の理論との対応、および N=4 での計算可否

> **方針：** 図の内容とリポジトリの実コードを照合し、事実のみを記述する。
> 推測・希望的観測・不確かな記述は一切含めない。

---

## 1. 図が示す時間発展演算子の分類

図（タイトル「時間発展演算子整理」）には 4 種類の演算子が示されている。

| 番号 | 名称 | 次元（図） | 式 |
|------|------|-----------|-----|
| ① | $U_\mathrm{onsite}(\Delta t/2)$ | $3\times3$（1 qudit） | $\exp\!\left(-i\begin{pmatrix}0&0&0\\0&E_T&0\\0&0&E_S\end{pmatrix}\dfrac{\Delta t}{2}\right)$ |
| ② | $U_\mathrm{pair}^{(i,j)}(\Delta t/2)$ | $9\times9$（2 qudit） | $\exp\!\left(-iV(|01\rangle\langle10|+|10\rangle\langle01|)\dfrac{\Delta t}{2}\right)$ |
| ③ | $U_\mathrm{Stine}^\mathrm{single}$ | $6\times6$（2 qudit） | $\exp\!\left(-i\sqrt{\Delta t/2}\begin{pmatrix}0_{3\times3}&L_\mathrm{local}^\dagger\\L_\mathrm{local}&0_{3\times3}\end{pmatrix}\right)\in\mathbb{C}^{6\times6}$ |
| ④ | $U_\mathrm{Stine}^\mathrm{pair}$ | $18\times18$（multi qudit） | $\exp\!\left(-i\sqrt{\Delta t/2}\begin{pmatrix}0_{9\times9}&L_\mathrm{pair}^\dagger\\L_\mathrm{pair}&0_{9\times9}\end{pmatrix}\right)\in\mathbb{C}^{18\times18}$ |

図の右欄には「全て MQT-qudits のカスタムゲートを使ってシミュレーションする」と記されている。

---

## 2. GKSL 方程式と Trotter 分解（理論背景）

### 2.1 GKSL 方程式

系の密度行列 $\rho$ の時間発展は GKSL（Lindblad）方程式に従う：

$$
\frac{d\rho}{dt}
= -\frac{i}{\hbar}[H,\rho]
+ \sum_\alpha \!\left(L_\alpha\rho L_\alpha^\dagger - \tfrac{1}{2}L_\alpha^\dagger L_\alpha\rho - \tfrac{1}{2}\rho L_\alpha^\dagger L_\alpha\right)
$$

ここで $H = H_0 + H_\mathrm{transfer}$（オンサイト + 移動項）、$L_\alpha$ は Lindblad 演算子（既に $\sqrt{\gamma}$ を含む）。

### 2.2 対称 Trotter 分解（Strang 分割）

1 ステップ分の時間発展（時間幅 $\Delta t$）を以下の順序で近似する：

$$
e^{\mathcal{L}\Delta t}
\approx
e^{\mathcal{L}_H \Delta t/2}
\cdot \prod_{\alpha=1}^{n} \mathcal{E}_\alpha(\Delta t/2)
\cdot \prod_{\alpha=n}^{1} \mathcal{E}_\alpha(\Delta t/2)
\cdot e^{\mathcal{L}_H \Delta t/2}
$$

- $e^{\mathcal{L}_H \Delta t/2}$ はハミルトニアン部分の半ステップ：$\rho \mapsto U_H \rho U_H^\dagger$、$U_H = e^{-iH\Delta t/2}$
- $\mathcal{E}_\alpha(\Delta t/2)$ は各 Lindblad チャンネル $\alpha$ の半ステップ写像
- 逆順（回文）適用により Lie-Trotter 積の先頭次 commutator 誤差を相殺

### 2.3 ハミルトニアンの対可換性

$H_0$ と $H_\mathrm{transfer}$ は交換可能である。具体的に：

$$
[H_0, H_\mathrm{transfer}]
= \sum_i h_i^{(0)} \cdot \sum_{\langle j,k\rangle} V(|01\rangle\langle10|+\mathrm{h.c.})_{jk}
$$

$H_0 = \sum_i (E_T|1\rangle\langle1| + E_S|2\rangle\langle2|)_i$ は各サイトの数演算子の線形結合であり、移動項 $|01\rangle\langle10|$ は $T_1$ 状態数を保存するため、$[H_0,H_\mathrm{transfer}]=0$ が成立する。

したがって：

$$
e^{-i(H_0+H_\mathrm{transfer})\Delta t/2}
= e^{-iH_0\Delta t/2}\cdot e^{-iH_\mathrm{transfer}\Delta t/2}
$$

が**厳密に**成立する（Trotter 誤差なし）。

---

## 3. Stinespring 処方と Kraus 演算子の関係（完全な数式）

### 3.1 Stinespring 処方の定義

任意の CPTP（完全正値・トレース保存）写像 $\mathcal{E}$ は、補助系（アンシラ）を導入して以下のように表せる（Stinespring 定理）：

$$
\mathcal{E}(\rho)
= \mathrm{tr}_\mathrm{anc}\!\left[U_\mathrm{Stine}(\rho\otimes|0\rangle\langle0|_\mathrm{anc})U_\mathrm{Stine}^\dagger\right]
$$

ここで $U_\mathrm{Stine}$ は系＋アンシラ全体に作用するユニタリ演算子。

### 3.2 図の生成子 $G$

図③④の Stinespring ユニタリの生成子は：

$$
G = \begin{pmatrix}0 & L^\dagger \\ L & 0\end{pmatrix}
$$

このとき：

$$
U_\mathrm{Stine} = e^{-i\sqrt{\Delta t/2}\,G}
$$

### 3.3 Kraus 演算子の導出

$U_\mathrm{Stine}$ を $d_\mathrm{anc} \times d_\mathrm{sys}$ のブロック構造で書き、アンシラの $k$ 番目の基底 $|k\rangle_\mathrm{anc}$ で射影すると：

$$
K_k = {}_\mathrm{anc}\langle k|\,U_\mathrm{Stine}\,|0\rangle_\mathrm{anc}
\quad (k = 0, 1, \ldots, d_\mathrm{anc}-1)
$$

CPTP 写像は：

$$
\mathcal{E}(\rho) = \sum_{k=0}^{d_\mathrm{anc}-1} K_k\,\rho\, K_k^\dagger
$$

このとき完全性条件 $\sum_k K_k^\dagger K_k = I$ が成立する（ユニタリ性から自動的に保証）。

### 3.4 $d_\mathrm{anc}=2$ の場合の解析的表式

図は $d_\mathrm{anc}=2$（クビット・アンシラ）を想定している（次元が $6=2\times3$、$18=2\times9$）。

行列関数として厳密に：

$$
U_\mathrm{Stine} = e^{-i\theta G}
= \cos\!\left(\theta\sqrt{G^2}\right) - i\sin\!\left(\theta\sqrt{G^2}\right)/\sqrt{G^2}\cdot G
$$

ここで $\theta = \sqrt{\Delta t/2}$、$G^2 = \begin{pmatrix}L^\dagger L&0\\0&LL^\dagger\end{pmatrix}$ である。

Kraus 演算子：

$$
K_0 = {}_\mathrm{anc}\langle 0|\,U_\mathrm{Stine}\,|0\rangle_\mathrm{anc}
= \cos\!\left(\sqrt{\tfrac{\Delta t}{2}}\sqrt{L^\dagger L}\right)
$$

$$
K_1 = {}_\mathrm{anc}\langle 1|\,U_\mathrm{Stine}\,|0\rangle_\mathrm{anc}
= -i\,(L^\dagger L)^{-1/2}\,L\,\sin\!\left(\sqrt{\tfrac{\Delta t}{2}}\sqrt{L^\dagger L}\right)
$$

$L$ が零固有値を持つ場合は $(L^\dagger L)^{-1/2}L$ を擬逆行列で定義する。

### 3.5 小さい $\Delta t$ での展開（GKSL との対応）

$\Delta t \ll 1$ で Taylor 展開すると：

$$
K_0 = I - \tfrac{\Delta t}{4}L^\dagger L + O(\Delta t^2)
$$

$$
K_1 = -i\sqrt{\tfrac{\Delta t}{2}}\,L + O(\Delta t^{3/2})
$$

CPTP 写像：

$$
\mathcal{E}(\rho)
= K_0\rho K_0^\dagger + K_1\rho K_1^\dagger
= \rho + \frac{\Delta t}{2}\!\left(L\rho L^\dagger - \tfrac{1}{2}L^\dagger L\rho - \tfrac{1}{2}\rho L^\dagger L\right)
+ O(\Delta t^{3/2})
$$

これは GKSL 方程式の半ステップ dissipator $\exp(\mathcal{D}_\alpha \cdot \Delta t/2)$ の**1 次近似**（$O(\Delta t)$ 精度）に対応する。残差は $O(\Delta t^{3/2})$ から始まり、全体の Trotter ステップ収束は $O(\Delta t)$（1 次）となる。

---

## 4. 図の方法は実装されているか

### 4.1 実装されている部分

**Stinespring 処方（`algorithm="stinespring"`）**は実装されている。

ファイル：`tutorials/stinespring_utils.py`

```python
def stinespring_unitary_from_lindblad(L, dt, d_anc=2):
    G = np.zeros((d_anc * d_sys, d_anc * d_sys), dtype=np.complex128)
    G[:d_sys, d_sys:2*d_sys] = L.conj().T
    G[d_sys:2*d_sys, :d_sys] = L
    theta = np.sqrt(dt)         # ← sqrt(dt) を使用
    U = expm(-1j * theta * G)   # ← 図と同じ構造
    return U
```

呼び出し時は `dt/2` を渡しているため実効的に $\theta = \sqrt{\Delta t/2}$ となり、**図③④の数式と完全一致**する。

適用部分（`tutorials/stinespring_utils.py`）：

```python
def apply_stinespring_to_density_matrix(rho, U, d_anc=2):
    rho_ext = np.kron(env0, rho)           # アンシラ |0><0| ⊗ ρ
    rho_prime = U @ rho_ext @ U.conj().T   # U(|0><0|⊗ρ)U†
    for k in range(d_anc):
        rho_out += rho_prime[k*d_sys:(k+1)*d_sys, k*d_sys:(k+1)*d_sys]  # 部分トレース
    return rho_out
```

これは $\sum_k K_k\rho K_k^\dagger$ の計算（$K_k = {}_\mathrm{anc}\langle k|U|0\rangle_\mathrm{anc}$）と**数学的に等価**である。

**対称 Trotter ステップ（回文順序）**も実装されている（`tutorials/qudit_gksl_simulator.py:206-250`）。

### 4.2 実装されているが図と次元が異なる部分

**アンシラ次元が異なる。**

| | 図の想定 | 実コードの実際 |
|--|---------|--------------|
| アンシラ次元 $d_\mathrm{anc}$ | 2（クビット） | **3**（クトリット）|
| $U_\mathrm{Stine}^\mathrm{single}$ の次元 | $6\times6$ | **$9\times9$** |
| $U_\mathrm{Stine}^\mathrm{pair}$ の次元 | $18\times18$ | **$27\times27$** |

コード（`tutorials/qudit_gksl_simulator.py:133`）：

```python
self.d_anc = params.d  # = 3（クトリット・アンシラ）
```

`d_anc=3` を使う理由：MQT-Qudits はクトリットネイティブ機であるため、アンシラもクトリットを使う設計になっている。`d_anc=3` の場合、$K_2 = {}_\mathrm{anc}\langle2|U|0\rangle_\mathrm{anc} \approx 0$（docstring に明記）であり、$K_0,K_1$ は $d_\mathrm{anc}=2$ と同一になる。

### 4.3 実装されていない部分（図との差異）

**ハミルトニアンの分解方法が異なる。**

- 図①：$U_\mathrm{onsite}(\Delta t/2) = \exp(-ih_\mathrm{loc}\Delta t/2)$ を**各サイト独立**の $3\times3$ カスタム 1-qudit ゲートとして適用
- 図②：$U_\mathrm{pair}^{(i,j)}(\Delta t/2)$ を**隣接ペア**の $9\times9$ カスタム 2-qudit ゲートとして適用

**実コード**（`tutorials/qudit_gksl_simulator.py:179`）は：

```python
self._U_H_half = expm(-1j * self.H_total * dt / 2)
# H_total = H_0 + H_transfer（全サイトを一括した dim×dim 行列）
```

として計算し、`cu_multi`（N-qudit 全体ユニタリ）として適用する。

$[H_0, H_\mathrm{transfer}]=0$ であるため図の分解は**理論的に等価**であるが、実コードは $U_\mathrm{onsite}$ と $U_\mathrm{pair}^{(i,j)}$ を分離せず一括 `expm` で計算している。

---

## 5. N=4 での計算可否

### 5.1 `algorithm="stinespring"` モード

**可能**。

- 系の次元：$d^N = 3^4 = 81$
- $U_H^\mathrm{half}$：$81\times81$ 行列（`scipy.linalg.expm` で計算）
- 各チャンネルの $U_\mathrm{Stine}$：単チャンネル $27\times27$、ペアチャンネル $27\times27$（$d_\mathrm{anc}=3$）

Lindblad チャンネル数（N=4 の場合）：
- TTA ペア：$(N-1)\times2 = 3\times2 = 6$ チャンネル
- 蛍光（Fluorescence）：$N=4$ チャンネル
- リン光（Phosphorescence）：$N=4$ チャンネル
- 内部変換（IC）：$N=4$ チャンネル
- ISC S→T：$N=4$ チャンネル
- ISC T→S：$N=4$ チャンネル
- **合計：26 チャンネル**（`gksl_physical_parameters.py` の `GKSLPhysicalParameters` デフォルト）

収束精度：**O($\Delta t$)（1 次）**。exact_local_channels モードより精度が低い。

### 5.2 `algorithm="exact_local_channels"` モード（推奨）

**可能、かつ高精度。**

各 Lindblad チャンネルを Stinespring 近似せず、局所超演算子を厳密に指数関数化する：

$$
\mathcal{E}_\alpha(\Delta t/2)(\rho) = e^{\mathcal{L}_{D_\alpha}^\mathrm{local}\,\Delta t/2}(\rho)
$$

局所超演算子：
- 単サイト：$\mathcal{L}_{D_\alpha}^\mathrm{local}$ は $9\times9$（$d^2\times d^2$、$d=3$）
- ペア：$\mathcal{L}_{D_\alpha}^\mathrm{local}$ は $81\times81$（$d^4\times d^4$）

全体の系密度行列（$81\times81$）には `numpy.einsum` による部分縮約で適用し、$81^2\times81^2$（= 6561×6561）の全超演算子を**構築しない**。

`exact_local_channels.py` の docstring（実測値）：

```
N=4 (dim=81), t_max=100: rate → 2.004 (n_steps=200 で漸近 2 次収束)
```

収束精度：**O($\Delta t^2$)（2 次）**。

### 5.3 `execute_on_backend="dmsim"` モード

**可能**（`algorithm="exact_local_channels"` と組み合わせ必須）。

このモードでは、各チャンネルの局所超演算子を**Choi-Jamiołkowski 同型**を通じて Kraus 演算子に変換し、MQT-Qudits の `KrausChannel` 命令としてバックエンドに渡す。

Choi 行列からの Kraus 変換（`tutorials/dmsim_kraus_helpers.py`）：

$$
C_{ij,kl} = M_{i+d\cdot j,\; k+d\cdot l} \quad\text{(column-major)}
$$

$$
C = \sum_\alpha \lambda_\alpha |v_\alpha\rangle\langle v_\alpha|
\quad\Rightarrow\quad
K_\alpha = \sqrt{\lambda_\alpha}\cdot \mathrm{reshape}(|v_\alpha\rangle,\;(d,d))
$$

$$
\sum_\alpha K_\alpha^\dagger K_\alpha = I \quad\text{（数値検証済み）}
$$

この Kraus 表現は Stinespring 由来ではなく、厳密な局所超演算子の Choi 分解から得られる。

---

## 6. 図の方法との対応まとめ（真実ベース）

| 図の要素 | 対応コード | 整合性 |
|----------|-----------|--------|
| ①$U_\mathrm{onsite}$（$3\times3$） | `expm(-1j*H_total*dt/2)` の一部として**暗黙に含まれる**が、独立した $3\times3$ ゲートとして**分離されない** | 理論等価（$[H_0,H_\mathrm{tr}]=0$），実装方法が異なる |
| ②$U_\mathrm{pair}^{(i,j)}$（$9\times9$） | 同上（一括 `cu_multi`） | 同上 |
| ③$U_\mathrm{Stine}^\mathrm{single}$（$6\times6$） | `stinespring_unitary_from_lindblad(L,dt/2,d_anc=3)`（**$9\times9$**） | 数式は一致、$d_\mathrm{anc}$ が異なる |
| ④$U_\mathrm{Stine}^\mathrm{pair}$（$18\times18$） | 同上（**$27\times27$**） | 数式は一致、$d_\mathrm{anc}$ が異なる |
| Kraus = $K_k=\langle k|U_\mathrm{Stine}|0\rangle$ | `apply_stinespring_to_density_matrix` が部分トレースで等価計算 | `algorithm="stinespring"` 時のみ |
| Kraus（dmsim バックエンド） | `kraus_from_local_superoperator`（Choi 分解） | Stinespring 由来ではなく、別ルート |

---

## 7. 結論（真実のみ）

1. **図の理論（Stinespring 処方 + 対称 Trotter）は原理的にこのリポジトリで考えることができる。**
   GKSL 方程式を Stinespring ユニタリで近似し、対称 Trotter（回文順序）で時間発展させる枠組みは `tutorials/stinespring_utils.py` および `tutorials/qudit_gksl_simulator.py` に実装されている。

2. **Kraus 演算子は Stinespring ユニタリのアンシラ射影として定義される。**
   具体的には $K_k = {}_\mathrm{anc}\langle k|U_\mathrm{Stine}|0\rangle_\mathrm{anc}$ であり、`apply_stinespring_to_density_matrix` 関数の部分トレースがこれと等価な計算を行っている。

3. **図と実コードには次元の相違がある（$d_\mathrm{anc}=2$ vs $d_\mathrm{anc}=3$）。**
   物理的な Kraus 演算子（$K_0, K_1$）は両者で同一であるが、アンシラ次元が異なるため行列サイズが異なる。$d_\mathrm{anc}=3$ では $K_2\approx0$ が追加されるだけで精度は同等。

4. **ハミルトニアン部分は図の分解方法（①②の分離）と実装方法（一括 `cu_multi`）が異なる。** ただし $[H_0,H_\mathrm{transfer}]=0$ のため両者は数学的に完全等価であり、Trotter 誤差は生じない。

5. **N=4 での計算は可能である。** `algorithm="stinespring"` は O($\Delta t$)（1 次）、`algorithm="exact_local_channels"` は O($\Delta t^2$)（2 次）で動作することが実測済み（`exact_local_channels.py` docstring に実測値記載）。

6. **dmsim バックエンド使用時は Stinespring ではなく Choi-Jamiołkowski 由来の Kraus が使われる。**
   この場合 Kraus は Stinespring の射影 $K_k=\langle k|U_\mathrm{Stine}|0\rangle$ ではなく、厳密局所超演算子の固有値分解から得られる。

---

## 付録 A：コードの参照箇所

| 内容 | ファイル | 行 |
|------|---------|-----|
| Stinespring ユニタリ構築 | `tutorials/stinespring_utils.py` | 14–51 |
| 部分トレース（= Kraus 適用） | `tutorials/stinespring_utils.py` | 54–83 |
| GKSL 超演算子構築 | `tutorials/stinespring_utils.py` | 86–124 |
| $d_\mathrm{anc} = d = 3$ の設定 | `tutorials/qudit_gksl_simulator.py` | 133 |
| 対称 Trotter ステップ | `tutorials/qudit_gksl_simulator.py` | 206–250 |
| dmsim バックエンド Trotter ステップ | `tutorials/qudit_gksl_simulator.py` | 256–312 |
| Choi → Kraus 変換 | `tutorials/dmsim_kraus_helpers.py` | 51–129 |
| 局所超演算子の厳密指数化 | `tutorials/exact_local_channels.py` | 159–202 |
| N=4 収束実測値 | `tutorials/exact_local_channels.py` | 54–55 |
| Lindblad 演算子の構築 | `tutorials/gksl_math_utils.py` | 66–128 |
| 物理パラメータのデフォルト（N=4） | `tutorials/gksl_physical_parameters.py` | 32 |
