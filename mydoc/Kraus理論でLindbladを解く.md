# Kraus 演算子で Lindblad 方程式の時間発展を解く

> **方針：** 数式を丁寧に導き、「なぜ Kraus 演算子が使えるのか」を
> ステップごとに示す。事実のみを記述し、曖昧さを含まない。

---

## 0. 前提：開放量子系と密度行列

量子系が環境（浴）と相互作用している場合、系の状態は純粋状態ベクトル
$|\psi\rangle$ では記述できない。代わりに**密度行列**を使う：

$$
\rho \in \mathbb{C}^{d \times d}, \quad
\rho^\dagger = \rho, \quad
\rho \ge 0, \quad
\mathrm{tr}(\rho) = 1
$$

目標は、時刻 $t=0$ の密度行列 $\rho_0$ から時刻 $t = n\Delta t$ の密度行列
$\rho(t)$ を数値的に求めることである。

---

## 1. GKSL（Lindblad）方程式

マルコフ近似・Lindblad 形式のもとで、密度行列の時間発展は次の微分方程式に従う：

$$
\frac{d\rho}{dt}
= \mathcal{L}[\rho]
= -i[H,\rho]
+ \sum_{\alpha} \left(
    L_\alpha \rho L_\alpha^\dagger
    - \frac{1}{2} L_\alpha^\dagger L_\alpha \rho
    - \frac{1}{2} \rho L_\alpha^\dagger L_\alpha
  \right)
$$

- $H$：系のハミルトニアン（エルミート行列）
- $L_\alpha$：Lindblad 演算子（減衰・脱位相などを記述）、係数 $\sqrt{\gamma_\alpha}$ を含む
- $\mathcal{L}$：Liouvillian 超演算子（密度行列に作用する線形写像）

形式的な解は $\rho(t) = e^{\mathcal{L} t}[\rho_0]$ だが、
$d \times d$ 密度行列に対する超演算子 $\mathcal{L}$ は $d^2 \times d^2$ 行列に相当し、
$d = 3^4 = 81$（N=4 qutrits）では $81^2 = 6561$ 次元となるため、
直接 $\exp(\mathcal{L} t)$ を計算するのではなく、
小さなステップ $\Delta t$ に分けて逐次計算する。

---

## 2. 1 ステップの時間発展を Kraus 演算子で表す

### 2.1 CPTP 写像とは

十分短い時間 $\Delta t$ での時間発展 $\rho \mapsto \rho(\Delta t)$ は、
**完全正値（CP）かつトレース保存（TP）な写像**（= CPTP 写像）である。

Kraus 表現定理によれば、任意の CPTP 写像 $\mathcal{E}$ は次の形で書ける：

$$
\mathcal{E}(\rho) = \sum_{k} K_k \rho K_k^\dagger
$$

ここで $K_k$ は**Kraus 演算子**であり、トレース保存条件

$$
\sum_{k} K_k^\dagger K_k = I
$$

を満たす。

**この形式が使える理由**：  
Lindblad 方程式の右辺 $\mathcal{L}[\rho]$ は $\rho$ の線形関数であり、
時間 $\Delta t$ だけの発展 $e^{\mathcal{L}\Delta t}$ は CPTP 写像になることが
数学的に保証されている（Gorini-Kossakowski-Sudarshan-Lindblad の定理）。

---

### 2.2 ハミルトニアン部分の Kraus 表現

ハミルトニアン項 $\mathcal{L}_H[\rho] = -i[H, \rho]$ の時間発展は：

$$
e^{\mathcal{L}_H \Delta t}[\rho]
= U_H(\Delta t)\, \rho\, U_H(\Delta t)^\dagger, \quad
U_H(\Delta t) = e^{-iH\Delta t}
$$

これは Kraus 演算子が $K_0 = U_H$ の 1 本だけの場合に相当する。ユニタリーなので
$K_0^\dagger K_0 = I$ が自動的に成立する。

---

### 2.3 各 Lindblad チャンネルの Kraus 表現

散逸項 1 本（演算子 $L$）のみを考える。
$\Delta t \ll 1$ の極限で、局所超演算子

$$
\mathcal{L}_D[\rho]
= L\rho L^\dagger - \frac{1}{2}L^\dagger L\rho - \frac{1}{2}\rho L^\dagger L
$$

の短時間発展は**1 次近似**で次の Kraus 演算子 2 本で表せる：

$$
K_0 = I - \frac{\Delta t}{2} L^\dagger L, \qquad
K_1 = \sqrt{\Delta t}\, L
$$

**検証（$\Delta t$ の 1 次まで）：**

$$
K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger
= \left(I - \frac{\Delta t}{2}L^\dagger L\right)\rho\left(I - \frac{\Delta t}{2}L^\dagger L\right)
+ \Delta t\, L\rho L^\dagger
$$

$(\Delta t)^2$ 項を無視すると：

$$
= \rho - \frac{\Delta t}{2}L^\dagger L\rho - \frac{\Delta t}{2}\rho L^\dagger L
+ \Delta t\, L\rho L^\dagger
= \rho + \Delta t\,\mathcal{L}_D[\rho]
$$

これは $e^{\mathcal{L}_D \Delta t}[\rho]$ の $\Delta t$ 展開と一致する。

**完全性の確認：**

$$
K_0^\dagger K_0 + K_1^\dagger K_1
= \left(I - \frac{\Delta t}{2}L^\dagger L\right)^2 + \Delta t\, L^\dagger L
$$

$(\Delta t)^2$ 項を無視すると：

$$
= I - \Delta t\, L^\dagger L + \Delta t\, L^\dagger L = I \quad \checkmark
$$

---

### 2.4 厳密な Kraus 表現（局所超演算子の行列指数）

1 次近似では Trotter 誤差が $O(\Delta t)$ となる。
より高精度を求める場合は、局所散逸超演算子を行列として直接指数化する：

$$
\mathcal{E}_\alpha^{\text{exact}}(\rho) = \sum_k K_k^{(\alpha)} \rho \left(K_k^{(\alpha)}\right)^\dagger
$$

この Kraus 演算子 $\{K_k^{(\alpha)}\}$ は、局所超演算子行列

$$
\mathbf{L}_D^{(\alpha)} \in \mathbb{C}^{d^2 \times d^2}
$$

（ここで $d$ は局所ヒルベルト空間の次元）を $\Delta t$ だけ指数化し、

$$
\mathbf{S} = e^{\mathbf{L}_D^{(\alpha)} \Delta t}
$$

として得た $d^2 \times d^2$ 行列 $\mathbf{S}$（完全正値写像の行列表現）を
Choi-Jamiołkowski 同型

$$
\mathbf{C} = \frac{1}{d}\sum_{i,j} |i\rangle\langle j| \otimes \mathbf{S}(|i\rangle\langle j|)
$$

経由でChoi行列に変換し、固有値分解 $\mathbf{C} = \sum_k \lambda_k |v_k\rangle\langle v_k|$
（$\lambda_k \ge 0$）から

$$
K_k = \sqrt{\lambda_k \cdot d}\; \mathrm{reshape}(v_k, [d, d])
$$

で抽出する。これが本リポジトリの `dmsim_kraus_helpers.kraus_from_local_superoperator()` が行う処理である。

---

## 3. 対称 Trotter 分解（1 ステップの全体構造）

1 ステップ $\Delta t$ の全体発展を次の順序で適用する：

$$
\rho(t + \Delta t)
\approx
e^{\mathcal{L}_H \frac{\Delta t}{2}} \circ
\prod_{\alpha=1}^{n} \mathcal{E}_\alpha\!\left(\frac{\Delta t}{2}\right) \circ
\prod_{\alpha=n}^{1} \mathcal{E}_\alpha\!\left(\frac{\Delta t}{2}\right) \circ
e^{\mathcal{L}_H \frac{\Delta t}{2}}
\left[\rho(t)\right]
$$

図式的に：

```
ρ(t)
  ↓  U_H(Δt/2)  ρ U_H†         ← ハミルトニアン半ステップ
  ↓  E_1(Δt/2)                  ← Lindblad チャンネル 1（前向き）
  ↓  E_2(Δt/2)                  ← Lindblad チャンネル 2（前向き）
     ...
  ↓  E_n(Δt/2)                  ← Lindblad チャンネル n（前向き）
  ↓  E_n(Δt/2)                  ← Lindblad チャンネル n（逆向き＝回文）
     ...
  ↓  E_2(Δt/2)
  ↓  E_1(Δt/2)
  ↓  U_H(Δt/2)  ρ U_H†         ← ハミルトニアン半ステップ
ρ(t + Δt)
```

**対称（回文）配置の理由：**  
$[H, L_\alpha] \ne 0$ の場合、ハミルトニアン項と散逸項の交換子

$$
[\mathcal{L}_H, \mathcal{L}_D]
$$

が非零となり、単純な Lie-Trotter 積 $e^{\mathcal{L}_H \Delta t} e^{\mathcal{L}_D \Delta t}$
は $O(\Delta t)$ の誤差を持つ。対称配置（Strang 分割）は次数 2 の BCH 展開で
先頭誤差項を相殺し、$O(\Delta t^2)$ の精度を実現する。

**注意：** ただし各 Lindblad チャンネル自体が Stinespring 近似で
$O(\Delta t)$ 精度しか持たない場合（`algorithm="stinespring"`）は、
対称 Trotter 全体の精度も $O(\Delta t)$ に制限される。
$O(\Delta t^2)$ を得るには `algorithm="exact_local_channels"` が必要。

---

## 4. Stinespring 拡張との関係

Kraus 表現と等価な別の見方が Stinespring 拡張（dilation）である。

### 4.1 Stinespring の定理

CPTP 写像 $\mathcal{E}$ は、アンシラ（補助系）$\mathcal{H}_\mathrm{anc}$ を導入すると、
必ずユニタリー発展

$$
\mathcal{E}(\rho)
= \mathrm{tr}_\mathrm{anc}\!\left(
    V (\rho \otimes |0\rangle\langle 0|) V^\dagger
  \right)
$$

として実現できる。ここで $V : \mathcal{H}_\mathrm{sys} \otimes \mathcal{H}_\mathrm{anc} \to \mathcal{H}_\mathrm{sys} \otimes \mathcal{H}_\mathrm{anc}$ はユニタリー行列。

### 4.2 Kraus 演算子との対応

アンシラの次元を $d_\mathrm{anc}$ とし、$|k\rangle$（$k = 0, \ldots, d_\mathrm{anc}-1$）
をアンシラの正規直交基底とすると：

$$
K_k = \langle k |_\mathrm{anc} V | 0 \rangle_\mathrm{anc}
$$

具体的に、散逸項 1 本・$d_\mathrm{anc} = 2$ の場合：

$$
V = \begin{pmatrix} K_0 & -K_1^\dagger \\ K_1 & K_0^\dagger \end{pmatrix}
= \begin{pmatrix}
    I - \frac{\Delta t}{2}L^\dagger L & -\sqrt{\Delta t}\, L^\dagger \\
    \sqrt{\Delta t}\, L & I - \frac{\Delta t}{2}LL^\dagger
  \end{pmatrix}
+ O(\Delta t^{3/2})
$$

これが「Stinespring ユニタリー」である。

$\mathrm{tr}_\mathrm{anc}$ を取ると：

$$
\mathrm{tr}_\mathrm{anc}\!\left(V(\rho\otimes|0\rangle\langle 0|)V^\dagger\right)
= K_0\rho K_0^\dagger + K_1\rho K_1^\dagger
= \mathcal{E}_D(\rho)
$$

つまり「アンシラを $|0\rangle$ に初期化してユニタリーを掛け、アンシラをトレースアウト」
する操作と、Kraus 演算子を直接適用する操作は数学的に等価である。

### 4.3 コードにおける扱い

| 操作 | コード実装 | バックエンド使用 |
|---|---|---|
| Stinespring ユニタリーをアンシラ付きで直接適用し partial trace | `stinespring_utils.apply_stinespring_to_density_matrix()` | NumPy のみ |
| Choi-Kraus 抽出した Kraus 演算子を KrausChannel 命令で適用 | `dmsim_kraus_helpers.kraus_from_local_superoperator()` + `DMSim` | MQT-Qudits DMSim |
| 上記 Kraus を Qiskit の `Kraus` instruction として適用 | `QiskitQubitGKSLSimulator` | Qiskit Aer density_matrix |

---

## 5. N=4 での数値コスト

N=4 qutrits（$d=3$）の場合、密度行列のサイズは：

$$
\dim(\rho) = 3^4 \times 3^4 = 81 \times 81 = 6561 \text{ 要素}
$$

各 Kraus 適用 $K_k \rho K_k^\dagger$ は $81 \times 81$ 行列積 2 回
（計算量 $O(d^{2N})$）。Lindblad チャンネル数 $n_\alpha = 26$（N=4）、
対称 Trotter で各チャンネルを 2 回（前後）適用するため、
1 ステップあたり $26 \times 2 = 52$ 回の Kraus 適用が発生する。

実測値（`pytest TestDMSimBackendExecution`, t_max=100, n_steps=100）：

- NumPy のみ（`algorithm="stinespring"`）：約 0.88 秒
- DMSim backend（`algorithm="exact_local_channels"`）：約 2.68 秒
- 数値誤差：$\|\rho_\mathrm{NumPy} - \rho_\mathrm{DMSim}\|_F = 8.55 \times 10^{-14}$

---

## 6. 各アルゴリズムの精度まとめ

| アルゴリズム | 1 ステップ誤差 | Trotter 全体精度 | バックエンド実行可否（N=4） |
|---|---|---|---|
| `stinespring`（デフォルト） | $O(\Delta t)$ | $O(\Delta t)$ | NumPy のみ（DMSim は ValueError） |
| `exact_local_channels` + NumPy | $O(\Delta t^2)$（Strang） | $O(\Delta t^2)$（実測 rate≈2.0） | NumPy のみ |
| `exact_local_channels` + DMSim | $O(\Delta t^2)$ | $O(\Delta t^2)$ | **DMSim で実行済み** |
| Stinespring + TNSim/MISim | — | — | **不可（コードで拒否＋メモリ不足）** |

---

## 7. まとめ：Kraus 演算子を使った時間発展アルゴリズム（擬似コード）

```python
rho = rho_0          # 初期密度行列 (d^N × d^N)
dt = T / n_steps

for step in range(n_steps):
    # 1. ハミルトニアン半ステップ
    U = expm(-1j * H * dt / 2)
    rho = U @ rho @ U.conj().T

    # 2. 各 Lindblad チャンネル（前向き、dt/2）
    for alpha in range(n_channels):
        Ks = kraus_ops[alpha]          # Kraus 演算子のリスト
        rho = sum(K @ rho @ K.conj().T for K in Ks)

    # 3. 逆順（回文 = 対称 Trotter）
    for alpha in reversed(range(n_channels)):
        Ks = kraus_ops[alpha]
        rho = sum(K @ rho @ K.conj().T for K in Ks)

    # 4. ハミルトニアン半ステップ
    rho = U @ rho @ U.conj().T
```

各ステップ後に $\mathrm{tr}(\rho) \approx 1$、$\rho \ge 0$ が保たれることが
トレース保存・完全正値性の数値的証拠となる。

---

## 参考：本リポジトリでの実装箇所

| 概念 | ファイル | 関数 / クラス |
|---|---|---|
| Stinespring ユニタリー構築 | `tutorials/stinespring_utils.py` | `stinespring_unitary_from_lindblad()` |
| NumPy で partial trace 適用 | `tutorials/stinespring_utils.py` | `apply_stinespring_to_density_matrix()` |
| 厳密 Kraus 抽出（Choi 経由） | `tutorials/dmsim_kraus_helpers.py` | `kraus_from_local_superoperator()` |
| 全 Trotter ループ（NumPy）  | `tutorials/qudit_gksl_simulator.py` | `QuditGKSLSimulator.simulate()` |
| 全 Trotter ループ（DMSim）  | `tutorials/qudit_gksl_simulator.py` | `QuditGKSLSimulator._trotter_step_dmsim()` |
| MQT-Qudits 密度行列 backend | `src/mqt/qudits/simulation/backends/dmsim.py` | `DMSim` |
| CPTP 命令 | `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py` | `KrausChannel` |
