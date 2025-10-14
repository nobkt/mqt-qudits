# 分子系における励起状態の量子ダイナミクス理論

## 1. はじめに

本文書では、基底１重項状態、励起３重項状態、励起１重項状態からなる分子系の量子ダイナミクスについて、省略無しに数式を用いて詳細に記述する。

## 2. 分子の電子状態

### 2.1 状態の定義

各分子は以下の3つの電子状態のいずれかをとる：

- **基底１重項状態** $|S_0\rangle$：スピン多重度が1の基底状態
- **励起３重項状態** $|T_1\rangle$：スピン多重度が3の励起状態
- **励起１重項状態** $|S_1\rangle$：スピン多重度が1の励起状態

### 2.2 エネルギー準位

各状態のエネルギーを以下のように定義する：

$$
E_{S_0} = 0
$$

$$
E_{T_1} = E_T
$$

$$
E_{S_1} = E_S
$$

ここで、一般的に以下のエネルギー順位関係が成り立つ：

$$
E_{S_0} < E_{T_1} < E_{S_1}
$$

すなわち、

$$
0 < E_T < E_S
$$

## 3. 量子状態の記述

### 3.1 単一分子のハミルトニアン

単一分子 $i$ のハミルトニアンは以下のように表される：

$$
\hat{H}_i = E_{S_0}|S_0\rangle_i\langle S_0| + E_{T_1}|T_1\rangle_i\langle T_1| + E_{S_1}|S_1\rangle_i\langle S_1|
$$

状態演算子を導入すると：

$$
\hat{n}_{S_0}^{(i)} = |S_0\rangle_i\langle S_0|
$$

$$
\hat{n}_{T_1}^{(i)} = |T_1\rangle_i\langle T_1|
$$

$$
\hat{n}_{S_1}^{(i)} = |S_1\rangle_i\langle S_1|
$$

これらは完全性関係を満たす：

$$
\hat{n}_{S_0}^{(i)} + \hat{n}_{T_1}^{(i)} + \hat{n}_{S_1}^{(i)} = \hat{1}
$$

### 3.2 N分子系の状態空間

N個の分子からなる系の状態は、テンソル積空間に属する：

$$
|\Psi\rangle \in \mathcal{H} = \bigotimes_{i=1}^{N} \mathcal{H}_i
$$

ここで、$\mathcal{H}_i = \text{span}\{|S_0\rangle_i, |T_1\rangle_i, |S_1\rangle_i\}$ は分子 $i$ の3次元ヒルベルト空間である。

### 3.3 多体基底状態

多体系の基底状態は：

$$
|\Phi_0\rangle = \bigotimes_{i=1}^{N} |S_0\rangle_i
$$

## 4. 励起エネルギー移動過程

### 4.1 三重項励起エネルギー移動

励起３重項状態の分子 $i$ が隣接した基底状態の分子 $j$ に励起エネルギーを移動する過程を考える。

#### 4.1.1 移動演算子

エネルギー移動演算子を以下のように定義する：

$$
\hat{T}_{ij}^{\text{transfer}} = |S_0\rangle_i\langle T_1| \otimes |T_1\rangle_j\langle S_0|
$$

これは以下の状態遷移を表す：

$$
|T_1\rangle_i \otimes |S_0\rangle_j \xrightarrow{\hat{T}_{ij}^{\text{transfer}}} |S_0\rangle_i \otimes |T_1\rangle_j
$$

#### 4.1.2 相互作用ハミルトニアン

隣接分子間の相互作用ハミルトニアンは：

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( \hat{T}_{ij}^{\text{transfer}} + \hat{T}_{ij}^{\text{transfer}\dagger} \right)
$$

ここで、$\langle i,j \rangle$ は隣接分子対を表し、$V_{ij}$ は移動積分（transfer integral）である。

エルミート演算子として明示的に書くと：

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( |S_0\rangle_i\langle T_1| \otimes |T_1\rangle_j\langle S_0| + |T_1\rangle_i\langle S_0| \otimes |S_0\rangle_j\langle T_1| \right)
$$

#### 4.1.3 フェルミの黄金律による移動速度

時間依存摂動論（フェルミの黄金律）により、遷移速度は：

$$
\Gamma_{ij}^{\text{transfer}} = \frac{2\pi}{\hbar} |V_{ij}|^2 \delta(E_{T_1} - E_{T_1})
$$

エネルギー保存により、三重項間のエネルギー移動速度は：

$$
\Gamma_{ij}^{\text{transfer}} = \frac{2\pi}{\hbar} |V_{ij}|^2 \rho(E_T)
$$

ここで、$\rho(E_T)$ は状態密度である。

## 5. 三重項-三重項消滅過程

### 5.1 TTA過程の定義

隣接した2つの励起３重項状態の分子 $i$ と $j$ が相互作用し、片方が励起１重項状態 $|S_1\rangle$ に励起され、もう片方が基底状態 $|S_0\rangle$ に戻る過程を三重項-三重項消滅（Triplet-Triplet Annihilation, TTA）という。

### 5.2 TTA演算子

TTA過程を表す演算子は：

$$
\hat{T}_{ij}^{\text{TTA}} = |S_1\rangle_i\langle T_1| \otimes |S_0\rangle_j\langle T_1| + |S_0\rangle_i\langle T_1| \otimes |S_1\rangle_j\langle T_1|
$$

これは以下の2つの遷移過程を含む：

1. 分子 $i$ が $|T_1\rangle \to |S_1\rangle$、分子 $j$ が $|T_1\rangle \to |S_0\rangle$
2. 分子 $i$ が $|T_1\rangle \to |S_0\rangle$、分子 $j$ が $|T_1\rangle \to |S_1\rangle$

### 5.3 TTA相互作用ハミルトニアン

TTA相互作用ハミルトニアンは：

$$
\hat{H}_{\text{TTA}} = \sum_{\langle i,j \rangle} J_{ij} \left( \hat{T}_{ij}^{\text{TTA}} + \hat{T}_{ij}^{\text{TTA}\dagger} \right)
$$

ここで、$J_{ij}$ はTTA相互作用定数である。

### 5.4 エネルギー保存則

TTA過程では以下のエネルギー保存が要求される：

$$
2E_{T_1} = E_{S_1} + E_{S_0}
$$

つまり：

$$
2E_T = E_S
$$

実際の系では、フォノンの関与により：

$$
2E_T \approx E_S + \Delta E_{\text{phonon}}
$$

ここで、$\Delta E_{\text{phonon}}$ は格子振動への散逸エネルギーである。

### 5.5 TTA速度定数

TTA過程の速度定数は、フェルミの黄金律により：

$$
k_{\text{TTA}} = \frac{2\pi}{\hbar} |J_{ij}|^2 g(E)
$$

ここで、$g(E)$ は状態密度とフランク-コンドン因子を含む有効状態密度である。

### 5.6 スピン選択則

３重項状態 $|T_1\rangle$ は3つのスピン副準位を持つ：

$$
|T_1, m_s\rangle, \quad m_s = -1, 0, +1
$$

TTA過程では、スピン角運動量保存則により：

$$
m_s^{(i)} + m_s^{(j)} = M_S
$$

ここで、$M_S$ は生成される一重項状態の全スピン量子数（$M_S = 0$）である。

スピン統計により、TTA効率は：

$$
\eta_{\text{TTA}} = \frac{1}{9} \quad (\text{ランダムスピン配向の場合})
$$

## 6. 蛍光発光過程

### 6.1 自然放出

励起１重項状態 $|S_1\rangle$ は自然放出により基底状態 $|S_0\rangle$ に遷移し、光子を放出する。

#### 6.1.1 自然放出演算子

自然放出過程は以下の演算子で記述される：

$$
\hat{A}_i = |S_0\rangle_i\langle S_1|
$$

#### 6.1.2 光子場との相互作用

光子場を含む完全なハミルトニアンは：

$$
\hat{H}_{\text{rad}} = \sum_i \sum_{\mathbf{k},\lambda} g_{\mathbf{k}\lambda} \left( \hat{A}_i \hat{a}_{\mathbf{k}\lambda}^\dagger + \hat{A}_i^\dagger \hat{a}_{\mathbf{k}\lambda} \right)
$$

ここで：
- $\hat{a}_{\mathbf{k}\lambda}^\dagger$ は波数 $\mathbf{k}$、偏光 $\lambda$ の光子生成演算子
- $g_{\mathbf{k}\lambda}$ は光-物質相互作用定数

$$
g_{\mathbf{k}\lambda} = \sqrt{\frac{\hbar \omega_k}{2\epsilon_0 V}} \mathbf{d}_{S_1 \to S_0} \cdot \boldsymbol{\epsilon}_{\mathbf{k}\lambda}
$$

ここで：
- $\omega_k = |\mathbf{k}|c$ は光子角周波数
- $V$ は量子化体積
- $\mathbf{d}_{S_1 \to S_0}$ は遷移双極子モーメント
- $\boldsymbol{\epsilon}_{\mathbf{k}\lambda}$ は偏光ベクトル

### 6.2 蛍光寿命と放出速度

#### 6.2.1 自然放出速度

フェルミの黄金律により、自然放出速度（Einstein A係数）は：

$$
A_{S_1 \to S_0} = \frac{\omega^3}{3\pi\epsilon_0\hbar c^3} |\mathbf{d}_{S_1 \to S_0}|^2
$$

ここで、$\omega = (E_S - E_{S_0})/\hbar = E_S/\hbar$ は遷移角周波数である。

#### 6.2.2 蛍光寿命

蛍光寿命 $\tau_{\text{fl}}$ は：

$$
\tau_{\text{fl}} = \frac{1}{A_{S_1 \to S_0}}
$$

典型的な一重項-一重項遷移では：

$$
\tau_{\text{fl}} \sim 1-10 \text{ ns}
$$

### 6.3 蛍光量子収率

蛍光量子収率 $\Phi_{\text{fl}}$ は、放射過程と無放射過程の競合により決まる：

$$
\Phi_{\text{fl}} = \frac{k_{\text{rad}}}{k_{\text{rad}} + k_{\text{nr}}}
$$

ここで：
- $k_{\text{rad}} = A_{S_1 \to S_0}$ は放射速度定数
- $k_{\text{nr}}$ は無放射失活速度定数

蛍光量子収率は以下のようにも表される：

$$
\Phi_{\text{fl}} = \frac{\tau_{\text{obs}}}{\tau_{\text{fl}}}
$$

ここで、$\tau_{\text{obs}}$ は観測される蛍光寿命である。

## 7. 完全な系のハミルトニアン

### 7.1 全ハミルトニアン

系全体のハミルトニアンは以下の和として表される：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}} + \hat{H}_{\text{rad}}
$$

ここで：

$$
\hat{H}_0 = \sum_{i=1}^{N} \left( E_T \hat{n}_{T_1}^{(i)} + E_S \hat{n}_{S_1}^{(i)} \right)
$$

### 7.2 明示的表現

完全なハミルトニアンを明示的に書くと：

$$
\begin{align}
\hat{H}_{\text{total}} = &\sum_{i=1}^{N} \left( E_T |T_1\rangle_i\langle T_1| + E_S |S_1\rangle_i\langle S_1| \right) \\
&+ \sum_{\langle i,j \rangle} V_{ij} \left( |S_0\rangle_i\langle T_1| \otimes |T_1\rangle_j\langle S_0| + \text{h.c.} \right) \\
&+ \sum_{\langle i,j \rangle} J_{ij} \left( |S_1\rangle_i\langle T_1| \otimes |S_0\rangle_j\langle T_1| + |S_0\rangle_i\langle T_1| \otimes |S_1\rangle_j\langle T_1| + \text{h.c.} \right) \\
&+ \sum_i \sum_{\mathbf{k},\lambda} g_{\mathbf{k}\lambda} \left( |S_0\rangle_i\langle S_1| \hat{a}_{\mathbf{k}\lambda}^\dagger + \text{h.c.} \right)
\end{align}
$$

ここで、h.c. はエルミート共役（Hermitian conjugate）を表す。

## 8. マスター方程式による動力学記述

### 8.1 密度行列の時間発展

系の状態を密度行列 $\hat{\rho}(t)$ で記述する。完全な時間発展は以下のリウヴィル-フォン・ノイマン方程式に従う：

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}_{\text{total}}, \hat{\rho}] + \hat{\mathcal{L}}_{\text{diss}}[\hat{\rho}]
$$

### 8.2 散逸項

散逸項 $\hat{\mathcal{L}}_{\text{diss}}$ は、リンドブラッド形式で表される：

$$
\hat{\mathcal{L}}_{\text{diss}}[\hat{\rho}] = \sum_\alpha \left( \hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger - \frac{1}{2}\{\hat{L}_\alpha^\dagger \hat{L}_\alpha, \hat{\rho}\} \right)
$$

ここで、$\hat{L}_\alpha$ はリンドブラッド演算子である。

### 8.3 蛍光放出のリンドブラッド演算子

蛍光放出に対するリンドブラッド演算子は：

$$
\hat{L}_i^{\text{fl}} = \sqrt{\Gamma_{\text{fl}}} |S_0\rangle_i\langle S_1|
$$

ここで、$\Gamma_{\text{fl}} = A_{S_1 \to S_0}$ は蛍光放出速度である。

### 8.4 個体群動力学

分子 $i$ の各状態の占有確率を：

$$
P_{S_0}^{(i)} = \langle S_0|_i \hat{\rho} |S_0\rangle_i, \quad P_{T_1}^{(i)} = \langle T_1|_i \hat{\rho} |T_1\rangle_i, \quad P_{S_1}^{(i)} = \langle S_1|_i \hat{\rho} |S_1\rangle_i
$$

と定義する。

## 9. 速度方程式

### 9.1 簡略化された速度方程式

平均場近似と非コヒーレント極限では、個体群の時間発展は速度方程式で記述される：

$$
\frac{dN_{S_0}}{dt} = -k_{\text{exc}} N_{S_0} + k_{\text{transfer}} N_{T_1} + \Gamma_{\text{fl}} N_{S_1} + \frac{1}{2}k_{\text{TTA}} N_{T_1}^2
$$

$$
\frac{dN_{T_1}}{dt} = k_{\text{exc}} N_{S_0} - k_{\text{transfer}} N_{T_1} - k_{\text{TTA}} N_{T_1}^2 - k_{\text{T-decay}} N_{T_1}
$$

$$
\frac{dN_{S_1}}{dt} = \frac{1}{2}k_{\text{TTA}} N_{T_1}^2 - \Gamma_{\text{fl}} N_{S_1} - k_{\text{S-decay}} N_{S_1}
$$

ここで：
- $N_{S_0}, N_{T_1}, N_{S_1}$ は各状態の分子数密度
- $k_{\text{exc}}$ は励起速度定数
- $k_{\text{transfer}}$ は三重項エネルギー移動速度定数
- $k_{\text{TTA}}$ はTTA速度定数
- $k_{\text{T-decay}}$ は三重項状態の無放射失活速度定数
- $k_{\text{S-decay}}$ は一重項状態の無放射失活速度定数
- $\Gamma_{\text{fl}}$ は蛍光放出速度

### 9.2 保存則

全分子数は保存される：

$$
N_{\text{total}} = N_{S_0} + N_{T_1} + N_{S_1} = \text{const.}
$$

### 9.3 TTA項の因子

TTA過程では2つの三重項分子が消滅するため、三重項個体数の減少率は $-k_{\text{TTA}} N_{T_1}^2$ となる。一方、生成される一重項分子数は2つの三重項から1つなので、$+\frac{1}{2}k_{\text{TTA}} N_{T_1}^2$ となる。

## 10. 空間依存性を考慮した拡散-反応方程式

### 10.1 拡散項を含む方程式

空間的に非一様な系では、拡散を考慮した反応-拡散方程式を用いる：

$$
\frac{\partial n_{S_0}(\mathbf{r},t)}{\partial t} = D_{S_0} \nabla^2 n_{S_0} + \text{(反応項)}
$$

$$
\frac{\partial n_{T_1}(\mathbf{r},t)}{\partial t} = D_{T_1} \nabla^2 n_{T_1} - k_{\text{TTA}}(\mathbf{r}) n_{T_1}^2 - k_{\text{T-decay}} n_{T_1}
$$

$$
\frac{\partial n_{S_1}(\mathbf{r},t)}{\partial t} = D_{S_1} \nabla^2 n_{S_1} + \frac{1}{2}k_{\text{TTA}}(\mathbf{r}) n_{T_1}^2 - \Gamma_{\text{fl}} n_{S_1}
$$

ここで：
- $n_{X}(\mathbf{r},t)$ は位置 $\mathbf{r}$、時刻 $t$ における状態 $X$ の数密度
- $D_X$ は状態 $X$ の拡散定数
- $\nabla^2$ はラプラシアン演算子

### 10.2 拡散制御TTA

TTA速度定数の空間依存性は、分子の拡散と反応の競合により決まる：

$$
k_{\text{TTA}}(\mathbf{r}) = k_{\text{TTA}}^0 \cdot f(n_{T_1}(\mathbf{r}))
$$

ここで、$f$ は三重項密度依存関数である。

## 11. 観測量

### 11.1 蛍光強度

時刻 $t$ における蛍光強度 $I_{\text{fl}}(t)$ は：

$$
I_{\text{fl}}(t) = \Gamma_{\text{fl}} \int d\mathbf{r} \, n_{S_1}(\mathbf{r},t)
$$

または、空間積分した系では：

$$
I_{\text{fl}}(t) = \Gamma_{\text{fl}} N_{S_1}(t)
$$

### 11.2 時間積分蛍光強度

全蛍光発光量は：

$$
\mathcal{I}_{\text{total}} = \int_0^\infty I_{\text{fl}}(t) dt = \Gamma_{\text{fl}} \int_0^\infty N_{S_1}(t) dt
$$

### 11.3 遅延蛍光

TTA過程により生成される一重項状態からの発光は遅延蛍光（delayed fluorescence）と呼ばれる。三重項寿命が通常マイクロ秒オーダーであるため、遅延蛍光は通常の蛍光（ナノ秒オーダー）より長い時間スケールで観測される。

遅延蛍光強度の時間依存性は：

$$
I_{\text{DF}}(t) \propto [N_{T_1}(t)]^2
$$

## 12. エネルギー図とプロセスの概要

### 12.1 エネルギー準位図

```
エネルギー
    ↑
    |
E_S |     ───── S₁ (励起一重項)
    |      ↑↓ Γ_fl (蛍光)
    |      |
E_T |     ───── T₁ (励起三重項)
    |      ⇄ V_ij (エネルギー移動)
    |     / \
    |    /   \ k_TTA (TTA)
    |   /     \
  0 |  ───── S₀ (基底一重項)
    |
    └─────────────→
```

### 12.2 プロセスのまとめ

1. **三重項エネルギー移動**: $|T_1\rangle_i |S_0\rangle_j \leftrightarrow |S_0\rangle_i |T_1\rangle_j$
   - 速度定数: $k_{\text{transfer}} \propto |V_{ij}|^2$

2. **三重項-三重項消滅**: $|T_1\rangle_i |T_1\rangle_j \to |S_1\rangle_i |S_0\rangle_j$ or $|S_0\rangle_i |S_1\rangle_j$
   - 速度定数: $k_{\text{TTA}} \propto |J_{ij}|^2$
   - エネルギー条件: $2E_T \approx E_S$

3. **蛍光発光**: $|S_1\rangle \to |S_0\rangle + h\nu$
   - 放出速度: $\Gamma_{\text{fl}} = A_{S_1 \to S_0}$
   - 蛍光寿命: $\tau_{\text{fl}} = 1/\Gamma_{\text{fl}}$

## 13. 結論

本文書では、基底一重項状態、励起三重項状態、励起一重項状態を持つ分子系における量子ダイナミクスを、完全な数式展開とともに記述した。特に以下の過程を詳細に解析した：

1. 隣接分子間の三重項励起エネルギー移動
2. 三重項-三重項消滅による一重項生成
3. 一重項状態からの蛍光発光

これらの過程は、有機フォトニクス、光アップコンバージョン、有機ELデバイスなど、多くの光物理・光化学現象において重要な役割を果たしている。

## 参考文献

本理論的枠組みは、以下の物理学・化学の基礎理論に基づいている：

- 量子力学の時間依存摂動論（フェルミの黄金律）
- 開放量子系の理論（リンドブラッド方程式）
- 光と物質の相互作用（量子電磁力学）
- 化学反応動力学（速度方程式論）
- 三重項-三重項消滅過程の理論
- 遅延蛍光の物理化学

---

**文書作成日**: 2025-10-14  
**分野**: 量子ダイナミクス、分子光物理学、励起状態ダイナミクス  
**対象**: MQT Qudits フレームワークにおける量子系シミュレーション
