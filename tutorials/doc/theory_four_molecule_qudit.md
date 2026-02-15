# 4分子直線配置モデルにおける分子三重項状態の量子ダイナミクス完全理論

## 1. はじめに

### 1.1 本文書の目的

本文書では、4つの分子が一次元直線上に配列した系における励起三重項状態（triplet state）の量子ダイナミクスを、Quditベースの量子アルゴリズムを用いて記述するための完全な理論的基礎を提供する。特に、両端の分子が励起三重項状態、中間の分子が基底状態にある初期状態からの時間発展を詳細に定式化する。

### 1.2 物理系の概要

考察する系は以下の特徴を持つ：

- **分子数**: $N = 4$
- **空間配置**: 一次元直線配列、最近接相互作用
- **各分子の電子状態**: 3つの状態（基底一重項、励起三重項、励起一重項）
- **エネルギー移動過程**: 三重項間エネルギー移動（Triplet Energy Transfer）、三重項-三重項消滅（Triplet-Triplet Annihilation）

### 1.3 Qudit表現の必然性

各分子が3つの電子状態を持つため、量子ビット（qubit, 2準位系）では表現できない。したがって、3準位系であるQutrit（Qudit, $d=3$）を用いる必要がある。4分子系は4つのQutritで表現され、ヒルベルト空間の次元は $3^4 = 81$ となる。

## 2. 分子の電子状態とエネルギー準位

### 2.1 単一分子の電子状態

各分子 $i$ （$i = 0, 1, 2, 3$）は以下の3つの電子状態を持つ：

#### 基底一重項状態（Ground Singlet State）

$$
|S_0\rangle_i
$$

- **エネルギー**: $E_{S_0} = 0$ （基準エネルギー）
- **スピン多重度**: 1（一重項）
- **電子配置**: 全ての電子が基底軌道に対占有

#### 励起三重項状態（Excited Triplet State）

$$
|T_1\rangle_i
$$

- **エネルギー**: $E_{T_1} = E_T = 1.5 \text{ eV}$
- **スピン多重度**: 3（三重項）
- **電子配置**: 1つの電子が励起軌道に昇位、スピンが平行

#### 励起一重項状態（Excited Singlet State）

$$
|S_1\rangle_i
$$

- **エネルギー**: $E_{S_1} = E_S = 3.0 \text{ eV}$
- **スピン多重度**: 1（一重項）
- **電子配置**: 2つの電子が励起軌道に対占有

### 2.2 エネルギー関係式

実験的に観測される典型的なエネルギー関係は：

$$
E_{S_1} \approx 2 E_{T_1}
$$

これは三重項-三重項消滅（TTA）過程において重要な役割を果たす：

$$
T_1 + T_1 \rightarrow S_1 + S_0
$$

この過程はエネルギー保存則：

$$
2 E_{T_1} = E_{S_1} + E_{S_0} = E_{S_1}
$$

を満たす。

## 3. Qutrit表現とヒルベルト空間

### 3.1 単一分子のQutrit表現

各分子の電子状態をQutrit基底に対応させる：

$$
\begin{align}
|S_0\rangle_i &\longleftrightarrow |0\rangle_i \\
|T_1\rangle_i &\longleftrightarrow |1\rangle_i \\
|S_1\rangle_i &\longleftrightarrow |2\rangle_i
\end{align}
$$

Qutrit基底は以下の性質を満たす：

$$
\langle m | n \rangle = \delta_{mn}, \quad m, n \in \{0, 1, 2\}
$$

$$
\sum_{n=0}^{2} |n\rangle \langle n| = \mathbb{I}_3
$$

### 3.2 4分子系の全ヒルベルト空間

4つの分子からなる系の全ヒルベルト空間は、各Qutritのテンソル積として構成される：

$$
\mathcal{H}_{\text{total}} = \mathcal{H}_0 \otimes \mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \mathcal{H}_3
$$

各 $\mathcal{H}_i$ は3次元ヒルベルト空間であり：

$$
\mathcal{H}_i = \text{span}\{|0\rangle_i, |1\rangle_i, |2\rangle_i\}
$$

全系の次元は：

$$
\dim(\mathcal{H}_{\text{total}}) = 3 \times 3 \times 3 \times 3 = 81
$$

### 3.3 基底状態の表記

全系の基底状態は：

$$
|n_0 n_1 n_2 n_3\rangle = |n_0\rangle \otimes |n_1\rangle \otimes |n_2\rangle \otimes |n_3\rangle
$$

ここで、$n_i \in \{0, 1, 2\}$ は分子 $i$ の状態を表す。

例えば：

- $|0000\rangle$: 全ての分子が基底状態
- $|1111\rangle$: 全ての分子が三重項状態
- $|1001\rangle$: 両端（分子0と分子3）が三重項状態、中間（分子1, 2）が基底状態

### 3.4 基底の完全性

81個の基底状態は完全系を形成する：

$$
\sum_{n_0, n_1, n_2, n_3 = 0}^{2} |n_0 n_1 n_2 n_3\rangle \langle n_0 n_1 n_2 n_3| = \mathbb{I}_{81}
$$

任意の状態ベクトル $|\Psi\rangle$ は：

$$
|\Psi\rangle = \sum_{n_0, n_1, n_2, n_3 = 0}^{2} c_{n_0 n_1 n_2 n_3} |n_0 n_1 n_2 n_3\rangle
$$

と展開される。規格化条件は：

$$
\sum_{n_0, n_1, n_2, n_3 = 0}^{2} |c_{n_0 n_1 n_2 n_3}|^2 = 1
$$

## 4. ハミルトニアンの定式化

### 4.1 全ハミルトニアン

系の全ハミルトニアンは3つの項の和として表される：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}}
$$

各項は物理的に異なるプロセスを記述する。

### 4.2 オンサイトエネルギー項 $\hat{H}_0$

#### 4.2.1 物理的意味

$\hat{H}_0$ は各分子の固有エネルギーを記述する対角項である。

#### 4.2.2 数学的表現

$$
\hat{H}_0 = \sum_{i=0}^{3} \left( E_{T_1} |1\rangle_i \langle 1|_i + E_{S_1} |2\rangle_i \langle 2|_i \right)
$$

基底状態 $|0\rangle$ のエネルギーは $E_{S_0} = 0$ と設定しているため、対応する項は省略されている。

#### 4.2.3 行列表現

単一分子 $i$ のオンサイトハミルトニアンは：

$$
\hat{h}_0^{(i)} = E_{T_1} |1\rangle_i \langle 1|_i + E_{S_1} |2\rangle_i \langle 2|_i
$$

Qutrit基底 $\{|0\rangle, |1\rangle, |2\rangle\}$ での行列表現は：

$$
\hat{h}_0^{(i)} = \begin{pmatrix}
0 & 0 & 0 \\
0 & E_{T_1} & 0 \\
0 & 0 & E_{S_1}
\end{pmatrix}
$$

全系のハミルトニアンは：

$$
\hat{H}_0 = \sum_{i=0}^{3} \mathbb{I}_0 \otimes \cdots \otimes \hat{h}_0^{(i)} \otimes \cdots \otimes \mathbb{I}_3
$$

ここで、$\hat{h}_0^{(i)}$ は分子 $i$ の位置にあり、他の位置には単位行列 $\mathbb{I}_3$ が配置される。

#### 4.2.4 具体例

分子0のエネルギー項：

$$
\hat{H}_0^{(0)} = \hat{h}_0^{(0)} \otimes \mathbb{I}_3 \otimes \mathbb{I}_3 \otimes \mathbb{I}_3
$$

これは $|1 n_1 n_2 n_3\rangle$ に作用すると：

$$
\hat{H}_0^{(0)} |1 n_1 n_2 n_3\rangle = E_{T_1} |1 n_1 n_2 n_3\rangle
$$

### 4.3 三重項エネルギー移動項 $\hat{H}_{\text{transfer}}$

#### 4.3.1 物理的過程

三重項エネルギー移動（Triplet Energy Transfer, TET）は、隣接する分子間で励起三重項状態が移動する過程である：

$$
T_1 + S_0 \leftrightarrow S_0 + T_1
$$

量子状態で表現すると：

$$
|1\rangle_i |0\rangle_{i+1} \leftrightarrow |0\rangle_i |1\rangle_{i+1}
$$

#### 4.3.2 デクスターメカニズム

この過程はデクスター電子交換機構（Dexter mechanism）によって媒介され、そのレートは分子間距離 $r$ に対して指数関数的に減衰する：

$$
J(r) = J_0 \exp(-\alpha r)
$$

ここで、$J_0$ は接触時の移動積分、$\alpha$ は減衰定数である。

最近接距離を $r_0 = 5$ Å とすると：

$$
J_{\text{transfer}} = J_0 \exp(-\alpha r_0)
$$

典型的な値として $J_{\text{transfer}} = 0.01$ eV を採用する。

#### 4.3.3 数学的表現

最近接相互作用として：

$$
\hat{H}_{\text{transfer}} = -J_{\text{transfer}} \sum_{i=0}^{2} \left( |0\rangle_i \langle 1|_i \otimes |1\rangle_{i+1} \langle 0|_{i+1} + |1\rangle_i \langle 0|_i \otimes |0\rangle_{i+1} \langle 1|_{i+1} \right)
$$

符号が負であることに注意。これは結合状態が非結合状態よりも低エネルギーであることを反映している。

#### 4.3.4 展開形式

各項を明示的に書くと：

$$
\begin{align}
\hat{H}_{\text{transfer}} = &-J_{\text{transfer}} \Big[
|0\rangle_0 \langle 1|_0 \otimes |1\rangle_1 \langle 0|_1 \otimes \mathbb{I}_2 \otimes \mathbb{I}_3 \\
&+ |1\rangle_0 \langle 0|_0 \otimes |0\rangle_1 \langle 1|_1 \otimes \mathbb{I}_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes |0\rangle_1 \langle 1|_1 \otimes |1\rangle_2 \langle 0|_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes |1\rangle_1 \langle 0|_1 \otimes |0\rangle_2 \langle 1|_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes \mathbb{I}_1 \otimes |0\rangle_2 \langle 1|_2 \otimes |1\rangle_3 \langle 0|_3 \\
&+ \mathbb{I}_0 \otimes \mathbb{I}_1 \otimes |1\rangle_2 \langle 0|_2 \otimes |0\rangle_3 \langle 1|_3 \Big]
\end{align}
$$

#### 4.3.5 部分空間での作用

$\hat{H}_{\text{transfer}}$ は、隣接する2つの分子が $(|01\rangle, |10\rangle)$ の形式をとる2次元部分空間でのみ非ゼロの行列要素を持つ：

$$
\begin{pmatrix}
|01\rangle \\
|10\rangle
\end{pmatrix}
\rightarrow
\begin{pmatrix}
0 & -J_{\text{transfer}} \\
-J_{\text{transfer}} & 0
\end{pmatrix}
\begin{pmatrix}
|01\rangle \\
|10\rangle
\end{pmatrix}
$$

この2×2行列は、実対称行列であり、固有値は：

$$
E_{\pm} = \pm J_{\text{transfer}}
$$

固有状態は：

$$
\begin{align}
|\phi_+\rangle &= \frac{1}{\sqrt{2}}(|01\rangle + |10\rangle) \quad \text{(bonding)} \\
|\phi_-\rangle &= \frac{1}{\sqrt{2}}(|01\rangle - |10\rangle) \quad \text{(antibonding)}
\end{align}
$$

### 4.4 三重項-三重項消滅項 $\hat{H}_{\text{TTA}}$

#### 4.4.1 物理的過程

三重項-三重項消滅（Triplet-Triplet Annihilation, TTA）は、2つの三重項励起が相互作用して、1つの励起一重項と1つの基底状態を生成する過程である：

$$
T_1 + T_1 \rightarrow S_1 + S_0
$$

量子状態で表現すると：

$$
|1\rangle_i |1\rangle_{i+1} \rightarrow |2\rangle_i |0\rangle_{i+1} + |0\rangle_i |2\rangle_{i+1}
$$

#### 4.4.2 エネルギー保存

エネルギー保存則：

$$
2 E_{T_1} = E_{S_1} + E_{S_0}
$$

を満たす共鳴過程であるため、高い反応確率を持つ。

#### 4.4.3 数学的表現

最近接相互作用として：

$$
\begin{align}
\hat{H}_{\text{TTA}} = -J_{\text{TTA}} \sum_{i=0}^{2} \Big[
&|0\rangle_i \langle 1|_i \otimes |2\rangle_{i+1} \langle 1|_{i+1} \\
&+ |2\rangle_i \langle 1|_i \otimes |0\rangle_{i+1} \langle 1|_{i+1} \\
&+ |1\rangle_i \langle 0|_i \otimes |1\rangle_{i+1} \langle 2|_{i+1} \\
&+ |1\rangle_i \langle 2|_i \otimes |1\rangle_{i+1} \langle 0|_{i+1} \Big]
\end{align}
$$

符号が負であることに注意。

#### 4.4.4 相互作用強度

TTA過程の相互作用強度 $J_{\text{TTA}}$ は、2電子過程であるため一般的にTETよりも強い：

$$
J_{\text{TTA}} = 0.05 \text{ eV} = 5 \times J_{\text{transfer}}
$$

#### 4.4.5 展開形式

各項を明示的に書くと：

$$
\begin{align}
\hat{H}_{\text{TTA}} = &-J_{\text{TTA}} \Big[
|0\rangle_0 \langle 1|_0 \otimes |2\rangle_1 \langle 1|_1 \otimes \mathbb{I}_2 \otimes \mathbb{I}_3 \\
&+ |2\rangle_0 \langle 1|_0 \otimes |0\rangle_1 \langle 1|_1 \otimes \mathbb{I}_2 \otimes \mathbb{I}_3 \\
&+ |1\rangle_0 \langle 0|_0 \otimes |1\rangle_1 \langle 2|_1 \otimes \mathbb{I}_2 \otimes \mathbb{I}_3 \\
&+ |1\rangle_0 \langle 2|_0 \otimes |1\rangle_{1} \langle 0|_1 \otimes \mathbb{I}_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes |0\rangle_1 \langle 1|_1 \otimes |2\rangle_2 \langle 1|_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes |2\rangle_1 \langle 1|_1 \otimes |0\rangle_2 \langle 1|_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes |1\rangle_1 \langle 0|_1 \otimes |1\rangle_2 \langle 2|_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes |1\rangle_1 \langle 2|_1 \otimes |1\rangle_2 \langle 0|_2 \otimes \mathbb{I}_3 \\
&+ \mathbb{I}_0 \otimes \mathbb{I}_1 \otimes |0\rangle_2 \langle 1|_2 \otimes |2\rangle_3 \langle 1|_3 \\
&+ \mathbb{I}_0 \otimes \mathbb{I}_1 \otimes |2\rangle_2 \langle 1|_2 \otimes |0\rangle_3 \langle 1|_3 \\
&+ \mathbb{I}_0 \otimes \mathbb{I}_1 \otimes |1\rangle_2 \langle 0|_2 \otimes |1\rangle_3 \langle 2|_3 \\
&+ \mathbb{I}_0 \otimes \mathbb{I}_1 \otimes |1\rangle_2 \langle 2|_2 \otimes |1\rangle_3 \langle 0|_3 \Big]
\end{align}
$$

#### 4.4.6 部分空間での作用

$\hat{H}_{\text{TTA}}$ は、隣接する2つの分子が $(|02\rangle, |11\rangle, |20\rangle)$ の形式をとる3次元部分空間でのみ非ゼロの行列要素を持つ：

$$
\begin{pmatrix}
|02\rangle \\
|11\rangle \\
|20\rangle
\end{pmatrix}
\rightarrow
\begin{pmatrix}
0 & -J_{\text{TTA}} & 0 \\
-J_{\text{TTA}} & 0 & -J_{\text{TTA}} \\
0 & -J_{\text{TTA}} & 0
\end{pmatrix}
\begin{pmatrix}
|02\rangle \\
|11\rangle \\
|20\rangle
\end{pmatrix}
$$

この3×3行列は実対称行列であり、固有値と固有ベクトルは解析的に求められる。

固有値：

$$
\begin{align}
E_0 &= 0 \\
E_{\pm} &= \pm \sqrt{2} J_{\text{TTA}}
\end{align}
$$

固有状態：

$$
\begin{align}
|\psi_0\rangle &= \frac{1}{\sqrt{2}}(|02\rangle - |20\rangle) \\
|\psi_+\rangle &= \frac{1}{2}(|02\rangle + \sqrt{2}|11\rangle + |20\rangle) \\
|\psi_-\rangle &= \frac{1}{2}(|02\rangle - \sqrt{2}|11\rangle + |20\rangle)
\end{align}
$$

### 4.5 ハミルトニアンの性質

#### 4.5.1 エルミート性

全ハミルトニアンはエルミート演算子である：

$$
\hat{H}_{\text{total}}^{\dagger} = \hat{H}_{\text{total}}
$$

各項もエルミート性を満たす：

$$
\hat{H}_0^{\dagger} = \hat{H}_0, \quad
\hat{H}_{\text{transfer}}^{\dagger} = \hat{H}_{\text{transfer}}, \quad
\hat{H}_{\text{TTA}}^{\dagger} = \hat{H}_{\text{TTA}}
$$

#### 4.5.2 実対称性

Qutrit基底での行列表現は実対称行列である。

#### 4.5.3 疎行列構造

全ハミルトニアンの $81 \times 81$ 行列表現は高度に疎である。各項の非ゼロ要素数：

- $\hat{H}_0$: 対角要素のみ、最大 $81$ 個
- $\hat{H}_{\text{transfer}}$: 各相互作用対につき2要素、最大 $6 \times 2 = 12$ 個
- $\hat{H}_{\text{TTA}}$: 各相互作用対につき4要素、最大 $12 \times 4 = 48$ 個

合計でも非ゼロ要素は $81 + 12 + 48 = 141$ 程度であり、全要素 $81^2 = 6561$ の約2%である。

## 5. 初期状態の設定

### 5.1 両端励起初期状態

本研究では、以下の初期状態を考察する：

$$
|\Psi(0)\rangle = |1001\rangle = |1\rangle_0 \otimes |0\rangle_1 \otimes |0\rangle_2 \otimes |1\rangle_3
$$

これは：

- 分子0（左端）: 励起三重項状態 $|T_1\rangle$
- 分子1（中央左）: 基底一重項状態 $|S_0\rangle$
- 分子2（中央右）: 基底一重項状態 $|S_0\rangle$
- 分子3（右端）: 励起三重項状態 $|T_1\rangle$

を表す。

### 5.2 初期状態の物理的意味

この初期状態は、系の両端に局在した励起を持つ非平衡状態である。時間発展により：

1. **三重項エネルギー移動**: 励起が中央に向かって移動
2. **三重項-三重項消滅**: 励起が出会うとTTAが発生
3. **エネルギー再分配**: 一重項励起の生成と緩和

という複雑なダイナミクスが展開される。

### 5.3 量子ゲートによる準備

MQT-Quditsフレームワークでは、初期状態を以下のように準備する：

1. **真空状態の準備**: $|0000\rangle$（自動的に準備される）
2. **Xゲートの適用**:
   - 分子0に $\hat{X}$ を適用: $|0000\rangle \rightarrow |1000\rangle$
   - 分子3に $\hat{X}$ を適用: $|1000\rangle \rightarrow |1001\rangle$

ここで、Xゲート（昇降演算子）は：

$$
\hat{X} |n\rangle = |n+1 \mod 3\rangle
$$

すなわち：

$$
\hat{X}|0\rangle = |1\rangle, \quad \hat{X}|1\rangle = |2\rangle, \quad \hat{X}|2\rangle = |0\rangle
$$

### 5.4 初期状態のエネルギー

初期状態のエネルギー期待値は：

$$
\langle \Psi(0) | \hat{H}_{\text{total}} | \Psi(0) \rangle = \langle 1001 | \hat{H}_0 | 1001 \rangle = 2 E_{T_1} = 3.0 \text{ eV}
$$

$\hat{H}_{\text{transfer}}$ と $\hat{H}_{\text{TTA}}$ は $|1001\rangle$ に対して対角要素を持たないため寄与しない。

### 5.5 初期個体数

各電子状態の個体数演算子を：

$$
\hat{N}_{S_0} = \sum_{i=0}^{3} |0\rangle_i \langle 0|_i, \quad
\hat{N}_{T_1} = \sum_{i=0}^{3} |1\rangle_i \langle 1|_i, \quad
\hat{N}_{S_1} = \sum_{i=0}^{3} |2\rangle_i \langle 2|_i
$$

と定義すると、初期状態での期待値は：

$$
\begin{align}
N_{S_0}(0) &= \langle \Psi(0) | \hat{N}_{S_0} | \Psi(0) \rangle = 2 \\
N_{T_1}(0) &= \langle \Psi(0) | \hat{N}_{T_1} | \Psi(0) \rangle = 2 \\
N_{S_1}(0) &= \langle \Psi(0) | \hat{N}_{S_1} | \Psi(0) \rangle = 0
\end{align}
$$

個体数保存則：

$$
N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = 4
$$

が全時刻で成立する。

## 6. 時間発展演算子とシュレーディンガー方程式

### 6.1 時間依存シュレーディンガー方程式

系の量子状態 $|\Psi(t)\rangle$ は時間依存シュレーディンガー方程式に従う：

$$
i\hbar \frac{\partial}{\partial t} |\Psi(t)\rangle = \hat{H}_{\text{total}} |\Psi(t)\rangle
$$

ハミルトニアンが時間に依存しない場合、形式解は：

$$
|\Psi(t)\rangle = \hat{U}(t) |\Psi(0)\rangle
$$

ここで、時間発展演算子は：

$$
\hat{U}(t) = \exp\left(-\frac{i}{\hbar} \hat{H}_{\text{total}} t\right)
$$

### 6.2 時間発展演算子の性質

#### 6.2.1 ユニタリ性

$$
\hat{U}^{\dagger}(t) \hat{U}(t) = \hat{U}(t) \hat{U}^{\dagger}(t) = \mathbb{I}
$$

これにより、状態ベクトルの規格化が保存される：

$$
\langle \Psi(t) | \Psi(t) \rangle = \langle \Psi(0) | \hat{U}^{\dagger}(t) \hat{U}(t) | \Psi(0) \rangle = \langle \Psi(0) | \Psi(0) \rangle = 1
$$

#### 6.2.2 群性質

$$
\hat{U}(t_1 + t_2) = \hat{U}(t_1) \hat{U}(t_2)
$$

$$
\hat{U}(0) = \mathbb{I}
$$

$$
\hat{U}^{-1}(t) = \hat{U}(-t) = \hat{U}^{\dagger}(t)
$$

### 6.3 単位系とプランク定数

本研究では、自然単位系を採用し：

$$
\hbar = 1
$$

と設定する。時間の単位はフェムト秒（fs, $10^{-15}$ s）、エネルギーの単位は電子ボルト（eV）を用いる。

換算関係：

$$
\hbar = 0.6582 \text{ eV} \cdot \text{fs}
$$

自然単位では時間発展演算子は：

$$
\hat{U}(t) = \exp(-i \hat{H}_{\text{total}} t)
$$

## 7. 鈴木トロッター分解による時間発展の数値計算

### 7.1 問題の本質

ハミルトニアンが複数の非可換な項の和である場合：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}}
$$

一般に以下は成立しない：

$$
e^{-i(\hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}})t} \neq e^{-i\hat{H}_0 t} e^{-i\hat{H}_{\text{transfer}} t} e^{-i\hat{H}_{\text{TTA}} t}
$$

なぜなら：

$$
[\hat{H}_0, \hat{H}_{\text{transfer}}] \neq 0, \quad
[\hat{H}_0, \hat{H}_{\text{TTA}}] \neq 0, \quad
[\hat{H}_{\text{transfer}}, \hat{H}_{\text{TTA}}] \neq 0
$$

### 7.2 Baker-Campbell-Hausdorff公式

2つの演算子 $\hat{A}$, $\hat{B}$ に対して：

$$
e^{\hat{A}} e^{\hat{B}} = \exp\left(\hat{A} + \hat{B} + \frac{1}{2}[\hat{A}, \hat{B}] + \frac{1}{12}[\hat{A}, [\hat{A}, \hat{B}]] - \frac{1}{12}[\hat{B}, [\hat{A}, \hat{B}]] + \cdots\right)
$$

交換子が非ゼロの場合、高次の補正項が無限に続く。

### 7.3 1次Lie-Trotter分解

時間区間を $N$ 個に分割し、各時間刻みを $\Delta t = T/N$ とする。

1次の分解公式：

$$
e^{-i(\hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}})\Delta t} \approx e^{-i\hat{H}_0 \Delta t} e^{-i\hat{H}_{\text{transfer}} \Delta t} e^{-i\hat{H}_{\text{TTA}} \Delta t}
$$

誤差は $\mathcal{O}(\Delta t^2)$ である。

全時間発展：

$$
\hat{U}(T) \approx \left(e^{-i\hat{H}_0 \Delta t} e^{-i\hat{H}_{\text{transfer}} \Delta t} e^{-i\hat{H}_{\text{TTA}} \Delta t}\right)^N
$$

累積誤差は $\mathcal{O}(T \Delta t) = \mathcal{O}(T^2/N)$ である。

### 7.4 2次Suzuki-Trotter分解（対称分解）

より高い精度を得るために、対称な分解を用いる：

$$
\begin{align}
\hat{U}_{\text{Trotter}}(\Delta t) = &\exp\left(-i\hat{H}_0 \frac{\Delta t}{2}\right)
\exp\left(-i\hat{H}_{\text{transfer}} \frac{\Delta t}{2}\right) \\
&\times \exp\left(-i\hat{H}_{\text{TTA}} \Delta t\right) \\
&\times \exp\left(-i\hat{H}_{\text{transfer}} \frac{\Delta t}{2}\right)
\exp\left(-i\hat{H}_0 \frac{\Delta t}{2}\right)
\end{align}
$$

この分解の局所誤差は $\mathcal{O}(\Delta t^3)$ であり、累積誤差は $\mathcal{O}(T \Delta t^2) = \mathcal{O}(T^3/N^2)$ となる。

### 7.5 具体的な実装手順

1ステップの時間発展 $\hat{U}_{\text{Trotter}}(\Delta t)$ を構成するために、以下の順序で演算子を適用する：

1. $\exp(-i\hat{H}_0 \Delta t/2)$
2. $\exp(-i\hat{H}_{\text{transfer}} \Delta t/2)$
3. $\exp(-i\hat{H}_{\text{TTA}} \Delta t)$
4. $\exp(-i\hat{H}_{\text{transfer}} \Delta t/2)$
5. $\exp(-i\hat{H}_0 \Delta t/2)$

全時間発展は：

$$
|\Psi(T)\rangle \approx \hat{U}_{\text{Trotter}}(\Delta t)^N |\Psi(0)\rangle
$$

### 7.6 誤差評価

2次Suzuki-Trotter分解の誤差は：

$$
\|\hat{U}(T) - \hat{U}_{\text{Trotter}}(\Delta t)^N\| = \mathcal{O}\left(\frac{T^3}{N^2}\right)
$$

精度を1桁上げるには、ステップ数を約3倍にする必要がある。

典型的には、$N = 20 \sim 100$ ステップで十分な精度が得られる。

## 8. MQT-Quditsによる量子ゲート実装

### 8.1 量子ゲート分解の原理

各時間発展演算子 $\exp(-i\hat{H}_j t)$ を、MQT-Quditsフレームワークで利用可能な基本量子ゲートの列に分解する。

### 8.2 オンサイト項 $\exp(-i\hat{H}_0 t)$ の実装

$\hat{H}_0$ は対角演算子であり、各Qutritに独立に作用する：

$$
\exp(-i\hat{H}_0 t) = \bigotimes_{i=0}^{3} \exp(-i\hat{h}_0^{(i)} t)
$$

単一Qutritの時間発展：

$$
\exp(-i\hat{h}_0^{(i)} t) = \mathbb{I} + (e^{-i E_{T_1} t} - 1)|1\rangle\langle 1| + (e^{-i E_{S_1} t} - 1)|2\rangle\langle 2|
$$

これは仮想Zゲート（VirtRz）で実装される：

$$
\text{VirtRz}(q, l, \theta) = e^{i\theta} |l\rangle\langle l|_q
$$

分子 $i$ に対して：

- `VirtRz(i, 1, -E_T1 * t)`: 準位1に位相 $-E_{T_1} t$ を付加
- `VirtRz(i, 2, -E_S1 * t)`: 準位2に位相 $-E_{S_1} t$ を付加

### 8.3 2体相互作用項の実装

#### 8.3.1 CustomTwoゲート

MQT-Quditsでは、2つのQutrit間の任意のユニタリ変換を記述する`CustomTwo`ゲートが提供される：

$$
\text{CustomTwo}(q_1, q_2, U) : U \in U(9)
$$

ここで、$U$ は $9 \times 9$ ユニタリ行列である。

#### 8.3.2 疎構造を利用した最適化

本研究の $\hat{H}_{\text{transfer}}$ と $\hat{H}_{\text{TTA}}$ は、2つのQutritの $3 \times 3 = 9$ 次元空間のうち、それぞれ2次元と3次元の部分空間でのみ非自明な作用を持つ。

疎構造認識コンパイラ（Sparse-Aware Compiler）は、この構造を自動検出し、効率的な分解を実行する：

- **$\hat{H}_{\text{transfer}}$**: 2次元部分空間 $\rightarrow$ 1個の基本ゲート
- **$\hat{H}_{\text{TTA}}$**: 3次元部分空間 $\rightarrow$ 6個の基本ゲート

従来の一般的な分解では各CustomTwoゲートあたり約1000個の基本ゲートを要していたため、劇的な効率化が達成されている。

#### 8.3.3 $\hat{H}_{\text{transfer}}$ の実装

隣接Qutrit対 $(i, i+1)$ に対して、2次元部分空間 $(|01\rangle, |10\rangle)$ での時間発展：

$$
U_{\text{transfer}}(t) = \begin{pmatrix}
\cos(J_{\text{transfer}} t) & -i\sin(J_{\text{transfer}} t) \\
-i\sin(J_{\text{transfer}} t) & \cos(J_{\text{transfer}} t)
\end{pmatrix}
$$

これは単一の制御励起ゲート（CEx）で実装できる：

$$
\text{CEx}(q_i, q_{i+1}, 0, 1, -J_{\text{transfer}} t)
$$

CExゲートは、制御Qudit $q_i$ が状態 $|c\rangle$ にあるとき、ターゲットQudit $q_{i+1}$ の準位 $|l\rangle$ と $|l+1\rangle$ の間で回転を実行する。

#### 8.3.4 $\hat{H}_{\text{TTA}}$ の実装

隣接Qutrit対 $(i, i+1)$ に対して、3次元部分空間 $(|02\rangle, |11\rangle, |20\rangle)$ での時間発展を記述するユニタリ行列を構成し、それを6個の基本ゲート（R, Rz, CEx）の列に分解する。

具体的な分解は、QR分解、Givens回転、ZYZ分解などの厳密な数学的手法を組み合わせて実行される。

## 9. 厳密対角化による解析解

### 9.1 ハミルトニアンの対角化

ハミルトニアン $\hat{H}_{\text{total}}$ の固有値問題を解く：

$$
\hat{H}_{\text{total}} |\phi_k\rangle = E_k |\phi_k\rangle, \quad k = 0, 1, \ldots, 80
$$

数値的には、$81 \times 81$ 実対称行列の固有値分解を実行する：

$$
\hat{H}_{\text{total}} = \sum_{k=0}^{80} E_k |\phi_k\rangle \langle \phi_k|
$$

固有状態は正規直交基底を形成する：

$$
\langle \phi_k | \phi_l \rangle = \delta_{kl}
$$

### 9.2 固有基底での時間発展

初期状態を固有基底で展開する：

$$
|\Psi(0)\rangle = \sum_{k=0}^{80} c_k |\phi_k\rangle
$$

ここで、展開係数は：

$$
c_k = \langle \phi_k | \Psi(0) \rangle
$$

時間発展は各固有状態に独立な位相因子を掛けるだけで実行される：

$$
|\Psi(t)\rangle = \sum_{k=0}^{80} c_k e^{-i E_k t} |\phi_k\rangle
$$

これは**厳密解**である。

### 9.3 計算複雑度

固有値分解の計算量は $\mathcal{O}(81^3) \approx 5 \times 10^5$ である。一度対角化すれば、任意の時刻での状態計算は $\mathcal{O}(81^2) \approx 6500$ の計算量で実行できる。

### 9.4 Quditシミュレーションとの比較

厳密対角化による解析解は、Quditシミュレーションの精度を検証するための基準として使用される。忠実度（fidelity）：

$$
F(t) = |\langle \Psi_{\text{exact}}(t) | \Psi_{\text{Qudit}}(t) \rangle|^2
$$

が1に近いほど、Quditシミュレーションが正確であることを示す。

## 10. 観測可能量の計算

### 10.1 個体数演算子

各電子状態の個体数演算子：

$$
\hat{N}_{S_0} = \sum_{i=0}^{3} |0\rangle_i \langle 0|_i
$$

$$
\hat{N}_{T_1} = \sum_{i=0}^{3} |1\rangle_i \langle 1|_i
$$

$$
\hat{N}_{S_1} = \sum_{i=0}^{3} |2\rangle_i \langle 2|_i
$$

### 10.2 期待値の時間発展

時刻 $t$ での期待値：

$$
N_{S_0}(t) = \langle \Psi(t) | \hat{N}_{S_0} | \Psi(t) \rangle
$$

$$
N_{T_1}(t) = \langle \Psi(t) | \hat{N}_{T_1} | \Psi(t) \rangle
$$

$$
N_{S_1}(t) = \langle \Psi(t) | \hat{N}_{S_1} | \Psi(t) \rangle
$$

### 10.3 保存則

個体数保存則：

$$
N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = 4 \quad \forall t
$$

これはハミルトニアンが各分子の占有数を変更せず、状態間での遷移のみを引き起こすことから保証される。

### 10.4 エネルギー保存

ハミルトニアンが時間に依存しないため、エネルギーは厳密に保存される：

$$
\langle \Psi(t) | \hat{H}_{\text{total}} | \Psi(t) \rangle = \langle \Psi(0) | \hat{H}_{\text{total}} | \Psi(0) \rangle = 2 E_{T_1} = 3.0 \text{ eV}
$$

## 11. 予想される物理的挙動

### 11.1 初期過程（$t \sim 0$-20 fs）

両端の三重項励起が、三重項エネルギー移動により中央に向かって拡散する。

$$
|1001\rangle \rightarrow \alpha|1001\rangle + \beta(|0101\rangle + |1010\rangle) + \cdots
$$

相互作用強度 $J_{\text{transfer}} = 0.01$ eV から、特性時間は：

$$
\tau_{\text{transfer}} \sim \frac{\hbar}{J_{\text{transfer}}} \sim \frac{0.66 \text{ eV} \cdot \text{fs}}{0.01 \text{ eV}} \sim 66 \text{ fs}
$$

### 11.2 中期過程（$t \sim$ 20-50 fs）

中央の分子（1, 2）に励起が到達すると、三重項-三重項消滅が発生する可能性が高まる。

$$
|0110\rangle \rightarrow |0020\rangle, |0200\rangle, \ldots
$$

TTA反応の特性時間は：

$$
\tau_{\text{TTA}} \sim \frac{\hbar}{J_{\text{TTA}}} \sim \frac{0.66 \text{ eV} \cdot \text{fs}}{0.05 \text{ eV}} \sim 13 \text{ fs}
$$

### 11.3 後期過程（$t >$ 50 fs）

生成された励起一重項 $|S_1\rangle$ は、周囲との相互作用により緩和する。ただし、本モデルでは一重項の直接的な緩和経路を含まないため、一重項個体数は準安定的に維持される。

複雑な量子干渉効果により、個体数は周期的または準周期的な振動を示す可能性がある。

### 11.4 対称性の破れ

初期状態 $|1001\rangle$ は左右対称であるが、量子力学的な重ね合わせと干渉により、時間発展は非対称なパターンを示すことがある。

## 12. まとめ

### 12.1 理論的枠組み

本文書では、4分子直線配置モデルにおける分子三重項状態の量子ダイナミクスを記述する完全な理論的枠組みを構築した。

主要な要素：

1. **ヒルベルト空間**: $3^4 = 81$ 次元Qutrit積空間
2. **ハミルトニアン**: オンサイト項、TET項、TTA項の和
3. **初期状態**: 両端励起状態 $|1001\rangle$
4. **時間発展**: 2次Suzuki-Trotter分解
5. **量子ゲート実装**: MQT-Qudits疎構造認識コンパイラ
6. **厳密解**: 固有値対角化による解析解

### 12.2 数学的厳密性

すべての定式化は数学的に厳密であり、近似を用いる場合（鈴木トロッター分解）は誤差評価を明示している。

### 12.3 計算可能性

本理論は完全に実装可能であり、MQT-Quditsフレームワークを用いた数値シミュレーションにより検証可能である。

### 12.4 物理的洞察

両端励起初期状態は、エネルギー移動とTTA過程の競合・協調を観察するための理想的な設定である。この系の詳細な解析により、分子集合体における励起ダイナミクスの本質的理解が深まることが期待される。

---

## 参考文献

1. M. Suzuki, "Fractal decomposition of exponential operators with applications to many-body theories and Monte Carlo simulations," Physics Letters A, vol. 146, pp. 319-323, 1990.

2. H. F. Trotter, "On the product of semi-groups of operators," Proceedings of the American Mathematical Society, vol. 10, pp. 545-551, 1959.

3. D. L. Dexter, "A Theory of Sensitized Luminescence in Solids," The Journal of Chemical Physics, vol. 21, pp. 836-850, 1953.

4. N. J. Turro, V. Ramamurthy, and J. C. Scaiano, "Modern Molecular Photochemistry of Organic Molecules," University Science Books, 2010.

5. S. M. Menke and R. J. Holmes, "Exciton diffusion in organic photovoltaic cells," Energy & Environmental Science, vol. 7, pp. 499-512, 2014.

6. MQT-Qudits Documentation: https://github.com/cda-tum/mqt-qudits

---

**作成日**: 2025年11月10日
**バージョン**: 1.0
**著者**: MQT-Qudits研究グループ
