# 4分子直線配置モデルにおける分子三重項状態の量子ダイナミクス完全理論：Qubit実装

## 1. はじめに

### 1.1 本文書の目的

本文書では、4つの分子が一次元直線上に配列した系における励起三重項状態（triplet state）の量子ダイナミクスを、Qubitベースの量子アルゴリズムを用いて記述するための完全な理論的基礎を提供する。特に、両端の分子が励起三重項状態、中間の分子が基底状態にある初期状態からの時間発展を詳細に定式化する。

### 1.2 物理系の概要

考察する系は以下の特徴を持つ：

- **分子数**: $N = 4$
- **空間配置**: 一次元直線配列、最近接相互作用
- **各分子の電子状態**: 3つの状態（基底一重項、励起三重項、励起一重項）
- **エネルギー移動過程**: 三重項間エネルギー移動（Triplet Energy Transfer）、三重項-三重項消滅（Triplet-Triplet Annihilation）

### 1.3 Qubitエンコーディングの必要性

各分子が3つの電子状態を持つため、これを量子コンピュータ上で表現するには工夫が必要である。理想的には3準位系であるQutrit（Qudit, $d=3$）を用いるべきであるが、現在の量子コンピュータの多くはQubit（2準位系）ベースである。したがって、本文書では**2個のQubitで1つの分子状態をエンコード**する手法を採用する。

4分子系は合計8個のQubitで表現され、全ヒルベルト空間の次元は $2^8 = 256$ となる。ただし、物理的に意味のある状態空間は $3^4 = 81$ 次元の部分空間に制限される。

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

## 3. Qubit表現とヒルベルト空間

### 3.1 2-Qubit エンコーディング方式

各分子の3つの電子状態を2個のQubitでエンコードする。分子 $i$ に対して、Qubit $2i$ と Qubit $2i+1$ の組を用いる。

#### エンコーディング規則（Big-endian表記）

$$
\begin{align}
|S_0\rangle_i &\longleftrightarrow |00\rangle_{2i+1,2i} = |0\rangle_{2i+1} \otimes |0\rangle_{2i} \\
|T_1\rangle_i &\longleftrightarrow |01\rangle_{2i+1,2i} = |0\rangle_{2i+1} \otimes |1\rangle_{2i} \\
|S_1\rangle_i &\longleftrightarrow |10\rangle_{2i+1,2i} = |1\rangle_{2i+1} \otimes |0\rangle_{2i}
\end{align}
$$

ここで、$|00\rangle_{2i+1,2i}$ は左が Qubit $2i+1$、右が Qubit $2i$ を表す（Big-endian）。

#### Qiskitのlittle-endian規約との関係

Qiskitは**little-endian**規約を採用しているため、ビット順序が逆になる点に注意が必要である：

- Qiskit little-endian: qubit indexが小さいほど右側（下位ビット）
- 計算基底状態のインデックスは little-endian で解釈される

具体例：

- Big-endian $|01\rangle$ = Qubit $2i+1$ が $|0\rangle$、Qubit $2i$ が $|1\rangle$
- Little-endian インデックス: Qubit $2i$ が右側なので $|10\rangle$ に対応
- 状態準備: Qubit $2i$ に $X$ ゲートを適用

### 3.2 4分子系の全ヒルベルト空間

4つの分子からなる系の全ヒルベルト空間は、8個のQubitのテンソル積として構成される：

$$
\mathcal{H}_{\text{total}} = \bigotimes_{q=0}^{7} \mathcal{H}_q
$$

各 $\mathcal{H}_q$ は2次元ヒルベルト空間（Qubit）であり：

$$
\mathcal{H}_q = \text{span}\{|0\rangle_q, |1\rangle_q\}
$$

全系の次元は：

$$
\dim(\mathcal{H}_{\text{total}}) = 2^8 = 256
$$

### 3.3 物理的部分空間

全256次元のヒルベルト空間のうち、物理的に許される状態は各分子が $\{|S_0\rangle, |T_1\rangle, |S_1\rangle\}$ のいずれかにある81次元の部分空間に限定される。

禁止された状態の例：

- $|11\rangle_{2i+1,2i}$：対応する分子状態が存在しない

### 3.4 基底状態の表記

全系の物理的基底状態は：

$$
|n_0 n_1 n_2 n_3\rangle = |n_0\rangle_{\text{mol}} \otimes |n_1\rangle_{\text{mol}} \otimes |n_2\rangle_{\text{mol}} \otimes |n_3\rangle_{\text{mol}}
$$

ここで、$n_i \in \{S_0, T_1, S_1\}$ は分子 $i$ の状態を表す。

Qubit表現では：

$$
|n_0 n_1 n_2 n_3\rangle_{\text{mol}} \longleftrightarrow |q_7 q_6 q_5 q_4 q_3 q_2 q_1 q_0\rangle_{\text{qubit}}
$$

例えば：

- $|S_0 S_0 S_0 S_0\rangle$: $|00\,00\,00\,00\rangle$ = 全Qubitが $|0\rangle$
- $|T_1 T_1 T_1 T_1\rangle$: $|01\,01\,01\,01\rangle$ = Qubits 0,2,4,6 が $|1\rangle$
- $|T_1 S_0 S_0 T_1\rangle$: $|01\,00\,00\,01\rangle$ = Qubits 0,6 が $|1\rangle$（両端励起）

### 3.5 基底の完全性

81個の物理的基底状態は物理的部分空間で完全系を形成する：

$$
\sum_{n_0, n_1, n_2, n_3 \in \{S_0, T_1, S_1\}} |n_0 n_1 n_2 n_3\rangle \langle n_0 n_1 n_2 n_3| = \mathbb{P}_{\text{phys}}
$$

ここで $\mathbb{P}_{\text{phys}}$ は物理的部分空間への射影演算子である。

任意の物理的状態ベクトル $|\Psi\rangle$ は：

$$
|\Psi\rangle = \sum_{n_0, n_1, n_2, n_3} c_{n_0 n_1 n_2 n_3} |n_0 n_1 n_2 n_3\rangle
$$

と展開される。規格化条件は：

$$
\sum_{n_0, n_1, n_2, n_3} |c_{n_0 n_1 n_2 n_3}|^2 = 1
$$

## 4. ハミルトニアンの定式化

### 4.1 全ハミルトニアン（分子表現）

系の全ハミルトニアンは3つの項の和として表される：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}}
$$

各項は物理的に異なるプロセスを記述する。

### 4.2 オンサイトエネルギー項 $\hat{H}_0$

#### 4.2.1 物理的意味

$\hat{H}_0$ は各分子の固有エネルギーを記述する対角項である。

#### 4.2.2 分子表現

$$
\hat{H}_0 = \sum_{i=0}^{3} \left( E_{T_1} |T_1\rangle_i \langle T_1|_i + E_{S_1} |S_1\rangle_i \langle S_1|_i \right)
$$

基底状態 $|S_0\rangle$ のエネルギーは $E_{S_0} = 0$ と設定しているため、対応する項は省略されている。

#### 4.2.3 Qubit表現への変換

分子演算子をQubit演算子に変換する。分子 $i$ に対して：

$$
|T_1\rangle_i \langle T_1|_i \longleftrightarrow |01\rangle_{2i+1,2i} \langle 01|_{2i+1,2i}
$$

これはQubit演算子で表現すると：

$$
|01\rangle \langle 01| = |0\rangle \langle 0| \otimes |1\rangle \langle 1| = \frac{1}{2}(I + Z) \otimes \frac{1}{2}(I - Z)
$$

ここで、

$$
\begin{align}
|0\rangle \langle 0| &= \frac{1}{2}(I + Z) \\
|1\rangle \langle 1| &= \frac{1}{2}(I - Z)
\end{align}
$$

を用いた。したがって：

$$
|T_1\rangle_i \langle T_1|_i = \frac{1}{4}(I_{2i+1} + Z_{2i+1})(I_{2i} - Z_{2i})
$$

同様に：

$$
|S_1\rangle_i \langle S_1|_i \longleftrightarrow |10\rangle_{2i+1,2i} \langle 10|_{2i+1,2i} = \frac{1}{4}(I_{2i+1} - Z_{2i+1})(I_{2i} + Z_{2i})
$$

#### 4.2.4 Qubit ハミルトニアン

分子 $i$ のオンサイトハミルトニアンは：

$$
\hat{h}_0^{(i)} = E_{T_1} \cdot \frac{1}{4}(I + Z_{2i+1})(I - Z_{2i}) + E_{S_1} \cdot \frac{1}{4}(I - Z_{2i+1})(I + Z_{2i})
$$

展開すると：

$$
\hat{h}_0^{(i)} = \frac{E_{T_1}}{4}(I - Z_{2i} + Z_{2i+1} - Z_{2i+1}Z_{2i}) + \frac{E_{S_1}}{4}(I + Z_{2i} - Z_{2i+1} - Z_{2i+1}Z_{2i})
$$

整理すると：

$$
\hat{h}_0^{(i)} = \frac{E_{T_1} + E_{S_1}}{4}I + \frac{E_{S_1} - E_{T_1}}{4}Z_{2i} + \frac{E_{T_1} - E_{S_1}}{4}Z_{2i+1} - \frac{E_{T_1} + E_{S_1}}{4}Z_{2i+1}Z_{2i}
$$

全系のオンサイトハミルトニアンは：

$$
\hat{H}_0 = \sum_{i=0}^{3} \hat{h}_0^{(i)}
$$

#### 4.2.5 時間発展演算子

オンサイト項の時間発展演算子は各分子で分離できる：

$$
\exp(-i\hat{H}_0 t/\hbar) = \prod_{i=0}^{3} \exp(-i\hat{h}_0^{(i)} t/\hbar)
$$

各分子の時間発展は、パウリ $Z$ ゲートと2-qubit $ZZ$ ゲートの積で実装できる：

$$
\exp(-i\hat{h}_0^{(i)} t/\hbar) = \exp(-i\alpha_0 I) \exp(-i\alpha_1 Z_{2i}) \exp(-i\alpha_2 Z_{2i+1}) \exp(-i\alpha_3 Z_{2i+1}Z_{2i})
$$

ここで：

$$
\begin{align}
\alpha_0 &= \frac{(E_{T_1} + E_{S_1})t}{4\hbar} \\
\alpha_1 &= \frac{(E_{S_1} - E_{T_1})t}{4\hbar} \\
\alpha_2 &= \frac{(E_{T_1} - E_{S_1})t}{4\hbar} \\
\alpha_3 &= -\frac{(E_{T_1} + E_{S_1})t}{4\hbar}
\end{align}
$$

### 4.3 三重項エネルギー移動項 $\hat{H}_{\text{transfer}}$

#### 4.3.1 物理的過程

三重項エネルギー移動（Triplet Energy Transfer, TET）は、隣接する分子間で励起エネルギーが移動する過程である：

$$
T_1 + S_0 \leftrightarrow S_0 + T_1
$$

これは**デクスターメカニズム**（Dexter mechanism）により電子交換を通じて起こる。

#### 4.3.2 分子表現

最近接分子 $(i, j)$ 間の相互作用：

$$
\hat{H}_{\text{transfer}}^{(i,j)} = V \left( |S_0\rangle_i \langle T_1|_i \otimes |T_1\rangle_j \langle S_0|_j + |T_1\rangle_i \langle S_0|_i \otimes |S_0\rangle_j \langle T_1|_j \right)
$$

ここで、$V$ は移動積分（transfer integral）である。

全系では：

$$
\hat{H}_{\text{transfer}} = \sum_{(i,j) \in \text{neighbors}} \hat{H}_{\text{transfer}}^{(i,j)}
$$

4分子直線配置では、隣接ペアは $(0,1), (1,2), (2,3)$ の3組。

#### 4.3.3 Qubit表現への変換

遷移演算子のQubit表現：

$$
|S_0\rangle_i \langle T_1|_i \longleftrightarrow |00\rangle \langle 01| = |0\rangle \langle 0| \otimes |0\rangle \langle 1|
$$

Qubit演算子では：

$$
|0\rangle \langle 1| = \frac{1}{2}(X - iY)
$$

したがって：

$$
|S_0\rangle_i \langle T_1|_i = \frac{1}{2}(I + Z_{2i+1}) \otimes \frac{1}{2}(X_{2i} - iY_{2i})
$$

同様に：

$$
|T_1\rangle_i \langle S_0|_i = \frac{1}{2}(I + Z_{2i+1}) \otimes \frac{1}{2}(X_{2i} + iY_{2i})
$$

#### 4.3.4 2分子間相互作用のQubit表現

分子 $i$ と $j$ の間の相互作用：

$$
\hat{H}_{\text{transfer}}^{(i,j)} = \frac{V}{4} \left[ (I + Z_{2i+1})(X_{2i} - iY_{2i}) \otimes (I + Z_{2j+1})(X_{2j} + iY_{2j}) + \text{h.c.} \right]
$$

これは4-qubitゲートで実装される複雑な演算子であるが、物理的部分空間では適切に作用する。

### 4.4 三重項-三重項消滅項 $\hat{H}_{\text{TTA}}$

#### 4.4.1 物理的過程

三重項-三重項消滅（Triplet-Triplet Annihilation, TTA）は、隣接する2つの三重項励起が相互作用して、1つが一重項励起に、もう1つが基底状態に変換される過程である：

$$
T_1 + T_1 \rightarrow S_1 + S_0
$$

#### 4.4.2 エネルギー保存

TTA過程はエネルギー保存則 $2E_{T_1} = E_{S_1}$ を満たす場合にのみ起こる。本系では $E_{T_1} = 1.5$ eV、$E_{S_1} = 3.0$ eV なので条件を満たす。

#### 4.4.3 分子表現

最近接分子 $(i, j)$ 間の相互作用：

$$
\hat{H}_{\text{TTA}}^{(i,j)} = J \left( |S_1\rangle_i \langle T_1|_i \otimes |S_0\rangle_j \langle T_1|_j + |S_0\rangle_i \langle T_1|_i \otimes |S_1\rangle_j \langle T_1|_j + \text{h.c.} \right)
$$

ここで、$J$ はTTA相互作用強度である。

#### 4.4.4 Qubit表現への変換

遷移演算子：

$$
|S_1\rangle_i \langle T_1|_i \longleftrightarrow |10\rangle \langle 01| = |1\rangle \langle 0| \otimes |0\rangle \langle 1|
$$

Qubit演算子では：

$$
\begin{align}
|1\rangle \langle 0| &= \frac{1}{2}(X + iY) \\
|0\rangle \langle 1| &= \frac{1}{2}(X - iY)
\end{align}
$$

したがって：

$$
|S_1\rangle_i \langle T_1|_i = \frac{1}{4}(X_{2i+1} + iY_{2i+1})(X_{2i} - iY_{2i})
$$

同様に：

$$
|S_0\rangle_j \langle T_1|_j = \frac{1}{4}(I + Z_{2j+1})(X_{2j} - iY_{2j})
$$

#### 4.4.5 2分子間TTA相互作用のQubit表現

分子 $i$ と $j$ の間のTTA相互作用は4-qubitゲートで表現される：

$$
\hat{H}_{\text{TTA}}^{(i,j)} = \frac{J}{16} \left[ (X_{2i+1} + iY_{2i+1})(X_{2i} - iY_{2i}) \otimes (I + Z_{2j+1})(X_{2j} - iY_{2j}) + \cdots \right]
$$

これは複数のパウリ演算子の積で構成される複雑な演算子である。

### 4.5 ハミルトニアンの性質

#### 4.5.1 エルミート性

全ハミルトニアン $\hat{H}_{\text{total}}$ はエルミート演算子である：

$$
\hat{H}_{\text{total}}^\dagger = \hat{H}_{\text{total}}
$$

これは各項がエルミートであることから保証される。

#### 4.5.2 実対称性

適切なゲージ選択により、ハミルトニアンは実対称行列として表現できる。

#### 4.5.3 疎行列構造

ハミルトニアンは最近接相互作用のみを含むため、$256 \times 256$ 行列として表現したとき疎行列となる。

## 5. 初期状態の設定

### 5.1 両端励起初期状態（Edge Triplet State）

本シミュレーションでは、以下の初期状態を考察する：

$$
|\Psi(0)\rangle = |T_1 S_0 S_0 T_1\rangle
$$

これは、分子0と分子3（両端）が励起三重項状態、分子1と分子2（中間）が基底状態にある状態である。

### 5.2 Qubit表現

この初期状態のQubit表現は：

$$
|\Psi(0)\rangle = |01\,00\,00\,01\rangle_{7,6,5,4,3,2,1,0}
$$

Big-endian表記では、各分子ペアは：

- 分子0: $|01\rangle$ → Qubits 1,0
- 分子1: $|00\rangle$ → Qubits 3,2
- 分子2: $|00\rangle$ → Qubits 5,4
- 分子3: $|01\rangle$ → Qubits 7,6

### 5.3 量子ゲートによる準備

初期状態は、全Qubitが $|0\rangle$ の状態から以下のゲート操作で準備される：

1. 全Qubitを $|0\rangle$ に初期化（デフォルト）
2. Qubit 0 に $X$ ゲートを適用（分子0を $|T_1\rangle$ に設定）
3. Qubit 6 に $X$ ゲートを適用（分子3を $|T_1\rangle$ に設定）

Qiskitのlittle-endian規約では：

- Qubit 0 が右端（最下位ビット）
- Qubit 6 も同様に対応する分子の右側qubit

したがって、状態準備コードは：

```python
circuit.x(0)  # 分子0をT1状態に
circuit.x(6)  # 分子3をT1状態に
```

### 5.4 初期状態の物理的意味

両端励起初期状態は以下の物理的状況を模擬する：

1. **光励起**: 系の両端に位置する分子が光吸収により励起される
2. **エネルギー移動の開始**: 両端から中央に向かってエネルギーが伝播する
3. **対称性**: 空間反転対称性を持つ初期配置

この初期状態から、以下のダイナミクスが期待される：

- **初期過程**: 三重項エネルギー移動により、励起が中央に拡散
- **中期過程**: 三重項が隣接すると TTA が発生し、一重項励起が生成
- **後期過程**: 一重項励起からの緩和過程

### 5.5 初期状態のエネルギー

初期状態のエネルギー期待値は：

$$
E(0) = \langle \Psi(0) | \hat{H}_0 | \Psi(0) \rangle = 2 E_{T_1} = 3.0 \text{ eV}
$$

これは2つの三重項励起のエネルギーに対応する。

### 5.6 初期個体数

初期状態における各状態の個体数は：

$$
\begin{align}
N_{S_0}(0) &= 2 \quad \text{(分子1と分子2)} \\
N_{T_1}(0) &= 2 \quad \text{(分子0と分子3)} \\
N_{S_1}(0) &= 0
\end{align}
$$

個体数保存則：

$$
N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = 4
$$

は時間発展中常に成立する。

## 6. 時間発展演算子とシュレーディンガー方程式

### 6.1 時間依存シュレーディンガー方程式

系の時間発展は時間依存シュレーディンガー方程式に従う：

$$
i\hbar \frac{\partial}{\partial t} |\Psi(t)\rangle = \hat{H}_{\text{total}} |\Psi(t)\rangle
$$

ハミルトニアンが時間に依存しない場合、形式解は：

$$
|\Psi(t)\rangle = \exp\left(-\frac{i}{\hbar}\hat{H}_{\text{total}} t\right) |\Psi(0)\rangle
$$

### 6.2 時間発展演算子

時間発展演算子は：

$$
\hat{U}(t) = \exp\left(-\frac{i}{\hbar}\hat{H}_{\text{total}} t\right)
$$

#### 6.2.1 ユニタリ性

$\hat{H}_{\text{total}}$ がエルミートであるため、$\hat{U}(t)$ はユニタリ演算子である：

$$
\hat{U}^\dagger(t) \hat{U}(t) = \hat{U}(t) \hat{U}^\dagger(t) = \mathbb{I}
$$

これは確率保存を保証する：

$$
\langle \Psi(t) | \Psi(t) \rangle = \langle \Psi(0) | \hat{U}^\dagger(t) \hat{U}(t) | \Psi(0) \rangle = \langle \Psi(0) | \Psi(0) \rangle = 1
$$

#### 6.2.2 群性質

時間発展演算子は連続群を形成する：

$$
\hat{U}(t_1 + t_2) = \hat{U}(t_1) \hat{U}(t_2)
$$

### 6.3 単位系とプランク定数

本シミュレーションでは以下の単位系を採用する：

- **エネルギー**: eV（電子ボルト）
- **時間**: fs（フェムト秒、$10^{-15}$ 秒）
- **換算プランク定数**: $\hbar = 0.6582$ eV·fs

この単位系では、時間発展の位相は：

$$
\phi = \frac{E \cdot t}{\hbar} = \frac{E \text{ [eV]} \times t \text{ [fs]}}{0.6582 \text{ [eV·fs]}}
$$

## 7. 鈴木トロッター分解による時間発展の数値計算

### 7.1 問題の本質

ハミルトニアンが複数の非可換項の和：

$$
\hat{H}_{\text{total}} = \hat{H}_0 + \hat{H}_{\text{transfer}} + \hat{H}_{\text{TTA}}
$$

である場合、時間発展演算子を直接計算することは困難である。なぜなら：

$$
\exp(-i(\hat{A} + \hat{B})t) \neq \exp(-i\hat{A}t) \exp(-i\hat{B}t) \quad \text{if } [\hat{A}, \hat{B}] \neq 0
$$

### 7.2 Baker-Campbell-Hausdorff公式

一般に、非可換演算子に対して：

$$
\exp(\hat{A}) \exp(\hat{B}) = \exp\left(\hat{A} + \hat{B} + \frac{1}{2}[\hat{A}, \hat{B}] + \frac{1}{12}[\hat{A}, [\hat{A}, \hat{B}]] + \cdots\right)
$$

この公式から、小さな時間刻み $\Delta t$ に対して近似が可能である。

### 7.3 1次Lie-Trotter分解

最も単純な近似は1次Lie-Trotter分解：

$$
\exp(-i\hat{H}_{\text{total}} \Delta t) \approx \exp(-i\hat{H}_0 \Delta t) \exp(-i\hat{H}_{\text{transfer}} \Delta t) \exp(-i\hat{H}_{\text{TTA}} \Delta t)
$$

誤差は $O(\Delta t^2)$ である。

### 7.4 2次Suzuki-Trotter分解（対称分解）

精度を向上させるため、2次の対称分解を用いる：

$$
\begin{align}
\exp(-i\hat{H}_{\text{total}} \Delta t) \approx &\exp(-i\hat{H}_0 \Delta t/2) \exp(-i\hat{H}_{\text{transfer}} \Delta t/2) \exp(-i\hat{H}_{\text{TTA}} \Delta t/2) \\
&\times \exp(-i\hat{H}_{\text{TTA}} \Delta t/2) \exp(-i\hat{H}_{\text{transfer}} \Delta t/2) \exp(-i\hat{H}_0 \Delta t/2)
\end{align}
$$

この分解の誤差は $O(\Delta t^3)$ であり、より高精度である。

### 7.5 具体的な実装手順

総時間 $T$ を $N$ ステップに分割し、時間刻み $\Delta t = T/N$ とする。

各ステップで以下の量子ゲート列を適用：

1. **前半のオンサイト項**: $\exp(-i\hat{H}_0 \Delta t/2)$

   - 各分子 $i = 0, 1, 2, 3$ に対して、$Z$ ゲートと $ZZ$ ゲートを適用

2. **前半の移動項**: $\exp(-i\hat{H}_{\text{transfer}} \Delta t/2)$

   - 各隣接ペア $(i, j)$ に対して、4-qubitゲートを適用

3. **前半のTTA項**: $\exp(-i\hat{H}_{\text{TTA}} \Delta t/2)$

   - 各隣接ペア $(i, j)$ に対して、4-qubitゲートを適用

4. **後半のTTA項**: $\exp(-i\hat{H}_{\text{TTA}} \Delta t/2)$（同じ）

5. **後半の移動項**: $\exp(-i\hat{H}_{\text{transfer}} \Delta t/2)$（同じ）

6. **後半のオンサイト項**: $\exp(-i\hat{H}_0 \Delta t/2)$（同じ）

### 7.6 誤差評価

トロッター分解の総誤差は：

$$
\left\| \hat{U}_{\text{exact}}(T) - \hat{U}_{\text{Trotter}}(T) \right\| = O(T^3 \Delta t^2 / N^2)
$$

したがって、ステップ数 $N$ を増やすことで精度を向上できる。

## 8. Qubit量子ゲート実装の詳細

### 8.1 量子ゲート分解の原理

各ハミルトニアン項 $\exp(-i\hat{H}_k t)$ を、基本的なQubitゲート（Pauli gates, CNOT, 回転ゲート）の列に分解する必要がある。

### 8.2 オンサイト項 $\exp(-i\hat{H}_0 t)$ の実装

分子 $i$ のオンサイト項は：

$$
\hat{h}_0^{(i)} = c_0 I + c_1 Z_{2i} + c_2 Z_{2i+1} + c_3 Z_{2i+1}Z_{2i}
$$

時間発展演算子は：

$$
\exp(-i\hat{h}_0^{(i)} t) = e^{-ic_0 t} \exp(-ic_1 Z_{2i} t) \exp(-ic_2 Z_{2i+1} t) \exp(-ic_3 Z_{2i+1}Z_{2i} t)
$$

（可換性により分離可能）

#### 実装手順

1. **単一qubit $Z$ 回転**: $\exp(-ic_1 Z_{2i} t) = R_Z(2c_1 t)$

   - Qiskitでは `rz(2*c1*t, qubit_2i)` で実装

2. **単一qubit $Z$ 回転**: $\exp(-ic_2 Z_{2i+1} t) = R_Z(2c_2 t)$

   - Qiskitでは `rz(2*c2*t, qubit_2i+1)` で実装

3. **2-qubit $ZZ$ 回転**: $\exp(-ic_3 Z_{2i+1}Z_{2i} t) = R_{ZZ}(2c_3 t)$

   - Qiskitでは `rzz(2*c3*t, qubit_2i+1, qubit_2i)` で実装

4. **グローバル位相**: $e^{-ic_0 t}$ は物理的観測量に影響しないため無視可能

### 8.3 2体相互作用項の実装

#### 8.3.1 一般的な2体演算子

移動項とTTA項は、2つの分子（4個のQubit）間の相互作用を記述する。

一般形：

$$
\hat{H}^{(i,j)} = \sum_{a,b,c,d} h_{abcd} \sigma_a^{(2i+1)} \sigma_b^{(2i)} \sigma_c^{(2j+1)} \sigma_d^{(2j)}
$$

ここで、$\sigma \in \{I, X, Y, Z\}$ はパウリ演算子。

#### 8.3.2 Qiskitでの実装

Qiskitでは、複雑な2体演算子は以下の方法で実装できる：

1. **行列表現**: $16 \times 16$ ユニタリ行列として計算
2. **カスタムゲート**: `UnitaryGate` を使用して回路に追加
3. **ゲート分解**: CNOTと単一qubitゲートに分解（KAK分解など）

本実装では、カスタムゲートとして直接適用する方法を採用する。

#### 8.3.3 $\hat{H}_{\text{transfer}}$ の実装

各隣接ペア $(i, j)$ に対して：

```python
# 行列表現を構築（物理的部分空間で）
H_transfer_matrix = build_transfer_matrix(i, j, V, dt)

# ユニタリ演算子
U_transfer = expm(-1j * H_transfer_matrix)

# カスタムゲートとして適用
circuit.unitary(U_transfer, [2 * i, 2 * i + 1, 2 * j, 2 * j + 1])
```

#### 8.3.4 $\hat{H}_{\text{TTA}}$ の実装

同様に、TTA項も4-qubitカスタムゲートとして実装：

```python
# 行列表現を構築
H_TTA_matrix = build_TTA_matrix(i, j, J, dt)

# ユニタリ演算子
U_TTA = expm(-1j * H_TTA_matrix)

# カスタムゲートとして適用
circuit.unitary(U_TTA, [2 * i, 2 * i + 1, 2 * j, 2 * j + 1])
```

### 8.4 物理的部分空間の維持

重要な点として、時間発展演算子は物理的部分空間を保存する必要がある。すなわち、初期状態が物理的であれば、時間発展後も物理的状態である必要がある。

これは、ハミルトニアンが物理的部分空間で閉じているため保証される。

## 9. 観測可能量の計算

### 9.1 個体数演算子

各電子状態の個体数演算子は、その状態にある分子の数をカウントする。

#### 分子表現

$$
\begin{align}
\hat{N}_{S_0} &= \sum_{i=0}^{3} |S_0\rangle_i \langle S_0|_i \\
\hat{N}_{T_1} &= \sum_{i=0}^{3} |T_1\rangle_i \langle T_1|_i \\
\hat{N}_{S_1} &= \sum_{i=0}^{3} |S_1\rangle_i \langle S_1|_i
\end{align}
$$

#### Qubit表現

各分子 $i$ に対して：

$$
\begin{align}
|S_0\rangle_i \langle S_0|_i &= |00\rangle \langle 00| = \frac{1}{4}(I + Z_{2i+1})(I + Z_{2i}) \\
|T_1\rangle_i \langle T_1|_i &= |01\rangle \langle 01| = \frac{1}{4}(I + Z_{2i+1})(I - Z_{2i}) \\
|S_1\rangle_i \langle S_1|_i &= |10\rangle \langle 10| = \frac{1}{4}(I - Z_{2i+1})(I + Z_{2i})
\end{align}
$$

### 9.2 状態ベクトルからの直接計算

量子回路シミュレーションでは、状態ベクトル $|\Psi(t)\rangle$ が得られる。

各時刻での個体数は：

$$
N_{S_0}(t) = \langle \Psi(t) | \hat{N}_{S_0} | \Psi(t) \rangle
$$

計算基底 $\{|q_7 \cdots q_0\rangle\}$ での展開：

$$
|\Psi(t)\rangle = \sum_{q_7, \ldots, q_0 = 0}^{1} c_{q_7 \cdots q_0} |q_7 \cdots q_0\rangle
$$

を用いて、各基底状態がどの分子状態に対応するかを判定し、個体数を計算する。

### 9.3 期待値の時間発展

個体数の時間発展をプロットすることで、以下の物理過程を可視化できる：

- **三重項個体数 $N_{T_1}(t)$**: 初期値2から減少
- **基底状態個体数 $N_{S_0}(t)$**: 初期値2から変化
- **一重項個体数 $N_{S_1}(t)$**: 初期値0から増加（TTA発生時）

### 9.4 保存則

個体数保存則：

$$
N_{S_0}(t) + N_{T_1}(t) + N_{S_1}(t) = 4
$$

は常に成立する。これは数値計算の正確性を検証する指標となる。

エネルギー保存則：

$$
E(t) = \langle \Psi(t) | \hat{H}_{\text{total}} | \Psi(t) \rangle = E(0) = 3.0 \text{ eV}
$$

も成立する（時間発展がユニタリであるため）。

## 10. 数値シミュレーションの詳細

### 10.1 Qiskitによる実装

本シミュレーションはQiskitフレームワークを用いて実装される。

#### 主要コンポーネント

1. **PhysicalParameters**: 物理パラメータの管理

   - 分子数、エネルギー、相互作用強度

2. **StateEncoder**: 初期状態の準備

   - `edge_triplet` モード：両端励起状態

3. **HamiltonianGates**: ハミルトニアン項のゲート実装

   - オンサイト項、移動項、TTA項の時間発展

4. **TrotterCircuitBuilder**: トロッター回路の構築

   - 2次対称分解の実装

5. **ObservableCalculator**: 観測量の計算

   - 個体数の計算

6. **QubitMolecularDynamicsSimulator**: メインシミュレータ
   - 時間発展の実行と結果の可視化

### 10.2 計算複雑度

#### メモリ使用量

8-qubit系の状態ベクトルは $2^8 = 256$ 個の複素振幅を持つ。

メモリ: $256 \times 16$ bytes (複素数) $= 4$ KB

これは非常に小さく、古典コンピュータで容易にシミュレート可能。

#### 計算時間

各トロッターステップは以下を含む：

- オンサイト項: 4分子 $\times$ (数個のゲート) $= O(10)$ ゲート
- 移動項: 3ペア $\times$ 1カスタムゲート $= 3$ カスタムゲート
- TTA項: 3ペア $\times$ 1カスタムゲート $= 3$ カスタムゲート

総ゲート数: $O(20)$ per step

20ステップの場合: $O(400)$ ゲート

状態ベクトルシミュレーションの計算量: $O(N_{\text{gates}} \times 2^8) = O(10^5)$ 演算

これは数秒以内で実行可能。

### 10.3 実装上の工夫

#### 物理的部分空間での演算

全256次元空間ではなく、物理的81次元部分空間のみで演算を行うことで効率化が可能（本実装では簡単のため全空間を使用）。

#### 疎行列利用

ハミルトニアンの疎性を利用した効率的な行列演算。

#### ゲート最適化

同じ種類のゲートをまとめて適用することで、回路深度を削減。

### 10.4 検証方法

#### 保存則のチェック

- 個体数保存: $N_{S_0} + N_{T_1} + N_{S_1} = 4$
- エネルギー保存（近似的）

#### トロッター誤差の評価

時間刻みを変えてシミュレーションを実行し、収束性を確認。

#### 物理的妥当性

- 初期過程でのエネルギー移動の方向
- TTA発生のタイミング
- 長時間での緩和挙動

## 11. 予想される物理的挙動

### 11.1 初期過程（$t \sim 0$-20 fs）

**三重項エネルギー移動の開始**

両端の励起三重項（分子0と分子3）から、隣接する基底状態分子（分子1と分子2）へエネルギーが移動する。

期待される変化：

- $N_{T_1}$: わずかに減少（または維持）
- $N_{S_0}$: わずかに変化
- $N_{S_1}$: ほぼ0（TTA未発生）

### 11.2 中期過程（$t \sim$ 20-50 fs）

**三重項-三重項消滅（TTA）の発生**

エネルギー移動により、隣接する分子が共に三重項状態になると、TTA過程が発生：

$$
T_1 + T_1 \rightarrow S_1 + S_0
$$

期待される変化：

- $N_{T_1}$: 急速に減少（2ずつ減少）
- $N_{S_0}$: 増加（1ずつ増加）
- $N_{S_1}$: 増加（1ずつ増加）

### 11.3 後期過程（$t >$ 50 fs）

**一重項励起からの緩和**

生成された一重項励起 $S_1$ は、蛍光放出や内部変換により基底状態に緩和する可能性がある（本モデルでは緩和項を含まない場合、振動する）。

期待される変化：

- 個体数の振動または準定常状態への到達

### 11.4 対称性

初期状態が空間反転対称性を持つため、ダイナミクスも対称性を維持する（近似的に）。

分子0と分子3、分子1と分子2の個体数分布は鏡像関係を保つ。

### 11.5 エネルギー移動の方向性

初期状態 $|T_1 S_0 S_0 T_1\rangle$ では、両端から中央に向かってエネルギーが伝播する。

中央（分子1と分子2の間）でエネルギーが集中し、TTA が優先的に発生する可能性がある。

## 12. Qubit実装とQudit実装の比較

### 12.1 表現の違い

| 項目                | Qudit (Qutrit) 実装 | Qubit 実装         |
| ------------------- | ------------------- | ------------------ |
| 1分子あたりの量子系 | 1 Qutrit (3準位)    | 2 Qubits           |
| 4分子系の量子系数   | 4 Qutrits           | 8 Qubits           |
| ヒルベルト空間次元  | $3^4 = 81$          | $2^8 = 256$        |
| 物理的状態空間      | 81次元（全空間）    | 81次元（部分空間） |

### 12.2 ゲート実装の違い

| 項目         | Qudit 実装              | Qubit 実装             |
| ------------ | ----------------------- | ---------------------- |
| オンサイト項 | 単一Qutrit回転          | 複数のPauliゲート      |
| 移動項       | 2-Qutrit カスタムゲート | 4-Qubit カスタムゲート |
| TTA項        | 2-Qutrit カスタムゲート | 4-Qubit カスタムゲート |

### 12.3 計算資源の比較

#### メモリ使用量

- Qudit: $81 \times 16$ bytes $= 1.3$ KB
- Qubit: $256 \times 16$ bytes $= 4$ KB

Qubit実装は約3倍のメモリを使用。

#### 計算時間

状態空間サイズの違いにより、Qudit実装の方が効率的（理論上）。

ただし、現実の量子ハードウェアはQubitベースであるため、Qubit実装の方が実機での実行に適している。

### 12.4 精度とエラー

#### 物理的部分空間からの逸脱

Qubit実装では、数値誤差により物理的部分空間から逸脱する可能性がある。

対策：

- 高精度演算の使用
- 各ステップでの部分空間への射影（必要に応じて）

#### トロッター誤差

両実装で同様のトロッター誤差が発生。時間刻みを小さくすることで改善。

### 12.5 実用的側面

#### 現在の量子コンピュータ

- IBM Quantum, Google Quantumなど：Qubitベース
- Qudit量子コンピュータ：研究段階

したがって、Qubit実装は現行ハードウェアで実行可能であり、実用的価値が高い。

#### スケーラビリティ

大規模系（多分子系）では：

- Qudit: $3^N$ 次元
- Qubit: $2^{2N}$ 次元

Qubit実装は指数的に資源を消費するが、量子もつれを利用した量子アルゴリズムにより効率化が期待される。

## 13. まとめ

### 13.1 理論的枠組み

本文書では、4分子直線配置モデルにおける分子三重項状態の量子ダイナミクスを、Qubitベースの量子アルゴリズムを用いて厳密に定式化した。

主要な成果：

1. **Qubitエンコーディング**: 3準位系（分子状態）を2-Qubitで表現する手法の確立
2. **ハミルトニアンの変換**: 分子演算子からQubit演算子（Pauli演算子）への変換公式
3. **初期状態の設定**: 両端励起状態（edge triplet）のQubit表現と準備方法
4. **時間発展の実装**: 鈴木トロッター分解による量子回路の構築

### 13.2 数学的厳密性

全ての数式は以下の点で厳密である：

- **ヒルベルト空間の定義**: テンソル積構造と次元の明示
- **ハミルトニアンの導出**: 物理的過程からの系統的な構築
- **演算子の変換**: 分子表現からQubit表現への明示的変換公式
- **誤差評価**: トロッター分解の誤差オーダーの明示

### 13.3 計算可能性

本理論は以下の点で計算可能である：

- **Qiskitによる実装**: 標準的な量子回路フレームワークでの実装が可能
- **古典シミュレーション**: 8-qubit系は古典コンピュータで容易にシミュレート可能
- **量子ハードウェア対応**: 現行のQubit量子コンピュータで実行可能（ノイズ対策は必要）

### 13.4 物理的洞察

両端励起初期状態からのダイナミクスは以下を示すと期待される：

- **エネルギー移動**: 両端から中央への三重項エネルギーの伝播
- **TTA発生**: 中央領域での三重項-三重項消滅と一重項生成
- **対称性**: 空間反転対称性の維持（近似的）
- **緩和過程**: 長時間での定常状態または振動状態への到達

### 13.5 今後の発展

本理論は以下の方向に拡張可能である：

1. **大規模系**: より多くの分子を含む系への拡張
2. **環境効果**: デコヒーレンスや散逸を含むモデル
3. **最適制御**: 外部場による励起状態の制御
4. **実験検証**: 実験データとの比較による理論の検証
5. **量子優位性**: 古典計算では困難な大規模系でのQubit量子アルゴリズムの優位性実証

### 13.6 結語

本文書は、分子三重項状態の量子ダイナミクスをQubitベースの量子コンピュータでシミュレートするための完全な理論的基盤を提供した。数学的厳密性と計算可能性を両立し、物理的洞察を与える本理論は、量子化学と量子情報科学の融合領域における重要な一歩である。

## 参考文献

1. Smith, M. B., & Michl, J. (2010). "Singlet fission." _Chemical Reviews_, 110(11), 6891-6936.

2. Casanova, D. (2018). "Theoretical modeling of singlet fission." _Chemical Reviews_, 118(15), 7164-7207.

3. Nielsen, M. A., & Chuang, I. L. (2010). _Quantum Computation and Quantum Information_. Cambridge University Press.

4. Suzuki, M. (1976). "Generalized Trotter's formula and systematic approximants of exponential operators and inner derivations with applications to many-body problems." _Communications in Mathematical Physics_, 51(2), 183-190.

5. Qiskit Development Team. (2021). _Qiskit: An Open-source Framework for Quantum Computing_. https://qiskit.org/

6. McClean, J. R., et al. (2016). "The theory of variational hybrid quantum-classical algorithms." _New Journal of Physics_, 18(2), 023023.

7. Dexter, D. L. (1953). "A theory of sensitized luminescence in solids." _The Journal of Chemical Physics_, 21(5), 836-850.

8. Pope, M., & Swenberg, C. E. (1999). _Electronic Processes in Organic Crystals and Polymers_. Oxford University Press.
