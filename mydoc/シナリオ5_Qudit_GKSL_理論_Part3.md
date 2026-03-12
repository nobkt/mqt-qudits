# シナリオ5: Qudit GKSL（ボソン無し）の省略無し詳細理論 — Part 3

## Trotter分解・シミュレーションアルゴリズム・誤差解析・ゲート数見積

本文書は Part 1（物理系・GKSL方程式）、Part 2（Stinespring拡張）に続き、シナリオ5の完全なシミュレーションアルゴリズムを省略無しに定式化する。Hamiltonian-Dissipator分割のStrang分割法、Lindbladチャネルの回文積順序、各分割の誤差解析、ゲート数見積、収束特性を扱う。

---

## 1. GKSL生成子の分割

### 1.1 Hamiltonian-Dissipator 分割

GKSL生成子 $\mathcal{L}$ を以下の2つの成分に分割する：

$$
\mathcal{L} = \mathcal{L}_H + \mathcal{L}_D
$$

ここで：

- $\mathcal{L}_H[\rho] = -i[H, \rho]$（ハミルトニアン部分：ユニタリ時間発展）
- $\mathcal{L}_D[\rho] = \sum_{\alpha=1}^{26} \mathcal{D}[\tilde{L}_\alpha](\rho)$（散逸部分：非ユニタリ過程）

### 1.2 散逸部分のさらなる分割

散逸生成子は26個の個別チャネル生成子に分割される：

$$
\mathcal{L}_D = \sum_{\alpha=1}^{26} \mathcal{L}_{D_\alpha}, \qquad \mathcal{L}_{D_\alpha}[\rho] = \mathcal{D}[\tilde{L}_\alpha](\rho)
$$

---

## 2. Strang分割（対称Trotter分割）

### 2.1 2次対称分割の定義

Hamiltonian部分 $\mathcal{L}_H$ と散逸部分 $\mathcal{L}_D$ の間の分割に、2次対称（Strang）分割を適用する。1タイムステップ $\Delta t$ の時間発展は以下の形式で近似される：

$$
e^{\mathcal{L}\Delta t} \approx e^{\mathcal{L}_H \Delta t/2}\, e^{\mathcal{L}_D \Delta t}\, e^{\mathcal{L}_H \Delta t/2}
$$

### 2.2 Strang分割の誤差

対称Trotter分割の誤差は：

$$
e^{\mathcal{L}_H \Delta t/2}\, e^{\mathcal{L}_D \Delta t}\, e^{\mathcal{L}_H \Delta t/2} = e^{\mathcal{L}\Delta t} + O(\Delta t^3)
$$

すなわち、1ステップあたりの局所誤差は $O(\Delta t^3)$ であり、$n_{\mathrm{steps}} = t_{\mathrm{max}}/\Delta t$ ステップの積算で全体（大域）誤差は $O(\Delta t^2)$ となる。

**Baker-Campbell-Hausdorff公式による導出**: 1次Trotter分割 $e^{A\tau}e^{B\tau} = e^{(A+B)\tau + \frac{\tau^2}{2}[A,B] + \cdots}$ に対して、対称分割は：

$$
e^{A\tau/2}e^{B\tau}e^{A\tau/2} = e^{(A+B)\tau + \frac{\tau^3}{24}([A,[A,B]] + 2[B,[B,A]]) + \cdots}
$$

であり、$\tau^2$ の項（交換子誤差）がキャンセルする。

### 2.3 ただし: 散逸部分内部の分割が必要

$e^{\mathcal{L}_D \Delta t}$ を直接計算することはできない（量子回路に対応しないため）。代わりに、$\mathcal{L}_D = \sum_{\alpha} \mathcal{L}_{D_\alpha}$ をさらに個別チャネルの積に分割する必要がある。

---

## 3. 回文積（Palindromic Product）によるLindblad チャネルの分割

### 3.1 Lie-Trotter積公式

$e^{\mathcal{L}_D \Delta t}$ を個別チャネルの積で近似する最も単純な方法は1次Lie-Trotter積：

$$
e^{\mathcal{L}_D \Delta t} \approx \prod_{\alpha=1}^{26} e^{\mathcal{L}_{D_\alpha} \Delta t}
$$

しかし、この近似の誤差は $O(\Delta t^2)$（1ステップ）であり、個別チャネル生成子間の交換子 $[\mathcal{L}_{D_\alpha}, \mathcal{L}_{D_\beta}]$ に起因する。

### 3.2 回文積（Palindromic Product）の定義

回文積は個別チャネルを前進・逆順の対称的な順序で適用する。$n_L = 26$ 個のチャネルに対して：

$$
e^{\mathcal{L}_D \Delta t} \approx \left(\prod_{\alpha=1}^{n_L} e^{\mathcal{L}_{D_\alpha} \Delta t/2}\right) \left(\prod_{\alpha=n_L}^{1} e^{\mathcal{L}_{D_\alpha} \Delta t/2}\right)
$$

すなわち、まず $\alpha = 1, 2, \ldots, 26$ の順で半ステップ $\Delta t/2$ を適用し、次に $\alpha = 26, 25, \ldots, 1$ の逆順で半ステップ $\Delta t/2$ を適用する。

### 3.3 回文積の誤差

回文積は2次対称Trotter分割の一般化であり、1次の交換子誤差がキャンセルする：

$$
\left(\prod_{\alpha=1}^{n_L} e^{\mathcal{L}_{D_\alpha} \Delta t/2}\right) \left(\prod_{\alpha=n_L}^{1} e^{\mathcal{L}_{D_\alpha} \Delta t/2}\right) = e^{\mathcal{L}_D \Delta t} + O(\Delta t^3)
$$

1ステップあたりの局所誤差は $O(\Delta t^3)$、大域誤差は $O(\Delta t^2)$。

**しかし**、シナリオ5では $e^{\mathcal{L}_{D_\alpha} \Delta t/2}$ の代わりにStinespring近似 $\mathcal{E}_\alpha(\Delta t/2)$ を使用する。Part 2 で示したように、Stinespring近似は $\mathcal{E}_\alpha(\Delta t/2) = e^{\mathcal{L}_{D_\alpha} \Delta t/2} + O(\Delta t^2)$（半ステップの場合は $O((\Delta t/2)^2) = O(\Delta t^2)$）であるため、回文積の2次精度はStinespring近似の1次精度に律速される。

---

## 4. 完全なTrotterステップの構成

### 4.1 1ステップの全体構造

ハミルトニアンのStrang分割と散逸部分の回文積を組み合わせた、シナリオ5の1タイムステップ $\Delta t$ の完全な構成は以下の通りである：

$$
\rho(t + \Delta t) \approx \mathcal{T}_{\Delta t}[\rho(t)]
$$

ここで $\mathcal{T}_{\Delta t}$ は以下の合成写像である（右から左に適用）：

$$
\mathcal{T}_{\Delta t} = \underbrace{\mathcal{U}_H^{\Delta t/2}}_{\text{(4)}} \circ \underbrace{\mathcal{E}_1^{\Delta t/2} \circ \mathcal{E}_2^{\Delta t/2} \circ \cdots \circ \mathcal{E}_{26}^{\Delta t/2}}_{\text{(3) 逆順}} \circ \underbrace{\mathcal{E}_{26}^{\Delta t/2} \circ \cdots \circ \mathcal{E}_2^{\Delta t/2} \circ \mathcal{E}_1^{\Delta t/2}}_{\text{(2) 順方向}} \circ \underbrace{\mathcal{U}_H^{\Delta t/2}}_{\text{(1)}}
$$

ここで：
- $\mathcal{U}_H^{\Delta t/2}[\rho] = U_H^{(1/2)}\, \rho\, U_H^{(1/2)\dagger}$、$U_H^{(1/2)} = e^{-iH\Delta t/2}$（$81 \times 81$ ユニタリ）
- $\mathcal{E}_\alpha^{\Delta t/2}$: 第 $\alpha$ Lindbladチャネルの Stinespring 近似（半ステップ $\Delta t/2$）

### 4.2 適用手順（左から右の計算順序）

$\mathcal{T}_{\Delta t}$ を $\rho$ に適用する計算手順は、以下の4段階に分かれる。**コードの実行順序に沿って（上から下に）記述する**：

**(1) 前半ハミルトニアン**:

$$
\rho^{(1)} = U_H^{(1/2)}\, \rho\, U_H^{(1/2)\dagger}
$$

**(2) 順方向Lindbladチャネル（$\alpha = 1, 2, \ldots, 26$）**:

$$
\rho^{(2,0)} = \rho^{(1)}
$$

$$
\rho^{(2,\alpha)} = \mathcal{E}_\alpha^{\Delta t/2}[\rho^{(2,\alpha-1)}], \qquad \alpha = 1, 2, \ldots, 26
$$

$$
\rho^{(2)} = \rho^{(2,26)}
$$

**(3) 逆順Lindbladチャネル（$\alpha = 26, 25, \ldots, 1$）**:

$$
\rho^{(3,0)} = \rho^{(2)}
$$

$$
\rho^{(3,k)} = \mathcal{E}_{27-k}^{\Delta t/2}[\rho^{(3,k-1)}], \qquad k = 1, 2, \ldots, 26
$$

$$
\rho^{(3)} = \rho^{(3,26)}
$$

**(4) 後半ハミルトニアン**:

$$
\rho(t + \Delta t) = U_H^{(1/2)}\, \rho^{(3)}\, U_H^{(1/2)\dagger}
$$

**コード対応** (`QuditGKSLSimulator._trotter_step`):

```python
def _trotter_step(self, rho):
    # (1) 前半ハミルトニアン
    rho = self._U_H_half @ rho @ self._U_H_half.conj().T
    # (2) 順方向Lindbladチャネル
    for U_stine_half in self._U_stines_half:
        rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
    # (3) 逆順Lindbladチャネル
    for U_stine_half in reversed(self._U_stines_half):
        rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
    # (4) 後半ハミルトニアン
    rho = self._U_H_half @ rho @ self._U_H_half.conj().T
    return rho
```

### 4.3 各操作の行列サイズ

| 操作 | 入力 $\rho$ | ユニタリ $U$ | 中間 $\rho_{\mathrm{ext}}$ | 出力 $\rho$ |
|---|---|---|---|---|
| ハミルトニアン | $81 \times 81$ | $81 \times 81$ | — | $81 \times 81$ |
| Stinespring | $81 \times 81$ | $243 \times 243$ | $243 \times 243$ | $81 \times 81$ |

1 Trotterステップ内でのStinespring操作の回数は $26 \times 2 = 52$ 回（順方向26回 + 逆順26回）であり、ハミルトニアン操作は2回（前半・後半）である。

---

## 5. 事前計算

### 5.1 事前計算される行列

`_precompute_unitaries(dt)` で以下の行列が事前計算され、全Trotterステップで再利用される：

1. **半ステップ・ハミルトニアン・ユニタリ**: $U_H^{(1/2)} = e^{-iH\,\Delta t/2} \in \mathbb{C}^{81 \times 81}$

2. **半ステップ・Stinespring ユニタリ**: $U_\alpha^{(1/2)} = \exp(-i\sqrt{\Delta t/2}\,G_\alpha) \in \mathbb{C}^{243 \times 243}$、$\alpha = 1, \ldots, 26$

合計で $1 + 26 = 27$ 個の行列が事前計算される。

### 5.2 計算コスト

| 行列 | サイズ | 計算法 | 計算量 |
|---|---|---|---|
| $U_H^{(1/2)}$ | $81 \times 81$ | `scipy.linalg.expm` | $O(81^3) \approx O(5.3 \times 10^5)$ |
| $U_\alpha^{(1/2)}$ | $243 \times 243$ | `scipy.linalg.expm` | $O(243^3) \approx O(1.4 \times 10^7)$ |
| 全 $U_\alpha$ (26個) | — | — | $O(26 \times 1.4 \times 10^7) \approx O(3.7 \times 10^8)$ |

事前計算のコストは全体の計算量の主要部分の1つであるが、1回のみ実行される。

---

## 6. 完全なシミュレーションアルゴリズム

### 6.1 アルゴリズムの全体フロー

```
入力: params, t_max, n_steps, initial_state
出力: times[], populations[], entropy[], purity[], trace[], rho_final

1. パラメータ検証
   - params.with_boson == False を確認
   
2. 演算子の構成
   - H_0 = build_onsite_hamiltonian(params)          (81×81)
   - H_transfer = build_transfer_hamiltonian(params)  (81×81)
   - H = H_0 + H_transfer                            (81×81)
   - lindblad_ops = build_lindblad_operators(params)   (26個の (81×81, float) タプル)

3. タイムステップの計算
   - dt = t_max / n_steps = 100.0 / 100 = 1.0

4. 事前計算
   - U_H_half = expm(-i * H * dt/2)                  (81×81)
   - U_stines_half[α] = stinespring_unitary_from_lindblad(L_α, dt/2, d_anc=3)
     for α = 1, ..., 26                              (各243×243)

5. 初期状態
   - psi = [0, 0, ..., 0, 1(位置28), 0, ..., 0]     (長さ81のベクトル)
   - rho = |psi><psi|                                 (81×81)

6. 初期物理量の記録
   - times[0] = 0
   - populations[0] = compute_populations(rho, params)
   - entropy[0] = compute_von_neumann_entropy(rho)
   - purity[0] = compute_purity(rho)
   - trace[0] = Tr[rho]

7. 時間発展ループ（step = 0, 1, ..., n_steps-1）
   - rho = trotter_step(rho)    （第4節のアルゴリズム）
   - times[step+1] = (step+1) * dt
   - 物理量の計算と記録

8. ゲート数見積（第8節参照）
```

### 6.2 シミュレーション設定

本ノートブックでは以下の設定が使用される：

| パラメータ | 値 |
|---|---|
| $t_{\mathrm{max}}$ | 100.0（自然単位） |
| $n_{\mathrm{steps}}$ | 100 |
| $\Delta t = t_{\mathrm{max}} / n_{\mathrm{steps}}$ | 1.0 |
| 初期状態 | `'edge_triplet'` = $\|1,0,0,1\rangle$ |

---

## 7. 誤差解析

### 7.1 誤差源の一覧

シナリオ5のTrotterステップには3つの独立した近似誤差源がある：

| 誤差源 | 1ステップ局所誤差 | 大域誤差 |
|---|---|---|
| (A) Strang分割（$\mathcal{L}_H$ vs $\mathcal{L}_D$） | $O(\Delta t^3)$ | $O(\Delta t^2)$ |
| (B) 回文積（個別$\mathcal{L}_{D_\alpha}$の分割） | $O(\Delta t^3)$ | $O(\Delta t^2)$ |
| (C) Stinespring近似（$\mathcal{E}_\alpha$ vs $e^{\mathcal{L}_{D_\alpha}}$） | $O(\Delta t^2)$ | $O(\Delta t)$ |

### 7.2 各誤差源の詳細

#### (A) Strang分割誤差

$$
e^{\mathcal{L}_H \Delta t/2}\, e^{\mathcal{L}_D \Delta t}\, e^{\mathcal{L}_H \Delta t/2} - e^{(\mathcal{L}_H + \mathcal{L}_D)\Delta t} = O(\Delta t^3)
$$

Baker-Campbell-Hausdorff公式により、$\Delta t^2$ の項（$[\mathcal{L}_H, \mathcal{L}_D]$ に比例する項）は対称性によりキャンセルされる。残る主要な誤差項は：

$$
\frac{\Delta t^3}{24}\left([\mathcal{L}_H, [\mathcal{L}_H, \mathcal{L}_D]] + 2[\mathcal{L}_D, [\mathcal{L}_D, \mathcal{L}_H]]\right)
$$

$n_{\mathrm{steps}}$ ステップの積算で $n_{\mathrm{steps}} \times O(\Delta t^3) = O(\Delta t^2)$。

#### (B) 回文積誤差

$$
\left(\prod_{\alpha=1}^{26} e^{\mathcal{L}_{D_\alpha} \Delta t/2}\right)\left(\prod_{\alpha=26}^{1} e^{\mathcal{L}_{D_\alpha} \Delta t/2}\right) - e^{\mathcal{L}_D \Delta t} = O(\Delta t^3)
$$

同様に、対称性により1次の交換子項 $[\mathcal{L}_{D_\alpha}, \mathcal{L}_{D_\beta}]$ がキャンセルされる。$n_{\mathrm{steps}}$ ステップの積算で $O(\Delta t^2)$。

#### (C) Stinespring近似誤差

Part 2 で示したように：

$$
\mathcal{E}_\alpha(\Delta t/2) - e^{\mathcal{L}_{D_\alpha} \Delta t/2} = O((\Delta t/2)^2) = O(\Delta t^2)
$$

26個のチャネルが各ステップで $2 \times 26 = 52$ 回適用されるため、1ステップあたりの累積誤差は $52 \times O(\Delta t^2) = O(\Delta t^2)$（定数は増加するが次数は変わらない）。

$n_{\mathrm{steps}} = t_{\mathrm{max}}/\Delta t$ ステップの積算で：

$$
n_{\mathrm{steps}} \times O(\Delta t^2) = \frac{t_{\mathrm{max}}}{\Delta t} \times O(\Delta t^2) = O(\Delta t)
$$

### 7.3 支配的誤差

3つの誤差源のうち、大域誤差が最も大きいのは **(C) Stinespring近似** であり、$O(\Delta t)$（1次）である。(A) と (B) はともに $O(\Delta t^2)$（2次）であり、(C) に支配される。

したがって、シナリオ5の**有効収束次数は1次**（$O(\Delta t)$）であり、トレース距離における全体誤差は：

$$
D\!\left(\rho_{\mathrm{exact}}(t_{\mathrm{max}}),\; \rho_{\mathrm{Trotter}}(t_{\mathrm{max}})\right) = O(\Delta t)
$$

### 7.4 収束の数値的確認

$\Delta t$ を半分にするとトレース距離がおよそ半分になることが期待される。コードのコメントでは：

> "The Stinespring approximation itself remains 1st-order. As a result, the effective convergence in trace distance is O(dt) (1st-order), not O(dt²)."

---

## 8. ゲート数見積

### 8.1 ゲート数の計算式

`QuditGKSLSimulator.simulate` で返される `estimated_gates_per_step` は、1 Trotterステップあたりの推定量子ゲート数であり、以下のように計算される：

$$
G_{\mathrm{step}} = 2(N + |\mathrm{neighbors}|) + 2 \times n_{\mathrm{Lindblad}}
$$

#### 各項の意味

| 項 | 値 | ゲート種類 | 物理的意味 |
|---|---|---|---|
| $2N$ | $2 \times 4 = 8$ | cu_one（1-qutritゲート） | $H_0$ 対角部分（2半ステップ） |
| $2\|\mathrm{neighbors}\|$ | $2 \times 3 = 6$ | cu_two（2-qutritゲート） | $H_{\mathrm{transfer}}$ ホッピング項（2半ステップ） |
| $2 n_{\mathrm{Lindblad}}$ | $2 \times 26 = 52$ | cu_two/cu_multi（Stinespring） | Lindbladチャネル（順方向+逆順） |

$$
G_{\mathrm{step}} = 8 + 6 + 52 = 66
$$

### 8.2 因数2の由来

ハミルトニアン部分（$2N + 2|\mathrm{neighbors}|$）の因数2は、Strang分割における**2つの半ステップ**（前半と後半のハミルトニアン適用）に対応する。

Lindblad部分（$2 n_{\mathrm{Lindblad}}$）の因数2は、回文積における**順方向と逆順の2回の走査**に対応する。

**コード対応**:

```python
n_lindblad = len(self.lindblad_ops)   # 26
gates_per_step = 2 * (self.n_system_qudits + len(self.params.neighbors)) + n_lindblad * 2
# = 2 * (4 + 3) + 26 * 2 = 14 + 52 = 66
```

### 8.3 全推定ゲート数

$$
G_{\mathrm{total}} = G_{\mathrm{step}} \times n_{\mathrm{steps}} = 66 \times 100 = 6600
$$

### 8.4 qubit エンコーディング（シナリオ3）との比較

qubit エンコーディング（シナリオ3）では、各分子を2 qubit（$d=4$）でエンコードする。ネイティブな2-qubitゲート（CXゲートなど）は1対のqubitにしか作用しないため、2-qutritゲートに相当する操作を実装するには複数のCXゲートが必要となる。そのため、シナリオ3のゲート数はシナリオ5より大幅に多くなる。

### 8.5 qudit量子回路の構成

シナリオ5の各操作は以下のqudit量子回路要素に対応する：

| 操作 | 量子回路要素 | 系qutrit | アンシラqutrit |
|---|---|---|---|
| $U_H^{(1/2)}$のうち$H_0$部分 | 各系qutritへの1体ユニタリ | 1 | 0 |
| $U_H^{(1/2)}$のうち$H_{\mathrm{transfer}}$部分 | 隣接系qutritペアへの2体ユニタリ | 2 | 0 |
| $\mathcal{E}_\alpha^{\Delta t/2}$ (Stinespring) | 系全体+アンシラqutritへの操作 | 4 | 1 |

系のqutrit数: $N = 4$

各Lindbladチャネルに1つのアンシラqutritが必要であり、アンシラのqutrit数は $n_{\mathrm{Lindblad}} = 26$（ただし、アンシラは各使用後に $|0\rangle$ にリセットすれば再利用可能）。

**コード対応**:

```python
self.n_system_qudits = params.N_molecules   # 4
self.n_ancilla_qudits = len(self.lindblad_ops)  # 26
```

---

## 9. 完全な1ステップの計算量

### 9.1 各操作の計算量

1 Trotterステップ内の計算量の内訳：

| 操作 | 回数 | 1回の計算量 | 合計 |
|---|---|---|---|
| $U_H \rho U_H^\dagger$ | 2 | $O(81^2 \times 81) = O(81^3)$ | $O(2 \times 81^3) \approx 10^6$ |
| Stinespring適用 | 52 | $O(243^2 \times 243 + 81^2) = O(243^3)$ | $O(52 \times 243^3) \approx 7.5 \times 10^8$ |

Stinespring適用が計算量の圧倒的大部分を占める。各Stinespring適用では：

1. $\rho_{\mathrm{ext}} = |0\rangle\langle 0| \otimes \rho$: $O(d_{\mathrm{anc}} \cdot d_{\mathrm{sys}}^2) = O(3 \times 81^2)$
2. $\rho' = U \rho_{\mathrm{ext}} U^\dagger$: $O((d_{\mathrm{anc}} d_{\mathrm{sys}})^3) = O(243^3)$
3. 部分トレース: $O(d_{\mathrm{anc}} \cdot d_{\mathrm{sys}}^2) = O(3 \times 81^2)$

### 9.2 全シミュレーションの計算量

$n_{\mathrm{steps}} = 100$ ステップで：

$$
\text{合計} \approx 100 \times 7.5 \times 10^8 \approx 7.5 \times 10^{10} \text{ 浮動小数点演算}
$$

---

## 10. 出力形式

### 10.1 `simulate` メソッドの返り値

`QuditGKSLSimulator.simulate` は以下の辞書を返す：

| キー | 型 | 内容 |
|---|---|---|
| `times` | `list[float]` (長さ101) | 時刻 $t_0 = 0, t_1 = 1, \ldots, t_{100} = 100$ |
| `populations` | `list[dict]` (長さ101) | 各時刻の$N_{\mathrm{S}_0}, N_{\mathrm{T}_1}, N_{\mathrm{S}_1}$と分子別占有数 |
| `entropy` | `list[float]` (長さ101) | フォン・ノイマンエントロピー $S(t_k)$ |
| `purity` | `list[float]` (長さ101) | 純度 $P(t_k)$ |
| `trace` | `list[float]` (長さ101) | トレース $\mathrm{Tr}[\rho(t_k)]$ |
| `rho_final` | `np.ndarray` ($81 \times 81$) | 最終時刻の密度行列 $\rho(t_{\mathrm{max}})$ |
| `elapsed_time` | `float` | 計算時間（秒） |
| `method` | `str` | `'qudit_gksl'` |
| `params` | `dict` | 物理パラメータの辞書表現 |
| `n_system_qudits` | `int` | 4 |
| `n_ancilla_qudits` | `int` | 26 |
| `d_anc` | `int` | 3 |
| `estimated_gates_per_step` | `int` | 66 |
| `total_estimated_gates` | `int` | 6600 |

### 10.2 ノートブックのセル6での出力

セル6では以下の情報が出力される：

```
Elapsed time: {result5["elapsed_time"]:.2f}s
System qutrits: 4
Ancilla qudits: 26
Estimated gates/step: 66
Total estimated gates: 6600
Trace conservation: max|Tr-1| = {max(abs(t-1) for t in result5["trace"]):.2e}
```

および `plot_population_dynamics` によるポピュレーション・ダイナミクスのグラフ。

---

## 11. シナリオ1（古典GKSL）との数学的関係

### 11.1 同一の物理モデル

シナリオ1とシナリオ5は**同一のGKSL方程式**を解く。具体的には：
- 同一のハミルトニアン $H = H_0 + H_{\mathrm{transfer}}$（$81 \times 81$）
- 同一のLindblad演算子 $\tilde{L}_\alpha$（26個、各 $81 \times 81$）
- 同一の初期状態 $\rho(0) = |1,0,0,1\rangle\langle 1,0,0,1|$

### 11.2 異なる数値解法

| 特徴 | シナリオ1（古典GKSL） | シナリオ5（Qudit GKSL） |
|---|---|---|
| 解法 | リウヴィリアン超演算子 $\mathcal{L}$ の行列指数関数 | Trotter分割 + Stinespring |
| 行列サイズ | $6561 \times 6561$（$\mathcal{L}$） | $81 \times 81$（$\rho$）+ $243 \times 243$（$U_\alpha$） |
| 精度 | 浮動小数点精度まで厳密 | $O(\Delta t)$ |
| 量子回路対応 | なし | あり |
| CPTP保証 | 数学的に厳密 | 各ステップで厳密（Kraus形式） |

### 11.3 $\Delta t \to 0$ の極限

$$
\lim_{\Delta t \to 0} \rho_{\mathrm{Trotter}}(t) = \rho_{\mathrm{exact}}(t)
$$

すなわち、ステップ数 $n_{\mathrm{steps}} \to \infty$ とすると、シナリオ5の結果はシナリオ1の結果に収束する。有効収束速度は $O(\Delta t)$（1次）である。

---

## 付録A. Trotter-Suzuki分解の一般論

### A.1 1次Lie-Trotter分割

$$
e^{(A+B)\tau} = e^{A\tau} e^{B\tau} + O(\tau^2)
$$

$$
= e^{(A+B)\tau + \frac{\tau^2}{2}[A,B] + O(\tau^3)}
$$

### A.2 2次対称（Strang）分割

$$
e^{(A+B)\tau} = e^{A\tau/2} e^{B\tau} e^{A\tau/2} + O(\tau^3)
$$

$$
= e^{(A+B)\tau + \frac{\tau^3}{24}([A,[A,B]] + 2[B,[B,A]]) + O(\tau^5)}
$$

$\tau^2$ の交換子項がキャンセルするため、2次精度が達成される。

### A.3 回文積の一般化

$n$ 個の演算子 $A_1, A_2, \ldots, A_n$ に対する回文積：

$$
\left(\prod_{k=1}^{n} e^{A_k \tau/2}\right)\left(\prod_{k=n}^{1} e^{A_k \tau/2}\right) = e^{(A_1 + A_2 + \cdots + A_n)\tau} + O(\tau^3)
$$

これは多体分割の対称化であり、2次精度を保証する。

## 付録B. Stinespring近似の高次項

$U_\alpha = e^{-i\sqrt{\Delta t}\,G_\alpha}$ のKraus演算子の3次までの展開：

$$
K_0 = I - \frac{\Delta t}{2}\tilde{L}_\alpha^\dagger \tilde{L}_\alpha + \frac{\Delta t^2}{8}(\tilde{L}_\alpha^\dagger \tilde{L}_\alpha)^2 + O(\Delta t^{5/2})
$$

$$
K_1 = -i\sqrt{\Delta t}\,\tilde{L}_\alpha + i\frac{\Delta t^{3/2}}{6}\tilde{L}_\alpha\,\tilde{L}_\alpha^\dagger \tilde{L}_\alpha + O(\Delta t^{5/2})
$$

$K_2 = 0$（厳密。$d_{\mathrm{anc}} = 3$ のブロック対角構造より）。

チャネル $\mathcal{E}_\alpha$ の2次項：

$$
\mathcal{E}_\alpha(\Delta t)(\rho) = \rho + \Delta t\,\mathcal{D}[\tilde{L}_\alpha](\rho) + \frac{\Delta t^2}{4}\left(\tilde{L}_\alpha^\dagger \tilde{L}_\alpha\,\rho\,\tilde{L}_\alpha^\dagger \tilde{L}_\alpha\right) + O(\Delta t^3)
$$

一方、厳密なLindblad時間発展の2次項は：

$$
e^{\mathcal{L}_{D_\alpha}\Delta t}(\rho) = \rho + \Delta t\,\mathcal{D}[\tilde{L}_\alpha](\rho) + \frac{\Delta t^2}{2}\mathcal{D}[\tilde{L}_\alpha]^2(\rho) + O(\Delta t^3)
$$

2次項の不一致がStinespring近似の1次誤差を決定する。

---

## 参考文献

- H. F. Trotter, "On the product of semi-groups of operators," Proc. Amer. Math. Soc. **10**, 545 (1959).
- G. Strang, "On the construction and comparison of difference schemes," SIAM J. Numer. Anal. **5**, 506 (1968).
- M. Suzuki, "Generalized Trotter's formula and systematic approximants of exponential operators and inner derivations with applications to many-body problems," Commun. Math. Phys. **51**, 183 (1976).
- A. M. Childs and Y. Su, "Nearly optimal lattice simulation by product formulas," Phys. Rev. Lett. **123**, 050503 (2019).

---

*本文書群（Part 1 ~ Part 3）により、シナリオ5: Qudit GKSL（ボソン無し）の理論全体が、プログラムで計算可能な解像度で省略無しに定式化された。*
