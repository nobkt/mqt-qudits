# なぜ Stinespring 処方 + 対称 Trotter は MQT-Qudits の backend で実装できないのか

> **方針：** このリポジトリの実コードを根拠とし、事実のみを記述する。
> 推測・希望的観測・曖昧な記述は一切含めない。

---

## 1. 結論（1 行）

**Stinespring 処方はアンシラ（補助 qudit）の追加と partial trace（偏トレース）が不可欠だが、
MQT-Qudits の circuit/backend API はどちらの操作も提供していない。
Kraus 演算子はその partial trace を事前に行列として吸収しているため、
`KrausChannel` 命令として DMSim backend に直接渡せる。**

---

## 2. 背景：2 種類の CPTP 写像の表現

### 2.1 Stinespring 処方

任意の CPTP 写像 $\mathcal{E}$ は補助系（アンシラ）を導入して

$$
\mathcal{E}(\rho)
= \mathrm{tr}_\mathrm{anc}\!\bigl[
    U_\mathrm{Stine}
    \,(\rho \otimes |0\rangle\langle0|_\mathrm{anc})\,
    U_\mathrm{Stine}^\dagger
  \bigr]
$$

と書ける（Stinespring 定理）。計算には**3 ステップ**が必要：

1. アンシラを $|0\rangle$ で初期化し、$\rho_\mathrm{ext} = |0\rangle\langle0|_\mathrm{anc} \otimes \rho$ を作る
2. 拡大空間のユニタリー $U_\mathrm{Stine}$ を適用する
3. アンシラを partial trace で**消去**する

### 2.2 Kraus 表現

同じ CPTP 写像は Kraus 演算子 $\{K_k\}$（$\sum_k K_k^\dagger K_k = I$）を使って

$$
\mathcal{E}(\rho) = \sum_{k} K_k \rho K_k^\dagger
$$

と書ける。Stinespring との関係は

$$
K_k = {}_\mathrm{anc}\langle k|\, U_\mathrm{Stine}\, |0\rangle_\mathrm{anc}
$$

であり、partial trace は $K_k$ を**定義する時点で数学的に完了**している。
$K_k$ は $d_\mathrm{sys} \times d_\mathrm{sys}$ の行列だけで記述できる。

---

## 3. MQT-Qudits backend が受け付ける命令

`src/mqt/qudits/simulation/backends/dmsim.py:_apply_instruction` のコード：

```python
if isinstance(instruction, KrausChannel):
    return apply_kraus_to_density(rho, instruction.kraus_operators, qudits, dims)

# ユニタリーゲートパス
u_matrix = instruction.to_matrix(identities=0)
return apply_unitary_to_density(rho, u_matrix, qudits, dims)
```

DMSim backend が受け付けるのは**2 種類だけ**：

| 命令の種類 | 内部処理 |
|-----------|---------|
| ユニタリーゲート（`Gate.to_matrix()` が返す行列） | $\rho \to U \rho U^\dagger$ |
| `KrausChannel`（Kraus 演算子リスト） | $\rho \to \sum_k K_k \rho K_k^\dagger$ |

**partial trace 命令は存在しない。ancilla qudit の追加・削除命令も存在しない。**

### 3.1 `NoiseModel` オプションも拒否される

`dmsim.py:220-224`:

```python
if self._options.get("noise_model", None) is not None:
    msg = (
        "DMSim does not accept a NoiseModel option. "
        "Add noise as KrausChannel instructions in the circuit instead."
    )
    raise ValueError(msg)
```

---

## 4. Stinespring 処方が backend で動かない理由（コードの事実）

### 4.1 回路次元が固定されている

`QuantumCircuit` の qudit 数と次元は回路生成時に確定し、途中で変更できない。

```python
qreg = QuantumRegister("sys", n, [d] * n)
circuit = QuantumCircuit(qreg)   # ← n qudits, dimension d each, 固定
```

Stinespring では各 Lindblad チャンネルごとに「アンシラ 1 本を追加→ユニタリー適用→アンシラを捨てる」を
繰り返さなければならない。1 Trotter ステップに 26 チャンネルあるなら、アンシラの追加・削除を 26 回行う必要がある。
この動的な次元変更は `QuantumCircuit` に実装されていない。

### 4.2 partial trace 命令がない

Stinespring ステップの最後には必ず

$$
\rho_\mathrm{sys}^{(k+1)} = \mathrm{tr}_\mathrm{anc}\!\bigl[U \rho_\mathrm{ext} U^\dagger\bigr]
$$

が必要だが、DMSim が実行できる命令は「行列をかける」だけであり、
「qudit インデックスを指定して偏トレースを取る」命令は存在しない。

### 4.3 コードに明示的なエラーが書かれている

`tutorials/qudit_gksl_simulator.py:119-126`:

```python
if algorithm != "exact_local_channels":
    msg = (
        "execute_on_backend='dmsim' requires "
        "algorithm='exact_local_channels' so that each Lindblad "
        "channel has an exact Kraus representation that can be "
        "fed to MQT-Qudits as a KrausChannel instruction.  "
        "Stinespring dilation has no Kraus representation on "
        "the system alone (it requires ancillas)."
    )
    raise ValueError(msg)
```

`algorithm="stinespring"` と `execute_on_backend="dmsim"` を同時に指定すると
**`ValueError` が上がって実行不可能**になる。これはアーキテクチャ上の制約を明示したものである。

---

## 5. 対称 Trotter（回文順序）は Stinespring の精度問題を救えない

ここが「対称 Trotter があれば 2 次精度になるはずでは？」という誤解の核心である。

### 5.1 対称 Trotter（Strang 分割）の効果

ハミルトニアン部と散逸部を **H-D-D-H** の順（回文）で適用すると、
ハミルトニアン–散逸間の交換子誤差（Lie-Trotter 積の 1 次誤差）は相殺され、
**ハミルトニアン–散逸の分割に関しては $O(\Delta t^2)$** の精度が得られる。

### 5.2 Stinespring 近似がかける精度限界

しかし、Stinespring 処方による各 Lindblad チャンネルの近似：

$$
U_\mathrm{Stine} = \exp\!\left(-i\sqrt{\tfrac{\Delta t}{2}}\,G\right),
\quad
G = \begin{pmatrix}0 & L^\dagger \\ L & 0\end{pmatrix}
$$

を展開すると、対応する Kraus 演算子は：

$$
K_0 = I - \frac{\Delta t}{4}L^\dagger L + O(\Delta t^2)
$$

$$
K_1 = -i\sqrt{\frac{\Delta t}{2}}\,L + O(\Delta t^{3/2})
$$

CPTP 写像を計算すると：

$$
\mathcal{E}(\rho)
= \rho
+ \frac{\Delta t}{2}\!\left(L\rho L^\dagger
  - \tfrac{1}{2}L^\dagger L\rho
  - \tfrac{1}{2}\rho L^\dagger L\right)
+ O(\Delta t^{3/2})
$$

厳密な Lindblad チャンネル $e^{\mathcal{D}_\alpha \Delta t/2}$ との差は $O(\Delta t^{3/2})$ から始まる。
これは対称 Trotter の回文順序がどれだけ精巧であっても**消えない誤差**である。

### 5.3 実効的な精度

`tutorials/qudit_gksl_simulator.py:24-29`（docstring）に明記されている：

```
the effective convergence in trace distance is O(dt) (1st-order), not O(dt²).
```

対称 Trotter は「ハミルトニアンと散逸の順序誤差」を相殺するが、
「各散逸チャンネルの Stinespring 近似誤差（$O(\Delta t^{3/2})$）」は相殺できない。
**Stinespring 近似が 1 次精度の壁を作っており、対称 Trotter はその壁を突破できない。**

---

## 6. Kraus 理論が backend で動く理由

### 6.1 exact_local_channels + Choi–Jamiołkowski 変換

`algorithm="exact_local_channels"` モードでは、各 Lindblad チャンネルの
局所超演算子を**厳密に指数関数化**する：

$$
\mathcal{E}_\alpha\!\left(\frac{\Delta t}{2}\right)(\rho)
= e^{\mathcal{L}_{D_\alpha}^\mathrm{local}\,\Delta t/2}(\rho)
$$

- 単サイトチャンネル：$\mathcal{L}_{D_\alpha}^\mathrm{local}$ は $9 \times 9$（$d^2 \times d^2$、$d=3$）
- ペアチャンネル：$\mathcal{L}_{D_\alpha}^\mathrm{local}$ は $81 \times 81$

この厳密指数関数を Choi–Jamiołkowski 同型経由で Kraus 演算子に変換する
（`tutorials/dmsim_kraus_helpers.py:kraus_from_local_superoperator`）：

$$
C = \sum_\alpha \lambda_\alpha |v_\alpha\rangle\langle v_\alpha|
\quad\Rightarrow\quad
K_\alpha = \sqrt{\lambda_\alpha}\cdot \mathrm{reshape}(|v_\alpha\rangle,\,(d,d))
$$

完全性条件 $\sum_\alpha K_\alpha^\dagger K_\alpha = I$ は数値検証済み（許容誤差 $10^{-10}$）。

### 6.2 backend への渡し方

`tutorials/qudit_gksl_simulator.py:_trotter_step_dmsim`:

```python
# ハミルトニアン半ステップ（全 N-qudit ユニタリー）
circuit.cu_multi(list(range(n)), self._U_H_half)

# 前向きパス（Lindblad チャンネル × n_channels）
for sites, kraus_ops, kind in self._dmsim_kraus_half:
    circuit.kraus_channel(target, kraus_ops)   # KrausChannel 命令

# 逆向きパス（回文順序）
for sites, kraus_ops, kind in reversed(self._dmsim_kraus_half):
    circuit.kraus_channel(target, kraus_ops)

# ハミルトニアン半ステップ（閉じ）
circuit.cu_multi(list(range(n)), self._U_H_half)

job = self._dmsim_backend.run(circuit, initial_density_matrix=rho)
```

ここで `kraus_ops` は $d_\mathrm{sys} \times d_\mathrm{sys}$ の行列リストであり、
アンシラも partial trace も一切不要。DMSim は `apply_kraus_to_density` を呼ぶだけで処理できる。

### 6.3 精度

Stinespring 近似を使わないため、各チャンネルは厳密。
回文順序により散逸チャンネル間の交換子誤差も相殺される。
実効精度は**$O(\Delta t^2)$（2 次）**（`tutorials/exact_local_channels.py:54-55` の実測値）。

---

## 7. 真実ベースの対比まとめ

| 項目 | Stinespring + 対称 Trotter | Kraus（exact_local_channels） |
|------|--------------------------|------------------------------|
| backend への命令 | アンシラ追加 + ユニタリー + partial trace | `KrausChannel`（行列リストのみ） |
| DMSim が partial trace を持つか | **持たない** | 不要 |
| 回路次元の動的変更 | 必要（チャンネルごとにアンシラ追加・削除） | 不要（固定次元） |
| `execute_on_backend="dmsim"` | **ValueError で拒否**（コード明記） | 動作する |
| 各チャンネルの精度 | $O(\Delta t^{3/2})$ 近似 | 厳密 |
| 対称 Trotter の効果 | H-D 分割誤差は相殺、Stinespring 近似誤差は**残る** | H-D + チャンネル間誤差、両方相殺 |
| 全体の精度 | $O(\Delta t)$（1 次） | $O(\Delta t^2)$（2 次） |
| アンシラ次元（チュートリアル実装） | $d_\mathrm{anc}=3$（クトリット） | 不要 |

---

## 8. 補足：なぜ「対称 Trotter があれば精度が上がるはず」は誤りか

対称 Trotter（Strang 分割）の精度向上は**分割する 2 つの生成子が非可換な場合の交換子誤差の相殺**に由来する。

$$
e^{(A+B)h} = e^{Ah/2} e^{Bh} e^{Ah/2} + O(h^3 [A,[A,B]] + h^3 [B,[B,A]])
$$

ここで $A = \mathcal{L}_H$、$B = \mathcal{L}_D$（散逸部）。

しかし Stinespring 処方では「$e^{Bh}$ を厳密に計算する」のではなく、
$U_\mathrm{Stine} = e^{-i\sqrt{h/2}\,G}$ という**近似写像**で代用している。
この近似誤差は Trotter の誤差解析の外側にあり、対称 Trotter の定理は適用できない。

「Stinespring 処方を使った場合の対称 Trotter が 2 次精度にならない」ことは
このリポジトリの docstring（`qudit_gksl_simulator.py:24-29`）に明確に記述されており、
実験的な収束測定（`tutorials/test_gksl_simulators.py`）でも確認されている。

---

## 9. コード参照箇所一覧

| 内容 | ファイル | 行（概略） |
|------|---------|-----------|
| Stinespring + DMSim の組合せを拒否する `ValueError` | `tutorials/qudit_gksl_simulator.py` | 119–126 |
| DMSim が受け付ける命令（2 種類のみ） | `src/mqt/qudits/simulation/backends/dmsim.py` | 272–288 |
| DMSim が `NoiseModel` を拒否する | `src/mqt/qudits/simulation/backends/dmsim.py` | 220–224 |
| `KrausChannel.to_matrix()` が `NotImplementedError` | `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py` | 147–155 |
| Stinespring ユニタリの実装 | `tutorials/stinespring_utils.py` | 14–51 |
| partial trace の実装（NumPy 専用） | `tutorials/stinespring_utils.py` | 54–83 |
| Stinespring 精度が O(dt) であることの明記 | `tutorials/qudit_gksl_simulator.py` | 24–29 |
| `_trotter_step_dmsim`（KrausChannel で backend 実行） | `tutorials/qudit_gksl_simulator.py` | 256–312 |
| Choi → Kraus 変換 | `tutorials/dmsim_kraus_helpers.py` | 51–129 |
| exact_local_channels の精度（O(dt²)、実測値） | `tutorials/exact_local_channels.py` | 54–55 |
