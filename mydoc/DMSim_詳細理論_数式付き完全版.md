# DMSim 詳細理論（省略なし・数式付き完全版）

> **本文書はコード実装の直接引用に基づく。数学的主張はすべてコード行番号で裏付けられている。**  
> コードから読めない事項には明示的に「コード外」と注記する。  
> 引用元ファイル:
> - `src/mqt/qudits/simulation/backends/dmsim.py`（全288行）  
> - `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`（全165行）  
> - `src/mqt/qudits/quantum_circuit/circuit.py`（`dimensions` プロパティ）

---

## 1. システムの表現：密度行列

### 1.1 Hilbert 空間の次元

MQT-Qudits の量子回路は N 個の qudit から成り、i 番目の qudit の局所次元を $d_i \in \mathbb{Z}_{>0}$ とする。
全体の Hilbert 空間の次元は

$$
D = \prod_{i=0}^{N-1} d_i
$$

である（`dmsim.py:237`: `d_total = int(reduce(operator.mul, dims, 1))`）。  
各 $d_i$ は独立に設定でき（`circuit.py:113`: `dims = num_qudits * [2] if args[1] is None else args[1]`）、  
qubit（$d=2$）、qutrit（$d=3$）、ququart（$d=4$）を混在させることができる。

### 1.2 密度行列の定義

系の状態を **密度演算子**（density operator）

$$
\rho \in \mathbb{C}^{D \times D}
$$

で表す。物理的な密度行列は以下の 3 条件を満たす（コードは構築時に CPTP 検証のみ行い、ρ 自体の Hermiticity・半正定値性は検証しない。`dmsim.py:197–200`）:

1. **Hermiticity**: $\rho = \rho^\dagger$  
2. **半正定値性**: $\rho \geq 0$（すべての固有値 $\geq 0$）  
3. **トレース規格化**: $\mathrm{tr}(\rho) = 1$

### 1.3 ストレージ形式

DMSim はρを `(D, D)` の `numpy.ndarray` (dtype=`complex128`) として保持する。  
ストレージ容量は

$$
\text{メモリ} = 16 D^2 \text{ bytes}
$$

（`complex128` = 16 bytes/要素、`dmsim.py:32–33`）。  
N=4 qutrit 系（$D=81$）では $16 \times 81^2 \approx 105\,\text{KB}$。

### 1.4 テンソル表現

DMSim は内部計算で ρ を **2N-way テンソル** として扱う。

$$
\rho_{\text{tensor}} \in \mathbb{C}^{d_0 \times d_1 \times \cdots \times d_{N-1} \times d_0 \times d_1 \times \cdots \times d_{N-1}}
$$

インデックス規約: 最初の $N$ 軸が「行レグ」、後の $N$ 軸が「列レグ」。

変換式: `rho_tensor = rho.reshape(*([d for d in dims] * 2))` (`dmsim.py:141,158`)

$$
(\rho_{\text{tensor}})_{i_0 i_1 \cdots i_{N-1},\; j_0 j_1 \cdots j_{N-1}} = \rho_{\,i_0 d_1 d_2 \cdots + \cdots,\;\, j_0 d_1 d_2 \cdots + \cdots}
$$

---

## 2. 初期状態の構築

### 2.1 デフォルト初期状態

引数 `initial_density_matrix=None` のとき、DMSim は計算基底の $|0\rangle$ 状態から純粋密度行列を構築する（`dmsim.py:239–242`）:

$$
|\psi_0\rangle = |0\rangle^{\otimes N} = e_0 \in \mathbb{C}^D \quad (e_0 = (1, 0, 0, \ldots, 0)^\top)
$$

$$
\rho_0 = |\psi_0\rangle\langle\psi_0| = e_0 e_0^\top
$$

実装:
```python
psi0 = np.zeros(d_total, dtype=np.complex128)
psi0[0] = 1.0
rho = np.outer(psi0, psi0.conj())
```

### 2.2 任意初期状態

`initial_density_matrix` として `(D, D)` 複素配列を渡すと、その配列（のコピー）が初期 ρ として使われる（`dmsim.py:243–252`）。  
**コードは shape チェックのみ行い、Hermiticity・半正定値性・トレース規格化は検証しない。** 不正な ρ を渡しても `ValueError` にならず、そのまま伝播する。

---

## 3. ユニタリゲートの適用：$\rho \to U_S \rho U_S^\dagger$

### 3.1 全体作用素としての表現

ゲート $U$ が qudit 部分集合 $S = \{q_0, q_1, \ldots, q_{k-1}\}$ に作用するとき（ただし各 $q_j \in \{0,\ldots,N-1\}$）、
Hilbert 空間全体への作用は

$$
U_{\text{total}} = I_{d_0} \otimes \cdots \otimes U_{d_{q_0}} \otimes \cdots \otimes I_{d_{N-1}}
$$

の形の局所ユニタリ拡張であるが、DMSim はこれを明示的に $D \times D$ に展開しない。  
代わりにテンソルレグ操作で実現する（下記参照）。

### 3.2 テンソルレグによる 2 段階計算

`apply_unitary_to_density(rho, u_matrix, qudits, dims)` (`dmsim.py:132–144`) は以下の 2 ステップを実行する:

**ステップ 1（行レグへの作用）**: $\rho_{\text{tensor}} \leftarrow (U \otimes I_{\bar S}) \cdot \rho_{\text{tensor}}$

テンソル縮約:

$$
(\rho'_{\text{tensor}})_{i_0 \cdots i_{N-1},\; j_0 \cdots j_{N-1}}
= \sum_{\ell_0,\ldots,\ell_{k-1}} U_{i_{q_0}\cdots i_{q_{k-1}},\; \ell_0\cdots\ell_{k-1}} 
  \cdot (\rho_{\text{tensor}})_{\underbrace{\cdots \ell_0 \cdots \ell_{k-1} \cdots}_{\text{位置}S},\; j_0 \cdots j_{N-1}}
$$

実装: `np.tensordot(op_t, rho_tensor, axes=(range(k, 2k), qudits))` (`dmsim.py:121`)

**ステップ 2（列レグへの作用）**: $\rho_{\text{tensor}} \leftarrow \rho'_{\text{tensor}} \cdot (U^* \otimes I_{\bar S})$

$$
(\rho''_{\text{tensor}})_{i_0 \cdots i_{N-1},\; j_0 \cdots j_{N-1}}
= \sum_{\ell_0,\ldots,\ell_{k-1}} U^*_{j_{q_0}\cdots j_{q_{k-1}},\; \ell_0\cdots\ell_{k-1}}
  \cdot (\rho'_{\text{tensor}})_{i_0 \cdots i_{N-1},\; \underbrace{\cdots \ell_0 \cdots \ell_{k-1} \cdots}_{\text{位置}S}}
$$

実装: `np.tensordot(op_t_conj, rho_tensor, axes=(range(k, 2k), [n+q for q in qudits]))` (`dmsim.py:121,113`)

最終結果を $D \times D$ に reshape して返す（`dmsim.py:144`）:

$$
\rho_{\text{new}} = U_S \,\rho\, U_S^\dagger
$$

### 3.3 テンソル縮約後の軸並び替え

`tensordot` の結果は "op の出力レグ" が先頭に出る順序になる。
これを元の `(d_0, d_1, …, d_{N-1}, d_0, …, d_{N-1})` の軸順に戻すために **逆置換** (`np.transpose`) を適用する（`dmsim.py:122–129`）。

置換の構築:
```
contract_axes = qudits  (または [n+q for q in qudits])
remaining = {0,…,2N-1} \ contract_axes  (順序を保つ)
tensordot 結果の軸順: [contract_axes[0], …, contract_axes[k-1], remaining[0], …]
逆置換 perm[label] = new_position
np.transpose(out, perm)
```

### 3.4 計算コスト（ユニタリゲート）

- ローカル次元: $D_{\text{local}} = \prod_{j} d_{q_j}$
- `tensordot` の縮約コスト: $O(D \cdot D_{\text{local}}^2)$（行レグ・列レグ各 1 回ずつ）
- 合計: $O(D \cdot D_{\text{local}}^2)$（`dmsim.py:35`）

単 qudit ゲート（$k=1$, $D_{\text{local}} = d_{q_0}$）の場合: $O(D \cdot d_{q_0}^2)$。

---

## 4. Kraus チャネルの適用：$\rho \to \sum_k K_k \rho K_k^\dagger$

### 4.1 CPTP 写像の定義

CPTP（完全正かつトレース保存）写像は Kraus 表現で

$$
\mathcal{E}(\rho) = \sum_{k=0}^{r-1} K_k \,\rho\, K_k^\dagger
$$

と書かれる。ここで $\{K_k\}_{k=0}^{r-1}$ は **Kraus 演算子** と呼ばれ、トレース保存条件

$$
\sum_{k=0}^{r-1} K_k^\dagger K_k = I_{D_{\text{local}}}
$$

を満たす。$r$ を Kraus **ランク**（Choi ランク）と呼ぶ。

### 4.2 CPTP 条件の検証（構築時・厳密・補正なし）

`KrausChannel.__init__` はコンストラクタ内で CPTP 条件を数値的に検証する（`kraus_channel.py:108–119`）:

$$
\Delta = \left\| \sum_{k=0}^{r-1} K_k^\dagger K_k - I \right\|_F
$$

- 許容値: `CPTP_TOLERANCE = 1e-9`（`kraus_channel.py:37`）
- $\Delta > 10^{-9}$ ならば `ValueError` を送出
- **補正は一切行わない**（明示的に「No heuristic correction is applied」と記される）

実装:
```python
completeness = sum(K.conj().T @ K for K in kraus_list)
defect = np.linalg.norm(completeness - identity, ord='fro')
if defect > CPTP_TOLERANCE:
    raise ValueError(f"||Σ K†K - I||_F = {defect:.3e}")
```

### 4.3 テンソルレグによる Kraus 適用

`apply_kraus_to_density(rho, kraus_ops, qudits, dims)` (`dmsim.py:147–162`) は以下を計算する:

$$
\rho_{\text{new}} = \sum_{k=0}^{r-1} K_k^{(S)} \,\rho\, (K_k^{(S)})^\dagger
$$

アルゴリズム:

```
accumulator = 0_{D×D}
for k in range(r):
    rho_tensor = rho.reshape(*dims, *dims)
    # ステップ1: 行レグに K_k を作用
    rho_tensor = _apply_local_op_one_side(rho_tensor, K_k, qudits, dims, 'row')
    # ステップ2: 列レグに K_k* を作用
    rho_tensor = _apply_local_op_one_side(rho_tensor, K_k.conj(), qudits, dims, 'col')
    accumulator += rho_tensor.reshape(d_total, d_total)
return accumulator
```

各 Kraus 項 $K_k \rho K_k^\dagger$ の計算はユニタリゲートの場合と同一の `_apply_local_op_one_side` を使い、$r$ 回繰り返して加算する。

**重要**: この処理はテンソル形式の「$K_k$ の行列版」を局所 qudit 部分集合 $S$ に埋め込んだ作用

$$
(\mathcal{E}_k \otimes \mathrm{id}_{\bar S})(\rho)
$$

であり、環境 qudit $\bar S$ は恒等作用を受ける。

### 4.4 計算コスト（Kraus チャネル）

- Kraus ランク: $r$  
- 1 項あたり: $O(D \cdot D_{\text{local}}^2)$  
- 合計: $O(r \cdot D \cdot D_{\text{local}}^2)$（`dmsim.py:36`）

---

## 5. `_apply_local_op_one_side` の完全なアルゴリズム

この関数がDMSim の核心計算を担う。完全な動作を記述する（`dmsim.py:75–129`）。

**入力:**
- `rho_tensor`: shape `(d_0, …, d_{N-1}, d_0, …, d_{N-1})` の 2N-way テンソル  
- `op_matrix`: shape `(D_local, D_local)`、`D_local = prod(dims[q] for q in qudits)`  
- `qudits`: 作用先 qudit のインデックスリスト `[q_0, …, q_{k-1}]`  
- `dims`: 全 qudit の次元リスト  
- `side`: `'row'`（行レグ）または `'col'`（列レグ）

**処理:**

1. `op_matrix` をテンソル形式に reshape:

$$
\text{op\_t} \leftarrow \text{op\_matrix.reshape}(d_{q_0}, \ldots, d_{q_{k-1}}, d_{q_0}, \ldots, d_{q_{k-1}})
$$

形状: $(d_{q_0}, \ldots, d_{q_{k-1}}, d_{q_0}, \ldots, d_{q_{k-1}})$

2. 縮約対象軸の決定:
   - `side='row'`: `contract_axes = [q_0, …, q_{k-1}]`（行レグの対応軸）
   - `side='col'`: `contract_axes = [N+q_0, …, N+q_{k-1}]`（列レグの対応軸）

3. `tensordot` による縮約:

$$
\text{out} = \text{tensordot}(\text{op\_t},\; \text{rho\_tensor},\; \text{axes}=(\{k, k+1, \ldots, 2k-1\},\; \text{contract\_axes}))
$$

`op_t` の後半 $k$ 軸（入力レグ）と `rho_tensor` の `contract_axes` を縮約する。  
結果の軸順: `(op_t の前半 k 軸 = 出力レグ, rho_tensor の残り 2N-k 軸)`

4. 軸並び替え（逆置換）:

`tensordot` 後の軸ラベル配置を元の `(0, 1, …, 2N-1)` 順に戻す逆置換 `perm` を計算し `np.transpose(out, perm)` を適用する。

$$
\text{perm}[i] = \text{tensordot 結果における元のラベル } i \text{ の位置}
$$

具体的には `contract_axes` が先頭 $k$ 位置に出るため、  
`label_to_pos[contract_axes[j]] = j`、`label_to_pos[remaining[j]] = k+j` として  
`perm = [label_to_pos[i] for i in range(2N)]` を構築する。

**出力:** 並び替え後のテンソル、shape `(d_0, …, d_{N-1}, d_0, …, d_{N-1})`（元と同じ）

---

## 6. 命令ループ：回路の逐次実行

`DMSim.execute` は回路の命令列を順番に処理する（`dmsim.py:254–257`）:

```python
for instruction in circuit.instructions:
    rho = self._apply_instruction(rho, instruction, dims, n)
```

各命令で `_apply_instruction` が呼ばれ（`dmsim.py:259–288`）:

- `isinstance(instruction, KrausChannel)` が `True` → `apply_kraus_to_density` を呼ぶ  
- それ以外 → `instruction.to_matrix(identities=0)` でユニタリ行列を取得し `apply_unitary_to_density` を呼ぶ

**注意**: `identities=0`（`dmsim.py:276`）は「qudit ごとの部分単位行列を埋め込まずに局所行列だけ返す」指定であり、Hilbert 空間全体への $U_{\text{total}}$ の展開は行われない。局所テンソル演算がその役割を担う。

---

## 7. トレース保存性の伝播（保証と非保証）

### 7.1 ユニタリゲート

ユニタリゲート $U_S$（$U_S^\dagger U_S = I$）のもとで:

$$
\mathrm{tr}(U_S \rho U_S^\dagger) = \mathrm{tr}(U_S^\dagger U_S \rho) = \mathrm{tr}(\rho)
$$

したがって初期 $\rho$ のトレースが 1 であれば、ユニタリゲート適用後も 1 を保つ（数学的保証、コードは検証なし）。

### 7.2 Kraus チャネル

Kraus 演算子が CPTP 条件 $\sum_k K_k^\dagger K_k = I$ を満たすとき:

$$
\mathrm{tr}\!\left(\sum_k K_k \rho K_k^\dagger\right) = \sum_k \mathrm{tr}(K_k \rho K_k^\dagger) = \sum_k \mathrm{tr}(K_k^\dagger K_k \rho) = \mathrm{tr}\!\left(\!\left(\sum_k K_k^\dagger K_k\right)\rho\right) = \mathrm{tr}(I \cdot \rho) = \mathrm{tr}(\rho)
$$

`KrausChannel` コンストラクタの CPTP 検証（Section 4.2）により、$\|\sum_k K_k^\dagger K_k - I\|_F \leq 10^{-9}$ が構築時に保証される。
ただし **DMSim はチャネル適用後の ρ のトレースを再検証しない**。浮動小数点誤差の蓄積は回路長とともに増大するが、コードはそれを追跡しない。

---

## 8. NoiseModel の拒否設計

`DMSim.run` は `noise_model` オプションが `None` でない場合に即座に `ValueError` を送出する（`dmsim.py:219–224`）:

```python
if self._options.get("noise_model", None) is not None:
    raise ValueError(
        "DMSim does not accept a NoiseModel option. "
        "Add noise as KrausChannel instructions in the circuit instead."
    )
```

**設計意図**（docstring より）: 「noise must be supplied as KrausChannel instructions so that the simulator and the user agree on what the physics is」

これは意図的な制約であり、ノイズは常に明示的な `KrausChannel` 命令として回路に組み込む必要がある。

---

## 9. KrausChannel 命令の制約

### 9.1 `to_matrix` は使用不可

`KrausChannel.to_matrix()` は `NotImplementedError` を送出する（`kraus_channel.py:154–159`）。  
同様に `__array__` も `NotImplementedError`（`kraus_channel.py:147–152`）。

これは CPTP 写像一般がユニタリ行列の形で表現できないことに起因する。  
**例**: 振幅減衰チャネル $K_0 = |0\rangle\langle 0| + \sqrt{1-\gamma}|1\rangle\langle 1|$、$K_1 = \sqrt{\gamma}|0\rangle\langle 1|$ はどちらも非ユニタリ。

### 9.2 QASM タグ

`KrausChannel` の QASM タグは `"kraus"` だが、QASM への完全なシリアライズは定義されていない（`kraus_channel.py:136`）。  
`params` には最初の Kraus 演算子 $K_0$ が渡されるが、これは QASM ブックキーピング用であり、「チャネル全体のユニタリ表現」ではない（`kraus_channel.py:126–127`）。

---

## 10. 結果の取り出し方

`DMSim.run` は `DensityMatrixJobResult` オブジェクトを含む `Job` を返す（`dmsim.py:226`）。

| 属性 / メソッド | 内容 | コード行 |
|----------------|------|---------|
| `job.result().density_matrix` | `(D, D)` 複素配列 `ρ` | `dmsim.py:179–180` |
| `job.result().state_vector` | `ρ` を flatten した `(1, D²)` 配列 | `dmsim.py:174–175` |
| `job.result().counts` | 空リスト `[]`（測定なし） | `dmsim.py:175` |

**注意**: `state_vector` 属性は既存 `JobResult` API との互換性のためだけに存在する。  
内容は **状態ベクトルではなく密度行列の flatten** である。

---

## 11. DMSim が実装していないこと（コードから確認）

以下の機能はコードに存在しない。

| 機能 | 状況 |
|------|------|
| GPU 演算 | コードに `cuda`・`cupy`・GPU 関連コードなし |
| スパース表現 | ρ は常に密（dense）な `(D,D)` 配列 |
| 適用後 CPTP 検証 | `apply_unitary_to_density` / `apply_kraus_to_density` は出力 ρ を検証しない |
| 測定・mid-circuit measurement | 測定命令を処理するコードなし |
| 状態ベクトルとの互換 API | `state_vector` は ρ の flatten（状態ベクトルではない） |
| 並列実行 | `for` ループによる逐次処理のみ |
| 初期 ρ の物理的検証 | shape チェックのみ |

---

## 12. 数値例：1-qutrit 位相ゲート

N=1、$d_0=3$（$D=3$）の系で位相ゲート $U = \mathrm{diag}(1, e^{i\phi}, e^{2i\phi})$ を適用する場合:

**初期状態**（デフォルト）:

$$
\rho_0 = |0\rangle\langle 0| = \begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix}
$$

**テンソル形式**: shape `(3, 3)` のまま（N=1 なので reshape は `(3,3)` → `(3,3)` = 恒等）

**行レグへの作用**: `_apply_local_op_one_side(rho_t, U, [0], [3], 'row')`

$$
\rho'_{ij} = \sum_\ell U_{i\ell} \rho_{0,\ell j} = U_{i0} \cdot \delta_{j0}
\implies \rho' = U \cdot \rho_0 = \begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix}
$$

**列レグへの作用**: `_apply_local_op_one_side(rho_t, U.conj(), [0], [3], 'col')`

$$
\rho''_{ij} = \sum_\ell \rho'_{i\ell} U^*_{j\ell} = \rho'_{i0} U^*_{j0} = \delta_{i0}\delta_{j0}
\implies \rho'' = \rho_0
$$

これは $U|0\rangle\langle 0|U^\dagger = |0\rangle\langle 0|$ の自明な例（位相は $|0\rangle$ に作用しない）と一致する。

---

## 参考：コード構造の対応表

| 数学的概念 | コード実装 | ファイル:行 |
|-----------|-----------|------------|
| $D = \prod d_i$ | `reduce(operator.mul, dims, 1)` | `dmsim.py:237` |
| $\rho_0 = \|0\cdots 0\rangle\langle 0\cdots 0\|$ | `np.outer(psi0, psi0.conj())` | `dmsim.py:242` |
| テンソル化 $\rho \to \rho_{\rm tensor}$ | `rho.reshape(*dims, *dims)` | `dmsim.py:141,158` |
| $U_S \rho U_S^\dagger$（行レグ） | `tensordot(op_t, rho_t, axes=(range(k,2k), qudits))` | `dmsim.py:121` |
| $U_S \rho U_S^\dagger$（列レグ） | `tensordot(op_t_conj, rho_t, axes=(range(k,2k), [n+q ...]))` | `dmsim.py:121,113` |
| 軸並び替え | `np.transpose(out, perm)` | `dmsim.py:129` |
| $\sum_k K_k \rho K_k^\dagger$ | `apply_kraus_to_density` ループ | `dmsim.py:147–162` |
| $\|\sum K^\dagger K - I\|_F$ | `np.linalg.norm(completeness - identity, ord='fro')` | `kraus_channel.py:112` |
| CPTP 許容値 $10^{-9}$ | `CPTP_TOLERANCE = 1e-9` | `kraus_channel.py:37` |
| $O(D \cdot D_{\rm local}^2)$ | docstring | `dmsim.py:35` |
| $O(r \cdot D \cdot D_{\rm local}^2)$ | docstring | `dmsim.py:36` |
