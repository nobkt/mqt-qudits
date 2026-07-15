# なぜ Kraus 演算子は DMSim が必要で TNSim へ適用できないのか

> **方針：** このリポジトリの実コードを根拠とし、事実のみを記述する。
> 推測・希望的観測・曖昧な記述は一切含めない。

---

## 1. 結論（1 行）

**TNSim は純粋状態ベクトル（complex 配列 1 本）を保持・更新するバックエンドであり、
`to_matrix()` が行列を返すことを前提に動作する。`KrausChannel` は `to_matrix()` を実装しておらず
`NotImplementedError` を送出するため、TNSim が `KrausChannel` を含む回路を実行しようとすると
必ずエラーになる。DMSim は密度行列 $\rho$（$D \times D$ 行列）を保持し、
`KrausChannel` を特別ケースとして識別・処理する唯一のバックエンドである。**

---

## 2. 各バックエンドが保持する状態の型

### 2.1 TNSim（テンソルネットワーク状態ベクトルシミュレータ）

ファイル：`src/mqt/qudits/simulation/backends/tnsim.py`

TNSim が保持するのは**純粋状態ベクトル** $|\psi\rangle \in \mathbb{C}^D$（$D = \prod_i d_i$）である。
内部的にはテンソルネットワーク（`tensornetwork` ライブラリ）で表現され、
最終的に `D` 次元の複素ベクトル 1 本として返す：

```python
result = self.__contract_circuit(self.system_sizes, self.circ_operations)
result = np.transpose(result.tensor, list(range(len(self.system_sizes))))
state_size = reduce(operator.mul, self.system_sizes, 1)
return result.reshape(1, state_size)   # ← 1×D の状態ベクトル
```

`JobResult.state_vector` は $1 \times D$ の配列である。

### 2.2 DMSim（密度行列シミュレータ）

ファイル：`src/mqt/qudits/simulation/backends/dmsim.py`

DMSim が保持するのは**密度行列** $\rho \in \mathbb{C}^{D \times D}$（$D = \prod_i d_i$）である：

```python
psi0 = np.zeros(d_total, dtype=np.complex128)
psi0[0] = 1.0
rho = np.outer(psi0, psi0.conj())   # ← D×D 行列
```

`JobResult.density_matrix` は $D \times D$ の配列である。
`JobResult.state_vector` は後方互換のために平坦化した $\rho$ を返すが、
これは状態ベクトルではなく密度行列の reshape である（docstring に明記）。

---

## 3. TNSim がゲートを処理する仕組み

`tnsim.py:__contract_circuit` より：

```python
for op in operations:
    op_matrix = op.to_matrix(identities=1)   # ← ここで必ず to_matrix() を呼ぶ
    lines = op.reference_lines

    if op.gate_type == GateTypes.SINGLE:
        op_matrix = op_matrix.T
        ...
    elif op.gate_type == GateTypes.TWO and not op.is_long_range:
        op_matrix = op_matrix.T
        op_matrix = op_matrix.reshape(...)
        ...
    elif op.is_long_range or op.gate_type == GateTypes.MULTI:
        op_matrix = op_matrix.T
        op_matrix = op_matrix.reshape(...)
        ...
    self.__apply_gate(qudits_legs, op_matrix, lines)
```

**すべてのゲートに対して `to_matrix()` を無条件に呼び出す。** TNSim はゲートの型を `KrausChannel` かどうかで分岐しない。

---

## 4. `KrausChannel.to_matrix()` の動作

ファイル：`src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`

```python
def to_matrix(self, identities: int = 0) -> NDArray[np.complex128]:
    msg = (
        "KrausChannel has no unitary matrix representation. "
        "Use a density-matrix backend (e.g. DMSim) instead."
    )
    raise NotImplementedError(msg)
```

`__array__` も同様に `NotImplementedError` を送出する：

```python
def __array__(self) -> NDArray[np.complex128]:
    msg = (
        "KrausChannel has no unitary matrix representation. "
        "Use a density-matrix backend (e.g. DMSim) instead."
    )
    raise NotImplementedError(msg)
```

**`KrausChannel` が `to_matrix()` を実装していない理由は数学的必然性**：
CPTP 写像 $\mathcal{E}(\rho) = \sum_k K_k \rho K_k^\dagger$ は、
`Kraus rank ≥ 2`（ノイズチャンネルの場合は常にそう）の場合、
系の純粋状態を必ず混合状態に変換する。
すなわち「$\mathcal{E}$」を実現する単一のユニタリー行列は（系の次元だけでは）存在しない。
行列表現は定義できないため `to_matrix()` は実装不可能である。

---

## 5. DMSim が `KrausChannel` を処理する仕組み

`dmsim.py:_apply_instruction` より：

```python
@staticmethod
def _apply_instruction(rho, instruction, dims, n_qudits):
    target = instruction.target_qudits
    qudits = (target,) if isinstance(target, int) else tuple(target)

    if isinstance(instruction, KrausChannel):
        # ← KrausChannel を識別し、to_matrix() を呼ばずに処理する
        return apply_kraus_to_density(rho, instruction.kraus_operators, qudits, dims)

    # ユニタリーゲートパス
    u_matrix = instruction.to_matrix(identities=0)
    return apply_unitary_to_density(rho, u_matrix, qudits, dims)
```

`isinstance(instruction, KrausChannel)` による分岐が存在するため、
`to_matrix()` を呼ばずに `instruction.kraus_operators`（$\{K_k\}$ のリスト）を直接取得できる。

`apply_kraus_to_density` の実装：

```python
def apply_kraus_to_density(rho, kraus_ops, qudits, dims):
    accumulator = np.zeros((d_total, d_total), dtype=np.complex128)
    for k_op in kraus_ops:
        rho_t = rho.reshape(...)
        rho_t = _apply_local_op_one_side(rho_t, k_op, qudits, dims, "row")   # K_k ρ
        rho_t = _apply_local_op_one_side(rho_t, k_op.conj(), qudits, dims, "col")  # K_k ρ K_k†
        accumulator += rho_t.reshape(d_total, d_total)
    return accumulator   # Σ_k K_k ρ K_k†
```

これは密度行列 $\rho$ が存在するから成立する。純粋状態ベクトル $|\psi\rangle$ に対しては
$\sum_k K_k |\psi\rangle \langle\psi| K_k^\dagger$ を状態ベクトル 1 本で表現できない（混合状態になるため）。

---

## 6. TNSim の「ノイズ」機能との混同について

TNSim は `noise_model` オプションと `stochastic_simulation` を受け付ける
（`stochastic_sim.py`）。しかしこれは `KrausChannel` とは**全く異なる仕組み**である。

### 6.1 TNSim のノイズ：確率的射影（Monte Carlo wave-function 法の一種）

`stochastic_sim.py:stochastic_simulation`:

```python
factory = NoisyCircuitFactory(noise_model, circuit)
args_tn = [(backend, factory) for _ in range(shots)]
results = pool.map(stochastic_execution_tn, args_tn)
```

`NoisyCircuitFactory.generate_circuit()` が行うのは：
1. 元の回路の各ゲートをコピーする
2. ノイズモデルに従って**確率的に** `NoiseX`（`X` ゲート）/ `NoiseY`（`Y` ゲート）/ `NoiseZ`（`Z` ゲート）を追加する
3. 追加されたゲートは全て `to_matrix()` が実装されたユニタリーゲートである

```python
# noisy_circuit_factory.py
noisy_circuit.noisex(dit)  # ← NoiseX は to_matrix() を持つ
noisy_circuit.noisez(dit)  # ← NoiseZ は to_matrix() を持つ
```

つまり TNSim のノイズは「ランダムなユニタリーを射影的に適用するショット」を `shots` 回繰り返し、
その平均として混合状態をサンプリングする。各ショットは純粋状態のままである。

### 6.2 KrausChannel のノイズ：決定論的 CPTP 写像

DMSim の `KrausChannel` は 1 回の実行で $\rho \to \sum_k K_k \rho K_k^\dagger$ を**決定論的**に実行する。
これは確率的サンプリングではなく、混合状態の密度行列を直接更新する。

| | TNSim + NoiseModel | DMSim + KrausChannel |
|--|---|---|
| 保持する状態 | 純粋状態ベクトル $|\psi\rangle$ | 密度行列 $\rho$ |
| ノイズの実現 | ランダムユニタリーを確率的に射影（Monte Carlo） | CPTP 写像を決定論的に $\rho$ に適用 |
| ショット数 | 複数ショットで平均を取る必要がある | 1 回で期待値を得る |
| `KrausChannel` 命令 | **使用不可**（`to_matrix()` が `NotImplementedError`） | **使用可能**（`isinstance` 分岐で識別） |
| ノイズに必要な API | `NoiseModel`（`Noise`/`SubspaceNoise`） | `KrausChannel` 命令を回路に追加 |

---

## 7. 技術的な数学的理由

### 7.1 状態ベクトルは混合状態を表現できない

純粋状態 $|\psi\rangle$ は任意の $U$ に対して $U|\psi\rangle$ として変換されるが、
CPTP 写像（Kraus rank ≥ 2）の出力は

$$
\mathcal{E}(|\psi\rangle\langle\psi|)
= \sum_{k} K_k |\psi\rangle\langle\psi| K_k^\dagger
$$

これは一般に**混合状態**であり、単一のベクトル $|\phi\rangle$ で $|\phi\rangle\langle\phi|$ とは書けない。
TNSim は $|\psi'\rangle$ を 1 本のベクトルとして返さなければならないが、
混合状態の $\rho'$ を 1 本のベクトルで表現する方法は存在しない。

### 7.2 Kraus 演算子に「等価な単一ユニタリー」は存在しない（系単独では）

$\mathcal{E}$ が Kraus rank $r \geq 2$ の CPTP 写像であるとき、
系（次元 $D$）だけの空間に $\mathcal{E}(\rho) = U \rho U^\dagger$ となる $U$ は存在しない。
Stinespring 定理により $U$ は系＋アンシラの拡大空間（次元 $r \cdot D$）に存在するが、
TNSim はアンシラを扱う機能を持たず、partial trace 命令も存在しない。

---

## 8. DMSim も「状態ベクトル」を返すか？

DMSim の `JobResult.state_vector` は `dmsim.py:175`:

```python
def __init__(self, density_matrix: NDArray[np.complex128]) -> None:
    flat = density_matrix.reshape(1, -1)   # D² 要素の配列
    super().__init__(state_vector=flat, counts=[])
```

これは **$\rho$ を平坦化したもの**であり、純粋状態ベクトルではない。
`JobResult` API との後方互換のためにこの名前になっているが、
物理的な状態ベクトル $|\psi\rangle$ とは別物である。
`JobResult.density_matrix` または `JobResult.get_density_matrix()` が正しいアクセス法。

---

## 9. qudit_gksl_simulator.py の明示的なエラーメッセージ

`tutorials/qudit_gksl_simulator.py:109-117`:

```python
if execute_on_backend != "dmsim":
    msg = (
        "execute_on_backend currently only supports 'dmsim' "
        "(MQT-Qudits density-matrix backend); got "
        f"{execute_on_backend!r}.  state-vector backends "
        "(tnsim/misim) cannot run the fresh-ancilla per-step "
        "circuit at N>=3 due to state-vector capacity limits "
        "(see STATUS_HONEST_2026-05.md A-1)."
    )
    raise ValueError(msg)
```

ここで言う「state-vector capacity limits」とは、もし Stinespring 処方を使う場合にアンシラを回路に含めると
状態ベクトルの次元が $d^{N_\mathrm{sys} + N_\mathrm{anc}}$ に膨張し TNSim では扱えなくなるという問題であり、
これは `KrausChannel` の問題とは独立した別の理由である（ただし根は同じ：TNSim は混合状態を扱えない）。

---

## 10. コード参照箇所一覧

| 内容 | ファイル | 行（概略） |
|------|---------|-----------|
| TNSim の状態ベクトル（$1 \times D$）返却 | `src/mqt/qudits/simulation/backends/tnsim.py` | 69–74 |
| TNSim が `to_matrix()` を無条件呼び出し | `src/mqt/qudits/simulation/backends/tnsim.py` | 80 |
| `KrausChannel.to_matrix()` の `NotImplementedError` | `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py` | 147–155 |
| `KrausChannel.__array__()` の `NotImplementedError` | `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py` | 139–145 |
| DMSim が `KrausChannel` を `isinstance` で識別 | `src/mqt/qudits/simulation/backends/dmsim.py` | 272–274 |
| `apply_kraus_to_density`（$\Sigma_k K_k \rho K_k^\dagger$） | `src/mqt/qudits/simulation/backends/dmsim.py` | 121–132 |
| DMSim の初期密度行列（$D \times D$ 行列） | `src/mqt/qudits/simulation/backends/dmsim.py` | 236–243 |
| DMSim docstring の「state-vector でない」注意書き | `src/mqt/qudits/simulation/backends/dmsim.py` | 38–42 |
| TNSim の確率的ノイズ（Monte Carlo ユニタリー） | `src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py` | 98–116 |
| NoiseX/NoiseZ（TNSim ノイズ、ユニタリーゲート） | `src/mqt/qudits/quantum_circuit/gates/noise_x.py` | 38–55 |
| tnsim/misim 制限の明示的エラー | `tutorials/qudit_gksl_simulator.py` | 109–117 |
