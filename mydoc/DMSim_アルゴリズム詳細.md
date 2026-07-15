# DMSim アルゴリズム詳細

> **本文書はコードの実装を直接引用して作成した。推測・希望的記述は含まない。**
> 引用元: `src/mqt/qudits/simulation/backends/dmsim.py`（288行）、
> `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`

---

## 1. DMSim とは何か

`DMSim`（`mqt.qudits.simulation.backends.DMSim`）は **密度行列シミュレータ**である。
状態ベクトル `|ψ⟩` を追跡するのではなく、密度行列 `ρ`（`D×D` 複素行列、`D = Π_i d_i`）を直接更新する。

状態ベクトルバックエンド `TNSim` / `MISim` とは別物であり、互換性はない（`JobResult.state_vector` には `ρ` の flatten が格納される）。

---

## 2. サポートする命令の種類

DMSim が認識する命令は **2種類のみ**。それ以外の命令は `to_matrix(identities=0)` でユニタリ行列を取り出してユニタリゲートとして処理される。

| 命令の種類 | 判定条件 | 更新則 |
|-----------|---------|-------|
| ユニタリゲート | `isinstance(instruction, KrausChannel)` が `False` | `ρ → U_S ρ U_S†` |
| Kraus チャネル | `isinstance(instruction, KrausChannel)` が `True` | `ρ → Σ_k K_k ρ K_k†` |

コード該当箇所（`dmsim.py:260–288`）:

```python
if isinstance(instruction, KrausChannel):
    return apply_kraus_to_density(rho, instruction.kraus_operators, qudits, dims)

u_matrix = instruction.to_matrix(identities=0)
return apply_unitary_to_density(rho, u_matrix, qudits, dims)
```

---

## 3. 密度行列の初期化

`dmsim.py:239–252`

```python
if initial_density_matrix is None:
    psi0 = np.zeros(d_total, dtype=np.complex128)
    psi0[0] = 1.0
    rho = np.outer(psi0, psi0.conj())   # = |0…0⟩⟨0…0|
else:
    rho = np.asarray(initial_density_matrix, dtype=np.complex128)
    if rho.shape != (d_total, d_total):
        raise ValueError(...)
    rho = rho.copy()
```

- デフォルトは `|0…0⟩⟨0…0|`（純粋状態）
- 任意の初期密度行列を `initial_density_matrix` オプションで渡せる
- **CPTP 検証・Hermiticity チェックは一切行わない**。物理的に正当な行列を渡す責任は呼び出し元にある

---

## 4. メインループ

`dmsim.py:254–256`

```python
for instruction in circuit.instructions:
    rho = self._apply_instruction(rho, instruction, dims, n)
```

`circuit.instructions` は `list[Gate]` であり、**回路に追加された順番に逐次適用**される。
並列化・並べ替えは行われない。

---

## 5. テンソル演算の実装

### 5.1 データ形状の変換

`ρ` は `(D, D)` の行列として保持されるが、内部演算では **テンソル形状** `(*dims, *dims)` に reshape して扱う。

- 前半 `N` 軸: 行インデックス（「行レグ」）
- 後半 `N` 軸: 列インデックス（「列レグ」）

例: qudit 次元が `[3, 3, 3, 3]`（4 qutrit）の場合、`ρ` は `(3,3,3,3,3,3,3,3)` テンソルとして扱われる。

### 5.2 局所演算子の片側作用: `_apply_local_op_one_side`

`dmsim.py:75–129`

```python
op_t = op_matrix.reshape(*local_dims, *local_dims)
out = np.tensordot(op_t, rho_tensor,
                   axes=(list(range(k, 2*k)), contract_axes))
```

- `side='row'` のとき: `contract_axes = list(qudits)`（行レグに作用 → `ρ_new = (op ⊗ I) ρ`）
- `side='col'` のとき: `contract_axes = [n+q for q in qudits]`（列レグに作用 → `ρ_new = ρ (op ⊗ I)†`）
- `tensordot` 後に軸の順序を復元するため `np.transpose` で置換する（`label_to_pos` 辞書による逆置換）

### 5.3 ユニタリゲートの適用: `apply_unitary_to_density`

`dmsim.py:132–144`

```python
rho_t = rho.reshape(*([d for d in dims] * 2))
rho_t = _apply_local_op_one_side(rho_t, u_matrix, qudits, dims, "row")
rho_t = _apply_local_op_one_side(rho_t, u_matrix.conj(), qudits, dims, "col")
return rho_t.reshape(d_total, d_total)
```

実質: `ρ → U_S ρ U_S†`  
- 行レグに `U` を作用させ、列レグに `U*`（複素共役、転置なし）を作用させる
- 転置なしの複素共役を列レグに作用させることで `U_S†` と等価になる（列方向に作用するため）

### 5.4 Kraus チャネルの適用: `apply_kraus_to_density`

`dmsim.py:147–162`

```python
accumulator = np.zeros((d_total, d_total), dtype=np.complex128)
for k_op in kraus_ops:
    rho_t = rho.reshape(*([d for d in dims] * 2))
    rho_t = _apply_local_op_one_side(rho_t, k_op, qudits, dims, "row")
    rho_t = _apply_local_op_one_side(rho_t, k_op.conj(), qudits, dims, "col")
    accumulator += rho_t.reshape(d_total, d_total)
return accumulator
```

実質: `ρ → Σ_k K_k ρ K_k†`  
各 Kraus 演算子について `ρ_k = K_k ρ K_k†` を計算してゼロ行列に足し込む。

---

## 6. KrausChannel 命令の構造

`src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`

### 6.1 CPTP 検証

コンストラクタで `Σ_k K_k† K_k` を計算し、単位行列からのフロベニウスノルム差を検証する:

```python
CPTP_TOLERANCE = 1e-9

completeness = sum(K.conj().T @ K for K in kraus_list)
defect = np.linalg.norm(completeness - identity, ord='fro')
if defect > CPTP_TOLERANCE:
    raise ValueError(...)
```

- 許容値: `1e-9`（フロベニウスノルム）
- ヒューリスティック補正は行わない。違反した場合は `ValueError` を送出して終了
- 検証は**構築時に1回だけ**行われる

### 6.2 `to_matrix` の動作

```python
def to_matrix(self, identities: int = 0):
    raise NotImplementedError(
        "KrausChannel has no unitary matrix representation. "
        "Use a density-matrix backend (e.g. DMSim) instead."
    )
```

`KrausChannel` にはユニタリ行列表現が存在しない。`DMSim` 以外のバックエンドに渡した場合、`to_matrix` 呼び出しで `NotImplementedError` が発生する。

### 6.3 Kraus 演算子の取り出し

`DMSim._apply_instruction` は `instruction.kraus_operators` プロパティ（防御的コピーなし）で演算子リストを取り出す。

---

## 7. `noise_model` オプションの扱い

`dmsim.py:219–224`

```python
if self._options.get("noise_model", None) is not None:
    raise ValueError(
        "DMSim does not accept a NoiseModel option. "
        "Add noise as KrausChannel instructions in the circuit instead."
    )
```

- `noise_model` を渡すと**即座に `ValueError`** を送出する
- これは設計上意図的な制約である（docstring に明記）
- ノイズを扱う唯一の方法は、回路に `KrausChannel` 命令を明示的に追加することである

---

## 8. 計算コスト

`dmsim.py` の module docstring より:

| 演算 | コスト |
|------|--------|
| ユニタリゲート（局所） | `O(D · D_local²)` |
| Kraus チャネル（局所） | `O(rank · D · D_local²)` |

ここで:
- `D = Π_i d_i`（全系の次元）
- `D_local = Π_{q ∈ target} d_q`（対象 qudit の局所次元）
- `rank` = Kraus 演算子の個数

メモリ使用量: `16 D²` バイト（`complex128` = 16 バイト/要素）  
例: `N=4` qutrit（`D=81`）→ `≈ 100 KB`

---

## 9. `execute` と `run` の関係

```
run(circuit, **options)
  ├── noise_model チェック（ValueError or スキップ）
  └── execute(circuit, initial_density_matrix=...)
        └── for instruction in circuit.instructions:
              rho = _apply_instruction(rho, instruction, dims, n)
        return rho
```

- `run` は `Job` を生成して結果をセットして返す（非同期実行なし）
- `execute` はループを回して最終 `ρ` を返す numpy 配列
- `DensityMatrixJobResult` は `ρ.reshape(1, -1)` を `state_vector` 属性に持ち、元の `ρ` を `density_matrix` プロパティで返す

---

## 10. 他のバックエンドとの対比

| 項目 | DMSim | TNSim / MISim |
|------|-------|---------------|
| 内部状態 | 密度行列 `ρ`（D×D） | 状態ベクトル `|ψ⟩`（D） |
| 混合状態 | 扱える | 扱えない |
| KrausChannel | 対応 | `to_matrix()` が `NotImplementedError` → 実行不可 |
| NoiseModel | 受け付けない（ValueError） | 受け付ける |
| メモリ | `O(D²)` | `O(D)` |

---

## 11. 事実として確認できないこと

以下は本コードから読み取れない（実装が存在しない）:

- GPU / 並列化による高速化
- スパース行列最適化
- Trotter 分解の内部実装（DMSim 自体は行わない。Trotter ステップは外部で組み立てた回路の命令列として渡される）
- CPTP 保存の積極的検証（適用後の `ρ` の trace や Hermiticity を検証するコードは存在しない）
