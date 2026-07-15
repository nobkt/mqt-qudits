# DMsim + Kraus表現：N=4の場合の状態空間はいくつか

> **真実ベースのみ。コード実測値に基づく。**
>
> 根拠ファイル：
> - `src/mqt/qudits/simulation/backends/dmsim.py`（DMSimバックエンド本体）
> - `tutorials/qudit_gksl_simulator.py:256–313`（`_trotter_step_dmsim`）
> - `tutorials/qudit_gksl_simulator.py:284–294`（QuantumRegisterの定義）

---

## 結論（先に述べる）

**YES。`execute_on_backend="dmsim"` で `algorithm="exact_local_channels"` を使うと、Lindblad項が何個あっても、状態空間は系（システム）の次元 $D = d^N$ の密度行列のみ。N=4, d=3なら $D = 3^4 = 81$、密度行列は 81×81。アンシラ（補助ビット）は一切追加されない。**

---

## 1. DMsimバックエンドが確保する行列のサイズ

`dmsim.py:235–237`：

```python
dims: list[int] = list(circuit.dimensions)   # 回路に登録された各量子ビットの次元
d_total = int(reduce(operator.mul, dims, 1)) # D = Π dims[i]
```

回路に登録された量子ビットの次元の積だけを使って $D \times D$ の密度行列を確保する。  
**Lindblad演算子の数・チャンネルの数は `d_total` の計算に一切影響しない。**

## 2. dmsim経路でのQuantumRegisterの定義

`qudit_gksl_simulator.py:289–293`：

```python
n = self.n_system_qudits          # N（分子数）
d = self.params.d                 # 局所次元（例：3）

qreg = QuantumRegister("sys", n, [d] * n)   # システムのみ：n個のd準位量子ビット
circuit = QuantumCircuit(qreg)
```

**`QuantumRegister` の引数は `"sys", n, [d]*n`（システムのみ）。アンシラの `QuantumRegister` は追加されていない。**

N=4, d=3 の場合：
- `dims = [3, 3, 3, 3]`
- `d_total = 3^4 = 81`
- 密度行列 ρ の形状：**81 × 81**

## 3. Lindblad項が複数あっても変わらない理由

DMsim + Kraus経路では、Lindblad各チャンネルは `KrausChannel` 命令として回路に追加される（`qudit_gksl_simulator.py:301, 304`）：

```python
for sites, kraus_ops, kind in self._dmsim_kraus_half:
    target = sites[0] if kind == "single" else list(sites)
    circuit.kraus_channel(target, kraus_ops)
```

`KrausChannel` の作用は（`dmsim.py:147–162`）：

```python
def apply_kraus_to_density(rho, kraus_ops, qudits, dims):
    d_total = int(reduce(operator.mul, dims, 1))
    accumulator = np.zeros((d_total, d_total), dtype=np.complex128)
    for k_op in kraus_ops:
        # rho（D×D）にテンソル縮約で局所的にK_op を作用
        ...
        accumulator += ...
    return accumulator   # D×D のまま
```

$\rho \to \sum_k K_k \rho K_k^\dagger$ を計算するが、$\rho$ のサイズは常に `d_total × d_total = 81 × 81`（N=4のとき）。チャンネルを何個通っても $\rho$ の形状は変わらない。

## 4. Stinespring経路（アンシラあり）との比較

| 比較項目 | Stinespring（algorithm="stinespring"） | Kraus + DMsim（algorithm="exact_local_channels" + dmsim） |
|---|---|---|
| アンシラ qutrit | 1チャンネルにつき1個使用（使用後に部分トレース） | **ゼロ**（追加なし） |
| 各ステップの最大行列サイズ | `d_anc × D = 3 × 81 = 243` 次元（N=4） | **D × D = 81 × 81**（N=4） |
| Lindblad数の影響 | チャンネルを逐次独立適用（ρはD×Dを保持）| チャンネルを回路命令として逐次適用（ρはD×D固定） |
| 状態空間の次元 | システムのみ（アンシラは一時使用後消去）= D | **システムのみ = D** |

注：Stinespring経路でも「アンシラを同時に全部テンソル積する」実装はしていない（mydoc/scipy_expm_26ch_なぜ計算できたか.md 参照）。

## 5. DMsim ドキュメント自身の記述

`dmsim.py:32–34`：

```
For ``D = prod(dimensions)`` the storage is ``16 D²`` bytes.  For the
Scenario 5 GKSL system (``N=4`` qutrits, ``D = 81``) this is
≈ 100 KB: trivially tractable.
```

コードの公式コメントで「N=4 qutrits → D=81 → ~100KB」と明記されている。

## 6. まとめ

| 問い | 答え |
|---|---|
| DMsim + Kraus で状態空間はいくつ？ | $D = d^N$（システムのみ） |
| N=4, d=3 なら？ | $D = 3^4 = 81$、密度行列は **81×81** |
| Lindblad項が26個あると空間は広がるか？ | **広がらない**。チャンネル数は回路命令の数であり、密度行列のサイズに無関係 |
| アンシラビットは使うか？ | **使わない**（QuantumRegister にシステムのみ登録） |

$3^{N+26}$ などという巨大空間になることは**ない**。
