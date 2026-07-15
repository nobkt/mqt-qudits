# N=3・26チャンネルでもscipyで計算できた理由

> **真実ベースのみ。コード実測値に基づく。疎行列かどうか含め、思い込みで書かない。**
>
> 根拠ファイル：
> - `tutorials/stinespring_utils.py`（`stinespring_unitary_from_lindblad`, `build_gksl_superoperator`）
> - `tutorials/classical_gksl_simulator.py`（`ClassicalGKSLSimulator`）
> - `tutorials/qudit_gksl_simulator.py`（`_precompute_unitaries`）

---

## 1. 前提：「26補助ビット」は何を指しているか

`qudit_gksl_simulator.py:67` のコメント：

```
4 qutrits for system + 26 ancilla qudits (d=3) for Lindblad channels.
```

- **26** = Lindblad チャンネルの総数（`len(self.lindblad_ops)`）
- **補助ビット** = 各 Stinespring チャンネルの環境アンシラ（1チャンネルに1アンシラ）

重要：この「26アンシラ」は**同時に1つのヒルベルト空間に連結されていない**。  
各チャンネルの Stinespring ユニタリは**チャンネルごとに独立して**計算される。

---

## 2. Stinespring 計算での行列サイズ（N=3 vs N=4）

`stinespring_utils.py:14–44`:

```python
def stinespring_unitary_from_lindblad(L, dt, d_anc=2):
    d_sys = L.shape[0]
    dim_total = d_anc * d_sys     # ← この行列に expm を呼ぶ
    ...
    U = expm(-1j * theta * G)    # scipy.linalg.expm（密行列）
```

| シミュレーション | system dim `d_sys` | `d_anc` | `dim_total = d_anc × d_sys` | expm 行列サイズ | expm 呼び出し回数 |
|---|---|---|---|---|---|
| N=3, d=3 | 3³ = **27** | 3 | 3 × 27 = **81** | 81 × 81 | 26 |
| N=4, d=3 | 3⁴ = **81** | 3 | 3 × 81 = **243** | 243 × 243 | 26 |

**重要**：26 チャンネル分の Stinespring ユニタリが「全部まとめて 3²⁶ 次元空間」で計算されることは**ない**。  
各チャンネルについて「81×81（N=3）」または「243×243（N=4）」の行列に対して、個別に `expm` を 1 回ずつ呼ぶ。

---

## 3. 疎行列を使っているかどうか

### stinespring_utils.py（Stinespring 計算）

```python
from scipy.linalg import expm   # ← scipy.sparse ではない
```

```python
G = np.zeros((dim_total, dim_total), dtype=np.complex128)   # 密行列
U = expm(-1j * theta * G)                                    # 密行列の expm
```

**疎行列は使っていない。** 81×81 または 243×243 の密行列に対して Padé 近似ベースの `scipy.linalg.expm` を呼ぶ。

### classical_gksl_simulator.py（古典リファレンス）

```python
from scipy.sparse.linalg import expm_multiply   # ← こちらは sparse.linalg
```

```python
self._L_super = build_gksl_superoperator(...)   # 729×729 (N=3) または 6561×6561 (N=4)
vecs = expm_multiply(self._L_super, vec_0, ...)
```

`build_gksl_superoperator` は `np.kron` で密行列として構築する（`stinespring_utils.py:86–131`）。  
`expm_multiply` は疎行列・密行列どちらも受け付けるが、**ここでは密行列を渡している**。  
`expm_multiply` の利点は「`exp(L)` を陽に作らずに `exp(L)·v` を計算する」点であり、疎性の活用ではない。

---

## 4. なぜ計算できたのか（真実）

### 理由1：各 expm の行列サイズが小さい

N=3, d=3 の場合、1チャンネルあたりの Stinespring expm 行列は **81×81**。  
81³ ≈ 53 万要素 の密行列演算であり、現代の PC では 1ms 以下で終わる。  
26 チャンネル分でも 26 ms 程度のオーダー（実測値ではなく概算）。

N=4 の場合でも **243×243**（= 243³ ≈ 1,430 万要素）であり、1 チャンネルあたり数十 ms 程度。

### 理由2：26チャンネルを独立計算するから「全体が 3²⁶ 次元」にならない

「26 アンシラ qutrit が全部同時に存在する空間」は 3²⁶ ≈ 2.5 兆次元になり計算不可能。  
しかし実装では各チャンネルの Stinespring ユニタリを**1チャンネルずつ独立に**計算し、直後に密度行列に作用（部分トレース）して**アンシラを消す**。  
アンシラが系の Hilbert 空間に累積することはない。

`qudit_gksl_simulator.py:226–234`（`_trotter_step`）:

```python
for U_stine_half in self._U_stines_half:
    rho = apply_stinespring_to_density_matrix(rho, U_stine_half, d_anc=self.d_anc)
    # ↑ 各チャンネルで 81×81 行列を rho に作用 → 部分トレース → rho は 27×27 のまま
```

ループのたびに rho は 27×27（N=3）に戻る。26 チャンネル全体を通じて、rho のサイズは変わらない。

### 理由3：古典シミュレータは expm_multiply で exp(L) を陽に作らない

`ClassicalGKSLSimulator` の Liouvillian は 729×729（N=3）。  
`expm_multiply(L, v)` は Krylov 部分空間法で `exp(L)·v` を計算し、exp(L) の 729×729 行列（= 729² ≈ 53 万要素）を RAM に展開しない。  
これが疎性ではなく「行列ベクトル積の繰り返し」による節約。  
ただし L 自体は密行列であり、各行列ベクトル積は密行列演算として実行される。

---

## 5. 疎行列仮説は正しいか

**正しくない。**

| 処理 | 疎行列を使うか | 実際の理由 |
|---|---|---|
| Stinespring expm（qudit_gksl_simulator） | **使わない**（scipy.linalg.expm / 密行列） | 各行列が小さい（81×81, 243×243）から |
| 古典 GKSL（classical_gksl_simulator） | scipy.sparse.linalg.expm_multiply を呼ぶが**入力は密行列** | exp(L) を陽に作らない（Krylov 法）から |
| exact_local_channels | 使わない（scipy.linalg.expm / 密行列） | 局所超演算子が小さい（9×9, 81×81）から |

疎行列の活用によって計算できていた、という解釈は**誤り**。  
計算できた理由は「各 expm 呼び出しが小さな密行列に対して独立に実行されるから」である。

---

## 6. 「26 補助ビット」という言葉が混乱を招く理由

コード中の `n_ancilla_qudits = 26` は「26 個のアンシラ qutrit が同時に系に結合している」という意味ではない。  
それは「26 チャンネルそれぞれが順番に自分のアンシラ（各 1 qutrit）と一時的に結合し、終わったら部分トレースでアンシラを捨てる」実装の累計数である。

全アンシラが同時に存在する Hilbert 空間のサイズを計算すると：

$$\dim = d^{N + 26} = 3^{4+26} = 3^{30} \approx 2.06 \times 10^{14}$$

これは当然計算不可能。しかし実装はそのような空間を一切構築していない。
