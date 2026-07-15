# Stinespring処方のアンシラ数・空間次元・scipy.linalg.expm使用実態

> **真実ベースのみで記述。コード実測値に基づく。推測・願望は一切含まない。**
>
> 根拠ファイル：
> - `tutorials/stinespring_utils.py`
> - `tutorials/gksl_math_utils.py`
> - `tutorials/exact_local_channels.py`
> - `tutorials/qudit_gksl_simulator.py`
> - `tutorials/dmsim_kraus_helpers.py`

---

## 1. N=4 でのアンシラ qutrit 数（Stinespring処方）

### 結論から

**各 Lindblad チャンネルにつき 1 個の qutrit アンシラ（d=3）を使う。**

アンシラを共有・使い回しはしない（コード上は逐次適用のため同一バッファを再利用しているが、
量子回路として考えると各チャンネルごとに新鮮なアンシラが必要）。

### N=4 の Lindblad チャンネル総数

コード（`gksl_math_utils.py:128`）：

```python
expected = 2 * len(params.neighbors) + 5 * N
# N=4: neighbors = [(0,1),(1,2),(2,3)] → 3ペア
# = 2 × 3 + 5 × 4 = 6 + 20 = 26
```

| チャンネル種別 | 式 | 個数 |
|---|---|---|
| TTA チャンネル1（\|2⟩_i\<1\| ⊗ \|0⟩_j\<1\|） | 2×3ペア の前半 | 3 |
| TTA チャンネル2（\|0⟩_i\<1\| ⊗ \|2⟩_j\<1\|） | 2×3ペア の後半 | 3 |
| 蛍光 Γ_fl（\|0⟩_i\<2\|） | 4サイト | 4 |
| 燐光 Γ_ph（\|0⟩_i\<1\|） | 4サイト | 4 |
| 内部転換 k_IC（\|0⟩_i\<2\|） | 4サイト | 4 |
| ISC S→T k_ISC_ST（\|1⟩_i\<2\|） | 4サイト | 4 |
| ISC T→S k_ISC_TS（\|0⟩_i\<1\|） | 4サイト | 4 |
| **合計** | | **26** |

### TTA・蛍光だけのアンシラ数

| 種別 | チャンネル数 | アンシラ qutrit 数 |
|---|---|---|
| TTA（3ペア × 2） | 6 | **6個** |
| 蛍光（4サイト × 1） | 4 | **4個** |
| 合計（TTA+蛍光のみ） | 10 | **10個** |

---

## 2. 空間次元の変化（N=4、d=3）

### システムだけの空間

$$
d^N = 3^4 = 81
$$

システムの密度行列は 81×81。これは **Stinespring の有無にかかわらず変化しない**。

### Stinespring 適用時の「一時的」な拡張

コード（`stinespring_utils.py:35-46`）：

```python
d_sys = L.shape[0]          # Lindblad 演算子の行次元
dim_total = d_anc * d_sys   # アンシラを付加した全空間
G = np.zeros((dim_total, dim_total), dtype=np.complex128)
```

**L は全系演算子（81×81）として `build_lindblad_operators` で構築される。**
したがって：

| パラメータ | 値 |
|---|---|
| `d_sys` | 81（= 3^4） |
| `d_anc` | 3（qutrit アンシラ） |
| `dim_total` | 243（= 81 × 3 = 3^5） |
| Stinespring ユニタリ U のサイズ | **243 × 243** |
| アンシラ付加後の ρ のサイズ | **243 × 243** |

```python
# apply_stinespring_to_density_matrix の実装（stinespring_utils.py:73-76）
env0 = np.zeros((d_anc, d_anc))   # 3×3 アンシラ |0><0|
env0[0, 0] = 1.0
rho_ext = np.kron(env0, rho)       # 3×3 ⊗ 81×81 = 243×243
rho_prime = U @ rho_ext @ U.conj().T  # 243×243 行列積
```

### 重要：拡張は逐次的・一時的

26チャンネルを**同時に**拡張するわけではない。各チャンネルに対して：

1. 81×81 の ρ に 3×3 アンシラを付加 → 243×243 の ρ_ext
2. 243×243 の U を適用
3. アンシラを部分トレースアウト → 81×81 に戻る
4. 次のチャンネルで繰り返す

したがって**同時に必要なメモリは 243×243 の行列 1 枚分**。
26 チャンネル分のアンシラを同時に保持するわけではない。

### 次元整理表（N=4、algorithm="stinespring"）

| 状況 | 空間次元 | サイズ |
|---|---|---|
| システム（平常時） | 3^4 | 81×81 |
| Stinespring 適用中（1 チャンネル） | 3^5 = 3 × 81 | 243×243 |
| Stinespring ユニタリ U | 3^5 × 3^5 | 243×243 |
| Stinespring ジェネレータ G | 3^5 × 3^5 | 243×243 |

---

## 3. exact_local_channels の場合との比較

`algorithm="exact_local_channels"` では Stinespring を**使わない**。

代わりに各チャンネルの**局所リンドブラッド超演算子**を直接指数関数化する：

```python
# exact_local_channels.py: 局所超演算子のサイズ
# 単一サイトチャンネル（蛍光など）: d × d = 3 × 3 の局所 ρ
# → 超演算子は (d²) × (d²) = 9 × 9

# ペアチャンネル（TTA）: d² × d² = 9 × 9 の局所 ρ
# → 超演算子は (d⁴) × (d⁴) = 81 × 81
```

| チャンネル種別 | 局所 ρ サイズ | 局所超演算子 M のサイズ | アンシラ |
|---|---|---|---|
| 蛍光（単一サイト） | 3×3 | 9×9 | **不要** |
| TTA（ペア） | 9×9 | 81×81 | **不要** |

`exact_local_channels` はアンシラ qutrit を一切使わない。代わりに局所超演算子を `scipy.linalg.expm` で指数関数化し、`numpy.einsum` で全系密度行列に作用させる。

---

## 4. scipy.linalg.expm の使用実態

### 核心の真実：アルゴリズム・バックエンドに関係なく expm は必ず使われる

コードでの expm 呼び出し箇所（`qudit_gksl_simulator.py:179`）：

```python
# _precompute_unitaries() — シミュレーション開始時に1回だけ呼ばれる
self._U_H_half = expm(-1j * self.H_total * dt / 2)  # 81×81 ハミルトニアン伝播子
```

これはどの `algorithm` を選んでも、`execute_on_backend` が何であっても**必ず実行される**。

### algorithm="stinespring" の場合

```python
# stinespring_utils.py:43
theta = np.sqrt(dt)
U = expm(-1j * theta * G)  # 243×243 ユニタリ — チャンネル数(26)回呼ばれる
```

expm の呼び出し：
- ハミルトニアン伝播子 U_H_half のため 1 回（81×81 行列）
- 各 Stinespring ユニタリのため 26 回（243×243 行列）
- **合計 27 回、事前計算時に呼ばれる**

### algorithm="exact_local_channels" の場合

```python
# exact_local_channels.py 内の precompute_exact_channels_half()
# 各局所超演算子に対して expm(L_D^local * dt/2) を計算
# 単一サイト: 9×9、ペア: 81×81
```

expm の呼び出し：
- ハミルトニアン伝播子 U_H_half のため 1 回（81×81）
- 各局所チャンネルの超演算子指数化のため 26 回
  - TTA（ペア）6 チャンネル: 81×81 行列の expm × 6
  - 蛍光等（単一）20 チャンネル: 9×9 行列の expm × 20
- **合計 27 回、事前計算時に呼ばれる**

### execute_on_backend="dmsim" の場合

コード（`qudit_gksl_simulator.py:256-311`）の docstring より：

```
no NumPy expm is called inside the simulation loop, only at
the once-per-simulation pre-compute stage (which is a property
of the Hamiltonian Trotter step itself, identical for any
backend choice).
```

DMSim を使っても：

| フェーズ | expm の使用 |
|---|---|
| 事前計算（`_precompute_unitaries`） | **あり**（scipy.linalg.expm × 27 回） |
| シミュレーションループ（各ステップ） | **なし**（MQT-Qudits DMSim が回路を実行） |

**DMSim を使う場合の時間発展の仕組み：**

1. 事前計算時（scipy.linalg.expm を使用）：
   - U_H_half（81×81 ユニタリ）を計算
   - 各 Lindblad チャンネルの Kraus 演算子を計算（Choi-Jamiolkowski 分解経由）
2. 各 Totter ステップ（expm は使わない）：
   - MQT-Qudits `QuantumCircuit` を構築
   - `cu_multi` ゲート（U_H_half）と `KrausChannel` 命令を登録
   - `dmsim.run(circuit, initial_density_matrix=rho)` でバックエンド実行

### まとめ：形式解による時間発展の本質

| 状況 | 時間発展方法 | expm 使用 |
|---|---|---|
| algorithm="stinespring"、backend なし | Stinespring + Trotter（ρ を NumPy で更新） | 事前計算時のみ |
| algorithm="exact_local_channels"、backend なし | 局所超演算子 + Trotter（ρ を NumPy で更新） | 事前計算時のみ |
| algorithm="exact_local_channels"、backend="dmsim" | MQT-Qudits 回路として実行（DMSim） | 事前計算時のみ |

**どの組み合わせでも `scipy.linalg.expm` を使って行列指数関数を計算することには変わりがない。**
違いは「expm の結果をどのように利用するか」である。

---

## 5. Trotterステップ内でのアンシラ付加・除去の流れ（algorithm="stinespring"）

1ステップで行われること（`qudit_gksl_simulator.py:220-250`）：

```
ρ (81×81)
  ↓ U_H_half @ ρ @ U_H_half†  (ハミルトニアン半ステップ)
ρ (81×81)
  ↓ channel 1/26: ρ → 243×243 → U_1 適用 → 部分トレース → 81×81
  ↓ channel 2/26: 同上
  ↓ ...
  ↓ channel 26/26: 同上
  ↓ channel 26/26（逆順） ...
  ↓ channel 1/26（逆順）
ρ (81×81)
  ↓ U_H_half @ ρ @ U_H_half†  (ハミルトニアン半ステップ)
ρ (81×81)
```

パリンドロームチャンネル順序（対称的順序）は Lie-Trotter 積の交換子誤差を消すためだが、
Stinespring 自体が1次近似のため、**全体の収束オーダーは O(dt)（1次）のまま**。

---

## 6. 操作しているコードの制約事項（事実のみ）

| 制約 | 内容 |
|---|---|
| execute_on_backend="dmsim" は algorithm="exact_local_channels" のみ対応 | `stinespring` は Kraus 表現をシステム単体で持てないため |
| state-vector バックエンド（tnsim/misim）は execute_on_backend 非対応 | N≥3 で状態ベクトルの容量制限 |
| DMSim バックエンド経由でも事前計算に expm 不可欠 | 回路パラメータ生成に必要 |
| Stinespring の全系 L（81×81）は `build_lindblad_operators` で構築される | 局所表現ではない全系演算子 |

---

## 7. 数値まとめ（N=4、d=3 の確定値）

| 項目 | 値 |
|---|---|
| システム次元 | 3^4 = **81** |
| 全チャンネル数 | **26**（TTA 6 + 蛍光/燐光/IC/ISC 20） |
| TTA アンシラ qutrit 数 | **6** |
| 蛍光アンシラ qutrit 数 | **4** |
| 全アンシラ qutrit 数（algorithm="stinespring"） | **26** |
| algorithm="exact_local_channels" のアンシラ数 | **0**（アンシラ不要） |
| 1 チャンネル Stinespring 適用中の全空間次元 | 3^5 = **243** |
| Stinespring ユニタリ U のサイズ | **243×243** |
| expm 呼び出し回数（事前計算時） | **27 回**（いずれのアルゴリズム・バックエンドでも） |
