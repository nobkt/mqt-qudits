# DMSim 事前計算（precompute）とは何をやっているのか

> **真実ベースのみ。コード実測値に基づく。推測・願望なし。**
>
> 根拠ファイル：
> - `tutorials/qudit_gksl_simulator.py` (`_precompute_unitaries`, `simulate`)
> - `tutorials/exact_local_channels.py` (`precompute_exact_channels_half`, `build_local_dissipator_super`)
> - `tutorials/stinespring_utils.py` (`stinespring_unitary_from_lindblad`)
> - `tutorials/dmsim_kraus_helpers.py` (`kraus_from_local_superoperator`)

---

## 1. 「事前計算」とはどこで何をするのか

`QuditGKSLSimulator.simulate()` の冒頭（`qudit_gksl_simulator.py:361`）：

```python
dt = t_max / n_steps
self._precompute_unitaries(dt)   # ← これが事前計算
rho = self.prepare_initial_state(initial_state)

for step in range(n_steps):
    rho = self._trotter_step(rho)   # ← ループ内でキャッシュを使う
```

「事前計算」= **`_precompute_unitaries(dt)` の1回呼び出し**。
シミュレーションループ（`n_steps` 回の `_trotter_step`）の**前に1度だけ**実行される。

---

## 2. 事前計算の中で何が起きているか（アルゴリズム別）

### 共通処理（すべてのアルゴリズム・バックエンドで実行される）

```python
# qudit_gksl_simulator.py:179
self._U_H_half = expm(-1j * self.H_total * dt / 2)
```

- **目的**：ハミルトニアン半ステップ伝播子 $U_H = e^{-i H_{total} \cdot dt/2}$ を計算する。
- **行列サイズ**：`H_total` は 81×81（N=4、d=3）。`_U_H_half` も 81×81。
- **計算コスト**：`scipy.linalg.expm`（行列指数関数）を 1 回。
- **必要理由**：各 Trotter ステップで `U_H_half @ rho @ U_H_half†` を 2 回計算するが、`U_H_half` 自体は dt が固定である限り変わらない。毎ステップ `expm` を呼ぶのは無駄なので事前に計算する。

---

### algorithm="stinespring" の場合

```python
# qudit_gksl_simulator.py:181-184
self._U_stines_half = [
    stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
    for L_op, _gamma in self.lindblad_ops  # 26 チャンネル
]
```

各チャンネル α に対して：

```python
# stinespring_utils.py:38-46
dim_total = d_anc * d_sys   # 3 × 81 = 243
G = np.zeros((243, 243), dtype=np.complex128)
G[:81, 81:162] = L.conj().T   # ブロック行列 G を組む
G[81:162, :81] = L
theta = np.sqrt(dt / 2)
U = expm(-1j * theta * G)   # 243×243 行列の指数関数
```

- **目的**：各 Lindblad チャンネルの Stinespring ユニタリ $U_\alpha = e^{-i\sqrt{dt/2} \cdot G_\alpha}$ を計算・キャッシュ。
- **行列サイズ**：1チャンネルあたり 243×243。
- **計算回数**：26 チャンネル分 → `expm` を 26 回（加えてハミルトニアン 1 回 = 合計 27 回）。
- **必要理由**：全 `n_steps` ステップで同じ `U_α` を繰り返し使う。dt が固定なので毎ステップ `expm` を呼ぶ必要はない。

---

### algorithm="exact_local_channels" の場合

```python
# qudit_gksl_simulator.py:189-191
self._exact_channels_half = precompute_exact_channels_half(self.params, dt)
```

内部（`exact_local_channels.py:318-339`）：

```python
for (sites, L_local, kind) in get_local_lindblad_ops(params):
    d_local = L_local.shape[0]
    L_D_local = build_local_dissipator_super(L_local, d_local)   # 局所超演算子
    exp_LD_half = expm(L_D_local * (dt / 2.0))                    # 局所超演算子の指数関数
    out.append((sites, exp_LD_half, kind))
```

各チャンネル α に対して：

1. **`build_local_dissipator_super(L_local, d_local)`**：  
   局所 GKSL 散逸超演算子 $\mathcal{L}_{D}^{local}$ を構築する。  
   - 単一サイト（蛍光等）：$d^2 \times d^2 = 9 \times 9$ 行列  
   - ペア（TTA）：$d^4 \times d^4 = 81 \times 81$ 行列

2. **`expm(L_D_local * (dt / 2))`**：  
   $e^{\mathcal{L}_{D,\alpha}^{local} \cdot dt/2}$ を計算。これが**正確なチャンネルの半ステップ**。  
   - 単一サイト：9×9 行列の `expm`（20 チャンネル）  
   - ペア：81×81 行列の `expm`（6 チャンネル）

- **合計**：ハミルトニアン 1 回 + 各チャンネル 26 回 = 27 回（同じ）。
- **必要理由**：各ステップで `expm` を呼ぶ代わりに、結果をキャッシュしてから全ステップで再利用。

---

### execute_on_backend="dmsim" の場合（追加処理）

`exact_local_channels` の事前計算に加えて：

```python
# qudit_gksl_simulator.py:193-204
for sites, M_half, kind in self._exact_channels_half:
    d_root = self.params.d if kind == "single" else self.params.d ** 2
    kraus = kraus_from_local_superoperator(M_half, d_root)
    kraus_list.append((sites, kraus, kind))
self._dmsim_kraus_half = kraus_list

self._dmsim_backend = MQTQuditProvider().get_backend("dmsim")
```

各チャンネルの局所超演算子 $M = e^{\mathcal{L}_{D,\alpha}^{local} \cdot dt/2}$（9×9 または 81×81）を  
**Kraus 演算子リストに変換**する。

変換手順（`dmsim_kraus_helpers.py:51-130`）：

1. **超演算子 → Choi 行列**：  
   $M$ のテンソル成分 $T[i,j,k,l] = M[i + d\cdot j, k + d\cdot l]$ を  
   添字転置して $C[i\cdot d + k,\ j\cdot d + l] = T[i,j,k,l]$ を構成。  
   Choi 行列 $C$ は 9×9（単一サイト）または 81×81（ペア）。

2. **Choi 行列の固有値分解**：  
   CPTP マップの Choi 行列は正半定値。  
   $C = \sum_\alpha \lambda_\alpha |v_\alpha\rangle\langle v_\alpha|$  
   ($\lambda_\alpha < 0$ で絶対値が `1e-12` を超えると `ValueError` を raise。ヒューリスティック補正なし。)

3. **Kraus 演算子の抽出**：  
   各固有値 $\lambda_\alpha > 0$ に対して  
   $K_\alpha = \sqrt{\lambda_\alpha} \cdot \text{reshape}(v_\alpha, (d, d))$

- **必要理由**：MQT-Qudits の `KrausChannel` 命令は Kraus 演算子リスト `{K_α}` を受け取る。超演算子形式では渡せない。各ステップで回路として発行するためにキャッシュが必要。

---

## 3. 事前計算が必要な理由（まとめ）

事前計算は**以下の 2 つの問題**を解決するために存在する：

### 問題1：`expm` は計算コストが高い

`scipy.linalg.expm` は Padé 近似ベースのアルゴリズムであり、  
$n \times n$ 行列に対して $O(n^3)$ のコストがかかる。

N=4 で 26 チャンネルの場合：

| 計算対象 | 行列サイズ | `expm` 回数 |
|---|---|---|
| $U_{H,half}$（ハミルトニアン） | 81×81 | 1 |
| stinespring: $U_\alpha$（各チャンネル） | 243×243 | 26 |
| exact_local: $e^{\mathcal{L}_{D,\alpha}^{local} \cdot dt/2}$（ペア） | 81×81 | 6 |
| exact_local: $e^{\mathcal{L}_{D,\alpha}^{local} \cdot dt/2}$（単一） | 9×9 | 20 |

`dt` が固定である限り、これらは全シミュレーション期間を通じて**変化しない**。  
`n_steps = 1000` であれば、事前計算がなければ `expm` を 27,000 回呼ぶことになる。  
事前計算によって 27 回で済む。

### 問題2：DMSim バックエンドに渡す形式が超演算子ではなく Kraus 形式

`KrausChannel` 命令は `{K_α}` リストを受け取るが、`exact_local_channels` が計算するのは  
超演算子 $M$（行列として）である。  
両者の間の変換（Choi-Jamiolkowski 分解）も dt が固定なら結果は変わらないため、  
事前計算で一度変換してキャッシュしておく。

---

## 4. 事前計算のタイミングと dt 依存性

```
simulate(t_max, n_steps) が呼ばれる
  ↓
dt = t_max / n_steps   ← dt が決まる
  ↓
_precompute_unitaries(dt) ← dt に依存した行列を全て計算・キャッシュ
  ↓
for step in range(n_steps):
    _trotter_step(rho)   ← キャッシュを使うだけ、expm は呼ばない
```

`n_steps` を変えて `simulate` を再実行すると、dt が変わるため事前計算も再実行される。  
同じ `dt` で複数回 `simulate` を呼ぶ場合でも、現在の実装は毎回 `_precompute_unitaries` を呼ぶ  
（キャッシュの無効化チェックなし）。

---

## 5. DMSim バックエンドを使う場合の事前計算のフロー（全体像）

```
dt が決まる
  ↓
expm(-i * H_total * dt/2)
  → _U_H_half (81×81) キャッシュ
  
各 Lindblad チャンネル α (26個) に対して:
  ↓
  build_local_dissipator_super(L_local, d_local)
    → L_D_local (9×9 または 81×81)
  ↓
  expm(L_D_local * dt/2)
    → exp_LD_half (9×9 または 81×81) キャッシュ
  ↓
  kraus_from_local_superoperator(exp_LD_half, d_local)
    → Choi 行列構築 (9×9 または 81×81)
    → 固有値分解
    → K_α = sqrt(λ_α) * reshape(v_α)
    → Kraus リスト キャッシュ

DMSim バックエンドの取得:
  MQTQuditProvider().get_backend("dmsim")

以上で事前計算完了。
```

以降のシミュレーションループでは：
- `_U_H_half` を `cu_multi` ゲートとして回路に追加
- 各 Kraus リストを `KrausChannel` 命令として回路に追加
- `dmsim.run(circuit, initial_density_matrix=rho)` を呼び出す

---

## 6. 事前計算で呼ばれる `expm` の総数（確定値）

| algorithm | execute_on_backend | `expm` 呼び出し総数 | 内訳 |
|---|---|---|---|
| stinespring | None | **27** | 1（H） + 26（Stinespring U） |
| exact_local_channels | None | **27** | 1（H） + 6（pair 81×81） + 20（single 9×9） |
| exact_local_channels | dmsim | **27** | 同上（Choi 分解は expm でなく固有値分解） |

Choi-Jamiolkowski 変換は **`expm` を使わない**（固有値分解 `np.linalg.eigh`）。
