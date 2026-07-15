# TNSim はアンシラを「扱えない」とは具体的に何を指すか

> **本文書はコード実装の直接引用に基づく。推測・希望的記述は含まない。**  
> 引用元ファイル:
> - `src/mqt/qudits/simulation/backends/tnsim.py`（全121行）
> - `tutorials/qudit_gksl_circuit_simulator.py`（モジュール docstring 他）
> - `STATUS_HONEST_2026-05.md`（A-1 項）

---

## 1. 前提：Stinespring 拡張と「アンシラ」

Lindblad GKSL 方程式の Trotter ステップを量子回路で実現する標準的な手法は **Stinespring 拡張** である。

Lindblad 演算子 $L_\alpha$（$d_{\rm sys} \times d_{\rm sys}$）に対し、Stinespring ユニタリ $U_\alpha$（$(d_{\rm anc} \cdot d_{\rm sys}) \times (d_{\rm anc} \cdot d_{\rm sys})$）を構築し、

$$
\mathcal{E}_\alpha(\rho_{\rm sys}) = \mathrm{tr}_{\rm anc}\!\bigl(U_\alpha \,(\rho_{\rm sys} \otimes |0\rangle\langle 0|_{\rm anc})\, U_\alpha^\dagger\bigr)
$$

として散逸チャネルを実現する。ここで「アンシラ（ancilla）」とは、この拡張に必要な補助 qudit のことである。

このリポジトリでは $d_{\rm anc} = d_{\rm sys} = 3$（qutrit）として実装されている（`qudit_gksl_circuit_simulator.py:94`: `self.d_anc = params.d`）。

---

## 2. TNSim の実行モデル

TNSim（`src/mqt/qudits/simulation/backends/tnsim.py`）は **純粋状態ベクトル** $|\psi\rangle$ をテンソルネットワーク上で伝播させるシミュレータである。

### 2.1 初期状態

```python
# tnsim.py:83-87
for s in system_sizes:
    z = [0] * s
    z[0] = 1
    state_nodes.append(tn.Node(np.array(z, dtype="complex")))
```

全 qudit を $|0\rangle$ に初期化し、各 qudit を独立なテンソルノードとして張る。

### 2.2 ゲート適用

```python
# tnsim.py:91
op_matrix = op.to_matrix(identities=1)
```

各命令について `to_matrix(identities=1)` でユニタリ行列を取り出し、対象 qudit の「テンソルレグ」に縮約する。

### 2.3 最終縮約

```python
# tnsim.py:120
return tn.contractors.auto(all_nodes, output_edge_order=qudits_legs)
```

全ノードを縮約して状態ベクトルを返す。

**TNSim の処理パイプラインは以上であり、これ以外の操作は一切存在しない。**

---

## 3. アンシラを「扱えない」ことの具体的内容

「TNSim はアンシラを扱えない」という主張は、以下の **3 つの独立した欠如** を意味する。それぞれコードで確認できる。

### 3.1 欠如①：KrausChannel 命令を認識しない

Stinespring 拡張後のアンシラを「トレースアウト」する操作は、状態ベクトルシミュレータでは **部分トレース（partial trace）** に相当し、非ユニタリ操作である。これは `KrausChannel` 命令で表現されるが：

- `KrausChannel.to_matrix()` は `NotImplementedError` を送出する（`kraus_channel.py:154–159`）
- TNSim は `op.to_matrix(identities=1)` を全命令に無条件で呼ぶ（`tnsim.py:91`）
- したがって TNSim の命令ループに `KrausChannel` が渡ると **即座に `NotImplementedError` で停止する**

これは DMSim が `isinstance(instruction, KrausChannel)` で分岐するのと対照的である（`dmsim.py:272`）。

### 3.2 欠如②：mid-circuit リセット命令が存在しない

Stinespring チャネルを複数回（チャネル数 $M$ 回）適用するには、チャネルごとにアンシラを $|0\rangle$ に戻す必要がある:

$$
\text{1 Trotter ステップ} = U_{H}(\Delta t/2) \circ \underbrace{\mathcal{E}_1 \circ \mathcal{E}_2 \circ \cdots \circ \mathcal{E}_M}_{\text{各チャネルの前後でアンシラを }\,|0\rangle\,\text{ にリセット}} \circ U_{H}(\Delta t/2)
$$

「アンシラリセット」とは回路の途中で特定の qudit を $|0\rangle$ に射影・初期化する操作である。

**コード確認**:
```
grep -n "reset\|measure" tnsim.py → 0件
grep -n "reset\|measure" misim.py → 0件
```

TNSim の命令ループ（`tnsim.py:90–118`）は `GateTypes.SINGLE`、`GateTypes.TWO`、`GateTypes.MULTI`/`is_long_range` の 3 分岐しか持たない。リセット・測定命令のクラスも存在しない。

これはコード中に明示的に記述されている（`qudit_gksl_circuit_simulator.py:19–21`）:

> MQT-Qudits backends (`tnsim`/`misim`) currently do **not** implement mid-circuit reset, so that combined circuit is a *visualization artefact* and is not executable on a MQT-Qudits backend in its current form.

### 3.3 欠如③：状態ベクトルの指数関数的なサイズ爆発

「mid-circuit リセットが使えない」問題を回避する唯一の方法は、**チャネルごとに独立したアンシラ qudit を割り当てる**「フレッシュアンシラ（fresh ancilla）」構成である。

このリポジトリには `build_backend_executable_per_step_circuit` メソッドがあり、これが唯一の backend-executable Stinespring 回路構築関数である。しかし：

| 条件 | アンシラ数 | 全 qudit 数 | 状態ベクトルサイズ |
|------|-----------|------------|-----------------|
| N=2 分子、palindromic=False | 12 | 14 | $3^{14} \approx 4.8 \times 10^6$（TNSim で実行可） |
| N=2 分子、palindromic=True | 24 | 26 | $3^{26} \approx 2.5 \times 10^{12}$（メモリ不可能） |
| N=4 分子、palindromic=False | 48 | 52 | $3^{52} \approx 1.7 \times 10^{24}$（物理的に不可能） |

（`qudit_gksl_circuit_simulator.py:727–732`, `STATUS_HONEST_2026-05.md:A-1 Part B`）

N=4 分子系（シナリオ 5）では palindromic の有無に関わらず TNSim で状態ベクトルを保持することは不可能である。これは実装の問題ではなく、**状態ベクトルシミュレーションの根本的な制約**（古典計算機のメモリが指数スケール）である。

A-1 Part B（`STATUS_HONEST_2026-05.md:108`）ではこれを実測で確認している：5 サブ回路中 2 件で `MemoryError (468 GiB)` が発生した。

---

## 4. 実際の回避策と残課題

### 4.1 N=2 + palindromic=False の場合のみ TNSim で実行可能

`STATUS_HONEST_2026-05.md` A-1 Part C:

> A-1 Part C (interleaved_n2 full per-step): 1.3s, tnsim_ok=True, ‖Δρ_sys‖_F = 0.0, Tr=1.0, worst gate range = 6, intermediate ≈ 8 MB

`ancilla_layout='interleaved_n2'`（`qudit_gksl_circuit_simulator.py:609`）でアンシラをシステム qudit の近傍に配置することで、MULTI ゲートの作用範囲（min〜max qudit インデックスの連続区間）を最小化し、TNSim の中間テンソルサイズを $3^6 \times 3^6 \approx 8\,\text{MB}$ に抑えている。

ただしこれは N=2、palindromic=False のみ動作し、N=4 では依然不可能（`qudit_gksl_circuit_simulator.py:638–648`）。

### 4.2 N=4 の本来の解決策：DMSim

N=4 に対する実際の解決策は TNSim ではなく `DMSim` バックエンドである。

DMSim は状態ベクトルを追わず密度行列 $\rho$（$D \times D$ 行列、$D = \prod d_i$）を直接更新し、`KrausChannel` 命令を `apply_kraus_to_density` で処理する。アンシラ qudit は不要で、Lindblad チャネルを直接 Kraus 演算子として与える。

`STATUS_HONEST_2026-05.md` A-1 DMSim:

> `QuditGKSLSimulator(execute_on_backend="dmsim")` で N=4 ボソン無し（Cell 9）を MQT-Qudits backend で実行成功（‖Δρ‖_F = 8.55e-14）

---

## 5. まとめ：3 層の制約

| 制約 | 内容 | コード根拠 |
|------|------|-----------|
| ① 命令非対応 | `KrausChannel` を渡すと `NotImplementedError`（`to_matrix` が未定義） | `tnsim.py:91`、`kraus_channel.py:154` |
| ② 操作欠如 | mid-circuit リセット・測定命令が TNSim に存在しない | `tnsim.py` 全体（分岐は GateTypes のみ） |
| ③ メモリ爆発 | フレッシュアンシラ構成でも N=4 では状態ベクトルが指数的に大きすぎる（$3^{52}$ 等） | `STATUS_HONEST_2026-05.md` A-1 Part B、`qudit_gksl_circuit_simulator.py:641–647` |

これらは独立した制約であり、①を回避しても②が残り、②を回避しても③が残る。

TNSim で「Stinespring 拡張付き Trotter ステップを N=4 全ステップ実行する」ことは、コードに手を加えて① ②を解決したとしても、③により古典コンピュータのメモリ容量の壁に衝突する。この壁は TNSim の設計問題ではなく、**状態ベクトルシミュレーションの指数コストという普遍的な事実**である。
