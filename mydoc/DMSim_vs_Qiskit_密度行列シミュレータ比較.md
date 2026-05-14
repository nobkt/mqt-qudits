# MQT-Qudits DMSim と Qiskit Aer 密度行列シミュレータの比較

> **本文書はコード実装の直接観察に基づく。**  
> 引用元:  
> - `src/mqt/qudits/simulation/backends/dmsim.py`  
> - `src/mqt/qudits/quantum_circuit/gates/kraus_channel.py`  
> - `tutorials/qiskit_qubit_gksl_simulator.py`（Qiskit Aer の使用箇所）  
> - `tutorials/dmsim_shot_simulator.py`（両バックエンド並走比較）
>
> **Qiskit Aer の内部実装については本リポジトリのコードから直接読めない部分がある。**  
> そのような箇所は「公開ドキュメント・API の観察」と明記する。

---

## 1. 前提：何と何を比較するか

| 対象 | 正確な名称 | 本リポジトリでの役割 |
|------|-----------|-------------------|
| **MQT-Qudits DMSim** | `mqt.qudits.simulation.backends.DMSim` | qudit（d 進量子ビット）の密度行列シミュレータ |
| **Qiskit Aer 密度行列** | `AerSimulator(method="density_matrix")` | qubit 専用の密度行列シミュレータ |

このリポジトリの `tutorials/qiskit_qubit_gksl_simulator.py` では、qutrit 系を qubit 対にエンコードしてから Qiskit Aer に渡している。DMSim はエンコードなしで直接 qutrit を扱う。

---

## 2. 類似点（コードから確認できる事実）

### 2.1 基本的な数学的定式

両者とも **同じ更新則** を実行する:

```
ρ → U ρ U†          （ユニタリゲート）
ρ → Σ_k K_k ρ K_k†  （Kraus チャネル）
```

`dmsim_shot_simulator.py` の docstring（行 14–35）は両バックエンドが同一の Trotter ステップ構造（`cu_multi` ユニタリ + `KrausChannel`）を用いていると明記している。

### 2.2 Kraus 命令の受け付け

- MQT-Qudits DMSim: `KrausChannel` 命令（`mqt.qudits.quantum_circuit.gates.KrausChannel`）
- Qiskit Aer: `qiskit.quantum_info.Kraus` 命令

どちらも `Σ_k K_k ρ K_k†` を直接実行する。Kraus 演算子のリストを渡すという API の概念は共通。

### 2.3 初期密度行列の入力

- MQT-Qudits DMSim: `execute(circuit, initial_density_matrix=rho)` オプション（`dmsim.py:233`）
- Qiskit Aer: `circuit.set_density_matrix(DensityMatrix(rho))` 命令

どちらも任意の初期密度行列を渡せる。

### 2.4 純粋古典シミュレーション

どちらも **CPU 上の古典数値シミュレーション** であり、実量子ハードウェアで動作しているわけではない。`qiskit_qubit_gksl_simulator.py` の docstring（行 59–65）に「This is a Qiskit Aer **simulator** of the qubit circuit, not a real hardware run.」と明記されている。

---

## 3. 相違点（コードから確認できる事実）

### 3.1 対象システムの次元

| 項目 | MQT-Qudits DMSim | Qiskit Aer |
|------|-----------------|------------|
| 1 qudit の次元 | **任意の正整数 `d`**（qubit, qutrit, ququart …） | **2 固定**（qubit のみ） |
| 次元の記述方法 | `circuit.dimensions` リスト（各 qudit が個別の `d_i` を持つ） | qubit 数（全 qudit が `d=2`） |
| D の定義 | `D = Π_i d_i` | `D = 2^N` |

このリポジトリで qutrit を Qiskit Aer に渡す際は `tutorials/qiskit_qubit_gksl_simulator.py` が qutrit（3 次元）を qubit 対（4 次元、うち 1 次元は使用禁止）にエンコードしてから渡している（`_embed_local_operator` 関数）。DMSim はそのエンコードを必要としない。

### 3.2 ノイズモデルの渡し方

| 項目 | MQT-Qudits DMSim | Qiskit Aer |
|------|-----------------|------------|
| `noise_model` オプション | **受け付けない。`ValueError` を送出**（`dmsim.py:219–224`） | `NoiseModel` オブジェクトを受け付ける（Qiskit Aer の標準 API） |
| ノイズの注入方法 | 回路内に `KrausChannel` 命令として明示的に追加する必要がある | `NoiseModel` で後付けに注入できる（ゲートに自動挿入） |

`dmsim.py` の設計コメント: 「noise must be supplied as KrausChannel instructions so that the simulator and the user agree on what the physics is」

### 3.3 テンソル演算の実装

| 項目 | MQT-Qudits DMSim | Qiskit Aer |
|------|-----------------|------------|
| 演算ライブラリ | NumPy `tensordot` + `transpose`（純 Python、`dmsim.py:121,129`） | C++/BLAS ベースの内部実装（本リポジトリのコードからは読めない） |
| GPU 対応 | **なし**（コードに該当なし） | あり（`AerSimulator(device="GPU")`、公開 API より） |
| スパース最適化 | **なし** | あり（Aer の自動最適化、本リポジトリのコードからは詳細不明） |

### 3.4 `JobResult` の返却形式

| 項目 | MQT-Qudits DMSim | Qiskit Aer |
|------|-----------------|------------|
| ρ の取り出し方 | `job.result().density_matrix`（`DensityMatrixJobResult.density_matrix` プロパティ） | `result.data(0)["density_matrix"]`（`save_density_matrix()` 命令が必要） |
| `state_vector` 属性 | `ρ.reshape(1, -1)`（flatten された密度行列。state vector ではない） | 存在しない（密度行列専用） |

### 3.5 CPTP 検証の責任の所在

| 項目 | MQT-Qudits DMSim | Qiskit Aer |
|------|-----------------|------------|
| 入力 Kraus の CPTP 検証 | `KrausChannel.__init__` が `‖Σ K†K − I‖_F ≤ 1e-9` を**構築時に強制**（`kraus_channel.py`） | Qiskit の `Kraus` クラスが検証するか否かは本リポジトリからは確認できない |
| 違反時の挙動 | `ValueError`、補正なし | 不明（本リポジトリのコードから判断不可） |

### 3.6 回路ごとの ρ 引き継ぎ方法

- **MQT-Qudits DMSim**: `execute(circuit, initial_density_matrix=rho)` で次ステップの ρ を渡す。これは `dmsim.py:233` の引数として実装されている。
- **Qiskit Aer**: 各ステップで新しい `QuantumCircuit` に `set_density_matrix(DensityMatrix(rho))` 命令を追加し、前ステップの ρ を埋め込む（`qiskit_qubit_gksl_simulator.py:396`）。

どちらも**同じ目的**（ステップ間で ρ を受け渡す）を達成しているが、API が異なる。

---

## 4. このリポジトリ内での実際の使われ方

`tutorials/dmsim_shot_simulator.py` および `tutorials/qiskit_qubit_gksl_simulator.py` では、**同一の物理系（qutrit 4 サイト GKSL 系）**を次のように並走させている:

```
qudit path:  MQT-Qudits DMSim (KrausChannel 命令, d=3 ネイティブ)
qubit path:  Qiskit Aer AerSimulator(method="density_matrix")
             （qutrit を qubit 対にエンコード後）
```

docstring（`dmsim_shot_simulator.py:60–65`）に明記: 両者の一致は「qubit エンコーディングの正しさを検証するもの」であり、「独立した開放系ダイナミクスを 2 つの方法で解いている」わけではない。

---

## 5. 事実として確認できないこと

以下は本リポジトリのコードからは読み取れない:

- Qiskit Aer の `density_matrix` メソッドの内部テンソル演算の詳細
- Qiskit の `Kraus` クラスが入力を検証するか否か、およびその許容値
- Qiskit Aer がスパース表現や GPU を自動的に選択する条件
- MQT-Qudits DMSim が将来 GPU 対応する予定があるか否か（未実装）

---

## 6. まとめ表

| 比較項目 | MQT-Qudits DMSim | Qiskit Aer density_matrix |
|---------|-----------------|--------------------------|
| 対象 | qudit（任意次元 d） | qubit（d=2 固定） |
| 更新則 | `ρ → UρU†` / `ρ → ΣKρK†` | 同じ |
| Kraus 命令 | `KrausChannel`（CPTP を構築時強制検証） | `qiskit.quantum_info.Kraus` |
| noise_model | **受け付けない**（ValueError） | 受け付ける |
| ノイズ注入 | 回路内の明示的 Kraus 命令のみ | NoiseModel で後付け可 |
| 初期 ρ 設定 | `initial_density_matrix` 引数 | `set_density_matrix()` 命令 |
| 実装言語 | Python（NumPy tensordot） | C++（詳細はコード外） |
| GPU 対応 | なし | あり（公開 API より） |
| 出力 API | `job.result().density_matrix` | `result.data(0)["density_matrix"]` |
| qudit ネイティブ対応 | **あり**（d=3, 4, … を直接扱う） | **なし**（qubit エンコードが必要） |
