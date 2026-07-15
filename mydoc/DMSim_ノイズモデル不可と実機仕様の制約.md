# DMSim はなぜ NoiseModel が効かないのか、なぜ実機仕様にならないのか

> **方針：** このリポジトリの実コードのみを根拠とする。推測・希望的観測は含めない。

---

## 1. 結論（1 行）

**DMSim は設計上 `NoiseModel` オプションを受け取ると `ValueError` で即時停止する。ノイズは `KrausChannel` 命令として回路に埋め込む必要がある。また、`FakeBackend`（実機仕様擬似バックエンド）は TNSim を継承しており DMSim とは独立した別系統である。DMSim に実機キャリブレーションデータを接続する機能はコードに存在しない。**

---

## 2. DMSim が NoiseModel を拒否するコード

ファイル：`src/mqt/qudits/simulation/backends/dmsim.py:219-224`

```python
def run(self, circuit: QuantumCircuit, **options: ...) -> Job:
    job = Job(self)
    self._options.update(options)
    if self._options.get("noise_model", None) is not None:
        msg = (
            "DMSim does not accept a NoiseModel option. "
            "Add noise as KrausChannel instructions in the circuit instead."
        )
        raise ValueError(msg)
    ...
```

**`noise_model` が `None` 以外の値であれば問答無用で `ValueError` を送出する。**
これはバグではなく意図的な設計である（docstring に "noise must be supplied as `KrausChannel` instructions" と明記）。

同ファイルの `execute()` シグネチャ（`dmsim.py:232`）でも：

```python
def execute(
    self,
    circuit: QuantumCircuit,
    noise_model: NoiseModel | None = None,  # noqa: ARG002
    ...
```

`ARG002`（未使用引数の警告を抑制）コメントが付いている。`noise_model` 引数は受け取るが**完全に無視**する。

---

## 3. なぜ DMSim は NoiseModel を拒否するのか（設計上の理由）

### 3.1 NoiseModel は TNSim の仕組みである

`NoiseModel` / `NoisyCircuitFactory` は **TNSim・MISim の確率的ユニタリー射影**のための API である（`stochastic_sim.py`）。
仕組みは以下の通り：

1. `NoisyCircuitFactory.generate_circuit()` が `NoiseModel` を読み込み、各ゲートの直後に確率的に `NoiseX` / `NoiseZ` などのユニタリーゲートを挿入した「ノイジー回路」を生成する
2. TNSim が生成された回路を実行する（各ゲートに `to_matrix()` を呼ぶ）
3. これを `shots` 回繰り返し、純粋状態ベクトルのサンプル平均でノイズを近似する

このアプローチは **DMSim の密度行列表現とは設計が根本的に異なる**：
- TNSim は各ショットで純粋状態 $|\psi\rangle$ を維持し、ランダムユニタリーを射影的に適用する
- DMSim は密度行列 $\rho$ を一発で更新し、混合状態を直接表現する

両者を混在させることに意味がない（DMSim はそもそも「1 回で期待値を計算する」ための仕組みであり、shots による近似を必要としない）。

### 3.2 DMSim の「ノイズ」は KrausChannel 命令として回路に書く

DMSim でノイズを加えるには、`KrausChannel` 命令を回路に明示的に挿入する：

```python
# 例：depolarizing channel を手動で回路に追加する
kraus_ops = depolarisation_kraus(d=3, p=0.01)
circuit.custom_operation(KrausChannel(circuit, "depol", target=0, kraus_operators=kraus_ops, dimensions=3))
```

DMSim はこれを `_apply_instruction` 内で `isinstance(instruction, KrausChannel)` として識別し、
$\rho \to \sum_k K_k \rho K_k^\dagger$ を実行する（`dmsim.py:272-274`）。

**ユーザーがどの Kraus 演算子を使うかを明示的に選択しなければならない。** NoiseModel は「ゲート名 → ノイズ確率」のマッピングを持つが、その変換（確率値 → Kraus 演算子）を DMSim は自動的には行わない。

---

## 4. FakeBackend（実機仕様バックエンド）の構造

### 4.1 FakeBackend は TNSim の継承クラス

ファイル：`src/mqt/qudits/simulation/backends/fake_backends/fake_traps2three.py`

```python
from ..tnsim import TNSim

class FakeIonTraps2Trits(TNSim):   # ← TNSim を継承（DMSim ではない）
    def __init__(self, ...):
        super().__init__(
            provider=provider,
            name="FakeTrap2",
            description="A Fake backend of an ion trap qudit machine",
            ...
        )
        self.options["noise_model"] = self.__noise_model()   # ← TNSim の noise_model
```

**全 3 種類の FakeBackend（`FakeIonTraps2Trits`, `FakeIonTraps2Six`, `FakeIonTraps3Six`）は TNSim を継承しており、DMSim とは無関係である。**

### 4.2 FakeBackend の「実機仕様」の実態

`__noise_model()` メソッド（`fake_traps2three.py:87-110`）：

```python
def __noise_model(self) -> NoiseModel:
    """Noise model coded in plain sight, just for prototyping reasons."""
    local_error = Noise(probability_depolarizing=0.001, probability_dephasing=0.001)
    local_error_rz = Noise(probability_depolarizing=0.03, probability_dephasing=0.03)
    entangling_error = Noise(probability_depolarizing=0.1, probability_dephasing=0.001)
    ...
    return noise_model
```

docstring に **"just for prototyping reasons"** と明記されている。

これらの確率値（0.001, 0.03, 0.1 など）はコードに固定値（ハードコード）で書かれており、
**実際の実機（イオントラップ量子コンピュータ）からのキャリブレーションデータを読み込む機能は存在しない**。

### 4.3 エネルギーレベルグラフは「接続トポロジー」のみ

`energy_level_graphs` プロパティはイオントラップのエネルギー準位グラフを定義するが、
これは「どの準位間の遷移が使えるか」というトポロジー情報のみであり、
実機の $T_1$（緩和時間）・$T_2$（コヒーレンス時間）・ゲート時間・クロストーク等は含まない。

---

## 5. なぜ DMSim は「実機仕様」にならないのか

| 実機仕様のために必要なもの | DMSim での状況 |
|--------------------------|---------------|
| `NoiseModel`（ゲートごとのエラー率） | **拒否される**（`ValueError`） |
| FakeBackend の実機 calibration データ | FakeBackend は TNSim 継承。DMSim と無関係 |
| $T_1$/$T_2$ コヒーレンス時間 | コードに存在しない |
| ゲート時間（pulse level） | コードに存在しない |
| 読み出しエラー（readout error） | コードに存在しない |
| クロストーク | コードに存在しない |

DMSim がノイズを扱うために必要なのは**ユーザーが Kraus 演算子を自分で計算して `KrausChannel` として回路に書くこと**だけである。
実機パラメータから Kraus 演算子を自動生成するパイプラインはこのリポジトリには存在しない。

---

## 6. TNSim + NoiseModel の「ノイズ近似」の問題

TNSim + FakeBackend が使う `NoiseModel` / `NoisyCircuitFactory` の確率的射影は：

1. **確率値は固定のハードコード**（実機データではない）
2. **適用されるノイズは離散的なユニタリー**（X, Y, Z ゲートのランダム射影）であり、連続的な Kraus チャンネルではない
3. コード内に `# TODO: ARE WE SURE THIS IS CORRECT?` というコメントが存在する（`noisy_circuit_factory.py:178`）

```python
prob_each = noise_info.probability_depolarizing / dim / dim  # TODO: ARE WE SURE THIS IS CORRECT?
```

これは開発者自身が確率の正確性に確信を持っていないことを示している。

---

## 7. 各バックエンドとノイズの対応まとめ

| バックエンド | ノイズ方式 | 実機仕様との関係 |
|------------|----------|---------------|
| TNSim | `NoiseModel` + `NoisyCircuitFactory` によるランダムユニタリー射影（Monte Carlo、shots 回） | FakeBackend がハードコードの確率値を注入（実機キャリブレーションなし） |
| MISim | 同上（`stochastic_simulation` 共有） | 同上 |
| DMSim | `KrausChannel` 命令を回路に明示挿入（決定論的、1 回） | FakeBackend が存在しない。Kraus 演算子はユーザーが手動で準備 |
| FakeIonTraps* | TNSim 継承 + ハードコードの `NoiseModel` | 実機からのデータ取り込みなし。prototyping 用と明記 |

---

## 8. コード参照箇所一覧

| 内容 | ファイル | 行 |
|------|---------|-----|
| DMSim の `noise_model` 拒否 | `src/mqt/qudits/simulation/backends/dmsim.py` | 219–224 |
| DMSim の `execute()` が noise_model を無視 | `src/mqt/qudits/simulation/backends/dmsim.py` | 232（ARG002コメント） |
| DMSim docstring「noise as KrausChannel」 | `src/mqt/qudits/simulation/backends/dmsim.py` | 201–203 |
| FakeBackend が TNSim を継承 | `src/mqt/qudits/simulation/backends/fake_backends/fake_traps2three.py` | 12,18 |
| FakeBackend の "prototyping reasons" docstring | `src/mqt/qudits/simulation/backends/fake_backends/fake_traps2three.py` | 87 |
| FakeBackend のハードコードノイズ確率 | `src/mqt/qudits/simulation/backends/fake_backends/fake_traps2three.py` | 88–110 |
| NoisyCircuitFactory のランダムユニタリー射影 | `src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py` | 173–220 |
| `# TODO: ARE WE SURE THIS IS CORRECT?` | `src/mqt/qudits/simulation/noise_tools/noisy_circuit_factory.py` | 178 |
| TNSim の stochastic_simulation 呼び出し | `src/mqt/qudits/simulation/backends/stochastic_sim.py` | 47–51 |
| KrausChannel 識別（DMSim） | `src/mqt/qudits/simulation/backends/dmsim.py` | 272–274 |
