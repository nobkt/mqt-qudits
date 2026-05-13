# TTA-UC GKSL検証反復運用手順と作業ログ（iteration6）

## 反復の概要

- **反復番号**: 6
- **日時**: 2026-02-22
- **前回反復**: iteration5（全11シナリオPASS、問題1-4の修正完了）
- **目的**: iteration5の検証結果とコード全体を再度詳細分析し、残存する問題を修正

## 手順1: 検証結果分析

### 入力

- `developing/verification_results/tta_uc_gksl_verification_20260222T143658Z.json`
- `developing/verification_results/tta_uc_gksl_verification_20260222T143658Z.md`
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration5.md`
- `developing/検証結果分析_20260222_iteration5.md`

### 分析範囲

iteration5 の検証結果だけでなく、コード全体（理論・設計・仕様）を
詳細に再分析した。対象ファイル:

- `tutorials/stinespring_utils.py` — Stinespring dilation + GKSL超演算子
- `tutorials/gksl_math_utils.py` — リンドブラッド演算子・ハミルトニアン構成
- `tutorials/gksl_physical_parameters.py` — 物理パラメータ
- `tutorials/gksl_validation.py` — 密度行列検証
- `tutorials/classical_gksl_simulator.py` — 古典参照シミュレータ
- `tutorials/classical_gksl_boson_simulator.py` — 古典ボソンシミュレータ
- `tutorials/qudit_gksl_simulator.py` — Qudit密度行列シミュレータ
- `tutorials/qubit_gksl_simulator.py` — Qubit密度行列シミュレータ
- `tutorials/qudit_gksl_boson_simulator.py` — Quditボソンシミュレータ
- `tutorials/qubit_gksl_boson_simulator.py` — Qubitボソンシミュレータ
- `tutorials/qudit_gksl_shot_simulator.py` — Quditショットシミュレータ
- `tutorials/qubit_gksl_shot_simulator.py` — Qubitショットシミュレータ
- `tutorials/run_tta_uc_gksl_verification.py` — 検証スクリプト

### 分析結果

iteration5 の修正（問題1-4）はすべて正しく実装されていることを確認。
新たに以下の2つの問題を発見:

1. **問題1（重要）**: ハミルトニアン転送ゲート後のディフェージングノイズ欠落
2. **問題2（中程度）**: 検証スクリプトのqubit/quditノイズモデル非対称（T1熱緩和）

詳細は `developing/検証結果分析_20260222_iteration6.md` に記載。

## 手順2: コード修正

### 修正1: ハミルトニアン転送ゲート後のディフェージングノイズ追加

**対象ファイル（2ファイル）:**
- `tutorials/qudit_gksl_shot_simulator.py` (QuditGKSLNoisyShotSimulator)
- `tutorials/qubit_gksl_shot_simulator.py` (QubitGKSLNoisyShotSimulator)

**変更内容:**

`_trotter_step_trajectory` メソッドの2箇所のハミルトニアン転送ゲートノイズブロック
（最初のhalf-Hamiltonian後と最後のhalf-Hamiltonian後）で、
脱分極ノイズに加えてディフェージングノイズを適用。

変更前:
```python
for i, j in self.params.neighbors:
    psi = depolarization_pair(psi, i, j, ...)
```

変更後:
```python
for i, j in self.params.neighbors:
    psi = depolarization_pair(psi, i, j, ...)
    if self.p_dephasing > 0.0:
        psi = dephasing_single(psi, i, ...)
        psi = dephasing_single(psi, j, ...)
```

**理論的根拠:**

ハミルトニアン転送ゲートとリンドブラッドStinespringゲートは
どちらも2体量子ゲートであり、同一のノイズモデルが適用されるべき。
現在のコードではStinespringゲートにのみディフェージングが適用されており、
転送ゲートのディフェージングが欠落していた。

**影響:**
- Trotterステップあたり追加ディフェージング: 2回 × 3対 × 2サイト = 12回
- qubit/qudit両方に対称的に適用（比較の公平性を維持）

### 修正2: 検証スクリプトのノイズモデル対称化

**対象ファイル:** `tutorials/run_tta_uc_gksl_verification.py`

**変更内容:**

`qubit_noisy_shot` シナリオから T1/t_gate パラメータを除去:

変更前:
```python
simulator = QubitGKSLNoisyShotSimulator(
    params, p_depol=0.01, p_dephasing=0.005,
    T1=50000.0, t_gate=300.0,
)
```

変更後:
```python
simulator = QubitGKSLNoisyShotSimulator(
    params, p_depol=0.01, p_dephasing=0.005,
)
```

**理論的根拠:**

T1 熱緩和は qubit ハードウェア固有のノイズであり、qudit には存在しない。
符号化方式（qubit vs qudit）の比較を行う場合、
ノイズモデルは同一でなければ比較の公平性が担保できない。

T1 熱緩和を含めた検証を行いたい場合は、
`QubitGKSLNoisyShotSimulator` を直接使用してカスタムパラメータで実行可能。

## 手順3: 検証

### ローカル検証結果

全シナリオの動作確認:
- classical: ✅ PASS
- qubit: ✅ PASS
- qudit: ✅ PASS
- classical_boson: ✅ PASS
- qubit_boson: ✅ PASS
- qudit_boson: ✅ PASS
- qudit_shot: ✅ PASS
- qubit_shot: ✅ PASS
- qudit_noisy_shot: ✅ PASS（dephasing on Hamiltonian gates）
- qubit_noisy_shot: ✅ PASS（T1=None, matched noise model）
- unitary: ✅ PASS

### ノイズモデルの対称性確認

修正後の qubit_noisy_shot:
- noise_params: {p_depol: 0.01, p_dephasing: 0.005, T1: None, p_reset: 0.0}
- qudit_noisy_shot と同一のノイズパラメータ

## 次のステップ

ユーザーによるフル検証スクリプトの実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification.py
```

全11シナリオの結果を `developing/verification_results/` にpushし、
iteration7 で結果を確認する。
