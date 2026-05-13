# TTA-UC GKSL検証反復運用手順と作業ログ（iteration7）

## 反復の概要

- **反復番号**: 7
- **日時**: 2026-02-23
- **前回反復**: iteration6（全11シナリオPASS、2問題の修正完了）
- **目的**: iteration6の検証結果（20260223T052721Z）とコード全体を再度詳細分析し、残存する問題を修正

## 手順1: 検証結果分析

### 入力

- `developing/verification_results/tta_uc_gksl_verification_20260223T052721Z.json`
- `developing/verification_results/tta_uc_gksl_verification_20260223T052721Z.md`
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration6.md`
- `developing/検証結果分析_20260222_iteration6.md`

### 分析範囲

iteration6の検証結果だけでなく、コード全体（理論・設計・仕様）を
網羅的に再分析した。対象ファイル:

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

iteration6の修正（問題1-2）はすべて正しく実装されていることを確認。
以下の項目すべてについて理論的正当性を検証し、問題なしと判断:

- Stinespring dilation (G=[[0,L†],[L,0]], U=exp(-i√dt·G))
- 部分トレース (ブロック構造: ρ_out = ρ'[:d,:d] + ρ'[d:,d:])
- 2次対称Trotter分解 (H-D分割O(dt³), チャネル間O(dt²))
- 古典参照シミュレータ (厳密Liouvillian + expm_multiply)
- Weyl-Heisenberg脱分極 (d=3, d²-1=8非恒等演算子)
- 2-qubitパウリ脱分極 (d=4, d²-1=15非恒等演算子)
- ペア脱分極 (qudit d⁴-1=80, qubit 256-1=255)
- ディフェージング (qudit d=3全物理的, qubit d=4 |11⟩リーケージ)
- qubit-qutritマッピング (kron規約整合)
- Lindblad演算子構成 (26個: TTA6+蛍光4+燐光4+IC4+ISC_ST4+ISC_TS4)

新たに以下の2つの問題を発見:

1. **問題1（重要）**: validate_density_matrix のトレース許容値ハック（ごまかしのheuristic）
2. **問題2（中程度）**: サブノーマライズ密度行列からのエントロピー・純度の比較不可能性

詳細は `developing/検証結果分析_20260222_iteration7.md` に記載。

## 手順2: コード修正

### 修正1: validate_density_matrix のサブノーマライズ対応

**対象ファイル:** `tutorials/gksl_validation.py`

**変更内容:**

`validate_density_matrix` に `allow_subnormalized: bool = False` パラメータを追加。

変更前:
```python
trace_val = float(np.real(np.trace(rho)))
if abs(trace_val - 1.0) > tol["trace"]:
    errors.append(...)
```

変更後:
```python
trace_val = float(np.real(np.trace(rho)))
if allow_subnormalized:
    if trace_val > 1.0 + tol["trace"]:
        errors.append(f"Step {step}: Trace = {trace_val} (exceeds 1.0)")
    if trace_val < 0.0:
        errors.append(f"Step {step}: Trace = {trace_val} (negative)")
else:
    if abs(trace_val - 1.0) > tol["trace"]:
        errors.append(f"Step {step}: Trace = {trace_val} (expected 1.0)")
```

**理論的根拠:**

qubit符号化の禁止状態リーケージにより、qutrit空間に射影された密度行列は
Tr(ρ) < 1 となる。これはサブノーマライズ密度行列であり、正しい検証条件は
|Tr(ρ) - 1| < ε ではなく 0 ≤ Tr(ρ) ≤ 1 + ε である。

### 修正2: 検証スクリプトのheuristic除去とノーマライズメトリクス追加

**対象ファイル:** `tutorials/run_tta_uc_gksl_verification.py`

**変更内容:**

1. 任意トレース許容値 `{"trace": 0.5}` を除去し、`allow_subnormalized=True` に置換:

変更前:
```python
dm_tolerance = {"trace": 0.5} if allow_leakage else None
density_validation = validate_density_matrix(
    result["rho_final"], step=n_steps, tolerance=dm_tolerance
)
```

変更後:
```python
density_validation = validate_density_matrix(
    result["rho_final"], step=n_steps, allow_subnormalized=allow_leakage
)
```

2. サブノーマライズ状態に対してノーマライズされたエントロピー/純度を追加:

```python
if allow_leakage and tr > 1e-10:
    normalized_entropy = (raw_entropy + tr * log(tr)) / tr
    normalized_purity = raw_purity / (tr * tr)
```

これにより、qubit_noisy_shot のメトリクスを他シナリオと定量的に比較可能になる。

## 手順3: 検証

### 単体テスト結果

- `test/python/tutorials/test_tta_uc_gksl_verification_script.py`: ✅ PASS

### validate_density_matrix エッジケース検証

| テスト | 入力 | 期待 | 結果 |
|--------|------|------|------|
| 正常 trace=1 | diag(0.5,0.3,0.2) | PASS | ✅ PASS |
| サブノーマライズ（許可） | diag(0.3,0.2,0.1), allow=True | PASS | ✅ PASS |
| サブノーマライズ（不許可） | diag(0.3,0.2,0.1), allow=False | FAIL | ✅ PhysicsViolationError |
| trace>1（許可あり） | diag(0.6,0.3,0.2), allow=True | FAIL | ✅ PhysicsViolationError |
| 負trace（許可あり） | diag(-0.1,0,0), allow=True | FAIL | ✅ PhysicsViolationError |

### qubit_noisy_shot + qudit_noisy_shot 統合テスト

- qubit_noisy_shot: trace=0.770, valid=True, normalized_entropy=1.507, normalized_purity=0.457
- qudit_noisy_shot: trace=1.000, valid=True, entropy=1.830, purity=0.327
- ノーマライズされたメトリクスにより定量比較が可能

## 次のステップ

ユーザーによるフル検証スクリプトの実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification.py
```

全11シナリオの結果を `developing/verification_results/` にpushし、
iteration8 で結果を確認する。
