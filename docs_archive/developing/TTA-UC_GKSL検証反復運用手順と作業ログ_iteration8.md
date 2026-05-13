# TTA-UC GKSL検証反復運用手順と作業ログ（iteration8）

## 反復の概要

- **反復番号**: 8
- **日時**: 2026-02-23
- **前回反復**: iteration7（全11シナリオPASS、2問題の修正完了）
- **目的**: iteration7の検証結果（20260223T063110Z）とコード全体を再度詳細分析し、残存する問題を修正

## 手順1: 検証結果分析

### 入力

- `developing/verification_results/tta_uc_gksl_verification_20260223T063110Z.json`
- `developing/verification_results/tta_uc_gksl_verification_20260223T063110Z.md`
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration7.md`
- `developing/検証結果分析_20260222_iteration7.md`

### 分析範囲

iteration7の検証結果だけでなく、コード全体（理論・設計・仕様）を
網羅的に再分析した。対象ファイル（全14ファイル）:

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
- `tutorials/gksl_visualization.py` — 可視化ユーティリティ

### 分析結果

iteration7の修正（問題1-2: allow_subnormalized、ノーマライズ済みメトリクス）は
すべて正しく実装されていることを確認。

以下の項目すべてについて理論的正当性を検証し、問題なしと判断:

- Stinespring dilation (G=[[0,L†],[L,0]], U=exp(-i√dt·G))
- 部分トレース (ブロック構造: ρ_out = ρ'[:d,:d] + ρ'[d:,d:])
- GKSL超演算子 (column-major vectorization)
- 2次対称Trotter分解 (H-D分割O(dt³), チャネル間O(dt²))
- 古典参照シミュレータ (厳密Liouvillian + expm_multiply)
- Weyl-Heisenberg脱分極 (d=3, d²-1=8非恒等演算子)
- 2-qubitパウリ脱分極 (d=4, d²-1=15非恒等演算子)
- ペア脱分極 (qudit d⁴-1=80, qubit 256-1=255)
- ディフェージング (qudit d=3全物理的, qubit d=4 |11⟩リーケージ)
- qubit-qutritマッピング (kron規約整合)
- Lindblad演算子構成 (26個: TTA6+蛍光4+燐光4+IC4+ISC_ST4+ISC_TS4)
- Holstein電子-フォノン結合
- ノーマライズ済みエントロピー S(ρ/tr) = (S_raw + tr·ln(tr))/tr
- サブノーマライズ密度行列検証

新たに以下の2つの問題を発見:

1. **問題1（中程度）**: ボソンシミュレータ間のトレース時系列の不整合
2. **問題2（軽微）**: quantum_fidelityの浮動小数点オーバーフロー（F > 1.0）

詳細は `developing/検証結果分析_20260222_iteration8.md` に記載。

## 手順2: コード修正

### 修正1: ボソンシミュレータのトレース時系列の一貫性修正

**対象ファイル:**
- `tutorials/qudit_gksl_boson_simulator.py`（行158, 166）
- `tutorials/qubit_gksl_boson_simulator.py`（行252, 260）

**変更内容:**

トレース記録を `Tr(ρ_total)`（電子+フォノン全空間）から `Tr(ρ_el)`（電子縮約密度行列）
に変更し、`classical_gksl_boson_simulator.py` と一致させた。

qudit_gksl_boson_simulator.py:
```python
# 変更前:
traces = [float(np.real(np.trace(rho)))]
traces.append(float(np.real(np.trace(rho))))
# 変更後:
traces = [float(np.real(np.trace(rho_el)))]
traces.append(float(np.real(np.trace(rho_el))))
```

qubit_gksl_boson_simulator.py:
```python
# 変更前:
traces = [float(np.real(np.trace(rho)))]
traces.append(float(np.real(np.trace(rho))))
# 変更後:
traces = [float(np.real(np.trace(rho_el)))]
traces.append(float(np.real(np.trace(rho_el))))
```

**理論的根拠:**

部分トレースの性質: `Tr(ρ_total) = Tr(Tr_ph(ρ_total)) = Tr(ρ_el)` であるため、
数学的に等価な変更である。しかし、entropy, purity, populations が全て ρ_el から
計算されているのに trace だけが ρ_total から計算されていたのは、意味の不一致であり、
これを解消する。

### 修正2: quantum_fidelity の F ≤ 1.0 クランプ

**対象ファイル:** `tutorials/run_tta_uc_gksl_verification.py`（quantum_fidelity関数）

**変更内容:**

```python
# 変更前:
return float(np.real(np.sum(np.sqrt(evals_m))) ** 2)
# 変更後:
return float(min(np.real(np.sum(np.sqrt(evals_m))) ** 2, 1.0))
```

**理論的根拠:**

Uhlmann忠実度は理論的に F(ρ,σ) ∈ [0, 1] であることが証明されている。
浮動小数点演算の誤差により F > 1 が生じる場合、理論が保証する上界 1.0 へのクランプは
数学的に正当化される（ごまかしのheuristicではない）。

## 手順3: 検証

### 単体テスト結果

- `test/python/tutorials/test_tta_uc_gksl_verification_script.py`: ✅ PASS (48.63s)

### ボソントレース一貫性検証

| 時点 | classical_boson Tr(ρ_el) | qudit_boson Tr(ρ_el) | qubit_boson Tr(ρ_el) |
|------|-------------------------|---------------------|---------------------|
| t=0 | 1.000000e+00 | 1.000000e+00 | 1.000000e+00 |
| t=1 | 1.000000e+00 | 1.000000e+00 | 1.000000e+00 |
| t=2 | 1.000000e+00 | 1.000000e+00 | 1.000000e+00 |

全3シミュレータが Tr(ρ_el) を報告し、一貫性が確保された。✅

### 忠実度クランプ検証

| テスト | F値 | F ≤ 1.0 |
|--------|-----|---------|
| F(ρ, ρ) | 1.000000 | ✅ |
| F(ρ, ρ+ε) | 1.000000 | ✅ |
| qudit_boson vs qubit_boson | 1.000000 | ✅ |

全ての忠実度が [0, 1] 範囲内。✅

### ボソンシナリオ統合テスト

- classical_boson: PASS ✅
- qudit_boson: PASS ✅
- qubit_boson: PASS ✅
- classical_boson vs qudit_boson: F=0.999913 ✅
- classical_boson vs qubit_boson: F=0.999913 ✅
- qudit_boson vs qubit_boson: F=1.000000 ✅

## 次のステップ

ユーザーによるフル検証スクリプトの実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification.py
```

全11シナリオの結果を `developing/verification_results/` にpushし、
iteration9 で結果を確認する。
