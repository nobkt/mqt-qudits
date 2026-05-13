# TTA-UC GKSL検証反復運用手順と作業ログ（iteration11）

## 反復の概要

- **反復番号**: 11
- **日時**: 2026-02-23
- **前回反復**: iteration10（ノートブック4箇所修正、検証済みコードとの整合）
- **目的**: iteration10で更新されたノートブック `quantum_dynamics_gksl_comparison.ipynb` の各セル実行出力、コード全体、理論・設計・仕様を詳細分析し、残存問題を発見・修正

## 手順1: ノートブック全セル分析

### 入力

- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration10.md`
- `developing/検証結果分析_20260222_iteration9.md`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`（iteration10更新版）
- `tutorials/qubit_gksl_shot_simulator.py`
- `tutorials/qudit_gksl_shot_simulator.py`
- `tutorials/run_tta_uc_gksl_verification.py`
- `developing/verification_results/tta_uc_gksl_verification_20260223T070836Z.md`（最新検証結果）

### 分析方法

全37セル（Cell 0〜36）の以下を分析:
1. ソースコードの正確性
2. 実行出力の物理的妥当性
3. マークダウン説明文とコード/出力の整合性
4. シミュレータ実装コードとの整合性
5. 理論的正当性

### 発見された問題点

#### 問題1: Cell 25 マークダウンのシナリオ3cノイズ記述不整合（テキストエラー）

- **対象**: Cell 25（ショットベースシミュレーション説明）のシナリオ一覧表
- **現状**: `| 3c | Qubit Shot-based | 脱分極+熱緩和 | QubitGKSLNoisyShotSimulator |`
- **正**: `| 3c | Qubit Shot-based | 脱分極+位相緩和 | QubitGKSLNoisyShotSimulator |`
- **影響**: iteration10でCell 32, 33, 36は修正されたが、Cell 25の表のみ修正漏れ
- **根拠**: Cell 33のコードは `QubitGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.005)` を使用しており、熱緩和（T1）ではなく位相緩和を使用

#### 問題2: Cell 33 禁止状態リーケージの診断出力欠如（重要）

- **対象**: Cell 33 シナリオ3c実行結果
- **現状の出力**:
  ```
  Final populations: N_S0=0.4786, N_T1=0.3669, N_S1=0.4030
  ```
  人口合計 = 0.4786 + 0.3669 + 0.4030 = **1.2485**（N_molecules=4に対して大幅に不足）
- **問題**: 禁止状態 |11⟩ へのリーケージにより人口の約69%が消失しているが、この情報が一切表示されない。読者に物理的挙動を誤解させる可能性が高い
- **必要な追加出力**:
  - 最終密度行列のトレース（Tr(ρ) < 1 はリーケージの指標）
  - 禁止状態への最終測定数（forbidden_count）
  - トレース欠損の時系列（リーケージの進行を可視化）
- **根拠**: `result3c` は `rho_final`, `trace`, `forbidden_count` を既に含んでいる（`qubit_gksl_shot_simulator.py` の `simulate()` メソッドが返す）。これらを表示するだけでよい

#### 問題3: Cell 35 忠実度計算がサブノーマライズ密度行列に対して不正確（重要）

- **対象**: Cell 35 のStatevector vs Shot-Based忠実度比較
- **現状のコード**:
  ```python
  F_qubit_noisy_shot = quantum_fidelity(result3['rho_final'], result3c['rho_final'])
  ```
- **問題**: `result3['rho_final']` はTr=1の正規な密度行列だが、`result3c['rho_final']` はTr≈0.31（100トロッターステップ後のリーケージにより大幅にサブノーマライズ）。Uhlmann忠実度 $F(\rho,\sigma) = (\mathrm{Tr}\sqrt{\sqrt{\rho}\sigma\sqrt{\rho}})^2$ は**両方の入力がTr=1の密度行列**であることを前提とする。サブノーマライズ行列を直接入力すると、忠実度は物理的意味を失う
- **現状の出力**:
  ```
  Qubit DM vs Qubit Shot (noisy): F = 0.033438
  ```
  これはリーケージによるトレース減少をそのまま反映した値であり、ノイズによる状態変化の指標としては不適切
- **検証スクリプトとの不整合**: `run_tta_uc_gksl_verification.py` の `_compute_fidelities()` 関数（行315-318）では正しくノーマライズしてから忠実度を計算:
  ```python
  tr_a = float(np.real(np.trace(rho_a)))
  tr_b = float(np.real(np.trace(rho_b)))
  rho_a_norm = rho_a / tr_a if tr_a > 1e-10 else rho_a
  rho_b_norm = rho_b / tr_b if tr_b > 1e-10 else rho_b
  f_val = quantum_fidelity(rho_a_norm, rho_b_norm)
  ```
- **修正方針**: 検証スクリプトと同一のノーマライズ手順を適用。これはヒューリスティックではなく、サブノーマライズ行列からの条件付き状態（物理部分空間に留まった場合の状態）に対する忠実度計算として数学的に正当
- **修正後の期待値**: 検証結果（n_steps=5）では、ノーマライズ後の忠実度は classical vs qubit_noisy = 0.393（classical vs qudit_noisy = 0.296より高い）。ノートブック（n_steps=100）ではリーケージがより大きいが、ノーマライズ後はqubit_noisyとqudit_noisyの忠実度差が大幅に縮小するはず

#### 問題4: Cell 35 qubit noisy vs qudit noisy の比較に関する説明不足

- **対象**: Cell 35 の出力解釈
- **現状**: 出力のみで、qubit noisyとqudit noisyの忠実度差の原因について説明がない
- **物理的背景**:
  - **Quditノイズ**: Weyl-Heisenberg演算子 $X^aZ^b$ はd=3空間内で完結。禁止状態が存在しないため、ノイズ後もTr(ρ)=1が厳密に保存
  - **Qubitノイズ**: 2-qubitパウリ演算子はd=4空間で作用。15個の非恒等パウリ中12個（80%）が物理状態→禁止状態 |11⟩ へのリーケージを引き起こす
  - この差異は符号化スキームの本質的な違いであり、実際のqubit量子コンピュータで2-qubitエンコーディングを使用した場合に必然的に発生する現象
- **修正方針**: Cell 35にノーマライズ前後の忠実度を両方表示し、トレース値も併記。Cell 34のマークダウンにリーケージの説明を追加

## 手順2: コード修正

### 修正1: Cell 25 マークダウン修正

**対象**: Cell 25 ショットベースシナリオ一覧表

```markdown
# 変更前:
| 3c | Qubit Shot-based | 脱分極+熱緩和 | QubitGKSLNoisyShotSimulator |
# 変更後:
| 3c | Qubit Shot-based | 脱分極+位相緩和 | QubitGKSLNoisyShotSimulator |
```

### 修正2: Cell 33 禁止状態リーケージ診断出力追加

**対象**: Cell 33 print文追加

以下の診断出力を追加:
- `Tr(rho_final)`: 最終密度行列トレース（1.0未満ならリーケージ発生）
- `forbidden_count`: 最終測定での禁止状態計数
- `Trace deficit (leakage)`: 最大トレース欠損

### 修正3: Cell 35 忠実度計算のノーマライズ対応

**対象**: Cell 35 忠実度計算

検証スクリプトと同一のノーマライズ手順を適用:
```python
def _normalize_rho(rho):
    """Normalize density matrix for fidelity computation."""
    tr = float(np.real(np.trace(rho)))
    return rho / tr if tr > 1e-10 else rho

# ノーマライズ忠実度（サブノーマライズ行列対応）
rho3_norm = _normalize_rho(result3['rho_final'])
rho3c_norm = _normalize_rho(result3c['rho_final'])
F_qubit_noisy_shot_norm = quantum_fidelity(rho3_norm, rho3c_norm)
```

### 修正4: Cell 34 マークダウン説明の拡充

qubit noisy特有の禁止状態リーケージに関する説明を追加。

## 手順3: 検証スクリプト作成

禁止状態リーケージの詳細診断を行う検証スクリプトを作成:
`tutorials/run_tta_uc_gksl_verification_iteration11.py`

このスクリプトは以下を出力:
1. qubit_noisy_shot のトレース欠損時系列
2. ノーマライズ前後の忠実度比較
3. 禁止状態リーケージ率の定量分析
4. qubit_noisy vs qudit_noisy の公平な比較

## 総合判定

iteration10でCells 32, 33, 36は正しく修正されたが、以下の問題が残存:
1. Cell 25のテキストエラー（修正漏れ）
2. Cell 33の禁止状態リーケージ診断出力の欠如
3. Cell 35のサブノーマライズ密度行列に対する不正確な忠実度計算
4. 比較の解釈に関する説明不足

これらはいずれもヒューリスティック/フォールバックではなく、
数学的に正当な修正と物理的事実に基づく診断出力の追加である。
