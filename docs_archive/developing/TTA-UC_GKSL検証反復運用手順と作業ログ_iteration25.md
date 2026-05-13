# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 25

## 日時
2026-02-25

## 目的
Iteration 24 の検証結果を詳細分析した結果、3件の分析改善点を発見。分析スクリプトに Frobenius 距離、合成 Rate(T) 予測、合成モデルフィッティングを追加し、ST vs Exact の収束挙動を完全に説明可能にする。

## Iteration 24 で確認された成果（正常動作の確認）

1. **全個別成分の収束次数が理論と完全一致** ✅
   - ST vs SCPT: Rate(T) = 1.00（全4ペア、最大偏差 0.001）
   - SCPT vs Exact: Rate(T) = 2.00（全4ペア、最大偏差 0.013）
   - CT vs Exact: Rate(T) = 2.00（全4ペア、最大偏差 0.013）
   - SCPT vs CT: Rate(T) = 2.00（全4ペア、最大偏差 0.018）

2. **高次フィッティングの高精度** ✅
   - 全成分で最大相対フィット誤差 < 0.23%

3. **Richardson 外挿との一致** ✅
   - 全成分で差 < 0.21%

4. **シミュレータコードにバグなし** ✅

## Iteration 24 で発見された分析上の問題点

### 問題点1（分析予測の不正確さ）：ST vs Exact の Rate(T) = 0.87
- 予測範囲：0.95〜1.05
- 実測値：0.87
- 原因：誤差相殺の数学的帰結（合成モデルで定量的に説明可能）
- 修正：合成モデルからの Rate(T) 予測を追加

### 問題点2（分析手法の限界）：トレース距離ベースの角度推定が不安定
- 角度変動：119°〜131°（12° 変動）
- 原因：L₁ ノルム（トレース距離）での余弦定理は近似
- 修正：Frobenius ノルム（L₂）ベースの厳密な角度推定を追加

### 問題点3（出力不足）：合成 Rate(T) の理論予測値が未出力
- 予測値と実測値の比較ができない
- 修正：合成モデルフィッティング + 予測 Rate(T) の出力

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260225_iteration25.md` を作成
- 問題点3件を数値的根拠とともに記載
- Rate(T) の定量的予測式を導出

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration25.py` を作成
- Test A-D: Iteration 24 から変更なし
- Test E: 以下の分析改善を追加
  - **frobenius_distance() 関数**: Frobenius 距離の計算
  - **Frobenius 距離の計算**: 主要ペアの Frobenius 距離を convergence_data に追加
  - **Part 6f**: 合成モデルフィッティング T(ST,ex)/dt = p₀ + p₁·dt
  - **Part 6g**: Frobenius ベースの誤差ベクトル角度推定（厳密な余弦定理）
  - **Part 6h**: 合成 Rate(T) 予測と実測の比較
  - Part 6a-e: Iteration 24 と同一

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration25.py
```
結果は `developing/verification_results/iteration25_frobenius_analysis_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 25 の結果を確認し、以下を検証:
  1. Frobenius ベースの角度 θ_F が dt に対して安定であること（変動 < 5°）
  2. 合成モデルの Rate(T) 予測が実測と ±0.03 以内で一致すること
  3. T(ST,ex)/dt の直接フィット p₀ が a₀（Stinespring 漸近値）と概ね一致すること
  4. Rate(T) > 0.95 に必要な dt の推定が妥当であること

## 期待される結果

### Frobenius ベース角度推定
- θ_F が dt に対して比較的安定（変動 < 5°）
- トレース距離ベース θ_T との差が定量化される
- L₁/L₂ ノルムの差を数値的に示す

### 合成 Rate(T) 予測
- 予測 Rate(T) と実測 Rate(T) の差が ±0.03 以内
- Rate(T) > 0.95 に必要な dt ≈ 0.05（n_steps=200）
- 合成モデルによる Rate(T) の説明が完全に一致

### 合成モデルフィッティング
- T(ST,ex)/dt = p₀ + p₁·dt で p₁ < 0（相殺効果）
- p₀ と a₀ の差 < 数%（大きい dt のデータによるバイアスのため完全一致はしない）

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論：Iteration 26 以降で実施

- Iteration 25 の結果を含めてからノートブックを更新する方が効率的
- Frobenius 距離と合成 Rate(T) の分析手法が安定した後に改修すべき
- 詳細は `developing/検証結果分析_20260225_iteration25.md` を参照

## 参照ファイル
- `developing/検証結果分析_20260225_iteration25.md` — 詳細分析
- `developing/検証結果分析_20260225_iteration24.md` — iteration 24 の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration24.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration25.py` — 検証スクリプト
- `tutorials/qudit_gksl_simulator.py` — シミュレータ（変更なし）
