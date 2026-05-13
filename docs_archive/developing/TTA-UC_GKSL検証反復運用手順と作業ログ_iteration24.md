# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 24

## 日時
2026-02-25

## 目的
Iteration 23 の検証結果の詳細分析に基づき、分析の精度をさらに改善する。具体的には、dt 範囲の拡張、高次フィッティング、誤差ベクトル幾何学解析、Richardson 外挿を追加する。

## Iteration 23 で確認された成果（正常動作の確認）

1. **全収束次数が理論予測と一致** ✅
   - 全ペアで Rate(T) が理論値と 0.02 以内の精度で一致

2. **成分別フィッティングの高精度** ✅
   - 全成分で相対フィット誤差 < 1.1%

3. **SCPT と CT 基準の相殺比の一致** ✅
   - 差 < 0.002

4. **シミュレータコードにバグなし** ✅

## Iteration 23 で発見された分析上の改善点

### 改善点1（軽微）：正規化係数の系統的ドリフト
- 平均値ベースフィッティングが大きい dt のデータにバイアスされる
- T(SCPT,ex)/dt² のドリフト幅 1.2%
- 修正：高次フィッティングと Richardson 外挿

### 改善点2（軽微）：誤差ベクトル幾何学の未定量化
- 相殺比の非単調性の原因が定性的説明のみ
- 修正：Frobenius ノルムベースの角度推定

### 改善点3（分析ギャップ）：漸近 Rate(T) の未確認
- ST vs Exact の Rate(T) = 0.91 < 1.0（dt=0.5→0.2）
- 理論的には dt→0 で Rate→1.0
- 修正：dt=0.1 の追加

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260225_iteration24.md` を作成
- 改善点3件を数値的根拠とともに記載
- 誤差ベクトル角度の推定値を計算（θ ≈ 119°〜131°）

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration24.py` を作成
- Test A-D: Iteration 23 から変更なし
- Test E: 以下の分析改善を追加
  - **dt 範囲拡張**: n_steps = [5, 10, 20, 50, 100]（dt=0.1 を追加）
  - **Part 6b'**: 高次フィッティング T = c₀·dt^n + c₁·dt^(n+1)
  - **Part 6d**: 誤差ベクトル角度推定
  - **Part 6e**: Richardson 外挿
  - Part 6a-c: Iteration 23 と同一

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration24.py
```
結果は `developing/verification_results/iteration24_asymptotic_analysis_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 24 の結果を確認し、以下を検証:
  1. dt=0.1 での ST vs Exact の Rate(T) が 0.95〜1.05 の範囲内であること
  2. 高次フィッティングの c₀ が平均値と近似的に一致すること
  3. 誤差ベクトル角度が dt に対してほぼ一定であること
  4. Richardson 外挿値が高次フィッティングの c₀ と一致すること

## 期待される結果

### dt=0.1 での追加データ
- T(ST,SCPT) ≈ 1.985e-05
- T(SCPT,ex) ≈ 3.069e-06
- T(SCPT,CT) ≈ 1.041e-08
- Stinespring/Strang 比 ≈ 6.47
- Rate(T) for ST vs Exact (dt=0.2→0.1): 0.95〜1.05

### 高次フィッティングの漸近係数
- a₀_stine ≈ 1.985e-04（平均値 1.984e-04 と < 0.1% の差）
- b₀_strang ≈ 3.067e-04（平均値 3.080e-04 と < 0.5% の差）
- c₀_palindromic ≈ 1.041e-06（平均値 1.035e-06 と < 0.6% の差）

### 誤差ベクトル角度
- θ ≈ 119°〜131°（大きい dt で小さく、小さい dt で安定）
- cos(θ) ≈ -0.49〜-0.66

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論：Iteration 25 以降で実施

- 技術的には可能
- Iteration 24 の結果を含めてからノートブックを更新する方が効率的
- 詳細は `developing/検証結果分析_20260225_iteration24.md` を参照

## 参照ファイル
- `developing/検証結果分析_20260225_iteration24.md` — 詳細分析
- `developing/検証結果分析_20260225_iteration23.md` — iteration 23 の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration23.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration24.py` — 検証スクリプト
- `tutorials/qudit_gksl_simulator.py` — シミュレータ（変更なし）
