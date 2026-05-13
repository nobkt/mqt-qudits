# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 23

## 日時
2026-02-25

## 目的
Iteration 22 の検証結果の詳細分析に基づき、誤差分析の精度を改善する。具体的には、成分別多項式フィッティング、SCPT 基準の相殺比、正規化誤差係数の追加を行う。

## Iteration 22 で確認された成果（正常動作の確認）

1. **全収束次数が理論予測と一致** ✅
   - ST vs SCPT: Rate(T) = 1.00（純粋 Stinespring 誤差、O(dt)）
   - ST vs CT: Rate(T) = 1.00（Stinespring 誤差の近似、O(dt)）
   - CT vs Exact: Rate(T) = 2.01（Strang 分割誤差、O(dt²)）
   - SCPT vs Exact: Rate(T) = 2.01（回文順残差 + Strang、O(dt²)）
   - SCPT vs CT: Rate(T) = 1.99（回文順 Lie-Trotter 残差、O(dt²)）

2. **SCPT 実装の正しさ確認** ✅
   - ST と同一の回文順構造を正確に再現
   - T(SCPT,CT) ≈ 10⁻⁶（Stinespring 誤差の 1/200 未満）

3. **シミュレータコードにバグなし** ✅

## Iteration 22 で発見された分析上の問題

### 問題1（中程度）：多項式フィッティングの不適切性
- 集約フィット T(ST,ex) ≈ a·dt + b·dt² の最大相対誤差が 29%
- 原因：誤差相殺により T(ST,ex) ≠ T(ST,SCPT) + T(SCPT,ex)
- 集約フィット係数 a = 7.15e-05 は真の Stinespring 係数 1.98e-04 の 2.8 分の1
- 詳細：`developing/検証結果分析_20260225_iteration23.md`

### 問題2（中程度）：SCPT 基準の相殺比の欠如
- 相殺比が CT 基準のみで計算されていた
- SCPT 基準の方が物理的に明確な分解を与える

### 問題3（軽微）：T(SCPT,CT) の主テーブル不記載
- 回文順残差の大きさが直接見えなかった

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260225_iteration23.md` を作成
- 成分別フィッティングの有効性を数値的に実証
- 集約フィッティングの不適切性を定量的に示した

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration23.py` を作成
- Test A-D: Iteration 22 から変更なし
- Test E: 以下の分析改善を追加
  - **Part 6a**: 集約フィッティング（参考値として維持）
  - **Part 6b**: 成分別フィッティング（新規）
    - T(ST,SCPT) = a_stine · dt（Stinespring 誤差）
    - T(SCPT,ex) = b_strang · dt²（Strang+回文順残差）
    - T(SCPT,CT) = c_palindromic · dt²（回文順残差）
    - 誤差大きさ比の表示
  - **Part 6c**: 正規化誤差係数テーブル
  - 主テーブルに T(SCPT,CT) 列を追加
  - SCPT 基準の相殺比 cancel_SCPT を追加

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration23.py
```
結果は `developing/verification_results/iteration23_component_analysis_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 23 の結果を確認し、以下を検証:
  1. 成分別フィッティングの相対誤差が 1% 未満であること
  2. SCPT 基準と CT 基準の相殺比が近似的に一致すること
  3. 正規化誤差係数が dt に対してほぼ一定であること
  4. 集約フィッティングとの乖離が定量的に確認できること

## 期待される結果

### 成分別フィッティングの係数
- a_stine ≈ 1.98e-04（Stinespring 誤差）
- b_strang ≈ 3.08e-04（Strang+回文順残差）
- c_palindromic ≈ 1.04e-06（回文順残差のみ、b_strang の 1/300）

### 正規化誤差係数の安定性
- T(ST,SCPT)/dt ≈ 1.98e-04（±0.3%）
- T(SCPT,ex)/dt² ≈ 3.08e-04（±0.8%）
- T(SCPT,CT)/dt² ≈ 1.04e-06（±1.7%）

### SCPT 基準と CT 基準の相殺比の一致
- cancel_SCPT ≈ cancel_CT（差 < 0.002）
- 両方が 0.43〜0.67 の範囲

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論：可能（ただし iteration 24 以降で実施すべき）

- 技術的には SCPT 比較の追加が可能
- ただし検証スクリプトの分析完成を待ってから反映すべき
- 詳細は `developing/検証結果分析_20260225_iteration23.md` を参照

## 参照ファイル
- `developing/検証結果分析_20260225_iteration23.md` — 詳細分析
- `developing/検証結果分析_20260224_iteration22.md` — iteration 22 の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration22.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration23.py` — 検証スクリプト
- `tutorials/qudit_gksl_simulator.py` — シミュレータ（変更なし）
