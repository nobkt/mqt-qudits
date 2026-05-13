# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 27

## 日時
2026-02-26

## 目的
Iteration 26 の検証結果を分析した結果、Frobenius 合成モデルの係数フィッティングに4つの問題を発見（精度改善の余地、一貫性欠如、検証不足、データ不足）。高次フィッティングの導入、非パラメトリック検証の追加、Frobenius Richardson 外挿の追加、および不足 Frobenius 距離ペアの追加により、モデル精度を 1.13% → 0.23% に改善する。

## Iteration 26 で確認された成果

1. **Frobenius 合成モデルの基本的妥当性** ✅
   - トレース距離モデル（最大誤差 8.69%）→ Frobenius モデル（最大誤差 1.13%）
   - Rate(d_F) 予測の最大 |Δ| = 0.017（トレース距離: 0.222）
2. **全個別成分の収束次数が理論と完全一致** ✅
3. **Frobenius 角度 θ_F が極めて安定** ✅（125.0°〜125.2°、変動 0.2°）
4. **T/d_F 比率分析が正常** ✅（個別成分: 変動 < 0.05%、合成: クロスオーバーで変動）
5. **シミュレータコードにバグなし** ✅

## Iteration 26 で発見された問題点

### 問題点1（精度改善）：Frobenius 係数が単純平均
- Part 6f が `a_F = mean(d_F(ST,SCPT)/dt)` を使用
- Part 6b' のトレース距離では高次フィッティング (a₀+a₁·dt) を使用
- 高次フィッティングで合成モデル最大誤差: 1.13% → 0.23%（4.8x 改善）

### 問題点2（一貫性欠如）：高次フィッティングと Richardson 外挿の不整合
- トレース距離: Part 6b'（高次）+ Part 6e（Richardson）✅
- Frobenius 距離: ❌ 未実装

### 問題点3（検証不足）：非パラメトリック余弦定理検証の欠如
- 余弦定理の数値的検証が未出力
- 独自分析で機械精度（~10⁻¹⁶）の成立を確認済み

### 問題点4（データ不足）：Frobenius 距離の計測ペアが不完全
- frob_CT_vs_exact, frob_SCPT_vs_CT が未計測
- Part 6j の Rate(d_F) 分析が不完全

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260226_iteration27.md` を作成
- 4問題点を数値的根拠とともに記載
- 高次フィッティングの改善効果を独自分析で定量化

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration27.py` を作成
- Test A-D: Iteration 26 から変更なし
- Test E の修正・追加部分:
  - **Frobenius 距離ペアの拡張**: frob_CT_vs_exact, frob_SCPT_vs_CT を追加
  - **Part 6f 改善: 高次 Frobenius 係数フィッティング**
    - a_F0 + a_F1·dt, b_F0 + b_F1·dt の最小二乗法
    - 高次モデル d_F_model_ho(dt) の精度評価
    - 単純モデルとの比較
  - **Part 6f 改善: 非パラメトリック余弦定理検証**
    - 測定値からの直接再構成
  - **Part 6f 改善: Frobenius Richardson 外挿**
    - a_F_Rich, b_F_Rich の計算
  - **Part 6h 改善: 高次モデルベースの Rate(d_F) 予測**
  - **Part 6j 拡張: CT ペアの Rate(d_F) 追加**

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration27.py
```
結果は `developing/verification_results/iteration27_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 27 の結果を確認し、以下を検証:
  1. 高次 Frobenius 合成モデルの全 dt で誤差 < 0.5% であること
  2. Rate(d_F) 予測の |Δ| が全区間で < 0.01 であること
  3. 非パラメトリック余弦定理検証が機械精度であること
  4. Frobenius Richardson 外挿が高次フィッティングと整合すること
  5. CT ペアの Rate(d_F) が理論値 2.0 に一致すること

## 期待される結果

### 高次 Frobenius 合成モデル
- 全 dt で相対誤差 < 0.3%（独自分析で確認済み: 最大 0.23%）
- 単純平均モデル（最大 1.13%）の 5x 改善

### 高次 Rate(d_F) 予測
- 全区間で |Δ| < 0.01（独自分析で確認済み: 最大 0.006）
- 単純モデル（最大 0.017）の 3x 改善

### 非パラメトリック検証
- ~10⁻¹⁶ の精度で余弦定理成立

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論：Iteration 27 の結果確認後に実施

- 高次 Frobenius モデルが確認されれば、ノートブックに以下を一括追加:
  - 個別成分の収束グラフ
  - Frobenius 合成モデルの可視化
  - T/d_F 比率の可視化
  - Rate(d_F) vs Rate(T) の比較

## 参照ファイル
- `developing/検証結果分析_20260226_iteration27.md` — 詳細分析
- `developing/検証結果分析_20260225_iteration26.md` — 前回の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration26.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration26.py` — 前回の検証スクリプト
- `tutorials/run_tta_uc_gksl_verification_iteration27.py` — 今回の検証スクリプト
