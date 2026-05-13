# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 28

## 日時
2026-02-26

## 目的
Iteration 27 の検証結果を分析した結果、高次 Frobenius 合成モデルに 3 つの問題を発見（HO モデルの小 dt 劣化、cos(θ_F) 定数近似、誤差分解診断の欠如）。cos(θ_F)(dt) = c₀ + c₁·dt² モデリングの導入、モデルバリアント体系比較の追加、および誤差分解診断の追加により、全 dt で最大誤差を 0.23% → 0.15% に改善する。

## Iteration 27 で確認された成果

1. **高次 Frobenius 係数フィッティング** ✅
   - Stinespring: a_F0 = 2.228760e-04, a_F1 = -2.289e-07（成分最大誤差 0.11%）
   - Strang: b_F0 = 3.359825e-04, b_F1 = 2.174e-06（成分最大誤差 0.18%）
2. **非パラメトリック余弦定理検証** ✅（最大誤差 1.86e-16）
3. **Frobenius Richardson 外挿** ✅（高次フィッティングと整合）
4. **CT ペアの Frobenius 距離と Rate(d_F)** ✅（理論値 2.0 と一致）
5. **シミュレータコードにバグなし** ✅

## Iteration 27 で発見された問題点

### 問題点1：HO モデルが小さい dt で単純モデルより悪化
- dt=0.1 で HO 誤差 0.099% vs 単純 0.006%（16.5x 悪い）
- dt=1.0 で HO 誤差 0.159% vs 単純 0.019%（8.5x 悪い）
- 原因: HO 係数フィットの残差と cos_mean 近似の非最適な相互作用

### 問題点2：cos(θ_F) が定数（平均値）で近似
- cos(θ_F) は dt=2.0 で -0.5730、dt→0 で -0.5767（変動 0.64%）
- cos(θ_F)(dt) = c₀ + c₁·dt² で正確にモデル化可能（最大誤差 9.12e-06）
- cos(θ_F) の dt² モデルで合成モデル最大誤差: 0.23% → 0.15%

### 問題点3：誤差分解診断の欠如
- 係数フィッティング誤差と cos(θ_F) 近似誤差の個別寄与が不明
- 分析結果: 係数誤差が主要な残存誤差源（最大 0.15%）
- cos(θ_F) 誤差は dt=2.0 で最大（0.12%）

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260226_iteration28.md` を作成
- 3 問題点を数値的根拠とともに記載
- cos(θ_F) dt² モデリングの有効性を独自分析で定量化
- 4 モデルバリアントの精度比較を実施

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration28.py` を作成
- Test A-D: Iteration 27 から変更なし
- Test E の修正・追加部分:
  - **Part 6f 改善: cos(θ_F)(dt) モデリング**
    - cos(θ_F) = c₀ + c₁·dt² の最小二乗フィッティング
    - 「完全 HO」モデル: HO 係数 + cos(θ_F)(dt²) の合成モデル
  - **Part 6f 改善: 4 モデルバリアント体系比較**
    - 単純, HO+cos_mean, HO+cos(dt²), HO+per-dt cos
  - **Part 6f 改善: 誤差分解診断**
    - 係数誤差のみ / cos 誤差のみ / 複合誤差の分解
  - **Part 6h 改善: cos(dt²) モデルの Rate(d_F) 予測**

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration28.py
```
結果は `developing/verification_results/iteration28_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 28 の結果を確認し、以下を検証:
  1. cos(θ_F) dt² モデルのフィット精度が ~10⁻⁵ であること
  2. HO+cos(dt²) 合成モデルの全 dt で誤差 < 0.2% であること
  3. 誤差分解で係数誤差が主要源であることの確認
  4. Rate(d_F) 予測の |Δ| が < 0.006 であること
  5. 4 モデルバリアント比較が正しく出力されること

## 期待される結果

### cos(θ_F) dt² モデル
- フィット最大誤差: ~10⁻⁵（独自分析で確認済み: 9.12e-06）

### 完全 HO 合成モデル（HO 係数 + cos(dt²)）
- 全 dt で相対誤差 < 0.16%（独自分析で確認済み: 最大 0.15%）
- 全 dt で均一に良好な精度

### 誤差分解
- 係数フィッティングが主要残存誤差源（0.08-0.15%）
- cos(θ_F) モデリングは dt=2.0 で効果最大（0.12%→~0%）

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論：Iteration 28 の結果確認後に実施

- cos(θ_F) dt² モデルが確認されれば、ノートブックに以下を一括追加:
  - 個別成分の収束グラフ
  - Frobenius 合成モデルの可視化（全バリアント比較）
  - cos(θ_F) の dt 依存性グラフ
  - 誤差分解の可視化
  - T/d_F 比率の可視化

## 参照ファイル
- `developing/検証結果分析_20260226_iteration28.md` — 詳細分析
- `developing/検証結果分析_20260226_iteration27.md` — 前回の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration27.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration27.py` — 前回の検証スクリプト
- `tutorials/run_tta_uc_gksl_verification_iteration28.py` — 今回の検証スクリプト
