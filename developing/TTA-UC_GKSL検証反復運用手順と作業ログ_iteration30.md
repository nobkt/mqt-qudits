# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 30

## 日時
2026-02-26

## 目的
Iteration 29 の検証結果分析により、残存誤差のボトルネックが Stinespring 係数モデル a_F(dt) = a_F0 + a_F1·dt の限界であることが判明した。d_F(ST,SCPT)/dt は dt ≈ 0.8 で極大を持つ非単調パターンを示し、2パラメータ線形モデルでは最大 0.113% の成分誤差が残る。3パラメータモデル a_F(dt) = a_F0 + a_F1·dt + a_F2·dt² は LOO 交差検証で過剰適合ではないことが確認され、合成モデル最大誤差を 0.074% → ~0.004%（20 倍改善）に削減できると予測される。

## Iteration 29 で確認された成果

1. **Strang dt² モデル** ✅ — 成分最大誤差 0.003%（58 倍改善）
2. **M5 合成モデル** ✅ — 最大誤差 0.074%（目標 < 0.08% 達成）
3. **cos(θ_F)(dt²) モデル** ✅ — cos 近似は本質的に解決
4. **Rate(d_F) 予測 M5** ✅ — max|Δ| = 0.0009

## Iteration 29 で発見された問題点

### 問題点: Stinespring 係数 a_F(dt) の非単調パターン

**症状**:
- M5 誤差が小さい dt で増大（dt=0.1 で 0.074% vs dt=2.0 で 0.008%）
- 誤差分解により、dt ≤ 0.5 では Stinespring 成分がほぼ 100% を支配

**原因**:
- a_F(dt) = d_F(ST,SCPT)/dt は dt ≈ 1.0 で極大（非単調）
- 変動幅: 0.279%
- 線形モデル a_F0 + a_F1·dt は単調関数のみ表現可能

**解決策**:
- 3パラメータモデル: a_F(dt) = a_F0 + a_F1·dt + a_F2·dt²
- 物理的動機: Stinespring 誤差 Taylor 展開の 3次項まで
- LOO 最大予測誤差: 0.0076%（過剰適合なし）

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_iteration29_結果分析.md` を作成
- Stinespring 非単調パターンの発見と分析
- 3パラメータモデルの LOO 交差検証
- M7 合成モデルの予測精度

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration30.py` を作成
- Iteration 29 からの変更点:
  - **追加 dt ポイント**: dt=0.25 (n=40), dt=0.125 (n=80)
  - **Part 6f 改善: 3パラメータ Stinespring モデル**
    - a_F(dt) = a_F0 + a_F1·dt + a_F2·dt² のフィッティング
    - LOO 交差検証の出力
  - **Part 6f 改善: M7 モデルバリアント追加**
    - M7: Stine-3param + Strang-dt² + cos(dt²)
    - 7 モデル比較表に拡張
  - **Part 6h 改善: M7 Rate(d_F) 予測**
  - Test A-D: 変更なし

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration30.py
```
結果は `developing/verification_results/iteration30_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 30 の結果を確認し、以下を検証:
  1. 3パラメータ Stinespring モデルの成分フィット最大誤差 < 0.002%
  2. M7 合成モデルの全 dt で誤差 < 0.01%
  3. LOO 交差検証の最大予測誤差 < 0.02%（新データ含む）
  4. M7 Rate(d_F) 予測の max|Δ| < 0.0005
  5. 追加 dt ポイントでの一貫性確認

## 期待される結果

### 3パラメータ Stinespring モデル
- 成分フィット最大誤差: ~0.001%
- LOO 最大予測誤差: ~0.01% 以下
- a_F0, a_F1, a_F2 の安定性（追加データで確認）

### 合成モデル M7（推奨）
- 全 dt で相対誤差 < 0.01%
- iteration 29 M5（最大 0.074%）の ~15-20 倍改善

### Rate(d_F) 予測
- max|Δ| の改善（M5 の 0.0009 以下）

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論: Iteration 30 の結果確認後に実施可能

- M7 モデルが確認されれば、ノートブックに一括追加:
  - SCPT 比較の追加
  - Frobenius 収束分析の可視化
  - cos(θ_F) の dt 依存性グラフ
  - 合成モデルバリアント比較のグラフ
  - 誤差分解の可視化（Stinespring 支配度クロスオーバー）

## 参照ファイル
- `developing/検証結果分析_iteration29_結果分析.md` — 今回の詳細分析
- `developing/検証結果分析_20260226_iteration29.md` — 前回の分析（iteration 29 準備）
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration29.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration29.py` — 前回の検証スクリプト
- `tutorials/run_tta_uc_gksl_verification_iteration30.py` — 今回の検証スクリプト
