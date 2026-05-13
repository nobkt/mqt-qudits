# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 29

## 日時
2026-02-26

## 目的
Iteration 28 の検証結果を分析した結果、Strang 係数モデルの形式に理論的誤りを発見。Strang（対称）分割の対称性から、正規化係数の補正は O(dt²)（偶数次）であるべきだが、現行モデルは O(dt)（奇数次含む）で近似している。b_F(dt) = b_F0 + b_F2·dt² への修正により、成分フィット誤差を 58 倍改善（0.176% → 0.003%）し、合成モデル最大誤差を 0.15% → 0.074% に半減させる。

## Iteration 28 で確認された成果

1. **cos(θ_F)(dt²) モデリング** ✅ — 最大フィット誤差 9.12e-06
2. **M3 合成モデル** ✅ — 最大誤差 0.1532%（M2 の 0.2343% から改善）
3. **cos 近似の解決** ✅ — M3 ≈ M4（差 < 0.001%）
4. **誤差分解** ✅ — 係数フィッティングが主要残存誤差源
5. **非パラメトリック余弦定理** ✅ — 機械精度で成立

## Iteration 28 で発見された問題点

### 問題点: Strang 係数モデル形式の理論的誤り

**現行**: d_F(SCPT,ex)/dt² = b_F0 + b_F1·dt（dt 線形補正）
**修正**: d_F(SCPT,ex)/dt² = b_F0 + b_F2·dt²（dt² 偶数次補正）

**根拠**: Strang 対称分割は奇数次誤差を消去するため、正規化係数の補正は偶数次のみ。

**数値的検証**:
- dt 線形モデル: 成分最大誤差 0.176%
- dt² 二乗モデル: 成分最大誤差 0.003%（58x 改善）
- 合成モデル M5（Stine-dt + Strang-dt² + cos(dt²)）: 最大 0.074%

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260226_iteration29.md` を作成
- Strang 係数モデルの理論的誤りを記載
- 代替モデル形式の数値比較を実施
- 6 モデルバリアント（M1-M6）の精度比較

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration29.py` を作成
- Test A-D: Iteration 28 から変更なし
- Test E の修正・追加部分:
  - **Part 6f 改善: Strang 係数モデルの修正**
    - d_F(SCPT,ex)/dt² = b_F0 + b_F2·dt² のフィッティング
    - 旧モデル（dt）との比較出力
  - **Part 6f 改善: Stinespring 係数の代替モデル検証**
    - a_F0+a_F1·dt vs a_F0+a_F2·dt² の比較
  - **Part 6f 改善: 6 モデルバリアント比較**
    - M1: 単純（iter26）
    - M2: HO-dt + cos_mean（iter27）
    - M3: HO-dt + cos(dt²)（iter28）
    - M4: HO-dt + per-dt cos（参照）
    - M5: Stine-dt + Strang-dt² + cos(dt²)（iter29 推奨）
    - M6: Stine-dt² + Strang-dt² + cos(dt²)（iter29 代替）
  - **Part 6h 改善: M5, M6 の Rate(d_F) 予測**

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration29.py
```
結果は `developing/verification_results/iteration29_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 29 の結果を確認し、以下を検証:
  1. Strang dt² モデルの成分フィット誤差が 0.005% 以下であること
  2. M5 合成モデルの全 dt で誤差 < 0.08% であること
  3. M5 が M3 より全 dt で改善していること
  4. Rate(d_F) 予測の max|Δ| が M3 以下であること
  5. Stinespring dt² モデルの結果が dt モデルとの比較で出力されること

## 期待される結果

### Strang dt² モデル
- 成分フィット最大誤差: ~0.003%
- 漸近値 b_F0 が Richardson 外挿と ~0.005% で一致

### 合成モデル M5（推奨）
- 全 dt で相対誤差 < 0.08%
- iteration 28 M3（最大 0.15%）の ~2x 改善

### Rate(d_F) 予測
- max|Δ| の改善（特に 2.0→1.0 範囲）

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論: Iteration 29 の結果確認後に実施可能

- Strang dt² モデルが確認されれば、ノートブックに以下を一括追加:
  - SCPT 比較の追加
  - Frobenius 収束分析の可視化
  - cos(θ_F) の dt 依存性グラフ
  - 合成モデルバリアント比較のグラフ
  - 誤差分解の可視化

## 参照ファイル
- `developing/検証結果分析_20260226_iteration29.md` — 今回の詳細分析
- `developing/検証結果分析_20260226_iteration28.md` — 前回の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration28.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration28.py` — 前回の検証スクリプト
- `tutorials/run_tta_uc_gksl_verification_iteration29.py` — 今回の検証スクリプト
