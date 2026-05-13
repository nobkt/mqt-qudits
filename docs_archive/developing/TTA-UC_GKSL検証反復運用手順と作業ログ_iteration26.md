# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 26

## 日時
2026-02-25

## 目的
Iteration 25 の検証結果を分析した結果、合成モデル（Part 6f-6h）にノルム選択に起因する精度問題を発見。トレース距離ベースの合成モデル（最大誤差 8.7%、Rate(T) 予測 Δ=-0.22）を Frobenius 距離ベース（最大誤差 1.2%、Rate(d_F) 予測 |Δ|<0.016）に再構築し、T/d_F 比率分析を追加する。

## Iteration 25 で確認された成果

1. **全個別成分の収束次数が理論と完全一致** ✅
2. **Frobenius 角度 θ_F が極めて安定** ✅
   - θ_F = 125.0°〜125.2°（変動 0.2°）
   - トレース距離角度 θ_T = 119.4°〜131.0°（変動 11.6°）
   - → Frobenius ノルムでの余弦定理が厳密に成立することが実証
3. **シミュレータコードにバグなし** ✅
4. **T/d_F 比率が個別成分で極めて安定** ✅（独自分析で発見）
   - Stinespring: T/d_F ≈ 0.8911（変動 0.031%）
   - Strang: T/d_F ≈ 0.9118（変動 0.047%）

## Iteration 25 で発見された分析上の問題点

### 問題点1（重大）：合成モデルがトレース距離角度を使用
- 原因: 余弦定理をトレース距離（L₁ ノルム）に適用 → cos(θ_T) が dt 依存
- 結果: T_model の最大誤差 8.69%、Rate(T) 予測の最大 |Δ| = 0.22
- 修正: Frobenius 距離ベースの合成モデルに再構築

### 問題点2（分析不足）：T/d_F 比率の dt 依存性が未分析
- 合成 T/d_F は 0.81〜0.96 の範囲で変動（18%）
- 原因: Stinespring (T/d_F=0.89) と Strang (T/d_F=0.91) の異なる固有値分布
- 修正: T/d_F 比率の分析を追加

### 問題点3（出力不足）：Frobenius 係数フィッティングが未出力
- a_F, b_F, r_F が未出力
- Frobenius 合成モデルとその精度が未出力
- Rate(d_F) の予測と実測の比較が未出力

### 問題点4（軽微）：変数存在チェックの実装
- `"a0_stine" in dir()` の不適切な使用
- 辞書チェックに置き換え

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260225_iteration26.md` を作成
- 問題点4件を数値的根拠とともに記載
- Frobenius 合成モデルの優位性を独自分析で実証

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration26.py` を作成
- Test A-D: Iteration 25 から変更なし
- Test E の Part 1-5, 6a-6e, 6g: 変更なし
- Test E の修正・追加部分:
  - **Part 6f 再構築**: Frobenius 距離ベースの合成モデル
    - Frobenius 係数フィッティング (a_F, b_F)
    - d_F_model(dt) = a_F·dt·√(1 + r_F²dt² + 2r_F·dt·cos(θ_F_mean))
    - トレース距離モデルとの精度比較
  - **Part 6h 再構築**: Frobenius ベース Rate(d_F) 予測
    - Rate(d_F) のモデル予測と実測の比較
    - Rate(T) との比較で T/d_F 効果を定量化
  - **Part 6i（新規）**: T/d_F 比率分析
    - 個別成分 T/d_F の安定性
    - 合成 T/d_F の dt 依存性
  - **Part 6j（新規）**: Frobenius 距離の収束次数 Rate(d_F)
  - **変数チェック修正**: `dir()` → 辞書キーチェック

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration26.py
```
結果は `developing/verification_results/iteration26_frobenius_model_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 26 の結果を確認し、以下を検証:
  1. Frobenius 合成モデルの精度が全 dt で < 2% であること
  2. Rate(d_F) 予測の |Δ| が全区間で < 0.02 であること
  3. T/d_F 比率の安定性分析が想定通りであること
  4. Rate(d_F) = 0.90 と Rate(T) = 0.87 の差が T/d_F 効果で説明されること

## 期待される結果

### Frobenius 合成モデル
- 全 dt で相対誤差 < 2%（独自分析で確認済み: 最大 1.17%）
- dt≤1.0 で相対誤差 < 0.1%

### Rate(d_F) 予測
- 全区間で |Δ| < 0.02（独自分析で確認済み: 最大 0.016）
- dt≤1.0 で |Δ| < 0.001

### T/d_F 比率
- 個別成分: 変動 < 0.05%（Stinespring: 0.031%, Strang: 0.047%）
- 合成: 変動 ~18%（Stinespring/Strang のクロスオーバーによる）
- Rate(T) と Rate(d_F) の差が T/d_F の対数微分で説明可能

## quantum_dynamics_gksl_comparison.ipynb の改修検討

### 結論：Iteration 26 の結果確認後に実施

- 個別成分の収束グラフは現時点でも追加可能だが、Frobenius 合成モデルの可視化を含めて一括改修する方が効率的
- Iteration 26 で Frobenius モデルの精度が確認された後にノートブックを改修
- 追加すべきコンテンツ:
  - 個別成分の収束グラフ
  - Frobenius 合成モデル vs トレース距離モデルの比較
  - T/d_F 比率の可視化
  - Rate(d_F) vs Rate(T) の比較

## 参照ファイル
- `developing/検証結果分析_20260225_iteration26.md` — 詳細分析
- `developing/検証結果分析_20260225_iteration25.md` — iteration 25 の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration25.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration26.py` — 検証スクリプト
- `tutorials/qudit_gksl_simulator.py` — シミュレータ（変更なし）
