# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 33

## 日時
2026-03-01

## 目的
Iteration 32 の検証結果分析で発見された **Cell 27 精度ガイダンステーブルの未修正問題** を実際に修正する。

### 前提条件（Iteration 32 で確認済み）
1. **全シミュレーターが正しく動作** — 問題なし
2. **収束特性が正しい** — O(dt) 収束は Stinespring dilation の数学的性質
3. **密度行列品質** — 全 n_steps で機械精度
4. **検証スクリプト（iteration32）は正しく動作** — 全 5 チェック PASS

### Iteration 32 で発見された問題（4 件）

1. **Cell 27 精度ガイダンステーブルが未修正**（深刻度: 高）
   - iteration 32 の作業ログでは「修正済み ✅」と記載されているが、実際のノートブックファイルには反映されていなかった
   - T 列の値が t_max=10.0 のデータのまま（~3e-04, ~8e-05, ~3e-05, ~2e-05）
   - 正しい値（t_max=100.0）: ~1.3e-02, ~4e-03, ~1.5e-03, ~7e-04

2. **Cell 27 の誤った記述が未修正**（深刻度: 高）
   - 「`n_steps=100` で `T < 2×10⁻⁵`」が残存（正しくは T ≈ 7.3e-04）

3. **Cell 27 に n_steps=200 の行が欠落**（深刻度: 中）
   - Cell 26 は n_steps=200 まで計算するが、テーブルに含まれていなかった

4. **iteration31 スクリプトの Markdown レポートラベル不整合**（深刻度: 低）
   - コードの閾値は 1e-3 に修正済みだが、Markdown レポートのラベルが「T < 5e-5」のまま

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_iteration32_結果分析.md` を作成
- 4 つの問題点を詳細に文書化

### ステップ 2: Cell 27 精度ガイダンステーブルの修正 ✅
- `tutorials/quantum_dynamics_gksl_comparison.ipynb` の Cell 27 を修正
- T 列の値を t_max=100.0 の実測値（Cell 26 出力）に更新
- n_steps=200 の行を追加
- 精度レベルの表記を実測値に合わせて修正
- `t_max=100.0` の条件を明示

### ステップ 3: Cell 27 の誤った記述の修正 ✅
- 「`T < 2×10⁻⁵`」→「`T < 1×10⁻³`（`t_max=100.0` の場合）」
- t_max 依存性の説明を新規追加

### ステップ 4: iteration31 スクリプトの Markdown ラベル修正 ✅
- 251 行目のラベルを「T < 5e-5」→「T < 1e-3」に修正

### ステップ 5: Iteration 33 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration33.py` を作成
- 既存の 5 チェックに加え、Cell 27 テーブル整合性チェック（Check 6）を追加
- Cell 27 の記載値と実測値の比率が 0.5〜2.0 の範囲内であることを検証

### ステップ 6: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration33.py
```
結果は `developing/verification_results/iteration33_*.json/md` に出力される。

### ステップ 7: 結果分析 ⬜
- iteration 33 の結果を確認し、全 6 チェックが PASS であることを検証

## 期待される結果

### 全チェック PASS
1. トレース距離単調減少: PASS（変更なし）
2. 収束次数 ≈ 1.0: PASS（変更なし）
3. 忠実度単調増加: PASS（変更なし）
4. 密度行列品質: PASS（変更なし）
5. n_steps=100 で T < 1e-3: PASS（変更なし）
6. Cell 27 テーブル整合性: PASS（新規チェック — 修正後のテーブル値と実測値が一致）

## 修正の正当性

### Cell 27 修正内容の根拠
- Cell 26 の実際の出力（t_max=100.0, Hilbert空間次元=81）に合わせた値
- 同一ノートブック内でコードの出力と文書の値が一致することは基本要件
- 修正後の値は計算再現性がある
- n_steps=200 の行を追加し、Cell 26 の全出力をカバー
- t_max 依存性の説明は Stinespring 近似誤差の理論的性質（T ∝ t_max · dt）に基づく

### Markdown ラベル修正の根拠
- コードの実際の閾値（1e-3）とレポートの記載が一致すべき

## 参照ファイル
- `developing/検証結果分析_iteration32_結果分析.md` — 分析レポート
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration32.md` — 前回の作業ログ
- `tutorials/quantum_dynamics_gksl_comparison.ipynb` — 修正対象（Cell 27）
- `tutorials/run_tta_uc_gksl_verification_iteration31.py` — Markdown ラベル修正
- `tutorials/run_tta_uc_gksl_verification_iteration33.py` — 修正後の検証スクリプト
