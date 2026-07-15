# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 31

## 日時
2026-02-28

## 目的
Iteration 30 の検証結果分析で確認された以下の結論に基づき、quantum_dynamics_gksl_comparison.ipynb の改修を実施する。

### 前提条件（Iteration 30 で確認済み）
1. **全シミュレーターが正しく動作** — ClassicalGKSLSimulator, QuditGKSLSimulator, QubitGKSLSimulator 等
2. **収束特性が完全に解明** — Stinespring 誤差 O(dt), Strang 誤差 O(dt²)
3. **密度行列品質** — トレース保存・エルミート性・正定値性すべて機械精度
4. **M7 合成モデル** — 誤差を 0.011% の精度で記述（ノートブックには不要）
5. **追加の検証イテレーション不要** — コードにバグは存在しない

## Iteration 30 結果分析からの改修方針

### 実施する改修
1. **Section 11b 追加**: Stinespring+Trotter 収束特性の分析セクション
   - 理論的背景（3段階近似と各次数）の Markdown ドキュメント
   - 複数 n_steps での収束解析コード（trace distance, fidelity, rate）
   - 精度ガイダンステーブル
2. **まとめセクション更新**: 収束特性の情報を追加

### 実施しない事項
- M1〜M7 のモデリング結果の追加（検証用副産物、ノートブックには不要）
- ヒューリスティックな補正・fallback の追加（コードは数学的に正しい）
- O(dt) 収束の「修正」（Stinespring dilation の正しい数学的性質）
- さらなるパラメトリックモデルの構築

## 作業内容

### ステップ 1: 分析レポート参照 ✅
- `developing/検証結果分析_iteration30_結果分析.md` を精読
- 改修方針を確定

### ステップ 2: ノートブック改修 ✅
- `tutorials/quantum_dynamics_gksl_comparison.ipynb` に以下を追加:
  - Cell 25 (Markdown): 収束特性の理論的背景
    - Strang 対称分割: O(dt³)/step → O(dt²) global
    - 回文順序積: Lie-Trotter 交換子誤差の2次消去
    - Stinespring dilation 近似: O(dt²)/step → O(dt) global（支配的）
  - Cell 26 (Code): 収束解析
    - ClassicalGKSLSimulator を厳密解として使用
    - QuditGKSLSimulator を n_steps=10,20,50,100,200 で実行
    - trace distance と fidelity を計算
    - 収束次数 Rate(T) を計算・表示
  - Cell 27 (Markdown): 精度ガイダンステーブルと収束次数の解釈
  - Cell 39 (Markdown): まとめセクションに収束特性を追加
- セル数: 37 → 40（3セル追加）
- 既存セルのコード変更: なし

### ステップ 3: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration31.py` を作成
- 検証内容:
  1. trace distance が n_steps 増加で単調減少すること
  2. 収束次数 Rate(T) ≈ 1.0（O(dt) 収束の確認）
  3. fidelity が n_steps 増加で単調増加すること
  4. 密度行列品質が全 n_steps で保持されること
  5. n_steps=100 で T < 5e-5 の実用精度が得られること

### ステップ 4: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration31.py
```
結果は `developing/verification_results/iteration31_*.json/md` に出力される。

### ステップ 5: 結果分析 ⬜
- iteration 31 の結果を確認し、以下を検証:
  1. 全5チェックが PASS であること
  2. 収束次数が理論値 ≈ 1.0 と一致すること
  3. 精度ガイダンステーブルの値が妥当であること

## 期待される結果

### 収束解析
- trace distance: n_steps に対して O(dt) = O(t_max/n_steps) で減少
- Rate(T) ≈ 1.0（許容範囲 ±0.3）
- n_steps=100 での trace distance: ~2e-05

### 密度行列品質
- 全 n_steps で:
  - |Tr(ρ)-1| < 1e-10
  - Hermiticity error < 1e-10
  - Min eigenvalue > -1e-10

## ノートブック改修の完了基準

1. 検証スクリプトの全チェックが PASS
2. ノートブックの全セル（40セル）が正常実行可能
3. 追加セクションの内容が理論的に正確
4. ヒューリスティックな処理・fallback が一切含まれていない

## 参照ファイル
- `developing/検証結果分析_iteration30_結果分析.md` — 改修方針の根拠
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration30.md` — 前回の作業ログ
- `tutorials/quantum_dynamics_gksl_comparison.ipynb` — 改修対象
- `tutorials/run_tta_uc_gksl_verification_iteration31.py` — 検証スクリプト
