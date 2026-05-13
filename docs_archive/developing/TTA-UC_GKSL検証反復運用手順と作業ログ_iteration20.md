# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 20

## 基本情報

- 日時: 2026-02-24
- 前回: Iteration 19（トレース距離の追加、古典 Trotter 参照解の追加）
- 目的: Iteration 19 の検証結果における3件の問題を修正

## 前回（Iteration 19）からの引き継ぎ

### 確認済み事項
- ✅ シミュレータコード（8ファイル）にバグなし
- ✅ Test A（ノイズ vs Trotter 分離）正常
- ✅ Test B（ノイズバジェット分解）正常
- ✅ Test C（ε(n) 追跡）正常
- ✅ Test D（Qubit忠実度追跡）正常
- ✅ 収束次数（トレース距離ベース）が理論と整合:
  - ST vs Exact: Rate(T) → 1.00（Stinespring O(dt) 累積誤差）
  - CT vs Exact: Rate(T) ≈ 2.01（2次 Strang 分割 O(dt²)）
  - ST vs CT: Rate(T) ≈ 0.99（Stinespring 支配的）
- ✅ N_S1 ポピュレーション収束も理論と整合

### 発見された問題
1. **重大**: `quantum_fidelity()` の数値精度喪失
   - CT vs Exact で F=1.0 にクリップ（n_steps≥20）
   - infid_CT_vs_exact = 0.0（真値は ~10⁻⁸〜10⁻⁹）
   - 収束次数が NaN に（3区間中2区間が計算不能）
2. **重大**: (1-F)/T² ≈ 1 の期待が混合状態に対して不正確
   - 実測値: 29.56〜179.16（1.0 とは程遠い）
   - 1-F = T² は純粋状態のみで成立する等式
   - 混合状態では 1-F > T²（Fuchs-van de Graaf 不等式の上界）
3. **中程度**: NaN 値の未処理
   - CT vs Exact の忠実度ベース収束次数で NaN 発生

### 詳細分析ドキュメント
- `developing/検証結果分析_20260224_iteration20.md`

## 作業手順

### ステップ 1: 分析ドキュメント作成 ✅
- `developing/検証結果分析_20260224_iteration20.md` を作成
- 3件の問題を詳細に記述（数値的証拠付き）

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration20.py` を作成
- 変更内容:
  - **Test A〜D**: Iteration 19 と同一（変更なし）
  - **Test E（修正）**:
    1. 忠実度精度限界の検出フラグ追加（F ≥ 1.0 にクリップされた場合）
    2. Bures 距離 d_B = √(2(1-√F)) の追加
    3. 各状態の純度（purity = Tr[ρ²]）の追加
    4. (1-F)/T² 列を削除し、d_B/T 比に置換
    5. 収束次数解析のトレース距離ベースへの一本化
    6. NaN 発生時の注釈を Markdown レポートに追加

### ステップ 3: ユーザーによるスクリプト実行
- ユーザーがローカルで `run_tta_uc_gksl_verification_iteration20.py` を実行
- 結果は `developing/verification_results/` に自動保存

### ステップ 4: 結果確認と修正
- Iteration 20 の結果を分析
- Bures 距離と純度の関係を確認
- 収束次数の最終確認

## コード変更方針

### 変更するファイル
- `tutorials/run_tta_uc_gksl_verification_iteration20.py`（新規作成）

### 変更しないファイル（シミュレータコア）
- `tutorials/qudit_gksl_simulator.py` — 変更なし
- `tutorials/qudit_gksl_noisy_simulator.py` — 変更なし
- `tutorials/qubit_gksl_simulator.py` — 変更なし
- `tutorials/qubit_gksl_noisy_simulator.py` — 変更なし
- `tutorials/stinespring_utils.py` — 変更なし
- `tutorials/gksl_math_utils.py` — 変更なし
- `tutorials/gksl_physical_parameters.py` — 変更なし
- `tutorials/classical_gksl_simulator.py` — 変更なし

## 理論的背景

### 忠実度とトレース距離の関係（混合状態の正しい理解）

Fuchs–van de Graaf 不等式:
```
T(ρ,σ)² ≤ 1 - F(ρ,σ) ≤ 2T(ρ,σ) - T(ρ,σ)²
```

- **純粋状態**（少なくとも一方が純粋）: 下界が等号 → T² = 1-F
- **混合状態**（両方が混合）: 下界は等号にならない → T² < 1-F
  - (1-F)/T² > 1 であり、状態の混合度に依存して大きな値をとる

### Bures 距離

Bures 距離: d_B(ρ,σ) = √(2(1-√F(ρ,σ)))

性質:
- 三角不等式を満たす距離メトリクス
- 小さい不忠実度で: d_B ≈ √(1-F)
- Fuchs-van de Graaf: d_B²/2 ≤ T ≤ d_B·√(1 - d_B²/4)
- 混合状態では d_B >> T となりうる

### 純度と (1-F)/T² 比の関係

参照状態 ρ の純度 Tr[ρ²] が低い（より混合的な）ほど、微小摂動 δ に対して:
- 忠実度は ρ の小さい固有値の逆数に比例して敏感になる（1-F が大きくなる）
- トレース距離は固有値構造に依存しない（T = ½||δ||₁）
- 結果として (1-F)/T² は混合状態で大きくなる
