# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 19

## 基本情報

- 日時: 2026-02-24
- 前回: Iteration 18（正確な Liouvillian 比較の追加）
- 目的: Iteration 18 の検証方法論の欠陥3件を修正

## 前回（Iteration 18）からの引き継ぎ

### 確認済み事項
- ✅ シミュレータコード（8ファイル）にバグなし
- ✅ Test A（ノイズ vs Trotter 分離）正常
- ✅ Test B（ノイズバジェット分解）正常
- ✅ Test C（ε(n) 追跡）正常 — dt 変動量のハードコード修正済み
- ✅ Test D（Qubit忠実度追跡）正常
- ✅ Test E の ClassicalGKSLSimulator 参照解は正確
- ✅ Test E の密度行列診断（Tr、正値性等）は正常

### 発見された問題
1. **重大**: Stinespring近似誤差の分離方法（`infid_stine = infid_total - infid_trotter`）が数学的に不正確
   - 忠実度不忠実度は加法的でない
   - "Trotter-only" 測定が Trotter 誤差だけを測定していない
2. **重大**: 忠実度ベースの収束次数 O(dt²) が真のトレース距離収束次数 O(dt) を隠蔽
   - 1-F ≈ T² の関係により見かけの収束次数が2倍
   - N_S1 差の収束次数 ≈ 1.0〜1.5 がこれを裏付け
3. **中程度**: Trotter 収束次数の異常値（3.51）が未分析

### 詳細分析ドキュメント
- `developing/検証結果分析_20260224_iteration19.md`

## 作業手順

### ステップ 1: 分析ドキュメント作成 ✅
- `developing/検証結果分析_20260224_iteration19.md` を作成
- 3件の問題を詳細に記述（数値的証拠付き）

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration19.py` を作成
- 変更内容:
  - **Test A〜D**: Iteration 18 と同一（変更なし）
  - **Test E（修正）**:
    1. `trace_distance()` 関数の追加（T(ρ,σ) = ½||ρ−σ||₁）
    2. 古典 Trotter 参照解の追加
       - `_build_dissipator_superoperator()`: 散逸子超演算子 L_D の構築
       - `_classical_trotter_simulate()`: 古典 Trotter シミュレーション
         - Hamiltonian 部分: 81×81 のユニタリ expm（高速）
         - 散逸子部分: expm_multiply による正確な exp(L_D·dt)
    3. 三方比較の実装:
       - ST vs Exact: 全体誤差
       - ClassicalTrotter vs Exact: H-D Trotter 分割誤差
       - ST vs ClassicalTrotter: Stinespring + 個別チャネル Trotter 誤差
    4. 各比較で忠実度不忠実度(1-F)とトレース距離(T)の両方を報告
    5. 収束次数を両メトリクスで算出
    6. 1-F ≈ T² 関係の検証（(1-F)/T² 比の計算）
    7. N_S1 収束分析の追加
    8. 不正確な減算ベース分離の廃止

### ステップ 3: ユーザーによるスクリプト実行
- ユーザーがローカルで `run_tta_uc_gksl_verification_iteration19.py` を実行
- 結果は `developing/verification_results/` に自動保存

### ステップ 4: 結果確認と修正
- Iteration 19 の結果を分析
- トレース距離と忠実度の収束次数比が ≈ 2 であることを確認
- 古典 Trotter と ST の誤差比較から支配的な誤差源を特定

## コード変更方針

### 変更するファイル
- `tutorials/run_tta_uc_gksl_verification_iteration19.py`（新規作成）

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

### トレース距離と忠実度の関係

Fuchs–van de Graaf 不等式:
```
1 − √F ≤ T ≤ √(1−F)

整理すると: T² ≤ 1−F ≤ 2T − T²
```

微小誤差（F → 1）の場合、下界が支配的:
```
1−F ≈ T²
```

これにより：
- トレース距離が O(dt^α) で収束 → 忠実度不忠実度は O(dt^{2α}) で収束
- α=1（Stinespring 近似）→ 忠実度は O(dt²)、トレース距離は O(dt)

### 古典 Trotter 分割

古典 Trotter ステップ:
```
exp(L·dt) ≈ exp(L_H · dt/2) · exp(L_D · dt) · exp(L_H · dt/2)
```

- L_H: Hamiltonian 部分の超演算子
- L_D: 散逸子部分の超演算子（全チャネル合算）

Stinespring+Trotter ステップ:
```
ρ → exp(-iH dt/2) ρ exp(iH dt/2)
  → E_1(ρ)  →  E_2(ρ)  → ... → E_26(ρ)
  → exp(-iH dt/2) ρ exp(iH dt/2)
```

違い:
- 古典 Trotter: 散逸子を一括処理（exp(L_D dt)）
- ST: 各チャネルを個別に Stinespring で逐次適用

ST vs ClassicalTrotter の差 = Stinespring 近似誤差 + 個別チャネル Trotter 誤差
