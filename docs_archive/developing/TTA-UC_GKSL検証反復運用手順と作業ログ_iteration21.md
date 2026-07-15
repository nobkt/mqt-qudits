# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 21

## 日時
2026-02-24

## 目的
Iteration 20 の検証結果から発見された重大な問題（ST シミュレータの収束次数 O(dt) vs 文書記載の O(dt²)）に対して：
1. 根本原因を分析し、分析レポートを作成
2. コード修正（対称 Lindblad 積 + ドキュメント修正）
3. 誤差分解のための新しい診断ツールを含む検証スクリプト（iteration 21）を作成

## Iteration 20 で発見された問題

### 問題1（重大）：収束次数の不一致
- **主張**: QuditGKSLSimulator / QubitGKSLSimulator は「2nd-order symmetric Trotter」
- **実測**: トレース距離の収束率 Rate(T) ≈ 1.0（O(dt) = 1次収束）
- **対比**: Classical Trotter は Rate(T) ≈ 2.0（O(dt²) = 2次収束）
- **根本原因**:
  (a) Stinespring 拡張は各 Lindblad チャネルの1次近似（O(dt²)/ステップ）
  (b) Lie-Trotter 積（逐次適用）は合計散逸子の1次近似（O(dt²)/ステップ）

詳細分析: `developing/検証結果分析_20260224_iteration21.md`

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260224_iteration21.md` を作成
- 3段階の誤差分解（Strang 分割 / Lie-Trotter 積 / Stinespring 近似）を詳述
- 数値データによる裏付け

### ステップ 2: コード修正 ✅

#### 修正 A: ドキュメント修正
- `tutorials/qudit_gksl_simulator.py`:
  - モジュール docstring: "2nd-order Trotter" → 実効 O(dt) 収束の説明追加
  - クラス docstring: 同上
  - `_trotter_step` docstring: 各誤差源と実効収束次数を明記
- `tutorials/qubit_gksl_simulator.py`: 同様の修正

#### 修正 B: 対称（回文順）Lindblad 積の実装
- `_precompute_unitaries(dt)`:
  - `_U_stines_half`: dt/2 の Stinespring ユニタリ（対称積用）
  - `_U_stines`: dt の Stinespring ユニタリ（ノイズ付きサブクラス互換用）
- `_trotter_step(rho)`:
  - `exp(L_H dt/2) · prod_{α=1..n} E_α(dt/2) · prod_{α=n..1} E_α(dt/2) · exp(L_H dt/2)`
  - Lie-Trotter 交換子誤差を除去（誤差源 (b) → O(dt³)/ステップ）
  - Stinespring 近似誤差は残存（誤差源 (a) → O(dt²)/ステップ）

### ステップ 3: 既存テストの確認 ✅
- `test_gksl_simulators.py`: 53 テスト合格、1 テスト不合格（既存のメモリ問題、無関係）
- 主要テスト合格:
  - `test_stinespring_fidelity_qudit`: F > 0.99 ✅
  - `test_classical_qudit_convergence`: diff < 1e-3 ✅
  - `test_classical_qubit_convergence`: diff < 1e-3 ✅
  - `test_convergence_improves_with_dt`: diff 減少 ✅

### ステップ 4: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration21.py` を作成
- Test A-D: Iteration 20 から変更なし
- Test E: 以下を追加
  - **Classical Product Trotter (CPT)** 比較対象の追加
  - `_build_individual_dissipator_superoperators()`: 個別 L_{D_α} 超演算子の構築
  - `_classical_product_trotter_simulate()`: `exp(L_H dt/2) · prod_α exp(L_{D_α} dt) · exp(L_H dt/2)`
  - 6組の比較ペア:
    | 比較 | 分離される誤差 |
    |------|-------------|
    | ST vs Exact | 全誤差 |
    | CT vs Exact | Strang 分割誤差のみ |
    | CPT vs Exact | Strang + Lie-Trotter 積誤差 |
    | ST vs CPT | **Stinespring 近似誤差のみ** |
    | CPT vs CT | **Lie-Trotter 積誤差のみ** |
    | ST vs CT | Stinespring + Lie-Trotter 積 |

### ステップ 5: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration21.py
```
結果は `developing/verification_results/iteration21_error_decomposition_*.json/md` に出力される。

### ステップ 6: 結果分析 ⬜
- iteration 21 の結果を確認し、以下を検証:
  1. 対称 Lindblad 積により ST の収束が改善されたか
  2. ST vs CPT（Stinespring 誤差のみ）の収束次数
  3. CPT vs CT（Lie-Trotter 積誤差のみ）の収束次数
  4. 改善の程度（定数係数の比較）

## 期待される結果

### 対称 Lindblad 積の効果
- **CPT vs CT**: Lie-Trotter 積誤差 → Rate(T) ≈ 1.0（1次）
- **ST vs CPT**: Stinespring 近似誤差 → Rate(T) ≈ 1.0（1次）
- **ST vs Exact**: 両方の合計 → Rate(T) ≈ 1.0（1次だが定数係数が改善される可能性）

対称積により ST vs CT（iteration 20 で Rate(T) ≈ 1.0）のうち Lie-Trotter 成分が除去され、残りは Stinespring 誤差のみ。改善の程度は Stinespring 誤差と Lie-Trotter 積誤差の比率に依存。

### 真の2次収束の達成について
Stinespring 拡張（2準位アンシラ）が本質的に1次近似であるため、真の O(dt²) 収束には高次の Stinespring 拡張（3準位以上のアンシラ）が必要。これは将来の研究課題。

## 参照ファイル
- `developing/検証結果分析_20260224_iteration21.md` — 詳細分析
- `developing/検証結果分析_20260224_iteration20.md` — iteration 20 の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration20.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration21.py` — 検証スクリプト
- `tutorials/qudit_gksl_simulator.py` — 修正済みシミュレータ
- `tutorials/qubit_gksl_simulator.py` — 修正済みシミュレータ
