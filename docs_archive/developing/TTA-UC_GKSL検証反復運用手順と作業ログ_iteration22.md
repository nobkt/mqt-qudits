# TTA-UC GKSL 検証反復運用手順と作業ログ — Iteration 22

## 日時
2026-02-24

## 目的
Iteration 21 の検証結果から発見された問題（CPT 比較の誤差分離欠陥）を修正し、正確な誤差分解のための新しい検証スクリプト（iteration 22）を作成する。

## Iteration 21 で発見された問題

### 問題1（重大）：CPT 比較が Stinespring 誤差を正しく分離できていない
- **主張**: ST vs CPT が「Stinespring 近似誤差のみ」を分離する
- **実際**: ST（回文順）と CPT（順方向のみ）で Lindblad チャネルの適用順序が異なるため、ST vs CPT には順序差と Stinespring 差の両方が混在
- **証拠**: dt=0.2 で T(ST,CPT) = 1.65e-04 >> T(ST,CT) = 3.97e-05。ST は CPT より CT に4倍以上近い。

### 問題2（中程度）：収束次数の報告が誤解を招く
- ST vs Exact の Rate(T) が 2.05→1.82→0.91 と変動し、平均 1.59 は物理的に無意味
- 原因: Stinespring 誤差 O(dt) と Strang 誤差 O(dt²) の遷移

### 問題3（軽微）：N_S1 収束の非単調性
- |Δ_ST| が dt=0.5→0.2 で微増（符号反転による）

詳細分析: `developing/検証結果分析_20260224_iteration22.md`

## 作業内容

### ステップ 1: 分析レポート作成 ✅
- `developing/検証結果分析_20260224_iteration22.md` を作成
- CPT 比較の設計欠陥を数値的に実証
- 誤差相殺の定量分析
- Iteration 20 → 21 の改善効果の確認

### ステップ 2: 検証スクリプト作成 ✅
- `tutorials/run_tta_uc_gksl_verification_iteration22.py` を作成
- Test A-D: Iteration 21 から変更なし
- Test E: 以下を修正・追加
  - **SCPT（Symmetric Classical Product Trotter）** の追加
    - `_symmetric_classical_product_trotter_simulate()`: 回文順の個別散逸子指数関数
    - ST と同じ順序構造 → ST vs SCPT で純粋 Stinespring 誤差を分離
  - 比較ペアを8組に拡張:
    | 比較 | 分離される誤差 | 期待 Rate(T) |
    |------|-------------|-------------|
    | ST vs Exact | 全誤差 | O(dt) 漸近 |
    | CT vs Exact | Strang 分割 | O(dt²) |
    | CPT vs Exact | 順方向 Lie-Trotter + Strang | O(dt) |
    | SCPT vs Exact | 回文順残差 + Strang | O(dt²) |
    | **ST vs SCPT** | **純粋 Stinespring 誤差** | **O(dt)** |
    | SCPT vs CT | 回文順 Lie-Trotter 残差 | O(dt²) |
    | CPT vs CT | 順方向 Lie-Trotter 積誤差 | O(dt) |
    | ST vs CT | Stinespring 誤差（近似） | O(dt) |
  - 符号付き N_S1 誤差の追加
  - 誤差相殺分析（三角不等式余裕度）
  - 多項式フィッティング `T ≈ a·dt + b·dt²`

### ステップ 3: ユーザーが検証スクリプトを実行 ⬜
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration22.py
```
結果は `developing/verification_results/iteration22_scpt_decomposition_*.json/md` に出力される。

### ステップ 4: 結果分析 ⬜
- iteration 22 の結果を確認し、以下を検証:
  1. SCPT vs Exact が O(dt²) を示すか（回文順残差 = O(dt²) の確認）
  2. ST vs SCPT が O(dt) を示すか（純粋 Stinespring 誤差の確認）
  3. SCPT vs CT の大きさ（回文順残差の定量化）
  4. 多項式フィッティングの妥当性

## 期待される結果

### SCPT の効果
- **SCPT vs Exact**: O(dt²) — 回文順 Lie-Trotter 残差は O(dt³)/ステップ、Strang は O(dt³)/ステップ → 合計 O(dt²) 大域
- **ST vs SCPT**: O(dt) — 純粋 Stinespring 近似誤差、2n チャネル × O(dt²)/チャネル
- **SCPT vs CT**: O(dt²) — 回文順残差は Strang と同次数

### 誤差の大きさの序列（小さい dt）
```
T(CT, exact) ≈ T(SCPT, exact) << T(ST, CT) ≈ T(ST, SCPT) << T(CPT, CT)
```

CT と SCPT はともに O(dt²) で近い大きさ。ST の追加誤差（Stinespring）は O(dt) で CT/SCPT より大きいが、CPT の順方向 Lie-Trotter 積誤差よりは小さい（iteration 21 のデータから推定）。

## 参照ファイル
- `developing/検証結果分析_20260224_iteration22.md` — 詳細分析
- `developing/検証結果分析_20260224_iteration21.md` — iteration 21 の分析
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration21.md` — 前回の作業ログ
- `tutorials/run_tta_uc_gksl_verification_iteration22.py` — 検証スクリプト
- `tutorials/qudit_gksl_simulator.py` — シミュレータ（変更なし）
