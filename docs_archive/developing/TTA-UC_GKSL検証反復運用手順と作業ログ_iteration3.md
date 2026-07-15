# TTA-UC GKSL検証 反復運用手順と作業ログ（第3回反復）

## 作業日時
2026-02-22

## 作業概要
第2回反復の検証結果（20260222T091949Z: 全11シナリオPASS）を受け、コード全体の
理論・設計・仕様を精査し、構造的問題を発見・修正した。

## 発見した問題点

### Issue 1: ClassicalGKSLSimulator ≡ QuditGKSLSimulator（独立参照解の消失）
- 第2回反復でRK45→Stinespring+Trotter化した結果、classical/qubit/quditが完全同一結果
- 独立したクロスバリデーションとして機能しなくなっていた

### Issue 2: `build_gksl_superoperator` のKronecker積引数順序バグ
- column-major (F-order) ベクトル化に対して、row-major (C-order) のKronecker積を使用
- 全ての `np.kron` の引数が逆転していた（潜在バグ、未使用コードパス）

### Issue 3: `build_trotter_step_classical` の同様のバグ
- Issue 2と同一の問題

### Issue 4: ClassicalGKSLBosonSimulator のBDF ODE積分
- solve_ivp(BDF) はCPTP構造を保存しない
- 小システムでは顕在化していなかったが、構造的に不適切

### Issue 5: ボソンシミュレータのrho_final不整合
- ClassicalGKSLBosonSimulator: ρ_el (9×9) を返却
- QuditGKSLBosonSimulator: ρ_total (36×36) を返却
- QubitGKSLBosonSimulator: ρ_total (36×36) を返却
- 忠実度比較不可能、検証エントロピー不一致の原因

## 実施した修正

### 1. `build_gksl_superoperator` 修正 (stinespring_utils.py)
- 全Kronecker積引数をcolumn-major規約に修正
- `np.kron(H, I)` → `np.kron(I, H)` 等
- 2×2系での検証で正確性を確認

### 2. `build_trotter_step_classical` 修正 (stinespring_utils.py)
- 同様にcolumn-major規約に修正

### 3. ClassicalGKSLSimulator → exact Liouvillian expm (classical_gksl_simulator.py)
- `build_gksl_superoperator` でGKSLリウビリアンL構築
- `scipy.sparse.linalg.expm_multiply(L, vec(ρ₀), ...)` で効率的に計算
- Stinespring+Trotterとは根本的に異なるアルゴリズム → 独立参照解
- CPTP保存はQDS理論により数学的に保証
- 性能: 6561×6561超演算子に対して約3秒（full expmの225秒から70倍高速化）

### 4. ClassicalGKSLBosonSimulator → exact Liouvillian expm (classical_gksl_boson_simulator.py)
- solve_ivp(BDF) → exact Liouvillian expm 方式に変更
- 拡張電子+フォノン空間(1296×1296超演算子)で計算、約1.7秒
- rho_finalはρ_el(電子系のみ)を返却

### 5. QuditGKSLBosonSimulator rho_final修正 (qudit_gksl_boson_simulator.py)
- `rho_final: rho` → `rho_final: partial_trace_phonon(rho, dim_el, dim_ph)`

### 6. QubitGKSLBosonSimulator rho_final修正 (qubit_gksl_boson_simulator.py)
- 同様にρ_el返却に統一

## 検証結果

### ClassicalGKSLSimulator (exact Liouvillian expm)
- Valid: True
- Trace: 1.000000000000001
- Min eigenvalue: -2.596e-16 (machine precision)
- Entropy: 0.030569978935907
- 旧Stinespring+Trotter結果(0.030047)との差: ~1.7% = Trotter誤差

### ClassicalGKSLBosonSimulator (exact Liouvillian expm)
- Valid: True
- rho_final shape: (9, 9) ← 修正前は9×9(classical)と36×36(qubit/qudit)で不整合
- 時系列/検証エントロピー一致: True

### Qudit/QubitGKSLBosonSimulator (rho_final修正)
- rho_final shape: (9, 9) ← 修正前は36×36
- 時系列/検証エントロピー一致: True

## 次回反復で必要な作業
1. ユーザーがローカルでrun_tta_uc_gksl_verification.pyを実行
2. 検証結果を確認して残存問題があれば修正

## 修正ファイル一覧
- `tutorials/stinespring_utils.py`: build_gksl_superoperator, build_trotter_step_classical
- `tutorials/classical_gksl_simulator.py`: exact Liouvillian expm方式
- `tutorials/classical_gksl_boson_simulator.py`: exact Liouvillian expm方式
- `tutorials/qudit_gksl_boson_simulator.py`: rho_final修正
- `tutorials/qubit_gksl_boson_simulator.py`: rho_final修正
- `developing/検証結果分析_20260222_iteration3.md`: 分析レポート
