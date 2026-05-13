# Iteration 34 検証結果分析

## 分析日時
2026-03-01

## 分析対象
1. `run_tta_uc_gksl_verification_iteration34.py` の実行結果
2. `quantum_dynamics_gksl_comparison.ipynb` の全セル実行結果
3. シミュレーターコード全体（classical_gksl_simulator.py, qudit_gksl_simulator.py, qubit_gksl_simulator.py, stinespring_utils.py, gksl_math_utils.py）
4. 詳細理論・詳細設計・詳細仕様

## 1. Iteration 34 検証スクリプト結果

### 実行結果サマリー

| チェック項目 | 結果 | 詳細 |
|---|---|---|
| トレース距離単調減少 | ✓ PASS | 全 n_steps で単調に減少 |
| 収束次数 ≈ 1.0 | ✓ PASS | 平均 1.2010（許容範囲 ±0.3） |
| 忠実度単調増加 | ✓ PASS | 全 n_steps で単調に増加 |
| 密度行列品質 | ✓ PASS | トレース誤差 < 1e-10, エルミート性 < 1e-10 |
| n_steps=100 で T < 1e-3 | ✓ PASS | T = 7.319e-04 |
| Cell 27 ノートブック内容検証 | ✓ PASS | テーブル値が実測値と一致（許容範囲内） |

### 収束解析データ

| n_steps | dt | T(ST, exact) | Fidelity | Rate(T) |
|--------:|-------:|--------------:|----------:|--------:|
| 10 | 10.0000 | 1.294305e-02 | 0.99966585 | --- |
| 20 | 5.0000 | 4.081287e-03 | 0.99996979 | 1.6651 |
| 50 | 2.0000 | 1.492513e-03 | 0.99999631 | 1.0979 |
| 100 | 1.0000 | 7.319238e-04 | 0.99999916 | 1.0280 |
| 200 | 0.5000 | 3.626913e-04 | 0.99999982 | 1.0130 |

## 2. ノートブック全セル分析

### 各セル実行結果の確認

| セル | 内容 | 結果 | 備考 |
|------|------|------|------|
| Cell 2 | パラメータ設定 | ✓ | Hilbert空間次元 81, バリデーション通過 |
| Cell 4 | Classical GKSL (ボソン無し) | ✓ | Trace conservation: 5.33e-15 |
| Cell 6 | Qudit GKSL (ボソン無し) | ✓ | Trace conservation: 2.02e-14 |
| Cell 8 | Qubit GKSL (ボソン無し) | ✓ | Trace conservation: 2.04e-14 |
| Cell 10 | Classical GKSL (ボソン有り) | ✓ | Trace conservation: 5.77e-15 |
| Cell 12 | Qudit GKSL (ボソン有り) | ✓ | Trace conservation: 1.33e-14 |
| Cell 14 | Qubit GKSL (ボソン有り) | ✓ | Trace conservation: 1.33e-14 |
| Cell 16 | 量子回路可視化 | ✓ | 全4シナリオの回路図生成成功 |
| Cell 18 | 包括的比較プロット | ✓ | 正常出力 |
| Cell 20 | ユニタリ vs GKSL 比較 | ✓ | ユニタリ: エントロピー ≈ 0（2.22e-16） |
| Cell 22 | 検証（保存則チェック） | ✓ | 全保存則が機械精度レベルで保存 |
| Cell 24 | 忠実度評価 | ✓ | Classical vs Qudit: F = 0.999999 |
| Cell 25 | 収束特性理論説明 | **△** | **「30回」が古い（後述）** |
| Cell 26 | 収束解析コード | ✓ | 結果は実測値と完全一致 |
| Cell 27 | 精度ガイダンステーブル | **△** | **n_steps=200 の近似値が若干粗い（後述）** |
| Cell 28-38 | ショットベース/ノイズ | ✓ | 物理的に妥当な結果 |
| Cell 39 | まとめ | **△** | **「30回」が古い（後述）** |

### エラー・警告

ノートブック全セルの出力にエラーおよび警告は**一切検出されなかった**。

## 3. コア物理実装の数学的検証

### 3.1 Stinespring Dilation の正当性

**結論: 数学的に正しい。**

生成子 $G = \begin{pmatrix} 0 & L^\dagger \\ L & 0 \end{pmatrix}$ からユニタリ $U = e^{-i\theta G}$（$\theta = \sqrt{\Delta t}$）を構成。
Kraus 演算子は:
- $K_0 = \cos(\theta\sqrt{L^\dagger L})$
- $K_1 = -iL \cdot \sin(\theta\sqrt{L^\dagger L}) / \sqrt{L^\dagger L}$

1次展開で:
$$\mathcal{E}(\rho) = K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger \approx \rho + \Delta t (L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L, \rho\}) + O(\Delta t^2)$$

- CPTP性: $K_0^\dagger K_0 + K_1^\dagger K_1 = I$（ユニタリ性より厳密成立）✓
- 近似次数: $O(\Delta t^2)$ /ステップ → $O(\Delta t)$ 全体 ✓

### 3.2 Trotter 分解の正当性

**結論: 正しい実装。**

対称 Strang 分割 + 回文順序積:
$$e^{\mathcal{L}\Delta t} \approx e^{\mathcal{L}_H \Delta t/2} \cdot \prod_{\alpha=1}^n \mathcal{E}_\alpha(\Delta t/2) \cdot \prod_{\alpha=n}^1 \mathcal{E}_\alpha(\Delta t/2) \cdot e^{\mathcal{L}_H \Delta t/2}$$

- Strang 分割: $O(\Delta t^3)$/ステップ ✓
- 回文順序積: Lie-Trotter 交換子誤差の2次消去 ✓
- Stinespring 近似: $O(\Delta t^2)$/ステップ（支配的）→ 全体 $O(\Delta t)$ ✓

### 3.3 超演算子（Liouvillian）の正当性

**結論: 正しい実装。**

列優先ベクトル化 $\text{vec}(A X B) = (B^T \otimes A)\text{vec}(X)$ に基づき:
$$\mathcal{L} = -\frac{i}{\hbar}(I \otimes H - H^T \otimes I) + \sum_\alpha \left[\bar{L}_\alpha \otimes L_\alpha - \frac{1}{2}(I \otimes L_\alpha^\dagger L_\alpha + (L_\alpha^\dagger L_\alpha)^T \otimes I)\right]$$

`build_gksl_superoperator` の実装がこの公式と完全に一致していることを確認 ✓

### 3.4 Lindblad 演算子の物理モデル

**結論: TTA-UC の物理を正しくモデル化。**

| チャネル | 数 | 遷移 | レート |
|---|---|---|---|
| TTA | 6 (3ペア×2方向) | $\|1,1\rangle \to \|2,0\rangle, \|0,2\rangle$ | $\gamma_{TTA}/2$ |
| 蛍光 | 4 | $S_1 \to S_0$ | $\Gamma_{fl}$ |
| りん光 | 4 | $T_1 \to S_0$ | $\Gamma_{ph}$ |
| 内部変換 | 4 | $S_1 \to S_0$ | $k_{IC}$ |
| ISC (S→T) | 4 | $S_1 \to T_1$ | $k_{ISC\_ST}$ |
| ISC (T→S) | 4 | $T_1 \to S_0$ | $k_{ISC\_TS}$ |
| **合計** | **26** | | |

全演算子に $\sqrt{\gamma}$ が含まれており、レートの二重カウントなし ✓

### 3.5 数値検証（2分子系での厳密比較）

2分子系（dim=9）で Stinespring+Trotter と厳密解（行列指数関数）を比較:

| n_steps | dt | T(Trotter, exact) | Rate(T) |
|--------:|-----:|--------------:|--------:|
| 10 | 1.0000 | 6.334e-04 | --- |
| 20 | 0.5000 | 3.163e-04 | 1.002 |
| 50 | 0.2000 | 1.264e-04 | 1.001 |
| 100 | 0.1000 | 6.318e-05 | 1.000 |
| 200 | 0.0500 | 3.159e-05 | 1.000 |
| 500 | 0.0200 | 1.263e-05 | 1.000 |

収束次数が **正確に 1.0** であることを確認。これは Stinespring dilation の数学的に予測される挙動と完全に一致する。

## 4. 発見された問題点

### 問題 1: Cell 25, 39 のイテレーション回数が古い（重要度: 低）

**内容**: Cell 25 と Cell 39 に「30回のイテレーション検証」と記載されているが、現在はイテレーション 34 まで完了している。

**対象箇所**:
- Cell 25: 「30回のイテレーション検証により確認された収束特性を定量的に示します」
- Cell 39: 「30回の検証イテレーションで確認済み」

**修正内容**: 「30回」→「34回」に更新。

### 問題 2: Cell 27 テーブルの n_steps=200 近似値の精度（重要度: 低）

**内容**: Cell 27 テーブルの n_steps=200 のトレース距離値が `~4e-04` と記載されているが、実測値は `3.627e-04` である。比率は 0.907 で、他の行と比較して最も精度が低い。

| n_steps | テーブル値 | 実測値 | 比率 | 精度 |
|--------:|:---------:|:-------:|:-----:|:----:|
| 10 | ~1.3e-02 | 1.294e-02 | 0.996 | ◎ |
| 20 | ~4e-03 | 4.081e-03 | 1.020 | ○ |
| 50 | ~1.5e-03 | 1.493e-03 | 0.995 | ◎ |
| 100 | ~7e-04 | 7.319e-04 | 1.046 | ○ |
| 200 | ~4e-04 | 3.627e-04 | 0.907 | △ |

**修正内容**: `~4e-04` → `~3.6e-04` に更新。

## 5. 結論

### コードの正当性

**全シミュレーターコードは数学的に正しく実装されている。** 具体的には:

1. **Stinespring dilation**: CPTP写像として厳密に正しい
2. **Trotter 分解**: 対称分割と回文順序積の組み合わせが正しい
3. **超演算子構成**: 列優先ベクトル化と一致
4. **Lindblad 演算子**: TTA-UC 物理モデルとして正確
5. **初期状態/可観測量**: 正しい tensor product 構造
6. **保存則**: トレース保存（< 1e-14）、粒子数保存（< 1e-13）

### 収束特性

$O(\Delta t)$ 収束は Stinespring dilation の数学的性質であり、バグではない。数値的に Rate(T) → 1.0 が厳密に確認された。

### 修正が必要な箇所

2つの軽微な文書整合性の問題のみ（上記 問題 1, 2）。コードロジックの修正は不要。
