# TTA-UC GKSL 検証反復運用手順と作業ログ: Iteration 36

## 実施日

2026-03-02

## 前回（Iteration 35）の概要

- 全9チェック PASS
- Cell 25/27/39 の文書修正を検証・確認
- コードバグなし、数学的正当性確認済み

## Iteration 36 の目的

Iteration 35 までの検証は主に最終時刻 t=t_max での収束特性の検証に限定されていた。Iteration 36 では、以下の**新規検証項目**を追加し、検証フレームワークの網羅性を大幅に向上させる。

### 新規チェック項目

| # | 項目 | 根拠 |
|:---:|:---|:---|
| 10 | 粒子数保存（全時刻ステップ） | N_S0 + N_T1 + N_S1 = N_molecules の検証 |
| 11 | 中間時刻での人口動態比較 | 最終時刻だけでなく中間時刻での一致確認 |
| 12 | 2分子系の収束次数 | Rate(T) → 1.0 をより小さなdt（200→500）で確認 |
| 13 | ユニタリ限界精度 | 散逸なしでハミルトニアン発展が機械精度であることの確認 |
| 14 | 固有値スペクトル比較 | 密度行列の固有値構造レベルでの一致 |

## 深層分析の結果

### 分析対象

1. コア実装コード全体（stinespring_utils.py, gksl_math_utils.py, classical_gksl_simulator.py, qudit_gksl_simulator.py, gksl_physical_parameters.py）
2. Stinespring dilation の数学的導出と検証
3. GKSL superoperator の構成と整合性
4. Lindblad 演算子の物理的正確性
5. ノートブック全40セルの実行結果

### 発見事項

#### コードバグ

**発見されなかった。**

以下の全項目について数学的正確性を検証した：

1. **Stinespring dilation**: G = [[0, L†], [L, 0]], U = exp(-iθG) で正しく GKSL 散逸子を O(dt²)/step で近似
2. **Kron 順序**: `np.kron(env0, rho)` と G 行列のブロック構造が整合
3. **部分トレース**: `rho_out = rho'[0:d,0:d] + rho'[d:,d:]` が正しい
4. **GKSL superoperator**: 列優先ベクトル化 vec(AXB) = (B^T ⊗ A) vec(X) と整合
5. **Lindblad 演算子**: 26 個（6 TTA + 4 蛍光 + 4 燐光 + 4 IC + 4 ISC S→T + 4 ISC T→S）、各演算子に √γ を含む
6. **Strang 分割**: exp(L_H dt/2) · [prod E_α(dt/2) · prod_rev E_α(dt/2)] · exp(L_H dt/2)
7. **回文順序**: Lie-Trotter 交換子誤差の2次消去

#### 数値検証結果（独自実施）

| 検証項目 | 結果 | 判定 |
|:---|:---|:---:|
| 粒子数保存（古典） | < 2.22e-14 | ✓ |
| 粒子数保存（Qudit） | < 1.24e-14 | ✓ |
| トレース保存（古典） | < 5.33e-15 | ✓ |
| トレース保存（Qudit） | < 2.89e-15 | ✓ |
| CPTP 性（全26チャネル） | Tr = 1.0, min_eig > 0 | ✓ |
| superoperator トレース保存 | ‖vec(I)^T L‖ < 3.5e-17 | ✓ |
| 単一チャネル誤差 O(dt²)/step | T/dt² → 1.04e-4 (定数) | ✓ |
| 全 Trotter ステップ誤差 O(dt²)/step | T/dt² → 1.04e-4 (定数) | ✓ |
| ユニタリ限界（散逸なし） | T ~ 1e-14 (機械精度) | ✓ |

#### 微小な文書記載の不正確性

iteration 34 作業ログの「Rate(T) = 1.0000 for 2-molecule systems」は漸近極限の表現としては正しいが、有限 dt での Rate はわずかに 1.0 より大きい（例: dt=0.5 で Rate ≈ 1.005）。O(dt) 収束の結論に影響しないため、コード修正は不要。

## 作業手順

### ステップ 1: 深層分析の実施

以下の数値検証を実施：
- 粒子数保存（全時刻ステップ）
- 中間時刻人口動態比較
- 2分子系収束検証
- ユニタリ限界検証
- 単一チャネルおよび全 Trotter ステップの誤差解析
- CPTP 性検証
- superoperator トレース保存性検証

### ステップ 2: 分析結果のドキュメント化

`developing/検証結果分析_iteration36_深層分析.md` を作成。

### ステップ 3: 検証スクリプトの作成

`tutorials/run_tta_uc_gksl_verification_iteration36.py` を作成。
- 既存の9チェック + 新規5チェック = 全14チェック

### ステップ 4: ユーザーによる検証スクリプトの実行

ユーザーがローカルで実行し、結果を `developing/verification_results/` に出力。

## 次回（Iteration 37）への引き継ぎ

- Iteration 36 の検証スクリプト実行結果の確認
- 全14チェックが PASS であれば、検証フレームワークは十分に網羅的であると判断
- 必要に応じて追加の検証項目を検討
