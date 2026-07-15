# TTA-UC GKSL 検証反復運用手順と作業ログ Iteration 14

## 概要

- 反復番号: 14
- 日時: 2026-02-24
- 前回: iteration 13 (developing/検証結果分析_20260224_iteration13.md)
- 目的: ペアのみ脱分極ノイズモードの実装とqubit/qudit精度比較

## 背景

Iteration 12の検証結果で、ノイズ有りシミュレーションの精度が極端に悪い
（F ≈ 0.09）ことが判明した。原因は1ステップあたり76回のノイズ適用。
単一サイトゲート後のノイズは物理的に過剰であるため、2-qubit/2-qudit以上の
相互作用のみに脱分極ノイズを適用するモードを追加する。

## 作業手順

### ステップ1: コード修正

以下の4つのノイジーシミュレータに `depol_pair_only` パラメータを追加:

1. `tutorials/qudit_gksl_noisy_simulator.py` - QuditGKSLNoisySimulator
2. `tutorials/qubit_gksl_noisy_simulator.py` - QubitGKSLNoisySimulator
3. `tutorials/qudit_gksl_shot_simulator.py` - QuditGKSLNoisyShotSimulator
4. `tutorials/qubit_gksl_shot_simulator.py` - QubitGKSLNoisyShotSimulator

修正内容:
- コンストラクタに `depol_pair_only: bool = False` パラメータ追加
- `_trotter_step` / `_trotter_step_trajectory` でsingle-siteノイズを条件分岐
- `simulate()` の返り値 `noise_params` に `depol_pair_only` を追加
- クラスdocstringを更新

### ステップ2: 検証スクリプト作成

`tutorials/run_tta_uc_gksl_verification_iteration14.py` を作成:
- Test 1: DM 1ステップ比較 (pair-only vs all-gates)
- Test 2: ショットベース 10ステップ忠実度追跡
- Test 3: エラーバジェット比較

### ステップ3: 分析ドキュメント作成

- `developing/検証結果分析_20260224_iteration14.md` - 問題分析と修正内容
- `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration14.md` - 本ファイル

### ステップ4: ユーザーによるローカル検証（次回）

ユーザーがローカルで以下を実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration14.py
```

結果が `developing/verification_results/` に出力される。

## 修正ファイル一覧

| ファイル | 変更種別 | 内容 |
|---------|---------|------|
| tutorials/qudit_gksl_noisy_simulator.py | 修正 | depol_pair_only パラメータ追加 |
| tutorials/qubit_gksl_noisy_simulator.py | 修正 | depol_pair_only パラメータ追加 |
| tutorials/qudit_gksl_shot_simulator.py | 修正 | depol_pair_only パラメータ追加 |
| tutorials/qubit_gksl_shot_simulator.py | 修正 | depol_pair_only パラメータ追加 |
| tutorials/run_tta_uc_gksl_verification_iteration14.py | 新規 | ペアのみ検証スクリプト |
| developing/検証結果分析_20260224_iteration14.md | 新規 | 問題分析ドキュメント |
| developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration14.md | 新規 | 本作業ログ |

## 次回の作業

1. ユーザーがローカルで iteration14 検証スクリプトを実行
2. 結果をリポジトリにpush
3. 検証結果を確認し、必要に応じて追加修正
