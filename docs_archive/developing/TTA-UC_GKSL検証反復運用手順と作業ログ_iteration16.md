# TTA-UC GKSL 検証反復運用手順と作業ログ Iteration 16

## 概要

- 反復番号: 16
- 日時: 2026-02-24
- 前回: iteration 15 (developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration15.md)
- 目的: iteration 15 検証結果の詳細分析、ノイズ-Trotterトレードオフ検証、
  ノイズバジェット分解、解析的ノイズ境界との比較

## 背景

Iteration 15 の検証結果を詳細に分析した結果:

1. **コードにバグは検出されなかった** — 8つのシミュレータファイルすべてにおいて
   数学的正確性を確認
2. **Trotter収束は優れている** — dt=1.0 でも F > 0.99997
3. **厳密GKSL解との一致は非常に高い** — F > 0.99996 at t=10
4. **DM-Shot整合性は高い** — F_norm > 0.99
5. **ノイズモデルは理論予測と定量的に一致** — N_S1 成長、禁止状態leakage ともに
   解析的推定と整合

しかし、以下の検証方法論上の改善点が特定された:

- ノイズ-Trotterトレードオフ（n_steps最適化）の未検証
- ノイズバジェット分解（ゲート種類別寄与）の不在
- Qubit忠実度指標の報告改善の必要性

## 作業手順

### ステップ1: 分析ドキュメント作成

`developing/検証結果分析_20260224_iteration16.md` を作成:
- iteration 15 結果の数学的検証
- コード全体の正確性確認
- ノイズ挙動の定量的説明
- 4つの問題点の特定（バグなし、方法論改善3件）

### ステップ2: 検証スクリプト作成

`tutorials/run_tta_uc_gksl_verification_iteration16.py` を作成:

- **Test A: ノイズ-Trotterトレードオフ** — n_steps=5,10,20,50,100 で
  ノイジー Qudit を実行し、最適 n_steps を探索
- **Test B: ノイズバジェット分解** — Hamiltonian ノイズのみ、Lindblad ノイズのみ、
  全ノイズの3条件で1ステップを実行し、各寄与を定量化
- **Test C: 解析的ノイズ境界** — 理論的忠実度減衰と実測値を比較
- **Test D: Qubit正規化忠実度追跡** — DMノイジー多ステップで raw F と
  F_norm の両方を各ステップで記録

### ステップ3: 作業ログ作成

本ファイル `developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration16.md`

### ステップ4: ユーザーによるローカル検証（次回）

ユーザーがローカルで以下を実行:
```bash
cd tutorials
python run_tta_uc_gksl_verification_iteration16.py
```

結果が `developing/verification_results/` に出力される。

## 修正ファイル一覧

| ファイル | 変更種別 | 内容 |
|---------|---------|------|
| tutorials/run_tta_uc_gksl_verification_iteration16.py | 新規 | iteration16 検証スクリプト |
| developing/検証結果分析_20260224_iteration16.md | 新規 | iteration15 結果分析 |
| developing/TTA-UC_GKSL検証反復運用手順と作業ログ_iteration16.md | 新規 | 本作業ログ |

## コード変更について

シミュレータ本体のコード（qudit_gksl_noisy_simulator.py 等）には変更なし。
iteration 15 のコード修正は正確であり、バグは検出されなかった。
今回は検証方法の改善（新しい検証スクリプト）のみ。

## 次回の作業

1. ユーザーがローカルで iteration16 検証スクリプトを実行
2. 結果をリポジトリにpush
3. 検証結果を確認し:
   - ノイズ-Trotterトレードオフの最適 n_steps を決定
   - 解析的ノイズ境界との一致を確認
   - 必要に応じてドキュメントに推奨パラメータを追記
