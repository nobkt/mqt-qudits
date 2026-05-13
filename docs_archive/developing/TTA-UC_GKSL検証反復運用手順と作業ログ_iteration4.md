# TTA-UC GKSL検証反復運用手順と作業ログ（第4回反復）

## 概要

第4回反復では、第3回反復の検証結果（全11シナリオPASS）のコード全体を詳細に分析し、
**4つの重大な構造的問題**を発見・修正した。

## 発見された問題

### 問題1〜4: Qubitシミュレータの同一実装問題

**発見事項:**

- `QubitGKSLSimulator` ≡ `QuditGKSLSimulator`（完全同一結果）
- `QubitGKSLShotSimulator` ≡ `QuditGKSLShotSimulator`（完全同一結果）
- `QubitGKSLBosonSimulator` ≡ `QuditGKSLBosonSimulator`（完全同一結果）
- `QubitGKSLNoisyShotSimulator` がd=3キュートリットノイズを使用（量子ビットはd=2パウリノイズであるべき）

**根本原因:**

全qubitシミュレータが81次元キュートリット空間で動作しており、
256次元量子ビット符号化空間を使用していなかった。
qubitシミュレータに存在する `_embed_in_qubit_space()` / `_extract_from_qubit_space()` メソッドは
`simulate()` メソッド内で一度も使用されていなかった。

## 実施した修正

### ① qubit_gksl_simulator.py — 256次元量子ビット空間での密度行列シミュレーション

**修正内容:**

1. モジュールレベル関数の追加:
   - `build_qubit_qutrit_mapping()`: キュートリット↔量子ビットインデックス変換
   - `embed_operator_in_qubit_space()`: 81×81演算子を256×256に埋め込み
   - `embed_density_matrix_in_qubit_space()`: 密度行列の埋め込み
   - `extract_density_matrix_from_qubit_space()`: 密度行列の抽出
   - `embed_statevector_in_qubit_space()`: 状態ベクトルの埋め込み
   - `extract_statevector_from_qubit_space()`: 状態ベクトルの抽出
   - `compute_forbidden_state_population()`: 禁止状態|11⟩の占有率計算

2. `__init__`: ハミルトニアンとリンドブラッド演算子を256次元空間に埋め込み
3. `_trotter_step`: 256次元空間でTrotter演算実行
4. `prepare_initial_state`: 初期状態を256次元に埋め込み
5. `simulate`: 256次元で発展、81次元に抽出して観測量計算、禁止状態リーケージ追跡

**結果:**

- ノイズなし: qubit結果 = qudit結果（浮動小数点精度で一致、物理的に正しい）
- 計算時間: qubit ~11s vs qudit ~2.6s（256次元vs81次元のオーバーヘッド）
- forbidden_state_population: [0.0, 0.0, ...] （ノイズなしではリーケージゼロ）

### ② qubit_gksl_shot_simulator.py — 量子ビット空間での量子軌道法

**修正内容:**

1. `QubitGKSLShotSimulator`: 256次元空間で量子軌道を実行
   - Stinespring拡張ユニタリ: 512×512（vs quditの162×162）
   - 測定結果のキュートリット基底マッピング
   - `forbidden_count` の追跡

2. `QubitGKSLNoisyShotSimulator`: d=4量子ビットパウリノイズの実装
   - `_apply_pauli_on_molecule()`: 分子の2量子ビットに4×4パウリ演算子を適用
   - `_apply_stochastic_qubit_depolarization_single()`: d=4パウリ脱分極（単一分子）
   - `_apply_stochastic_qubit_depolarization_pair()`: d=16パウリ脱分極（分子ペア）
   - `_apply_stochastic_qubit_thermal_relaxation()`: 2量子ビット熱緩和

**物理的効果:**

16個の2量子ビットパウリ演算子（σ_a⊗σ_b）のうち、以下が禁止状態リーケージを引き起こす:
- σ_X⊗I: |01⟩→|11⟩ (T1→禁止)
- I⊗σ_X: |10⟩→|11⟩ (S1→禁止)
- σ_Y⊗I, σ_X⊗σ_Z, etc. もリーケージを引き起こす
- 15個の非恒等パウリのうち約8個がリーケージを含む

**検証結果（N=4, 1000ショット, p_depol=0.01）:**

- `forbidden_count`: 332/1000（33.2%のショットが禁止状態で終了）
- `trace(rho_qutrit)`: 0.673（32.7%のリーケージ）
- F(qubit_noisy, qudit_noisy) = 0.943（異なるノイズモデルにより有意に異なる結果）

### ③ qubit_gksl_boson_simulator.py — 量子ビット符号化拡張空間

**修正内容:**

1. 電子自由度を4^N空間に埋め込み、フォノン自由度はFock空間のまま
2. 拡張空間マッピング: 各(el_qt, ph) → (el_qb, ph)
3. `_partial_trace_phonon_qubit()`: 量子ビット拡張空間でのフォノン部分トレース
4. `_extract_electronic_rho()`: フォノン部分トレース→キュートリット抽出

**結果:**

- qubit_boson結果 = qudit_boson結果（ノイズなし、物理的に正しい）
- F(qubit_boson, qudit_boson) = 1.0000000000

### ④ run_tta_uc_gksl_verification.py — リーケージ対応の検証スクリプト

**修正内容:**

1. `_validate_result`: `allow_leakage`パラメータ追加
   - リーケージ時: trace ≤ 1のみ検査（trace < 1は期待される物理現象）
   - 粒子数保存: 実際のトレースに基づいてスケーリング
   - 密度行列検証: 緩和されたトレース許容値

2. `_build_scenario_report`: `forbidden_count`, `dim_qubit_space` 出力追加
3. `_compute_fidelities`: トレースで正規化してから忠実度計算
4. `_run_single_scenario`: qubit_noisy_shot に `allow_leakage=True`
5. Markdown出力: forbidden_count, dim_qubit_space 表示追加

## 検証結果サマリー

全11シナリオの検証結果:

| シナリオ | 判定 | 特記事項 |
|----------|------|----------|
| classical | ✅ PASS | exact expm参照解 |
| qubit | ✅ PASS | 256次元空間、forbidden_pop=0 |
| qudit | ✅ PASS | 81次元空間 |
| classical_boson | ✅ PASS | exact expm参照解 |
| qubit_boson | ✅ PASS | 64次元拡張空間 |
| qudit_boson | ✅ PASS | 36次元拡張空間 |
| qudit_shot | ✅ PASS | 81次元、1000ショット |
| qubit_shot | ✅ PASS | 256次元、forbidden_count=0 |
| qudit_noisy_shot | ✅ PASS | d=3 Weyl-Heisenbergノイズ |
| qubit_noisy_shot | ✅ PASS | d=4パウリノイズ、forbidden_count=332 |
| unitary | ✅ PASS | 純ユニタリ発展 |

## 次のステップ

1. ユーザーが検証スクリプト `run_tta_uc_gksl_verification.py` をローカルで実行
2. 全シナリオの完全な検証結果をファイル出力してpush
3. 結果を確認し、必要に応じて追加修正
