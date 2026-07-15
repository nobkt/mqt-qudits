# 量子ダイナミクス完全比較ノートブック 検証結果分析レポート（Iteration 5）

## 作成日: 2026-03-10
## 対象ファイル:
- `tutorials/quantum_dynamics_complete_comparison.ipynb`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`
## 前回参照: `developing/検証結果分析_quantum_dynamics_comparison_iteration4.md`

---

## 1. 問題の特定

### 問題1: GKSLノートブックのノイズモデルが不統一（重大度: 高）

**原因**: `quantum_dynamics_gksl_comparison.ipynb` のノイズ有りショットシミュレータ（Cell 34, Cell 38）が、
`quantum_dynamics_complete_comparison.ipynb` のノイズモデルと異なる設定を使用していた。

**具体的な不整合**:

| 項目 | complete_comparison | gksl_comparison（修正前） | gksl_comparison（修正後） |
|------|-------------------|-------------------------|-------------------------|
| 脱分極率 | 1.0% (depol_2q=0.01) | 1.0% (p_depol=0.01) | 1.0% (p_depol=0.01) |
| 位相緩和率 | 0% (なし) | **0.5% (p_dephasing=0.005)** | 0% (p_dephasing=0.0) |
| ノイズ適用対象 | 2-qubit/2-quditゲートのみ | **全ゲート（pair+Lindblad）** | ペアゲートのみ (depol_pair_only=True) |

**影響**:
- `p_dephasing=0.005` による過剰なノイズ蓄積
- 単一サイトLindbladチャネルへのノイズ適用による物理的に不正確な結果
- Qubit noisy（シナリオ3c）での過大な禁止状態リーケージ
- 両ノートブック間の比較が不公平

### 問題2: complete_comparison.ipynb のノイズモデル（確認結果: 問題なし）

**確認結果**: `quantum_dynamics_complete_comparison.ipynb` のノイズモデルは**既に正しく設定されていた**。

- Cell 15 (Qubit noisy): `QubitMolecularDynamicsSimulatorNoisy` は `depol_2q=0.01` のみ使用、
  2-qubitゲート（CX, CZ）にのみノイズ適用。1-qubitゲートは理想的。
- Cell 23 (Qudit noisy): `NoisyQuditMolecularDynamicsSimulator` は `depol_2q=0.01` のみ使用、
  2-quditゲート（cx, csum, ls, ms）にのみノイズ適用。1-quditゲートは理想的。

両セルとも脱分極ノイズのみ、2-qubit/2-quditゲートのみに適用されており、修正不要。

---

## 2. 修正内容の一覧

| ファイル | セル | 修正内容 | 理由 |
|---------|------|---------|------|
| gksl_comparison | Cell 30 (markdown) | ノイズ列を「脱分極+位相緩和」→「脱分極のみ（ペアゲートのみ）」に修正 | ドキュメント整合性 |
| gksl_comparison | Cell 33 (markdown) | ノイズ説明を脱分極のみ・ペアゲートのみに更新 | ドキュメント整合性 |
| gksl_comparison | Cell 34 (code) | `p_dephasing=0.005` → `p_dephasing=0.0`, `depol_pair_only=True` 追加 | **核心的修正** |
| gksl_comparison | Cell 37 (markdown) | ノイズ説明を脱分極のみ・ペアゲートのみに更新 | ドキュメント整合性 |
| gksl_comparison | Cell 38 (code) | `p_dephasing=0.005` → `p_dephasing=0.0`, `depol_pair_only=True` 追加 | **核心的修正** |
| gksl_comparison | Cell 39 (markdown) | ノイズモデル統一の説明を追加 | ドキュメント整合性 |
| gksl_comparison | Cell 41 (markdown) | 比較表を「脱分極のみ（ペアゲートのみ）」に更新 | ドキュメント整合性 |
| gksl_comparison | Cell 43 (markdown) | まとめのノイズモデル記述を更新 | ドキュメント整合性 |

---

## 3. 技術的詳細

### 3.1 `depol_pair_only=True` の意味

GKSLシミュレーションのTrotterステップは以下の構成：

```
1. ハミルトニアン半ステップ (NN結合ペアゲート) → 2-qudit/2-qubitゲート
2. Lindbladチャネル (放射減衰など) → 単一サイト操作
3. Lindbladチャネル (逆順) → 単一サイト操作
4. ハミルトニアン半ステップ (NN結合ペアゲート) → 2-qudit/2-qubitゲート
```

`depol_pair_only=True` の場合：
- ステップ1, 4（ペアゲート）の後にのみ脱分極ノイズを適用 ✅
- ステップ2, 3（Lindbladチャネル）の後にはノイズを適用しない ✅

`depol_pair_only=False` の場合（旧設定）：
- 全ステップの後にノイズを適用 ❌（過剰なノイズ蓄積）

### 3.2 `p_dephasing=0.0` の根拠

`quantum_dynamics_complete_comparison.ipynb` では：
- Qubit: `QubitMolecularDynamicsSimulatorNoisy` は脱分極ノイズのみ（Qiskit Aer `depolarizing_error`）
- Qudit: `NoisyQuditMolecularDynamicsSimulator` は脱分極ノイズのみ（密度行列チャネル）

両方とも位相緩和（dephasing）ノイズを使用していない。
GKSLノートブックでも同一条件にすることで、公平な比較が可能になる。

### 3.3 ノイズモデル統一のまとめ

| ノートブック | Qubit noisy | Qudit noisy |
|------------|------------|------------|
| complete_comparison | depol_2q=0.01, 2-qubitゲートのみ | depol_2q=0.01, 2-quditゲートのみ |
| gksl_comparison | p_depol=0.01, p_dephasing=0.0, depol_pair_only=True | p_depol=0.01, p_dephasing=0.0, depol_pair_only=True |

**共通条件**: 脱分極ノイズのみ (1.0%)、2-qubit/2-quditペアゲートのみに適用。

### 3.4 ヒューリスティック・フォールバックの不使用

本修正では以下の原則を厳守：
- ❌ 密度行列の強制正規化は使用しない
- ❌ 固有値クリッピングは使用しない
- ❌ ノイズパラメータの人為的調整（見た目を良くするための値の変更）は行わない
- ✅ 物理的に正当な理由に基づくパラメータ設定のみ
- ✅ 両ノートブック間のノイズモデル統一という明確な根拠

---

## 4. Iteration 4からの変更点

| 項目 | Iteration 4 | Iteration 5 |
|------|------------|------------|
| GKSL qudit noisy | p_dephasing=0.005 | **p_dephasing=0.0, depol_pair_only=True** |
| GKSL qubit noisy | p_dephasing=0.005 | **p_dephasing=0.0, depol_pair_only=True** |
| complete qubit noisy | depol_2q=0.01のみ（変更なし） | 変更なし ✅ |
| complete qudit noisy | depol_2q=0.01のみ（変更なし） | 変更なし ✅ |
| ノイズモデル統一 | 不統一 | **統一完了** |
| ドキュメント整合性 | 不整合あり | **整合性確認済み** |

---

## 5. 検証結果

### 静的検証（run_iteration5_verification.py）

```
検証結果: 38/38 チェック通過
✓ 全チェック通過
```

検証項目：
1. GKSLノートブックのコードセル（Cell 34, 38）のノイズパラメータ修正 ✅
2. GKSLノートブックのマークダウンセル（Cell 30, 33, 37, 39, 41, 43）の記述更新 ✅
3. complete_comparisonノートブックのノイズモデル確認（既に正しい） ✅
4. シミュレータクラスの`depol_pair_only`パラメータ存在確認 ✅
5. 両ノートブック間のノイズモデル一貫性確認 ✅

---

## 6. 次回の検証で確認すべき事項

1. 両ノートブックの全セルが正常に実行されること（ユーザーによるローカル実行）
2. `depol_pair_only=True` によりLindbladチャネル後のノイズが適用されないこと
3. Qubit noisy（シナリオ3c）の禁止状態リーケージが旧設定より改善されること
4. ノイズ有りの結果がノイズ無しと比較して物理的に妥当な乖離を示すこと
5. `complete_comparison` と `gksl_comparison` のノイズ影響の比較が定量的に妥当であること

---

## 7. 検証スクリプト

```bash
cd tutorials && python run_iteration5_verification.py
```

結果は `developing/verification_results/iteration5_noise_model_*.json` に保存される。
