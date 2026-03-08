# 量子ダイナミクス完全比較ノートブック 検証結果分析レポート（Iteration 2）

## 作成日: 2026-03-08
## 対象ファイル:
- `tutorials/quantum_dynamics_complete_comparison.ipynb`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`
- `tutorials/run_quantum_dynamics_verification.py`
## 前回参照: `developing/検証結果分析_quantum_dynamics_comparison_iteration1.md`

---

## 1. Iteration 1 からの修正状況確認

### 1.1 修正済み問題

Iteration 1 で報告された3つの問題の修正状況を確認した。

| 問題 | ステータス | 確認方法 |
|------|-----------|---------|
| 問題1: Quditシミュレータの放射減衰 | ✅ 修正済み | `simulate()`と`simulate_shot_based()`のコードに放射減衰の呼び出しがないことを確認 |
| 問題1b: 累積時間バグ | ✅ 修正済み | 上記と同時に解消（放射減衰自体を除去） |
| 問題2: 量子リソース表示 | ✅ 修正済み | Cell 35の出力で「8 qubits」「4 qutrits」と正しく表示されることを確認 |
| 問題3: 回路可視化の不整合 | ✅ 修正済み | Cell 13, 20, 26に`_draw_qiskit_subcircuits()`と`_draw_qudit_circuit_mpl()`を導入済み |

### 1.2 独立検証スクリプトの結果

`run_quantum_dynamics_verification.py` の最新結果（`quantum_dynamics_verification_20260306T061845Z.json`）:

| テスト | ステータス | 主要数値 |
|--------|-----------|---------|
| H0位相検証 | PASS | max_error = 2.22e-16 |
| Trotter比較 | PASS | max_error = 1.08e-14 |
| 脱分極チャネル | PASS | 全6サブテスト合格 |
| ノイズ蓄積分析 | PASS | 88/12 gates, 7.33×削減 |

### 1.3 ノートブック実行結果の分析

Cell 32の精度比較出力を確認:

```
精度比較: Qubit実装 vs 古典シミュレーション
  N_T1の最大誤差: 0.015810（ショットノイズレベル）

精度比較: Qudit実装 vs 古典シミュレーション
  N_T1の最大誤差: 0.015744（ショットノイズレベル）
```

両者とも ~0.016 のショットノイズレベル（10000ショットで 1/√10000 ≈ 0.01 のオーダー）であり、Iteration 1で報告された0.169のQudit誤差は完全に解消されている。

---

## 2. 新たに発見された問題

### 2.1 問題A（軽微）: `_draw_qiskit_subcircuits`の型注釈の不整合

#### 問題の詳細

`quantum_dynamics_complete_comparison.ipynb` の Cell 13:
```python
def _draw_qiskit_subcircuits(
    title: str,
    named_circuits: list,     # ← 型注釈が不完全
    prefix: str,
) -> None:
```

`quantum_dynamics_gksl_comparison.ipynb` の Cell 16:
```python
def _draw_qiskit_subcircuits(
    title: str,
    named_circuits: list[tuple[str, QuantumCircuit]],  # ← 正確な型注釈
    prefix: str,
) -> None:
```

#### 影響
動作に影響はないが、型情報の一貫性の観点から修正すべき。

#### 修正方針
Cell 13の型注釈をGKSLノートブックと同一にする。

---

## 3. 2つのノートブック間の整合性検証

### 3.1 物理モデルの違い

2つのノートブックは**意図的に異なる物理モデル**を使用している：

| 項目 | `complete_comparison` | `gksl_comparison` |
|------|----------------------|-------------------|
| **TTA模型** | ハミルトニアン結合 (H_TTA = J[...+h.c.]) | Lindblad散逸子 (L_TTA = √γ_TTA[...]) |
| **全ハミルトニアン** | H = H₀ + H_transfer + H_TTA | H = H₀ + H_transfer |
| **TTA過程の性質** | 可逆（ユニタリ） | 不可逆（散逸） |
| **S1生成** | コヒーレント振動: \|T1,T1⟩ ↔ \|S0,S1⟩ | 不可逆遷移: \|T1,T1⟩ → \|S0,S1⟩ |
| **エネルギー保存** | 厳密に保存 | 環境への散逸あり |
| **目的** | 量子コンピューティングプラットフォーム比較 | 開放量子系の散逸過程シミュレーション |

### 3.2 パラメータの対応関係

| パラメータ | `complete_comparison` | `gksl_comparison` | 一致 |
|-----------|----------------------|-------------------|------|
| E_T | 1.5 eV | 1.5 eV | ✅ |
| E_S | 3.0 eV | 3.0 eV | ✅ |
| V | 0.1 eV | 0.1 eV | ✅ |
| J (TTA結合定数) | 0.05 eV | N/A (散逸率として表現) | 異なる物理量 |
| γ_TTA | N/A | 0.05 eV/ℏ | 異なる物理量 |
| N_molecules | 4 | 4 | ✅ |
| t_max | 100 fs | 100 fs | ✅ |
| N_steps | 20 | 100 | 異なる（精度要件の違い） |
| ℏ | 0.6582 eV·fs | 1.0（自然単位系） | **注意**: 単位系が異なる |

### 3.3 単位系の整合性

- `complete_comparison`: ℏ = 0.6582119569 eV·fs（SI単位系を採用）
- `gksl_comparison`: ℏ = 1.0（自然単位系を採用、GKSLPhysicalParameters.hbar プロパティ）

この違いは設計上のものであり、各ノートブック内で一貫している。ただし、2つのノートブック間で数値を直接比較する場合は単位系の変換が必要。

### 3.4 GKSL ユニタリ極限との比較

GKSLノートブック（Cell 20）のユニタリ極限（全散逸率=0）では：
- H = H₀ + H_transfer のみ（H_TTAなし）
- N_T1 = 2.0, N_S1 = 0.0 のまま変化なし（三重項が分子間を移動するのみ）
- エントロピー: max = 2.22e-16（マシンイプシロン、純粋状態を維持）

complete_comparisonノートブック（Cell 5）の結果：
- H = H₀ + H_transfer + H_TTA
- N_T1 ≈ 1.19, N_S1 ≈ 0.40（100 fs後）
- H_TTAによるコヒーレントな T1+T1 ↔ S0+S1 振動

**この差異はバグではなく、異なる物理モデルの正しい結果である。**

### 3.5 結論

2つのノートブックは補完的な役割を果たしており、整合性に問題はない：
- `complete_comparison`: 量子コンピューティングプラットフォーム（qubit vs qudit）の公平な比較
- `gksl_comparison`: 開放量子系における散逸過程の忠実なシミュレーション

---

## 4. TTAモデル比較の検討

### 4.1 比較の動機

ユーザーの要求: 「quantum_dynamics_gksl_comparison.ipynbで散逸項をTTAだけにして、quantum_dynamics_complete_comparison.ipynbと比較してTTAがユニタリ相互作用から散逸項になった場合のポピュレーション変化の差異を検証」

### 4.2 比較対象

| ケース | ハミルトニアン | 散逸項 | 物理的意味 |
|--------|-------------|--------|-----------|
| A: ユニタリTTA | H = H₀ + H_transfer + H_TTA | なし | TTA過程が可逆的なコヒーレント相互作用 |
| B: 散逸TTA only | H = H₀ + H_transfer | L_TTA のみ | TTA過程のみが不可逆な散逸チャネル |
| C: 全散逸GKSL | H = H₀ + H_transfer | 全L_α | TTA + 蛍光 + 燐光 + IC + ISC の全散逸 |

### 4.3 実装方法

**ケースA**（ユニタリTTA）は`complete_comparison`ノートブックの古典シミュレータと同一の物理。GKSLフレームワーク内で再現するには、`gksl_math_utils.py`の`build_onsite_hamiltonian`と`build_transfer_hamiltonian`に加えて、新たにH_TTAハミルトニアンを構築し、H = H₀ + H_transfer + H_TTA で散逸なしの時間発展を計算する。

ただし**重要な注意**: 2つのノートブックでは単位系が異なる。
- `complete_comparison`: ℏ = 0.6582119569 eV·fs で、時間発展は exp(-iHt/ℏ)
- `gksl_comparison`: ℏ = 1.0（自然単位系）で、時間発展は exp(-iHt)

同じ100 fs の時間発展でも、H の係数に ℏ の違いが反映される。`complete_comparison`の PhysicalParameters では J = 0.05 eV だが、ℏ=1 の自然単位系では J/ℏ に相当する値を使う必要がある。実際には `complete_comparison` の exp(-iHt/ℏ) = exp(-iH·100/0.6582) に対し、`gksl_comparison` では exp(-iHt) = exp(-iH·100) となる。

これを正しく比較するために、GKSLノートブック内にユニタリTTAシミュレーションを追加する場合は：
1. GKSLの単位系（ℏ=1）で H_TTA を構築
2. GKSLと同じパラメータ（E_T, E_S, V）を使用
3. J_TTA パラメータを GKSLの γ_TTA と比較可能な形で設定

ただし、J（ハミルトニアン結合定数）と γ_TTA（散逸率）は物理的に異なる量であり、「同じ値」にすることに物理的根拠はない。J_TTA = γ_TTA と設定するのは便宜上の選択であり、定性的な振る舞いの違いを示すためのものである。

### 4.4 実装の実現可能性

**実装可能**。以下の手順で実装する：

1. GKSLノートブックに新セルを追加
2. `gksl_math_utils.py`の既存関数を使って H₀, H_transfer を構築
3. H_TTA ハミルトニアンを構築する関数を追加
4. `scipy.linalg.expm`で exp(-iHt) を計算
5. ケースA（ユニタリTTA）、B（散逸TTA only）、C（全散逸）を比較プロット

### 4.5 注意事項

- この比較は**定性的なもの**であり、J_TTA = γ_TTA は物理的な等価性を意味しない
- ユニタリTTAではコヒーレント振動が見られ、散逸TTAでは指数的減衰が見られるはず
- この違いはTTA過程のコヒーレンス/デコヒーレンスの物理的差異を反映する

---

## 5. 実施する修正一覧

| ファイル | 修正内容 | 理由 |
|----------|----------|------|
| `tutorials/quantum_dynamics_complete_comparison.ipynb` Cell 13 | `_draw_qiskit_subcircuits`の型注釈を修正 | GKSL notebookとの一貫性 |
| `tutorials/quantum_dynamics_gksl_comparison.ipynb` 新セル追加 | ユニタリTTA vs 散逸TTAの比較セルを追加 | ユーザー要求 |

---

## 6. 作業ログ

| 手順 | ステータス | 説明 |
|------|-----------|------|
| ① 前回修正状況の確認 | ✅ 完了 | 3つの問題すべて修正済み |
| ② 検証結果分析 | ✅ 完了 | 全4テストPASS、精度はショットノイズレベル |
| ③ ノートブック間整合性分析 | ✅ 完了 | 異なる物理モデル（設計上正しい） |
| ④ 型注釈修正 | 🔄 実施中 | Cell 13の型注釈をGKSL notebookに合わせる |
| ⑤ TTA比較セル追加 | 🔄 実施中 | GKSLノートブックに新セル追加 |
| ⑥ 修正後の検証 | ⏳ 未実施 | 修正後に検証実施予定 |
