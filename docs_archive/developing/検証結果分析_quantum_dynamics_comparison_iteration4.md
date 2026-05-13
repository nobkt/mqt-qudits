# 量子ダイナミクス完全比較ノートブック 検証結果分析レポート（Iteration 4）

## 作成日: 2026-03-09
## 対象ファイル:
- `tutorials/quantum_dynamics_complete_comparison.ipynb`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`
## 前回参照: `developing/検証結果分析_quantum_dynamics_comparison_iteration3.md`

---

## 1. 問題の特定

### 問題1: 基本ゲート分解後の回路可視化が膨大（重大度: 高）

**原因**: `quantum_dynamics_complete_comparison.ipynb` において、基本ゲート分解後の量子回路を可視化すると、
数百〜数千ゲートの膨大な回路図が出力される。これはノートブックの可読性を著しく損なう。

**影響を受けるセル**:
- Cell 9: `qubit_results['step_circuit']`（基本ゲート分解版）の回路可視化
  - `circuit_drawer(step_circuit, output='mpl', fold=100)` により膨大な回路図が出力される
  - 88+ CNOTゲートを含む回路は1画面に収まらない
- Cell 13: `("基本ゲート分解版", step_circuit_basic)` の回路可視化
  - 基本ゲート分解版とUnitaryGate版の両方を描画していた

**修正**:
- Cell 9: `circuit_drawer` による回路可視化を完全に削除。ゲート統計のテキスト出力は維持
- Cell 13: `基本ゲート分解版` の可視化を削除。`UnitaryGate版` のみを表示
- Cell 26: 既に`CustomTwoゲート版`のみ表示（前回Iteration 3で修正済み）

**根拠**: 
- ユニタリゲートやカスタムゲートの状態での可視化は、高レベルの回路構造を理解するのに有用
- 基本ゲート分解後の可視化は、数百ゲートの膨大な図となり実用的でない
- ゲート数の比較分析（Cell 12, Cell 24）はテキスト出力として保持し、定量的情報は失わない

### 問題2: 9b. ユニタリTTA vs 散逸TTA 比較にショットベースシミュレーション結果が欠如（重大度: 中）

**原因**: `quantum_dynamics_gksl_comparison.ipynb` のセクション9b（Cell 22）では、
密度行列レベルの厳密シミュレーション（ケースA, B, C）のみで比較を行っていた。
ショットベースシミュレータ（ノイズ無し・ノイズ有り）の結果が含まれていなかった。

**影響**: 
- 量子回路実装レベルでの比較ができない
- ショットノイズやハードウェアノイズの影響が不明
- `quantum_dynamics_complete_comparison.ipynb`（ユニタリTTA）と`quantum_dynamics_gksl_comparison.ipynb`（散逸TTA）の
  実装レベルでの比較が不完全

**修正**:
- Cell 40の後に新セクション「9c. ユニタリTTA vs 散逸TTA の比較（ショットベースシミュレーション付き）」を追加
- 既に計算済みのショットベース結果（result5b, result5c, result3b, result3c）を使用
- 4パネル構成の比較プロット:
  - (1) Qudit Shot（ノイズ無し）vs ユニタリTTA vs GKSL DM
  - (2) Qudit Shot（ノイズ有り）vs ユニタリTTA vs GKSL DM
  - (3) Qubit Shot（ノイズ無し）vs ユニタリTTA vs GKSL DM
  - (4) Qubit Shot（ノイズ有り）vs ユニタリTTA vs GKSL DM
- ノイズ影響比較プロット:
  - Qudit: DM vs Shot (No Noise) vs Shot (Noisy)
  - Qubit: DM vs Shot (No Noise) vs Shot (Noisy)
- 定量的サマリーテーブル

**注意事項**:
- ヒューリスティックな処理やfallbackは一切使用していない
- すべてのシミュレーション結果は既に他のセルで正しく計算されたものを再利用
- セクション9cはセクション12b（ショットベースシミュレーション）の後に配置し、
  すべてのshot結果変数が利用可能な状態で実行

---

## 2. 修正内容の一覧

| ファイル | セル | 修正内容 | 理由 |
|---------|------|---------|------|
| complete_comparison | Cell 9 | `circuit_drawer`による回路可視化を削除 | 基本ゲート分解後の膨大な回路出力を抑制 |
| complete_comparison | Cell 13 | `基本ゲート分解版`の可視化を削除、`UnitaryGate版`のみ表示 | 同上 |
| gksl_comparison | Cell 41 (新規) | セクション9cのMarkdown解説を追加 | ショットベース比較の背景説明 |
| gksl_comparison | Cell 42 (新規) | ショットベースTTA比較のプロット・サマリーコードを追加 | ノイズ有無での比較を追加 |

---

## 3. 技術的詳細

### 3.1 回路可視化の方針

| 回路タイプ | 可視化 | 理由 |
|-----------|--------|------|
| UnitaryGate版（Qubit） | ✅ 表示 | 高レベル構造が把握しやすい |
| CustomTwoゲート版（Qudit） | ✅ 表示 | 同上 |
| 基本ゲート分解版（Qubit） | ❌ テキストのみ | 88+ CNOTで膨大 |
| 基本ゲート分解版（Qudit） | ❌ テキストのみ | 分解後ゲート数が多い |

### 3.2 ショットベース比較の構成

**新セクション9cのデータフロー**:
```
Cell 22 (9b): result_unitary_tta, result_tta_only → 密度行列レベル比較
Cell 4:  result1 (Classical GKSL DM)
Cell 6:  result5 (Qudit GKSL DM)
Cell 8:  result3 (Qubit GKSL DM)
Cell 32: result5b (Qudit Shot, no noise)
Cell 34: result5c (Qudit Shot, noisy)
Cell 36: result3b (Qubit Shot, no noise)
Cell 38: result3c (Qubit Shot, noisy)
    ↓
Cell 42 (9c): ショットベース比較プロット + サマリーテーブル
```

### 3.3 期待される比較結果

1. **ユニタリTTA vs 散逸GKSL**: コヒーレント振動 vs 不可逆減衰の定性的差異
2. **Shot (no noise) vs DM**: ショットノイズによるわずかな揺らぎのみ
3. **Shot (noisy) vs DM**: ハードウェアノイズによる追加的デコヒーレンス
4. **Qudit vs Qubit**: ゲート数の違いによるノイズ蓄積差

---

## 4. Iteration 3からの変更点

| 項目 | Iteration 3 | Iteration 4 |
|------|------------|------------|
| 回路可視化 | 分解版を一部削除 | 基本ゲート分解版の可視化を全面削除 |
| ショットベース比較 | なし | 新セクション9cを追加 |
| N_steps | 100（修正済み） | 変更なし |
| matplotlibバックエンド | Agg + savefig（修正済み） | 変更なし |

---

## 5. 検証手順

修正後のノートブックを検証するには：

### quantum_dynamics_complete_comparison.ipynb:
1. 全セルを順次実行
2. 以下を確認：
   - Cell 9: 回路図が表示されないこと（ゲート統計テキストのみ）
   - Cell 13: UnitaryGate版のみ1つの回路図が表示されること
   - Cell 26: CustomTwoゲート版のみ1つの回路図が表示されること
   - その他のプロット・テーブルが正常に表示されること

### quantum_dynamics_gksl_comparison.ipynb:
1. 全セルを順次実行
2. 以下を確認：
   - Cell 22 (9b): 従来の3ケース比較が正常に表示されること
   - Cell 42 (9c 新規): 4パネル比較プロットが表示されること
   - Cell 42 (9c 新規): ノイズ影響比較プロットが表示されること
   - Cell 42 (9c 新規): サマリーテーブルが正しく出力されること
   - Cell 43 (まとめ): 従来のまとめセクションが正常であること

---

## 6. 次回の検証で確認すべき事項

1. 両ノートブックの全セルが正常に実行されること
2. 新セクション9cのプロットが正しいデータを表示していること
3. ショットベース結果と密度行列結果の一致度が妥当であること
4. ノイズ有りの結果がノイズ無しと比較して物理的に妥当な乖離を示すこと
5. Qubit noisy結果の禁止状態リーケージが適切に扱われていること
