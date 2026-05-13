# 量子ダイナミクス完全比較ノートブック 検証結果分析レポート（Iteration 3）

## 作成日: 2026-03-09
## 対象ファイル:
- `tutorials/quantum_dynamics_complete_comparison.ipynb`
- `tutorials/quantum_dynamics_gksl_comparison.ipynb`
## 前回参照: `developing/検証結果分析_quantum_dynamics_comparison_iteration2.md`

---

## 1. 問題の特定

### 問題1: グラフのプロットがギザギザ（重大度: 中）

**原因**: `quantum_dynamics_complete_comparison.ipynb` の `N_steps = 20` により、プロットデータ点が21点しかない。一方 `quantum_dynamics_gksl_comparison.ipynb` は `n_steps = 100` で101データ点を使用しており、滑らかなプロットとなっている。

**影響を受けるセル**:
- Cell 6: `plot_population_dynamics` / `plot_per_molecule_populations` 関数定義
- Cell 9: Qubit結果の可視化
- Cell 18: Qudit結果の可視化
- Cell 29: ノイズモデルの影響可視化
- Cell 33: 誤差の時間発展プロット
- Cell 34: 3手法の個体数動態比較

**修正**: `N_steps` を 20 から 100 に変更。これにより：
- Trotter近似の精度も向上（dt = 5.0 fs → 1.0 fs）
- プロットのデータ点が21 → 101に増加
- 物理的に正当な改善であり、ヒューリスティックではない

加えて、100データ点でのマーカー密度を適切にするため：
- `markersize` を 5 → 3 に縮小
- `markevery=max(1, len(times)//20)` を追加（約20点ごとにマーカー表示）

### 問題2: 一部セルの出力が生成されない（重大度: 高）

**原因**: Cell 13（Qubit回路可視化）で `matplotlib.use("Agg")` が設定される。Aggバックエンドは非インタラクティブなため、以降のセルで `plt.show()` を呼んでも画面出力が生成されない。

**影響を受けるセル**（すべてCell 13以降）:
- Cell 29: ノイズモデルの影響可視化 → **出力なし**
- Cell 33: 誤差の時間発展プロット → **出力なし**
- Cell 34: 3手法の個体数動態比較 → **出力なし**
- Cell 35: 包括的比較表 → テーブル図の出力なし（テキスト出力は正常）

**注記**: Cell 6の `plot_population_dynamics` はCell 13より前に実行されるため、`plt.show()` が機能する。しかしCell 9以降で再呼び出しされた場合は出力されない。

**修正**: 影響を受けるすべてのセルを以下のパターンに統一：
```python
fig.savefig(filename, dpi=150, bbox_inches="tight")
plt.close(fig)
display(Image(filename=filename))
```
これはGKSLノートブック（`gksl_visualization.py`）と同じパターンであり、Aggバックエンドでも正しく動作する。

**追加修正**: Cell 6 と Cell 9 も同様に `matplotlib.use("Agg")` + `savefig` + `display(Image())` パターンに統一。これにより、ノートブック全体でmatplotlibバックエンドの一貫性を確保。

### 問題3: 基本ゲート分解後の回路可視化が膨大（重大度: 低）

**原因**: Cell 13 で `分解後UnitaryGate版`（KAK分解による）、Cell 26 で `基本ゲート分解版`（CustomTwo→基本ゲート分解）の回路を可視化すると、数百～数千ゲートの膨大な回路図が出力される。

**修正**:
- Cell 13: `分解後UnitaryGate版` の可視化を削除。`UnitaryGate版` と `基本ゲート分解版` のみ表示
- Cell 26: `基本ゲート分解版` の可視化を削除。`CustomTwoゲート版` のみ表示

**注記**: ゲート数の比較分析（Cell 12, Cell 24）はそのまま保持。分解後のゲート数統計は重要な情報であるため、テキスト出力は維持。

---

## 2. 修正内容の一覧

| セル | 修正内容 | 理由 |
|------|---------|------|
| Cell 3 | `N_steps = 20` → `N_steps = 100` | プロットの滑らかさ + Trotter精度向上 |
| Cell 6 | `matplotlib.use("Agg")` + `savefig` + `display(Image())` + `markevery` | 出力確保 + 滑らかプロット |
| Cell 9 | `matplotlib.use("Agg")` + `savefig` + `display(Image())` | 出力確保 |
| Cell 13 | `分解後UnitaryGate版` の可視化を削除 | 膨大な回路出力の抑制 |
| Cell 26 | `基本ゲート分解版` の可視化を削除 | 膨大な回路出力の抑制 |
| Cell 29 | `plt.show()` → `savefig` + `display(Image())` + `markevery` | Agg非互換修正 |
| Cell 33 | `plt.show()` → `savefig` + `display(Image())` + `markevery` | Agg非互換修正 |
| Cell 34 | `plt.show()` → `savefig` + `display(Image())` + `markevery` | Agg非互換修正 |
| Cell 35 | `plt.show()` → `savefig` + `display(Image())` | Agg非互換修正 |

---

## 3. PR216との関連

PR216（#216）で行われた変更:
1. Cell 13 に `from qiskit import QuantumCircuit` を追加
2. `_draw_qiskit_subcircuits` の型アノテーションを `list` → `list[tuple[str, QuantumCircuit]]` に修正

これらの変更自体は直接的に出力を壊すものではないが、PR216のマージに伴いCell 13の `matplotlib.use("Agg")` が既にセットされている状態でノートブックが再実行された結果、問題2が顕在化した。

**根本原因**: Cell 13での `matplotlib.use("Agg")` は以前から存在していたが、以前はCell 13以降のプロットセル（29, 33, 34）が正常に実行されない状態が見過ごされていた可能性がある。PR216でノートブックを再実行した際に問題が明確になった。

---

## 4. GKSLノートブックとの整合性

`quantum_dynamics_gksl_comparison.ipynb` の可視化パターンとの整合性を確認：

| 項目 | GKSL notebook | Complete comparison (修正後) |
|------|--------------|---------------------------|
| matplotlibバックエンド | Agg | Agg（統一） |
| プロット出力方法 | `savefig` + `plt.close` + `display(Image())` | 同一パターン（統一） |
| データ点数 | 100 (`n_steps=100`) | 100 (`N_steps=100`)（統一） |
| マーカーサイズ | markersize=4 | markersize=3 + markevery（適切） |

---

## 5. 検証手順

修正後のノートブックを検証するには：

1. `cd tutorials`
2. ノートブックを全セル実行
3. 以下を確認：
   - Cell 6: 2つのプロットが滑らかに表示される（display_data出力あり）
   - Cell 9: 回路図が表示される（display_data出力あり）
   - Cell 13: 2つの回路図のみ表示（UnitaryGate版 + 基本ゲート分解版）
   - Cell 20: Qudit回路図が表示される
   - Cell 26: 1つの回路図のみ表示（CustomTwoゲート版のみ）
   - Cell 29: ノイズ比較プロットが表示される（display_data出力あり）
   - Cell 33: 誤差プロットが表示される（display_data出力あり）
   - Cell 34: 3手法比較プロットが表示される（display_data出力あり）
   - Cell 35: 比較テーブルが表示される（display_data出力あり）

---

## 6. 次回の検証で確認すべき事項

1. ノートブックの全セルが正常に実行されること
2. すべてのプロットが滑らかに描画されること
3. 回路可視化が適切なサイズで表示されること
4. N_steps=100での計算精度がN_steps=20と比較して改善されていること
5. ノイズモデルの比較が正しく表示されること
