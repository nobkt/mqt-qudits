# 検証結果分析: Quantum Dynamics Comparison - Iteration 10

## 1. 背景

### 1.1 Iteration 9の状態

Iteration 9で以下の修正が行われた：

- `qubit_noisy_simulator.py`のdocstringで`cx_per_pair_gate: int`→`float`に修正
- GKSLノートブック Cell 40のF_rawに関する記述を「misleadingly low」→「correctly includes leakage penalty」に修正
- 検証: 27/27チェック全てPASS

### 1.2 Iteration 10で発見された問題点

Iteration 9の検証結果とノートブック実行結果を注意深く分析した結果、以下の重要な物理的問題を発見した。

**問題: GKSLノートブック Cell 38での`cx_per_pair_gate`にH_TTAのCXゲート数が含まれていた**

GKSLノートブックのCell 38で、`cx_per_pair_gate`に`cx_info['cx_avg']`（H_transferとH_TTAの平均CXゲート数）を使用していた。しかし、GKSLモデルのハミルトニアンは以下の構造である：

```python
self.H_0 = build_onsite_hamiltonian(params)
self.H_transfer = build_transfer_hamiltonian(params)
self.H_total_qutrit = self.H_0 + self.H_transfer  # H_TTAなし！
```

**GKSLモデルではH_TTAはハミルトニアンに含まれない。** TTAプロセスはLindblad演算子として実装される。
したがって、ハミルトニアンのCXゲートオーバーヘッドにH_TTAのCXゲート数を含めるのは物理的に不正確である。

### 1.3 問題の定量的影響

| パラメータ | 修正前（cx_avg使用） | 修正後（cx_transfer使用） |
|-----------|-------------------|----------------------|
| cx_per_pair_gate | 66.0 (avg of 46, 86) | 46 (H_transfer only) |
| p_eff | 6.39% | 4.50% |
| ノイズ過大評価率 | 42%過大 | 正確 |

修正前は、GKSLモデルに存在しないH_TTA（86 CXゲート）のCXゲート数を含めることで、
有効脱分極率が42%過大に評価されていた。

### 1.4 2つのノートブックの正しい使い分け

| 項目 | Complete Comparison | GKSL Comparison |
|------|-------------------|-----------------|
| ハミルトニアン | H_0 + H_transfer + H_TTA | H_0 + H_transfer |
| TTA実装 | ユニタリ（ハミルトニアン） | 散逸的（Lindblad演算子） |
| cx_per_pair_gate | cx_avg（transfer+TTAの平均） ✓ | cx_transfer（H_transferのみ） ✓ |
| 物理的根拠 | 両方のペア相互作用がゲート | H_transferのみがゲート |

## 2. Iteration 10での修正内容

### 2.1 GKSLノートブック Cell 38の修正

**修正前:**
```python
cx_per_pair = cx_info['cx_avg']
```

**修正後:**
```python
# GKSLモデルのハミルトニアンは H_0 + H_transfer のみ（H_TTAなし）。
# TTAはLindblad演算子で実装されるため、ハミルトニアンのCXオーバーヘッドは
# H_transferの分解結果のみを使用する。
cx_per_pair = cx_info['cx_transfer']
```

### 2.2 GKSLノートブック Cell 37のマークダウン更新

GKSLモデルにおけるCXゲート数の選択理由を追記：

- GKSLハミルトニアンがH_0 + H_transferのみを含むこと
- TTAプロセスがLindblad演算子で実装されること
- Lindbladペアチャネルの正確なCXゲート数はStinespring拡張の分解から推定が必要であること
- 現時点ではH_transferのCXゲート数を近似値として使用すること

### 2.3 Cell 38の出力形式の改善

H_TTAのCXゲート数は参考値として表示し、実際に使用する値（cx_transfer）を明示する：

```
CXゲート数推定結果:
  H_transfer: 46 CXゲート
  H_TTA: 86 CXゲート (参考値: GKSLではLindblad実装)
  使用値: 46 CXゲート/ペア (H_transferのみ、GKSLハミルトニアンに基づく)
```

## 3. 技術的分析

### 3.1 GKSLシミュレータのノイズイベント構造

`_trotter_step_trajectory`メソッドの各ステップにおけるノイズ適用箇所：

1. **半ステップハミルトニアン** (`_U_H_half = expm(-i*(H_0+H_transfer)*dt/2)`)
   - 3ペアノイズイベント（H_transferのみ）
   - 使用率: p_eff = 1-(1-p_phys)^cx_transfer

2. **順方向Lindbladチャネル**（6ペアチャネル: 3近傍×2 TTA演算子）
   - 6ペアノイズイベント
   - 使用率: p_eff = 1-(1-p_phys)^cx_transfer（近似値）

3. **逆方向Lindbladチャネル**（回文的順序）
   - 6ペアノイズイベント

4. **半ステップハミルトニアン**
   - 3ペアノイズイベント

**合計: 18ペアノイズイベント/ステップ**

### 3.2 Lindbladペアチャネルのノイズに関する議論

LindbladのTTAペアチャネルは、Stinespring拡張により実装される：

```
U_Stinespring = stinespring_unitary_from_lindblad(L_TTA, dt/2)
```

このStinespring拡張ユニタリはシステム+補助ビットの空間で動作し、
qubitハードウェア上ではH_transferとは異なる数のCXゲートが必要になる可能性がある。

現時点では以下の理由によりH_transferのCXゲート数を近似値として使用：
- Stinespring拡張のCXゲート数推定機能が未実装
- H_transferのCXゲート数はGKSLモデルの文脈から得られた妥当な値
- 以前のcx_avg使用（H_TTAを含む）よりも物理的に正確

### 3.3 Complete Comparisonノートブックへの影響

Complete Comparisonノートブック（Cell 15）は変更なし。理由：
- Complete Comparisonのハミルトニアンは H_0 + H_transfer + H_TTA を含む
- H_transferとH_TTAの両方がペアゲートとしてノイズを受ける
- cx_avg使用は物理的に正確

## 4. 検証結果

### 4.1 Iteration 10検証スクリプト結果

- 検証スクリプト: `run_iteration10_verification.py`
- 結果: 22/22チェック全てPASS（CXゲート推定はQiskit依存のためSKIP）
- 回帰テスト: Iteration 9のチェックも全てPASS

### 4.2 予想される実行結果の変化

ノートブック再実行後、以下の変化が予想される：

| 指標 | 修正前（cx=66.0） | 修正後（cx≈46） |
|------|-----------------|---------------|
| p_eff | 6.39% | ~4.50% |
| 最終N_T1 | ~0.41 | やや増加（ノイズ減少） |
| Tr(rho_final) | ~0.31 | やや増加 |
| F_raw | ~0.034 | やや増加 |
| F_norm | ~0.110 | やや増加 |

ノイズが減少するため、物理的状態の保持率が改善され、
qubitとquditの性能差は依然として大きいものの、以前の比較よりもやや縮小する見込み。

## 5. 残存する検討事項

### 5.1 Lindbladペアチャネルの正確なCXゲート数推定

Stinespring拡張ユニタリのCXゲート数を推定する機能の実装が望ましい。
これにより、Lindbladペアチャネルのノイズ率をより正確にモデル化できる。

### 5.2 QubitGKSLNoisyShotSimulatorクラスの拡張可能性

現在の`cx_per_pair_gate`パラメータは全ペア相互作用に同じ値を使用する。
将来的に、ハミルトニアンペアゲートとLindbladペアゲートに異なるCX数を
指定できるよう拡張することが考えられる。

## 6. 次のステップ

ユーザーが以下を実行する必要がある：
1. `cd tutorials && python run_iteration10_verification.py` で検証
2. `quantum_dynamics_gksl_comparison.ipynb` を再実行（修正されたcx_per_pair_gateで）
3. 検証結果をリポジトリにpush
4. 実行結果を確認し、さらなる修正が必要かどうか判断
