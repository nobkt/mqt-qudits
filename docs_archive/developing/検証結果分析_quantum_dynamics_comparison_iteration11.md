# 検証結果分析: Quantum Dynamics Comparison - Iteration 11

## 1. 背景

### 1.1 Iteration 10の修正内容

Iteration 10で以下の修正が行われた：
- GKSLノートブック Cell 38の`cx_per_pair_gate`を`cx_info['cx_avg']`（66.0）から`cx_info['cx_transfer']`（46）に変更
- GKSLノートブック Cell 37のドキュメント更新（GKSLハミルトニアンがH_0 + H_transferのみであることの説明追記）
- 検証結果: 26/26チェック全てPASS（iteration10_gksl_cx_fix_20260311T054951Z.json）

### 1.2 Iteration 11の目的

Iteration 10の修正後にノートブックを再実行した結果を注意深く分析し、
コード・ノートブック・物理モデルに問題がないか包括的に調査する。

## 2. 検証結果の分析

### 2.1 GKSLノートブック実行結果

#### Cell 34: Qudit Shot+Noise (p_depol=0.001)
```
Final populations: N_S0=2.6132, N_T1=0.7746, N_S1=0.6122
Sum = 4.0000 ✓（トレース保存）
```

#### Cell 38: Qubit Shot+Noise (p_depol=0.001, cx_per_pair_gate=46)
```
p_depol_pair_eff = 0.0450
Final populations: N_S0=0.4440, N_T1=0.4287, N_S1=0.4297
Tr(rho_final) = 0.325596（67.4%がforbidden状態にリーク）
Forbidden state count: 669/1000
```

#### Cell 40: フィデリティ比較
```
Qudit DM vs Qudit Shot (no noise): F = 0.998056 ✓
Qubit DM vs Qubit Shot (no noise): F = 0.998056 ✓
Qudit noisy: Tr=1.000000, F_raw=0.534629, F_norm=0.534629
Qubit noisy: Tr=0.325596, F_raw=0.032488, F_norm=0.099780
```

#### Cell 42: 全シナリオ比較
```
| Scenario               | N_T1   | N_S1   | N_S0   | Sum    |
|------------------------|--------|--------|--------|--------|
| GKSL DM (Classical)    | 0.4184 | 0.3054 | 3.2762 | 4.0000 |
| Qudit DM (Stinespring) | 0.4172 | 0.3053 | 3.2774 | 3.9999 |
| Qubit DM (Stinespring) | 0.4172 | 0.3053 | 3.2774 | 3.9999 |
| Qudit Shot (1000)      | 0.4600 | 0.3000 | 3.2400 | 4.0000 |
| Qubit Shot (1000)      | 0.4600 | 0.3000 | 3.2400 | 4.0000 |
| Qudit Shot+Noise       | 0.7746 | 0.6122 | 2.6132 | 4.0000 |
| Qubit Shot+Noise       | 0.4287 | 0.4297 | 0.4440 | 1.3024 |
```

### 2.2 Complete Comparisonノートブック実行結果

#### Cell 15: Qubit Noisy (密度行列, cx_avg=66.5)
```
p_eff = 6.4368%
ペアゲート数/ステップ: 12
Final: N_S0=0.4340, N_T1=0.4324, N_S1=0.4264
非物理的状態: 0.6768 (67.7%)
```

#### Cell 23: Qudit Noisy (密度行列, p_depol=0.001)
```
2-quditゲート数/ステップ: 12
Final: N_S0=1.9142, N_T1=1.2790, N_S1=0.8068
Sum = 4.0000 ✓
```

### 2.3 クロスノートブック整合性検証

| 項目 | Complete | GKSL | 整合性 |
|------|----------|------|--------|
| Qubit cx_per_pair_gate | 66.5 (cx_avg) | 46 (cx_transfer) | ✓ 各モデルに適切 |
| Qubit p_eff | 6.44% | 4.50% | ✓ CX数の差に対応 |
| ペアノイズイベント/ステップ | 12 | 18 | ✓ モデル構造の差 |
| 総ノイズ/ステップ | 0.773 | 0.810 | ✓ 同程度 |
| Qubit forbidden率 | 67.7% | 67.4% | ✓ 一致 |
| Qudit トレース保存 | ✓ (Sum=4.0) | ✓ (Sum=4.0) | ✓ |

### 2.4 ノイズイベント構造の詳細検証

#### Complete Comparison (qubit_noisy_simulator.py)
```
Forward half:
  H0: 4 single-molecule unitaries (no noise)
  H_transfer: 3 pair gates + 3 noise events (p_eff=6.44%)
  H_TTA: 3 pair gates + 3 noise events (p_eff=6.44%)
Backward half:
  H_TTA: 3 pair gates + 3 noise events
  H_transfer: 3 pair gates + 3 noise events
  H0: 4 single-molecule unitaries (no noise)
Total: 12 pair noise events/step
```

#### GKSL Comparison (qubit_gksl_shot_simulator.py)
```
Half Hamiltonian (H_0 + H_transfer):
  1 full unitary + 3 pair noise events (p_eff=4.50%)
Forward Lindblad (26 channels, depol_pair_only=True):
  6 TTA pair channels: 6 noise events (p_eff=4.50%)
  20 single-site channels: noise skipped
Reverse Lindblad (palindromic):
  6 TTA pair channels: 6 noise events (p_eff=4.50%)
  20 single-site channels: noise skipped
Half Hamiltonian:
  1 full unitary + 3 pair noise events (p_eff=4.50%)
Total: 18 pair noise events/step
```

### 2.5 物理的整合性検証

1. **個体数保存** (Qudit):
   - Complete: Sum=4.0000 ✓
   - GKSL: Sum=4.0000 ✓
   - トレース保存が正しく機能している

2. **Forbidden状態リーク** (Qubit):
   - Complete: 67.7%（12 events × p_eff=6.44% → 77.3 expected errors/100steps）
   - GKSL: 67.4%（18 events × p_eff=4.50% → 81.0 expected errors/100steps）
   - 類似した総ノイズレベル → 類似したリーク率 ✓

3. **フィデリティ計算**:
   - F_norm = F_raw / Tr(ρ) = 0.032488 / 0.325596 = 0.099780 ✓
   - Qudit: F_raw = F_norm = 0.534629 (Tr=1なので同一) ✓

4. **無ノイズ結果の一致**:
   - Qudit DM = Qubit DM (no noise): N_T1=0.4172 ✓
   - Qudit Shot = Qubit Shot (no noise): N_T1=0.4600 ✓（ショットノイズのみ）

5. **DM vs Shot (no noise) フィデリティ**:
   - F = 0.998056 for both qudit and qubit ✓（1000 shots → ショットノイズ ~0.002）

## 3. 分析結果

### 3.1 コードバグの有無

**コードバグは発見されなかった。** 以下の項目を全て検証した：

1. ✓ GKSL Cell 38: `cx_per_pair = cx_info['cx_transfer']` を正しく使用
2. ✓ Complete Cell 15: `cx_per_pair = cx_info['cx_avg']` を正しく使用
3. ✓ 各モデルのCX数選択が物理的に正当
4. ✓ ノイズイベント数がモデル構造と一致（GKSL: 18, Complete: 12）
5. ✓ 個体数保存（qudit）とトレース欠損（qubit）が正しく計算
6. ✓ フィデリティ計算が正確
7. ✓ Lindblad演算子の順序と`_compute_lindblad_sites`のマッピングが一致
8. ✓ `extract_statevector_from_qubit_space`が正しくphysical subspaceのみ抽出
9. ✓ ペア脱分極チャネルの数式実装が正確
10. ✓ Weyl-Heisenberg演算子（qudit）と4-qubit Pauli（qubit）の実装が正確

### 3.2 既知の近似とその影響

以下は意図的な近似であり、バグではない：

1. **Lindblad ペアチャネルのCXゲート数近似**
   - 現状: H_transferのCX数（46）をLindblad TTAペアチャネルにも使用
   - 実際: Stinespring拡張ユニタリは異なるCX数が必要（推定100-200以上）
   - 影響: Qubitのノイズを過小評価（実際はさらに悪化する）
   - ドキュメント: Iteration 10分析 Section 3.2で記載済み

2. **Qudit Lindbladチャネルのネイティブゲート数近似**
   - 現状: 各Lindbladペアチャネルを1ネイティブゲートとして扱う
   - 実際: Stinespring拡張ユニタリは複数のネイティブゲートが必要
   - 影響: Quditのノイズも過小評価（ただしqubitほど大きくない）
   - 注: quditのネイティブゲート優位性は保持される

3. **Complete Comparisonの平均CXゲート数**
   - 現状: cx_avg = (cx_transfer + cx_tta) / 2 = 66.5を全ペアに使用
   - 実際: H_transfer=47, H_TTA=86と異なる
   - 影響: H_transferは過大ノイズ、H_TTAは過小ノイズだが平均は正確

### 3.3 結果の物理的妥当性

全ての結果が物理的に妥当である：

- **Qubit 70%リーク**: p_eff=4.5%, 18イベント/ステップ × 100ステップ = 1800イベントでは
  事実上全トラジェクトリにエラーが発生し、大量のforbidden状態リークが予想される
- **Qubit nearly maximally mixed**: 物理空間内の個体数がほぼ均等（~0.43）は、
  十分なノイズで完全混合状態に収束する挙動と一致
- **Qudit F_raw=0.535**: 1800イベント × 0.1%で約1.8エラー。
  約17%のトラジェクトリがエラーフリー（F=1に寄与）、残りはエラーあり（F<1）。
  平均F_raw≈0.535は妥当

## 4. 今後の改善提案

### 4.1 Stinespring拡張ユニタリのCXゲート数推定（優先度: 中）

Lindbladペアチャネルの正確なCXゲート数を推定する機能の実装：
- Stinespring拡張ユニタリは2分子+アンシラ = 5 qubit (32×32) ユニタリ
- Qiskitの`transpile`で分解してCXゲート数をカウント可能
- これにより、ハミルトニアンペアゲートとLindbladペアチャネルに異なるp_effを適用可能

### 4.2 ペアごとの個別CXゲート数使用（優先度: 低）

Complete ComparisonでH_transferとH_TTAに異なるp_effを適用：
- 現状: cx_avg=66.5を全ペアに使用
- 改善: H_transferには47、H_TTAには86を個別適用
- シミュレータの`simulate()`にdict形式のcx_per_pair_gateを受け取る拡張が必要

## 5. Iteration 11検証スクリプト

`tutorials/run_iteration11_verification.py`を作成。以下を検証：
1. Iteration 10のチェック回帰テスト（21項目）
2. ノートブック出力の数値整合性チェック
3. ノイズイベント数の検証
4. シミュレータコードの深い整合性チェック
5. 脱分極チャネル数式の正確性テスト

## 6. 結論

**Iteration 10の修正は正しく適用されており、コードバグは発見されなかった。**
シミュレーション結果は物理的に妥当であり、全ての検証チェックがPASSしている。

既知の近似（Lindbladチャネルのノイズ率）はドキュメントされており、
qubitとquditの比較の定性的結論に影響を与えない。
むしろ、より正確なモデルではqubitの性能がさらに悪化するため、
quditのネイティブ量子ゲート優位性はさらに強くなる。

## 7. 次のステップ

現時点でコードバグは発見されておらず、即座の修正は不要。
以下の改善を将来の反復で検討可能：

1. Stinespring CX推定の実装（Section 4.1）
2. ペア個別CX数の実装（Section 4.2）
3. ノートブックの追加ドキュメント（近似の詳細説明）

ユーザーが上記改善のいずれかを希望する場合、次の反復で対応する。
