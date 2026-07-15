# 検証結果分析: Quantum Dynamics Comparison - Iteration 9

## 1. 背景

### 1.1 Iteration 8の状態

Iteration 8で以下の修正が行われた：

- `qubit_noisy_simulator.py`の`cx_per_pair_gate`型ヒントを`int`→`float`に修正
- `quantum_dynamics_complete_comparison.ipynb` Cell 15に`cx_per_pair_gate`を正しく渡す修正
- `qubit_gksl_shot_simulator.py`の`QubitGKSLNoisyShotSimulator`に`cx_per_pair_gate`サポートを追加
- `quantum_dynamics_gksl_comparison.ipynb` Cell 37-38にCXオーバーヘッドの説明・実装を追加

### 1.2 Iteration 8の検証結果

- 検証スクリプト: 33/33チェック全てPASS
- CXゲート数: H_transfer=47, H_TTA=86, 平均66.5（complete）/ H_transfer=46, H_TTA=86, 平均66.0（GKSL）
- 有効ペアエラー率: p_eff = 6.44%（p_phys=0.001, cx=66.5の場合）

### 1.3 Iteration 9で発見された問題点

Iteration 8の検証結果とノートブック実行結果を注意深く分析した結果、以下の問題を発見：

**問題1: `qubit_noisy_simulator.py`のdocstringに残存する型不整合**

`simulate()`メソッドの型ヒントは正しく`cx_per_pair_gate: float = 1`（line 351）に修正されていたが、
同メソッドのdocstring（line 382）では`cx_per_pair_gate : int`のままだった。
`estimate_cx_per_pair_gate()`が返す`cx_avg`はfloat（例: 66.5）であるため、
docstringも`float`に修正する必要がある。

**問題2: GKSL比較ノートブック Cell 40のF_raw忠実度に関する不正確な記述**

Cell 40の出力で以下のように表示されていた：

```
Note: F_raw for qubit noisy is misleadingly low due to forbidden-state leakage (Tr<1).
```

この記述は事実と異なる。F_rawは「misleadingly low（誤解を招くほど低い）」のではなく、
禁止状態への漏洩によるペナルティを正確に含んだ物理的に正しい忠実度指標である。

- F_raw = 0.034: 禁止状態漏洩（69%）+ 物理部分空間内の脱分極の両方を正しく反映
- F_norm = 0.110: 物理部分空間内に残った状態のみの条件付き忠実度

F_rawが低いのは「誤解を招く」のではなく、qubitエンコーディングの実際の性能劣化を正確に示している。

## 2. Iteration 9での修正内容

### 2.1 `qubit_noisy_simulator.py`のdocstring修正

**修正前（line 382）:**
```
cx_per_pair_gate : int
    Number of CX gates per pair interaction on qubit hardware.
    Default 1 gives same noise level as qudit.
    Use estimate_cx_per_pair_gate() to compute the actual value.
```

**修正後:**
```
cx_per_pair_gate : float
    Number of CX gates per pair interaction on qubit hardware.
    Accepts float values for averaged counts (e.g. 66.5).
    Default 1 gives same noise level as qudit.
    Use estimate_cx_per_pair_gate() to compute the actual value.
```

### 2.2 GKSLノートブック Cell 40の記述修正

**修正前:**
```
Note: F_raw for qubit noisy is misleadingly low due to forbidden-state leakage (Tr<1).
F_norm normalizes the sub-normalized matrix first, giving the conditional fidelity
(fidelity of the state conditioned on remaining in the physical subspace).
```

**修正後:**
```
Note: For qubit noisy, F_raw correctly includes the forbidden-state leakage penalty (Tr<1).
F_norm normalizes the sub-normalized matrix first, giving the conditional fidelity
(fidelity conditioned on remaining in the physical subspace, removing the leakage penalty).
```

## 3. ノイズモデルの検証結果

### 3.1 脱分極チャネルの正確性

両シミュレータの脱分極チャネルを数学的に検証した：

**密度行列シミュレータ（qubit_noisy_simulator.py）:**
```
E(ρ) = (1-p)ρ + p·Tr_{ij}(ρ) ⊗ I_{16}/16
```
ペア部分空間の完全脱分極チャネル。トレース保存性を満たす。

**ショットベースシミュレータ（qubit_gksl_shot_simulator.py）:**
```
- 確率 (1-p+p/256): 恒等操作（エラーなし）
- 確率 p/256: 各非恒等4-qubitパウリ演算子を適用
```
256個のパウリ演算子の確率的適用。期待値は密度行列版と一致：
(1/256)Σ_k P_k ρ P_k† = Tr_{pair}(ρ) ⊗ I_{16}/16

### 3.2 トレース保存性の検証

**Quditシミュレータ:**
- Tr(ρ_final) = 1.000000（完全保存）
- 理由: 全81状態が物理的。脱分極ノイズは物理部分空間内に留まる。

**Qubitシミュレータ:**
- Tr(ρ_final in qutrit space) = 0.311811
- 理由: 全256状態のうち81状態のみが物理的。最大混合状態 I/256 の物理部分空間射影のトレースは 81/256 = 0.3164。
- 観測値0.312は理論値0.316に非常に近く、系がほぼ完全に最大混合状態に達したことを示す。
- **注意**: 完全4^N空間でのTr(ρ) = 1.0は保存されている。射影によるTr < 1は禁止状態への漏洩を正しく反映。

### 3.3 ノイズイベント数の分析

**Complete comparison (DM simulator):**
- 前半: 3(transfer) + 3(TTA) = 6ペアノイズイベント
- 後半: 3(transfer) + 3(TTA) = 6ペアノイズイベント
- 合計: 12ペアノイズイベント/ステップ
- 1ステップあたりの有効ノイズ: 1-(1-0.0644)^12 ≈ 55.0%

**GKSL comparison (Shot simulator, depol_pair_only=True):**
- H_half後: 3ペアノイズイベント
- Lindblad前半: 6ペアチャネルノイズイベント（3近傍×2チャネル）
- Lindblad後半: 6ペアチャネルノイズイベント（回文的）
- H_half後: 3ペアノイズイベント
- 合計: 18ペアノイズイベント/ステップ

GKSLシミュレータはLindblad散逸チャネルの追加ペア相互作用により、
DMシミュレータより多くのノイズイベントを持つ。これは異なる物理モデルの
正しい反映であり、バグではない。

## 4. 全体の正確性評価

### 4.1 正しく動作している項目

1. **CXゲートオーバーヘッド計算**: p_eff = 1-(1-p_phys)^cx が正しく実装
2. **脱分極チャネル**: 密度行列版とショットベース版で数学的に等価
3. **禁止状態追跡**: qubitの非物理的状態が正確にカウントされている
4. **後方互換性**: cx_per_pair_gate=1でp_eff=p_physに還元
5. **Quditのネイティブゲート**: cx_per_pair_gateパラメータなし（正しい）
6. **CXゲート数推定**: Qiskitトランスパイラによる実測値使用
7. **パラメータ整合性**: 両ノートブックでdepol=0.001を使用

### 4.2 修正された項目

1. **docstring型**: `int` → `float`（機能影響なし、ドキュメント修正のみ）
2. **F_raw記述**: 「misleadingly low」→ 正確な物理的解釈に修正

### 4.3 残存するリスクなし

- コードの物理モデルは正確
- ノイズチャネルの数学的実装は正しい
- ノートブック実行結果は期待される物理に一致

## 5. 次のステップ

ユーザーが以下を実行する必要がある：
1. `cd tutorials && python run_iteration9_verification.py` で検証
2. `quantum_dynamics_gksl_comparison.ipynb` を再実行（Cell 40の修正テキスト確認）
3. 検証結果をリポジトリにpush
4. 実行結果を確認し、さらなる修正が必要かどうか判断
