# 検証結果分析: Quantum Dynamics Comparison - Iteration 12

## 1. 背景

### 1.1 Iteration 11の結果

Iteration 11の検証（35/35チェック全PASS）の後、iteration 11分析ドキュメントでは
「コードバグは発見されなかった」と結論していた。しかし、以下の詳細なコードレビューを
実施したところ、quditアンシラ修正（PR#228）の適用漏れが1ファイルで発見された。

### 1.2 Iteration 12で発見された問題

`qudit_gksl_circuit_boson_simulator.py`のStinespring回路構築メソッドで、
アンシラ次元が**d=2（qubit）にハードコード**されていた。

PR#228（quditアンシラ修正）では以下の修正が行われたが：
- `stinespring_utils.py`: d_ancパラメータ追加 ✓
- `qudit_gksl_simulator.py`: d_anc=params.d ✓
- `qudit_gksl_shot_simulator.py`: d_anc=params.d ✓
- `qudit_gksl_noisy_simulator.py`: 親クラスから継承 ✓
- `qudit_gksl_boson_simulator.py`: d_anc=params.d ✓
- `qudit_gksl_circuit_simulator.py`: d_anc=params.d, [d, d_anc], [d, d, d_anc] ✓
- **`qudit_gksl_circuit_boson_simulator.py`: d_anc=params.d ✓ だが回路構築でd_anc未使用 ✗**

回路構築部分が修正漏れとなっていた。

## 2. 発見されたバグの詳細

### 2.0 バグ #0: _trotter_step の非回文構造（Trotter分割不一致）

`qudit_gksl_circuit_boson_simulator.py`の`_trotter_step`が非回文（non-palindromic）
1次Trotter分割を使用していた。これは`QuditGKSLBosonSimulator`の回文2次分割と異なる。

**修正前:**
```python
H(dt/2) → D_1...D_26(dt) → H(dt/2)  # 非回文、full dt
```

**修正後:**
```python
H(dt/2) → D_1...D_26(dt/2) → D_26...D_1(dt/2) → H(dt/2)  # 回文、half dt
```

この修正により`test_matches_matrix_boson_simulator`テストが通過するようになった
（修正前: diff=3.6e-4でFAIL、修正後: diff<1e-10でPASS）。

### 2.1 バグ #1: build_stinespring_circuit_single

**修正前:**
```python
circuit = QuantumCircuit(2, [d, 2], 0)  # アンシラがd=2にハードコード
```

**修正後:**
```python
d_anc = self.d_anc
circuit = QuantumCircuit(2, [d, d_anc], 0)  # アンシラがd_ancを使用
```

### 2.2 バグ #2: build_stinespring_circuit_pair

**修正前:**
```python
circuit = QuantumCircuit(3, [d, d, 2], 0)  # アンシラがd=2にハードコード
```

**修正後:**
```python
d_anc = self.d_anc
circuit = QuantumCircuit(3, [d, d, d_anc], 0)  # アンシラがd_ancを使用
```

### 2.3 バグの影響

`self.d_anc = params.d` で d_anc=3 が設定されているにもかかわらず、
`_build_local_stinespring_unitary` は d_anc=3 を使用して 9×9（single）または
27×27（pair）のユニタリを生成する。一方、回路は [d, 2]=[3, 2] → 6次元 または
[d, d, 2]=[3, 3, 2] → 18次元 で構築されるため、**ユニタリと回路の次元が不一致**する。

| ケース | ユニタリ次元 | 回路次元（修正前） | 回路次元（修正後） |
|--------|-------------|-------------------|-------------------|
| single-site | 9×9 (d_anc×d) | 6 (d×2) ✗ | 9 (d×d_anc) ✓ |
| pair | 27×27 (d_anc×d²) | 18 (d²×2) ✗ | 27 (d²×d_anc) ✓ |

この次元不一致は、`cu_two` または `cu_multi` ゲート適用時にランタイムエラーを
引き起こすか、もし何らかの理由で暗黙の処理が行われた場合は不正な結果を生む。

### 2.4 非ボソン版との比較

`qudit_gksl_circuit_simulator.py`（非ボソン版）は正しく修正されていた：
```python
# 非ボソン版（正しい）
d_anc = self.d_anc
circuit = QuantumCircuit(2, [d, d_anc], 0)    # single: line 322
circuit = QuantumCircuit(3, [d, d, d_anc], 0)  # pair: line 348
```

ボソン版のみが修正漏れとなっていた。

## 3. 修正内容

### 3.1 変更ファイル
- `tutorials/qudit_gksl_circuit_boson_simulator.py`
  - `build_stinespring_circuit_single`: `[d, 2]` → `[d, d_anc]`, `d_anc = self.d_anc` 追加
  - `build_stinespring_circuit_pair`: `[d, d, 2]` → `[d, d, d_anc]`, `d_anc = self.d_anc` 追加

### 3.2 変更なし
- `_build_local_stinespring_unitary`: 既にd_anc使用（変更不要）
- `_extract_kraus_from_local_stinespring`: 既にd_anc使用（変更不要）
- 他の全シミュレータ: 既に正しく修正済み

## 4. 検証結果

`tutorials/run_iteration12_verification.py`を実行: **32/32チェック全てPASS**

主要検証項目：
1. ✓ build_stinespring_circuit_single が [d, d_anc] を使用
2. ✓ build_stinespring_circuit_pair が [d, d, d_anc] を使用
3. ✓ d_anc = self.d_anc がメソッド内で抽出されている
4. ✓ 全quditシミュレータで d_anc = params.d が設定されている
5. ✓ 全qubitシミュレータで d_anc が未設定（デフォルトd_anc=2を使用）
6. ✓ 非ボソン版との整合性確認
7. ✓ 次元計算の数値検証（single: 9×9, pair: 27×27）
8. ✓ 旧バグでの次元不一致の確認（single: 6≠9, pair: 18≠27）
9. ✓ Iteration 11回帰テスト（cx_transfer, depol_001等）

## 5. Iteration 11分析との関係

Iteration 11の検証スクリプトは `qudit_gksl_circuit_boson_simulator.py` の
回路構築部分を直接検証していなかったため、このバグを検出できなかった。
Iteration 12では、全quditシミュレータのd_anc使用を網羅的に検証するチェックを追加した。

## 6. 次のステップ

コードバグは修正済み。ユーザーに以下を推奨：
1. `tutorials/run_iteration12_verification.py` をローカルで実行して検証結果を確認
2. 必要に応じてノートブックを再実行して結果を確認
3. 結果に問題がなければ、次のiteration改善提案（Section 4.1, 4.2 from iteration 11）を検討
