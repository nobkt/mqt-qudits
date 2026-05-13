# 検証結果分析: Quantum Dynamics Comparison - Iteration 13

## 1. 背景

### 1.1 Iteration 12の結果

Iteration 12では`qudit_gksl_circuit_boson_simulator.py`の3つのバグが修正された：
- Bug #0: `_trotter_step`の非回文構造（回文2次Trotterに修正）
- Bug #1: `build_stinespring_circuit_single`のアンシラd=2ハードコード
- Bug #2: `build_stinespring_circuit_pair`のアンシラd=2ハードコード

しかし、**同一のTrotterステップバグが非ボソン版の回路シミュレータにも存在していた**
ことが見落とされていた。

### 1.2 Iteration 13で発見された問題

Iteration 12のPR#229（closed, not merged）および検証結果を詳細に分析した結果、
以下の6つのバグが発見された。

## 2. 発見されたバグの詳細

### 2.1 Bug A: `qudit_gksl_circuit_simulator.py` の非回文Trotterステップ

**問題**: `_trotter_step_circuit`が非回文の1次Trotterを使用

```python
# 修正前（非回文、1次精度）
H(dt/2) → D_1...D_n(dt) → H(dt/2)
```

行列版の`QuditGKSLSimulator._trotter_step`は回文2次Trotterを使用：
```python
# 行列版（回文、2次精度）
H(dt/2) → D_1...D_n(dt/2) → D_n...D_1(dt/2) → H(dt/2)
```

**修正後**:
```python
# 回文、2次精度（行列版と一致）
H(dt/2) → D_1...D_n(dt/2) → D_n...D_1(dt/2) → H(dt/2)
```

**影響**: Iteration 12のBug #0（ボソン版）と完全に同一のバグパターン。
テストが`1e-3`の緩い許容誤差を使用していたため、ハミルトニアンTrotter誤差
（O(dt²) ≈ 3.5e-5）に隠れて見落とされていた。

### 2.2 Bug B: `qubit_gksl_circuit_simulator.py` の非回文Trotterステップ

Bug Aと同一のパターン。qubit版の`_trotter_step_circuit`も非回文だった。
行列版`QubitGKSLSimulator._trotter_step`は回文を使用。

### 2.3 Bug C: `qudit_gksl_circuit_boson_simulator.py` のゲートカウント

**問題**: 回文Trotterは各Lindbladチャネルを2回適用（forward + reverse）するが、
ゲートカウントは1回分しか計上していなかった。

```python
# 修正前
gates_per_step = 2 * gates_per_half_ham + n_stinespring_single + n_stinespring_pair

# 修正後
gates_per_step = 2 * gates_per_half_ham + 2 * n_stinespring_single + 2 * n_stinespring_pair
```

**具体例**（N=2, boson, g_eph=0.005）：
- 修正前: gates_per_step = 26（Lindblad 12ゲート）
- 修正後: gates_per_step = 38（Lindblad 24ゲート）

### 2.4 Bug D/E: qudit/qubit回路シミュレータのゲートカウント

Bug Cと同一パターン。Bug A/Bの回文修正後、Lindbladゲート数も2倍に修正。

**具体例**（N=4, default）：
- 修正前: gates_per_step = 40（Lindblad 26ゲート）
- 修正後: gates_per_step = 66（Lindblad 52ゲート）

### 2.5 Bug F: `build_full_trotter_step_circuit` の非回文構造

両回路シミュレータの`build_full_trotter_step_circuit`メソッドが非回文構造で
回路を構築していた。実際のシミュレーション（`_trotter_step_circuit`）とは
異なる構造の回路が生成されていた。

修正: 回文構造（forward + reverse circuits）に更新。
docstringも `H(dt/2) → D_1...D_n(dt/2) → D_n...D_1(dt/2) → H(dt/2)` に更新。

## 3. バグ発見の経緯

### 3.1 なぜIteration 12で見落とされたか

1. Iteration 12の修正はボソン版（`qudit_gksl_circuit_boson_simulator.py`）のみ対象
2. 非ボソン版のテスト`test_circuit_matches_matrix_simulator`は`1e-3`許容誤差を使用
3. 実際の誤差は`~3.5e-5`で`1e-3`内に収まっていたため、バグが検出されなかった
4. ゲートカウントのバグは、ボソン版のテストがStinespringゲート数を
   直接検証していなかったため検出されなかった

### 3.2 N=2での特殊事情

N=2の場合、ハミルトニアンの各項（on-site, transfer）が可換であるため、
ハミルトニアンTrotterの追加誤差がゼロになる。このため：

- ボソン版テスト（N=2）: 回路と行列の差が`1e-17`（機械精度）
- 非ボソン版テスト（N=4, default）: 差が`3.5e-5`（ハミルトニアンTrotter誤差）

回文修正後は、両方とも（同じNであれば）同一の精度を示す。

## 4. 修正内容

### 4.1 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `qudit_gksl_circuit_simulator.py` | `_trotter_step_circuit`: dt→dt/2 + 回文化; `build_full_trotter_step_circuit`: 回文化; `simulate()`: ゲートカウント2× |
| `qubit_gksl_circuit_simulator.py` | 同上 |
| `qudit_gksl_circuit_boson_simulator.py` | `simulate()`: ゲートカウント2×のみ（Trotterは既に回文） |
| `test_gksl_simulators.py` | 許容誤差1e-3→1e-4, ゲートカウント検証値更新 |

### 4.2 変更なし
- `qudit_gksl_simulator.py`: 既に回文（変更不要）
- `qubit_gksl_simulator.py`: 既に回文（変更不要）
- `qudit_gksl_boson_simulator.py`: 既に回文（変更不要）
- `qudit_gksl_circuit_boson_simulator.py` `_trotter_step`: Iteration 12で修正済み
- `stinespring_utils.py`: 変更不要

## 5. 検証結果

`tutorials/run_iteration13_verification.py`を実行: **33/33チェック全てPASS**

### 5.1 主要検証項目

| カテゴリ | チェック数 | 結果 |
|---------|-----------|------|
| Qudit回路回文Trotter | 4 | 全PASS |
| Qubit回路回文Trotter | 4 | 全PASS |
| ゲートカウント修正 | 6 | 全PASS |
| build_full回文構造 | 4 | 全PASS |
| ボソンゲートカウント | 4 | 全PASS |
| Iteration 12回帰 | 3 | 全PASS |
| 数値検証 | 8 | 全PASS |

### 5.2 数値検証結果

| テスト | パラメータ | 差分 | 判定 |
|--------|-----------|------|------|
| Qudit回路 vs 行列 (N=2) | t_max=5.0, n_steps=5 | 2.50e-16 | 機械精度 ✓ |
| Qudit回路 vs 行列 (N=4) | t_max=5.0, n_steps=5 | 3.55e-05 | O(dt²) ✓ |
| Qubit回路 vs 行列 (N=2) | t_max=5.0, n_steps=5 | 2.50e-16 | 機械精度 ✓ |
| ボソン回路 vs 行列 (N=2) | t_max=2.0, n_steps=5 | 2.78e-17 | 機械精度 ✓ |

N=4での`3.55e-05`差はハミルトニアンTrotter分割（回路分解固有の近似）に起因し、
バグではない。dt→0で収束する。

## 6. テスト更新

| テスト名 | 変更前 | 変更後 |
|---------|--------|--------|
| test_circuit_matches_matrix_simulator (qudit) | 1e-3 | 1e-4 |
| test_circuit_matches_matrix_simulator (qubit) | 1e-3 | 1e-4 |
| test_gate_breakdown (qudit) | n_single, n_pair | 2*n_single, 2*n_pair |
| test_gate_breakdown (qubit) | n_single, n_pair | 2*n_single, 2*n_pair |
| test_build_full_trotter_step_circuit (qudit) | 40, 26 | 66, 52 |
| test_build_full_trotter_step_circuit (qubit) | 40, 26 | 66, 52 |

## 7. 次のステップ

1. ユーザーが`tutorials/run_iteration13_verification.py`をローカルで実行して結果を確認
2. 必要に応じてノートブックを再実行して結果を確認
3. ノートブック内のゲートカウント関連のテキストが更新が必要かどうかを確認
