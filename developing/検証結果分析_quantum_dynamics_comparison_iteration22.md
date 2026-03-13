# 検証結果分析: quantum_dynamics_comparison iteration 22

## 1. 概要

Iteration 21の検証結果(42/42 PASS)およびノートブック実行結果を詳細に分析した。
シミュレータコードに致命的なバグは発見されなかったが、**検証スクリプトのボソンパラメータがノートブックと不一致**である問題（Issue AA）を発見・修正した。

### 1.1 分析対象

1. `run_iteration21_verification.py` の実行結果 (`iteration21_cell36_display_20260313T074320Z.json`: 42/42 PASS)
2. `quantum_dynamics_gksl_comparison.ipynb` の全セル実行結果
3. `quantum_dynamics_complete_comparison.ipynb` の全セル実行結果
4. 全6シミュレータのソースコード詳細レビュー

### 1.2 分析結論

- **シミュレータコードにバグなし**: 全6シミュレータ（qudit/qubit × DM/shot/boson）の物理的正確性を確認
- **Issue AA発見**: 検証スクリプトがノートブックと異なるボソンパラメータ（n_max）を使用していた

## 2. 詳細分析

### 2.1 シミュレータコードの物理的正確性

以下の項目を詳細にレビューし、すべて正確であることを確認した：

| 項目 | 結果 | 詳細 |
|------|------|------|
| Stinespring拡張 | ✓ 正確 | d_anc=2(qubit)とd_anc=3(qudit)の両方で正しく実装 |
| 個体群計算 | ✓ 正確 | Kronecker積順序と一致する分解ループ |
| トレース保存 | ✓ 正確 | max|Tr-1| ≈ 1e-14（機械精度内） |
| 禁止状態検出 | ✓ 正確 | noiseless qubitでforbidden_count=0 |
| ゲート数推定 | ✓ 正確 | palindromic 2×因子を含む |
| 物理部分空間保存 | ✓ 正確 | 埋め込み演算子が禁止部分空間を混合しない |
| ノイズモデル | ✓ 正確 | Pauli脱分極、CXオーバーヘッド |

### 2.2 ノートブック実行結果の検証

#### quantum_dynamics_gksl_comparison.ipynb

| セル | 内容 | 結果 | 値 |
|------|------|------|-----|
| Cell 4 | ClassicalGKSL | ✓ | max\|Tr-1\|=5.33e-15 |
| Cell 6 | QuditGKSL DM | ✓ | n_total=30, GPS=66, max\|Tr-1\|=2.02e-14 |
| Cell 8 | QubitGKSL DM | ✓ | n_total=34, GPS=388, max\|Tr-1\|=2.04e-14 |
| Cell 12 | QuditGKSLBoson | ✓ | n_total=16, GPS=38, max\|Tr-1\|=9.99e-15 |
| Cell 14 | QubitGKSLBoson | ✓ | n_total=18, GPS=216, max\|Tr-1\|=1.01e-14 |
| Cell 26 | Fidelity | ✓ | Classical vs Qudit: F=0.999999, Qubit vs Qudit: F=1.000000 |
| Cell 28 | 収束解析 | ✓ | Rate(T)≈1.0 → O(dt)収束（Stinespring近似による） |
| Cell 32 | QuditShot | ✓ | 1000 shots, 結果=DM一致、counts={0:426,...} |
| Cell 34 | QuditNoisyShot | ✓ | p_depol=0.001, 個体群合計=4.0（漏洩なし） |
| Cell 36 | QubitShot | ✓ | forbidden_count=0, counts=Cell 32と完全一致 |
| Cell 38 | QubitNoisyShot | ✓ | forbidden_count=669/1000, Tr=0.326（物理的に正しい漏洩） |

#### quantum_dynamics_complete_comparison.ipynb

| セル | 内容 | 結果 | 値 |
|------|------|------|-----|
| Cell 8 | QubitシミュレーションQiskit | ✓ | GPS=3484, 深さ=1989 |
| Cell 17 | QuditシミュレーションMQT-Qudits | ✓ | GPS=88 |
| Cell 30 | ノイズ影響定量化 | ✓ | Qubit非物理的状態=67.88%, Qudit漏洩なし |
| Cell 35 | 包括的比較表 | ✓ | 全3手法の一貫した比較 |

### 2.3 Issue AA: ボソンパラメータ不一致

**問題**: Iteration 21の検証スクリプトがノートブックと異なるn_maxを使用

| パラメータ | Iteration 21検証 | ノートブック（Cell 2） |
|-----------|----------------|---------------------|
| n_max | 2（デフォルト） | 1 |
| qubit_boson n_total | 20 | 18 |
| qubit_boson GPS | 220 | 216 |
| qubit_boson_N4 n_total | 42 | 38 |
| qubit_boson_N4 GPS | 484 | 476 |

**根本原因**: ノートブックは計算コスト削減のために`n_max=1`を使用しているが、検証スクリプトはGKSLPhysicalParametersのデフォルト`n_max=2`を使用していた。

**影響**: シミュレータコード自体にはバグなし。検証スクリプトは独自の内部整合性チェックとしては正確だったが、ノートブックの設定と不一致なため、ノートブック実行結果を正確に検証できていなかった。

**修正**: `run_iteration22_verification.py`でノートブック一致パラメータ（n_max=1, omega_ph=0.15, g_eph=0.02）を使用。

### 2.4 qudit bosonのn_max非依存性

注目すべき特性として、qudit bosonシミュレータのn_totalとGPSはn_maxに依存しない：

| 計量 | n_max=1 | n_max=2 | 理由 |
|------|---------|---------|------|
| n_phonon_qudits | 2 (=N) | 2 (=N) | 1 qutrit/phonon mode |
| n_total_qudits | 16 | 16 | n_phonon固定 |
| GPS | 38 | 38 | ゲート数はレジスタ数のみに依存 |

一方、qubit bosonではn_ph_qubits = ceil(log2(n_max+1)) * Nがn_maxに依存するため、n_totalとGPSも変化する。

## 3. 修正内容

### 3.1 新規作成ファイル

- `tutorials/run_iteration22_verification.py`: ノートブック一致パラメータでの検証スクリプト（49項目）

### 3.2 検証項目一覧（49項目）

| セクション | 項目数 | 内容 |
|-----------|-------|------|
| Section 1: n_total | 8 | 全シミュレータのn_total値（ボソンはn_max=1） |
| Section 2: 交差整合性 | 4 | ancilla一致、DM/Shot一致 |
| Section 3: GPS回帰 | 8 | 全GPS値（qubit_boson: 216/476に更新） |
| Section 4: ノートブックセル | 13 | Cell 6/8/12/14/32/36のソース検査 |
| Section 5: Result dict | 6 | 戻り値の完全性 |
| Section 6: forbidden_count | 3 | noiseless qubitで0確認 |
| Section 7: n_max非依存性 | 4 | qudit bosonのn_max非依存性チェック（新規） |
| Section 8: Cell 14出力 | 3 | ノートブック出力値の直接検証（新規） |

## 4. 残課題

### 4.1 Qubit回路シミュレータの抽象レベル不一致

既知問題（iteration 17以降）: qubit回路シミュレータは66個のcompositeゲート（UnitaryGate）をカウントするが、qubit DM/shotは388個のbasicゲートを推定する。これらは異なる抽象レベルの計測であり、公正な比較のための統一が必要。

### 4.2 次の検証手順

1. ユーザーが`run_iteration22_verification.py`をローカルで実行し、結果をpush
2. 結果を確認し、必要に応じて追加修正

## 5. 技術メモ

### 5.1 CXゲート数の差異

GKSL比較ノートブック（hbar=1.0）とcomplete比較ノートブック（hbar=0.6582）でCXゲート数が異なる：
- GKSL: H_transfer=46 CXゲート（回転角=V*dt/1.0=0.1 rad）
- Complete: H_transfer=47 CXゲート（回転角=V*dt/0.6582=0.152 rad）

これは異なる単位系によるものであり、バグではない。GKSLシミュレータは自然単位系（hbar=1）を使用し、completeシミュレータは物理単位系（hbar=0.6582 eV·fs）を使用している。

### 5.2 Stinespring拡張のd_anc依存性

quditシミュレータはd_anc=3（qutrit ancilla）を使用し、qubitシミュレータはd_anc=2（qubit ancilla）を使用する。d_anc=3の場合、第3Kraus演算子K_2=0であるため、物理的には同等。これは設計上の意図であり、それぞれのハードウェアアーキテクチャを正しくモデル化している。
