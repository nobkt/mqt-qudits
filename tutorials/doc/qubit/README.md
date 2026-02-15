# Qubitによる分子三重項状態量子ダイナミクス ドキュメント

## 概要

本ディレクトリには、`tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`で実装されたQudit（Qutrit）ベースの分子三重項状態量子ダイナミクスシミュレーションと同等の計算を、**Qubit（2準位系）とQiskitフレームワーク**を用いて実施するための完全なドキュメントが含まれています。

## 文書構成

### 1. 理論書 (Theory)

**ファイル**: `qubit_quantum_dynamics_molecular_triplet_states_theory.md`

**内容**:
- 3準位分子系のQubit表現（2-qubitエンコーディング）
- N分子系の状態空間
- ハミルトニアンのQubit表現（Pauli演算子展開）
- Qiskitゲートによる実装理論
- 鈴木トロッター分解
- 観測量の計算
- 理論的正当性と制約

**特徴**:
- すべての数式を省略無しに展開
- 物理的部分空間の保存を厳密に証明
- Qiskitの標準ゲートのみで実装

**行数**: 1,079行

### 2. 仕様書 (Specification)

**ファイル**: `qubit_implementation_specification.md`

**内容**:
- システム要件とパラメータ仕様
- Qiskitゲートカタログ（20種類以上の完全な定義）
- 状態エンコーディング仕様
- ハミルトニアン項の実装仕様（詳細なゲート列）
- 鈴木トロッター回路仕様
- 観測量計算仕様
- エラーハンドリングと検証
- 性能仕様とベンチマーク

**特徴**:
- 実装レベルの詳細仕様
- 各ゲートの行列表現とコード例
- ゲート数の見積もり
- 計算量とメモリ使用量の分析

**行数**: 1,409行

### 3. 設計書 (Detailed Design)

**ファイル**: `qubit_detailed_design.md`

**内容**:
- システムアーキテクチャ
- クラス設計（6つの主要クラス）
- 完全なゲート分解アルゴリズム
- 実装アルゴリズム（疑似コード付き）
- 完全なPython実装コード（約500行）
- 収束性とエラー解析
- 実行例とチュートリアル

**特徴**:
- 即座に実行可能なPythonコード
- アルゴリズムの詳細な説明
- 期待される出力例
- 収束性の数値検証手法

**行数**: 1,447行

## 読み方

### 初学者向け

1. **理論書**から始める
   - 基本的な概念とQubitエンコーディングを理解
   - ハミルトニアンの表現を学ぶ

2. **仕様書**で詳細を確認
   - 使用するゲートの具体的な定義
   - 実装の細部を理解

3. **設計書**で実装を学ぶ
   - 完全なコード例を確認
   - 実際に実行してみる

### 実装者向け

1. **設計書**から始める
   - クラス設計とアーキテクチャを理解
   - サンプルコードを実行

2. **仕様書**で詳細を確認
   - 必要な仕様を参照
   - ゲート分解の詳細を確認

3. **理論書**で理論的背景を確認
   - 実装の正当性を検証
   - 数式の導出を理解

## 主要な特徴

### ✅ 数学的厳密性

- すべての数式が省略無しに展開されている
- 物理的整合性が理論的に保証されている
- 近似やヒューリスティックを一切使用していない

### ✅ ヒューリスティック排除

以下の手法を**使用禁止**としている：
- `scipy.linalg.expm` による行列指数関数の直接計算
- 近似的なfallback処理
- 非物理的な状態への遷移

### ✅ Qiskitネイティブ

- Qiskitの標準ゲートのみを使用
- 実機（IBMQ等）での実行に対応
- 広く利用可能なハードウェアプラットフォームに対応

### ✅ 完全な実装

- 即座に実行可能なPythonコード
- 単体テストと検証機能
- エラーハンドリングと物理的整合性チェック

## QuditとQubitの比較

| 項目 | Qutrit (MQT-Qudits) | Qubit (本実装) |
|------|-------------------|---------------|
| 1分子の表現 | 1 Qutrit（3次元） | 2 Qubit（4次元、1次元未使用） |
| 4分子系の次元 | $3^4 = 81$ | $2^8 = 256$（物理的: 81） |
| Qubit/Qutrit数 | 4 | 8 |
| ゲート数（1ステップ） | 約55個 | 約430個 |
| 実装の自然性 | 高い | 中程度（エンコーディング必要） |
| ハードウェア可用性 | 実験段階 | 広く利用可能 |
| 開発フレームワーク | MQT-Qudits | Qiskit |

**結論**: Qubit版はQudit版の約8倍のゲート数を必要とするが、現在のハードウェアで広く利用可能。

## 参照文書

本ドキュメントは以下の既存文書を参照・引用している：

- `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb` - Qudit実装のリファレンス
- `tutorials/doc/quantum_dynamics_molecular_triplet_states.md` - 基礎理論
- `tutorials/doc/suzuki_trotter_decomposition_theory.md` - 数値計算理論
- `tutorials/doc/qudit_quantum_algorithm_for_molecular_triplet_dynamics.md` - Qudit完全実装
- `tutorials/doc/mqt_qudits_gates_and_bases_reference.md` - ゲートリファレンス
- `tutorials/doc/n_molecule_triplet_dynamics_basic_gates.md` - N分子系への一般化

## 実装例

```python
from qiskit import QuantumCircuit
import numpy as np

# パラメータ設定
params = PhysicalParameters(
    N_molecules=4,
    E_T=1.5,  # eV
    E_S=3.0,  # eV
    V=0.1,    # eV
    J=0.05,   # eV
    Gamma_fl=0.01  # fs^-1
)

# シミュレータの初期化
simulator = QubitMolecularDynamicsSimulator(params)

# シミュレーション実行
results = simulator.simulate(
    T_total=100.0,  # fs
    N_steps=20,
    initial_state='all_triplet',
    track_dynamics=True
)

# 結果の可視化
simulator.plot_results(results)
```

## システム要件

- Python ≥ 3.8
- Qiskit ≥ 0.40.0
- NumPy ≥ 1.20.0
- Matplotlib ≥ 3.5 (可視化用)
- メモリ: 4分子系で最低 2GB

## 今後の展望

### 実装の最適化

- ゲート分解の効率化
- 回路transpilationの最適化
- 並列化による高速化

### 機能拡張

- 2次元格子系への対応
- 不均一系の実装
- 時間依存ハミルトニアンの導入

### 実機での検証

- IBMQなどの実機での実行
- ノイズモデルの導入
- 誤り訂正の適用

## ライセンスと引用

本ドキュメントは、MQT Quditsプロジェクトの一部として提供されています。

**作成日**: 2025-10-19  
**バージョン**: 1.0.0  
**対象フレームワーク**: Qiskit

---

**問い合わせ**: MQT Qudits プロジェクト  
**リポジトリ**: https://github.com/nobkt/mqt-qudits
