# Qubit版分子三重項状態量子ダイナミクスチュートリアル

## 📚 概要

このディレクトリには、`tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`（Qudit版）と同等の分子三重項状態量子ダイナミクスシミュレーションを、**Qubit（2準位系）とQiskitフレームワーク**を用いて実施するための資料が含まれています。

---

## 📁 ディレクトリ構成

```
tutorials/qubit/
├── README.md                                                  # このファイル
├── four_molecule_linear_chain_quantum_dynamics_qubit.ipynb   # ✅ Qubit版チュートリアル（ドキュメント専用版）
├── NOTEBOOK_COMPLETION_REPORT.md                              # ✅ ノートブック作成完了報告
├── IMPLEMENTATION_STATUS.md                                   # 実装状況の詳細報告
├── IMPLEMENTATION_GUIDE.md                                    # 実装ガイド
└── ../doc/qubit/                                              # 詳細ドキュメント（理論・仕様・設計）
    ├── README.md
    ├── COMPLETION_SUMMARY.md
    ├── CONTINUATION_PLAN.md
    ├── qubit_quantum_dynamics_molecular_triplet_states_theory.md
    ├── qubit_implementation_specification.md
    └── qubit_detailed_design.md
```

---

## 🎯 目的

### 主な目標

1. **QubitとQuditの比較**: 
   - 同一の物理系に対する異なる実装方式の比較
   - ゲート数、回路深さ、実装の自然性の評価
   
2. **理論的基盤の提供**:
   - 3準位分子系のQubit表現理論
   - 数学的に厳密な実装方法
   - ヒューリスティックを排除した厳密な手法

3. **広く利用可能なハードウェアへの対応**:
   - Qiskitを用いた実装
   - IBMQ、Rigetti、IonQなどのクラウド量子コンピュータで実行可能

---

## 📊 現在の状況

### ✅ 完成しているもの

#### 包括的なドキュメント（4,159行、120,000文字）

`tutorials/doc/qubit/`に以下の完全なドキュメントが用意されています：

1. **理論書** (1,079行)
   - 3準位分子系の2-Qubitエンコーディング
   - ハミルトニアンのPauli演算子表現
   - 鈴木トロッター分解の数学的定式化

2. **仕様書** (1,409行)
   - Qiskitゲートカタログ（20種類以上）
   - ハミルトニアン項の完全なゲート分解
   - 性能仕様とベンチマーク

3. **設計書** (1,447行)
   - 6つの主要クラスの完全設計
   - 実装可能なPythonコード例（約500行）
   - ゲート分解アルゴリズムの詳細

### ✅ 新規追加（2025-10-20）

4. **Jupyter Notebookチュートリアル** ✅ **NEW**
   - [`four_molecule_linear_chain_quantum_dynamics_qubit.ipynb`](./four_molecule_linear_chain_quantum_dynamics_qubit.ipynb)
   - ドキュメント専用版として完成
   - 理論的背景とQubit表現
   - QubitとQuditの詳細比較
   - 実装アーキテクチャと疑似コード
   - 完全実装への道筋

5. **ノートブック作成完了報告** ✅ **NEW**
   - [`NOTEBOOK_COMPLETION_REPORT.md`](./NOTEBOOK_COMPLETION_REPORT.md)
   - 成果物の詳細説明
   - ヒューリスティック排除の方針
   - 実装シナリオと推奨事項

### ⏳ 実装中/計画中

- **完全なPython実装**: 設計完了、コーディング待ち（Qiskit依存の追加が必要）

詳細は [`IMPLEMENTATION_STATUS.md`](./IMPLEMENTATION_STATUS.md) および [`NOTEBOOK_COMPLETION_REPORT.md`](./NOTEBOOK_COMPLETION_REPORT.md) を参照してください。

---

## 🔬 QubitとQuditの比較

### 定量的比較

| 項目 | Qutrit (MQT-Qudits) | Qubit (本実装) | 比率 |
|------|-------------------|---------------|------|
| **1分子の表現** | 1 Qutrit（3次元） | 2 Qubit（4次元） | 2倍 |
| **4分子系の次元** | 3⁴ = 81 | 2⁸ = 256 | 3.16倍 |
| **物理的次元** | 81 | 81 | 同じ |
| **未使用状態** | 0 | 175 (68%) | - |
| **Qubit/Qutrit数** | 4 | 8 | 2倍 |
| **ゲート数/ステップ** | 約55個 | 約430個 | **約8倍** |

### 定性的比較

| 側面 | Qutrit | Qubit |
|------|--------|-------|
| **実装の自然性** | ⭐⭐⭐⭐⭐ 高い | ⭐⭐⭐ 中程度 |
| **ハードウェア可用性** | ⭐⭐ 実験段階 | ⭐⭐⭐⭐⭐ 広く利用可能 |
| **ゲート効率** | ⭐⭐⭐⭐⭐ 高効率 | ⭐⭐ 低効率（8倍） |
| **実機での実行** | ⭐⭐ 限定的 | ⭐⭐⭐⭐⭐ 可能 |
| **数学的厳密性** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

### 結論

- **Qutrit版**: より自然で効率的だが、ハードウェアが限定的
- **Qubit版**: ゲート数は多いが、現在広く利用可能なハードウェアで実行可能

---

## 🚀 使い方（将来の実装）

### 前提条件

```bash
# Python環境
python >= 3.8

# 必要なパッケージ
pip install qiskit>=0.40.0
pip install qiskit-aer>=0.11.0
pip install numpy>=1.20.0
pip install matplotlib>=3.5.0
```

### 基本的な使用例（計画）

```python
from qubit_molecular_dynamics import (
    PhysicalParameters,
    QubitMolecularDynamicsSimulator
)

# パラメータ設定
params = PhysicalParameters(
    N_molecules=4,
    E_T=1.5,  # 三重項エネルギー (eV)
    E_S=3.0,  # 一重項エネルギー (eV)
    V=0.1,    # エネルギー移動積分 (eV)
    J=0.05,   # TTA相互作用定数 (eV)
    Gamma_fl=0.01  # 蛍光放出速度 (fs^-1)
)

# シミュレータの初期化
simulator = QubitMolecularDynamicsSimulator(params)

# シミュレーション実行
results = simulator.simulate(
    T_total=100.0,  # 総時間 (fs)
    N_steps=20,     # ステップ数
    initial_state='all_triplet',  # 初期状態
    track_dynamics=True
)

# 結果の可視化
simulator.plot_results(results)
```

---

## 📖 ドキュメントの読み方

### 初学者向け

1. **まずREADMEを読む** (このファイル)
   - プロジェクト全体の概要を理解

2. **理論書を読む** (`../doc/qubit/qubit_quantum_dynamics_molecular_triplet_states_theory.md`)
   - 3準位分子系のQubit表現
   - ハミルトニアンの数学的表現
   - 鈴木トロッター分解

3. **仕様書を確認** (`../doc/qubit/qubit_implementation_specification.md`)
   - 使用するゲートの定義
   - 実装の詳細仕様

4. **設計書で実装を学ぶ** (`../doc/qubit/qubit_detailed_design.md`)
   - クラス設計
   - コード例

### 実装者向け

1. **設計書から始める** (`../doc/qubit/qubit_detailed_design.md`)
   - クラス設計とアーキテクチャ
   - サンプルコード

2. **仕様書で詳細を確認** (`../doc/qubit/qubit_implementation_specification.md`)
   - ゲートの行列表現
   - ゲート分解の詳細

3. **理論書で検証** (`../doc/qubit/qubit_quantum_dynamics_molecular_triplet_states_theory.md`)
   - 実装の正当性
   - 数式の導出

---

## ⚠️ 重要な実装方針

### ✅ 許可されているもの

- **Qiskitの標準ゲートのみ**
  - 単一qubitゲート: X, Y, Z, H, RX, RY, RZ, S, T, P
  - 2-qubitゲート: CNOT, CZ, SWAP, CRX, CRY, CRZ, RXX, RYY, RZZ
  - 多qubitゲート: Toffoli, Fredkin

- **数学的に厳密な演算**
  - 鈴木トロッター分解
  - ゲートの組み合わせによる正確な演算子実装

### ❌ 禁止されているもの

- **ヒューリスティックな手法**
  - `scipy.linalg.expm` による行列指数関数の直接計算
  - 近似的なfallback処理
  - 非物理的な状態への遷移

---

## 🔗 関連リンク

### ドキュメント

- [実装状況報告](./IMPLEMENTATION_STATUS.md)
- [理論書](../doc/qubit/qubit_quantum_dynamics_molecular_triplet_states_theory.md)
- [仕様書](../doc/qubit/qubit_implementation_specification.md)
- [設計書](../doc/qubit/qubit_detailed_design.md)
- [完了報告](../doc/qubit/COMPLETION_SUMMARY.md)

### 参照実装

- [Qudit版チュートリアル](../four_molecule_linear_chain_quantum_dynamics.ipynb)
- [Qudit版実装コード](../mqt_qudits_four_molecule_implementation.py)

### 理論的背景

- [分子三重項状態の量子ダイナミクス理論](../doc/quantum_dynamics_molecular_triplet_states.md)
- [鈴木トロッター分解理論](../doc/suzuki_trotter_decomposition_theory.md)
- [厳密対角化理論](../doc/exact_diagonalization_theory.md)

---

## 🛠️ 開発ロードマップ

### 短期（1-3ヶ月）

- [ ] Qiskitの依存関係追加
- [ ] 完全なPython実装の作成
- [ ] Jupyter Notebookチュートリアルの作成
- [ ] 単体テストの追加
- [ ] Qudit版との比較ベンチマーク

### 中期（3-6ヶ月）

- [ ] ゲート分解の最適化
- [ ] 回路transpilationの最適化
- [ ] 実機（IBMQ等）での実行テスト
- [ ] ノイズモデルの導入

### 長期（6ヶ月以降）

- [ ] 変分量子アルゴリズム（VQE）の適用
- [ ] 大規模系への拡張
- [ ] 実験データとの比較
- [ ] 他の物理系への応用

---

## 🤝 コントリビューション

### 実装に貢献したい方

1. [実装状況報告](./IMPLEMENTATION_STATUS.md)を確認
2. [設計書](../doc/qubit/qubit_detailed_design.md)を参照
3. issueを立てて実装計画を共有
4. プルリクエストを作成

### ドキュメント改善

- 誤字・脱字の修正
- 説明の追加・改善
- コード例の追加

---

## 📞 サポート

### 質問・相談

- **ドキュメントの内容**: MQT Qudits issue tracker
- **実装の詳細**: 設計書のコメント欄
- **Qubitの理論**: 理論書を参照

### 連絡先

- **プロジェクト**: MQT Qudits
- **リポジトリ**: https://github.com/nobkt/mqt-qudits

---

## 📜 ライセンス

本ドキュメントおよび将来の実装は、MQT Quditsプロジェクトのライセンス（MIT License）に従います。

---

## 📝 変更履歴

| 日付 | バージョン | 変更内容 |
|------|----------|---------|
| 2025-10-19 | 1.0.0 | 初版作成、ドキュメント完成 |
| 2025-10-20 | 1.1.0 | Jupyter Notebookチュートリアル（ドキュメント専用版）追加 |

---

**作成日**: 2025-10-19  
**最終更新**: 2025-10-20  
**ステータス**: ✅ **Jupyter Notebook完成、ドキュメント専用版として提供**  
**バージョン**: 1.1.0
