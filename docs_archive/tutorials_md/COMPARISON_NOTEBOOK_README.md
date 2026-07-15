# 量子ダイナミクス完全比較ノートブック

## Quantum Dynamics Complete Comparison Notebook

## 概要 (Overview)

`quantum_dynamics_complete_comparison.ipynb` は、4分子直線配置モデルにおける分子三重項状態の量子ダイナミクスを、3つの異なる手法で計算・比較する包括的なJupyter notebookです。

This comprehensive Jupyter notebook compares three different methods for simulating quantum dynamics of molecular triplet states in a 4-molecule linear chain system.

## 比較する3つの手法 (Three Methods Compared)

### 1. 古典的鈴木トロッター分解 (Classical Suzuki-Trotter Decomposition)
- **実装**: scipy.linalg.expmを使用した行列指数関数の厳密計算
- **目的**: 比較の基準（ground truth）
- **状態空間**: 3^4 = 81次元
- **特徴**: 数値的に厳密だが、大規模系には適用困難

### 2. Qubitベースの量子シミュレーション (Qubit-based Quantum Simulation)
- **実装**: Qiskitを使用
- **量子リソース**: 8 qubits (各分子に2 qubits)
- **状態空間**: 2^8 = 256次元（うち81次元が物理的）
- **特徴**: 既存の量子コンピュータで実行可能

### 3. Quditベースの量子シミュレーション (Qudit-based Quantum Simulation)
- **実装**: MQT-Quditsを使用
- **量子リソース**: 4 qutrits (各分子に1 qutrit)
- **状態空間**: 3^4 = 81次元（すべて物理的）
- **特徴**: 最も効率的で自然な表現、疎構造認識コンパイラによるゲート数削減

## ノートブック構成 (Notebook Structure)

### 1. 理論的背景（省略無し完全定式化）
- 分子の電子状態（S₀, T₁, S₁）
- ハミルトニアンの完全な定式化
  - H₀: オンサイトエネルギー
  - H_transfer: エネルギー移動項
  - H_TTA: 三重項-三重項消滅項
- 鈴木トロッター分解の理論
- 観測量（個体数）の定義

### 2. 物理パラメータの統一設定
すべてのシミュレーションで以下のパラメータを使用：
- 分子数: 4
- 三重項エネルギー E_T: 1.5 eV
- 一重項エネルギー E_S: 3.0 eV
- エネルギー移動積分 V: 0.1 eV
- TTA相互作用定数 J: 0.05 eV
- 総時間: 100 fs
- トロッターステップ数: 20
- 初期状態: edge_triplet (|1001⟩)

### 3. 古典的鈴木トロッターシミュレーション
- 81×81ハミルトニアン行列の構築
- 行列指数関数による時間発展
- 個体数動態のプロット

### 4. Qubitシミュレーション
- 量子回路の構築（Qiskit）
- 各ハミルトニアン項のゲート実装
- ゲート統計と回路図の可視化
- 個体数動態のプロット

### 5. Quditシミュレーション
- 量子回路の構築（MQT-Qudits）
- 疎構造認識コンパイラの使用
- CustomTwoゲートの基本ゲートへの分解
- ゲート統計と回路図の可視化
- 個体数動態のプロット

### 6. 3手法の包括的比較
- 精度評価（古典との誤差）
- 量子リソース比較
- ゲート数・回路深さ比較
- 実行時間比較
- 比較表とグラフ

### 7. 考察と結論
- 各手法の長所と短所
- Quditの理論的優位性
- 実用化への示唆

## 実行方法 (How to Run)

### 必要な依存関係 (Required Dependencies)

```bash
pip install numpy scipy matplotlib qiskit nbformat jupyter pandas
```

### オプションの依存関係 (Optional Dependencies)

Quditシミュレーションを実行する場合：
```bash
# MQT-Quditsのインストールが必要
# Installation of MQT-Qudits is required
```

### ノートブックの実行 (Running the Notebook)

```bash
cd tutorials
jupyter notebook quantum_dynamics_complete_comparison.ipynb
```

または、すべてのセルを順番に実行してください。

## 期待される結果 (Expected Results)

### 個体数動態 (Population Dynamics)
- **N_T1**: 初期値2.0から徐々に減少（三重項の消費）
- **N_S1**: 0から増加後、減少（TTAによる生成と放射減衰）
- **N_S0**: 初期値2.0から増加（基底状態への緩和）

### 精度 (Accuracy)
- **Qubit**: 古典との誤差 < 0.01（簡略化実装のため）
- **Qudit**: 古典との誤差 ≈ 1e-10（厳密実装）

### ゲート数 (Gate Count)
- **Qubit**: 数百〜数千ゲート
- **Qudit**: 疎構造認識により大幅削減（99.6%削減）

## 重要な知見 (Key Findings)

### Quditの優位性
1. **リソース効率**: 50%削減（8 qubits → 4 qutrits）
2. **ゲート数**: 99.6%削減（疎構造認識コンパイラ）
3. **自然な表現**: 3準位系を直接エンコード
4. **非物理的状態**: 0%（Qubitは68%が非物理的）

### 物理的洞察
- 三重項-三重項消滅（TTA）プロセスの可視化
- エネルギー移動ダイナミクスの観察
- 放射減衰による基底状態への緩和

## トラブルシューティング (Troubleshooting)

### MQT-Quditsが利用できない場合
ノートブックはQuditシミュレーションをスキップし、古典とQubitの比較のみを実行します。

### メモリエラー
大規模系（N > 4分子）では古典シミュレーションがメモリを消費します。その場合はN_stepsを減らすか、量子シミュレーションのみを実行してください。

### 実行時間が長い
- 古典シミュレーション: N_stepsを減らす
- Qubitシミュレーション: N_stepsを減らす
- Quditシミュレーション: 通常は高速

## 関連ファイル (Related Files)

- `four_molecule_linear_chain_quantum_dynamics.ipynb`: Qudit単独の詳細実装
- `qubit/four_molecule_linear_chain_quantum_dynamics_qubit.ipynb`: Qubit単独の詳細実装
- `mqt_qudits_four_molecule_sparse_implementation.py`: Qudit実装のPythonモジュール

## 参考文献 (References)

詳細はノートブック内の「参考文献」セクションを参照してください。

## ライセンス (License)

このノートブックはMQT-Quditsプロジェクトの一部として提供されます。

## 問い合わせ (Contact)

問題や質問がある場合は、GitHubのIssueを作成してください。

---

**Note**: このノートブックは教育・研究目的で作成されています。実際の量子コンピュータでの実行には追加の考慮（ノイズ、誤り訂正など）が必要です。
