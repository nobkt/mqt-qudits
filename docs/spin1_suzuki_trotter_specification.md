# スピンS=1量子ダイナミクス - 鈴木トロッター分解実装仕様書

## 1. 概要

本仕様書は、スピンS=1の量子系の時間発展を鈴木トロッター分解を用いてMQT Quditsライブラリで実装するための要件と仕様を定義します。スピンS=1系は3準位（qutrit）量子系として表現され、一般的な量子ビット（2準位系）よりも豊かな量子状態空間を持ちます。

## 2. 目的

- スピンS=1系のハミルトニアン時間発展を効率的にシミュレートする
- 鈴木トロッター分解を用いて、連続的な時間発展を離散的な量子ゲート操作に近似する
- MQT Quditsライブラリの既存機能を活用した実装を提供する
- 物理的に意味のあるスピン相互作用（イジング型、ハイゼンベルク型など）をサポートする

## 3. システム要件

### 3.1 ハードウェア要件
- Python 3.9以上が動作する環境
- メモリ：最小2GB（推奨4GB以上）
- CPU：x86_64またはARM64アーキテクチャ

### 3.2 ソフトウェア要件
- Python 3.9+
- MQT Qudits ライブラリ
- NumPy >= 1.20
- SciPy >= 1.7
- Matplotlib（可視化用、オプション）

### 3.3 依存ライブラリ
```python
mqt.qudits
numpy
scipy
matplotlib  # オプション
```

## 4. 機能要件

### 4.1 スピンS=1演算子

スピンS=1系では、3つのスピン演算子 $S_x$, $S_y$, $S_z$ を実装する必要があります。

#### 4.1.1 スピン演算子の行列表現

**$S_z$ 演算子（対角演算子）：**
```
S_z = ℏ × |1⟩⟨1| - |(-1)⟩⟨(-1)|
```

**$S_x$ 演算子：**
```
S_x = (ℏ/√2) × (|1⟩⟨0| + |0⟩⟨1| + |0⟩⟨-1| + |-1⟩⟨0|)
```

**$S_y$ 演算子：**
```
S_y = (ℏ/√2i) × (|1⟩⟨0| - |0⟩⟨1| + |0⟩⟨-1| - |-1⟩⟨0|)
```

### 4.2 ハミルトニアンの種類

#### 4.2.1 単一スピンハミルトニアン
磁場中のスピン：
```
H_single = -B_z S_z - B_x S_x - B_y S_y
```

#### 4.2.2 二体相互作用ハミルトニアン

**イジング型相互作用：**
```
H_Ising = J S_z^(i) ⊗ S_z^(j)
```

**ハイゼンベルク型相互作用：**
```
H_Heisenberg = J(S_x^(i) ⊗ S_x^(j) + S_y^(i) ⊗ S_y^(j) + S_z^(i) ⊗ S_z^(j))
```

**XXZ型相互作用：**
```
H_XXZ = J_⊥(S_x^(i) ⊗ S_x^(j) + S_y^(i) ⊗ S_y^(j)) + J_z S_z^(i) ⊗ S_z^(j)
```

### 4.3 鈴木トロッター分解

#### 4.3.1 一次鈴木トロッター分解
ハミルトニアン $H = H_A + H_B$ の時間発展演算子を近似：
```
exp(-iHt) ≈ [exp(-iH_A δt) exp(-iH_B δt)]^n
```
ここで、$t = n × δt$、$δt$ はトロッターステップサイズ。

誤差：$O(δt^2)$

#### 4.3.2 二次鈴木トロッター分解（対称分解）
より高精度な近似：
```
exp(-iHt) ≈ [exp(-iH_A δt/2) exp(-iH_B δt) exp(-iH_A δt/2)]^n
```

誤差：$O(δt^3)$

#### 4.3.3 高次鈴木トロッター分解
さらに高精度が必要な場合、4次以上の分解を使用可能：
```
S_4(t) = S_2(p_1 t) S_2(p_2 t) S_2(p_1 t)
```
ここで、$p_1 = 1/(4-4^{1/3})$、$p_2 = 1 - 4p_1$

## 5. 性能要件

### 5.1 精度要件
- トロッターステップサイズに応じた近似誤差の制御
- 一次分解：$O(δt^2)$ の誤差
- 二次分解：$O(δt^3)$ の誤差
- ユニタリ性の保持：$||U^\dagger U - I|| < 10^{-10}$

### 5.2 計算量要件
- N-qudit系の場合、状態ベクトルのサイズ：$3^N$
- 単一ゲート操作：$O(3^N)$
- 二体ゲート操作：$O(3^N)$
- 全体の時間複雑度：$O(n_{steps} × n_{terms} × 3^N)$

## 6. インターフェース仕様

### 6.1 クラス設計

```python
class Spin1Operator:
    """スピンS=1演算子の基底クラス"""
    
    def __init__(self, dimension: int = 3):
        """
        Parameters:
            dimension: qudit次元（スピンS=1の場合は3）
        """
        pass
    
    def Sx(self) -> np.ndarray:
        """S_x演算子の行列を返す"""
        pass
    
    def Sy(self) -> np.ndarray:
        """S_y演算子の行列を返す"""
        pass
    
    def Sz(self) -> np.ndarray:
        """S_z演算子の行列を返す"""
        pass

class Spin1Hamiltonian:
    """スピンS=1系のハミルトニアン"""
    
    def __init__(self, n_spins: int):
        """
        Parameters:
            n_spins: スピンの数
        """
        pass
    
    def add_single_spin_term(self, site: int, Bx: float, By: float, Bz: float):
        """単一スピン項を追加"""
        pass
    
    def add_ising_interaction(self, site_i: int, site_j: int, J: float):
        """イジング相互作用を追加"""
        pass
    
    def add_heisenberg_interaction(self, site_i: int, site_j: int, J: float):
        """ハイゼンベルク相互作用を追加"""
        pass
    
    def to_matrix(self) -> np.ndarray:
        """ハミルトニアンの完全行列表現を返す"""
        pass

class SuzukiTrotterEvolution:
    """鈴木トロッター分解による時間発展"""
    
    def __init__(self, hamiltonian: Spin1Hamiltonian, order: int = 1):
        """
        Parameters:
            hamiltonian: 時間発展させるハミルトニアン
            order: トロッター分解の次数（1, 2, 4をサポート）
        """
        pass
    
    def build_circuit(self, time: float, n_steps: int) -> QuantumCircuit:
        """
        トロッター分解された量子回路を構築
        
        Parameters:
            time: 総時間発展時間
            n_steps: トロッターステップ数
        
        Returns:
            構築された量子回路
        """
        pass
    
    def evolve_state(self, initial_state: np.ndarray, time: float, 
                     n_steps: int) -> np.ndarray:
        """
        初期状態を時間発展させる
        
        Parameters:
            initial_state: 初期量子状態
            time: 総時間発展時間
            n_steps: トロッターステップ数
        
        Returns:
            時間発展後の量子状態
        """
        pass
```

### 6.2 MQT Quditsライブラリとの統合

MQT Quditsライブラリの既存機能を使用：
- `QuantumCircuit`: 量子回路の構築
- `QuantumRegister`: qutritレジスタの定義
- `Rz`, `R`, `Rh`: 単一qudit回転ゲート
- `LS`, `MS`: 多体ゲート（Lorch-Schack, Mølmer-Sørensen）
- `MQTQuditProvider`: シミュレーションバックエンド

## 7. 出力仕様

### 7.1 時間発展の結果
- 各時刻における状態ベクトル
- 物理量の期待値（磁化、エネルギーなど）
- 忠実度（理論値との比較）

### 7.2 可視化出力
- 時間発展のアニメーション
- 物理量の時間依存性グラフ
- 量子状態のブロッホ球表現（単一スピンの場合）

## 8. エラーハンドリング

### 8.1 入力検証
- スピン数が正の整数であることを確認
- トロッター次数が実装済み（1, 2, 4）であることを確認
- 時間パラメータが正であることを確認
- ハミルトニアン項が有効なサイトインデックスを持つことを確認

### 8.2 計算時のエラー
- メモリ不足エラーのハンドリング
- 数値不安定性の検出（ユニタリ性の破れ）
- 行列次元の不一致の検出

## 9. テストとバリデーション

### 9.1 単体テスト
- スピン演算子の交換関係の確認：$[S_i, S_j] = i\epsilon_{ijk} S_k$
- ハミルトニアンのエルミート性の確認
- 時間発展演算子のユニタリ性の確認

### 9.2 統合テスト
- 既知の解析解との比較（単一スピン、二スピン系）
- トロッター次数と精度の関係の検証
- ステップ数と計算時間のスケーリング

### 9.3 物理的妥当性の検証
- エネルギー保存（時間非依存ハミルトニアンの場合）
- 全スピンの保存（対称ハミルトニアンの場合）
- カオス的挙動の再現（特定のパラメータ領域）

## 10. 拡張性

### 10.1 将来の拡張
- 任意のスピン量子数Sへの一般化
- 時間依存ハミルトニアンのサポート
- 開放系ダイナミクス（マスター方程式）
- GPU加速による大規模系のシミュレーション

### 10.2 カスタマイズポイント
- ユーザー定義ハミルトニアン項の追加
- カスタム観測量の定義
- 異なる初期状態の準備方法

## 11. 参考文献

1. Suzuki, M. (1976). "Generalized Trotter's formula and systematic approximants of exponential operators and inner derivations with applications to many-body problems." Communications in Mathematical Physics, 51(2), 183-190.

2. Lloyd, S. (1996). "Universal quantum simulators." Science, 273(5278), 1073-1078.

3. Hatano, N., & Suzuki, M. (2005). "Finding exponential product formulas of higher orders." In Quantum Annealing and Other Optimization Methods (pp. 37-68). Springer, Berlin, Heidelberg.

4. Childs, A. M., et al. (2021). "Theory of Trotter error with commutator scaling." Physical Review X, 11(1), 011020.

## 12. 変更履歴

| 版 | 日付 | 変更内容 | 作成者 |
|----|------|---------|--------|
| 1.0 | 2025-10-14 | 初版作成 | MQT Qudits Team |

