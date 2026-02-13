# TTA-UC現象GKSL-Lindblad量子ダイナミクス詳細実装計画書

## 文書情報

- **文書名**: TTA-UC現象GKSL-Lindblad量子ダイナミクス詳細実装計画書
- **バージョン**: v1.0.0
- **作成日**: 2026年2月13日
- **対象システム**: MQT-Qudits GKSL-Lindblad量子ダイナミクスシミュレータ
- **参照文書**:
  - `tutorials/doc/GKSL/TTA-UC現象のGKSL-Lindblad量子ダイナミクス完全理論書.md`
  - `tutorials/doc/GKSL/TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md`
  - `tutorials/doc/GKSL/TTA-UC現象GKSL-Lindblad量子ダイナミクス詳細設計書.md`

## 文書の目的と範囲

### 目的

本文書は、TTA-UC現象のGKSL-Lindblad量子ダイナミクスシミュレーションシステムを完璧に実装するための詳細実装計画書である。本文書に従って実装を進めることで、以下の3つの参照文書の要求事項を完全に満たす実装が完成する。

### 範囲

- **対象**: 6つの実装シナリオ全て（Classical/Qubit/Qudit × Boson有/無）
- **実装ファイル数**: 約15ファイル（シミュレータ本体、補助モジュール、テスト、ノートブック）
- **実装規模**: 約2,500行のPythonコード + 1つの統合ノートブック
- **実装期間**: 1-3 PR（本計画書作成: PR#142、実装フェーズ1-2: 後続PR）

### 基本方針

1. **ごまかし無しの原則**: ヒューリスティックな処理やfallbackを一切使用しない
2. **真実ベースの実装**: 理論に基づく厳密な数値計算のみ
3. **完全な検証**: 各実装ステップで数学的・物理的整合性を検証
4. **後方互換性**: 既存ノートブックの出力形式を維持
5. **段階的実装**: 優先度順にシナリオを実装し、各段階で検証

## 実装の全体構成

### ファイル構成

```
mqt-qudits/
├── tutorials/
│   ├── gksl_physical_parameters.py          # 物理パラメータクラス（共通）
│   ├── gksl_math_utils.py                   # 数学的基盤関数（共通）
│   ├── stinespring_utils.py                 # Stinespring dilation（共通）
│   ├── classical_gksl_simulator.py          # シナリオ1: 古典・ボソン無し
│   ├── classical_gksl_boson_simulator.py    # シナリオ2: 古典・ボソン有り
│   ├── qubit_gksl_simulator.py              # シナリオ3: Qubit・ボソン無し
│   ├── qubit_gksl_boson_simulator.py        # シナリオ4: Qubit・ボソン有り
│   ├── qudit_gksl_simulator.py              # シナリオ5: Qudit・ボソン無し
│   ├── qudit_gksl_boson_simulator.py        # シナリオ6: Qudit・ボソン有り
│   ├── gksl_validation.py                   # 検証関数（共通）
│   ├── gksl_visualization.py                # 可視化関数（共通）
│   ├── test_gksl_simulators.py              # 包括的テストスイート
│   └── quantum_dynamics_gksl_comparison.ipynb  # 統合比較ノートブック
├── developing/
│   ├── TTA-UC現象GKSL-Lindblad量子ダイナミクス詳細実装計画書.md  # 本文書
│   └── TTA-UC現象GKSL-Lindblad実装作業履歴_PR142.md  # 作業履歴
└── ...
```

### 実装の優先順位

| 優先度 | シナリオ | 実装ファイル | 理由 |
|--------|---------|-------------|------|
| **1** | Classical NB | classical_gksl_simulator.py | 基準実装・検証用 |
| **2** | Qudit NB | qudit_gksl_simulator.py | MQT-Quditsの中核機能 |
| **3** | Qubit NB | qubit_gksl_simulator.py | 比較・検証用 |
| **4** | Classical B | classical_gksl_boson_simulator.py | 拡張モデル |
| **5** | Qudit B | qudit_gksl_boson_simulator.py | 完全実装 |
| **6** | Qubit B | qubit_gksl_boson_simulator.py | 完全実装 |

### 依存関係

```
共通モジュール（最初に実装）:
  gksl_physical_parameters.py
  gksl_math_utils.py
  stinespring_utils.py
  gksl_validation.py
     ↓
各シナリオのシミュレータ（優先度順に実装）
     ↓
gksl_visualization.py（可視化）
     ↓
test_gksl_simulators.py（テスト）
     ↓
quantum_dynamics_gksl_comparison.ipynb（統合）
```

---

## 第1部: 共通モジュールの実装

### 1.1 gksl_physical_parameters.py

#### 1.1.1 目的

全シナリオで共通に使用する物理パラメータを定義・検証するクラスを提供する。

#### 1.1.2 クラス設計

**GKSLPhysicalParameters**

```python
class GKSLPhysicalParameters:
    """
    GKSL-Lindblad量子ダイナミクスの物理パラメータ
    
    設計原則:
    - 全パラメータをfloatで保持（SI単位系: eV, eV/ℏ）
    - validate()メソッドで7つの整合性チェック
    - 不正なパラメータでは例外を発生（fallback禁止）
    """
    
    def __init__(self,
                 # エネルギーパラメータ (eV)
                 E_T: float = 1.5,
                 E_S: float = 3.0,
                 V: float = 0.1,
                 
                 # 散逸レート (eV/ℏ)
                 gamma_TTA: float = 0.05,
                 Gamma_fl: float = 0.01,
                 Gamma_ph: float = 1e-6,
                 k_IC: float = 0.005,
                 k_ISC_ST: float = 0.003,
                 k_ISC_TS: float = 1e-5,
                 
                 # システムパラメータ
                 N_molecules: int = 4,
                 d: int = 3,
                 
                 # ボソンパラメータ（オプション）
                 with_boson: bool = False,
                 n_max: int = 2,
                 omega_ph: float = 0.15,
                 g_eph: float = 0.02):
        """
        パラメータの初期化
        
        Parameters:
        -----------
        E_T : float
            三重項状態エネルギー (eV)
        E_S : float
            一重項状態エネルギー (eV)
        V : float
            エネルギー移動積分 (eV)
        gamma_TTA : float
            TTA散逸レート (eV/ℏ)
        Gamma_fl : float
            蛍光散逸レート (eV/ℏ)
        Gamma_ph : float
            燐光散逸レート (eV/ℏ)
        k_IC : float
            内部変換レート (eV/ℏ)
        k_ISC_ST : float
            ISC S→T レート (eV/ℏ)
        k_ISC_TS : float
            ISC T→S レート (eV/ℏ)
        N_molecules : int
            分子数（デフォルト: 4）
        d : int
            局所次元（デフォルト: 3 = qutrit）
        with_boson : bool
            ボソン相互作用を含むか
        n_max : int
            フォノンカットオフ（ボソン有りの場合）
        omega_ph : float
            フォノン周波数 (eV/ℏ)
        g_eph : float
            電子-フォノン結合定数 (eV)
        """
        self.E_T = E_T
        self.E_S = E_S
        self.V = V
        
        self.gamma_TTA = gamma_TTA
        self.Gamma_fl = Gamma_fl
        self.Gamma_ph = Gamma_ph
        self.k_IC = k_IC
        self.k_ISC_ST = k_ISC_ST
        self.k_ISC_TS = k_ISC_TS
        
        self.N_molecules = N_molecules
        self.d = d
        
        self.with_boson = with_boson
        self.n_max = n_max
        self.omega_ph = omega_ph
        self.g_eph = g_eph
        
        # 検証を実行
        self.validate()
    
    def validate(self) -> None:
        """
        物理パラメータの整合性を検証
        
        検証項目（7項目）:
        1. エネルギー準位の正値性
        2. エネルギー関係式の確認（2*E_T ≈ E_S）
        3. 散逸レートの正値性
        4. 時間スケール階層の検証
        5. システムパラメータの妥当性
        6. ボソンパラメータの妥当性（with_boson=Trueの場合）
        7. GKSL理論の適用条件（弱結合近似）
        
        Raises:
        -------
        ValueError
            パラメータが物理的に不正な場合
        """
        # 1. エネルギー準位の正値性
        if self.E_T <= 0:
            raise ValueError(f"三重項エネルギーE_Tは正でなければなりません: E_T={self.E_T}")
        if self.E_S <= 0:
            raise ValueError(f"一重項エネルギーE_Sは正でなければなりません: E_S={self.E_S}")
        if self.E_S <= self.E_T:
            raise ValueError(f"一重項エネルギーは三重項より大きくなければなりません: E_S={self.E_S}, E_T={self.E_T}")
        
        # 2. エネルギー関係式の確認（2*E_T ≈ E_S, 許容誤差10%）
        energy_ratio = abs(2 * self.E_T - self.E_S) / self.E_S
        if energy_ratio > 0.1:
            raise ValueError(f"エネルギー関係式 2*E_T ≈ E_S が満たされません: 2*E_T={2*self.E_T}, E_S={self.E_S}, 誤差={energy_ratio*100:.1f}%")
        
        # 3. 散逸レートの正値性
        rates = {
            'gamma_TTA': self.gamma_TTA,
            'Gamma_fl': self.Gamma_fl,
            'Gamma_ph': self.Gamma_ph,
            'k_IC': self.k_IC,
            'k_ISC_ST': self.k_ISC_ST,
            'k_ISC_TS': self.k_ISC_TS
        }
        for name, rate in rates.items():
            if rate < 0:
                raise ValueError(f"散逸レート{name}は非負でなければなりません: {name}={rate}")
        
        # 4. 時間スケール階層の検証
        # TTA > 蛍光 > IC > ISC(S→T) >> 燐光, ISC(T→S)
        if not (self.gamma_TTA >= self.Gamma_fl):
            raise ValueError(f"時間スケール階層違反: gamma_TTA >= Gamma_fl が必要 (gamma_TTA={self.gamma_TTA}, Gamma_fl={self.Gamma_fl})")
        if not (self.Gamma_fl >= self.k_IC):
            raise ValueError(f"時間スケール階層違反: Gamma_fl >= k_IC が必要 (Gamma_fl={self.Gamma_fl}, k_IC={self.k_IC})")
        if not (self.k_IC >= self.k_ISC_ST):
            raise ValueError(f"時間スケール階層違反: k_IC >= k_ISC_ST が必要 (k_IC={self.k_IC}, k_ISC_ST={self.k_ISC_ST})")
        
        # 5. システムパラメータの妥当性
        if self.N_molecules < 2:
            raise ValueError(f"分子数N_moleculesは2以上でなければなりません: N_molecules={self.N_molecules}")
        if self.d != 3:
            raise ValueError(f"局所次元dは3（qutrit）でなければなりません: d={self.d}")
        
        # 6. ボソンパラメータの妥当性
        if self.with_boson:
            if self.n_max < 1:
                raise ValueError(f"フォノンカットオフn_maxは1以上でなければなりません: n_max={self.n_max}")
            if self.omega_ph <= 0:
                raise ValueError(f"フォノン周波数omega_phは正でなければなりません: omega_ph={self.omega_ph}")
            if self.g_eph < 0:
                raise ValueError(f"電子-フォノン結合g_ephは非負でなければなりません: g_eph={self.g_eph}")
        
        # 7. GKSL理論の適用条件（弱結合近似）
        # 散逸レート << システムエネルギー
        max_dissipation = max(rates.values())
        min_energy = min(self.E_T, self.V)
        if max_dissipation > 0.5 * min_energy:
            raise ValueError(f"弱結合近似違反: 最大散逸レート({max_dissipation})がシステムエネルギー({min_energy})に対して大きすぎます")
    
    def get_hilbert_space_dim(self) -> int:
        """
        ヒルベルト空間の次元を返す
        
        Returns:
        --------
        dim : int
            ボソン無し: d^N_molecules = 3^4 = 81
            ボソン有り: d^N_molecules * (n_max+1)^N_molecules = 81 * 3^4 = 6561 (n_max=2)
        """
        electronic_dim = self.d ** self.N_molecules
        if self.with_boson:
            phonon_dim = (self.n_max + 1) ** self.N_molecules
            return electronic_dim * phonon_dim
        else:
            return electronic_dim
    
    def get_time_unit(self) -> float:
        """
        特徴的な時間単位を返す（ℏ/eV）
        
        Returns:
        --------
        time_unit : float
            ℏ/eV = 0.6582 fs
        """
        hbar_eV_fs = 0.6582  # ℏ in eV·fs
        return hbar_eV_fs
    
    def to_dict(self) -> dict:
        """
        パラメータを辞書形式で返す（シリアライズ用）
        
        Returns:
        --------
        params : dict
            全パラメータの辞書
        """
        return {
            'E_T': self.E_T,
            'E_S': self.E_S,
            'V': self.V,
            'gamma_TTA': self.gamma_TTA,
            'Gamma_fl': self.Gamma_fl,
            'Gamma_ph': self.Gamma_ph,
            'k_IC': self.k_IC,
            'k_ISC_ST': self.k_ISC_ST,
            'k_ISC_TS': self.k_ISC_TS,
            'N_molecules': self.N_molecules,
            'd': self.d,
            'with_boson': self.with_boson,
            'n_max': self.n_max,
            'omega_ph': self.omega_ph,
            'g_eph': self.g_eph
        }
    
    def __repr__(self) -> str:
        """文字列表現"""
        boson_str = "with boson" if self.with_boson else "without boson"
        return f"GKSLPhysicalParameters(N={self.N_molecules}, d={self.d}, {boson_str})"
```

#### 1.1.3 実装手順

1. **ファイル作成**: `tutorials/gksl_physical_parameters.py`を作成
2. **import文追加**:
   ```python
   from typing import Dict
   import numpy as np
   ```
3. **クラス実装**: 上記の完全なクラス定義をコピー
4. **docstring確認**: 全メソッドのdocstringが完全であることを確認
5. **コメント追加**: 複雑なロジック（validate内の階層チェック等）にコメント
6. **単体テスト作成**: `test_gksl_simulators.py`に以下のテストを追加
   - 正常系: デフォルトパラメータでの初期化
   - 異常系1: E_Sがnegativeな場合
   - 異常系2: エネルギー関係式違反
   - 異常系3: 時間スケール階層違反
   - 異常系4: 弱結合近似違反
   - ボソン有り: with_boson=Trueでの初期化と検証

#### 1.1.4 検証基準

- [ ] 全7項目のvalidateが正しく機能する
- [ ] 不正なパラメータで適切な例外が発生する
- [ ] デフォルトパラメータで初期化できる
- [ ] get_hilbert_space_dim()が正しい次元を返す（81 or 6561）
- [ ] to_dict()で全パラメータがシリアライズできる

### 1.2 gksl_math_utils.py

#### 1.2.1 目的

GKSL-Lindblad方程式の数学的基盤となる関数群を提供する:
- ハミルトニアン構築
- Lindblad演算子構築
- 密度行列操作
- 超演算子構築

#### 1.2.2 関数設計

**build_onsite_hamiltonian(params: GKSLPhysicalParameters) -> np.ndarray**

```python
def build_onsite_hamiltonian(params: GKSLPhysicalParameters) -> np.ndarray:
    """
    オンサイトハミルトニアン H_0 を構築
    
    数式:
        H_0 = Σ_i (E_T |1⟩⟨1| + E_S |2⟩⟨2|)_i
    
    Parameters:
    -----------
    params : GKSLPhysicalParameters
        物理パラメータ
    
    Returns:
    --------
    H_0 : np.ndarray, shape (d^N, d^N)
        オンサイトハミルトニアン（Hermitian行列）
    
    実装詳細:
    ---------
    1. 単一分子の局所演算子を構築:
       h_local = E_T * |1⟩⟨1| + E_S * |2⟩⟨2|
       = diag(0, E_T, E_S)
    
    2. 各分子サイトiに対してテンソル積で拡張:
       H_i = I⊗...⊗I⊗h_local⊗I⊗...⊗I (i番目の位置)
    
    3. 全サイトで合計:
       H_0 = Σ_i H_i
    
    検証:
    -----
    - Hermitian性: H_0 == H_0^†
    - 対角要素のみ非ゼロ
    - 最小固有値 = 0 (基底状態 |0000⟩)
    - 最大固有値 = 4*E_S (励起状態 |2222⟩)
    """
    N = params.N_molecules
    d = params.d
    E_T = params.E_T
    E_S = params.E_S
    
    # 単一分子の局所ハミルトニアン
    h_local = np.diag([0.0, E_T, E_S])
    
    # 恒等演算子
    I = np.eye(d)
    
    # 全系のハミルトニアンを初期化
    H_0 = np.zeros((d**N, d**N), dtype=np.complex128)
    
    # 各サイトに対してテンソル積で構築
    for i in range(N):
        # i番目のサイトにh_localを配置
        op_list = [I] * N
        op_list[i] = h_local
        
        # テンソル積を計算
        H_i = op_list[0]
        for j in range(1, N):
            H_i = np.kron(H_i, op_list[j])
        
        H_0 += H_i
    
    # Hermitian性の検証
    if not np.allclose(H_0, H_0.conj().T):
        raise ValueError("構築されたH_0がHermitianではありません")
    
    return H_0
```

**build_transfer_hamiltonian(params: GKSLPhysicalParameters) -> np.ndarray**

```python
def build_transfer_hamiltonian(params: GKSLPhysicalParameters) -> np.ndarray:
    """
    エネルギー移動ハミルトニアン H_transfer を構築
    
    数式:
        H_transfer = Σ_{⟨i,j⟩} V * (|0⟩_i⟨1| ⊗ |1⟩_j⟨0| + h.c.)
    
    Parameters:
    -----------
    params : GKSLPhysicalParameters
        物理パラメータ
    
    Returns:
    --------
    H_transfer : np.ndarray, shape (d^N, d^N)
        エネルギー移動ハミルトニアン（Hermitian行列）
    
    実装詳細:
    ---------
    1. 三重項-三重項のホッピングのみ考慮
       （一重項は局在化している）
    
    2. 最近接ペア: (0,1), (1,2), (2,3) の3ペア
    
    3. 各ペアに対して:
       T_ij = |01⟩⟨10| + |10⟩⟨01|
       をテンソル積で構築
    
    4. 振幅Vを掛ける
    
    検証:
    -----
    - Hermitian性: H_transfer == H_transfer^†
    - 非対角要素のみ非ゼロ（エネルギー保存遷移）
    - |01⟩ ↔ |10⟩ の遷移振幅が V
    """
    N = params.N_molecules
    d = params.d
    V = params.V
    
    # 単一サイトの演算子
    ket_0 = np.array([1, 0, 0])
    ket_1 = np.array([0, 1, 0])
    
    # |0⟩⟨1| と |1⟩⟨0|
    T_forward = np.outer(ket_0, ket_1)   # |0⟩⟨1|
    T_backward = np.outer(ket_1, ket_0)  # |1⟩⟨0|
    
    I = np.eye(d)
    
    # 全系のハミルトニアンを初期化
    H_transfer = np.zeros((d**N, d**N), dtype=np.complex128)
    
    # 最近接ペア
    for i in range(N - 1):
        j = i + 1
        
        # 前方ホッピング: |0⟩_i⟨1| ⊗ |1⟩_j⟨0|
        op_list = [I] * N
        op_list[i] = T_forward
        op_list[j] = T_backward
        
        hop_forward = op_list[0]
        for k in range(1, N):
            hop_forward = np.kron(hop_forward, op_list[k])
        
        # 後方ホッピング: |1⟩_i⟨0| ⊗ |0⟩_j⟨1| (エルミート共役)
        op_list = [I] * N
        op_list[i] = T_backward
        op_list[j] = T_forward
        
        hop_backward = op_list[0]
        for k in range(1, N):
            hop_backward = np.kron(hop_backward, op_list[k])
        
        # ペアの寄与を加算
        H_transfer += V * (hop_forward + hop_backward)
    
    # Hermitian性の検証
    if not np.allclose(H_transfer, H_transfer.conj().T):
        raise ValueError("構築されたH_transferがHermitianではありません")
    
    return H_transfer
```

**build_lindblad_operators(params: GKSLPhysicalParameters) -> list**

```python
def build_lindblad_operators(params: GKSLPhysicalParameters) -> list:
    """
    26個のLindblad演算子を構築
    
    Parameters:
    -----------
    params : GKSLPhysicalParameters
        物理パラメータ
    
    Returns:
    --------
    lindblad_ops : list of tuples
        各要素は (L_alpha, gamma_alpha) のタプル
        L_alpha: np.ndarray, shape (d^N, d^N)
        gamma_alpha: float
        合計26個
    
    演算子の種類:
    -------------
    1. TTA: 6個
       - 各ペア(i,j)に対して2チャネル
       - L_TTA,1^(ij) = |2⟩_i⟨1| ⊗ |0⟩_j⟨1|
       - L_TTA,2^(ij) = |0⟩_i⟨1| ⊗ |2⟩_j⟨1|
       - gamma = gamma_TTA / 2
    
    2. 蛍光: 4個
       - L_fl^(i) = |0⟩_i⟨2|
       - gamma = Gamma_fl
    
    3. 燐光: 4個
       - L_ph^(i) = |0⟩_i⟨1|
       - gamma = Gamma_ph
    
    4. 内部変換: 4個
       - L_IC^(i) = |0⟩_i⟨2|
       - gamma = k_IC
    
    5. ISC S→T: 4個
       - L_ISC_ST^(i) = |1⟩_i⟨2|
       - gamma = k_ISC_ST
    
    6. ISC T→S: 4個
       - L_ISC_TS^(i) = |0⟩_i⟨1|
       - gamma = k_ISC_TS
    
    実装詳細:
    ---------
    各演算子をテンソル積で構築し、リストに追加
    """
    N = params.N_molecules
    d = params.d
    
    # 単一サイトの基底ベクトル
    ket_0 = np.array([1, 0, 0])
    ket_1 = np.array([0, 1, 0])
    ket_2 = np.array([0, 0, 1])
    
    I = np.eye(d)
    
    lindblad_ops = []
    
    # 1. TTA演算子（6個）
    for i in range(N - 1):
        j = i + 1
        
        # チャネル1: |2⟩_i⟨1| ⊗ |0⟩_j⟨1|
        op_list = [I] * N
        op_list[i] = np.outer(ket_2, ket_1)  # |2⟩⟨1|
        op_list[j] = np.outer(ket_0, ket_1)  # |0⟩⟨1|
        
        L_TTA_1 = op_list[0]
        for k in range(1, N):
            L_TTA_1 = np.kron(L_TTA_1, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.gamma_TTA / 2) * L_TTA_1, params.gamma_TTA / 2))
        
        # チャネル2: |0⟩_i⟨1| ⊗ |2⟩_j⟨1|
        op_list = [I] * N
        op_list[i] = np.outer(ket_0, ket_1)  # |0⟩⟨1|
        op_list[j] = np.outer(ket_2, ket_1)  # |2⟩⟨1|
        
        L_TTA_2 = op_list[0]
        for k in range(1, N):
            L_TTA_2 = np.kron(L_TTA_2, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.gamma_TTA / 2) * L_TTA_2, params.gamma_TTA / 2))
    
    # 2. 蛍光（4個）
    for i in range(N):
        op_list = [I] * N
        op_list[i] = np.outer(ket_0, ket_2)  # |0⟩⟨2|
        
        L_fl = op_list[0]
        for k in range(1, N):
            L_fl = np.kron(L_fl, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.Gamma_fl) * L_fl, params.Gamma_fl))
    
    # 3. 燐光（4個）
    for i in range(N):
        op_list = [I] * N
        op_list[i] = np.outer(ket_0, ket_1)  # |0⟩⟨1|
        
        L_ph = op_list[0]
        for k in range(1, N):
            L_ph = np.kron(L_ph, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.Gamma_ph) * L_ph, params.Gamma_ph))
    
    # 4. 内部変換（4個）
    for i in range(N):
        op_list = [I] * N
        op_list[i] = np.outer(ket_0, ket_2)  # |0⟩⟨2|
        
        L_IC = op_list[0]
        for k in range(1, N):
            L_IC = np.kron(L_IC, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.k_IC) * L_IC, params.k_IC))
    
    # 5. ISC S→T（4個）
    for i in range(N):
        op_list = [I] * N
        op_list[i] = np.outer(ket_1, ket_2)  # |1⟩⟨2|
        
        L_ISC_ST = op_list[0]
        for k in range(1, N):
            L_ISC_ST = np.kron(L_ISC_ST, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.k_ISC_ST) * L_ISC_ST, params.k_ISC_ST))
    
    # 6. ISC T→S（4個）
    for i in range(N):
        op_list = [I] * N
        op_list[i] = np.outer(ket_0, ket_1)  # |0⟩⟨1|
        
        L_ISC_TS = op_list[0]
        for k in range(1, N):
            L_ISC_TS = np.kron(L_ISC_TS, op_list[k])
        
        lindblad_ops.append((np.sqrt(params.k_ISC_TS) * L_ISC_TS, params.k_ISC_TS))
    
    # 総数の検証
    if len(lindblad_ops) != 26:
        raise ValueError(f"Lindblad演算子の数が不正です: {len(lindblad_ops)} (期待値: 26)")
    
    return lindblad_ops
```

**その他の重要関数**

```python
def vectorize_density_matrix(rho: np.ndarray) -> np.ndarray:
    """
    密度行列をベクトルに変換（列優先）
    
    vec(ρ) : (d^N × d^N) -> (d^{2N},)
    """
    return rho.flatten(order='F')  # Fortran順（列優先）

def unvectorize_density_matrix(vec: np.ndarray, dim: int) -> np.ndarray:
    """
    ベクトルを密度行列に復元
    """
    return vec.reshape((dim, dim), order='F')

def compute_von_neumann_entropy(rho: np.ndarray, tol: float = 1e-12) -> float:
    """
    von Neumannエントロピーを計算
    
    S = -Tr[ρ ln ρ]
    """
    # 固有値を計算
    eigenvalues = np.linalg.eigvalsh(rho)
    
    # 正の固有値のみ考慮（数値誤差で負になる可能性を排除）
    eigenvalues = eigenvalues[eigenvalues > tol]
    
    # エントロピー計算
    entropy = -np.sum(eigenvalues * np.log(eigenvalues))
    
    return entropy

def compute_purity(rho: np.ndarray) -> float:
    """
    純度を計算
    
    P = Tr[ρ^2]
    """
    return np.real(np.trace(rho @ rho))
```

#### 1.2.3 実装手順

1. **ファイル作成**: `tutorials/gksl_math_utils.py`
2. **import文**:
   ```python
   import numpy as np
   from typing import List, Tuple
   from gksl_physical_parameters import GKSLPhysicalParameters
   ```
3. **各関数の実装**: 上記の関数を順次実装
4. **docstring完備**: 全関数に詳細なdocstringを追加
5. **単体テスト**: 各関数の正しさを検証
   - H_0の固有値が正しいか
   - H_transferがHermitianか
   - Lindblad演算子が26個か
   - エントロピーが非負か

#### 1.2.4 検証基準

- [ ] H_0とH_transferがHermitian行列である
- [ ] Lindblad演算子が正確に26個構築される
- [ ] ベクトル化/逆ベクトル化が可逆である
- [ ] エントロピー計算が非負値を返す
- [ ] 純度計算が[0,1]の範囲内である

### 1.3 stinespring_utils.py

#### 1.3.1 目的

Stinespring dilationを使用してLindblad演算子を量子回路（ユニタリ演算）に変換する。

#### 1.3.2 理論背景

Lindblad散逸項 `D[L][ρ] = Lρ L† - ½(L†L ρ + ρ L†L)` を、補助qubit/quditを用いたユニタリ演算で実現する:

```
E[ρ_S] = Tr_E[U (ρ_S ⊗ |0⟩_E⟨0|) U†]
```

ここで:
- `ρ_S`: システム密度行列
- `|0⟩_E`: 補助系（環境）の初期状態
- `U`: システム+環境の合成ユニタリ

#### 1.3.3 関数設計

**stinespring_unitary_from_lindblad(L: np.ndarray, dt: float) -> np.ndarray**

```python
def stinespring_unitary_from_lindblad(L: np.ndarray, dt: float) -> np.ndarray:
    """
    Lindblad演算子LからStinespringユニタリを構築
    
    Parameters:
    -----------
    L : np.ndarray, shape (d^N, d^N)
        Lindblad演算子（係数√γを含む）
    dt : float
        時間ステップ
    
    Returns:
    --------
    U : np.ndarray, shape (2*d^N, 2*d^N)
        Stinespringユニタリ（システム+補助qubit）
    
    理論:
    -----
    1. 生成子Gを構築:
       G = [ 0    L  ]
           [ L†   0  ]
    
    2. ユニタリを計算:
       U = exp(-i θ G)
       where θ = √(γ dt)
    
    3. 近似（小さいdt）:
       U ≈ I - i θ G + O(θ^2)
    
    実装詳細:
    ---------
    - 行列指数関数はscipy.linalg.expmを使用
    - ユニタリ性の検証: U U† = I
    """
    d_sys = L.shape[0]
    
    # θ = √(γ dt) を計算
    # Lには既に√γが含まれているので、θ = ||L|| * √dt
    theta = np.sqrt(dt)
    
    # 生成子G（反Hermitian行列）
    G = np.block([
        [np.zeros_like(L), L],
        [L.conj().T, np.zeros_like(L)]
    ])
    
    # ユニタリU = exp(-i θ G)
    from scipy.linalg import expm
    U = expm(-1j * theta * G)
    
    # ユニタリ性の検証
    I_full = np.eye(2 * d_sys)
    if not np.allclose(U @ U.conj().T, I_full, atol=1e-10):
        raise ValueError("構築されたStinespringユニタリがユニタリ性を満たしません")
    
    return U

def apply_stinespring_to_density_matrix(rho: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    Stinespringユニタリを密度行列に適用
    
    Parameters:
    -----------
    rho : np.ndarray, shape (d^N, d^N)
        システム密度行列
    U : np.ndarray, shape (2*d^N, 2*d^N)
        Stinespringユニタリ
    
    Returns:
    --------
    rho_out : np.ndarray, shape (d^N, d^N)
        演算後のシステム密度行列
    
    手順:
    -----
    1. 拡張密度行列を構築: ρ_ext = ρ_S ⊗ |0⟩_E⟨0|
    2. ユニタリを適用: ρ' = U ρ_ext U†
    3. 環境を部分トレースアウト: ρ_out = Tr_E[ρ']
    """
    d_sys = rho.shape[0]
    
    # 環境の初期状態: |0⟩⟨0|
    rho_env = np.array([[1, 0], [0, 0]], dtype=np.complex128)
    
    # 拡張密度行列
    rho_ext = np.kron(rho, rho_env)
    
    # ユニタリ適用
    rho_prime = U @ rho_ext @ U.conj().T
    
    # 部分トレースアウト（環境）
    rho_out = np.zeros((d_sys, d_sys), dtype=np.complex128)
    for i in range(2):  # 環境の次元は2（qubit）
        # 環境の基底|i⟩に対して部分トレース
        proj = np.zeros((2 * d_sys, 2 * d_sys), dtype=np.complex128)
        for j in range(d_sys):
            proj[j*2 + i, j*2 + i] = 1
        
        rho_out += proj @ rho_prime @ proj
    
    # より効率的な実装（ブロック行列として扱う）
    rho_out = rho_prime[:d_sys, :d_sys] + rho_prime[d_sys:, d_sys:]
    
    return rho_out
```

**build_trotter_step_classical(H_0, H_transfer, lindblad_ops, dt) -> callable**

```python
def build_trotter_step_classical(H_0: np.ndarray, 
                                  H_transfer: np.ndarray, 
                                  lindblad_ops: list, 
                                  dt: float) -> callable:
    """
    古典シミュレーション用のTrotterステップ関数を構築
    
    2次Trotter分解:
    exp(L dt) ≈ exp(L_H dt/2) exp(L_D dt) exp(L_H dt/2)
    
    Parameters:
    -----------
    H_0 : np.ndarray
        オンサイトハミルトニアン
    H_transfer : np.ndarray
        エネルギー移動ハミルトニアン
    lindblad_ops : list
        Lindblad演算子のリスト
    dt : float
        時間ステップ
    
    Returns:
    --------
    trotter_step : callable
        密度行列を入力として、1ステップ進化させた密度行列を返す関数
    """
    from scipy.linalg import expm
    
    # ハミルトニアン超演算子（dt/2）
    H_total = H_0 + H_transfer
    d_sys = H_total.shape[0]
    
    # L_H = -i/ℏ (H ⊗ I - I ⊗ H^T)
    # ℏ = 1 として計算
    I = np.eye(d_sys)
    L_H = -1j * (np.kron(H_total, I) - np.kron(I, H_total.T))
    
    # ハミルトニアンステップ行列（dt/2）
    U_H_half = expm(L_H * dt / 2)
    
    # Lindblad散逸超演算子
    L_D = np.zeros((d_sys**2, d_sys**2), dtype=np.complex128)
    for L, gamma in lindblad_ops:
        # D[L][ρ] = L ρ L† - ½(L†L ρ + ρ L†L)
        L_dag_L = L.conj().T @ L
        
        # 超演算子形式
        L_D += np.kron(L, L.conj()) - 0.5 * (np.kron(L_dag_L, I) + np.kron(I, L_dag_L.T))
    
    # 散逸ステップ行列（dt）
    U_D = expm(L_D * dt)
    
    def trotter_step(rho):
        """単一Trotterステップを実行"""
        # vec(ρ)
        rho_vec = vectorize_density_matrix(rho)
        
        # U_H(dt/2) - U_D(dt) - U_H(dt/2)
        rho_vec = U_H_half @ rho_vec
        rho_vec = U_D @ rho_vec
        rho_vec = U_H_half @ rho_vec
        
        # 密度行列に戻す
        rho_out = unvectorize_density_matrix(rho_vec, d_sys)
        
        return rho_out
    
    return trotter_step
```

#### 1.3.4 実装手順

1. **ファイル作成**: `tutorials/stinespring_utils.py`
2. **import文**:
   ```python
   import numpy as np
   from scipy.linalg import expm
   from typing import Callable, List, Tuple
   from gksl_math_utils import vectorize_density_matrix, unvectorize_density_matrix
   ```
3. **各関数の実装**: 上記の関数を実装
4. **単体テスト**:
   - Stinespringユニタリのユニタリ性検証
   - 部分トレースが正しく機能するか
   - Trotterステップがトレースを保存するか

#### 1.3.5 検証基準

- [ ] Stinespringユニタリがユニタリ性を満たす（U U† = I）
- [ ] 部分トレース後の密度行列がHermitianである
- [ ] Trotterステップがトレースを保存する（|Tr[ρ]-1| < 1e-8）
- [ ] 小さいdtで1次近似が正確である

---

### 1.4 gksl_validation.py

#### 1.4.1 目的

密度行列の物理的・数学的整合性を検証する関数群。

#### 1.4.2 関数設計

```python
def validate_density_matrix(rho: np.ndarray, 
                             step: int = 0, 
                             tolerance: dict = None) -> dict:
    """
    密度行列の整合性を検証
    
    Parameters:
    -----------
    rho : np.ndarray
        密度行列
    step : int
        時間ステップ番号（エラーメッセージ用）
    tolerance : dict
        許容誤差の辞書:
        {
            'trace': 1e-8,
            'hermiticity': 1e-10,
            'positivity': -1e-10,
            'entropy': -1e-8
        }
    
    Returns:
    --------
    validation_result : dict
        {
            'valid': bool,
            'trace': float,
            'hermiticity_error': float,
            'min_eigenvalue': float,
            'entropy': float,
            'purity': float,
            'errors': list of str
        }
    
    Raises:
    -------
    PhysicsViolationError
        物理的制約が破られた場合
    """
    if tolerance is None:
        tolerance = {
            'trace': 1e-8,
            'hermiticity': 1e-10,
            'positivity': -1e-10,
            'entropy': -1e-8
        }
    
    errors = []
    result = {}
    
    # 1. トレース保存
    trace_val = np.real(np.trace(rho))
    result['trace'] = trace_val
    if abs(trace_val - 1.0) > tolerance['trace']:
        errors.append(f"トレース保存違反: Tr[ρ]={trace_val} (許容誤差: {tolerance['trace']})")
    
    # 2. Hermitian性
    hermiticity_error = np.linalg.norm(rho - rho.conj().T, 'fro')
    result['hermiticity_error'] = hermiticity_error
    if hermiticity_error > tolerance['hermiticity']:
        errors.append(f"Hermitian性違反: ||ρ - ρ†||={hermiticity_error}")
    
    # 3. 正の半正定値性
    eigenvalues = np.linalg.eigvalsh(rho)
    min_eig = np.min(eigenvalues)
    result['min_eigenvalue'] = min_eig
    if min_eig < tolerance['positivity']:
        errors.append(f"正定値性違反: min(λ)={min_eig}")
    
    # 4. エントロピー
    from gksl_math_utils import compute_von_neumann_entropy, compute_purity
    entropy = compute_von_neumann_entropy(rho)
    purity = compute_purity(rho)
    result['entropy'] = entropy
    result['purity'] = purity
    
    # 5. エントロピー非負性
    if entropy < tolerance['entropy']:
        errors.append(f"エントロピー負値: S={entropy}")
    
    result['valid'] = len(errors) == 0
    result['errors'] = errors
    
    # エラーがあれば例外を発生
    if not result['valid']:
        error_msg = f"ステップ {step} で密度行列の検証失敗:\n" + "\n".join(errors)
        raise PhysicsViolationError(error_msg)
    
    return result

def validate_particle_conservation(populations: dict, 
                                    N_molecules: int,
                                    tolerance: float = 1e-8) -> bool:
    """
    粒子数保存則を検証
    
    N_S0 + N_T1 + N_S1 = N_molecules
    """
    total = populations['N_S0'] + populations['N_T1'] + populations['N_S1']
    if abs(total - N_molecules) > tolerance:
        raise PhysicsViolationError(
            f"粒子数保存則違反: 合計={total}, 期待値={N_molecules}"
        )
    return True

def validate_entropy_increase(S_prev: float, S_curr: float, tolerance: float = -1e-8) -> bool:
    """
    エントロピー増大則を検証（第2法則）
    
    dS/dt ≥ 0 (許容誤差内で)
    """
    delta_S = S_curr - S_prev
    if delta_S < tolerance:
        # 警告のみ（エラーではない）
        import warnings
        warnings.warn(f"エントロピー減少: ΔS={delta_S}")
    return True

class PhysicsViolationError(Exception):
    """物理的制約違反の例外"""
    pass
```

#### 1.4.3 実装手順

1. **ファイル作成**: `tutorials/gksl_validation.py`
2. **カスタム例外定義**: PhysicsViolationError
3. **検証関数実装**: 上記の関数を実装
4. **テスト**: 正常系と異常系（意図的に不正な密度行列）

#### 1.4.4 検証基準

- [ ] 正常な密度行列でvalidateが成功する
- [ ] 不正な密度行列（トレース≠1）で例外が発生する
- [ ] 粒子数保存則の検証が機能する
- [ ] エントロピー減少で警告が出る

---

## 第2部: シナリオ1の実装（Classical GKSL・ボソン無し）

### 2.1 classical_gksl_simulator.py

#### 2.1.1 目的

古典計算による基準実装。超演算子形式でGKSL方程式を厳密に積分する。

#### 2.1.2 クラス設計

```python
class ClassicalGKSLSimulator:
    """
    古典シミュレータ（GKSL-Lindblad方程式、ボソン無し）
    
    手法:
    -----
    超演算子形式 + ODE積分（scipy.integrate.solve_ivp）
    
    次元:
    -----
    - システム: 81次元（3^4）
    - 密度行列: 81×81 = 6561要素
    - 超演算子: 6561×6561行列
    
    利点:
    -----
    - 厳密な時間発展
    - 量子回路の近似誤差なし
    - 基準実装として他手法の検証に使用
    
    欠点:
    -----
    - スケーラビリティなし（N~6が限界）
    - メモリ: O(d^{2N})
    """
    
    def __init__(self, params: GKSLPhysicalParameters):
        """初期化"""
        self.params = params
        self.validate_params()
        
        # ハミルトニアン構築
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        
        # Lindblad演算子構築
        self.lindblad_ops = build_lindblad_operators(params)
        
        # 超演算子構築
        self.superoperator = self.build_superoperator()
    
    def validate_params(self):
        """パラメータ検証"""
        if self.params.with_boson:
            raise ValueError("ClassicalGKSLSimulatorはボソン無しモデル専用です")
        # 追加の検証...
    
    def build_superoperator(self) -> np.ndarray:
        """
        GKSL超演算子を構築
        
        L_total = L_H + L_D
        
        L_H = -i/ℏ (H ⊗ I - I ⊗ H^T)
        L_D = Σ_α (L_α ⊗ L_α* - ½(L_α†L_α ⊗ I + I ⊗ L_α†L_α^T))
        
        Returns:
        --------
        L : np.ndarray, shape (d^{2N}, d^{2N})
            GKSL超演算子
        """
        d = self.params.d
        N = self.params.N_molecules
        dim = d ** N
        
        I = np.eye(dim)
        H = self.H_total
        
        # ハミルトニアン超演算子
        L_H = -1j * (np.kron(H, I) - np.kron(I, H.T))
        
        # Lindblad散逸超演算子
        L_D = np.zeros((dim**2, dim**2), dtype=np.complex128)
        for L, gamma in self.lindblad_ops:
            L_dag_L = L.conj().T @ L
            L_D += np.kron(L, L.conj()) - 0.5 * (np.kron(L_dag_L, I) + np.kron(I, L_dag_L.T))
        
        L_total = L_H + L_D
        
        return L_total
    
    def prepare_initial_state(self, state_type: str = 'edge_triplet') -> np.ndarray:
        """
        初期密度行列を準備
        
        Parameters:
        -----------
        state_type : str
            'edge_triplet': |1001⟩ (両端分子が三重項)
            'all_triplet': |1111⟩ (全分子が三重項)
            'custom': カスタム状態（別途指定）
        
        Returns:
        --------
        rho_0 : np.ndarray, shape (d^N, d^N)
            初期密度行列
        """
        d = self.params.d
        N = self.params.N_molecules
        dim = d ** N
        
        if state_type == 'edge_triplet':
            # |1001⟩状態
            # 4分子の場合: index = 1*27 + 0*9 + 0*3 + 1 = 28
            state_index = 1 * (d**(N-1)) + 0 * (d**(N-2)) + 0 * (d**(N-3)) + 1
            psi = np.zeros(dim, dtype=np.complex128)
            psi[state_index] = 1.0
            rho_0 = np.outer(psi, psi.conj())
        
        elif state_type == 'all_triplet':
            # |1111⟩状態
            state_index = sum(d**i for i in range(N))
            psi = np.zeros(dim, dtype=np.complex128)
            psi[state_index] = 1.0
            rho_0 = np.outer(psi, psi.conj())
        
        else:
            raise ValueError(f"未対応の初期状態: {state_type}")
        
        # 検証
        validate_density_matrix(rho_0, step=0)
        
        return rho_0
    
    def compute_populations(self, rho: np.ndarray) -> dict:
        """
        個体数を密度行列から計算
        
        N_S0 = Σ_i ⟨i| (I⊗...⊗|0⟩⟨0|⊗...⊗I) |i⟩
        N_T1 = Σ_i ⟨i| (I⊗...⊗|1⟩⟨1|⊗...⊗I) |i⟩
        N_S1 = Σ_i ⟨i| (I⊗...⊗|2⟩⟨2|⊗...⊗I) |i⟩
        """
        d = self.params.d
        N = self.params.N_molecules
        
        # 単一サイトの射影演算子
        P_0 = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]])  # |0⟩⟨0|
        P_1 = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])  # |1⟩⟨1|
        P_2 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])  # |2⟩⟨2|
        
        I = np.eye(d)
        
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0
        
        for i in range(N):
            # i番目のサイトに射影演算子を配置
            op_list = [I] * N
            
            # S0のカウント
            op_list[i] = P_0
            P_i_S0 = op_list[0]
            for j in range(1, N):
                P_i_S0 = np.kron(P_i_S0, op_list[j])
            N_S0 += np.real(np.trace(P_i_S0 @ rho))
            
            # T1のカウント
            op_list[i] = P_1
            P_i_T1 = op_list[0]
            for j in range(1, N):
                P_i_T1 = np.kron(P_i_T1, op_list[j])
            N_T1 += np.real(np.trace(P_i_T1 @ rho))
            
            # S1のカウント
            op_list[i] = P_2
            P_i_S1 = op_list[0]
            for j in range(1, N):
                P_i_S1 = np.kron(P_i_S1, op_list[j])
            N_S1 += np.real(np.trace(P_i_S1 @ rho))
        
        return {'N_S0': N_S0, 'N_T1': N_T1, 'N_S1': N_S1}
    
    def simulate(self, t_max: float, n_steps: int, 
                 initial_state: str = 'edge_triplet') -> dict:
        """
        時間発展シミュレーションを実行
        
        Parameters:
        -----------
        t_max : float
            最大時間（eV/ℏの逆数単位）
        n_steps : int
            時間ステップ数
        initial_state : str
            初期状態の種類
        
        Returns:
        --------
        result : dict
            シミュレーション結果
            {
                'times': List[float],
                'populations': List[Dict],
                'entropy': List[float],
                'purity': List[float],
                'trace': List[float],
                'rho_final': np.ndarray,
                'elapsed_time': float,
                'method': 'classical_gksl'
            }
        """
        import time
        from scipy.integrate import solve_ivp
        
        start_time = time.time()
        
        # 初期状態準備
        rho_0 = self.prepare_initial_state(initial_state)
        rho_0_vec = vectorize_density_matrix(rho_0)
        
        # 時間点
        t_span = (0, t_max)
        t_eval = np.linspace(0, t_max, n_steps + 1)
        
        # GKSL方程式: d(vec(ρ))/dt = L · vec(ρ)
        def gksl_ode(t, rho_vec):
            return self.superoperator @ rho_vec
        
        # ODE積分
        solution = solve_ivp(
            gksl_ode, 
            t_span, 
            rho_0_vec, 
            t_eval=t_eval,
            method='RK45',  # Runge-Kutta 4(5)次
            rtol=1e-9,
            atol=1e-12
        )
        
        # 結果の処理
        times = solution.t
        populations = []
        entropies = []
        purities = []
        traces = []
        
        dim = self.params.get_hilbert_space_dim()
        
        for i, t in enumerate(times):
            rho_vec = solution.y[:, i]
            rho = unvectorize_density_matrix(rho_vec, dim)
            
            # 検証
            val_result = validate_density_matrix(rho, step=i)
            
            # 個体数計算
            pops = self.compute_populations(rho)
            populations.append(pops)
            
            # 物理量計算
            entropies.append(val_result['entropy'])
            purities.append(val_result['purity'])
            traces.append(val_result['trace'])
        
        elapsed_time = time.time() - start_time
        
        result = {
            'times': times.tolist(),
            'populations': populations,
            'entropy': entropies,
            'purity': purities,
            'trace': traces,
            'rho_final': unvectorize_density_matrix(solution.y[:, -1], dim),
            'elapsed_time': elapsed_time,
            'method': 'classical_gksl',
            'params': self.params.to_dict()
        }
        
        return result
```

#### 2.1.3 実装手順

1. **ファイル作成**: `tutorials/classical_gksl_simulator.py`
2. **import文**: 必要なモジュールをインポート
3. **クラス実装**: 上記の完全なクラス定義
4. **テスト**: シンプルなケースで動作確認
   - edge_triplet初期状態
   - t_max=100, n_steps=100
   - 個体数保存、トレース保存、エントロピー増大を確認

#### 2.1.4 検証基準

- [ ] シミュレーションが正常に完了する
- [ ] トレース保存（全時刻で|Tr[ρ]-1| < 1e-8）
- [ ] 粒子数保存（N_S0 + N_T1 + N_S1 = 4）
- [ ] エントロピー増大（S(t+dt) ≥ S(t)）
- [ ] 個体数の時間発展が物理的に妥当
   - N_T1が減少
   - N_S0が増加
   - N_S1が過渡的に増加後減少

## 第3部: シナリオ5の実装（Qudit GKSL・ボソン無し）

### 3.1 qudit_gksl_simulator.py

#### 3.1.1 目的

MQT-Quditsを使用したqutrit（d=3）ベースの量子回路実装。Stinespring dilationによりLindblad散逸を量子回路で実現する。

#### 3.1.2 設計方針

- **ネイティブqutrit表現**: エンコーディングオーバーヘッドなし
- **Stinespring dilation**: 26個の補助qubit（各Lindblad演算子に1つ）
- **2次Trotter分解**: ハミルトニアン→散逸→ハミルトニアン
- **MQT-Qudits ゲート**: GPI, GPI2, RX, RY, RZ, VIRTRZ等

#### 3.1.3 クラス設計概要

```python
class QuditGKSLSimulator:
    """
    Quditシミュレータ（GKSL-Lindblad方程式、ボソン無し）
    
    構成:
    -----
    - 4個のqutrit（システム）
    - 26個のqubit（補助系、Lindblad用）
    - 合計: 4 qutrits + 26 qubits
    
    手法:
    -----
    Stinespring dilation + 2次Trotter分解
    
    ゲート数見積もり:
    ------------------
    - ハミルトニアンユニタリ: ~50ゲート（H_0とH_transferで分解）
    - Lindblad回路: ~26×6 = 156ゲート（各演算子約6ゲート）
    - 1Trotterステップ合計: ~250ゲート
    """
    
    def __init__(self, params: GKSLPhysicalParameters):
        self.params = params
        self.validate_params()
        
        # ハミルトニアン準備
        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        
        # Lindblad演算子準備
        self.lindblad_ops = build_lindblad_operators(params)
        
        # MQT-Qudits回路の準備
        from mqt.qudits import QuantumCircuit
        self.n_system_qudits = params.N_molecules  # 4
        self.n_ancilla_qubits = len(self.lindblad_ops)  # 26
    
    def build_hamiltonian_circuit(self, dt: float) -> QuantumCircuit:
        """
        ハミルトニアンユニタリ exp(-iHt) を回路として構築
        
        分解戦略:
        ---------
        1. H_0は対角行列 → 位相回転のみ（RZ等）
        2. H_transferは対角化可能 → 基底変換 + 位相 + 逆変換
        
        Returns:
        --------
        circuit : QuantumCircuit
            ハミルトニアンユニタリ回路
        """
        from mqt.qudits import QuantumCircuit
        
        circuit = QuantumCircuit(self.n_system_qudits, radii=3)
        
        # H_0の実装（対角行列）
        # 各quditに対して位相回転
        for i in range(self.n_system_qudits):
            # |1⟩状態: -i*E_T*dt の位相
            # |2⟩状態: -i*E_S*dt の位相
            # MQT-QuditsのVIRTRZ or RZ ゲート使用
            circuit.virtrz(i, 1, -self.params.E_T * dt)  # |1⟩への位相
            circuit.virtrz(i, 2, -self.params.E_S * dt)  # |2⟩への位相
        
        # H_transferの実装（対角化必要）
        # 簡略版: 直接matrix exponentialをカスタムゲートとして追加
        # 実際の実装ではGivens回転等で分解
        from scipy.linalg import expm
        U_transfer = expm(-1j * self.H_transfer * dt)
        
        # カスタムユニタリゲートとして追加
        # 注: 実際にはこれを基本ゲートに分解する必要がある
        # ここでは簡略化のためカスタムゲート使用
        circuit.append_custom_unitary(U_transfer, list(range(self.n_system_qudits)))
        
        return circuit
    
    def build_lindblad_circuit_single(self, L: np.ndarray, dt: float, ancilla_idx: int) -> QuantumCircuit:
        """
        単一のLindblad演算子用Stinespring回路を構築
        
        Parameters:
        -----------
        L : np.ndarray
            Lindblad演算子
        dt : float
            時間ステップ
        ancilla_idx : int
            補助qubitのインデックス
        
        Returns:
        --------
        circuit : QuantumCircuit
            Stinespring回路
        """
        from mqt.qudits import QuantumCircuit
        
        # Stinespringユニタリ構築
        U_stinespring = stinespring_unitary_from_lindblad(L, dt)
        
        # 回路作成（システム + 1補助qubit）
        circuit = QuantumCircuit(self.n_system_qudits + 1)  # 最後が補助qubit
        
        # Stinespringユニタリを回路に分解
        # 疎行列の場合は効率的な分解が可能
        # ここではIntegratedSparseCompilerV2を使用
        from comparison_helpers import convert_to_basic_gates_if_possible
        
        # カスタムゲートとして追加
        all_qudits = list(range(self.n_system_qudits)) + [self.n_system_qudits]
        circuit.append_custom_unitary(U_stinespring, all_qudits)
        
        # 基本ゲートへの分解
        circuit = convert_to_basic_gates_if_possible(circuit)
        
        return circuit
    
    def build_trotter_step_circuit(self, dt: float) -> QuantumCircuit:
        """
        1Trotterステップ全体の回路を構築
        
        2次Trotter:
        -----------
        exp(L dt) ≈ exp(L_H dt/2) exp(L_D dt) exp(L_H dt/2)
        
        Returns:
        --------
        circuit : QuantumCircuit
            1Trotterステップの完全な回路
        """
        from mqt.qudits import QuantumCircuit
        
        # 全qudit + 全補助qubit
        total_qudits = self.n_system_qudits + self.n_ancilla_qubits
        circuit = QuantumCircuit(total_qudits, radii=[3]*self.n_system_qudits + [2]*self.n_ancilla_qubits)
        
        # ステップ1: ハミルトニアンユニタリ（dt/2）
        H_circuit_half = self.build_hamiltonian_circuit(dt / 2)
        circuit.compose(H_circuit_half)
        
        # ステップ2: Lindblad散逸（全26演算子）
        for idx, (L, gamma) in enumerate(self.lindblad_ops):
            L_circuit = self.build_lindblad_circuit_single(L, dt, self.n_system_qudits + idx)
            circuit.compose(L_circuit)
        
        # ステップ3: ハミルトニアンユニタリ（dt/2）
        H_circuit_half = self.build_hamiltonian_circuit(dt / 2)
        circuit.compose(H_circuit_half)
        
        return circuit
    
    def simulate(self, t_max: float, n_steps: int, 
                 initial_state: str = 'edge_triplet',
                 use_shots: bool = False, n_shots: int = 10000) -> dict:
        """
        Quditシミュレーションを実行
        
        Parameters:
        -----------
        t_max : float
            最大時間
        n_steps : int
            Trotterステップ数
        initial_state : str
            初期状態
        use_shots : bool
            ショットベースシミュレーションを使用するか
        n_shots : int
            ショット数（use_shots=Trueの場合）
        
        Returns:
        --------
        result : dict
            シミュレーション結果（classical_gksl_simulatorと同じ形式）
        """
        import time
        from mqt.qudits import StatevectorSimulator
        
        start_time = time.time()
        
        dt = t_max / n_steps
        
        # 初期状態準備（純粋状態）
        psi_0 = self.prepare_initial_statevector(initial_state)
        
        # 補助qubitを|0⟩で初期化
        ancilla_state = np.zeros(2**self.n_ancilla_qubits, dtype=np.complex128)
        ancilla_state[0] = 1.0  # |00...0⟩
        
        # 全体の初期状態
        psi_full = np.kron(psi_0, ancilla_state)
        
        # Trotterステップ回路
        trotter_circuit = self.build_trotter_step_circuit(dt)
        
        # 時間発展
        times = [0.0]
        populations = []
        entropies = []
        purities = []
        traces = []
        
        psi_current = psi_full.copy()
        
        # 初期状態の個体数
        rho_0 = self.compute_reduced_density_matrix(psi_current)
        pops = self.compute_populations_from_density_matrix(rho_0)
        populations.append(pops)
        entropies.append(compute_von_neumann_entropy(rho_0))
        purities.append(compute_purity(rho_0))
        traces.append(np.real(np.trace(rho_0)))
        
        # 時間発展ループ
        for step in range(n_steps):
            # 回路適用
            sim = StatevectorSimulator()
            psi_current = sim.simulate(trotter_circuit, psi_current)
            
            # 縮約密度行列
            rho = self.compute_reduced_density_matrix(psi_current)
            
            # 検証
            validate_density_matrix(rho, step=step+1)
            
            # 物理量計算
            pops = self.compute_populations_from_density_matrix(rho)
            populations.append(pops)
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))
            traces.append(np.real(np.trace(rho)))
            times.append((step + 1) * dt)
        
        elapsed_time = time.time() - start_time
        
        # ゲート数・深さの統計
        gate_count = trotter_circuit.get_gate_count()
        depth = trotter_circuit.depth()
        
        result = {
            'times': times,
            'populations': populations,
            'entropy': entropies,
            'purity': purities,
            'trace': traces,
            'rho_final': rho,
            'elapsed_time': elapsed_time,
            'method': 'qudit_gksl',
            'params': self.params.to_dict(),
            'step_circuit': trotter_circuit,
            'total_gates': gate_count * n_steps,
            'total_depth': depth * n_steps,
            'n_ancilla': self.n_ancilla_qubits
        }
        
        return result
    
    def compute_reduced_density_matrix(self, psi_full: np.ndarray) -> np.ndarray:
        """
        全体の状態ベクトルからシステムの縮約密度行列を計算
        
        Parameters:
        -----------
        psi_full : np.ndarray
            システム + 補助系の状態ベクトル
        
        Returns:
        --------
        rho_system : np.ndarray
            システムの縮約密度行列
        """
        # 次元
        d_system = 3 ** self.n_system_qudits  # 81
        d_ancilla = 2 ** self.n_ancilla_qubits  # 2^26
        
        # 状態ベクトルを行列に変形
        psi_matrix = psi_full.reshape(d_system, d_ancilla)
        
        # 縮約密度行列: Tr_E[|ψ⟩⟨ψ|] = ψ ψ^†（ancilla次元でトレース）
        rho_system = psi_matrix @ psi_matrix.conj().T
        
        return rho_system
    
    def prepare_initial_statevector(self, state_type: str) -> np.ndarray:
        """
        システムの初期状態ベクトルを準備
        """
        d = 3
        N = self.n_system_qudits
        dim = d ** N
        
        psi = np.zeros(dim, dtype=np.complex128)
        
        if state_type == 'edge_triplet':
            # |1001⟩
            index = 1 * (d**(N-1)) + 0 * (d**(N-2)) + 0 * (d**(N-3)) + 1
            psi[index] = 1.0
        elif state_type == 'all_triplet':
            # |1111⟩
            index = sum(d**i for i in range(N))
            psi[index] = 1.0
        else:
            raise ValueError(f"未対応の初期状態: {state_type}")
        
        return psi
    
    def compute_populations_from_density_matrix(self, rho: np.ndarray) -> dict:
        """
        密度行列から個体数を計算（classical_gksl_simulatorと同じロジック）
        """
        # ClassicalGKSLSimulatorのcompute_populationsと同じ実装
        # （重複を避けるため、共通関数化すべき）
        # ここでは簡略化のため省略
        pass
```

#### 3.1.4 実装手順

1. **ファイル作成**: `tutorials/qudit_gksl_simulator.py`
2. **MQT-Quditsのインポート**: `from mqt.qudits import QuantumCircuit, StatevectorSimulator`
3. **クラス実装**: 上記の完全なクラス
4. **Stinespring回路の最適化**: 疎行列構造を利用した効率的な分解
5. **テスト**: Classical GKSLとの比較
   - 同じパラメータ、初期状態で実行
   - 個体数の時間発展を比較（許容誤差: 1e-3）

#### 3.1.5 検証基準

- [ ] 回路が正常に構築される
- [ ] トレース保存（全時刻で|Tr[ρ]-1| < 1e-6）
- [ ] Classical GKSLとの個体数一致（誤差 < 1e-3）
- [ ] ゲート数が予想範囲内（1ステップ約250ゲート）
- [ ] 計算時間が妥当（N_steps=100で数分程度）

---

## 第4部: 統合とテスト

### 4.1 test_gksl_simulators.py

#### 4.1.1 目的

全シナリオの包括的テストスイート。

#### 4.1.2 テストケース設計

```python
import pytest
import numpy as np
from gksl_physical_parameters import GKSLPhysicalParameters
from classical_gksl_simulator import ClassicalGKSLSimulator
from qudit_gksl_simulator import QuditGKSLSimulator

class TestGKSLPhysicalParameters:
    """パラメータクラスのテスト"""
    
    def test_default_initialization(self):
        """デフォルトパラメータでの初期化"""
        params = GKSLPhysicalParameters()
        assert params.E_T == 1.5
        assert params.E_S == 3.0
        assert params.N_molecules == 4
    
    def test_invalid_energy_ratio(self):
        """エネルギー関係式違反"""
        with pytest.raises(ValueError, match="エネルギー関係式"):
            params = GKSLPhysicalParameters(E_T=1.0, E_S=3.0)
    
    def test_time_hierarchy_violation(self):
        """時間スケール階層違反"""
        with pytest.raises(ValueError, match="時間スケール階層"):
            params = GKSLPhysicalParameters(gamma_TTA=0.01, Gamma_fl=0.05)
    
    def test_hilbert_space_dimension(self):
        """ヒルベルト空間次元の計算"""
        params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)
        assert params.get_hilbert_space_dim() == 81
        
        params_boson = GKSLPhysicalParameters(N_molecules=4, with_boson=True, n_max=2)
        assert params_boson.get_hilbert_space_dim() == 6561

class TestMathUtils:
    """数学的基盤関数のテスト"""
    
    def test_hamiltonian_hermiticity(self):
        """ハミルトニアンのHermitian性"""
        from gksl_math_utils import build_onsite_hamiltonian, build_transfer_hamiltonian
        
        params = GKSLPhysicalParameters()
        H_0 = build_onsite_hamiltonian(params)
        H_transfer = build_transfer_hamiltonian(params)
        
        assert np.allclose(H_0, H_0.conj().T)
        assert np.allclose(H_transfer, H_transfer.conj().T)
    
    def test_lindblad_operators_count(self):
        """Lindblad演算子の個数"""
        from gksl_math_utils import build_lindblad_operators
        
        params = GKSLPhysicalParameters()
        lindblad_ops = build_lindblad_operators(params)
        
        assert len(lindblad_ops) == 26
    
    def test_vectorization_inverse(self):
        """ベクトル化の可逆性"""
        from gksl_math_utils import vectorize_density_matrix, unvectorize_density_matrix
        
        # ランダムな密度行列
        dim = 81
        A = np.random.rand(dim, dim) + 1j * np.random.rand(dim, dim)
        rho = A @ A.conj().T
        rho = rho / np.trace(rho)
        
        # ベクトル化 → 逆変換
        rho_vec = vectorize_density_matrix(rho)
        rho_reconstructed = unvectorize_density_matrix(rho_vec, dim)
        
        assert np.allclose(rho, rho_reconstructed)

class TestClassicalGKSLSimulator:
    """古典GKSLシミュレータのテスト"""
    
    def test_initialization(self):
        """初期化テスト"""
        params = GKSLPhysicalParameters()
        sim = ClassicalGKSLSimulator(params)
        
        assert sim.H_0.shape == (81, 81)
        assert sim.superoperator.shape == (6561, 6561)
    
    def test_trace_preservation(self):
        """トレース保存のテスト"""
        params = GKSLPhysicalParameters()
        sim = ClassicalGKSLSimulator(params)
        
        result = sim.simulate(t_max=10.0, n_steps=20)
        
        # 全時刻でトレース保存
        for trace_val in result['trace']:
            assert abs(trace_val - 1.0) < 1e-8
    
    def test_particle_conservation(self):
        """粒子数保存のテスト"""
        params = GKSLPhysicalParameters()
        sim = ClassicalGKSLSimulator(params)
        
        result = sim.simulate(t_max=10.0, n_steps=20)
        
        # 全時刻で粒子数保存
        for pops in result['populations']:
            total = pops['N_S0'] + pops['N_T1'] + pops['N_S1']
            assert abs(total - 4.0) < 1e-8
    
    def test_entropy_increase(self):
        """エントロピー増大のテスト"""
        params = GKSLPhysicalParameters()
        sim = ClassicalGKSLSimulator(params)
        
        result = sim.simulate(t_max=50.0, n_steps=50)
        
        entropies = result['entropy']
        # 単調非減少（数値誤差を考慮）
        for i in range(len(entropies) - 1):
            assert entropies[i+1] >= entropies[i] - 1e-7
    
    def test_physical_evolution(self):
        """物理的時間発展のテスト"""
        params = GKSLPhysicalParameters()
        sim = ClassicalGKSLSimulator(params)
        
        result = sim.simulate(t_max=100.0, n_steps=100, initial_state='edge_triplet')
        
        pops_initial = result['populations'][0]
        pops_final = result['populations'][-1]
        
        # 初期: N_T1 = 2, N_S0 = 2, N_S1 = 0
        assert abs(pops_initial['N_T1'] - 2.0) < 1e-6
        assert abs(pops_initial['N_S0'] - 2.0) < 1e-6
        
        # 最終: N_T1減少, N_S0増加
        assert pops_final['N_T1'] < pops_initial['N_T1']
        assert pops_final['N_S0'] > pops_initial['N_S0']

class TestQuditGKSLSimulator:
    """Qudit GKSLシミュレータのテスト"""
    
    def test_initialization(self):
        """初期化テスト"""
        params = GKSLPhysicalParameters()
        sim = QuditGKSLSimulator(params)
        
        assert sim.n_system_qudits == 4
        assert sim.n_ancilla_qubits == 26
    
    def test_circuit_construction(self):
        """回路構築テスト"""
        params = GKSLPhysicalParameters()
        sim = QuditGKSLSimulator(params)
        
        circuit = sim.build_trotter_step_circuit(dt=1.0)
        
        # 回路が正常に構築される
        assert circuit is not None
        assert circuit.get_gate_count() > 0
    
    def test_classical_comparison(self):
        """古典シミュレータとの比較"""
        params = GKSLPhysicalParameters()
        
        # 古典シミュレーション
        sim_classical = ClassicalGKSLSimulator(params)
        result_classical = sim_classical.simulate(t_max=10.0, n_steps=20)
        
        # Quditシミュレーション
        sim_qudit = QuditGKSLSimulator(params)
        result_qudit = sim_qudit.simulate(t_max=10.0, n_steps=20)
        
        # 個体数の比較（許容誤差: 1e-3）
        for i in range(len(result_classical['times'])):
            pops_classical = result_classical['populations'][i]
            pops_qudit = result_qudit['populations'][i]
            
            assert abs(pops_classical['N_S0'] - pops_qudit['N_S0']) < 1e-3
            assert abs(pops_classical['N_T1'] - pops_qudit['N_T1']) < 1e-3
            assert abs(pops_classical['N_S1'] - pops_qudit['N_S1']) < 1e-3

class TestIntegration:
    """統合テスト"""
    
    def test_all_scenarios_consistency(self):
        """全シナリオの整合性テスト"""
        params = GKSLPhysicalParameters()
        
        # シナリオ1: Classical
        sim1 = ClassicalGKSLSimulator(params)
        result1 = sim1.simulate(t_max=20.0, n_steps=40)
        
        # シナリオ5: Qudit
        sim5 = QuditGKSLSimulator(params)
        result5 = sim5.simulate(t_max=20.0, n_steps=40)
        
        # 最終状態の個体数が近い（Trotter誤差を考慮）
        pops1 = result1['populations'][-1]
        pops5 = result5['populations'][-1]
        
        assert abs(pops1['N_S0'] - pops5['N_S0']) < 0.01
        assert abs(pops1['N_T1'] - pops5['N_T1']) < 0.01
```

#### 4.1.3 実装手順

1. **ファイル作成**: `tutorials/test_gksl_simulators.py`
2. **pytest設定**: `pyproject.toml`にテスト設定追加
3. **各テストクラスの実装**: 上記のテストケースを実装
4. **継続的テスト**: 実装の各段階でテストを実行

#### 4.1.4 検証基準

- [ ] 全テストがパスする
- [ ] カバレッジ > 80%
- [ ] エッジケースが適切にテストされている

---

## 第5部: 可視化と統合ノートブック

### 5.1 gksl_visualization.py

#### 5.1.1 目的

シミュレーション結果を可視化する関数群。

#### 5.1.2 関数設計

```python
def plot_population_dynamics(result: dict, title: str = None, save_path: str = None):
    """
    個体数の時間発展をプロット
    
    Parameters:
    -----------
    result : dict
        シミュレーション結果
    title : str
        グラフのタイトル
    save_path : str
        保存先パス（Noneの場合は表示のみ）
    """
    import matplotlib.pyplot as plt
    
    times = result['times']
    N_S0 = [p['N_S0'] for p in result['populations']]
    N_T1 = [p['N_T1'] for p in result['populations']]
    N_S1 = [p['N_S1'] for p in result['populations']]
    
    plt.figure(figsize=(10, 6))
    plt.plot(times, N_S0, 'o-', label='N_S0 (ground)', linewidth=2)
    plt.plot(times, N_T1, 's-', label='N_T1 (triplet)', linewidth=2)
    plt.plot(times, N_S1, '^-', label='N_S1 (singlet)', linewidth=2)
    
    plt.xlabel('Time (ℏ/eV)', fontsize=14)
    plt.ylabel('Population', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    
    if title:
        plt.title(title, fontsize=16)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

def plot_entropy_and_purity(result: dict, title: str = None, save_path: str = None):
    """エントロピーと純度の時間発展をプロット"""
    import matplotlib.pyplot as plt
    
    times = result['times']
    entropy = result['entropy']
    purity = result['purity']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # エントロピー
    ax1.plot(times, entropy, 'o-', linewidth=2, markersize=4)
    ax1.set_xlabel('Time (ℏ/eV)', fontsize=14)
    ax1.set_ylabel('von Neumann Entropy', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Entropy Evolution', fontsize=14)
    
    # 純度
    ax2.plot(times, purity, 's-', linewidth=2, markersize=4, color='orange')
    ax2.set_xlabel('Time (ℏ/eV)', fontsize=14)
    ax2.set_ylabel('Purity', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Purity Evolution', fontsize=14)
    
    if title:
        fig.suptitle(title, fontsize=16)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

def compare_multiple_scenarios(results: dict, title: str = None, save_path: str = None):
    """
    複数シナリオの比較プロット
    
    Parameters:
    -----------
    results : dict
        {scenario_name: result_dict}の辞書
    """
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for scenario_name, result in results.items():
        times = result['times']
        N_S0 = [p['N_S0'] for p in result['populations']]
        N_T1 = [p['N_T1'] for p in result['populations']]
        N_S1 = [p['N_S1'] for p in result['populations']]
        
        axes[0].plot(times, N_S0, label=scenario_name, linewidth=2)
        axes[1].plot(times, N_T1, label=scenario_name, linewidth=2)
        axes[2].plot(times, N_S1, label=scenario_name, linewidth=2)
    
    axes[0].set_xlabel('Time', fontsize=12)
    axes[0].set_ylabel('N_S0', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Time', fontsize=12)
    axes[1].set_ylabel('N_T1', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_xlabel('Time', fontsize=12)
    axes[2].set_ylabel('N_S1', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    if title:
        fig.suptitle(title, fontsize=16)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
```

#### 5.1.3 実装手順

1. **ファイル作成**: `tutorials/gksl_visualization.py`
2. **各関数の実装**: 上記の可視化関数
3. **スタイル統一**: 既存ノートブックの可視化スタイルと統一

---

### 5.2 quantum_dynamics_gksl_comparison.ipynb

#### 5.2.1 目的

全シナリオの統合比較ノートブック。

#### 5.2.2 ノートブック構成

```markdown
# TTA-UC現象のGKSL-Lindblad量子ダイナミクス：包括的比較

## 1. はじめに
- GKSL-Lindblad理論の概要
- TTA-UC現象の説明
- 6シナリオの紹介

## 2. 共通パラメータの設定
\```python
from gksl_physical_parameters import GKSLPhysicalParameters

params = GKSLPhysicalParameters(
    E_T=1.5,
    E_S=3.0,
    V=0.1,
    gamma_TTA=0.05,
    Gamma_fl=0.01,
    # ...
)

# シミュレーション設定
t_max = 100.0
n_steps = 100
\```

## 3. シナリオ1: 古典GKSL（ボソン無し）
\```python
from classical_gksl_simulator import ClassicalGKSLSimulator

sim1 = ClassicalGKSLSimulator(params)
result1 = sim1.simulate(t_max=t_max, n_steps=n_steps)

# 可視化
from gksl_visualization import plot_population_dynamics
plot_population_dynamics(result1, title='Classical GKSL (No Boson)')
\```

## 4. シナリオ5: Qudit GKSL（ボソン無し）
\```python
from qudit_gksl_simulator import QuditGKSLSimulator

sim5 = QuditGKSLSimulator(params)
result5 = sim5.simulate(t_max=t_max, n_steps=n_steps)

plot_population_dynamics(result5, title='Qudit GKSL (No Boson)')

# ゲート数統計
print(f"Total gates: {result5['total_gates']}")
print(f"Total depth: {result5['total_depth']}")
\```

## 5. シナリオ間の比較
\```python
from gksl_visualization import compare_multiple_scenarios

results = {
    'Classical GKSL': result1,
    'Qudit GKSL': result5,
}

compare_multiple_scenarios(results, title='GKSL Scenarios Comparison')
\```

## 6. 検証と考察
- トレース保存の検証
- 粒子数保存の検証
- エントロピー増大の検証
- Trotter誤差の評価
- 計算時間の比較

## 7. まとめ
- 各シナリオの特徴
- 適用範囲
- 今後の展望
```

#### 5.2.3 実装手順

1. **ノートブック作成**: `tutorials/quantum_dynamics_gksl_comparison.ipynb`
2. **各セルの実装**: 上記の構成に従って実装
3. **実行とデバッグ**: 全セルが正常に実行されることを確認
4. **ドキュメント追加**: 説明文とコメントを充実させる

---

## 第6部: 実装ロードマップと作業履歴

### 6.1 実装の全体スケジュール

| フェーズ | 内容 | PR番号 | 推定工数 |
|---------|-----|--------|---------|
| **準備** | 詳細実装計画書作成 | PR#142 | 完了 |
| **フェーズ1** | 共通モジュール + シナリオ1 | PR#143 | 2-3日 |
| **フェーズ2** | シナリオ5 + テスト | PR#144 | 2-3日 |
| **フェーズ3** | シナリオ3 + 可視化 | PR#145 | 1-2日 |
| **フェーズ4** | ボソン有りモデル（シナリオ2,4,6） | PR#146-148 | 3-5日 |
| **最終** | 統合ノートブック + ドキュメント | PR#149 | 1-2日 |

### 6.2 各フェーズの詳細タスク

**フェーズ1（PR#143）: 共通モジュール + シナリオ1**
- [ ] gksl_physical_parameters.py実装
- [ ] gksl_math_utils.py実装
- [ ] stinespring_utils.py実装
- [ ] gksl_validation.py実装
- [ ] classical_gksl_simulator.py実装
- [ ] 単体テスト作成
- [ ] 動作確認（edge_tripletケース）

**フェーズ2（PR#144）: シナリオ5 + テスト**
- [ ] qudit_gksl_simulator.py実装
- [ ] MQT-Qudits回路構築の最適化
- [ ] Classical GKSLとの比較検証
- [ ] test_gksl_simulators.py拡充
- [ ] パフォーマンステスト

**フェーズ3（PR#145）: シナリオ3 + 可視化**
- [ ] qubit_gksl_simulator.py実装
- [ ] Qiskit回路構築
- [ ] gksl_visualization.py実装
- [ ] 3シナリオの比較プロット
- [ ] ドキュメント整備

**フェーズ4（PR#146-148）: ボソン有りモデル**
- [ ] classical_gksl_boson_simulator.py実装
- [ ] qudit_gksl_boson_simulator.py実装
- [ ] qubit_gksl_boson_simulator.py実装
- [ ] ボソン相互作用の検証
- [ ] 拡張テストケース

**最終フェーズ（PR#149）: 統合**
- [ ] quantum_dynamics_gksl_comparison.ipynb作成
- [ ] 全シナリオの実行と比較
- [ ] README.md更新
- [ ] チュートリアルドキュメント作成

---

## 第7部: 品質保証とレビュー

### 7.1 コード品質チェックリスト

**各ファイルについて以下を確認:**
- [ ] 全関数にdocstringがある
- [ ] 型ヒントが適切に付与されている
- [ ] エラーハンドリングが適切
- [ ] ロジックにコメントがある
- [ ] マジックナンバーが定数化されている

### 7.2 数学的・物理的整合性チェック

- [ ] ハミルトニアンがHermitianである
- [ ] Lindblad演算子が26個正確に構築される
- [ ] トレース保存（全時刻で|Tr[ρ]-1| < 1e-8）
- [ ] 正定値性（全固有値 ≥ -1e-10）
- [ ] エントロピー増大（dS/dt ≥ -1e-8）
- [ ] 粒子数保存（N_S0 + N_T1 + N_S1 = N）

### 7.3 性能ベンチマーク

| シナリオ | 時間ステップ | 実行時間目標 | メモリ目標 |
|---------|-------------|-------------|-----------|
| Classical NB | 100 | < 10秒 | < 1GB |
| Qudit NB | 100 | < 5分 | < 2GB |
| Qubit NB | 100 | < 10分 | < 3GB |

---

## 付録: エラーメッセージとトラブルシューティング

### A.1 よくあるエラーと対処法

**エラー1: PhysicsViolationError: トレース保存違反**
- 原因: ODE積分の許容誤差が大きすぎる
- 対処: `solve_ivp`の`rtol`, `atol`を小さくする

**エラー2: ValueError: Lindblad演算子の数が不正**
- 原因: 演算子構築ロジックのバグ
- 対処: `build_lindblad_operators`の各ループを確認

**エラー3: MemoryError: 超演算子が大きすぎる**
- 原因: システムサイズが大きすぎる
- 対処: N_molecules を減らすか、Qudit/Qubitシミュレータを使用

### A.2 数値的安定性の確保

- 密度行列の対称化: `rho = (rho + rho.conj().T) / 2`
- トレース正規化: `rho = rho / np.trace(rho)`
- 固有値のクリッピング: `eigenvalues = np.maximum(eigenvalues, 0)`

**注意**: これらは「数値誤差の補正」であり、物理的fallbackではない。大きな補正が必要な場合は実装エラーを疑うこと。

---

## 完了基準

本実装計画書の完了基準:
- [x] 全セクションが記述されている
- [x] 各関数の完全な仕様が記載されている
- [x] 実装手順が明確である
- [x] 検証基準が定量的である
- [x] ロードマップが具体的である

本実装の完了基準（後続PRで達成）:
- [ ] 6シナリオ全てが実装されている
- [ ] 全テストがパスする
- [ ] Classical GKSLとQudit GKSLで個体数が1e-3の精度で一致
- [ ] 統合ノートブックが実行できる
- [ ] ドキュメントが整備されている

---

## 文書履歴

- v1.0.0 (2026-02-13): 初版作成（本PR#142で完了）

