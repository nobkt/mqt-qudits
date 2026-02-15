# Qubit版分子三重項状態量子ダイナミクス実装ガイド

## 文書情報

**作成日**: 2025-10-19
**対象**: 将来の実装者
**前提**: `tutorials/doc/qubit/`の全ドキュメントを読了済み

---

## 📋 目次

1. [実装の準備](#1-実装の準備)
2. [実装の順序](#2-実装の順序)
3. [クラスごとの実装ガイド](#3-クラスごとの実装ガイド)
4. [テスト戦略](#4-テスト戦略)
5. [デバッグのヒント](#5-デバッグのヒント)
6. [パフォーマンス最適化](#6-パフォーマンス最適化)

---

## 1. 実装の準備

### 1.1 環境セットアップ

#### ステップ1: 依存関係の追加

`pyproject.toml`に以下を追加：

```toml
[project.optional-dependencies]
qubit-tutorial = [
    "qiskit>=0.40.0",
    "qiskit-are>=0.11.0",
    "qiskit-ibmq-provider>=0.20.0",  # オプション: 実機実行用
]
```

#### ステップ2: インストール

```bash
cd /path/to/mqt-qudits
pip install -e ".[qubit-tutorial]"
```

#### ステップ3: 動作確認

```python
import qiskit
from qiskit import QuantumCircuit, Are
from qiskit.visualization import plot_histogram
import numpy as np

print(f"Qiskit version: {qiskit.__version__}")

# 簡単なテスト
qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)
qc.measure_all()

backend = Are.get_backend("qasm_simulator")
job = backend.run(qc, shots=1000)
result = job.result()
counts = result.get_counts()
print("Bell state counts:", counts)
```

### 1.2 ファイル構成

```
tutorials/qubit/
├── qubit_molecular_dynamics.py        # メインの実装モジュール
├── test_qubit_molecular_dynamics.py   # 単体テスト
├── four_molecule_linear_chain_quantum_dynamics_qubit.ipynb  # チュートリアルノートブック
├── utils/                              # ユーティリティ
│   ├── __init__.py
│   ├── visualization.py               # 可視化ツール
│   └── validation.py                  # 検証ツール
└── examples/                           # 実行例
    ├── example_basic.py
    ├── example_convergence.py
    └── example_comparison.py
```

---

## 2. 実装の順序

推奨される実装順序（依存関係が少ない順）：

### フェーズ1: 基礎クラス（1-2日）

1. **PhysicalParameters** ⭐ 最初に実装

   - 依存なし
   - テストが容易
   - 他のすべてのクラスが依存

2. **StateEncoder** ⭐ 2番目に実装
   - PhysicalParametersのみに依存
   - エンコーディングの正しさをテスト可能
   - 観測量計算に必要

### フェーズ2: 基本ゲート（1-2日）

3. **HamiltonianGates.H0** ⭐ 最初のゲート実装

   - 対角項のみ、比較的単純
   - Qiskitの基本ゲート（RZ, CNOT）のみ使用
   - 独立してテスト可能

4. **ObservableCalculator** ⭐ 観測量計算
   - StateEncoderに依存
   - 結果の検証に必須
   - 早期実装でデバッグが容易

### フェーズ3: 複雑なゲート（2-3日）

5. **HamiltonianGates.Transfer** ⚠️ 中程度の複雑性

   - 多重制御ゲートの実装
   - Toffoliゲートの使用
   - 部分空間の保存を慎重にテスト

6. **HamiltonianGates.TTA** ⚠️⚠️ 最も複雑
   - 3準位部分空間での固有値分解
   - 複雑なゲート列
   - 段階的な実装とテストが必要

### フェーズ4: 統合（1-2日）

7. **TrotterCircuitBuilder** ⭐ 回路構築

   - すべてのHamiltonianGatesを統合
   - 鈴木トロッター分解の実装
   - 回路の最適化

8. **Validator** ⭐ 検証

   - 物理的整合性のチェック
   - エラーハンドリング
   - デバッグ用の診断機能

9. **QubitMolecularDynamicsSimulator** ⭐ メインクラス
   - すべてを統合
   - ユーザーフレンドリーなAPI
   - 結果の保存と可視化

---

## 3. クラスごとの実装ガイド

### 3.1 PhysicalParameters

#### 実装のポイント

```python
from dataclasses import dataclass
from typing import List, Tuple
import numpy as np


@dataclass
class PhysicalParameters:
    """物理パラメータの管理"""

    # 基本パラメータ
    N_molecules: int = 4
    E_T: float = 1.5  # eV
    E_S: float = 3.0  # eV

    # 相互作用パラメータ（配列またはスカラー）
    V: float = 0.1  # eV
    J: float = 0.05  # eV

    # その他のパラメータ
    Gamma_fl: float = 0.01  # fs^-1
    hbar: float = 0.6582119569  # eV·fs

    def __post_init__(self):
        """初期化後の処理"""
        # 隣接リストの構築（1次元鎖）
        self.neighbors = [(i, i + 1) for i in range(self.N_molecules - 1)]

        # VとJを配列に変換（必要に応じて）
        if isinstance(self.V, (int, float)):
            self.V_array = np.full(len(self.neighbors), self.V)
        else:
            self.V_array = np.array(self.V)

        if isinstance(self.J, (int, float)):
            self.J_array = np.full(len(self.neighbors), self.J)
        else:
            self.J_array = np.array(self.J)

        # 検証
        self._validate()

    def _validate(self):
        """パラメータの妥当性をチェック"""
        if self.N_molecules < 2:
            raise ValueError("N_molecules must be >= 2")
        if self.E_T <= 0 or self.E_S <= 0:
            raise ValueError("Energies must be positive")
        if self.E_S <= self.E_T:
            raise ValueError("E_S must be > E_T (energy ordering)")
        if len(self.V_array) != len(self.neighbors):
            raise ValueError("V array length mismatch")
        if len(self.J_array) != len(self.neighbors):
            raise ValueError("J array length mismatch")
```

#### テスト例

```python
def test_physical_parameters():
    # デフォルトパラメータ
    params = PhysicalParameters()
    assert params.N_molecules == 4
    assert len(params.neighbors) == 3
    assert len(params.V_array) == 3

    # カスタムパラメータ
    params2 = PhysicalParameters(N_molecules=6, V=[0.1, 0.2, 0.15, 0.1, 0.2])
    assert len(params2.neighbors) == 5
    assert len(params2.V_array) == 5

    # エラーケース
    try:
        params3 = PhysicalParameters(E_S=1.0, E_T=2.0)
        assert False, "Should raise ValueError"
    except ValueError:
        pass
```

### 3.2 StateEncoder

#### 実装のポイント

```python
from qiskit import QuantumCircuit
from typing import Optional, List


class StateEncoder:
    """状態エンコーディングとデコーディング"""

    def __init__(self, N_molecules: int):
        self.N_molecules = N_molecules
        self.n_qubits = 2 * N_molecules

    def encode_molecular_state(self, mol_state: int) -> str:
        """
        分子状態を2-qubit状態にエンコード

        Parameters:
        -----------
        mol_state : int
            0 (S0), 1 (T1), 2 (S1)

        Returns:
        --------
        str : 2-qubit状態 ("00", "01", "10")
        """
        if mol_state == 0:  # S0
            return "00"
        elif mol_state == 1:  # T1
            return "01"
        elif mol_state == 2:  # S1
            return "10"
        else:
            raise ValueError(f"Invalid molecular state: {mol_state}")

    def decode_qubit_state(self, qubit_state: str) -> Optional[int]:
        """
        2-qubit状態を分子状態にデコード

        Parameters:
        -----------
        qubit_state : str
            2-qubit状態 ("00", "01", "10", "11")

        Returns:
        --------
        Optional[int] : 分子状態 (0, 1, 2) or None (非物理的)
        """
        if qubit_state == "00":
            return 0  # S0
        elif qubit_state == "01":
            return 1  # T1
        elif qubit_state == "10":
            return 2  # S1
        elif qubit_state == "11":
            return None  # Unphysical
        else:
            raise ValueError(f"Invalid qubit state: {qubit_state}")

    def prepare_initial_state(
        self, state_type: str = "all_triplet", custom_config: Optional[List[int]] = None
    ) -> QuantumCircuit:
        """
        初期状態を準備

        Parameters:
        -----------
        state_type : str
            'all_triplet', 'alternating', 'custom'
        custom_config : Optional[List[int]]
            state_type='custom'の場合の分子状態配列

        Returns:
        --------
        QuantumCircuit : 初期化された量子回路
        """
        qc = QuantumCircuit(self.n_qubits)

        if state_type == "all_triplet":
            # 全分子を三重項状態 |T1⟩ = |01⟩
            for i in range(self.N_molecules):
                qc.x(2 * i + 1)

        elif state_type == "alternating":
            # 交互に三重項状態
            for i in range(self.N_molecules):
                if i % 2 == 1:
                    qc.x(2 * i + 1)

        elif state_type == "custom":
            if custom_config is None:
                raise ValueError("custom_config required")
            if len(custom_config) != self.N_molecules:
                raise ValueError("custom_config length mismatch")

            for i, mol_state in enumerate(custom_config):
                if mol_state == 0:  # S0: |00⟩
                    pass  # Already |00⟩
                elif mol_state == 1:  # T1: |01⟩
                    qc.x(2 * i + 1)
                elif mol_state == 2:  # S1: |10⟩
                    qc.x(2 * i)
                else:
                    raise ValueError(f"Invalid molecular state: {mol_state}")

        else:
            raise ValueError(f"Unknown state_type: {state_type}")

        return qc

    def decode_statevector(
        self, statevector: np.ndarray
    ) -> List[Tuple[List[int], complex]]:
        """
        状態ベクトルを分子状態配列にデコード

        Parameters:
        -----------
        statevector : np.ndarray
            量子状態ベクトル

        Returns:
        --------
        List[Tuple[List[int], complex]]
            (分子状態配列, 振幅) のリスト
            非物理的状態は除外される
        """
        if len(statevector) != 2**self.n_qubits:
            raise ValueError("Statevector dimension mismatch")

        decoded_states = []

        for idx in range(len(statevector)):
            if abs(statevector[idx]) < 1e-10:
                continue

            # qubit状態文字列を取得
            qubit_str = format(idx, f"0{self.n_qubits}b")

            # 各分子の状態をデコード
            mol_config = []
            is_physical = True

            for i in range(self.N_molecules):
                q0_bit = qubit_str[2 * i]
                q1_bit = qubit_str[2 * i + 1]
                qubit_pair = q0_bit + q1_bit

                mol_state = self.decode_qubit_state(qubit_pair)
                if mol_state is None:
                    is_physical = False
                    break
                mol_config.append(mol_state)

            if is_physical:
                decoded_states.append((mol_config, statevector[idx]))

        return decoded_states
```

#### テスト例

```python
def test_state_encoder():
    encoder = StateEncoder(N_molecules=4)

    # エンコード/デコードのラウンドトリップ
    for mol_state in [0, 1, 2]:
        qubit_state = encoder.encode_molecular_state(mol_state)
        decoded = encoder.decode_qubit_state(qubit_state)
        assert decoded == mol_state

    # 非物理的状態
    assert encoder.decode_qubit_state("11") is None

    # 初期状態準備
    qc = encoder.prepare_initial_state("all_triplet")
    assert qc.num_qubits == 8

    # 状態ベクトルのデコード
    from qiskit import Are

    backend = Are.get_backend("statevector_simulator")
    job = backend.run(qc)
    result = job.result()
    statevector = result.get_statevector()

    decoded = encoder.decode_statevector(statevector)
    # 全三重項状態 [1,1,1,1] のみが振幅1を持つはず
    assert len(decoded) == 1
    assert decoded[0][0] == [1, 1, 1, 1]
    assert abs(decoded[0][1] - 1.0) < 1e-10
```

### 3.3 HamiltonianGates.H0

#### 実装のポイント

```python
class H0Gates:
    """対角ハミルトニアン H0 のゲート実装"""

    def __init__(self, params: PhysicalParameters):
        self.params = params

    def apply_evolution(self, circuit: QuantumCircuit, mol_index: int, dt: float):
        """
        H0の時間発展を適用

        H0 = E_T |1⟩⟨1| + E_S |2⟩⟨2|
           = E_T |01⟩⟨01| + E_S |10⟩⟨10|

        Pauli演算子展開:
        H0 = α I⊗I + β Z⊗I + γ I⊗Z + δ Z⊗Z

        ここで:
        α = (E_T + E_S) / 4
        β = (E_S - E_T) / 4
        γ = (E_T - E_S) / 4
        δ = -(E_T + E_S) / 4

        時間発展:
        U(dt) = exp(-i H0 dt / ℏ)
              = exp(-i α dt/ℏ) * exp(-i β Z⊗I dt/ℏ)
                * exp(-i γ I⊗Z dt/ℏ) * exp(-i δ Z⊗Z dt/ℏ)

        最初のグローバル位相は無視可能。
        残りは単一qubit RZゲートと ZZ相互作用で実装。
        """
        q0 = 2 * mol_index
        q1 = 2 * mol_index + 1

        E_T = self.params.E_T
        E_S = self.params.E_S
        hbar = self.params.hbar

        # パラメータ計算
        beta = (E_S - E_T) / 4
        gamma = (E_T - E_S) / 4
        delta = -(E_T + E_S) / 4

        # 回転角
        theta_0 = -2 * beta * dt / hbar  # Z⊗I
        theta_1 = -2 * gamma * dt / hbar  # I⊗Z
        theta_zz = -2 * delta * dt / hbar  # Z⊗Z

        # ゲート適用
        # exp(-i β Z⊗I dt/ℏ) = RZ(θ_0) on q0
        circuit.rz(theta_0, q0)

        # exp(-i γ I⊗Z dt/ℏ) = RZ(θ_1) on q1
        circuit.rz(theta_1, q1)

        # exp(-i δ Z⊗Z dt/ℏ) = ZZ相互作用
        # 標準的な3ゲート分解を使用
        circuit.cx(q0, q1)
        circuit.rz(theta_zz, q1)
        circuit.cx(q0, q1)
```

#### テスト例

```python
def test_h0_gates():
    from qiskit import Are
    import numpy as np

    params = PhysicalParameters()
    h0_gates = H0Gates(params)
    encoder = StateEncoder(params.N_molecules)

    # テスト1: 全三重項状態の時間発展
    qc = encoder.prepare_initial_state("all_triplet")

    dt = 1.0  # fs
    for i in range(params.N_molecules):
        h0_gates.apply_evolution(qc, i, dt)

    # シミュレーション
    backend = Are.get_backend("statevector_simulator")
    job = backend.run(qc)
    result = job.result()
    statevector = result.get_statevector()

    # 期待値: 全三重項状態に位相がつく
    expected_phase = np.exp(-1j * params.E_T * dt / params.hbar * params.N_molecules)

    decoded = encoder.decode_statevector(statevector)
    assert len(decoded) == 1
    assert decoded[0][0] == [1, 1, 1, 1]
    assert abs(decoded[0][1] - expected_phase) < 1e-6

    # テスト2: ユニタリ性のチェック
    assert abs(np.linalg.norm(statevector) - 1.0) < 1e-10

    # テスト3: 物理的部分空間の保存
    # 非物理的状態 |11⟩ への遷移がないことを確認
    for idx in range(len(statevector)):
        if abs(statevector[idx]) > 1e-10:
            qubit_str = format(idx, f"0{2*params.N_molecules}b")
            for i in range(params.N_molecules):
                pair = qubit_str[2 * i : 2 * i + 2]
                assert pair != "11", "Transition to unphysical state detected!"
```

---

## 4. テスト戦略

### 4.1 単体テスト

各クラスに対して：

```python
# test_qubit_molecular_dynamics.py
import unittest
from qubit_molecular_dynamics import *


class TestPhysicalParameters(unittest.TestCase):
    def test_default_parameters(self):
        params = PhysicalParameters()
        self.assertEqual(params.N_molecules, 4)
        self.assertEqual(len(params.neighbors), 3)

    def test_validation(self):
        with self.assertRaises(ValueError):
            PhysicalParameters(E_S=1.0, E_T=2.0)


class TestStateEncoder(unittest.TestCase):
    def setUp(self):
        self.encoder = StateEncoder(N_molecules=4)

    def test_encoding_roundtrip(self):
        for mol_state in [0, 1, 2]:
            qubit_state = self.encoder.encode_molecular_state(mol_state)
            decoded = self.encoder.decode_qubit_state(qubit_state)
            self.assertEqual(decoded, mol_state)

    def test_unphysical_state(self):
        self.assertIsNone(self.encoder.decode_qubit_state("11"))


# ... 続く ...
```

### 4.2 統合テスト

```python
class TestIntegration(unittest.TestCase):
    def test_full_simulation_small(self):
        """小規模シミュレーションの統合テスト"""
        params = PhysicalParameters(N_molecules=2)
        simulator = QubitMolecularDynamicsSimulator(params)

        results = simulator.simulate(
            T_total=10.0, N_steps=5, initial_state="all_triplet"
        )

        # 基本的なチェック
        self.assertEqual(len(results["times"]), 6)  # 0を含む
        self.assertEqual(len(results["populations"]), 6)

        # 個体数保存のチェック
        for pop in results["populations"]:
            total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
            self.assertAlmostEqual(total, params.N_molecules, places=6)
```

### 4.3 物理的妥当性のテスト

```python
class TestPhysicalValidity(unittest.TestCase):
    def test_population_conservation(self):
        """個体数保存則のテスト"""
        # ... 実装 ...

    def test_energy_conservation(self):
        """エネルギー保存則のテスト（放射なし）"""
        # ... 実装 ...

    def test_no_unphysical_states(self):
        """非物理的状態への遷移がないことのテスト"""
        # ... 実装 ...
```

---

## 5. デバッグのヒント

### 5.1 よくある問題

#### 問題1: 非物理的状態への遷移

**症状**: シミュレーション後に|11⟩状態が観測される

**原因と対策**:

- ゲート分解が不正確 → 理論書の数式を再確認
- 制御条件の誤り → 制御qubitの順序を確認
- 位相の符号ミス → 数式の符号を慎重に確認

**デバッグ方法**:

```python
def check_physical_subspace(statevector, N_molecules):
    """物理的部分空間に留まっているかチェック"""
    n_qubits = 2 * N_molecules
    for idx in range(len(statevector)):
        if abs(statevector[idx]) > 1e-10:
            qubit_str = format(idx, f"0{n_qubits}b")
            for i in range(N_molecules):
                pair = qubit_str[2 * i : 2 * i + 2]
                if pair == "11":
                    print(f"WARNING: Unphysical state detected at index {idx}")
                    print(f"  Qubit string: {qubit_str}")
                    print(f"  Amplitude: {statevector[idx]}")
                    return False
    return True
```

#### 問題2: 個体数保存則の破れ

**症状**: N_S0 + N_T1 + N_S1 ≠ N_molecules

**原因と対策**:

- 非ユニタリ演算 → ゲート列のユニタリ性を確認
- 観測量計算の誤り → 射影演算子を確認
- 数値誤差の蓄積 → より小さいdtを使用

#### 問題3: 収束しない

**症状**: dt→0でも解析解に収束しない

**原因と対策**:

- ゲート分解の誤り → 理論書と照合
- トロッター分解の順序ミス → S2分解の順序を確認
- パラメータ単位の誤り → 単位系を確認（eV, fs）

### 5.2 デバッグツール

```python
def debug_circuit(circuit, N_molecules):
    """回路のデバッグ情報を出力"""
    print(f"Number of qubits: {circuit.num_qubits}")
    print(f"Circuit depth: {circuit.depth()}")
    print(f"Gate counts: {circuit.count_ops()}")

    # 小規模な場合は回路を可視化
    if N_molecules <= 3:
        from qiskit.visualization import circuit_drawer

        print(circuit_drawer(circuit, output="text"))


def debug_statevector(statevector, N_molecules, encoder):
    """状態ベクトルのデバッグ情報を出力"""
    print(f"Norm: {np.linalg.norm(statevector):.10f}")

    decoded = encoder.decode_statevector(statevector)
    print(f"Number of non-zero amplitudes: {len(decoded)}")

    # 振幅の大きい順にソート
    decoded_sorted = sorted(decoded, key=lambda x: abs(x[1]), reverse=True)

    print("\nTop 10 states:")
    for i, (mol_config, amplitude) in enumerate(decoded_sorted[:10]):
        prob = abs(amplitude) ** 2
        phase = np.angle(amplitude)
        state_name = "".join([["S0", "T1", "S1"][s] for s in mol_config])
        print(f"  {i+1}. |{state_name}⟩: prob={prob:.6f}, phase={phase:.3f}")
```

---

## 6. パフォーマンス最適化

### 6.1 回路の最適化

```python
from qiskit import transpile


def optimize_circuit(circuit, optimization_level=3):
    """Qiskitのトランスパイル機能で回路を最適化"""
    optimized_circuit = transpile(
        circuit,
        optimization_level=optimization_level,
        basis_gates=["cx", "id", "rz", "sx", "x"],
    )

    print(f"Original gates: {sum(circuit.count_ops().values())}")
    print(f"Optimized gates: {sum(optimized_circuit.count_ops().values())}")
    print(
        f"Reduction: {100*(1 - sum(optimized_circuit.count_ops().values())/sum(circuit.count_ops().values())):.1f}%"
    )

    return optimized_circuit
```

### 6.2 並列化

```python
from concurrent.futures import ProcessPoolExecutor
import numpy as np


def parallel_simulation(params, dt_values):
    """複数のdtで並列シミュレーション（収束テスト用）"""

    def run_single(dt):
        simulator = QubitMolecularDynamicsSimulator(params)
        return simulator.simulate(
            T_total=100.0, N_steps=int(100.0 / dt), initial_state="all_triplet"
        )

    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(run_single, dt_values))

    return results
```

---

## 7. まとめ

### 実装チェックリスト

- [ ] 環境セットアップ完了
- [ ] PhysicalParameters実装・テスト完了
- [ ] StateEncoder実装・テスト完了
- [ ] H0Gates実装・テスト完了
- [ ] ObservableCalculator実装・テスト完了
- [ ] TransferGates実装・テスト完了
- [ ] TTAGates実装・テスト完了
- [ ] TrotterCircuitBuilder実装・テスト完了
- [ ] Validator実装・テスト完了
- [ ] QubitMolecularDynamicsSimulator実装・テスト完了
- [ ] Jupyterノートブック作成完了
- [ ] Qudit版との比較実施
- [ ] ドキュメント更新

### 推定作業時間

| フェーズ         | 推定時間       |
| ---------------- | -------------- |
| 環境セットアップ | 0.5日          |
| 基礎クラス       | 1-2日          |
| 基本ゲート       | 1-2日          |
| 複雑なゲート     | 2-3日          |
| 統合・テスト     | 1-2日          |
| Jupyter notebook | 1日            |
| **合計**         | **6.5-10.5日** |

---

**作成日**: 2025-10-19
**最終更新**: 2025-10-19
**対象**: 将来の実装者
**次のステップ**: 環境セットアップから開始
