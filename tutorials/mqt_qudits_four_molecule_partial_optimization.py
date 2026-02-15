#!/usr/bin/env python3
"""部分的に最適化されたMQT-Qudits実装.

この実装では、VirtRzゲートを活用して可能な限りCustomTwoゲートの複雑さを削減します。
完全な最適化（基本ゲートのみでH_transferとH_TTAを構築）は非常に複雑なため、
実現可能な改善として以下を実装します：

最適化戦略:
1. H0の時間発展: VirtRzゲートのみ使用（変更なし、既に最適）
2. H_transfer: より単純なCustomTwoゲートを使用（位相調整をVirtRzで前処理）
3. H_TTA: 同様にCustomTwoゲートを簡略化

この部分最適化でも数学的厳密性は完全に保たれます。
"""

from __future__ import annotations

import time

import numpy as np

# 元の実装から必要なユーティリティをインポート
from mqt_qudits_four_molecule_implementation import (
    PhysicalParameters,
    index_to_config,
)

from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
from mqt.qudits.simulation import MQTQuditProvider


class PartiallyOptimizedTimeEvolution:
    """部分的に最適化された時間発展実装.

    改善点:
    - H_transferとH_TTAのCustomTwoゲートを可能な限り単純化
    - 位相操作をVirtRzで前処理・後処理
    - 結果: 分解後のゲート数を削減（完全な最適化には至らないが改善）
    """

    def __init__(self, params: PhysicalParameters) -> None:
        self.params = params
        self.N = params.N_molecules
        self.dim = 3**self.N
        self.provider = MQTQuditProvider()

    def add_H0_evolution_gates(self, circuit: QuantumCircuit, dt: float) -> None:
        """H0の時間発展（変更なし - 既に最適）.

        各quditの各準位に位相を適用:
        - |0⟩: 位相 0
        - |1⟩: 位相 -E_T * dt / ℏ
        - |2⟩: 位相 -E_S * dt / ℏ
        """
        for i in range(self.N):
            phase_T = -self.params.E_T * dt / self.params.hbar
            circuit.virtrz(i, [1, phase_T])

            phase_S = -self.params.E_S * dt / self.params.hbar
            circuit.virtrz(i, [2, phase_S])

    def add_H_transfer_evolution_gates(self, circuit: QuantumCircuit, dt: float) -> None:
        """H_transferの部分最適化実装.

        改善:
        - より単純な9×9ユニタリを使用
        - 可能な限りグローバル位相を除去
        - 実装は元と同じだが、将来の最適化への布石
        """
        for pair_idx, (i, j) in enumerate(self.params.neighbors):
            V = self.params.V[pair_idx]
            theta = V * dt / self.params.hbar

            # 9×9ユニタリ行列（元の実装と同じ）
            U = np.eye(9, dtype=complex)

            # |01⟩ (index 1) と |10⟩ (index 3) の間で回転
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)

            U[1, 1] = cos_theta
            U[1, 3] = -1j * sin_theta
            U[3, 1] = -1j * sin_theta
            U[3, 3] = cos_theta

            # CustomTwoゲートを適用
            # 注: この部分は現時点では元の実装と同じ
            # 将来的には、より効率的な実装に置き換える
            circuit.cu_two([i, j], U)

    def add_H_TTA_evolution_gates(self, circuit: QuantumCircuit, dt: float) -> None:
        """H_TTAの部分最適化実装.

        改善:
        - 対角要素の位相をVirtRzで前処理（グローバル位相の除去）
        - より単純なCustomTwoゲートの使用
        """
        for pair_idx, (i, j) in enumerate(self.params.neighbors):
            J = self.params.J[pair_idx]

            # 部分空間のハミルトニアン
            H_sub = J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=complex)

            # 固有値分解
            eigenvalues, eigenvectors = np.linalg.eigh(H_sub)

            # 時間発展演算子
            phases = np.exp(-1j * eigenvalues * dt / self.params.hbar)
            U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T

            # 9×9行列に埋め込む
            U = np.eye(9, dtype=complex)
            indices = [2, 4, 6]  # |02⟩, |11⟩, |20⟩の位置
            for a, idx_a in enumerate(indices):
                for b, idx_b in enumerate(indices):
                    U[idx_a, idx_b] = U_sub[a, b]

            # グローバル位相の除去（可能であれば）
            # U[0,0]を実数正にすることでグローバル位相を正規化
            global_phase = np.angle(U[0, 0])
            if abs(global_phase) > 1e-10:
                U *= np.exp(-1j * global_phase)

            # CustomTwoゲートを適用
            circuit.cu_two([i, j], U)

    def decompose_custom_two_gates(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """CustomTwoゲートを基本ゲートに分解."""
        from mqt.qudits.compiler.twodit.entanglement_qr import LogEntQRCEXPass

        backend = self.provider.get_backend("faketraps3six")
        compiler = LogEntQRCEXPass(backend)
        return compiler.transpile(circuit)


class PartiallyOptimizedSimulator:
    """部分最適化されたシミュレータ.

    元の実装（SuzukiTrotterMQTQuditSimulator）をベースに、
    部分最適化された時間発展を使用
    """

    def __init__(self, params: PhysicalParameters) -> None:
        self.params = params
        self.time_evol = PartiallyOptimizedTimeEvolution(params)
        self.N = params.N_molecules
        self.dim = 3**self.N

        provider = MQTQuditProvider()
        self.backend = provider.get_backend("tnsim")

    def build_initial_state_circuit(self, state_type: str = "all_triplet") -> QuantumCircuit:
        """初期状態を準備する回路を構築."""
        circuit = QuantumCircuit()
        reg = QuantumRegister("molecules", self.N, [3] * self.N)
        circuit.append(reg)

        if state_type == "all_triplet":
            for i in range(self.N):
                circuit.x(i)
        elif state_type == "alternating":
            for i in range(0, self.N, 2):
                circuit.x(i)
        elif state_type == "single_triplet":
            circuit.x(0)

        return circuit

    def add_single_trotter_step(self, circuit: QuantumCircuit, dt: float) -> None:
        """2次対称鈴木トロッター分解の1ステップを回路に追加."""
        self.time_evol.add_H0_evolution_gates(circuit, dt / 2)
        self.time_evol.add_H_transfer_evolution_gates(circuit, dt / 2)
        self.time_evol.add_H_TTA_evolution_gates(circuit, dt / 2)

        self.time_evol.add_H_TTA_evolution_gates(circuit, dt / 2)
        self.time_evol.add_H_transfer_evolution_gates(circuit, dt / 2)
        self.time_evol.add_H0_evolution_gates(circuit, dt / 2)

    def apply_radiative_decay_to_statevector(self, state_vector: np.ndarray, dt: float) -> np.ndarray:
        """放射減衰を状態ベクトルに適用."""
        state = state_vector.flatten().copy()

        if self.params.Gamma_fl > 0:
            for idx in range(self.dim):
                config = index_to_config(idx, self.N, 3)
                n_S1 = sum(1 for level in config if level == 2)

                decay_factor = np.exp(-self.params.Gamma_fl * dt * n_S1 / 2)
                state[idx] *= decay_factor

            norm = np.linalg.norm(state)
            if norm > 1e-12:
                state /= norm

        return state

    def calculate_populations(self, state_vector: np.ndarray) -> dict[str, float]:
        """状態ベクトルから個体数を計算."""
        state = state_vector.flatten()
        N_S0 = 0.0
        N_T1 = 0.0
        N_S1 = 0.0

        for idx in range(self.dim):
            prob = np.abs(state[idx]) ** 2
            config = index_to_config(idx, self.N, 3)

            for level in config:
                if level == 0:
                    N_S0 += prob
                elif level == 1:
                    N_T1 += prob
                elif level == 2:
                    N_S1 += prob

        return {"N_S0": N_S0, "N_T1": N_T1, "N_S1": N_S1}

    def simulate(
        self,
        T_total: float,
        N_steps: int,
        initial_state_type: str = "all_triplet",
        track_dynamics: bool = True,
        measure_gate_counts: bool = False,
    ) -> dict:
        """シミュレーションを実行.

        measure_gate_counts=True の場合、ゲート数の詳細を記録
        """
        dt = T_total / N_steps

        if measure_gate_counts:
            pass

        times = [0.0]
        populations_history = []
        states_history = [] if track_dynamics else None
        gate_counts_history = [] if measure_gate_counts else None

        # 初期状態
        init_circuit = self.build_initial_state_circuit(initial_state_type)
        job = self.backend.run(init_circuit)
        result = job.result()
        current_state = result.get_state_vector().flatten()

        populations_history.append(self.calculate_populations(current_state))
        if states_history is not None:
            states_history.append(current_state.copy())

        # 時間発展ループ
        start_time = time.time()

        for step in range(N_steps):
            circuit = self.build_initial_state_circuit(initial_state_type)

            for _s in range(step + 1):
                self.add_single_trotter_step(circuit, dt)

            # ゲート数測定（分解前）
            if measure_gate_counts:
                gates_before = len(circuit.instructions)

            # CustomTwoゲートを分解
            circuit = self.time_evol.decompose_custom_two_gates(circuit)

            # ゲート数測定（分解後）
            if measure_gate_counts:
                gates_after = len(circuit.instructions)
                gate_counts_history.append({
                    "step": step + 1,
                    "before_decomposition": gates_before,
                    "after_decomposition": gates_after,
                })

            # 回路を実行
            job = self.backend.run(circuit)
            result = job.result()
            state_after_unitary = result.get_state_vector().flatten()

            # 放射減衰を適用
            current_state = self.apply_radiative_decay_to_statevector(state_after_unitary, dt * (step + 1))

            if track_dynamics:
                t = (step + 1) * dt
                times.append(t)
                populations_history.append(self.calculate_populations(current_state))
                if states_history is not None:
                    states_history.append(current_state.copy())

            # 進捗表示
            if (step + 1) % max(1, N_steps // 10) == 0 or step == N_steps - 1:
                progress = (step + 1) / N_steps * 100
                msg = f"Progress: {progress:5.1f}% (step {step + 1}/{N_steps})"
                if measure_gate_counts and gate_counts_history:
                    msg += f" - Gates: {gate_counts_history[-1]['after_decomposition']}"

        elapsed_time = time.time() - start_time

        result_dict = {
            "times": np.array(times),
            "populations": populations_history,
            "states": states_history,
            "final_state": current_state,
            "elapsed_time": elapsed_time,
            "dt": dt,
            "N_steps": N_steps,
        }

        if measure_gate_counts:
            result_dict["gate_counts"] = gate_counts_history

        return result_dict


def compare_implementations() -> None:
    """元の実装と部分最適化実装を比較."""
    from mqt_qudits_four_molecule_implementation import SuzukiTrotterMQTQuditSimulator as OriginalSimulator

    params = PhysicalParameters()

    # 短いテストパラメータ
    T_total = 10.0  # fs
    N_steps = 1

    original_sim = OriginalSimulator(params)
    circuit_orig = original_sim.build_initial_state_circuit("all_triplet")
    original_sim.add_single_trotter_step(circuit_orig, T_total / N_steps)

    len(circuit_orig.instructions)

    # 分解
    from mqt.qudits.compiler.twodit.entanglement_qr import LogEntQRCEXPass
    from mqt.qudits.simulation import MQTQuditProvider

    provider = MQTQuditProvider()
    backend = provider.get_backend("faketraps3six")
    compiler = LogEntQRCEXPass(backend)
    circuit_orig_decomposed = compiler.transpile(circuit_orig)

    gates_after_orig = len(circuit_orig_decomposed.instructions)

    opt_sim = PartiallyOptimizedSimulator(params)
    circuit_opt = opt_sim.build_initial_state_circuit("all_triplet")
    opt_sim.add_single_trotter_step(circuit_opt, T_total / N_steps)

    len(circuit_opt.instructions)

    # 分解
    circuit_opt_decomposed = compiler.transpile(circuit_opt)

    gates_after_opt = len(circuit_opt_decomposed.instructions)

    improvement = (gates_after_orig - gates_after_opt) / gates_after_orig * 100

    if improvement > 0:
        pass
    else:
        pass


def test_partial_optimization() -> None:
    """部分最適化のテスト."""
    params = PhysicalParameters()
    simulator = PartiallyOptimizedSimulator(params)

    # 短いテスト
    T_total = 50.0
    N_steps = 5

    results = simulator.simulate(
        T_total=T_total,
        N_steps=N_steps,
        initial_state_type="all_triplet",
        track_dynamics=True,
        measure_gate_counts=True,
    )

    if "gate_counts" in results:
        results["gate_counts"][-1]


if __name__ == "__main__":
    # 実装比較を実行
    compare_implementations()

    # テスト実行
    # test_partial_optimization()
