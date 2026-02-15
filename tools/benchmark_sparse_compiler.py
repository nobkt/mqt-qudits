#!/usr/bin/env python3
"""疎構造認識コンパイラの性能ベンチマーク.

測定項目:
1. ゲート数削減率
2. コンパイル時間
3. 忠実度
4. メモリ使用量
"""

from __future__ import annotations

import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "tutorials"))
sys.path.insert(0, str(Path(__file__).parent))

from mqt_qudits_four_molecule_sparse_implementation import PhysicalParameters, SparseAwareMQTQuditTimeEvolution


def benchmark_compilation_performance() -> None:
    """コンパイル性能のベンチマーク."""
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)

    # トレースメモリ開始
    tracemalloc.start()

    # タイミング開始
    start_time = time.time()

    # 1トロッターステップのゲート生成をシミュレート
    dt = 10.0  # fs

    # H_transfer（疎構造認識）

    transfer_gates = []
    for pair_idx, (i, j) in enumerate(params.neighbors):
        V = params.V[pair_idx]
        theta = V * dt / params.hbar

        U = np.eye(9, dtype=complex)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        U[1, 1] = cos_theta
        U[1, 3] = -1j * sin_theta
        U[3, 1] = -1j * sin_theta
        U[3, 3] = cos_theta

        gate_info = time_evol.gate_generator.compile_unitary_to_gates(U, [i, j])
        transfer_gates.append(gate_info)

    # H_TTA（疎構造認識）
    tta_gates = []
    for pair_idx, (i, j) in enumerate(params.neighbors):
        J = params.J[pair_idx]

        U = np.eye(9, dtype=complex)

        # 部分空間のハミルトニアン
        H_sub = J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=complex)

        # 固有値分解
        eigenvalues, eigenvectors = np.linalg.eigh(H_sub)

        # 時間発展演算子
        phases = np.exp(-1j * eigenvalues * dt / params.hbar)
        U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T

        # 9×9行列の該当部分に埋め込む
        indices = [2, 4, 6]
        for a, idx_a in enumerate(indices):
            for b, idx_b in enumerate(indices):
                U[idx_a, idx_b] = U_sub[a, b]

        gate_info = time_evol.gate_generator.compile_unitary_to_gates(U, [i, j])
        tta_gates.append(gate_info)

    # タイミング終了
    end_time = time.time()
    end_time - start_time

    # メモリ使用量取得
    _current, _peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # 結果集計
    total_transfer_gates = sum(g["gate_count"] for g in transfer_gates)
    total_tta_gates = sum(g["gate_count"] for g in tta_gates)
    h0_gates = 8  # VirtRzゲート

    total_gates_sparse = h0_gates + total_transfer_gates + total_tta_gates
    total_gates_traditional = 8 + 6 * 1000  # H0: 8, CustomTwo×6: 6000
    (1 - total_gates_sparse / total_gates_traditional) * 100

    # 比較表の作成

    qubit_gates = 112
    qubit_gates / total_gates_sparse
    6182 / qubit_gates
    total_gates_sparse / qubit_gates


def benchmark_exact_diagonalization() -> None:
    """厳密対角化のベンチマーク."""
    from mqt_qudits_four_molecule_sparse_implementation import ExactDiagonalizationSolver

    params = PhysicalParameters()
    solver = ExactDiagonalizationSolver(params)

    # ハミルトニアン構築
    start_time = time.time()
    solver.build_hamiltonian()
    time.time() - start_time

    # 対角化
    start_time = time.time()
    solver.diagonalize()
    time.time() - start_time

    # 時間発展
    initial_state = np.zeros(81, dtype=complex)
    initial_state[40] = 1.0  # |1,1,0,0⟩

    start_time = time.time()
    solver.time_evolve(initial_state, 200.0)
    time.time() - start_time


if __name__ == "__main__":
    try:
        benchmark_compilation_performance()
        benchmark_exact_diagonalization()

    except Exception:
        import traceback

        traceback.print_exc()
        sys.exit(1)
