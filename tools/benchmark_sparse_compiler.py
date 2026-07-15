#!/usr/bin/env python3
"""
疎構造認識コンパイラの性能ベンチマーク

測定項目:
1. ゲート数削減率
2. コンパイル時間
3. 忠実度
4. メモリ使用量
"""

import sys
from pathlib import Path
import numpy as np
import time
import tracemalloc

sys.path.insert(0, str(Path(__file__).parent.parent / 'tutorials'))
sys.path.insert(0, str(Path(__file__).parent))

from mqt_qudits_four_molecule_sparse_implementation import (
    PhysicalParameters,
    SparseAwareMQTQuditTimeEvolution
)


def benchmark_compilation_performance():
    """コンパイル性能のベンチマーク"""
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)
    
    # トレースメモリ開始
    tracemalloc.start()
    
    # タイミング開始
    start_time = time.time()
    
    # 1トロッターステップのゲート生成をシミュレート
    dt = 10.0  # fs
    
    # H_transfer（疎構造認識）
    print("=" * 70)
    print("H_transfer コンパイル中...")
    print("=" * 70)
    
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
        print(f"  ペア {i}-{j}: {gate_info['structure_type']}, "
              f"{gate_info['gate_count']} ゲート, "
              f"忠実度={gate_info['fidelity']:.10f}")
    
    print()
    print("=" * 70)
    print("H_TTA コンパイル中...")
    print("=" * 70)
    
    # H_TTA（疎構造認識）
    tta_gates = []
    for pair_idx, (i, j) in enumerate(params.neighbors):
        J = params.J[pair_idx]
        
        U = np.eye(9, dtype=complex)
        
        # 部分空間のハミルトニアン
        H_sub = J * np.array([
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 0]
        ], dtype=complex)
        
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
        print(f"  ペア {i}-{j}: {gate_info['structure_type']}, "
              f"{gate_info['gate_count']} ゲート, "
              f"忠実度={gate_info['fidelity']:.10f}")
    
    # タイミング終了
    end_time = time.time()
    compilation_time = end_time - start_time
    
    # メモリ使用量取得
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # 結果集計
    total_transfer_gates = sum(g['gate_count'] for g in transfer_gates)
    total_tta_gates = sum(g['gate_count'] for g in tta_gates)
    h0_gates = 8  # VirtRzゲート
    
    total_gates_sparse = h0_gates + total_transfer_gates + total_tta_gates
    total_gates_traditional = 8 + 6 * 1000  # H0: 8, CustomTwo×6: 6000
    reduction_rate = (1 - total_gates_sparse / total_gates_traditional) * 100
    
    # 比較表の作成
    print()
    print("=" * 70)
    print("ベンチマーク結果")
    print("=" * 70)
    print()
    print("【ゲート数比較】")
    print(f"  H0: {h0_gates} ゲート")
    print(f"  H_transfer: {total_transfer_gates} ゲート (3ペア)")
    print(f"  H_TTA: {total_tta_gates} ゲート (3ペア)")
    print(f"  合計（疎構造認識）: {total_gates_sparse} ゲート")
    print(f"  合計（従来方式推定）: {total_gates_traditional} ゲート")
    print(f"  削減率: {reduction_rate:.1f}%")
    print()
    print("【性能メトリクス】")
    print(f"  コンパイル時間: {compilation_time*1000:.2f} ms")
    print(f"  ピークメモリ使用量: {peak / 1024 / 1024:.2f} MB")
    print()
    print("【Qubit vs Qudit 比較】")
    print()
    print("| 実装方式 | Qubit/Qudit数 | ゲート数/ステップ | 20ステップ総数 | 備考 |")
    print("|---------|--------------|-----------------|---------------|------|")
    print("| Qubit | 8 qubits | 112 | 2,240 | 標準的な実装 |")
    print(f"| Qudit（従来） | 4 qutrits | ~6,182 | ~123,640 | LogEntQRCEXPass |")
    print(f"| **Qudit（改良）** | 4 qutrits | **{total_gates_sparse}** | **{total_gates_sparse*20}** | **疎構造認識** |")
    print()
    
    qubit_gates = 112
    improvement_over_qubit = qubit_gates / total_gates_sparse
    traditional_ratio = 6182 / qubit_gates
    sparse_ratio = total_gates_sparse / qubit_gates
    
    print("【性能向上】")
    print(f"  従来Qudit vs Qubit: {traditional_ratio:.1f}倍遅い")
    print(f"  改良Qudit vs Qubit: {improvement_over_qubit:.1f}倍高速")
    print(f"  改良による改善: {traditional_ratio / sparse_ratio:.0f}倍の性能向上")
    print()
    print("=" * 70)
    print(time_evol.get_compilation_report())
    print("=" * 70)


def benchmark_exact_diagonalization():
    """厳密対角化のベンチマーク"""
    from mqt_qudits_four_molecule_sparse_implementation import ExactDiagonalizationSolver
    
    print()
    print("=" * 70)
    print("厳密対角化ベンチマーク")
    print("=" * 70)
    
    params = PhysicalParameters()
    solver = ExactDiagonalizationSolver(params)
    
    # ハミルトニアン構築
    start_time = time.time()
    H = solver.build_hamiltonian()
    build_time = time.time() - start_time
    
    # 対角化
    start_time = time.time()
    solver.diagonalize()
    diag_time = time.time() - start_time
    
    # 時間発展
    initial_state = np.zeros(81, dtype=complex)
    initial_state[40] = 1.0  # |1,1,0,0⟩
    
    start_time = time.time()
    final_state = solver.time_evolve(initial_state, 200.0)
    evolve_time = time.time() - start_time
    
    print()
    print(f"ハミルトニアンサイズ: {H.shape}")
    print(f"ハミルトニアン構築: {build_time*1000:.2f} ms")
    print(f"対角化: {diag_time*1000:.2f} ms")
    print(f"時間発展: {evolve_time*1000:.2f} ms")
    print()
    print(f"固有値範囲: [{solver.eigenvalues[0]:.4f}, {solver.eigenvalues[-1]:.4f}] eV")
    print(f"エネルギーギャップ: {solver.eigenvalues[1] - solver.eigenvalues[0]:.6f} eV")


if __name__ == "__main__":
    print("=" * 70)
    print("疎構造認識コンパイラ 性能ベンチマーク")
    print("=" * 70)
    print()
    
    try:
        benchmark_compilation_performance()
        benchmark_exact_diagonalization()
        
        print()
        print("=" * 70)
        print("✓✓✓ ベンチマーク完了")
        print("=" * 70)
    except Exception as e:
        print(f"✗ ベンチマーク失敗: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
