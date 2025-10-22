#!/usr/bin/env python3
"""
疎構造認識実装のテスト

テスト項目:
1. ゲート数削減の検証
2. 忠実度の検証（fidelity = 1.0）
3. 疎構造検出の精度
4. 統計レポートの生成
"""

import sys
from pathlib import Path
import numpy as np
import pytest

# tutorials/からのインポート
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / 'tutorials'))
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / 'tools'))

from mqt_qudits_four_molecule_sparse_implementation import (
    PhysicalParameters,
    SparseAwareMQTQuditTimeEvolution,
    SparseAwareMQTGateGenerator
)


def test_gate_count_reduction():
    """ゲート数削減のテスト"""
    params = PhysicalParameters()
    gate_generator = SparseAwareMQTGateGenerator()
    
    # H_transfer型の2×2ユニタリ
    theta = 0.1
    U_transfer = np.eye(9, dtype=complex)
    U_transfer[1, 1] = np.cos(theta)
    U_transfer[1, 3] = -1j * np.sin(theta)
    U_transfer[3, 1] = -1j * np.sin(theta)
    U_transfer[3, 3] = np.cos(theta)
    
    # コンパイル
    gate_info = gate_generator.compile_unitary_to_gates(U_transfer, [0, 1])
    
    # 2×2部分空間として検出されるべき
    assert gate_info['structure_type'] == 'sparse_2x2', \
        f"Expected sparse_2x2 but got {gate_info['structure_type']}"
    
    # ゲート数は小さいはず（1-5ゲート程度）
    assert gate_info['gate_count'] <= 5, \
        f"Gate count too high: {gate_info['gate_count']}"
    
    print(f"✓ ゲート数削減テスト合格")
    print(f"  構造タイプ: {gate_info['structure_type']}")
    print(f"  ゲート数: {gate_info['gate_count']}")


def test_fidelity_preservation():
    """忠実度保存のテスト"""
    from integrated_sparse_compiler_v2 import IntegratedSparseCompilerV2
    
    # H_transfer型の2×2ユニタリ
    theta = 0.1
    U_transfer = np.eye(9, dtype=complex)
    U_transfer[1, 1] = np.cos(theta)
    U_transfer[1, 3] = -1j * np.sin(theta)
    U_transfer[3, 1] = -1j * np.sin(theta)
    U_transfer[3, 3] = np.cos(theta)
    
    # コンパイル
    compiler = IntegratedSparseCompilerV2()
    result = compiler.compile(U_transfer)
    
    # 忠実度チェック
    assert result.gate_sequence.fidelity > 0.9999, \
        f"忠実度が低すぎます: {result.gate_sequence.fidelity}"
    
    print(f"✓ 忠実度保存テスト合格")
    print(f"  忠実度: {result.gate_sequence.fidelity:.10f}")


def test_sparse_structure_detection():
    """疎構造検出のテスト"""
    from integrated_sparse_compiler_v2 import IntegratedSparseCompilerV2
    
    compiler = IntegratedSparseCompilerV2()
    
    # Test 1: 2×2部分空間
    theta = 0.1
    U_2x2 = np.eye(9, dtype=complex)
    U_2x2[1, 1] = np.cos(theta)
    U_2x2[1, 3] = -1j * np.sin(theta)
    U_2x2[3, 1] = -1j * np.sin(theta)
    U_2x2[3, 3] = np.cos(theta)
    
    result_2x2 = compiler.compile(U_2x2)
    assert result_2x2.structure_info.active_dimension == 2, \
        "2×2部分空間が検出されませんでした"
    
    # Test 2: 3×3部分空間
    J = 0.05
    H_sub = J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=complex)
    eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
    phases = np.exp(-1j * eigenvalues * 0.1)
    U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T
    
    U_3x3 = np.eye(9, dtype=complex)
    indices = [2, 4, 6]
    for a, idx_a in enumerate(indices):
        for b, idx_b in enumerate(indices):
            U_3x3[idx_a, idx_b] = U_sub[a, b]
    
    result_3x3 = compiler.compile(U_3x3)
    assert result_3x3.structure_info.active_dimension == 3, \
        "3×3部分空間が検出されませんでした"
    
    print(f"✓ 疎構造検出テスト合格")
    print(f"  2×2部分空間: 検出成功")
    print(f"  3×3部分空間: 検出成功")


def test_statistics_report():
    """統計レポート生成のテスト"""
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)
    
    # レポート生成
    report = time_evol.get_compilation_report()
    
    # レポートが空でないことを確認
    assert len(report) > 0
    assert "疎構造認識コンパイラ統計" in report
    
    print(f"✓ 統計レポート生成テスト合格")


def test_time_evolution_instantiation():
    """時間発展演算子のインスタンス化テスト"""
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)
    
    # 基本的なプロパティを確認
    assert time_evol.N == 4
    assert time_evol.dim == 81  # 3^4
    assert time_evol.gate_generator is not None
    
    print(f"✓ 時間発展演算子インスタンス化テスト合格")


def test_hamiltonian_construction():
    """ハミルトニアン構築のテスト"""
    from mqt_qudits_four_molecule_sparse_implementation import ExactDiagonalizationSolver
    
    params = PhysicalParameters()
    solver = ExactDiagonalizationSolver(params)
    
    # ハミルトニアン構築
    H = solver.build_hamiltonian()
    
    # エルミート性の確認
    assert np.allclose(H, H.conj().T), "ハミルトニアンがエルミートではありません"
    
    # サイズの確認
    assert H.shape == (81, 81), f"Expected (81, 81) but got {H.shape}"
    
    print(f"✓ ハミルトニアン構築テスト合格")
    print(f"  サイズ: {H.shape}")
    print(f"  エルミート性: OK")


if __name__ == "__main__":
    print("=" * 70)
    print("疎構造認識実装テストスイート")
    print("=" * 70)
    print()
    
    try:
        test_time_evolution_instantiation()
        print()
        test_hamiltonian_construction()
        print()
        test_fidelity_preservation()
        print()
        test_sparse_structure_detection()
        print()
        test_gate_count_reduction()
        print()
        test_statistics_report()
        print()
    except Exception as e:
        print(f"✗ テスト失敗: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("=" * 70)
    print("✓✓✓ すべてのテストに合格")
    print("=" * 70)
