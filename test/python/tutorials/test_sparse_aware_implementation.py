#!/usr/bin/env python3
"""
疎構造認識実装のテスト

テスト項目:
1. ゲート数削減の検証
2. 忠実度の検証（fidelity = 1.0）
3. 疎構造検出の精度
4. 統計レポートの生成
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# tutorials/ ディレクトリへのパスを追加
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "tutorials"))
# tools/ ディレクトリへのパスを追加
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "tools"))

from integrated_sparse_compiler_v2 import IntegratedSparseCompilerV2
from mqt_qudits_four_molecule_sparse_implementation import (
    PhysicalParameters,
    SparseAwareMQTQuditTimeEvolution,
)


def test_gate_count_reduction():
    """ゲート数削減のテスト"""
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)

    # 実際にゲートを生成して統計を収集
    dt = 10.0

    # H_transfer型のユニタリをコンパイル
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

        time_evol.gate_generator.compile_unitary_to_gates(U, [i, j])

    # H_TTA型のユニタリをコンパイル
    for pair_idx, (i, j) in enumerate(params.neighbors):
        J = params.J[pair_idx]

        U = np.eye(9, dtype=complex)
        H_sub = J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=complex)
        eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
        phases = np.exp(-1j * eigenvalues * dt / params.hbar)
        U_sub = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T

        indices = [2, 4, 6]
        for a, idx_a in enumerate(indices):
            for b, idx_b in enumerate(indices):
                U[idx_a, idx_b] = U_sub[a, b]

        time_evol.gate_generator.compile_unitary_to_gates(U, [i, j])

    # 統計情報を確認
    stats = time_evol.gate_generator.compilation_stats

    # 少なくとも疎構造が検出されるべき
    assert stats["sparse_2x2"] > 0 or stats["sparse_3x3"] > 0, "疎構造が検出されませんでした"

    # 密構造として扱われてはいけない
    assert stats["dense"] == 0, "密構造として扱われてはいけません"

    print("✓ ゲート数削減テスト合格")
    print(f"  検出された2×2部分空間: {stats['sparse_2x2']}")
    print(f"  検出された3×3部分空間: {stats['sparse_3x3']}")


def test_fidelity_preservation():
    """忠実度保存のテスト"""
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
    assert result.gate_sequence.fidelity > 0.9999, f"忠実度が低すぎます: {result.gate_sequence.fidelity}"

    print("✓ 忠実度保存テスト合格")
    print(f"  忠実度: {result.gate_sequence.fidelity:.10f}")


def test_sparse_structure_detection():
    """疎構造検出のテスト"""
    compiler = IntegratedSparseCompilerV2()

    # Test 1: 2×2部分空間
    theta = 0.1
    U_2x2 = np.eye(9, dtype=complex)
    U_2x2[1, 1] = np.cos(theta)
    U_2x2[1, 3] = -1j * np.sin(theta)
    U_2x2[3, 1] = -1j * np.sin(theta)
    U_2x2[3, 3] = np.cos(theta)

    result_2x2 = compiler.compile(U_2x2)
    assert result_2x2.structure_info.active_dimension == 2, "2×2部分空間が検出されませんでした"

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
    assert result_3x3.structure_info.active_dimension == 3, "3×3部分空間が検出されませんでした"

    print("✓ 疎構造検出テスト合格")
    print("  2×2部分空間: 検出成功")
    print("  3×3部分空間: 検出成功")


def test_statistics_report():
    """統計レポート生成のテスト"""
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)

    # レポート生成
    report = time_evol.get_compilation_report()

    # レポートが空でないことを確認
    assert len(report) > 0, "レポートが空です"
    assert "疎構造認識コンパイラ統計" in report, "レポートに必要な情報が含まれていません"

    print("✓ 統計レポート生成テスト合格")


if __name__ == "__main__":
    print("=" * 70)
    print("疎構造認識実装テストスイート")
    print("=" * 70)
    print()

    test_fidelity_preservation()
    test_sparse_structure_detection()
    test_gate_count_reduction()
    test_statistics_report()

    print()
    print("=" * 70)
    print("✓✓✓ すべてのテストに合格")
    print("=" * 70)
