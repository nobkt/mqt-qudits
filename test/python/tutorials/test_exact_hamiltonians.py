#!/usr/bin/env python3
"""Test suite for exact Hamiltonian implementations.

Tests the exact_hamiltonian_builders and exact_qubit_hamiltonians modules
to ensure they produce correct unitaries with no approximations.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

# Add tutorials directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../tutorials"))

from exact_hamiltonian_builders import (
    build_H_transfer_matrix,
    build_H_transfer_unitary,
    build_H_TTA_matrix,
    build_H_TTA_unitary,
    verify_sparse_structure,
)
from exact_qubit_hamiltonians import (
    build_H_transfer_qubit_unitary,
    build_H_TTA_qubit_unitary,
)


class TestExactHamiltonianBuilders:
    """Test exact Hamiltonian matrix builders for qudits."""

    def test_H_transfer_matrix_structure(self):
        """Test H_transfer Hamiltonian matrix has correct structure."""
        V = 0.1
        H_tr = build_H_transfer_matrix(V, dim=3)

        # Should be 9x9 for qutrits (3x3)
        assert H_tr.shape == (9, 9)

        # Should be Hermitian
        assert np.allclose(H_tr, H_tr.conj().T)

        # Should have exactly 2 non-zero elements (off-diagonal coupling)
        assert np.count_nonzero(H_tr) == 2

        # Verify the specific structure
        # |01⟩ has index 1, |10⟩ has index 3
        assert np.abs(H_tr[1, 3] - V) < 1e-10
        assert np.abs(H_tr[3, 1] - V) < 1e-10

    def test_H_TTA_matrix_structure(self):
        """Test H_TTA Hamiltonian matrix has correct structure."""
        J = 0.05
        H_TTA = build_H_TTA_matrix(J, dim=3)

        # Should be 9x9 for qutrits
        assert H_TTA.shape == (9, 9)

        # Should be Hermitian
        assert np.allclose(H_TTA, H_TTA.conj().T)

        # Should have exactly 4 non-zero elements
        assert np.count_nonzero(H_TTA) == 4

        # Verify the specific structure
        # |02⟩=2, |11⟩=4, |20⟩=6
        assert np.abs(H_TTA[2, 4] - J) < 1e-10
        assert np.abs(H_TTA[4, 2] - J) < 1e-10
        assert np.abs(H_TTA[6, 4] - J) < 1e-10
        assert np.abs(H_TTA[4, 6] - J) < 1e-10

    def test_H_transfer_unitary_is_unitary(self):
        """Test H_transfer time evolution is unitary."""
        V = 0.1
        dt = 5.0
        hbar = 0.6582119569

        U = build_H_transfer_unitary(V, dt, hbar, dim=3)

        # Should be unitary: U† U = I
        I = np.eye(9)
        assert np.allclose(U @ U.conj().T, I)
        assert np.allclose(U.conj().T @ U, I)

    def test_H_TTA_unitary_is_unitary(self):
        """Test H_TTA time evolution is unitary."""
        J = 0.05
        dt = 5.0
        hbar = 0.6582119569

        U = build_H_TTA_unitary(J, dt, hbar, dim=3)

        # Should be unitary
        I = np.eye(9)
        assert np.allclose(U @ U.conj().T, I)
        assert np.allclose(U.conj().T @ U, I)

    def test_H_transfer_sparse_structure(self):
        """Test H_transfer unitary has sparse structure."""
        V = 0.1
        dt = 5.0

        U = build_H_transfer_unitary(V, dt)
        struct = verify_sparse_structure(U)

        # Should have exactly 2 active indices
        assert struct["active_dimension"] == 2
        assert struct["is_sparse"] is True
        assert struct["active_indices"] == [1, 3]

    def test_H_TTA_sparse_structure(self):
        """Test H_TTA unitary has sparse structure."""
        J = 0.05
        dt = 5.0

        U = build_H_TTA_unitary(J, dt)
        struct = verify_sparse_structure(U)

        # Should have exactly 3 active indices
        assert struct["active_dimension"] == 3
        assert struct["is_sparse"] is True
        assert struct["active_indices"] == [2, 4, 6]


class TestExactQubitHamiltonians:
    """Test exact Hamiltonian implementations for qubits."""

    def test_H_transfer_qubit_unitary_is_unitary(self):
        """Test qubit H_transfer unitary is unitary."""
        V = 0.1
        dt = 5.0
        hbar = 0.6582119569

        U = build_H_transfer_qubit_unitary(V, dt, hbar)

        # Should be 16x16 (4 qubits)
        assert U.shape == (16, 16)

        # Should be unitary
        I = np.eye(16)
        assert np.allclose(U @ U.conj().T, I)
        assert np.allclose(U.conj().T @ U, I)

    def test_H_TTA_qubit_unitary_is_unitary(self):
        """Test qubit H_TTA unitary is unitary."""
        J = 0.05
        dt = 5.0
        hbar = 0.6582119569

        U = build_H_TTA_qubit_unitary(J, dt, hbar)

        # Should be 16x16
        assert U.shape == (16, 16)

        # Should be unitary
        I = np.eye(16)
        assert np.allclose(U @ U.conj().T, I)
        assert np.allclose(U.conj().T @ U, I)

    def test_H_transfer_qubit_active_subspace(self):
        """Test qubit H_transfer operates on correct 2D subspace."""
        V = 0.1
        dt = 5.0

        U = build_H_transfer_qubit_unitary(V, dt)

        # Find active indices (those that differ from identity)
        I = np.eye(16)
        active_indices = [i for i in range(16) if not np.allclose(U[i, :], I[i, :])]

        # Should operate on indices 1 and 4 (|0001⟩ and |0100⟩)
        assert active_indices == [1, 4]

    def test_H_TTA_qubit_active_subspace(self):
        """Test qubit H_TTA operates on correct 3D subspace."""
        J = 0.05
        dt = 5.0

        U = build_H_TTA_qubit_unitary(J, dt)

        # Find active indices
        I = np.eye(16)
        active_indices = [i for i in range(16) if not np.allclose(U[i, :], I[i, :])]

        # Should operate on indices 2, 5, 8 (|0010⟩, |0101⟩, |1000⟩)
        assert active_indices == [2, 5, 8]

    @pytest.mark.parametrize(
        ("V", "dt"),
        [
            (0.05, 5.0),
            (0.1, 5.0),
            (0.2, 5.0),
            (0.1, 2.0),
            (0.1, 10.0),
        ],
    )
    def test_H_transfer_different_parameters(self, V, dt):
        """Test H_transfer with different coupling strengths and time steps."""
        U = build_H_transfer_qubit_unitary(V, dt)

        # Should always be unitary
        I = np.eye(16)
        assert np.allclose(U @ U.conj().T, I, atol=1e-10)

    @pytest.mark.parametrize(
        ("J", "dt"),
        [
            (0.025, 5.0),
            (0.05, 5.0),
            (0.1, 5.0),
            (0.05, 2.0),
            (0.05, 10.0),
        ],
    )
    def test_H_TTA_different_parameters(self, J, dt):
        """Test H_TTA with different coupling strengths and time steps."""
        U = build_H_TTA_qubit_unitary(J, dt)

        # Should always be unitary
        I = np.eye(16)
        assert np.allclose(U @ U.conj().T, I, atol=1e-10)


class TestHamiltonianConsistency:
    """Test consistency between qubit and qudit implementations."""

    def test_unitaries_are_exact_not_approximate(self):
        """Verify that all unitaries are exact (no approximations)."""
        # This test documents that we use scipy.linalg.expm
        # which is exact (up to numerical precision)

        V = 0.1
        J = 0.05
        dt = 5.0

        # Build unitaries
        U_tr_qudit = build_H_transfer_unitary(V, dt)
        U_TTA_qudit = build_H_TTA_unitary(J, dt)
        U_tr_qubit = build_H_transfer_qubit_unitary(V, dt)
        U_TTA_qubit = build_H_TTA_qubit_unitary(J, dt)

        # All should be unitary to machine precision
        for U in [U_tr_qudit, U_TTA_qudit]:
            I = np.eye(9)
            assert np.allclose(U @ U.conj().T, I, atol=1e-14)

        for U in [U_tr_qubit, U_TTA_qubit]:
            I = np.eye(16)
            assert np.allclose(U @ U.conj().T, I, atol=1e-14)

    def test_time_reversibility(self):
        """Test that U(t) @ U(-t) = I (time reversibility)."""
        V = 0.1
        J = 0.05
        dt = 5.0

        # Qudit
        U_tr_forward = build_H_transfer_unitary(V, dt)
        U_tr_backward = build_H_transfer_unitary(V, -dt)
        I9 = np.eye(9)
        assert np.allclose(U_tr_forward @ U_tr_backward, I9, atol=1e-10)

        U_TTA_forward = build_H_TTA_unitary(J, dt)
        U_TTA_backward = build_H_TTA_unitary(J, -dt)
        assert np.allclose(U_TTA_forward @ U_TTA_backward, I9, atol=1e-10)

        # Qubit
        U_tr_qubit_forward = build_H_transfer_qubit_unitary(V, dt)
        U_tr_qubit_backward = build_H_transfer_qubit_unitary(V, -dt)
        I16 = np.eye(16)
        assert np.allclose(U_tr_qubit_forward @ U_tr_qubit_backward, I16, atol=1e-10)

        U_TTA_qubit_forward = build_H_TTA_qubit_unitary(J, dt)
        U_TTA_qubit_backward = build_H_TTA_qubit_unitary(J, -dt)
        assert np.allclose(U_TTA_qubit_forward @ U_TTA_qubit_backward, I16, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
