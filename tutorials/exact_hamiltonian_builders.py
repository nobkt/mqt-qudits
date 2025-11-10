#!/usr/bin/env python3
"""
Exact Hamiltonian Matrix Builders

This module provides functions to build exact unitary matrices for H_transfer and H_TTA
time evolution, which can then be compiled into quantum gates.

No approximations or heuristics are used - all implementations are mathematically exact.
"""

import numpy as np
import scipy.linalg
from typing import Tuple


def build_H_transfer_matrix(V: float, dim: int = 3) -> np.ndarray:
    """
    Build the exact H_transfer Hamiltonian matrix for a pair of qudits.
    
    For two d-level systems (default d=3 for qutrits):
    H_transfer = V (|01⟩⟨10| + |10⟩⟨01|)
    
    This operates only on the 2×2 subspace {|01⟩, |10⟩}.
    
    Args:
        V: Energy transfer coupling strength (eV)
        dim: Dimension of each qudit (default 3 for qutrits)
    
    Returns:
        H_transfer: (dim²×dim²) Hermitian matrix
    """
    total_dim = dim * dim
    H = np.zeros((total_dim, total_dim), dtype=complex)
    
    # |01⟩ has index 0*dim + 1 = 1
    # |10⟩ has index 1*dim + 0 = dim
    idx_01 = 0 * dim + 1  # = 1 for dim=3
    idx_10 = 1 * dim + 0  # = 3 for dim=3
    
    # H_transfer = V(|01⟩⟨10| + |10⟩⟨01|)
    H[idx_01, idx_10] = V
    H[idx_10, idx_01] = V
    
    return H


def build_H_TTA_matrix(J: float, dim: int = 3) -> np.ndarray:
    """
    Build the exact H_TTA Hamiltonian matrix for a pair of qudits.
    
    For two d-level systems (default d=3 for qutrits):
    H_TTA = J (|02⟩⟨11| + |11⟩⟨02| + |20⟩⟨11| + |11⟩⟨20|)
    
    This operates only on the 3×3 subspace {|02⟩, |11⟩, |20⟩}.
    
    Args:
        J: TTA coupling strength (eV)
        dim: Dimension of each qudit (default 3 for qutrits)
    
    Returns:
        H_TTA: (dim²×dim²) Hermitian matrix
    """
    total_dim = dim * dim
    H = np.zeros((total_dim, total_dim), dtype=complex)
    
    # State indices in |ij⟩ notation (i,j ∈ {0,1,2} for qutrits)
    # Index = i*dim + j
    idx_02 = 0 * dim + 2  # = 2 for dim=3
    idx_11 = 1 * dim + 1  # = 4 for dim=3
    idx_20 = 2 * dim + 0  # = 6 for dim=3
    
    # H_TTA = J(|02⟩⟨11| + |11⟩⟨02| + |20⟩⟨11| + |11⟩⟨20|)
    H[idx_02, idx_11] = J
    H[idx_11, idx_02] = J
    H[idx_20, idx_11] = J
    H[idx_11, idx_20] = J
    
    return H


def build_time_evolution_unitary(H: np.ndarray, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """
    Build exact time evolution unitary from Hamiltonian.
    
    U(t) = exp(-i * H * t / ℏ)
    
    Args:
        H: Hermitian Hamiltonian matrix
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    
    Returns:
        U: Unitary time evolution operator
    """
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    return U


def build_H_transfer_unitary(V: float, dt: float, hbar: float = 0.6582119569, dim: int = 3) -> np.ndarray:
    """
    Build exact time evolution unitary for H_transfer.
    
    Args:
        V: Energy transfer coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
        dim: Qudit dimension (default 3)
    
    Returns:
        U_transfer: (dim²×dim²) unitary matrix
    """
    H = build_H_transfer_matrix(V, dim)
    U = build_time_evolution_unitary(H, dt, hbar)
    return U


def build_H_TTA_unitary(J: float, dt: float, hbar: float = 0.6582119569, dim: int = 3) -> np.ndarray:
    """
    Build exact time evolution unitary for H_TTA.
    
    Args:
        J: TTA coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
        dim: Qudit dimension (default 3)
    
    Returns:
        U_TTA: (dim²×dim²) unitary matrix
    """
    H = build_H_TTA_matrix(J, dim)
    U = build_time_evolution_unitary(H, dt, hbar)
    return U


def verify_sparse_structure(U: np.ndarray, tolerance: float = 1e-10) -> dict:
    """
    Verify and analyze the sparse structure of a unitary matrix.
    
    Args:
        U: Unitary matrix to analyze
        tolerance: Threshold for considering elements as zero
    
    Returns:
        Dictionary with structure information:
        - active_indices: List of indices with non-identity behavior
        - active_dimension: Size of active subspace
        - is_sparse: Whether matrix has significant sparse structure
        - identity_fraction: Fraction of rows/cols that are identity
    """
    dim = U.shape[0]
    I = np.eye(dim, dtype=complex)
    
    # Find rows/columns that differ from identity
    active_indices = []
    for i in range(dim):
        # Check if row i differs from identity row
        row_diff = np.linalg.norm(U[i, :] - I[i, :])
        # Check if column i differs from identity column
        col_diff = np.linalg.norm(U[:, i] - I[:, i])
        
        if row_diff > tolerance or col_diff > tolerance:
            active_indices.append(i)
    
    active_dim = len(active_indices)
    is_sparse = active_dim < dim / 2
    identity_fraction = 1.0 - (active_dim / dim)
    
    return {
        'active_indices': active_indices,
        'active_dimension': active_dim,
        'is_sparse': is_sparse,
        'identity_fraction': identity_fraction
    }


def extract_active_subspace_unitary(U: np.ndarray, active_indices: list) -> np.ndarray:
    """
    Extract the active subspace unitary from a larger sparse unitary.
    
    Args:
        U: Full unitary matrix
        active_indices: List of active indices
    
    Returns:
        U_sub: Reduced unitary acting only on active subspace
    """
    active_dim = len(active_indices)
    U_sub = np.zeros((active_dim, active_dim), dtype=complex)
    
    for i, idx_i in enumerate(active_indices):
        for j, idx_j in enumerate(active_indices):
            U_sub[i, j] = U[idx_i, idx_j]
    
    return U_sub


# ===================================================================
# Validation and Testing
# ===================================================================

def test_hamiltonians():
    """Test Hamiltonian construction"""
    print("="*70)
    print("Testing Exact Hamiltonian Builders")
    print("="*70)
    
    # Test H_transfer
    print("\n1. H_transfer Matrix")
    V = 0.1
    H_tr = build_H_transfer_matrix(V, dim=3)
    print(f"   V = {V} eV")
    print(f"   Matrix shape: {H_tr.shape}")
    print(f"   Non-zero elements: {np.count_nonzero(H_tr)}")
    print(f"   Is Hermitian: {np.allclose(H_tr, H_tr.conj().T)}")
    
    # Check structure
    expected_nnz = 2  # Only H[1,3] and H[3,1]
    actual_nnz = np.count_nonzero(H_tr)
    print(f"   Expected non-zeros: {expected_nnz}, Actual: {actual_nnz}")
    assert actual_nnz == expected_nnz, "H_transfer structure incorrect"
    print("   ✓ Structure correct")
    
    # Test H_TTA
    print("\n2. H_TTA Matrix")
    J = 0.05
    H_TTA = build_H_TTA_matrix(J, dim=3)
    print(f"   J = {J} eV")
    print(f"   Matrix shape: {H_TTA.shape}")
    print(f"   Non-zero elements: {np.count_nonzero(H_TTA)}")
    print(f"   Is Hermitian: {np.allclose(H_TTA, H_TTA.conj().T)}")
    
    # Check structure
    expected_nnz = 4  # H[2,4], H[4,2], H[6,4], H[4,6]
    actual_nnz = np.count_nonzero(H_TTA)
    print(f"   Expected non-zeros: {expected_nnz}, Actual: {actual_nnz}")
    assert actual_nnz == expected_nnz, "H_TTA structure incorrect"
    print("   ✓ Structure correct")
    
    # Test time evolution unitaries
    print("\n3. Time Evolution Unitaries")
    dt = 5.0
    hbar = 0.6582119569
    
    U_tr = build_H_transfer_unitary(V, dt, hbar, dim=3)
    print(f"   U_transfer shape: {U_tr.shape}")
    print(f"   Is unitary: {np.allclose(U_tr @ U_tr.conj().T, np.eye(9))}")
    
    # Verify sparse structure
    struct_tr = verify_sparse_structure(U_tr)
    print(f"   Active indices: {struct_tr['active_indices']}")
    print(f"   Active dimension: {struct_tr['active_dimension']}")
    print(f"   Is sparse: {struct_tr['is_sparse']}")
    print(f"   Identity fraction: {struct_tr['identity_fraction']:.2%}")
    
    U_TTA = build_H_TTA_unitary(J, dt, hbar, dim=3)
    print(f"\n   U_TTA shape: {U_TTA.shape}")
    print(f"   Is unitary: {np.allclose(U_TTA @ U_TTA.conj().T, np.eye(9))}")
    
    # Verify sparse structure
    struct_TTA = verify_sparse_structure(U_TTA)
    print(f"   Active indices: {struct_TTA['active_indices']}")
    print(f"   Active dimension: {struct_TTA['active_dimension']}")
    print(f"   Is sparse: {struct_TTA['is_sparse']}")
    print(f"   Identity fraction: {struct_TTA['identity_fraction']:.2%}")
    
    print("\n" + "="*70)
    print("All tests passed ✓")
    print("="*70)


if __name__ == "__main__":
    test_hamiltonians()
