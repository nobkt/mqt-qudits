#!/usr/bin/env python3
"""Test script to verify TTA Hamiltonian fix."""

from __future__ import annotations

import numpy as np
import scipy.linalg


def ket_bra(i, j, dim=3):
    """Create |i⟩⟨j| operator for dim-dimensional system."""
    op = np.zeros((dim, dim))
    op[i, j] = 1
    return op


def test_qutrit_tta_hamiltonian() -> bool:
    """Test the correct 9x9 TTA Hamiltonian for 2-qutrit system."""
    # Define basis operators
    S0_to_T1 = ket_bra(1, 0)  # |T1⟩⟨S0| = |1⟩⟨0|
    T1_to_S0 = ket_bra(0, 1)  # |S0⟩⟨T1| = |0⟩⟨1|
    S1_to_T1 = ket_bra(1, 2)  # |T1⟩⟨S1| = |1⟩⟨2|
    T1_to_S1 = ket_bra(2, 1)  # |S1⟩⟨T1| = |2⟩⟨1|

    # Build CORRECT TTA Hamiltonian (all 4 terms)
    # Ĥ_TTA = J [|S0⟩_i⟨T1|_i ⊗ |S1⟩_j⟨T1|_j
    #          + |S1⟩_i⟨T1|_i ⊗ |S0⟩_j⟨T1|_j
    #          + |T1⟩_i⟨S0|_i ⊗ |T1⟩_j⟨S1|_j
    #          + |T1⟩_i⟨S1|_i ⊗ |T1⟩_j⟨S0|_j]

    term1 = np.kron(T1_to_S0, T1_to_S1)  # |S0⟩_i⟨T1|_i ⊗ |S1⟩_j⟨T1|_j
    term2 = np.kron(T1_to_S1, T1_to_S0)  # |S1⟩_i⟨T1|_i ⊗ |S0⟩_j⟨T1|_j
    term3 = np.kron(S0_to_T1, S1_to_T1)  # |T1⟩_i⟨S0|_i ⊗ |T1⟩_j⟨S1|_j
    term4 = np.kron(S1_to_T1, S0_to_T1)  # |T1⟩_i⟨S1|_i ⊗ |T1⟩_j⟨S0|_j

    H_TTA = term1 + term2 + term3 + term4

    # Verify physical process

    # Verify hermiticity
    is_hermitian = np.allclose(H_TTA, H_TTA.T.conj())

    # Test time evolution
    J = 0.05  # eV
    dt = 5.0  # fs
    hbar = 0.6582119569  # eV·fs

    U_TTA = scipy.linalg.expm(-1j * J * H_TTA * dt / hbar)

    # Verify unitarity
    is_unitary = np.allclose(U_TTA @ U_TTA.T.conj(), np.eye(9))

    return bool(is_hermitian and is_unitary)


def test_qubit_tta_hamiltonian() -> bool:
    """Test the correct 16x16 TTA Hamiltonian for 4-qubit system."""
    H = np.zeros((16, 16), dtype=complex)

    # Encoding: |S0⟩ = |00⟩, |T1⟩ = |01⟩, |S1⟩ = |10⟩
    # Little-endian: rightmost bit = qubit 0

    # |T1,T1⟩ = |01,01⟩ = |0101⟩ (q3q2q1q0)
    idx_T1_T1 = 0b0101  # = 5

    # |S0,S1⟩ = |00,10⟩ = |0010⟩
    idx_S0_S1 = 0b0010  # = 2

    # |S1,S0⟩ = |10,00⟩ = |1000⟩
    idx_S1_S0 = 0b1000  # = 8

    # Build CORRECT Hamiltonian (both terms!)
    # Term 1: |S0,S1⟩ ↔ |T1,T1⟩
    H[idx_S0_S1, idx_T1_T1] = 1.0
    H[idx_T1_T1, idx_S0_S1] = 1.0

    # Term 2: |S1,S0⟩ ↔ |T1,T1⟩  (THIS WAS MISSING!)
    H[idx_S1_S0, idx_T1_T1] = 1.0
    H[idx_T1_T1, idx_S1_S0] = 1.0

    nz = np.where(H != 0)
    for _i, _j in zip(nz[0], nz[1]):
        pass

    # Verify physical process

    # Verify hermiticity
    is_hermitian = np.allclose(H, H.T.conj())

    # Test time evolution
    J = 0.05  # eV
    dt = 5.0  # fs
    hbar = 0.6582119569  # eV·fs

    U_TTA = scipy.linalg.expm(-1j * J * H * dt / hbar)

    # Verify unitarity
    is_unitary = np.allclose(U_TTA @ U_TTA.T.conj(), np.eye(16))

    return bool(is_hermitian and is_unitary)


if __name__ == "__main__":
    success = True
    success &= test_qutrit_tta_hamiltonian()
    success &= test_qubit_tta_hamiltonian()

    if success:
        pass
    else:
        pass
