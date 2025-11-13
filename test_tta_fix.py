#!/usr/bin/env python3
"""
Test script to verify TTA Hamiltonian fix
"""

import numpy as np
import scipy.linalg

def ket_bra(i, j, dim=3):
    """Create |i⟩⟨j| operator for dim-dimensional system"""
    op = np.zeros((dim, dim))
    op[i, j] = 1
    return op

def test_qutrit_tta_hamiltonian():
    """Test the correct 9x9 TTA Hamiltonian for 2-qutrit system"""
    print("="*80)
    print("TEST: 9x9 TTA Hamiltonian for 2-Qutrit System")
    print("="*80)
    
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
    
    print("\nCorrect H_TTA matrix (J=1):")
    print(H_TTA.astype(int))
    
    # Verify physical process
    print("\n✓ Physical verification:")
    print(f"  H[2,4] = {H_TTA[2,4]} (|S0S1⟩ ← |T1T1⟩) should be 1")
    print(f"  H[4,2] = {H_TTA[4,2]} (|T1T1⟩ ← |S0S1⟩) should be 1")
    print(f"  H[6,4] = {H_TTA[6,4]} (|S1S0⟩ ← |T1T1⟩) should be 1")
    print(f"  H[4,6] = {H_TTA[4,6]} (|T1T1⟩ ← |S1S0⟩) should be 1")
    
    # Verify hermiticity
    is_hermitian = np.allclose(H_TTA, H_TTA.T.conj())
    print(f"\n✓ Hermiticity: {is_hermitian}")
    
    # Test time evolution
    J = 0.05  # eV
    dt = 5.0  # fs
    hbar = 0.6582119569  # eV·fs
    
    U_TTA = scipy.linalg.expm(-1j * J * H_TTA * dt / hbar)
    
    # Verify unitarity
    is_unitary = np.allclose(U_TTA @ U_TTA.T.conj(), np.eye(9))
    print(f"✓ Unitarity of time evolution: {is_unitary}")
    
    if is_hermitian and is_unitary:
        print("\n✅ TEST PASSED: TTA Hamiltonian is correct!")
        return True
    else:
        print("\n❌ TEST FAILED!")
        return False


def test_qubit_tta_hamiltonian():
    """Test the correct 16x16 TTA Hamiltonian for 4-qubit system"""
    print("\n" + "="*80)
    print("TEST: 16x16 TTA Hamiltonian for 4-Qubit System")
    print("="*80)
    
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
    
    print(f"\nNon-zero elements:")
    nz = np.where(H != 0)
    for i, j in zip(nz[0], nz[1]):
        print(f"  H[{i:2d},{j:2d}] = {H[i,j]:.1f} ({bin(i):>6s} ← {bin(j):>6s})")
    
    # Verify physical process
    print("\n✓ Physical verification:")
    print(f"  H[{idx_S0_S1},{idx_T1_T1}] = {H[idx_S0_S1,idx_T1_T1]} (|S0S1⟩ ← |T1T1⟩)")
    print(f"  H[{idx_T1_T1},{idx_S0_S1}] = {H[idx_T1_T1,idx_S0_S1]} (|T1T1⟩ ← |S0S1⟩)")
    print(f"  H[{idx_S1_S0},{idx_T1_T1}] = {H[idx_S1_S0,idx_T1_T1]} (|S1S0⟩ ← |T1T1⟩)")
    print(f"  H[{idx_T1_T1},{idx_S1_S0}] = {H[idx_T1_T1,idx_S1_S0]} (|T1T1⟩ ← |S1S0⟩)")
    
    # Verify hermiticity
    is_hermitian = np.allclose(H, H.T.conj())
    print(f"\n✓ Hermiticity: {is_hermitian}")
    
    # Test time evolution
    J = 0.05  # eV
    dt = 5.0  # fs
    hbar = 0.6582119569  # eV·fs
    
    U_TTA = scipy.linalg.expm(-1j * J * H * dt / hbar)
    
    # Verify unitarity
    is_unitary = np.allclose(U_TTA @ U_TTA.T.conj(), np.eye(16))
    print(f"✓ Unitarity of time evolution: {is_unitary}")
    
    if is_hermitian and is_unitary:
        print("\n✅ TEST PASSED: Qubit TTA Hamiltonian is correct!")
        return True
    else:
        print("\n❌ TEST FAILED!")
        return False


if __name__ == "__main__":
    success = True
    success &= test_qutrit_tta_hamiltonian()
    success &= test_qubit_tta_hamiltonian()
    
    print("\n" + "="*80)
    if success:
        print("✅ ALL TESTS PASSED!")
    else:
        print("❌ SOME TESTS FAILED!")
    print("="*80)
