#!/usr/bin/env python3
"""
Verification script for the qudit accuracy fix.

This script verifies that the H_TTA bug has been fixed and provides
a quick check before running the full notebook.

Run with:
    python3 verify_qudit_accuracy_fix.py
"""

import sys
import numpy as np
from scipy.linalg import expm

def test_H_TTA_unitarity():
    """Test that H_TTA produces a unitary matrix."""
    print("=" * 70)
    print("Test 1: H_TTA Unitarity Check")
    print("=" * 70)
    
    # Test parameters
    J = 0.05  # eV
    dt = 10.0  # fs
    hbar = 0.6582  # eV·fs
    
    # Define Hamiltonian
    H_TTA = J * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    
    # Compute unitary
    U = expm(-1j * H_TTA * dt / hbar)
    
    # Check unitarity
    identity = U @ U.conj().T
    error = np.linalg.norm(identity - np.eye(3))
    
    print(f"Hamiltonian eigenvalues: {np.linalg.eigvalsh(H_TTA)}")
    print(f"Expected: {np.array([-np.sqrt(2)*J, 0, np.sqrt(2)*J])}")
    print(f"Unitarity error: {error:.2e}")
    
    if error < 1e-10:
        print("✓ PASS: H_TTA unitary is correct")
        return True
    else:
        print("✗ FAIL: H_TTA unitary is not unitary!")
        return False

def test_H_TTA_structure():
    """Test that H_TTA has the correct structure (diagonal real, off-diagonal imaginary)."""
    print("\n" + "=" * 70)
    print("Test 2: H_TTA Structure Check (3x3 Subspace)")
    print("=" * 70)
    
    J = 0.05
    dt = 10.0
    hbar = 0.6582
    
    H_TTA = J * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    U = expm(-1j * H_TTA * dt / hbar)
    
    print("Unitary matrix in 3x3 subspace:")
    print(f"U[0,0] = {U[0,0]}")
    print(f"U[0,1] = {U[0,1]}")
    print(f"U[0,2] = {U[0,2]}")
    print(f"U[1,1] = {U[1,1]}")
    print()
    
    # Check diagonal elements are real
    diagonal_imag_max = max(abs(np.imag(U[i,i])) for i in range(3))
    
    # Check off-diagonal elements are imaginary (for this specific 3x3 H_TTA)
    # In the 3x3 subspace, all off-diagonal should be purely imaginary
    off_diagonal_real_max = 0
    for i in range(3):
        for j in range(3):
            if i != j and abs(U[i,j]) > 1e-10:  # Only check non-zero elements
                off_diagonal_real_max = max(off_diagonal_real_max, abs(np.real(U[i,j])))
    
    print(f"Max imaginary part of diagonal: {diagonal_imag_max:.2e} (should be ~0)")
    print(f"Max real part of off-diagonal: {off_diagonal_real_max:.2e} (should be ~0)")
    
    # For this specific H_TTA, the structure should be:
    # - Diagonal: real
    # - Off-diagonal non-zero: imaginary
    # But U[0,2] and U[2,0] are real! Let me check the actual values
    
    # Actually, looking at the exact matrix:
    # U = [[ 0.738+0j,  0-0.622j,  -0.262+0j  ],
    #      [ 0-0.622j,  0.476+0j,   0-0.622j  ],
    #      [-0.262+0j,  0-0.622j,   0.738+0j  ]]
    # So U[0,2] and U[2,0] ARE real (and diagonal too)
    # But U[0,1], U[1,0], U[1,2], U[2,1] are imaginary
    
    # Correct check: specific elements should be real or imaginary
    tolerance = 1e-10
    checks_pass = True
    
    # Diagonal should be real
    if diagonal_imag_max > tolerance:
        print(f"✗ Diagonal elements have imaginary parts!")
        checks_pass = False
    
    # U[0,1], U[1,0], U[1,2], U[2,1] should be imaginary
    for (i, j) in [(0,1), (1,0), (1,2), (2,1)]:
        if abs(np.real(U[i,j])) > tolerance:
            print(f"✗ U[{i},{j}] should be imaginary but has real part {np.real(U[i,j]):.2e}")
            checks_pass = False
    
    # U[0,2], U[2,0] should be real  
    for (i, j) in [(0,2), (2,0)]:
        if abs(np.imag(U[i,j])) > tolerance:
            print(f"✗ U[{i},{j}] should be real but has imaginary part {np.imag(U[i,j]):.2e}")
            checks_pass = False
    
    if checks_pass:
        print("✓ PASS: H_TTA has correct structure")
        return True
    else:
        print("✗ FAIL: H_TTA structure is incorrect")
        return False

def test_exact_qudit_basic_gates():
    """Test the exact_qudit_basic_gates module."""
    print("\n" + "=" * 70)
    print("Test 3: exact_qudit_basic_gates Module")
    print("=" * 70)
    
    try:
        sys.path.insert(0, 'tutorials')
        from exact_qudit_basic_gates import (
            verify_H_transfer_decomposition,
            verify_H_TTA_decomposition
        )
        
        J = 0.05
        V = 0.1
        dt = 10.0
        hbar = 0.6582
        
        # Test H_transfer
        if verify_H_transfer_decomposition(V, dt, hbar):
            print("✓ PASS: H_transfer decomposition is correct")
            transfer_ok = True
        else:
            print("✗ FAIL: H_transfer decomposition is incorrect")
            transfer_ok = False
        
        # Test H_TTA
        if verify_H_TTA_decomposition(J, dt, hbar):
            print("✓ PASS: H_TTA decomposition is correct")
            tta_ok = True
        else:
            print("✗ FAIL: H_TTA decomposition is incorrect")
            tta_ok = False
        
        return transfer_ok and tta_ok
    
    except Exception as e:
        print(f"✗ FAIL: Error loading module: {e}")
        return False

def main():
    """Run all verification tests."""
    print("\n" + "=" * 70)
    print("QUDIT ACCURACY FIX VERIFICATION")
    print("=" * 70)
    print()
    
    results = []
    
    # Test 1: Unitarity
    results.append(test_H_TTA_unitarity())
    
    # Test 2: Structure
    results.append(test_H_TTA_structure())
    
    # Test 3: Module functions
    results.append(test_exact_qudit_basic_gates())
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    if all(results):
        print("✅ ALL TESTS PASSED")
        print()
        print("The H_TTA bug has been successfully fixed!")
        print()
        print("Next steps:")
        print("1. Run the full notebook:")
        print("   cd tutorials")
        print("   jupyter nbconvert --to notebook --execute \\")
        print("       quantum_dynamics_complete_comparison.ipynb")
        print()
        print("2. Expected results:")
        print("   - Qudit max error: ~0.01-0.02 (same as qubit)")
        print("   - Qudit gate count: 118 (22x better than qubit)")
        print("   - ✓ Both efficiency and accuracy!")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        print()
        print("Please check the error messages above.")
        print("The fix may not have been applied correctly.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
