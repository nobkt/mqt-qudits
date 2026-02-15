#!/usr/bin/env python3
"""Verification script for the qudit accuracy fix.

This script verifies that the H_TTA bug has been fixed and provides
a quick check before running the full notebook.

Run with:
    python3 verify_qudit_accuracy_fix.py
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.linalg import expm


def test_H_TTA_unitarity() -> bool:
    """Test that H_TTA produces a unitary matrix."""
    # Test parameters
    J = 0.05  # eV
    dt = 10.0  # fs
    hbar = 0.6582  # eV·fs

    # Define Hamiltonian
    H_TTA = J * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])

    # Compute unitary
    U = expm(-1j * H_TTA * dt / hbar)

    # Check unitarity
    identity = U @ U.conj().T
    error = np.linalg.norm(identity - np.eye(3))

    return error < 1e-10


def test_H_TTA_structure() -> bool:
    """Test that H_TTA has the correct structure (diagonal real, off-diagonal imaginary)."""
    J = 0.05
    dt = 10.0
    hbar = 0.6582

    H_TTA = J * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    U = expm(-1j * H_TTA * dt / hbar)

    # Check diagonal elements are real
    diagonal_imag_max = max(abs(np.imag(U[i, i])) for i in range(3))

    # Check off-diagonal elements are imaginary (for this specific 3x3 H_TTA)
    # In the 3x3 subspace, all off-diagonal should be purely imaginary
    off_diagonal_real_max = 0
    for i in range(3):
        for j in range(3):
            if i != j and abs(U[i, j]) > 1e-10:  # Only check non-zero elements
                off_diagonal_real_max = max(off_diagonal_real_max, abs(np.real(U[i, j])))

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
        checks_pass = False

    # U[0,1], U[1,0], U[1,2], U[2,1] should be imaginary
    for i, j in [(0, 1), (1, 0), (1, 2), (2, 1)]:
        if abs(np.real(U[i, j])) > tolerance:
            checks_pass = False

    # U[0,2], U[2,0] should be real
    for i, j in [(0, 2), (2, 0)]:
        if abs(np.imag(U[i, j])) > tolerance:
            checks_pass = False

    return bool(checks_pass)


def test_exact_qudit_basic_gates():
    """Test the exact_qudit_basic_gates module."""
    try:
        sys.path.insert(0, "tutorials")
        from exact_qudit_basic_gates import verify_H_transfer_decomposition, verify_H_TTA_decomposition

        J = 0.05
        V = 0.1
        dt = 10.0
        hbar = 0.6582

        # Test H_transfer
        transfer_ok = bool(verify_H_transfer_decomposition(V, dt, hbar))

        # Test H_TTA
        tta_ok = bool(verify_H_TTA_decomposition(J, dt, hbar))

        return transfer_ok and tta_ok

    except Exception:
        return False


def main() -> int:
    """Run all verification tests."""
    results = []

    # Test 1: Unitarity
    results.append(test_H_TTA_unitarity())

    # Test 2: Structure
    results.append(test_H_TTA_structure())

    # Test 3: Module functions
    results.append(test_exact_qudit_basic_gates())

    # Summary

    if all(results):
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
