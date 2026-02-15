#!/usr/bin/env python3
"""
Analyze if the current H_TTA gate decomposition could possibly be correct.

Since we can't easily extract the unitary from a circuit, we'll analyze
the structure of the gates and compare with what we know must be true.
"""

import numpy as np
from scipy.linalg import expm

def analyze_correct_unitary():
    """Analyze the structure of the correct unitary matrix."""
    J = 0.05
    dt = 10.0
    hbar = 0.6582
    
    # Correct Hamiltonian
    H_TTA = J * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    
    # Correct unitary
    U = expm(-1j * H_TTA * dt / hbar)
    
    print("=" * 70)
    print("ANALYSIS OF CORRECT H_TTA UNITARY MATRIX")
    print("=" * 70)
    print()
    
    print("Parameters:")
    print(f"  J = {J} eV")
    print(f"  dt = {dt} fs")
    print(f"  ℏ = {hbar} eV·fs")
    print(f"  ω = √2·J·dt/ℏ = {np.sqrt(2)*J*dt/hbar:.6f}")
    print()
    
    print("Correct Unitary Matrix:")
    print(U)
    print()
    
    print("Matrix Properties:")
    print(f"  Unitarity error: {np.linalg.norm(U @ U.conj().T - np.eye(3)):.2e}")
    print(f"  Determinant: {np.linalg.det(U):.6f}")
    print()
    
    print("Element Analysis:")
    for i in range(3):
        for j in range(3):
            real_part = np.real(U[i,j])
            imag_part = np.imag(U[i,j])
            mag = np.abs(U[i,j])
            
            if abs(real_part) > 1e-10 and abs(imag_part) > 1e-10:
                print(f"  U[{i},{j}] = {U[i,j]:.6f} - BOTH real and imaginary!")
            elif abs(real_part) > 1e-10:
                print(f"  U[{i},{j}] = {real_part:.6f} (real)")
            elif abs(imag_part) > 1e-10:
                print(f"  U[{i},{j}] = {imag_part:.6f}j (imaginary)")
            else:
                print(f"  U[{i},{j}] = 0 (zero)")
    
    print()
    print("Structure Summary:")
    print("  Diagonal elements: ALL REAL")
    print("  Off-diagonal (0,1), (1,0), (1,2), (2,1): ALL PURELY IMAGINARY")
    print("  Off-diagonal (0,2), (2,0): ALL REAL")
    print()
    
    omega = np.sqrt(2) * J * dt / hbar
    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)
    
    print("Analytical Formula Verification:")
    print(f"  ω = {omega:.6f}")
    print(f"  cos(ω) = {cos_omega:.6f}")
    print(f"  sin(ω) = {sin_omega:.6f}")
    print()
    
    U_analytical = np.array([
        [0.5*(1+cos_omega),  -1j*sin_omega/np.sqrt(2),  -0.5*(1-cos_omega)],
        [-1j*sin_omega/np.sqrt(2),  cos_omega,  -1j*sin_omega/np.sqrt(2)],
        [-0.5*(1-cos_omega),  -1j*sin_omega/np.sqrt(2),  0.5*(1+cos_omega)]
    ])
    
    formula_error = np.linalg.norm(U - U_analytical)
    print(f"  Error vs analytical formula: {formula_error:.2e}")
    
    if formula_error < 1e-10:
        print("  ✓ Analytical formula MATCHES scipy.linalg.expm")
    else:
        print("  ✗ Analytical formula DOES NOT MATCH")
    
    print()
    return U

def analyze_gate_sequence_structure():
    """Analyze what the current gate sequence is trying to do."""
    print("=" * 70)
    print("ANALYSIS OF CURRENT GATE DECOMPOSITION")
    print("=" * 70)
    print()
    
    print("Current gate sequence from exact_qudit_basic_gates.py:")
    print()
    print("  1-6. VirtRz gates (phase preparation)")
    print("    - qudit_i: levels 0, 1, 2")
    print("    - qudit_j: levels 0, 1, 2")
    print("    - Level 1 gets phase -π/2 (prepare for imaginary coupling)")
    print()
    print("  7-8. R gates (single-qudit rotations)")
    print("    - qudit_i: levels 0-1, angle ω/2")
    print("    - qudit_j: levels 0-1, angle ω/2")
    print()
    print("  9. CEx gate (two-qudit coupling)")
    print("    - Control: qudit_i, level 1")
    print("    - Target: qudit_j, level 2")
    print("    - Angle: ω/√2")
    print()
    print("  10. CEx gate (two-qudit coupling)")
    print("    - Control: qudit_j, level 0")
    print("    - Target: qudit_i, level 2")
    print("    - Angle: ω/2")
    print()
    print("  11-12. VirtRz gates (phase correction)")
    print("    - qudit_i, level 1: +π/2")
    print("    - qudit_j, level 1: +π/2")
    print()
    
    print("CRITICAL ISSUE:")
    print("  The code admits this is a 'simplified version'")
    print("  Quote: 'Note: This is a simplified version - a full decomposition'")
    print("          'would require more sophisticated gate sequence optimization'")
    print()
    print("  This suggests the decomposition may NOT be mathematically exact!")
    print()

def check_if_verification_exists():
    """Check if there's any verification of the gate sequence."""
    print("=" * 70)
    print("VERIFICATION STATUS")
    print("=" * 70)
    print()
    
    print("Current verification (verify_H_TTA_decomposition):")
    print("  ✓ Verifies scipy.linalg.expm produces unitary matrix")
    print("  ✓ Verifies eigenvalues are correct")
    print("  ✓ Verifies analytical formula matches scipy.linalg.expm")
    print()
    print("  ✗ Does NOT verify gate sequence produces correct unitary!")
    print("  ✗ Does NOT check ||U_gates - U_target|| < ε")
    print()
    
    print("CONCLUSION:")
    print("  The gate decomposition has NEVER been verified to be correct.")
    print("  Given that the code admits it's 'simplified', it is likely WRONG.")
    print()

if __name__ == '__main__':
    U_correct = analyze_correct_unitary()
    print()
    analyze_gate_sequence_structure()
    print()
    check_if_verification_exists()
    
    print("=" * 70)
    print("RECOMMENDATION")
    print("=" * 70)
    print()
    print("Based on PR#86's requirements:")
    print("  1. NO heuristics or fallbacks allowed")
    print("  2. Must be mathematically exact")
    print("  3. Must completely fix the issue")
    print()
    print("The current 'simplified' gate decomposition does NOT meet these")
    print("requirements. We need to either:")
    print()
    print("  Option A: Implement and verify a correct decomposition")
    print("  Option B: Use CustomTwo gate with the exact unitary from scipy.linalg.expm")
    print()
    print("However, the problem statement says NOT to degrade the stable tutorial,")
    print("and the comments claim CustomTwo is not used. So Option A is required.")
    print()
    print("NEXT STEPS:")
    print("  1. Implement a provably correct gate decomposition")
    print("  2. Add verification that checks ||U_gates - U_target|| < 1e-10")
    print("  3. Ensure no 'simplified' or approximate methods")
    print()
