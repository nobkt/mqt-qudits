#!/usr/bin/env python3
"""
Comprehensive test to verify the H_TTA fix is complete and correct.
This test validates that the fix resolves the accuracy issue reported in the problem statement.
"""

import sys
sys.path.insert(0, '/home/runner/work/mqt-qudits/mqt-qudits/tutorials')

import numpy as np
from scipy.linalg import expm
from mqt.qudits.quantum_circuit import QuantumCircuit
from exact_qudit_basic_gates import apply_H_TTA_basic_gates, verify_H_TTA_decomposition

print("="*70)
print("COMPREHENSIVE H_TTA FIX VERIFICATION")
print("="*70)
print()

# Test parameters (matching notebook)
J = 0.05  # eV
dt = 10.0  # fs
hbar = 0.6582  # eV·fs

# Test 1: Verify the verification function passes
print("Test 1: Verification function")
print("-" * 70)
result = verify_H_TTA_decomposition(J, dt, hbar)
if result:
    print("✅ PASS: verify_H_TTA_decomposition() returns True")
else:
    print("❌ FAIL: verify_H_TTA_decomposition() returns False")
    sys.exit(1)
print()

# Test 2: Verify the exact unitary is computed correctly
print("Test 2: Exact unitary computation")
print("-" * 70)
H_TTA = J * np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])
U_classical = expm(-1j * H_TTA * dt / hbar)

# Check analytical formula
omega = np.sqrt(2) * J * dt / hbar
cos_omega = np.cos(omega)
sin_omega = np.sin(omega)

U_analytical = np.array([
    [0.5*(1+cos_omega),  -1j*sin_omega/np.sqrt(2),  -0.5*(1-cos_omega)],
    [-1j*sin_omega/np.sqrt(2),  cos_omega,  -1j*sin_omega/np.sqrt(2)],
    [-0.5*(1-cos_omega),  -1j*sin_omega/np.sqrt(2),  0.5*(1+cos_omega)]
])

formula_error = np.linalg.norm(U_classical - U_analytical)
print(f"||U_expm - U_analytical|| = {formula_error:.2e}")
if formula_error < 1e-10:
    print("✅ PASS: Analytical formula matches scipy.linalg.expm")
else:
    print(f"❌ FAIL: Formula error too large: {formula_error:.2e}")
    sys.exit(1)
print()

# Test 3: Verify CustomTwo gate creation
print("Test 3: CustomTwo gate creation")
print("-" * 70)
circuit = QuantumCircuit(2, [3, 3], 0)
apply_H_TTA_basic_gates(circuit, 0, 1, J, dt, hbar)

if len(circuit.instructions) != 1:
    print(f"❌ FAIL: Expected 1 gate, got {len(circuit.instructions)}")
    sys.exit(1)

gate = circuit.instructions[0]
if type(gate).__name__ != 'CustomTwo':
    print(f"❌ FAIL: Expected CustomTwo gate, got {type(gate).__name__}")
    sys.exit(1)

print(f"✅ PASS: Created 1 CustomTwo gate")
print()

# Test 4: Verify the CustomTwo gate contains the exact unitary
print("Test 4: CustomTwo gate unitary correctness")
print("-" * 70)
U_9x9 = gate.to_matrix(identities=0)

# Check unitarity of 9×9 matrix
unitarity_error = np.linalg.norm(U_9x9 @ U_9x9.conj().T - np.eye(9))
print(f"9×9 unitarity error: {unitarity_error:.2e}")
if unitarity_error > 1e-10:
    print(f"❌ FAIL: 9×9 matrix is not unitary")
    sys.exit(1)

# Extract active subspace
active_indices = [2, 4, 6]  # |02⟩, |11⟩, |20⟩
U_active = U_9x9[np.ix_(active_indices, active_indices)]

# Compare with classical
active_error = np.linalg.norm(U_active - U_classical)
print(f"Active subspace error: {active_error:.2e}")

if active_error < 1e-10:
    print("✅ PASS: CustomTwo gate contains exact H_TTA unitary")
else:
    print(f"❌ FAIL: Active subspace error too large: {active_error:.2e}")
    sys.exit(1)
print()

# Test 5: Verify inactive subspace is identity
print("Test 5: Inactive subspace is identity")
print("-" * 70)
all_indices = list(range(9))
inactive_indices = [i for i in all_indices if i not in active_indices]

max_inactive_error = 0
for i in inactive_indices:
    for j in range(9):
        if i == j:
            expected = 1.0
        else:
            expected = 0.0
        error = abs(U_9x9[i, j] - expected)
        max_inactive_error = max(max_inactive_error, error)

print(f"Max inactive element error: {max_inactive_error:.2e}")
if max_inactive_error < 1e-10:
    print("✅ PASS: Inactive subspace is identity")
else:
    print(f"❌ FAIL: Inactive subspace error: {max_inactive_error:.2e}")
    sys.exit(1)
print()

# Summary
print("="*70)
print("FINAL RESULT: ALL TESTS PASSED ✅")
print("="*70)
print()
print("Summary:")
print("  - H_TTA unitary is mathematically exact (scipy.linalg.expm)")
print("  - Analytical formula verified")
print("  - CustomTwo gate created with exact 9×9 unitary")
print("  - Active subspace matches classical: error < 1e-10")
print("  - Inactive subspace is identity: error < 1e-10")
print()
print("The fix should resolve the accuracy problem:")
print("  Before: Max error 0.615, Avg error 0.237")
print("  After:  Max error ~0.01, Avg error ~0.002 (expected)")
print()
