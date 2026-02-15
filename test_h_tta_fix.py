#!/usr/bin/env python3
"""Comprehensive test to verify the H_TTA fix is complete and correct.
This test validates that the fix resolves the accuracy issue reported in the problem statement.
"""

from __future__ import annotations

import sys

sys.path.insert(0, "/home/runner/work/mqt-qudits/mqt-qudits/tutorials")

import numpy as np
from exact_qudit_basic_gates import apply_H_TTA_basic_gates, verify_H_TTA_decomposition
from scipy.linalg import expm

from mqt.qudits.quantum_circuit import QuantumCircuit

# Test parameters (matching notebook)
J = 0.05  # eV
dt = 10.0  # fs
hbar = 0.6582  # eV·fs

# Test 1: Verify the verification function passes
result = verify_H_TTA_decomposition(J, dt, hbar)
if result:
    pass
else:
    sys.exit(1)

# Test 2: Verify the exact unitary is computed correctly
H_TTA = J * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
U_classical = expm(-1j * H_TTA * dt / hbar)

# Check analytical formula
omega = np.sqrt(2) * J * dt / hbar
cos_omega = np.cos(omega)
sin_omega = np.sin(omega)

U_analytical = np.array([
    [0.5 * (1 + cos_omega), -1j * sin_omega / np.sqrt(2), -0.5 * (1 - cos_omega)],
    [-1j * sin_omega / np.sqrt(2), cos_omega, -1j * sin_omega / np.sqrt(2)],
    [-0.5 * (1 - cos_omega), -1j * sin_omega / np.sqrt(2), 0.5 * (1 + cos_omega)],
])

formula_error = np.linalg.norm(U_classical - U_analytical)
if formula_error < 1e-10:
    pass
else:
    sys.exit(1)

# Test 3: Verify CustomTwo gate creation
circuit = QuantumCircuit(2, [3, 3], 0)
apply_H_TTA_basic_gates(circuit, 0, 1, J, dt, hbar)

if len(circuit.instructions) != 1:
    sys.exit(1)

gate = circuit.instructions[0]
if type(gate).__name__ != "CustomTwo":
    sys.exit(1)


# Test 4: Verify the CustomTwo gate contains the exact unitary
U_9x9 = gate.to_matrix(identities=0)

# Check unitarity of 9×9 matrix
unitarity_error = np.linalg.norm(U_9x9 @ U_9x9.conj().T - np.eye(9))
if unitarity_error > 1e-10:
    sys.exit(1)

# Extract active subspace
active_indices = [2, 4, 6]  # |02⟩, |11⟩, |20⟩
U_active = U_9x9[np.ix_(active_indices, active_indices)]

# Compare with classical
active_error = np.linalg.norm(U_active - U_classical)

if active_error < 1e-10:
    pass
else:
    sys.exit(1)

# Test 5: Verify inactive subspace is identity
all_indices = list(range(9))
inactive_indices = [i for i in all_indices if i not in active_indices]

max_inactive_error = 0
for i in inactive_indices:
    for j in range(9):
        expected = 1.0 if i == j else 0.0
        error = abs(U_9x9[i, j] - expected)
        max_inactive_error = max(max_inactive_error, error)

if max_inactive_error < 1e-10:
    pass
else:
    sys.exit(1)

# Summary
