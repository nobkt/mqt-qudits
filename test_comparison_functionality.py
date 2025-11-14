#!/usr/bin/env python3
"""
Quick test to verify the new comparison functionality works correctly.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'tutorials'))

import numpy as np
from qiskit import QuantumCircuit

# Test 1: comparison_helpers
print("="*70)
print("Test 1: comparison_helpers module")
print("="*70)

from comparison_helpers import count_gates_by_type, compare_gate_counts

# Create a simple Qiskit circuit
qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)
qc.rz(0.5, 0)
qc.ry(0.3, 1)

gates = count_gates_by_type(qc, is_qiskit=True)
print(f"\nGate counts: {gates}")
print("✓ count_gates_by_type works")

# Test comparison
gates2 = {'h': 2, 'cx': 2, 'rz': 1}
compare_gate_counts(gates, gates2, "Circuit 1", "Circuit 2")
print("✓ compare_gate_counts works")

# Test 2: qubit_unitary_simulator
print("\n" + "="*70)
print("Test 2: qubit_unitary_simulator module")
print("="*70)

from qubit_unitary_simulator import QubitMolecularDynamicsSimulatorUnitary

# Simple parameters for testing
class TestParams:
    def __init__(self):
        self.N_molecules = 2  # Smaller for quick test
        self.E_T = 1.5
        self.E_S = 3.0
        self.V = 0.1
        self.J = 0.05
        self.hbar = 0.6582119569
        self.neighbors = [(0, 1)]

params = TestParams()
simulator = QubitMolecularDynamicsSimulatorUnitary(params)
print("✓ QubitMolecularDynamicsSimulatorUnitary initialized")

# Build a single step circuit
circuit = simulator.build_single_trotter_step(dt=5.0)
print(f"✓ Built circuit with {len(circuit.data)} gates")

# Test 3: exact_qubit_hamiltonians
print("\n" + "="*70)
print("Test 3: exact_qubit_hamiltonians module")
print("="*70)

from exact_qubit_hamiltonians import (
    build_H_transfer_qubit_unitary,
    build_H_TTA_qubit_unitary,
    apply_exact_H_transfer_qubit,
    apply_exact_H_TTA_qubit
)

# Test unitary builders
U_transfer = build_H_transfer_qubit_unitary(V=0.1, dt=5.0, hbar=0.6582)
print(f"✓ Built H_transfer unitary: {U_transfer.shape}")

# Check unitarity
identity = U_transfer @ U_transfer.conj().T
error = np.linalg.norm(identity - np.eye(16))
print(f"  Unitarity error: {error:.2e}")
assert error < 1e-10, "H_transfer is not unitary!"

U_TTA = build_H_TTA_qubit_unitary(J=0.05, dt=5.0, hbar=0.6582)
print(f"✓ Built H_TTA unitary: {U_TTA.shape}")

# Check unitarity
identity = U_TTA @ U_TTA.conj().T
error = np.linalg.norm(identity - np.eye(16))
print(f"  Unitarity error: {error:.2e}")
assert error < 1e-10, "H_TTA is not unitary!"

# Test gate application
test_circuit = QuantumCircuit(4)
apply_exact_H_transfer_qubit(test_circuit, 0, 1, V=0.1, dt=5.0, hbar=0.6582)
print(f"✓ Applied H_transfer gate: {len(test_circuit.data)} gates added")

apply_exact_H_TTA_qubit(test_circuit, 0, 1, J=0.05, dt=5.0, hbar=0.6582)
print(f"✓ Applied H_TTA gate: {len(test_circuit.data)} total gates")

# Final summary
print("\n" + "="*70)
print("All tests passed!")
print("="*70)
print("\nThe new comparison functionality is working correctly:")
print("  ✓ comparison_helpers: gate counting and comparison")
print("  ✓ qubit_unitary_simulator: UnitaryGate-based simulation")
print("  ✓ exact_qubit_hamiltonians: exact unitary construction")
print("\nThe notebook should execute successfully.")
