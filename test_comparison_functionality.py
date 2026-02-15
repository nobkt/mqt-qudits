#!/usr/bin/env python3
"""Quick test to verify the new comparison functionality works correctly."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "tutorials"))

import numpy as np

# Test 1: comparison_helpers
from comparison_helpers import compare_gate_counts, count_gates_by_type
from qiskit import QuantumCircuit

# Create a simple Qiskit circuit
qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)
qc.rz(0.5, 0)
qc.ry(0.3, 1)

gates = count_gates_by_type(qc, is_qiskit=True)

# Test comparison
gates2 = {"h": 2, "cx": 2, "rz": 1}
compare_gate_counts(gates, gates2, "Circuit 1", "Circuit 2")

# Test 2: qubit_unitary_simulator

from qubit_unitary_simulator import QubitMolecularDynamicsSimulatorUnitary


# Simple parameters for testing
class TestParams:
    def __init__(self) -> None:
        self.N_molecules = 2  # Smaller for quick test
        self.E_T = 1.5
        self.E_S = 3.0
        self.V = 0.1
        self.J = 0.05
        self.hbar = 0.6582119569
        self.neighbors = [(0, 1)]


params = TestParams()
simulator = QubitMolecularDynamicsSimulatorUnitary(params)

# Build a single step circuit
circuit = simulator.build_single_trotter_step(dt=5.0)

# Test 3: exact_qubit_hamiltonians

from exact_qubit_hamiltonians import (
    apply_exact_H_transfer_qubit,
    apply_exact_H_TTA_qubit,
    build_H_transfer_qubit_unitary,
    build_H_TTA_qubit_unitary,
)

# Test unitary builders
U_transfer = build_H_transfer_qubit_unitary(V=0.1, dt=5.0, hbar=0.6582)

# Check unitarity
identity = U_transfer @ U_transfer.conj().T
error = np.linalg.norm(identity - np.eye(16))
assert error < 1e-10, "H_transfer is not unitary!"

U_TTA = build_H_TTA_qubit_unitary(J=0.05, dt=5.0, hbar=0.6582)

# Check unitarity
identity = U_TTA @ U_TTA.conj().T
error = np.linalg.norm(identity - np.eye(16))
assert error < 1e-10, "H_TTA is not unitary!"

# Test gate application
test_circuit = QuantumCircuit(4)
apply_exact_H_transfer_qubit(test_circuit, 0, 1, V=0.1, dt=5.0, hbar=0.6582)

apply_exact_H_TTA_qubit(test_circuit, 0, 1, J=0.05, dt=5.0, hbar=0.6582)

# Final summary
