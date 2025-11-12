#!/usr/bin/env python3
"""
Test if the H_TTA gate decomposition actually produces the correct unitary.

This is the missing verification - PR#86 verified that the TARGET unitary is correct,
but didn't verify that the GATE SEQUENCE produces that unitary.
"""

import sys
import numpy as np
from scipy.linalg import expm

# Add tutorials to path
sys.path.insert(0, 'tutorials')

from mqt.qudits.quantum_circuit import QuantumCircuit

def test_H_TTA_gate_decomposition():
    """
    Test if apply_H_TTA_basic_gates actually produces the correct unitary.
    """
    from exact_qudit_basic_gates import apply_H_TTA_basic_gates
    
    # Test parameters
    J = 0.05  # eV
    dt = 10.0  # fs
    hbar = 0.6582  # eV·fs
    
    # Expected unitary from scipy.linalg.expm
    H_TTA = J * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    U_expected = expm(-1j * H_TTA * dt / hbar)
    
    print("Expected unitary (from scipy.linalg.expm):")
    print(U_expected)
    print()
    
    # Create a circuit and apply the gates
    circuit = QuantumCircuit(2, [3, 3], 0)
    apply_H_TTA_basic_gates(circuit, 0, 1, J, dt, hbar)
    
    print(f"Number of gates applied: {len(circuit.instructions)}")
    print("Gate sequence:")
    for i, gate in enumerate(circuit.instructions):
        print(f"  {i+1}. {type(gate).__name__} on qudits {gate.target_qudits}")
    print()
    
    # Get the unitary matrix from the circuit
    # This would require executing the circuit and comparing the result
    # For now, we can at least check that gates were applied
    
    print("WARNING: Cannot directly verify the gate sequence produces correct unitary")
    print("without executing the circuit on actual MQT-Qudits backend.")
    print()
    print("The gate decomposition code says:")
    print('  "Note: This is a simplified version - a full decomposition would require"')
    print('  "more sophisticated gate sequence optimization"')
    print()
    print("This suggests the current decomposition may NOT be exact!")
    
    return circuit

if __name__ == '__main__':
    circuit = test_H_TTA_gate_decomposition()
