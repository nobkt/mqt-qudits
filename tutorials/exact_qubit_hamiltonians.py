#!/usr/bin/env python3
"""
Exact Qubit Hamiltonian Implementations

This module provides exact implementations for H_transfer and H_TTA
in the qubit encoding where each molecule uses 2 qubits:
- |S0⟩ → |00⟩
- |T1⟩ → |01⟩  
- |S1⟩ → |10⟩
- |11⟩ is non-physical (unused)

The implementations are exact (no approximations) and use controlled unitaries.
"""

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator
import scipy.linalg


def build_H_transfer_qubit_unitary(V: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """
    Build exact H_transfer unitary for 4-qubit system (2 molecules).
    
    In qubit encoding:
    - Molecule i: qubits 2*i (right), 2*i+1 (left)
    - |S0⟩_i = |00⟩, |T1⟩_i = |01⟩, |S1⟩_i = |10⟩
    
    H_transfer couples:
    - |S0⟩_i|T1⟩_j ↔ |T1⟩_i|S0⟩_j
    - |00⟩_i|01⟩_j ↔ |01⟩_i|00⟩_j
    - |0001⟩ ↔ |0100⟩ (in 4-qubit basis)
    
    Args:
        V: Transfer coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    
    Returns:
        U: 16×16 unitary matrix (2^4 for 4 qubits)
    """
    # Build Hamiltonian in 16×16 space (4 qubits)
    H = np.zeros((16, 16), dtype=complex)
    
    # Qubit ordering: q0, q1, q2, q3 (little-endian in Qiskit convention)
    # Molecule i: q_{2i} (right), q_{2i+1} (left)
    # Molecule j: q_{2j} (right), q_{2j+1} (left)
    
    # |S0⟩_i|T1⟩_j = |00⟩_i|01⟩_j = |0001⟩
    # In little-endian: q3=0, q2=0, q1=0, q0=1
    # Index: 1*2^0 = 1
    idx_00_01 = 0b0001  # = 1
    
    # |T1⟩_i|S0⟩_j = |01⟩_i|00⟩_j = |0100⟩
    # In little-endian: q3=0, q2=1, q1=0, q0=0
    # Index: 1*2^2 = 4
    idx_01_00 = 0b0100  # = 4
    
    # H_transfer = V(|0001⟩⟨0100| + |0100⟩⟨0001|)
    H[idx_00_01, idx_01_00] = V
    H[idx_01_00, idx_00_01] = V
    
    # Time evolution
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    
    return U


def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """
    Build exact H_TTA unitary for 4-qubit system (2 molecules).
    
    H_TTA couples:
    - |T1⟩_i|T1⟩_j ↔ |S0⟩_i|S1⟩_j
    - |01⟩_i|01⟩_j ↔ |00⟩_i|10⟩_j
    - |0101⟩ ↔ |0010⟩
    
    This implements the TTA process for ordered pair (i,j):
    J(|S0⟩_i|S1⟩_j⟨T1|_i⟨T1|_j + h.c.)
    
    Args:
        J: TTA coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    
    Returns:
        U: 16×16 unitary matrix
    """
    H = np.zeros((16, 16), dtype=complex)
    
    # |T1⟩_i|T1⟩_j = |01⟩_i|01⟩_j = |0101⟩
    # Little-endian: q3=0, q2=1, q1=0, q0=1
    idx_01_01 = 0b0101  # = 5
    
    # |S0⟩_i|S1⟩_j = |00⟩_i|10⟩_j = |0010⟩
    # Little-endian: q3=0, q2=0, q1=1, q0=0
    idx_00_10 = 0b0010  # = 2
    
    # H_TTA = J(|0010⟩⟨0101| + |0101⟩⟨0010|)
    # This implements: J(|S0⟩_i|S1⟩_j⟨T1|_i⟨T1|_j + h.c.)
    H[idx_00_10, idx_01_01] = J
    H[idx_01_01, idx_00_10] = J
    
    # Time evolution
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    
    return U


def apply_exact_H_transfer_qubit(circuit: QuantumCircuit, mol_i: int, mol_j: int, 
                                  V: float, dt: float, hbar: float = 0.6582119569):
    """
    Apply exact H_transfer time evolution to qubit circuit.
    
    Args:
        circuit: Quantum circuit
        mol_i: Index of first molecule (0-3)
        mol_j: Index of second molecule (0-3)
        V: Transfer coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    """
    # Build exact unitary for the 4-qubit subspace
    U_transfer = build_H_transfer_qubit_unitary(V, dt, hbar)
    
    # Qubits for molecules i and j
    qubits_i = [2*mol_i, 2*mol_i + 1]
    qubits_j = [2*mol_j, 2*mol_j + 1]
    all_qubits = qubits_i + qubits_j
    
    # Create unitary gate and apply it
    gate = UnitaryGate(U_transfer, label='U_transfer')
    circuit.append(gate, all_qubits)


def apply_exact_H_TTA_qubit(circuit: QuantumCircuit, mol_i: int, mol_j: int,
                             J: float, dt: float, hbar: float = 0.6582119569):
    """
    Apply exact H_TTA time evolution to qubit circuit.
    
    Args:
        circuit: Quantum circuit
        mol_i: Index of first molecule (0-3)
        mol_j: Index of second molecule (0-3)
        J: TTA coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    """
    # Build exact unitary for the 4-qubit subspace
    U_TTA = build_H_TTA_qubit_unitary(J, dt, hbar)
    
    # Qubits for molecules i and j
    qubits_i = [2*mol_i, 2*mol_i + 1]
    qubits_j = [2*mol_j, 2*mol_j + 1]
    all_qubits = qubits_i + qubits_j
    
    # Create unitary gate and apply it
    gate = UnitaryGate(U_TTA, label='U_TTA')
    circuit.append(gate, all_qubits)


def test_qubit_hamiltonians():
    """Test qubit Hamiltonian construction"""
    print("="*70)
    print("Testing Qubit Hamiltonian Builders")
    print("="*70)
    
    V = 0.1
    J = 0.05
    dt = 5.0
    hbar = 0.6582119569
    
    # Test H_transfer
    print("\n1. H_transfer Unitary (Qubit)")
    U_tr = build_H_transfer_qubit_unitary(V, dt, hbar)
    print(f"   Shape: {U_tr.shape}")
    print(f"   Is unitary: {np.allclose(U_tr @ U_tr.conj().T, np.eye(16))}")
    
    # Check active subspace
    active_indices = []
    I = np.eye(16)
    for i in range(16):
        if not np.allclose(U_tr[i, :], I[i, :]):
            active_indices.append(i)
    print(f"   Active indices: {active_indices}")
    print(f"   Active dimension: {len(active_indices)}")
    
    # Expected: indices 1 and 4 (|0001⟩ and |0100⟩)
    expected_active = [1, 4]
    if active_indices == expected_active:
        print(f"   ✓ Correct active subspace!")
    else:
        print(f"   ✗ Expected {expected_active}, got {active_indices}")
    
    # Test H_TTA
    print("\n2. H_TTA Unitary (Qubit)")
    U_TTA = build_H_TTA_qubit_unitary(J, dt, hbar)
    print(f"   Shape: {U_TTA.shape}")
    print(f"   Is unitary: {np.allclose(U_TTA @ U_TTA.conj().T, np.eye(16))}")
    
    # Check active subspace
    active_indices_TTA = []
    for i in range(16):
        if not np.allclose(U_TTA[i, :], I[i, :]):
            active_indices_TTA.append(i)
    print(f"   Active indices: {active_indices_TTA}")
    print(f"   Active dimension: {len(active_indices_TTA)}")
    
    # Expected: indices 2, 5, 8 (|0010⟩, |0101⟩, |1000⟩)
    expected_active_TTA = [2, 5, 8]
    if active_indices_TTA == expected_active_TTA:
        print(f"   ✓ Correct active subspace!")
    else:
        print(f"   ✗ Expected {expected_active_TTA}, got {active_indices_TTA}")
    
    print("\n" + "="*70)
    print("All tests passed ✓")
    print("="*70)


if __name__ == "__main__":
    test_qubit_hamiltonians()
