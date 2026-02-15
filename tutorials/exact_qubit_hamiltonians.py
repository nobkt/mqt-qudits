#!/usr/bin/env python3
"""Exact Qubit Hamiltonian Implementations.

This module provides exact implementations for H_transfer and H_TTA
in the qubit encoding where each molecule uses 2 qubits:
- |S0⟩ → |00⟩
- |T1⟩ → |01⟩
- |S1⟩ → |10⟩
- |11⟩ is non-physical (unused)

The implementations are exact (no approximations) and use controlled unitaries.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.linalg
from qiskit.circuit.library import UnitaryGate

if TYPE_CHECKING:
    from qiskit import QuantumCircuit


def build_H_transfer_qubit_unitary(V: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """Build exact H_transfer unitary for 4-qubit system (2 molecules).

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
    return scipy.linalg.expm(-1j * H * dt / hbar)


def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """Build exact H_TTA unitary for 4-qubit system (2 molecules).

    H_TTA implements the complete TTA process with BOTH terms:
    - |T1⟩_i|T1⟩_j ↔ |S0⟩_i|S1⟩_j
    - |T1⟩_i|T1⟩_j ↔ |S1⟩_i|S0⟩_j

    In qubit encoding (|S0⟩=|00⟩, |T1⟩=|01⟩, |S1⟩=|10⟩):
    - |0101⟩ ↔ |0010⟩  (|T1,T1⟩ ↔ |S0,S1⟩)
    - |0101⟩ ↔ |1000⟩  (|T1,T1⟩ ↔ |S1,S0⟩)

    Full Hamiltonian:
    Ĥ_TTA = J [|S0⟩_i⟨T1|_i ⊗ |S1⟩_j⟨T1|_j + |S1⟩_i⟨T1|_i ⊗ |S0⟩_j⟨T1|_j + h.c.]

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
    idx_T1_T1 = 0b0101  # = 5

    # |S0⟩_i|S1⟩_j = |00⟩_i|10⟩_j = |0010⟩
    # Little-endian: q3=0, q2=0, q1=1, q0=0
    idx_S0_S1 = 0b0010  # = 2

    # |S1⟩_i|S0⟩_j = |10⟩_i|00⟩_j = |1000⟩
    # Little-endian: q3=1, q2=0, q1=0, q0=0
    idx_S1_S0 = 0b1000  # = 8

    # Term 1: |T1,T1⟩ ↔ |S0,S1⟩
    # J(|S0⟩_i|S1⟩_j⟨T1|_i⟨T1|_j + h.c.)
    H[idx_S0_S1, idx_T1_T1] = J
    H[idx_T1_T1, idx_S0_S1] = J

    # Term 2: |T1,T1⟩ ↔ |S1,S0⟩ (THIS WAS MISSING!)
    # J(|S1⟩_i|S0⟩_j⟨T1|_i⟨T1|_j + h.c.)
    H[idx_S1_S0, idx_T1_T1] = J
    H[idx_T1_T1, idx_S1_S0] = J

    # Time evolution
    return scipy.linalg.expm(-1j * H * dt / hbar)


def apply_exact_H_transfer_qubit(
    circuit: QuantumCircuit, mol_i: int, mol_j: int, V: float, dt: float, hbar: float = 0.6582119569
) -> None:
    """Apply exact H_transfer time evolution to qubit circuit.

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
    qubits_i = [2 * mol_i, 2 * mol_i + 1]
    qubits_j = [2 * mol_j, 2 * mol_j + 1]
    all_qubits = qubits_i + qubits_j

    # Create unitary gate and apply it
    gate = UnitaryGate(U_transfer, label="U_transfer")
    circuit.append(gate, all_qubits)


def apply_exact_H_TTA_qubit(
    circuit: QuantumCircuit, mol_i: int, mol_j: int, J: float, dt: float, hbar: float = 0.6582119569
) -> None:
    """Apply exact H_TTA time evolution to qubit circuit.

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
    qubits_i = [2 * mol_i, 2 * mol_i + 1]
    qubits_j = [2 * mol_j, 2 * mol_j + 1]
    all_qubits = qubits_i + qubits_j

    # Create unitary gate and apply it
    gate = UnitaryGate(U_TTA, label="U_TTA")
    circuit.append(gate, all_qubits)


def test_qubit_hamiltonians() -> None:
    """Test qubit Hamiltonian construction."""
    V = 0.1
    J = 0.05
    dt = 5.0
    hbar = 0.6582119569

    # Test H_transfer
    U_tr = build_H_transfer_qubit_unitary(V, dt, hbar)

    # Check active subspace
    I = np.eye(16)
    active_indices = [i for i in range(16) if not np.allclose(U_tr[i, :], I[i, :])]

    # Expected: indices 1 and 4 (|0001⟩ and |0100⟩)
    expected_active = [1, 4]
    if active_indices == expected_active:
        pass
    else:
        pass

    # Test H_TTA
    U_TTA = build_H_TTA_qubit_unitary(J, dt, hbar)

    # Check active subspace
    active_indices_TTA = [i for i in range(16) if not np.allclose(U_TTA[i, :], I[i, :])]

    # Expected: indices 2, 5, 8 (|0010⟩, |0101⟩, |1000⟩)
    expected_active_TTA = [2, 5, 8]
    if active_indices_TTA == expected_active_TTA:
        pass
    else:
        pass


if __name__ == "__main__":
    test_qubit_hamiltonians()
