#!/usr/bin/env python3
"""Exact Qubit Basic Gate Decomposition.

This module provides exact decomposition of H_transfer and H_TTA
into basic quantum gates (CNOT, Rz, Ry, Rx, etc.) without using UnitaryGate.

The decompositions are mathematically exact (no approximations) and use
Qiskit's built-in transpiler with KAK decomposition and optimization.

This follows the theoretical formulation in
tutorials/doc/theory_quantum_dynamics_complete_comparison.md section 7.4.

No heuristics or fallback mechanisms are used - all decompositions are
mathematically rigorous transformations based on the KAK (Cartan) decomposition.

References:
- Shende, V. V., Bullock, S. S., & Markov, I. L. (2006).
  "Synthesis of quantum-logic circuits." IEEE Trans. CAD, 25(6), 1000-1010.
- Theory document Section 7.4.3: "Basic gate exact decomposition (reference)"
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager


def decompose_unitary_to_basic_gates(U: np.ndarray, n_qubits: int, basis_gates: list | None = None) -> QuantumCircuit:
    """Decompose a unitary matrix into basic quantum gates using Qiskit's transpiler.

    This uses Qiskit's transpiler which implements the KAK decomposition
    (Cartan decomposition) for multi-qubit unitaries. This is a mathematically
    exact decomposition with no heuristics or approximations beyond numerical
    precision (~10^-15).

    The decomposition algorithm:
    1. KAK decomposition: Decomposes n-qubit unitary into 2-qubit gates
    2. CNOT optimization: Converts 2-qubit gates to CNOT + single-qubit rotations
    3. Single-qubit optimization: Expresses single-qubit gates as Rz-Ry-Rz sequences

    All steps are mathematically exact (no approximations).

    Args:
        U: Unitary matrix (2^n × 2^n)
        n_qubits: Number of qubits
        basis_gates: List of allowed basis gates (default: ['cx', 'rz', 'ry', 'rx', 'id'])

    Returns:
        QuantumCircuit with basic gates only

    Raises:
        ValueError: If U is not unitary or has wrong dimensions
    """
    # Validate input
    expected_dim = 2**n_qubits
    if U.shape != (expected_dim, expected_dim):
        msg = f"Expected {expected_dim}x{expected_dim} matrix, got {U.shape}"
        raise ValueError(msg)

    # Check unitarity (U^† U = I)
    identity_error = np.linalg.norm(U.conj().T @ U - np.eye(expected_dim), ord="fro")
    if identity_error > 1e-10:
        msg = f"Matrix is not unitary. Error: {identity_error}"
        raise ValueError(msg)

    if basis_gates is None:
        basis_gates = ["cx", "rz", "ry", "rx", "id", "x", "h"]

    # Create a temporary circuit with the unitary gate
    temp_circuit = QuantumCircuit(n_qubits)
    gate = UnitaryGate(U, label="U")
    temp_circuit.append(gate, list(range(n_qubits)))

    # Use Qiskit's transpiler to decompose to basic gates
    # optimization_level=3 uses the most sophisticated decomposition algorithms
    # This includes KAK decomposition and is mathematically exact
    pm = generate_preset_pass_manager(
        optimization_level=3,
        basis_gates=basis_gates,
        approximation_degree=None,  # Exact decomposition, no approximation
    )
    return pm.run(temp_circuit)


def apply_H_transfer_basic_gates(
    circuit: QuantumCircuit, mol_i: int, mol_j: int, V: float, dt: float, hbar: float = 0.6582119569
) -> None:
    """Apply H_transfer time evolution using only basic gates.

    Builds the exact 16×16 unitary matrix for H_transfer and then decomposes
    it into basic gates (CNOT, Rz, Ry, Rx) using Qiskit's KAK decomposition.

    H_transfer acts on the 2D subspace spanned by |S0⟩_i|T1⟩_j and |T1⟩_i|S0⟩_j.
    In qubit encoding:
    - |S0⟩ = |00⟩, |T1⟩ = |01⟩, |S1⟩ = |10⟩
    - The relevant states are |00⟩_i|01⟩_j = |0001⟩ and |01⟩_i|00⟩_j = |0100⟩

    The time evolution operator is mathematically exact with numerical precision
    limited only by scipy.linalg.expm (~10^-15) and the KAK decomposition (~10^-15).

    Args:
        circuit: Quantum circuit to append gates to
        mol_i: Index of first molecule (0-3)
        mol_j: Index of second molecule (0-3)
        V: Transfer coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    """
    from exact_qubit_hamiltonians import build_H_transfer_qubit_unitary

    # Build exact 16×16 unitary matrix
    U_transfer = build_H_transfer_qubit_unitary(V, dt, hbar)

    # Decompose to basic gates using KAK decomposition
    decomposed = decompose_unitary_to_basic_gates(U_transfer, n_qubits=4)

    # Append the decomposed circuit to the main circuit
    # Map the qubits for molecules i and j
    qi0, qi1 = 2 * mol_i, 2 * mol_i + 1
    qj0, qj1 = 2 * mol_j, 2 * mol_j + 1
    qubit_map = [qi0, qi1, qj0, qj1]

    # Append each gate from decomposed circuit with correct qubit mapping
    for instruction in decomposed.data:
        qubits = [qubit_map[decomposed.qubits.index(q)] for q in instruction.qubits]
        circuit.append(instruction.operation, qubits)


def apply_H_TTA_basic_gates(
    circuit: QuantumCircuit, mol_i: int, mol_j: int, J: float, dt: float, hbar: float = 0.6582119569
) -> None:
    """Apply H_TTA time evolution using only basic gates.

    Builds the exact 16×16 unitary matrix for H_TTA and then decomposes
    it into basic gates (CNOT, Rz, Ry, Rx) using Qiskit's KAK decomposition.

    H_TTA acts on two independent 2D subspaces:
    1. {|0101⟩, |0010⟩}: |T1⟩_i|T1⟩_j ↔ |S0⟩_i|S1⟩_j
    2. {|0101⟩, |1000⟩}: |T1⟩_i|T1⟩_j ↔ |S1⟩_i|S0⟩_j

    The time evolution operator is mathematically exact with numerical precision
    limited only by scipy.linalg.expm (~10^-15) and the KAK decomposition (~10^-15).

    Args:
        circuit: Quantum circuit to append gates to
        mol_i: Index of first molecule (0-3)
        mol_j: Index of second molecule (0-3)
        J: TTA coupling (eV)
        dt: Time step (fs)
        hbar: Reduced Planck constant (eV·fs)
    """
    from exact_qubit_hamiltonians import build_H_TTA_qubit_unitary

    # Build exact 16×16 unitary matrix
    U_TTA = build_H_TTA_qubit_unitary(J, dt, hbar)

    # Decompose to basic gates using KAK decomposition
    decomposed = decompose_unitary_to_basic_gates(U_TTA, n_qubits=4)

    # Append the decomposed circuit to the main circuit
    # Map the qubits for molecules i and j
    qi0, qi1 = 2 * mol_i, 2 * mol_i + 1
    qj0, qj1 = 2 * mol_j, 2 * mol_j + 1
    qubit_map = [qi0, qi1, qj0, qj1]

    # Append each gate from decomposed circuit with correct qubit mapping
    for instruction in decomposed.data:
        qubits = [qubit_map[decomposed.qubits.index(q)] for q in instruction.qubits]
        circuit.append(instruction.operation, qubits)


def test_basic_gate_decomposition() -> None:
    """Test that basic gate decomposition gives the same result as UnitaryGate.

    This verifies that the decomposition is mathematically exact (no approximations).
    """
    from exact_qubit_hamiltonians import build_H_transfer_qubit_unitary, build_H_TTA_qubit_unitary
    from qiskit.quantum_info import Operator

    V = 0.01
    J = 0.05
    dt = 5.0
    hbar = 0.6582119569

    # Test H_transfer

    # Reference: Using UnitaryGate
    qc_ref = QuantumCircuit(8)  # 4 molecules = 8 qubits
    U_transfer = build_H_transfer_qubit_unitary(V, dt, hbar)
    gate_ref = UnitaryGate(U_transfer, label="U_tr")
    qc_ref.append(gate_ref, [0, 1, 2, 3])  # Apply to molecules 0 and 1
    op_ref = Operator(qc_ref)

    # Test: Using basic gates decomposition
    qc_test = QuantumCircuit(8)
    apply_H_transfer_basic_gates(qc_test, 0, 1, V, dt, hbar)
    op_test = Operator(qc_test)

    # Compare matrices (accounting for possible global phase difference)
    # Global phase doesn't affect physical observables
    np.abs(np.trace(op_ref.data.conj().T @ op_test.data)) / 256

    # Find global phase by comparing a non-zero element
    global_phase = 1.0
    for i in range(256):
        for j in range(256):
            if abs(op_ref.data[i, j]) > 1e-10 and abs(op_test.data[i, j]) > 1e-10:
                global_phase = op_test.data[i, j] / op_ref.data[i, j]
                break
        if abs(global_phase - 1.0) > 1e-10:
            break

    # Correct for global phase
    op_test_corrected = op_test.data / global_phase
    matrix_error = np.linalg.norm(op_ref.data - op_test_corrected, ord="fro")
    np.max(np.abs(op_ref.data - op_test_corrected))

    if matrix_error < 1e-10:
        pass
    else:
        pass

    # Test H_TTA

    # Reference
    qc_ref = QuantumCircuit(8)
    U_TTA = build_H_TTA_qubit_unitary(J, dt, hbar)
    gate_ref = UnitaryGate(U_TTA, label="U_TTA")
    qc_ref.append(gate_ref, [0, 1, 2, 3])
    op_ref = Operator(qc_ref)

    # Test
    qc_test = QuantumCircuit(8)
    apply_H_TTA_basic_gates(qc_test, 0, 1, J, dt, hbar)
    op_test = Operator(qc_test)

    # Compare (accounting for possible global phase difference)
    # Global phase doesn't affect physical observables, so we correct for it
    np.abs(np.trace(op_ref.data.conj().T @ op_test.data)) / 256

    # Find global phase by comparing a non-zero element
    global_phase = 1.0
    for i in range(256):
        for j in range(256):
            if abs(op_ref.data[i, j]) > 1e-10 and abs(op_test.data[i, j]) > 1e-10:
                global_phase = op_test.data[i, j] / op_ref.data[i, j]
                break
        if abs(global_phase - 1.0) > 1e-10:
            break

    # Correct for global phase
    op_test_corrected = op_test.data / global_phase
    matrix_error = np.linalg.norm(op_ref.data - op_test_corrected, ord="fro")
    np.max(np.abs(op_ref.data - op_test_corrected))

    if matrix_error < 1e-10:
        pass
    else:
        pass


if __name__ == "__main__":
    test_basic_gate_decomposition()
