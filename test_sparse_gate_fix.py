#!/usr/bin/env python3
"""Test script to verify the sparse gate decomposition fix.

This test verifies that:
1. Gate counting is dramatically improved (from ~7252 to ~88 per Trotter step)
2. The sparse structure recognition works correctly
3. Mathematical exactness is maintained (fidelity = 1.0)
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "tutorials")
sys.path.insert(0, "tools")

from integrated_sparse_compiler_v2 import IntegratedSparseCompilerV2
from mqt_qudits_four_molecule_sparse_implementation import PhysicalParameters, SparseAwareMQTQuditTimeEvolution

from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister


def test_h_tta_sparse_structure() -> bool:
    """Test that H_TTA unitary has the expected 3×3 sparse structure."""
    # Build H_TTA unitary
    J = 0.05
    dt = 10.0
    hbar = 0.6582

    omega = np.sqrt(2) * J * dt / hbar
    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)

    # 3×3 H_TTA unitary
    U_3x3 = np.array(
        [
            [0.5 * (1 + cos_omega), -1j * sin_omega / np.sqrt(2), -0.5 * (1 - cos_omega)],
            [-1j * sin_omega / np.sqrt(2), cos_omega, -1j * sin_omega / np.sqrt(2)],
            [-0.5 * (1 - cos_omega), -1j * sin_omega / np.sqrt(2), 0.5 * (1 + cos_omega)],
        ],
        dtype=complex,
    )

    # Embed into 9×9 at indices [2, 4, 6] (|02⟩, |11⟩, |20⟩)
    U_9x9 = np.eye(9, dtype=complex)
    active_indices = [2, 4, 6]
    for i, idx_i in enumerate(active_indices):
        for j, idx_j in enumerate(active_indices):
            U_9x9[idx_i, idx_j] = U_3x3[i, j]

    # Compile with sparse compiler
    compiler = IntegratedSparseCompilerV2(tolerance=1e-10, optimize_gates=True)
    result = compiler.compile(U_9x9)

    # Verify results
    assert result.structure_info.structure_type == "sparse_subspace", (
        f"Expected sparse_subspace, got {result.structure_info.structure_type}"
    )
    assert result.structure_info.active_dimension == 3, (
        f"Expected active dimension 3, got {result.structure_info.active_dimension}"
    )
    assert result.structure_info.active_subspace == [2, 4, 6], (
        f"Expected active subspace [2, 4, 6], got {result.structure_info.active_subspace}"
    )
    assert result.gate_count_estimate == 6, f"Expected 6 gates, got {result.gate_count_estimate}"
    assert result.fidelity > 0.9999999999, f"Expected fidelity ~1.0, got {result.fidelity}"

    return True


def test_gate_count_per_trotter_step() -> bool:
    """Test that a full Trotter step has the expected gate count."""
    # Create parameters
    params = PhysicalParameters()
    time_evol = SparseAwareMQTQuditTimeEvolution(params)

    # Build a circuit with one full symmetric Trotter step
    circuit = QuantumCircuit()
    reg = QuantumRegister("molecules", params.N_molecules, [3] * params.N_molecules)
    circuit.append(reg)

    dt = 5.0

    # Add one full symmetric Trotter step:
    # Forward: H0/2, H_transfer/2, H_TTA/2
    # Backward: H_TTA/2, H_transfer/2, H0/2

    # Forward
    time_evol.add_H0_evolution_gates(circuit, dt / 2)
    time_evol.add_H_transfer_evolution_gates(circuit, dt / 2)
    time_evol.add_H_TTA_evolution_gates(circuit, dt / 2)

    # Backward
    time_evol.add_H_TTA_evolution_gates(circuit, dt / 2)
    time_evol.add_H_transfer_evolution_gates(circuit, dt / 2)
    time_evol.add_H0_evolution_gates(circuit, dt / 2)

    # Count gates
    custom_two_count = sum(1 for gate in circuit.instructions if gate.__class__.__name__ == "CustomTwo")
    other_gates = sum(1 for gate in circuit.instructions if gate.__class__.__name__ != "CustomTwo")

    # Estimate gate count
    estimated_gates = time_evol.decompose_custom_two_gates(circuit)

    # Verify results
    # Expected:
    # - H0: 8 gates per half-step × 2 = 16
    # - H_transfer: 18 gates per half-step × 2 = 36
    # - H_TTA: 18 gates per half-step × 2 = 36
    # Total: 88 gates

    expected_range = (85, 95)  # Allow some tolerance
    assert expected_range[0] <= estimated_gates <= expected_range[1], (
        f"Expected {expected_range[0]}-{expected_range[1]} gates, got {estimated_gates}"
    )

    # Old implementation would have been ~7252 gates
    old_estimate = other_gates + custom_two_count * 1200
    reduction = (old_estimate - estimated_gates) / old_estimate * 100

    assert reduction > 95, f"Expected >95% reduction, got {reduction:.1f}%"

    return True


def test_mathematical_exactness() -> bool:
    """Test that the sparse compiler maintains mathematical exactness."""
    # Test with random sparse unitaries
    compiler = IntegratedSparseCompilerV2(tolerance=1e-10, optimize_gates=True)

    test_cases = []

    # Test case 1: 2×2 sparse structure
    U1 = np.eye(9, dtype=complex)
    angle = 0.3
    U1[1, 1] = np.cos(angle)
    U1[1, 2] = -1j * np.sin(angle)
    U1[2, 1] = -1j * np.sin(angle)
    U1[2, 2] = np.cos(angle)
    test_cases.append(("2×2 sparse", U1))

    # Test case 2: 3×3 sparse structure (H_TTA-like)
    J = 0.05
    dt = 10.0
    hbar = 0.6582
    omega = np.sqrt(2) * J * dt / hbar
    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)
    U_3x3 = np.array(
        [
            [0.5 * (1 + cos_omega), -1j * sin_omega / np.sqrt(2), -0.5 * (1 - cos_omega)],
            [-1j * sin_omega / np.sqrt(2), cos_omega, -1j * sin_omega / np.sqrt(2)],
            [-0.5 * (1 - cos_omega), -1j * sin_omega / np.sqrt(2), 0.5 * (1 + cos_omega)],
        ],
        dtype=complex,
    )
    U2 = np.eye(9, dtype=complex)
    active_indices = [2, 4, 6]
    for i, idx_i in enumerate(active_indices):
        for j, idx_j in enumerate(active_indices):
            U2[idx_i, idx_j] = U_3x3[i, j]
    test_cases.append(("3×3 sparse (H_TTA)", U2))

    all_passed = True
    for _name, U in test_cases:
        result = compiler.compile(U)
        fidelity = result.fidelity

        if fidelity < 0.9999999999:
            all_passed = False
        else:
            pass

    assert all_passed, "Some test cases failed fidelity check"

    return True


def main() -> int:
    """Run all tests."""
    tests = [test_h_tta_sparse_structure, test_gate_count_per_trotter_step, test_mathematical_exactness]

    passed = 0
    failed = 0

    for test in tests:
        try:
            if test():
                passed += 1
        except AssertionError:
            failed += 1
        except Exception:
            import traceback

            traceback.print_exc()
            failed += 1

    if failed == 0:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
