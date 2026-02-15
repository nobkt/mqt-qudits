#!/usr/bin/env python3
"""Test that verifies the fix for the AssertionError.

This test ensures that:
1. Multi-qudit subspaces generate CustomTwo gates (not R gates with invalid indices)
2. Single-qudit subspaces generate R gates with proper local indices
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "tutorials"))
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "tools"))

try:
    import numpy as np

    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("NumPy not available, skipping tests")
    sys.exit(0)

from mqt_qudits_four_molecule_sparse_implementation import SparseAwareMQTGateGenerator


def test_h_transfer_generates_customtwo():
    """Test that H_transfer (multi-qudit |01⟩↔|10⟩) generates CustomTwo gate.

    This is the case that was causing the AssertionError before the fix.
    """
    print("Test: H_transfer generates CustomTwo gate")
    print("-" * 70)

    gate_generator = SparseAwareMQTGateGenerator()

    # H_transfer: rotation between |01⟩ and |10⟩
    theta = 0.1
    U_transfer = np.eye(9, dtype=complex)
    U_transfer[1, 1] = np.cos(theta)
    U_transfer[1, 3] = -1j * np.sin(theta)
    U_transfer[3, 1] = -1j * np.sin(theta)
    U_transfer[3, 3] = np.cos(theta)

    gate_info = gate_generator.compile_unitary_to_gates(U_transfer, [0, 1])

    # Should detect as sparse_2x2 structure
    assert gate_info["structure_type"] == "sparse_2x2", f"Expected sparse_2x2, got {gate_info['structure_type']}"

    # Should generate CustomTwo gate (not R gates with invalid indices)
    assert len(gate_info["gates"]) == 1, f"Expected 1 CustomTwo gate, got {len(gate_info['gates'])} gates"

    gate = gate_info["gates"][0]
    assert gate["type"] == "CustomTwo", f"Expected CustomTwo gate, got {gate['type']}"

    # The unitary should be preserved
    assert "unitary" in gate["params"], "CustomTwo should have unitary parameter"
    assert np.allclose(gate["params"]["unitary"], U_transfer), "Unitary matrix should be preserved"

    print(f"  ✓ Structure type: {gate_info['structure_type']}")
    print(f"  ✓ Gate type: {gate['type']}")
    print(f"  ✓ Fidelity: {gate_info['fidelity']:.10f}")
    print("  ✓ No R gates with invalid indices (AssertionError avoided!)")
    print()


def test_h_tta_generates_customtwo():
    """Test that H_TTA (multi-qudit 3x3 subspace) generates CustomTwo gate."""
    print("Test: H_TTA generates CustomTwo gate")
    print("-" * 70)

    gate_generator = SparseAwareMQTGateGenerator()

    # H_TTA: rotation in |02⟩, |11⟩, |20⟩ subspace
    U_tta = np.eye(9, dtype=complex)
    # Simplified version - just make it non-trivial in the 3x3 subspace
    for i, j in [(2, 4), (4, 6)]:
        U_tta[i, j] = 0.1j
        U_tta[j, i] = -0.1j

    gate_info = gate_generator.compile_unitary_to_gates(U_tta, [0, 1])

    # Should detect as sparse_3x3 structure
    assert gate_info["structure_type"] == "sparse_3x3", f"Expected sparse_3x3, got {gate_info['structure_type']}"

    # Should generate CustomTwo gate for multi-qudit subspace
    # (Even though it's 3x3, it still spans multiple qudits)
    found_customtwo = any(g["type"] == "CustomTwo" for g in gate_info["gates"])
    assert found_customtwo, "Expected at least one CustomTwo gate for multi-qudit 3x3 subspace"

    # Verify no R gates have invalid indices
    for gate in gate_info["gates"]:
        if gate["type"] in {"R", "Rz", "Rh"} and "level1" in gate["params"]:
            level1 = gate["params"]["level1"]
            level2 = gate["params"]["level2"]
            assert level1 < 3 and level2 < 3, (
                f"R gate has invalid indices: level1={level1}, level2={level2} (should be < 3)"
            )

    print(f"  ✓ Structure type: {gate_info['structure_type']}")
    print(f"  ✓ Contains CustomTwo: {found_customtwo}")
    print(f"  ✓ Fidelity: {gate_info['fidelity']:.10f}")
    print("  ✓ All R gates have valid indices (< 3)")
    print()


def test_single_qudit_subspace():
    """Test that single-qudit subspaces generate R gates with proper local indices.

    This verifies the local index conversion works correctly.
    """
    print("Test: Single-qudit subspace generates proper R gates")
    print("-" * 70)

    gate_generator = SparseAwareMQTGateGenerator()

    # Create a 2x2 subspace within a single qudit
    # E.g., rotation in |10⟩↔|11⟩ (only qudit 1 varies, qudit 0 stays at 1)
    theta = 0.1
    U_single = np.eye(9, dtype=complex)
    U_single[3, 3] = np.cos(theta)
    U_single[3, 4] = -1j * np.sin(theta)
    U_single[4, 3] = -1j * np.sin(theta)
    U_single[4, 4] = np.cos(theta)

    gate_info = gate_generator.compile_unitary_to_gates(U_single, [0, 1])

    # Should detect as sparse_2x2
    assert gate_info["structure_type"] == "sparse_2x2", f"Expected sparse_2x2, got {gate_info['structure_type']}"

    # Should NOT use CustomTwo for single-qudit subspace
    customtwo_count = sum(1 for g in gate_info["gates"] if g["type"] == "CustomTwo")
    assert customtwo_count == 0, f"Single-qudit subspace should not use CustomTwo, but found {customtwo_count}"

    # All R gates should have valid local indices (< 3)
    for gate in gate_info["gates"]:
        if gate["type"] in {"R", "Rz", "Rh"} and "level1" in gate["params"]:
            level1 = gate["params"]["level1"]
            level2 = gate["params"]["level2"]
            assert 0 <= level1 < 3, f"level1={level1} is invalid (should be 0-2)"
            assert 0 <= level2 < 3, f"level2={level2} is invalid (should be 0-2)"
            assert level1 != level2, "level1 and level2 should be different"
            print(f"  ✓ R gate with local indices: [{level1}, {level2}]")

    print(f"  ✓ Structure type: {gate_info['structure_type']}")
    print("  ✓ No CustomTwo gates (single-qudit optimization)")
    print("  ✓ All indices are valid local qudit levels")
    print()


if __name__ == "__main__":
    if not NUMPY_AVAILABLE:
        sys.exit(0)

    print("=" * 70)
    print("Testing gate generation fix for AssertionError")
    print("=" * 70)
    print()

    try:
        test_h_transfer_generates_customtwo()
        test_h_tta_generates_customtwo()
        test_single_qudit_subspace()

        print("=" * 70)
        print("All gate generation tests passed! ✓")
        print("=" * 70)
        print()
        print("Summary:")
        print("  ✓ Multi-qudit subspaces use CustomTwo (avoids AssertionError)")
        print("  ✓ Single-qudit subspaces use R gates with local indices")
        print("  ✓ Mathematical rigor preserved (exact unitary matrix)")

    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
