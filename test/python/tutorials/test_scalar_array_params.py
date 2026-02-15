#!/usr/bin/env python3
"""Test that the sparse implementation handles both scalar and array parameters.

This addresses the bug in quantum_dynamics_complete_comparison.ipynb where
the notebook uses scalar V and J, but the sparse implementation expected arrays.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Add tutorials to path
tutorials_path = Path(__file__).parent.parent.parent.parent / "tutorials"
sys.path.insert(0, str(tutorials_path))

from mqt_qudits_four_molecule_sparse_implementation import (
    ExactDiagonalizationSolver,
    PhysicalParameters,
    SuzukiTrotterMQTQuditSimulator,
)


def test_scalar_parameters():
    """Test that scalar V and J parameters work (notebook style)."""

    class ScalarParams:
        """Parameters with scalar V and J like in the notebook."""

        def __init__(self):
            self.N_molecules = 4
            self.E_T = 1.5
            self.E_S = 3.0
            self.V = 0.1  # scalar!
            self.J = 0.05  # scalar!
            self.Gamma_fl = 0.01
            self.hbar = 0.6582119569
            self.neighbors = [(0, 1), (1, 2), (2, 3)]

    params = ScalarParams()

    # This should not raise TypeError about subscripting a float
    # Test with Suzuki-Trotter simulator
    sim = SuzukiTrotterMQTQuditSimulator(params)

    # Test building Trotter step unitary (this uses V and J subscripting)
    dt = 10.0
    U = sim.build_trotter_step_unitary_direct(dt)
    assert U.shape == (81, 81), "Unitary should be 81x81"

    # Verify unitarity
    unitarity_error = np.linalg.norm(U @ U.conj().T - np.eye(81))
    assert unitarity_error < 1e-10, f"Unitary error too large: {unitarity_error}"

    # Test with exact diagonalization solver
    solver = ExactDiagonalizationSolver(params)
    H = solver.build_hamiltonian()
    assert H.shape == (81, 81), "Hamiltonian should be 81x81 for 4 qutrits"


def test_array_parameters():
    """Test that array V and J parameters still work (backward compatibility)."""
    # Use the standard PhysicalParameters class which has array V and J
    params = PhysicalParameters()

    # This should work as before
    sim = SuzukiTrotterMQTQuditSimulator(params)

    # Test building Trotter step unitary (this uses V and J subscripting)
    dt = 10.0
    U = sim.build_trotter_step_unitary_direct(dt)
    assert U.shape == (81, 81), "Unitary should be 81x81"

    # Verify unitarity
    unitarity_error = np.linalg.norm(U @ U.conj().T - np.eye(81))
    assert unitarity_error < 1e-10, f"Unitary error too large: {unitarity_error}"

    # Test with exact diagonalization solver
    solver = ExactDiagonalizationSolver(params)
    H = solver.build_hamiltonian()
    assert H.shape == (81, 81), "Hamiltonian should be 81x81 for 4 qutrits"


def test_scalar_and_array_give_same_results():
    """Test that scalar and array params with same values give identical results."""

    # Scalar params
    class ScalarParams:
        def __init__(self):
            self.N_molecules = 4
            self.E_T = 1.5
            self.E_S = 3.0
            self.V = 0.1
            self.J = 0.05
            self.Gamma_fl = 0.01
            self.hbar = 0.6582119569
            self.neighbors = [(0, 1), (1, 2), (2, 3)]

    # Array params with same values
    class ArrayParams:
        def __init__(self):
            self.N_molecules = 4
            self.E_T = 1.5
            self.E_S = 3.0
            self.V = np.array([0.1, 0.1, 0.1])
            self.J = np.array([0.05, 0.05, 0.05])
            self.Gamma_fl = 0.01
            self.hbar = 0.6582119569
            self.neighbors = [(0, 1), (1, 2), (2, 3)]

    scalar_params = ScalarParams()
    array_params = ArrayParams()

    # Test with exact diagonalization solver
    solver_scalar = ExactDiagonalizationSolver(scalar_params)
    solver_array = ExactDiagonalizationSolver(array_params)

    # Build Hamiltonians
    H_scalar = solver_scalar.build_hamiltonian()
    H_array = solver_array.build_hamiltonian()

    # They should be identical
    diff = np.linalg.norm(H_scalar - H_array)
    assert diff < 1e-14, f"Hamiltonians differ: {diff}"

    # Test with Suzuki-Trotter simulator
    sim_scalar = SuzukiTrotterMQTQuditSimulator(scalar_params)
    sim_array = SuzukiTrotterMQTQuditSimulator(array_params)

    # Build Trotter step unitaries
    dt = 10.0
    U_scalar = sim_scalar.build_trotter_step_unitary_direct(dt)
    U_array = sim_array.build_trotter_step_unitary_direct(dt)

    # They should be identical
    diff = np.linalg.norm(U_scalar - U_array)
    assert diff < 1e-14, f"Unitaries differ: {diff}"


if __name__ == "__main__":
    print("Testing scalar parameters...")
    test_scalar_parameters()
    print("✓ Scalar parameters work correctly")

    print("\nTesting array parameters...")
    test_array_parameters()
    print("✓ Array parameters work correctly")

    print("\nTesting that scalar and array give same results...")
    test_scalar_and_array_give_same_results()
    print("✓ Scalar and array parameters give identical results")

    print("\n" + "=" * 70)
    print("All tests passed!")
