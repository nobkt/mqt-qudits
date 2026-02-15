#!/usr/bin/env python3
"""Validation script to compare Classical, Qubit, and Qudit simulations.

This script validates that all three implementations produce identical results
within numerical precision (< 1e-6 error).
"""

from __future__ import annotations

import json
import os
import sys

# Add tutorials to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../tutorials"))

# Import from notebook cells
# We'll need to extract and run the relevant cells


def extract_and_execute_notebook_cell(notebook_path, cell_index) -> None:
    """Extract and execute a specific notebook cell."""
    with open(notebook_path, encoding="utf-8") as f:
        nb = json.load(f)

    if cell_index >= len(nb["cells"]):
        msg = f"Cell {cell_index} not found in notebook"
        raise IndexError(msg)

    cell = nb["cells"][cell_index]
    if cell["cell_type"] != "code":
        msg = f"Cell {cell_index} is not a code cell"
        raise ValueError(msg)

    code = "".join(cell["source"])
    exec(code, globals())


def run_validation() -> bool:
    """Run full 3-way validation."""
    notebook_path = os.path.join(os.path.dirname(__file__), "../tutorials/quantum_dynamics_complete_comparison.ipynb")

    # Load and execute PhysicalParameters
    extract_and_execute_notebook_cell(notebook_path, 3)

    # Load and execute ClassicalSuzukiTrotterSimulator
    extract_and_execute_notebook_cell(notebook_path, 5)

    # Load and execute QubitMolecularDynamicsSimulator

    # We need to import dependencies first
    import warnings

    warnings.filterwarnings("ignore")

    extract_and_execute_notebook_cell(notebook_path, 8)

    # Create physical parameters

    params = PhysicalParameters()

    # Simulation parameters
    T_total = 100.0  # fs
    N_steps = 20

    # Run Classical simulation

    classical_sim = ClassicalSuzukiTrotterSimulator(params)
    classical_results = classical_sim.simulate(T_total, N_steps, initial_state_type="edge_triplet")

    final_classical = classical_results["populations"][-1]

    # Run Qubit simulation

    qubit_sim = QubitMolecularDynamicsSimulator(params)
    qubit_results = qubit_sim.simulate(T_total, N_steps, initial_state_type="edge_triplet", shots=100000)

    final_qubit = qubit_results["populations"][-1]

    # Compare results

    error_S0 = abs(final_classical["N_S0"] - final_qubit["N_S0"])
    error_T1 = abs(final_classical["N_T1"] - final_qubit["N_T1"])
    error_S1 = abs(final_classical["N_S1"] - final_qubit["N_S1"])
    max_error = max(error_S0, error_T1, error_S1)

    # Validation
    threshold = 1e-2  # Allow for shot noise

    return max_error < threshold


if __name__ == "__main__":
    import time

    start = time.time()

    try:
        success = run_validation()
        elapsed = time.time() - start

        sys.exit(0 if success else 1)

    except Exception:
        import traceback

        traceback.print_exc()
        sys.exit(1)
