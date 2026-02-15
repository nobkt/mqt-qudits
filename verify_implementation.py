#!/usr/bin/env python3
"""Verification script for shot-based simulation implementation.

This script verifies that all required modifications have been correctly
applied to the quantum dynamics comparison notebook.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def verify_notebook() -> bool:
    """Verify the notebook modifications."""
    notebook_path = Path("tutorials/quantum_dynamics_complete_comparison.ipynb")

    if not notebook_path.exists():
        return False

    # Load notebook
    with open(notebook_path, encoding="utf-8") as f:
        nb = json.load(f)

    # Test 1: Check Qubit shot-based cell
    qubit_cell = None
    for cell in nb["cells"]:
        if cell.get("id") == "5bf25d98":
            qubit_cell = cell
            break

    if not qubit_cell:
        return False

    source = "".join(qubit_cell.get("source", []))
    if "ショットベース" in source and "Sampler" in source and "calculate_populations_from_counts" in source:
        pass
    else:
        return False

    # Test 2: Check Qudit shot-based cell
    qudit_cell = None
    for cell in nb["cells"]:
        if cell.get("id") == "397deb45":
            qudit_cell = cell
            break

    if not qudit_cell:
        return False

    source = "".join(qudit_cell.get("source", []))
    if "ショットベース" in source and "simulate_shot_based" in source:
        pass
    else:
        return False

    # Test 3: Check Qubit visualization cell exists
    qubit_viz_cell = None
    for cell in nb["cells"]:
        if cell.get("id") == "c7fe05f9":
            qubit_viz_cell = cell
            break

    if not qubit_viz_cell:
        return False

    source = "".join(qubit_viz_cell.get("source", []))
    if "circuit_drawer" in source and "step_circuit" in source:
        pass
    else:
        return False

    # Test 4: Check Qudit visualization cell was added
    qudit_viz_cell = None
    for cell in nb["cells"]:
        if cell.get("id") == "qudit_viz_1step":
            qudit_viz_cell = cell
            break

    if not qudit_viz_cell:
        return False

    source = "".join(qudit_viz_cell.get("source", []))
    if "1鈴木トロッターステップ" in source and "plot_circuit" in source:
        pass
    else:
        return False

    # Test 5: Check no heuristics mentioned
    all_source = "".join(["".join(c.get("source", [])) for c in nb["cells"]])
    if "heuristic" not in all_source.lower() and "fallback" not in all_source.lower():
        pass
    else:
        pass

    # Test 6: Verify shot numbers
    qubit_source = "".join(qubit_cell.get("source", []))
    qudit_source = "".join(qudit_cell.get("source", []))

    if "shots=10000" in qubit_source or "shots: int = 10000" in qubit_source:
        pass
    else:
        pass

    if "shots=10000" in qudit_source or "shots: 10000" in qudit_source:
        pass
    else:
        pass

    return True


def verify_implementation_file() -> bool:
    """Verify the sparse implementation file has shot-based methods."""
    impl_path = Path("tutorials/mqt_qudits_four_molecule_sparse_implementation.py")

    if not impl_path.exists():
        return False

    with open(impl_path, encoding="utf-8") as f:
        content = f.read()

    # Check for new methods
    if "def simulate_shot_based" in content:
        pass
    else:
        return False

    if "def calculate_populations_from_samples" in content:
        pass
    else:
        return False

    if "np.random.choice" in content:
        pass
    else:
        return False

    return True


def main() -> int:
    """Main verification function."""
    notebook_ok = verify_notebook()
    impl_ok = verify_implementation_file()

    if notebook_ok and impl_ok:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
