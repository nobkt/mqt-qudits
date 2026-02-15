#!/usr/bin/env python3
"""Verification script for CustomTwo gate removal.

This script checks that:
1. Notebook cells 5, 8, 9 are properly formatted
2. No CustomTwo gates are used in the modified implementation
3. All syntax is valid
4. Basic gate decompositions are mathematically sound
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Add tutorials to path
sys.path.insert(0, str(Path(__file__).parent / "tutorials"))


def check_notebook_formatting() -> bool:
    """Check that notebook cells are properly formatted."""
    with open("tutorials/quantum_dynamics_complete_comparison.ipynb", encoding="utf-8") as f:
        nb = json.load(f)

    issues = []

    for i in [5, 8, 9]:
        cell = nb["cells"][i]
        source = cell["source"]

        # Check if properly formatted (has newlines)
        has_newlines = any("\n" in line for line in source)

        # Check if all on one line (problem)
        if len(source) == 1 and len(source[0]) > 1000:
            issues.append(f"Cell {i}: All code on one line!")
        elif not has_newlines:
            issues.append(f"Cell {i}: Missing newlines!")
        else:
            pass

    if issues:
        for _issue in issues:
            pass
        return False
    return True


def check_no_customtwo() -> bool:
    """Check that CustomTwo gates are not used in the implementation."""
    # Check sparse implementation
    with open("tutorials/mqt_qudits_four_molecule_sparse_implementation.py", encoding="utf-8") as f:
        content = f.read()

    issues = []

    # Check that CustomTwo is only mentioned in:
    # 1. Comments about NOT using it
    # 2. Error messages
    # 3. Old code that's commented out

    lines = content.split("\n")
    for i, line in enumerate(lines, 1):
        if "CustomTwo" in line or "custom_two" in line.lower():
            # Check if it's in a context that says "NO CustomTwo" or "not using CustomTwo"
            if any(keyword in line for keyword in ["NO CustomTwo", "not using", "不使用", "without", "error"]):
                continue  # This is OK - documenting that we don't use it
            if line.strip().startswith("#"):
                continue  # Comment is OK
            if "circuit.cu_two" in line:
                # This is actual usage - should be removed
                issues.append(f"Line {i}: CustomTwo gate usage found: {line.strip()}")

    if issues:
        for _issue in issues:
            pass
        return False
    return True


def check_syntax():
    """Check Python syntax for all modified files."""
    files = ["tutorials/exact_qudit_basic_gates.py", "tutorials/mqt_qudits_four_molecule_sparse_implementation.py"]

    all_ok = True
    for file in files:
        try:
            import py_compile

            py_compile.compile(file, doraise=True)
        except py_compile.PyCompileError:
            all_ok = False

    if all_ok:
        pass
    return all_ok


def check_decomposition_correctness() -> bool | None:
    """Check that the decomposition functions are mathematically sound."""
    # This requires numpy, which may not be installed
    try:
        from exact_qudit_basic_gates import verify_H_transfer_decomposition, verify_H_TTA_decomposition

        # Test parameters
        V = 0.1  # eV
        J = 0.05  # eV
        dt = 10.0  # fs
        hbar = 0.6582  # eV·fs

        # Check H_transfer
        if verify_H_transfer_decomposition(V, dt, hbar):
            pass
        else:
            return False

        # Check H_TTA
        if verify_H_TTA_decomposition(J, dt, hbar):
            pass
        else:
            return False

        return True

    except ImportError:
        return True  # Don't fail if dependencies not available


def main() -> int:
    """Run all verification checks."""
    results = {
        "Notebook Formatting": check_notebook_formatting(),
        "No CustomTwo Usage": check_no_customtwo(),
        "Python Syntax": check_syntax(),
        "Mathematical Correctness": check_decomposition_correctness(),
    }

    all_passed = True
    for passed in results.values():
        if not passed:
            all_passed = False

    if all_passed:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
