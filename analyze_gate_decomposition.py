#!/usr/bin/env python3
"""Analyze if the current H_TTA gate decomposition could possibly be correct.

Since we can't easily extract the unitary from a circuit, we'll analyze
the structure of the gates and compare with what we know must be true.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm


def analyze_correct_unitary():
    """Analyze the structure of the correct unitary matrix."""
    J = 0.05
    dt = 10.0
    hbar = 0.6582

    # Correct Hamiltonian
    H_TTA = J * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])

    # Correct unitary
    U = expm(-1j * H_TTA * dt / hbar)

    for i in range(3):
        for j in range(3):
            real_part = np.real(U[i, j])
            imag_part = np.imag(U[i, j])
            np.abs(U[i, j])

            if (abs(real_part) > 1e-10 and abs(imag_part) > 1e-10) or abs(real_part) > 1e-10 or abs(imag_part) > 1e-10:
                pass
            else:
                pass

    omega = np.sqrt(2) * J * dt / hbar
    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)

    U_analytical = np.array([
        [0.5 * (1 + cos_omega), -1j * sin_omega / np.sqrt(2), -0.5 * (1 - cos_omega)],
        [-1j * sin_omega / np.sqrt(2), cos_omega, -1j * sin_omega / np.sqrt(2)],
        [-0.5 * (1 - cos_omega), -1j * sin_omega / np.sqrt(2), 0.5 * (1 + cos_omega)],
    ])

    formula_error = np.linalg.norm(U - U_analytical)

    if formula_error < 1e-10:
        pass
    else:
        pass

    return U


def analyze_gate_sequence_structure() -> None:
    """Analyze what the current gate sequence is trying to do."""


def check_if_verification_exists() -> None:
    """Check if there's any verification of the gate sequence."""


if __name__ == "__main__":
    U_correct = analyze_correct_unitary()
    analyze_gate_sequence_structure()
    check_if_verification_exists()
