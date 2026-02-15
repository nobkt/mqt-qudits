"""Validation functions for GKSL-Lindblad quantum dynamics."""

from __future__ import annotations

import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import compute_purity, compute_von_neumann_entropy


class PhysicsViolationError(Exception):
    """Raised when a physical constraint is violated."""


def validate_density_matrix(
    rho: np.ndarray,
    step: int = 0,
    tolerance: dict[str, float] | None = None,
) -> dict:
    """Validate physical properties of a density matrix.

    Args:
        rho: Density matrix to validate.
        step: Time-step index (for error messages).
        tolerance: Override default tolerances for 'trace', 'hermiticity',
                   'positivity', and 'entropy'.

    Returns:
        Dict with keys: valid, trace, hermiticity_error, min_eigenvalue,
        entropy, purity, errors.

    Raises:
        PhysicsViolationError: If any physical constraint is violated.
    """
    default_tol = {
        "trace": 1e-8,
        "hermiticity": 1e-10,
        "positivity": -1e-10,
        "entropy": -1e-8,
    }
    if tolerance is not None:
        default_tol.update(tolerance)
    tol = default_tol

    errors: list[str] = []

    trace_val = float(np.real(np.trace(rho)))
    if abs(trace_val - 1.0) > tol["trace"]:
        errors.append(f"Step {step}: Trace = {trace_val} (expected 1.0)")

    herm_err = float(np.linalg.norm(rho - rho.conj().T, "fro"))
    if herm_err > tol["hermiticity"]:
        errors.append(f"Step {step}: Hermiticity error = {herm_err}")

    eigenvalues = np.linalg.eigvalsh(rho)
    min_eig = float(eigenvalues[0])
    if min_eig < tol["positivity"]:
        errors.append(f"Step {step}: Min eigenvalue = {min_eig} (negative)")

    entropy = compute_von_neumann_entropy(rho)
    if entropy < tol["entropy"]:
        errors.append(f"Step {step}: Entropy = {entropy} (negative)")

    purity = compute_purity(rho)

    valid = len(errors) == 0
    result = {
        "valid": valid,
        "trace": trace_val,
        "hermiticity_error": herm_err,
        "min_eigenvalue": min_eig,
        "entropy": entropy,
        "purity": purity,
        "errors": errors,
    }

    if not valid:
        raise PhysicsViolationError("; ".join(errors))

    return result


def validate_particle_conservation(
    populations: dict[str, float],
    N_molecules: float,
    tolerance: float = 1e-8,
) -> bool:
    """Check that N_S0 + N_T1 + N_S1 = N_molecules.

    Args:
        populations: Dict with keys 'N_S0', 'N_T1', 'N_S1'.
        N_molecules: Expected total number of molecules.
        tolerance: Allowed deviation.

    Returns:
        True if conservation holds.

    Raises:
        PhysicsViolationError: If particle number is not conserved.
    """
    total = populations["N_S0"] + populations["N_T1"] + populations["N_S1"]
    if abs(total - N_molecules) > tolerance:
        msg = f"Particle conservation violated: N_S0 + N_T1 + N_S1 = {total}, expected {N_molecules}"
        raise PhysicsViolationError(msg)
    return True


def validate_entropy_increase(
    S_prev: float,
    S_curr: float,
    tolerance: float = -1e-8,
) -> bool:
    """Warn if entropy decreased beyond tolerance.

    Args:
        S_prev: Entropy at previous time step.
        S_curr: Entropy at current time step.
        tolerance: Allowed decrease (negative value).

    Returns:
        True if entropy did not decrease significantly.
    """
    diff = S_curr - S_prev
    if diff < tolerance:
        warnings.warn(
            f"Entropy decreased: S_prev={S_prev:.6e}, S_curr={S_curr:.6e}, dS={diff:.6e}",
            stacklevel=2,
        )
        return False
    return True
