"""Stinespring dilation utilities for GKSL-Lindblad quantum dynamics."""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def stinespring_unitary_from_lindblad(L: np.ndarray, dt: float) -> np.ndarray:
    """Build Stinespring unitary from Lindblad operator L (which already includes sqrt(gamma)).

    The generator is the block matrix G = [[0, L†], [L, 0]] and the unitary is
    U = expm(-i * theta * G) with theta = sqrt(dt).

    At leading order this implements the GKSL dissipator:
    E(ρ) ≈ ρ + dt * (L ρ L† − ½{L†L, ρ})
    """
    d_sys = L.shape[0]
    L = np.asarray(L, dtype=np.complex128)

    G = np.zeros((2 * d_sys, 2 * d_sys), dtype=np.complex128)
    G[:d_sys, d_sys:] = L.conj().T
    G[d_sys:, :d_sys] = L

    theta = np.sqrt(dt)
    U = expm(-1j * theta * G)

    residual = np.linalg.norm(U.conj().T @ U - np.eye(2 * d_sys), ord="fro")
    if residual >= 1e-10:
        msg = f"Stinespring unitary failed unitarity check: ||U†U - I||_F = {residual}"
        raise ValueError(msg)

    return U


def apply_stinespring_to_density_matrix(rho: np.ndarray, U: np.ndarray) -> np.ndarray:
    """Apply Stinespring unitary to a system density matrix.

    Extends rho with an environment qubit in |0>, applies U, then partial-traces
    over the environment.
    """
    d_sys = rho.shape[0]
    rho = np.asarray(rho, dtype=np.complex128)

    # Environment: |0><0| for a 2-level ancilla
    env0 = np.zeros((2, 2), dtype=np.complex128)
    env0[0, 0] = 1.0

    # Use kron(env, system) ordering to match Stinespring unitary block structure
    rho_ext = np.kron(env0, rho)
    rho_prime = U @ rho_ext @ U.conj().T

    # Partial trace over environment (2-level) using block structure
    rho_out = rho_prime[:d_sys, :d_sys] + rho_prime[d_sys:, d_sys:]

    return rho_out


def build_gksl_superoperator(
    H_total: np.ndarray,
    lindblad_ops: list,
    hbar: float = 1.0,
) -> np.ndarray:
    """Build the full GKSL Liouvillian superoperator in column-major vectorized form.

    Uses the identity vec(AXB) = (B^T ⊗ A) vec(X) for column-major (Fortran-order)
    vectorization, consistent with vectorize_density_matrix / unvectorize_density_matrix.

    L_total = L_H + L_D where:
      L_H = -i/hbar * (I ⊗ H - H^T ⊗ I)
      L_D = sum_alpha [ conj(L_alpha) ⊗ L_alpha
                        - 0.5*(I ⊗ L†_alpha L_alpha + (L†_alpha L_alpha)^T ⊗ I) ]

    lindblad_ops: list of (L_alpha, gamma_alpha) tuples.
    L_alpha already contains sqrt(gamma) factor.
    """
    dim = H_total.shape[0]
    I = np.eye(dim, dtype=np.complex128)

    # Hamiltonian part: vec(-i[H,rho]/hbar) = (-i/hbar)(I⊗H - H^T⊗I) vec(rho)
    L_H = (-1j / hbar) * (np.kron(I, H_total) - np.kron(H_total.T, I))

    # Dissipator part
    L_D = np.zeros((dim * dim, dim * dim), dtype=np.complex128)
    for item in lindblad_ops:
        if isinstance(item, tuple):
            L_op = np.asarray(item[0], dtype=np.complex128)
        else:
            L_op = np.asarray(item, dtype=np.complex128)
        LdL = L_op.conj().T @ L_op
        L_D += (
            np.kron(L_op.conj(), L_op)
            - 0.5 * np.kron(I, LdL)
            - 0.5 * np.kron(LdL.T, I)
        )

    return L_H + L_D


