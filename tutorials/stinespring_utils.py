"""Stinespring dilation utilities for GKSL-Lindblad quantum dynamics."""

from __future__ import annotations

import os
import sys
from collections.abc import Callable

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import unvectorize_density_matrix, vectorize_density_matrix


def stinespring_unitary_from_lindblad(L: np.ndarray, dt: float) -> np.ndarray:
    """Build Stinespring unitary from Lindblad operator L (which already includes sqrt(gamma)).

    The generator is the block matrix G = [[0, L], [L†, 0]] and the unitary is
    U = expm(-i * theta * G) with theta = sqrt(dt).
    """
    d_sys = L.shape[0]
    L = np.asarray(L, dtype=np.complex128)

    G = np.zeros((2 * d_sys, 2 * d_sys), dtype=np.complex128)
    G[:d_sys, d_sys:] = L
    G[d_sys:, :d_sys] = L.conj().T

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

    rho_ext = np.kron(rho, env0)
    rho_prime = U @ rho_ext @ U.conj().T

    # Partial trace over environment (2-level) using block structure
    rho_out = rho_prime[:d_sys, :d_sys] + rho_prime[d_sys:, d_sys:]

    return rho_out


def build_gksl_superoperator(
    H_total: np.ndarray,
    lindblad_ops: list,
    hbar: float = 1.0,
) -> np.ndarray:
    """Build the full GKSL Liouvillian superoperator in vectorized (column-major) form.

    L_total = L_H + L_D where:
      L_H = -i/hbar * (H ⊗ I - I ⊗ H^T)
      L_D = sum_alpha [ L_alpha ⊗ conj(L_alpha)
                        - 0.5*(L†_alpha L_alpha ⊗ I + I ⊗ (L†_alpha L_alpha)^T) ]

    lindblad_ops: list of (L_alpha, gamma_alpha) tuples.
    L_alpha already contains sqrt(gamma) factor.
    """
    dim = H_total.shape[0]
    I = np.eye(dim, dtype=np.complex128)

    # Hamiltonian part
    L_H = (-1j / hbar) * (np.kron(H_total, I) - np.kron(I, H_total.T))

    # Dissipator part
    L_D = np.zeros((dim * dim, dim * dim), dtype=np.complex128)
    for item in lindblad_ops:
        if isinstance(item, tuple):
            L_op = np.asarray(item[0], dtype=np.complex128)
        else:
            L_op = np.asarray(item, dtype=np.complex128)
        LdL = L_op.conj().T @ L_op
        L_D += (
            np.kron(L_op, L_op.conj())
            - 0.5 * np.kron(LdL, I)
            - 0.5 * np.kron(I, LdL.T)
        )

    return L_H + L_D


def build_trotter_step_classical(
    H_0: np.ndarray,
    H_transfer: np.ndarray,
    lindblad_ops: list,
    dt: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """Build a 2nd-order symmetric Trotter step function.

    exp(L dt) ≈ exp(L_H dt/2) exp(L_D dt) exp(L_H dt/2)

    Returns a function trotter_step(rho) -> rho'.
    lindblad_ops: list of (L_alpha, gamma_alpha) tuples or plain arrays.
    """
    dim = H_0.shape[0]
    I = np.eye(dim, dtype=np.complex128)
    H_total = np.asarray(H_0, dtype=np.complex128) + np.asarray(H_transfer, dtype=np.complex128)

    # Hamiltonian superoperator
    L_H = (-1j) * (np.kron(H_total, I) - np.kron(I, H_total.T))

    # Dissipator superoperator
    L_D = np.zeros((dim * dim, dim * dim), dtype=np.complex128)
    for item in lindblad_ops:
        if isinstance(item, tuple):
            L_op = np.asarray(item[0], dtype=np.complex128)
        else:
            L_op = np.asarray(item, dtype=np.complex128)
        LdL = L_op.conj().T @ L_op
        L_D += (
            np.kron(L_op, L_op.conj())
            - 0.5 * np.kron(LdL, I)
            - 0.5 * np.kron(I, LdL.T)
        )

    # Precompute matrix exponentials
    U_H_half = expm(L_H * (dt / 2.0))
    U_D = expm(L_D * dt)

    def trotter_step(rho: np.ndarray) -> np.ndarray:
        vec = vectorize_density_matrix(rho)
        vec = U_H_half @ vec
        vec = U_D @ vec
        vec = U_H_half @ vec
        return unvectorize_density_matrix(vec, dim)

    return trotter_step
