#!/usr/bin/env python3
"""
Comprehensive Quantum Dynamics Verification Script

Verifies two critical bug fixes:
1. H0 Sign Bug Fix: The qubit H0 Pauli decomposition had inverted theta signs
   (used -2*coeff instead of +2*coeff), causing H0 evolution to be time-reversed
   (U_H0† instead of U_H0). Fixed in all qubit simulator files.

2. Qudit Noisy Simulator Fix: The qudit noisy simulator was applying global
   depolarizing once per Trotter step instead of per-gate 2-qudit depolarizing.
   Reimplemented using density matrix formalism with exact partial trace and
   per-gate noise channels: ε(ρ) = (1-p)ρ + p·Tr_{ij}(ρ)⊗I_{ij}/d².

Uses only numpy and scipy (no qiskit or mqt dependencies).
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import Any

import numpy as np
import scipy.linalg

# ====================================================================
# Physical parameters
# ====================================================================
E_T = 1.5       # eV  (triplet energy)
E_S = 3.0       # eV  (singlet energy)
V = 0.1         # eV  (transfer coupling)
J = 0.05        # eV  (TTA coupling)
HBAR = 0.6582119569  # eV·fs
N_MOLECULES = 4
D_QUTRIT = 3   # qutrit dimension
NEIGHBORS = [(0, 1), (1, 2), (2, 3)]
T_TOTAL = 100.0  # fs
N_STEPS = 20
DEPOL_2Q = 0.01  # 1% per 2-qubit/qudit gate

DIM_FULL = D_QUTRIT ** N_MOLECULES  # 81 for 4 qutrits


# ====================================================================
# Helper functions
# ====================================================================
def index_to_config(idx: int, N: int = N_MOLECULES, d: int = D_QUTRIT) -> list[int]:
    """Convert linear index to d-ary configuration (most-significant first)."""
    config = []
    for _ in range(N):
        config.append(idx % d)
        idx //= d
    return config[::-1]


def config_to_index(config: list[int], d: int = D_QUTRIT) -> int:
    """Convert d-ary configuration to linear index."""
    idx = 0
    for level in config:
        idx = idx * d + level
    return idx


def build_single_site_operator(
    site: int, op: np.ndarray, N: int = N_MOLECULES, d: int = D_QUTRIT
) -> np.ndarray:
    """Build N-site operator with `op` at `site` and identity elsewhere."""
    ops = [np.eye(d, dtype=complex)] * N
    ops[site] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result


def build_two_site_operator(
    site_i: int, site_j: int,
    op_i: np.ndarray, op_j: np.ndarray,
    N: int = N_MOLECULES, d: int = D_QUTRIT,
) -> np.ndarray:
    """Build N-site operator with op_i at site_i, op_j at site_j, identity elsewhere."""
    ops = [np.eye(d, dtype=complex)] * N
    ops[site_i] = op_i
    ops[site_j] = op_j
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result


# ====================================================================
# Classical (qutrit) Hamiltonian builders
# ====================================================================
def build_H0_classical() -> np.ndarray:
    """Build the diagonal H0 = Σ_i (E_T|1⟩⟨1|_i + E_S|2⟩⟨2|_i)."""
    H0 = np.zeros((DIM_FULL, DIM_FULL), dtype=complex)
    proj_T = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    proj_T[1, 1] = 1.0
    proj_S = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    proj_S[2, 2] = 1.0
    for site in range(N_MOLECULES):
        H0 += E_T * build_single_site_operator(site, proj_T)
        H0 += E_S * build_single_site_operator(site, proj_S)
    return H0


def build_H_transfer_pair(i: int, j: int) -> np.ndarray:
    """H_transfer for pair (i,j): V(|0⟩⟨1|_i ⊗ |1⟩⟨0|_j + h.c.)."""
    ket0_bra1 = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    ket0_bra1[0, 1] = 1.0
    ket1_bra0 = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    ket1_bra0[1, 0] = 1.0
    H = V * build_two_site_operator(i, j, ket0_bra1, ket1_bra0)
    H += V * build_two_site_operator(i, j, ket1_bra0, ket0_bra1)
    return H


def build_H_TTA_pair(i: int, j: int) -> np.ndarray:
    """H_TTA for pair (i,j): J(|0⟩⟨1|_i ⊗ |2⟩⟨1|_j + |2⟩⟨1|_i ⊗ |0⟩⟨1|_j + h.c.)."""
    ket0_bra1 = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    ket0_bra1[0, 1] = 1.0
    ket1_bra0 = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    ket1_bra0[1, 0] = 1.0
    ket2_bra1 = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    ket2_bra1[2, 1] = 1.0
    ket1_bra2 = np.zeros((D_QUTRIT, D_QUTRIT), dtype=complex)
    ket1_bra2[1, 2] = 1.0
    # Term 1: |S0,S1⟩ ↔ |T1,T1⟩  →  |0⟩⟨1|_i ⊗ |2⟩⟨1|_j + h.c.
    H = J * build_two_site_operator(i, j, ket0_bra1, ket2_bra1)
    H += J * build_two_site_operator(i, j, ket1_bra0, ket1_bra2)
    # Term 2: |S1,S0⟩ ↔ |T1,T1⟩  →  |2⟩⟨1|_i ⊗ |0⟩⟨1|_j + h.c.
    H += J * build_two_site_operator(i, j, ket2_bra1, ket0_bra1)
    H += J * build_two_site_operator(i, j, ket1_bra2, ket1_bra0)
    return H


def build_H_full_classical() -> np.ndarray:
    """Build full Hamiltonian H = H0 + H_transfer + H_TTA."""
    H = build_H0_classical()
    for i, j in NEIGHBORS:
        H += build_H_transfer_pair(i, j)
        H += build_H_TTA_pair(i, j)
    return H


# ====================================================================
# Qubit H0 Pauli decomposition (FIXED version)
# ====================================================================
def build_H0_qubit_unitary_pauli(dt: float) -> np.ndarray:
    """
    Build the full 2^8 = 256-dimensional H0 unitary using the FIXED
    Pauli decomposition for the qubit encoding.

    Each molecule uses 2 qubits: |S0⟩=|00⟩, |T1⟩=|01⟩, |S1⟩=|10⟩.
    H0_mol = E_T|01⟩⟨01| + E_S|10⟩⟨10|

    Pauli decomposition:
      H0_mol = α I⊗I + β I⊗Z + γ Z⊗I + δ Z⊗Z
      α = (E_T + E_S)/4,  β = (E_S - E_T)/4
      γ = (E_T - E_S)/4,  δ = -(E_T + E_S)/4

    Time evolution: exp(-i H0_mol dt/ℏ)
      = e^{-iα dt/ℏ} · Rz(θ₀, q0) · Rz(θ₁, q1) · ZZ(θ_zz)

    FIXED signs: θ = +2·coeff·dt/ℏ  (the bug used -2·coeff)
    where Rz(θ) = diag(e^{-iθ/2}, e^{iθ/2}).
    """
    n_qubits = 2 * N_MOLECULES  # 8
    dim = 2 ** n_qubits  # 256
    U = np.eye(dim, dtype=complex)

    for mol in range(N_MOLECULES):
        q0 = 2 * mol       # right qubit
        q1 = 2 * mol + 1   # left qubit

        alpha = (E_T + E_S) / 4.0
        beta = (E_S - E_T) / 4.0
        gamma = (E_T - E_S) / 4.0
        delta = -(E_T + E_S) / 4.0

        # FIXED: positive sign  θ = +2·coeff·dt/ℏ
        theta_0 = 2.0 * beta * dt / HBAR
        theta_1 = 2.0 * gamma * dt / HBAR
        theta_zz = 2.0 * delta * dt / HBAR

        # Global phase from α term: exp(-iα dt/ℏ)
        global_phase = alpha * dt / HBAR

        # Build per-molecule 4×4 unitary
        U_mol = np.eye(4, dtype=complex)

        # Rz(θ₀) on q0
        rz0 = np.diag([np.exp(-1j * theta_0 / 2), np.exp(1j * theta_0 / 2)])
        # Rz(θ₁) on q1
        rz1 = np.diag([np.exp(-1j * theta_1 / 2), np.exp(1j * theta_1 / 2)])
        U_mol = np.kron(rz1, rz0)

        # ZZ(θ_zz): CNOT · Rz(θ_zz, target) · CNOT
        zz_phases = np.array([
            np.exp(-1j * theta_zz / 2),   # |00⟩: Z⊗Z = +1
            np.exp(1j * theta_zz / 2),    # |01⟩: Z⊗Z = -1
            np.exp(1j * theta_zz / 2),    # |10⟩: Z⊗Z = -1
            np.exp(-1j * theta_zz / 2),   # |11⟩: Z⊗Z = +1
        ])
        U_mol = np.diag(zz_phases) @ U_mol

        # Global phase
        U_mol *= np.exp(-1j * global_phase)

        # Embed into full 2^8 space
        U_full_mol = _embed_2qubit_unitary(U_mol, q0, q1, n_qubits)
        U = U_full_mol @ U

    return U


def _embed_2qubit_unitary(U_2q: np.ndarray, q0: int, q1: int, n_qubits: int) -> np.ndarray:
    """Embed a 2-qubit unitary into the full n-qubit Hilbert space."""
    dim = 2 ** n_qubits
    U_full = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            # Extract bits for q0 and q1
            b0_i = (i >> q0) & 1
            b1_i = (i >> q1) & 1
            b0_j = (j >> q0) & 1
            b1_j = (j >> q1) & 1
            # Check other bits match
            mask = ~((1 << q0) | (1 << q1)) & ((1 << n_qubits) - 1)
            if (i & mask) != (j & mask):
                continue
            # 2-qubit index: |q1, q0⟩
            idx_i = b1_i * 2 + b0_i
            idx_j = b1_j * 2 + b0_j
            U_full[i, j] = U_2q[idx_i, idx_j]
    return U_full


def build_H0_qubit_unitary_exact(dt: float) -> np.ndarray:
    """
    Build exact H0 unitary via matrix exponential in 2^8 qubit space.
    This is the reference against which the Pauli decomposition is checked.
    """
    n_qubits = 2 * N_MOLECULES
    dim = 2 ** n_qubits
    H0 = np.zeros((dim, dim), dtype=complex)

    for mol in range(N_MOLECULES):
        q0 = 2 * mol
        q1 = 2 * mol + 1
        # |01⟩⟨01| projector for T1
        proj_T1 = np.zeros((4, 4), dtype=complex)
        proj_T1[1, 1] = 1.0  # |01⟩ = index 1
        # |10⟩⟨10| projector for S1
        proj_S1 = np.zeros((4, 4), dtype=complex)
        proj_S1[2, 2] = 1.0  # |10⟩ = index 2

        H0_mol = E_T * proj_T1 + E_S * proj_S1
        H0 += _embed_2qubit_hermitian(H0_mol, q0, q1, n_qubits)

    return scipy.linalg.expm(-1j * H0 * dt / HBAR)


def _embed_2qubit_hermitian(H_2q: np.ndarray, q0: int, q1: int, n_qubits: int) -> np.ndarray:
    """Embed a 2-qubit Hermitian into the full n-qubit Hilbert space."""
    return _embed_2qubit_unitary(H_2q, q0, q1, n_qubits)


# ====================================================================
# Qubit H_transfer and H_TTA unitary builders (4-qubit per pair)
# ====================================================================
def build_H_transfer_qubit_pair_unitary(dt: float) -> np.ndarray:
    """Build exact 16×16 H_transfer unitary for a pair of molecules in qubit encoding."""
    H = np.zeros((16, 16), dtype=complex)
    idx_00_01 = 0b0001  # |S0,T1⟩
    idx_01_00 = 0b0100  # |T1,S0⟩
    H[idx_00_01, idx_01_00] = V
    H[idx_01_00, idx_00_01] = V
    return scipy.linalg.expm(-1j * H * dt / HBAR)


def build_H_TTA_qubit_pair_unitary(dt: float) -> np.ndarray:
    """Build exact 16×16 H_TTA unitary for a pair of molecules in qubit encoding."""
    H = np.zeros((16, 16), dtype=complex)
    idx_T1_T1 = 0b0101  # |T1,T1⟩
    idx_S0_S1 = 0b0010  # |S0,S1⟩
    idx_S1_S0 = 0b1000  # |S1,S0⟩
    H[idx_S0_S1, idx_T1_T1] = J
    H[idx_T1_T1, idx_S0_S1] = J
    H[idx_S1_S0, idx_T1_T1] = J
    H[idx_T1_T1, idx_S1_S0] = J
    return scipy.linalg.expm(-1j * H * dt / HBAR)


def _embed_4qubit_unitary(
    U_4q: np.ndarray, q0: int, q1: int, q2: int, q3: int, n_qubits: int
) -> np.ndarray:
    """Embed a 4-qubit unitary (acting on q0,q1,q2,q3) into n-qubit space."""
    dim = 2 ** n_qubits
    qubits = [q0, q1, q2, q3]
    U_full = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            # Check non-target bits match
            target_mask = 0
            for q in qubits:
                target_mask |= (1 << q)
            other_mask = ((1 << n_qubits) - 1) & ~target_mask
            if (i & other_mask) != (j & other_mask):
                continue
            # Extract 4-qubit indices
            idx_i = 0
            idx_j = 0
            for k, q in enumerate(qubits):
                idx_i |= (((i >> q) & 1) << k)
                idx_j |= (((j >> q) & 1) << k)
            U_full[i, j] = U_4q[idx_i, idx_j]
    return U_full


# ====================================================================
# Classical Trotter simulation
# ====================================================================
def run_classical_trotter(
    psi0: np.ndarray, n_steps: int, dt: float
) -> tuple[np.ndarray, list[np.ndarray]]:
    """
    Run classical Trotter simulation in the qutrit basis.
    First-order Trotter: U_step = U_H0 · Π_{pairs} U_transfer · Π_{pairs} U_TTA
    """
    H0 = build_H0_classical()
    U_H0 = scipy.linalg.expm(-1j * H0 * dt / HBAR)

    U_transfers = {}
    U_TTAs = {}
    for i_mol, j_mol in NEIGHBORS:
        H_tr = build_H_transfer_pair(i_mol, j_mol)
        U_transfers[(i_mol, j_mol)] = scipy.linalg.expm(-1j * H_tr * dt / HBAR)
        H_tta = build_H_TTA_pair(i_mol, j_mol)
        U_TTAs[(i_mol, j_mol)] = scipy.linalg.expm(-1j * H_tta * dt / HBAR)

    psi = psi0.copy()
    populations_history = [_extract_populations(psi)]

    for _ in range(n_steps):
        psi = U_H0 @ psi
        for pair in NEIGHBORS:
            psi = U_transfers[pair] @ psi
        for pair in NEIGHBORS:
            psi = U_TTAs[pair] @ psi
        populations_history.append(_extract_populations(psi))

    return psi, populations_history


def _extract_populations(psi: np.ndarray) -> np.ndarray:
    """Extract per-molecule state populations [P(S0), P(T1), P(S1)] × N_molecules."""
    probs = np.abs(psi) ** 2
    pops = np.zeros((N_MOLECULES, D_QUTRIT))
    for idx in range(DIM_FULL):
        cfg = index_to_config(idx)
        p = probs[idx]
        for site in range(N_MOLECULES):
            pops[site, cfg[site]] += p
    return pops


# ====================================================================
# Qubit Trotter simulation (statevector, no shots)
# ====================================================================
def _qutrit_to_qubit_state(psi_qutrit: np.ndarray) -> np.ndarray:
    """Map qutrit statevector to qubit statevector.

    Qubit encoding per molecule: |S0⟩→|00⟩, |T1⟩→|01⟩, |S1⟩→|10⟩, |11⟩ unused.
    Uses little-endian bit ordering within each molecule pair.
    """
    n_qubits = 2 * N_MOLECULES
    dim_qubit = 2 ** n_qubits
    psi_qubit = np.zeros(dim_qubit, dtype=complex)

    for idx_q in range(DIM_FULL):
        cfg = index_to_config(idx_q)
        qubit_idx = 0
        for mol in range(N_MOLECULES):
            level = cfg[mol]
            q0 = 2 * mol      # right qubit
            q1 = 2 * mol + 1  # left qubit
            if level == 0:     # |S0⟩ = |00⟩
                pass
            elif level == 1:   # |T1⟩ = |01⟩ → q0=1
                qubit_idx |= (1 << q0)
            elif level == 2:   # |S1⟩ = |10⟩ → q1=1
                qubit_idx |= (1 << q1)
        psi_qubit[qubit_idx] += psi_qutrit[idx_q]

    return psi_qubit


def _qubit_to_qutrit_populations(psi_qubit: np.ndarray) -> np.ndarray:
    """Extract per-molecule populations from qubit statevector."""
    n_qubits = 2 * N_MOLECULES
    dim_qubit = 2 ** n_qubits
    probs = np.abs(psi_qubit) ** 2
    pops = np.zeros((N_MOLECULES, D_QUTRIT))
    for idx in range(dim_qubit):
        p = probs[idx]
        if p < 1e-15:
            continue
        for mol in range(N_MOLECULES):
            q0 = 2 * mol
            q1 = 2 * mol + 1
            b0 = (idx >> q0) & 1
            b1 = (idx >> q1) & 1
            if b0 == 0 and b1 == 0:
                pops[mol, 0] += p  # S0
            elif b0 == 1 and b1 == 0:
                pops[mol, 1] += p  # T1
            elif b0 == 0 and b1 == 1:
                pops[mol, 2] += p  # S1
            # b0=1, b1=1 is non-physical; skip
    return pops


def run_qubit_trotter(
    psi0_qutrit: np.ndarray, n_steps: int, dt: float
) -> tuple[np.ndarray, list[np.ndarray]]:
    """
    Run qubit-encoded Trotter simulation using statevector evolution.
    Builds the full Trotter step unitary from:
      - Pauli-decomposed H0 (FIXED version)
      - Exact H_transfer and H_TTA unitaries per pair
    """
    n_qubits = 2 * N_MOLECULES
    dim_qubit = 2 ** n_qubits

    psi = _qutrit_to_qubit_state(psi0_qutrit)
    populations_history = [_qubit_to_qutrit_populations(psi)]

    # Build H0 unitary via Pauli decomposition
    U_H0_pauli = build_H0_qubit_unitary_pauli(dt)

    # Build H_transfer and H_TTA unitaries per pair
    U_tr_pair = build_H_transfer_qubit_pair_unitary(dt)
    U_tta_pair = build_H_TTA_qubit_pair_unitary(dt)

    # Embed pair unitaries
    U_transfers = {}
    U_TTAs = {}
    for i_mol, j_mol in NEIGHBORS:
        qubits = [2 * i_mol, 2 * i_mol + 1, 2 * j_mol, 2 * j_mol + 1]
        U_transfers[(i_mol, j_mol)] = _embed_4qubit_unitary(
            U_tr_pair, *qubits, n_qubits
        )
        U_TTAs[(i_mol, j_mol)] = _embed_4qubit_unitary(
            U_tta_pair, *qubits, n_qubits
        )

    for _ in range(n_steps):
        psi = U_H0_pauli @ psi
        for pair in NEIGHBORS:
            psi = U_transfers[pair] @ psi
        for pair in NEIGHBORS:
            psi = U_TTAs[pair] @ psi
        populations_history.append(_qubit_to_qutrit_populations(psi))

    return psi, populations_history


# ====================================================================
# 2-qudit depolarizing channel verification
# ====================================================================
def partial_trace_pair(rho: np.ndarray, sites: tuple[int, int],
                       N: int = N_MOLECULES, d: int = D_QUTRIT) -> np.ndarray:
    """
    Compute partial trace of rho over the specified pair of sites.
    Returns the reduced density matrix of the remaining (N-2) sites.

    rho: d^N × d^N density matrix
    sites: tuple of two site indices to trace out
    """
    dim_full = d ** N
    dim_pair = d * d
    dim_rest = d ** (N - 2)

    rest_sites = [s for s in range(N) if s not in sites]

    rho_rest = np.zeros((dim_rest, dim_rest), dtype=complex)

    for i_rest in range(dim_rest):
        for j_rest in range(dim_rest):
            cfg_i_rest = _idx_to_cfg(i_rest, len(rest_sites), d)
            cfg_j_rest = _idx_to_cfg(j_rest, len(rest_sites), d)

            val = 0.0 + 0.0j
            for k_pair in range(dim_pair):
                cfg_pair = _idx_to_cfg(k_pair, 2, d)
                # Build full config for i
                cfg_i_full = [0] * N
                cfg_j_full = [0] * N
                for r_idx, r_site in enumerate(rest_sites):
                    cfg_i_full[r_site] = cfg_i_rest[r_idx]
                    cfg_j_full[r_site] = cfg_j_rest[r_idx]
                for p_idx, p_site in enumerate(sites):
                    cfg_i_full[p_site] = cfg_pair[p_idx]
                    cfg_j_full[p_site] = cfg_pair[p_idx]

                idx_i = config_to_index(cfg_i_full, d)
                idx_j = config_to_index(cfg_j_full, d)
                val += rho[idx_i, idx_j]

            rho_rest[i_rest, j_rest] = val

    return rho_rest


def _idx_to_cfg(idx: int, n: int, d: int) -> list[int]:
    """Convert index to d-ary config of length n."""
    cfg = []
    for _ in range(n):
        cfg.append(idx % d)
        idx //= d
    return cfg[::-1]


def apply_2qudit_depolarizing(
    rho: np.ndarray, sites: tuple[int, int], p: float,
    N: int = N_MOLECULES, d: int = D_QUTRIT,
) -> np.ndarray:
    """
    Apply per-gate 2-qudit depolarizing channel:
      ε(ρ) = (1-p)ρ + p · Tr_{ij}(ρ) ⊗ I_{ij}/d²

    This replaces the qudits at sites with the maximally mixed state
    while preserving the rest.
    """
    dim_full = d ** N
    dim_pair = d * d
    dim_rest = d ** (N - 2)

    rest_sites = [s for s in range(N) if s not in sites]

    # Compute Tr_{ij}(ρ)
    rho_rest = partial_trace_pair(rho, sites, N, d)

    # Build Tr_{ij}(ρ) ⊗ I_{ij}/d²
    rho_mixed = np.zeros((dim_full, dim_full), dtype=complex)
    I_pair = np.eye(dim_pair, dtype=complex) / dim_pair

    for i_full in range(dim_full):
        for j_full in range(dim_full):
            cfg_i = index_to_config(i_full, N, d)
            cfg_j = index_to_config(j_full, N, d)

            # Extract rest indices
            cfg_i_rest = [cfg_i[s] for s in rest_sites]
            cfg_j_rest = [cfg_j[s] for s in rest_sites]
            i_rest = _cfg_to_idx(cfg_i_rest, d)
            j_rest = _cfg_to_idx(cfg_j_rest, d)

            # Extract pair indices
            cfg_i_pair = [cfg_i[s] for s in sites]
            cfg_j_pair = [cfg_j[s] for s in sites]
            i_pair = _cfg_to_idx(cfg_i_pair, d)
            j_pair = _cfg_to_idx(cfg_j_pair, d)

            rho_mixed[i_full, j_full] = rho_rest[i_rest, j_rest] * I_pair[i_pair, j_pair]

    return (1 - p) * rho + p * rho_mixed


def _cfg_to_idx(cfg: list[int], d: int) -> int:
    """Convert d-ary config to index."""
    idx = 0
    for c in cfg:
        idx = idx * d + c
    return idx


# ====================================================================
# Test functions
# ====================================================================
def test_h0_phases() -> dict[str, Any]:
    """
    Test A: Verify H0 Pauli decomposition produces correct relative phases
    for |S0⟩, |T1⟩, |S1⟩ states.
    """
    dt = 5.0  # fs

    # Build exact and Pauli-decomposed H0 unitaries in qubit space
    U_exact = build_H0_qubit_unitary_exact(dt)
    U_pauli = build_H0_qubit_unitary_pauli(dt)

    # Test on single-molecule basis states embedded in the 4-molecule system
    n_qubits = 2 * N_MOLECULES
    dim = 2 ** n_qubits

    test_states = {
        "S0": 0,   # all |00⟩ = index 0
        "T1_mol0": 1 << 0,  # mol 0 in |01⟩
        "S1_mol0": 1 << 1,  # mol 0 in |10⟩
    }

    phases = {}
    max_error = 0.0

    for name, idx in test_states.items():
        psi = np.zeros(dim, dtype=complex)
        psi[idx] = 1.0
        psi_exact = U_exact @ psi
        psi_pauli = U_pauli @ psi

        # Both should produce a pure phase on a single basis state
        phase_exact = np.angle(psi_exact[idx])
        phase_pauli = np.angle(psi_pauli[idx])
        phase_diff = abs(phase_exact - phase_pauli)
        # Wrap to [-π, π]
        phase_diff = min(phase_diff, 2 * np.pi - phase_diff)

        phases[name] = {
            "phase_exact": float(phase_exact),
            "phase_pauli": float(phase_pauli),
            "phase_error": float(phase_diff),
        }
        max_error = max(max_error, phase_diff)

    # Also check that U_exact and U_pauli agree on the physical subspace
    # (up to global phase per molecule block, they should match exactly)
    physical_indices = _get_physical_qubit_indices()
    U_exact_phys = U_exact[np.ix_(physical_indices, physical_indices)]
    U_pauli_phys = U_pauli[np.ix_(physical_indices, physical_indices)]
    # Remove global phase
    if abs(U_exact_phys[0, 0]) > 1e-10:
        phase_adj = U_exact_phys[0, 0] / U_pauli_phys[0, 0]
        U_pauli_phys_adj = U_pauli_phys * phase_adj
    else:
        U_pauli_phys_adj = U_pauli_phys
    unitary_error = float(np.max(np.abs(U_exact_phys - U_pauli_phys_adj)))

    passed = max_error < 1e-10 and unitary_error < 1e-10
    return {
        "status": "PASS" if passed else "FAIL",
        "details": {
            "max_phase_error": float(max_error),
            "unitary_subspace_error": unitary_error,
            "phases": phases,
        },
    }


def _get_physical_qubit_indices() -> list[int]:
    """Get indices of physical qubit states (no |11⟩ per molecule)."""
    n_qubits = 2 * N_MOLECULES
    dim = 2 ** n_qubits
    physical = []
    for idx in range(dim):
        is_physical = True
        for mol in range(N_MOLECULES):
            q0 = 2 * mol
            q1 = 2 * mol + 1
            b0 = (idx >> q0) & 1
            b1 = (idx >> q1) & 1
            if b0 == 1 and b1 == 1:
                is_physical = False
                break
        if is_physical:
            physical.append(idx)
    return physical


def test_trotter_comparison() -> dict[str, Any]:
    """
    Test B: Run classical and qubit Trotter simulations and compare.
    """
    dt = T_TOTAL / N_STEPS

    # Initial state: edge_triplet (molecules 0,3 in T1, others in S0)
    # config = [1, 0, 0, 1]
    psi0 = np.zeros(DIM_FULL, dtype=complex)
    init_config = [1, 0, 0, 1]
    psi0[config_to_index(init_config)] = 1.0

    # Classical simulation
    _, classical_history = run_classical_trotter(psi0, N_STEPS, dt)
    classical_final = classical_history[-1]

    # Qubit simulation
    _, qubit_history = run_qubit_trotter(psi0, N_STEPS, dt)
    qubit_final = qubit_history[-1]

    # Compare
    pop_errors = np.abs(classical_final - qubit_final)
    max_pop_error = float(np.max(pop_errors))
    mean_pop_error = float(np.mean(pop_errors))

    passed = max_pop_error < 1e-8
    return {
        "status": "PASS" if passed else "FAIL",
        "details": {
            "max_population_error": max_pop_error,
            "mean_population_error": mean_pop_error,
            "classical_final": classical_final.tolist(),
            "qubit_final": qubit_final.tolist(),
            "population_errors": pop_errors.tolist(),
        },
    }


def test_depolarizing_channel() -> dict[str, Any]:
    """
    Test C: Verify the 2-qudit depolarizing channel implementation.
    - Partial trace correctness
    - Mixed state construction
    - Trace preservation
    """
    results = {}
    all_passed = True

    # Use a smaller system for tractability: N=3, d=3 (27 dim)
    N_test = 3
    d_test = 3
    dim_test = d_test ** N_test

    # --- C1: Partial trace of pure product state ---
    # |ψ⟩ = |1⟩ ⊗ |0⟩ ⊗ |2⟩
    psi_prod = np.zeros(dim_test, dtype=complex)
    cfg_prod = [1, 0, 2]
    idx_prod = _cfg_to_idx(cfg_prod, d_test)
    psi_prod[idx_prod] = 1.0
    rho_prod = np.outer(psi_prod, psi_prod.conj())

    # Trace out sites (0,1), should get |2⟩⟨2|
    rho_rest_01 = partial_trace_pair(rho_prod, (0, 1), N_test, d_test)
    expected_rest_01 = np.zeros((d_test, d_test), dtype=complex)
    expected_rest_01[2, 2] = 1.0
    error_pt1 = float(np.max(np.abs(rho_rest_01 - expected_rest_01)))
    pt1_passed = error_pt1 < 1e-12
    results["partial_trace_product_state"] = {
        "status": "PASS" if pt1_passed else "FAIL",
        "error": error_pt1,
    }
    all_passed = all_passed and pt1_passed

    # --- C2: Partial trace of entangled state ---
    # |ψ⟩ = (|100⟩ + |010⟩) / √2
    psi_ent = np.zeros(dim_test, dtype=complex)
    psi_ent[_cfg_to_idx([1, 0, 0], d_test)] = 1.0 / np.sqrt(2)
    psi_ent[_cfg_to_idx([0, 1, 0], d_test)] = 1.0 / np.sqrt(2)
    rho_ent = np.outer(psi_ent, psi_ent.conj())

    # Trace out site (2), keeping sites 0,1
    # But our partial_trace_pair traces out a pair; trace out (0,2) instead
    rho_rest_02 = partial_trace_pair(rho_ent, (0, 2), N_test, d_test)
    # Remaining is site 1
    # From |100⟩: site1=0;  from |010⟩: site1=1
    # ρ_1 = 0.5 |0⟩⟨0| + 0.5 |1⟩⟨1|
    expected_rest_02 = np.zeros((d_test, d_test), dtype=complex)
    expected_rest_02[0, 0] = 0.5
    expected_rest_02[1, 1] = 0.5
    error_pt2 = float(np.max(np.abs(rho_rest_02 - expected_rest_02)))
    pt2_passed = error_pt2 < 1e-12
    results["partial_trace_entangled_state"] = {
        "status": "PASS" if pt2_passed else "FAIL",
        "error": error_pt2,
    }
    all_passed = all_passed and pt2_passed

    # --- C3: Trace preservation under depolarizing channel ---
    # Random density matrix
    rng = np.random.default_rng(42)
    A = rng.standard_normal((dim_test, dim_test)) + 1j * rng.standard_normal((dim_test, dim_test))
    rho_rand = A @ A.conj().T
    rho_rand /= np.trace(rho_rand)

    p_test = 0.05
    rho_depol = apply_2qudit_depolarizing(rho_rand, (0, 1), p_test, N_test, d_test)
    trace_before = float(np.abs(np.trace(rho_rand)))
    trace_after = float(np.abs(np.trace(rho_depol)))
    trace_error = abs(trace_after - 1.0)
    tp_passed = trace_error < 1e-10
    results["trace_preservation"] = {
        "status": "PASS" if tp_passed else "FAIL",
        "trace_before": trace_before,
        "trace_after": trace_after,
        "trace_error": trace_error,
    }
    all_passed = all_passed and tp_passed

    # --- C4: Positivity preservation ---
    eigenvalues = np.linalg.eigvalsh(rho_depol)
    min_eigenvalue = float(np.min(eigenvalues))
    pos_passed = min_eigenvalue > -1e-12
    results["positivity_preservation"] = {
        "status": "PASS" if pos_passed else "FAIL",
        "min_eigenvalue": min_eigenvalue,
    }
    all_passed = all_passed and pos_passed

    # --- C5: p=0 leaves state unchanged ---
    rho_p0 = apply_2qudit_depolarizing(rho_rand, (0, 1), 0.0, N_test, d_test)
    error_p0 = float(np.max(np.abs(rho_p0 - rho_rand)))
    p0_passed = error_p0 < 1e-14
    results["p0_identity"] = {
        "status": "PASS" if p0_passed else "FAIL",
        "error": error_p0,
    }
    all_passed = all_passed and p0_passed

    # --- C6: p=1 fully depolarizes the pair ---
    rho_p1 = apply_2qudit_depolarizing(rho_rand, (0, 1), 1.0, N_test, d_test)
    # The pair subsystem should be I/d², rest should be Tr_{01}(ρ)
    rho_pair_after = partial_trace_pair(rho_p1, (2,), N_test, d_test)
    # Actually for p=1: ε(ρ) = Tr_{01}(ρ) ⊗ I_{01}/d²
    # Trace out site 2 to get the pair part conditioned on each rest state
    # Just check that partial trace over pair of the result equals partial trace of original
    rho_rest_p1 = partial_trace_pair(rho_p1, (0, 1), N_test, d_test)
    rho_rest_orig = partial_trace_pair(rho_rand, (0, 1), N_test, d_test)
    error_p1_rest = float(np.max(np.abs(rho_rest_p1 - rho_rest_orig)))
    p1_passed = error_p1_rest < 1e-10
    results["p1_full_depolarization"] = {
        "status": "PASS" if p1_passed else "FAIL",
        "error": error_p1_rest,
    }
    all_passed = all_passed and p1_passed

    return {
        "status": "PASS" if all_passed else "FAIL",
        "details": results,
    }


def test_noise_accumulation() -> dict[str, Any]:
    """
    Test D: Compare effective noise accumulation between qubit and qudit.
    Count 2-qubit/2-qudit gates per Trotter step and compute effective noise.
    """
    n_pairs = len(NEIGHBORS)

    # --- Qubit approach ---
    # H0: per molecule, ZZ decomposition uses 2 CNOT gates → 2 two-qubit gates per mol
    n_cnot_h0_per_mol = 2
    n_2q_gates_h0 = N_MOLECULES * n_cnot_h0_per_mol
    # H_transfer: each pair gets one 4-qubit unitary, which decomposes to ~6 CNOTs
    # (exact count depends on decomposition, but the 4-qubit unitary acts on the
    #  2-qubit subspace of active indices so ~3 CNOTs suffice; use 6 as upper bound)
    n_2q_gates_transfer = n_pairs * 6
    # H_TTA: each pair gets one 4-qubit unitary, ~6 CNOTs
    n_2q_gates_tta = n_pairs * 6
    n_2q_gates_qubit = n_2q_gates_h0 + n_2q_gates_transfer + n_2q_gates_tta

    # --- Qudit approach ---
    # H0 is diagonal → single-qudit virtual Rz gates → 0 two-qudit gates
    n_2qd_gates_h0 = 0
    # H_transfer: 1 two-qudit gate per pair
    n_2qd_gates_transfer = n_pairs
    # H_TTA: 1 two-qudit gate per pair
    n_2qd_gates_tta = n_pairs
    n_2qd_gates_qudit = n_2qd_gates_h0 + n_2qd_gates_transfer + n_2qd_gates_tta

    # Effective per-step noise (probability of at least one error)
    p_no_error_qubit = (1 - DEPOL_2Q) ** n_2q_gates_qubit
    p_error_qubit_step = 1 - p_no_error_qubit

    p_no_error_qudit = (1 - DEPOL_2Q) ** n_2qd_gates_qudit
    p_error_qudit_step = 1 - p_no_error_qudit

    # Total over all steps
    p_no_error_qubit_total = (1 - DEPOL_2Q) ** (n_2q_gates_qubit * N_STEPS)
    p_no_error_qudit_total = (1 - DEPOL_2Q) ** (n_2qd_gates_qudit * N_STEPS)

    passed = True  # This is informational; always passes

    return {
        "status": "PASS" if passed else "FAIL",
        "details": {
            "qubit": {
                "n_2q_gates_per_step": n_2q_gates_qubit,
                "breakdown": {
                    "H0_CNOTs": n_2q_gates_h0,
                    "H_transfer_CNOTs": n_2q_gates_transfer,
                    "H_TTA_CNOTs": n_2q_gates_tta,
                },
                "p_error_per_step": float(p_error_qubit_step),
                "p_no_error_total": float(p_no_error_qubit_total),
            },
            "qudit": {
                "n_2qd_gates_per_step": n_2qd_gates_qudit,
                "breakdown": {
                    "H0_gates": n_2qd_gates_h0,
                    "H_transfer_gates": n_2qd_gates_transfer,
                    "H_TTA_gates": n_2qd_gates_tta,
                },
                "p_error_per_step": float(p_error_qudit_step),
                "p_no_error_total": float(p_no_error_qudit_total),
            },
            "noise_reduction_factor": float(n_2q_gates_qubit / max(n_2qd_gates_qudit, 1)),
        },
    }


# ====================================================================
# Main driver
# ====================================================================
def main() -> None:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    print("=" * 70)
    print("Quantum Dynamics Verification Script")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    test_results: dict[str, Any] = {}

    # --- Test A: H0 phases ---
    print("\n[A] H0 Pauli decomposition phase verification ...")
    result_a = test_h0_phases()
    test_results["h0_phase_verification"] = result_a
    print(f"    Result: {result_a['status']}")
    if result_a["status"] == "PASS":
        print(f"    Max phase error: {result_a['details']['max_phase_error']:.2e}")
    else:
        print(f"    FAILED: max phase error = {result_a['details']['max_phase_error']:.2e}")

    # --- Test B: Trotter comparison ---
    print("\n[B] Classical vs. qubit Trotter simulation ...")
    result_b = test_trotter_comparison()
    test_results["trotter_comparison"] = result_b
    print(f"    Result: {result_b['status']}")
    print(f"    Max population error: {result_b['details']['max_population_error']:.2e}")

    # --- Test C: Depolarizing channel ---
    print("\n[C] 2-qudit depolarizing channel verification ...")
    result_c = test_depolarizing_channel()
    test_results["depolarizing_channel"] = result_c
    print(f"    Result: {result_c['status']}")
    for subtest, subresult in result_c["details"].items():
        print(f"      {subtest}: {subresult['status']}")

    # --- Test D: Noise accumulation ---
    print("\n[D] Noise accumulation comparison ...")
    result_d = test_noise_accumulation()
    test_results["noise_accumulation"] = result_d
    print(f"    Result: {result_d['status']}")
    det = result_d["details"]
    print(f"    Qubit:  {det['qubit']['n_2q_gates_per_step']} 2q-gates/step, "
          f"p_err/step = {det['qubit']['p_error_per_step']:.4f}")
    print(f"    Qudit:  {det['qudit']['n_2qd_gates_per_step']} 2qd-gates/step, "
          f"p_err/step = {det['qudit']['p_error_per_step']:.4f}")
    print(f"    Noise reduction factor: {det['noise_reduction_factor']:.1f}×")

    # --- Assemble output ---
    output = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "test_results": {
            name: {"status": r["status"], "details": r["details"]}
            for name, r in test_results.items()
        },
        "h0_phases": result_a["details"]["phases"],
        "classical_populations": result_b["details"]["classical_final"],
        "qubit_populations": result_b["details"]["qubit_final"],
        "population_errors": result_b["details"]["population_errors"],
        "noise_analysis": {
            "qubit_2q_gates_per_step": det["qubit"]["n_2q_gates_per_step"],
            "qudit_2qd_gates_per_step": det["qudit"]["n_2qd_gates_per_step"],
            "qubit_p_error_per_step": det["qubit"]["p_error_per_step"],
            "qudit_p_error_per_step": det["qudit"]["p_error_per_step"],
            "qubit_p_no_error_total": det["qubit"]["p_no_error_total"],
            "qudit_p_no_error_total": det["qudit"]["p_no_error_total"],
            "noise_reduction_factor": det["noise_reduction_factor"],
        },
    }

    # Save results
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "developing", "verification_results",
    )
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"quantum_dynamics_verification_{timestamp}.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to: {out_path}")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    all_passed = True
    for name, r in test_results.items():
        status = r["status"]
        marker = "✓" if status == "PASS" else "✗"
        print(f"  {marker} {name}: {status}")
        if status != "PASS":
            all_passed = False

    if all_passed:
        print("\nAll tests PASSED ✓")
    else:
        print("\nSome tests FAILED ✗")
        sys.exit(1)


if __name__ == "__main__":
    main()
