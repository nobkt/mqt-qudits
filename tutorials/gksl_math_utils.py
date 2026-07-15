"""Mathematical utility functions for GKSL-Lindblad quantum dynamics."""

from __future__ import annotations

import os
import sys
from functools import reduce

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters


def build_single_site_operator(
    op: np.ndarray, site: int, N: int, d: int
) -> np.ndarray:
    """Build full-system operator from a single-site operator using tensor products."""
    eye = np.eye(d, dtype=np.complex128)
    op_list = [eye] * N
    op_list[site] = np.asarray(op, dtype=np.complex128)
    return reduce(np.kron, op_list)


def build_onsite_hamiltonian(params: GKSLPhysicalParameters) -> np.ndarray:
    """H_0 = sum_i (E_T |1><1| + E_S |2><2|)_i."""
    N = params.N_molecules
    d = params.d
    dim = d**N
    h_local = np.diag(np.array([0.0, params.E_T, params.E_S], dtype=np.complex128))
    H0 = np.zeros((dim, dim), dtype=np.complex128)
    for i in range(N):
        H0 += build_single_site_operator(h_local, i, N, d)
    if not np.allclose(H0, H0.conj().T):
        raise ValueError("Onsite Hamiltonian is not Hermitian")
    return H0


def build_transfer_hamiltonian(params: GKSLPhysicalParameters) -> np.ndarray:
    """H_transfer = sum_{<i,j>} V * (|0>_i<1| x |1>_j<0| + h.c.)."""
    N = params.N_molecules
    d = params.d
    dim = d**N
    H_t = np.zeros((dim, dim), dtype=np.complex128)
    eye = np.eye(d, dtype=np.complex128)

    # |0><1| and |1><0| single-site operators
    ket0_bra1 = np.zeros((d, d), dtype=np.complex128)
    ket0_bra1[0, 1] = 1.0
    ket1_bra0 = np.zeros((d, d), dtype=np.complex128)
    ket1_bra0[1, 0] = 1.0

    for i, j in params.neighbors:
        # |0>_i<1| x |1>_j<0|
        op_list_fwd = [eye] * N
        op_list_fwd[i] = ket0_bra1
        op_list_fwd[j] = ket1_bra0
        fwd = reduce(np.kron, op_list_fwd)
        H_t += params.V * (fwd + fwd.conj().T)

    if not np.allclose(H_t, H_t.conj().T):
        raise ValueError("Transfer Hamiltonian is not Hermitian")
    return H_t


def build_lindblad_operators(
    params: GKSLPhysicalParameters,
) -> list[tuple[np.ndarray, float]]:
    """Build all Lindblad operators. Returns list of (L_alpha, gamma_alpha).

    Each L_alpha already includes the sqrt(gamma) factor.
    """
    N = params.N_molecules
    d = params.d
    eye = np.eye(d, dtype=np.complex128)
    ops: list[tuple[np.ndarray, float]] = []

    def _ket_bra(a: int, b: int) -> np.ndarray:
        m = np.zeros((d, d), dtype=np.complex128)
        m[a, b] = 1.0
        return m

    # TTA: 3 pairs x 2 channels = 6
    for i, j in params.neighbors:
        gamma = params.gamma_TTA / 2.0
        sq = np.sqrt(gamma)
        # Channel 1: |2>_i<1| x |0>_j<1|
        op_list = [eye] * N
        op_list[i] = _ket_bra(2, 1)
        op_list[j] = _ket_bra(0, 1)
        ops.append((sq * reduce(np.kron, op_list), gamma))
        # Channel 2: |0>_i<1| x |2>_j<1|
        op_list = [eye] * N
        op_list[i] = _ket_bra(0, 1)
        op_list[j] = _ket_bra(2, 1)
        ops.append((sq * reduce(np.kron, op_list), gamma))

    # Fluorescence: |0>_i<2|, gamma = Gamma_fl
    for i in range(N):
        gamma = params.Gamma_fl
        L = np.sqrt(gamma) * build_single_site_operator(_ket_bra(0, 2), i, N, d)
        ops.append((L, gamma))

    # Phosphorescence: |0>_i<1|, gamma = Gamma_ph
    for i in range(N):
        gamma = params.Gamma_ph
        L = np.sqrt(gamma) * build_single_site_operator(_ket_bra(0, 1), i, N, d)
        ops.append((L, gamma))

    # Internal conversion: |0>_i<2|, gamma = k_IC
    for i in range(N):
        gamma = params.k_IC
        L = np.sqrt(gamma) * build_single_site_operator(_ket_bra(0, 2), i, N, d)
        ops.append((L, gamma))

    # ISC S->T: |1>_i<2|, gamma = k_ISC_ST
    for i in range(N):
        gamma = params.k_ISC_ST
        L = np.sqrt(gamma) * build_single_site_operator(_ket_bra(1, 2), i, N, d)
        ops.append((L, gamma))

    # ISC T->S: |0>_i<1|, gamma = k_ISC_TS
    for i in range(N):
        gamma = params.k_ISC_TS
        L = np.sqrt(gamma) * build_single_site_operator(_ket_bra(0, 1), i, N, d)
        ops.append((L, gamma))

    expected = 2 * len(params.neighbors) + 5 * N
    assert len(ops) == expected, f"Expected {expected} Lindblad operators, got {len(ops)}"
    return ops


def vectorize_density_matrix(rho: np.ndarray) -> np.ndarray:
    """Column-major vectorization of a density matrix."""
    return rho.flatten(order="F")


def unvectorize_density_matrix(vec: np.ndarray, dim: int) -> np.ndarray:
    """Reshape vector back to density matrix (column-major)."""
    return vec.reshape((dim, dim), order="F")


def compute_von_neumann_entropy(rho: np.ndarray, tol: float = 1e-12) -> float:
    """S = -Tr[rho ln rho] = -sum lambda_i ln(lambda_i) for lambda_i > tol."""
    eigenvalues = np.linalg.eigvalsh(rho)
    pos = eigenvalues[eigenvalues > tol]
    return -float(np.sum(pos * np.log(pos)))


def compute_purity(rho: np.ndarray) -> float:
    """P = Tr[rho^2]."""
    return float(np.real(np.trace(rho @ rho)))


def compute_populations_from_density_matrix(
    rho: np.ndarray, params: GKSLPhysicalParameters
) -> dict:
    """Compute N_S0, N_T1, N_S1 populations from density matrix diagonal.

    For each computational basis state, decompose into per-molecule local states
    and accumulate diagonal probability into the corresponding population counter.
    Also returns per_molecule_populations[mol_idx][state].
    """
    N = params.N_molecules
    d = params.d
    dim = d**N
    diag = np.real(np.diag(rho))

    N_S0 = 0.0
    N_T1 = 0.0
    N_S1 = 0.0
    per_molecule: list[dict[int, float]] = [{0: 0.0, 1: 0.0, 2: 0.0} for _ in range(N)]

    for idx in range(dim):
        p = diag[idx]
        remainder = idx
        for mol in range(N - 1, -1, -1):
            local_state = remainder % d
            remainder //= d
            if local_state == 0:
                N_S0 += p
            elif local_state == 1:
                N_T1 += p
            else:
                N_S1 += p
            per_molecule[mol][local_state] += p

    return {
        "N_S0": N_S0,
        "N_T1": N_T1,
        "N_S1": N_S1,
        "per_molecule_populations": per_molecule,
    }


def build_phonon_operators(
    n_max: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (a, a_dag, n_op) for phonon Fock space with dimension n_max+1."""
    dim = n_max + 1
    a = np.zeros((dim, dim), dtype=np.complex128)
    for n in range(1, dim):
        a[n - 1, n] = np.sqrt(n)
    a_dag = a.conj().T
    n_op = a_dag @ a
    return a, a_dag, n_op


def build_H_phonon(params: GKSLPhysicalParameters) -> np.ndarray:
    """H_phonon = I_el x sum_i omega_ph * n_hat_i (in full extended space)."""
    N = params.N_molecules
    d = params.d
    n_max = params.n_max
    dim_el = d**N
    dim_ph_single = n_max + 1
    dim_ph = dim_ph_single**N

    _, _, n_op = build_phonon_operators(n_max)
    eye_ph = np.eye(dim_ph_single, dtype=np.complex128)

    H_ph = np.zeros((dim_ph, dim_ph), dtype=np.complex128)
    for i in range(N):
        op_list = [eye_ph] * N
        op_list[i] = n_op
        H_ph += params.omega_ph * reduce(np.kron, op_list)

    return np.kron(np.eye(dim_el, dtype=np.complex128), H_ph)


def build_H_eph(params: GKSLPhysicalParameters) -> np.ndarray:
    """Holstein coupling: g_eph * sum_i |1>_i<1| x (a_i + a_dag_i)."""
    N = params.N_molecules
    d = params.d
    n_max = params.n_max
    dim_ph_single = n_max + 1

    a, a_dag, _ = build_phonon_operators(n_max)
    x_op = a + a_dag  # (a + a†)
    eye_el = np.eye(d, dtype=np.complex128)
    eye_ph = np.eye(dim_ph_single, dtype=np.complex128)

    # |1><1| for single site
    proj_T = np.zeros((d, d), dtype=np.complex128)
    proj_T[1, 1] = 1.0

    dim_el = d**N
    dim_ph = dim_ph_single**N
    H_eph = np.zeros((dim_el * dim_ph, dim_el * dim_ph), dtype=np.complex128)

    for i in range(N):
        el_list = [eye_el] * N
        el_list[i] = proj_T
        el_op = reduce(np.kron, el_list)

        ph_list = [eye_ph] * N
        ph_list[i] = x_op
        ph_op = reduce(np.kron, ph_list)

        H_eph += params.g_eph * np.kron(el_op, ph_op)

    return H_eph


def build_H_total_boson(params: GKSLPhysicalParameters) -> np.ndarray:
    """Full Hamiltonian: H_el_ext + H_phonon + H_eph."""
    N = params.N_molecules
    d = params.d
    n_max = params.n_max
    dim_el = d**N
    dim_ph_single = n_max + 1
    dim_ph = dim_ph_single**N

    H_el = build_onsite_hamiltonian(params) + build_transfer_hamiltonian(params)
    H_el_ext = np.kron(H_el, np.eye(dim_ph, dtype=np.complex128))

    return H_el_ext + build_H_phonon(params) + build_H_eph(params)


def extend_lindblad_operators(
    lindblad_ops_el: list[tuple[np.ndarray, float]], dim_phonon: int
) -> list[tuple[np.ndarray, float]]:
    """L_ext = L_el x I_phonon for each operator."""
    eye_ph = np.eye(dim_phonon, dtype=np.complex128)
    return [(np.kron(L, eye_ph), gamma) for L, gamma in lindblad_ops_el]


def partial_trace_phonon(
    rho_total: np.ndarray, dim_el: int, dim_ph: int
) -> np.ndarray:
    """Trace out phonon degrees of freedom."""
    rho_reshaped = rho_total.reshape((dim_el, dim_ph, dim_el, dim_ph))
    return np.trace(rho_reshaped, axis1=1, axis2=3).astype(np.complex128)
