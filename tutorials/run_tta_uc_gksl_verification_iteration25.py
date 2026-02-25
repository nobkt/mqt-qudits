#!/usr/bin/env python3
"""TTA-UC GKSL iteration 25 verification: Frobenius analysis and composite Rate(T) prediction.

Changes from iteration 24:
  - Test A-D: unchanged
  - Test E (extended analysis):
    * Added frobenius_distance() for L2-norm distance computation
    * Added Frobenius distance computation for key comparison pairs
    * Part 6f (new): Composite model fitting T(ST,ex)/dt = p0 + p1*dt
    * Part 6g (new): Frobenius-based error vector angle estimation (exact law of cosines)
    * Part 6h (new): Analytical Rate(T) prediction and comparison with observed values
    * Parts 6a-6e retained from iteration 24

  Key issues addressed from iteration 24:
    1. ST vs Exact Rate(T) = 0.87 at dt=0.2->0.1 (predicted 0.95-1.05):
       -> Not a bug; mathematical consequence of error cancellation.
       -> Part 6h adds analytical prediction to verify.
    2. Trace distance-based angle estimation unstable (119-131 deg):
       -> Part 6g adds Frobenius-based angle (exact L2 law of cosines).
    3. No way to predict composite Rate(T):
       -> Part 6f adds composite model fitting.

Results are written to developing/verification_results/ as JSON and Markdown.
"""

from __future__ import annotations

import json
import os
import sys
import time as time_module
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_populations_from_density_matrix,
    compute_purity,
    unvectorize_density_matrix,
    vectorize_density_matrix,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator
from qudit_gksl_simulator import QuditGKSLSimulator


def quantum_fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Compute quantum state fidelity F(rho, sigma) via eigendecomposition."""
    evals, evecs = np.linalg.eigh(rho)
    evals = np.maximum(evals, 0.0)
    sqrt_rho = evecs @ np.diag(np.sqrt(evals)) @ evecs.conj().T
    m = sqrt_rho @ sigma @ sqrt_rho
    m = (m + m.conj().T) / 2
    evals_m = np.linalg.eigvalsh(m)
    evals_m = np.maximum(evals_m, 0.0)
    return float(min(np.real(np.sum(np.sqrt(evals_m))) ** 2, 1.0))


def trace_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Compute trace distance T(rho, sigma) = 0.5 * ||rho - sigma||_1.

    For Hermitian delta = rho - sigma, ||delta||_1 = sum |eigenvalues|.
    """
    delta = rho - sigma
    delta = (delta + delta.conj().T) / 2  # enforce Hermiticity
    eigvals = np.linalg.eigvalsh(delta)
    return 0.5 * float(np.sum(np.abs(eigvals)))


def bures_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Compute Bures distance d_B(rho, sigma) = sqrt(2(1 - sqrt(F(rho, sigma))))."""
    f = quantum_fidelity(rho, sigma)
    val = 2.0 * (1.0 - np.sqrt(f))
    return float(np.sqrt(max(val, 0.0)))


def frobenius_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Compute Frobenius (Hilbert-Schmidt) distance ||rho - sigma||_F.

    Unlike trace distance (L1 norm), the Frobenius norm is an inner-product
    norm (L2), so the law of cosines applies exactly for angle estimation.
    """
    delta = rho - sigma
    return float(np.linalg.norm(delta, "fro"))


def normalize_rho(rho: np.ndarray) -> np.ndarray:
    """Normalize density matrix to trace 1."""
    tr = float(np.real(np.trace(rho)))
    return rho / tr if tr > 1e-10 else rho


# =========================================================================
# Test A: Noise vs Trotter separation (unchanged from iteration 17)
# =========================================================================
def test_noise_vs_trotter_separation(params: GKSLPhysicalParameters) -> dict:
    """Separately measure noise infidelity and Trotter infidelity."""
    print("\n" + "=" * 70)
    print("Test A: Noise vs Trotter Separation (Qudit, pair-only)")
    print("=" * 70)

    t_max = 10.0
    step_counts = [5, 10, 20, 50]
    p_depol_values = [1e-4, 1e-3, 1e-2]

    print("  Running noiseless reference (n_steps=100)...", flush=True)
    sim_ref = QuditGKSLSimulator(params)
    r_ref = sim_ref.simulate(t_max=t_max, n_steps=100)
    rho_ref_finest = r_ref["rho_final"]

    all_results = []

    for p_depol in p_depol_values:
        print(f"\n  --- p_depol = {p_depol} ---")
        print(
            f"  {'n_steps':>7} | {'Trotter_infid':>14} | {'Noise_infid':>14} | "
            f"{'Total_infid':>14} | {'Noise/Trotter':>14}"
        )
        print("  " + "-" * 70)

        for ns in step_counts:
            t0 = time_module.time()

            sim_noiseless = QuditGKSLSimulator(params)
            r_noiseless = sim_noiseless.simulate(t_max=t_max, n_steps=ns)
            rho_noiseless = r_noiseless["rho_final"]

            sim_noisy = QuditGKSLNoisySimulator(
                params, p_depol=p_depol, depol_pair_only=True
            )
            r_noisy = sim_noisy.simulate(t_max=t_max, n_steps=ns)
            rho_noisy = r_noisy["rho_final"]

            elapsed = time_module.time() - t0

            f_trotter = quantum_fidelity(rho_ref_finest, rho_noiseless)
            f_noise = quantum_fidelity(rho_noiseless, rho_noisy)
            f_total = quantum_fidelity(rho_ref_finest, rho_noisy)

            trotter_infid = 1.0 - f_trotter
            noise_infid = 1.0 - f_noise
            total_infid = 1.0 - f_total
            ratio = noise_infid / trotter_infid if trotter_infid > 1e-15 else float("inf")

            pop_noisy = compute_populations_from_density_matrix(rho_noisy, params)

            entry = {
                "p_depol": p_depol,
                "n_steps": ns,
                "dt": t_max / ns,
                "total_noise_events": 12 * ns,
                "F_trotter": f_trotter,
                "F_noise": f_noise,
                "F_total": f_total,
                "trotter_infidelity": trotter_infid,
                "noise_infidelity": noise_infid,
                "total_infidelity": total_infid,
                "noise_to_trotter_ratio": ratio,
                "noisy_N_S1": pop_noisy["N_S1"],
                "elapsed": elapsed,
            }
            all_results.append(entry)

            print(
                f"  {ns:7d} | {trotter_infid:14.8f} | {noise_infid:14.8f} | "
                f"{total_infid:14.8f} | {ratio:14.1f}x"
            )

    return {
        "test": "noise_vs_trotter_separation",
        "t_max": t_max,
        "depol_pair_only": True,
        "p_depol_values": p_depol_values,
        "n_steps_values": step_counts,
        "results": all_results,
    }


# =========================================================================
# Test B: Noise budget decomposition (unchanged from iteration 16)
# =========================================================================
def test_noise_budget(params: GKSLPhysicalParameters) -> dict:
    """Decompose noise budget: Hamiltonian gates vs Lindblad channel gates."""
    print("\n" + "=" * 70)
    print("Test B: Noise Budget Decomposition (Qudit, 1 step)")
    print("=" * 70)

    t_max_step = 1.0
    p_depol = 0.01

    sim_noiseless = QuditGKSLSimulator(params)
    sim_noiseless._precompute_unitaries(t_max_step)
    rho_0 = sim_noiseless.prepare_initial_state("edge_triplet")
    rho_noiseless = sim_noiseless._trotter_step(rho_0.copy())

    sim_full = QuditGKSLNoisySimulator(
        params, p_depol=p_depol, depol_pair_only=True
    )
    sim_full._precompute_unitaries(t_max_step)
    rho_full_noise = sim_full._trotter_step(rho_0.copy())

    f_full = quantum_fidelity(rho_noiseless, rho_full_noise)

    from qudit_gksl_noisy_simulator import _apply_local_depolarization_pair
    from stinespring_utils import apply_stinespring_to_density_matrix

    d = params.d
    N = params.N_molecules

    rho_h_only = sim_noiseless._U_H_half @ rho_0.copy() @ sim_noiseless._U_H_half.conj().T
    for i, j in params.neighbors:
        rho_h_only = _apply_local_depolarization_pair(rho_h_only, i, j, d, N, p_depol)
    for U_stine in sim_noiseless._U_stines:
        rho_h_only = apply_stinespring_to_density_matrix(rho_h_only, U_stine)
    rho_h_only = sim_noiseless._U_H_half @ rho_h_only @ sim_noiseless._U_H_half.conj().T
    for i, j in params.neighbors:
        rho_h_only = _apply_local_depolarization_pair(rho_h_only, i, j, d, N, p_depol)

    f_h_only = quantum_fidelity(rho_noiseless, rho_h_only)

    rho_l_only = sim_noiseless._U_H_half @ rho_0.copy() @ sim_noiseless._U_H_half.conj().T
    lindblad_sites = sim_full._lindblad_sites
    for k, U_stine in enumerate(sim_noiseless._U_stines):
        rho_l_only = apply_stinespring_to_density_matrix(rho_l_only, U_stine)
        sites = lindblad_sites[k]
        if len(sites) == 2:
            rho_l_only = _apply_local_depolarization_pair(
                rho_l_only, sites[0], sites[1], d, N, p_depol
            )
    rho_l_only = sim_noiseless._U_H_half @ rho_l_only @ sim_noiseless._U_H_half.conj().T

    f_l_only = quantum_fidelity(rho_noiseless, rho_l_only)

    pop_noiseless = compute_populations_from_density_matrix(rho_noiseless, params)
    pop_full = compute_populations_from_density_matrix(rho_full_noise, params)
    pop_h_only = compute_populations_from_density_matrix(rho_h_only, params)
    pop_l_only = compute_populations_from_density_matrix(rho_l_only, params)

    infid_h = 1.0 - f_h_only
    infid_l = 1.0 - f_l_only
    f_composed = (1.0 - infid_h) * (1.0 - infid_l)

    print(f"\n  === 1-Step Fidelity vs Noiseless ===")
    print(f"  Full noise:              F = {f_full:.8f}  (12 events)")
    print(f"  H-gate noise only:       F = {f_h_only:.8f}  (6 events)")
    print(f"  L-channel noise only:    F = {f_l_only:.8f}  (6 events)")
    print(f"  Multiplicative estimate: F = {f_composed:.8f}  ((1-infid_H)(1-infid_L))")
    print(f"\n  === N_S1 Populations ===")
    print(f"  Noiseless:    N_S1 = {pop_noiseless['N_S1']:.8f}")
    print(f"  Full noise:   N_S1 = {pop_full['N_S1']:.8f}")
    print(f"  H-gate only:  N_S1 = {pop_h_only['N_S1']:.8f}")
    print(f"  L-ch only:    N_S1 = {pop_l_only['N_S1']:.8f}")

    return {
        "test": "noise_budget",
        "dt": t_max_step,
        "p_depol": p_depol,
        "depol_pair_only": True,
        "F_full_noise": f_full,
        "F_hamiltonian_noise_only": f_h_only,
        "F_lindblad_noise_only": f_l_only,
        "F_multiplicative_estimate": f_composed,
        "N_S1_noiseless": pop_noiseless["N_S1"],
        "N_S1_full_noise": pop_full["N_S1"],
        "N_S1_hamiltonian_only": pop_h_only["N_S1"],
        "N_S1_lindblad_only": pop_l_only["N_S1"],
        "noise_events_hamiltonian": 6,
        "noise_events_lindblad": 6,
        "noise_events_total": 12,
    }


# =========================================================================
# Test C: Per-step epsilon tracking (fixed dt-variation report)
# =========================================================================
def test_per_step_epsilon_tracking(params: GKSLPhysicalParameters) -> dict:
    """Track per-step fidelity ratio to demonstrate diminishing-epsilon effect.

    Fixed from iteration 17: dt-variation description is computed from data
    instead of hardcoded.
    """
    print("\n" + "=" * 70)
    print("Test C: Per-Step Epsilon Tracking (Qudit, pair-only, p_depol=0.01)")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    dt = t_max / n_steps
    p_depol = 0.01
    d = params.d

    infidelity_per_event = p_depol * (1.0 - 1.0 / (d * d))
    events_per_step = 12
    eps_theoretical = 1.0 - (1.0 - infidelity_per_event) ** events_per_step

    # --- Part 1: Per-step tracking at dt=1.0 ---
    sim_noiseless = QuditGKSLSimulator(params)
    sim_noiseless._precompute_unitaries(dt)
    sim_noisy = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    sim_noisy._precompute_unitaries(dt)

    rho_noiseless = sim_noiseless.prepare_initial_state("edge_triplet")
    rho_noisy = sim_noisy.prepare_initial_state("edge_triplet")

    print(f"\n  Part 1: Per-step epsilon at dt={dt}")
    print(f"  Theoretical eps (constant model): {eps_theoretical:.8f}")
    print(f"\n  {'Step':>4} | {'F(n)':>10} | {'eps(n)':>10} | {'F_ratio':>10} | "
          f"{'Simple_pred':>12} | {'Delta':>10}")
    print("  " + "-" * 70)

    step_data = []
    f_prev = 1.0
    for s in range(n_steps):
        rho_noiseless = sim_noiseless._trotter_step(rho_noiseless)
        rho_noisy = sim_noisy._trotter_step(rho_noisy)
        f_n = quantum_fidelity(rho_noiseless, rho_noisy)
        eps_n = 1.0 - f_n / f_prev if f_prev > 1e-10 else float("nan")
        f_ratio = f_n / f_prev if f_prev > 1e-10 else float("nan")

        if s == 0:
            eps_step1 = eps_n
        simple_pred = (1.0 - eps_step1) ** (s + 1)
        delta = f_n - simple_pred

        sd = {
            "step": s + 1,
            "time": (s + 1) * dt,
            "F_n": f_n,
            "eps_n": eps_n,
            "F_ratio": f_ratio,
            "simple_model_prediction": simple_pred,
            "delta_vs_simple": delta,
        }
        step_data.append(sd)

        print(
            f"  {s+1:4d} | {f_n:10.8f} | {eps_n:10.8f} | {f_ratio:10.8f} | "
            f"{simple_pred:12.8f} | {delta:+10.6f}"
        )
        f_prev = f_n

    # --- Part 2: Per-dt epsilon measurement ---
    print(f"\n  Part 2: Epsilon at different dt (1-step measurement, same initial state)")
    print(f"  {'dt':>8} | {'F_1step':>10} | {'eps':>10}")
    print("  " + "-" * 35)

    dt_values = [0.2, 0.5, 1.0, 2.0]
    per_dt_data = []

    rho_0 = QuditGKSLSimulator(params).prepare_initial_state("edge_triplet")

    for dt_val in dt_values:
        sim_nl = QuditGKSLSimulator(params)
        sim_nl._precompute_unitaries(dt_val)
        rho_nl = sim_nl._trotter_step(rho_0.copy())

        sim_ny = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
        sim_ny._precompute_unitaries(dt_val)
        rho_ny = sim_ny._trotter_step(rho_0.copy())

        f = quantum_fidelity(rho_nl, rho_ny)
        eps = 1.0 - f

        per_dt_data.append({"dt": dt_val, "F_1step": f, "eps": eps})
        print(f"  {dt_val:8.2f} | {f:10.8f} | {eps:10.8f}")

    # Summary — computed from data, not hardcoded
    eps_values = [d_item["eps"] for d_item in per_dt_data]
    eps_mean = sum(eps_values) / len(eps_values)
    dt_abs_variation = max(eps_values) - min(eps_values)
    dt_rel_variation = dt_abs_variation / eps_mean * 100 if eps_mean > 1e-15 else 0.0

    print(f"\n  === Summary ===")
    print(f"  eps(step 1) = {step_data[0]['eps_n']:.8f}")
    print(f"  eps(step {n_steps}) = {step_data[-1]['eps_n']:.8f}")
    eps_decrease = (1.0 - step_data[-1]["eps_n"] / step_data[0]["eps_n"]) * 100
    print(f"  Decrease: {eps_decrease:.1f}%")
    print(f"  Theoretical constant eps: {eps_theoretical:.8f}")
    print(f"  dt absolute variation: {dt_abs_variation:.8f}")
    print(f"  dt relative variation: {dt_rel_variation:.2f}%")
    print(f"  CONCLUSION: eps(n) is NOT constant — it decreases by ~{eps_decrease:.0f}% over {n_steps} steps")

    return {
        "test": "per_step_epsilon_tracking",
        "t_max": t_max,
        "n_steps": n_steps,
        "dt": dt,
        "p_depol": p_depol,
        "depol_pair_only": True,
        "eps_theoretical": eps_theoretical,
        "infidelity_per_event": infidelity_per_event,
        "events_per_step": events_per_step,
        "per_step_data": step_data,
        "per_dt_data": per_dt_data,
        "dt_abs_variation": dt_abs_variation,
        "dt_rel_variation_pct": dt_rel_variation,
    }


# =========================================================================
# Test D: Improved qubit fidelity tracking (unchanged from iteration 16)
# =========================================================================
def test_qubit_improved_fidelity(params: GKSLPhysicalParameters) -> dict:
    """Improved qubit fidelity tracking with both raw and normalized metrics."""
    print("\n" + "=" * 70)
    print("Test D: Improved Qubit Fidelity Tracking (pair-only, p_depol=0.01)")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    dt = t_max / n_steps
    p_depol = 0.01

    sim_qd_ref = QuditGKSLSimulator(params)
    sim_qd_ref._precompute_unitaries(dt)
    rho_qd_ref = sim_qd_ref.prepare_initial_state("edge_triplet")

    sim_qb_ref = QubitGKSLSimulator(params)
    sim_qb_ref._precompute_unitaries(dt)
    rho_qb_ref = sim_qb_ref.prepare_initial_state("edge_triplet")

    sim_qd_noisy = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    sim_qd_noisy._precompute_unitaries(dt)
    rho_qd_noisy = sim_qd_noisy.prepare_initial_state("edge_triplet")

    sim_qb_noisy = QubitGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    sim_qb_noisy._precompute_unitaries(dt)
    rho_qb_noisy = sim_qb_noisy.prepare_initial_state("edge_triplet")

    print(f"\n  {'Step':>4} | {'F_qd':>8} | {'F_qb_raw':>9} | {'F_qb_norm':>10} | "
          f"{'Tr_qb':>8} | {'P_forb':>8} | {'qd_N_S1':>9} | {'qb_N_S1_cond':>13}")
    print("  " + "-" * 90)

    step_data = []

    for s in range(n_steps + 1):
        if s > 0:
            rho_qd_ref = sim_qd_ref._trotter_step(rho_qd_ref)
            rho_qb_ref = sim_qb_ref._trotter_step(rho_qb_ref)
            rho_qd_noisy = sim_qd_noisy._trotter_step(rho_qd_noisy)
            rho_qb_noisy = sim_qb_noisy._trotter_step(rho_qb_noisy)

        f_qd = quantum_fidelity(rho_qd_ref, rho_qd_noisy)
        rho_qb_ref_qt = sim_qb_ref._extract_from_qubit_space(rho_qb_ref)
        rho_qb_noisy_qt = sim_qb_noisy._extract_from_qubit_space(rho_qb_noisy)

        f_qb_raw = quantum_fidelity(rho_qb_ref_qt, rho_qb_noisy_qt)
        f_qb_norm = quantum_fidelity(
            normalize_rho(rho_qb_ref_qt), normalize_rho(rho_qb_noisy_qt)
        )

        tr_qb = float(np.real(np.trace(rho_qb_noisy_qt)))
        p_forb = 1.0 - tr_qb

        pop_qd_noisy = compute_populations_from_density_matrix(rho_qd_noisy, params)
        pop_qb_noisy = compute_populations_from_density_matrix(rho_qb_noisy_qt, params)
        n_s1_qb_conditional = pop_qb_noisy["N_S1"] / tr_qb if tr_qb > 1e-10 else 0.0

        sd = {
            "step": s,
            "time": s * dt,
            "F_qudit": f_qd,
            "F_qubit_raw": f_qb_raw,
            "F_qubit_norm": f_qb_norm,
            "trace_qubit": tr_qb,
            "forbidden_pop": p_forb,
            "qudit_noisy_N_S1": pop_qd_noisy["N_S1"],
            "qubit_noisy_N_S1_raw": pop_qb_noisy["N_S1"],
            "qubit_noisy_N_S1_conditional": n_s1_qb_conditional,
        }
        step_data.append(sd)

        print(
            f"  {s:4d} | {f_qd:.6f} | {f_qb_raw:.7f} | {f_qb_norm:.8f} | "
            f"{tr_qb:.6f} | {p_forb:.6f} | {pop_qd_noisy['N_S1']:.7f} | "
            f"{n_s1_qb_conditional:.7f}"
        )

    return {
        "test": "qubit_improved_fidelity",
        "t_max": t_max,
        "n_steps": n_steps,
        "p_depol": p_depol,
        "depol_pair_only": True,
        "step_data": step_data,
    }


# =========================================================================
# Test E: Exact Liouvillian comparison (MODIFIED in iteration 21)
# =========================================================================
def _build_dissipator_superoperator(
    H_total: np.ndarray, lindblad_ops: list
) -> np.ndarray:
    """Build the dissipator-only superoperator L_D (no Hamiltonian part).

    L_D = sum_alpha [ conj(L_alpha) ⊗ L_alpha
                      - 0.5*(I ⊗ L†_alpha L_alpha + (L†_alpha L_alpha)^T ⊗ I) ]
    """
    dim = H_total.shape[0]
    I_dim = np.eye(dim, dtype=np.complex128)
    L_D = np.zeros((dim * dim, dim * dim), dtype=np.complex128)
    for item in lindblad_ops:
        if isinstance(item, tuple):
            L_op = np.asarray(item[0], dtype=np.complex128)
        else:
            L_op = np.asarray(item, dtype=np.complex128)
        LdL = L_op.conj().T @ L_op
        L_D += (
            np.kron(L_op.conj(), L_op)
            - 0.5 * np.kron(I_dim, LdL)
            - 0.5 * np.kron(LdL.T, I_dim)
        )
    return L_D


def _build_individual_dissipator_superoperators(
    H_total: np.ndarray, lindblad_ops: list
) -> list:
    """Build individual dissipator superoperators L_{D_alpha} for each Lindblad channel.

    Returns a list of superoperator matrices, one per Lindblad operator.
    Each L_{D_alpha} = conj(L_alpha) ⊗ L_alpha
                       - 0.5*(I ⊗ L†_alpha L_alpha + (L†_alpha L_alpha)^T ⊗ I)
    """
    dim = H_total.shape[0]
    I_dim = np.eye(dim, dtype=np.complex128)
    L_D_list = []
    for item in lindblad_ops:
        if isinstance(item, tuple):
            L_op = np.asarray(item[0], dtype=np.complex128)
        else:
            L_op = np.asarray(item, dtype=np.complex128)
        LdL = L_op.conj().T @ L_op
        L_D_alpha = (
            np.kron(L_op.conj(), L_op)
            - 0.5 * np.kron(I_dim, LdL)
            - 0.5 * np.kron(LdL.T, I_dim)
        )
        L_D_list.append(L_D_alpha)
    return L_D_list


def _classical_trotter_simulate(
    H_total: np.ndarray,
    L_D_super: np.ndarray,
    rho_0: np.ndarray,
    t_max: float,
    n_steps: int,
    dim: int,
) -> np.ndarray:
    """Run classical Trotter simulation: exp(L_H dt/2) exp(L_D dt) exp(L_H dt/2).

    Uses exact unitary for the Hamiltonian part (81x81 expm) and
    expm_multiply for the dissipator superoperator (6561-vector).
    This avoids computing full 6561x6561 matrix exponentials.
    """
    dt = t_max / n_steps
    U_H_half = expm(-1j * H_total * dt / 2)
    U_H_half_dag = U_H_half.conj().T

    rho = rho_0.copy()
    for _ in range(n_steps):
        # Half Hamiltonian (exact unitary, fast)
        rho = U_H_half @ rho @ U_H_half_dag
        # Full dissipator (exact exp(L_D dt), via expm_multiply)
        vec = vectorize_density_matrix(rho)
        vec = expm_multiply(L_D_super * dt, vec)
        rho = unvectorize_density_matrix(vec, dim)
        # Half Hamiltonian (exact unitary, fast)
        rho = U_H_half @ rho @ U_H_half_dag
    return rho


def _classical_product_trotter_simulate(
    H_total: np.ndarray,
    L_D_individual: list,
    rho_0: np.ndarray,
    t_max: float,
    n_steps: int,
    dim: int,
) -> np.ndarray:
    """Run classical product Trotter simulation with Strang splitting for H-D.

    Applies exp(L_H dt/2) · prod_alpha exp(L_{D_alpha} dt) · exp(L_H dt/2)
    where individual dissipator channels are applied sequentially in forward
    order within each dissipator step (not symmetrized/palindromic), while
    the overall H-D decomposition uses Strang splitting. This forward-only
    product of dissipator channels isolates the Lie-Trotter product error.
    """
    dt = t_max / n_steps
    U_H_half = expm(-1j * H_total * dt / 2)
    U_H_half_dag = U_H_half.conj().T

    rho = rho_0.copy()
    for _ in range(n_steps):
        # Half Hamiltonian (exact unitary, fast)
        rho = U_H_half @ rho @ U_H_half_dag
        # Product of individual dissipator channels (forward order)
        vec = vectorize_density_matrix(rho)
        for L_D_alpha in L_D_individual:
            vec = expm_multiply(L_D_alpha * dt, vec)
        rho = unvectorize_density_matrix(vec, dim)
        # Half Hamiltonian (exact unitary, fast)
        rho = U_H_half @ rho @ U_H_half_dag
    return rho


def _symmetric_classical_product_trotter_simulate(
    H_total: np.ndarray,
    L_D_individual: list,
    rho_0: np.ndarray,
    t_max: float,
    n_steps: int,
    dim: int,
) -> np.ndarray:
    """Symmetric (palindromic) classical product Trotter simulation.

    Applies exp(L_H dt/2) · [prod_{α=1..n} exp(L_{D_α} dt/2)]
                           · [prod_{α=n..1} exp(L_{D_α} dt/2)]
                           · exp(L_H dt/2)

    This matches the ST simulator's palindromic Lindblad channel ordering but
    uses exact channel exponentials instead of Stinespring approximations.
    Therefore ST vs SCPT cleanly isolates the Stinespring approximation error.
    """
    dt = t_max / n_steps
    U_H_half = expm(-1j * H_total * dt / 2)
    U_H_half_dag = U_H_half.conj().T

    rho = rho_0.copy()
    for _ in range(n_steps):
        # Half Hamiltonian
        rho = U_H_half @ rho @ U_H_half_dag
        # Forward half-step for all individual dissipators
        vec = vectorize_density_matrix(rho)
        for L_D_alpha in L_D_individual:
            vec = expm_multiply(L_D_alpha * (dt / 2), vec)
        # Reverse half-step (palindromic)
        for L_D_alpha in reversed(L_D_individual):
            vec = expm_multiply(L_D_alpha * (dt / 2), vec)
        rho = unvectorize_density_matrix(vec, dim)
        # Half Hamiltonian
        rho = U_H_half @ rho @ U_H_half_dag
    return rho


def test_exact_liouvillian_comparison(params: GKSLPhysicalParameters) -> dict:
    """Compare ST against exact GKSL, CT, CPT, and SCPT.

    MODIFIED from iteration 24:
    - Added Frobenius distance computation for key pairs
    - Part 6f (new): Composite model fitting T(ST,ex)/dt = p0 + p1*dt
    - Part 6g (new): Frobenius-based error vector angle estimation
    - Part 6h (new): Analytical Rate(T) prediction and comparison
    - Parts 6a-6e retained from iteration 24
    """
    print("\n" + "=" * 70)
    print("Test E: Exact Liouvillian Comparison (ST vs Exact vs CT vs CPT vs SCPT)")
    print("       [iter 24: extended dt range + higher-order analysis]")
    print("       [NOTE: ST uses symmetric Lindblad product ordering]")
    print("=" * 70)

    t_max = 10.0
    step_counts = [5, 10, 20, 50, 100]
    dim = params.d ** params.N_molecules

    # --- Part 1: Compute exact reference ---
    print("\n  Part 1: Computing exact GKSL Liouvillian solution...", flush=True)
    t0 = time_module.time()
    sim_exact = ClassicalGKSLSimulator(params)
    r_exact = sim_exact.simulate(t_max=t_max, n_steps=100)
    rho_exact_final = r_exact["rho_final"]
    exact_elapsed = time_module.time() - t0
    print(f"  Exact solution computed in {exact_elapsed:.2f}s")

    pop_exact = compute_populations_from_density_matrix(rho_exact_final, params)
    purity_exact = compute_purity(rho_exact_final)
    print(f"  Exact N_S1 = {pop_exact['N_S1']:.10f}")
    print(f"  Exact purity = {purity_exact:.10f}")

    # --- Part 1b: Build dissipator superoperators ---
    print("  Building dissipator superoperators for CT, CPT, and SCPT...", flush=True)
    H_total = sim_exact.H_total
    lindblad_ops = sim_exact.lindblad_ops
    L_D_super = _build_dissipator_superoperator(H_total, lindblad_ops)
    L_D_individual = _build_individual_dissipator_superoperators(H_total, lindblad_ops)
    rho_0 = sim_exact.prepare_initial_state("edge_triplet")
    print(f"  Done. {len(L_D_individual)} individual dissipator channels built.")

    # --- Part 2: Five-way convergence comparison ---
    print(f"\n  Part 2: Five-way convergence (primary metric: trace distance)")
    print(
        f"  {'n_steps':>7} | {'dt':>5} | "
        f"{'T(ST,ex)':>11} | {'T(CT,ex)':>11} | {'T(CPT,ex)':>11} | "
        f"{'T(SCPT,ex)':>11} | {'T(ST,SCPT)':>11} | {'T(ST,CT)':>11} | "
        f"{'T(SCPT,CT)':>11}"
    )
    print("  " + "-" * 115)

    convergence_data = []

    for ns in step_counts:
        dt = t_max / ns
        t0 = time_module.time()

        # Stinespring+Trotter
        sim_st = QuditGKSLSimulator(params)
        r_st = sim_st.simulate(t_max=t_max, n_steps=ns)
        rho_st = r_st["rho_final"]
        elapsed_st = time_module.time() - t0

        # Classical Trotter
        t0 = time_module.time()
        rho_ct = _classical_trotter_simulate(
            H_total, L_D_super, rho_0, t_max, ns, dim
        )
        elapsed_ct = time_module.time() - t0

        # Classical Product Trotter
        t0 = time_module.time()
        rho_cpt = _classical_product_trotter_simulate(
            H_total, L_D_individual, rho_0, t_max, ns, dim
        )
        elapsed_cpt = time_module.time() - t0

        # Symmetric Classical Product Trotter (new in iteration 22)
        t0 = time_module.time()
        rho_scpt = _symmetric_classical_product_trotter_simulate(
            H_total, L_D_individual, rho_0, t_max, ns, dim
        )
        elapsed_scpt = time_module.time() - t0

        # Purity tracking
        purity_st = compute_purity(rho_st)
        purity_ct = compute_purity(rho_ct)
        purity_cpt = compute_purity(rho_cpt)
        purity_scpt = compute_purity(rho_scpt)

        # --- Metrics: ST vs Exact ---
        f_st_ex = quantum_fidelity(rho_exact_final, rho_st)
        fidelity_clipped_st_ex = f_st_ex >= 1.0
        infid_st_ex = 0.0 if fidelity_clipped_st_ex else 1.0 - f_st_ex
        t_st_ex = trace_distance(rho_exact_final, rho_st)
        db_st_ex = bures_distance(rho_exact_final, rho_st)
        ratio_db_over_T_st_ex = db_st_ex / t_st_ex if t_st_ex > 1e-15 else float("nan")

        # --- Metrics: CT vs Exact ---
        f_ct_ex = quantum_fidelity(rho_exact_final, rho_ct)
        fidelity_clipped_ct_ex = f_ct_ex >= 1.0
        infid_ct_ex = 0.0 if fidelity_clipped_ct_ex else 1.0 - f_ct_ex
        t_ct_ex = trace_distance(rho_exact_final, rho_ct)
        db_ct_ex = bures_distance(rho_exact_final, rho_ct)

        # --- Metrics: CPT vs Exact ---
        f_cpt_ex = quantum_fidelity(rho_exact_final, rho_cpt)
        fidelity_clipped_cpt_ex = f_cpt_ex >= 1.0
        infid_cpt_ex = 0.0 if fidelity_clipped_cpt_ex else 1.0 - f_cpt_ex
        t_cpt_ex = trace_distance(rho_exact_final, rho_cpt)
        db_cpt_ex = bures_distance(rho_exact_final, rho_cpt)

        # --- Metrics: SCPT vs Exact (new in iteration 22) ---
        f_scpt_ex = quantum_fidelity(rho_exact_final, rho_scpt)
        fidelity_clipped_scpt_ex = f_scpt_ex >= 1.0
        infid_scpt_ex = 0.0 if fidelity_clipped_scpt_ex else 1.0 - f_scpt_ex
        t_scpt_ex = trace_distance(rho_exact_final, rho_scpt)
        db_scpt_ex = bures_distance(rho_exact_final, rho_scpt)

        # --- Metrics: ST vs CT ---
        f_st_ct = quantum_fidelity(rho_ct, rho_st)
        fidelity_clipped_st_ct = f_st_ct >= 1.0
        infid_st_ct = 0.0 if fidelity_clipped_st_ct else 1.0 - f_st_ct
        t_st_ct = trace_distance(rho_ct, rho_st)
        db_st_ct = bures_distance(rho_ct, rho_st)

        # --- Metrics: ST vs CPT ---
        f_st_cpt = quantum_fidelity(rho_cpt, rho_st)
        fidelity_clipped_st_cpt = f_st_cpt >= 1.0
        infid_st_cpt = 0.0 if fidelity_clipped_st_cpt else 1.0 - f_st_cpt
        t_st_cpt = trace_distance(rho_cpt, rho_st)
        db_st_cpt = bures_distance(rho_cpt, rho_st)

        # --- Metrics: ST vs SCPT (isolates Stinespring error — new in iter 22) ---
        f_st_scpt = quantum_fidelity(rho_scpt, rho_st)
        fidelity_clipped_st_scpt = f_st_scpt >= 1.0
        infid_st_scpt = 0.0 if fidelity_clipped_st_scpt else 1.0 - f_st_scpt
        t_st_scpt = trace_distance(rho_scpt, rho_st)
        db_st_scpt = bures_distance(rho_scpt, rho_st)

        # --- Metrics: CPT vs CT ---
        f_cpt_ct = quantum_fidelity(rho_ct, rho_cpt)
        fidelity_clipped_cpt_ct = f_cpt_ct >= 1.0
        infid_cpt_ct = 0.0 if fidelity_clipped_cpt_ct else 1.0 - f_cpt_ct
        t_cpt_ct = trace_distance(rho_ct, rho_cpt)
        db_cpt_ct = bures_distance(rho_ct, rho_cpt)

        # --- Metrics: SCPT vs CT (new in iteration 22) ---
        f_scpt_ct = quantum_fidelity(rho_ct, rho_scpt)
        fidelity_clipped_scpt_ct = f_scpt_ct >= 1.0
        infid_scpt_ct = 0.0 if fidelity_clipped_scpt_ct else 1.0 - f_scpt_ct
        t_scpt_ct = trace_distance(rho_ct, rho_scpt)
        db_scpt_ct = bures_distance(rho_ct, rho_scpt)

        # --- Frobenius distances (new in iteration 25) ---
        frob_st_ex = frobenius_distance(rho_exact_final, rho_st)
        frob_st_scpt = frobenius_distance(rho_scpt, rho_st)
        frob_scpt_ex = frobenius_distance(rho_exact_final, rho_scpt)

        fidelity_clipped = (
            fidelity_clipped_st_ex or fidelity_clipped_ct_ex or fidelity_clipped_cpt_ex
            or fidelity_clipped_scpt_ex or fidelity_clipped_st_ct
            or fidelity_clipped_st_cpt or fidelity_clipped_st_scpt
            or fidelity_clipped_cpt_ct or fidelity_clipped_scpt_ct
        )

        # Cancellation ratio: how much ST benefits from error cancellation
        # CT-based (original)
        sum_separate = t_st_ct + t_ct_ex
        cancel_ratio = t_st_ex / sum_separate if sum_separate > 1e-18 else float("nan")
        # SCPT-based (more principled: clean Stinespring + Strang decomposition)
        sum_via_scpt = t_st_scpt + t_scpt_ex
        cancel_ratio_scpt = t_st_ex / sum_via_scpt if sum_via_scpt > 1e-18 else float("nan")

        # Density matrix diagnostics for ST
        tr_val = float(np.real(np.trace(rho_st)))
        herm_err = float(np.linalg.norm(rho_st - rho_st.conj().T, "fro"))
        min_eig = float(np.linalg.eigvalsh(rho_st)[0])

        pop_st = compute_populations_from_density_matrix(rho_st, params)
        pop_ct = compute_populations_from_density_matrix(rho_ct, params)
        pop_cpt = compute_populations_from_density_matrix(rho_cpt, params)
        pop_scpt = compute_populations_from_density_matrix(rho_scpt, params)

        entry = {
            "n_steps": ns,
            "dt": dt,
            # ST vs Exact
            "F_ST_vs_exact": f_st_ex,
            "T_ST_vs_exact": t_st_ex,
            "infid_ST_vs_exact": infid_st_ex,
            "bures_ST_vs_exact": db_st_ex,
            "ratio_dB_over_T_ST_vs_exact": ratio_db_over_T_st_ex,
            "fidelity_clipped_ST_vs_exact": fidelity_clipped_st_ex,
            # CT vs Exact
            "F_CT_vs_exact": f_ct_ex,
            "T_CT_vs_exact": t_ct_ex,
            "infid_CT_vs_exact": infid_ct_ex,
            "bures_CT_vs_exact": db_ct_ex,
            "fidelity_clipped_CT_vs_exact": fidelity_clipped_ct_ex,
            # CPT vs Exact
            "F_CPT_vs_exact": f_cpt_ex,
            "T_CPT_vs_exact": t_cpt_ex,
            "infid_CPT_vs_exact": infid_cpt_ex,
            "bures_CPT_vs_exact": db_cpt_ex,
            "fidelity_clipped_CPT_vs_exact": fidelity_clipped_cpt_ex,
            # SCPT vs Exact (new)
            "F_SCPT_vs_exact": f_scpt_ex,
            "T_SCPT_vs_exact": t_scpt_ex,
            "infid_SCPT_vs_exact": infid_scpt_ex,
            "bures_SCPT_vs_exact": db_scpt_ex,
            "fidelity_clipped_SCPT_vs_exact": fidelity_clipped_scpt_ex,
            # ST vs CT
            "F_ST_vs_CT": f_st_ct,
            "T_ST_vs_CT": t_st_ct,
            "infid_ST_vs_CT": infid_st_ct,
            "bures_ST_vs_CT": db_st_ct,
            "fidelity_clipped_ST_vs_CT": fidelity_clipped_st_ct,
            # ST vs CPT
            "F_ST_vs_CPT": f_st_cpt,
            "T_ST_vs_CPT": t_st_cpt,
            "infid_ST_vs_CPT": infid_st_cpt,
            "bures_ST_vs_CPT": db_st_cpt,
            "fidelity_clipped_ST_vs_CPT": fidelity_clipped_st_cpt,
            # ST vs SCPT (isolates Stinespring error — new)
            "F_ST_vs_SCPT": f_st_scpt,
            "T_ST_vs_SCPT": t_st_scpt,
            "infid_ST_vs_SCPT": infid_st_scpt,
            "bures_ST_vs_SCPT": db_st_scpt,
            "fidelity_clipped_ST_vs_SCPT": fidelity_clipped_st_scpt,
            # CPT vs CT (forward Lie-Trotter product error)
            "F_CPT_vs_CT": f_cpt_ct,
            "T_CPT_vs_CT": t_cpt_ct,
            "infid_CPT_vs_CT": infid_cpt_ct,
            "bures_CPT_vs_CT": db_cpt_ct,
            "fidelity_clipped_CPT_vs_CT": fidelity_clipped_cpt_ct,
            # SCPT vs CT (palindromic residual — new)
            "F_SCPT_vs_CT": f_scpt_ct,
            "T_SCPT_vs_CT": t_scpt_ct,
            "infid_SCPT_vs_CT": infid_scpt_ct,
            "bures_SCPT_vs_CT": db_scpt_ct,
            "fidelity_clipped_SCPT_vs_CT": fidelity_clipped_scpt_ct,
            # Fidelity clipped (any)
            "fidelity_clipped": fidelity_clipped,
            # Cancellation ratio
            "cancellation_ratio": cancel_ratio,
            "cancellation_ratio_scpt": cancel_ratio_scpt,
            # Normalized error coefficients
            "T_ST_vs_SCPT_over_dt": t_st_scpt / dt,
            "T_SCPT_vs_exact_over_dt2": t_scpt_ex / (dt * dt),
            "T_CT_vs_exact_over_dt2": t_ct_ex / (dt * dt),
            "T_SCPT_vs_CT_over_dt2": t_scpt_ct / (dt * dt),
            # Purity
            "purity_exact": purity_exact,
            "purity_ST": purity_st,
            "purity_CT": purity_ct,
            "purity_CPT": purity_cpt,
            "purity_SCPT": purity_scpt,
            # Populations
            "N_S1_ST": pop_st["N_S1"],
            "N_S1_CT": pop_ct["N_S1"],
            "N_S1_CPT": pop_cpt["N_S1"],
            "N_S1_SCPT": pop_scpt["N_S1"],
            "N_S1_exact": pop_exact["N_S1"],
            # Signed N_S1 errors
            "delta_N_S1_ST": pop_st["N_S1"] - pop_exact["N_S1"],
            "delta_N_S1_CT": pop_ct["N_S1"] - pop_exact["N_S1"],
            "delta_N_S1_CPT": pop_cpt["N_S1"] - pop_exact["N_S1"],
            "delta_N_S1_SCPT": pop_scpt["N_S1"] - pop_exact["N_S1"],
            # Density matrix diagnostics
            "trace_ST": tr_val,
            "hermiticity_error_ST": herm_err,
            "min_eigenvalue_ST": min_eig,
            # Frobenius distances (new in iteration 25)
            "frob_ST_vs_exact": frob_st_ex,
            "frob_ST_vs_SCPT": frob_st_scpt,
            "frob_SCPT_vs_exact": frob_scpt_ex,
            # Timing
            "elapsed_ST": elapsed_st,
            "elapsed_CT": elapsed_ct,
            "elapsed_CPT": elapsed_cpt,
            "elapsed_SCPT": elapsed_scpt,
        }
        convergence_data.append(entry)

        print(
            f"  {ns:7d} | {dt:5.2f} | "
            f"{t_st_ex:11.4e} | {t_ct_ex:11.4e} | {t_cpt_ex:11.4e} | "
            f"{t_scpt_ex:11.4e} | {t_st_scpt:11.4e} | {t_st_ct:11.4e} | "
            f"{t_scpt_ct:11.4e}"
        )

    # --- Part 3: Convergence rate analysis ---
    print(f"\n  Part 3: Convergence rate analysis (primary metric: trace distance)")

    # 9 comparison pairs × 3 metrics each
    rates = {
        key: [] for key in [
            "infid_st_ex", "tdist_st_ex", "bures_st_ex",
            "infid_ct_ex", "tdist_ct_ex", "bures_ct_ex",
            "infid_cpt_ex", "tdist_cpt_ex", "bures_cpt_ex",
            "infid_scpt_ex", "tdist_scpt_ex", "bures_scpt_ex",
            "infid_st_ct", "tdist_st_ct", "bures_st_ct",
            "infid_st_cpt", "tdist_st_cpt", "bures_st_cpt",
            "infid_st_scpt", "tdist_st_scpt", "bures_st_scpt",
            "infid_cpt_ct", "tdist_cpt_ct", "bures_cpt_ct",
            "infid_scpt_ct", "tdist_scpt_ct", "bures_scpt_ct",
        ]
    }

    def _compute_rate(val_prev: float, val_curr: float, dt_ratio: float) -> float:
        if val_curr > 1e-18 and val_prev > 1e-18:
            return np.log(val_prev / val_curr) / np.log(dt_ratio)
        return float("nan")

    # Mapping from rate key to convergence_data key suffix
    pair_map = {
        "st_ex": ("ST_vs_exact", "ST_vs_exact", "ST_vs_exact"),
        "ct_ex": ("CT_vs_exact", "CT_vs_exact", "CT_vs_exact"),
        "cpt_ex": ("CPT_vs_exact", "CPT_vs_exact", "CPT_vs_exact"),
        "scpt_ex": ("SCPT_vs_exact", "SCPT_vs_exact", "SCPT_vs_exact"),
        "st_ct": ("ST_vs_CT", "ST_vs_CT", "ST_vs_CT"),
        "st_cpt": ("ST_vs_CPT", "ST_vs_CPT", "ST_vs_CPT"),
        "st_scpt": ("ST_vs_SCPT", "ST_vs_SCPT", "ST_vs_SCPT"),
        "cpt_ct": ("CPT_vs_CT", "CPT_vs_CT", "CPT_vs_CT"),
        "scpt_ct": ("SCPT_vs_CT", "SCPT_vs_CT", "SCPT_vs_CT"),
    }

    if len(convergence_data) >= 2:
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k - 1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr

            for pair_short, (infid_suffix, tdist_suffix, bures_suffix) in pair_map.items():
                rates[f"infid_{pair_short}"].append(_compute_rate(
                    convergence_data[k-1][f"infid_{infid_suffix}"],
                    convergence_data[k][f"infid_{infid_suffix}"], dt_ratio))
                rates[f"tdist_{pair_short}"].append(_compute_rate(
                    convergence_data[k-1][f"T_{tdist_suffix}"],
                    convergence_data[k][f"T_{tdist_suffix}"], dt_ratio))
                rates[f"bures_{pair_short}"].append(_compute_rate(
                    convergence_data[k-1][f"bures_{bures_suffix}"],
                    convergence_data[k][f"bures_{bures_suffix}"], dt_ratio))

        def _clip_marker(k_idx: int) -> str:
            """Return '*' if fidelity was clipped at either endpoint of the rate interval."""
            if convergence_data[k_idx]["fidelity_clipped"] or convergence_data[k_idx + 1]["fidelity_clipped"]:
                return "*"
            return ""

        pair_labels = [
            ("st_ex", "ST vs Exact"),
            ("ct_ex", "Classical Trotter vs Exact"),
            ("cpt_ex", "Classical Product Trotter vs Exact"),
            ("scpt_ex", "Symmetric Classical Product Trotter vs Exact"),
            ("st_ct", "ST vs Classical Trotter"),
            ("st_cpt", "ST vs CPT (Stinespring + ordering diff)"),
            ("st_scpt", "ST vs SCPT (pure Stinespring error)"),
            ("cpt_ct", "CPT vs CT (forward Lie-Trotter product error)"),
            ("scpt_ct", "SCPT vs CT (palindromic Lie-Trotter residual)"),
        ]

        for pair_short, pair_label in pair_labels:
            print(f"\n  === {pair_label} ===")
            print(f"  {'dt_range':>15} | {'Rate(1-F)':>10} | {'Rate(T)':>10} | {'Rate(d_B)':>10} | {'clip':>4}")
            print("  " + "-" * 60)
            for k in range(len(rates[f"infid_{pair_short}"])):
                dt_prev = convergence_data[k]["dt"]
                dt_curr = convergence_data[k + 1]["dt"]
                cm = _clip_marker(k)
                print(
                    f"  {dt_prev:.1f} → {dt_curr:.1f}      | "
                    f"{rates[f'infid_{pair_short}'][k]:10.2f}{cm} | "
                    f"{rates[f'tdist_{pair_short}'][k]:10.2f} | "
                    f"{rates[f'bures_{pair_short}'][k]:10.2f} | {cm:>4s}"
                )
            avg_infid = np.nanmean(rates[f"infid_{pair_short}"])
            avg_tdist = np.nanmean(rates[f"tdist_{pair_short}"])
            avg_bures = np.nanmean(rates[f"bures_{pair_short}"])
            print(f"  Average: {avg_infid:10.2f} | {avg_tdist:10.2f} | {avg_bures:10.2f}")

    # Compute averages dict for return value
    avg_rates = {}
    for key in rates:
        avg_rates[f"avg_rate_{key}"] = float(np.nanmean(rates[key])) if rates[key] else float("nan")

    # --- Part 4: Step-by-step comparison at n_steps=10 ---
    print(f"\n  Part 4: Step-by-step ST vs Exact at n_steps=10 (with Bures distance & purity)")

    dt_detail = t_max / 10
    sim_st_detail = QuditGKSLSimulator(params)
    sim_st_detail._precompute_unitaries(dt_detail)
    rho_st_detail = sim_st_detail.prepare_initial_state("edge_triplet")

    print(
        f"  {'Step':>4} | {'1-F':>10} | {'T_dist':>10} | {'d_B':>10} | {'d_B/T':>8} | "
        f"{'P_ST':>8} | {'P_ex':>8} | {'N_S1_ST':>10} | {'N_S1_exact':>10} | {'clip':>4}"
    )
    print("  " + "-" * 110)

    step_detail = []
    for s in range(11):
        if s > 0:
            rho_st_detail = sim_st_detail._trotter_step(rho_st_detail)

        if s == 0:
            rho_exact_step = sim_exact.prepare_initial_state("edge_triplet")
        else:
            vec_0 = vectorize_density_matrix(
                sim_exact.prepare_initial_state("edge_triplet")
            )
            vecs = expm_multiply(
                sim_exact._L_super, vec_0,
                start=0.0, stop=s * dt_detail,
                num=2, endpoint=True,
            )
            rho_exact_step = unvectorize_density_matrix(vecs[-1], sim_exact.dim)

        f_vs_exact_s = quantum_fidelity(rho_exact_step, rho_st_detail)
        fidelity_clipped_s = f_vs_exact_s >= 1.0
        infid_s = 0.0 if fidelity_clipped_s else 1.0 - f_vs_exact_s
        t_vs_exact_s = trace_distance(rho_exact_step, rho_st_detail)
        db_s = bures_distance(rho_exact_step, rho_st_detail)
        ratio_db_T_s = db_s / t_vs_exact_s if t_vs_exact_s > 1e-15 else float("nan")
        purity_st_s = compute_purity(rho_st_detail)
        purity_exact_s = compute_purity(rho_exact_step)
        pop_st_s = compute_populations_from_density_matrix(rho_st_detail, params)
        pop_exact_s = compute_populations_from_density_matrix(rho_exact_step, params)

        clip_s = "*" if fidelity_clipped_s else ""

        sd = {
            "step": s,
            "time": s * dt_detail,
            "F_vs_exact": f_vs_exact_s,
            "trace_dist_vs_exact": t_vs_exact_s,
            "infidelity": infid_s,
            "bures_dist": db_s,
            "ratio_dB_over_T": ratio_db_T_s,
            "fidelity_clipped": fidelity_clipped_s,
            "purity_ST": purity_st_s,
            "purity_exact": purity_exact_s,
            "N_S1_ST": pop_st_s["N_S1"],
            "N_S1_exact": pop_exact_s["N_S1"],
        }
        step_detail.append(sd)

        print(
            f"  {s:4d} | {infid_s:10.4e} | {t_vs_exact_s:10.4e} | "
            f"{db_s:10.4e} | {ratio_db_T_s:8.4f} | "
            f"{purity_st_s:8.6f} | {purity_exact_s:8.6f} | "
            f"{pop_st_s['N_S1']:10.8f} | {pop_exact_s['N_S1']:10.8f} | {clip_s:>4s}"
        )

    # --- Part 5: N_S1 convergence analysis (with signed errors) ---
    print(f"\n  Part 5: N_S1 convergence analysis (signed errors)")
    print(f"  {'n_steps':>7} | {'dt':>5} | {'Δ_ST':>12} | {'Δ_CT':>12} | "
          f"{'Δ_SCPT':>12} | {'Δ_CPT':>12} | {'cancel':>7}")
    print("  " + "-" * 85)

    ns1_diffs_st = []
    ns1_diffs_ct = []
    ns1_diffs_cpt = []
    ns1_diffs_scpt = []
    for cd in convergence_data:
        diff_st = abs(cd["N_S1_ST"] - cd["N_S1_exact"])
        diff_ct = abs(cd["N_S1_CT"] - cd["N_S1_exact"])
        diff_cpt = abs(cd["N_S1_CPT"] - cd["N_S1_exact"])
        diff_scpt = abs(cd["N_S1_SCPT"] - cd["N_S1_exact"])
        ns1_diffs_st.append(diff_st)
        ns1_diffs_ct.append(diff_ct)
        ns1_diffs_cpt.append(diff_cpt)
        ns1_diffs_scpt.append(diff_scpt)
        print(
            f"  {cd['n_steps']:7d} | {cd['dt']:5.2f} | "
            f"{cd['delta_N_S1_ST']:+12.4e} | {cd['delta_N_S1_CT']:+12.4e} | "
            f"{cd['delta_N_S1_SCPT']:+12.4e} | {cd['delta_N_S1_CPT']:+12.4e} | "
            f"{cd['cancellation_ratio']:7.3f}"
        )

    rates_ns1_st = []
    rates_ns1_ct = []
    rates_ns1_cpt = []
    rates_ns1_scpt = []
    if len(convergence_data) >= 2:
        print(f"\n  N_S1 convergence rates (|Δ|):")
        print(f"  {'dt_range':>15} | {'Rate(ST)':>10} | {'Rate(CT)':>10} | "
              f"{'Rate(SCPT)':>10} | {'Rate(CPT)':>10}")
        print("  " + "-" * 70)
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k-1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr
            rate_st = _compute_rate(ns1_diffs_st[k-1], ns1_diffs_st[k], dt_ratio)
            rate_ct = _compute_rate(ns1_diffs_ct[k-1], ns1_diffs_ct[k], dt_ratio)
            rate_cpt = _compute_rate(ns1_diffs_cpt[k-1], ns1_diffs_cpt[k], dt_ratio)
            rate_scpt = _compute_rate(ns1_diffs_scpt[k-1], ns1_diffs_scpt[k], dt_ratio)
            rates_ns1_st.append(rate_st)
            rates_ns1_ct.append(rate_ct)
            rates_ns1_cpt.append(rate_cpt)
            rates_ns1_scpt.append(rate_scpt)
            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | {rate_st:10.2f} | "
                f"{rate_ct:10.2f} | {rate_scpt:10.2f} | {rate_cpt:10.2f}"
            )
        print(
            f"  Average: {np.nanmean(rates_ns1_st):10.2f} | "
            f"{np.nanmean(rates_ns1_ct):10.2f} | "
            f"{np.nanmean(rates_ns1_scpt):10.2f} | "
            f"{np.nanmean(rates_ns1_cpt):10.2f}"
        )

    # --- Part 6: Polynomial fitting ---
    # Part 6a: Aggregate fitting T(ST,exact) ≈ a*dt + b*dt² (reference only)
    print(f"\n  Part 6a: Aggregate fitting T(ST,exact) ≈ a·dt + b·dt² (reference — see 6b for accurate fits)")
    dts = np.array([cd["dt"] for cd in convergence_data])
    t_st_vals = np.array([cd["T_ST_vs_exact"] for cd in convergence_data])
    # Fit T = a*dt + b*dt² using least squares: T/dt = a + b*dt
    if len(dts) >= 2:
        A_mat = np.column_stack([np.ones_like(dts), dts])
        b_vec = t_st_vals / dts
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, b_vec, rcond=None)
        a_fit, b_fit = coeffs
        t_fitted = a_fit * dts + b_fit * dts**2
        fit_errors = np.abs(t_st_vals - t_fitted) / t_st_vals
        print(f"  Aggregate fit: T ≈ {a_fit:.6e} · dt + {b_fit:.6e} · dt²")
        print(f"  Relative fit errors: {', '.join(f'{e:.4f}' for e in fit_errors)}")
        print(f"  ⚠ Large errors due to error cancellation — use component-wise fits below")
        polyfit_result = {
            "a_coefficient": float(a_fit),
            "b_coefficient": float(b_fit),
            "relative_fit_errors": [float(e) for e in fit_errors],
        }
    else:
        polyfit_result = {}

    # Part 6b: Component-wise fitting (accurate)
    print(f"\n  Part 6b: Component-wise fitting (accurate decomposition via SCPT)")
    t_stine_vals = np.array([cd["T_ST_vs_SCPT"] for cd in convergence_data])
    t_strang_vals = np.array([cd["T_SCPT_vs_exact"] for cd in convergence_data])
    t_palindromic_vals = np.array([cd["T_SCPT_vs_CT"] for cd in convergence_data])

    # Stinespring error: T(ST,SCPT) ≈ a_stine * dt
    a_stine = float(np.mean(t_stine_vals / dts))
    stine_fitted = a_stine * dts
    stine_fit_errors = np.abs(t_stine_vals - stine_fitted) / t_stine_vals
    print(f"  Stinespring: T(ST,SCPT) ≈ {a_stine:.6e} · dt")
    print(f"  Relative errors: {', '.join(f'{e:.4f}' for e in stine_fit_errors)}")

    # Strang error: T(SCPT,exact) ≈ b_strang * dt^2
    b_strang = float(np.mean(t_strang_vals / dts**2))
    strang_fitted = b_strang * dts**2
    strang_fit_errors = np.abs(t_strang_vals - strang_fitted) / t_strang_vals
    print(f"  Strang+palindromic: T(SCPT,ex) ≈ {b_strang:.6e} · dt²")
    print(f"  Relative errors: {', '.join(f'{e:.4f}' for e in strang_fit_errors)}")

    # Palindromic residual: T(SCPT,CT) ≈ c_palindromic * dt^2
    c_palindromic = float(np.mean(t_palindromic_vals / dts**2))
    palindromic_fitted = c_palindromic * dts**2
    palindromic_fit_errors = np.abs(t_palindromic_vals - palindromic_fitted) / t_palindromic_vals
    print(f"  Palindromic residual: T(SCPT,CT) ≈ {c_palindromic:.6e} · dt²")
    print(f"  Relative errors: {', '.join(f'{e:.4f}' for e in palindromic_fit_errors)}")

    print(f"\n  Error magnitude ratios:")
    print(f"    a_stine / b_strang = {a_stine / b_strang:.2f} "
          f"(Stinespring / Strang)")
    print(f"    b_strang / c_palindromic = {b_strang / c_palindromic:.0f} "
          f"(Strang / palindromic residual)")
    print(f"    Stinespring dominates for dt < {b_strang / a_stine:.2f}")

    component_fit_result = {
        "a_stine_coefficient": a_stine,
        "a_stine_relative_errors": [float(e) for e in stine_fit_errors],
        "b_strang_coefficient": b_strang,
        "b_strang_relative_errors": [float(e) for e in strang_fit_errors],
        "c_palindromic_coefficient": c_palindromic,
        "c_palindromic_relative_errors": [float(e) for e in palindromic_fit_errors],
    }

    # Part 6b': Higher-order fitting T = c0*dt^n + c1*dt^(n+1) (new in iteration 24)
    print(f"\n  Part 6b': Higher-order fitting T = c0·dt^n + c1·dt^(n+1)")
    higher_order_fit_result = {}

    if len(dts) >= 3:
        # Stinespring: T(ST,SCPT)/dt = a0 + a1*dt (linear regression)
        A_ho = np.column_stack([np.ones_like(dts), dts])
        coeffs_stine_ho, _, _, _ = np.linalg.lstsq(A_ho, t_stine_vals / dts, rcond=None)
        a0_stine, a1_stine = coeffs_stine_ho
        stine_ho_fitted = (a0_stine + a1_stine * dts) * dts
        stine_ho_fit_errors = np.abs(t_stine_vals - stine_ho_fitted) / t_stine_vals
        print(f"  Stinespring: T(ST,SCPT) = ({a0_stine:.6e} + {a1_stine:.6e}·dt) · dt")
        print(f"    a0 (asymptotic) = {a0_stine:.6e}, a1 (correction) = {a1_stine:.6e}")
        print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in stine_ho_fit_errors)}")

        # Strang: T(SCPT,ex)/dt^2 = b0 + b1*dt (linear regression)
        coeffs_strang_ho, _, _, _ = np.linalg.lstsq(A_ho, t_strang_vals / dts**2, rcond=None)
        b0_strang, b1_strang = coeffs_strang_ho
        strang_ho_fitted = (b0_strang + b1_strang * dts) * dts**2
        strang_ho_fit_errors = np.abs(t_strang_vals - strang_ho_fitted) / t_strang_vals
        print(f"  Strang: T(SCPT,ex) = ({b0_strang:.6e} + {b1_strang:.6e}·dt) · dt²")
        print(f"    b0 (asymptotic) = {b0_strang:.6e}, b1 (correction) = {b1_strang:.6e}")
        print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in strang_ho_fit_errors)}")

        # Palindromic: T(SCPT,CT)/dt^2 = c0 + c1*dt (linear regression)
        coeffs_pal_ho, _, _, _ = np.linalg.lstsq(A_ho, t_palindromic_vals / dts**2, rcond=None)
        c0_palindromic, c1_palindromic = coeffs_pal_ho
        pal_ho_fitted = (c0_palindromic + c1_palindromic * dts) * dts**2
        pal_ho_fit_errors = np.abs(t_palindromic_vals - pal_ho_fitted) / t_palindromic_vals
        print(f"  Palindromic: T(SCPT,CT) = ({c0_palindromic:.6e} + {c1_palindromic:.6e}·dt) · dt²")
        print(f"    c0 (asymptotic) = {c0_palindromic:.6e}, c1 (correction) = {c1_palindromic:.6e}")
        print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in pal_ho_fit_errors)}")

        print(f"\n  Comparison: mean vs asymptotic coefficients")
        print(f"    Stinespring: mean={a_stine:.6e}, asymptotic c0={a0_stine:.6e}, diff={abs(a_stine-a0_stine)/a_stine*100:.3f}%")
        print(f"    Strang: mean={b_strang:.6e}, asymptotic c0={b0_strang:.6e}, diff={abs(b_strang-b0_strang)/b_strang*100:.3f}%")
        print(f"    Palindromic: mean={c_palindromic:.6e}, asymptotic c0={c0_palindromic:.6e}, diff={abs(c_palindromic-c0_palindromic)/c_palindromic*100:.3f}%")

        higher_order_fit_result = {
            "stinespring_a0": float(a0_stine),
            "stinespring_a1": float(a1_stine),
            "stinespring_ho_relative_errors": [float(e) for e in stine_ho_fit_errors],
            "strang_b0": float(b0_strang),
            "strang_b1": float(b1_strang),
            "strang_ho_relative_errors": [float(e) for e in strang_ho_fit_errors],
            "palindromic_c0": float(c0_palindromic),
            "palindromic_c1": float(c1_palindromic),
            "palindromic_ho_relative_errors": [float(e) for e in pal_ho_fit_errors],
        }

    # Part 6c: Normalized error coefficients table
    print(f"\n  Part 6c: Normalized error coefficients")
    print(f"  {'dt':>5} | {'T(ST,SCPT)/dt':>14} | {'T(SCPT,ex)/dt²':>15} | "
          f"{'T(CT,ex)/dt²':>13} | {'T(SCPT,CT)/dt²':>15} | "
          f"{'cancel_CT':>10} | {'cancel_SCPT':>12}")
    print("  " + "-" * 105)
    for cd_entry in convergence_data:
        print(
            f"  {cd_entry['dt']:5.2f} | "
            f"{cd_entry['T_ST_vs_SCPT_over_dt']:14.6e} | "
            f"{cd_entry['T_SCPT_vs_exact_over_dt2']:15.6e} | "
            f"{cd_entry['T_CT_vs_exact_over_dt2']:13.6e} | "
            f"{cd_entry['T_SCPT_vs_CT_over_dt2']:15.6e} | "
            f"{cd_entry['cancellation_ratio']:10.4f} | "
            f"{cd_entry['cancellation_ratio_scpt']:12.4f}"
        )

    # Part 6d: Error vector angle estimation (new in iteration 24)
    print(f"\n  Part 6d: Error vector angle estimation (trace distance law of cosines)")
    print(f"  Using cos(θ) ≈ (T²(ST,ex) - T²(ST,SCPT) - T²(SCPT,ex)) / (2·T(ST,SCPT)·T(SCPT,ex))")
    print(f"  {'dt':>5} | {'T(ST,ex)':>11} | {'T(ST,SCPT)':>11} | {'T(SCPT,ex)':>11} | "
          f"{'cos(θ)':>8} | {'θ (deg)':>8} | {'Stine/Strang':>12}")
    print("  " + "-" * 85)

    error_angle_data = []
    for cd_entry in convergence_data:
        t_ab = cd_entry["T_ST_vs_exact"]
        t_a = cd_entry["T_ST_vs_SCPT"]
        t_b = cd_entry["T_SCPT_vs_exact"]
        # Approximate angle using law of cosines with trace distances
        denominator = 2.0 * t_a * t_b
        if denominator > 1e-30:
            cos_theta = (t_ab**2 - t_a**2 - t_b**2) / denominator
            cos_theta = max(-1.0, min(1.0, cos_theta))
            theta_deg = float(np.degrees(np.arccos(cos_theta)))
        else:
            cos_theta = float("nan")
            theta_deg = float("nan")
        stine_over_strang = t_a / t_b if t_b > 1e-30 else float("nan")

        error_angle_data.append({
            "dt": cd_entry["dt"],
            "cos_theta_approx": cos_theta,
            "theta_deg_approx": theta_deg,
            "stinespring_over_strang_ratio": stine_over_strang,
        })
        print(
            f"  {cd_entry['dt']:5.2f} | {t_ab:11.4e} | {t_a:11.4e} | {t_b:11.4e} | "
            f"{cos_theta:8.4f} | {theta_deg:8.1f} | {stine_over_strang:12.2f}"
        )

    # Part 6e: Richardson extrapolation for asymptotic coefficients (new in iteration 24)
    print(f"\n  Part 6e: Richardson extrapolation (dt→0 limit from two smallest dt values)")
    richardson_result = {}
    if len(dts) >= 2:
        # Use the two smallest dt values
        dt_small = dts[-1]
        dt_next = dts[-2]
        r0_stine_small = t_stine_vals[-1] / dt_small
        r0_stine_next = t_stine_vals[-2] / dt_next
        a0_rich = r0_stine_small - dt_small * (r0_stine_next - r0_stine_small) / (dt_next - dt_small)

        r0_strang_small = t_strang_vals[-1] / dt_small**2
        r0_strang_next = t_strang_vals[-2] / dt_next**2
        b0_rich = r0_strang_small - dt_small * (r0_strang_next - r0_strang_small) / (dt_next - dt_small)

        r0_pal_small = t_palindromic_vals[-1] / dt_small**2
        r0_pal_next = t_palindromic_vals[-2] / dt_next**2
        c0_rich = r0_pal_small - dt_small * (r0_pal_next - r0_pal_small) / (dt_next - dt_small)

        print(f"  Using dt={dt_small:.2f} and dt={dt_next:.2f}:")
        print(f"    Stinespring a0_Richardson = {a0_rich:.6e} (mean = {a_stine:.6e}, diff = {abs(a0_rich-a_stine)/a_stine*100:.3f}%)")
        print(f"    Strang b0_Richardson = {b0_rich:.6e} (mean = {b_strang:.6e}, diff = {abs(b0_rich-b_strang)/b_strang*100:.3f}%)")
        print(f"    Palindromic c0_Richardson = {c0_rich:.6e} (mean = {c_palindromic:.6e}, diff = {abs(c0_rich-c_palindromic)/c_palindromic*100:.3f}%)")

        richardson_result = {
            "dt_pair": [float(dt_small), float(dt_next)],
            "a0_stinespring_richardson": float(a0_rich),
            "b0_strang_richardson": float(b0_rich),
            "c0_palindromic_richardson": float(c0_rich),
        }

    # Part 6f: Composite model fitting T(ST,ex) = p0*dt + p1*dt^2 (new in iteration 25)
    print(f"\n  Part 6f: Composite model fitting T(ST,ex)/dt = p₀ + p₁·dt")
    composite_fit_result = {}
    if len(dts) >= 3:
        A_comp = np.column_stack([np.ones_like(dts), dts])
        coeffs_comp, _, _, _ = np.linalg.lstsq(A_comp, t_st_vals / dts, rcond=None)
        p0_comp, p1_comp = coeffs_comp
        t_comp_fitted = (p0_comp + p1_comp * dts) * dts
        comp_fit_errors = np.abs(t_st_vals - t_comp_fitted) / t_st_vals
        print(f"  T(ST,ex) = ({p0_comp:.6e} + {p1_comp:.6e}·dt) · dt")
        print(f"  p₀ (effective leading coeff) = {p0_comp:.6e}")
        print(f"  p₁ (cancellation effect)     = {p1_comp:.6e} {'(< 0 → cancellation)' if p1_comp < 0 else '(> 0 → no cancellation at large dt)'}")
        print(f"  Relative fit errors: {', '.join(f'{e:.6f}' for e in comp_fit_errors)}")
        print(f"  Comparison with Stinespring a₀:")
        print(f"    a₀ (Stinespring asymptotic) = {a0_stine:.6e}")
        print(f"    p₀ (composite leading)      = {p0_comp:.6e}")
        print(f"    diff = {abs(p0_comp - a0_stine) / a0_stine * 100:.2f}%")
        if p1_comp < 0:
            print(f"  p₁ < 0 confirms error cancellation between Stinespring and Strang components")
        composite_fit_result = {
            "p0_leading": float(p0_comp),
            "p1_cancellation": float(p1_comp),
            "relative_fit_errors": [float(e) for e in comp_fit_errors],
        }

    # Part 6g: Frobenius-based error vector angle estimation (new in iteration 25)
    print(f"\n  Part 6g: Frobenius-based error vector angle estimation (exact L₂ law of cosines)")
    print(f"  cos(θ_F) = (d_F²(ST,ex) - d_F²(ST,SCPT) - d_F²(SCPT,ex)) / (2·d_F(ST,SCPT)·d_F(SCPT,ex))")
    print(f"  {'dt':>5} | {'d_F(ST,ex)':>11} | {'d_F(ST,SCPT)':>13} | {'d_F(SCPT,ex)':>13} | "
          f"{'cos(θ_F)':>9} | {'θ_F(deg)':>9} | {'cos(θ_T)':>9} | {'θ_T(deg)':>9} | {'Δθ':>5}")
    print("  " + "-" * 115)

    frobenius_angle_data = []
    for k_idx, cd_entry in enumerate(convergence_data):
        df_ab = cd_entry["frob_ST_vs_exact"]
        df_a = cd_entry["frob_ST_vs_SCPT"]
        df_b = cd_entry["frob_SCPT_vs_exact"]
        # Exact law of cosines (Frobenius norm is Hilbert-Schmidt inner product norm)
        denom_f = 2.0 * df_a * df_b
        if denom_f > 1e-30:
            cos_theta_f = (df_ab**2 - df_a**2 - df_b**2) / denom_f
            cos_theta_f = max(-1.0, min(1.0, cos_theta_f))
            theta_f_deg = float(np.degrees(np.arccos(cos_theta_f)))
        else:
            cos_theta_f = float("nan")
            theta_f_deg = float("nan")

        # Trace distance angle for comparison
        cos_theta_t = error_angle_data[k_idx]["cos_theta_approx"]
        theta_t_deg = error_angle_data[k_idx]["theta_deg_approx"]
        delta_theta = theta_f_deg - theta_t_deg if not (np.isnan(theta_f_deg) or np.isnan(theta_t_deg)) else float("nan")

        frobenius_angle_data.append({
            "dt": cd_entry["dt"],
            "frob_ST_ex": df_ab,
            "frob_ST_SCPT": df_a,
            "frob_SCPT_ex": df_b,
            "cos_theta_frobenius": cos_theta_f,
            "theta_frobenius_deg": theta_f_deg,
            "cos_theta_trace": cos_theta_t,
            "theta_trace_deg": theta_t_deg,
            "delta_theta_deg": delta_theta,
        })
        print(
            f"  {cd_entry['dt']:5.2f} | {df_ab:11.4e} | {df_a:13.4e} | {df_b:13.4e} | "
            f"{cos_theta_f:9.4f} | {theta_f_deg:9.1f} | {cos_theta_t:9.4f} | {theta_t_deg:9.1f} | "
            f"{delta_theta:+5.1f}"
        )

    # Stability analysis for Frobenius angle
    theta_f_vals = [a["theta_frobenius_deg"] for a in frobenius_angle_data if not np.isnan(a["theta_frobenius_deg"])]
    theta_t_vals = [a["theta_trace_deg"] for a in frobenius_angle_data if not np.isnan(a["theta_trace_deg"])]
    if theta_f_vals:
        print(f"\n  Frobenius angle range: {min(theta_f_vals):.1f}° – {max(theta_f_vals):.1f}° "
              f"(variation: {max(theta_f_vals)-min(theta_f_vals):.1f}°)")
    if theta_t_vals:
        print(f"  Trace dist angle range: {min(theta_t_vals):.1f}° – {max(theta_t_vals):.1f}° "
              f"(variation: {max(theta_t_vals)-min(theta_t_vals):.1f}°)")

    # Part 6h: Analytical Rate(T) prediction (new in iteration 25)
    print(f"\n  Part 6h: Analytical Rate(T) prediction for ST vs Exact")
    print(f"  Using composite model T(ST,ex) ≈ p₀·dt + p₁·dt² to predict Rate(T)")
    rate_prediction_data = []
    if composite_fit_result and len(convergence_data) >= 2:
        p0_val = composite_fit_result["p0_leading"]
        p1_val = composite_fit_result["p1_cancellation"]
        print(f"  Model: T(ST,ex) ≈ {p0_val:.6e}·dt + {p1_val:.6e}·dt²")
        print(f"  {'dt range':>15} | {'T_pred(dt₁)':>12} | {'T_pred(dt₂)':>12} | "
              f"{'Rate_pred':>10} | {'Rate_obs':>10} | {'Δ_Rate':>8}")
        print("  " + "-" * 85)

        observed_rates = rates.get("tdist_st_ex", [])
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k - 1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr

            # Predicted T from composite model
            t_pred_prev = p0_val * dt_prev + p1_val * dt_prev**2
            t_pred_curr = p0_val * dt_curr + p1_val * dt_curr**2

            if t_pred_prev > 0 and t_pred_curr > 0:
                rate_predicted = np.log(t_pred_prev / t_pred_curr) / np.log(dt_ratio)
            else:
                rate_predicted = float("nan")

            rate_observed = observed_rates[k - 1] if k - 1 < len(observed_rates) else float("nan")
            delta_rate = rate_predicted - rate_observed if not (np.isnan(rate_predicted) or np.isnan(rate_observed)) else float("nan")

            rate_prediction_data.append({
                "dt_prev": dt_prev,
                "dt_curr": dt_curr,
                "T_predicted_prev": t_pred_prev,
                "T_predicted_curr": t_pred_curr,
                "rate_predicted": float(rate_predicted),
                "rate_observed": float(rate_observed),
                "delta_rate": float(delta_rate),
            })

            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | {t_pred_prev:12.4e} | {t_pred_curr:12.4e} | "
                f"{rate_predicted:10.4f} | {rate_observed:10.4f} | {delta_rate:+8.4f}"
            )

        # Predict dt threshold for Rate > 0.95
        if p0_val > 0 and p1_val < 0:
            # Rate ≈ 1 + p1/p0 * dt (approximate for small dt)
            # Rate > 0.95  →  p1/p0 * dt > -0.05  →  dt < 0.05 * p0 / |p1|
            dt_threshold_095 = 0.05 * p0_val / abs(p1_val)
            # Rate > 0.99  →  dt < 0.01 * p0 / |p1|
            dt_threshold_099 = 0.01 * p0_val / abs(p1_val)
            print(f"\n  Predicted dt thresholds (using Rate ≈ 1 + (p₁/p₀)·dt):")
            print(f"    Rate(T) > 0.95 requires dt < {dt_threshold_095:.4f} (n_steps > {10.0/dt_threshold_095:.0f})")
            print(f"    Rate(T) > 0.99 requires dt < {dt_threshold_099:.4f} (n_steps > {10.0/dt_threshold_099:.0f})")
            print(f"    p₁/p₀ = {p1_val/p0_val:.4f}")
        else:
            dt_threshold_095 = float("nan")
            dt_threshold_099 = float("nan")

        rate_prediction_data_result = {
            "model_p0": float(p0_val),
            "model_p1": float(p1_val),
            "predictions": rate_prediction_data,
            "dt_threshold_rate_095": float(dt_threshold_095) if not np.isnan(dt_threshold_095) else None,
            "dt_threshold_rate_099": float(dt_threshold_099) if not np.isnan(dt_threshold_099) else None,
        }
    else:
        rate_prediction_data_result = {}


    return {
        "test": "exact_liouvillian_comparison",
        "t_max": t_max,
        "n_lindblad_channels": len(L_D_individual),
        "exact_N_S1": pop_exact["N_S1"],
        "exact_purity": purity_exact,
        "exact_elapsed": exact_elapsed,
        "convergence_data": convergence_data,
        # All convergence rates (flat)
        "rates": {k: v for k, v in rates.items()},
        # Average rates
        **avg_rates,
        # N_S1 convergence
        "rates_N_S1_ST": rates_ns1_st,
        "rates_N_S1_CT": rates_ns1_ct,
        "rates_N_S1_CPT": rates_ns1_cpt,
        "rates_N_S1_SCPT": rates_ns1_scpt,
        # Polynomial fitting (aggregate — reference only)
        "polyfit_aggregate": polyfit_result,
        # Component-wise fitting (accurate)
        "polyfit_component": component_fit_result,
        # Higher-order fitting (new in iteration 24)
        "polyfit_higher_order": higher_order_fit_result,
        # Error vector angle estimation (new in iteration 24)
        "error_angle_data": error_angle_data,
        # Richardson extrapolation (new in iteration 24)
        "richardson_extrapolation": richardson_result,
        # Composite model fitting (new in iteration 25)
        "composite_fit": composite_fit_result,
        # Frobenius angle estimation (new in iteration 25)
        "frobenius_angle_data": frobenius_angle_data,
        # Rate(T) prediction (new in iteration 25)
        "rate_prediction": rate_prediction_data_result,
        # Step detail
        "step_detail_n10": step_detail,
    }


# =========================================================================
# Main
# =========================================================================
def main() -> int:
    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    output_dir = Path(__file__).resolve().parents[1] / "developing" / "verification_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    params = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
    )

    print("=" * 70)
    print("TTA-UC GKSL Iteration 25: Frobenius Analysis and Composite Rate(T) Prediction")
    print(f"Timestamp: {timestamp}")
    print("NOTE: ST simulator uses symmetric Lindblad product ordering")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 25,
        "purpose": "Frobenius analysis and composite Rate(T) prediction. "
                   "Iteration 24 showed Rate(T)=0.87 for ST vs Exact at dt=0.2->0.1 "
                   "(predicted 0.95-1.05). This is NOT a bug but a mathematical consequence "
                   "of error cancellation. This iteration adds: "
                   "(1) Frobenius distance for exact L2 angle estimation, "
                   "(2) composite model fitting T(ST,ex) = p0*dt + p1*dt^2, "
                   "(3) analytical Rate(T) prediction.",
        "simulator_note": "ST uses symmetric (palindromic) Lindblad product ordering. "
                          "SCPT uses the same palindromic ordering with exact channel exponentials.",
        "tests": [],
    }

    # Test A: Noise vs Trotter separation
    ta = test_noise_vs_trotter_separation(params)
    results["tests"].append(ta)

    # Test B: Noise budget decomposition
    tb = test_noise_budget(params)
    results["tests"].append(tb)

    # Test C: Per-step epsilon tracking (fixed reporting)
    tc = test_per_step_epsilon_tracking(params)
    results["tests"].append(tc)

    # Test D: Improved qubit fidelity tracking
    td = test_qubit_improved_fidelity(params)
    results["tests"].append(td)

    # Test E: Exact Liouvillian comparison (MODIFIED with SCPT)
    te = test_exact_liouvillian_comparison(params)
    results["tests"].append(te)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration25_frobenius_analysis_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 25: Frobenius 解析と合成 Rate(T) 予測",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "- iteration 24 からの主要修正:",
        "  - Frobenius 距離の追加計算（L₂ ノルム）",
        "  - Part 6f: 合成モデルフィッティング T(ST,ex)/dt = p₀ + p₁·dt",
        "  - Part 6g: Frobenius ベースの誤差ベクトル角度推定（厳密な余弦定理）",
        "  - Part 6h: 合成 Rate(T) 予測と実測の比較",
        "  - Parts 6a-6e: iteration 24 と同一",
        "",
    ]

    # Test A summary
    md_lines.extend([
        "## Test A: ノイズ vs Trotter 分離",
        "",
    ])
    for p_depol in ta["p_depol_values"]:
        md_lines.extend([
            f"### p_depol = {p_depol}",
            "",
            "| n_steps | Trotter不忠実度 | ノイズ不忠実度 | 合計不忠実度 | ノイズ/Trotter比 |",
            "|---------|----------------|--------------|------------|-----------------|",
        ])
        for r in ta["results"]:
            if r["p_depol"] == p_depol:
                ratio_str = f"{r['noise_to_trotter_ratio']:.1f}x" if r["noise_to_trotter_ratio"] < 1e10 else "∞"
                md_lines.append(
                    f"| {r['n_steps']} | {r['trotter_infidelity']:.8f} | "
                    f"{r['noise_infidelity']:.8f} | {r['total_infidelity']:.8f} | {ratio_str} |"
                )
        md_lines.append("")

    md_lines.extend([
        "**結論**: すべての実用的 p_depol 値でノイズ >> Trotter 誤差。",
        "最小 n_steps が常に最適。",
        "",
    ])

    # Test B summary
    md_lines.extend([
        "## Test B: ノイズバジェット分解（1ステップ、dt=1.0）",
        "",
        f"- p_depol: {tb['p_depol']}",
        "",
        "| 条件 | イベント数 | F(vs noiseless) | N_S1 |",
        "|------|----------|-----------------|------|",
        f"| ノイズなし | 0 | 1.00000000 | {tb['N_S1_noiseless']:.8f} |",
        f"| 全ノイズ | {tb['noise_events_total']} | "
        f"{tb['F_full_noise']:.8f} | {tb['N_S1_full_noise']:.8f} |",
        f"| Hゲートのみ | {tb['noise_events_hamiltonian']} | "
        f"{tb['F_hamiltonian_noise_only']:.8f} | {tb['N_S1_hamiltonian_only']:.8f} |",
        f"| Lチャネルのみ | {tb['noise_events_lindblad']} | "
        f"{tb['F_lindblad_noise_only']:.8f} | {tb['N_S1_lindblad_only']:.8f} |",
        f"| 乗算推定 | — | {tb['F_multiplicative_estimate']:.8f} | — |",
        "",
    ])

    # Test C summary
    md_lines.extend([
        "## Test C: ステップごとの ε(n) 追跡",
        "",
        f"- dt: {tc['dt']}, p_depol: {tc['p_depol']}",
        f"- 理論 ε (一定モデル): {tc['eps_theoretical']:.8f}",
        "",
        "### ステップごとの ε(n) — 収穫逓減効果の実証",
        "",
        "| Step | F(n) | ε(n) | F_ratio | 簡易予測 | Δ |",
        "|------|------|------|---------|---------|---|",
    ])
    for sd in tc["per_step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['F_n']:.8f} | {sd['eps_n']:.8f} | "
            f"{sd['F_ratio']:.8f} | {sd['simple_model_prediction']:.8f} | "
            f"{sd['delta_vs_simple']:+.6f} |"
        )

    eps_decrease_pct = (1.0 - tc['per_step_data'][-1]['eps_n'] / tc['per_step_data'][0]['eps_n']) * 100

    md_lines.extend([
        "",
        f"**ε(n) は step 1 → step {tc['n_steps']} で約{eps_decrease_pct:.0f}%減少。"
        f"F=(1-ε)^n モデルは不正確。**",
        "",
        "### dt 依存性検証（1ステップ測定）",
        "",
        "| dt | F_1step | ε |",
        "|----|---------|---|",
    ])
    for pd in tc["per_dt_data"]:
        md_lines.append(f"| {pd['dt']:.2f} | {pd['F_1step']:.8f} | {pd['eps']:.8f} |")
    md_lines.extend([
        "",
        f"**ε の dt 依存性: 絶対変動 {tc['dt_abs_variation']:.8f}、"
        f"相対変動 {tc['dt_rel_variation_pct']:.2f}%。ノイズは本質的にゲートノイズ。**",
        "",
    ])

    # Test D summary
    md_lines.extend([
        "## Test D: 改善版Qubit忠実度追跡（pair-only, p_depol=0.01）",
        "",
        "| Step | F_qudit | F_qubit_raw | F_qubit_norm | Tr_qubit | P_forbidden | "
        "Qudit N_S1 | Qubit N_S1_cond |",
        "|------|---------|-------------|--------------|----------|-------------|"
        "-----------|-----------------|",
    ])
    for sd in td["step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['F_qudit']:.6f} | {sd['F_qubit_raw']:.7f} | "
            f"{sd['F_qubit_norm']:.8f} | {sd['trace_qubit']:.6f} | "
            f"{sd['forbidden_pop']:.6f} | {sd['qudit_noisy_N_S1']:.7f} | "
            f"{sd['qubit_noisy_N_S1_conditional']:.7f} |"
        )
    md_lines.append("")

    # Test E summary (MODIFIED for iteration 24)
    md_lines.extend([
        "## Test E: 正確な GKSL Liouvillian 解との比較（iteration 25 Frobenius 解析）",
        "",
        "- 参照: ClassicalGKSLSimulator（exp(L·t) による厳密解）",
        f"- 正確解の N_S1(t=10): {te['exact_N_S1']:.10f}",
        f"- 正確解の純度: {te['exact_purity']:.10f}",
        f"- Lindblad チャネル数: {te['n_lindblad_channels']}",
        "- **Frobenius 解析（iteration 25 追加）**: Frobenius 距離 + 合成モデルフィッティング + 合成 Rate(T) 予測",
        "- **漸近解析（iteration 24）**: 拡張 dt 範囲 + 高次フィッティング + Richardson 外挿",
        "- **改善版分析**: 成分別多項式フィッティング + 正規化誤差係数",
        "- **修正された誤差分解**:",
        "  - **ST vs SCPT → 純粋な Stinespring 近似誤差**",
        "  - **T(SCPT,CT) → 回文順 Lie-Trotter 残差の直接測定**",
        "  - SCPT 基準の相殺比分析",
        "- **主要指標**: トレース距離 T(ρ,σ) = ½||ρ−σ||₁",
        "",
        "### 五方比較: ST vs Exact vs CT vs CPT vs SCPT",
        "",
        "| n_steps | dt | T(ST,ex) | T(CT,ex) | T(CPT,ex) | T(SCPT,ex) | T(ST,SCPT) | T(ST,CT) | T(SCPT,CT) | cancel_CT | cancel_SCPT |",
        "|---------|-----|---------|---------|----------|-----------|----------|---------|----------|----------|------------|",
    ])
    for cd in te["convergence_data"]:
        md_lines.append(
            f"| {cd['n_steps']} | {cd['dt']:.2f} | "
            f"{cd['T_ST_vs_exact']:.4e} | {cd['T_CT_vs_exact']:.4e} | "
            f"{cd['T_CPT_vs_exact']:.4e} | {cd['T_SCPT_vs_exact']:.4e} | "
            f"{cd['T_ST_vs_SCPT']:.4e} | {cd['T_ST_vs_CT']:.4e} | "
            f"{cd['T_SCPT_vs_CT']:.4e} | "
            f"{cd['cancellation_ratio']:.3f} | "
            f"{cd['cancellation_ratio_scpt']:.3f} |"
        )

    md_lines.extend([
        "",
        "### 誤差分解の解釈（iteration 24）",
        "",
        "| 比較ペア | 分離される誤差 | 期待 Rate(T) | 説明 |",
        "|----------|--------------|-------------|------|",
        "| ST vs Exact | 全誤差 | O(dt) 漸近 | Stinespring + 回文順残差 + Strang |",
        "| CT vs Exact | Strang H-D 分割 | O(dt²) | 厳密散逸子 + Strang 分割 |",
        "| CPT vs Exact | 順方向 Lie-Trotter + Strang | O(dt) | 個別散逸子の順方向積 + Strang |",
        "| SCPT vs Exact | 回文順残差 + Strang | O(dt²) | 個別散逸子の回文順積 + Strang |",
        "| **ST vs SCPT** | **純粋 Stinespring 近似誤差** | **O(dt)** | **同一順序での Stinespring vs 厳密チャネルの差** |",
        "| ST vs CT | Stinespring + 残差 | O(dt) | 小さい dt で ≈ Stinespring 誤差 |",
        "| CPT vs CT | 順方向 Lie-Trotter 積誤差 | O(dt) | 順方向積 vs 合計散逸子 |",
        "| **SCPT vs CT** | **回文順 Lie-Trotter 残差** | **O(dt²)** | **回文順積 vs 合計散逸子（無視可能なサイズ）** |",
        "| ST vs CPT | Stinespring + 順序差 | O(dt) | **注意**: 順序差を含むため分離不完全 |",
        "",
    ])

    # Convergence rates
    md_lines.extend([
        "### 収束次数（主要指標: トレース距離 T）",
        "",
    ])

    pair_labels_md = [
        ("st_ex", "ST vs Exact"),
        ("ct_ex", "Classical Trotter vs Exact"),
        ("cpt_ex", "Classical Product Trotter vs Exact"),
        ("scpt_ex", "Symmetric Classical Product Trotter vs Exact"),
        ("st_ct", "ST vs Classical Trotter"),
        ("st_cpt", "ST vs CPT（Stinespring + 順序差）"),
        ("st_scpt", "ST vs SCPT（純粋 Stinespring 誤差）"),
        ("cpt_ct", "CPT vs CT（順方向 Lie-Trotter 積誤差）"),
        ("scpt_ct", "SCPT vs CT（回文順 Lie-Trotter 残差）"),
    ]

    for pair_short, pair_label in pair_labels_md:
        md_lines.extend([
            f"#### {pair_label}",
            "",
            "| dt 範囲 | Rate(1-F) | Rate(T) | Rate(d_B) | clip |",
            "|---------|----------|---------|----------|------|",
        ])
        rates_infid = te["rates"].get(f"infid_{pair_short}", [])
        rates_tdist = te["rates"].get(f"tdist_{pair_short}", [])
        rates_bures = te["rates"].get(f"bures_{pair_short}", [])
        for k in range(len(rates_infid)):
            dt_prev = te["convergence_data"][k]["dt"]
            dt_curr = te["convergence_data"][k + 1]["dt"]
            clip_str = "`*`" if (te["convergence_data"][k].get("fidelity_clipped", False) or te["convergence_data"][k+1].get("fidelity_clipped", False)) else ""
            md_lines.append(
                f"| {dt_prev:.1f} → {dt_curr:.1f} | "
                f"{rates_infid[k]:.2f} | "
                f"{rates_tdist[k]:.2f} | "
                f"{rates_bures[k]:.2f} | {clip_str} |"
            )
        avg_infid = te.get(f"avg_rate_infid_{pair_short}", float("nan"))
        avg_tdist = te.get(f"avg_rate_tdist_{pair_short}", float("nan"))
        avg_bures = te.get(f"avg_rate_bures_{pair_short}", float("nan"))
        md_lines.extend([
            "",
            f"**平均 Rate(1-F): {avg_infid:.2f}**",
            f"**平均 Rate(T): {avg_tdist:.2f}** (← 主要指標)",
            f"**平均 Rate(d_B): {avg_bures:.2f}**",
            "",
        ])

    md_lines.extend([
        "### 混合状態における 1-F と T² の関係",
        "",
        "Fuchs–van de Graaf 不等式: T² ≤ 1-F ≤ 2T−T²。",
        "**1-F = T² は純粋状態でのみ成立。混合状態では (1-F)/T² >> 1 が想定される。**",
        "Bures距離 d_B はこの中間的な指標で、d_B² = 2(1−√F) の関係を持つ。",
        "",
    ])

    # Step detail
    md_lines.extend([
        "### ステップごとの精度（n_steps=10, dt=1.0）",
        "",
        "| Step | 1-F | T_dist | d_B | d_B/T | P(ST) | P(exact) | N_S1(ST) | N_S1(exact) | clip |",
        "|------|-----|--------|-----|-------|-------|----------|----------|-------------|------|",
    ])
    for sd in te["step_detail_n10"]:
        ratio_str = f"{sd['ratio_dB_over_T']:.4f}" if not np.isnan(sd.get("ratio_dB_over_T", float("nan"))) else "—"
        clip_str = "`*`" if sd.get("fidelity_clipped", False) else ""
        md_lines.append(
            f"| {sd['step']} | {sd['infidelity']:.4e} | {sd['trace_dist_vs_exact']:.4e} | "
            f"{sd['bures_dist']:.4e} | {ratio_str} | "
            f"{sd['purity_ST']:.6f} | {sd['purity_exact']:.6f} | "
            f"{sd['N_S1_ST']:.8f} | {sd['N_S1_exact']:.8f} | {clip_str} |"
        )

    # N_S1 convergence (with signed errors and SCPT)
    md_lines.extend([
        "",
        "### N_S1 ポピュレーション収束（符号付き誤差）",
        "",
        "| n_steps | dt | Δ_ST | Δ_CT | Δ_SCPT | Δ_CPT | 相殺比 |",
        "|---------|-----|------|------|--------|-------|--------|",
    ])
    for cd in te["convergence_data"]:
        md_lines.append(
            f"| {cd['n_steps']} | {cd['dt']:.2f} | "
            f"{cd['delta_N_S1_ST']:+.4e} | {cd['delta_N_S1_CT']:+.4e} | "
            f"{cd['delta_N_S1_SCPT']:+.4e} | {cd['delta_N_S1_CPT']:+.4e} | "
            f"{cd['cancellation_ratio']:.3f} |"
        )

    md_lines.append("")

    # Polynomial fitting — aggregate (reference)
    if te.get("polyfit_aggregate"):
        pf = te["polyfit_aggregate"]
        md_lines.extend([
            "### 多項式フィッティング（集約 — 参考値）: T(ST,exact) ≈ a·dt + b·dt²",
            "",
            f"- a = {pf['a_coefficient']:.6e}",
            f"- b = {pf['b_coefficient']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.4f}' for e in pf['relative_fit_errors'])}",
            "- ⚠ 誤差相殺により最大 29% の誤差。成分別フィッティングを参照。",
            "",
        ])

    # Polynomial fitting — component-wise (accurate)
    if te.get("polyfit_component"):
        cf = te["polyfit_component"]
        md_lines.extend([
            "### 成分別多項式フィッティング（正確な誤差分解）",
            "",
            f"#### Stinespring 誤差: T(ST,SCPT) ≈ a_stine · dt",
            "",
            f"- a_stine = {cf['a_stine_coefficient']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.4f}' for e in cf['a_stine_relative_errors'])}",
            "",
            f"#### Strang+回文順残差: T(SCPT,exact) ≈ b_strang · dt²",
            "",
            f"- b_strang = {cf['b_strang_coefficient']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.4f}' for e in cf['b_strang_relative_errors'])}",
            "",
            f"#### 回文順 Lie-Trotter 残差: T(SCPT,CT) ≈ c_palindromic · dt²",
            "",
            f"- c_palindromic = {cf['c_palindromic_coefficient']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.4f}' for e in cf['c_palindromic_relative_errors'])}",
            "",
        ])

    # Normalized error coefficients table
    md_lines.extend([
        "### 正規化誤差係数（T/dt^n）",
        "",
        "| dt | T(ST,SCPT)/dt | T(SCPT,ex)/dt² | T(CT,ex)/dt² | T(SCPT,CT)/dt² | cancel_CT | cancel_SCPT |",
        "|----|-------------|---------------|-------------|---------------|----------|------------|",
    ])
    for cd_entry in te["convergence_data"]:
        md_lines.append(
            f"| {cd_entry['dt']:.2f} | "
            f"{cd_entry['T_ST_vs_SCPT_over_dt']:.6e} | "
            f"{cd_entry['T_SCPT_vs_exact_over_dt2']:.6e} | "
            f"{cd_entry['T_CT_vs_exact_over_dt2']:.6e} | "
            f"{cd_entry['T_SCPT_vs_CT_over_dt2']:.6e} | "
            f"{cd_entry['cancellation_ratio']:.4f} | "
            f"{cd_entry['cancellation_ratio_scpt']:.4f} |"
        )
    md_lines.append("")

    # Higher-order fitting (new in iteration 24)
    if te.get("polyfit_higher_order"):
        hf = te["polyfit_higher_order"]
        md_lines.extend([
            "### 高次フィッティング T = c₀·dt^n + c₁·dt^(n+1)（iteration 24 追加）",
            "",
            "#### Stinespring 誤差: T(ST,SCPT) = (a₀ + a₁·dt) · dt",
            "",
            f"- a₀（漸近値）= {hf['stinespring_a0']:.6e}",
            f"- a₁（補正項）= {hf['stinespring_a1']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.6f}' for e in hf['stinespring_ho_relative_errors'])}",
            "",
            "#### Strang+回文順残差: T(SCPT,ex) = (b₀ + b₁·dt) · dt²",
            "",
            f"- b₀（漸近値）= {hf['strang_b0']:.6e}",
            f"- b₁（補正項）= {hf['strang_b1']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.6f}' for e in hf['strang_ho_relative_errors'])}",
            "",
            "#### 回文順残差: T(SCPT,CT) = (c₀ + c₁·dt) · dt²",
            "",
            f"- c₀（漸近値）= {hf['palindromic_c0']:.6e}",
            f"- c₁（補正項）= {hf['palindromic_c1']:.6e}",
            f"- 相対フィット誤差: {', '.join(f'{e:.6f}' for e in hf['palindromic_ho_relative_errors'])}",
            "",
        ])

    # Error vector angle estimation (new in iteration 24)
    if te.get("error_angle_data"):
        md_lines.extend([
            "### 誤差ベクトル角度推定（iteration 24 追加）",
            "",
            "Stinespring 誤差と Strang 誤差のなす角度を余弦定理（トレース距離による近似）で推定。",
            "cos(θ) ≈ (T²(ST,ex) - T²(ST,SCPT) - T²(SCPT,ex)) / (2·T(ST,SCPT)·T(SCPT,ex))",
            "",
            "| dt | cos(θ) | θ (deg) | Stine/Strang 比 |",
            "|----|--------|---------|----------------|",
        ])
        for ad in te["error_angle_data"]:
            md_lines.append(
                f"| {ad['dt']:.2f} | {ad['cos_theta_approx']:.4f} | "
                f"{ad['theta_deg_approx']:.1f} | "
                f"{ad['stinespring_over_strang_ratio']:.2f} |"
            )
        md_lines.append("")

    # Richardson extrapolation (new in iteration 24)
    if te.get("richardson_extrapolation"):
        re_data = te["richardson_extrapolation"]
        md_lines.extend([
            "### Richardson 外挿（iteration 24 追加）",
            "",
            f"- 使用 dt ペア: {re_data['dt_pair'][0]:.2f}, {re_data['dt_pair'][1]:.2f}",
            f"- a₀_Richardson（Stinespring）= {re_data['a0_stinespring_richardson']:.6e}",
            f"- b₀_Richardson（Strang）= {re_data['b0_strang_richardson']:.6e}",
            f"- c₀_Richardson（回文順残差）= {re_data['c0_palindromic_richardson']:.6e}",
            "",
        ])

    # Composite model fitting (new in iteration 25)
    if te.get("composite_fit"):
        cf_data = te["composite_fit"]
        md_lines.extend([
            "### 合成モデルフィッティング（iteration 25 追加）",
            "",
            f"T(ST,ex) = (p₀ + p₁·dt) · dt",
            "",
            f"- p₀（実効先導係数）= {cf_data['p0_leading']:.6e}",
            f"- p₁（相殺効果）= {cf_data['p1_cancellation']:.6e}"
            f" {'(< 0 → 相殺確認)' if cf_data['p1_cancellation'] < 0 else ''}",
            f"- 相対フィット誤差: {', '.join(f'{e:.6f}' for e in cf_data['relative_fit_errors'])}",
            "",
        ])

    # Frobenius angle estimation (new in iteration 25)
    if te.get("frobenius_angle_data"):
        md_lines.extend([
            "### Frobenius ベース誤差ベクトル角度推定（iteration 25 追加）",
            "",
            "Frobenius ノルム（L₂）は内積空間ノルムであるため、余弦定理が厳密に成立する。",
            "cos(θ_F) = (d_F²(ST,ex) - d_F²(ST,SCPT) - d_F²(SCPT,ex)) / (2·d_F(ST,SCPT)·d_F(SCPT,ex))",
            "",
            "| dt | cos(θ_F) | θ_F (deg) | cos(θ_T) | θ_T (deg) | Δθ |",
            "|----|----------|-----------|----------|-----------|-----|",
        ])
        for ad in te["frobenius_angle_data"]:
            delta_str = f"{ad['delta_theta_deg']:+.1f}" if not np.isnan(ad.get("delta_theta_deg", float("nan"))) else "—"
            md_lines.append(
                f"| {ad['dt']:.2f} | {ad['cos_theta_frobenius']:.4f} | "
                f"{ad['theta_frobenius_deg']:.1f} | "
                f"{ad['cos_theta_trace']:.4f} | "
                f"{ad['theta_trace_deg']:.1f} | "
                f"{delta_str} |"
            )
        md_lines.append("")

    # Rate(T) prediction (new in iteration 25)
    if te.get("rate_prediction") and te["rate_prediction"].get("predictions"):
        rp = te["rate_prediction"]
        md_lines.extend([
            "### 合成 Rate(T) 予測と実測の比較（iteration 25 追加）",
            "",
            f"合成モデル: T(ST,ex) ≈ {rp['model_p0']:.6e}·dt + {rp['model_p1']:.6e}·dt²",
            "",
            "| dt 範囲 | Rate 予測 | Rate 実測 | 差 |",
            "|---------|----------|----------|-----|",
        ])
        for pred in rp["predictions"]:
            delta_str = f"{pred['delta_rate']:+.4f}" if not np.isnan(pred.get("delta_rate", float("nan"))) else "—"
            md_lines.append(
                f"| {pred['dt_prev']:.1f} → {pred['dt_curr']:.1f} | "
                f"{pred['rate_predicted']:.4f} | "
                f"{pred['rate_observed']:.4f} | "
                f"{delta_str} |"
            )
        md_lines.append("")
        if rp.get("dt_threshold_rate_095") is not None:
            md_lines.extend([
                f"- Rate(T) > 0.95 に必要な dt < {rp['dt_threshold_rate_095']:.4f} "
                f"(n_steps > {10.0/rp['dt_threshold_rate_095']:.0f})",
                f"- Rate(T) > 0.99 に必要な dt < {rp['dt_threshold_rate_099']:.4f} "
                f"(n_steps > {10.0/rp['dt_threshold_rate_099']:.0f})",
                "",
            ])

    md_path = output_dir / f"iteration25_frobenius_analysis_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY (primary convergence metric: trace distance)")
    print("NOTE: ST uses symmetric Lindblad product ordering")
    print("=" * 70)

    # Test A
    for p_depol in ta["p_depol_values"]:
        entries = [r for r in ta["results"] if r["p_depol"] == p_depol]
        if entries:
            best = min(entries, key=lambda r: r["total_infidelity"])
            worst_ratio = min(entries, key=lambda r: r["noise_to_trotter_ratio"])
            print(
                f"  Test A (p={p_depol}): Best n_steps={best['n_steps']} "
                f"(total infid={best['total_infidelity']:.6f}), "
                f"min noise/Trotter ratio={worst_ratio['noise_to_trotter_ratio']:.1f}x"
            )

    # Test B
    print(f"  Test B: Full F={tb['F_full_noise']:.6f}, "
          f"H-only F={tb['F_hamiltonian_noise_only']:.6f}, "
          f"L-only F={tb['F_lindblad_noise_only']:.6f}, "
          f"composed={tb['F_multiplicative_estimate']:.6f}")

    # Test C
    n_steps_c = tc["n_steps"]
    print(f"  Test C: eps(1)={tc['per_step_data'][0]['eps_n']:.6f}, "
          f"eps({n_steps_c})={tc['per_step_data'][-1]['eps_n']:.6f} "
          f"(decrease: {eps_decrease_pct:.1f}%), "
          f"dt rel. variation: {tc['dt_rel_variation_pct']:.2f}%")

    # Test D
    if td["step_data"]:
        final = td["step_data"][-1]
        print(f"  Test D: F_qd={final['F_qudit']:.6f}, "
              f"F_qb_raw={final['F_qubit_raw']:.7f}, "
              f"F_qb_norm={final['F_qubit_norm']:.8f}")

    # Test E
    print(f"  Test E: avg convergence rate ST vs Exact: "
          f"T(primary)={te.get('avg_rate_tdist_st_ex', float('nan')):.2f}, "
          f"1-F={te.get('avg_rate_infid_st_ex', float('nan')):.2f}")
    print(f"  Test E: avg convergence rate CT vs Exact: "
          f"T(primary)={te.get('avg_rate_tdist_ct_ex', float('nan')):.2f}")
    print(f"  Test E: avg convergence rate SCPT vs Exact: "
          f"T(primary)={te.get('avg_rate_tdist_scpt_ex', float('nan')):.2f}")
    print(f"  Test E: avg convergence rate CPT vs Exact: "
          f"T(primary)={te.get('avg_rate_tdist_cpt_ex', float('nan')):.2f}")
    print(f"  Test E: --- Error decomposition (corrected) ---")
    print(f"  Test E: avg rate ST vs SCPT (pure Stinespring error): "
          f"T={te.get('avg_rate_tdist_st_scpt', float('nan')):.2f}")
    print(f"  Test E: avg rate ST vs CT (Stinespring approx at small dt): "
          f"T={te.get('avg_rate_tdist_st_ct', float('nan')):.2f}")
    print(f"  Test E: avg rate CPT vs CT (forward Lie-Trotter error): "
          f"T={te.get('avg_rate_tdist_cpt_ct', float('nan')):.2f}")
    print(f"  Test E: avg rate SCPT vs CT (palindromic residual): "
          f"T={te.get('avg_rate_tdist_scpt_ct', float('nan')):.2f}")
    if te.get("polyfit_aggregate"):
        pf = te["polyfit_aggregate"]
        print(f"  Test E: aggregate fit T(ST,ex) ≈ {pf['a_coefficient']:.2e}·dt + {pf['b_coefficient']:.2e}·dt² (reference)")
    if te.get("polyfit_component"):
        cf = te["polyfit_component"]
        print(f"  Test E: component fit T(ST,SCPT) = {cf['a_stine_coefficient']:.2e}·dt (Stinespring)")
        print(f"  Test E: component fit T(SCPT,ex) = {cf['b_strang_coefficient']:.2e}·dt² (Strang)")
        print(f"  Test E: component fit T(SCPT,CT) = {cf['c_palindromic_coefficient']:.2e}·dt² (palindromic)")
    if te.get("polyfit_higher_order"):
        hf = te["polyfit_higher_order"]
        print(f"  Test E: higher-order Stinespring a0={hf['stinespring_a0']:.2e} (asymptotic)")
        print(f"  Test E: higher-order Strang b0={hf['strang_b0']:.2e} (asymptotic)")
        print(f"  Test E: higher-order Palindromic c0={hf['palindromic_c0']:.2e} (asymptotic)")
    if te.get("error_angle_data"):
        angles = te["error_angle_data"]
        cos_values = [a["cos_theta_approx"] for a in angles if not np.isnan(a["cos_theta_approx"])]
        if cos_values:
            print(f"  Test E: error vector angle range: {min(a['theta_deg_approx'] for a in angles if not np.isnan(a['theta_deg_approx'])):.1f}° – "
                  f"{max(a['theta_deg_approx'] for a in angles if not np.isnan(a['theta_deg_approx'])):.1f}°")
    if te.get("richardson_extrapolation"):
        re_data = te["richardson_extrapolation"]
        print(f"  Test E: Richardson a0={re_data['a0_stinespring_richardson']:.2e}, "
              f"b0={re_data['b0_strang_richardson']:.2e}, "
              f"c0={re_data['c0_palindromic_richardson']:.2e}")
    if te.get("composite_fit"):
        cf_data = te["composite_fit"]
        print(f"  Test E: composite fit T(ST,ex) = ({cf_data['p0_leading']:.2e} + {cf_data['p1_cancellation']:.2e}·dt)·dt")
    if te.get("frobenius_angle_data"):
        f_angles = te["frobenius_angle_data"]
        theta_f_vals_summary = [a["theta_frobenius_deg"] for a in f_angles if not np.isnan(a.get("theta_frobenius_deg", float("nan")))]
        if theta_f_vals_summary:
            print(f"  Test E: Frobenius angle range: {min(theta_f_vals_summary):.1f}° – {max(theta_f_vals_summary):.1f}° "
                  f"(variation: {max(theta_f_vals_summary)-min(theta_f_vals_summary):.1f}°)")
    if te.get("rate_prediction") and te["rate_prediction"].get("predictions"):
        rp = te["rate_prediction"]
        max_delta = max(abs(p["delta_rate"]) for p in rp["predictions"] if not np.isnan(p.get("delta_rate", float("nan"))))
        print(f"  Test E: Rate(T) prediction max |Δ| = {max_delta:.4f}")
        if rp.get("dt_threshold_rate_095") is not None:
            print(f"  Test E: Rate(T) > 0.95 requires dt < {rp['dt_threshold_rate_095']:.4f}")
    print(f"  Test E: exact purity={te['exact_purity']:.6f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
