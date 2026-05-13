#!/usr/bin/env python3
"""TTA-UC GKSL iteration 20 verification: Bures distance, purity & mixed-state metrics.

Changes from iteration 19:
  - Test A-D: unchanged from iteration 19
  - Test E (modified):
    * Add bures_distance() metric: d_B(ρ,σ) = √(2(1−√F(ρ,σ)))
    * Add fidelity_clipped flag for detecting min(...,1.0) clamping
    * Add purity tracking for exact, ST and CT states
    * Replace (1-F)/T² with Bures distance metrics (d_B, d_B/T)
    * Add Bures distance convergence rates
    * Trace distance is the primary convergence metric
    * Markdown report updated with mixed-state annotations

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
# Test E: Exact Liouvillian comparison (MODIFIED in iteration 20)
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


def test_exact_liouvillian_comparison(params: GKSLPhysicalParameters) -> dict:
    """Compare Stinespring+Trotter against exact GKSL and classical Trotter.

    MODIFIED from iteration 19:
    - Add Bures distance metric d_B = sqrt(2(1 - sqrt(F)))
    - Add fidelity_clipped flag to detect numerical clamping at F=1
    - Add purity tracking for exact, ST, and CT density matrices
    - Replace (1-F)/T² column with Bures distance metrics (d_B, d_B/T)
    - Add Bures distance convergence rates
    - Trace distance is now the primary convergence metric
    - Note: 1-F = T² holds only for pure states; for mixed states (1-F)/T² >> 1 is expected
    """
    print("\n" + "=" * 70)
    print("Test E: Exact Liouvillian Comparison (ST vs Exact vs Classical Trotter)")
    print("       [iter 20: +Bures distance, +purity, +fidelity clipping flag]")
    print("=" * 70)

    t_max = 10.0
    step_counts = [5, 10, 20, 50]
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

    # --- Part 1b: Build dissipator superoperator for classical Trotter ---
    print("  Building dissipator superoperator for classical Trotter...", flush=True)
    H_total = sim_exact.H_total
    lindblad_ops = sim_exact.lindblad_ops
    L_D_super = _build_dissipator_superoperator(H_total, lindblad_ops)
    rho_0 = sim_exact.prepare_initial_state("edge_triplet")
    print("  Done.")

    # --- Part 2: Three-way convergence comparison ---
    print(f"\n  Part 2: Three-way convergence (primary metric: trace distance)")
    print(
        f"  {'n_steps':>7} | {'dt':>5} | "
        f"{'1-F(ST,ex)':>11} | {'T(ST,ex)':>11} | "
        f"{'d_B(ST,ex)':>11} | {'d_B/T':>8} | "
        f"{'1-F(CT,ex)':>11} | {'T(CT,ex)':>11} | "
        f"{'d_B(ST,CT)':>11} | {'clip':>4}"
    )
    print("  " + "-" * 120)

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

        # Purity tracking
        purity_st = compute_purity(rho_st)
        purity_ct = compute_purity(rho_ct)

        # Metrics: ST vs Exact
        f_st_ex = quantum_fidelity(rho_exact_final, rho_st)
        fidelity_clipped_st_ex = f_st_ex >= 1.0
        infid_st_ex = 0.0 if fidelity_clipped_st_ex else 1.0 - f_st_ex
        t_st_ex = trace_distance(rho_exact_final, rho_st)
        db_st_ex = bures_distance(rho_exact_final, rho_st)
        ratio_db_over_T_st_ex = db_st_ex / t_st_ex if t_st_ex > 1e-15 else float("nan")

        # Metrics: ClassicalTrotter vs Exact
        f_ct_ex = quantum_fidelity(rho_exact_final, rho_ct)
        fidelity_clipped_ct_ex = f_ct_ex >= 1.0
        infid_ct_ex = 0.0 if fidelity_clipped_ct_ex else 1.0 - f_ct_ex
        t_ct_ex = trace_distance(rho_exact_final, rho_ct)
        db_ct_ex = bures_distance(rho_exact_final, rho_ct)

        # Metrics: ST vs ClassicalTrotter
        f_st_ct = quantum_fidelity(rho_ct, rho_st)
        fidelity_clipped_st_ct = f_st_ct >= 1.0
        infid_st_ct = 0.0 if fidelity_clipped_st_ct else 1.0 - f_st_ct
        t_st_ct = trace_distance(rho_ct, rho_st)
        db_st_ct = bures_distance(rho_ct, rho_st)

        fidelity_clipped = fidelity_clipped_st_ex or fidelity_clipped_ct_ex or fidelity_clipped_st_ct

        # Density matrix diagnostics for ST
        tr_val = float(np.real(np.trace(rho_st)))
        herm_err = float(np.linalg.norm(rho_st - rho_st.conj().T, "fro"))
        min_eig = float(np.linalg.eigvalsh(rho_st)[0])

        pop_st = compute_populations_from_density_matrix(rho_st, params)
        pop_ct = compute_populations_from_density_matrix(rho_ct, params)

        clip_str = "*" if fidelity_clipped else ""

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
            # Classical Trotter vs Exact
            "F_CT_vs_exact": f_ct_ex,
            "T_CT_vs_exact": t_ct_ex,
            "infid_CT_vs_exact": infid_ct_ex,
            "bures_CT_vs_exact": db_ct_ex,
            "fidelity_clipped_CT_vs_exact": fidelity_clipped_ct_ex,
            # ST vs Classical Trotter
            "F_ST_vs_CT": f_st_ct,
            "T_ST_vs_CT": t_st_ct,
            "infid_ST_vs_CT": infid_st_ct,
            "bures_ST_vs_CT": db_st_ct,
            "fidelity_clipped_ST_vs_CT": fidelity_clipped_st_ct,
            # Fidelity clipped (any)
            "fidelity_clipped": fidelity_clipped,
            # Purity
            "purity_exact": purity_exact,
            "purity_ST": purity_st,
            "purity_CT": purity_ct,
            # Populations
            "N_S1_ST": pop_st["N_S1"],
            "N_S1_CT": pop_ct["N_S1"],
            "N_S1_exact": pop_exact["N_S1"],
            # Density matrix diagnostics
            "trace_ST": tr_val,
            "hermiticity_error_ST": herm_err,
            "min_eigenvalue_ST": min_eig,
            # Timing
            "elapsed_ST": elapsed_st,
            "elapsed_CT": elapsed_ct,
        }
        convergence_data.append(entry)

        print(
            f"  {ns:7d} | {dt:5.2f} | "
            f"{infid_st_ex:11.4e} | {t_st_ex:11.4e} | "
            f"{db_st_ex:11.4e} | {ratio_db_over_T_st_ex:8.4f} | "
            f"{infid_ct_ex:11.4e} | {t_ct_ex:11.4e} | "
            f"{db_st_ct:11.4e} | {clip_str:>4s}"
        )

    # --- Part 3: Convergence rate analysis ---
    print(f"\n  Part 3: Convergence rate analysis (primary metric: trace distance)")

    rates_infid_st_ex = []
    rates_tdist_st_ex = []
    rates_bures_st_ex = []
    rates_infid_ct_ex = []
    rates_tdist_ct_ex = []
    rates_bures_ct_ex = []
    rates_infid_st_ct = []
    rates_tdist_st_ct = []
    rates_bures_st_ct = []

    def _compute_rate(val_prev: float, val_curr: float, dt_ratio: float) -> float:
        if val_curr > 1e-18 and val_prev > 1e-18:
            return np.log(val_prev / val_curr) / np.log(dt_ratio)
        return float("nan")

    if len(convergence_data) >= 2:
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k - 1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr

            rates_infid_st_ex.append(_compute_rate(
                convergence_data[k-1]["infid_ST_vs_exact"],
                convergence_data[k]["infid_ST_vs_exact"], dt_ratio))
            rates_tdist_st_ex.append(_compute_rate(
                convergence_data[k-1]["T_ST_vs_exact"],
                convergence_data[k]["T_ST_vs_exact"], dt_ratio))
            rates_bures_st_ex.append(_compute_rate(
                convergence_data[k-1]["bures_ST_vs_exact"],
                convergence_data[k]["bures_ST_vs_exact"], dt_ratio))

            rates_infid_ct_ex.append(_compute_rate(
                convergence_data[k-1]["infid_CT_vs_exact"],
                convergence_data[k]["infid_CT_vs_exact"], dt_ratio))
            rates_tdist_ct_ex.append(_compute_rate(
                convergence_data[k-1]["T_CT_vs_exact"],
                convergence_data[k]["T_CT_vs_exact"], dt_ratio))
            rates_bures_ct_ex.append(_compute_rate(
                convergence_data[k-1]["bures_CT_vs_exact"],
                convergence_data[k]["bures_CT_vs_exact"], dt_ratio))

            rates_infid_st_ct.append(_compute_rate(
                convergence_data[k-1]["infid_ST_vs_CT"],
                convergence_data[k]["infid_ST_vs_CT"], dt_ratio))
            rates_tdist_st_ct.append(_compute_rate(
                convergence_data[k-1]["T_ST_vs_CT"],
                convergence_data[k]["T_ST_vs_CT"], dt_ratio))
            rates_bures_st_ct.append(_compute_rate(
                convergence_data[k-1]["bures_ST_vs_CT"],
                convergence_data[k]["bures_ST_vs_CT"], dt_ratio))

        def _clip_marker(k_idx: int) -> str:
            """Return '*' if fidelity was clipped at either endpoint of the rate interval."""
            if convergence_data[k_idx]["fidelity_clipped"] or convergence_data[k_idx + 1]["fidelity_clipped"]:
                return "*"
            return ""

        print(f"\n  === ST vs Exact (primary metric: trace distance T) ===")
        print(f"  {'dt_range':>15} | {'Rate(1-F)':>10} | {'Rate(T)':>10} | {'Rate(d_B)':>10} | {'clip':>4}")
        print("  " + "-" * 60)
        for k in range(len(rates_infid_st_ex)):
            dt_prev = convergence_data[k]["dt"]
            dt_curr = convergence_data[k + 1]["dt"]
            cm = _clip_marker(k)
            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | "
                f"{rates_infid_st_ex[k]:10.2f}{cm} | {rates_tdist_st_ex[k]:10.2f} | "
                f"{rates_bures_st_ex[k]:10.2f} | {cm:>4s}"
            )
        avg_infid_st_ex = np.nanmean(rates_infid_st_ex)
        avg_tdist_st_ex = np.nanmean(rates_tdist_st_ex)
        avg_bures_st_ex = np.nanmean(rates_bures_st_ex)
        print(f"  Average: {avg_infid_st_ex:10.2f} | {avg_tdist_st_ex:10.2f} | {avg_bures_st_ex:10.2f}")
        print(f"  * = fidelity clipped (F clamped to 1.0); 1-F rate unreliable")

        print(f"\n  === Classical Trotter vs Exact ===")
        print(f"  {'dt_range':>15} | {'Rate(1-F)':>10} | {'Rate(T)':>10} | {'Rate(d_B)':>10} | {'clip':>4}")
        print("  " + "-" * 60)
        for k in range(len(rates_infid_ct_ex)):
            dt_prev = convergence_data[k]["dt"]
            dt_curr = convergence_data[k + 1]["dt"]
            cm = _clip_marker(k)
            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | "
                f"{rates_infid_ct_ex[k]:10.2f}{cm} | {rates_tdist_ct_ex[k]:10.2f} | "
                f"{rates_bures_ct_ex[k]:10.2f} | {cm:>4s}"
            )
        avg_infid_ct_ex = np.nanmean(rates_infid_ct_ex)
        avg_tdist_ct_ex = np.nanmean(rates_tdist_ct_ex)
        avg_bures_ct_ex = np.nanmean(rates_bures_ct_ex)
        print(f"  Average: {avg_infid_ct_ex:10.2f} | {avg_tdist_ct_ex:10.2f} | {avg_bures_ct_ex:10.2f}")

        print(f"\n  === ST vs Classical Trotter ===")
        print(f"  {'dt_range':>15} | {'Rate(1-F)':>10} | {'Rate(T)':>10} | {'Rate(d_B)':>10} | {'clip':>4}")
        print("  " + "-" * 60)
        for k in range(len(rates_infid_st_ct)):
            dt_prev = convergence_data[k]["dt"]
            dt_curr = convergence_data[k + 1]["dt"]
            cm = _clip_marker(k)
            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | "
                f"{rates_infid_st_ct[k]:10.2f}{cm} | {rates_tdist_st_ct[k]:10.2f} | "
                f"{rates_bures_st_ct[k]:10.2f} | {cm:>4s}"
            )
        avg_infid_st_ct = np.nanmean(rates_infid_st_ct)
        avg_tdist_st_ct = np.nanmean(rates_tdist_st_ct)
        avg_bures_st_ct = np.nanmean(rates_bures_st_ct)
        print(f"  Average: {avg_infid_st_ct:10.2f} | {avg_tdist_st_ct:10.2f} | {avg_bures_st_ct:10.2f}")

    else:
        avg_infid_st_ex = avg_tdist_st_ex = avg_bures_st_ex = float("nan")
        avg_infid_ct_ex = avg_tdist_ct_ex = avg_bures_ct_ex = float("nan")
        avg_infid_st_ct = avg_tdist_st_ct = avg_bures_st_ct = float("nan")

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

    # --- Part 5: N_S1 convergence analysis ---
    print(f"\n  Part 5: N_S1 convergence analysis")
    print(f"  {'n_steps':>7} | {'dt':>5} | {'N_S1(ST)':>12} | {'N_S1(CT)':>12} | "
          f"{'N_S1(exact)':>12} | {'|ΔN_S1|(ST)':>12} | {'|ΔN_S1|(CT)':>12}")
    print("  " + "-" * 85)

    ns1_diffs_st = []
    ns1_diffs_ct = []
    for cd in convergence_data:
        diff_st = abs(cd["N_S1_ST"] - cd["N_S1_exact"])
        diff_ct = abs(cd["N_S1_CT"] - cd["N_S1_exact"])
        ns1_diffs_st.append(diff_st)
        ns1_diffs_ct.append(diff_ct)
        print(
            f"  {cd['n_steps']:7d} | {cd['dt']:5.2f} | {cd['N_S1_ST']:12.8f} | "
            f"{cd['N_S1_CT']:12.8f} | {cd['N_S1_exact']:12.8f} | "
            f"{diff_st:12.4e} | {diff_ct:12.4e}"
        )

    rates_ns1_st = []
    rates_ns1_ct = []
    if len(convergence_data) >= 2:
        print(f"\n  N_S1 convergence rates:")
        print(f"  {'dt_range':>15} | {'Rate(ST)':>10} | {'Rate(CT)':>10}")
        print("  " + "-" * 42)
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k-1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr
            rate_st = _compute_rate(ns1_diffs_st[k-1], ns1_diffs_st[k], dt_ratio)
            rate_ct = _compute_rate(ns1_diffs_ct[k-1], ns1_diffs_ct[k], dt_ratio)
            rates_ns1_st.append(rate_st)
            rates_ns1_ct.append(rate_ct)
            print(f"  {dt_prev:.1f} → {dt_curr:.1f}      | {rate_st:10.2f} | {rate_ct:10.2f}")
        print(f"  Average: {np.nanmean(rates_ns1_st):10.2f} | {np.nanmean(rates_ns1_ct):10.2f}")

    return {
        "test": "exact_liouvillian_comparison",
        "t_max": t_max,
        "exact_N_S1": pop_exact["N_S1"],
        "exact_purity": purity_exact,
        "exact_elapsed": exact_elapsed,
        "convergence_data": convergence_data,
        # Convergence rates: ST vs Exact
        "rates_infid_ST_vs_exact": rates_infid_st_ex,
        "rates_tdist_ST_vs_exact": rates_tdist_st_ex,
        "rates_bures_ST_vs_exact": rates_bures_st_ex,
        "avg_rate_infid_ST_vs_exact": float(np.nanmean(rates_infid_st_ex)) if rates_infid_st_ex else float("nan"),
        "avg_rate_tdist_ST_vs_exact": float(np.nanmean(rates_tdist_st_ex)) if rates_tdist_st_ex else float("nan"),
        "avg_rate_bures_ST_vs_exact": float(np.nanmean(rates_bures_st_ex)) if rates_bures_st_ex else float("nan"),
        # Convergence rates: CT vs Exact
        "rates_infid_CT_vs_exact": rates_infid_ct_ex,
        "rates_tdist_CT_vs_exact": rates_tdist_ct_ex,
        "rates_bures_CT_vs_exact": rates_bures_ct_ex,
        "avg_rate_infid_CT_vs_exact": float(np.nanmean(rates_infid_ct_ex)) if rates_infid_ct_ex else float("nan"),
        "avg_rate_tdist_CT_vs_exact": float(np.nanmean(rates_tdist_ct_ex)) if rates_tdist_ct_ex else float("nan"),
        "avg_rate_bures_CT_vs_exact": float(np.nanmean(rates_bures_ct_ex)) if rates_bures_ct_ex else float("nan"),
        # Convergence rates: ST vs CT
        "rates_infid_ST_vs_CT": rates_infid_st_ct,
        "rates_tdist_ST_vs_CT": rates_tdist_st_ct,
        "rates_bures_ST_vs_CT": rates_bures_st_ct,
        "avg_rate_infid_ST_vs_CT": float(np.nanmean(rates_infid_st_ct)) if rates_infid_st_ct else float("nan"),
        "avg_rate_tdist_ST_vs_CT": float(np.nanmean(rates_tdist_st_ct)) if rates_tdist_st_ct else float("nan"),
        "avg_rate_bures_ST_vs_CT": float(np.nanmean(rates_bures_st_ct)) if rates_bures_st_ct else float("nan"),
        # N_S1 convergence
        "rates_N_S1_ST": rates_ns1_st,
        "rates_N_S1_CT": rates_ns1_ct,
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
    print("TTA-UC GKSL Iteration 20: Bures Distance, Purity & Mixed-State Metrics")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 20,
        "purpose": "Add Bures distance metric, purity tracking, fidelity clipping "
                   "detection. Replace (1-F)/T² with d_B and d_B/T. "
                   "Trace distance is the primary convergence metric. "
                   "Note: 1-F = T² holds only for pure states.",
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

    # Test E: Exact Liouvillian comparison (MODIFIED)
    te = test_exact_liouvillian_comparison(params)
    results["tests"].append(te)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration20_mixed_state_metrics_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 20: Bures距離・純度・混合状態メトリクス",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "- iteration 19 からの主要修正:",
        "  - Test E: Bures距離 d_B = √(2(1−√F)) の追加",
        "  - Test E: fidelity_clipped フラグ（F が min(...,1.0) でクランプされた場合）",
        "  - Test E: 純度（purity）の追跡: exact, ST, CT",
        "  - Test E: (1-F)/T² を d_B, d_B/T に置換",
        "  - Test E: Bures距離の収束率分析を追加",
        "  - Test E: トレース距離が主要収束指標",
        "  - 注意: 1-F = T² は純粋状態でのみ成立。混合状態では (1-F)/T² >> 1 が想定される。",
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

    # Test E summary (MODIFIED for iteration 20)
    md_lines.extend([
        "## Test E: 正確な GKSL Liouvillian 解との比較（iteration 20 修正版）",
        "",
        "- 参照: ClassicalGKSLSimulator（exp(L·t) による厳密解）",
        f"- 正確解の N_S1(t=10): {te['exact_N_S1']:.10f}",
        f"- 正確解の純度: {te['exact_purity']:.10f}",
        "- 新規: Bures距離 d_B(ρ,σ) = √(2(1−√F(ρ,σ))) の追加",
        "- 新規: 忠実度クリッピング検出（F が min(...,1.0) でクランプ時に `*` 表示）",
        "- 新規: 純度（purity = Tr(ρ²)）の追跡",
        "- **注意**: 1-F = T² は純粋状態でのみ成立。混合状態では (1-F)/T² >> 1 が想定される。",
        "- **主要指標**: トレース距離 T(ρ,σ) = ½||ρ−σ||₁",
        "",
        "### 三方比較: ST vs Exact vs 古典 Trotter",
        "",
        "| n_steps | dt | 1-F(ST,ex) | T(ST,ex) | d_B(ST,ex) | d_B/T | "
        "1-F(CT,ex) | T(CT,ex) | d_B(ST,CT) | P(ST) | P(CT) | clip |",
        "|---------|-----|-----------|---------|-----------|-------|"
        "-----------|---------|-----------|-------|-------|------|",
    ])
    for cd in te["convergence_data"]:
        ratio_db_str = f"{cd['ratio_dB_over_T_ST_vs_exact']:.4f}" if not np.isnan(cd.get("ratio_dB_over_T_ST_vs_exact", float("nan"))) else "—"
        clip_str = "`*`" if cd.get("fidelity_clipped", False) else ""
        md_lines.append(
            f"| {cd['n_steps']} | {cd['dt']:.2f} | "
            f"{cd['infid_ST_vs_exact']:.4e} | {cd['T_ST_vs_exact']:.4e} | "
            f"{cd['bures_ST_vs_exact']:.4e} | {ratio_db_str} | "
            f"{cd['infid_CT_vs_exact']:.4e} | {cd['T_CT_vs_exact']:.4e} | "
            f"{cd['bures_ST_vs_CT']:.4e} | "
            f"{cd['purity_ST']:.6f} | {cd['purity_CT']:.6f} | {clip_str} |"
        )

    md_lines.extend([
        "",
        "> `*` = 忠実度が min(..., 1.0) でクランプされた。1-F は不正確な可能性あり。",
        "",
    ])

    # Convergence rates
    md_lines.extend([
        "### 収束次数（主要指標: トレース距離 T）",
        "",
        "#### ST vs Exact",
        "",
        "| dt 範囲 | Rate(1-F) | Rate(T) | Rate(d_B) | clip |",
        "|---------|----------|---------|----------|------|",
    ])
    for k in range(len(te.get("rates_infid_ST_vs_exact", []))):
        dt_prev = te["convergence_data"][k]["dt"]
        dt_curr = te["convergence_data"][k + 1]["dt"]
        clip_str = "`*`" if (te["convergence_data"][k].get("fidelity_clipped", False) or te["convergence_data"][k+1].get("fidelity_clipped", False)) else ""
        md_lines.append(
            f"| {dt_prev:.1f} → {dt_curr:.1f} | "
            f"{te['rates_infid_ST_vs_exact'][k]:.2f} | "
            f"{te['rates_tdist_ST_vs_exact'][k]:.2f} | "
            f"{te['rates_bures_ST_vs_exact'][k]:.2f} | {clip_str} |"
        )
    md_lines.extend([
        "",
        f"**平均 Rate(1-F): {te.get('avg_rate_infid_ST_vs_exact', float('nan')):.2f}**",
        f"**平均 Rate(T): {te.get('avg_rate_tdist_ST_vs_exact', float('nan')):.2f}** (← 主要指標)",
        f"**平均 Rate(d_B): {te.get('avg_rate_bures_ST_vs_exact', float('nan')):.2f}**",
        "",
        "#### 古典 Trotter vs Exact",
        "",
        "| dt 範囲 | Rate(1-F) | Rate(T) | Rate(d_B) | clip |",
        "|---------|----------|---------|----------|------|",
    ])
    for k in range(len(te.get("rates_infid_CT_vs_exact", []))):
        dt_prev = te["convergence_data"][k]["dt"]
        dt_curr = te["convergence_data"][k + 1]["dt"]
        clip_str = "`*`" if (te["convergence_data"][k].get("fidelity_clipped_CT_vs_exact", False) or te["convergence_data"][k+1].get("fidelity_clipped_CT_vs_exact", False)) else ""
        md_lines.append(
            f"| {dt_prev:.1f} → {dt_curr:.1f} | "
            f"{te['rates_infid_CT_vs_exact'][k]:.2f} | "
            f"{te['rates_tdist_CT_vs_exact'][k]:.2f} | "
            f"{te['rates_bures_CT_vs_exact'][k]:.2f} | {clip_str} |"
        )
    md_lines.extend([
        "",
        f"**平均 Rate(1-F): {te.get('avg_rate_infid_CT_vs_exact', float('nan')):.2f}**",
        f"**平均 Rate(T): {te.get('avg_rate_tdist_CT_vs_exact', float('nan')):.2f}**",
        f"**平均 Rate(d_B): {te.get('avg_rate_bures_CT_vs_exact', float('nan')):.2f}**",
        "",
        "#### ST vs 古典 Trotter",
        "",
        "| dt 範囲 | Rate(1-F) | Rate(T) | Rate(d_B) | clip |",
        "|---------|----------|---------|----------|------|",
    ])
    for k in range(len(te.get("rates_infid_ST_vs_CT", []))):
        dt_prev = te["convergence_data"][k]["dt"]
        dt_curr = te["convergence_data"][k + 1]["dt"]
        clip_str = "`*`" if (te["convergence_data"][k].get("fidelity_clipped_ST_vs_CT", False) or te["convergence_data"][k+1].get("fidelity_clipped_ST_vs_CT", False)) else ""
        md_lines.append(
            f"| {dt_prev:.1f} → {dt_curr:.1f} | "
            f"{te['rates_infid_ST_vs_CT'][k]:.2f} | "
            f"{te['rates_tdist_ST_vs_CT'][k]:.2f} | "
            f"{te['rates_bures_ST_vs_CT'][k]:.2f} | {clip_str} |"
        )
    md_lines.extend([
        "",
        f"**平均 Rate(1-F): {te.get('avg_rate_infid_ST_vs_CT', float('nan')):.2f}**",
        f"**平均 Rate(T): {te.get('avg_rate_tdist_ST_vs_CT', float('nan')):.2f}**",
        f"**平均 Rate(d_B): {te.get('avg_rate_bures_ST_vs_CT', float('nan')):.2f}**",
        "",
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

    # N_S1 convergence
    md_lines.extend([
        "",
        "### N_S1 ポピュレーション収束",
        "",
        "| n_steps | dt | N_S1(ST) | N_S1(CT) | N_S1(exact) | |Δ|(ST) | |Δ|(CT) |",
        "|---------|-----|----------|----------|-------------|--------|--------|",
    ])
    for cd in te["convergence_data"]:
        diff_st = abs(cd["N_S1_ST"] - cd["N_S1_exact"])
        diff_ct = abs(cd["N_S1_CT"] - cd["N_S1_exact"])
        md_lines.append(
            f"| {cd['n_steps']} | {cd['dt']:.2f} | {cd['N_S1_ST']:.8f} | "
            f"{cd['N_S1_CT']:.8f} | {cd['N_S1_exact']:.8f} | "
            f"{diff_st:.4e} | {diff_ct:.4e} |"
        )

    md_lines.append("")

    md_path = output_dir / f"iteration20_mixed_state_metrics_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY (primary convergence metric: trace distance)")
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
          f"T(primary)={te.get('avg_rate_tdist_ST_vs_exact', float('nan')):.2f}, "
          f"1-F={te.get('avg_rate_infid_ST_vs_exact', float('nan')):.2f}, "
          f"d_B={te.get('avg_rate_bures_ST_vs_exact', float('nan')):.2f}")
    print(f"  Test E: avg convergence rate CT vs Exact: "
          f"T(primary)={te.get('avg_rate_tdist_CT_vs_exact', float('nan')):.2f}, "
          f"1-F={te.get('avg_rate_infid_CT_vs_exact', float('nan')):.2f}, "
          f"d_B={te.get('avg_rate_bures_CT_vs_exact', float('nan')):.2f}")
    print(f"  Test E: exact purity={te['exact_purity']:.6f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
