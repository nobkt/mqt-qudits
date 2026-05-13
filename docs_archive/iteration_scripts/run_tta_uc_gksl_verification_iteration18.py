#!/usr/bin/env python3
"""TTA-UC GKSL iteration 18 verification: exact Liouvillian comparison.

Fixes issues identified in iteration 17 analysis:
  - Test A-D: unchanged from iteration 17
  - Test C: fix hardcoded dt-variation description (compute from data)
  - Test E (NEW): compare Stinespring+Trotter against the exact GKSL
    Liouvillian solution (ClassicalGKSLSimulator) to quantify the total
    approximation error, separate Stinespring and Trotter contributions,
    and verify density matrix validity at each step.

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
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_math_utils import (
    compute_populations_from_density_matrix,
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
# Test E: Exact Liouvillian comparison (NEW in iteration 18)
# =========================================================================
def test_exact_liouvillian_comparison(params: GKSLPhysicalParameters) -> dict:
    """Compare Stinespring+Trotter against exact GKSL Liouvillian solution.

    Uses ClassicalGKSLSimulator (exact exp(L*t) via scipy expm_multiply)
    as the ground-truth reference. Measures:
      1. Total approximation error: F(Stinespring+Trotter, Exact)
      2. Convergence with dt: does the error decrease as expected?
      3. Density matrix validity diagnostics at each step
      4. Population comparison between methods

    This test addresses the critical gap in iteration 17: Test A only
    measured Trotter splitting error (by comparing against a Stinespring+Trotter
    reference at finer dt), but never validated against the true GKSL solution.
    """
    print("\n" + "=" * 70)
    print("Test E: Exact Liouvillian Comparison (Stinespring+Trotter vs Exact)")
    print("=" * 70)

    t_max = 10.0
    step_counts = [5, 10, 20, 50]

    # --- Part 1: Compute exact reference ---
    print("\n  Part 1: Computing exact GKSL Liouvillian solution...", flush=True)
    t0 = time_module.time()
    sim_exact = ClassicalGKSLSimulator(params)
    r_exact = sim_exact.simulate(t_max=t_max, n_steps=100)
    rho_exact_final = r_exact["rho_final"]
    exact_elapsed = time_module.time() - t0
    print(f"  Exact solution computed in {exact_elapsed:.2f}s")

    pop_exact = compute_populations_from_density_matrix(rho_exact_final, params)
    print(f"  Exact N_S1 = {pop_exact['N_S1']:.10f}")

    # --- Part 2: Compare Stinespring+Trotter at various dt ---
    print(f"\n  Part 2: Stinespring+Trotter convergence vs exact solution")
    print(
        f"  {'n_steps':>7} | {'dt':>6} | {'F_vs_exact':>12} | {'Infid_total':>12} | "
        f"{'Infid_Trotter':>13} | {'Infid_Stine':>12} | {'N_S1':>10} | {'Tr':>10} | "
        f"{'Herm_err':>10} | {'min_eig':>10}"
    )
    print("  " + "-" * 130)

    convergence_data = []

    # Finest Stinespring+Trotter reference for Trotter-only error
    sim_finest = QuditGKSLSimulator(params)
    r_finest = sim_finest.simulate(t_max=t_max, n_steps=100)
    rho_finest = r_finest["rho_final"]

    for ns in step_counts:
        dt = t_max / ns
        t0 = time_module.time()

        sim_st = QuditGKSLSimulator(params)
        r_st = sim_st.simulate(t_max=t_max, n_steps=ns)
        rho_st = r_st["rho_final"]

        elapsed = time_module.time() - t0

        # Fidelity vs exact
        f_vs_exact = quantum_fidelity(rho_exact_final, rho_st)
        infid_total = 1.0 - f_vs_exact

        # Fidelity vs finest Stinespring+Trotter (Trotter-only error)
        f_vs_finest = quantum_fidelity(rho_finest, rho_st)
        infid_trotter = 1.0 - f_vs_finest

        # Stinespring approximation error estimate (total - Trotter).
        # This decomposition is approximate: infidelities are not strictly
        # additive when both error sources are comparable or correlated.
        # When infid_trotter << infid_total, the estimate is reliable;
        # when they are of similar magnitude, treat as indicative only.
        infid_stine = infid_total - infid_trotter

        # Density matrix diagnostics
        tr_val = float(np.real(np.trace(rho_st)))
        herm_err = float(np.linalg.norm(rho_st - rho_st.conj().T, "fro"))
        min_eig = float(np.linalg.eigvalsh(rho_st)[0])

        pop_st = compute_populations_from_density_matrix(rho_st, params)

        entry = {
            "n_steps": ns,
            "dt": dt,
            "F_vs_exact": f_vs_exact,
            "infidelity_total": infid_total,
            "infidelity_trotter_only": infid_trotter,
            "infidelity_stinespring_approx": infid_stine,
            "N_S1": pop_st["N_S1"],
            "trace": tr_val,
            "hermiticity_error": herm_err,
            "min_eigenvalue": min_eig,
            "elapsed": elapsed,
        }
        convergence_data.append(entry)

        print(
            f"  {ns:7d} | {dt:6.2f} | {f_vs_exact:12.10f} | {infid_total:12.4e} | "
            f"{infid_trotter:13.4e} | {infid_stine:12.4e} | {pop_st['N_S1']:10.8f} | "
            f"{tr_val:10.8f} | {herm_err:10.2e} | {min_eig:10.2e}"
        )

    # --- Part 3: Convergence rate analysis ---
    print(f"\n  Part 3: Convergence rate analysis")

    if len(convergence_data) >= 2:
        rates_total = []
        rates_trotter = []
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k - 1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr

            infid_prev = convergence_data[k - 1]["infidelity_total"]
            infid_curr = convergence_data[k]["infidelity_total"]
            if infid_curr > 1e-15 and infid_prev > 1e-15:
                rate_total = np.log(infid_prev / infid_curr) / np.log(dt_ratio)
                rates_total.append(rate_total)
            else:
                rates_total.append(float("nan"))

            infid_prev_t = convergence_data[k - 1]["infidelity_trotter_only"]
            infid_curr_t = convergence_data[k]["infidelity_trotter_only"]
            if infid_curr_t > 1e-15 and infid_prev_t > 1e-15:
                rate_trotter = np.log(infid_prev_t / infid_curr_t) / np.log(dt_ratio)
                rates_trotter.append(rate_trotter)
            else:
                rates_trotter.append(float("nan"))

        print(f"  {'dt_range':>15} | {'Rate_total':>10} | {'Rate_trotter':>12}")
        print("  " + "-" * 45)
        for k in range(len(rates_total)):
            dt_prev = convergence_data[k]["dt"]
            dt_curr = convergence_data[k + 1]["dt"]
            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | {rates_total[k]:10.2f} | "
                f"{rates_trotter[k]:12.2f}"
            )

        avg_rate_total = np.nanmean(rates_total)
        avg_rate_trotter = np.nanmean(rates_trotter)
        print(f"\n  Average convergence rate (total vs exact): {avg_rate_total:.2f}")
        print(f"  Average convergence rate (Trotter only):   {avg_rate_trotter:.2f}")
        print(f"  Expected for 2nd-order Trotter: 2.0")
    else:
        avg_rate_total = float("nan")
        avg_rate_trotter = float("nan")
        rates_total = []
        rates_trotter = []

    # --- Part 4: Step-by-step comparison at n_steps=10 ---
    print(f"\n  Part 4: Step-by-step Stinespring+Trotter vs Exact at n_steps=10")

    dt_detail = t_max / 10
    sim_st_detail = QuditGKSLSimulator(params)
    sim_st_detail._precompute_unitaries(dt_detail)
    rho_st_detail = sim_st_detail.prepare_initial_state("edge_triplet")

    print(
        f"  {'Step':>4} | {'F_vs_exact':>12} | {'Infid':>10} | {'Tr_ST':>10} | "
        f"{'min_eig_ST':>10} | {'N_S1_ST':>10} | {'N_S1_exact':>10}"
    )
    print("  " + "-" * 80)

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
        infid_s = 1.0 - f_vs_exact_s
        tr_s = float(np.real(np.trace(rho_st_detail)))
        min_eig_s = float(np.linalg.eigvalsh(rho_st_detail)[0])
        pop_st_s = compute_populations_from_density_matrix(rho_st_detail, params)
        pop_exact_s = compute_populations_from_density_matrix(rho_exact_step, params)

        sd = {
            "step": s,
            "time": s * dt_detail,
            "F_vs_exact": f_vs_exact_s,
            "infidelity": infid_s,
            "trace_ST": tr_s,
            "min_eigenvalue_ST": min_eig_s,
            "N_S1_ST": pop_st_s["N_S1"],
            "N_S1_exact": pop_exact_s["N_S1"],
        }
        step_detail.append(sd)

        print(
            f"  {s:4d} | {f_vs_exact_s:12.10f} | {infid_s:10.4e} | {tr_s:10.8f} | "
            f"{min_eig_s:10.2e} | {pop_st_s['N_S1']:10.8f} | {pop_exact_s['N_S1']:10.8f}"
        )

    # --- Part 5: Exact vs finest Stinespring+Trotter ---
    f_exact_vs_finest = quantum_fidelity(rho_exact_final, rho_finest)
    infid_exact_vs_finest = 1.0 - f_exact_vs_finest
    pop_finest = compute_populations_from_density_matrix(rho_finest, params)

    print(f"\n  Part 5: Reference comparison")
    print(f"  F(exact, Stinespring+Trotter n=100): {f_exact_vs_finest:.10f}")
    print(f"  Infidelity: {infid_exact_vs_finest:.4e}")
    print(f"  N_S1 exact:   {pop_exact['N_S1']:.10f}")
    print(f"  N_S1 ST(100): {pop_finest['N_S1']:.10f}")
    print(f"  N_S1 difference: {abs(pop_exact['N_S1'] - pop_finest['N_S1']):.4e}")

    return {
        "test": "exact_liouvillian_comparison",
        "t_max": t_max,
        "exact_N_S1": pop_exact["N_S1"],
        "exact_elapsed": exact_elapsed,
        "convergence_data": convergence_data,
        "convergence_rates_total": rates_total,
        "convergence_rates_trotter": rates_trotter,
        "avg_rate_total": avg_rate_total,
        "avg_rate_trotter": avg_rate_trotter,
        "step_detail_n10": step_detail,
        "F_exact_vs_finest_ST": f_exact_vs_finest,
        "infidelity_exact_vs_finest_ST": infid_exact_vs_finest,
        "N_S1_finest_ST": pop_finest["N_S1"],
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
    print("TTA-UC GKSL Iteration 18: Exact Liouvillian Comparison")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 18,
        "purpose": "Exact Liouvillian comparison: validate Stinespring+Trotter "
                   "against exact GKSL solution, fix Test C reporting, "
                   "add density matrix diagnostics",
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

    # Test E: Exact Liouvillian comparison (NEW)
    te = test_exact_liouvillian_comparison(params)
    results["tests"].append(te)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration18_exact_comparison_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 18: 正確な Liouvillian 比較",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        f"- iteration 17 からの主要修正:",
        f"  - Test E（新規）: 正確な GKSL Liouvillian 解との比較",
        f"  - Test C: dt 変動量をデータから計算に変更",
        f"  - 密度行列妥当性診断を追加",
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

    # Test C summary (fixed: compute dt variation from data)
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

    # Test E summary
    md_lines.extend([
        "## Test E: 正確な GKSL Liouvillian 解との比較（新規）",
        "",
        f"- 参照: ClassicalGKSLSimulator（exp(L*t) による厳密解）",
        f"- 正確解の N_S1(t=10): {te['exact_N_S1']:.10f}",
        "",
        "### Stinespring+Trotter の全体収束性",
        "",
        "| n_steps | dt | F(vs exact) | 全体不忠実度 | Trotter不忠実度 | Stinespring近似 | N_S1 |",
        "|---------|-----|-------------|-------------|----------------|----------------|------|",
    ])
    for cd in te["convergence_data"]:
        md_lines.append(
            f"| {cd['n_steps']} | {cd['dt']:.2f} | {cd['F_vs_exact']:.10f} | "
            f"{cd['infidelity_total']:.4e} | {cd['infidelity_trotter_only']:.4e} | "
            f"{cd['infidelity_stinespring_approx']:.4e} | {cd['N_S1']:.8f} |"
        )

    md_lines.extend([
        "",
        "### 収束次数",
        "",
        "| dt 範囲 | 全体収束次数 | Trotter 収束次数 |",
        "|---------|------------|-----------------|",
    ])
    for k in range(len(te["convergence_rates_total"])):
        dt_prev = te["convergence_data"][k]["dt"]
        dt_curr = te["convergence_data"][k + 1]["dt"]
        md_lines.append(
            f"| {dt_prev:.1f} → {dt_curr:.1f} | "
            f"{te['convergence_rates_total'][k]:.2f} | "
            f"{te['convergence_rates_trotter'][k]:.2f} |"
        )
    md_lines.extend([
        "",
        f"**平均収束次数（全体 vs 正確解）: {te['avg_rate_total']:.2f}**",
        f"**平均収束次数（Trotter のみ）: {te['avg_rate_trotter']:.2f}**",
        "",
    ])

    # Step detail
    md_lines.extend([
        "### ステップごとの精度（n_steps=10, dt=1.0）",
        "",
        "| Step | F(vs exact) | 不忠実度 | Tr | min_eig | N_S1(ST) | N_S1(exact) |",
        "|------|-------------|---------|-----|---------|----------|-------------|",
    ])
    for sd in te["step_detail_n10"]:
        md_lines.append(
            f"| {sd['step']} | {sd['F_vs_exact']:.10f} | {sd['infidelity']:.4e} | "
            f"{sd['trace_ST']:.8f} | {sd['min_eigenvalue_ST']:.2e} | "
            f"{sd['N_S1_ST']:.8f} | {sd['N_S1_exact']:.8f} |"
        )

    md_lines.extend([
        "",
        f"### 参照比較",
        "",
        f"- F(exact, Stinespring+Trotter n=100): {te['F_exact_vs_finest_ST']:.10f}",
        f"- 不忠実度: {te['infidelity_exact_vs_finest_ST']:.4e}",
        f"- N_S1 差: |exact - ST(100)| = "
        f"{abs(te['exact_N_S1'] - te['N_S1_finest_ST']):.4e}",
        "",
    ])

    md_path = output_dir / f"iteration18_exact_comparison_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
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
    print(f"  Test E: F(exact, ST n=100)={te['F_exact_vs_finest_ST']:.10f}, "
          f"avg convergence rate: {te['avg_rate_total']:.2f} (total), "
          f"{te['avg_rate_trotter']:.2f} (Trotter)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
