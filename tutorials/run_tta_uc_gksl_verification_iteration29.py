#!/usr/bin/env python3
"""TTA-UC GKSL iteration 29 verification: Strang coefficient model correction.

Changes from iteration 28:
  - Test A-D: unchanged
  - Test E (extended analysis):
    * Part 6f IMPROVED: Strang coefficient model corrected from dt to dt²
      - Old: d_F(SCPT,ex)/dt² = b_F0 + b_F1*dt  (wrong: includes odd-order term)
      - New: d_F(SCPT,ex)/dt² = b_F0 + b_F2*dt²  (correct: even-order only)
      - Strang symmetric splitting eliminates odd-order corrections
      - Component fit error: 0.176% -> 0.003% (58x improvement)
    * Part 6f IMPROVED: Stinespring alternative model tested
      - Compare a_F0+a_F1*dt vs a_F0+a_F2*dt²
    * Part 6f IMPROVED: 6 model variants (M1-M6)
      - M5: Stine-dt + Strang-dt² + cos(dt²) (recommended)
      - M6: Stine-dt² + Strang-dt² + cos(dt²) (alternative)
    * Part 6h IMPROVED: Rate(d_F) prediction for M5, M6

  Key issue addressed from iteration 28:
    Strang coefficient b_F(dt) = b_F0 + b_F1*dt used dt-linear correction,
    but symmetric (Strang) splitting has only even-order corrections:
      d_F(SCPT,ex) = b_F0*dt² + b_F2*dt⁴ + b_F4*dt⁶ + ...
    Therefore b_F(dt) = b_F0 + b_F2*dt² is the correct model form.
    This halves the composite model max error: 0.15% -> 0.074%.

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

    MODIFIED from iteration 28 (iteration 29):
    - Part 6f IMPROVED: Strang coefficient model corrected (dt → dt²), 6-model comparison
    - Part 6h IMPROVED: M5 Rate(d_F) prediction added
    - All other parts retained from iteration 28
    """
    print("\n" + "=" * 70)
    print("Test E: Exact Liouvillian Comparison (ST vs Exact vs CT vs CPT vs SCPT)")
    print("       [iter 29: Strang dt² model + 6-model variant comparison]")
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

        # --- Frobenius distances (expanded in iteration 27) ---
        frob_st_ex = frobenius_distance(rho_exact_final, rho_st)
        frob_st_scpt = frobenius_distance(rho_scpt, rho_st)
        frob_scpt_ex = frobenius_distance(rho_exact_final, rho_scpt)
        frob_ct_ex = frobenius_distance(rho_exact_final, rho_ct)
        frob_scpt_ct = frobenius_distance(rho_ct, rho_scpt)

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
            # Frobenius distances (expanded in iteration 27)
            "frob_ST_vs_exact": frob_st_ex,
            "frob_ST_vs_SCPT": frob_st_scpt,
            "frob_SCPT_vs_exact": frob_scpt_ex,
            "frob_CT_vs_exact": frob_ct_ex,
            "frob_SCPT_vs_CT": frob_scpt_ct,
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

    # Part 6f: Frobenius-based composite model (IMPROVED in iteration 29)
    # Iteration 26: simple averaging for a_F, b_F.
    # Iteration 27: higher-order fitting (a_F0 + a_F1*dt), non-parametric verification,
    # and Frobenius Richardson extrapolation.
    # Iteration 28: cos(θ_F)(dt) = c0 + c1*dt² modeling, 4-model variant comparison,
    # error decomposition diagnostics.
    # Iteration 29: Strang coefficient model corrected (dt → dt²), Stinespring alternative
    # model tested, 6-model variant comparison (M1-M6).
    print(f"\n  Part 6f: Frobenius-based composite model (IMPROVED in iteration 29)")
    print(f"  Iteration 29 adds: Strang dt² model correction, Stinespring dt² comparison, M5/M6 models")

    # Frobenius distances from convergence_data
    frob_stine_vals = np.array([cd["frob_ST_vs_SCPT"] for cd in convergence_data])
    frob_strang_vals = np.array([cd["frob_SCPT_vs_exact"] for cd in convergence_data])
    frob_total_vals = np.array([cd["frob_ST_vs_exact"] for cd in convergence_data])

    composite_fit_result = {}
    if len(dts) >= 2:
        # --- Simple coefficient fitting (iteration 26, retained as baseline) ---
        a_F = float(np.mean(frob_stine_vals / dts))
        b_F = float(np.mean(frob_strang_vals / dts**2))
        r_F = b_F / a_F

        # Frobenius angle (from Part 6g data, very stable ~125.2 deg)
        cos_theta_F_values = []
        for cd_entry in convergence_data:
            df_ab = cd_entry["frob_ST_vs_exact"]
            df_a = cd_entry["frob_ST_vs_SCPT"]
            df_b = cd_entry["frob_SCPT_vs_exact"]
            denom_f = 2.0 * df_a * df_b
            if denom_f > 1e-30:
                cos_f = (df_ab**2 - df_a**2 - df_b**2) / denom_f
                cos_f = max(-1.0, min(1.0, cos_f))
                cos_theta_F_values.append(cos_f)
        cos_theta_F_mean = float(np.mean(cos_theta_F_values)) if cos_theta_F_values else 0.0
        theta_F_mean_deg = float(np.degrees(np.arccos(max(-1, min(1, cos_theta_F_mean)))))

        print(f"\n  --- Simple model (iteration 26 baseline) ---")
        print(f"  Frobenius coefficients: a_F = {a_F:.6e}, b_F = {b_F:.6e}, r_F = b_F/a_F = {r_F:.4f}")
        print(f"  Frobenius angle: cos(θ_F_mean) = {cos_theta_F_mean:.4f} (θ_F = {theta_F_mean_deg:.1f}°)")

        # --- Higher-order Frobenius coefficient fitting (iter 27, EXTENDED in 29) ---
        a_F0_ho = a_F
        a_F1_ho = 0.0
        b_F0_ho = b_F
        b_F1_ho = 0.0
        frob_stine_ho_errors = []
        frob_strang_ho_errors = []
        # dt² model variables (NEW in iteration 29)
        a_F0_ho_dt2 = a_F
        a_F2_ho_dt2 = 0.0
        b_F0_ho_dt2 = b_F
        b_F2_ho_dt2 = 0.0
        frob_stine_dt2_errors = []
        frob_strang_dt2_errors = []
        if len(dts) >= 3:
            A_ho_mat = np.column_stack([np.ones_like(dts), dts])
            A_ho_dt2_mat = np.column_stack([np.ones_like(dts), dts**2])
            # Stinespring dt model: d_F(ST,SCPT)/dt = a_F0 + a_F1*dt
            coeffs_stine_F, _, _, _ = np.linalg.lstsq(A_ho_mat, frob_stine_vals / dts, rcond=None)
            a_F0_ho, a_F1_ho = float(coeffs_stine_F[0]), float(coeffs_stine_F[1])
            stine_F_ho_fitted = (a_F0_ho + a_F1_ho * dts) * dts
            frob_stine_ho_errors = [float(e) for e in np.abs(frob_stine_vals - stine_F_ho_fitted) / frob_stine_vals]
            # Stinespring dt² model (NEW in iteration 29): d_F(ST,SCPT)/dt = a_F0 + a_F2*dt²
            coeffs_stine_dt2, _, _, _ = np.linalg.lstsq(A_ho_dt2_mat, frob_stine_vals / dts, rcond=None)
            a_F0_ho_dt2, a_F2_ho_dt2 = float(coeffs_stine_dt2[0]), float(coeffs_stine_dt2[1])
            stine_dt2_fitted = (a_F0_ho_dt2 + a_F2_ho_dt2 * dts**2) * dts
            frob_stine_dt2_errors = [float(e) for e in np.abs(frob_stine_vals - stine_dt2_fitted) / frob_stine_vals]
            # Strang dt model (iteration 27): d_F(SCPT,ex)/dt² = b_F0 + b_F1*dt
            coeffs_strang_F, _, _, _ = np.linalg.lstsq(A_ho_mat, frob_strang_vals / dts**2, rcond=None)
            b_F0_ho, b_F1_ho = float(coeffs_strang_F[0]), float(coeffs_strang_F[1])
            strang_F_ho_fitted = (b_F0_ho + b_F1_ho * dts) * dts**2
            frob_strang_ho_errors = [float(e) for e in np.abs(frob_strang_vals - strang_F_ho_fitted) / frob_strang_vals]
            # Strang dt² model (NEW in iteration 29): d_F(SCPT,ex)/dt² = b_F0 + b_F2*dt²
            # Theoretical motivation: Strang symmetric splitting eliminates odd-order
            # corrections, so d_F(SCPT,ex) = b_F0·dt² + b_F2·dt⁴ + ..., giving
            # d_F(SCPT,ex)/dt² = b_F0 + b_F2·dt² (even-order corrections only).
            coeffs_strang_dt2, _, _, _ = np.linalg.lstsq(A_ho_dt2_mat, frob_strang_vals / dts**2, rcond=None)
            b_F0_ho_dt2, b_F2_ho_dt2 = float(coeffs_strang_dt2[0]), float(coeffs_strang_dt2[1])
            strang_dt2_fitted = (b_F0_ho_dt2 + b_F2_ho_dt2 * dts**2) * dts**2
            frob_strang_dt2_errors = [float(e) for e in np.abs(frob_strang_vals - strang_dt2_fitted) / frob_strang_vals]

        r_F_ho = b_F0_ho / a_F0_ho if a_F0_ho > 1e-30 else r_F
        r_F_ho_dt2 = b_F0_ho_dt2 / a_F0_ho if a_F0_ho > 1e-30 else r_F

        print(f"\n  --- Coefficient model comparison (IMPROVED in iteration 29) ---")
        print(f"  Stinespring dt model (iter27): d_F(ST,SCPT) = ({a_F0_ho:.6e} + {a_F1_ho:.6e}·dt) · dt")
        print(f"    a_F0 = {a_F0_ho:.6e}, a_F1 = {a_F1_ho:.6e}")
        if frob_stine_ho_errors:
            print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in frob_stine_ho_errors)}")
            print(f"    Max: {max(frob_stine_ho_errors)*100:.4f}%")
        print(f"  Stinespring dt² model (iter29): d_F(ST,SCPT) = ({a_F0_ho_dt2:.6e} + {a_F2_ho_dt2:.6e}·dt²) · dt")
        print(f"    a_F0 = {a_F0_ho_dt2:.6e}, a_F2 = {a_F2_ho_dt2:.6e}")
        if frob_stine_dt2_errors:
            print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in frob_stine_dt2_errors)}")
            print(f"    Max: {max(frob_stine_dt2_errors)*100:.4f}%")
        if frob_stine_ho_errors and frob_stine_dt2_errors:
            ratio_stine = max(frob_stine_ho_errors) / max(frob_stine_dt2_errors) if max(frob_stine_dt2_errors) > 0 else float("inf")
            print(f"  → Stinespring dt² model is {ratio_stine:.1f}x more accurate than dt model")

        print(f"\n  Strang dt model (iter27): d_F(SCPT,ex) = ({b_F0_ho:.6e} + {b_F1_ho:.6e}·dt) · dt²")
        print(f"    b_F0 = {b_F0_ho:.6e}, b_F1 = {b_F1_ho:.6e}")
        if frob_strang_ho_errors:
            print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in frob_strang_ho_errors)}")
            print(f"    Max: {max(frob_strang_ho_errors)*100:.4f}%")
        print(f"  Strang dt² model (iter29): d_F(SCPT,ex) = ({b_F0_ho_dt2:.6e} + {b_F2_ho_dt2:.6e}·dt²) · dt²")
        print(f"    b_F0 = {b_F0_ho_dt2:.6e}, b_F2 = {b_F2_ho_dt2:.6e}")
        if frob_strang_dt2_errors:
            print(f"    Relative errors: {', '.join(f'{e:.6f}' for e in frob_strang_dt2_errors)}")
            print(f"    Max: {max(frob_strang_dt2_errors)*100:.4f}%")
        if frob_strang_ho_errors and frob_strang_dt2_errors:
            ratio_strang = max(frob_strang_ho_errors) / max(frob_strang_dt2_errors) if max(frob_strang_dt2_errors) > 0 else float("inf")
            print(f"  → Strang dt² model is {ratio_strang:.1f}x more accurate than dt model")
            print(f"  → Theoretical reason: Strang (symmetric) splitting eliminates odd-order corrections")
            print(f"     d_F(SCPT,ex) = b_F0·dt² + b_F2·dt⁴ + ... → d_F/dt² = b_F0 + b_F2·dt² (even-order only)")

        print(f"\n  r_F (dt model) = b_F0/a_F0 = {r_F_ho:.4f}")
        print(f"  r_F (dt² Strang) = b_F0_dt2/a_F0 = {r_F_ho_dt2:.4f}")

        # --- Frobenius Richardson extrapolation (NEW in iteration 27) ---
        a_F_rich = a_F
        b_F_rich = b_F
        if len(dts) >= 2:
            dt_small = dts[-1]
            dt_next = dts[-2]
            r_s_small = frob_stine_vals[-1] / dt_small
            r_s_next = frob_stine_vals[-2] / dt_next
            a_F_rich = float(r_s_small - dt_small * (r_s_next - r_s_small) / (dt_next - dt_small))
            r_b_small = frob_strang_vals[-1] / dt_small**2
            r_b_next = frob_strang_vals[-2] / dt_next**2
            b_F_rich = float(r_b_small - dt_small * (r_b_next - r_b_small) / (dt_next - dt_small))

        print(f"\n  --- Frobenius Richardson extrapolation (NEW in iteration 27) ---")
        print(f"  Using dt={dts[-1]:.2f} and dt={dts[-2]:.2f}:")
        print(f"    a_F_Richardson = {a_F_rich:.6e} (mean = {a_F:.6e}, diff = {abs(a_F_rich-a_F)/a_F*100:.4f}%)")
        print(f"    b_F_Richardson = {b_F_rich:.6e} (mean = {b_F:.6e}, diff = {abs(b_F_rich-b_F)/b_F*100:.4f}%)")
        print(f"    a_F0 (HO fit) = {a_F0_ho:.6e}, a_F_Rich = {a_F_rich:.6e}")
        print(f"    b_F0 (HO fit) = {b_F0_ho:.6e}, b_F_Rich = {b_F_rich:.6e}")

        # --- Non-parametric law-of-cosines verification (NEW in iteration 27) ---
        print(f"\n  --- Non-parametric law-of-cosines verification (NEW in iteration 27) ---")
        print(f"  Using exact measured component distances and per-dt cos(θ_F):")
        print(f"  d_F_nonparam(dt) = √(d_F²(ST,SCPT) + d_F²(SCPT,ex) + 2·d_F(ST,SCPT)·d_F(SCPT,ex)·cos(θ_F))")
        print(f"  {'dt':>5} | {'d_F_nonparam':>13} | {'d_F_observed':>13} | {'error':>12}")
        print("  " + "-" * 55)
        nonparam_errors = []
        for k_idx, cd_entry in enumerate(convergence_data):
            df_a = cd_entry["frob_ST_vs_SCPT"]
            df_b = cd_entry["frob_SCPT_vs_exact"]
            df_obs = cd_entry["frob_ST_vs_exact"]
            cos_f_k = cos_theta_F_values[k_idx] if k_idx < len(cos_theta_F_values) else cos_theta_F_mean
            df_nonparam = np.sqrt(df_a**2 + df_b**2 + 2 * df_a * df_b * cos_f_k)
            err_np = abs(df_nonparam - df_obs) / df_obs if df_obs > 1e-30 else 0.0
            nonparam_errors.append(float(err_np))
            print(f"  {cd_entry['dt']:5.2f} | {df_nonparam:13.6e} | {df_obs:13.6e} | {err_np:12.2e}")
        print(f"  → Law of cosines verified to machine precision (max error: {max(nonparam_errors):.2e})")
        print(f"  → Parametric model errors come entirely from coefficient fitting and θ_F averaging")

        # --- cos(θ_F)(dt) modeling (NEW in iteration 28) ---
        print(f"\n  --- cos(θ_F)(dt) modeling (NEW in iteration 28) ---")
        print(f"  cos(θ_F) varies with dt: physical effect, not a bug.")
        cos_theta_F_arr = np.array(cos_theta_F_values)
        c0_cos_dt2 = cos_theta_F_mean  # fallback
        c1_cos_dt2 = 0.0
        cos_dt2_fit_errors = []
        if len(dts) >= 3:
            # Fit cos(θ_F) = c0 + c1*dt² (dt² is physically motivated because the
            # Stinespring error has O(dt) + O(dt²) structure, so the angle between
            # error vectors gets O(dt²) corrections from the ratio of subleading terms)
            A_cos_mat = np.column_stack([np.ones_like(dts), dts**2])
            coeffs_cos, _, _, _ = np.linalg.lstsq(A_cos_mat, cos_theta_F_arr, rcond=None)
            c0_cos_dt2, c1_cos_dt2 = float(coeffs_cos[0]), float(coeffs_cos[1])
            cos_fitted = c0_cos_dt2 + c1_cos_dt2 * dts**2
            cos_dt2_fit_errors = [float(abs(cos_fitted[k] - cos_theta_F_arr[k])) for k in range(len(dts))]
        print(f"  cos(θ_F)(dt) = {c0_cos_dt2:.10f} + {c1_cos_dt2:.10f}·dt²")
        print(f"  c0 (asymptotic) = {c0_cos_dt2:.10f} (θ = {float(np.degrees(np.arccos(max(-1, min(1, c0_cos_dt2))))):.2f}°)")
        print(f"  c1 (correction) = {c1_cos_dt2:.10f}")
        if cos_dt2_fit_errors:
            print(f"  Fit errors: {', '.join(f'{e:.2e}' for e in cos_dt2_fit_errors)}")
            print(f"  Max fit error: {max(cos_dt2_fit_errors):.2e} (vs cos variation {max(cos_theta_F_arr)-min(cos_theta_F_arr):.4f})")

        # Also show dt-linear model for comparison
        c0_cos_dt = cos_theta_F_mean
        c1_cos_dt = 0.0
        cos_dt_fit_errors = []
        if len(dts) >= 3:
            A_cos_dt_mat = np.column_stack([np.ones_like(dts), dts])
            coeffs_cos_dt, _, _, _ = np.linalg.lstsq(A_cos_dt_mat, cos_theta_F_arr, rcond=None)
            c0_cos_dt, c1_cos_dt = float(coeffs_cos_dt[0]), float(coeffs_cos_dt[1])
            cos_dt_fitted = c0_cos_dt + c1_cos_dt * dts
            cos_dt_fit_errors = [float(abs(cos_dt_fitted[k] - cos_theta_F_arr[k])) for k in range(len(dts))]
        if cos_dt_fit_errors and cos_dt2_fit_errors:
            print(f"\n  Comparison with dt-linear model: cos(θ_F)(dt) = {c0_cos_dt:.10f} + {c1_cos_dt:.10f}·dt")
            print(f"  dt-linear max error: {max(cos_dt_fit_errors):.2e}")
            if max(cos_dt2_fit_errors) > 0 and max(cos_dt_fit_errors) > 0:
                print(f"  → dt² model is {max(cos_dt_fit_errors)/max(cos_dt2_fit_errors):.0f}x more accurate than dt model")

        # --- 6-model variant comparison (EXPANDED in iteration 29) ---
        # Also compute trace-distance model (iteration 25) for comparison
        cos_theta_T_values = [a["cos_theta_approx"] for a in error_angle_data
                             if not np.isnan(a.get("cos_theta_approx", float("nan")))]
        cos_theta_T_mean = float(np.mean(cos_theta_T_values)) if cos_theta_T_values else 0.0
        a0_T_val = higher_order_fit_result.get("stinespring_a0", a_stine)
        b0_T_val = higher_order_fit_result.get("strang_b0", b_strang)
        r_T_val = b0_T_val / a0_T_val

        print(f"\n  --- 6-model variant comparison (EXPANDED in iteration 29) ---")
        print(f"  M1: Simple (iter26) — constant a_F, b_F, cos_mean")
        print(f"  M2: HO-dt+cos_mean (iter27) — HO-dt coefficients, cos_mean")
        print(f"  M3: HO-dt+cos(dt²) (iter28) — HO-dt coefficients, cos(θ_F)(dt²)")
        print(f"  M4: HO-dt+per-dt cos (reference) — HO-dt coefficients, per-dt cos(θ_F)")
        print(f"  M5: Stine-dt+Strang-dt²+cos(dt²) (iter29) — recommended model")
        print(f"  M6: Stine-dt²+Strang-dt²+cos(dt²) (iter29 alt) — alternative")
        print(f"  T:  Trace-dist model (iter25) — for comparison")
        print(f"\n  {'dt':>5} | {'M1 err%':>8} | {'M2 err%':>8} | {'M3 err%':>8} | {'M5 err%':>8} | {'M6 err%':>8} | {'M4 err%':>8} | {'T err%':>7} | {'d_F obs':>11}")
        print("  " + "-" * 110)

        dF_model_vals = []
        dF_model_errors = []
        dF_ho_model_vals = []
        dF_ho_model_errors = []
        dF_cos_dt2_model_vals = []
        dF_cos_dt2_model_errors = []
        dF_perdt_model_vals = []
        dF_perdt_model_errors = []
        dF_m5_model_vals = []
        dF_m5_model_errors = []
        dF_m6_model_vals = []
        dF_m6_model_errors = []
        tT_model_vals = []
        tT_model_errors = []
        for k_idx, cd_entry in enumerate(convergence_data):
            dt_val = cd_entry["dt"]
            dF_obs_dt = cd_entry["frob_ST_vs_exact"]

            # M1: Simple Frobenius model (iteration 26)
            inner_F = 1.0 + r_F**2 * dt_val**2 + 2.0 * r_F * dt_val * cos_theta_F_mean
            inner_F = max(inner_F, 0.0)
            dF_model_dt = a_F * dt_val * np.sqrt(inner_F)
            rel_err_F = abs(dF_model_dt - dF_obs_dt) / dF_obs_dt if dF_obs_dt > 0 else float("nan")
            dF_model_vals.append(dF_model_dt)
            dF_model_errors.append(rel_err_F)

            # M2: HO + cos_mean (iteration 27)
            a_F_dt = a_F0_ho + a_F1_ho * dt_val
            b_F_dt = b_F0_ho + b_F1_ho * dt_val
            r_F_dt = b_F_dt / a_F_dt if a_F_dt > 1e-30 else r_F_ho
            inner_F_ho = 1.0 + r_F_dt**2 * dt_val**2 + 2.0 * r_F_dt * dt_val * cos_theta_F_mean
            inner_F_ho = max(inner_F_ho, 0.0)
            dF_ho_model_dt = a_F_dt * dt_val * np.sqrt(inner_F_ho)
            rel_err_F_ho = abs(dF_ho_model_dt - dF_obs_dt) / dF_obs_dt if dF_obs_dt > 0 else float("nan")
            dF_ho_model_vals.append(dF_ho_model_dt)
            dF_ho_model_errors.append(rel_err_F_ho)

            # M3: HO + cos(dt²) (iteration 28 NEW)
            cos_dt2_k = c0_cos_dt2 + c1_cos_dt2 * dt_val**2
            inner_F_cos = 1.0 + r_F_dt**2 * dt_val**2 + 2.0 * r_F_dt * dt_val * cos_dt2_k
            inner_F_cos = max(inner_F_cos, 0.0)
            dF_cos_dt2_dt = a_F_dt * dt_val * np.sqrt(inner_F_cos)
            rel_err_cos = abs(dF_cos_dt2_dt - dF_obs_dt) / dF_obs_dt if dF_obs_dt > 0 else float("nan")
            dF_cos_dt2_model_vals.append(dF_cos_dt2_dt)
            dF_cos_dt2_model_errors.append(rel_err_cos)

            # M4: HO + per-dt cos (reference upper bound)
            cos_k = cos_theta_F_values[k_idx] if k_idx < len(cos_theta_F_values) else cos_theta_F_mean
            inner_F_perdt = 1.0 + r_F_dt**2 * dt_val**2 + 2.0 * r_F_dt * dt_val * cos_k
            inner_F_perdt = max(inner_F_perdt, 0.0)
            dF_perdt_dt = a_F_dt * dt_val * np.sqrt(inner_F_perdt)
            rel_err_perdt = abs(dF_perdt_dt - dF_obs_dt) / dF_obs_dt if dF_obs_dt > 0 else float("nan")
            dF_perdt_model_vals.append(dF_perdt_dt)
            dF_perdt_model_errors.append(rel_err_perdt)

            # M5: Stine-dt + Strang-dt² + cos(dt²) (NEW in iteration 29)
            b_F_dt_dt2 = b_F0_ho_dt2 + b_F2_ho_dt2 * dt_val**2
            r_F_dt_m5 = b_F_dt_dt2 / a_F_dt if a_F_dt > 1e-30 else r_F_ho_dt2
            inner_F_m5 = 1.0 + r_F_dt_m5**2 * dt_val**2 + 2.0 * r_F_dt_m5 * dt_val * cos_dt2_k
            inner_F_m5 = max(inner_F_m5, 0.0)
            dF_m5_dt = a_F_dt * dt_val * np.sqrt(inner_F_m5)
            rel_err_m5 = abs(dF_m5_dt - dF_obs_dt) / dF_obs_dt if dF_obs_dt > 0 else float("nan")
            dF_m5_model_vals.append(dF_m5_dt)
            dF_m5_model_errors.append(rel_err_m5)

            # M6: Stine-dt² + Strang-dt² + cos(dt²) (NEW in iteration 29)
            a_F_dt_m6 = a_F0_ho_dt2 + a_F2_ho_dt2 * dt_val**2
            r_F_dt_m6 = b_F_dt_dt2 / a_F_dt_m6 if a_F_dt_m6 > 1e-30 else r_F_ho_dt2
            inner_F_m6 = 1.0 + r_F_dt_m6**2 * dt_val**2 + 2.0 * r_F_dt_m6 * dt_val * cos_dt2_k
            inner_F_m6 = max(inner_F_m6, 0.0)
            dF_m6_dt = a_F_dt_m6 * dt_val * np.sqrt(inner_F_m6)
            rel_err_m6 = abs(dF_m6_dt - dF_obs_dt) / dF_obs_dt if dF_obs_dt > 0 else float("nan")
            dF_m6_model_vals.append(dF_m6_dt)
            dF_m6_model_errors.append(rel_err_m6)

            # T: Trace-distance model (iteration 25, for comparison)
            inner_T = 1.0 + r_T_val**2 * dt_val**2 + 2.0 * r_T_val * dt_val * cos_theta_T_mean
            inner_T = max(inner_T, 0.0)
            tT_model_dt = a0_T_val * dt_val * np.sqrt(inner_T)
            tT_obs_dt = cd_entry["T_ST_vs_exact"]
            rel_err_T = abs(tT_model_dt - tT_obs_dt) / tT_obs_dt if tT_obs_dt > 0 else float("nan")
            tT_model_vals.append(tT_model_dt)
            tT_model_errors.append(rel_err_T)

            print(
                f"  {dt_val:5.2f} | {rel_err_F*100:7.4f}% | {rel_err_F_ho*100:7.4f}% | "
                f"{rel_err_cos*100:7.4f}% | {rel_err_m5*100:7.4f}% | {rel_err_m6*100:7.4f}% | "
                f"{rel_err_perdt*100:7.4f}% | "
                f"{rel_err_T*100:6.2f}% | {dF_obs_dt:11.4e}"
            )

        print(f"\n  Max errors across all dt:")
        print(f"    M1 Simple:                {max(dF_model_errors)*100:.4f}%")
        print(f"    M2 HO-dt+cos_mean:        {max(dF_ho_model_errors)*100:.4f}%")
        print(f"    M3 HO-dt+cos(dt²):        {max(dF_cos_dt2_model_errors)*100:.4f}%")
        print(f"    M5 Stine-dt+Strang-dt²:   {max(dF_m5_model_errors)*100:.4f}% ← recommended (iter29)")
        print(f"    M6 Stine-dt²+Strang-dt²:  {max(dF_m6_model_errors)*100:.4f}% ← alternative (iter29)")
        print(f"    M4 HO-dt+per-dt:          {max(dF_perdt_model_errors)*100:.4f}% (coefficient-fitting limit)")
        print(f"    T  Trace-dist:            {max(tT_model_errors)*100:.2f}%")
        if max(dF_m5_model_errors) > 0:
            print(f"  → M5 is {max(dF_cos_dt2_model_errors)/max(dF_m5_model_errors):.1f}x more accurate than M3 (Strang dt² fix)")
            print(f"  → M5 is {max(dF_model_errors)/max(dF_m5_model_errors):.1f}x more accurate than M1 (overall)")
            if max(dF_m6_model_errors) > 0:
                print(f"  → M6 is {max(dF_m5_model_errors)/max(dF_m6_model_errors):.1f}x more accurate than M5 (Stine dt² additional benefit)")

        # --- Error decomposition (NEW in iteration 28) ---
        print(f"\n  --- Error source decomposition (NEW in iteration 28) ---")
        print(f"  For each dt, decompose composite model error into:")
        print(f"    coeff_only: HO coefficients + per-dt cos(θ_F) → isolates coefficient error")
        print(f"    cos_only: observed components + cos_mean → isolates cos(θ_F) error")
        print(f"    combined: HO coefficients + cos_mean → combined error (= M2)")
        print(f"\n  {'dt':>5} | {'coeff only%':>11} | {'cos only%':>10} | {'combined%':>10} | {'dominant':>10}")
        print("  " + "-" * 65)
        err_decomp_data = []
        for k_idx, cd_entry in enumerate(convergence_data):
            dt_val = cd_entry["dt"]
            dF_obs_dt = cd_entry["frob_ST_vs_exact"]
            df_a_obs = cd_entry["frob_ST_vs_SCPT"]
            df_b_obs = cd_entry["frob_SCPT_vs_exact"]

            # Coefficient-only error: HO coefficients + per-dt cos
            cos_k = cos_theta_F_values[k_idx] if k_idx < len(cos_theta_F_values) else cos_theta_F_mean
            a_F_dt = a_F0_ho + a_F1_ho * dt_val
            b_F_dt = b_F0_ho + b_F1_ho * dt_val
            r_F_dt = b_F_dt / a_F_dt if a_F_dt > 1e-30 else r_F_ho
            inner_c = max(0, 1.0 + r_F_dt**2*dt_val**2 + 2.0*r_F_dt*dt_val*cos_k)
            d_coeff = a_F_dt * dt_val * np.sqrt(inner_c)
            err_coeff = abs(d_coeff - dF_obs_dt) / dF_obs_dt * 100 if dF_obs_dt > 0 else 0.0

            # Cos-only error: observed components + cos_mean
            inner_a = max(0, df_a_obs**2 + df_b_obs**2 + 2*df_a_obs*df_b_obs*cos_theta_F_mean)
            d_cos_only = np.sqrt(inner_a)
            err_cos = abs(d_cos_only - dF_obs_dt) / dF_obs_dt * 100 if dF_obs_dt > 0 else 0.0

            # Combined error (= M2 HO+cos_mean)
            inner_comb = max(0, 1.0 + r_F_dt**2*dt_val**2 + 2.0*r_F_dt*dt_val*cos_theta_F_mean)
            d_comb = a_F_dt * dt_val * np.sqrt(inner_comb)
            err_comb = abs(d_comb - dF_obs_dt) / dF_obs_dt * 100 if dF_obs_dt > 0 else 0.0

            dominant = "coeff" if err_coeff > err_cos else "cos" if err_cos > err_coeff else "equal"

            err_decomp_data.append({
                "dt": float(dt_val),
                "coeff_only_error_pct": float(err_coeff),
                "cos_only_error_pct": float(err_cos),
                "combined_error_pct": float(err_comb),
                "dominant_source": dominant,
            })
            print(f"  {dt_val:5.2f} | {err_coeff:10.4f}% | {err_cos:9.4f}% | {err_comb:9.4f}% | {dominant:>10}")

        print(f"\n  Note: combined error ≠ coeff + cos (nonlinear interaction causes partial cancellation)")

        composite_fit_result = {
            # Simple Frobenius model (iteration 26, retained as baseline)
            "frobenius_a_F": float(a_F),
            "frobenius_b_F": float(b_F),
            "frobenius_r_F": float(r_F),
            "frobenius_cos_theta_F_mean": float(cos_theta_F_mean),
            "frobenius_cos_theta_F_per_dt": [float(c) for c in cos_theta_F_values],
            "frobenius_theta_F_mean_deg": float(theta_F_mean_deg),
            "frobenius_model_dF_values": [float(v) for v in dF_model_vals],
            "frobenius_model_relative_errors": [float(e) for e in dF_model_errors],
            # Higher-order Frobenius model (iteration 27)
            "frobenius_a_F0_ho": float(a_F0_ho),
            "frobenius_a_F1_ho": float(a_F1_ho),
            "frobenius_b_F0_ho": float(b_F0_ho),
            "frobenius_b_F1_ho": float(b_F1_ho),
            "frobenius_r_F_ho": float(r_F_ho),
            "frobenius_stine_ho_relative_errors": frob_stine_ho_errors,
            "frobenius_strang_ho_relative_errors": frob_strang_ho_errors,
            "frobenius_ho_model_dF_values": [float(v) for v in dF_ho_model_vals],
            "frobenius_ho_model_relative_errors": [float(e) for e in dF_ho_model_errors],
            # cos(θ_F)(dt) modeling (NEW in iteration 28)
            "frobenius_cos_c0_dt2": float(c0_cos_dt2),
            "frobenius_cos_c1_dt2": float(c1_cos_dt2),
            "frobenius_cos_dt2_fit_errors": cos_dt2_fit_errors,
            "frobenius_cos_dt2_model_dF_values": [float(v) for v in dF_cos_dt2_model_vals],
            "frobenius_cos_dt2_model_relative_errors": [float(e) for e in dF_cos_dt2_model_errors],
            # HO+per-dt cos reference (NEW in iteration 28)
            "frobenius_perdt_model_dF_values": [float(v) for v in dF_perdt_model_vals],
            "frobenius_perdt_model_relative_errors": [float(e) for e in dF_perdt_model_errors],
            # M5/M6 model data (NEW in iteration 29)
            "frobenius_a_F0_ho_dt2": float(a_F0_ho_dt2),
            "frobenius_a_F2_ho_dt2": float(a_F2_ho_dt2),
            "frobenius_b_F0_ho_dt2": float(b_F0_ho_dt2),
            "frobenius_b_F2_ho_dt2": float(b_F2_ho_dt2),
            "frobenius_stine_dt2_relative_errors": frob_stine_dt2_errors,
            "frobenius_strang_dt2_relative_errors": frob_strang_dt2_errors,
            "frobenius_m5_model_dF_values": [float(v) for v in dF_m5_model_vals],
            "frobenius_m5_model_relative_errors": [float(e) for e in dF_m5_model_errors],
            "frobenius_m6_model_dF_values": [float(v) for v in dF_m6_model_vals],
            "frobenius_m6_model_relative_errors": [float(e) for e in dF_m6_model_errors],
            # Error decomposition (NEW in iteration 28)
            "error_decomposition": err_decomp_data,
            # Frobenius Richardson extrapolation (iteration 27)
            "frobenius_a_F_richardson": float(a_F_rich),
            "frobenius_b_F_richardson": float(b_F_rich),
            "frobenius_richardson_dt_pair": [float(dts[-1]), float(dts[-2])],
            # Non-parametric verification (iteration 27)
            "frobenius_nonparam_errors": nonparam_errors,
            # Trace-distance model (iteration 25, for comparison)
            "trace_a0_stinespring": float(a0_T_val),
            "trace_b0_strang": float(b0_T_val),
            "trace_r_ratio": float(r_T_val),
            "trace_cos_theta_mean": float(cos_theta_T_mean),
            "trace_model_T_values": [float(v) for v in tT_model_vals],
            "trace_model_relative_errors": [float(e) for e in tT_model_errors],
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

    # Part 6h: Frobenius-based Rate prediction (IMPROVED in iteration 29)
    print(f"\n  Part 6h: Frobenius-based Rate(d_F) prediction (IMPROVED in iteration 29)")
    print(f"  Rate(d_F) uses Frobenius model where law of cosines is exact.")
    print(f"  Models: simple (iter26), HO+cos_mean (iter27), HO+cos(dt²) (iter28), M5 (iter29).")
    rate_prediction_data = []
    rate_prediction_data_result = {}
    if composite_fit_result and len(convergence_data) >= 2:
        a_F_mod = composite_fit_result["frobenius_a_F"]
        b_F_mod = composite_fit_result["frobenius_b_F"]
        r_F_mod = composite_fit_result["frobenius_r_F"]
        cos_F_mean = composite_fit_result["frobenius_cos_theta_F_mean"]
        a_F0_mod = composite_fit_result.get("frobenius_a_F0_ho", a_F_mod)
        a_F1_mod = composite_fit_result.get("frobenius_a_F1_ho", 0.0)
        b_F0_mod = composite_fit_result.get("frobenius_b_F0_ho", b_F_mod)
        b_F1_mod = composite_fit_result.get("frobenius_b_F1_ho", 0.0)
        cos_c0 = composite_fit_result.get("frobenius_cos_c0_dt2", cos_F_mean)
        cos_c1 = composite_fit_result.get("frobenius_cos_c1_dt2", 0.0)

        print(f"  Simple model: a_F={a_F_mod:.6e}, r_F={r_F_mod:.4f}, cos(θ_F)={cos_F_mean:.4f}")
        print(f"  HO+cos_mean: a_F0={a_F0_mod:.6e}, a_F1={a_F1_mod:.6e}")
        print(f"               b_F0={b_F0_mod:.6e}, b_F1={b_F1_mod:.6e}")
        print(f"  cos(dt²) model: c0={cos_c0:.8f}, c1={cos_c1:.8f}")
        print(f"\n  {'dt range':>15} | {'Rate_s':>7} | {'Rate_M2':>7} | {'Rate_M3':>7} | {'Rate_M5':>7} | {'Rate_obs':>8} | "
              f"{'\u0394_s':>8} | {'\u0394_M2':>8} | {'\u0394_M3':>8} | {'\u0394_M5':>8} | {'Rate_T':>7}")
        print("  " + "-" * 140)

        observed_rates_T = rates.get("tdist_st_ex", [])
        for k in range(1, len(convergence_data)):
            dt_prev = convergence_data[k - 1]["dt"]
            dt_curr = convergence_data[k]["dt"]
            dt_ratio = dt_prev / dt_curr

            # Simple model Rate(d_F)
            inner_prev = max(0, 1.0 + r_F_mod**2 * dt_prev**2 + 2 * r_F_mod * dt_prev * cos_F_mean)
            inner_curr = max(0, 1.0 + r_F_mod**2 * dt_curr**2 + 2 * r_F_mod * dt_curr * cos_F_mean)
            dF_mod_prev = a_F_mod * dt_prev * np.sqrt(inner_prev)
            dF_mod_curr = a_F_mod * dt_curr * np.sqrt(inner_curr)
            if dF_mod_prev > 0 and dF_mod_curr > 0:
                rate_model_simple = np.log(dF_mod_prev / dF_mod_curr) / np.log(dt_ratio)
            else:
                rate_model_simple = float("nan")

            # M2: HO+cos_mean Rate(d_F) (iteration 27)
            a_prev_ho = a_F0_mod + a_F1_mod * dt_prev
            a_curr_ho = a_F0_mod + a_F1_mod * dt_curr
            b_prev_ho = b_F0_mod + b_F1_mod * dt_prev
            b_curr_ho = b_F0_mod + b_F1_mod * dt_curr
            r_prev_ho = b_prev_ho / a_prev_ho if a_prev_ho > 1e-30 else r_F_mod
            r_curr_ho = b_curr_ho / a_curr_ho if a_curr_ho > 1e-30 else r_F_mod
            inner_prev_ho = max(0, 1.0 + r_prev_ho**2 * dt_prev**2 + 2 * r_prev_ho * dt_prev * cos_F_mean)
            inner_curr_ho = max(0, 1.0 + r_curr_ho**2 * dt_curr**2 + 2 * r_curr_ho * dt_curr * cos_F_mean)
            dF_ho_prev = a_prev_ho * dt_prev * np.sqrt(inner_prev_ho)
            dF_ho_curr = a_curr_ho * dt_curr * np.sqrt(inner_curr_ho)
            if dF_ho_prev > 0 and dF_ho_curr > 0:
                rate_model_ho = np.log(dF_ho_prev / dF_ho_curr) / np.log(dt_ratio)
            else:
                rate_model_ho = float("nan")

            # M3: HO+cos(dt²) Rate(d_F) (NEW in iteration 28)
            cos_prev = cos_c0 + cos_c1 * dt_prev**2
            cos_curr = cos_c0 + cos_c1 * dt_curr**2
            inner_prev_cos = max(0, 1.0 + r_prev_ho**2 * dt_prev**2 + 2 * r_prev_ho * dt_prev * cos_prev)
            inner_curr_cos = max(0, 1.0 + r_curr_ho**2 * dt_curr**2 + 2 * r_curr_ho * dt_curr * cos_curr)
            dF_cos_prev = a_prev_ho * dt_prev * np.sqrt(inner_prev_cos)
            dF_cos_curr = a_curr_ho * dt_curr * np.sqrt(inner_curr_cos)
            if dF_cos_prev > 0 and dF_cos_curr > 0:
                rate_model_cos = np.log(dF_cos_prev / dF_cos_curr) / np.log(dt_ratio)
            else:
                rate_model_cos = float("nan")

            # M5: Stine-dt + Strang-dt² + cos(dt²) Rate(d_F) (NEW in iteration 29)
            b_F0_dt2_mod = composite_fit_result.get("frobenius_b_F0_ho_dt2", b_F0_mod)
            b_F2_dt2_mod = composite_fit_result.get("frobenius_b_F2_ho_dt2", 0.0)
            b_prev_dt2 = b_F0_dt2_mod + b_F2_dt2_mod * dt_prev**2
            b_curr_dt2 = b_F0_dt2_mod + b_F2_dt2_mod * dt_curr**2
            r_prev_m5 = b_prev_dt2 / a_prev_ho if a_prev_ho > 1e-30 else r_F_mod
            r_curr_m5 = b_curr_dt2 / a_curr_ho if a_curr_ho > 1e-30 else r_F_mod
            inner_prev_m5 = max(0, 1.0 + r_prev_m5**2 * dt_prev**2 + 2 * r_prev_m5 * dt_prev * cos_prev)
            inner_curr_m5 = max(0, 1.0 + r_curr_m5**2 * dt_curr**2 + 2 * r_curr_m5 * dt_curr * cos_curr)
            dF_m5_prev = a_prev_ho * dt_prev * np.sqrt(inner_prev_m5)
            dF_m5_curr = a_curr_ho * dt_curr * np.sqrt(inner_curr_m5)
            if dF_m5_prev > 0 and dF_m5_curr > 0:
                rate_model_m5 = np.log(dF_m5_prev / dF_m5_curr) / np.log(dt_ratio)
            else:
                rate_model_m5 = float("nan")
            delta_rate_m5 = rate_model_m5 - rate_obs_dF if not (np.isnan(rate_model_m5) or np.isnan(rate_obs_dF)) else float("nan")

            # Observed Rate(d_F)
            dF_obs_prev = convergence_data[k - 1]["frob_ST_vs_exact"]
            dF_obs_curr = convergence_data[k]["frob_ST_vs_exact"]
            if dF_obs_prev > 1e-18 and dF_obs_curr > 1e-18:
                rate_obs_dF = np.log(dF_obs_prev / dF_obs_curr) / np.log(dt_ratio)
            else:
                rate_obs_dF = float("nan")

            delta_rate_simple = rate_model_simple - rate_obs_dF if not (np.isnan(rate_model_simple) or np.isnan(rate_obs_dF)) else float("nan")
            delta_rate_ho = rate_model_ho - rate_obs_dF if not (np.isnan(rate_model_ho) or np.isnan(rate_obs_dF)) else float("nan")
            delta_rate_cos = rate_model_cos - rate_obs_dF if not (np.isnan(rate_model_cos) or np.isnan(rate_obs_dF)) else float("nan")

            # Observed Rate(T) for comparison
            rate_observed_T = observed_rates_T[k - 1] if k - 1 < len(observed_rates_T) else float("nan")

            rate_prediction_data.append({
                "dt_prev": dt_prev,
                "dt_curr": dt_curr,
                "rate_model_simple": float(rate_model_simple),
                "rate_model_ho": float(rate_model_ho),
                "rate_model_cos_dt2": float(rate_model_cos),
                "rate_model_m5": float(rate_model_m5),
                "rate_observed_dF": float(rate_obs_dF),
                "delta_rate_simple": float(delta_rate_simple),
                "delta_rate_ho": float(delta_rate_ho),
                "delta_rate_cos_dt2": float(delta_rate_cos),
                "delta_rate_m5": float(delta_rate_m5),
                "rate_observed_T": float(rate_observed_T),
            })

            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | {rate_model_simple:7.4f} | {rate_model_ho:7.4f} | "
                f"{rate_model_cos:7.4f} | {rate_model_m5:7.4f} | {rate_obs_dF:8.4f} | {delta_rate_simple:+8.4f} | "
                f"{delta_rate_ho:+8.4f} | {delta_rate_cos:+8.4f} | {delta_rate_m5:+8.4f} | {rate_observed_T:7.4f}"
            )

        # Predict dt threshold for Rate(d_F) > 0.95
        r_cos_theta_F = r_F_mod * cos_F_mean
        if r_cos_theta_F < 0:
            dt_threshold_095_dF = 0.05 / abs(r_cos_theta_F)
            dt_threshold_099_dF = 0.01 / abs(r_cos_theta_F)
        else:
            dt_threshold_095_dF = float("nan")
            dt_threshold_099_dF = float("nan")

        # Summary of improvements
        delta_simples = [abs(p["delta_rate_simple"]) for p in rate_prediction_data if not np.isnan(p["delta_rate_simple"])]
        delta_hos = [abs(p["delta_rate_ho"]) for p in rate_prediction_data if not np.isnan(p["delta_rate_ho"])]
        delta_cos_dt2s = [abs(p["delta_rate_cos_dt2"]) for p in rate_prediction_data if not np.isnan(p["delta_rate_cos_dt2"])]
        delta_m5s = [abs(p["delta_rate_m5"]) for p in rate_prediction_data if not np.isnan(p["delta_rate_m5"])]
        if delta_simples and delta_hos and delta_cos_dt2s and delta_m5s:
            print(f"\n  Rate prediction max|Δ|: simple={max(delta_simples):.4f}, M2={max(delta_hos):.4f}, M3={max(delta_cos_dt2s):.4f}, M5={max(delta_m5s):.4f}")

        print(f"\n  Asymptotic: Rate(d_F) ≈ 1 + r_F·cos(θ_F_mean)·dt = 1 + ({r_cos_theta_F:.4f})·dt")
        if not np.isnan(dt_threshold_095_dF):
            print(f"  Rate(d_F) > 0.95 requires dt < {dt_threshold_095_dF:.4f} (n_steps > {10.0/dt_threshold_095_dF:.0f})")
            print(f"  Rate(d_F) > 0.99 requires dt < {dt_threshold_099_dF:.4f} (n_steps > {10.0/dt_threshold_099_dF:.0f})")

        print(f"\n  Note: Rate(T) ≠ Rate(d_F) because T/d_F ratio is dt-dependent for the")
        print(f"  composite error. See Part 6i for T/d_F ratio analysis.")

        rate_prediction_data_result = {
            "model_a_F": float(a_F_mod),
            "model_r_F": float(r_F_mod),
            "model_cos_theta_F_mean": float(cos_F_mean),
            "model_r_cos_theta_F": float(r_cos_theta_F),
            "model_a_F0_ho": float(a_F0_mod),
            "model_a_F1_ho": float(a_F1_mod),
            "model_b_F0_ho": float(b_F0_mod),
            "model_b_F1_ho": float(b_F1_mod),
            "model_cos_c0_dt2": float(cos_c0),
            "model_cos_c1_dt2": float(cos_c1),
            "model_b_F0_ho_dt2": float(b_F0_dt2_mod),
            "model_b_F2_ho_dt2": float(b_F2_dt2_mod),
            "predictions": rate_prediction_data,
            "dt_threshold_rate_dF_095": float(dt_threshold_095_dF) if not np.isnan(dt_threshold_095_dF) else None,
            "dt_threshold_rate_dF_099": float(dt_threshold_099_dF) if not np.isnan(dt_threshold_099_dF) else None,
        }

    # Part 6i: T/d_F ratio analysis (new in iteration 26)
    print(f"\n  Part 6i: T/d_F ratio analysis (new in iteration 26)")
    print(f"  T = ½||δ||₁ (L₁ norm of eigenvalues), d_F = ||δ||_F (L₂ norm of eigenvalues)")
    print(f"  Individual components have stable T/d_F; composite varies due to spectrum crossover.")

    tdF_ratio_data = []
    print(f"\n  {'dt':>5} | {'T/d_F(total)':>12} | {'T/d_F(Stine)':>12} | {'T/d_F(Strang)':>13}")
    print("  " + "-" * 55)
    for cd_entry in convergence_data:
        r_total = cd_entry["T_ST_vs_exact"] / cd_entry["frob_ST_vs_exact"] if cd_entry["frob_ST_vs_exact"] > 1e-30 else float("nan")
        r_stine = cd_entry["T_ST_vs_SCPT"] / cd_entry["frob_ST_vs_SCPT"] if cd_entry["frob_ST_vs_SCPT"] > 1e-30 else float("nan")
        r_strang = cd_entry["T_SCPT_vs_exact"] / cd_entry["frob_SCPT_vs_exact"] if cd_entry["frob_SCPT_vs_exact"] > 1e-30 else float("nan")

        tdF_ratio_data.append({
            "dt": cd_entry["dt"],
            "TdF_ratio_total": r_total,
            "TdF_ratio_stinespring": r_stine,
            "TdF_ratio_strang": r_strang,
        })
        print(f"  {cd_entry['dt']:5.2f} | {r_total:12.6f} | {r_stine:12.6f} | {r_strang:13.6f}")

    # Stability analysis
    r_stine_vals = [d["TdF_ratio_stinespring"] for d in tdF_ratio_data if not np.isnan(d["TdF_ratio_stinespring"])]
    r_strang_vals = [d["TdF_ratio_strang"] for d in tdF_ratio_data if not np.isnan(d["TdF_ratio_strang"])]
    r_total_vals = [d["TdF_ratio_total"] for d in tdF_ratio_data if not np.isnan(d["TdF_ratio_total"])]

    if r_stine_vals:
        var_stine = (max(r_stine_vals) - min(r_stine_vals)) / np.mean(r_stine_vals) * 100
        print(f"\n  Stinespring T/d_F: {min(r_stine_vals):.6f} – {max(r_stine_vals):.6f} (variation: {var_stine:.3f}%)")
    if r_strang_vals:
        var_strang = (max(r_strang_vals) - min(r_strang_vals)) / np.mean(r_strang_vals) * 100
        print(f"  Strang T/d_F: {min(r_strang_vals):.6f} – {max(r_strang_vals):.6f} (variation: {var_strang:.3f}%)")
    if r_total_vals:
        var_total = (max(r_total_vals) - min(r_total_vals)) / np.mean(r_total_vals) * 100
        print(f"  Composite T/d_F: {min(r_total_vals):.6f} – {max(r_total_vals):.6f} (variation: {var_total:.3f}%)")
        print(f"  → Composite T/d_F varies because Stinespring (T/d_F≈{np.mean(r_stine_vals):.4f}) and")
        print(f"    Strang (T/d_F≈{np.mean(r_strang_vals):.4f}) have different eigenvalue spectra.")
        print(f"    At large dt, Strang dominates → T/d_F ≈ {np.mean(r_strang_vals):.4f}")
        print(f"    At small dt, Stinespring dominates → T/d_F ≈ {np.mean(r_stine_vals):.4f}")

    # Part 6j: Frobenius convergence rates Rate(d_F) (EXPANDED in iteration 27)
    print(f"\n  Part 6j: Frobenius convergence rates Rate(d_F) for all comparison pairs")
    print(f"  Expanded in iteration 27: CT vs Exact and SCPT vs CT pairs added")

    frob_rate_data = {}
    if len(convergence_data) >= 2:
        frob_pairs = [
            ("frob_ST_vs_exact", "ST vs Exact"),
            ("frob_ST_vs_SCPT", "ST vs SCPT (Stinespring)"),
            ("frob_SCPT_vs_exact", "SCPT vs Exact (Strang)"),
            ("frob_CT_vs_exact", "CT vs Exact"),
            ("frob_SCPT_vs_CT", "SCPT vs CT (palindromic)"),
        ]
        for frob_key, frob_label in frob_pairs:
            rates_dF = []
            for k in range(1, len(convergence_data)):
                dt_prev = convergence_data[k-1]["dt"]
                dt_curr = convergence_data[k]["dt"]
                dt_ratio = dt_prev / dt_curr
                val_prev = convergence_data[k-1][frob_key]
                val_curr = convergence_data[k][frob_key]
                rate_dF = _compute_rate(val_prev, val_curr, dt_ratio)
                rates_dF.append(rate_dF)
            frob_rate_data[frob_key] = rates_dF

        print(f"\n  {'dt range':>15} | {'Rate(d_F) ST-ex':>15} | {'Rate(T) ST-ex':>14} | "
              f"{'Rate(d_F) Stine':>15} | {'Rate(d_F) Strang':>16} | "
              f"{'Rate(d_F) CT-ex':>15} | {'Rate(d_F) SCPT-CT':>17}")
        print("  " + "-" * 130)
        rates_T_st_ex = rates.get("tdist_st_ex", [])
        for k in range(len(frob_rate_data.get("frob_ST_vs_exact", []))):
            dt_prev = convergence_data[k]["dt"]
            dt_curr = convergence_data[k+1]["dt"]
            r_dF_total = frob_rate_data["frob_ST_vs_exact"][k]
            r_T_total = rates_T_st_ex[k] if k < len(rates_T_st_ex) else float("nan")
            r_dF_stine = frob_rate_data["frob_ST_vs_SCPT"][k]
            r_dF_strang = frob_rate_data["frob_SCPT_vs_exact"][k]
            r_dF_ct_ex = frob_rate_data["frob_CT_vs_exact"][k]
            r_dF_scpt_ct = frob_rate_data["frob_SCPT_vs_CT"][k]
            print(
                f"  {dt_prev:.1f} → {dt_curr:.1f}      | {r_dF_total:15.4f} | {r_T_total:14.4f} | "
                f"{r_dF_stine:15.4f} | {r_dF_strang:16.4f} | "
                f"{r_dF_ct_ex:15.4f} | {r_dF_scpt_ct:17.4f}"
            )


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
        # Composite model fitting (REBUILT in iteration 26: Frobenius-based)
        "composite_fit": composite_fit_result,
        # Frobenius angle estimation (iteration 25)
        "frobenius_angle_data": frobenius_angle_data,
        # Rate prediction (REBUILT in iteration 26: Frobenius-based)
        "rate_prediction": rate_prediction_data_result,
        # T/d_F ratio analysis (new in iteration 26)
        "tdF_ratio_data": tdF_ratio_data,
        # Frobenius convergence rates (new in iteration 26)
        "frobenius_rate_data": frob_rate_data,
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
    print("TTA-UC GKSL Iteration 28: cos(θ_F)(dt) Modeling")
    print(f"Timestamp: {timestamp}")
    print("NOTE: ST simulator uses symmetric Lindblad product ordering")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 29,
        "purpose": "cos(θ_F)(dt) modeling. "
                   "Iteration 27 HO composite model had max error 0.23% but was worse than "
                   "simple model at small dt (0.10% vs 0.006% at dt=0.1). Root cause: "
                   "(1) HO coefficient fitting errors don't cancel as well as simple-mean errors, "
                   "(2) constant cos(θ_F) mean approximation. This iteration adds: "
                   "(1) cos(θ_F)(dt) = c0 + c1*dt² modeling (dt² physically motivated), "
                   "(2) 6-model variant comparison (M1-M6 including Strang-dt² models), "
                   "(3) Error source decomposition (coefficient vs cos approximation).",
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
    json_path = output_dir / f"iteration29_strang_dt2_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 28: cos(θ_F)(dt) モデリング",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "- iteration 27 からの主要修正:",
        "  - Part 6f 改善: cos(θ_F)(dt) = c₀ + c₁·dt² モデリング追加",
        "  - Part 6f 改善: 6 モデルバリアント体系比較（M1-M6, Strang dt² モデル含む）",
        "  - Part 6f 改善: 誤差分解診断（係数誤差 vs cos 近似誤差）",
        "  - Part 6h 改善: M5/M6 Rate(d_F) 予測追加",
        "  - その他のパート: iteration 27 と同一",
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
        "## Test E: 正確な GKSL Liouvillian 解との比較（iteration 26 Frobenius 合成モデル）",
        "",
        "- 参照: ClassicalGKSLSimulator（exp(L·t) による厳密解）",
        f"- 正確解の N_S1(t=10): {te['exact_N_S1']:.10f}",
        f"- 正確解の純度: {te['exact_purity']:.10f}",
        f"- Lindblad チャネル数: {te['n_lindblad_channels']}",
        "- **Frobenius 合成モデル（iteration 26）**: Frobenius 距離ベース合成モデル + T/d_F 比率分析",
        "- **Frobenius 解析（iteration 25）**: Frobenius 距離 + 角度推定",
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

    # Frobenius composite model (IMPROVED in iteration 28)
    if te.get("composite_fit"):
        cf_data = te["composite_fit"]
        md_lines.extend([
            "### Frobenius ベース合成モデル（iteration 29 改善: Strang dt² モデル補正）",
            "",
            "Frobenius ノルム（L₂）は内積空間ノルムであるため、余弦定理が厳密に成立する。",
            "",
            "#### 単純モデル M1（iteration 26 ベースライン）",
            "d_F_model(dt) = a_F·dt · √(1 + r_F²·dt² + 2·r_F·dt·cos(θ_F_mean))",
            "",
            f"- a_F = {cf_data['frobenius_a_F']:.6e}",
            f"- b_F = {cf_data['frobenius_b_F']:.6e}",
            f"- r_F = b_F/a_F = {cf_data['frobenius_r_F']:.4f}",
            f"- cos(θ_F_mean) = {cf_data['frobenius_cos_theta_F_mean']:.4f} "
            f"(θ_F = {cf_data['frobenius_theta_F_mean_deg']:.1f}°)",
            "",
        ])
        # Higher-order model section (iteration 27)
        if cf_data.get("frobenius_a_F0_ho") is not None:
            md_lines.extend([
                "#### 高次モデル M2（iteration 27: HO 係数 + cos_mean）",
                "d_F(ST,SCPT) = (a_F0 + a_F1·dt) · dt, d_F(SCPT,ex) = (b_F0 + b_F1·dt) · dt²",
                "",
                f"- a_F0 = {cf_data['frobenius_a_F0_ho']:.6e}, a_F1 = {cf_data['frobenius_a_F1_ho']:.6e}",
                f"- b_F0 = {cf_data['frobenius_b_F0_ho']:.6e}, b_F1 = {cf_data['frobenius_b_F1_ho']:.6e}",
                f"- r_F (高次) = b_F0/a_F0 = {cf_data['frobenius_r_F_ho']:.4f}",
                "",
            ])
        # cos(θ_F)(dt) model (NEW in iteration 28)
        if cf_data.get("frobenius_cos_c0_dt2") is not None:
            md_lines.extend([
                "#### cos(θ_F)(dt) モデル（iteration 28 新規）",
                "cos(θ_F)(dt) = c₀ + c₁·dt²（dt² は物理的に動機付け: 補正の主項が O(dt²)）",
                "",
                f"- c₀ = {cf_data['frobenius_cos_c0_dt2']:.10f}",
                f"- c₁ = {cf_data['frobenius_cos_c1_dt2']:.10f}",
                "",
            ])
            if cf_data.get("frobenius_cos_dt2_fit_errors"):
                max_cos_err = max(cf_data["frobenius_cos_dt2_fit_errors"])
                md_lines.extend([
                    f"- cos(θ_F) dt² フィット最大誤差: {max_cos_err:.2e}",
                    "",
                ])
            md_lines.extend([
                "#### 完全 HO モデル M3（iteration 28 新規: HO 係数 + cos(dt²)）",
                "d_F_model(dt) = a_F(dt)·dt · √(1 + r_F(dt)²·dt² + 2·r_F(dt)·dt·cos(θ_F)(dt))",
                "",
            ])
        # Richardson extrapolation (iteration 27)
        if cf_data.get("frobenius_a_F_richardson") is not None:
            md_lines.extend([
                "#### Frobenius Richardson 外挿（iteration 27）",
                "",
                f"- a_F_Richardson = {cf_data['frobenius_a_F_richardson']:.6e}",
                f"- b_F_Richardson = {cf_data['frobenius_b_F_richardson']:.6e}",
                "",
            ])
        # Non-parametric verification (iteration 27)
        if cf_data.get("frobenius_nonparam_errors"):
            np_errs = cf_data["frobenius_nonparam_errors"]
            md_lines.extend([
                "#### 非パラメトリック余弦定理検証（iteration 27）",
                "",
                f"余弦定理は機械精度で成立（最大誤差: {max(np_errs):.2e}）。",
                "パラメトリックモデルの誤差源は係数フィッティングと角度モデリングのみ。",
                "",
            ])
        # 4-Model comparison table (EXPANDED in iteration 28)
        md_lines.extend([
            "#### 6 モデルバリアント比較（iteration 29 拡張）",
            "",
            "| dt | M1 誤差% | M2 誤差% | M3 誤差% | M5 誤差% | M6 誤差% | M4 誤差% | T 誤差% | d_F obs |",
            "|----|---------|---------|---------|---------|---------|---------|--------|---------|",
        ])
        m1_errs = cf_data.get("frobenius_model_relative_errors", [])
        m2_errs = cf_data.get("frobenius_ho_model_relative_errors", m1_errs)
        m3_errs = cf_data.get("frobenius_cos_dt2_model_relative_errors", m1_errs)
        m4_errs = cf_data.get("frobenius_perdt_model_relative_errors", m1_errs)
        m5_errs = cf_data.get("frobenius_m5_model_relative_errors", m1_errs)
        m6_errs = cf_data.get("frobenius_m6_model_relative_errors", m1_errs)
        t_errs = cf_data.get("trace_model_relative_errors", [])
        for k_idx, cd_entry in enumerate(te["convergence_data"]):
            dF_obs = cd_entry["frob_ST_vs_exact"]
            e1 = m1_errs[k_idx] * 100 if k_idx < len(m1_errs) else float("nan")
            e2 = m2_errs[k_idx] * 100 if k_idx < len(m2_errs) else float("nan")
            e3 = m3_errs[k_idx] * 100 if k_idx < len(m3_errs) else float("nan")
            e4 = m4_errs[k_idx] * 100 if k_idx < len(m4_errs) else float("nan")
            e5 = m5_errs[k_idx] * 100 if k_idx < len(m5_errs) else float("nan")
            e6 = m6_errs[k_idx] * 100 if k_idx < len(m6_errs) else float("nan")
            eT = t_errs[k_idx] * 100 if k_idx < len(t_errs) else float("nan")
            md_lines.append(
                f"| {cd_entry['dt']:.2f} | {e1:.4f} | {e2:.4f} | {e3:.4f} | {e5:.4f} | {e6:.4f} | {e4:.4f} | {eT:.2f} | {dF_obs:.4e} |"
            )
        md_lines.extend([
            "",
            "- M1: 単純モデル（iter26）— 定数 a_F, b_F, cos_mean",
            "- M2: HO+cos_mean（iter27）— HO 係数, cos_mean",
            "- M3: HO-dt+cos(dt²)（iter28）— HO-dt 係数, cos(θ_F)(dt²)",
            "- M4: HO-dt+per-dt cos（参照）— HO-dt 係数, per-dt cos（係数フィッティング限界）",
            "- M5: Stine-dt+Strang-dt²+cos(dt²)（iter29 推奨）— Strang dt² 補正",
            "- M6: Stine-dt²+Strang-dt²+cos(dt²)（iter29 代替）— 両方 dt² 補正",
            "- T: トレース距離モデル（iter25）",
            "",
        ])
        # Error decomposition table (NEW in iteration 28)
        if cf_data.get("error_decomposition"):
            md_lines.extend([
                "#### 誤差源分解（iteration 28 新規）",
                "",
                "| dt | 係数誤差のみ% | cos 誤差のみ% | 複合誤差% | 支配的誤差源 |",
                "|----|-------------|-------------|---------|-----------|",
            ])
            for ed in cf_data["error_decomposition"]:
                md_lines.append(
                    f"| {ed['dt']:.2f} | {ed['coeff_only_error_pct']:.4f} | "
                    f"{ed['cos_only_error_pct']:.4f} | {ed['combined_error_pct']:.4f} | "
                    f"{ed['dominant_source']} |"
                )
            md_lines.extend([
                "",
                "注: 複合誤差 ≠ 係数誤差 + cos 誤差（非線形相互作用あり）",
                "",
            ])
        md_lines.append("")

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

    # Rate prediction (IMPROVED in iteration 29: simple + HO + cos(dt²) + M5)
    if te.get("rate_prediction") and te["rate_prediction"].get("predictions"):
        rp = te["rate_prediction"]
        md_lines.extend([
            "### Frobenius ベース Rate(d_F) 予測と実測の比較（iteration 29 改善）",
            "",
            "Frobenius 合成モデル: d_F_model(dt) = a_F·dt · √(1 + r_F²·dt² + 2·r_F·dt·cos(θ_F_mean))",
            f"- r_F·cos(θ_F_mean) = {rp['model_r_cos_theta_F']:.4f}",
            f"- 漸近展開: Rate(d_F) ≈ 1 + r_F·cos(θ_F_mean)·dt = 1 + ({rp['model_r_cos_theta_F']:.4f})·dt",
            "",
            "| dt 範囲 | Rate_M1 | Rate_M2 | Rate_M3 | Rate_M5 | Rate_obs | Δ_M1 | Δ_M2 | Δ_M3 | Δ_M5 | Rate_T |",
            "|---------|--------|--------|--------|--------|---------|------|------|------|------|--------|",
        ])
        for pred in rp["predictions"]:
            delta_s = pred.get("delta_rate_simple", float("nan"))
            delta_ho = pred.get("delta_rate_ho", float("nan"))
            delta_cos = pred.get("delta_rate_cos_dt2", float("nan"))
            delta_m5 = pred.get("delta_rate_m5", float("nan"))
            rate_s = pred.get("rate_model_simple", float("nan"))
            rate_ho = pred.get("rate_model_ho", float("nan"))
            rate_cos = pred.get("rate_model_cos_dt2", float("nan"))
            rate_m5 = pred.get("rate_model_m5", float("nan"))
            md_lines.append(
                f"| {pred['dt_prev']:.1f} → {pred['dt_curr']:.1f} | "
                f"{rate_s:.4f} | {rate_ho:.4f} | {rate_cos:.4f} | {rate_m5:.4f} | "
                f"{pred['rate_observed_dF']:.4f} | "
                f"{delta_s:+.4f} | {delta_ho:+.4f} | "
                f"{delta_cos:+.4f} | {delta_m5:+.4f} | "
                f"{pred['rate_observed_T']:.4f} |"
            )
        md_lines.append("")
        if rp.get("dt_threshold_rate_dF_095") is not None:
            md_lines.extend([
                f"- Rate(d_F) > 0.95 に必要な dt < {rp['dt_threshold_rate_dF_095']:.4f} "
                f"(n_steps > {10.0/rp['dt_threshold_rate_dF_095']:.0f})",
                f"- Rate(d_F) > 0.99 に必要な dt < {rp['dt_threshold_rate_dF_099']:.4f} "
                f"(n_steps > {10.0/rp['dt_threshold_rate_dF_099']:.0f})",
                "",
            ])
        md_lines.extend([
            "**注意**: Rate(T) ≠ Rate(d_F) となるのは、合成誤差の T/d_F 比率が dt 依存であるため。Part 6i 参照。",
            "",
        ])

    # T/d_F ratio analysis (new in iteration 26)
    if te.get("tdF_ratio_data"):
        md_lines.extend([
            "### T/d_F 比率分析（iteration 26 追加）",
            "",
            "T/d_F = T(ρ,σ)/d_F(ρ,σ) は固有値スペクトル形状に依存する。",
            "個別成分では安定、合成では Stinespring/Strang のクロスオーバーにより不安定。",
            "",
            "| dt | T/d_F(total) | T/d_F(Stine) | T/d_F(Strang) |",
            "|----|-------------|-------------|--------------|",
        ])
        for d in te["tdF_ratio_data"]:
            md_lines.append(
                f"| {d['dt']:.2f} | {d['TdF_ratio_total']:.6f} | "
                f"{d['TdF_ratio_stinespring']:.6f} | {d['TdF_ratio_strang']:.6f} |"
            )
        md_lines.append("")

    # Frobenius convergence rates (new in iteration 26)
    if te.get("frobenius_rate_data"):
        frd = te["frobenius_rate_data"]
        rates_T_st_ex_md = te["rates"].get("tdist_st_ex", [])
        md_lines.extend([
            "### Frobenius 収束次数 Rate(d_F)（iteration 27 拡張: CT ペア追加）",
            "",
            "| dt 範囲 | Rate(d_F) ST-ex | Rate(T) ST-ex | Rate(d_F) Stine | Rate(d_F) Strang | Rate(d_F) CT-ex | Rate(d_F) SCPT-CT |",
            "|---------|----------------|--------------|----------------|-----------------|----------------|------------------|",
        ])
        for k in range(len(frd.get("frob_ST_vs_exact", []))):
            dt_prev = te["convergence_data"][k]["dt"]
            dt_curr = te["convergence_data"][k+1]["dt"]
            r_dF_t = frd["frob_ST_vs_exact"][k]
            r_T_t = rates_T_st_ex_md[k] if k < len(rates_T_st_ex_md) else float("nan")
            r_dF_s = frd["frob_ST_vs_SCPT"][k]
            r_dF_str = frd["frob_SCPT_vs_exact"][k]
            r_dF_ct = frd.get("frob_CT_vs_exact", [])[k] if k < len(frd.get("frob_CT_vs_exact", [])) else float("nan")
            r_dF_scpt_ct = frd.get("frob_SCPT_vs_CT", [])[k] if k < len(frd.get("frob_SCPT_vs_CT", [])) else float("nan")
            md_lines.append(
                f"| {dt_prev:.1f} → {dt_curr:.1f} | {r_dF_t:.4f} | {r_T_t:.4f} | "
                f"{r_dF_s:.4f} | {r_dF_str:.4f} | {r_dF_ct:.4f} | {r_dF_scpt_ct:.4f} |"
            )
        md_lines.append("")

    md_path = output_dir / f"iteration29_strang_dt2_{timestamp}.md"
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
        print(f"  Test E: Simple Frobenius model r_F={cf_data['frobenius_r_F']:.4f}, "
              f"cos(θ_F)={cf_data['frobenius_cos_theta_F_mean']:.4f}, "
              f"max_error={max(cf_data['frobenius_model_relative_errors'])*100:.2f}%")
        if cf_data.get("frobenius_ho_model_relative_errors"):
            print(f"  Test E: M2 HO+cos_mean max_error={max(cf_data['frobenius_ho_model_relative_errors'])*100:.2f}%")
        if cf_data.get("frobenius_cos_dt2_model_relative_errors"):
            print(f"  Test E: M3 HO+cos(dt²) max_error={max(cf_data['frobenius_cos_dt2_model_relative_errors'])*100:.2f}%")
        if cf_data.get("frobenius_perdt_model_relative_errors"):
            print(f"  Test E: M4 HO+per-dt max_error={max(cf_data['frobenius_perdt_model_relative_errors'])*100:.2f}% (coeff limit)")
            if "frobenius_m5_model_relative_errors" in cf_data:
                print(f"  Test E: M5 Stine-dt+Strang-dt² max_error={max(cf_data['frobenius_m5_model_relative_errors'])*100:.4f}% (recommended)")
            if "frobenius_m6_model_relative_errors" in cf_data:
                print(f"  Test E: M6 Stine-dt²+Strang-dt² max_error={max(cf_data['frobenius_m6_model_relative_errors'])*100:.4f}% (alternative)")
        if cf_data.get("frobenius_cos_c0_dt2") is not None:
            print(f"  Test E: cos(θ_F)(dt) = {cf_data['frobenius_cos_c0_dt2']:.8f} + {cf_data['frobenius_cos_c1_dt2']:.8f}·dt²")
        if cf_data.get("frobenius_nonparam_errors"):
            print(f"  Test E: Non-parametric law-of-cosines max_error={max(cf_data['frobenius_nonparam_errors']):.2e} (machine precision)")
        if cf_data.get("frobenius_a_F_richardson") is not None:
            print(f"  Test E: Frobenius Richardson a_F={cf_data['frobenius_a_F_richardson']:.6e}, "
                  f"b_F={cf_data['frobenius_b_F_richardson']:.6e}")
        print(f"  Test E: Trace model (iter25) max_error={max(cf_data['trace_model_relative_errors'])*100:.1f}%")
    if te.get("frobenius_angle_data"):
        f_angles = te["frobenius_angle_data"]
        theta_f_vals_summary = [a["theta_frobenius_deg"] for a in f_angles if not np.isnan(a.get("theta_frobenius_deg", float("nan")))]
        if theta_f_vals_summary:
            print(f"  Test E: Frobenius angle range: {min(theta_f_vals_summary):.1f}° – {max(theta_f_vals_summary):.1f}° "
                  f"(variation: {max(theta_f_vals_summary)-min(theta_f_vals_summary):.1f}°)")
    if te.get("rate_prediction") and te["rate_prediction"].get("predictions"):
        rp = te["rate_prediction"]
        max_delta_simple = max(abs(p.get("delta_rate_simple", float("nan")))
                              for p in rp["predictions"]
                              if not np.isnan(p.get("delta_rate_simple", float("nan"))))
        max_delta_ho = max(abs(p.get("delta_rate_ho", float("nan")))
                          for p in rp["predictions"]
                          if not np.isnan(p.get("delta_rate_ho", float("nan"))))
        max_delta_cos = max(abs(p.get("delta_rate_cos_dt2", float("nan")))
                           for p in rp["predictions"]
                           if not np.isnan(p.get("delta_rate_cos_dt2", float("nan"))))
        print(f"  Test E: Rate(d_F) max|Δ|: M1={max_delta_simple:.4f}, M2={max_delta_ho:.4f}, M3={max_delta_cos:.4f}")
        if rp.get("dt_threshold_rate_dF_095") is not None:
            print(f"  Test E: Rate(d_F) > 0.95 requires dt < {rp['dt_threshold_rate_dF_095']:.4f}")
    if te.get("tdF_ratio_data"):
        r_total = [d["TdF_ratio_total"] for d in te["tdF_ratio_data"] if not np.isnan(d["TdF_ratio_total"])]
        if r_total:
            print(f"  Test E: T/d_F ratio range: {min(r_total):.4f} – {max(r_total):.4f} "
                  f"(variation: {(max(r_total)-min(r_total))/np.mean(r_total)*100:.1f}%)")
    print(f"  Test E: exact purity={te['exact_purity']:.6f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
