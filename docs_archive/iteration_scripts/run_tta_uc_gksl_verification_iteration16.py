#!/usr/bin/env python3
"""TTA-UC GKSL iteration 16 verification: noise-Trotter tradeoff & noise budget.

Diagnostic tests motivated by iteration 15 analysis (code confirmed correct):
  A. Noise-Trotter tradeoff: vary n_steps for noisy sim to find optimal point
  B. Noise budget decomposition: separate Hamiltonian-gate vs Lindblad-channel noise
  C. Analytical noise bound comparison
  D. Improved qubit fidelity tracking with normalized metrics

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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from gksl_math_utils import (
    compute_populations_from_density_matrix,
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
# Test A: Noise-Trotter tradeoff
# =========================================================================
def test_noise_trotter_tradeoff(params: GKSLPhysicalParameters) -> dict:
    """Vary n_steps for noisy qudit simulator to find the optimal tradeoff.

    Increasing n_steps improves Trotter accuracy O(dt^2) but increases total
    noise events (12 * n_steps for pair-only mode).

    For each n_steps, runs both noiseless and noisy simulators and computes
    the fidelity F(noiseless, noisy).
    """
    print("\n" + "=" * 70)
    print("Test A: Noise-Trotter Tradeoff (Qudit, pair-only, p_depol=0.01)")
    print("=" * 70)

    t_max = 10.0
    p_depol = 0.01
    step_counts = [5, 10, 20, 50, 100]

    # Run noiseless at finest resolution for reference
    print("  Running noiseless reference (n_steps=100)...", flush=True)
    sim_ref = QuditGKSLSimulator(params)
    r_ref = sim_ref.simulate(t_max=t_max, n_steps=100)
    rho_ref_finest = r_ref["rho_final"]

    results = []
    for ns in step_counts:
        t0 = time_module.time()

        # Noiseless
        sim_noiseless = QuditGKSLSimulator(params)
        r_noiseless = sim_noiseless.simulate(t_max=t_max, n_steps=ns)
        rho_noiseless = r_noiseless["rho_final"]

        # Noisy
        sim_noisy = QuditGKSLNoisySimulator(
            params, p_depol=p_depol, depol_pair_only=True
        )
        r_noisy = sim_noisy.simulate(t_max=t_max, n_steps=ns)
        rho_noisy = r_noisy["rho_final"]

        elapsed = time_module.time() - t0

        # Fidelities
        f_noisy_vs_noiseless = quantum_fidelity(rho_noiseless, rho_noisy)
        f_noisy_vs_finest = quantum_fidelity(rho_ref_finest, rho_noisy)
        f_noiseless_vs_finest = quantum_fidelity(rho_ref_finest, rho_noiseless)

        # Populations
        pop_noiseless = compute_populations_from_density_matrix(rho_noiseless, params)
        pop_noisy = compute_populations_from_density_matrix(rho_noisy, params)

        total_noise_events = 12 * ns  # pair-only mode

        entry = {
            "n_steps": ns,
            "dt": t_max / ns,
            "total_noise_events": total_noise_events,
            "F_noisy_vs_noiseless": f_noisy_vs_noiseless,
            "F_noisy_vs_finest_ref": f_noisy_vs_finest,
            "F_noiseless_vs_finest_ref": f_noiseless_vs_finest,
            "noiseless_N_S1": pop_noiseless["N_S1"],
            "noisy_N_S1": pop_noisy["N_S1"],
            "noisy_trace": r_noisy["trace"][-1],
            "elapsed": elapsed,
        }
        results.append(entry)

        print(
            f"  n_steps={ns:4d}  dt={t_max/ns:.4f}  "
            f"noise_events={total_noise_events:5d}  "
            f"F(noisy,noiseless)={f_noisy_vs_noiseless:.6f}  "
            f"F(noisy,finest)={f_noisy_vs_finest:.6f}  "
            f"({elapsed:.1f}s)"
        )

    return {
        "test": "noise_trotter_tradeoff",
        "t_max": t_max,
        "p_depol": p_depol,
        "depol_pair_only": True,
        "results": results,
    }


# =========================================================================
# Test B: Noise budget decomposition
# =========================================================================
def test_noise_budget(params: GKSLPhysicalParameters) -> dict:
    """Decompose noise budget: Hamiltonian gates vs Lindblad channel gates.

    Runs 1 Trotter step with:
      1. Full noise (both Hamiltonian and Lindblad gate noise)
      2. A custom decomposition separating contributions

    To separate contributions, we run step-by-step with partial noise:
      - Apply half-H with noise + Lindblad without noise + half-H with noise
        → Hamiltonian gate noise only
      - Apply half-H without noise + Lindblad with noise + half-H without noise
        → Lindblad channel noise only
    """
    print("\n" + "=" * 70)
    print("Test B: Noise Budget Decomposition (Qudit, 1 step)")
    print("=" * 70)

    t_max_step = 1.0  # Single step
    p_depol = 0.01

    # Reference: noiseless
    sim_noiseless = QuditGKSLSimulator(params)
    sim_noiseless._precompute_unitaries(t_max_step)
    rho_0 = sim_noiseless.prepare_initial_state("edge_triplet")
    rho_noiseless = sim_noiseless._trotter_step(rho_0.copy())

    # Full noise
    sim_full = QuditGKSLNoisySimulator(
        params, p_depol=p_depol, depol_pair_only=True
    )
    sim_full._precompute_unitaries(t_max_step)
    rho_full_noise = sim_full._trotter_step(rho_0.copy())

    f_full = quantum_fidelity(rho_noiseless, rho_full_noise)

    # Hamiltonian-gate noise only:
    # Manually execute the Trotter step applying noise only at Hamiltonian gates
    from qudit_gksl_noisy_simulator import (
        _apply_local_depolarization_pair,
    )
    from stinespring_utils import apply_stinespring_to_density_matrix

    d = params.d
    N = params.N_molecules

    # Half H with noise
    rho_h_only = sim_noiseless._U_H_half @ rho_0.copy() @ sim_noiseless._U_H_half.conj().T
    for i, j in params.neighbors:
        rho_h_only = _apply_local_depolarization_pair(rho_h_only, i, j, d, N, p_depol)
    # Lindblad channels WITHOUT noise
    for U_stine in sim_noiseless._U_stines:
        rho_h_only = apply_stinespring_to_density_matrix(rho_h_only, U_stine)
    # Half H with noise
    rho_h_only = sim_noiseless._U_H_half @ rho_h_only @ sim_noiseless._U_H_half.conj().T
    for i, j in params.neighbors:
        rho_h_only = _apply_local_depolarization_pair(rho_h_only, i, j, d, N, p_depol)

    f_h_only = quantum_fidelity(rho_noiseless, rho_h_only)

    # Lindblad-channel noise only:
    # Half H without noise
    rho_l_only = sim_noiseless._U_H_half @ rho_0.copy() @ sim_noiseless._U_H_half.conj().T
    # Lindblad channels WITH noise (pair channels only)
    lindblad_sites = sim_full._lindblad_sites
    for k, U_stine in enumerate(sim_noiseless._U_stines):
        rho_l_only = apply_stinespring_to_density_matrix(rho_l_only, U_stine)
        sites = lindblad_sites[k]
        if len(sites) == 2:
            rho_l_only = _apply_local_depolarization_pair(
                rho_l_only, sites[0], sites[1], d, N, p_depol
            )
    # Half H without noise
    rho_l_only = sim_noiseless._U_H_half @ rho_l_only @ sim_noiseless._U_H_half.conj().T

    f_l_only = quantum_fidelity(rho_noiseless, rho_l_only)

    # Population comparisons
    pop_noiseless = compute_populations_from_density_matrix(rho_noiseless, params)
    pop_full = compute_populations_from_density_matrix(rho_full_noise, params)
    pop_h_only = compute_populations_from_density_matrix(rho_h_only, params)
    pop_l_only = compute_populations_from_density_matrix(rho_l_only, params)

    print(f"\n  === 1-Step Fidelity vs Noiseless ===")
    print(f"  Full noise:              F = {f_full:.8f}  (6 H-gate + 6 L-gate noise events)")
    print(f"  H-gate noise only:       F = {f_h_only:.8f}  (6 H-gate noise events)")
    print(f"  L-channel noise only:    F = {f_l_only:.8f}  (6 L-gate noise events)")
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
        "N_S1_noiseless": pop_noiseless["N_S1"],
        "N_S1_full_noise": pop_full["N_S1"],
        "N_S1_hamiltonian_only": pop_h_only["N_S1"],
        "N_S1_lindblad_only": pop_l_only["N_S1"],
        "noise_events_hamiltonian": 6,  # 3 NN pairs × 2 half-H
        "noise_events_lindblad": 6,  # 3 NN pairs × 2 TTA channels
        "noise_events_total": 12,
    }


# =========================================================================
# Test C: Analytical noise bound comparison
# =========================================================================
def test_analytical_noise_bound(params: GKSLPhysicalParameters) -> dict:
    """Compare observed fidelity decay with analytical bounds.

    For depolarizing noise with effective rate ε per step:
      F(n_steps) ≈ (1 - ε)^n_steps

    We estimate ε from 1-step simulation and predict multi-step fidelity,
    then compare with actual multi-step results.
    """
    print("\n" + "=" * 70)
    print("Test C: Analytical Noise Bound Comparison (Qudit)")
    print("=" * 70)

    t_max = 10.0
    p_depol = 0.01

    # 1-step measurement for ε estimation
    sim_noiseless_1 = QuditGKSLSimulator(params)
    sim_noiseless_1._precompute_unitaries(t_max / 10)  # dt=1.0 for 10-step case
    rho_0 = sim_noiseless_1.prepare_initial_state("edge_triplet")
    rho_noiseless_1 = sim_noiseless_1._trotter_step(rho_0.copy())

    sim_noisy_1 = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    sim_noisy_1._precompute_unitaries(t_max / 10)
    rho_noisy_1 = sim_noisy_1._trotter_step(rho_0.copy())

    f_1step = quantum_fidelity(rho_noiseless_1, rho_noisy_1)
    epsilon_measured = 1.0 - f_1step

    print(f"  1-step fidelity: F = {f_1step:.8f}")
    print(f"  Measured ε = 1 - F = {epsilon_measured:.8f}")

    # Multi-step: predict vs actual
    print(f"\n  === Multi-Step Fidelity: Predicted vs Actual ===")
    print(f"  {'n_steps':>7} | {'F_predicted':>12} | {'F_actual':>12} | {'Δ':>10}")
    print("  " + "-" * 50)

    step_counts = [5, 10, 20, 50]
    step_comparison = []

    for ns in step_counts:
        # Predicted using ε from 1-step
        f_predicted = (1.0 - epsilon_measured) ** ns

        # Actual
        dt = t_max / ns
        sim_noiseless = QuditGKSLSimulator(params)
        sim_noiseless._precompute_unitaries(dt)
        rho_ref = rho_0.copy()
        for _ in range(ns):
            rho_ref = sim_noiseless._trotter_step(rho_ref)

        sim_noisy = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
        sim_noisy._precompute_unitaries(dt)
        rho_noisy = rho_0.copy()
        for _ in range(ns):
            rho_noisy = sim_noisy._trotter_step(rho_noisy)

        f_actual = quantum_fidelity(rho_ref, rho_noisy)

        delta = f_actual - f_predicted
        step_comparison.append({
            "n_steps": ns,
            "dt": dt,
            "F_predicted": f_predicted,
            "F_actual": f_actual,
            "delta": delta,
        })
        print(f"  {ns:7d} | {f_predicted:12.8f} | {f_actual:12.8f} | {delta:+10.6f}")

    # Theoretical bound: per-event infidelity for d=3 pair depol
    # For pair depol with p_depol on d=3 pair: effective infidelity ≈ p * (1 - 1/d^2)
    d = params.d
    infidelity_per_event = p_depol * (1.0 - 1.0 / (d * d))
    # With 12 events per step (pair-only, dt=1.0)
    n_events_per_step = 12
    epsilon_theoretical = 1.0 - (1.0 - infidelity_per_event) ** n_events_per_step

    print(f"\n  === Theoretical Estimates ===")
    print(f"  Per-event infidelity (d={d}): p*(1-1/d²) = {infidelity_per_event:.6f}")
    print(f"  Events per step: {n_events_per_step}")
    print(f"  Theoretical ε per step: {epsilon_theoretical:.6f}")
    print(f"  Measured ε per step: {epsilon_measured:.6f}")

    return {
        "test": "analytical_noise_bound",
        "t_max": t_max,
        "p_depol": p_depol,
        "depol_pair_only": True,
        "F_1step": f_1step,
        "epsilon_measured": epsilon_measured,
        "epsilon_theoretical": epsilon_theoretical,
        "infidelity_per_event": infidelity_per_event,
        "events_per_step": n_events_per_step,
        "step_comparison": step_comparison,
    }


# =========================================================================
# Test D: Improved qubit fidelity tracking
# =========================================================================
def test_qubit_improved_fidelity(params: GKSLPhysicalParameters) -> dict:
    """Improved qubit fidelity tracking with both raw and normalized metrics.

    Runs DM noisy qubit and qudit simulators for 10 steps with detailed
    per-step reporting of:
      - F_raw: raw fidelity (affected by trace mismatch)
      - F_norm: normalized fidelity (physical subspace structure)
      - Trace, forbidden population
      - Conditional N_S1 (normalized by trace)
    """
    print("\n" + "=" * 70)
    print("Test D: Improved Qubit Fidelity Tracking (pair-only, p_depol=0.01)")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    dt = t_max / n_steps
    p_depol = 0.01

    # Noiseless references (step-by-step)
    sim_qd_ref = QuditGKSLSimulator(params)
    sim_qd_ref._precompute_unitaries(dt)
    rho_qd_ref = sim_qd_ref.prepare_initial_state("edge_triplet")

    sim_qb_ref = QubitGKSLSimulator(params)
    sim_qb_ref._precompute_unitaries(dt)
    rho_qb_ref = sim_qb_ref.prepare_initial_state("edge_triplet")

    # Noisy simulators
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

        # Qudit fidelity
        f_qd = quantum_fidelity(rho_qd_ref, rho_qd_noisy)

        # Qubit: extract qutrit-space
        rho_qb_ref_qt = sim_qb_ref._extract_from_qubit_space(rho_qb_ref)
        rho_qb_noisy_qt = sim_qb_noisy._extract_from_qubit_space(rho_qb_noisy)

        # Raw and normalised fidelity
        f_qb_raw = quantum_fidelity(rho_qb_ref_qt, rho_qb_noisy_qt)
        f_qb_norm = quantum_fidelity(
            normalize_rho(rho_qb_ref_qt), normalize_rho(rho_qb_noisy_qt)
        )

        tr_qb = float(np.real(np.trace(rho_qb_noisy_qt)))
        p_forb = 1.0 - tr_qb

        # Populations
        pop_qd_noisy = compute_populations_from_density_matrix(rho_qd_noisy, params)
        pop_qb_noisy = compute_populations_from_density_matrix(rho_qb_noisy_qt, params)

        # Conditional N_S1 (normalised by trace)
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
    print("TTA-UC GKSL Iteration 16: Noise-Trotter Tradeoff & Noise Budget")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 16,
        "purpose": "Noise-Trotter tradeoff, noise budget decomposition, "
                   "analytical bounds, improved qubit fidelity",
        "tests": [],
    }

    # Test A: Noise-Trotter tradeoff
    ta = test_noise_trotter_tradeoff(params)
    results["tests"].append(ta)

    # Test B: Noise budget decomposition
    tb = test_noise_budget(params)
    results["tests"].append(tb)

    # Test C: Analytical noise bound comparison
    tc = test_analytical_noise_bound(params)
    results["tests"].append(tc)

    # Test D: Improved qubit fidelity tracking
    td = test_qubit_improved_fidelity(params)
    results["tests"].append(td)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration16_noise_tradeoff_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 16: ノイズ-Trotterトレードオフ & ノイズバジェット",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "",
    ]

    # Test A summary
    md_lines.extend([
        "## Test A: ノイズ-Trotterトレードオフ（Qudit, pair-only, p_depol=0.01）",
        "",
        f"- t_max: {ta['t_max']}, p_depol: {ta['p_depol']}",
        "",
        "| n_steps | dt | noise_events | F(noisy, noiseless) | F(noisy, finest_ref) |",
        "|---------|------|-------|---------------------|----------------------|",
    ])
    for r in ta["results"]:
        md_lines.append(
            f"| {r['n_steps']} | {r['dt']:.4f} | {r['total_noise_events']} | "
            f"{r['F_noisy_vs_noiseless']:.8f} | {r['F_noisy_vs_finest_ref']:.8f} |"
        )
    md_lines.append("")

    # Test B summary
    md_lines.extend([
        "## Test B: ノイズバジェット分解（1ステップ、dt=1.0）",
        "",
        f"- p_depol: {tb['p_depol']}",
        "",
        "| 条件 | ノイズイベント数 | F(vs noiseless) | N_S1 |",
        "|------|-----------------|-----------------|------|",
        f"| ノイズなし | 0 | 1.00000000 | {tb['N_S1_noiseless']:.8f} |",
        f"| 全ノイズ | {tb['noise_events_total']} | "
        f"{tb['F_full_noise']:.8f} | {tb['N_S1_full_noise']:.8f} |",
        f"| Hゲートノイズのみ | {tb['noise_events_hamiltonian']} | "
        f"{tb['F_hamiltonian_noise_only']:.8f} | {tb['N_S1_hamiltonian_only']:.8f} |",
        f"| Lチャネルノイズのみ | {tb['noise_events_lindblad']} | "
        f"{tb['F_lindblad_noise_only']:.8f} | {tb['N_S1_lindblad_only']:.8f} |",
        "",
    ])

    # Test C summary
    md_lines.extend([
        "## Test C: 解析的ノイズ境界との比較",
        "",
        f"- 1ステップ忠実度: F = {tc['F_1step']:.8f}",
        f"- 実測 ε: {tc['epsilon_measured']:.8f}",
        f"- 理論 ε: {tc['epsilon_theoretical']:.8f}",
        "",
        "### 多ステップ予測 vs 実測",
        "",
        "| n_steps | F_predicted | F_actual | Δ |",
        "|---------|-------------|----------|---|",
    ])
    for sc in tc["step_comparison"]:
        md_lines.append(
            f"| {sc['n_steps']} | {sc['F_predicted']:.8f} | "
            f"{sc['F_actual']:.8f} | {sc['delta']:+.6f} |"
        )
    md_lines.append("")

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

    md_path = output_dir / f"iteration16_noise_tradeoff_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    # Test A: Find optimal n_steps
    if ta["results"]:
        best = max(ta["results"], key=lambda r: r["F_noisy_vs_noiseless"])
        print(f"  Test A: Best F(noisy,noiseless) = {best['F_noisy_vs_noiseless']:.6f} "
              f"at n_steps={best['n_steps']}")

    # Test B
    print(f"  Test B: Full noise F = {tb['F_full_noise']:.6f}, "
          f"H-only F = {tb['F_hamiltonian_noise_only']:.6f}, "
          f"L-only F = {tb['F_lindblad_noise_only']:.6f}")

    # Test C
    print(f"  Test C: ε_measured = {tc['epsilon_measured']:.6f}, "
          f"ε_theoretical = {tc['epsilon_theoretical']:.6f}")

    # Test D
    if td["step_data"]:
        final = td["step_data"][-1]
        print(f"  Test D: F_qd={final['F_qudit']:.6f}, "
              f"F_qb_raw={final['F_qubit_raw']:.7f}, "
              f"F_qb_norm={final['F_qubit_norm']:.8f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
