#!/usr/bin/env python3
"""TTA-UC GKSL iteration 15 verification: Trotter convergence & DM-noisy consistency.

Diagnostic tests motivated by iteration 14 analysis:
  1. Trotter convergence (noiseless): vary n_steps to quantify Trotter error
  2. Exact GKSL reference vs Trotter (qudit): compare with full Liouvillian expm
  3. DM noisy multi-step fidelity with per-step tracking
  4. DM noisy vs shot noisy consistency check
  5. Forbidden-state leakage analysis per step

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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_populations_from_density_matrix,
    vectorize_density_matrix,
    unvectorize_density_matrix,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator
from qubit_gksl_shot_simulator import QubitGKSLNoisyShotSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator
from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator
from qudit_gksl_simulator import QuditGKSLSimulator
from stinespring_utils import build_gksl_superoperator


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
# Test 1: Trotter convergence
# =========================================================================
def test_trotter_convergence(params: GKSLPhysicalParameters) -> dict:
    """Compare noiseless qudit results at different Trotter step counts.

    Runs the same total time t_max=10 with n_steps=10,20,50,100
    and checks how the populations converge.
    """
    print("\n" + "=" * 70)
    print("Test 1: Trotter Convergence (Noiseless Qudit)")
    print("=" * 70)

    t_max = 10.0
    step_counts = [10, 20, 50, 100]

    results_by_steps: dict[int, dict] = {}
    rho_finals: dict[int, np.ndarray] = {}

    for ns in step_counts:
        t0 = time_module.time()
        sim = QuditGKSLSimulator(params)
        r = sim.simulate(t_max=t_max, n_steps=ns)
        elapsed = time_module.time() - t0
        pop_final = r["populations"][-1]
        rho_finals[ns] = r["rho_final"]
        results_by_steps[ns] = {
            "n_steps": ns,
            "dt": t_max / ns,
            "N_S0": pop_final["N_S0"],
            "N_T1": pop_final["N_T1"],
            "N_S1": pop_final["N_S1"],
            "trace": r["trace"][-1],
            "elapsed": elapsed,
        }
        print(f"  n_steps={ns:4d}  dt={t_max/ns:.4f}  "
              f"N_S0={pop_final['N_S0']:.8f}  N_T1={pop_final['N_T1']:.8f}  "
              f"N_S1={pop_final['N_S1']:.8f}  Tr={r['trace'][-1]:.10f}  "
              f"({elapsed:.2f}s)")

    # Fidelity between consecutive refinements
    print("\n  === Fidelity Between Refinements ===")
    fidelity_pairs = []
    for i in range(len(step_counts) - 1):
        ns_a = step_counts[i]
        ns_b = step_counts[i + 1]
        f = quantum_fidelity(rho_finals[ns_a], rho_finals[ns_b])
        fidelity_pairs.append({
            "n_steps_a": ns_a, "n_steps_b": ns_b, "fidelity": f,
        })
        print(f"  F(n={ns_a}, n={ns_b}) = {f:.10f}")

    # Fidelity against finest (n_steps=100)
    ns_ref = step_counts[-1]
    print(f"\n  === Fidelity vs Finest (n_steps={ns_ref}) ===")
    fidelity_vs_finest = []
    for ns in step_counts[:-1]:
        f = quantum_fidelity(rho_finals[ns], rho_finals[ns_ref])
        fidelity_vs_finest.append({"n_steps": ns, "fidelity": f})
        print(f"  F(n={ns}, n={ns_ref}) = {f:.10f}")

    return {
        "test": "trotter_convergence",
        "t_max": t_max,
        "results_by_steps": results_by_steps,
        "fidelity_pairs": fidelity_pairs,
        "fidelity_vs_finest": fidelity_vs_finest,
    }


# =========================================================================
# Test 2: Exact GKSL vs Trotter
# =========================================================================
def test_exact_gksl_vs_trotter(params: GKSLPhysicalParameters) -> dict:
    """Compare exact GKSL evolution (full Liouvillian expm) with Trotter.

    Computes exact evolution for the qudit (81-dim) system using:
      rho(t) = expm(L_total * t) @ vec(rho_0)
    and compares with the Trotter simulation at each time step.

    Uses manual step-by-step Trotter evolution to capture intermediate rhos.
    """
    print("\n" + "=" * 70)
    print("Test 2: Exact GKSL vs Trotter (Qudit)")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    dt = t_max / n_steps

    # Build exact Liouvillian
    H_total = build_onsite_hamiltonian(params) + build_transfer_hamiltonian(params)
    lindblad_ops = build_lindblad_operators(params)
    L_super = build_gksl_superoperator(H_total, lindblad_ops)

    # Prepare initial state and Trotter simulator
    sim = QuditGKSLSimulator(params)
    rho_0 = sim.prepare_initial_state("edge_triplet")
    vec_0 = vectorize_density_matrix(rho_0)

    # Exact evolution at each time step
    dim = rho_0.shape[0]
    print(f"  Liouvillian dimension: {L_super.shape[0]}x{L_super.shape[1]}")
    print(f"  Computing exact evolutions...", flush=True)

    exact_rhos = [rho_0]
    t0 = time_module.time()
    for s in range(1, n_steps + 1):
        t = s * dt
        prop = expm(L_super * t)
        vec_t = prop @ vec_0
        rho_t = unvectorize_density_matrix(vec_t, dim)
        exact_rhos.append(rho_t)
    exact_time = time_module.time() - t0
    print(f"  Exact computation: {exact_time:.2f}s")

    # Trotter evolution: step-by-step to capture intermediate rhos
    print(f"  Computing Trotter evolution (step-by-step)...", flush=True)
    t0 = time_module.time()
    sim._precompute_unitaries(dt)
    rho_tr = rho_0.copy()
    trotter_rhos = [rho_0]
    for _step in range(n_steps):
        rho_tr = sim._trotter_step(rho_tr)
        trotter_rhos.append(rho_tr.copy())
    trotter_time = time_module.time() - t0
    print(f"  Trotter computation: {trotter_time:.2f}s")

    # Per-step comparison
    print(f"\n  === Per-Step Comparison (Exact vs Trotter, n_steps={n_steps}) ===")
    print(f"  {'Step':>4} | {'t':>5} | {'F(exact,trotter)':>18} | "
          f"{'||Δ||_F':>10} | {'Δ(N_S1)':>12}")
    print("  " + "-" * 65)

    step_data = []
    for s in range(n_steps + 1):
        rho_ex = exact_rhos[s]
        rho_trot = trotter_rhos[s]
        f = quantum_fidelity(rho_ex, rho_trot)
        frob = float(np.linalg.norm(rho_ex - rho_trot, "fro"))
        pop_ex = compute_populations_from_density_matrix(rho_ex, params)
        pop_tr = compute_populations_from_density_matrix(rho_trot, params)
        d_ns1 = pop_tr["N_S1"] - pop_ex["N_S1"]
        tr_exact = float(np.real(np.trace(rho_ex)))

        sd = {
            "step": s,
            "time": s * dt,
            "fidelity": f,
            "frobenius": frob,
            "exact_N_S0": pop_ex["N_S0"],
            "exact_N_T1": pop_ex["N_T1"],
            "exact_N_S1": pop_ex["N_S1"],
            "trotter_N_S0": pop_tr["N_S0"],
            "trotter_N_T1": pop_tr["N_T1"],
            "trotter_N_S1": pop_tr["N_S1"],
            "delta_N_S1": d_ns1,
            "exact_trace": tr_exact,
        }
        step_data.append(sd)
        print(f"  {s:4d} | {s*dt:5.1f} | {f:18.12f} | {frob:10.6f} | {d_ns1:+12.8f}")

    return {
        "test": "exact_gksl_vs_trotter",
        "t_max": t_max,
        "n_steps": n_steps,
        "dt": dt,
        "step_data": step_data,
        "exact_computation_time": exact_time,
        "trotter_computation_time": trotter_time,
    }


# =========================================================================
# Test 3: DM noisy multi-step fidelity
# =========================================================================
def test_dm_noisy_multistep(params: GKSLPhysicalParameters) -> dict:
    """DM noisy simulator multi-step fidelity tracking (pair-only).

    Runs DM noisy simulators for 10 steps and records per-step fidelity,
    populations, trace, and forbidden-state population (qubit).
    No shot noise — this isolates hardware noise effects.

    Uses manual step-by-step evolution to capture intermediate density matrices.
    """
    print("\n" + "=" * 70)
    print("Test 3: DM Noisy Multi-Step Fidelity (Pair-Only)")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    dt = t_max / n_steps
    p_depol = 0.01

    # --- Noiseless references (step-by-step) ---
    print("\n  Running noiseless references (step-by-step)...", flush=True)
    sim_qd_ref = QuditGKSLSimulator(params)
    sim_qd_ref._precompute_unitaries(dt)
    rho_qd_ref = sim_qd_ref.prepare_initial_state("edge_triplet")

    sim_qb_ref = QubitGKSLSimulator(params)
    sim_qb_ref._precompute_unitaries(dt)
    rho_qb_ref = sim_qb_ref.prepare_initial_state("edge_triplet")

    # --- Noisy simulators (step-by-step) ---
    print("  Running DM noisy (step-by-step)...", flush=True)
    sim_qd_noisy = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    sim_qd_noisy._precompute_unitaries(dt)
    rho_qd_noisy = sim_qd_noisy.prepare_initial_state("edge_triplet")

    sim_qb_noisy = QubitGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    sim_qb_noisy._precompute_unitaries(dt)
    rho_qb_noisy = sim_qb_noisy.prepare_initial_state("edge_triplet")

    step_data = []

    for s in range(n_steps + 1):
        if s > 0:
            rho_qd_ref = sim_qd_ref._trotter_step(rho_qd_ref)
            rho_qb_ref = sim_qb_ref._trotter_step(rho_qb_ref)
            rho_qd_noisy = sim_qd_noisy._trotter_step(rho_qd_noisy)
            rho_qb_noisy = sim_qb_noisy._trotter_step(rho_qb_noisy)

        # Qudit fidelity (both in 81-dim space)
        f_qd = quantum_fidelity(rho_qd_ref, rho_qd_noisy)
        frob_qd = float(np.linalg.norm(rho_qd_ref - rho_qd_noisy, "fro"))

        # Qubit: extract qutrit-space density matrices for comparison
        rho_qb_ref_qt = sim_qb_ref._extract_from_qubit_space(rho_qb_ref)
        rho_qb_noisy_qt = sim_qb_noisy._extract_from_qubit_space(rho_qb_noisy)
        f_qb = quantum_fidelity(rho_qb_ref_qt, rho_qb_noisy_qt)
        frob_qb = float(np.linalg.norm(rho_qb_ref_qt - rho_qb_noisy_qt, "fro"))

        tr_qb = float(np.real(np.trace(rho_qb_noisy_qt)))
        p_forb = 1.0 - tr_qb

        pop_qd_ref = compute_populations_from_density_matrix(rho_qd_ref, params)
        pop_qd_noisy = compute_populations_from_density_matrix(rho_qd_noisy, params)
        pop_qb_noisy = compute_populations_from_density_matrix(rho_qb_noisy_qt, params)

        sd = {
            "step": s,
            "time": s * dt,
            "F_qudit": f_qd,
            "F_qubit": f_qb,
            "frob_qudit": frob_qd,
            "frob_qubit": frob_qb,
            "trace_qubit": tr_qb,
            "forbidden_pop_qubit": p_forb,
            "populations_qudit_ref": {k: pop_qd_ref[k] for k in ("N_S0", "N_T1", "N_S1")},
            "populations_qudit_noisy": {k: pop_qd_noisy[k] for k in ("N_S0", "N_T1", "N_S1")},
            "populations_qubit_noisy": {k: pop_qb_noisy[k] for k in ("N_S0", "N_T1", "N_S1")},
        }
        step_data.append(sd)

        print(f"  Step {s:2d} | t={s*dt:5.1f} | "
              f"F_qd={f_qd:.6f} | F_qb={f_qb:.6f} | "
              f"Tr_qb={tr_qb:.8f} | P_forb={p_forb:.8f}")

    return {
        "test": "dm_noisy_multistep",
        "t_max": t_max,
        "n_steps": n_steps,
        "p_depol": p_depol,
        "depol_pair_only": True,
        "step_data": step_data,
    }


# =========================================================================
# Test 4: DM noisy vs shot noisy consistency
# =========================================================================
def test_dm_vs_shot_consistency(params: GKSLPhysicalParameters) -> dict:
    """Compare DM noisy and shot noisy results for pair-only mode.

    Both should produce the same physics (up to statistical error in shots).
    """
    print("\n" + "=" * 70)
    print("Test 4: DM Noisy vs Shot Noisy Consistency (Pair-Only)")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    p_depol = 0.01
    n_shots = 2000
    seed = 42

    # --- DM noisy ---
    print("\n  Running DM noisy (qudit)...", flush=True)
    sim_qd_dm = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_qd_dm = sim_qd_dm.simulate(t_max=t_max, n_steps=n_steps)

    print("  Running DM noisy (qubit)...", flush=True)
    sim_qb_dm = QubitGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_qb_dm = sim_qb_dm.simulate(t_max=t_max, n_steps=n_steps)

    # --- Shot noisy ---
    print(f"  Running shot noisy (qudit, {n_shots} shots)...", flush=True)
    sim_qd_shot = QuditGKSLNoisyShotSimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_qd_shot = sim_qd_shot.simulate(
        t_max=t_max, n_steps=n_steps, n_shots=n_shots, seed=seed,
    )

    print(f"  Running shot noisy (qubit, {n_shots} shots)...", flush=True)
    sim_qb_shot = QubitGKSLNoisyShotSimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_qb_shot = sim_qb_shot.simulate(
        t_max=t_max, n_steps=n_steps, n_shots=n_shots, seed=seed,
    )

    # --- Final state fidelity ---
    rho_qd_dm_f = r_qd_dm["rho_final"]
    rho_qd_shot_f = r_qd_shot["rho_final"]
    rho_qb_dm_f = r_qb_dm["rho_final"]
    rho_qb_shot_f = r_qb_shot["rho_final"]

    f_qd = quantum_fidelity(rho_qd_dm_f, rho_qd_shot_f)
    f_qb = quantum_fidelity(rho_qb_dm_f, rho_qb_shot_f)
    f_qb_norm = quantum_fidelity(normalize_rho(rho_qb_dm_f), normalize_rho(rho_qb_shot_f))

    frob_qd = float(np.linalg.norm(rho_qd_dm_f - rho_qd_shot_f, "fro"))
    frob_qb = float(np.linalg.norm(rho_qb_dm_f - rho_qb_shot_f, "fro"))

    print(f"\n  === Final State Comparison (DM vs Shot, 10 steps) ===")
    print(f"  Qudit: F(DM, shot) = {f_qd:.6f}  ||Δ||_F = {frob_qd:.6f}")
    print(f"  Qubit: F(DM, shot) = {f_qb:.6f}  ||Δ||_F = {frob_qb:.6f}")
    print(f"  Qubit: F_norm(DM, shot) = {f_qb_norm:.6f}")

    # --- Per-step population comparison ---
    print(f"\n  === Per-Step Population Comparison ===")
    print(f"  {'Step':>4} | {'Qudit DM N_S1':>14} | {'Qudit Shot N_S1':>15} | "
          f"{'Δ':>10} | {'Qubit DM Tr':>12} | {'Qubit Shot Tr':>14} | {'Δ':>10}")
    print("  " + "-" * 95)

    step_comparison = []
    for s in range(n_steps + 1):
        pop_qd_dm = r_qd_dm["populations"][s]
        pop_qd_shot = r_qd_shot["populations"][s]
        pop_qb_dm = r_qb_dm["populations"][s]
        pop_qb_shot = r_qb_shot["populations"][s]
        tr_qb_dm = r_qb_dm["trace"][s]
        tr_qb_shot = r_qb_shot["trace"][s]

        sc = {
            "step": s,
            "time": r_qd_dm["times"][s],
            "qudit_dm_N_S1": pop_qd_dm["N_S1"],
            "qudit_shot_N_S1": pop_qd_shot["N_S1"],
            "delta_qudit_N_S1": pop_qd_shot["N_S1"] - pop_qd_dm["N_S1"],
            "qubit_dm_trace": tr_qb_dm,
            "qubit_shot_trace": tr_qb_shot,
            "delta_qubit_trace": tr_qb_shot - tr_qb_dm,
            "qudit_dm_N_S0": pop_qd_dm["N_S0"],
            "qudit_shot_N_S0": pop_qd_shot["N_S0"],
            "qubit_dm_N_S1": pop_qb_dm["N_S1"],
            "qubit_shot_N_S1": pop_qb_shot["N_S1"],
        }
        step_comparison.append(sc)

        print(f"  {s:4d} | {pop_qd_dm['N_S1']:14.8f} | {pop_qd_shot['N_S1']:15.8f} | "
              f"{pop_qd_shot['N_S1']-pop_qd_dm['N_S1']:+10.6f} | "
              f"{tr_qb_dm:12.8f} | {tr_qb_shot:14.8f} | "
              f"{tr_qb_shot-tr_qb_dm:+10.6f}")

    return {
        "test": "dm_vs_shot_consistency",
        "t_max": t_max,
        "n_steps": n_steps,
        "p_depol": p_depol,
        "n_shots": n_shots,
        "seed": seed,
        "depol_pair_only": True,
        "final_fidelity": {
            "qudit": f_qd,
            "qubit_raw": f_qb,
            "qubit_norm": f_qb_norm,
        },
        "final_frobenius": {
            "qudit": frob_qd,
            "qubit": frob_qb,
        },
        "step_comparison": step_comparison,
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
    print("TTA-UC GKSL Iteration 15: Trotter Convergence & DM-Noisy Consistency")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 15,
        "purpose": "Trotter convergence, exact GKSL reference, DM-noisy consistency",
        "tests": [],
    }

    # Test 1: Trotter convergence
    t1 = test_trotter_convergence(params)
    results["tests"].append(t1)

    # Test 2: Exact GKSL vs Trotter
    t2 = test_exact_gksl_vs_trotter(params)
    results["tests"].append(t2)

    # Test 3: DM noisy multi-step
    t3 = test_dm_noisy_multistep(params)
    results["tests"].append(t3)

    # Test 4: DM vs shot consistency
    t4 = test_dm_vs_shot_consistency(params)
    results["tests"].append(t4)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration15_trotter_consistency_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 15: Trotter収束性 & DMノイジー整合性検証",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "",
    ]

    # Test 1 summary
    md_lines.extend([
        "## Test 1: Trotter収束テスト（ノイズレス Qudit）",
        "",
        f"- t_max: {t1['t_max']}",
        "",
        "### 最終状態 populations (t=10)",
        "",
        "| n_steps | dt | N_S0 | N_T1 | N_S1 | Trace |",
        "|---------|-----|------|------|------|-------|",
    ])
    for ns_key, r in t1["results_by_steps"].items():
        md_lines.append(
            f"| {r['n_steps']} | {r['dt']:.4f} | {r['N_S0']:.8f} | "
            f"{r['N_T1']:.8f} | {r['N_S1']:.8f} | {r['trace']:.10f} |"
        )
    md_lines.append("")

    if t1["fidelity_vs_finest"]:
        md_lines.extend([
            f"### 最精細(n_steps={t1['results_by_steps'][max(t1['results_by_steps'])]['n_steps']})との忠実度",
            "",
        ])
        for fp in t1["fidelity_vs_finest"]:
            md_lines.append(f"- F(n={fp['n_steps']}, n=100) = {fp['fidelity']:.10f}")
        md_lines.append("")

    # Test 2 summary
    md_lines.extend([
        "## Test 2: 厳密GKSL解 vs Trotter (Qudit)",
        "",
        f"- n_steps: {t2['n_steps']}, dt: {t2['dt']}",
        "",
        "### ステップ毎の比較",
        "",
        "| Step | Time | F(exact, trotter) | ||Δ||_F | Exact N_S1 | Trotter N_S1 | Δ(N_S1) |",
        "|------|------|-------------------|---------|------------|--------------|---------|",
    ])
    for sd in t2["step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['time']:.1f} | {sd['fidelity']:.10f} | "
            f"{sd['frobenius']:.6f} | {sd['exact_N_S1']:.8f} | "
            f"{sd['trotter_N_S1']:.8f} | {sd['delta_N_S1']:+.8f} |"
        )
    md_lines.append("")

    # Test 3 summary
    md_lines.extend([
        "## Test 3: DMノイジー10ステップ忠実度（pair-only, p_depol=0.01）",
        "",
        "### ステップ毎の忠実度",
        "",
        "| Step | Time | F_qudit | F_qubit | Tr_qubit | P_forbidden |",
        "|------|------|---------|---------|----------|-------------|",
    ])
    for sd in t3["step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['time']:.1f} | {sd['F_qudit']:.6f} | "
            f"{sd['F_qubit']:.6f} | {sd['trace_qubit']:.8f} | "
            f"{sd['forbidden_pop_qubit']:.8f} |"
        )
    md_lines.extend(["", "### ステップ毎の populations", ""])
    md_lines.extend([
        "| Step | Ref N_S1 | Qudit noisy N_S1 | Qubit noisy N_S1 |",
        "|------|----------|------------------|------------------|",
    ])
    for sd in t3["step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['populations_qudit_ref']['N_S1']:.8f} | "
            f"{sd['populations_qudit_noisy']['N_S1']:.8f} | "
            f"{sd['populations_qubit_noisy']['N_S1']:.8f} |"
        )
    md_lines.append("")

    # Test 4 summary
    md_lines.extend([
        "## Test 4: DMノイジー vs ショットノイジー整合性（pair-only）",
        "",
        f"- n_shots: {t4['n_shots']}, seed: {t4['seed']}",
        "",
        "### 最終状態忠実度",
        "",
        f"- Qudit F(DM, shot) = {t4['final_fidelity']['qudit']:.6f}  ||Δ||_F = {t4['final_frobenius']['qudit']:.6f}",
        f"- Qubit F(DM, shot) = {t4['final_fidelity']['qubit_raw']:.6f}  ||Δ||_F = {t4['final_frobenius']['qubit']:.6f}",
        f"- Qubit F_norm(DM, shot) = {t4['final_fidelity']['qubit_norm']:.6f}",
        "",
        "### ステップ毎の N_S1 比較",
        "",
        "| Step | Qudit DM | Qudit Shot | Δ | Qubit DM Tr | Qubit Shot Tr | Δ |",
        "|------|----------|------------|---|-------------|---------------|---|",
    ])
    for sc in t4["step_comparison"]:
        md_lines.append(
            f"| {sc['step']} | {sc['qudit_dm_N_S1']:.8f} | {sc['qudit_shot_N_S1']:.8f} | "
            f"{sc['delta_qudit_N_S1']:+.6f} | {sc['qubit_dm_trace']:.8f} | "
            f"{sc['qubit_shot_trace']:.8f} | {sc['delta_qubit_trace']:+.6f} |"
        )
    md_lines.append("")

    md_path = output_dir / f"iteration15_trotter_consistency_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    if t1["fidelity_vs_finest"]:
        print(f"  Test 1: F(n=10, n=100) = {t1['fidelity_vs_finest'][0]['fidelity']:.10f}")
    if t2["step_data"]:
        final_t2 = t2["step_data"][-1]
        print(f"  Test 2: F(exact, trotter) at t=10 = {final_t2['fidelity']:.10f}")
    if t3["step_data"]:
        final_t3 = [s for s in t3["step_data"] if s["step"] == t3["n_steps"]]
        if final_t3:
            print(f"  Test 3: F_qudit={final_t3[0]['F_qudit']:.6f}, "
                  f"F_qubit={final_t3[0]['F_qubit']:.6f}, "
                  f"P_forb={final_t3[0]['forbidden_pop_qubit']:.6f}")
    print(f"  Test 4: F(DM,shot) qudit={t4['final_fidelity']['qudit']:.6f}, "
          f"qubit={t4['final_fidelity']['qubit_raw']:.6f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
