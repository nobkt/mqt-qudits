#!/usr/bin/env python3
"""TTA-UC GKSL iteration 14 verification: pair-only depolarization comparison.

This script verifies the depol_pair_only feature and compares qubit vs qudit
encodings under pair-only depolarization noise, where noise is applied only
to 2-qubit/2-qudit (or higher) interactions.

Physical motivation:
  - Single-qubit/qudit gates have much lower error rates on real hardware
  - Applying noise only to multi-body interactions is a more realistic model
  - This allows a fair comparison of qubit vs qudit noise resilience

Verification tests:
  1. DM noisy simulator: pair-only vs all-gates noise comparison (1 step)
  2. Shot-based: pair-only noise with per-step fidelity tracking (10 steps)
  3. Error budget comparison: all-gates vs pair-only
  4. Qubit vs Qudit fair comparison under pair-only noise

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

from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator
from qubit_gksl_shot_simulator import QubitGKSLNoisyShotSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator
from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator
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


def test_dm_pair_only_comparison(params: GKSLPhysicalParameters) -> dict:
    """Test 1: DM noisy simulator pair-only vs all-gates noise (1 step).

    Compares density matrix evolution with noise on all gates vs only pair gates
    to verify the depol_pair_only feature works correctly.
    """
    print("\n" + "=" * 70)
    print("Test 1: DM Noisy Pair-Only vs All-Gates (1 Trotter Step)")
    print("=" * 70)

    p_depol = 0.01
    t_max = 1.0
    n_steps = 1

    # --- Qudit ---
    print("\n  --- Qudit DM ---")
    sim_ref = QuditGKSLSimulator(params)
    r_ref = sim_ref.simulate(t_max=t_max, n_steps=n_steps)

    sim_all = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=False)
    r_all = sim_all.simulate(t_max=t_max, n_steps=n_steps)

    sim_pair = QuditGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_pair = sim_pair.simulate(t_max=t_max, n_steps=n_steps)

    rho_ref = r_ref["rho_final"]
    f_all = quantum_fidelity(rho_ref, r_all["rho_final"])
    f_pair = quantum_fidelity(rho_ref, r_pair["rho_final"])
    frob_all = float(np.linalg.norm(rho_ref - r_all["rho_final"], "fro"))
    frob_pair = float(np.linalg.norm(rho_ref - r_pair["rho_final"], "fro"))

    print(f"  Qudit F(noiseless, all-gates):   {f_all:.6f}  ||Δ||_F = {frob_all:.6f}")
    print(f"  Qudit F(noiseless, pair-only):    {f_pair:.6f}  ||Δ||_F = {frob_pair:.6f}")
    print(f"  Qudit improvement (pair vs all):  F: +{f_pair - f_all:.6f}")

    # --- Qubit ---
    print("\n  --- Qubit DM ---")
    sim_ref_qb = QubitGKSLSimulator(params)
    r_ref_qb = sim_ref_qb.simulate(t_max=t_max, n_steps=n_steps)

    sim_all_qb = QubitGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=False)
    r_all_qb = sim_all_qb.simulate(t_max=t_max, n_steps=n_steps)

    sim_pair_qb = QubitGKSLNoisySimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_pair_qb = sim_pair_qb.simulate(t_max=t_max, n_steps=n_steps)

    rho_ref_qb = r_ref_qb["rho_final"]
    f_all_qb = quantum_fidelity(rho_ref_qb, r_all_qb["rho_final"])
    f_pair_qb = quantum_fidelity(rho_ref_qb, r_pair_qb["rho_final"])
    frob_all_qb = float(np.linalg.norm(rho_ref_qb - r_all_qb["rho_final"], "fro"))
    frob_pair_qb = float(np.linalg.norm(rho_ref_qb - r_pair_qb["rho_final"], "fro"))

    # Trace values (qubit noisy trace decreases due to forbidden-state leakage)
    tr_all_qb = r_all_qb["trace"][-1]
    tr_pair_qb = r_pair_qb["trace"][-1]

    print(f"  Qubit F(noiseless, all-gates):   {f_all_qb:.6f}  ||Δ||_F = {frob_all_qb:.6f}  Tr = {tr_all_qb:.8f}")
    print(f"  Qubit F(noiseless, pair-only):    {f_pair_qb:.6f}  ||Δ||_F = {frob_pair_qb:.6f}  Tr = {tr_pair_qb:.8f}")
    print(f"  Qubit improvement (pair vs all):  F: +{f_pair_qb - f_all_qb:.6f}")

    # --- Cross comparison ---
    f_cross_pair = quantum_fidelity(r_pair["rho_final"], r_pair_qb["rho_final"])
    f_cross_pair_norm = quantum_fidelity(
        normalize_rho(r_pair["rho_final"]), normalize_rho(r_pair_qb["rho_final"])
    )
    print(f"\n  Cross (qudit pair vs qubit pair): F_raw = {f_cross_pair:.6f}, F_norm = {f_cross_pair_norm:.6f}")

    return {
        "test": "dm_pair_only_comparison",
        "p_depol": p_depol,
        "qudit": {
            "F_all_gates": f_all,
            "F_pair_only": f_pair,
            "frob_all_gates": frob_all,
            "frob_pair_only": frob_pair,
        },
        "qubit": {
            "F_all_gates": f_all_qb,
            "F_pair_only": f_pair_qb,
            "frob_all_gates": frob_all_qb,
            "frob_pair_only": frob_pair_qb,
            "trace_all_gates": tr_all_qb,
            "trace_pair_only": tr_pair_qb,
        },
        "cross_pair_only": {
            "F_raw": f_cross_pair,
            "F_norm": f_cross_pair_norm,
        },
    }


def test_shot_pair_only_fidelity(params: GKSLPhysicalParameters) -> dict:
    """Test 2: Shot-based pair-only noise with per-step fidelity tracking.

    Runs both qudit and qubit shot-based simulators with pair-only depolarization
    and tracks the per-step population deviation from the noiseless reference.
    """
    print("\n" + "=" * 70)
    print("Test 2: Shot-Based Pair-Only Noise Per-Step Fidelity")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    n_shots = 1000
    seed = 42
    p_depol = 0.01

    # --- Reference simulations ---
    print("\n  Running noiseless references...", flush=True)
    sim_qd_ref = QuditGKSLSimulator(params)
    r_qd_ref = sim_qd_ref.simulate(t_max=t_max, n_steps=n_steps)

    sim_qb_ref = QubitGKSLSimulator(params)
    r_qb_ref = sim_qb_ref.simulate(t_max=t_max, n_steps=n_steps)

    # --- Pair-only noise simulations ---
    print("  Running qudit pair-only noisy shot...", flush=True)
    sim_qd_pair = QuditGKSLNoisyShotSimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_qd_pair = sim_qd_pair.simulate(
        t_max=t_max, n_steps=n_steps, n_shots=n_shots, seed=seed,
    )

    print("  Running qubit pair-only noisy shot...", flush=True)
    sim_qb_pair = QubitGKSLNoisyShotSimulator(params, p_depol=p_depol, depol_pair_only=True)
    r_qb_pair = sim_qb_pair.simulate(
        t_max=t_max, n_steps=n_steps, n_shots=n_shots, seed=seed,
    )

    # --- All-gates noise for comparison ---
    print("  Running qudit all-gates noisy shot...", flush=True)
    sim_qd_all = QuditGKSLNoisyShotSimulator(params, p_depol=p_depol, depol_pair_only=False)
    r_qd_all = sim_qd_all.simulate(
        t_max=t_max, n_steps=n_steps, n_shots=n_shots, seed=seed,
    )

    print("  Running qubit all-gates noisy shot...", flush=True)
    sim_qb_all = QubitGKSLNoisyShotSimulator(params, p_depol=p_depol, depol_pair_only=False)
    r_qb_all = sim_qb_all.simulate(
        t_max=t_max, n_steps=n_steps, n_shots=n_shots, seed=seed,
    )

    # --- Per-step comparison ---
    print("\n  === Per-Step Population Deviation (vs Noiseless DM Reference) ===")
    print(f"  {'Step':>4} | {'t':>5} | {'Qdit pair |ΔN|':>15} | {'Qdit all |ΔN|':>15} | "
          f"{'Qbit pair |ΔN|':>15} | {'Qbit all |ΔN|':>15} | "
          f"{'Qbit Tr pair':>13} | {'Qbit Tr all':>13}")
    print("  " + "-" * 115)

    step_data = []
    for s in range(n_steps + 1):
        pop_qd_ref = r_qd_ref["populations"][s]
        pop_qd_pair = r_qd_pair["populations"][s]
        pop_qd_all = r_qd_all["populations"][s]
        pop_qb_ref = r_qb_ref["populations"][s]
        pop_qb_pair = r_qb_pair["populations"][s]
        pop_qb_all = r_qb_all["populations"][s]

        # Sum of absolute population deviations (all 3 species)
        dn_qd_pair = sum(
            abs(pop_qd_pair[k] - pop_qd_ref[k]) for k in ("N_S0", "N_T1", "N_S1")
        )
        dn_qd_all = sum(
            abs(pop_qd_all[k] - pop_qd_ref[k]) for k in ("N_S0", "N_T1", "N_S1")
        )
        dn_qb_pair = sum(
            abs(pop_qb_pair[k] - pop_qb_ref[k]) for k in ("N_S0", "N_T1", "N_S1")
        )
        dn_qb_all = sum(
            abs(pop_qb_all[k] - pop_qb_ref[k]) for k in ("N_S0", "N_T1", "N_S1")
        )
        tr_qb_pair = r_qb_pair["trace"][s]
        tr_qb_all = r_qb_all["trace"][s]

        sd = {
            "step": s,
            "time": r_qd_ref["times"][s],
            "delta_N_qudit_pair": float(dn_qd_pair),
            "delta_N_qudit_all": float(dn_qd_all),
            "delta_N_qubit_pair": float(dn_qb_pair),
            "delta_N_qubit_all": float(dn_qb_all),
            "trace_qubit_pair": float(tr_qb_pair),
            "trace_qubit_all": float(tr_qb_all),
            "populations_qudit_pair": {k: pop_qd_pair[k] for k in ("N_S0", "N_T1", "N_S1")},
            "populations_qubit_pair": {k: pop_qb_pair[k] for k in ("N_S0", "N_T1", "N_S1")},
            "populations_ref": {k: pop_qd_ref[k] for k in ("N_S0", "N_T1", "N_S1")},
        }
        step_data.append(sd)
        print(
            f"  {s:4d} | {sd['time']:5.1f} | {dn_qd_pair:15.6f} | {dn_qd_all:15.6f} | "
            f"{dn_qb_pair:15.6f} | {dn_qb_all:15.6f} | "
            f"{tr_qb_pair:13.8f} | {tr_qb_all:13.8f}"
        )

    # --- Final fidelity ---
    rho_qd_ref_f = r_qd_ref["rho_final"]
    rho_qd_pair_f = r_qd_pair["rho_final"]
    rho_qd_all_f = r_qd_all["rho_final"]
    rho_qb_ref_f = r_qb_ref["rho_final"]
    rho_qb_pair_f = r_qb_pair["rho_final"]
    rho_qb_all_f = r_qb_all["rho_final"]

    f_qd_pair = quantum_fidelity(rho_qd_ref_f, rho_qd_pair_f)
    f_qd_all = quantum_fidelity(rho_qd_ref_f, rho_qd_all_f)
    f_qb_pair_raw = quantum_fidelity(rho_qb_ref_f, rho_qb_pair_f)
    f_qb_all_raw = quantum_fidelity(rho_qb_ref_f, rho_qb_all_f)
    f_qb_pair_norm = quantum_fidelity(
        normalize_rho(rho_qb_ref_f), normalize_rho(rho_qb_pair_f)
    )
    f_qb_all_norm = quantum_fidelity(
        normalize_rho(rho_qb_ref_f), normalize_rho(rho_qb_all_f)
    )

    # Cross comparison (qudit pair vs qubit pair)
    f_cross_raw = quantum_fidelity(rho_qd_pair_f, rho_qb_pair_f)
    f_cross_norm = quantum_fidelity(
        normalize_rho(rho_qd_pair_f), normalize_rho(rho_qb_pair_f)
    )

    print(f"\n  === Final State Fidelity (10 steps) ===")
    print(f"  Qudit pair-only:  F = {f_qd_pair:.6f}   (all-gates: {f_qd_all:.6f})")
    print(f"  Qubit pair-only:  F_raw = {f_qb_pair_raw:.6f}, F_norm = {f_qb_pair_norm:.6f}")
    print(f"  Qubit all-gates:  F_raw = {f_qb_all_raw:.6f}, F_norm = {f_qb_all_norm:.6f}")
    print(f"  Cross (qudit pair vs qubit pair): F_raw = {f_cross_raw:.6f}, F_norm = {f_cross_norm:.6f}")
    print(f"  Qudit pair advantage over qubit pair (F): {f_qd_pair - f_qb_pair_norm:+.6f}")

    fidelity_results = {
        "qudit_pair_only": f_qd_pair,
        "qudit_all_gates": f_qd_all,
        "qubit_pair_only_raw": f_qb_pair_raw,
        "qubit_pair_only_norm": f_qb_pair_norm,
        "qubit_all_gates_raw": f_qb_all_raw,
        "qubit_all_gates_norm": f_qb_all_norm,
        "cross_pair_raw": f_cross_raw,
        "cross_pair_norm": f_cross_norm,
    }

    return {
        "test": "shot_pair_only_fidelity",
        "config": {
            "t_max": t_max,
            "n_steps": n_steps,
            "n_shots": n_shots,
            "seed": seed,
            "p_depol": p_depol,
        },
        "step_data": step_data,
        "fidelity_results": fidelity_results,
    }


def test_error_budget_comparison(params: GKSLPhysicalParameters) -> dict:
    """Test 3: Error budget comparison between all-gates and pair-only modes."""
    print("\n" + "=" * 70)
    print("Test 3: Error Budget Comparison (All-Gates vs Pair-Only)")
    print("=" * 70)

    d = params.d
    N = params.N_molecules
    n_neighbors = len(params.neighbors)
    p_depol = 0.01

    # All-gates mode: noise on every gate
    all_pair_depol = 2 * n_neighbors + 2 * n_neighbors  # Hamiltonian(2×n_neighbors) + TTA(2×n_neighbors)
    all_single_depol = 5 * N  # 20 single-site Lindblad channels
    all_total = all_pair_depol + all_single_depol  # 32 depol (excluding dephasing)

    # Pair-only mode: noise only on pair interactions
    pair_pair_depol = 2 * n_neighbors + 2 * n_neighbors  # Same: 12
    pair_single_depol = 0  # Skipped!
    pair_total = pair_pair_depol  # 12 depol only

    # Expected error events per step
    p_error_pair = p_depol * (1 - 1 / (d ** 4))  # prob of non-identity in pair depol
    p_error_single = p_depol * (1 - 1 / (d ** 2))  # prob of non-identity in single depol

    exp_all = all_pair_depol * p_error_pair + all_single_depol * p_error_single
    exp_pair = pair_pair_depol * p_error_pair

    noise_reduction = 1.0 - pair_total / all_total if all_total > 0 else 0.0

    print(f"\n  === Noise Gates Per Trotter Step (depolarization only) ===")
    print(f"  {'':30s} | {'All-gates':>12} | {'Pair-only':>12} |")
    print(f"  {'-' * 60}")
    print(f"  {'Pair depolarization':30s} | {all_pair_depol:12d} | {pair_pair_depol:12d} |")
    print(f"  {'Single depolarization':30s} | {all_single_depol:12d} | {pair_single_depol:12d} |")
    print(f"  {'Total depolarization':30s} | {all_total:12d} | {pair_total:12d} |")
    print(f"  {'Noise reduction':30s} | {'-':>12s} | {noise_reduction:11.1%} |")

    print(f"\n  === Expected Errors Per Step (p_depol={p_depol}) ===")
    print(f"  All-gates:  {exp_all:.4f} errors/step")
    print(f"  Pair-only:  {exp_pair:.4f} errors/step")
    print(f"  Reduction:  {1.0 - exp_pair / exp_all:.1%}")

    # For 100 steps
    print(f"\n  === Cumulative (100 steps) ===")
    print(f"  All-gates:  {exp_all * 100:.1f} errors")
    print(f"  Pair-only:  {exp_pair * 100:.1f} errors")

    # Qubit-specific leakage analysis for pair-only mode
    n_paulis_pair = 256  # 4^4 Pauli operators for 4-qubit pair
    # In pair-only mode, only pair depolarization can cause leakage
    # Each pair depol applies a random 4-qubit Pauli (2 molecules × 2 qubits each)
    # Leakage occurs if any Pauli component at a molecule maps physical → forbidden
    print(f"\n  === Qubit Pair-Only Leakage Analysis ===")
    print(f"  Only pair depolarization can cause leakage in pair-only mode")
    print(f"  Pair depol gates per step: {pair_pair_depol}")
    print(f"  Per molecule: 4×4 Pauli, 12/15 non-identity cause |11> leakage")

    return {
        "test": "error_budget_comparison",
        "all_gates": {
            "pair_depol": all_pair_depol,
            "single_depol": all_single_depol,
            "total": all_total,
            "expected_errors_per_step": exp_all,
        },
        "pair_only": {
            "pair_depol": pair_pair_depol,
            "single_depol": pair_single_depol,
            "total": pair_total,
            "expected_errors_per_step": exp_pair,
        },
        "noise_reduction_fraction": noise_reduction,
        "error_reduction_fraction": 1.0 - exp_pair / exp_all if exp_all > 0 else 0.0,
    }


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
    print("TTA-UC GKSL Iteration 14: Pair-Only Depolarization Comparison")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 14,
        "purpose": "Pair-only depolarization: qubit vs qudit accuracy comparison",
        "tests": [],
    }

    # Test 1: DM pair-only comparison
    t1 = test_dm_pair_only_comparison(params)
    results["tests"].append(t1)

    # Test 2: Shot-based pair-only fidelity
    t2 = test_shot_pair_only_fidelity(params)
    results["tests"].append(t2)

    # Test 3: Error budget comparison
    t3 = test_error_budget_comparison(params)
    results["tests"].append(t3)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration14_pair_only_noise_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 14: ペアのみ脱分極ノイズ比較結果",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "",
        "## Test 1: DM ノイジーシミュレータ ペアのみ vs 全ゲート (1 Trotterステップ)",
        "",
        f"- 脱分極確率: p_depol = {t1['p_depol']}",
        "",
        "### Qudit (d=3, dim=81)",
        "",
        f"- 全ゲートノイズ: F = {t1['qudit']['F_all_gates']:.6f}, ||Δ||_F = {t1['qudit']['frob_all_gates']:.6f}",
        f"- ペアのみノイズ: F = {t1['qudit']['F_pair_only']:.6f}, ||Δ||_F = {t1['qudit']['frob_pair_only']:.6f}",
        "",
        "### Qubit (d=4, dim=256)",
        "",
        f"- 全ゲートノイズ: F = {t1['qubit']['F_all_gates']:.6f}, Tr = {t1['qubit']['trace_all_gates']:.8f}",
        f"- ペアのみノイズ: F = {t1['qubit']['F_pair_only']:.6f}, Tr = {t1['qubit']['trace_pair_only']:.8f}",
        "",
        "### Qudit vs Qubit (ペアのみ)",
        "",
        f"- F_raw = {t1['cross_pair_only']['F_raw']:.6f}",
        f"- F_norm = {t1['cross_pair_only']['F_norm']:.6f}",
        "",
    ]

    # Test 2
    t2_fid = t2["fidelity_results"]
    t2_cfg = t2["config"]
    md_lines.extend([
        "## Test 2: ショットベース ペアのみノイズ ステップ毎忠実度 (10ステップ)",
        "",
        f"- t_max: {t2_cfg['t_max']}, n_steps: {t2_cfg['n_steps']}, n_shots: {t2_cfg['n_shots']}",
        f"- p_depol: {t2_cfg['p_depol']}, p_dephasing: 0.0",
        "",
        "### ステップ毎の総占有数偏差 |ΔN|",
        "",
        "| Step | Time | Qudit pair |ΔN| | Qudit all |ΔN| | Qubit pair |ΔN| | Qubit all |ΔN| | Qubit Tr pair | Qubit Tr all |",
        "|------|------|------------|------------|------------|------------|------------|------------|",
    ])
    for sd in t2["step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['time']:.1f} "
            f"| {sd['delta_N_qudit_pair']:.6f} | {sd['delta_N_qudit_all']:.6f} "
            f"| {sd['delta_N_qubit_pair']:.6f} | {sd['delta_N_qubit_all']:.6f} "
            f"| {sd['trace_qubit_pair']:.8f} | {sd['trace_qubit_all']:.8f} |"
        )
    md_lines.extend([
        "",
        "### 最終状態忠実度 (10ステップ後)",
        "",
        f"- **Qudit ペアのみ: F = {t2_fid['qudit_pair_only']:.6f}** (全ゲート: {t2_fid['qudit_all_gates']:.6f})",
        f"- **Qubit ペアのみ: F_norm = {t2_fid['qubit_pair_only_norm']:.6f}** (全ゲート: {t2_fid['qubit_all_gates_norm']:.6f})",
        f"- Cross (qudit pair vs qubit pair): F_raw = {t2_fid['cross_pair_raw']:.6f}, F_norm = {t2_fid['cross_pair_norm']:.6f}",
        "",
    ])

    # Test 3
    md_lines.extend([
        "## Test 3: エラーバジェット比較",
        "",
        "### 1ステップあたりの脱分極ノイズ適用回数",
        "",
        "| 種別 | 全ゲート | ペアのみ |",
        "|------|---------|---------|",
        f"| ペア脱分極 | {t3['all_gates']['pair_depol']} | {t3['pair_only']['pair_depol']} |",
        f"| 単一サイト脱分極 | {t3['all_gates']['single_depol']} | {t3['pair_only']['single_depol']} |",
        f"| **合計** | **{t3['all_gates']['total']}** | **{t3['pair_only']['total']}** |",
        "",
        f"- ノイズゲート削減率: **{t3['noise_reduction_fraction']:.1%}**",
        f"- 期待エラー削減率: **{t3['error_reduction_fraction']:.1%}**",
        "",
        "### 結論",
        "",
        "ペアのみモードでは、20個の単一サイトLindblad後のノイズを省略し、",
        "ペア相互作用（Hamiltonian転送+TTA Lindblad）のみにノイズを適用する。",
        "これは物理的により現実的なモデルであり、qubit vs qudit比較に適している。",
        "",
    ])

    md_path = output_dir / f"iteration14_pair_only_noise_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Test 1 (DM 1-step):")
    print(f"    Qudit F(pair-only)={t1['qudit']['F_pair_only']:.6f} vs F(all)={t1['qudit']['F_all_gates']:.6f}")
    print(f"    Qubit F(pair-only)={t1['qubit']['F_pair_only']:.6f} vs F(all)={t1['qubit']['F_all_gates']:.6f}")
    print(f"  Test 2 (Shot 10-step):")
    print(f"    Qudit pair-only: F={t2_fid['qudit_pair_only']:.6f}")
    print(f"    Qubit pair-only: F_norm={t2_fid['qubit_pair_only_norm']:.6f}")
    print(f"    Cross pair-only: F_norm={t2_fid['cross_pair_norm']:.6f}")
    print(f"  Test 3 (Error budget):")
    print(f"    Noise gate reduction: {t3['noise_reduction_fraction']:.1%}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
