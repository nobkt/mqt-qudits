#!/usr/bin/env python3
"""TTA-UC GKSL stable regression test script.

This is the stable (non-iteration-specific) regression test script that validates
the mathematical and physical correctness of the GKSL simulator implementation.

It runs all regression checks (convergence, density matrix quality, particle
conservation, simulator consistency, etc.) WITHOUT self-referential iteration
counter checks that caused the infinite verification loop (checks 7, 8, 31).

See developing/検証結果分析_iteration45_根本原因分析.md for root cause analysis.

Results are written to developing/verification_results/ as JSON and Markdown.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time as time_module
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_simulator import QubitGKSLSimulator
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
    """Compute trace distance T(rho, sigma) = 0.5 * ||rho - sigma||_1."""
    delta = rho - sigma
    delta = (delta + delta.conj().T) / 2
    eigenvalues = np.linalg.eigvalsh(delta)
    return float(0.5 * np.sum(np.abs(eigenvalues)))


def check_density_matrix_quality(rho: np.ndarray) -> dict:
    """Check trace, Hermiticity, and positivity of a density matrix."""
    tr = float(np.real(np.trace(rho)))
    herm_err = float(np.max(np.abs(rho - rho.conj().T)))
    eigenvalues = np.linalg.eigvalsh(rho)
    min_eig = float(np.min(eigenvalues))
    return {
        "trace": tr,
        "trace_error": abs(tr - 1.0),
        "hermiticity_error": herm_err,
        "min_eigenvalue": min_eig,
    }


def parse_cell27_table_from_notebook(notebook_path: str) -> dict:
    """Read the notebook ipynb file and parse Cell 27 markdown table values."""
    with open(notebook_path) as f:
        nb = json.load(f)

    cell27_source = "".join(nb["cells"][27]["source"])

    table_values = {}
    for line in cell27_source.split("\n"):
        m = re.match(
            r"\|\s*(\d+)\s*\|\s*[\d.]+\s*\|\s*~?([\d.eE+-]+)\s*\|",
            line,
        )
        if m:
            n_steps = int(m.group(1))
            t_value = float(m.group(2))
            table_values[n_steps] = t_value

    return {
        "table_values": table_values,
        "full_text": cell27_source,
    }


def check_documentation_fixes(notebook_path: str) -> dict:
    """Check that documentation content is consistent with computed values.

    Note: Iteration count checks (cell 25/39 "N回") have been removed because
    they created a self-referential infinite loop. See developing/検証結果分析_iteration45_根本原因分析.md.
    """
    with open(notebook_path) as f:
        nb = json.load(f)

    results = {}

    # Cell 25/39 iteration count checks removed (self-referential loop)
    results["cell25_ok"] = True
    results["cell39_ok"] = True

    cell27_text = "".join(nb["cells"][27]["source"])
    has_36e04 = "~3.6e-04" in cell27_text
    has_old_4e04 = "~4e-04" in cell27_text
    results["cell27_has_36e04"] = has_36e04
    results["cell27_has_old_4e04"] = has_old_4e04
    results["cell27_n200_ok"] = has_36e04 and not has_old_4e04

    results["all_doc_fixes_ok"] = results["cell27_n200_ok"]

    return results


def main() -> None:
    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    results_dir = Path(__file__).resolve().parent.parent / "developing" / "verification_results"
    results_dir.mkdir(parents=True, exist_ok=True)

    notebook_path = str(Path(__file__).resolve().parent / "quantum_dynamics_gksl_comparison.ipynb")

    params = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
    )
    t_max = 100.0

    results: dict = {
        "iteration": "stable",
        "timestamp": timestamp,
        "description": "Stable regression tests (self-referential iteration counter checks removed)",
        "params": {
            "t_max": t_max,
            "hilbert_dim": params.get_hilbert_space_dim(),
        },
    }

    print("=" * 70)
    print("TTA-UC GKSL Stable Regression Test")
    print("=" * 70)

    # ===================================================================
    # Part 1: Exact reference solution
    # ===================================================================
    print("\n--- Part 1: Exact reference (ClassicalGKSLSimulator) ---")
    t0 = time_module.time()
    sim_exact = ClassicalGKSLSimulator(params)
    result_exact = sim_exact.simulate(t_max=t_max, n_steps=1000, initial_state="edge_triplet")
    t_exact = time_module.time() - t0
    rho_exact = result_exact["rho_final"]
    quality_exact = check_density_matrix_quality(rho_exact)
    print(f"  Elapsed: {t_exact:.2f}s")
    print(f"  Trace: {quality_exact['trace']:.10f}")
    results["exact_reference"] = {
        "elapsed_time": t_exact,
        "quality": quality_exact,
    }

    # ===================================================================
    # Part 2: Convergence analysis (4-molecule system)
    # ===================================================================
    print("\n--- Part 2: Convergence analysis (4-molecule) ---")
    n_steps_list = [10, 20, 50, 100, 200]
    convergence_data = []

    result_qudit_100 = None

    for ns in n_steps_list:
        dt = t_max / ns
        t0 = time_module.time()
        sim_st = QuditGKSLSimulator(params)
        result_st = sim_st.simulate(t_max=t_max, n_steps=ns, initial_state="edge_triplet")
        elapsed = time_module.time() - t0
        rho_st = result_st["rho_final"]

        td = trace_distance(rho_exact, rho_st)
        fid = quantum_fidelity(rho_exact, rho_st)
        quality = check_density_matrix_quality(rho_st)

        entry = {
            "n_steps": ns,
            "dt": dt,
            "trace_distance": td,
            "fidelity": fid,
            "elapsed_time": elapsed,
            "quality": quality,
        }
        convergence_data.append(entry)
        print(f"  n_steps={ns:>4}, dt={dt:>8.4f}, T={td:.6e}, F={fid:.8f}")

        if ns == 100:
            result_qudit_100 = result_st

    for i in range(1, len(convergence_data)):
        prev = convergence_data[i - 1]
        curr = convergence_data[i]
        ratio_T = prev["trace_distance"] / curr["trace_distance"]
        ratio_dt = prev["dt"] / curr["dt"]
        if ratio_dt > 1 and ratio_T > 0:
            rate = float(np.log(ratio_T) / np.log(ratio_dt))
        else:
            rate = float("nan")
        convergence_data[i]["rate_T"] = rate
        print(f"  Rate(T) [{prev['n_steps']}->{curr['n_steps']}]: {rate:.4f}")

    results["convergence_data"] = convergence_data

    # ===================================================================
    # Part 3: Original validation checks (1–6)
    # ===================================================================
    print("\n--- Part 3: Original validation checks (1-6) ---")
    checks = {}

    tds = [d["trace_distance"] for d in convergence_data]
    monotonic = all(tds[i] > tds[i + 1] for i in range(len(tds) - 1))
    checks["check1_trace_distance_monotonic"] = monotonic
    print(f"  Check 1 - Trace distance monotonically decreasing: {monotonic}")

    rates = [d.get("rate_T", float("nan")) for d in convergence_data if "rate_T" in d]
    avg_rate = float(np.mean(rates)) if rates else float("nan")
    rate_near_1 = abs(avg_rate - 1.0) < 0.3
    checks["check2_convergence_rate"] = rate_near_1
    checks["check2_average_rate"] = avg_rate
    print(f"  Check 2 - Average Rate(T): {avg_rate:.4f} (near 1.0: {rate_near_1})")

    fids = [d["fidelity"] for d in convergence_data]
    fid_monotonic = all(fids[i] < fids[i + 1] for i in range(len(fids) - 1))
    checks["check3_fidelity_monotonic"] = fid_monotonic
    print(f"  Check 3 - Fidelity monotonically increasing: {fid_monotonic}")

    all_quality_ok = all(
        d["quality"]["trace_error"] < 1e-10
        and d["quality"]["hermiticity_error"] < 1e-10
        and d["quality"]["min_eigenvalue"] > -1e-10
        for d in convergence_data
    )
    checks["check4_density_matrix_quality"] = all_quality_ok
    print(f"  Check 4 - Density matrix quality OK: {all_quality_ok}")

    t100 = next(d["trace_distance"] for d in convergence_data if d["n_steps"] == 100)
    practical_accuracy = t100 < 1e-3
    checks["check5_n100_accuracy"] = practical_accuracy
    checks["check5_n100_trace_distance"] = t100
    print(f"  Check 5 - n_steps=100 T < 1e-3: {practical_accuracy} (T={t100:.6e})")

    # Check 6: Cell 27 content validation
    print("\n--- Cell 27 notebook content validation ---")
    parsed = parse_cell27_table_from_notebook(notebook_path)
    table_values = parsed["table_values"]
    full_text = parsed["full_text"]

    missing = [ns for ns in [10, 20, 50, 100, 200] if ns not in table_values]
    all_present = len(missing) == 0

    comparison = {}
    all_match = True
    for d in convergence_data:
        ns = d["n_steps"]
        if ns in table_values:
            nv = table_values[ns]
            av = d["trace_distance"]
            ratio = av / nv
            ok = 0.5 < ratio < 2.0
            comparison[f"n_steps_{ns}"] = {
                "notebook_value": nv,
                "actual_value": av,
                "ratio": ratio,
                "within_tolerance": ok,
            }
            status = "✓" if ok else "✗"
            print(f"  {status} n_steps={ns}: notebook={nv:.2e}, actual={av:.4e}, ratio={ratio:.3f}")
            if not ok:
                all_match = False

    no_incorrect_2e5 = "2×10⁻⁵" not in full_text and "2e-5" not in full_text and "2e-05" not in full_text
    has_tmax = "100.0" in full_text or "100" in full_text

    cell27_ok = all_present and all_match and no_incorrect_2e5 and has_tmax
    checks["check6_cell27_validation"] = {
        "all_present": all_present,
        "all_match": all_match,
        "no_incorrect_2e5": no_incorrect_2e5,
        "has_tmax": has_tmax,
        "comparison": comparison,
        "overall": cell27_ok,
    }
    print(f"  Cell 27 overall: {cell27_ok}")

    # ===================================================================
    # Part 4: Documentation fixes validation (checks 7-9)
    # ===================================================================
    print("\n--- Part 4: Documentation fixes validation (7-9) ---")
    doc_fixes = check_documentation_fixes(notebook_path)
    checks["check7_cell25_iteration_count"] = doc_fixes["cell25_ok"]
    checks["check8_cell39_iteration_count"] = doc_fixes["cell39_ok"]
    checks["check9_cell27_n200_value"] = doc_fixes["cell27_n200_ok"]

    print(f"  Check 7 - Cell 25 (iteration count check removed): {doc_fixes['cell25_ok']}")
    print(f"  Check 8 - Cell 39 (iteration count check removed): {doc_fixes['cell39_ok']}")
    print(f"  Check 9 - Cell 27 n_steps=200 is '~3.6e-04': {doc_fixes['cell27_n200_ok']}")

    # ===================================================================
    # Part 5: Particle number conservation (check 10)
    # ===================================================================
    print("\n--- Part 5: Particle number conservation (check 10) ---")
    sim_exact_pn = ClassicalGKSLSimulator(params)
    result_exact_pn = sim_exact_pn.simulate(t_max=t_max, n_steps=100, initial_state="edge_triplet")

    max_pnum_err_cl = 0.0
    for pop in result_exact_pn["populations"]:
        total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
        max_pnum_err_cl = max(max_pnum_err_cl, abs(total - 4.0))

    max_pnum_err_qd = 0.0
    assert result_qudit_100 is not None
    for pop in result_qudit_100["populations"]:
        total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
        max_pnum_err_qd = max(max_pnum_err_qd, abs(total - 4.0))

    pnum_ok = bool(max_pnum_err_cl < 1e-12 and max_pnum_err_qd < 1e-12)
    checks["check10_particle_conservation"] = {
        "classical_max_error": float(max_pnum_err_cl),
        "qudit_max_error": float(max_pnum_err_qd),
        "ok": pnum_ok,
    }
    print(f"  Classical max |N_total - 4|: {max_pnum_err_cl:.2e}")
    print(f"  Qudit max |N_total - 4|: {max_pnum_err_qd:.2e}")
    print(f"  Check 10 - Particle conservation: {pnum_ok}")

    # ===================================================================
    # Part 6: Intermediate-time population comparison (check 11)
    # ===================================================================
    print("\n--- Part 6: Intermediate-time population comparison (check 11) ---")

    interp_data = []
    max_pop_diff = 0.0

    for t_idx in range(0, 101, 10):
        pop_cl = result_exact_pn["populations"][t_idx]
        pop_qd = result_qudit_100["populations"][t_idx]
        t = result_exact_pn["times"][t_idx]

        delta_S0 = float(abs(pop_cl["N_S0"] - pop_qd["N_S0"]))
        delta_T1 = float(abs(pop_cl["N_T1"] - pop_qd["N_T1"]))
        delta_S1 = float(abs(pop_cl["N_S1"] - pop_qd["N_S1"]))
        max_delta = max(delta_S0, delta_T1, delta_S1)
        max_pop_diff = max(max_pop_diff, max_delta)

        interp_data.append({
            "time": t,
            "delta_N_S0": delta_S0,
            "delta_N_T1": delta_T1,
            "delta_N_S1": delta_S1,
            "max_delta": max_delta,
        })
        print(f"  t={t:6.1f}: max|ΔN|={max_delta:.4e}")

    interp_ok = bool(max_pop_diff < 5e-3)
    checks["check11_intermediate_populations"] = {
        "max_population_difference": float(max_pop_diff),
        "data": interp_data,
        "ok": interp_ok,
    }
    print(f"  Max pop difference over all times: {max_pop_diff:.4e}")
    print(f"  Check 11 - Intermediate populations: {interp_ok}")

    # ===================================================================
    # Part 7: 2-molecule system convergence (check 12)
    # ===================================================================
    print("\n--- Part 7: 2-molecule system convergence (check 12) ---")

    params2 = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
        N_molecules=2,
    )

    sim_exact2 = ClassicalGKSLSimulator(params2)
    result_exact2 = sim_exact2.simulate(t_max=100.0, n_steps=1000, initial_state="edge_triplet")
    rho_exact2 = result_exact2["rho_final"]

    n_steps_2mol = [10, 20, 50, 100, 200, 500]
    conv2_data = []
    for ns in n_steps_2mol:
        sim_st2 = QuditGKSLSimulator(params2)
        result_st2 = sim_st2.simulate(t_max=100.0, n_steps=ns, initial_state="edge_triplet")
        td = trace_distance(rho_exact2, result_st2["rho_final"])
        conv2_data.append({"n_steps": ns, "dt": 100.0 / ns, "trace_distance": td})
        print(f"  2-mol: n_steps={ns:>4}, dt={100.0/ns:>8.4f}, T={td:.6e}")

    rates2 = []
    for i in range(1, len(conv2_data)):
        prev = conv2_data[i - 1]
        curr = conv2_data[i]
        ratio_T = prev["trace_distance"] / curr["trace_distance"]
        ratio_dt = prev["dt"] / curr["dt"]
        if ratio_dt > 1 and ratio_T > 0:
            rate2 = float(np.log(ratio_T) / np.log(ratio_dt))
        else:
            rate2 = float("nan")
        conv2_data[i]["rate_T"] = rate2
        rates2.append(rate2)
        print(f"  Rate(T) [{prev['n_steps']}->{curr['n_steps']}]: {rate2:.4f}")

    last_rate_2mol = rates2[-1] if rates2 else float("nan")
    conv2_ok = bool(abs(last_rate_2mol - 1.0) < 0.01)
    checks["check12_2mol_convergence"] = {
        "convergence_data": conv2_data,
        "last_rate": last_rate_2mol,
        "ok": conv2_ok,
    }
    print(f"  Last Rate(T) [2-mol]: {last_rate_2mol:.4f}")
    print(f"  Check 12 - 2-mol convergence: {conv2_ok}")

    # ===================================================================
    # Part 8: Unitary-only accuracy (check 13)
    # ===================================================================
    print("\n--- Part 8: Unitary-only (zero dissipation) accuracy (check 13) ---")

    params_u = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.0, Gamma_fl=0.0, Gamma_ph=0.0,
        k_IC=0.0, k_ISC_ST=0.0, k_ISC_TS=0.0,
        N_molecules=2,
    )

    sim_exact_u = ClassicalGKSLSimulator(params_u)
    result_exact_u = sim_exact_u.simulate(t_max=100.0, n_steps=1000, initial_state="edge_triplet")
    rho_exact_u = result_exact_u["rho_final"]

    sim_st_u = QuditGKSLSimulator(params_u)
    result_st_u = sim_st_u.simulate(t_max=100.0, n_steps=100, initial_state="edge_triplet")
    td_u = trace_distance(rho_exact_u, result_st_u["rho_final"])

    unitary_ok = bool(td_u < 1e-12)
    checks["check13_unitary_only"] = {
        "trace_distance": td_u,
        "ok": unitary_ok,
    }
    print(f"  Unitary-only T(exact, trotter): {td_u:.2e}")
    print(f"  Check 13 - Unitary-only accuracy: {unitary_ok}")

    # ===================================================================
    # Part 9: Eigenvalue spectrum comparison (check 14)
    # ===================================================================
    print("\n--- Part 9: Eigenvalue spectrum comparison (check 14) ---")

    evals_exact = np.sort(np.linalg.eigvalsh(rho_exact))[::-1]
    rho_st_100 = result_qudit_100["rho_final"]
    evals_st = np.sort(np.linalg.eigvalsh(rho_st_100))[::-1]

    max_eval_diff = float(np.max(np.abs(evals_exact - evals_st)))
    top5_diffs = [float(abs(evals_exact[i] - evals_st[i])) for i in range(5)]

    print(f"  Top 5 eigenvalue differences:")
    for i in range(5):
        print(f"    λ_{i}: exact={evals_exact[i]:.10f}, stine={evals_st[i]:.10f}, diff={top5_diffs[i]:.6e}")

    eval_ok = bool(max_eval_diff < 2e-3)
    checks["check14_eigenvalue_spectrum"] = {
        "max_eigenvalue_difference": float(max_eval_diff),
        "top5_differences": [float(d) for d in top5_diffs],
        "ok": eval_ok,
    }
    print(f"  Max eigenvalue difference: {max_eval_diff:.6e}")
    print(f"  Check 14 - Eigenvalue spectrum: {eval_ok}")

    # ===================================================================
    # Part 10: Shot-based qudit vs DM qudit consistency (check 15)
    # ===================================================================
    print("\n--- Part 10: Shot-based qudit vs DM qudit consistency (check 15) ---")

    from qudit_gksl_shot_simulator import QuditGKSLShotSimulator

    sim_dm2 = QuditGKSLSimulator(params2)
    result_dm2 = sim_dm2.simulate(t_max=10.0, n_steps=10, initial_state="edge_triplet")
    rho_dm2 = result_dm2["rho_final"]

    sim_shot2 = QuditGKSLShotSimulator(params2)
    result_shot2 = sim_shot2.simulate(
        t_max=10.0, n_steps=10, initial_state="edge_triplet",
        n_shots=5000, seed=42,
    )
    rho_shot2 = result_shot2["rho_final"]

    td_shot_dm = trace_distance(rho_dm2, rho_shot2)
    fid_shot_dm = quantum_fidelity(rho_dm2, rho_shot2)

    shot_dm_ok = bool(td_shot_dm < 0.05)
    checks["check15_qudit_shot_dm_consistency"] = {
        "trace_distance": td_shot_dm,
        "fidelity": fid_shot_dm,
        "n_shots": 5000,
        "ok": shot_dm_ok,
    }
    print(f"  T(qudit_shot, qudit_dm): {td_shot_dm:.6e}")
    print(f"  F(qudit_shot, qudit_dm): {fid_shot_dm:.8f}")
    print(f"  Check 15 - Qudit shot-DM consistency: {shot_dm_ok}")

    # ===================================================================
    # Part 11: Shot-based qubit vs DM qubit consistency (check 16)
    # ===================================================================
    print("\n--- Part 11: Shot-based qubit vs DM qubit consistency (check 16) ---")

    from qubit_gksl_shot_simulator import QubitGKSLShotSimulator

    sim_qb_dm2 = QubitGKSLSimulator(params2)
    result_qb_dm2 = sim_qb_dm2.simulate(t_max=10.0, n_steps=10, initial_state="edge_triplet")
    rho_qb_dm2 = result_qb_dm2["rho_final"]

    sim_qb_shot2 = QubitGKSLShotSimulator(params2)
    result_qb_shot2 = sim_qb_shot2.simulate(
        t_max=10.0, n_steps=10, initial_state="edge_triplet",
        n_shots=5000, seed=42,
    )
    rho_qb_shot2 = result_qb_shot2["rho_final"]

    td_qb_shot_dm = trace_distance(rho_qb_dm2, rho_qb_shot2)
    fid_qb_shot_dm = quantum_fidelity(rho_qb_dm2, rho_qb_shot2)

    qb_shot_dm_ok = bool(td_qb_shot_dm < 0.05)
    checks["check16_qubit_shot_dm_consistency"] = {
        "trace_distance": td_qb_shot_dm,
        "fidelity": fid_qb_shot_dm,
        "n_shots": 5000,
        "ok": qb_shot_dm_ok,
    }
    print(f"  T(qubit_shot, qubit_dm): {td_qb_shot_dm:.6e}")
    print(f"  F(qubit_shot, qubit_dm): {fid_qb_shot_dm:.8f}")
    print(f"  Check 16 - Qubit shot-DM consistency: {qb_shot_dm_ok}")

    # ===================================================================
    # Part 12: Noisy DM qudit zero-noise matches ideal (check 17)
    # ===================================================================
    print("\n--- Part 12: Noisy DM qudit zero-noise matches ideal (check 17) ---")

    from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

    sim_noisy_zero = QuditGKSLNoisySimulator(params2, p_depol=0.0, p_dephasing=0.0)
    result_noisy_zero = sim_noisy_zero.simulate(
        t_max=10.0, n_steps=10, initial_state="edge_triplet",
    )
    rho_noisy_zero = result_noisy_zero["rho_final"]

    td_noisy_zero = trace_distance(rho_dm2, rho_noisy_zero)
    noisy_zero_ok = bool(td_noisy_zero < 1e-12)
    checks["check17_qudit_noisy_zero_matches_ideal"] = {
        "trace_distance": td_noisy_zero,
        "ok": noisy_zero_ok,
    }
    print(f"  T(qudit_noisy_zero, qudit_dm): {td_noisy_zero:.2e}")
    print(f"  Check 17 - Qudit noisy zero-noise matches ideal: {noisy_zero_ok}")

    # ===================================================================
    # Part 13: Noisy DM qubit zero-noise matches ideal (check 18)
    # ===================================================================
    print("\n--- Part 13: Noisy DM qubit zero-noise matches ideal (check 18) ---")

    from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

    sim_qb_noisy_zero = QubitGKSLNoisySimulator(params2, p_depol=0.0, p_dephasing=0.0)
    result_qb_noisy_zero = sim_qb_noisy_zero.simulate(
        t_max=10.0, n_steps=10, initial_state="edge_triplet",
    )
    rho_qb_noisy_zero = result_qb_noisy_zero["rho_final"]

    td_qb_noisy_zero = trace_distance(rho_qb_dm2, rho_qb_noisy_zero)
    qb_noisy_zero_ok = bool(td_qb_noisy_zero < 1e-12)
    checks["check18_qubit_noisy_zero_matches_ideal"] = {
        "trace_distance": td_qb_noisy_zero,
        "ok": qb_noisy_zero_ok,
    }
    print(f"  T(qubit_noisy_zero, qubit_dm): {td_qb_noisy_zero:.2e}")
    print(f"  Check 18 - Qubit noisy zero-noise matches ideal: {qb_noisy_zero_ok}")

    # ===================================================================
    # Part 14: NEW — Qudit DM ancilla count (4-molecule, check 19)
    # ===================================================================
    print("\n--- Part 14: Qudit DM ancilla count (4-molecule) (check 19) ---")

    sim_qd4 = QuditGKSLSimulator(params)
    n_lindblad_4 = len(sim_qd4.lindblad_ops)
    ancilla_match_4 = sim_qd4.n_ancilla_qubits == n_lindblad_4
    # Expected: 2*3(TTA) + 5*4(single-site) = 6+20 = 26
    expected_4 = 2 * len(params.neighbors) + 5 * params.N_molecules
    ancilla_expected_4 = sim_qd4.n_ancilla_qubits == expected_4

    checks["check19_qudit_dm_ancilla_4mol"] = {
        "n_ancilla_qubits": sim_qd4.n_ancilla_qubits,
        "n_lindblad_ops": n_lindblad_4,
        "expected": expected_4,
        "matches_lindblad": ancilla_match_4,
        "matches_expected": ancilla_expected_4,
        "ok": ancilla_match_4 and ancilla_expected_4,
    }
    print(f"  n_ancilla_qubits={sim_qd4.n_ancilla_qubits}, n_lindblad={n_lindblad_4}, expected={expected_4}")
    print(f"  Check 19 - Qudit DM ancilla (4-mol): {ancilla_match_4 and ancilla_expected_4}")

    # ===================================================================
    # Part 15: NEW — Qudit DM ancilla count (2-molecule, check 20)
    # ===================================================================
    print("\n--- Part 15: Qudit DM ancilla count (2-molecule) (check 20) ---")

    sim_qd2 = QuditGKSLSimulator(params2)
    n_lindblad_2 = len(sim_qd2.lindblad_ops)
    ancilla_match_2 = sim_qd2.n_ancilla_qubits == n_lindblad_2
    # Expected: 2*1(TTA) + 5*2(single-site) = 2+10 = 12
    expected_2 = 2 * len(params2.neighbors) + 5 * params2.N_molecules
    ancilla_expected_2 = sim_qd2.n_ancilla_qubits == expected_2

    checks["check20_qudit_dm_ancilla_2mol"] = {
        "n_ancilla_qubits": sim_qd2.n_ancilla_qubits,
        "n_lindblad_ops": n_lindblad_2,
        "expected": expected_2,
        "matches_lindblad": ancilla_match_2,
        "matches_expected": ancilla_expected_2,
        "ok": ancilla_match_2 and ancilla_expected_2,
    }
    print(f"  n_ancilla_qubits={sim_qd2.n_ancilla_qubits}, n_lindblad={n_lindblad_2}, expected={expected_2}")
    print(f"  Check 20 - Qudit DM ancilla (2-mol): {ancilla_match_2 and ancilla_expected_2}")

    # ===================================================================
    # Part 16: NEW — Qudit DM gate count (4-molecule, check 21)
    # ===================================================================
    print("\n--- Part 16: Qudit DM gate count (4-molecule) (check 21) ---")

    result_qd4 = sim_qd4.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
    n_lindblad4 = len(sim_qd4.lindblad_ops)
    expected_gates_4 = params.N_molecules + len(params.neighbors) + n_lindblad4 * 2
    actual_gates_4 = result_qd4["estimated_gates_per_step"]
    gates_ok_4 = actual_gates_4 == expected_gates_4

    checks["check21_qudit_dm_gate_count_4mol"] = {
        "actual_gates_per_step": actual_gates_4,
        "expected_gates_per_step": expected_gates_4,
        "breakdown": {
            "VirtRz": params.N_molecules,
            "CustomTwo_transfer": len(params.neighbors),
            "Stinespring_palindromic": n_lindblad4 * 2,
        },
        "ok": gates_ok_4,
    }
    print(f"  gates_per_step={actual_gates_4}, expected={expected_gates_4}")
    print(f"  Check 21 - Qudit DM gates (4-mol): {gates_ok_4}")

    # ===================================================================
    # Part 17: NEW — Qudit DM gate count (2-molecule, check 22)
    # ===================================================================
    print("\n--- Part 17: Qudit DM gate count (2-molecule) (check 22) ---")

    result_qd2 = sim_qd2.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
    n_lindblad2 = len(sim_qd2.lindblad_ops)
    expected_gates_2 = params2.N_molecules + len(params2.neighbors) + n_lindblad2 * 2
    actual_gates_2 = result_qd2["estimated_gates_per_step"]
    gates_ok_2 = actual_gates_2 == expected_gates_2

    checks["check22_qudit_dm_gate_count_2mol"] = {
        "actual_gates_per_step": actual_gates_2,
        "expected_gates_per_step": expected_gates_2,
        "breakdown": {
            "VirtRz": params2.N_molecules,
            "CustomTwo_transfer": len(params2.neighbors),
            "Stinespring_palindromic": n_lindblad2 * 2,
        },
        "ok": gates_ok_2,
    }
    print(f"  gates_per_step={actual_gates_2}, expected={expected_gates_2}")
    print(f"  Check 22 - Qudit DM gates (2-mol): {gates_ok_2}")

    # ===================================================================
    # Part 18: NEW — Qudit Shot gate count (4-molecule, check 23)
    # ===================================================================
    print("\n--- Part 18: Qudit Shot gate count (4-molecule) (check 23) ---")

    sim_shot4 = QuditGKSLShotSimulator(params)
    result_shot4 = sim_shot4.simulate(
        t_max=5.0, n_steps=5, initial_state="edge_triplet",
        n_shots=10, seed=0,
    )
    actual_shot_gates_4 = result_shot4["estimated_gates_per_step"]
    shot_gates_ok_4 = actual_shot_gates_4 == expected_gates_4

    checks["check23_qudit_shot_gate_count_4mol"] = {
        "actual_gates_per_step": actual_shot_gates_4,
        "expected_gates_per_step": expected_gates_4,
        "ok": shot_gates_ok_4,
    }
    print(f"  shot gates_per_step={actual_shot_gates_4}, expected={expected_gates_4}")
    print(f"  Check 23 - Qudit Shot gates (4-mol): {shot_gates_ok_4}")

    # ===================================================================
    # Part 19: NEW — Qudit Shot gate count (2-molecule, check 24)
    # ===================================================================
    print("\n--- Part 19: Qudit Shot gate count (2-molecule) (check 24) ---")

    sim_shot2_meta = QuditGKSLShotSimulator(params2)
    result_shot2_meta = sim_shot2_meta.simulate(
        t_max=5.0, n_steps=5, initial_state="edge_triplet",
        n_shots=10, seed=0,
    )
    actual_shot_gates_2 = result_shot2_meta["estimated_gates_per_step"]
    shot_gates_ok_2 = actual_shot_gates_2 == expected_gates_2

    checks["check24_qudit_shot_gate_count_2mol"] = {
        "actual_gates_per_step": actual_shot_gates_2,
        "expected_gates_per_step": expected_gates_2,
        "ok": shot_gates_ok_2,
    }
    print(f"  shot gates_per_step={actual_shot_gates_2}, expected={expected_gates_2}")
    print(f"  Check 24 - Qudit Shot gates (2-mol): {shot_gates_ok_2}")

    # ===================================================================
    # Part 20: NEW — Qubit Shot includes gate count (check 25)
    # ===================================================================
    print("\n--- Part 20: Qubit Shot includes gate count (check 25) ---")

    has_gate_key = "estimated_gates_per_step" in result_qb_shot2
    has_total_key = "total_estimated_gates" in result_qb_shot2
    if has_gate_key:
        qb_shot_gates = result_qb_shot2["estimated_gates_per_step"]
        # Expected: n_sys_qubits + (N-1)*10 + n_ancilla*6*2
        n_sys_qb = 2 * params2.N_molecules
        n_ancilla_qb = len(sim_qb_shot2.lindblad_ops)
        expected_qb_gates = n_sys_qb + len(params2.neighbors) * 10 + n_ancilla_qb * 6 * 2
        qb_gates_match = qb_shot_gates == expected_qb_gates
    else:
        qb_shot_gates = None
        expected_qb_gates = None
        qb_gates_match = False

    qb_shot_gate_ok = has_gate_key and has_total_key and qb_gates_match
    checks["check25_qubit_shot_has_gate_count"] = {
        "has_estimated_gates_per_step": has_gate_key,
        "has_total_estimated_gates": has_total_key,
        "actual_gates_per_step": qb_shot_gates,
        "expected_gates_per_step": expected_qb_gates,
        "gates_match": qb_gates_match,
        "ok": qb_shot_gate_ok,
    }
    print(f"  has estimated_gates_per_step: {has_gate_key}")
    print(f"  has total_estimated_gates: {has_total_key}")
    if has_gate_key:
        print(f"  gates_per_step={qb_shot_gates}, expected={expected_qb_gates}")
    print(f"  Check 25 - Qubit Shot has gate count: {qb_shot_gate_ok}")

    # ===================================================================
    # Part 21: NEW — Qudit-Qubit Shot identical with same seed (check 26)
    # ===================================================================
    print("\n--- Part 21: Qudit-Qubit Shot identical results (2-mol) (check 26) ---")

    td_qudit_qubit_shot = trace_distance(rho_shot2, rho_qb_shot2)
    identical_ok = bool(td_qudit_qubit_shot < 1e-14)

    checks["check26_qudit_qubit_shot_identical"] = {
        "trace_distance": td_qudit_qubit_shot,
        "ok": identical_ok,
        "note": "Qubit and Qudit shot simulators produce identical results "
                "with same seed because physical subspace evolution is identical",
    }
    print(f"  T(qudit_shot, qubit_shot): {td_qudit_qubit_shot:.2e}")
    print(f"  Check 26 - Qudit-Qubit Shot identical: {identical_ok}")

    # ===================================================================
    # Part 22: NEW — Boson Qudit convergence test (check 27)
    # ===================================================================
    print("\n--- Part 22: Boson Qudit convergence test (check 27) ---")

    from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator
    from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
    from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator

    params_boson = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
        N_molecules=2,
        with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02,
    )
    t_max_boson = 10.0

    sim_boson_exact = ClassicalGKSLBosonSimulator(params_boson)
    result_boson_exact = sim_boson_exact.simulate(
        t_max=t_max_boson, n_steps=1000, initial_state="edge_triplet",
    )
    rho_boson_exact = result_boson_exact["rho_final"]

    n_steps_boson = [5, 10, 20, 50]
    conv_boson_qudit = []
    rho_qudit_boson_20 = None
    for ns in n_steps_boson:
        dt_b = t_max_boson / ns
        sim_qd_boson = QuditGKSLBosonSimulator(params_boson)
        result_qd_boson = sim_qd_boson.simulate(
            t_max=t_max_boson, n_steps=ns, initial_state="edge_triplet",
        )
        rho_qd_boson = result_qd_boson["rho_final"]
        td_b = trace_distance(rho_boson_exact, rho_qd_boson)
        conv_boson_qudit.append({"n_steps": ns, "dt": dt_b, "trace_distance": td_b})
        print(f"  Boson Qudit: n_steps={ns:>4}, dt={dt_b:>8.4f}, T={td_b:.6e}")
        if ns == 20:
            rho_qudit_boson_20 = rho_qd_boson

    rates_boson_qudit = []
    for i in range(1, len(conv_boson_qudit)):
        prev = conv_boson_qudit[i - 1]
        curr = conv_boson_qudit[i]
        ratio_T = prev["trace_distance"] / curr["trace_distance"]
        ratio_dt = prev["dt"] / curr["dt"]
        if ratio_dt > 1 and ratio_T > 0:
            rate_b = float(np.log(ratio_T) / np.log(ratio_dt))
        else:
            rate_b = float("nan")
        conv_boson_qudit[i]["rate_T"] = rate_b
        rates_boson_qudit.append(rate_b)
        print(f"  Rate(T) [{prev['n_steps']}->{curr['n_steps']}]: {rate_b:.4f}")

    last_rate_boson_qudit = rates_boson_qudit[-1] if rates_boson_qudit else float("nan")
    boson_qudit_ok = bool(abs(last_rate_boson_qudit - 1.0) < 0.3)
    checks["check27_boson_qudit_convergence"] = {
        "convergence_data": conv_boson_qudit,
        "last_rate": last_rate_boson_qudit,
        "ok": boson_qudit_ok,
    }
    print(f"  Last Rate(T) [Boson Qudit]: {last_rate_boson_qudit:.4f}")
    print(f"  Check 27 - Boson Qudit convergence: {boson_qudit_ok}")

    # ===================================================================
    # Part 23: NEW — Boson Qubit convergence test (check 28)
    # ===================================================================
    print("\n--- Part 23: Boson Qubit convergence test (check 28) ---")

    conv_boson_qubit = []
    rho_qubit_boson_20 = None
    for ns in n_steps_boson:
        dt_b = t_max_boson / ns
        sim_qb_boson = QubitGKSLBosonSimulator(params_boson)
        result_qb_boson = sim_qb_boson.simulate(
            t_max=t_max_boson, n_steps=ns, initial_state="edge_triplet",
        )
        rho_qb_boson = result_qb_boson["rho_final"]
        td_b = trace_distance(rho_boson_exact, rho_qb_boson)
        conv_boson_qubit.append({"n_steps": ns, "dt": dt_b, "trace_distance": td_b})
        print(f"  Boson Qubit: n_steps={ns:>4}, dt={dt_b:>8.4f}, T={td_b:.6e}")
        if ns == 20:
            rho_qubit_boson_20 = rho_qb_boson

    rates_boson_qubit = []
    for i in range(1, len(conv_boson_qubit)):
        prev = conv_boson_qubit[i - 1]
        curr = conv_boson_qubit[i]
        ratio_T = prev["trace_distance"] / curr["trace_distance"]
        ratio_dt = prev["dt"] / curr["dt"]
        if ratio_dt > 1 and ratio_T > 0:
            rate_b = float(np.log(ratio_T) / np.log(ratio_dt))
        else:
            rate_b = float("nan")
        conv_boson_qubit[i]["rate_T"] = rate_b
        rates_boson_qubit.append(rate_b)
        print(f"  Rate(T) [{prev['n_steps']}->{curr['n_steps']}]: {rate_b:.4f}")

    last_rate_boson_qubit = rates_boson_qubit[-1] if rates_boson_qubit else float("nan")
    boson_qubit_ok = bool(abs(last_rate_boson_qubit - 1.0) < 0.3)
    checks["check28_boson_qubit_convergence"] = {
        "convergence_data": conv_boson_qubit,
        "last_rate": last_rate_boson_qubit,
        "ok": boson_qubit_ok,
    }
    print(f"  Last Rate(T) [Boson Qubit]: {last_rate_boson_qubit:.4f}")
    print(f"  Check 28 - Boson Qubit convergence: {boson_qubit_ok}")

    # ===================================================================
    # Part 24: NEW — Boson Qudit-Qubit consistency (check 29)
    # ===================================================================
    print("\n--- Part 24: Boson Qudit-Qubit consistency (check 29) ---")

    assert rho_qudit_boson_20 is not None
    assert rho_qubit_boson_20 is not None
    td_boson_qd_qb = trace_distance(rho_qudit_boson_20, rho_qubit_boson_20)
    boson_qd_qb_ok = bool(td_boson_qd_qb < 1e-14)

    checks["check29_boson_qudit_qubit_consistency"] = {
        "trace_distance": td_boson_qd_qb,
        "n_steps": 20,
        "t_max": t_max_boson,
        "ok": boson_qd_qb_ok,
        "note": "Qudit and Qubit boson simulators use the same physical subspace",
    }
    print(f"  T(qudit_boson, qubit_boson) [n_steps=20]: {td_boson_qd_qb:.2e}")
    print(f"  Check 29 - Boson Qudit-Qubit consistency: {boson_qd_qb_ok}")

    # ===================================================================
    # Part 25: NEW — Boson no-noise baseline (check 30)
    # ===================================================================
    print("\n--- Part 25: Boson no-noise baseline (check 30) ---")

    quality_boson_qudit = check_density_matrix_quality(rho_qudit_boson_20)
    boson_trace_ok = quality_boson_qudit["trace_error"] < 1e-10
    boson_herm_ok = quality_boson_qudit["hermiticity_error"] < 1e-10
    boson_pos_ok = quality_boson_qudit["min_eigenvalue"] > -1e-10
    boson_baseline_ok = bool(boson_trace_ok and boson_herm_ok and boson_pos_ok)

    checks["check30_boson_no_noise_baseline"] = {
        "trace": quality_boson_qudit["trace"],
        "trace_error": quality_boson_qudit["trace_error"],
        "hermiticity_error": quality_boson_qudit["hermiticity_error"],
        "min_eigenvalue": quality_boson_qudit["min_eigenvalue"],
        "ok": boson_baseline_ok,
    }
    print(f"  Trace: {quality_boson_qudit['trace']:.10f}")
    print(f"  Hermiticity error: {quality_boson_qudit['hermiticity_error']:.2e}")
    print(f"  Min eigenvalue: {quality_boson_qudit['min_eigenvalue']:.2e}")
    print(f"  Check 30 - Boson no-noise baseline: {boson_baseline_ok}")

    # ===================================================================
    # Part 26: Iteration count verification (check 31 - removed)
    # ===================================================================
    print("\n--- Part 26: Iteration count verification (check 31 - REMOVED) ---")

    # Self-referential iteration counter check removed.
    # See developing/検証結果分析_iteration45_根本原因分析.md for details.
    cell_stable_ok = True

    checks["check31_iteration_count_removed"] = {
        "note": "Self-referential iteration counter check removed to break infinite loop",
        "ok": cell_stable_ok,
    }
    print(f"  Check 31 - Iteration count check removed (was self-referential): {cell_stable_ok}")

    # ===================================================================
    # Part 27: NEW — Non-boson simulators no _U_stines attribute (check 32)
    # ===================================================================
    print("\n--- Part 27: Non-boson simulators no _U_stines attribute (check 32) ---")

    stines_results = {}

    sim_qd_check = QuditGKSLSimulator(params2)
    sim_qd_check.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
    has_stines_qd = hasattr(sim_qd_check, "_U_stines")
    stines_results["QuditGKSLSimulator"] = not has_stines_qd
    print(f"  QuditGKSLSimulator has _U_stines: {has_stines_qd} (want False)")

    sim_qb_check = QubitGKSLSimulator(params2)
    sim_qb_check.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
    has_stines_qb = hasattr(sim_qb_check, "_U_stines")
    stines_results["QubitGKSLSimulator"] = not has_stines_qb
    print(f"  QubitGKSLSimulator has _U_stines: {has_stines_qb} (want False)")

    sim_qd_shot_check = QuditGKSLShotSimulator(params2)
    sim_qd_shot_check.simulate(
        t_max=5.0, n_steps=5, initial_state="edge_triplet",
        n_shots=10, seed=0,
    )
    has_stines_qd_shot = hasattr(sim_qd_shot_check, "_U_stines")
    stines_results["QuditGKSLShotSimulator"] = not has_stines_qd_shot
    print(f"  QuditGKSLShotSimulator has _U_stines: {has_stines_qd_shot} (want False)")

    sim_qb_shot_check = QubitGKSLShotSimulator(params2)
    sim_qb_shot_check.simulate(
        t_max=5.0, n_steps=5, initial_state="edge_triplet",
        n_shots=10, seed=0,
    )
    has_stines_qb_shot = hasattr(sim_qb_shot_check, "_U_stines")
    stines_results["QubitGKSLShotSimulator"] = not has_stines_qb_shot
    print(f"  QubitGKSLShotSimulator has _U_stines: {has_stines_qb_shot} (want False)")

    no_stines_ok = bool(all(stines_results.values()))
    checks["check32_no_U_stines_attribute"] = {
        "simulator_results": stines_results,
        "ok": no_stines_ok,
    }
    print(f"  Check 32 - No _U_stines in non-boson simulators: {no_stines_ok}")

    # ===================================================================
    # Part 28: NEW — Stinespring error coefficient stability (check 33)
    # ===================================================================
    print("\n--- Part 28: Stinespring error coefficient stability (check 33) ---")

    # Verify that T/dt converges to a constant (confirming O(dt) convergence)
    coefficients = []
    for d in convergence_data:
        coeff = d["trace_distance"] / d["dt"]
        coefficients.append(coeff)
    # Last two coefficients should be within 2% of each other
    if len(coefficients) >= 2:
        last = coefficients[-1]
        second_last = coefficients[-2]
        coeff_ratio = abs(last - second_last) / last if last > 0 else float("inf")
        coeff_stable = bool(coeff_ratio < 0.02)
    else:
        coeff_ratio = float("nan")
        coeff_stable = False

    checks["check33_error_coefficient_stability"] = {
        "coefficients": [{"n_steps": convergence_data[i]["n_steps"], "T_over_dt": coefficients[i]}
                         for i in range(len(coefficients))],
        "last_two_ratio": coeff_ratio,
        "ok": coeff_stable,
    }
    print(f"  T/dt values: {[f'{c:.6e}' for c in coefficients]}")
    print(f"  Last two T/dt ratio deviation: {coeff_ratio:.4e}")
    print(f"  Check 33 - Error coefficient stability: {coeff_stable}")

    # ===================================================================
    # Part 30: NEW — Stinespring single-channel error O(dt²) (check 34)
    # ===================================================================
    print("\n--- Part 30: Stinespring single-channel error O(dt²) (check 34) ---")

    from stinespring_utils import (
        apply_stinespring_to_density_matrix as apply_stine,
        build_gksl_superoperator,
        stinespring_unitary_from_lindblad,
    )

    # Use first Lindblad operator from 2-molecule system
    L_test = sim_qd2.lindblad_ops[0][0]  # L (includes sqrt(gamma))
    d_test = L_test.shape[0]  # 9 for 2-molecule qutrit

    # Build exact single-channel dissipator superoperator
    L_D_single = build_gksl_superoperator(
        np.zeros((d_test, d_test), dtype=np.complex128),  # zero Hamiltonian
        [(L_test, 0.0)],  # single Lindblad channel
    )

    # Test state: the initial edge_triplet state
    rho_test = sim_qd2.prepare_initial_state("edge_triplet")

    # Measure Stinespring error vs exact channel for different dt values
    from gksl_math_utils import vectorize_density_matrix, unvectorize_density_matrix

    dt_values = [1.0, 0.1, 0.01, 0.001]
    stine_errors = []
    for dt_sc in dt_values:
        # Exact: exp(L_D dt) rho
        vec_rho = vectorize_density_matrix(rho_test)
        vec_exact = expm_multiply(L_D_single, vec_rho, start=0.0, stop=dt_sc, num=2, endpoint=True)[-1]
        rho_exact_sc = unvectorize_density_matrix(vec_exact, d_test)

        # Stinespring: E(rho)
        U_stine = stinespring_unitary_from_lindblad(L_test, dt_sc)
        rho_stine = apply_stine(rho_test, U_stine)

        err = trace_distance(rho_exact_sc, rho_stine)
        stine_errors.append({"dt": dt_sc, "error": err})
        print(f"  dt={dt_sc:.4f}: ||E(rho) - exp(D dt) rho|| = {err:.6e}")

    # Check O(dt²) convergence: error/dt² should approach a constant
    # Only use dt ≤ 0.1 for rate check: larger dt values have significant
    # higher-order corrections that distort the O(dt²) rate estimate.
    if len(stine_errors) >= 2:
        ratios_dt2 = [e["error"] / (e["dt"] ** 2) for e in stine_errors if e["dt"] <= 0.1]
        if len(ratios_dt2) >= 2:
            last_ratio = ratios_dt2[-1]
            second_last_ratio = ratios_dt2[-2]
            ratio_deviation = abs(last_ratio - second_last_ratio) / last_ratio if last_ratio > 0 else float("inf")
            stine_channel_ok = bool(ratio_deviation < 0.05)  # 5% tolerance
        else:
            ratio_deviation = float("nan")
            stine_channel_ok = False
    else:
        ratio_deviation = float("nan")
        stine_channel_ok = False

    checks["check34_stinespring_channel_error"] = {
        "errors": stine_errors,
        "error_over_dt2": [{"dt": e["dt"], "E_over_dt2": e["error"] / (e["dt"] ** 2)}
                           for e in stine_errors],
        "ratio_deviation": ratio_deviation,
        "ok": stine_channel_ok,
    }
    e_dt2_values = [f"{e['error']/(e['dt']**2):.6e}" for e in stine_errors]
    print(f"  E/dt² values: {e_dt2_values}")
    print(f"  Ratio deviation (last two, dt≤0.1): {ratio_deviation:.4e}")
    print(f"  Check 34 - Stinespring channel error O(dt²): {stine_channel_ok}")

    results["checks"] = checks

    # ===================================================================
    # Overall verdict
    # ===================================================================
    all_pass = all([
        monotonic, rate_near_1, fid_monotonic, all_quality_ok,
        practical_accuracy, cell27_ok,
        doc_fixes["cell25_ok"], doc_fixes["cell39_ok"], doc_fixes["cell27_n200_ok"],
        pnum_ok, interp_ok, conv2_ok, unitary_ok, eval_ok,
        shot_dm_ok, qb_shot_dm_ok, noisy_zero_ok, qb_noisy_zero_ok,
        ancilla_match_4 and ancilla_expected_4,
        ancilla_match_2 and ancilla_expected_2,
        gates_ok_4, gates_ok_2,
        shot_gates_ok_4, shot_gates_ok_2,
        qb_shot_gate_ok,
        identical_ok,
        # Checks 27-32
        boson_qudit_ok,
        boson_qubit_ok,
        boson_qd_qb_ok,
        boson_baseline_ok,
        cell_stable_ok,
        no_stines_ok,
        # Check 33-34
        coeff_stable,
        stine_channel_ok,
    ])
    results["all_checks_passed"] = all_pass
    print(f"\n{'=' * 70}")
    print(f"ALL CHECKS PASSED: {all_pass}")
    print(f"{'=' * 70}")

    # ===================================================================
    # Save results
    # ===================================================================
    json_path = results_dir / f"stable_regression_{timestamp}.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nJSON results saved to: {json_path}")

    md_path = results_dir / f"stable_regression_{timestamp}.md"
    md_lines = [
        f"# TTA-UC GKSL 安定回帰テスト結果\n",
        f"\n",
        f"- 実行日時: {timestamp}\n",
        f"- 検証対象: 全チェック（自己参照的イテレーションカウンタチェック除去済み）\n",
        f"\n",
        f"## 収束解析結果（4分子系）\n",
        f"\n",
        f"| n_steps | dt | T(ST,exact) | Fidelity | Rate(T) |\n",
        f"|--------:|-------:|--------------:|----------:|--------:|\n",
    ]
    for d in convergence_data:
        rate_str = f"{d.get('rate_T', float('nan')):.4f}" if "rate_T" in d else "---"
        md_lines.append(
            f"| {d['n_steps']:>7} | {d['dt']:>8.4f} | {d['trace_distance']:.6e} "
            f"| {d['fidelity']:.8f} | {rate_str} |\n"
        )
    md_lines.extend([
        f"\n",
        f"## 全検証チェック\n",
        f"\n",
        f"| # | チェック項目 | 結果 |\n",
        f"|---|---|---|\n",
        f"| 1 | トレース距離単調減少 | {'✓' if monotonic else '✗'} |\n",
        f"| 2 | 収束次数 ≈ 1.0 (平均: {avg_rate:.4f}) | {'✓' if rate_near_1 else '✗'} |\n",
        f"| 3 | 忠実度単調増加 | {'✓' if fid_monotonic else '✗'} |\n",
        f"| 4 | 密度行列品質 | {'✓' if all_quality_ok else '✗'} |\n",
        f"| 5 | n_steps=100 で T < 1e-3 | {'✓' if practical_accuracy else '✗'} |\n",
        f"| 6 | Cell 27 テーブル値の整合性 | {'✓' if cell27_ok else '✗'} |\n",
        f"| 7 | Cell 25 イテレーションカウンタ（チェック除去済み） | ✓ |\n",
        f"| 8 | Cell 39 イテレーションカウンタ（チェック除去済み） | ✓ |\n",
        f"| 9 | Cell 27 n_steps=200 値 (~3.6e-04) | {'✓' if doc_fixes['cell27_n200_ok'] else '✗'} |\n",
        f"| 10 | 粒子数保存 (Cl: {max_pnum_err_cl:.2e}, Qd: {max_pnum_err_qd:.2e}) | {'✓' if pnum_ok else '✗'} |\n",
        f"| 11 | 中間時刻人口動態 (max diff: {max_pop_diff:.4e}) | {'✓' if interp_ok else '✗'} |\n",
        f"| 12 | 2分子系収束次数 (Rate: {last_rate_2mol:.4f}) | {'✓' if conv2_ok else '✗'} |\n",
        f"| 13 | ユニタリ限界精度 (T: {td_u:.2e}) | {'✓' if unitary_ok else '✗'} |\n",
        f"| 14 | 固有値スペクトル一致 (max diff: {max_eval_diff:.6e}) | {'✓' if eval_ok else '✗'} |\n",
        f"| 15 | Qudit Shot-DM一致 (T: {td_shot_dm:.4e}, F: {fid_shot_dm:.6f}) | {'✓' if shot_dm_ok else '✗'} |\n",
        f"| 16 | Qubit Shot-DM一致 (T: {td_qb_shot_dm:.4e}, F: {fid_qb_shot_dm:.6f}) | {'✓' if qb_shot_dm_ok else '✗'} |\n",
        f"| 17 | Qudit Noisy(p=0)=DM (T: {td_noisy_zero:.2e}) | {'✓' if noisy_zero_ok else '✗'} |\n",
        f"| 18 | Qubit Noisy(p=0)=DM (T: {td_qb_noisy_zero:.2e}) | {'✓' if qb_noisy_zero_ok else '✗'} |\n",
        f"| 19 | Qudit DM アンシラ数 (4分子: {sim_qd4.n_ancilla_qubits}={expected_4}) | {'✓' if ancilla_match_4 and ancilla_expected_4 else '✗'} |\n",
        f"| 20 | Qudit DM アンシラ数 (2分子: {sim_qd2.n_ancilla_qubits}={expected_2}) | {'✓' if ancilla_match_2 and ancilla_expected_2 else '✗'} |\n",
        f"| 21 | Qudit DM ゲート数 (4分子: {actual_gates_4}={expected_gates_4}) | {'✓' if gates_ok_4 else '✗'} |\n",
        f"| 22 | Qudit DM ゲート数 (2分子: {actual_gates_2}={expected_gates_2}) | {'✓' if gates_ok_2 else '✗'} |\n",
        f"| 23 | Qudit Shot ゲート数 (4分子: {actual_shot_gates_4}={expected_gates_4}) | {'✓' if shot_gates_ok_4 else '✗'} |\n",
        f"| 24 | Qudit Shot ゲート数 (2分子: {actual_shot_gates_2}={expected_gates_2}) | {'✓' if shot_gates_ok_2 else '✗'} |\n",
        f"| 25 | Qubit Shot ゲート数推定あり | {'✓' if qb_shot_gate_ok else '✗'} |\n",
        f"| 26 | Qudit-Qubit Shot 同一結果 (T: {td_qudit_qubit_shot:.2e}) | {'✓' if identical_ok else '✗'} |\n",
        f"| 27 | Boson Qudit 収束次数 (Rate: {last_rate_boson_qudit:.4f}) | {'✓' if boson_qudit_ok else '✗'} |\n",
        f"| 28 | Boson Qubit 収束次数 (Rate: {last_rate_boson_qubit:.4f}) | {'✓' if boson_qubit_ok else '✗'} |\n",
        f"| 29 | Boson Qudit-Qubit一致 (T: {td_boson_qd_qb:.2e}) | {'✓' if boson_qd_qb_ok else '✗'} |\n",
        f"| 30 | Boson 密度行列品質 (tr_err: {quality_boson_qudit['trace_error']:.2e}) | {'✓' if boson_baseline_ok else '✗'} |\n",
        f"| 31 | イテレーションカウンタ（チェック除去済み - 自己参照ループ防止） | ✓ |\n",
        f"| 32 | 非Bosonシミュレータに_U_stinesなし | {'✓' if no_stines_ok else '✗'} |\n",
        f"| 33 | 誤差係数 T/dt の安定性 (偏差: {coeff_ratio:.4e}) | {'✓' if coeff_stable else '✗'} |\n",
        f"| 34 | Stinespring チャネル誤差 O(dt²) (偏差: {ratio_deviation:.4e}) | {'✓' if stine_channel_ok else '✗'} |\n",
        f"\n",
        f"## 総合判定: {'PASS ✓' if all_pass else 'FAIL ✗'}\n",
    ])
    with open(md_path, "w") as f:
        f.writelines(md_lines)
    print(f"Markdown report saved to: {md_path}")


if __name__ == "__main__":
    main()
