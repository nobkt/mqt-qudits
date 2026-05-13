#!/usr/bin/env python3
"""TTA-UC GKSL iteration 35 verification: confirm Cell 25/27/39 documentation fixes.

This iteration validates that the documentation fixes from the iteration 34
analysis have been correctly applied to the notebook:

  1. Cell 25: "34回のイテレーション検証" (was "30回")
  2. Cell 27: n_steps=200 value updated to "~3.6e-04" (was "~4e-04")
  3. Cell 39: "34回の検証イテレーション" (was "30回")

Additionally, all iteration 34 checks are re-run to confirm no regression.

Verification targets:
  1. Trace distance T(ST, exact) decreases with increasing n_steps
  2. Convergence rate Rate(T) ≈ 1.0 (O(dt) from Stinespring dilation)
  3. Fidelity F(ST, exact) increases with increasing n_steps
  4. Density matrix quality preserved at all n_steps
  5. n_steps=100 (dt=1.0) gives T < 1e-3 (threshold for t_max=100.0)
  6. Cell 27 table values (read from notebook file) match actual computed values
  7. Cell 25 contains "34回" (not "30回")
  8. Cell 39 contains "34回" (not "30回")
  9. Cell 27 n_steps=200 value is "~3.6e-04" (not "~4e-04")

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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_physical_parameters import GKSLPhysicalParameters
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
    """Check that documentation fixes from iteration 34 analysis are applied."""
    with open(notebook_path) as f:
        nb = json.load(f)

    results = {}

    # Check 7: Cell 25 contains "34回"
    cell25_text = "".join(nb["cells"][25]["source"])
    has_34_cell25 = "34回のイテレーション検証" in cell25_text
    has_30_cell25 = "30回のイテレーション検証" in cell25_text
    results["cell25_has_34"] = has_34_cell25
    results["cell25_has_old_30"] = has_30_cell25
    results["cell25_ok"] = has_34_cell25 and not has_30_cell25

    # Check 8: Cell 39 contains "34回"
    cell39_text = "".join(nb["cells"][39]["source"])
    has_34_cell39 = "34回の検証イテレーション" in cell39_text
    has_30_cell39 = "30回の検証イテレーション" in cell39_text
    results["cell39_has_34"] = has_34_cell39
    results["cell39_has_old_30"] = has_30_cell39
    results["cell39_ok"] = has_34_cell39 and not has_30_cell39

    # Check 9: Cell 27 n_steps=200 value is "~3.6e-04"
    cell27_text = "".join(nb["cells"][27]["source"])
    has_36e04 = "~3.6e-04" in cell27_text
    has_old_4e04 = "~4e-04" in cell27_text
    results["cell27_has_36e04"] = has_36e04
    results["cell27_has_old_4e04"] = has_old_4e04
    results["cell27_n200_ok"] = has_36e04 and not has_old_4e04

    results["all_doc_fixes_ok"] = (
        results["cell25_ok"]
        and results["cell39_ok"]
        and results["cell27_n200_ok"]
    )

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
        "iteration": 35,
        "timestamp": timestamp,
        "description": "Documentation fixes validation (Cell 25/27/39 updates)",
        "params": {
            "t_max": t_max,
            "hilbert_dim": params.get_hilbert_space_dim(),
        },
    }

    print("=" * 70)
    print("Iteration 35: Documentation Fixes Validation")
    print("=" * 70)

    # Part 1: Exact reference solution
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

    # Part 2: Convergence analysis
    print("\n--- Part 2: Convergence analysis ---")
    n_steps_list = [10, 20, 50, 100, 200]
    convergence_data = []

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

    # Part 3: Original validation checks (1-5)
    print("\n--- Part 3: Original validation checks ---")
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

    # Part 4: Cell 27 content validation (check 6)
    print("\n--- Part 4: Cell 27 notebook content validation ---")
    parsed = parse_cell27_table_from_notebook(notebook_path)
    table_values = parsed["table_values"]
    full_text = parsed["full_text"]

    missing = [ns for ns in [10, 20, 50, 100, 200] if ns not in table_values]
    all_present = len(missing) == 0
    print(f"  All n_steps present: {all_present}")

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

    # Part 5: Documentation fixes validation (checks 7-9)
    print("\n--- Part 5: Documentation fixes validation ---")
    doc_fixes = check_documentation_fixes(notebook_path)
    checks["check7_cell25_iteration_count"] = doc_fixes["cell25_ok"]
    checks["check8_cell39_iteration_count"] = doc_fixes["cell39_ok"]
    checks["check9_cell27_n200_value"] = doc_fixes["cell27_n200_ok"]

    print(f"  Check 7 - Cell 25 has '34回': {doc_fixes['cell25_ok']}")
    if not doc_fixes["cell25_ok"]:
        print(f"    has_34: {doc_fixes['cell25_has_34']}, has_old_30: {doc_fixes['cell25_has_old_30']}")
    print(f"  Check 8 - Cell 39 has '34回': {doc_fixes['cell39_ok']}")
    if not doc_fixes["cell39_ok"]:
        print(f"    has_34: {doc_fixes['cell39_has_34']}, has_old_30: {doc_fixes['cell39_has_old_30']}")
    print(f"  Check 9 - Cell 27 n_steps=200 is '~3.6e-04': {doc_fixes['cell27_n200_ok']}")
    if not doc_fixes["cell27_n200_ok"]:
        print(f"    has_36e04: {doc_fixes['cell27_has_36e04']}, has_old_4e04: {doc_fixes['cell27_has_old_4e04']}")

    results["checks"] = checks

    # Overall verdict
    all_pass = all([
        monotonic, rate_near_1, fid_monotonic, all_quality_ok,
        practical_accuracy, cell27_ok,
        doc_fixes["cell25_ok"], doc_fixes["cell39_ok"], doc_fixes["cell27_n200_ok"],
    ])
    results["all_checks_passed"] = all_pass
    print(f"\n{'=' * 70}")
    print(f"ALL CHECKS PASSED: {all_pass}")
    print(f"{'=' * 70}")

    # Save results
    json_path = results_dir / f"iteration35_docfix_validation_{timestamp}.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nJSON results saved to: {json_path}")

    md_path = results_dir / f"iteration35_docfix_validation_{timestamp}.md"
    md_lines = [
        f"# Iteration 35 検証結果: 文書修正の検証\n",
        f"\n",
        f"- 実行日時: {timestamp}\n",
        f"- 検証対象: Cell 25/27/39 の文書修正 + 全既存チェックの回帰テスト\n",
        f"\n",
        f"## 収束解析結果\n",
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
        f"## 検証チェック\n",
        f"\n",
        f"| # | チェック項目 | 結果 |\n",
        f"|---|---|---|\n",
        f"| 1 | トレース距離単調減少 | {'✓' if monotonic else '✗'} |\n",
        f"| 2 | 収束次数 ≈ 1.0 (平均: {avg_rate:.4f}) | {'✓' if rate_near_1 else '✗'} |\n",
        f"| 3 | 忠実度単調増加 | {'✓' if fid_monotonic else '✗'} |\n",
        f"| 4 | 密度行列品質 | {'✓' if all_quality_ok else '✗'} |\n",
        f"| 5 | n_steps=100 で T < 1e-3 | {'✓' if practical_accuracy else '✗'} |\n",
        f"| 6 | Cell 27 テーブル値の整合性 | {'✓' if cell27_ok else '✗'} |\n",
        f"| 7 | Cell 25 「34回」表記 | {'✓' if doc_fixes['cell25_ok'] else '✗'} |\n",
        f"| 8 | Cell 39 「34回」表記 | {'✓' if doc_fixes['cell39_ok'] else '✗'} |\n",
        f"| 9 | Cell 27 n_steps=200 値 (~3.6e-04) | {'✓' if doc_fixes['cell27_n200_ok'] else '✗'} |\n",
        f"\n",
        f"## Cell 27 テーブル値比較\n",
        f"\n",
        f"| n_steps | ノートブック記載値 | 実測値 | 比率 | 判定 |\n",
        f"|--------:|:-----------------:|:---------:|:-----:|:----:|\n",
    ])
    for ns in [10, 20, 50, 100, 200]:
        key = f"n_steps_{ns}"
        if key in comparison:
            info = comparison[key]
            status = "✓" if info["within_tolerance"] else "✗"
            md_lines.append(
                f"| {ns:>7} | {info['notebook_value']:.2e} | {info['actual_value']:.4e} "
                f"| {info['ratio']:.3f} | {status} |\n"
            )
    md_lines.extend([
        f"\n",
        f"## 総合判定: {'PASS ✓' if all_pass else 'FAIL ✗'}\n",
    ])
    with open(md_path, "w") as f:
        f.writelines(md_lines)
    print(f"Markdown report saved to: {md_path}")


if __name__ == "__main__":
    main()
