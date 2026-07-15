#!/usr/bin/env python3
"""TTA-UC GKSL iteration 33 verification: Cell 27 table consistency validation.

This iteration validates that Cell 27 precision guidance table values in
quantum_dynamics_gksl_comparison.ipynb are consistent with the actual
computed values from Cell 26.

Changes from iteration 32:
  - Cell 27 markdown: trace distance values updated to match t_max=100.0 actual output
  - Cell 27 markdown: n_steps=200 row added
  - Cell 27 markdown: incorrect "T < 2e-5" statement corrected to "T < 1e-3"
  - Cell 27 markdown: t_max dependency explanation added
  - iteration31 script: markdown report label corrected (5e-5 → 1e-3)

Verification targets:
  1. Trace distance T(ST, exact) decreases with increasing n_steps
  2. Convergence rate Rate(T) ≈ 1.0 (O(dt) from Stinespring dilation)
  3. Density matrix quality preserved at all n_steps
  4. Fidelity F(ST, exact) increases with increasing n_steps
  5. n_steps=100 (dt=1.0) gives T < 1e-3 (threshold for t_max=100.0)
  6. Cell 27 table values match actual computed values within tolerance

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


def check_cell27_table_consistency(convergence_data: list[dict]) -> dict:
    """Check that Cell 27 table values match actual computed values.

    Cell 27 expected table values (t_max=100.0):
      n_steps=10:  ~1.3e-02
      n_steps=20:  ~4e-03
      n_steps=50:  ~1.5e-03
      n_steps=100: ~7e-04
      n_steps=200: ~4e-04
    """
    cell27_expected = {
        10: 1.3e-02,
        20: 4e-03,
        50: 1.5e-03,
        100: 7e-04,
        200: 4e-04,
    }

    results = {}
    all_ok = True
    for d in convergence_data:
        ns = d["n_steps"]
        if ns in cell27_expected:
            expected = cell27_expected[ns]
            actual = d["trace_distance"]
            # Allow ±30% relative tolerance for "~" approximate values in Cell 27.
            # Actual ratios observed: 0.91–1.05, so 0.7–1.3 provides adequate margin.
            ratio = actual / expected
            ok = 0.7 < ratio < 1.3
            results[f"n_steps_{ns}"] = {
                "cell27_value": expected,
                "actual_value": actual,
                "ratio": ratio,
                "within_tolerance": ok,
            }
            if not ok:
                all_ok = False

    results["all_consistent"] = all_ok
    return results


def main() -> None:
    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    results_dir = Path(__file__).resolve().parent.parent / "developing" / "verification_results"
    results_dir.mkdir(parents=True, exist_ok=True)

    params = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
    )
    t_max = 100.0

    results: dict = {
        "iteration": 33,
        "timestamp": timestamp,
        "description": "Cell 27 table consistency validation",
        "params": {
            "t_max": t_max,
            "hilbert_dim": params.get_hilbert_space_dim(),
        },
    }

    print("=" * 70)
    print("Iteration 33: Cell 27 Table Consistency Validation")
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
    print(f"  Hermiticity error: {quality_exact['hermiticity_error']:.2e}")
    print(f"  Min eigenvalue: {quality_exact['min_eigenvalue']:.2e}")
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
        max_trace_dev = max(abs(t - 1.0) for t in result_st["trace"])

        entry = {
            "n_steps": ns,
            "dt": dt,
            "trace_distance": td,
            "fidelity": fid,
            "elapsed_time": elapsed,
            "quality": quality,
            "max_trace_deviation": max_trace_dev,
        }
        convergence_data.append(entry)
        print(f"  n_steps={ns:>4}, dt={dt:>8.4f}, T={td:.6e}, F={fid:.8f}, time={elapsed:.2f}s")

    # Compute convergence rates
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

    # Part 3: Validation checks
    print("\n--- Part 3: Validation checks ---")
    checks = {}

    # Check 1: Monotonic decrease of trace distance
    tds = [d["trace_distance"] for d in convergence_data]
    monotonic = all(tds[i] > tds[i + 1] for i in range(len(tds) - 1))
    checks["trace_distance_monotonic_decrease"] = monotonic
    print(f"  Check 1 - Trace distance monotonically decreasing: {monotonic}")

    # Check 2: Convergence rate ≈ 1.0
    rates = [d.get("rate_T", float("nan")) for d in convergence_data if "rate_T" in d]
    avg_rate = float(np.mean(rates)) if rates else float("nan")
    rate_near_1 = abs(avg_rate - 1.0) < 0.3  # Allow ±0.3 tolerance
    checks["convergence_rate_near_1"] = rate_near_1
    checks["average_rate"] = avg_rate
    print(f"  Check 2 - Average Rate(T): {avg_rate:.4f} (expected ≈ 1.0, tolerance ±0.3)")
    print(f"  Rate near 1.0: {rate_near_1}")

    # Check 3: Monotonic increase of fidelity
    fids = [d["fidelity"] for d in convergence_data]
    fid_monotonic = all(fids[i] < fids[i + 1] for i in range(len(fids) - 1))
    checks["fidelity_monotonic_increase"] = fid_monotonic
    print(f"  Check 3 - Fidelity monotonically increasing: {fid_monotonic}")

    # Check 4: Density matrix quality at all n_steps
    all_quality_ok = all(
        d["quality"]["trace_error"] < 1e-10
        and d["quality"]["hermiticity_error"] < 1e-10
        and d["quality"]["min_eigenvalue"] > -1e-10
        for d in convergence_data
    )
    checks["density_matrix_quality_all_ok"] = all_quality_ok
    print(f"  Check 4 - Density matrix quality OK at all n_steps: {all_quality_ok}")

    # Check 5: n_steps=100 gives T < 1e-3
    t100 = next(d["trace_distance"] for d in convergence_data if d["n_steps"] == 100)
    practical_accuracy = t100 < 1e-3
    checks["n100_trace_distance"] = t100
    checks["n100_practical_accuracy"] = practical_accuracy
    print(f"  Check 5 - n_steps=100 trace distance: {t100:.6e} (< 1e-3: {practical_accuracy})")

    # Check 6: Cell 27 table values match actual computed values
    print("\n--- Part 4: Cell 27 table consistency ---")
    cell27_check = check_cell27_table_consistency(convergence_data)
    cell27_consistent = cell27_check["all_consistent"]
    checks["cell27_table_consistency"] = cell27_check
    for ns in [10, 20, 50, 100, 200]:
        key = f"n_steps_{ns}"
        if key in cell27_check:
            info = cell27_check[key]
            status = "✓" if info["within_tolerance"] else "✗"
            print(f"  {status} n_steps={ns}: Cell 27={info['cell27_value']:.1e}, "
                  f"actual={info['actual_value']:.4e}, ratio={info['ratio']:.3f}")

    results["checks"] = checks

    # Overall verdict
    all_pass = all([
        monotonic, rate_near_1, fid_monotonic, all_quality_ok,
        practical_accuracy, cell27_consistent,
    ])
    results["all_checks_passed"] = all_pass
    print(f"\n{'=' * 70}")
    print(f"ALL CHECKS PASSED: {all_pass}")
    print(f"{'=' * 70}")

    # Save results
    json_path = results_dir / f"iteration33_cell27_consistency_{timestamp}.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nJSON results saved to: {json_path}")

    # Save markdown report
    md_path = results_dir / f"iteration33_cell27_consistency_{timestamp}.md"
    md_lines = [
        f"# Iteration 33 検証結果: Cell 27 テーブル整合性検証\n",
        f"\n",
        f"- 実行日時: {timestamp}\n",
        f"- 検証対象: quantum_dynamics_gksl_comparison.ipynb Cell 27 修正後の整合性\n",
        f"- 変更点: Cell 27 精度ガイダンステーブルの値を t_max=100.0 の実測値に修正\n",
        f"- 変更点: 誤った記述「T < 2×10⁻⁵」を「T < 1×10⁻³」に修正\n",
        f"- 変更点: n_steps=200 の行を追加\n",
        f"- 変更点: t_max 依存性の説明を追加\n",
        f"- 変更点: iteration31 スクリプトの Markdown ラベル修正（5e-5 → 1e-3）\n",
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
        f"| チェック項目 | 結果 |\n",
        f"|---|---|\n",
        f"| トレース距離単調減少 | {'✓' if monotonic else '✗'} |\n",
        f"| 収束次数 ≈ 1.0 (平均: {avg_rate:.4f}) | {'✓' if rate_near_1 else '✗'} |\n",
        f"| 忠実度単調増加 | {'✓' if fid_monotonic else '✗'} |\n",
        f"| 密度行列品質 (全 n_steps) | {'✓' if all_quality_ok else '✗'} |\n",
        f"| n_steps=100 で T < 1e-3 | {'✓' if practical_accuracy else '✗'} |\n",
        f"| Cell 27 テーブル整合性 | {'✓' if cell27_consistent else '✗'} |\n",
        f"\n",
        f"## Cell 27 テーブル整合性の詳細\n",
        f"\n",
        f"| n_steps | Cell 27 記載値 | 実測値 | 比率 | 判定 |\n",
        f"|--------:|:--------------:|:---------:|:-----:|:----:|\n",
    ])
    for ns in [10, 20, 50, 100, 200]:
        key = f"n_steps_{ns}"
        if key in cell27_check:
            info = cell27_check[key]
            status = "✓" if info["within_tolerance"] else "✗"
            md_lines.append(
                f"| {ns:>7} | {info['cell27_value']:.1e} | {info['actual_value']:.4e} "
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
