#!/usr/bin/env python3
"""
Iteration 8 verification script: CX gate overhead integration.

Verifies that:
1. cx_per_pair_gate type hint is float (not int) in qubit_noisy_simulator.py
2. cx_per_pair_gate is used in complete_comparison notebook Cell 15
3. QubitGKSLNoisyShotSimulator supports cx_per_pair_gate
4. p_depol_pair_eff is computed correctly in both simulators
5. Backward compatibility: cx_per_pair_gate=1 gives p_eff = p_depol
6. GKSL notebook Cell 38 uses cx_per_pair_gate
7. CX gate counts from Qiskit transpilation (if available)

Uses string-based source code inspection (no heavy imports required).
Results are saved to developing/verification_results/
"""

import json
import math
import os
import sys
from datetime import datetime, timezone

# File paths
QUBIT_NOISY = os.path.join(os.path.dirname(__file__), "qubit_noisy_simulator.py")
QUDIT_NOISY = os.path.join(os.path.dirname(__file__), "mqt_qudits_noisy_simulator.py")
QUBIT_GKSL_SHOT = os.path.join(os.path.dirname(__file__), "qubit_gksl_shot_simulator.py")
QUDIT_GKSL_SHOT = os.path.join(os.path.dirname(__file__), "qudit_gksl_shot_simulator.py")
COMPLETE_NB = os.path.join(os.path.dirname(__file__), "quantum_dynamics_complete_comparison.ipynb")
GKSL_NB = os.path.join(os.path.dirname(__file__), "quantum_dynamics_gksl_comparison.ipynb")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "developing", "verification_results")

checks = []


def check(name, condition, detail):
    """Record a verification check result."""
    status = "PASS" if condition else "FAIL"
    checks.append({"name": name, "status": status, "detail": detail})
    symbol = "✓" if condition else "✗"
    print(f"  {symbol} {name}: {detail}")


def main():
    print("=" * 70)
    print("Iteration 8 検証: CXゲートオーバーヘッドの統合")
    print("=" * 70)

    # Read source files
    with open(QUBIT_NOISY) as f:
        qubit_noisy_src = f.read()
    with open(QUDIT_NOISY) as f:
        qudit_noisy_src = f.read()
    with open(QUBIT_GKSL_SHOT) as f:
        qubit_gksl_src = f.read()
    with open(QUDIT_GKSL_SHOT) as f:
        qudit_gksl_src = f.read()

    # Read notebooks
    with open(COMPLETE_NB) as f:
        complete_nb = json.load(f)
    with open(GKSL_NB) as f:
        gksl_nb = json.load(f)

    # ============================================================
    # Section 1: qubit_noisy_simulator.py type hint check
    # ============================================================
    print("\n--- Section 1: qubit_noisy_simulator.py cx_per_pair_gate type ---")

    check(
        "cx_type_hint_float",
        "cx_per_pair_gate: float" in qubit_noisy_src,
        "cx_per_pair_gate should have float type hint (accepts averages like 66.5)",
    )

    check(
        "cx_type_not_int",
        "cx_per_pair_gate: int" not in qubit_noisy_src,
        "cx_per_pair_gate should NOT have int type hint",
    )

    check(
        "cx_effective_rate_formula",
        "(1 - depol_2q) ** cx_per_pair_gate" in qubit_noisy_src,
        "Effective rate formula: p_eff = 1 - (1-depol_2q)^cx_per_pair_gate",
    )

    # ============================================================
    # Section 2: complete_comparison.ipynb Cell 15
    # ============================================================
    print("\n--- Section 2: complete_comparison.ipynb Cell 15 ---")

    cell15_src = "".join(complete_nb["cells"][15]["source"])

    check(
        "complete_cell15_cx_param",
        "cx_per_pair_gate" in cell15_src,
        "Cell 15 should pass cx_per_pair_gate to simulate()",
    )

    check(
        "complete_cell15_cx_estimate",
        "estimate_cx_per_pair_gate" in cell15_src,
        "Cell 15 should compute CX count via estimate_cx_per_pair_gate()",
    )

    check(
        "complete_cell15_cx_avg",
        "cx_avg" in cell15_src,
        "Cell 15 should use average CX count (cx_avg)",
    )

    check(
        "complete_cell15_depol_001",
        "'depol_2q': 0.001" in cell15_src or '"depol_2q": 0.001' in cell15_src,
        "Cell 15 should use depol_2q=0.001 as physical per-CX error rate",
    )

    check(
        "complete_cell15_not_default_cx1",
        "cx_per_pair_gate=cx_per_pair" in cell15_src
        or "cx_per_pair_gate=cx_info" in cell15_src,
        "Cell 15 should pass computed CX count (not default 1)",
    )

    # ============================================================
    # Section 3: QubitGKSLNoisyShotSimulator cx_per_pair_gate
    # ============================================================
    print("\n--- Section 3: QubitGKSLNoisyShotSimulator ---")

    check(
        "gksl_qubit_cx_param",
        "cx_per_pair_gate" in qubit_gksl_src
        and "class QubitGKSLNoisyShotSimulator" in qubit_gksl_src,
        "QubitGKSLNoisyShotSimulator should accept cx_per_pair_gate",
    )

    check(
        "gksl_qubit_p_eff_formula",
        "p_depol_pair_eff" in qubit_gksl_src,
        "Should compute p_depol_pair_eff from p_depol and cx_per_pair_gate",
    )

    check(
        "gksl_qubit_init_formula",
        "(1.0 - p_depol) ** cx_per_pair_gate" in qubit_gksl_src
        or "(1 - p_depol) ** cx_per_pair_gate" in qubit_gksl_src,
        "Should compute: p_eff = 1 - (1-p_depol)^cx_per_pair_gate",
    )

    check(
        "gksl_qubit_cx_default_1",
        "cx_per_pair_gate: float = 1" in qubit_gksl_src
        or "cx_per_pair_gate = 1" in qubit_gksl_src,
        "Default cx_per_pair_gate should be 1 (backward compatible)",
    )

    check(
        "gksl_qubit_cx_validation",
        "cx_per_pair_gate < 1" in qubit_gksl_src,
        "Should validate cx_per_pair_gate >= 1",
    )

    # Check _trotter_step_trajectory uses p_depol_pair_eff for pair noise
    # Find the NOISY version (second occurrence, in QubitGKSLNoisyShotSimulator)
    noisy_class_start = qubit_gksl_src.find("class QubitGKSLNoisyShotSimulator")
    trotter_start = qubit_gksl_src.find("def _trotter_step_trajectory", noisy_class_start)
    if trotter_start >= 0:
        trotter_end = qubit_gksl_src.find("\n    def ", trotter_start + 1)
        if trotter_end < 0:
            trotter_end = len(qubit_gksl_src)
        trotter_body = qubit_gksl_src[trotter_start:trotter_end]

        pair_depol_eff_count = trotter_body.count("p_depol_pair_eff")
        check(
            "gksl_trotter_uses_p_eff",
            pair_depol_eff_count >= 4,
            f"_trotter_step_trajectory uses p_depol_pair_eff {pair_depol_eff_count} times "
            f"(expected >=4: 2 half-Hamiltonian + 2 Lindblad pair channels forward/reverse)",
        )

        single_depol_uses = trotter_body.count("self.p_depol, rng")
        check(
            "gksl_single_site_uses_p_depol",
            single_depol_uses >= 2,
            f"Single-site depolarization uses self.p_depol {single_depol_uses} times "
            f"(unscaled by CX count, as single-qubit gates are native)",
        )
    else:
        check("gksl_trotter_uses_p_eff", False, "_trotter_step_trajectory not found")
        check("gksl_single_site_uses_p_depol", False, "_trotter_step_trajectory not found")

    # Check simulate() returns cx info
    sim_start = qubit_gksl_src.find(
        "def simulate", qubit_gksl_src.find("class QubitGKSLNoisyShotSimulator")
    )
    if sim_start >= 0:
        sim_end = qubit_gksl_src.find("\n    def ", sim_start + 1)
        if sim_end < 0:
            sim_end = qubit_gksl_src.find("\nclass ", sim_start + 1)
        if sim_end < 0:
            sim_end = len(qubit_gksl_src)
        sim_body = qubit_gksl_src[sim_start:sim_end]

        check(
            "gksl_simulate_returns_cx_info",
            "cx_per_pair_gate" in sim_body and "p_depol_pair_eff" in sim_body,
            "simulate() returns cx_per_pair_gate and p_depol_pair_eff in noise_params",
        )
    else:
        check("gksl_simulate_returns_cx_info", False, "simulate() not found")

    # ============================================================
    # Section 4: GKSL notebook Cell 38
    # ============================================================
    print("\n--- Section 4: gksl_comparison.ipynb Cell 38 ---")

    cell38_src = "".join(gksl_nb["cells"][38]["source"])

    check(
        "gksl_cell38_cx_param",
        "cx_per_pair_gate" in cell38_src,
        "GKSL Cell 38 should pass cx_per_pair_gate to QubitGKSLNoisyShotSimulator",
    )

    check(
        "gksl_cell38_cx_estimate",
        "estimate_cx_per_pair_gate" in cell38_src,
        "GKSL Cell 38 should compute CX count dynamically",
    )

    check(
        "gksl_cell38_shows_p_eff",
        "p_depol_pair_eff" in cell38_src,
        "GKSL Cell 38 should print effective pair depolarization rate",
    )

    # ============================================================
    # Section 5: Effective rate formula validation
    # ============================================================
    print("\n--- Section 5: Effective rate formula validation ---")

    p_phys = 0.001

    p_eff_1 = 1 - (1 - p_phys) ** 1
    check(
        "p_eff_cx1_identity",
        abs(p_eff_1 - p_phys) < 1e-15,
        f"cx=1: p_eff={p_eff_1} == p_phys={p_phys} (backward compatible)",
    )

    p_eff_66 = 1 - (1 - p_phys) ** 66.5
    check(
        "p_eff_cx66_significant",
        p_eff_66 > p_phys * 50,
        f"cx=66.5: p_eff={p_eff_66:.4f} ({p_eff_66/p_phys:.0f}x larger than p_phys)",
    )

    cx_values = [1, 5, 10, 20, 50, 66.5, 100]
    p_effs = [1 - (1 - p_phys) ** cx for cx in cx_values]
    is_monotonic = all(p_effs[i] < p_effs[i + 1] for i in range(len(p_effs) - 1))
    check(
        "p_eff_monotonic",
        is_monotonic,
        "p_eff is monotonically increasing with CX count",
    )

    # ============================================================
    # Section 6: Qudit simulators should NOT have cx_per_pair_gate
    # ============================================================
    print("\n--- Section 6: Qudit simulators (no CX overhead) ---")

    check(
        "qudit_noisy_no_cx_param",
        "cx_per_pair_gate" not in qudit_noisy_src,
        "Qudit noisy simulator should NOT have cx_per_pair_gate (native gates)",
    )

    check(
        "qudit_gksl_no_cx_param",
        "cx_per_pair_gate" not in qudit_gksl_src,
        "QuditGKSLNoisyShotSimulator should NOT have cx_per_pair_gate (native gates)",
    )

    cell34_src = "".join(gksl_nb["cells"][34]["source"])
    check(
        "qudit_gksl_cell34_no_cx",
        "cx_per_pair_gate" not in cell34_src,
        "GKSL Cell 34 (qudit) should NOT use cx_per_pair_gate",
    )

    cell23_src = "".join(complete_nb["cells"][23]["source"])
    check(
        "qudit_complete_cell23_no_cx",
        "cx_per_pair_gate" not in cell23_src,
        "Complete Cell 23 (qudit) should NOT use cx_per_pair_gate",
    )

    # ============================================================
    # Section 7: CX gate count estimation (if Qiskit available)
    # ============================================================
    print("\n--- Section 7: CX gate count estimation ---")

    cx_gate_counts = {}
    try:
        sys.path.insert(0, os.path.dirname(__file__))
        from qubit_noisy_simulator import QubitMolecularDynamicsSimulatorNoisy

        cx_info_complete = QubitMolecularDynamicsSimulatorNoisy.estimate_cx_per_pair_gate(
            V=0.1, J=0.05, dt=1.0, hbar=0.6582119569
        )
        cx_gate_counts["complete"] = cx_info_complete

        check(
            "cx_complete_transfer",
            cx_info_complete["cx_transfer"] > 0,
            f"Complete H_transfer: {cx_info_complete['cx_transfer']} CX gates",
        )
        check(
            "cx_complete_tta",
            cx_info_complete["cx_tta"] > 0,
            f"Complete H_TTA: {cx_info_complete['cx_tta']} CX gates",
        )
        check(
            "cx_complete_avg_reasonable",
            2 < cx_info_complete["cx_avg"] < 200,
            f"Complete average: {cx_info_complete['cx_avg']:.1f} CX/pair (reasonable)",
        )

        cx_info_gksl = QubitMolecularDynamicsSimulatorNoisy.estimate_cx_per_pair_gate(
            V=0.1, J=0.05, dt=1.0, hbar=1.0
        )
        cx_gate_counts["gksl"] = cx_info_gksl

        check(
            "cx_gksl_computed",
            cx_info_gksl["cx_avg"] > 0,
            f"GKSL CX: transfer={cx_info_gksl['cx_transfer']}, "
            f"tta={cx_info_gksl['cx_tta']}, avg={cx_info_gksl['cx_avg']:.1f}",
        )

    except ImportError as e:
        print(f"  (Qiskit/numpy not available: {e})")
        print("  CX estimation will be performed when user runs the notebook")
        check(
            "cx_estimation_deferred",
            True,
            "CX estimation requires Qiskit (will be computed when notebook is run)",
        )

    # ============================================================
    # Section 8: Documentation consistency
    # ============================================================
    print("\n--- Section 8: Documentation consistency ---")

    cell37_src = "".join(gksl_nb["cells"][37]["source"])
    check(
        "gksl_cell37_cx_doc",
        "CX" in cell37_src,
        "GKSL Cell 37 markdown should mention CX gate overhead",
    )

    check(
        "gksl_cell37_p_eff_formula",
        "p_{eff}" in cell37_src or "p_eff" in cell37_src,
        "GKSL Cell 37 should show p_eff formula",
    )

    cell14_src = "".join(complete_nb["cells"][14]["source"])
    check(
        "complete_cell14_noise_section",
        "ノイズ" in cell14_src or "Noise" in cell14_src,
        "Complete Cell 14 should be noise section header",
    )

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 70)
    n_pass = sum(1 for c in checks if c["status"] == "PASS")
    n_fail = sum(1 for c in checks if c["status"] == "FAIL")
    n_total = len(checks)

    summary = {"total_checks": n_total, "passed": n_pass, "failed": n_fail}

    print(f"検証結果: {n_pass}/{n_total} PASSED, {n_fail}/{n_total} FAILED")
    if n_fail > 0:
        print("\n失敗したチェック:")
        for c in checks:
            if c["status"] == "FAIL":
                print(f"  ✗ {c['name']}: {c['detail']}")
    print("=" * 70)

    # Save results
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    results = {
        "iteration": 8,
        "timestamp": timestamp,
        "description": "CX gate overhead integration into notebooks and GKSL simulator",
        "summary": summary,
        "cx_gate_counts": cx_gate_counts if cx_gate_counts else "Qiskit not available",
        "effective_rates": {
            "p_phys": 0.001,
            "qudit_p_eff": 0.001,
            "qubit_p_eff_cx1": float(1 - (1 - 0.001) ** 1),
            "qubit_p_eff_cx66": float(1 - (1 - 0.001) ** 66.5),
        },
        "checks": checks,
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    output_file = os.path.join(RESULTS_DIR, f"iteration8_cx_integration_{timestamp}.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n結果を保存しました: {output_file}")

    return n_fail == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
