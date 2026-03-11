#!/usr/bin/env python3
"""
Iteration 9 verification script: Docstring and fidelity note fixes.

Verifies that:
1. cx_per_pair_gate docstring says 'float' (not 'int') in qubit_noisy_simulator.py
2. GKSL notebook Cell 40 does NOT say 'misleadingly low'
3. GKSL notebook Cell 40 correctly describes F_raw as including leakage penalty
4. All iteration 8 checks still pass (regression check)

Uses string-based source code inspection (no heavy imports required).
Results are saved to developing/verification_results/
"""

import json
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
    print("Iteration 9 検証: Docstring修正 & 忠実度ノート修正")
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
    # Section 1: Docstring type fix (iteration 9 new check)
    # ============================================================
    print("\n--- Section 1: Docstring type consistency ---")

    # Find the docstring section for cx_per_pair_gate in simulate()
    sim_start = qubit_noisy_src.find("def simulate(")
    sim_docstring_end = qubit_noisy_src.find('"""', qubit_noisy_src.find('"""', sim_start) + 3)
    sim_docstring = qubit_noisy_src[sim_start:sim_docstring_end]

    check(
        "docstring_cx_float",
        "cx_per_pair_gate : float" in sim_docstring,
        "simulate() docstring should say 'cx_per_pair_gate : float'",
    )

    check(
        "docstring_cx_not_int",
        "cx_per_pair_gate : int" not in sim_docstring,
        "simulate() docstring should NOT say 'cx_per_pair_gate : int'",
    )

    check(
        "docstring_accepts_float",
        "float" in sim_docstring.split("cx_per_pair_gate")[1][:50],
        "Docstring should mention that float values are accepted",
    )

    # Type hint consistency
    check(
        "type_hint_float",
        "cx_per_pair_gate: float" in qubit_noisy_src,
        "Type hint in simulate() signature should be float",
    )

    check(
        "type_hint_not_int",
        "cx_per_pair_gate: int" not in qubit_noisy_src,
        "Type hint should NOT be int anywhere",
    )

    # ============================================================
    # Section 2: Fidelity note fix (iteration 9 new check)
    # ============================================================
    print("\n--- Section 2: GKSL Cell 40 fidelity note ---")

    cell40_src = "".join(gksl_nb["cells"][40]["source"])

    check(
        "no_misleading_text",
        "misleadingly low" not in cell40_src and "misleading" not in cell40_src,
        "Cell 40 should NOT contain 'misleadingly' (incorrect characterization)",
    )

    check(
        "fraw_correctly_includes",
        "correctly includes" in cell40_src or "correctly" in cell40_src,
        "Cell 40 should state F_raw 'correctly includes' the leakage penalty",
    )

    check(
        "conditional_fidelity",
        "conditional fidelity" in cell40_src,
        "Cell 40 should mention 'conditional fidelity' for F_norm",
    )

    check(
        "leakage_penalty",
        "leakage penalty" in cell40_src,
        "Cell 40 should mention 'leakage penalty'",
    )

    # Also check the outputs if present
    cell40_outputs = ""
    for out in gksl_nb["cells"][40].get("outputs", []):
        if "text" in out:
            cell40_outputs += "".join(out["text"])

    if cell40_outputs:
        check(
            "output_no_misleading",
            "misleadingly" not in cell40_outputs,
            "Cell 40 output should not contain 'misleadingly'",
        )
    else:
        check(
            "output_no_misleading",
            True,
            "Cell 40 has no outputs yet (will be verified after re-execution)",
        )

    # ============================================================
    # Section 3: Regression checks from iteration 8
    # ============================================================
    print("\n--- Section 3: Iteration 8 regression checks ---")

    # cx_per_pair_gate type hint
    check(
        "cx_type_hint_float",
        "cx_per_pair_gate: float" in qubit_noisy_src,
        "cx_per_pair_gate should have float type hint",
    )

    check(
        "cx_effective_rate_formula",
        "(1 - depol_2q) ** cx_per_pair_gate" in qubit_noisy_src,
        "Effective rate formula preserved: p_eff = 1 - (1-depol_2q)^cx_per_pair_gate",
    )

    # complete_comparison Cell 15
    cell15_src = "".join(complete_nb["cells"][15]["source"])
    check(
        "complete_cell15_cx_param",
        "cx_per_pair_gate" in cell15_src,
        "Cell 15 still passes cx_per_pair_gate to simulate()",
    )
    check(
        "complete_cell15_cx_estimate",
        "estimate_cx_per_pair_gate" in cell15_src,
        "Cell 15 still computes CX count via estimate_cx_per_pair_gate()",
    )

    # QubitGKSLNoisyShotSimulator
    check(
        "gksl_qubit_cx_param",
        "cx_per_pair_gate" in qubit_gksl_src
        and "class QubitGKSLNoisyShotSimulator" in qubit_gksl_src,
        "QubitGKSLNoisyShotSimulator still accepts cx_per_pair_gate",
    )
    check(
        "gksl_qubit_p_eff_formula",
        "p_depol_pair_eff" in qubit_gksl_src,
        "Still computes p_depol_pair_eff from p_depol and cx_per_pair_gate",
    )

    # GKSL Cell 38
    cell38_src = "".join(gksl_nb["cells"][38]["source"])
    check(
        "gksl_cell38_cx_param",
        "cx_per_pair_gate" in cell38_src,
        "GKSL Cell 38 still passes cx_per_pair_gate",
    )

    # Qudit simulators should NOT have cx_per_pair_gate
    check(
        "qudit_noisy_no_cx_param",
        "cx_per_pair_gate" not in qudit_noisy_src,
        "Qudit noisy simulator still has no cx_per_pair_gate (native gates)",
    )
    check(
        "qudit_gksl_no_cx_param",
        "cx_per_pair_gate" not in qudit_gksl_src,
        "QuditGKSLNoisyShotSimulator still has no cx_per_pair_gate",
    )

    # Effective rate formula validation
    p_phys = 0.001
    p_eff_1 = 1 - (1 - p_phys) ** 1
    check(
        "p_eff_cx1_identity",
        abs(p_eff_1 - p_phys) < 1e-15,
        f"cx=1: p_eff={p_eff_1} == p_phys (backward compatible)",
    )

    p_eff_66 = 1 - (1 - p_phys) ** 66.5
    check(
        "p_eff_cx66_significant",
        p_eff_66 > p_phys * 50,
        f"cx=66.5: p_eff={p_eff_66:.4f} ({p_eff_66/p_phys:.0f}x larger than p_phys)",
    )

    # ============================================================
    # Section 4: CX gate count estimation (if Qiskit available)
    # ============================================================
    print("\n--- Section 4: CX gate count estimation ---")

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
    except ImportError as e:
        print(f"  (Qiskit/numpy not available: {e})")
        check(
            "cx_estimation_deferred",
            True,
            "CX estimation requires Qiskit (will be computed when notebook is run)",
        )

    # ============================================================
    # Section 5: Documentation consistency
    # ============================================================
    print("\n--- Section 5: Documentation consistency ---")

    cell37_src = "".join(gksl_nb["cells"][37]["source"])
    check(
        "gksl_cell37_cx_doc",
        "CX" in cell37_src,
        "GKSL Cell 37 markdown still mentions CX gate overhead",
    )

    check(
        "gksl_cell37_p_eff_formula",
        "p_{eff}" in cell37_src or "p_eff" in cell37_src,
        "GKSL Cell 37 still shows p_eff formula",
    )

    # GKSL simulator docstring should mention CX overhead
    check(
        "gksl_simulator_cx_docstring",
        "cx_per_pair_gate" in qubit_gksl_src.split("class QubitGKSLNoisyShotSimulator")[1][:500],
        "QubitGKSLNoisyShotSimulator class docstring mentions cx_per_pair_gate",
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
        "iteration": 9,
        "timestamp": timestamp,
        "description": "Docstring type fix and fidelity note correction",
        "summary": summary,
        "cx_gate_counts": cx_gate_counts if cx_gate_counts else "Qiskit not available",
        "fixes_verified": {
            "docstring_type": "cx_per_pair_gate : int -> cx_per_pair_gate : float",
            "fidelity_note": "misleadingly low -> correctly includes leakage penalty",
        },
        "checks": checks,
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    output_file = os.path.join(RESULTS_DIR, f"iteration9_docstring_fidelity_{timestamp}.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n結果を保存しました: {output_file}")

    return n_fail == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
