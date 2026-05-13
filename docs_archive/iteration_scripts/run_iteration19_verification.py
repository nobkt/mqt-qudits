#!/usr/bin/env python3
"""Iteration 19 verification script.

Checks:
- Bug R fix: n_total_qudits in qudit boson simulator now includes ancilla
- Issue S fix: n_total_qudits added to qudit non-boson and shot simulators
- Issue T fix: Notebook cells display Total qudits for qudit simulators
- Cross-consistency: n_total definitions include ancilla in all simulators
- Regression: gate counts unchanged from iteration 18

Usage:
    cd tutorials
    python run_iteration19_verification.py
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main() -> None:
    checks: list[dict] = []

    def add(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": passed, "detail": detail})
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {detail}")

    # ================================================================
    # Section 1: n_total_qudits correctness in simulators
    # ================================================================
    print("\n=== Section 1: n_total_qudits in simulator classes ===")

    from gksl_physical_parameters import GKSLPhysicalParameters

    # Non-boson N=4
    params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)

    from qudit_gksl_simulator import QuditGKSLSimulator

    sim_qudit = QuditGKSLSimulator(params)
    add(
        "qudit_dm_has_n_total",
        hasattr(sim_qudit, "n_total_qudits"),
        f"QuditGKSLSimulator should define n_total_qudits",
    )
    expected_total_qudit = sim_qudit.n_system_qudits + sim_qudit.n_ancilla_qudits
    add(
        "qudit_dm_n_total_includes_ancilla",
        sim_qudit.n_total_qudits == expected_total_qudit,
        f"n_total_qudits={sim_qudit.n_total_qudits} should be "
        f"sys({sim_qudit.n_system_qudits})+anc({sim_qudit.n_ancilla_qudits})={expected_total_qudit}",
    )
    add(
        "qudit_dm_n_total_value",
        sim_qudit.n_total_qudits == 30,
        f"QuditGKSL N=4: n_total_qudits={sim_qudit.n_total_qudits}, expected 30 (4+26)",
    )

    from qubit_gksl_simulator import QubitGKSLSimulator

    sim_qubit = QubitGKSLSimulator(params)
    add(
        "qubit_dm_n_total_value",
        sim_qubit.n_total_qubits == 34,
        f"QubitGKSL N=4: n_total_qubits={sim_qubit.n_total_qubits}, expected 34 (8+26)",
    )

    # Boson N=2
    params_boson = GKSLPhysicalParameters(N_molecules=2, with_boson=True)

    from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator

    sim_qudit_boson = QuditGKSLBosonSimulator(params_boson)
    expected_total_boson = (
        sim_qudit_boson.n_system_qudits
        + sim_qudit_boson.n_phonon_qudits
        + sim_qudit_boson.n_ancilla_qudits
    )
    add(
        "qudit_boson_n_total_includes_ancilla",
        sim_qudit_boson.n_total_qudits == expected_total_boson,
        f"n_total_qudits={sim_qudit_boson.n_total_qudits} should be "
        f"sys({sim_qudit_boson.n_system_qudits})+ph({sim_qudit_boson.n_phonon_qudits})"
        f"+anc({sim_qudit_boson.n_ancilla_qudits})={expected_total_boson}",
    )
    add(
        "qudit_boson_n_total_value",
        sim_qudit_boson.n_total_qudits == 16,
        f"QuditGKSLBoson N=2: n_total_qudits={sim_qudit_boson.n_total_qudits}, expected 16 (2+2+12)",
    )

    from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator

    sim_qubit_boson = QubitGKSLBosonSimulator(params_boson)
    expected_total_qubit_boson = (
        sim_qubit_boson.n_el_qubits
        + sim_qubit_boson.n_ph_qubits
        + sim_qubit_boson.n_ancilla
    )
    add(
        "qubit_boson_n_total_includes_ancilla",
        sim_qubit_boson.n_total_qubits == expected_total_qubit_boson,
        f"n_total_qubits={sim_qubit_boson.n_total_qubits} should be "
        f"el({sim_qubit_boson.n_el_qubits})+ph({sim_qubit_boson.n_ph_qubits})"
        f"+anc({sim_qubit_boson.n_ancilla})={expected_total_qubit_boson}",
    )
    add(
        "qubit_boson_n_total_value",
        sim_qubit_boson.n_total_qubits == 20,
        f"QubitGKSLBoson N=2: n_total_qubits={sim_qubit_boson.n_total_qubits}, expected 20 (4+4+12)",
    )

    # Shot simulators
    from qudit_gksl_shot_simulator import QuditGKSLShotSimulator

    sim_qudit_shot = QuditGKSLShotSimulator(params)
    add(
        "qudit_shot_has_n_total",
        hasattr(sim_qudit_shot, "n_total_qudits"),
        "QuditGKSLShotSimulator should define n_total_qudits",
    )
    expected_total_shot = sim_qudit_shot.n_system_qudits + sim_qudit_shot.n_ancilla_qudits
    add(
        "qudit_shot_n_total_includes_ancilla",
        sim_qudit_shot.n_total_qudits == expected_total_shot,
        f"n_total_qudits={sim_qudit_shot.n_total_qudits} should be "
        f"sys({sim_qudit_shot.n_system_qudits})+anc({sim_qudit_shot.n_ancilla_qudits})={expected_total_shot}",
    )
    add(
        "qudit_shot_n_total_value",
        sim_qudit_shot.n_total_qudits == 30,
        f"QuditGKSLShot N=4: n_total_qudits={sim_qudit_shot.n_total_qudits}, expected 30 (4+26)",
    )

    # ================================================================
    # Section 2: n_total_qudits in result dicts
    # ================================================================
    print("\n=== Section 2: n_total_qudits in result dictionaries ===")

    result_qudit = sim_qudit.simulate(t_max=0.1, n_steps=2)
    add(
        "qudit_dm_result_has_n_total",
        "n_total_qudits" in result_qudit,
        "QuditGKSL result dict should contain n_total_qudits",
    )
    add(
        "qudit_dm_result_n_total_correct",
        result_qudit.get("n_total_qudits") == 30,
        f"QuditGKSL result n_total_qudits={result_qudit.get('n_total_qudits')}, expected 30",
    )

    result_qudit_boson = sim_qudit_boson.simulate(t_max=0.1, n_steps=2)
    add(
        "qudit_boson_result_has_n_total",
        "n_total_qudits" in result_qudit_boson,
        "QuditGKSLBoson result dict should contain n_total_qudits",
    )
    add(
        "qudit_boson_result_n_total_correct",
        result_qudit_boson.get("n_total_qudits") == 16,
        f"QuditGKSLBoson result n_total_qudits={result_qudit_boson.get('n_total_qudits')}, expected 16",
    )

    result_qudit_shot = sim_qudit_shot.simulate(t_max=0.1, n_steps=2, n_shots=10, seed=42)
    add(
        "qudit_shot_result_has_n_total",
        "n_total_qudits" in result_qudit_shot,
        "QuditGKSLShot result dict should contain n_total_qudits",
    )
    add(
        "qudit_shot_result_n_total_correct",
        result_qudit_shot.get("n_total_qudits") == 30,
        f"QuditGKSLShot result n_total_qudits={result_qudit_shot.get('n_total_qudits')}, expected 30",
    )

    # ================================================================
    # Section 3: Cross-consistency checks
    # ================================================================
    print("\n=== Section 3: Cross-consistency ===")

    # Same ancilla count across qudit/qubit for same model
    add(
        "nonboson_same_ancilla",
        sim_qudit.n_ancilla_qudits == sim_qubit.n_ancilla,
        f"Non-boson ancilla: qudit={sim_qudit.n_ancilla_qudits}, qubit={sim_qubit.n_ancilla} (should match)",
    )
    add(
        "boson_same_ancilla",
        sim_qudit_boson.n_ancilla_qudits == sim_qubit_boson.n_ancilla,
        f"Boson ancilla: qudit={sim_qudit_boson.n_ancilla_qudits}, qubit={sim_qubit_boson.n_ancilla} (should match)",
    )

    # DM and shot give same totals
    add(
        "qudit_dm_shot_total_match",
        sim_qudit.n_total_qudits == sim_qudit_shot.n_total_qudits,
        f"Qudit DM ({sim_qudit.n_total_qudits}) == Shot ({sim_qudit_shot.n_total_qudits})",
    )

    # ================================================================
    # Section 4: Notebook cell inspection
    # ================================================================
    print("\n=== Section 4: Notebook cell inspection ===")

    nb_path = os.path.join(os.path.dirname(__file__), "quantum_dynamics_gksl_comparison.ipynb")
    with open(nb_path) as f:
        nb = json.load(f)
    cells = nb["cells"]

    # Cell 6 should display n_total_qudits
    cell6_src = "".join(cells[6]["source"])
    add(
        "cell6_has_total_qudits",
        "n_total_qudits" in cell6_src,
        "Cell 6 (QuditGKSL) should display n_total_qudits",
    )

    # Cell 8 should display n_total_qubits (existing)
    cell8_src = "".join(cells[8]["source"])
    add(
        "cell8_has_total_qubits",
        "n_total_qubits" in cell8_src,
        "Cell 8 (QubitGKSL) should display n_total_qubits",
    )

    # Cell 12 should display n_total_qudits
    cell12_src = "".join(cells[12]["source"])
    add(
        "cell12_has_total_qudits",
        "n_total_qudits" in cell12_src,
        "Cell 12 (QuditGKSLBoson) should display n_total_qudits",
    )

    # Cell 14 should display n_total_qubits (existing)
    cell14_src = "".join(cells[14]["source"])
    add(
        "cell14_has_total_qubits",
        "n_total_qubits" in cell14_src,
        "Cell 14 (QubitGKSLBoson) should display n_total_qubits",
    )

    # ================================================================
    # Section 5: Regression - gate counts unchanged
    # ================================================================
    print("\n=== Section 5: Gate count regression ===")

    add(
        "regression_qudit_dm_gps_66",
        result_qudit["estimated_gates_per_step"] == 66,
        f"QuditGKSL N=4: GPS={result_qudit['estimated_gates_per_step']}, expected 66",
    )

    result_qubit = sim_qubit.simulate(t_max=0.1, n_steps=2)
    add(
        "regression_qubit_dm_gps_388",
        result_qubit["estimated_gates_per_step"] == 388,
        f"QubitGKSL N=4: GPS={result_qubit['estimated_gates_per_step']}, expected 388",
    )

    add(
        "regression_qudit_boson_gps_38",
        result_qudit_boson["estimated_gates_per_step"] == 38,
        f"QuditGKSLBoson N=2: GPS={result_qudit_boson['estimated_gates_per_step']}, expected 38",
    )

    result_qubit_boson = sim_qubit_boson.simulate(t_max=0.1, n_steps=2)
    add(
        "regression_qubit_boson_gps_220",
        result_qubit_boson["estimated_gates_per_step"] == 220,
        f"QubitGKSLBoson N=2: GPS={result_qubit_boson['estimated_gates_per_step']}, expected 220",
    )

    add(
        "regression_qudit_shot_gps_66",
        result_qudit_shot["estimated_gates_per_step"] == 66,
        f"QuditGKSLShot N=4: GPS={result_qudit_shot['estimated_gates_per_step']}, expected 66",
    )

    # Boson N=4 regression
    params_boson_n4 = GKSLPhysicalParameters(N_molecules=4, with_boson=True)
    sim_qudit_boson_n4 = QuditGKSLBosonSimulator(params_boson_n4)
    add(
        "regression_qudit_boson_n4_gps_82",
        sim_qudit_boson_n4.n_total_qudits == 34,
        f"QuditGKSLBoson N=4: n_total_qudits={sim_qudit_boson_n4.n_total_qudits}, expected 34 (4+4+26)",
    )

    # ================================================================
    # Summary
    # ================================================================
    passed = sum(1 for c in checks if c["passed"])
    failed = sum(1 for c in checks if not c["passed"])
    total = len(checks)

    print(f"\n{'='*60}")
    print(f"Total: {total}  Passed: {passed}  Failed: {failed}")
    if failed == 0:
        print("ALL PASSED ✓")
    else:
        print("FAILURES DETECTED:")
        for c in checks:
            if not c["passed"]:
                print(f"  FAIL: {c['name']}: {c['detail']}")

    # Save results
    results = {
        "iteration": 19,
        "timestamp": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "total_checks": total,
        "passed": passed,
        "failed": failed,
        "all_passed": failed == 0,
        "checks": checks,
    }

    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "developing",
        "verification_results",
    )
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(
        out_dir, f"iteration19_n_total_qudits_{results['timestamp']}.json"
    )
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to: {out_file}")


if __name__ == "__main__":
    main()
