#!/usr/bin/env python3
"""Iteration 21 verification script.

Checks:
- Issue Y fix: Cell 36 (Qubit Shot) now displays measurement counts + forbidden_count
- Issue Z fix: Section 5 result dict completeness expanded to include qubit DM/boson
- Regression: All iteration 20 checks still pass (34 checks)

Usage:
    cd tutorials
    python run_iteration21_verification.py
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
    # Section 1: n_total correctness in simulators (regression)
    # ================================================================
    print("\n=== Section 1: n_total in simulator classes ===")

    from gksl_physical_parameters import GKSLPhysicalParameters

    # Non-boson N=4
    params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)

    from qudit_gksl_simulator import QuditGKSLSimulator

    sim_qudit = QuditGKSLSimulator(params)
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
    add(
        "qudit_boson_n_total_value",
        sim_qudit_boson.n_total_qudits == 16,
        f"QuditGKSLBoson N=2: n_total_qudits={sim_qudit_boson.n_total_qudits}, expected 16 (2+2+12)",
    )

    from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator

    sim_qubit_boson = QubitGKSLBosonSimulator(params_boson)
    add(
        "qubit_boson_n_total_value",
        sim_qubit_boson.n_total_qubits == 20,
        f"QubitGKSLBoson N=2: n_total_qubits={sim_qubit_boson.n_total_qubits}, expected 20 (4+4+12)",
    )

    # Shot simulators
    from qudit_gksl_shot_simulator import QuditGKSLShotSimulator

    sim_qudit_shot = QuditGKSLShotSimulator(params)
    add(
        "qudit_shot_n_total_value",
        sim_qudit_shot.n_total_qudits == 30,
        f"QuditGKSLShot N=4: n_total_qudits={sim_qudit_shot.n_total_qudits}, expected 30 (4+26)",
    )

    from qubit_gksl_shot_simulator import QubitGKSLShotSimulator

    sim_qubit_shot = QubitGKSLShotSimulator(params)
    add(
        "qubit_shot_n_total_value",
        sim_qubit_shot.n_total_qubits == 34,
        f"QubitGKSLShot N=4: n_total_qubits={sim_qubit_shot.n_total_qubits}, expected 34 (8+26)",
    )

    # Boson N=4
    params_boson_n4 = GKSLPhysicalParameters(N_molecules=4, with_boson=True)
    sim_qudit_boson_n4 = QuditGKSLBosonSimulator(params_boson_n4)
    add(
        "qudit_boson_n4_n_total_value",
        sim_qudit_boson_n4.n_total_qudits == 34,
        f"QuditGKSLBoson N=4: n_total_qudits={sim_qudit_boson_n4.n_total_qudits}, expected 34 (4+4+26)",
    )

    # N=4 boson qubit: constructor requires ~7 GB (20736-dim Hilbert space),
    # so compute expected values from the formulas directly.
    import numpy as np

    n_el_qubits_n4 = 2 * params_boson_n4.N_molecules
    n_ph_qubits_n4 = int(np.ceil(np.log2(params_boson_n4.n_max + 1))) * params_boson_n4.N_molecules
    n_ancilla_qubit_boson_n4 = 2 * len(params_boson_n4.neighbors) + 5 * params_boson_n4.N_molecules
    n_total_qubit_boson_n4 = n_el_qubits_n4 + n_ph_qubits_n4 + n_ancilla_qubit_boson_n4
    add(
        "qubit_boson_n4_n_total_value",
        n_total_qubit_boson_n4 == 42,
        f"QubitGKSLBoson N=4 (formula): n_total_qubits={n_total_qubit_boson_n4}, expected 42 (8+8+26)",
    )

    # ================================================================
    # Section 2: Cross-consistency checks
    # ================================================================
    print("\n=== Section 2: Cross-consistency ===")

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
    add(
        "qudit_dm_shot_total_match",
        sim_qudit.n_total_qudits == sim_qudit_shot.n_total_qudits,
        f"Qudit DM ({sim_qudit.n_total_qudits}) == Shot ({sim_qudit_shot.n_total_qudits})",
    )
    add(
        "qubit_dm_shot_total_match",
        sim_qubit.n_total_qubits == sim_qubit_shot.n_total_qubits,
        f"Qubit DM ({sim_qubit.n_total_qubits}) == Shot ({sim_qubit_shot.n_total_qubits})",
    )

    # ================================================================
    # Section 3: GPS regression (run simulations)
    # ================================================================
    print("\n=== Section 3: Gate count regression ===")

    result_qudit = sim_qudit.simulate(t_max=0.1, n_steps=2)
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

    result_qudit_boson = sim_qudit_boson.simulate(t_max=0.1, n_steps=2)
    add(
        "regression_qudit_boson_n2_gps_38",
        result_qudit_boson["estimated_gates_per_step"] == 38,
        f"QuditGKSLBoson N=2: GPS={result_qudit_boson['estimated_gates_per_step']}, expected 38",
    )

    result_qubit_boson = sim_qubit_boson.simulate(t_max=0.1, n_steps=2)
    add(
        "regression_qubit_boson_n2_gps_220",
        result_qubit_boson["estimated_gates_per_step"] == 220,
        f"QubitGKSLBoson N=2: GPS={result_qubit_boson['estimated_gates_per_step']}, expected 220",
    )

    result_qudit_shot = sim_qudit_shot.simulate(t_max=0.1, n_steps=2, n_shots=10, seed=42)
    add(
        "regression_qudit_shot_gps_66",
        result_qudit_shot["estimated_gates_per_step"] == 66,
        f"QuditGKSLShot N=4: GPS={result_qudit_shot['estimated_gates_per_step']}, expected 66",
    )

    result_qubit_shot = sim_qubit_shot.simulate(
        t_max=0.1, n_steps=2, n_shots=10, seed=42
    )
    add(
        "regression_qubit_shot_gps_388",
        result_qubit_shot["estimated_gates_per_step"] == 388,
        f"QubitGKSLShot N=4: GPS={result_qubit_shot['estimated_gates_per_step']}, expected 388",
    )

    # Boson N=4 GPS regression (formula-based to avoid 7GB memory)
    gps_qudit_boson_n4 = (
        2 * (sim_qudit_boson_n4.n_system_qudits + sim_qudit_boson_n4.n_phonon_qudits
             + len(params_boson_n4.neighbors) + params_boson_n4.N_molecules)
        + sim_qudit_boson_n4.n_ancilla_qudits * 2
    )
    add(
        "regression_qudit_boson_n4_gps_82",
        gps_qudit_boson_n4 == 82,
        f"QuditGKSLBoson N=4: GPS={gps_qudit_boson_n4}, expected 82",
    )

    gps_qubit_boson_n4 = (
        2 * (n_el_qubits_n4 + n_ph_qubits_n4
             + len(params_boson_n4.neighbors) * 10 + params_boson_n4.N_molecules * 10)
        + n_ancilla_qubit_boson_n4 * 6 * 2
    )
    add(
        "regression_qubit_boson_n4_gps_484",
        gps_qubit_boson_n4 == 484,
        f"QubitGKSLBoson N=4: GPS={gps_qubit_boson_n4}, expected 484",
    )

    # ================================================================
    # Section 4: Notebook cell inspection
    # ================================================================
    print("\n=== Section 4: Notebook cell inspection ===")

    nb_path = os.path.join(
        os.path.dirname(__file__), "quantum_dynamics_gksl_comparison.ipynb"
    )
    with open(nb_path) as f:
        nb = json.load(f)
    cells = nb["cells"]

    # Cells 6, 8, 12, 14 (DM and boson - from iteration 19)
    cell6_src = "".join(cells[6]["source"])
    add(
        "cell6_has_total_qudits",
        "n_total_qudits" in cell6_src,
        "Cell 6 (QuditGKSL) should display n_total_qudits",
    )

    cell8_src = "".join(cells[8]["source"])
    add(
        "cell8_has_total_qubits",
        "n_total_qubits" in cell8_src,
        "Cell 8 (QubitGKSL) should display n_total_qubits",
    )

    cell12_src = "".join(cells[12]["source"])
    add(
        "cell12_has_total_qudits",
        "n_total_qudits" in cell12_src,
        "Cell 12 (QuditGKSLBoson) should display n_total_qudits",
    )

    cell14_src = "".join(cells[14]["source"])
    add(
        "cell14_has_total_qubits",
        "n_total_qubits" in cell14_src,
        "Cell 14 (QubitGKSLBoson) should display n_total_qubits",
    )

    # Cells 32, 36 (Shot simulators)
    cell32_src = "".join(cells[32]["source"])
    add(
        "cell32_has_total_qudits",
        "n_total_qudits" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display n_total_qudits",
    )
    add(
        "cell32_has_system_qudits",
        "n_system_qudits" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display n_system_qudits",
    )
    add(
        "cell32_has_ancilla_qudits",
        "n_ancilla_qudits" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display n_ancilla_qudits",
    )
    add(
        "cell32_has_counts",
        "counts" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display measurement counts",
    )

    cell36_src = "".join(cells[36]["source"])
    add(
        "cell36_has_total_qubits",
        "n_total_qubits" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display n_total_qubits",
    )
    add(
        "cell36_has_ancilla",
        "n_ancilla" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display n_ancilla",
    )
    add(
        "cell36_has_trace",
        "Trace conservation" in cell36_src or "trace" in cell36_src.lower(),
        "Cell 36 (QubitGKSLShot) should display trace conservation",
    )
    # Issue Y fix: Cell 36 should now display measurement counts and forbidden_count
    add(
        "cell36_has_counts",
        "counts" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display measurement counts",
    )
    add(
        "cell36_has_forbidden_count",
        "forbidden_count" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display forbidden_count",
    )

    # ================================================================
    # Section 5: Result dict completeness (expanded with qubit checks)
    # ================================================================
    print("\n=== Section 5: Result dict completeness ===")

    add(
        "qudit_dm_result_n_total",
        result_qudit.get("n_total_qudits") == 30,
        f"QuditGKSL result n_total_qudits={result_qudit.get('n_total_qudits')}, expected 30",
    )
    add(
        "qubit_dm_result_n_total",
        result_qubit.get("n_total_qubits") == 34,
        f"QubitGKSL result n_total_qubits={result_qubit.get('n_total_qubits')}, expected 34",
    )
    add(
        "qudit_boson_result_n_total",
        result_qudit_boson.get("n_total_qudits") == 16,
        f"QuditGKSLBoson result n_total_qudits={result_qudit_boson.get('n_total_qudits')}, expected 16",
    )
    add(
        "qubit_boson_result_n_total",
        result_qubit_boson.get("n_total_qubits") == 20,
        f"QubitGKSLBoson result n_total_qubits={result_qubit_boson.get('n_total_qubits')}, expected 20",
    )
    add(
        "qudit_shot_result_n_total",
        result_qudit_shot.get("n_total_qudits") == 30,
        f"QuditGKSLShot result n_total_qudits={result_qudit_shot.get('n_total_qudits')}, expected 30",
    )
    add(
        "qubit_shot_result_n_total",
        result_qubit_shot.get("n_total_qubits") == 34,
        f"QubitGKSLShot result n_total_qubits={result_qubit_shot.get('n_total_qubits')}, expected 34",
    )

    # ================================================================
    # Section 6: Qubit shot forbidden_count (noiseless should be 0)
    # ================================================================
    print("\n=== Section 6: Qubit shot forbidden_count ===")

    add(
        "qubit_shot_forbidden_count_zero",
        result_qubit_shot.get("forbidden_count") == 0,
        f"QubitGKSLShot noiseless: forbidden_count={result_qubit_shot.get('forbidden_count')}, expected 0",
    )
    add(
        "qubit_shot_has_counts_key",
        "counts" in result_qubit_shot,
        "QubitGKSLShot result dict should contain 'counts' key",
    )
    add(
        "qubit_shot_has_forbidden_key",
        "forbidden_count" in result_qubit_shot,
        "QubitGKSLShot result dict should contain 'forbidden_count' key",
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
        "iteration": 21,
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
        out_dir, f"iteration21_cell36_display_{results['timestamp']}.json"
    )
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to: {out_file}")


if __name__ == "__main__":
    main()
