"""Iteration 18 verification: Notebook gate count display completeness.

Verifies that:
1. All simulator cells in the GKSL comparison notebook display estimated_gates_per_step
2. Boson simulator cells display qudit/qubit register counts
3. Gate count values match expected values for each simulator type
4. Cross-consistency between DM and shot simulator gate counts
5. Boson qubit/qudit ratio is consistent

Run:
    cd tutorials
    python run_iteration18_verification.py
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

    # =========================================================
    # Section 1: Notebook cell source code verification
    # =========================================================
    print("\n=== Section 1: GKSL notebook gate count display completeness ===")

    nb_path = os.path.join(os.path.dirname(__file__),
                           "quantum_dynamics_gksl_comparison.ipynb")
    with open(nb_path) as f:
        nb = json.load(f)

    cells = nb["cells"]

    # Cell 6: QuditGKSLSimulator (non-boson)
    src6 = "".join(cells[6]["source"])
    add("cell6_has_gates_per_step",
        'estimated_gates_per_step' in src6,
        "Cell 6 (QuditGKSL) should display estimated_gates_per_step")

    # Cell 8: QubitGKSLSimulator (non-boson)
    src8 = "".join(cells[8]["source"])
    add("cell8_has_gates_per_step",
        'estimated_gates_per_step' in src8,
        "Cell 8 (QubitGKSL) should display estimated_gates_per_step")
    add("cell8_has_total_gates",
        'total_estimated_gates' in src8,
        "Cell 8 (QubitGKSL) should display total_estimated_gates")

    # Cell 12: QuditGKSLBosonSimulator
    src12 = "".join(cells[12]["source"])
    add("cell12_has_gates_per_step",
        'estimated_gates_per_step' in src12,
        "Cell 12 (QuditGKSLBoson) should display estimated_gates_per_step")
    add("cell12_has_system_qudits",
        'n_system_qudits' in src12,
        "Cell 12 (QuditGKSLBoson) should display system qudit count")
    add("cell12_has_phonon_qudits",
        'n_phonon_qudits' in src12,
        "Cell 12 (QuditGKSLBoson) should display phonon qudit count")
    add("cell12_has_ancilla",
        'n_ancilla_qudits' in src12,
        "Cell 12 (QuditGKSLBoson) should display ancilla qudit count")

    # Cell 14: QubitGKSLBosonSimulator
    src14 = "".join(cells[14]["source"])
    add("cell14_has_gates_per_step",
        'estimated_gates_per_step' in src14,
        "Cell 14 (QubitGKSLBoson) should display estimated_gates_per_step")
    add("cell14_has_el_qubits",
        'n_el_qubits' in src14,
        "Cell 14 (QubitGKSLBoson) should display electronic qubit count")
    add("cell14_has_ph_qubits",
        'n_ph_qubits' in src14,
        "Cell 14 (QubitGKSLBoson) should display phonon qubit count")
    add("cell14_has_total_qubits",
        'n_total_qubits' in src14,
        "Cell 14 (QubitGKSLBoson) should display total qubit count")

    # Cell 32: QuditGKSLShotSimulator
    src32 = "".join(cells[32]["source"])
    add("cell32_has_gates_per_step",
        'estimated_gates_per_step' in src32,
        "Cell 32 (QuditGKSLShot) should display estimated_gates_per_step")

    # Cell 36: QubitGKSLShotSimulator
    src36 = "".join(cells[36]["source"])
    add("cell36_has_gates_per_step",
        'estimated_gates_per_step' in src36,
        "Cell 36 (QubitGKSLShot) should display estimated_gates_per_step")

    # =========================================================
    # Section 2: Gate count value verification (runtime)
    # =========================================================
    print("\n=== Section 2: Gate count value verification ===")

    from gksl_physical_parameters import GKSLPhysicalParameters

    # Non-boson N=4
    params = GKSLPhysicalParameters()

    from qudit_gksl_simulator import QuditGKSLSimulator
    sim5 = QuditGKSLSimulator(params)
    r5 = sim5.simulate(t_max=100.0, n_steps=10, initial_state="edge_triplet")
    add("qudit_dm_gps_66",
        r5["estimated_gates_per_step"] == 66,
        f"QuditGKSL N=4: expected 66, got {r5['estimated_gates_per_step']}")

    from qubit_gksl_simulator import QubitGKSLSimulator
    sim3 = QubitGKSLSimulator(params)
    r3 = sim3.simulate(t_max=100.0, n_steps=10, initial_state="edge_triplet")
    add("qubit_dm_gps_388",
        r3["estimated_gates_per_step"] == 388,
        f"QubitGKSL N=4: expected 388, got {r3['estimated_gates_per_step']}")

    # Boson N=2
    params_boson = GKSLPhysicalParameters(N_molecules=2, with_boson=True, n_max=2)

    from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
    sim6 = QuditGKSLBosonSimulator(params_boson)
    r6 = sim6.simulate(t_max=100.0, n_steps=10, initial_state="edge_triplet")
    add("qudit_boson_gps_38",
        r6["estimated_gates_per_step"] == 38,
        f"QuditGKSLBoson N=2: expected 38, got {r6['estimated_gates_per_step']}")

    from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator
    sim4 = QubitGKSLBosonSimulator(params_boson)
    r4 = sim4.simulate(t_max=100.0, n_steps=10, initial_state="edge_triplet")
    add("qubit_boson_gps_220",
        r4["estimated_gates_per_step"] == 220,
        f"QubitGKSLBoson N=2: expected 220, got {r4['estimated_gates_per_step']}")

    # Shot simulators
    from qudit_gksl_shot_simulator import QuditGKSLShotSimulator
    sim5b = QuditGKSLShotSimulator(params)
    r5b = sim5b.simulate(t_max=100.0, n_steps=10, initial_state="edge_triplet",
                          n_shots=100, seed=42)
    add("qudit_shot_gps_66",
        r5b["estimated_gates_per_step"] == 66,
        f"QuditGKSLShot N=4: expected 66, got {r5b['estimated_gates_per_step']}")

    from qubit_gksl_shot_simulator import QubitGKSLShotSimulator
    sim3b = QubitGKSLShotSimulator(params)
    r3b = sim3b.simulate(t_max=100.0, n_steps=10, initial_state="edge_triplet",
                          n_shots=100, seed=42)
    add("qubit_shot_gps_388",
        r3b["estimated_gates_per_step"] == 388,
        f"QubitGKSLShot N=4: expected 388, got {r3b['estimated_gates_per_step']}")

    # =========================================================
    # Section 3: Cross-consistency checks
    # =========================================================
    print("\n=== Section 3: Cross-consistency checks ===")

    # DM and shot should match for same simulator type
    add("qudit_dm_shot_match",
        r5["estimated_gates_per_step"] == r5b["estimated_gates_per_step"],
        f"QuditGKSL DM ({r5['estimated_gates_per_step']}) == "
        f"Shot ({r5b['estimated_gates_per_step']})")

    add("qubit_dm_shot_match",
        r3["estimated_gates_per_step"] == r3b["estimated_gates_per_step"],
        f"QubitGKSL DM ({r3['estimated_gates_per_step']}) == "
        f"Shot ({r3b['estimated_gates_per_step']})")

    # Boson qubit/qudit ratio
    ratio = r4["estimated_gates_per_step"] / r6["estimated_gates_per_step"]
    add("boson_ratio_n2",
        abs(ratio - 220 / 38) < 0.01,
        f"Boson N=2: qubit/qudit = {ratio:.2f}× ({r4['estimated_gates_per_step']}/{r6['estimated_gates_per_step']})")

    # Non-boson qubit/qudit ratio
    ratio_nb = r3["estimated_gates_per_step"] / r5["estimated_gates_per_step"]
    add("nonboson_ratio_n4",
        abs(ratio_nb - 388 / 66) < 0.01,
        f"Non-boson N=4: qubit/qudit = {ratio_nb:.2f}× ({r3['estimated_gates_per_step']}/{r5['estimated_gates_per_step']})")

    # =========================================================
    # Section 4: Boson register count verification
    # =========================================================
    print("\n=== Section 4: Boson register count verification ===")

    # QuditGKSLBoson register counts
    add("qudit_boson_sys_qudits",
        r6["n_system_qudits"] == 2,
        f"QuditGKSLBoson N=2: n_system_qudits = {r6['n_system_qudits']} (expected 2)")
    add("qudit_boson_ph_qudits",
        r6["n_phonon_qudits"] == 2,
        f"QuditGKSLBoson N=2: n_phonon_qudits = {r6['n_phonon_qudits']} (expected 2)")
    add("qudit_boson_ancilla",
        r6["n_ancilla_qudits"] == 12,
        f"QuditGKSLBoson N=2: n_ancilla_qudits = {r6['n_ancilla_qudits']} (expected 12)")

    # QubitGKSLBoson register counts
    add("qubit_boson_el_qubits",
        r4["n_el_qubits"] == 4,
        f"QubitGKSLBoson N=2: n_el_qubits = {r4['n_el_qubits']} (expected 4)")
    add("qubit_boson_ph_qubits",
        r4["n_ph_qubits"] == 4,
        f"QubitGKSLBoson N=2: n_ph_qubits = {r4['n_ph_qubits']} (expected 4)")
    add("qubit_boson_n_total",
        r4["n_total_qubits"] == 4 + 4 + 12,
        f"QubitGKSLBoson N=2: n_total_qubits = {r4['n_total_qubits']} (expected 20)")

    # =========================================================
    # Section 5: Regression — iteration 17 boson gate counts (N=4)
    # Compute from simulator attributes (no full simulation needed —
    # N=4 boson Hilbert space dim=6561 is too large for CI memory).
    # =========================================================
    print("\n=== Section 5: Regression — iteration 17 boson gate counts (N=4) ===")

    params_boson4 = GKSLPhysicalParameters(N_molecules=4, with_boson=True, n_max=2)

    # Qudit boson N=4: compute gates_per_step from attributes
    sim6_n4 = QuditGKSLBosonSimulator(params_boson4)
    gps_qudit_n4 = (
        2 * (sim6_n4.n_system_qudits + sim6_n4.n_phonon_qudits
             + len(sim6_n4.params.neighbors) + sim6_n4.params.N_molecules)
        + sim6_n4.n_ancilla_qudits * 2
    )
    add("regression_qudit_boson_n4_82",
        gps_qudit_n4 == 82,
        f"QuditGKSLBoson N=4: expected 82, computed {gps_qudit_n4}")

    # Qubit boson N=4: compute gates_per_step from attributes
    sim4_n4 = QubitGKSLBosonSimulator(params_boson4)
    gps_qubit_n4 = (
        2 * (sim4_n4.n_el_qubits + sim4_n4.n_ph_qubits
             + len(sim4_n4.params.neighbors) * 10 + sim4_n4.params.N_molecules * 10)
        + sim4_n4.n_ancilla * 6 * 2
    )
    add("regression_qubit_boson_n4_484",
        gps_qubit_n4 == 484,
        f"QubitGKSLBoson N=4: expected 484, computed {gps_qubit_n4}")

    # =========================================================
    # Summary
    # =========================================================
    n_pass = sum(1 for c in checks if c["passed"])
    n_fail = sum(1 for c in checks if not c["passed"])
    n_total = len(checks)

    print(f"\n{'=' * 60}")
    print(f"Total: {n_total} checks, {n_pass} passed, {n_fail} failed")
    if n_fail > 0:
        print("FAILED checks:")
        for c in checks:
            if not c["passed"]:
                print(f"  - {c['name']}: {c['detail']}")
    print(f"{'=' * 60}")

    # Save results
    results = {
        "iteration": 18,
        "timestamp": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "total_checks": n_total,
        "passed": n_pass,
        "failed": n_fail,
        "all_passed": n_fail == 0,
        "checks": checks,
    }

    out_dir = os.path.join(os.path.dirname(__file__),
                           "..", "developing", "verification_results")
    os.makedirs(out_dir, exist_ok=True)
    ts = results["timestamp"]
    out_path = os.path.join(out_dir,
                            f"iteration18_notebook_gate_display_{ts}.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
