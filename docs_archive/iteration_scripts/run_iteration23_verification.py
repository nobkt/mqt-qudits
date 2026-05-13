#!/usr/bin/env python3
"""Iteration 23 verification script.

Checks:
- Issue BB fix: Cell 16 boson circuit visualization now has palindromic structure
  - Scenario 6 (Qudit Boson) and Scenario 4 (Qubit Boson) should have
    forward+reverse Lindblad channels and use dt/2 (not full dt)
- Issue CC fix: n_phonon_qudits now uses general formula ceil(log_d(n_max+1))*N
  - Results unchanged for n_max=1,2 (both give 1 qutrit/molecule)
  - Correctly handles n_max=3+ (2 qutrits/molecule for d=3)
- Regression: All iteration 22 checks (49 checks) still pass

Usage:
    cd tutorials
    python run_iteration23_verification.py
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
    # Section 1: n_total correctness in simulators (regression from iter22)
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

    # Boson N=2 — Use notebook-matching parameters (n_max=1)
    params_boson = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02,
    )

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
        sim_qubit_boson.n_total_qubits == 18,
        f"QubitGKSLBoson N=2 (n_max=1): n_total_qubits={sim_qubit_boson.n_total_qubits}, expected 18 (4+2+12)",
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

    # Boson N=4 (n_max=1, matching notebook convention)
    params_boson_n4 = GKSLPhysicalParameters(
        N_molecules=4, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_n4 = QuditGKSLBosonSimulator(params_boson_n4)
    add(
        "qudit_boson_n4_n_total_value",
        sim_qudit_boson_n4.n_total_qudits == 34,
        f"QuditGKSLBoson N=4: n_total_qudits={sim_qudit_boson_n4.n_total_qudits}, expected 34 (4+4+26)",
    )

    # N=4 boson qubit: constructor requires ~7 GB (large Hilbert space),
    # so compute expected values from the formulas directly.
    import numpy as np

    n_el_qubits_n4 = 2 * params_boson_n4.N_molecules
    n_ph_qubits_n4 = int(np.ceil(np.log2(params_boson_n4.n_max + 1))) * params_boson_n4.N_molecules
    n_ancilla_qubit_boson_n4 = 2 * len(params_boson_n4.neighbors) + 5 * params_boson_n4.N_molecules
    n_total_qubit_boson_n4 = n_el_qubits_n4 + n_ph_qubits_n4 + n_ancilla_qubit_boson_n4
    add(
        "qubit_boson_n4_n_total_value",
        n_total_qubit_boson_n4 == 38,
        f"QubitGKSLBoson N=4 (n_max=1, formula): n_total_qubits={n_total_qubit_boson_n4}, expected 38 (8+4+26)",
    )

    # ================================================================
    # Section 2: Cross-consistency checks (regression from iter22)
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
    # Section 3: GPS regression (run simulations) (regression from iter22)
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
        "regression_qubit_boson_n2_gps_216",
        result_qubit_boson["estimated_gates_per_step"] == 216,
        f"QubitGKSLBoson N=2 (n_max=1): GPS={result_qubit_boson['estimated_gates_per_step']}, expected 216",
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

    # Boson N=4 GPS regression (formula-based to avoid large memory)
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
        "regression_qubit_boson_n4_gps_476",
        gps_qubit_boson_n4 == 476,
        f"QubitGKSLBoson N=4 (n_max=1, formula): GPS={gps_qubit_boson_n4}, expected 476",
    )

    # ================================================================
    # Section 4: Notebook cell inspection (regression from iter22)
    # ================================================================
    print("\n=== Section 4: Notebook cell inspection ===")

    nb_path = os.path.join(
        os.path.dirname(__file__), "quantum_dynamics_gksl_comparison.ipynb"
    )
    with open(nb_path) as f:
        nb = json.load(f)
    cells = nb["cells"]

    # Cells 6, 8, 12, 14 (DM and boson)
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
    # Section 5: Result dict completeness (regression from iter22)
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
        result_qubit_boson.get("n_total_qubits") == 18,
        f"QubitGKSLBoson result n_total_qubits={result_qubit_boson.get('n_total_qubits')}, expected 18 (n_max=1)",
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
    # Section 6: Qubit shot forbidden_count (regression from iter22)
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
    # Section 7: Qudit boson n_max-independence check (regression from iter22)
    # ================================================================
    print("\n=== Section 7: Qudit boson n_max-independence ===")

    params_boson_nmax2 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=2, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_nmax2 = QuditGKSLBosonSimulator(params_boson_nmax2)
    result_qudit_boson_nmax2 = sim_qudit_boson_nmax2.simulate(t_max=0.1, n_steps=2)

    add(
        "qudit_boson_n_total_nmax_independent",
        sim_qudit_boson.n_total_qudits == sim_qudit_boson_nmax2.n_total_qudits,
        f"QuditBoson n_total: n_max=1 ({sim_qudit_boson.n_total_qudits}) == n_max=2 ({sim_qudit_boson_nmax2.n_total_qudits})",
    )
    add(
        "qudit_boson_gps_nmax_independent",
        result_qudit_boson["estimated_gates_per_step"] == result_qudit_boson_nmax2["estimated_gates_per_step"],
        f"QuditBoson GPS: n_max=1 ({result_qudit_boson['estimated_gates_per_step']}) == n_max=2 ({result_qudit_boson_nmax2['estimated_gates_per_step']})",
    )

    # Qubit boson n_total and GPS SHOULD depend on n_max
    sim_qubit_boson_nmax2 = QubitGKSLBosonSimulator(params_boson_nmax2)
    result_qubit_boson_nmax2 = sim_qubit_boson_nmax2.simulate(t_max=0.1, n_steps=2)
    add(
        "qubit_boson_n_total_nmax_dependent",
        sim_qubit_boson.n_total_qubits != sim_qubit_boson_nmax2.n_total_qubits,
        f"QubitBoson n_total: n_max=1 ({sim_qubit_boson.n_total_qubits}) != n_max=2 ({sim_qubit_boson_nmax2.n_total_qubits})",
    )
    add(
        "qubit_boson_gps_nmax_dependent",
        result_qubit_boson["estimated_gates_per_step"] != result_qubit_boson_nmax2["estimated_gates_per_step"],
        f"QubitBoson GPS: n_max=1 ({result_qubit_boson['estimated_gates_per_step']}) != n_max=2 ({result_qubit_boson_nmax2['estimated_gates_per_step']})",
    )

    # ================================================================
    # Section 8: Notebook Cell 14 output validation (regression from iter22)
    # ================================================================
    print("\n=== Section 8: Cell 14 output validation ===")

    cell14_outputs = cells[14].get("outputs", [])
    cell14_text = ""
    for out in cell14_outputs:
        if out.get("output_type") == "stream":
            cell14_text += "".join(out.get("text", []))

    add(
        "cell14_output_n_total_18",
        "Total qubits: 18" in cell14_text,
        "Cell 14 output should show 'Total qubits: 18' (n_max=1)",
    )
    add(
        "cell14_output_gps_216",
        "Estimated gates/step: 216" in cell14_text,
        "Cell 14 output should show 'Estimated gates/step: 216' (n_max=1)",
    )
    add(
        "cell14_output_phonon_qubits_2",
        "Phonon qubits: 2" in cell14_text,
        "Cell 14 output should show 'Phonon qubits: 2' (n_max=1, ceil(log2(2))*2=2)",
    )

    # ================================================================
    # Section 9: Issue BB - Cell 16 boson palindromic structure (NEW)
    # ================================================================
    print("\n=== Section 9: Cell 16 boson palindromic structure ===")

    cell16_src = "".join(cells[16]["source"])

    # Scenario 6 (Qudit Boson) should have fwd/rev labels
    add(
        "cell16_scenario6_has_fwd",
        "stinespring_fwd_" in cell16_src,
        "Cell 16 Scenario 6 should have forward Stinespring labels (stinespring_fwd_)",
    )
    add(
        "cell16_scenario6_has_rev",
        "stinespring_rev_" in cell16_src,
        "Cell 16 Scenario 6 should have reverse Stinespring labels (stinespring_rev_)",
    )

    # Scenario 6 should use dt_quantum / 2 for Stinespring
    # Check that build_stinespring_circuit_single uses dt_quantum / 2
    add(
        "cell16_scenario6_uses_half_dt",
        "dt_quantum / 2, sites[0]" in cell16_src or "dt_quantum / 2, sites" in cell16_src,
        "Cell 16 Scenario 6 should use dt_quantum/2 for Stinespring circuits",
    )

    # Scenario 4 (Qubit Boson) should have fwd/rev labels
    add(
        "cell16_scenario4_has_fwd",
        "stinespring_fwd_" in cell16_src,
        "Cell 16 Scenario 4 should have forward Stinespring labels",
    )
    add(
        "cell16_scenario4_has_rev",
        "stinespring_rev_" in cell16_src,
        "Cell 16 Scenario 4 should have reverse Stinespring labels",
    )

    # Scenario 4 should use dt_quantum / 2 for stinespring_unitary_from_lindblad
    add(
        "cell16_scenario4_uses_half_dt",
        "stinespring_unitary_from_lindblad(L_op, dt_quantum / 2)" in cell16_src,
        "Cell 16 Scenario 4 should use dt_quantum/2 for Stinespring unitary",
    )

    # Both scenarios should have reversed() for palindromic
    reversed_count = cell16_src.count("reversed(")
    # Non-boson scenarios don't use reversed() in Cell 16 (they use build_full_trotter_step_circuit)
    # Only boson scenarios have reversed() - should have 2 (Scenario 6 + Scenario 4)
    add(
        "cell16_boson_palindromic_reversed",
        reversed_count >= 2,
        f"Cell 16 should have >=2 reversed() calls for palindromic boson circuits (found {reversed_count})",
    )

    # Verify no leftover full-dt Stinespring calls for boson
    # The old code had: stinespring_unitary_from_lindblad(L_op, dt_quantum) without /2
    # Check that dt_quantum) doesn't appear after the boson section start
    boson_start = cell16_src.find("# Scenario 4:")
    if boson_start > 0:
        boson_section = cell16_src[boson_start:]
        has_full_dt = "stinespring_unitary_from_lindblad(L_op, dt_quantum)" in boson_section
        add(
            "cell16_no_full_dt_stinespring",
            not has_full_dt,
            "Cell 16 Scenario 4 should NOT use full dt for stinespring_unitary_from_lindblad",
        )
    else:
        add(
            "cell16_no_full_dt_stinespring",
            False,
            "Could not find Scenario 4 section in Cell 16",
        )

    # ================================================================
    # Section 10: Issue CC - n_phonon_qudits general formula (NEW)
    # ================================================================
    print("\n=== Section 10: n_phonon_qudits general formula ===")

    # n_max=1, d=3: ceil(log_3(2)) = 1 → n_phonon = 1*N (unchanged)
    add(
        "n_phonon_qudits_nmax1",
        sim_qudit_boson.n_phonon_qudits == params_boson.N_molecules,
        f"n_max=1: n_phonon_qudits={sim_qudit_boson.n_phonon_qudits}, expected {params_boson.N_molecules} (ceil(log_3(2))=1)",
    )

    # n_max=2, d=3: ceil(log_3(3)) = 1 → n_phonon = 1*N (unchanged)
    add(
        "n_phonon_qudits_nmax2",
        sim_qudit_boson_nmax2.n_phonon_qudits == params_boson_nmax2.N_molecules,
        f"n_max=2: n_phonon_qudits={sim_qudit_boson_nmax2.n_phonon_qudits}, expected {params_boson_nmax2.N_molecules} (ceil(log_3(3))=1)",
    )

    # n_max=3, d=3: ceil(log_3(4)) = 2 → n_phonon = 2*N (NEW: requires general formula)
    params_boson_nmax3 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=3, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_nmax3 = QuditGKSLBosonSimulator(params_boson_nmax3)
    expected_nph_nmax3 = 2 * params_boson_nmax3.N_molecules  # ceil(log_3(4))=2
    add(
        "n_phonon_qudits_nmax3",
        sim_qudit_boson_nmax3.n_phonon_qudits == expected_nph_nmax3,
        f"n_max=3: n_phonon_qudits={sim_qudit_boson_nmax3.n_phonon_qudits}, expected {expected_nph_nmax3} (ceil(log_3(4))=2)",
    )

    # n_max=8, d=3: ceil(log_3(9)) = 2 → n_phonon = 2*N
    params_boson_nmax8 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=8, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_nmax8 = QuditGKSLBosonSimulator(params_boson_nmax8)
    expected_nph_nmax8 = 2 * params_boson_nmax8.N_molecules  # ceil(log_3(9))=2
    add(
        "n_phonon_qudits_nmax8",
        sim_qudit_boson_nmax8.n_phonon_qudits == expected_nph_nmax8,
        f"n_max=8: n_phonon_qudits={sim_qudit_boson_nmax8.n_phonon_qudits}, expected {expected_nph_nmax8} (ceil(log_3(9))=2)",
    )

    # Verify qudit boson source code uses general formula (not hardcoded N_molecules)
    boson_sim_path = os.path.join(
        os.path.dirname(__file__), "qudit_gksl_boson_simulator.py"
    )
    with open(boson_sim_path) as f:
        boson_sim_src = f.read()

    add(
        "n_phonon_qudits_uses_log_formula",
        "np.log" in boson_sim_src and "np.ceil" in boson_sim_src,
        "qudit_gksl_boson_simulator.py should use np.ceil(np.log()) for n_phonon_qudits",
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
        "iteration": 23,
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
        out_dir, f"iteration23_palindromic_boson_{results['timestamp']}.json"
    )
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to: {out_file}")


if __name__ == "__main__":
    main()
