#!/usr/bin/env python3
"""Iteration 25 verification script.

Checks:
- Issue FF fix: QuditGKSLCircuitBosonSimulator.traces now uses rho_el
  (electronic-only partial trace) instead of full rho, matching
  QuditGKSLBosonSimulator.
- Issue GG fix: QuditGKSLCircuitBosonSimulator precomputes unitaries
  once per simulation (like DM simulator), not every Trotter step.
- Regression: All iteration 24 checks (74 checks) still pass.

Usage:
    cd tutorials
    python run_iteration25_verification.py
"""

from __future__ import annotations

import json
import os
import sys
import time as time_module
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main() -> None:
    checks: list[dict] = []

    def add(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {detail}")

    # ================================================================
    # Section 1-12: Iteration 24 regression checks (74 checks)
    # ================================================================
    import numpy as np

    from gksl_physical_parameters import GKSLPhysicalParameters

    # --- Section 1: n_total correctness ---
    print("\n=== Section 1: n_total in simulator classes ===")

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

    params_boson_n4 = GKSLPhysicalParameters(
        N_molecules=4, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_n4 = QuditGKSLBosonSimulator(params_boson_n4)
    add(
        "qudit_boson_n4_n_total_value",
        sim_qudit_boson_n4.n_total_qudits == 34,
        f"QuditGKSLBoson N=4: n_total_qudits={sim_qudit_boson_n4.n_total_qudits}, expected 34 (4+4+26)",
    )

    n_el_qubits_n4 = 2 * params_boson_n4.N_molecules
    n_ph_qubits_n4 = int(np.ceil(np.log2(params_boson_n4.n_max + 1))) * params_boson_n4.N_molecules
    n_ancilla_qubit_boson_n4 = 2 * len(params_boson_n4.neighbors) + 5 * params_boson_n4.N_molecules
    n_total_qubit_boson_n4 = n_el_qubits_n4 + n_ph_qubits_n4 + n_ancilla_qubit_boson_n4
    add(
        "qubit_boson_n4_n_total_value",
        n_total_qubit_boson_n4 == 38,
        f"QubitGKSLBoson N=4 (n_max=1, formula): n_total_qubits={n_total_qubit_boson_n4}, expected 38 (8+4+26)",
    )

    # --- Section 2: Cross-consistency ---
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

    # --- Section 3: GPS regression ---
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

    # --- Section 4: Notebook cell inspection ---
    print("\n=== Section 4: Notebook cell inspection ===")

    nb_path = os.path.join(
        os.path.dirname(__file__), "quantum_dynamics_gksl_comparison.ipynb"
    )
    with open(nb_path) as f:
        nb = json.load(f)
    cells = nb["cells"]

    cell6_src = "".join(cells[6]["source"])
    add("cell6_has_total_qudits", "n_total_qudits" in cell6_src,
        "Cell 6 (QuditGKSL) should display n_total_qudits")
    cell8_src = "".join(cells[8]["source"])
    add("cell8_has_total_qubits", "n_total_qubits" in cell8_src,
        "Cell 8 (QubitGKSL) should display n_total_qubits")
    cell12_src = "".join(cells[12]["source"])
    add("cell12_has_total_qudits", "n_total_qudits" in cell12_src,
        "Cell 12 (QuditGKSLBoson) should display n_total_qudits")
    cell14_src = "".join(cells[14]["source"])
    add("cell14_has_total_qubits", "n_total_qubits" in cell14_src,
        "Cell 14 (QubitGKSLBoson) should display n_total_qubits")

    cell32_src = "".join(cells[32]["source"])
    add("cell32_has_total_qudits", "n_total_qudits" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display n_total_qudits")
    add("cell32_has_system_qudits", "n_system_qudits" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display n_system_qudits")
    add("cell32_has_ancilla_qudits", "n_ancilla_qudits" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display n_ancilla_qudits")
    add("cell32_has_counts", "counts" in cell32_src,
        "Cell 32 (QuditGKSLShot) should display measurement counts")

    cell36_src = "".join(cells[36]["source"])
    add("cell36_has_total_qubits", "n_total_qubits" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display n_total_qubits")
    add("cell36_has_ancilla", "n_ancilla" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display n_ancilla")
    add("cell36_has_trace", "Trace conservation" in cell36_src or "trace" in cell36_src.lower(),
        "Cell 36 (QubitGKSLShot) should display trace conservation")
    add("cell36_has_counts", "counts" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display measurement counts")
    add("cell36_has_forbidden_count", "forbidden_count" in cell36_src,
        "Cell 36 (QubitGKSLShot) should display forbidden_count")

    # --- Section 5: Result dict completeness ---
    print("\n=== Section 5: Result dict completeness ===")

    add("qudit_dm_result_n_total", result_qudit.get("n_total_qudits") == 30,
        f"QuditGKSL result n_total_qudits={result_qudit.get('n_total_qudits')}, expected 30")
    add("qubit_dm_result_n_total", result_qubit.get("n_total_qubits") == 34,
        f"QubitGKSL result n_total_qubits={result_qubit.get('n_total_qubits')}, expected 34")
    add("qudit_boson_result_n_total", result_qudit_boson.get("n_total_qudits") == 16,
        f"QuditGKSLBoson result n_total_qudits={result_qudit_boson.get('n_total_qudits')}, expected 16")
    add("qubit_boson_result_n_total", result_qubit_boson.get("n_total_qubits") == 18,
        f"QubitGKSLBoson result n_total_qubits={result_qubit_boson.get('n_total_qubits')}, expected 18 (n_max=1)")
    add("qudit_shot_result_n_total", result_qudit_shot.get("n_total_qudits") == 30,
        f"QuditGKSLShot result n_total_qudits={result_qudit_shot.get('n_total_qudits')}, expected 30")
    add("qubit_shot_result_n_total", result_qubit_shot.get("n_total_qubits") == 34,
        f"QubitGKSLShot result n_total_qubits={result_qubit_shot.get('n_total_qubits')}, expected 34")

    # --- Section 6: Qubit shot forbidden_count ---
    print("\n=== Section 6: Qubit shot forbidden_count ===")

    add("qubit_shot_forbidden_count_zero", result_qubit_shot.get("forbidden_count") == 0,
        f"QubitGKSLShot noiseless: forbidden_count={result_qubit_shot.get('forbidden_count')}, expected 0")
    add("qubit_shot_has_counts_key", "counts" in result_qubit_shot,
        "QubitGKSLShot result dict should contain 'counts' key")
    add("qubit_shot_has_forbidden_key", "forbidden_count" in result_qubit_shot,
        "QubitGKSLShot result dict should contain 'forbidden_count' key")

    # --- Section 7: n_max-independence ---
    print("\n=== Section 7: Qudit boson n_max-independence ===")

    params_boson_nmax2 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=2, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_nmax2 = QuditGKSLBosonSimulator(params_boson_nmax2)
    result_qudit_boson_nmax2 = sim_qudit_boson_nmax2.simulate(t_max=0.1, n_steps=2)

    add("qudit_boson_n_total_nmax_independent",
        sim_qudit_boson.n_total_qudits == sim_qudit_boson_nmax2.n_total_qudits,
        f"QuditBoson n_total: n_max=1 ({sim_qudit_boson.n_total_qudits}) == n_max=2 ({sim_qudit_boson_nmax2.n_total_qudits})")
    add("qudit_boson_gps_nmax_independent",
        result_qudit_boson["estimated_gates_per_step"] == result_qudit_boson_nmax2["estimated_gates_per_step"],
        f"QuditBoson GPS: n_max=1 ({result_qudit_boson['estimated_gates_per_step']}) == n_max=2 ({result_qudit_boson_nmax2['estimated_gates_per_step']})")

    sim_qubit_boson_nmax2 = QubitGKSLBosonSimulator(params_boson_nmax2)
    result_qubit_boson_nmax2 = sim_qubit_boson_nmax2.simulate(t_max=0.1, n_steps=2)
    add("qubit_boson_n_total_nmax_dependent",
        sim_qubit_boson.n_total_qubits != sim_qubit_boson_nmax2.n_total_qubits,
        f"QubitBoson n_total: n_max=1 ({sim_qubit_boson.n_total_qubits}) != n_max=2 ({sim_qubit_boson_nmax2.n_total_qubits})")
    add("qubit_boson_gps_nmax_dependent",
        result_qubit_boson["estimated_gates_per_step"] != result_qubit_boson_nmax2["estimated_gates_per_step"],
        f"QubitBoson GPS: n_max=1 ({result_qubit_boson['estimated_gates_per_step']}) != n_max=2 ({result_qubit_boson_nmax2['estimated_gates_per_step']})")

    # --- Section 8: Cell 14 output ---
    print("\n=== Section 8: Cell 14 output validation ===")

    cell14_outputs = cells[14].get("outputs", [])
    cell14_text = ""
    for out in cell14_outputs:
        if out.get("output_type") == "stream":
            cell14_text += "".join(out.get("text", []))

    add("cell14_output_n_total_18", "Total qubits: 18" in cell14_text,
        "Cell 14 output should show 'Total qubits: 18' (n_max=1)")
    add("cell14_output_gps_216", "Estimated gates/step: 216" in cell14_text,
        "Cell 14 output should show 'Estimated gates/step: 216' (n_max=1)")
    add("cell14_output_phonon_qubits_2", "Phonon qubits: 2" in cell14_text,
        "Cell 14 output should show 'Phonon qubits: 2' (n_max=1, ceil(log2(2))*2=2)")

    # --- Section 9: Cell 16 palindromic (iteration 23) ---
    print("\n=== Section 9: Cell 16 boson palindromic structure ===")

    cell16_src = "".join(cells[16]["source"])

    add("cell16_scenario6_has_fwd", "stinespring_fwd_" in cell16_src,
        "Cell 16 Scenario 6 should have forward Stinespring labels (stinespring_fwd_)")
    add("cell16_scenario6_has_rev", "stinespring_rev_" in cell16_src,
        "Cell 16 Scenario 6 should have reverse Stinespring labels (stinespring_rev_)")
    add("cell16_scenario6_uses_half_dt",
        "dt_quantum / 2, sites[0]" in cell16_src or "dt_quantum / 2, sites" in cell16_src,
        "Cell 16 Scenario 6 should use dt_quantum/2 for Stinespring circuits")
    add("cell16_scenario4_has_fwd", "stinespring_fwd_" in cell16_src,
        "Cell 16 Scenario 4 should have forward Stinespring labels")
    add("cell16_scenario4_has_rev", "stinespring_rev_" in cell16_src,
        "Cell 16 Scenario 4 should have reverse Stinespring labels")
    add("cell16_scenario4_uses_half_dt",
        "stinespring_unitary_from_lindblad(L_op, dt_quantum / 2)" in cell16_src,
        "Cell 16 Scenario 4 should use dt_quantum/2 for Stinespring unitary")

    reversed_count = cell16_src.count("reversed(")
    add("cell16_boson_palindromic_reversed", reversed_count >= 2,
        f"Cell 16 should have >=2 reversed() calls for palindromic boson circuits (found {reversed_count})")

    boson_start = cell16_src.find("# Scenario 4:")
    if boson_start > 0:
        boson_section = cell16_src[boson_start:]
        has_full_dt = "stinespring_unitary_from_lindblad(L_op, dt_quantum)" in boson_section
        add("cell16_no_full_dt_stinespring", not has_full_dt,
            "Cell 16 Scenario 4 should NOT use full dt for stinespring_unitary_from_lindblad")
    else:
        add("cell16_no_full_dt_stinespring", False,
            "Could not find Scenario 4 section in Cell 16")

    # --- Section 10: n_phonon_qudits formula (iteration 23) ---
    print("\n=== Section 10: n_phonon_qudits general formula ===")

    add("n_phonon_qudits_nmax1",
        sim_qudit_boson.n_phonon_qudits == params_boson.N_molecules,
        f"n_max=1: n_phonon_qudits={sim_qudit_boson.n_phonon_qudits}, expected {params_boson.N_molecules} (ceil(log_3(2))=1)")
    add("n_phonon_qudits_nmax2",
        sim_qudit_boson_nmax2.n_phonon_qudits == params_boson_nmax2.N_molecules,
        f"n_max=2: n_phonon_qudits={sim_qudit_boson_nmax2.n_phonon_qudits}, expected {params_boson_nmax2.N_molecules} (ceil(log_3(3))=1)")

    params_boson_nmax3 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=3, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_nmax3 = QuditGKSLBosonSimulator(params_boson_nmax3)
    expected_nph_nmax3 = 2 * params_boson_nmax3.N_molecules
    add("n_phonon_qudits_nmax3",
        sim_qudit_boson_nmax3.n_phonon_qudits == expected_nph_nmax3,
        f"n_max=3: n_phonon_qudits={sim_qudit_boson_nmax3.n_phonon_qudits}, expected {expected_nph_nmax3} (ceil(log_3(4))=2)")

    params_boson_nmax8 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=8, omega_ph=0.15, g_eph=0.02,
    )
    sim_qudit_boson_nmax8 = QuditGKSLBosonSimulator(params_boson_nmax8)
    expected_nph_nmax8 = 2 * params_boson_nmax8.N_molecules
    add("n_phonon_qudits_nmax8",
        sim_qudit_boson_nmax8.n_phonon_qudits == expected_nph_nmax8,
        f"n_max=8: n_phonon_qudits={sim_qudit_boson_nmax8.n_phonon_qudits}, expected {expected_nph_nmax8} (ceil(log_3(9))=2)")

    boson_sim_path = os.path.join(
        os.path.dirname(__file__), "qudit_gksl_boson_simulator.py"
    )
    with open(boson_sim_path) as f:
        boson_sim_src = f.read()
    add("n_phonon_qudits_uses_log_formula",
        "np.log" in boson_sim_src and "np.ceil" in boson_sim_src,
        "qudit_gksl_boson_simulator.py should use np.ceil(np.log()) for n_phonon_qudits")

    # --- Section 11: Circuit boson rho_final consistency (iteration 24) ---
    print("\n=== Section 11: Circuit boson rho_final consistency ===")

    from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

    sim_circ_boson = QuditGKSLCircuitBosonSimulator(params_boson)
    result_circ_boson = sim_circ_boson.simulate(
        t_max=0.1, n_steps=2, initial_state="edge_triplet"
    )

    dim_el = params_boson.d ** params_boson.N_molecules
    rho_circ = result_circ_boson["rho_final"]
    add(
        "circuit_boson_rho_final_shape",
        rho_circ.shape == (dim_el, dim_el),
        f"Circuit boson rho_final shape should be ({dim_el},{dim_el}), got {rho_circ.shape}",
    )

    result_dm_boson_check = sim_qudit_boson.simulate(t_max=0.1, n_steps=2)
    rho_dm = result_dm_boson_check["rho_final"]
    diff_norm = float(np.linalg.norm(rho_circ - rho_dm))
    add(
        "circuit_boson_rho_matches_dm",
        diff_norm < 1e-12,
        f"||rho_circuit - rho_dm||_F = {diff_norm:.2e} (should be < 1e-12)",
    )

    circ_trace = float(np.real(np.trace(rho_circ)))
    add(
        "circuit_boson_rho_final_trace",
        abs(circ_trace - 1.0) < 1e-10,
        f"Circuit boson rho_final trace = {circ_trace:.10f} (should be ~1.0)",
    )

    # --- Section 12: Circuit boson result dict key consistency (iteration 24) ---
    print("\n=== Section 12: Circuit boson result dict key consistency ===")

    add(
        "circuit_boson_has_n_system_qudits",
        "n_system_qudits" in result_circ_boson,
        f"Circuit boson should have 'n_system_qudits' key (value={result_circ_boson.get('n_system_qudits')})",
    )
    add(
        "circuit_boson_has_n_phonon_qudits",
        "n_phonon_qudits" in result_circ_boson,
        f"Circuit boson should have 'n_phonon_qudits' key (value={result_circ_boson.get('n_phonon_qudits')})",
    )
    add(
        "circuit_boson_has_n_total_qudits",
        "n_total_qudits" in result_circ_boson,
        f"Circuit boson should have 'n_total_qudits' key (value={result_circ_boson.get('n_total_qudits')})",
    )
    add(
        "circuit_boson_has_estimated_gps",
        "estimated_gates_per_step" in result_circ_boson,
        f"Circuit boson should have 'estimated_gates_per_step' key (value={result_circ_boson.get('estimated_gates_per_step')})",
    )

    add(
        "circuit_boson_n_system_matches_dm",
        result_circ_boson.get("n_system_qudits") == result_dm_boson_check.get("n_system_qudits"),
        f"n_system_qudits: circuit={result_circ_boson.get('n_system_qudits')}, dm={result_dm_boson_check.get('n_system_qudits')}",
    )
    add(
        "circuit_boson_n_phonon_matches_dm",
        result_circ_boson.get("n_phonon_qudits") == result_dm_boson_check.get("n_phonon_qudits"),
        f"n_phonon_qudits: circuit={result_circ_boson.get('n_phonon_qudits')}, dm={result_dm_boson_check.get('n_phonon_qudits')}",
    )
    add(
        "circuit_boson_n_total_matches_dm",
        result_circ_boson.get("n_total_qudits") == result_dm_boson_check.get("n_total_qudits"),
        f"n_total_qudits: circuit={result_circ_boson.get('n_total_qudits')}, dm={result_dm_boson_check.get('n_total_qudits')}",
    )
    add(
        "circuit_boson_gps_matches_dm",
        result_circ_boson.get("estimated_gates_per_step") == result_dm_boson_check.get("estimated_gates_per_step"),
        f"GPS: circuit={result_circ_boson.get('estimated_gates_per_step')}, dm={result_dm_boson_check.get('estimated_gates_per_step')}",
    )

    add(
        "circuit_boson_no_old_keys",
        "n_el_qutrits" not in result_circ_boson and "n_ph_qutrits" not in result_circ_boson,
        "Circuit boson should not have old keys 'n_el_qutrits'/'n_ph_qutrits'",
    )

    # ================================================================
    # Section 13: Issue FF - Circuit boson traces consistency (NEW)
    # ================================================================
    print("\n=== Section 13: Issue FF - Circuit boson traces consistency ===")

    # Verify traces are computed from rho_el (electronic-only), not full rho
    # Both DM and circuit should give trace ~1.0 from rho_el
    dm_traces = result_dm_boson_check["trace"]
    circ_traces = result_circ_boson["trace"]

    add(
        "circuit_boson_trace_length_matches_dm",
        len(circ_traces) == len(dm_traces),
        f"Circuit traces length ({len(circ_traces)}) should match DM ({len(dm_traces)})",
    )

    # All traces should be very close to 1.0 (since Tr(rho_el) = 1.0 for valid density matrix)
    max_trace_dev_circ = max(abs(t - 1.0) for t in circ_traces)
    max_trace_dev_dm = max(abs(t - 1.0) for t in dm_traces)
    add(
        "circuit_boson_trace_near_unity",
        max_trace_dev_circ < 1e-10,
        f"Circuit max |Tr-1| = {max_trace_dev_circ:.2e} (should be < 1e-10)",
    )

    # Cross-compare: circuit traces should match DM traces closely
    max_trace_diff = max(abs(c - d) for c, d in zip(circ_traces, dm_traces))
    add(
        "circuit_boson_traces_match_dm",
        max_trace_diff < 1e-12,
        f"max |trace_circuit - trace_dm| = {max_trace_diff:.2e} (should be < 1e-12)",
    )

    # Source code verification: traces should use rho_el, not full rho
    circ_boson_path = os.path.join(
        os.path.dirname(__file__), "qudit_gksl_circuit_boson_simulator.py"
    )
    with open(circ_boson_path) as f:
        circ_boson_src = f.read()

    # Find the simulate method's trace computation lines
    # Should use np.trace(rho_el), NOT np.trace(rho) in the simulation loop
    simulate_start = circ_boson_src.find("def simulate(")
    if simulate_start >= 0:
        simulate_body = circ_boson_src[simulate_start:]
        # Count np.trace(rho_el) vs np.trace(rho) in simulate method
        trace_rho_el_count = simulate_body.count("np.trace(rho_el)")
        trace_rho_only_count = simulate_body.count("np.trace(rho))")
        add(
            "circuit_boson_traces_use_rho_el",
            trace_rho_el_count >= 2 and trace_rho_only_count == 0,
            f"simulate() should use np.trace(rho_el) (found {trace_rho_el_count}) "
            f"and NOT np.trace(rho) (found {trace_rho_only_count})",
        )
    else:
        add(
            "circuit_boson_traces_use_rho_el",
            False,
            "Could not find simulate() method in circuit boson simulator",
        )

    # ================================================================
    # Section 14: Issue GG - Circuit boson precomputation (NEW)
    # ================================================================
    print("\n=== Section 14: Issue GG - Circuit boson precomputation ===")

    # Verify that _precompute_unitaries method exists
    add(
        "circuit_boson_has_precompute",
        hasattr(sim_circ_boson, "_precompute_unitaries"),
        "Circuit boson should have _precompute_unitaries method",
    )

    # Verify that precomputed values are set after simulate()
    add(
        "circuit_boson_has_U_H_half",
        hasattr(sim_circ_boson, "_U_H_half"),
        "After simulate(), circuit boson should have _U_H_half attribute",
    )
    add(
        "circuit_boson_has_kraus_list",
        hasattr(sim_circ_boson, "_kraus_list"),
        "After simulate(), circuit boson should have _kraus_list attribute",
    )

    # Verify _U_H_half is a unitary matrix of correct dimension
    if hasattr(sim_circ_boson, "_U_H_half"):
        U_H = sim_circ_boson._U_H_half
        dim_total = sim_circ_boson.dim_total
        add(
            "circuit_boson_U_H_half_shape",
            U_H.shape == (dim_total, dim_total),
            f"_U_H_half shape should be ({dim_total},{dim_total}), got {U_H.shape}",
        )
        unitarity_err = np.linalg.norm(U_H.conj().T @ U_H - np.eye(dim_total))
        add(
            "circuit_boson_U_H_half_unitary",
            unitarity_err < 1e-12,
            f"_U_H_half unitarity error: {unitarity_err:.2e} (should be < 1e-12)",
        )
    else:
        add("circuit_boson_U_H_half_shape", False, "_U_H_half not found")
        add("circuit_boson_U_H_half_unitary", False, "_U_H_half not found")

    # Verify _kraus_list has correct length (= number of Lindblad channels)
    if hasattr(sim_circ_boson, "_kraus_list"):
        n_lindblad = len(sim_circ_boson.lindblad_local_info)
        add(
            "circuit_boson_kraus_list_length",
            len(sim_circ_boson._kraus_list) == n_lindblad,
            f"_kraus_list length should be {n_lindblad}, got {len(sim_circ_boson._kraus_list)}",
        )
    else:
        add("circuit_boson_kraus_list_length", False, "_kraus_list not found")

    # Source code verification: _trotter_step should NOT call
    # _build_local_stinespring_unitary or _compute_hamiltonian_unitary
    trotter_start = circ_boson_src.find("def _trotter_step(")
    if trotter_start >= 0:
        # Find the end of the method (next def at same indentation)
        next_def = circ_boson_src.find("\n    def ", trotter_start + 1)
        if next_def > 0:
            trotter_body = circ_boson_src[trotter_start:next_def]
        else:
            trotter_body = circ_boson_src[trotter_start:]

        has_stinespring_call = "_build_local_stinespring_unitary" in trotter_body
        has_hamiltonian_call = "_compute_hamiltonian_unitary" in trotter_body
        add(
            "circuit_boson_trotter_no_recompute",
            not has_stinespring_call and not has_hamiltonian_call,
            f"_trotter_step should NOT recompute unitaries "
            f"(stinespring={has_stinespring_call}, hamiltonian={has_hamiltonian_call})",
        )

        # Verify that _trotter_step uses precomputed _U_H_half
        uses_precomputed = "_U_H_half" in trotter_body and "_kraus_list" in trotter_body
        add(
            "circuit_boson_trotter_uses_precomputed",
            uses_precomputed,
            "_trotter_step should use precomputed _U_H_half and _kraus_list",
        )
    else:
        add("circuit_boson_trotter_no_recompute", False,
            "Could not find _trotter_step method")
        add("circuit_boson_trotter_uses_precomputed", False,
            "Could not find _trotter_step method")

    # Verify precompute is called in simulate, not in _trotter_step
    if simulate_start >= 0:
        simulate_section = circ_boson_src[simulate_start:]
        has_precompute_call = "_precompute_unitaries" in simulate_section
        add(
            "circuit_boson_simulate_calls_precompute",
            has_precompute_call,
            "simulate() should call _precompute_unitaries(dt)",
        )
    else:
        add("circuit_boson_simulate_calls_precompute", False,
            "Could not find simulate() method")

    # ================================================================
    # Section 15: Performance validation (NEW)
    # ================================================================
    print("\n=== Section 15: Performance comparison ===")

    # Run both simulators with more steps and compare timing
    n_perf_steps = 20
    t_perf = 1.0

    t0 = time_module.time()
    result_dm_perf = sim_qudit_boson.simulate(t_max=t_perf, n_steps=n_perf_steps)
    t_dm = time_module.time() - t0

    t0 = time_module.time()
    result_circ_perf = sim_circ_boson.simulate(t_max=t_perf, n_steps=n_perf_steps)
    t_circ = time_module.time() - t0

    # Verify numerical agreement for longer simulation
    rho_dm_perf = result_dm_perf["rho_final"]
    rho_circ_perf = result_circ_perf["rho_final"]
    diff_norm_perf = float(np.linalg.norm(rho_circ_perf - rho_dm_perf))
    add(
        "circuit_boson_long_run_matches_dm",
        diff_norm_perf < 1e-10,
        f"||rho_circuit - rho_dm||_F after {n_perf_steps} steps = {diff_norm_perf:.2e}",
    )

    # Log timing (not a pass/fail, just informational)
    add(
        "circuit_boson_timing_logged",
        True,
        f"DM: {t_dm:.3f}s, Circuit: {t_circ:.3f}s for {n_perf_steps} steps",
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
        "iteration": 25,
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
        out_dir, f"iteration25_trace_precompute_{results['timestamp']}.json"
    )
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to: {out_file}")


if __name__ == "__main__":
    main()
