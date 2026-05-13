"""Iteration 16 verification: estimated_gates_per_step consistency across ALL simulators.

Checks that shot-based and boson simulators' estimated_gates_per_step correctly
reflects the palindromic 2nd-order Trotter step (2× for Hamiltonian half-steps),
matching the density matrix simulators fixed in iteration 15.

Bugs fixed in iteration 16:
  Bug L: qudit_gksl_shot_simulator.py - estimated_gates_per_step was 59 (should be 66)
  Bug M: qubit_gksl_shot_simulator.py - estimated_gates_per_step was 350 (should be 388)
  Bug N: qudit_gksl_boson_simulator.py - missing 2× for Hamiltonian AND Stinespring
  Bug O: qubit_gksl_boson_simulator.py - missing 2× for Hamiltonian AND Stinespring

Run: cd tutorials && python run_iteration16_verification.py
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_checks() -> list[dict]:
    checks: list[dict] = []

    def add(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"name": name, "passed": passed, "detail": detail})

    from gksl_physical_parameters import GKSLPhysicalParameters

    params = GKSLPhysicalParameters()
    N = params.N_molecules  # 4
    n_pairs = len(params.neighbors)  # 3

    n_lindblad = 2 * n_pairs + 5 * N  # 6 + 20 = 26

    # ================================================================
    # Section 1: Shot-based qudit simulator (Bug L fix)
    # ================================================================
    from qudit_gksl_shot_simulator import QuditGKSLShotSimulator

    sim_shot_qudit = QuditGKSLShotSimulator(params)
    r_shot_qudit = sim_shot_qudit.simulate(
        t_max=1.0, n_steps=2, n_shots=10, seed=42
    )

    expected_qudit_gps = 2 * (N + n_pairs) + n_lindblad * 2  # 14 + 52 = 66
    actual_shot_qudit_gps = r_shot_qudit["estimated_gates_per_step"]

    add(
        "qudit_shot_gps_value",
        actual_shot_qudit_gps == expected_qudit_gps,
        f"Expected {expected_qudit_gps}, got {actual_shot_qudit_gps}",
    )

    old_broken_qudit = N + n_pairs + n_lindblad * 2  # 59
    add(
        "qudit_shot_not_old_59",
        actual_shot_qudit_gps != old_broken_qudit,
        f"Should not be old value {old_broken_qudit}",
    )

    # ================================================================
    # Section 2: Shot-based qubit simulator (Bug M fix)
    # ================================================================
    from qubit_gksl_shot_simulator import QubitGKSLShotSimulator

    sim_shot_qubit = QubitGKSLShotSimulator(params)
    r_shot_qubit = sim_shot_qubit.simulate(
        t_max=1.0, n_steps=2, n_shots=10, seed=42
    )

    n_sys_qubits = 2 * N  # 8
    n_ancilla = n_lindblad  # 26
    expected_qubit_gps = 2 * (n_sys_qubits + n_pairs * 10) + n_ancilla * 6 * 2  # 388
    actual_shot_qubit_gps = r_shot_qubit["estimated_gates_per_step"]

    add(
        "qubit_shot_gps_value",
        actual_shot_qubit_gps == expected_qubit_gps,
        f"Expected {expected_qubit_gps}, got {actual_shot_qubit_gps}",
    )

    old_broken_qubit = n_sys_qubits + n_pairs * 10 + n_ancilla * 6 * 2  # 350
    add(
        "qubit_shot_not_old_350",
        actual_shot_qubit_gps != old_broken_qubit,
        f"Should not be old value {old_broken_qubit}",
    )

    # ================================================================
    # Section 3: Cross-check shot vs density matrix simulators
    # ================================================================
    from qudit_gksl_simulator import QuditGKSLSimulator
    from qubit_gksl_simulator import QubitGKSLSimulator

    sim_dm_qudit = QuditGKSLSimulator(params)
    r_dm_qudit = sim_dm_qudit.simulate(t_max=1.0, n_steps=2)

    sim_dm_qubit = QubitGKSLSimulator(params)
    r_dm_qubit = sim_dm_qubit.simulate(t_max=1.0, n_steps=2)

    add(
        "qudit_shot_matches_dm",
        actual_shot_qudit_gps == r_dm_qudit["estimated_gates_per_step"],
        f"Shot {actual_shot_qudit_gps} vs DM {r_dm_qudit['estimated_gates_per_step']}",
    )

    add(
        "qubit_shot_matches_dm",
        actual_shot_qubit_gps == r_dm_qubit["estimated_gates_per_step"],
        f"Shot {actual_shot_qubit_gps} vs DM {r_dm_qubit['estimated_gates_per_step']}",
    )

    # ================================================================
    # Section 4: Noisy shot simulators inherit correct values
    # ================================================================
    from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator
    from qubit_gksl_shot_simulator import QubitGKSLNoisyShotSimulator

    sim_noisy_qudit = QuditGKSLNoisyShotSimulator(
        params, p_depol=0.001, p_dephasing=0.0, depol_pair_only=True
    )
    r_noisy_qudit = sim_noisy_qudit.simulate(
        t_max=1.0, n_steps=2, n_shots=10, seed=42
    )

    add(
        "qudit_noisy_shot_gps",
        r_noisy_qudit["estimated_gates_per_step"] == expected_qudit_gps,
        f"Expected {expected_qudit_gps}, got {r_noisy_qudit['estimated_gates_per_step']}",
    )

    sim_noisy_qubit = QubitGKSLNoisyShotSimulator(
        params, p_depol=0.001, p_dephasing=0.0, depol_pair_only=True
    )
    r_noisy_qubit = sim_noisy_qubit.simulate(
        t_max=1.0, n_steps=2, n_shots=10, seed=42
    )

    add(
        "qubit_noisy_shot_gps",
        r_noisy_qubit["estimated_gates_per_step"] == expected_qubit_gps,
        f"Expected {expected_qubit_gps}, got {r_noisy_qubit['estimated_gates_per_step']}",
    )

    # ================================================================
    # Section 5: Boson simulators formula verification (static check)
    # (Cannot run full simulation due to memory; verify formula only)
    # ================================================================
    from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
    from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator

    import inspect

    # Check qudit boson simulator has 2× factor in source
    src_qudit_boson = inspect.getsource(QuditGKSLBosonSimulator.simulate)
    has_2x_qudit_boson = "2 * (" in src_qudit_boson and "n_ancilla_qudits * 2" in src_qudit_boson
    add(
        "qudit_boson_formula_has_2x",
        has_2x_qudit_boson,
        "Source should contain '2 * (' and 'n_ancilla_qudits * 2' for palindromic Trotter",
    )

    # Check qubit boson simulator has 2× factor in source
    src_qubit_boson = inspect.getsource(QubitGKSLBosonSimulator.simulate)
    has_2x_qubit_boson = "2 * (" in src_qubit_boson and "n_ancilla * 6 * 2" in src_qubit_boson
    add(
        "qubit_boson_formula_has_2x",
        has_2x_qubit_boson,
        "Source should contain '2 * (' and 'n_ancilla * 6 * 2' for palindromic Trotter",
    )

    # ================================================================
    # Section 6: total_estimated_gates consistency
    # ================================================================
    n_steps = 2
    add(
        "qudit_shot_total_gates",
        r_shot_qudit["total_estimated_gates"] == expected_qudit_gps * n_steps,
        f"Expected {expected_qudit_gps * n_steps}, got {r_shot_qudit['total_estimated_gates']}",
    )
    add(
        "qubit_shot_total_gates",
        r_shot_qubit["total_estimated_gates"] == expected_qubit_gps * n_steps,
        f"Expected {expected_qubit_gps * n_steps}, got {r_shot_qubit['total_estimated_gates']}",
    )

    # ================================================================
    # Section 7: Iteration 15 regression checks (density matrix sims)
    # ================================================================
    add(
        "regression_qudit_dm_gps_66",
        r_dm_qudit["estimated_gates_per_step"] == 66,
        f"Expected 66, got {r_dm_qudit['estimated_gates_per_step']}",
    )
    add(
        "regression_qubit_dm_gps_388",
        r_dm_qubit["estimated_gates_per_step"] == 388,
        f"Expected 388, got {r_dm_qubit['estimated_gates_per_step']}",
    )

    # ================================================================
    # Section 8: Circuit simulator regression
    # ================================================================
    from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator
    from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

    sim_circ_qudit = QuditGKSLCircuitSimulator(params)
    r_circ_qudit = sim_circ_qudit.simulate(t_max=1.0, n_steps=2)

    add(
        "regression_qudit_circuit_gps_66",
        r_circ_qudit["gates_per_step"] == 66,
        f"Expected 66, got {r_circ_qudit['gates_per_step']}",
    )

    sim_circ_qubit = QubitGKSLCircuitSimulator(params)
    r_circ_qubit = sim_circ_qubit.simulate(t_max=1.0, n_steps=2)

    add(
        "regression_qubit_circuit_gps_66",
        r_circ_qubit["gates_per_step"] == 66,
        f"Expected 66, got {r_circ_qubit['gates_per_step']}",
    )

    # ================================================================
    # Section 9: Trace preservation regression
    # ================================================================
    add(
        "regression_qudit_shot_trace",
        abs(r_shot_qudit["trace"][-1] - 1.0) < 1e-6,
        f"Final trace = {r_shot_qudit['trace'][-1]}",
    )
    add(
        "regression_qubit_shot_trace",
        abs(r_shot_qubit["trace"][-1] - 1.0) < 1e-6,
        f"Final trace = {r_shot_qubit['trace'][-1]}",
    )

    return checks


def main() -> None:
    checks = run_checks()
    n_pass = sum(1 for c in checks if c["passed"])
    n_total = len(checks)

    print(f"\n{'='*60}")
    print(f"Iteration 16 Verification: {n_pass}/{n_total} checks passed")
    print(f"{'='*60}")
    for c in checks:
        status = "PASS" if c["passed"] else "FAIL"
        print(f"  [{status}] {c['name']}: {c['detail']}")

    # Save results
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "developing",
        "verification_results",
    )
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"iteration16_gate_count_all_sims_{ts}.json")
    with open(out_path, "w") as f:
        json.dump(
            {
                "iteration": 16,
                "timestamp": ts,
                "total_checks": n_total,
                "passed": n_pass,
                "failed": n_total - n_pass,
                "all_passed": n_pass == n_total,
                "checks": checks,
            },
            f,
            indent=2,
        )
    print(f"\nResults saved to: {out_path}")

    if n_pass < n_total:
        sys.exit(1)


if __name__ == "__main__":
    main()
