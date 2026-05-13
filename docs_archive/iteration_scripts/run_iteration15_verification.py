"""Iteration 15 verification: estimated_gates_per_step consistency.

Checks that non-circuit simulators' estimated_gates_per_step correctly
reflects the palindromic 2nd-order Trotter step (2× for Hamiltonian half-steps),
and matches the circuit simulators' exact gate counts where applicable.

Bugs fixed in iteration 15:
  Bug J: qudit_gksl_simulator.py - estimated_gates_per_step was 59 (should be 66)
  Bug K: qubit_gksl_simulator.py - estimated_gates_per_step was 350 (should be 388)

Run: cd tutorials && python run_iteration15_verification.py
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

    # ================================================================
    # Section 1: Non-circuit qudit simulator gate count (Bug J fix)
    # ================================================================
    from qudit_gksl_simulator import QuditGKSLSimulator

    sim_qudit = QuditGKSLSimulator(params)
    r_qudit = sim_qudit.simulate(t_max=1.0, n_steps=2)

    n_lindblad = 2 * n_pairs + 5 * N  # 6 + 20 = 26
    # Palindromic: 2 × (N + n_pairs) Hamiltonian + 2 × n_lindblad Stinespring
    expected_qudit_gps = 2 * (N + n_pairs) + n_lindblad * 2  # 14 + 52 = 66
    actual_qudit_gps = r_qudit["estimated_gates_per_step"]

    add(
        "qudit_noncircuit_gps_value",
        actual_qudit_gps == expected_qudit_gps,
        f"Expected {expected_qudit_gps}, got {actual_qudit_gps}",
    )

    # Check it's NOT the old broken value 59
    old_broken = N + n_pairs + n_lindblad * 2  # 4 + 3 + 52 = 59
    add(
        "qudit_noncircuit_not_old_59",
        actual_qudit_gps != old_broken,
        f"Should not be old value {old_broken}",
    )

    # ================================================================
    # Section 2: Non-circuit qubit simulator gate count (Bug K fix)
    # ================================================================
    from qubit_gksl_simulator import QubitGKSLSimulator

    sim_qubit = QubitGKSLSimulator(params)
    r_qubit = sim_qubit.simulate(t_max=1.0, n_steps=2)

    n_sys_qubits = 2 * N  # 8
    n_ancilla = n_lindblad  # 26
    # Palindromic: 2 × Hamiltonian_basic + 2 × Stinespring_basic
    expected_qubit_gps = 2 * (n_sys_qubits + n_pairs * 10) + n_ancilla * 6 * 2  # 76+312=388
    actual_qubit_gps = r_qubit["estimated_gates_per_step"]

    add(
        "qubit_noncircuit_gps_value",
        actual_qubit_gps == expected_qubit_gps,
        f"Expected {expected_qubit_gps}, got {actual_qubit_gps}",
    )

    # Check it's NOT the old broken value 350
    old_broken_qb = n_sys_qubits + n_pairs * 10 + n_ancilla * 6 * 2  # 8+30+312=350
    add(
        "qubit_noncircuit_not_old_350",
        actual_qubit_gps != old_broken_qb,
        f"Should not be old value {old_broken_qb}",
    )

    # ================================================================
    # Section 3: Cross-check with circuit simulators
    # ================================================================
    from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

    sim_circ_qudit = QuditGKSLCircuitSimulator(params)
    r_circ_qudit = sim_circ_qudit.simulate(t_max=1.0, n_steps=2)

    # Qudit non-circuit estimated should match circuit exact
    add(
        "qudit_noncircuit_matches_circuit",
        actual_qudit_gps == r_circ_qudit["gates_per_step"],
        f"Non-circuit {actual_qudit_gps} vs circuit {r_circ_qudit['gates_per_step']}",
    )

    from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

    sim_circ_qubit = QubitGKSLCircuitSimulator(params)
    r_circ_qubit = sim_circ_qubit.simulate(t_max=1.0, n_steps=2)

    # Qubit circuit should still be 66 (high-level gates)
    add(
        "qubit_circuit_gps_66",
        r_circ_qubit["gates_per_step"] == 66,
        f"Expected 66, got {r_circ_qubit['gates_per_step']}",
    )

    # Qudit circuit should still be 66
    add(
        "qudit_circuit_gps_66",
        r_circ_qudit["gates_per_step"] == 66,
        f"Expected 66, got {r_circ_qudit['gates_per_step']}",
    )

    # ================================================================
    # Section 4: total_estimated_gates consistency
    # ================================================================
    n_steps = 2
    add(
        "qudit_total_estimated_gates",
        r_qudit["total_estimated_gates"] == expected_qudit_gps * n_steps,
        f"Expected {expected_qudit_gps * n_steps}, got {r_qudit['total_estimated_gates']}",
    )
    add(
        "qubit_total_estimated_gates",
        r_qubit["total_estimated_gates"] == expected_qubit_gps * n_steps,
        f"Expected {expected_qubit_gps * n_steps}, got {r_qubit['total_estimated_gates']}",
    )

    # ================================================================
    # Section 5: Simulation physics regression (trace preservation)
    # ================================================================
    # Verify non-circuit qudit trace is preserved
    qudit_trace_final = r_qudit["trace"][-1]
    add(
        "qudit_noncircuit_trace_preserved",
        abs(qudit_trace_final - 1.0) < 1e-10,
        f"Final trace = {qudit_trace_final}",
    )

    # Verify non-circuit qubit trace is preserved
    qubit_trace_final = r_qubit["trace"][-1]
    add(
        "qubit_noncircuit_trace_preserved",
        abs(qubit_trace_final - 1.0) < 1e-10,
        f"Final trace = {qubit_trace_final}",
    )

    # ================================================================
    # Section 6: Iteration 14 regression checks
    # ================================================================
    # Docstrings still correct
    doc_qudit_circ = QuditGKSLCircuitSimulator.__doc__ or ""
    add(
        "regression_qudit_docstring_66",
        "66" in doc_qudit_circ,
        "Circuit qudit docstring should contain '66'",
    )

    doc_qubit_circ = QubitGKSLCircuitSimulator.__doc__ or ""
    add(
        "regression_qubit_docstring_66",
        "66" in doc_qubit_circ,
        "Circuit qubit docstring should contain '66'",
    )

    # Circuit simulator gate breakdowns still correct
    add(
        "regression_qudit_circuit_breakdown_single",
        r_circ_qudit["gate_breakdown"]["cu_two_stinespring_single"] == 40,
        f"Expected 40, got {r_circ_qudit['gate_breakdown']['cu_two_stinespring_single']}",
    )
    add(
        "regression_qudit_circuit_breakdown_pair",
        r_circ_qudit["gate_breakdown"]["cu_multi_stinespring_pair"] == 12,
        f"Expected 12, got {r_circ_qudit['gate_breakdown']['cu_multi_stinespring_pair']}",
    )

    return checks


def main() -> None:
    checks = run_checks()
    n_pass = sum(1 for c in checks if c["passed"])
    n_total = len(checks)

    print(f"\n{'='*60}")
    print(f"Iteration 15 Verification: {n_pass}/{n_total} checks passed")
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
    out_path = os.path.join(out_dir, f"iteration15_gate_count_consistency_{ts}.json")
    with open(out_path, "w") as f:
        json.dump(
            {
                "iteration": 15,
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
