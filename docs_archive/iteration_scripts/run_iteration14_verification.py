"""Iteration 14 verification: docstring and per_step_summary consistency.

Checks that class docstrings, compile_to_native_gates(), and
transpile_to_basic_gates() correctly reflect palindromic 2x Stinespring
gate counts introduced in Iteration 13.

Run: cd tutorials && python run_iteration14_verification.py
"""

from __future__ import annotations

import inspect
import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_checks() -> list[dict]:
    checks: list[dict] = []

    def add(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"name": name, "passed": passed, "detail": detail})

    # ================================================================
    # Section 1: Docstring gate count consistency
    # ================================================================
    from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator
    from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator
    from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

    # 1a: qudit docstring should mention 66
    doc_qudit = QuditGKSLCircuitSimulator.__doc__ or ""
    add(
        "qudit_docstring_total_66",
        "66" in doc_qudit,
        f"Found '66' in docstring: {'66' in doc_qudit}",
    )
    add(
        "qudit_docstring_no_40_total",
        "Total: 40" not in doc_qudit,
        "Should not contain 'Total: 40'",
    )

    # 1b: qubit docstring should mention 66
    doc_qubit = QubitGKSLCircuitSimulator.__doc__ or ""
    add(
        "qubit_docstring_total_66",
        "66" in doc_qubit,
        f"Found '66' in docstring: {'66' in doc_qubit}",
    )
    add(
        "qubit_docstring_no_40_total",
        "Total: 40" not in doc_qubit,
        "Should not contain 'Total: 40'",
    )

    # 1c: boson docstring should mention fwd+rev or x 2
    doc_boson = QuditGKSLCircuitBosonSimulator.__doc__ or ""
    add(
        "boson_docstring_palindromic",
        "fwd+rev" in doc_boson or "x 2" in doc_boson or "palindromic" in doc_boson,
        "Should mention palindromic/fwd+rev/x2",
    )
    add(
        "boson_docstring_stinespring_40",
        "40" in doc_boson,
        "Should show 40 single-site Stinespring (20 x 2)",
    )
    add(
        "boson_docstring_stinespring_12",
        "12" in doc_boson,
        "Should show 12 pair Stinespring (6 x 2)",
    )

    # ================================================================
    # Section 2: compile_to_native_gates per_step_summary (qudit)
    # ================================================================
    src_compile = inspect.getsource(QuditGKSLCircuitSimulator.compile_to_native_gates)
    add(
        "compile_native_2x_single",
        "2 * single_native_total" in src_compile,
        "per_step_summary should use 2 * single_native_total",
    )
    add(
        "compile_native_2x_pair",
        "2 * pair_uncompiled_total" in src_compile,
        "per_step_summary should use 2 * pair_uncompiled_total",
    )

    # ================================================================
    # Section 3: transpile_to_basic_gates per_step_summary (qubit)
    # ================================================================
    src_transpile = inspect.getsource(QubitGKSLCircuitSimulator.transpile_to_basic_gates)
    add(
        "transpile_2x_single",
        "* 2" in src_transpile.split("total_single")[1].split("\n")[0],
        "total_single should be multiplied by 2",
    )
    add(
        "transpile_2x_pair",
        "* 2" in src_transpile.split("total_pair")[1].split("\n")[0],
        "total_pair should be multiplied by 2",
    )

    # ================================================================
    # Section 4: simulate() gate counts still correct (regression)
    # ================================================================
    from gksl_physical_parameters import GKSLPhysicalParameters

    params = GKSLPhysicalParameters()
    N = params.N_molecules  # 4
    n_pairs = len(params.neighbors)  # 3

    # 4a: qudit simulate gate counts
    sim_qudit = QuditGKSLCircuitSimulator(params)
    r_qudit = sim_qudit.simulate(t_max=1.0, n_steps=2)
    n_single = 5 * N  # 20
    n_pair = 2 * n_pairs  # 6
    expected_gps = 2 * (N + n_pairs) + 2 * n_single + 2 * n_pair  # 66
    add(
        "qudit_simulate_gates_per_step",
        r_qudit["gates_per_step"] == expected_gps,
        f"Expected {expected_gps}, got {r_qudit['gates_per_step']}",
    )
    add(
        "qudit_simulate_stinespring_single_2x",
        r_qudit["gate_breakdown"]["cu_two_stinespring_single"] == 2 * n_single,
        f"Expected {2 * n_single}, got {r_qudit['gate_breakdown']['cu_two_stinespring_single']}",
    )
    add(
        "qudit_simulate_stinespring_pair_2x",
        r_qudit["gate_breakdown"]["cu_multi_stinespring_pair"] == 2 * n_pair,
        f"Expected {2 * n_pair}, got {r_qudit['gate_breakdown']['cu_multi_stinespring_pair']}",
    )

    # 4b: qubit simulate gate counts
    sim_qubit = QubitGKSLCircuitSimulator(params)
    r_qubit = sim_qubit.simulate(t_max=1.0, n_steps=2)
    add(
        "qubit_simulate_gates_per_step",
        r_qubit["gates_per_step"] == expected_gps,
        f"Expected {expected_gps}, got {r_qubit['gates_per_step']}",
    )

    # ================================================================
    # Section 5: build_full regression
    # ================================================================
    info_qudit = sim_qudit.build_full_trotter_step_circuit(dt=1.0)
    add(
        "qudit_build_full_total_66",
        info_qudit["total_gates"] == 66,
        f"Expected 66, got {info_qudit['total_gates']}",
    )
    add(
        "qudit_build_full_stinespring_52",
        info_qudit["n_stinespring_gates"] == 52,
        f"Expected 52, got {info_qudit['n_stinespring_gates']}",
    )

    try:
        info_qubit = sim_qubit.build_full_trotter_step_circuit(dt=1.0)
        add(
            "qubit_build_full_total_66",
            info_qubit["total_gates"] == 66,
            f"Expected 66, got {info_qubit['total_gates']}",
        )
        add(
            "qubit_build_full_stinespring_52",
            info_qubit["n_stinespring_gates"] == 52,
            f"Expected 52, got {info_qubit['n_stinespring_gates']}",
        )
    except ImportError:
        add("qubit_build_full_total_66", True, "SKIP: qiskit not available")
        add("qubit_build_full_stinespring_52", True, "SKIP: qiskit not available")

    return checks


def main() -> None:
    checks = run_checks()
    n_pass = sum(1 for c in checks if c["passed"])
    n_total = len(checks)

    print(f"\n{'='*60}")
    print(f"Iteration 14 Verification: {n_pass}/{n_total} checks passed")
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
    out_path = os.path.join(out_dir, f"iteration14_docstring_per_step_{ts}.json")
    with open(out_path, "w") as f:
        json.dump(
            {
                "iteration": 14,
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
