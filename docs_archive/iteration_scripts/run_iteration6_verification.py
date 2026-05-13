#!/usr/bin/env python3
"""
Iteration 6 Verification Script

Verifies:
1. qubit_noisy_simulator.py uses density matrix simulation (not Qiskit Aer)
2. Depolarization rates are consistent (0.001 = 0.1%) across all notebooks
3. Noise model documentation is consistent
4. GKSL notebook noise parameters are correct
5. Simulator classes have expected interface
"""

import json
import os
import re
import sys
from datetime import datetime, timezone

GKSL_NB = os.path.join(os.path.dirname(__file__), "quantum_dynamics_gksl_comparison.ipynb")
COMPLETE_NB = os.path.join(os.path.dirname(__file__), "quantum_dynamics_complete_comparison.ipynb")
QUBIT_NOISY = os.path.join(os.path.dirname(__file__), "qubit_noisy_simulator.py")
QUDIT_NOISY = os.path.join(os.path.dirname(__file__), "mqt_qudits_noisy_simulator.py")

checks = []


def check(name, condition, detail):
    status = "PASS" if condition else "FAIL"
    checks.append({"name": name, "status": status, "detail": detail})
    symbol = "✓" if condition else "✗"
    print(f"  {symbol} {name}: {detail}")


def main():
    print("=" * 70)
    print("Iteration 6 Verification")
    print("=" * 70)

    # Load notebooks
    with open(GKSL_NB) as f:
        gksl = json.load(f)
    with open(COMPLETE_NB) as f:
        complete = json.load(f)
    with open(QUBIT_NOISY) as f:
        qubit_noisy_src = f.read()
    with open(QUDIT_NOISY) as f:
        qudit_noisy_src = f.read()

    # ============================================================
    # Section 1: qubit_noisy_simulator.py Architecture
    # ============================================================
    print("\n--- Section 1: qubit_noisy_simulator.py Architecture ---")

    check(
        "no_qiskit_aer_import",
        "AerSimulator" not in qubit_noisy_src and "qiskit_aer" not in qubit_noisy_src,
        "qubit_noisy_simulator should NOT import Qiskit Aer",
    )

    check(
        "no_transpile_import",
        "transpile" not in qubit_noisy_src,
        "qubit_noisy_simulator should NOT use transpile",
    )

    check(
        "density_matrix_method",
        "_apply_2qubit_pair_depolarizing_dm" in qubit_noisy_src,
        "Should have density matrix depolarizing channel method",
    )

    check(
        "pair_unitary_builder",
        "_build_per_pair_unitaries" in qubit_noisy_src,
        "Should have per-pair unitary builder",
    )

    check(
        "uses_build_qubit_unitary",
        "build_H_transfer_qubit_unitary" in qubit_noisy_src
        and "build_H_TTA_qubit_unitary" in qubit_noisy_src,
        "Should import 16x16 qubit unitary builders",
    )

    check(
        "d_equals_4",
        "self.d = 4" in qubit_noisy_src or "self.d=4" in qubit_noisy_src,
        "Local dimension should be d=4 (2 qubits per molecule)",
    )

    check(
        "dim_256",
        "self.d ** self.N" in qubit_noisy_src,
        "System dimension should be 4^N = 256",
    )

    check(
        "pair_level_noise_in_simulate",
        "_apply_2qubit_pair_depolarizing_dm" in qubit_noisy_src
        and "U_tr @ rho @" in qubit_noisy_src,
        "Simulate should use pair-level noise with density matrix evolution",
    )

    check(
        "no_circuit_compose",
        "circuit.compose" not in qubit_noisy_src,
        "Should NOT compose circuits for time evolution",
    )

    # ============================================================
    # Section 2: Depolarization Rate Consistency (0.001)
    # ============================================================
    print("\n--- Section 2: Depolarization Rate Consistency ---")

    # GKSL notebook code cells
    gksl_cell34_src = "".join(gksl["cells"][34]["source"])
    gksl_cell38_src = "".join(gksl["cells"][38]["source"])

    check(
        "gksl_cell34_depol_rate",
        "p_depol=0.001" in gksl_cell34_src,
        "GKSL Cell 34 should use p_depol=0.001",
    )

    check(
        "gksl_cell38_depol_rate",
        "p_depol=0.001" in gksl_cell38_src,
        "GKSL Cell 38 should use p_depol=0.001",
    )

    # Complete notebook code cells
    complete_cell15_src = "".join(complete["cells"][15]["source"])
    complete_cell23_src = "".join(complete["cells"][23]["source"])

    check(
        "complete_cell15_depol_rate",
        "'depol_2q': 0.001" in complete_cell15_src
        or '"depol_2q": 0.001' in complete_cell15_src,
        "Complete Cell 15 should use depol_2q=0.001",
    )

    check(
        "complete_cell23_depol_rate",
        "'depol_2q': 0.001" in complete_cell23_src
        or '"depol_2q": 0.001' in complete_cell23_src,
        "Complete Cell 23 should use depol_2q=0.001",
    )

    # ============================================================
    # Section 3: Documentation Consistency
    # ============================================================
    print("\n--- Section 3: Documentation Consistency ---")

    # GKSL markdown cells should say 0.001
    gksl_cell33_src = "".join(gksl["cells"][33]["source"])
    gksl_cell37_src = "".join(gksl["cells"][37]["source"])
    gksl_cell39_src = "".join(gksl["cells"][39]["source"])
    gksl_cell43_src = "".join(gksl["cells"][43]["source"])

    check(
        "gksl_cell33_doc",
        "0.001" in gksl_cell33_src and "0.01" not in gksl_cell33_src.replace("0.001", ""),
        "GKSL Cell 33 markdown should reference p_depol=0.001",
    )

    check(
        "gksl_cell37_doc",
        "0.001" in gksl_cell37_src and "0.01" not in gksl_cell37_src.replace("0.001", ""),
        "GKSL Cell 37 markdown should reference p_depol=0.001",
    )

    check(
        "gksl_cell39_doc",
        "0.001" in gksl_cell39_src,
        "GKSL Cell 39 should reference p_depol=0.001",
    )

    check(
        "gksl_cell43_doc",
        "0.001" in gksl_cell43_src,
        "GKSL Cell 43 should reference p_depol=0.001",
    )

    # Complete notebook markdown cells
    complete_cell14_src = "".join(complete["cells"][14]["source"])
    complete_cell22_src = "".join(complete["cells"][22]["source"])

    check(
        "complete_cell14_doc",
        "0.1%" in complete_cell14_src,
        "Complete Cell 14 should reference 0.1% noise rate",
    )

    check(
        "complete_cell22_doc",
        "0.1%" in complete_cell22_src,
        "Complete Cell 22 should reference 0.1% noise rate",
    )

    # Code comments should say 0.1%
    check(
        "complete_cell15_comment",
        "0.1%" in complete_cell15_src,
        "Complete Cell 15 comment should say 0.1%",
    )

    check(
        "complete_cell23_comment",
        "0.1%" in complete_cell23_src,
        "Complete Cell 23 comment should say 0.1%",
    )

    # ============================================================
    # Section 4: GKSL Notebook Noise Parameters
    # ============================================================
    print("\n--- Section 4: GKSL Noise Parameters ---")

    check(
        "gksl_cell34_p_dephasing_0",
        "p_dephasing=0.0" in gksl_cell34_src,
        "GKSL Cell 34 should have p_dephasing=0.0",
    )

    check(
        "gksl_cell34_depol_pair_only",
        "depol_pair_only=True" in gksl_cell34_src,
        "GKSL Cell 34 should have depol_pair_only=True",
    )

    check(
        "gksl_cell38_p_dephasing_0",
        "p_dephasing=0.0" in gksl_cell38_src,
        "GKSL Cell 38 should have p_dephasing=0.0",
    )

    check(
        "gksl_cell38_depol_pair_only",
        "depol_pair_only=True" in gksl_cell38_src,
        "GKSL Cell 38 should have depol_pair_only=True",
    )

    # ============================================================
    # Section 5: Noise Granularity Match
    # ============================================================
    print("\n--- Section 5: Noise Granularity Match ---")

    check(
        "qubit_no_gate_level_noise",
        "add_all_qubit_quantum_error" not in qubit_noisy_src,
        "Qubit simulator should NOT use gate-level noise model",
    )

    check(
        "qubit_partial_trace",
        "einsum" in qubit_noisy_src and ("Tr_{ij}" in qubit_noisy_src or "partial trace" in qubit_noisy_src.lower()),
        "Qubit simulator should use partial trace for depolarizing channel",
    )

    check(
        "qudit_partial_trace",
        "einsum" in qudit_noisy_src,
        "Qudit simulator should use partial trace for depolarizing channel",
    )

    # Both should have 12 noise events per step (3 pairs x 2 half-steps x 2 types)
    check(
        "qubit_pair_noise_structure",
        "U_transfer_half_list" in qubit_noisy_src and "U_TTA_half_list" in qubit_noisy_src,
        "Qubit should have separate transfer and TTA unitary lists",
    )

    # ============================================================
    # Summary
    # ============================================================
    passed = sum(1 for c in checks if c["status"] == "PASS")
    failed = sum(1 for c in checks if c["status"] == "FAIL")

    print("\n" + "=" * 70)
    print(f"Results: {passed}/{len(checks)} passed, {failed} failed")
    if failed == 0:
        print("✓ All checks passed!")
    else:
        print("✗ Some checks failed:")
        for c in checks:
            if c["status"] == "FAIL":
                print(f"  FAIL: {c['name']}: {c['detail']}")
    print("=" * 70)

    # Save results
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    results_dir = os.path.join(
        os.path.dirname(__file__), "..", "developing", "verification_results"
    )
    os.makedirs(results_dir, exist_ok=True)
    results_file = os.path.join(
        results_dir, f"iteration6_qubit_noise_arch_{timestamp}.json"
    )

    with open(results_file, "w") as f:
        json.dump(
            {
                "timestamp": timestamp,
                "iteration": 6,
                "checks": checks,
                "passed": passed,
                "failed": failed,
            },
            f,
            indent=2,
        )

    print(f"\nResults saved to: {results_file}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
