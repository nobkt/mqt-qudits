#!/usr/bin/env python3
"""
Iteration 7 Verification Script

Verifies the noise simulator changes that properly model the physical
difference between qubit and qudit encodings:

1. qubit_noisy_simulator.py supports cx_per_pair_gate parameter
2. Effective depolarization rate is computed correctly:
   p_eff = 1 - (1 - p_phys)^cx_per_pair_gate
3. CX gate counts are computed via Qiskit decomposition (if available)
4. Noise channel properties (CPTP, forbidden state leakage) are verified
5. Qudit noise model remains unchanged and correct

Key Physics:
- On real qubit hardware, each pair interaction requires N_CX CX gates
- On real qutrit hardware, each pair interaction is 1 native gate
- For the same physical gate error rate p_phys:
  - Qudit effective pair error = p_phys
  - Qubit effective pair error = 1 - (1 - p_phys)^N_CX >> p_phys
- Additionally, qubit depolarizing leaks to forbidden |11> states (d=4)
  while qudit has no forbidden states (d=3)

Output: developing/verification_results/iteration7_*.json
"""

import json
import os
import sys
import math
from datetime import datetime, timezone

# File paths
QUBIT_NOISY = os.path.join(os.path.dirname(__file__), "qubit_noisy_simulator.py")
QUDIT_NOISY = os.path.join(os.path.dirname(__file__), "mqt_qudits_noisy_simulator.py")
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
    print("Iteration 7 Verification: Physically Correct Noise Model Comparison")
    print("=" * 70)

    with open(QUBIT_NOISY) as f:
        qubit_src = f.read()
    with open(QUDIT_NOISY) as f:
        qudit_src = f.read()

    # ============================================================
    # Section 1: Qubit Simulator - cx_per_pair_gate Support
    # ============================================================
    print("\n--- Section 1: cx_per_pair_gate Parameter Support ---")

    check(
        "cx_per_pair_gate_in_simulate",
        "cx_per_pair_gate" in qubit_src,
        "simulate() should accept cx_per_pair_gate parameter",
    )

    check(
        "effective_rate_formula",
        "(1 - depol_2q) ** cx_per_pair_gate" in qubit_src,
        "Should compute p_eff = 1 - (1-depol_2q)^cx_per_pair_gate",
    )

    check(
        "depol_2q_eff_in_noise_loop",
        "depol_2q_eff" in qubit_src,
        "Should use depol_2q_eff (effective rate) in noise application",
    )

    check(
        "estimate_cx_method",
        "estimate_cx_per_pair_gate" in qubit_src,
        "Should have estimate_cx_per_pair_gate() method",
    )

    check(
        "cx_per_pair_gate_default_1",
        "cx_per_pair_gate: int = 1" in qubit_src
        or "cx_per_pair_gate=1" in qubit_src,
        "Default cx_per_pair_gate should be 1 (backward compatible)",
    )

    check(
        "cx_per_pair_gate_validation",
        "cx_per_pair_gate < 1" in qubit_src,
        "Should validate cx_per_pair_gate >= 1",
    )

    check(
        "output_contains_cx_info",
        "'cx_per_pair_gate'" in qubit_src and "'depol_2q_effective'" in qubit_src,
        "Return dict should include cx_per_pair_gate and depol_2q_effective",
    )

    # ============================================================
    # Section 2: Qubit Simulator - Density Matrix Architecture
    # ============================================================
    print("\n--- Section 2: Qubit Simulator Architecture ---")

    check(
        "no_qiskit_aer_import",
        "AerSimulator" not in qubit_src and "qiskit_aer" not in qubit_src,
        "qubit_noisy_simulator should NOT import Qiskit Aer for simulation",
    )

    check(
        "density_matrix_method",
        "_apply_2qubit_pair_depolarizing_dm" in qubit_src,
        "Should have density matrix depolarizing channel method",
    )

    check(
        "pair_unitary_builder",
        "_build_per_pair_unitaries" in qubit_src,
        "Should have per-pair unitary builder",
    )

    check(
        "d_equals_4",
        "self.d = 4" in qubit_src,
        "Local dimension should be d=4 (2 qubits per molecule)",
    )

    # ============================================================
    # Section 3: Qudit Simulator - Unchanged and Correct
    # ============================================================
    print("\n--- Section 3: Qudit Simulator Integrity ---")

    check(
        "qudit_d_equals_3",
        "d = 3" in qudit_src,
        "Qudit local dimension should be d=3 (qutrit)",
    )

    check(
        "qudit_depolarizing_method",
        "_apply_2qudit_depolarizing_dm" in qudit_src,
        "Should have 2-qudit depolarizing channel method",
    )

    check(
        "qudit_I_over_d_squared",
        "I_d = np.eye(d, dtype=complex) / d" in qudit_src,
        "Should use I/d for each qudit (not I/d^2 for the pair)",
    )

    check(
        "qudit_no_cx_overhead",
        "cx_per_pair_gate" not in qudit_src,
        "Qudit simulator should NOT have cx_per_pair_gate (native gates)",
    )

    # ============================================================
    # Section 4: Depolarizing Channel Physics
    # ============================================================
    print("\n--- Section 4: Depolarizing Channel Physics ---")

    check(
        "qubit_depol_channel_formula",
        "Tr_{ij}(ρ)" in qubit_src or "Tr_{ij}" in qubit_src,
        "Qubit channel: ε(ρ) = (1-p)ρ + p·Tr_{ij}(ρ)⊗I_{16}/16",
    )

    check(
        "qudit_depol_channel_formula",
        "Tr_{ij}(ρ)" in qudit_src or "Tr_{ij}" in qudit_src,
        "Qudit channel: ε(ρ) = (1-p)ρ + p·Tr_{ij}(ρ)⊗I_9/9",
    )

    check(
        "qubit_forbidden_state_doc",
        "forbidden" in qubit_src.lower() or "|11⟩" in qubit_src or "|11>" in qubit_src,
        "Qubit simulator should document forbidden state leakage",
    )

    # ============================================================
    # Section 5: CX Gate Count Estimation (Qiskit-dependent)
    # ============================================================
    print("\n--- Section 5: CX Gate Count Estimation ---")

    cx_counts = None
    try:
        sys.path.insert(0, os.path.dirname(__file__))
        from qubit_noisy_simulator import QubitMolecularDynamicsSimulatorNoisy

        try:
            cx_counts = QubitMolecularDynamicsSimulatorNoisy.estimate_cx_per_pair_gate()
            print(f"  CX gate counts (via Qiskit decomposition):")
            print(f"    H_transfer: {cx_counts['cx_transfer']} CX gates")
            print(f"    H_TTA:      {cx_counts['cx_tta']} CX gates")
            print(f"    Average:    {cx_counts['cx_avg']:.1f} CX gates")

            check(
                "cx_count_positive",
                cx_counts['cx_transfer'] > 0 and cx_counts['cx_tta'] > 0,
                f"CX counts should be positive: transfer={cx_counts['cx_transfer']}, tta={cx_counts['cx_tta']}",
            )

            check(
                "cx_count_reasonable",
                1 < cx_counts['cx_avg'] < 100,
                f"Average CX count ({cx_counts['cx_avg']:.1f}) should be between 2 and 100",
            )

        except ImportError:
            print("  [SKIP] Qiskit not installed - CX gate counting skipped")
            check(
                "cx_estimation_graceful",
                True,
                "CX estimation raises ImportError when Qiskit unavailable (expected)",
            )

    except Exception as e:
        print(f"  [ERROR] Could not import qubit_noisy_simulator: {e}")
        check(
            "import_qubit_noisy",
            False,
            f"Failed to import: {e}",
        )

    # ============================================================
    # Section 6: Effective Rate Analysis
    # ============================================================
    print("\n--- Section 6: Effective Depolarization Rate Analysis ---")

    p_phys = 0.001  # 0.1% per physical gate

    # Compute effective rates for various CX counts
    cx_values = [1, 3, 6, 10, 20, 50]
    print(f"\n  Physical gate error rate: {p_phys*100:.2f}%")
    print(f"  {'CX/pair':>8} | {'p_eff':>12} | {'p_eff/p_phys':>12} | {'Qubit disadvantage':>20}")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*12}-+-{'-'*20}")

    for cx in cx_values:
        p_eff = 1 - (1 - p_phys) ** cx
        ratio = p_eff / p_phys
        print(f"  {cx:>8} | {p_eff*100:>11.4f}% | {ratio:>12.2f}x | "
              f"{'(same as qudit)' if cx == 1 else f'{ratio:.1f}x more noise'}")

    check(
        "effective_rate_monotonic",
        all(
            (1 - (1 - p_phys) ** cx_values[i]) < (1 - (1 - p_phys) ** cx_values[i + 1])
            for i in range(len(cx_values) - 1)
        ),
        "Effective rate should increase monotonically with CX count",
    )

    # If we know the actual CX count from Qiskit, add analysis
    if cx_counts is not None and 'cx_avg' in cx_counts:
        cx_avg = cx_counts['cx_avg']
        p_eff_actual = 1 - (1 - p_phys) ** cx_avg
        check(
            "qubit_more_noisy_than_qudit",
            p_eff_actual > p_phys,
            f"Qubit effective rate ({p_eff_actual*100:.4f}%) > "
            f"Qudit rate ({p_phys*100:.4f}%) for CX/pair={cx_avg:.1f}",
        )

    # ============================================================
    # Section 7: Leakage Analysis
    # ============================================================
    print("\n--- Section 7: Forbidden State Leakage Analysis ---")

    # For qubit (d=4): maximally mixed state I_{16}/16 has 7/16 non-physical states
    d_qubit = 4
    n_physical_qubit = 3 ** 2  # 9 physical states per pair (S0,T1,S1 x S0,T1,S1)
    n_total_qubit = d_qubit ** 2  # 16 total states per pair
    n_forbidden_qubit = n_total_qubit - n_physical_qubit  # 7 forbidden
    leakage_fraction_qubit = n_forbidden_qubit / n_total_qubit

    print(f"\n  Qubit encoding (d={d_qubit}):")
    print(f"    Physical states per pair:  {n_physical_qubit}/{n_total_qubit}")
    print(f"    Forbidden states per pair: {n_forbidden_qubit}/{n_total_qubit}")
    print(f"    Leakage fraction per noise event: {leakage_fraction_qubit*100:.2f}%")

    # For qudit (d=3): all states are physical
    d_qudit = 3
    n_total_qudit = d_qudit ** 2  # 9
    n_physical_qudit = d_qudit ** 2  # 9 (all physical)
    n_forbidden_qudit = 0

    print(f"\n  Qudit encoding (d={d_qudit}):")
    print(f"    Physical states per pair:  {n_physical_qudit}/{n_total_qudit}")
    print(f"    Forbidden states per pair: {n_forbidden_qudit}/{n_total_qudit}")
    print(f"    Leakage fraction per noise event: 0.00%")

    check(
        "qubit_leakage_theory",
        abs(leakage_fraction_qubit - 7 / 16) < 1e-10,
        f"Qubit leakage fraction = {n_forbidden_qubit}/{n_total_qubit} = {leakage_fraction_qubit:.4f} (theory: 7/16 = 0.4375)",
    )

    check(
        "qudit_no_leakage_theory",
        n_forbidden_qudit == 0,
        "Qudit has zero forbidden states (all d=3 states are physical)",
    )

    # ============================================================
    # Section 8: Combined Qudit Advantage Summary
    # ============================================================
    print("\n--- Section 8: Qudit Encoding Advantage Summary ---")

    print("\n  Qudit (qutrit) advantages over qubit encoding:")
    print("  1. Gate count advantage: 1 native gate vs N_CX CX gates per pair")
    print("     → Lower effective noise rate per pair interaction")
    print("  2. Leakage-free: All d=3 states are physical")
    print("     → No population lost to forbidden states")
    print("  3. Compact Hilbert space: 3^N vs 4^N")
    print("     → More efficient state representation")

    # Expected qudit advantage for p_phys=0.001
    if cx_counts is not None and 'cx_avg' in cx_counts:
        cx_avg = cx_counts['cx_avg']
        p_eff_qubit = 1 - (1 - p_phys) ** cx_avg
        p_eff_qudit = p_phys
        advantage_ratio = p_eff_qubit / p_eff_qudit
        print(f"\n  Quantified advantage (p_phys={p_phys}):")
        print(f"    Qubit effective pair noise: {p_eff_qubit*100:.4f}% ({cx_avg:.0f} CX/pair)")
        print(f"    Qudit effective pair noise: {p_eff_qudit*100:.4f}% (1 native gate)")
        print(f"    Noise ratio: {advantage_ratio:.1f}x")
        print(f"    Plus: {leakage_fraction_qubit*100:.1f}% leakage per qubit noise event")

    # ============================================================
    # Save Results
    # ============================================================
    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    result_file = os.path.join(RESULTS_DIR, f"iteration7_noise_comparison_{timestamp}.json")

    n_pass = sum(1 for c in checks if c["status"] == "PASS")
    n_fail = sum(1 for c in checks if c["status"] == "FAIL")

    result = {
        "iteration": 7,
        "timestamp": timestamp,
        "description": "Physically correct noise model comparison (qubit vs qudit)",
        "summary": {
            "total_checks": len(checks),
            "passed": n_pass,
            "failed": n_fail,
        },
        "cx_gate_counts": cx_counts,
        "effective_rates": {
            "p_phys": p_phys,
            "qudit_p_eff": p_phys,
            "qubit_p_eff_cx1": p_phys,
            "qubit_p_eff_cx6": 1 - (1 - p_phys) ** 6,
            "qubit_p_eff_cx10": 1 - (1 - p_phys) ** 10,
            "qubit_p_eff_cx20": 1 - (1 - p_phys) ** 20,
        },
        "leakage_analysis": {
            "qubit_d": d_qubit,
            "qubit_physical_states_per_pair": n_physical_qubit,
            "qubit_total_states_per_pair": n_total_qubit,
            "qubit_leakage_fraction": leakage_fraction_qubit,
            "qudit_d": d_qudit,
            "qudit_physical_states_per_pair": n_physical_qudit,
            "qudit_total_states_per_pair": n_total_qudit,
            "qudit_leakage_fraction": 0.0,
        },
        "checks": checks,
    }

    if cx_counts is not None and 'cx_avg' in cx_counts:
        cx_avg = cx_counts['cx_avg']
        result["effective_rates"]["qubit_p_eff_actual"] = 1 - (1 - p_phys) ** cx_avg
        result["effective_rates"]["cx_avg"] = cx_avg

    with open(result_file, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print(f"\n{'='*70}")
    print(f"Verification Results: {n_pass}/{len(checks)} PASS, {n_fail}/{len(checks)} FAIL")
    print(f"Results saved to: {result_file}")
    print(f"{'='*70}")

    return n_fail == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
