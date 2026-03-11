"""Iteration 13 verification: palindromic Trotter step in circuit simulators.

Checks:
  Section 1: Qudit circuit simulator palindromic Trotter fix
  Section 2: Qubit circuit simulator palindromic Trotter fix
  Section 3: Gate count corrections (all circuit simulators)
  Section 4: build_full_trotter_step_circuit palindromic structure
  Section 5: Boson circuit gate count fix
  Section 6: Iteration 12 regression (d_anc checks)
  Section 7: N=2 exact match (Hamiltonian terms commute)

Usage:
  cd tutorials
  python run_iteration13_verification.py
"""

from __future__ import annotations

import json
import os
import re
import sys
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

results: list[dict] = []


def check(name: str, passed: bool, detail: str) -> None:
    status = "PASS" if passed else "FAIL"
    results.append({"name": name, "status": status, "detail": detail})
    mark = "✓" if passed else "✗"
    print(f"  [{mark}] {name}: {detail}")


# ---------------------------------------------------------------------------
# Section 1: Qudit circuit simulator palindromic Trotter
# ---------------------------------------------------------------------------
print("\n=== Section 1: Qudit circuit simulator palindromic Trotter ===")

src = open("qudit_gksl_circuit_simulator.py").read()

# Check _trotter_step_circuit uses dt/2 (not full dt) for Stinespring
check(
    "qudit_circuit_dt_half",
    "self._build_local_stinespring_unitary(L_local, dt / 2)" in src,
    "Uses dt/2 for Stinespring unitary (not full dt)",
)

# Check palindromic structure (reversed loop)
check(
    "qudit_circuit_palindromic",
    "for op_type, sites, kraus_ops in reversed(kraus_list):" in src,
    "Has reversed() loop for palindromic Lindblad ordering",
)

# Check docstring mentions palindromic
trotter_method_match = re.search(
    r'def _trotter_step_circuit.*?"""(.*?)"""', src, re.DOTALL
)
if trotter_method_match:
    docstring = trotter_method_match.group(1)
    check(
        "qudit_circuit_docstring",
        "palindromic" in docstring.lower(),
        "Docstring mentions palindromic",
    )
else:
    check("qudit_circuit_docstring", False, "Could not find docstring")

# Check precompute pattern
check(
    "qudit_circuit_precompute",
    "kraus_list.append((op_type, sites, kraus_ops))" in src,
    "Precomputes Kraus operators before forward/reverse passes",
)

# ---------------------------------------------------------------------------
# Section 2: Qubit circuit simulator palindromic Trotter
# ---------------------------------------------------------------------------
print("\n=== Section 2: Qubit circuit simulator palindromic Trotter ===")

src_qubit = open("qubit_gksl_circuit_simulator.py").read()

check(
    "qubit_circuit_dt_half",
    "self._build_local_stinespring_unitary(L_local, dt / 2)" in src_qubit,
    "Uses dt/2 for Stinespring unitary",
)

check(
    "qubit_circuit_palindromic",
    "for op_type, sites, K0, K1 in reversed(kraus_list):" in src_qubit,
    "Has reversed() loop for palindromic Lindblad ordering",
)

trotter_qubit_match = re.search(
    r'def _trotter_step_circuit.*?"""(.*?)"""', src_qubit, re.DOTALL
)
if trotter_qubit_match:
    docstring_qubit = trotter_qubit_match.group(1)
    check(
        "qubit_circuit_docstring",
        "palindromic" in docstring_qubit.lower(),
        "Docstring mentions palindromic",
    )
else:
    check("qubit_circuit_docstring", False, "Could not find docstring")

check(
    "qubit_circuit_precompute",
    "kraus_list.append((op_type, sites, K0, K1))" in src_qubit,
    "Precomputes Kraus operators before forward/reverse passes",
)

# ---------------------------------------------------------------------------
# Section 3: Gate count corrections
# ---------------------------------------------------------------------------
print("\n=== Section 3: Gate count corrections ===")

# Qudit circuit: Lindblad gates should be 2x
check(
    "qudit_circuit_gate_single_2x",
    '"cu_two_stinespring_single": 2 * n_stinespring_single' in src,
    "Qudit circuit: Stinespring single gate count ×2",
)

check(
    "qudit_circuit_gate_pair_2x",
    '"cu_multi_stinespring_pair": 2 * n_stinespring_pair' in src,
    "Qudit circuit: Stinespring pair gate count ×2",
)

# Qubit circuit: same
check(
    "qubit_circuit_gate_single_2x",
    '"unitary_8x8_stinespring_single": 2 * n_stinespring_single' in src_qubit,
    "Qubit circuit: Stinespring single gate count ×2",
)

check(
    "qubit_circuit_gate_pair_2x",
    '"unitary_32x32_stinespring_pair": 2 * n_stinespring_pair' in src_qubit,
    "Qubit circuit: Stinespring pair gate count ×2",
)

# gates_per_step formula
check(
    "qudit_circuit_gates_per_step",
    "2 * n_stinespring_single" in src and "2 * n_stinespring_pair" in src,
    "Qudit circuit: gates_per_step uses 2× Lindblad counts",
)

check(
    "qubit_circuit_gates_per_step",
    "2 * n_stinespring_single" in src_qubit and "2 * n_stinespring_pair" in src_qubit,
    "Qubit circuit: gates_per_step uses 2× Lindblad counts",
)

# ---------------------------------------------------------------------------
# Section 4: build_full_trotter_step_circuit palindromic
# ---------------------------------------------------------------------------
print("\n=== Section 4: build_full_trotter_step_circuit palindromic ===")

# Qudit: check for dt/2 in Stinespring circuit construction
check(
    "qudit_build_full_palindromic",
    'L_local, dt / 2, sites[0]' in src
    and '"n_stinespring_gates": 2 * len(self.lindblad_local_info)' in src,
    "Qudit: build_full uses dt/2 and 2× stinespring count",
)

check(
    "qudit_build_full_reverse",
    "stinespring_rev_single_" in src or "stinespring_rev_pair_" in src,
    "Qudit: build_full has reverse pass circuits",
)

check(
    "qubit_build_full_palindromic",
    'L_local, dt / 2, sites[0]' in src_qubit
    and '"n_stinespring_gates": 2 * len(self.lindblad_local_info)' in src_qubit,
    "Qubit: build_full uses dt/2 and 2× stinespring count",
)

check(
    "qubit_build_full_reverse",
    "stinespring_rev_single_" in src_qubit or "stinespring_rev_pair_" in src_qubit,
    "Qubit: build_full has reverse pass circuits",
)

# ---------------------------------------------------------------------------
# Section 5: Boson circuit gate count fix
# ---------------------------------------------------------------------------
print("\n=== Section 5: Boson circuit gate count fix ===")

src_boson = open("qudit_gksl_circuit_boson_simulator.py").read()

check(
    "boson_gate_single_2x",
    '"cu_two_stinespring_single": 2 * n_stinespring_single' in src_boson,
    "Boson circuit: Stinespring single gate count ×2",
)

check(
    "boson_gate_pair_2x",
    '"cu_multi_stinespring_pair": 2 * n_stinespring_pair' in src_boson,
    "Boson circuit: Stinespring pair gate count ×2",
)

check(
    "boson_gates_per_step_formula",
    "2 * n_stinespring_single + 2 * n_stinespring_pair" in src_boson,
    "Boson circuit: gates_per_step uses 2× Lindblad counts",
)

# Boson circuit _trotter_step should still be palindromic (unchanged from iter 12)
check(
    "boson_trotter_palindromic",
    "for op_type, sites, kraus_ops in reversed(kraus_list):" in src_boson,
    "Boson circuit: _trotter_step still palindromic",
)

# ---------------------------------------------------------------------------
# Section 6: Iteration 12 regression (d_anc checks)
# ---------------------------------------------------------------------------
print("\n=== Section 6: Iteration 12 regression ===")

check(
    "boson_d_anc_single",
    "circuit = QuantumCircuit(2, [d, d_anc], 0)" in src_boson,
    "Boson build_stinespring_circuit_single uses [d, d_anc]",
)

check(
    "boson_d_anc_pair",
    "circuit = QuantumCircuit(3, [d, d, d_anc], 0)" in src_boson,
    "Boson build_stinespring_circuit_pair uses [d, d, d_anc]",
)

check(
    "boson_d_anc_from_params",
    "self.d_anc = params.d" in src_boson,
    "Boson: d_anc set from params.d",
)

# ---------------------------------------------------------------------------
# Section 7: Numerical verification (N=2 exact match)
# ---------------------------------------------------------------------------
print("\n=== Section 7: Numerical verification ===")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_simulator import QuditGKSLSimulator
from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator
from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

# N=2: Hamiltonian terms commute, so circuit matches matrix exactly
params_n2 = GKSLPhysicalParameters(N_molecules=2)
r_m2 = QuditGKSLSimulator(params_n2).simulate(t_max=5.0, n_steps=5)
r_c2 = QuditGKSLCircuitSimulator(params_n2).simulate(t_max=5.0, n_steps=5)
max_diff_n2 = max(
    abs(r_m2["populations"][-1][k] - r_c2["populations"][-1][k])
    for k in ["N_S0", "N_T1", "N_S1"]
)
check(
    "qudit_N2_exact_match",
    max_diff_n2 < 1e-14,
    f"Qudit circuit N=2 matches matrix: max_diff={max_diff_n2:.2e}",
)

# N=4: Hamiltonian Trotter error remains (O(dt²))
params_n4 = GKSLPhysicalParameters()
r_m4 = QuditGKSLSimulator(params_n4).simulate(t_max=5.0, n_steps=5)
r_c4 = QuditGKSLCircuitSimulator(params_n4).simulate(t_max=5.0, n_steps=5)
max_diff_n4 = max(
    abs(r_m4["populations"][-1][k] - r_c4["populations"][-1][k])
    for k in ["N_S0", "N_T1", "N_S1"]
)
check(
    "qudit_N4_trotter_error",
    max_diff_n4 < 1e-4,
    f"Qudit circuit N=4 Hamiltonian Trotter error: max_diff={max_diff_n4:.2e}",
)

# Qubit N=2 exact match
r_mq2 = QubitGKSLSimulator(params_n2).simulate(t_max=5.0, n_steps=5)
r_cq2 = QubitGKSLCircuitSimulator(params_n2).simulate(t_max=5.0, n_steps=5)
max_diff_qubit_n2 = max(
    abs(r_mq2["populations"][-1][k] - r_cq2["populations"][-1][k])
    for k in ["N_S0", "N_T1", "N_S1"]
)
check(
    "qubit_N2_exact_match",
    max_diff_qubit_n2 < 1e-14,
    f"Qubit circuit N=2 matches matrix: max_diff={max_diff_qubit_n2:.2e}",
)

# Boson circuit still matches boson matrix
params_boson = GKSLPhysicalParameters(
    N_molecules=2, with_boson=True, n_max=1, g_eph=0.005
)
r_cb = QuditGKSLCircuitBosonSimulator(params_boson).simulate(t_max=2.0, n_steps=5)
r_mb = QuditGKSLBosonSimulator(params_boson).simulate(t_max=2.0, n_steps=5)
max_diff_boson = max(
    abs(r_cb["populations"][-1][k] - r_mb["populations"][-1][k])
    for k in ["N_S0", "N_T1", "N_S1"]
)
check(
    "boson_circuit_matrix_match",
    max_diff_boson < 1e-10,
    f"Boson circuit matches matrix: max_diff={max_diff_boson:.2e}",
)

# Gate count verification
check(
    "qudit_N4_gate_count",
    r_c4["gates_per_step"] == 66,
    f"Qudit N=4 gates_per_step={r_c4['gates_per_step']} (expected 66)",
)

check(
    "qubit_N2_gate_count",
    r_cq2["gates_per_step"] == 30,
    f"Qubit N=2 gates_per_step={r_cq2['gates_per_step']} "
    f"(expected 30 = 2*(N+n_pairs)+2*(n_single+n_pair) = 2*(2+1)+2*(10+2))",
)

boson_gb = r_cb["gate_breakdown"]
check(
    "boson_stinespring_single_count",
    boson_gb["cu_two_stinespring_single"] == 20,
    f"Boson Stinespring single={boson_gb['cu_two_stinespring_single']} (expected 20 = 2×10)",
)

check(
    "boson_stinespring_pair_count",
    boson_gb["cu_multi_stinespring_pair"] == 4,
    f"Boson Stinespring pair={boson_gb['cu_multi_stinespring_pair']} (expected 4 = 2×2)",
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
passed = sum(1 for r in results if r["status"] == "PASS")
failed = sum(1 for r in results if r["status"] == "FAIL")
total = len(results)
print(f"Results: {passed}/{total} PASS, {failed}/{total} FAIL")

if failed > 0:
    print("\nFailed checks:")
    for r in results:
        if r["status"] == "FAIL":
            print(f"  ✗ {r['name']}: {r['detail']}")

# Save to file
timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
output_dir = os.path.join(os.path.dirname(__file__), "..", "developing", "verification_results")
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, f"iteration13_palindromic_trotter_{timestamp}.json")

output = {
    "iteration": 13,
    "timestamp": timestamp,
    "description": "Palindromic Trotter step fix for circuit simulators + gate count corrections",
    "summary": {"total_checks": total, "passed": passed, "failed": failed},
    "bugs_fixed": [
        {
            "id": "A",
            "file": "qudit_gksl_circuit_simulator.py",
            "method": "_trotter_step_circuit",
            "old": "Non-palindromic: H(dt/2) → D_1...D_n(dt) → H(dt/2)",
            "new": "Palindromic: H(dt/2) → D_1...D_n(dt/2) → D_n...D_1(dt/2) → H(dt/2)",
            "impact": "1st-order Lindblad accuracy → 2nd-order, matches matrix simulator",
        },
        {
            "id": "B",
            "file": "qubit_gksl_circuit_simulator.py",
            "method": "_trotter_step_circuit",
            "old": "Same non-palindromic issue as Bug A",
            "new": "Same palindromic fix as Bug A",
            "impact": "Same as Bug A for qubit encoding",
        },
        {
            "id": "C",
            "file": "qudit_gksl_circuit_boson_simulator.py",
            "method": "simulate (gate count)",
            "old": "Lindblad gates counted once despite palindromic Trotter",
            "new": "Lindblad gates counted twice (2× for forward + reverse)",
            "impact": "Gate count underreported by factor of 2",
        },
        {
            "id": "D",
            "file": "qudit_gksl_circuit_simulator.py",
            "method": "simulate (gate count)",
            "old": "Same gate count issue as Bug C (after Bug A fix)",
            "new": "Same fix as Bug C",
            "impact": "Same as Bug C",
        },
        {
            "id": "E",
            "file": "qubit_gksl_circuit_simulator.py",
            "method": "simulate (gate count)",
            "old": "Same gate count issue as Bug C",
            "new": "Same fix as Bug C",
            "impact": "Same as Bug C",
        },
        {
            "id": "F",
            "file": "qudit_gksl_circuit_simulator.py, qubit_gksl_circuit_simulator.py",
            "method": "build_full_trotter_step_circuit",
            "old": "Non-palindromic structure with full dt",
            "new": "Palindromic structure with dt/2, forward + reverse circuits",
            "impact": "Circuit visualization now matches actual simulation",
        },
    ],
    "checks": results,
}

with open(output_file, "w") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

print(f"\nResults saved to: {output_file}")
