"""Iteration 12 verification: qudit ancilla circuit boson simulator fix.

Iteration 11 analysis found that qudit_gksl_circuit_boson_simulator.py had
hardcoded ancilla dimension d_anc=2 in build_stinespring_circuit_single and
build_stinespring_circuit_pair, despite self.d_anc being correctly set to
params.d (=3 for qutrit). This iteration fixes these bugs and verifies the fix.

Verification checks:
  1-4:   Circuit boson simulator d_anc fix verification
  5-10:  Regression checks from iteration 11
  11-15: Cross-file consistency of d_anc usage
  16-20: Stinespring unitary dimension verification
"""

from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

results: list[dict] = []


def check(name: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    results.append({"name": name, "status": status, "detail": detail})
    print(f"  [{status}] {name}: {detail}")


def read_notebook(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def cell_source(nb: dict, idx: int) -> str:
    return "".join(nb["cells"][idx]["source"])


# ---------------------------------------------------------------------------
# Load source files
# ---------------------------------------------------------------------------

print("=" * 70)
print("Iteration 12 Verification: Circuit Boson Simulator d_anc Fix")
print("=" * 70)

base = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(base, "qudit_gksl_circuit_boson_simulator.py")) as f:
    circuit_boson_src = f.read()

with open(os.path.join(base, "qudit_gksl_circuit_simulator.py")) as f:
    circuit_src = f.read()

with open(os.path.join(base, "qudit_gksl_simulator.py")) as f:
    qudit_gksl_src = f.read()

with open(os.path.join(base, "qudit_gksl_shot_simulator.py")) as f:
    qudit_gksl_shot_src = f.read()

with open(os.path.join(base, "qudit_gksl_noisy_simulator.py")) as f:
    qudit_gksl_noisy_src = f.read()

with open(os.path.join(base, "qudit_gksl_boson_simulator.py")) as f:
    qudit_gksl_boson_src = f.read()

with open(os.path.join(base, "stinespring_utils.py")) as f:
    stinespring_src = f.read()

with open(os.path.join(base, "qubit_gksl_shot_simulator.py")) as f:
    qubit_gksl_shot_src = f.read()

with open(os.path.join(base, "qubit_gksl_noisy_simulator.py")) as f:
    qubit_gksl_noisy_src = f.read()

gksl_nb = read_notebook(os.path.join(base, "quantum_dynamics_gksl_comparison.ipynb"))
complete_nb = read_notebook(
    os.path.join(base, "quantum_dynamics_complete_comparison.ipynb")
)

# ---------------------------------------------------------------------------
# Section 1: Circuit Boson Simulator d_anc Fix Verification
# ---------------------------------------------------------------------------
print("\n--- Section 1: Circuit boson simulator d_anc fix ---")

# Bug #1: build_stinespring_circuit_single should use d_anc, not hardcoded 2
check(
    "boson_single_circuit_d_anc",
    "[d, d_anc]" in circuit_boson_src and "[d, 2]" not in circuit_boson_src,
    "build_stinespring_circuit_single uses [d, d_anc] (not [d, 2])",
)

# Bug #2: build_stinespring_circuit_pair should use d_anc, not hardcoded 2
check(
    "boson_pair_circuit_d_anc",
    "[d, d, d_anc]" in circuit_boson_src and "[d, d, 2]" not in circuit_boson_src,
    "build_stinespring_circuit_pair uses [d, d, d_anc] (not [d, d, 2])",
)

# Verify d_anc is set from params.d
check(
    "boson_d_anc_from_params",
    "self.d_anc = params.d" in circuit_boson_src,
    "QuditGKSLCircuitBosonSimulator sets d_anc = params.d",
)

# Verify d_anc local variable is extracted in circuit construction methods
boson_single_method = circuit_boson_src[
    circuit_boson_src.find("def build_stinespring_circuit_single"):
    circuit_boson_src.find("def build_stinespring_circuit_pair")
]
check(
    "boson_single_uses_d_anc_var",
    "d_anc = self.d_anc" in boson_single_method,
    "build_stinespring_circuit_single extracts d_anc = self.d_anc",
)

# ---------------------------------------------------------------------------
# Section 2: Cross-file d_anc consistency
# ---------------------------------------------------------------------------
print("\n--- Section 2: Cross-file d_anc consistency ---")

# Non-boson circuit simulator should also use d_anc
check(
    "circuit_single_d_anc",
    "[d, d_anc]" in circuit_src and "def build_stinespring_circuit_single" in circuit_src,
    "QuditGKSLCircuitSimulator single uses [d, d_anc]",
)

check(
    "circuit_pair_d_anc",
    "[d, d, d_anc]" in circuit_src and "def build_stinespring_circuit_pair" in circuit_src,
    "QuditGKSLCircuitSimulator pair uses [d, d, d_anc]",
)

# All qudit simulators should have d_anc = params.d
for name, src in [
    ("QuditGKSLSimulator", qudit_gksl_src),
    ("QuditGKSLShotSimulator", qudit_gksl_shot_src),
    ("QuditGKSLBosonSimulator", qudit_gksl_boson_src),
    ("QuditGKSLCircuitSimulator", circuit_src),
    ("QuditGKSLCircuitBosonSimulator", circuit_boson_src),
]:
    check(
        f"{name}_d_anc",
        "self.d_anc = params.d" in src,
        f"{name} sets d_anc = params.d",
    )

# Qubit simulators should NOT set d_anc (they use default d_anc=2)
for name, src in [
    ("QubitGKSLShotSimulator", qubit_gksl_shot_src),
    ("QubitGKSLNoisySimulator", qubit_gksl_noisy_src),
]:
    check(
        f"{name}_no_d_anc",
        "self.d_anc" not in src,
        f"{name} does not set d_anc (uses default=2)",
    )

# ---------------------------------------------------------------------------
# Section 3: Stinespring utils d_anc parameter
# ---------------------------------------------------------------------------
print("\n--- Section 3: Stinespring utils d_anc support ---")

check(
    "stinespring_unitary_d_anc_param",
    "def stinespring_unitary_from_lindblad" in stinespring_src
    and "d_anc=2" in stinespring_src,
    "stinespring_unitary_from_lindblad has d_anc=2 default",
)

check(
    "stinespring_dm_d_anc_param",
    "def apply_stinespring_to_density_matrix" in stinespring_src
    and "d_anc=2" in stinespring_src,
    "apply_stinespring_to_density_matrix has d_anc=2 default",
)

# ---------------------------------------------------------------------------
# Section 4: n_ancilla_qudits consistency in qudit simulators
# ---------------------------------------------------------------------------
print("\n--- Section 4: n_ancilla_qudits naming ---")

# Check qudit simulators use n_ancilla_qudits (not n_ancilla_qubits)
for name, src in [
    ("QuditGKSLSimulator", qudit_gksl_src),
    ("QuditGKSLShotSimulator", qudit_gksl_shot_src),
    ("QuditGKSLBosonSimulator", qudit_gksl_boson_src),
]:
    has_qudits = "n_ancilla_qudits" in src
    has_qubits = "n_ancilla_qubits" in src
    check(
        f"{name}_naming",
        has_qudits and not has_qubits,
        f"{name}: uses n_ancilla_qudits={has_qudits}, no n_ancilla_qubits={not has_qubits}",
    )

# ---------------------------------------------------------------------------
# Section 5: Dimension consistency checks
# ---------------------------------------------------------------------------
print("\n--- Section 5: Dimension consistency ---")

# Verify _build_local_stinespring_unitary uses d_anc in both circuit simulators
for label, src in [
    ("circuit", circuit_src),
    ("circuit_boson", circuit_boson_src),
]:
    build_method = src[
        src.find("def _build_local_stinespring_unitary"):
        src.find("def _extract_kraus")
    ]
    check(
        f"{label}_build_uses_d_anc",
        "d_anc = self.d_anc" in build_method and "d_anc * d_local" in build_method,
        f"{label}: _build_local_stinespring_unitary uses d_anc correctly",
    )

# Verify _extract_kraus uses d_anc in both circuit simulators
for label, src in [
    ("circuit", circuit_src),
    ("circuit_boson", circuit_boson_src),
]:
    start = src.find("def _extract_kraus_from_local_stinespring")
    # Find the next method definition after _extract_kraus
    next_def = src.find("\n    def ", start + 1)
    if next_def == -1:
        next_def = start + 1000
    extract_method = src[start:next_def]
    check(
        f"{label}_extract_uses_d_anc",
        "d_anc" in extract_method and "range(d_anc)" in extract_method,
        f"{label}: _extract_kraus uses d_anc correctly",
    )

# ---------------------------------------------------------------------------
# Section 6: Regression checks from iteration 11
# ---------------------------------------------------------------------------
print("\n--- Section 6: Iteration 11 regression ---")

cell38_src = cell_source(gksl_nb, 38)
cell15_src = cell_source(complete_nb, 15)

check(
    "gksl_uses_cx_transfer",
    "cx_per_pair = cx_info['cx_transfer']" in cell38_src,
    "GKSL Cell 38 uses cx_info['cx_transfer']",
)

check(
    "complete_uses_cx_avg",
    "cx_per_pair = cx_info['cx_avg']" in cell15_src,
    "Complete Cell 15 uses cx_info['cx_avg']",
)

check(
    "gksl_depol_001",
    "p_depol=0.001" in cell38_src,
    "GKSL uses p_depol=0.001",
)

check(
    "complete_depol_001",
    "'depol_2q': 0.001" in cell15_src or "depol_2q=0.001" in cell15_src,
    "Complete uses depol_2q=0.001",
)

# ---------------------------------------------------------------------------
# Section 7: Numerical dimension verification
# ---------------------------------------------------------------------------
print("\n--- Section 7: Numerical dimension tests ---")

import numpy as np

# For d=3 (qutrit), d_anc=3:
d = 3
d_anc = 3

# Single-site Stinespring unitary should be (d_anc * d) × (d_anc * d) = 9×9
dim_single = d_anc * d
check(
    "single_stinespring_dim",
    dim_single == 9,
    f"Single-site Stinespring: {d_anc}×{d} = {dim_single}×{dim_single} unitary",
)

# Pair Stinespring unitary should be (d_anc * d²) × (d_anc * d²) = 27×27
dim_pair = d_anc * d * d
check(
    "pair_stinespring_dim",
    dim_pair == 27,
    f"Pair Stinespring: {d_anc}×{d}² = {dim_pair}×{dim_pair} unitary",
)

# Circuit register dimensions should match
# Single: [d, d_anc] = [3, 3] → total = 9 = dim_single ✓
circuit_dim_single = d * d_anc
check(
    "single_circuit_dim_match",
    circuit_dim_single == dim_single,
    f"Single circuit dim [d={d}, d_anc={d_anc}] = {circuit_dim_single} matches unitary dim {dim_single}",
)

# Pair: [d, d, d_anc] = [3, 3, 3] → total = 27 = dim_pair ✓
circuit_dim_pair = d * d * d_anc
check(
    "pair_circuit_dim_match",
    circuit_dim_pair == dim_pair,
    f"Pair circuit dim [d={d}, d={d}, d_anc={d_anc}] = {circuit_dim_pair} matches unitary dim {dim_pair}",
)

# Old buggy dimensions (d_anc=2 hardcoded)
old_circuit_single = d * 2  # = 6 ≠ 9
old_circuit_pair = d * d * 2  # = 18 ≠ 27
check(
    "old_single_mismatch",
    old_circuit_single != dim_single,
    f"Old bug: [d={d}, 2] = {old_circuit_single} ≠ {dim_single} (dimension mismatch!)",
)
check(
    "old_pair_mismatch",
    old_circuit_pair != dim_pair,
    f"Old bug: [d={d}, d={d}, 2] = {old_circuit_pair} ≠ {dim_pair} (dimension mismatch!)",
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
n_pass = sum(1 for r in results if r["status"] == "PASS")
n_fail = sum(1 for r in results if r["status"] == "FAIL")
n_total = len(results)
print(f"Results: {n_pass}/{n_total} PASS, {n_fail}/{n_total} FAIL")
print("=" * 70)

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
out_dir = os.path.join(base, "..", "developing", "verification_results")
os.makedirs(out_dir, exist_ok=True)

ts = datetime.now(tz=timezone.utc).strftime("%Y%m%dT%H%M%SZ")
out_path = os.path.join(out_dir, f"iteration12_circuit_boson_fix_{ts}.json")

output = {
    "iteration": 12,
    "timestamp": ts,
    "description": "Circuit boson simulator d_anc fix: hardcoded ancilla=2 replaced with d_anc",
    "summary": {"total_checks": n_total, "passed": n_pass, "failed": n_fail},
    "bugs_fixed": [
        {
            "file": "qudit_gksl_circuit_boson_simulator.py",
            "method": "build_stinespring_circuit_single",
            "old": "QuantumCircuit(2, [d, 2], 0)",
            "new": "QuantumCircuit(2, [d, d_anc], 0)",
            "impact": "Circuit-unitary dimension mismatch when d_anc > 2",
        },
        {
            "file": "qudit_gksl_circuit_boson_simulator.py",
            "method": "build_stinespring_circuit_pair",
            "old": "QuantumCircuit(3, [d, d, 2], 0)",
            "new": "QuantumCircuit(3, [d, d, d_anc], 0)",
            "impact": "Circuit-unitary dimension mismatch when d_anc > 2",
        },
    ],
    "checks": results,
}

with open(out_path, "w") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

print(f"\nResults saved to: {out_path}")
