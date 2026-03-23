#!/usr/bin/env python3
"""Iteration 48 verification script: Combined Trotter step circuit visualization.

Tests the new build_combined_trotter_step_circuit() method added to:
  - QubitGKSLCircuitSimulator (Scenario 3)
  - QuditGKSLCircuitSimulator (Scenario 5)
  - QuditGKSLCircuitBosonSimulator (Scenario 6)
And the manual combined circuit construction for Scenario 4 (Qubit Boson).

Outputs results to developing/verification_results/iteration48_results.txt
"""

import os
import sys
import traceback
from datetime import datetime

import numpy as np
from scipy.linalg import expm

# Ensure tutorials directory is on path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_physical_parameters import GKSLPhysicalParameters
from stinespring_utils import stinespring_unitary_from_lindblad

# Output setup
os.makedirs("../developing/verification_results", exist_ok=True)
output_path = "../developing/verification_results/iteration48_results.txt"
results = []
pass_count = 0
fail_count = 0


def check(name, condition, detail=""):
    global pass_count, fail_count
    status = "PASS" if condition else "FAIL"
    if condition:
        pass_count += 1
    else:
        fail_count += 1
    msg = f"[{status}] {name}"
    if detail:
        msg += f" -- {detail}"
    results.append(msg)
    print(msg)


def section(title):
    sep = "=" * 70
    results.append("")
    results.append(sep)
    results.append(title)
    results.append(sep)
    print(f"\n{sep}\n{title}\n{sep}")


# ======================================================================
# Parameters (matching notebook Cell 2)
# ======================================================================
params = GKSLPhysicalParameters(N_molecules=4)
params_boson = GKSLPhysicalParameters(
    N_molecules=4, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02
)
t_max = 100.0
n_steps_quantum = 200
dt = t_max / n_steps_quantum

# ======================================================================
# Section 1: Scenario 3 - Qubit GKSL (No Boson) Combined Circuit
# ======================================================================
section("Section 1: Scenario 3 - Qubit GKSL (No Boson)")

try:
    from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

    sim3 = QubitGKSLCircuitSimulator(params)

    # Build combined circuit
    combined3 = sim3.build_combined_trotter_step_circuit(dt)

    check("S3-01: build_combined_trotter_step_circuit returns dict",
          isinstance(combined3, dict))

    check("S3-02: dict has 'circuit' key",
          "circuit" in combined3)

    check("S3-03: dict has 'total_gates' key",
          "total_gates" in combined3)

    check("S3-04: dict has 'n_hamiltonian_gates' key",
          "n_hamiltonian_gates" in combined3)

    check("S3-05: dict has 'n_stinespring_gates' key",
          "n_stinespring_gates" in combined3)

    check("S3-06: dict has 'n_system_qubits' key",
          "n_system_qubits" in combined3)

    check("S3-07: dict has 'n_ancilla_qubits' key",
          "n_ancilla_qubits" in combined3)

    check("S3-08: dict has 'n_total_qubits' key",
          "n_total_qubits" in combined3)

    # Build separate circuits for comparison
    separate3 = sim3.build_full_trotter_step_circuit(dt)

    check("S3-09: total_gates matches separate circuits",
          combined3["total_gates"] == separate3["total_gates"],
          f"combined={combined3['total_gates']}, separate={separate3['total_gates']}")

    check("S3-10: n_hamiltonian_gates matches",
          combined3["n_hamiltonian_gates"] == separate3["n_hamiltonian_gates"],
          f"combined={combined3['n_hamiltonian_gates']}, separate={separate3['n_hamiltonian_gates']}")

    check("S3-11: n_stinespring_gates matches",
          combined3["n_stinespring_gates"] == separate3["n_stinespring_gates"],
          f"combined={combined3['n_stinespring_gates']}, separate={separate3['n_stinespring_gates']}")

    N = params.N_molecules
    expected_sys = 2 * N
    expected_total = expected_sys + 1

    check("S3-12: n_system_qubits = 2*N",
          combined3["n_system_qubits"] == expected_sys,
          f"got {combined3['n_system_qubits']}, expected {expected_sys}")

    check("S3-13: n_ancilla_qubits = 1",
          combined3["n_ancilla_qubits"] == 1)

    check("S3-14: n_total_qubits = 2*N + 1",
          combined3["n_total_qubits"] == expected_total,
          f"got {combined3['n_total_qubits']}, expected {expected_total}")

    # Circuit object checks
    circ3 = combined3["circuit"]
    check("S3-15: circuit is Qiskit QuantumCircuit",
          type(circ3).__name__ == "QuantumCircuit",
          f"type={type(circ3).__name__}")

    check("S3-16: circuit has correct num_qubits",
          circ3.num_qubits == expected_total,
          f"got {circ3.num_qubits}, expected {expected_total}")

    # Gate count from Qiskit circuit
    # Note: reset operations are also in the circuit but not counted as gates
    n_unitary_gates = sum(1 for inst in circ3.data if inst.operation.name != "reset")
    n_reset_ops = sum(1 for inst in circ3.data if inst.operation.name == "reset")
    n_lindblad_channels = len(sim3.lindblad_local_info)

    check("S3-17: number of unitary gates matches total_gates",
          n_unitary_gates == combined3["total_gates"],
          f"unitary_gates={n_unitary_gates}, total_gates={combined3['total_gates']}")

    check("S3-18: number of reset operations = 2 * n_lindblad_channels",
          n_reset_ops == 2 * n_lindblad_channels,
          f"resets={n_reset_ops}, expected={2 * n_lindblad_channels}")

    # Regression: expected values
    check("S3-19: total gates = 66",
          combined3["total_gates"] == 66,
          f"got {combined3['total_gates']}")

    check("S3-20: n_hamiltonian_gates = 14",
          combined3["n_hamiltonian_gates"] == 14,
          f"got {combined3['n_hamiltonian_gates']}")

    check("S3-21: n_stinespring_gates = 52",
          combined3["n_stinespring_gates"] == 52,
          f"got {combined3['n_stinespring_gates']}")

except Exception as e:
    check(f"S3-ERROR: Exception during Scenario 3 tests", False,
          f"{type(e).__name__}: {e}")
    traceback.print_exc()

# ======================================================================
# Section 2: Scenario 5 - Qudit GKSL (No Boson) Combined Circuit
# ======================================================================
section("Section 2: Scenario 5 - Qudit GKSL (No Boson)")

try:
    from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

    sim5 = QuditGKSLCircuitSimulator(params)

    # Build combined circuit
    combined5 = sim5.build_combined_trotter_step_circuit(dt)

    check("S5-01: build_combined_trotter_step_circuit returns dict",
          isinstance(combined5, dict))

    check("S5-02: dict has 'circuit' key",
          "circuit" in combined5)

    # Build separate circuits for comparison
    separate5 = sim5.build_full_trotter_step_circuit(dt)

    check("S5-03: total_gates matches separate circuits",
          combined5["total_gates"] == separate5["total_gates"],
          f"combined={combined5['total_gates']}, separate={separate5['total_gates']}")

    check("S5-04: n_hamiltonian_gates matches",
          combined5["n_hamiltonian_gates"] == separate5["n_hamiltonian_gates"],
          f"combined={combined5['n_hamiltonian_gates']}, separate={separate5['n_hamiltonian_gates']}")

    check("S5-05: n_stinespring_gates matches",
          combined5["n_stinespring_gates"] == separate5["n_stinespring_gates"],
          f"combined={combined5['n_stinespring_gates']}, separate={separate5['n_stinespring_gates']}")

    N = params.N_molecules
    expected_sys = N
    expected_total = N + 1

    check("S5-06: n_system_qudits = N",
          combined5["n_system_qudits"] == expected_sys,
          f"got {combined5['n_system_qudits']}, expected {expected_sys}")

    check("S5-07: n_ancilla_qudits = 1",
          combined5["n_ancilla_qudits"] == 1)

    check("S5-08: n_total_qudits = N + 1",
          combined5["n_total_qudits"] == expected_total,
          f"got {combined5['n_total_qudits']}, expected {expected_total}")

    # Circuit object checks
    circ5 = combined5["circuit"]
    check("S5-09: circuit is MQT-Qudits QuantumCircuit",
          type(circ5).__name__ == "QuantumCircuit",
          f"type={type(circ5).__name__}")

    check("S5-10: circuit has correct num_qudits",
          circ5.num_qudits == expected_total,
          f"got {circ5.num_qudits}, expected {expected_total}")

    expected_dims = [3] * N + [3]  # system d=3 + ancilla d=3
    check("S5-11: circuit dimensions correct",
          circ5.dimensions == expected_dims,
          f"got {circ5.dimensions}, expected {expected_dims}")

    # Gate count from circuit instructions
    n_instructions = len(circ5.instructions)
    check("S5-12: number of circuit instructions matches total_gates",
          n_instructions == combined5["total_gates"],
          f"instructions={n_instructions}, total_gates={combined5['total_gates']}")

    # Regression
    check("S5-13: total gates = 66",
          combined5["total_gates"] == 66,
          f"got {combined5['total_gates']}")

    check("S5-14: n_hamiltonian_gates = 14",
          combined5["n_hamiltonian_gates"] == 14,
          f"got {combined5['n_hamiltonian_gates']}")

    check("S5-15: n_stinespring_gates = 52",
          combined5["n_stinespring_gates"] == 52,
          f"got {combined5['n_stinespring_gates']}")

    # Gate target checks: Hamiltonian gates should NOT touch ancilla
    from mqt.qudits.quantum_circuit.components.extensions.gate_types import GateTypes

    ancilla_idx = N
    ham_gates_touch_ancilla = False
    for gate in circ5.instructions[:7]:  # First 7 = hamiltonian_half_1
        targets = gate.target_qudits if isinstance(gate.target_qudits, list) else [gate.target_qudits]
        if ancilla_idx in targets:
            ham_gates_touch_ancilla = True
            break
    check("S5-16: Hamiltonian gates do not touch ancilla qudit",
          not ham_gates_touch_ancilla)

    # Stinespring gates should touch ancilla
    st_gates_touch_ancilla = True
    for gate in circ5.instructions[7:33]:  # Gates 7-32 = forward Stinespring
        targets = gate.target_qudits if isinstance(gate.target_qudits, list) else [gate.target_qudits]
        if ancilla_idx not in targets:
            st_gates_touch_ancilla = False
            break
    check("S5-17: All Stinespring gates touch ancilla qudit",
          st_gates_touch_ancilla)

except Exception as e:
    check(f"S5-ERROR: Exception during Scenario 5 tests", False,
          f"{type(e).__name__}: {e}")
    traceback.print_exc()

# ======================================================================
# Section 3: Scenario 6 - Qudit GKSL (With Boson) Combined Circuit
# ======================================================================
section("Section 3: Scenario 6 - Qudit GKSL (With Boson)")

try:
    from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

    sim6 = QuditGKSLCircuitBosonSimulator(params_boson)

    # Build combined circuit
    combined6 = sim6.build_combined_trotter_step_circuit(dt)

    check("S6-01: build_combined_trotter_step_circuit returns dict",
          isinstance(combined6, dict))

    check("S6-02: dict has required keys",
          all(k in combined6 for k in [
              "circuit", "total_gates", "n_hamiltonian_gates",
              "n_stinespring_gates", "n_system_qudits",
              "n_ancilla_qudits", "n_total_qudits"
          ]))

    N = params_boson.N_molecules
    n_max = params_boson.n_max
    d = 3
    d_ph = n_max + 1  # d_ph = 2 for n_max=1

    expected_sys = 2 * N  # electronic + phonon
    expected_total = expected_sys + 1

    check("S6-03: n_system_qudits = 2*N",
          combined6["n_system_qudits"] == expected_sys,
          f"got {combined6['n_system_qudits']}, expected {expected_sys}")

    check("S6-04: n_total_qudits = 2*N + 1",
          combined6["n_total_qudits"] == expected_total,
          f"got {combined6['n_total_qudits']}, expected {expected_total}")

    # Circuit object checks
    circ6 = combined6["circuit"]
    check("S6-05: circuit has correct num_qudits",
          circ6.num_qudits == expected_total,
          f"got {circ6.num_qudits}, expected {expected_total}")

    expected_dims = [d] * N + [d_ph] * N + [d]  # el + ph + ancilla(d=3)
    check("S6-06: circuit dimensions correct",
          circ6.dimensions == expected_dims,
          f"got {circ6.dimensions}, expected {expected_dims}")

    # Gate count checks
    n_instructions = len(circ6.instructions)
    check("S6-07: instructions match total_gates",
          n_instructions == combined6["total_gates"],
          f"instructions={n_instructions}, total_gates={combined6['total_gates']}")

    # Expected Hamiltonian gate count per half-step:
    # electronic on-site: N=4, electronic transfer: 3 (chain), phonon on-site: N=4, eph coupling: N=4
    n_ham_per_half = N + len(params_boson.neighbors) + N + N  # 4+3+4+4=15
    expected_ham = 2 * n_ham_per_half

    check("S6-08: n_hamiltonian_gates correct",
          combined6["n_hamiltonian_gates"] == expected_ham,
          f"got {combined6['n_hamiltonian_gates']}, expected {expected_ham}")

    # Expected Stinespring: 2 * n_lindblad_channels
    n_lindblad = len(sim6.lindblad_local_info)
    expected_st = 2 * n_lindblad

    check("S6-09: n_stinespring_gates correct",
          combined6["n_stinespring_gates"] == expected_st,
          f"got {combined6['n_stinespring_gates']}, expected {expected_st}")

    expected_total_gates = expected_ham + expected_st
    check("S6-10: total_gates = n_ham + n_st",
          combined6["total_gates"] == expected_total_gates,
          f"got {combined6['total_gates']}, expected {expected_total_gates}")

    # Ancilla index check
    ancilla_idx_6 = 2 * N
    # Stinespring gates should reference ancilla
    st_start = n_ham_per_half  # After first Hamiltonian half
    st_end = st_start + n_lindblad
    st_touch_anc = True
    for gate in circ6.instructions[st_start:st_end]:
        targets = gate.target_qudits if isinstance(gate.target_qudits, list) else [gate.target_qudits]
        if ancilla_idx_6 not in targets:
            st_touch_anc = False
            break
    check("S6-11: Forward Stinespring gates touch ancilla",
          st_touch_anc)

except Exception as e:
    check(f"S6-ERROR: Exception during Scenario 6 tests", False,
          f"{type(e).__name__}: {e}")
    traceback.print_exc()

# ======================================================================
# Section 4: Scenario 4 - Qubit GKSL (With Boson) Combined Circuit
# ======================================================================
section("Section 4: Scenario 4 - Qubit GKSL (With Boson)")

try:
    from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator

    # Use N=2 for CI feasibility (N=4 requires ~7GB RAM)
    params_boson_n2 = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02
    )
    sim4 = QubitGKSLBosonSimulator(params_boson_n2)
    n_sys_qubits = sim4.n_el_qubits + sim4.n_ph_qubits
    dim_sys = 2 ** n_sys_qubits
    n_total_qubits = n_sys_qubits + 1
    ancilla_idx_4 = n_sys_qubits

    from qiskit import QuantumCircuit
    from qiskit.circuit.library import UnitaryGate

    def _embed_unitary_to_power_of_two(U, target_dim):
        embedded = np.eye(target_dim, dtype=np.complex128)
        embedded[:U.shape[0], :U.shape[1]] = U
        return embedded

    combined_s4 = QuantumCircuit(n_total_qubits)

    # Hamiltonian half-step 1
    U_half = expm(-1j * sim4.H_total * (dt / 2))
    U_half_emb = _embed_unitary_to_power_of_two(U_half, dim_sys)
    combined_s4.append(UnitaryGate(U_half_emb, label="H_half"), range(n_sys_qubits))

    n_ham = 1
    n_st = 0

    # Forward Lindblad
    for idx, (L_op, _gamma) in enumerate(sim4.lindblad_ops):
        combined_s4.reset(ancilla_idx_4)
        U_st = stinespring_unitary_from_lindblad(L_op, dt / 2)
        U_st_emb = _embed_unitary_to_power_of_two(U_st, 2 * dim_sys)
        combined_s4.append(
            UnitaryGate(U_st_emb, label=f"D_{idx}"),
            range(n_total_qubits),
        )
        n_st += 1

    # Reverse Lindblad (palindromic)
    for idx, (L_op, _gamma) in reversed(list(enumerate(sim4.lindblad_ops))):
        combined_s4.reset(ancilla_idx_4)
        U_st = stinespring_unitary_from_lindblad(L_op, dt / 2)
        U_st_emb = _embed_unitary_to_power_of_two(U_st, 2 * dim_sys)
        combined_s4.append(
            UnitaryGate(U_st_emb, label=f"D_{idx}"),
            range(n_total_qubits),
        )
        n_st += 1

    # Hamiltonian half-step 2
    combined_s4.append(UnitaryGate(U_half_emb, label="H_half"), range(n_sys_qubits))
    n_ham += 1

    check("S4-01: Scenario 4 combined circuit built successfully", True)

    check("S4-02: circuit has correct num_qubits",
          combined_s4.num_qubits == n_total_qubits,
          f"got {combined_s4.num_qubits}, expected {n_total_qubits}")

    n_unitary = sum(1 for inst in combined_s4.data if inst.operation.name != "reset")
    n_resets = sum(1 for inst in combined_s4.data if inst.operation.name == "reset")
    total_gates_s4 = n_ham + n_st
    n_lindblad_s4 = len(sim4.lindblad_ops)

    check("S4-03: number of unitary gates correct",
          n_unitary == total_gates_s4,
          f"unitary={n_unitary}, expected={total_gates_s4}")

    check("S4-04: number of reset operations = 2*n_lindblad",
          n_resets == 2 * n_lindblad_s4,
          f"resets={n_resets}, expected={2 * n_lindblad_s4}")

    check("S4-05: n_sys_qubits > 0",
          n_sys_qubits > 0,
          f"n_sys_qubits={n_sys_qubits}")

    check("S4-06: Hamiltonian gates = 2",
          n_ham == 2)

    check("S4-07: Stinespring gates = 2*n_lindblad",
          n_st == 2 * n_lindblad_s4,
          f"n_st={n_st}, expected={2 * n_lindblad_s4}")

except Exception as e:
    check(f"S4-ERROR: Exception during Scenario 4 tests", False,
          f"{type(e).__name__}: {e}")
    traceback.print_exc()

# ======================================================================
# Section 5: Cross-scenario consistency checks
# ======================================================================
section("Section 5: Cross-scenario consistency")

try:
    # Ensure both combined3 and combined5 exist
    _c3 = combined3  # noqa: F841 - check existence
    _c5 = combined5  # noqa: F841

    # Scenario 3 vs 5: same physics, different encoding
    # Gate counts should match
    check("X-01: S3 and S5 have same total_gates",
          combined3["total_gates"] == combined5["total_gates"],
          f"S3={combined3['total_gates']}, S5={combined5['total_gates']}")

    check("X-02: S3 and S5 have same n_hamiltonian_gates",
          combined3["n_hamiltonian_gates"] == combined5["n_hamiltonian_gates"],
          f"S3={combined3['n_hamiltonian_gates']}, S5={combined5['n_hamiltonian_gates']}")

    check("X-03: S3 and S5 have same n_stinespring_gates",
          combined3["n_stinespring_gates"] == combined5["n_stinespring_gates"],
          f"S3={combined3['n_stinespring_gates']}, S5={combined5['n_stinespring_gates']}")

    # Qudit combined should have fewer qudits than qubit combined has qubits
    check("X-04: Qudit uses fewer registers than qubit",
          combined5["n_total_qudits"] < combined3["n_total_qubits"],
          f"qudit={combined5['n_total_qudits']}, qubit={combined3['n_total_qubits']}")

except NameError as e:
    check("X-SKIP: Cross-checks skipped (prerequisite test failed)", False,
          f"Missing variable: {e}")
except Exception as e:
    check(f"X-ERROR: Exception during cross-checks", False,
          f"{type(e).__name__}: {e}")
    traceback.print_exc()

# ======================================================================
# Section 6: Visualization feasibility (no display, just structure checks)
# ======================================================================
section("Section 6: Visualization feasibility")

try:
    # Check that Qiskit circuits can be drawn (no display)
    from qiskit.visualization import circuit_drawer
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Scenario 3 combined
    fig3 = circuit_drawer(combined3["circuit"], output="mpl", fold=40)
    check("V-01: Scenario 3 combined circuit drawable",
          fig3 is not None)
    plt.close(fig3)

    # Scenario 4 combined
    fig4 = circuit_drawer(combined_s4, output="mpl", fold=40)
    check("V-02: Scenario 4 combined circuit drawable",
          fig4 is not None)
    plt.close(fig4)

    # MQT-Qudits circuits: check structure for custom renderer
    circ5 = combined5["circuit"]
    check("V-03: Scenario 5 circuit has instructions list",
          hasattr(circ5, "instructions") and len(circ5.instructions) > 0,
          f"n_instructions={len(circ5.instructions)}")

    check("V-04: Scenario 5 circuit has dimensions",
          hasattr(circ5, "dimensions") and len(circ5.dimensions) > 0,
          f"dims={circ5.dimensions}")

    circ6 = combined6["circuit"]
    check("V-05: Scenario 6 circuit has instructions list",
          hasattr(circ6, "instructions") and len(circ6.instructions) > 0,
          f"n_instructions={len(circ6.instructions)}")

    check("V-06: Scenario 6 circuit has dimensions",
          hasattr(circ6, "dimensions") and len(circ6.dimensions) > 0,
          f"dims={circ6.dimensions}")

except Exception as e:
    check(f"V-ERROR: Exception during visualization checks", False,
          f"{type(e).__name__}: {e}")
    traceback.print_exc()

# ======================================================================
# Summary
# ======================================================================
section("SUMMARY")

total = pass_count + fail_count
summary = f"Total: {total} checks, {pass_count} PASS, {fail_count} FAIL"
results.append(summary)
print(summary)

timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
results.insert(0, f"Iteration 48 Verification Results - {timestamp}")
results.insert(1, "=" * 70)

with open(output_path, "w") as f:
    f.write("\n".join(results))

print(f"\nResults saved to: {output_path}")
sys.exit(0 if fail_count == 0 else 1)
