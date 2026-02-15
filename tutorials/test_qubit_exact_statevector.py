#!/usr/bin/env python3
"""Direct statevector-based validation (no shot noise).

This validates the exact qubit implementation against classical simulation
using exact statevector evolution (not shot-based sampling).
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "."))

from exact_qubit_hamiltonians import build_H_transfer_qubit_unitary, build_H_TTA_qubit_unitary
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector


class PhysicalParameters:
    """Physical parameters for the system."""

    def __init__(self) -> None:
        self.N_molecules = 4
        self.E_T = 1.5  # eV
        self.E_S = 3.0  # eV
        self.V = 0.1  # eV
        self.J = 0.05  # eV
        self.hbar = 0.6582119569  # eV·fs
        self.neighbors = [(0, 1), (1, 2), (2, 3)]


def apply_H0_evolution_exact(circuit, mol_idx, params, dt) -> None:
    """Apply H0 evolution to a single molecule."""
    q0 = 2 * mol_idx
    q1 = 2 * mol_idx + 1

    E_T = params.E_T
    E_S = params.E_S
    hbar = params.hbar

    # Pauli decomposition
    (E_T + E_S) / 4
    beta = (E_S - E_T) / 4
    gamma = (E_T - E_S) / 4
    delta = -(E_T + E_S) / 4

    theta_0 = -2 * beta * dt / hbar
    theta_1 = -2 * gamma * dt / hbar
    theta_zz = -2 * delta * dt / hbar

    circuit.rz(theta_0, q0)
    circuit.rz(theta_1, q1)

    # Z⊗Z interaction
    circuit.cx(q0, q1)
    circuit.rz(theta_zz, q1)
    circuit.cx(q0, q1)


def build_qubit_trotter_step_exact(params, dt):
    """Build exact Trotter step using exact unitaries."""
    N = params.N_molecules
    n_qubits = 2 * N
    circuit = QuantumCircuit(n_qubits)

    # Forward
    for i in range(N):
        apply_H0_evolution_exact(circuit, i, params, dt / 2)

    for i, j in params.neighbors:
        # H_transfer
        U_tr = build_H_transfer_qubit_unitary(params.V, dt / 2, params.hbar)
        qi0, qi1 = 2 * i, 2 * i + 1
        qj0, qj1 = 2 * j, 2 * j + 1
        qubits = [qi0, qi1, qj0, qj1]
        gate_tr = UnitaryGate(U_tr, label="U_tr")
        circuit.append(gate_tr, qubits)

    for i, j in params.neighbors:
        # H_TTA
        U_tta = build_H_TTA_qubit_unitary(params.J, dt / 2, params.hbar)
        qi0, qi1 = 2 * i, 2 * i + 1
        qj0, qj1 = 2 * j, 2 * j + 1
        qubits = [qi0, qi1, qj0, qj1]
        gate_tta = UnitaryGate(U_tta, label="U_TTA")
        circuit.append(gate_tta, qubits)

    # Backward
    for i, j in reversed(params.neighbors):
        U_tta = build_H_TTA_qubit_unitary(params.J, dt / 2, params.hbar)
        qi0, qi1 = 2 * i, 2 * i + 1
        qj0, qj1 = 2 * j, 2 * j + 1
        qubits = [qi0, qi1, qj0, qj1]
        gate_tta = UnitaryGate(U_tta, label="U_TTA")
        circuit.append(gate_tta, qubits)

    for i, j in reversed(params.neighbors):
        U_tr = build_H_transfer_qubit_unitary(params.V, dt / 2, params.hbar)
        qi0, qi1 = 2 * i, 2 * i + 1
        qj0, qj1 = 2 * j, 2 * j + 1
        qubits = [qi0, qi1, qj0, qj1]
        gate_tr = UnitaryGate(U_tr, label="U_tr")
        circuit.append(gate_tr, qubits)

    for i in reversed(range(N)):
        apply_H0_evolution_exact(circuit, i, params, dt / 2)

    return circuit


def calculate_populations_from_statevector(statevector, N_molecules):
    """Calculate populations from statevector."""
    probabilities = statevector.probabilities_dict()

    N_S0 = N_T1 = N_S1 = 0.0
    unphysical = 0.0

    for bitstring, prob in probabilities.items():
        if prob < 1e-15:
            continue

        # Qiskit uses big-endian for display, but we need little-endian for indexing
        bits = bitstring[::-1]  # Reverse to get little-endian

        is_unphysical = False
        for mol in range(N_molecules):
            q0_bit = int(bits[2 * mol])  # Right qubit
            q1_bit = int(bits[2 * mol + 1])  # Left qubit

            # Check for |11⟩ (unphysical)
            if q0_bit == 1 and q1_bit == 1:
                is_unphysical = True
                break

        if is_unphysical:
            unphysical += prob
        else:
            for mol in range(N_molecules):
                q0_bit = int(bits[2 * mol])
                q1_bit = int(bits[2 * mol + 1])

                if q1_bit == 0 and q0_bit == 0:  # |00⟩ = S0
                    N_S0 += prob
                elif q1_bit == 0 and q0_bit == 1:  # |01⟩ = T1
                    N_T1 += prob
                elif q1_bit == 1 and q0_bit == 0:  # |10⟩ = S1
                    N_S1 += prob

    return {"N_S0": N_S0, "N_T1": N_T1, "N_S1": N_S1, "unphysical": unphysical}


def run_qubit_simulation_exact(params, T_total, N_steps):
    """Run qubit simulation with exact statevector evolution."""
    dt = T_total / N_steps
    N = params.N_molecules
    n_qubits = 2 * N

    # Prepare initial state: edge_triplet
    init_circuit = QuantumCircuit(n_qubits)
    init_circuit.x(0)  # Molecule 0, right qubit (T1)
    init_circuit.x(2 * (N - 1))  # Molecule N-1, right qubit (T1)

    # Build single Trotter step
    trotter_step = build_qubit_trotter_step_exact(params, dt)

    # Initial populations
    state = Statevector(init_circuit)
    pop_0 = calculate_populations_from_statevector(state, N)

    times = [0.0]
    populations = [pop_0]

    # Time evolution
    for step in range(1, N_steps + 1):
        circuit = init_circuit.copy()
        for _ in range(step):
            circuit.compose(trotter_step, inplace=True)

        state = Statevector(circuit)
        pop = calculate_populations_from_statevector(state, N)

        t = step * dt
        times.append(t)
        populations.append(pop)

        if step == N_steps:
            pass

    return {"times": times, "populations": populations}


if __name__ == "__main__":
    params = PhysicalParameters()
    T_total = 100.0
    N_steps = 20

    result = run_qubit_simulation_exact(params, T_total, N_steps)

    final = result["populations"][-1]
