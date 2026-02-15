#!/usr/bin/env python3
"""Standalone exact qubit simulator - complete implementation from scratch.

This is a reference implementation to validate our exact Hamiltonians work correctly.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from exact_qubit_hamiltonians import build_H_transfer_qubit_unitary, build_H_TTA_qubit_unitary
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector


class PhysicalParameters:
    def __init__(self) -> None:
        self.N_molecules = 4
        self.E_T = 1.5
        self.E_S = 3.0
        self.V = 0.1
        self.J = 0.05
        self.hbar = 0.6582119569
        self.neighbors = [(0, 1), (1, 2), (2, 3)]


class QubitSimulatorExact:
    def __init__(self, params) -> None:
        self.params = params
        self.N = params.N_molecules
        self.n_qubits = 2 * self.N

    def apply_H0_evolution(self, circuit, mol_idx, dt) -> None:
        """Apply H0 evolution to a single molecule."""
        q0 = 2 * mol_idx
        q1 = 2 * mol_idx + 1

        E_T = self.params.E_T
        E_S = self.params.E_S
        hbar = self.params.hbar

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

    def apply_H_transfer_exact(self, circuit, mol_i, mol_j, dt) -> None:
        """Apply exact H_transfer evolution."""
        V = self.params.V
        hbar = self.params.hbar

        # Build exact unitary
        U_tr = build_H_transfer_qubit_unitary(V, dt, hbar)

        # Apply to qubits
        qi0, qi1 = 2 * mol_i, 2 * mol_i + 1
        qj0, qj1 = 2 * mol_j, 2 * mol_j + 1
        qubits = [qi0, qi1, qj0, qj1]

        gate = UnitaryGate(U_tr, label="U_tr")
        circuit.append(gate, qubits)

    def apply_H_TTA_exact(self, circuit, mol_i, mol_j, dt) -> None:
        """Apply exact H_TTA evolution."""
        J = self.params.J
        hbar = self.params.hbar

        # Build exact unitary
        U_tta = build_H_TTA_qubit_unitary(J, dt, hbar)

        # Apply to qubits
        qi0, qi1 = 2 * mol_i, 2 * mol_i + 1
        qj0, qj1 = 2 * mol_j, 2 * mol_j + 1
        qubits = [qi0, qi1, qj0, qj1]

        gate = UnitaryGate(U_tta, label="U_TTA")
        circuit.append(gate, qubits)

    def build_trotter_step(self, dt):
        """Build one Trotter step."""
        circuit = QuantumCircuit(self.n_qubits)

        # Forward
        for i in range(self.N):
            self.apply_H0_evolution(circuit, i, dt / 2)

        for i, j in self.params.neighbors:
            self.apply_H_transfer_exact(circuit, i, j, dt / 2)

        for i, j in self.params.neighbors:
            self.apply_H_TTA_exact(circuit, i, j, dt / 2)

        # Backward
        for i, j in reversed(self.params.neighbors):
            self.apply_H_TTA_exact(circuit, i, j, dt / 2)

        for i, j in reversed(self.params.neighbors):
            self.apply_H_transfer_exact(circuit, i, j, dt / 2)

        for i in reversed(range(self.N)):
            self.apply_H0_evolution(circuit, i, dt / 2)

        return circuit

    def calculate_populations(self, statevector):
        """Calculate populations from statevector."""
        probabilities = statevector.probabilities_dict()

        N_S0 = N_T1 = N_S1 = 0.0

        for bitstring, prob in probabilities.items():
            if prob < 1e-15:
                continue

            bits = bitstring[::-1]  # Convert to little-endian

            for mol in range(self.N):
                q0_bit = int(bits[2 * mol])
                q1_bit = int(bits[2 * mol + 1])

                if q1_bit == 0 and q0_bit == 0:
                    N_S0 += prob
                elif q1_bit == 0 and q0_bit == 1:
                    N_T1 += prob
                elif q1_bit == 1 and q0_bit == 0:
                    N_S1 += prob

        return {"N_S0": N_S0, "N_T1": N_T1, "N_S1": N_S1}

    def simulate(self, T_total, N_steps):
        """Run exact statevector simulation."""
        dt = T_total / N_steps

        # Prepare initial state: edge_triplet
        init_circuit = QuantumCircuit(self.n_qubits)
        init_circuit.x(0)  # Molecule 0 to T1
        init_circuit.x(2 * (self.N - 1))  # Molecule N-1 to T1

        # Build Trotter step
        trotter_step = self.build_trotter_step(dt)

        # Initial populations
        state = Statevector(init_circuit)
        self.calculate_populations(state)

        # Time evolution
        for step in range(1, N_steps + 1):
            circuit = init_circuit.copy()
            for _ in range(step):
                circuit.compose(trotter_step, inplace=True)

            state = Statevector(circuit)

        return self.calculate_populations(state)


if __name__ == "__main__":
    params = PhysicalParameters()
    sim = QubitSimulatorExact(params)
    result = sim.simulate(T_total=100.0, N_steps=20)
