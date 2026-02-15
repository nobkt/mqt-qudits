#!/usr/bin/env python3
"""Test for the MQT-Qudits gate-based 4-molecule implementation."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "tutorials"))

import numpy as np
from mqt_qudits_four_molecule_implementation import (
    MQTQuditTimeEvolution,
    PhysicalParameters,
    SuzukiTrotterMQTQuditSimulator,
)


def test_parameters():
    """Test parameter initialization."""
    params = PhysicalParameters()
    assert params.N_molecules == 4
    assert params.E_T == 1.5
    assert params.E_S == 3.0
    assert len(params.neighbors) == 3
    print("✓ Parameters test passed")


def test_gate_construction():
    """Test quantum gate construction."""
    params = PhysicalParameters()
    time_evol = MQTQuditTimeEvolution(params)

    from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister

    circuit = QuantumCircuit()
    reg = QuantumRegister("molecules", 4, [3, 3, 3, 3])
    circuit.append(reg)

    dt = 1.0

    # Test H0 gates
    initial_count = len(circuit.instructions)
    time_evol.add_H0_evolution_gates(circuit, dt)
    assert len(circuit.instructions) == initial_count + 8  # 4 qudits × 2 levels

    # Test H_transfer gates
    initial_count = len(circuit.instructions)
    time_evol.add_H_transfer_evolution_gates(circuit, dt)
    assert len(circuit.instructions) == initial_count + 3  # 3 neighbor pairs

    # Test H_TTA gates
    initial_count = len(circuit.instructions)
    time_evol.add_H_TTA_evolution_gates(circuit, dt)
    assert len(circuit.instructions) == initial_count + 3  # 3 neighbor pairs

    print("✓ Gate construction test passed")


def test_initial_state():
    """Test initial state preparation."""
    params = PhysicalParameters()
    simulator = SuzukiTrotterMQTQuditSimulator(params)

    # Test all_triplet state
    circuit = simulator.build_initial_state_circuit("all_triplet")
    job = simulator.backend.run(circuit)
    result = job.result()
    state = result.get_state_vector().flatten()

    # Should be in state |1111⟩ which is index 40
    non_zero = np.where(np.abs(state) > 0.1)[0]
    assert len(non_zero) == 1
    assert non_zero[0] == 40  # |1111⟩ = 1*27 + 1*9 + 1*3 + 1*1 = 40

    print("✓ Initial state test passed")


def test_short_simulation():
    """Test a short simulation run."""
    params = PhysicalParameters()
    simulator = SuzukiTrotterMQTQuditSimulator(params)

    results = simulator.simulate(T_total=20.0, N_steps=4, initial_state_type="all_triplet", track_dynamics=True)

    # Check initial state
    init_pops = results["populations"][0]
    assert abs(init_pops["N_S0"] - 0.0) < 1e-6
    assert abs(init_pops["N_T1"] - 4.0) < 1e-6
    assert abs(init_pops["N_S1"] - 0.0) < 1e-6

    # Check that dynamics occurred
    final_pops = results["populations"][-1]
    assert final_pops["N_T1"] < init_pops["N_T1"]  # Triplets should decrease
    assert final_pops["N_S1"] > init_pops["N_S1"]  # Singlets should increase
    assert final_pops["N_S0"] > init_pops["N_S0"]  # Ground state should increase

    print("✓ Short simulation test passed")


def test_population_conservation():
    """Test that total population is approximately conserved."""
    params = PhysicalParameters()
    params.Gamma_fl = 0.0  # Turn off decay for this test
    simulator = SuzukiTrotterMQTQuditSimulator(params)

    results = simulator.simulate(T_total=10.0, N_steps=5, initial_state_type="all_triplet", track_dynamics=True)

    # Check conservation
    for pops in results["populations"]:
        total = sum(pops.values())
        assert abs(total - 4.0) < 1e-6, f"Total population not conserved: {total}"

    print("✓ Population conservation test passed")


if __name__ == "__main__":
    print("Running MQT-Qudits Implementation Tests\n")
    print("=" * 60)

    test_parameters()
    test_gate_construction()
    test_initial_state()
    test_short_simulation()
    test_population_conservation()

    print("=" * 60)
    print("\nAll tests passed! ✓")
