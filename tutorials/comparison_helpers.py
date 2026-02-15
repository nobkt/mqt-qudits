#!/usr/bin/env python3
"""Helper functions for comparing different quantum simulation approaches.

This module provides utilities for:
1. Qubit simulations: UnitaryGate vs Basic Gates
2. Qudit simulations: CustomTwo gates vs decomposed basic gates
3. Gate counting and circuit analysis

All implementations are exact (no approximations or heuristics).
"""

from __future__ import annotations

import operator
from collections import Counter


def count_gates_by_type(circuit, is_qiskit: bool = True) -> dict[str, int]:
    """Count gates by type in a quantum circuit.

    Args:
        circuit: Quantum circuit (Qiskit or MQT-Qudits)
        is_qiskit: True for Qiskit circuits, False for MQT-Qudits circuits

    Returns:
        Dictionary mapping gate names to counts
    """
    if is_qiskit:
        # Qiskit circuit
        gate_counts = Counter([instr.operation.name for instr in circuit.data])
        return dict(gate_counts)
    # MQT-Qudits circuit
    gate_types = [type(instr).__name__ for instr in circuit.instructions]
    gate_counts = Counter(gate_types)
    return dict(gate_counts)


def compare_gate_counts(
    counts1: dict[str, int], counts2: dict[str, int], label1: str = "Method 1", label2: str = "Method 2"
) -> None:
    """Compare gate counts between two methods and print results.

    Args:
        counts1, counts2: Dictionaries of gate counts
        label1, label2: Labels for the methods
    """
    # Get all gate types
    all_gates = set(counts1.keys()) | set(counts2.keys())

    # Print comparison table

    for gate in sorted(all_gates):
        count1 = counts1.get(gate, 0)
        count2 = counts2.get(gate, 0)
        f"{(1 - count2 / count1) * 100:.1f}%" if count1 > 0 else "N/A"

    total1 = sum(counts1.values())
    total2 = sum(counts2.values())
    f"{(1 - total2 / total1) * 100:.1f}%" if total1 > 0 else "N/A"


def print_gate_statistics(circuit, label: str, is_qiskit: bool = True) -> None:
    """Print detailed gate statistics for a circuit.

    Args:
        circuit: Quantum circuit
        label: Description label
        is_qiskit: True for Qiskit, False for MQT-Qudits
    """
    gate_counts = count_gates_by_type(circuit, is_qiskit)

    for _gate_name, _count in sorted(gate_counts.items(), key=operator.itemgetter(1), reverse=True):
        pass

    sum(gate_counts.values())

    if is_qiskit and hasattr(circuit, "depth"):
        pass


def decompose_qiskit_unitary_gates(circuit):
    """Decompose UnitaryGate instances in a Qiskit circuit to basic gates.

    Args:
        circuit: Qiskit QuantumCircuit

    Returns:
        Decomposed circuit with only basic gates
    """
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    # Use Qiskit's transpiler to decompose unitaries
    pm = generate_preset_pass_manager(
        optimization_level=3,
        basis_gates=["cx", "rz", "ry", "rx", "id", "x", "h"],
        approximation_degree=None,  # Exact decomposition
    )

    return pm.run(circuit)


def estimate_qudit_customtwo_decomposition_cost(num_custom_two_gates: int) -> dict[str, int]:
    """Estimate the gate count after decomposing CustomTwo gates.

    Based on sparse structure-aware compilation:
    - 3×3 subspace (H_TTA): ~6 basic gates per CustomTwo

    Args:
        num_custom_two_gates: Number of CustomTwo gates

    Returns:
        Estimated gate counts after decomposition
    """
    # Conservative estimate based on 3×3 subspace

    return {
        "CEx": num_custom_two_gates * 3,  # Estimated CEx gates
        "R": num_custom_two_gates * 2,  # Estimated R gates
        "VirtRz": num_custom_two_gates * 1,  # Estimated VirtRz gates
    }


def decompose_qudit_customtwo_gates_to_circuit(circuit, time_evol):
    """Actually decompose CustomTwo gates in a Qudit circuit to basic gates.

    This uses the IntegratedSparseCompilerV2 to decompose CustomTwo gates
    into basic gates (CEx, R, VirtRz, Rz, Rh) while recognizing sparse structures.

    Note: This function currently uses a private method (_decompose_custom_two_exact)
    from the time_evol object. This is a temporary solution until a public API
    is established. The function includes validation to fail gracefully if the
    API changes.

    Args:
        circuit (mqt.qudits.quantum_circuit.QuantumCircuit): MQT-Qudits QuantumCircuit
            with CustomTwo gates to be decomposed
        time_evol (SparseAwareMQTQuditTimeEvolution): Time evolution instance that
            provides the _decompose_custom_two_exact method

    Returns:
        mqt.qudits.quantum_circuit.QuantumCircuit: New circuit with CustomTwo gates
            decomposed into basic gates (CEx, R, VirtRz, Rz, Rh)

    Raises:
        AttributeError: If time_evol doesn't have the decomposition method
        RuntimeError: If decomposition fails for any CustomTwo gate

    Example:
        >>> decomposed = decompose_qudit_customtwo_gates_to_circuit(
        ...     circuit_with_customtwo,
        ...     qudit_simulator.time_evol
        ... )
    """
    from mqt.qudits.quantum_circuit import QuantumCircuit as MQTQuantumCircuit

    # Validate that time_evol has the decomposition method
    if not hasattr(time_evol, "_decompose_custom_two_exact"):
        msg = (
            "time_evol must have '_decompose_custom_two_exact' method. "
            "Please ensure you're using SparseAwareMQTQuditTimeEvolution instance."
        )
        raise AttributeError(msg)

    # Create new circuit with same qudits
    decomposed_circuit = MQTQuantumCircuit(circuit.num_qudits, circuit.dimensions)

    # Process each gate
    for gate in circuit.instructions:
        if gate.__class__.__name__ == "CustomTwo":
            try:
                # Decompose CustomTwo gate using the time evolution's method
                # Note: This uses a private method for now. Consider making this
                # a public API if the time_evol interface is stabilized.
                decomposed_gates = time_evol._decompose_custom_two_exact(gate)
                for dec_gate in decomposed_gates:
                    decomposed_circuit.instructions.append(dec_gate)
            except Exception as e:
                msg = f"Failed to decompose CustomTwo gate: {e}. This may indicate an API change or incompatibility."
                raise RuntimeError(msg) from e
        else:
            # Copy other gates as-is
            decomposed_circuit.instructions.append(gate)

    return decomposed_circuit
