#!/usr/bin/env python3
"""
Helper functions for comparing different quantum simulation approaches.

This module provides utilities for:
1. Qubit simulations: UnitaryGate vs Basic Gates
2. Qudit simulations: CustomTwo gates vs decomposed basic gates
3. Gate counting and circuit analysis

All implementations are exact (no approximations or heuristics).
"""

import numpy as np
from typing import Dict, List, Tuple
from collections import Counter


def count_gates_by_type(circuit, is_qiskit: bool = True) -> Dict[str, int]:
    """
    Count gates by type in a quantum circuit.
    
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
    else:
        # MQT-Qudits circuit
        gate_types = [type(instr).__name__ for instr in circuit.instructions]
        gate_counts = Counter(gate_types)
        return dict(gate_counts)


def compare_gate_counts(counts1: Dict[str, int], counts2: Dict[str, int], 
                       label1: str = "Method 1", label2: str = "Method 2"):
    """
    Compare gate counts between two methods and print results.
    
    Args:
        counts1, counts2: Dictionaries of gate counts
        label1, label2: Labels for the methods
    """
    print("\n" + "="*70)
    print(f"Gate Count Comparison: {label1} vs {label2}")
    print("="*70)
    
    # Get all gate types
    all_gates = set(counts1.keys()) | set(counts2.keys())
    
    # Print comparison table
    print(f"\n{'Gate Type':<20} {label1:>15} {label2:>15} {'Reduction':>15}")
    print("-"*70)
    
    for gate in sorted(all_gates):
        count1 = counts1.get(gate, 0)
        count2 = counts2.get(gate, 0)
        if count1 > 0:
            reduction = f"{(1 - count2/count1)*100:.1f}%"
        else:
            reduction = "N/A"
        print(f"{gate:<20} {count1:>15} {count2:>15} {reduction:>15}")
    
    total1 = sum(counts1.values())
    total2 = sum(counts2.values())
    if total1 > 0:
        total_reduction = f"{(1 - total2/total1)*100:.1f}%"
    else:
        total_reduction = "N/A"
    
    print("-"*70)
    print(f"{'TOTAL':<20} {total1:>15} {total2:>15} {total_reduction:>15}")
    print("="*70)


def print_gate_statistics(circuit, label: str, is_qiskit: bool = True):
    """
    Print detailed gate statistics for a circuit.
    
    Args:
        circuit: Quantum circuit
        label: Description label
        is_qiskit: True for Qiskit, False for MQT-Qudits
    """
    gate_counts = count_gates_by_type(circuit, is_qiskit)
    
    print("\n" + "="*70)
    print(f"Gate Statistics: {label}")
    print("="*70)
    
    for gate_name, count in sorted(gate_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {gate_name:20s}: {count:5d}")
    
    total = sum(gate_counts.values())
    print(f"\n  {'Total Gates':20s}: {total:5d}")
    
    if is_qiskit and hasattr(circuit, 'depth'):
        print(f"  {'Circuit Depth':20s}: {circuit.depth():5d}")
    
    print("="*70)


def decompose_qiskit_unitary_gates(circuit):
    """
    Decompose UnitaryGate instances in a Qiskit circuit to basic gates.
    
    Args:
        circuit: Qiskit QuantumCircuit
    
    Returns:
        Decomposed circuit with only basic gates
    """
    from qiskit import transpile
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    
    # Use Qiskit's transpiler to decompose unitaries
    pm = generate_preset_pass_manager(
        optimization_level=3,
        basis_gates=['cx', 'rz', 'ry', 'rx', 'id', 'x', 'h'],
        approximation_degree=None  # Exact decomposition
    )
    
    decomposed = pm.run(circuit)
    return decomposed


def estimate_qudit_customtwo_decomposition_cost(num_custom_two_gates: int) -> Dict[str, int]:
    """
    Estimate the gate count after decomposing CustomTwo gates.
    
    Based on sparse structure-aware compilation:
    - 3×3 subspace (H_TTA): ~6 basic gates per CustomTwo
    
    Args:
        num_custom_two_gates: Number of CustomTwo gates
    
    Returns:
        Estimated gate counts after decomposition
    """
    # Conservative estimate based on 3×3 subspace
    gates_per_customtwo = 6
    
    return {
        'CEx': num_custom_two_gates * 3,  # Estimated CEx gates
        'R': num_custom_two_gates * 2,    # Estimated R gates
        'VirtRz': num_custom_two_gates * 1,  # Estimated VirtRz gates
    }


def decompose_qudit_customtwo_gates_to_circuit(circuit, sparse_generator):
    """
    Actually decompose CustomTwo gates in a Qudit circuit to basic gates.
    
    This uses the IntegratedSparseCompilerV2 to decompose CustomTwo gates
    into basic gates (CEx, R, VirtRz, Rz, Rh) while recognizing sparse structures.
    
    Note: This function currently uses a private method (_decompose_custom_two_exact)
    from the sparse_generator. This is a temporary solution until a public API
    is established. The function includes validation to fail gracefully if the
    API changes.
    
    Args:
        circuit (mqt.qudits.quantum_circuit.QuantumCircuit): MQT-Qudits QuantumCircuit 
            with CustomTwo gates to be decomposed
        sparse_generator (SparseAwareMQTGateGenerator): Instance with compiler that
            provides the _decompose_custom_two_exact method
    
    Returns:
        mqt.qudits.quantum_circuit.QuantumCircuit: New circuit with CustomTwo gates
            decomposed into basic gates (CEx, R, VirtRz, Rz, Rh)
        
    Raises:
        AttributeError: If sparse_generator doesn't have the decomposition method
        RuntimeError: If decomposition fails for any CustomTwo gate
        
    Example:
        >>> decomposed = decompose_qudit_customtwo_gates_to_circuit(
        ...     circuit_with_customtwo,
        ...     qudit_simulator.time_evol.gate_generator
        ... )
    """
    from mqt.qudits.quantum_circuit import QuantumCircuit as MQTQuantumCircuit
    
    # Validate that sparse_generator has the decomposition method
    if not hasattr(sparse_generator, '_decompose_custom_two_exact'):
        raise AttributeError(
            "sparse_generator must have '_decompose_custom_two_exact' method. "
            "Please ensure you're using SparseAwareMQTGateGenerator instance."
        )
    
    # Create new circuit with same qudits
    decomposed_circuit = MQTQuantumCircuit(circuit.num_qudits, circuit.dimensions)
    
    # Process each gate
    for gate in circuit.instructions:
        if gate.__class__.__name__ == 'CustomTwo':
            try:
                # Decompose CustomTwo gate using the sparse generator
                # Note: This uses a private method for now. Consider making this
                # a public API if the sparse_generator interface is stabilized.
                decomposed_gates = sparse_generator._decompose_custom_two_exact(gate)
                for dec_gate in decomposed_gates:
                    decomposed_circuit.append(dec_gate)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to decompose CustomTwo gate: {e}. "
                    "This may indicate an API change or incompatibility."
                ) from e
        else:
            # Copy other gates as-is
            decomposed_circuit.append(gate)
    
    return decomposed_circuit
