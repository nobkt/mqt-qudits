from __future__ import annotations

from ..quantum_circuit import QuantumCircuit, gates
from ..quantum_circuit.components.extensions.gate_types import GateTypes


def draw_qudit_local(circuit: QuantumCircuit, max_gates: int | None = None) -> None:
    """
    Draw a text-based representation of a quantum circuit.
    
    Args:
        circuit: The quantum circuit to visualize
        max_gates: Maximum number of gates to display (None for all gates)
    """
    instructions = circuit.instructions[:max_gates] if max_gates else circuit.instructions
    
    for line in range(circuit.num_qudits):
        print(f"q{line}|0>---", end="")
        for gate in instructions:
            if gate.gate_type == GateTypes.SINGLE and line == gate.target_qudits:
                if isinstance(gate, gates.VirtRz):
                    print(f"--[VRz{gate.lev_a}({gate.phi:.2f})]--", end="")
                elif isinstance(gate, gates.R):
                    print(f"--[R{gate.lev_a}{gate.lev_b}({gate.theta:.2f},{gate.phi:.2f})]--", end="")
                elif isinstance(gate, gates.Rh):
                    print(f"--[Rh{gate.lev_a}{gate.lev_b}({gate.theta:.2f})]--", end="")
                elif isinstance(gate, gates.Rz):
                    print(f"--[Rz{gate.lev_a}{gate.lev_b}({gate.phi:.2f})]--", end="")
                elif isinstance(gate, gates.X):
                    print("--[X]--", end="")
                elif isinstance(gate, gates.CustomOne):
                    print("--[CuOne]--", end="")
                else:
                    print(f"--[{gate.__class__.__name__}]--", end="")
            elif gate.gate_type == GateTypes.TWO:
                # Two-qudit gates
                # target_qudits can be either int or list[int], normalize to list
                if isinstance(gate.target_qudits, list):
                    target_qudits = gate.target_qudits
                else:
                    target_qudits = [gate.target_qudits]
                if line in target_qudits:
                    if isinstance(gate, gates.CustomTwo):
                        print("--[CustomTwo]--", end="")
                    elif isinstance(gate, gates.CEx):
                        print("--[CEx]--", end="")
                    else:
                        print(f"--[{gate.__class__.__name__}]--", end="")
                else:
                    print("----", end="")
            else:
                # Gate on different qudit, draw empty wire
                print("----", end="")
        if max_gates and len(circuit.instructions) > max_gates:
            print("--...--", end="")
        print("---=||")


def draw_circuit_summary(circuit: QuantumCircuit) -> None:
    """
    Print a summary of the quantum circuit.
    
    Args:
        circuit: The quantum circuit to summarize
    """
    print(f"Circuit Summary:")
    print(f"  Number of qudits: {circuit.num_qudits}")
    print(f"  Qudit dimensions: {circuit.dimensions}")
    print(f"  Total gates: {len(circuit.instructions)}")
    
    # Count gate types
    gate_counts = {}
    for gate in circuit.instructions:
        gate_name = gate.__class__.__name__
        gate_counts[gate_name] = gate_counts.get(gate_name, 0) + 1
    
    print(f"  Gate composition:")
    for gate_name, count in sorted(gate_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"    {gate_name}: {count}")
