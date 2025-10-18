#!/usr/bin/env python3
"""
Quantum Circuit Visualization Tool for MQT-Qudits

This script provides comprehensive visualization for quantum circuits,
with special support for CustomTwo gates and their decomposition.

References src/ code without modifying it:
- mqt.qudits.quantum_circuit for circuit structures
- mqt.qudits.quantum_circuit.gates for gate types
- mqt.qudits.quantum_circuit.components.extensions.gate_types for GateTypes enum
"""

from __future__ import annotations
from typing import TYPE_CHECKING
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

if TYPE_CHECKING:
    from mqt.qudits.quantum_circuit import QuantumCircuit

# Import gate types from src/
from mqt.qudits.quantum_circuit import gates
from mqt.qudits.quantum_circuit.components.extensions.gate_types import GateTypes


class CircuitVisualizer:
    """
    Visualize quantum circuits with proper handling of all gate types,
    including CustomTwo gates and their decomposition to basic gates.
    """
    
    def __init__(self, circuit: QuantumCircuit):
        """
        Initialize the visualizer with a quantum circuit.
        
        Args:
            circuit: The QuantumCircuit to visualize
        """
        self.circuit = circuit
        self.num_qudits = circuit.num_qudits
        self.instructions = circuit.instructions
        
        # Visual styling parameters
        self.qudit_spacing = 1.0  # Vertical spacing between qudit lines
        self.gate_width = 0.6
        self.gate_height = 0.4
        self.wire_color = 'black'
        self.gate_colors = {
            'VirtRz': '#FFE5B4',      # Peach for single-qudit phase gates
            'R': '#B4D7FF',            # Light blue for rotation gates
            'Rh': '#B4FFD7',           # Light green for hermitian rotation
            'Rz': '#D7B4FF',           # Light purple for Z rotation
            'CEx': '#FFB4D7',          # Pink for controlled exchange
            'CustomTwo': '#FFD700',    # Gold for custom two-qudit gates
            'CustomOne': '#FFA500',    # Orange for custom single-qudit gates
            'X': '#87CEEB',            # Sky blue for X gates
            'default': '#E0E0E0'       # Gray for other gates
        }
    
    def _get_gate_label(self, gate) -> str:
        """
        Generate a descriptive label for a gate.
        
        Args:
            gate: The gate instruction
            
        Returns:
            String label for the gate
        """
        gate_type = type(gate).__name__
        
        if isinstance(gate, gates.VirtRz):
            return f"VRz_{gate.lev_a}\n({gate.phi:.3f})"
        elif isinstance(gate, gates.R):
            return f"R_{gate.lev_a}{gate.lev_b}\n(θ={gate.theta:.2f},\nφ={gate.phi:.2f})"
        elif isinstance(gate, gates.Rh):
            return f"Rh_{gate.lev_a}{gate.lev_b}\n(θ={gate.theta:.2f})"
        elif isinstance(gate, gates.Rz):
            return f"Rz_{gate.lev_a}{gate.lev_b}\n(θ={gate.theta:.2f})"
        elif isinstance(gate, gates.CEx):
            # CEx is a two-qudit gate
            return "CEx"
        elif isinstance(gate, gates.CustomTwo):
            # CustomTwo gate with dimension info
            dims = gate.dimensions if hasattr(gate, 'dimensions') else 'N/A'
            return f"CustomTwo\n{dims[0]}×{dims[1]}"
        elif isinstance(gate, gates.CustomOne):
            return "CustomOne"
        elif isinstance(gate, gates.X):
            return f"X_{gate.lev_a}{gate.lev_b}"
        else:
            return gate_type
    
    def _get_gate_color(self, gate) -> str:
        """
        Get the color for a gate type.
        
        Args:
            gate: The gate instruction
            
        Returns:
            Color string for the gate
        """
        gate_type = type(gate).__name__
        return self.gate_colors.get(gate_type, self.gate_colors['default'])
    
    def draw_circuit(self, title: str = "Quantum Circuit", figsize: tuple = (16, 10)):
        """
        Draw the quantum circuit with matplotlib.
        
        Args:
            title: Title for the circuit diagram
            figsize: Figure size as (width, height)
            
        Returns:
            matplotlib Figure and Axes objects
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Calculate circuit width based on number of gates
        circuit_width = max(len(self.instructions) * 1.2 + 2, 10)
        
        # Set up the axes
        ax.set_xlim(0, circuit_width)
        ax.set_ylim(-0.5, self.num_qudits - 0.5)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw qudit lines
        for qudit_idx in range(self.num_qudits):
            y = self.num_qudits - 1 - qudit_idx  # Invert to show qudit 0 at top
            ax.plot([0.5, circuit_width - 0.5], [y, y], 
                   color=self.wire_color, linewidth=1.5, zorder=1)
            # Label qudit lines
            ax.text(0.2, y, f'q{qudit_idx}', ha='right', va='center', 
                   fontsize=12, fontweight='bold')
            ax.text(0.3, y, '|0⟩', ha='left', va='center', fontsize=10)
        
        # Draw gates
        gate_x_position = 1.0
        x_positions = []
        
        for gate_idx, gate in enumerate(self.instructions):
            gate_type_enum = gate.gate_type
            targets = gate.target_qudits
            
            # Ensure targets is a list
            if isinstance(targets, int):
                targets = [targets]
            
            # Convert target indices to y-coordinates (inverted)
            target_ys = [self.num_qudits - 1 - t for t in targets]
            
            if gate_type_enum == GateTypes.SINGLE:
                # Single qudit gate
                self._draw_single_gate(ax, gate, gate_x_position, target_ys[0])
            elif gate_type_enum == GateTypes.TWO:
                # Two qudit gate
                self._draw_two_gate(ax, gate, gate_x_position, target_ys)
            elif gate_type_enum == GateTypes.MULTI:
                # Multi qudit gate
                self._draw_multi_gate(ax, gate, gate_x_position, target_ys)
            
            x_positions.append(gate_x_position)
            gate_x_position += 1.2  # Move to next gate position
        
        # Add title
        ax.text(circuit_width / 2, self.num_qudits + 0.3, title,
               ha='center', va='bottom', fontsize=16, fontweight='bold')
        
        # Add gate count information
        gate_counts = self._count_gates()
        info_text = f"Total gates: {len(self.instructions)}"
        ax.text(circuit_width / 2, -0.8, info_text,
               ha='center', va='top', fontsize=12)
        
        # Add legend
        self._add_legend(ax, circuit_width)
        
        plt.tight_layout()
        return fig, ax
    
    def _draw_single_gate(self, ax, gate, x: float, y: float):
        """
        Draw a single-qudit gate.
        
        Args:
            ax: Matplotlib axes
            gate: Gate instruction
            x: X position
            y: Y position
        """
        color = self._get_gate_color(gate)
        label = self._get_gate_label(gate)
        
        # Draw gate box
        rect = mpatches.FancyBboxPatch(
            (x - self.gate_width/2, y - self.gate_height/2),
            self.gate_width, self.gate_height,
            boxstyle="round,pad=0.05",
            facecolor=color,
            edgecolor='black',
            linewidth=1.5,
            zorder=3
        )
        ax.add_patch(rect)
        
        # Add gate label
        ax.text(x, y, label, ha='center', va='center',
               fontsize=8, fontweight='bold', zorder=4)
    
    def _draw_two_gate(self, ax, gate, x: float, target_ys: list):
        """
        Draw a two-qudit gate.
        
        Args:
            ax: Matplotlib axes
            gate: Gate instruction
            x: X position
            target_ys: List of Y positions for target qudits
        """
        color = self._get_gate_color(gate)
        label = self._get_gate_label(gate)
        
        # Sort target positions
        y_min = min(target_ys)
        y_max = max(target_ys)
        
        # Draw connecting line between qudits
        ax.plot([x, x], [y_min, y_max], 
               color='black', linewidth=2, zorder=2)
        
        # Draw gate boxes on each qudit
        for y in target_ys:
            rect = mpatches.FancyBboxPatch(
                (x - self.gate_width/2, y - self.gate_height/2),
                self.gate_width, self.gate_height,
                boxstyle="round,pad=0.05",
                facecolor=color,
                edgecolor='black',
                linewidth=1.5,
                zorder=3
            )
            ax.add_patch(rect)
        
        # Add label in the middle of the gate
        mid_y = (y_min + y_max) / 2
        
        # For CustomTwo gates, show label near the connecting line
        if isinstance(gate, gates.CustomTwo):
            ax.text(x + 0.4, mid_y, label, ha='left', va='center',
                   fontsize=8, fontweight='bold', zorder=4,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor=color, 
                            edgecolor='black', alpha=0.9))
        else:
            # For other two-qudit gates, show label in boxes
            for y in target_ys:
                display_label = label if len(label) < 10 else label.split('\n')[0]
                ax.text(x, y, display_label, ha='center', va='center',
                       fontsize=7, fontweight='bold', zorder=4)
    
    def _draw_multi_gate(self, ax, gate, x: float, target_ys: list):
        """
        Draw a multi-qudit gate.
        
        Args:
            ax: Matplotlib axes
            gate: Gate instruction
            x: X position
            target_ys: List of Y positions for target qudits
        """
        color = self._get_gate_color(gate)
        label = self._get_gate_label(gate)
        
        # Sort target positions
        y_min = min(target_ys)
        y_max = max(target_ys)
        
        # Draw connecting line between qudits
        ax.plot([x, x], [y_min, y_max], 
               color='black', linewidth=2, zorder=2)
        
        # Draw a single large box spanning all qudits
        height = y_max - y_min + self.gate_height
        rect = mpatches.FancyBboxPatch(
            (x - self.gate_width/2, y_min - self.gate_height/2),
            self.gate_width, height,
            boxstyle="round,pad=0.05",
            facecolor=color,
            edgecolor='black',
            linewidth=1.5,
            zorder=3
        )
        ax.add_patch(rect)
        
        # Add label in the middle
        mid_y = (y_min + y_max) / 2
        ax.text(x, mid_y, label, ha='center', va='center',
               fontsize=8, fontweight='bold', zorder=4)
    
    def _count_gates(self) -> dict:
        """
        Count gates by type.
        
        Returns:
            Dictionary mapping gate type names to counts
        """
        counts = {}
        for gate in self.instructions:
            gate_type = type(gate).__name__
            counts[gate_type] = counts.get(gate_type, 0) + 1
        return counts
    
    def _add_legend(self, ax, circuit_width: float):
        """
        Add a legend showing gate types and colors.
        
        Args:
            ax: Matplotlib axes
            circuit_width: Width of the circuit
        """
        gate_counts = self._count_gates()
        if not gate_counts:
            return
        
        # Create legend entries for gates that appear in the circuit
        legend_elements = []
        for gate_type, count in sorted(gate_counts.items()):
            color = self.gate_colors.get(gate_type, self.gate_colors['default'])
            legend_elements.append(
                mpatches.Patch(facecolor=color, edgecolor='black', 
                             label=f'{gate_type}: {count}')
            )
        
        # Place legend below the circuit
        ax.legend(handles=legend_elements, loc='upper right', 
                 bbox_to_anchor=(1.0, -0.05), ncol=3, fontsize=10,
                 framealpha=0.9)
    
    def print_circuit_summary(self):
        """
        Print a text summary of the circuit.
        """
        print("=" * 70)
        print("CIRCUIT SUMMARY")
        print("=" * 70)
        print(f"Number of qudits: {self.num_qudits}")
        print(f"Total gates: {len(self.instructions)}")
        print()
        
        print("Gate composition:")
        gate_counts = self._count_gates()
        for gate_type, count in sorted(gate_counts.items()):
            print(f"  {gate_type:20s}: {count:4d}")
        print()
        
        print("Gate sequence:")
        for i, gate in enumerate(self.instructions):
            gate_type = type(gate).__name__
            targets = gate.target_qudits
            if isinstance(targets, int):
                targets = [targets]
            print(f"  {i:3d}. {gate_type:15s} on qudits {targets}")
        print("=" * 70)


def visualize_circuit_with_decomposition(
    circuit_before: QuantumCircuit,
    circuit_after: QuantumCircuit,
    title_before: str = "Circuit with CustomTwo Gates",
    title_after: str = "Circuit with Decomposed Basic Gates",
    save_path: str | None = None
) -> tuple:
    """
    Visualize a circuit before and after CustomTwo gate decomposition.
    
    Args:
        circuit_before: Circuit with CustomTwo gates
        circuit_after: Circuit after decomposing CustomTwo gates
        title_before: Title for the first circuit
        title_after: Title for the second circuit
        save_path: Optional path to save the figure
        
    Returns:
        Tuple of (fig1, ax1, fig2, ax2) for both visualizations
    """
    # Visualize circuit before decomposition
    viz_before = CircuitVisualizer(circuit_before)
    print("\n" + "=" * 70)
    print("BEFORE DECOMPOSITION")
    viz_before.print_circuit_summary()
    fig1, ax1 = viz_before.draw_circuit(title=title_before)
    
    # Visualize circuit after decomposition
    viz_after = CircuitVisualizer(circuit_after)
    print("\n" + "=" * 70)
    print("AFTER DECOMPOSITION")
    viz_after.print_circuit_summary()
    fig2, ax2 = viz_after.draw_circuit(title=title_after)
    
    if save_path:
        fig1.savefig(f"{save_path}_before.png", dpi=150, bbox_inches='tight')
        fig2.savefig(f"{save_path}_after.png", dpi=150, bbox_inches='tight')
        print(f"\nSaved visualizations to {save_path}_before.png and {save_path}_after.png")
    
    return fig1, ax1, fig2, ax2


def visualize_circuit(
    circuit: QuantumCircuit,
    title: str = "Quantum Circuit",
    save_path: str | None = None
) -> tuple:
    """
    Visualize a single quantum circuit.
    
    Args:
        circuit: The quantum circuit to visualize
        title: Title for the circuit diagram
        save_path: Optional path to save the figure
        
    Returns:
        Tuple of (fig, ax) for the visualization
    """
    viz = CircuitVisualizer(circuit)
    viz.print_circuit_summary()
    fig, ax = viz.draw_circuit(title=title)
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved visualization to {save_path}")
    
    return fig, ax


# Example usage (when run as a script)
if __name__ == "__main__":
    print("Circuit Visualization Tool for MQT-Qudits")
    print("==========================================")
    print()
    print("This module provides functions to visualize quantum circuits.")
    print("Import this module in your notebook or script:")
    print()
    print("  from tools.visualize_circuit import visualize_circuit_with_decomposition")
    print("  visualize_circuit_with_decomposition(circuit_before, circuit_after)")
    print()
    print("Or for a single circuit:")
    print()
    print("  from tools.visualize_circuit import visualize_circuit")
    print("  visualize_circuit(circuit)")
