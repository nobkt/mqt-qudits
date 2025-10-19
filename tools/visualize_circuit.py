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
    
    Supports automatic circuit splitting into multiple rows when the circuit
    becomes too large to display in a single figure.
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
        self.gate_width = 0.6  # Gate width in plot units
        self.gate_height = 0.4  # Gate height in plot units
        self.gate_spacing = 1.2  # Horizontal spacing between gates in plot units
        self.wire_color = 'black'
        
        # Size constraints for proper visualization
        self.min_gate_width_inches = 0.15  # Minimum gate width in inches for readability
        self.max_figure_width_inches = 30  # Maximum figure width before wrapping
        self.row_height_inches = 2.5  # Height per qudit in inches (for multi-row layouts)
        
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
    
    def _calculate_layout(self) -> dict:
        """
        Calculate optimal layout for visualizing the circuit.
        
        Determines whether to use single-row or multi-row layout,
        and calculates appropriate figure dimensions.
        
        Returns:
            Dictionary with layout parameters:
            - 'num_rows': Number of rows needed
            - 'gates_per_row': List of gate counts for each row
            - 'figsize': Tuple (width, height) for the figure
            - 'circuit_width_units': Width in plot units per row
        """
        num_gates = len(self.instructions)
        
        if num_gates == 0:
            return {
                'num_rows': 1,
                'gates_per_row': [0],
                'figsize': (10, max(self.num_qudits * self.row_height_inches, 6)),
                'circuit_width_units': 10
            }
        
        # Calculate required width in plot units for all gates
        circuit_width_units = num_gates * self.gate_spacing + 2
        
        # Calculate required width in inches to maintain minimum gate width
        # Each gate occupies gate_spacing units, and we need min_gate_width_inches per gate
        # Formula: required_width = circuit_width_units * (min_gate_width_inches / gate_width)
        required_width_inches = circuit_width_units * (self.min_gate_width_inches / self.gate_width)
        
        # If circuit fits in maximum width, use single row
        if required_width_inches <= self.max_figure_width_inches:
            figsize = (
                max(required_width_inches, 10),  # At least 10 inches
                max(self.num_qudits * self.row_height_inches, 6)  # At least 6 inches
            )
            return {
                'num_rows': 1,
                'gates_per_row': [num_gates],
                'figsize': figsize,
                'circuit_width_units': circuit_width_units
            }
        
        # Circuit is too large - split into multiple rows
        # Calculate gates per row to fit in max width
        # We want: (gates_per_row * gate_spacing + 2) * (min_gate_width_inches / gate_width) = max_figure_width_inches
        # Solving for gates_per_row:
        gates_per_row = int((self.max_figure_width_inches * self.gate_width / self.min_gate_width_inches - 2) / self.gate_spacing)
        
        # Ensure at least 1 gate per row
        gates_per_row = max(1, gates_per_row)
        
        # Calculate number of rows needed
        num_rows = (num_gates + gates_per_row - 1) // gates_per_row
        
        # Distribute gates across rows
        gates_per_row_list = []
        remaining_gates = num_gates
        for _ in range(num_rows):
            gates_in_this_row = min(gates_per_row, remaining_gates)
            gates_per_row_list.append(gates_in_this_row)
            remaining_gates -= gates_in_this_row
        
        # Calculate figure size for multi-row layout
        # Width: based on gates per row (use maximum to keep consistent width)
        max_gates_in_row = max(gates_per_row_list)
        width_per_row = max_gates_in_row * self.gate_spacing + 2
        figsize_width = width_per_row * (self.min_gate_width_inches / self.gate_width)
        
        # Height: number of rows * (qudits + spacing between rows)
        row_spacing = 1.5  # Extra space between rows
        figsize_height = num_rows * self.num_qudits * self.row_height_inches + (num_rows - 1) * row_spacing
        
        return {
            'num_rows': num_rows,
            'gates_per_row': gates_per_row_list,
            'figsize': (figsize_width, figsize_height),
            'circuit_width_units': width_per_row
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
    
    def draw_circuit(self, title: str = "Quantum Circuit", figsize: tuple = None):
        """
        Draw the quantum circuit with matplotlib.
        
        Automatically adjusts layout based on circuit size:
        - Small circuits: single row with appropriate width
        - Large circuits: multiple rows to maintain readability
        
        Args:
            title: Title for the circuit diagram
            figsize: Optional figure size as (width, height). If None, calculated automatically.
            
        Returns:
            matplotlib Figure and Axes objects
        """
        # Calculate optimal layout
        layout = self._calculate_layout()
        
        if figsize is None:
            figsize = layout['figsize']
        
        fig, ax = plt.subplots(figsize=figsize)
        
        num_rows = layout['num_rows']
        gates_per_row = layout['gates_per_row']
        circuit_width = layout['circuit_width_units']
        
        # Set up the axes
        if num_rows == 1:
            # Single row layout
            ax.set_xlim(0, circuit_width)
            ax.set_ylim(-0.5, self.num_qudits - 0.5)
        else:
            # Multi-row layout
            ax.set_xlim(0, circuit_width)
            row_spacing = 1.5
            total_height = num_rows * self.num_qudits + (num_rows - 1) * row_spacing
            ax.set_ylim(-0.5, total_height - 0.5)
        
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw circuit rows
        gate_idx = 0
        for row_idx in range(num_rows):
            num_gates_in_row = gates_per_row[row_idx]
            
            # Calculate vertical offset for this row
            if num_rows == 1:
                y_offset = 0
            else:
                row_spacing = 1.5
                y_offset = row_idx * (self.num_qudits + row_spacing)
            
            # Draw qudit lines for this row
            for qudit_idx in range(self.num_qudits):
                y = self.num_qudits - 1 - qudit_idx + y_offset  # Invert to show qudit 0 at top
                ax.plot([0.5, circuit_width - 0.5], [y, y], 
                       color=self.wire_color, linewidth=1.5, zorder=1)
                
                # Label qudit lines (only for first row or at start of each row)
                ax.text(0.2, y, f'q{qudit_idx}', ha='right', va='center', 
                       fontsize=12, fontweight='bold')
                if row_idx == 0:
                    ax.text(0.3, y, '|0⟩', ha='left', va='center', fontsize=10)
            
            # Add row indicator for multi-row layouts
            if num_rows > 1:
                row_label = f"Gates {gate_idx+1}-{gate_idx+num_gates_in_row}"
                ax.text(circuit_width / 2, self.num_qudits + y_offset + 0.5, 
                       row_label, ha='center', va='bottom',
                       fontsize=10, style='italic', color='gray')
            
            # Draw gates in this row
            gate_x_position = 1.0
            for _ in range(num_gates_in_row):
                if gate_idx >= len(self.instructions):
                    break
                    
                gate = self.instructions[gate_idx]
                gate_type_enum = gate.gate_type
                targets = gate.target_qudits
                
                # Ensure targets is a list
                if isinstance(targets, int):
                    targets = [targets]
                
                # Convert target indices to y-coordinates (inverted) with offset
                target_ys = [self.num_qudits - 1 - t + y_offset for t in targets]
                
                if gate_type_enum == GateTypes.SINGLE:
                    # Single qudit gate
                    self._draw_single_gate(ax, gate, gate_x_position, target_ys[0])
                elif gate_type_enum == GateTypes.TWO:
                    # Two qudit gate
                    self._draw_two_gate(ax, gate, gate_x_position, target_ys)
                elif gate_type_enum == GateTypes.MULTI:
                    # Multi qudit gate
                    self._draw_multi_gate(ax, gate, gate_x_position, target_ys)
                
                gate_x_position += self.gate_spacing
                gate_idx += 1
        
        # Add title
        if num_rows == 1:
            title_y = self.num_qudits + 0.3
        else:
            row_spacing = 1.5
            total_height = num_rows * self.num_qudits + (num_rows - 1) * row_spacing
            title_y = total_height + 0.3
        
        ax.text(circuit_width / 2, title_y, title,
               ha='center', va='bottom', fontsize=16, fontweight='bold')
        
        # Add gate count information
        gate_counts = self._count_gates()
        info_text = f"Total gates: {len(self.instructions)}"
        if num_rows > 1:
            info_text += f" (displayed in {num_rows} rows)"
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
    
    Automatically adjusts layout for large circuits after decomposition,
    splitting into multiple rows if needed for readability.
    
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
    
    # Get layout info for the decomposed circuit
    layout_after = viz_after._calculate_layout()
    
    # Print layout information
    print("\nVisualization Layout:")
    print(f"  Figure size: {layout_after['figsize'][0]:.1f} x {layout_after['figsize'][1]:.1f} inches")
    print(f"  Number of rows: {layout_after['num_rows']}")
    if layout_after['num_rows'] > 1:
        print(f"  Gates per row: {layout_after['gates_per_row']}")
        print("  (Circuit wrapped into multiple rows for better readability)")
    
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
    
    Automatically adjusts layout based on circuit size, splitting into
    multiple rows if needed for large circuits.
    
    Args:
        circuit: The quantum circuit to visualize
        title: Title for the circuit diagram
        save_path: Optional path to save the figure
        
    Returns:
        Tuple of (fig, ax) for the visualization
    """
    viz = CircuitVisualizer(circuit)
    viz.print_circuit_summary()
    
    # Get layout info
    layout = viz._calculate_layout()
    
    # Print layout information
    print("\nVisualization Layout:")
    print(f"  Figure size: {layout['figsize'][0]:.1f} x {layout['figsize'][1]:.1f} inches")
    print(f"  Number of rows: {layout['num_rows']}")
    if layout['num_rows'] > 1:
        print(f"  Gates per row: {layout['gates_per_row']}")
        print("  (Circuit wrapped into multiple rows for better readability)")
    
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
