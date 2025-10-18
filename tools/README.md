# Tools Directory

This directory contains utility scripts for working with MQT-Qudits.

## visualize_circuit.py

A comprehensive quantum circuit visualization tool for MQT-Qudits with special support for CustomTwo gates and their decomposition.

### Features

- **CustomTwo Gate Visualization**: Displays CustomTwo gates before decomposition
- **Basic Gate Visualization**: Shows circuits decomposed into basic gates (VirtRz, R, Rh, Rz, CEx)
- **Side-by-Side Comparison**: Compare circuits before and after decomposition
- **No Heuristics**: Uses exact implementations from MQT-Qudits without any approximations
- **No src/ Modifications**: References code from `src/` without modifying it

### Usage

```python
import sys
sys.path.append('.')
from tools.visualize_circuit import visualize_circuit_with_decomposition

# Visualize circuit before and after decomposition
fig1, ax1, fig2, ax2 = visualize_circuit_with_decomposition(
    circuit_before,
    circuit_after,
    title_before="Circuit with CustomTwo Gates",
    title_after="Circuit with Basic Gates Only"
)
```

For a single circuit:

```python
from tools.visualize_circuit import visualize_circuit

fig, ax = visualize_circuit(circuit, title="My Quantum Circuit")
```

### Requirements

- matplotlib
- numpy
- mqt.qudits (installed package)

### Example

See `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb` for a complete example of using the visualization tool.

### Implementation Details

The visualization tool:

1. **Imports from src/**: Uses `mqt.qudits.quantum_circuit.gates` and related modules
2. **Handles All Gate Types**: Single-qudit, two-qudit, and multi-qudit gates
3. **Color Coding**: Different colors for different gate types
4. **Detailed Labels**: Shows gate parameters (angles, levels, etc.)
5. **Circuit Summary**: Prints gate composition and sequence

### Gate Colors

- VirtRz: Peach (#FFE5B4)
- R: Light blue (#B4D7FF)
- Rh: Light green (#B4FFD7)
- Rz: Light purple (#D7B4FF)
- CEx: Pink (#FFB4D7)
- CustomTwo: Gold (#FFD700)
- CustomOne: Orange (#FFA500)
- X: Sky blue (#87CEEB)
- Default: Gray (#E0E0E0)

### No Compilation Required

This is a pure Python script that works without any compilation steps. Simply import and use.
