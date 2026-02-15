# Implementation Summary: Quantum Circuit Visualization

## Objective

Add quantum circuit visualization to `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb` with the following requirements:

1. Do NOT modify any code under `src/`
2. Create a Python script for circuit visualization under `./tools/`
3. The script should work without compilation
4. Visualize two versions:
   - Quantum circuit with CustomTwo gates
   - Quantum circuit after decomposing CustomTwo gates to basic gates
5. NO heuristics or fallback workarounds

## Implementation

### 1. Created `tools/visualize_circuit.py`

A comprehensive circuit visualization tool with the following features:

**Key Functions:**

- `visualize_circuit(circuit, title, save_path)`: Visualize a single circuit
- `visualize_circuit_with_decomposition(circuit_before, circuit_after, ...)`: Compare circuits before and after decomposition
- `CircuitVisualizer` class: Main visualization engine

**Features:**

- Handles all gate types (single-qudit, two-qudit, multi-qudit)
- Color-coded gates for easy identification
- Detailed gate labels showing parameters (angles, levels)
- Circuit summary with gate counts and sequence
- Side-by-side comparison for CustomTwo vs decomposed circuits

**Gate Support:**

- VirtRz (Virtual Z rotation)
- R (Rotation gate)
- Rh (Hermitian rotation)
- Rz (Z rotation)
- CEx (Controlled Exchange)
- CustomTwo (Custom two-qudit gate)
- CustomOne (Custom single-qudit gate)
- X (Pauli-X)

**Design Principles:**

- References code from `src/mqt/qudits` without modifying it
- No heuristics or approximations
- Pure Python implementation (no compilation required)
- Uses matplotlib for visualization

### 2. Updated `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`

Added a new cell (5.5) for circuit visualization:

- Imports the visualization tool from `tools/visualize_circuit`
- Visualizes `test_circuit` (with CustomTwo gates)
- Visualizes `decomposed_circuit` (with basic gates only)
- Provides informative output explaining the decomposition

### 3. Added Documentation

Created `tools/README.md` documenting:

- Tool features
- Usage examples
- Requirements
- Implementation details
- Gate color scheme

### 4. Updated `.gitignore`

Added `tools/__pycache__/` to prevent committing Python cache files

## Testing & Validation

### Test Results

```
Circuit before decomposition: 14 gates
  - VirtRz: 8 gates
  - CustomTwo: 6 gates

Circuit after decomposition: 6182 gates
  - VirtRz: 662 gates
  - R: 1560 gates
  - Rh: 1632 gates
  - Rz: 1266 gates
  - CEx: 1062 gates
```

### Requirements Validation

✅ **No src/ modifications**: Verified no files under `src/` were modified
✅ **Tools directory created**: `./tools/visualize_circuit.py` exists
✅ **Works without compilation**: Pure Python script
✅ **CustomTwo gates visualized**: Shows 6 CustomTwo gates before decomposition
✅ **Decomposed gates visualized**: Shows 6182 basic gates after decomposition
✅ **No heuristics**: All gates are exact implementations from MQT-Qudits

## Files Changed

1. `tools/visualize_circuit.py` (new) - 453 lines
2. `tools/README.md` (new) - 75 lines
3. `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb` (modified) - Added 1 cell
4. `.gitignore` (modified) - Added 1 line

## Usage Example

```python
import sys

sys.path.append(".")
from tools.visualize_circuit import visualize_circuit_with_decomposition

# Visualize circuit before and after decomposition
fig1, ax1, fig2, ax2 = visualize_circuit_with_decomposition(
    test_circuit,
    decomposed_circuit,
    title_before="分解前: CustomTwoゲートを含む量子回路",
    title_after="分解後: 基本ゲート（VirtRz, R, Rh, Rz, CEx）のみの量子回路",
)

plt.show()
```

## Technical Details

### Visualization Approach

1. Parse circuit instructions to extract gate types and targets
2. Create matplotlib figure with qudit lines
3. Draw gates as boxes with appropriate labels
4. Connect two-qudit gates with vertical lines
5. Add legend and statistics

### Gate Decomposition

- Uses `LogEntQRCEXPass` compiler from MQT-Qudits
- Decomposes CustomTwo gates into basic gates automatically
- No manual intervention or heuristics required
- Exact decomposition guaranteed by the compiler

### Color Scheme

- VirtRz: Peach (#FFE5B4)
- R: Light blue (#B4D7FF)
- Rh: Light green (#B4FFD7)
- Rz: Light purple (#D7B4FF)
- CEx: Pink (#FFB4D7)
- CustomTwo: Gold (#FFD700)
- CustomOne: Orange (#FFA500)
- X: Sky blue (#87CEEB)

## Compliance with Requirements

### Requirement 1: No src/ modifications ✅

Confirmed by: `git diff --name-only | grep "^src/" returns no results`

### Requirement 2: Tools directory ✅

Created `./tools/visualize_circuit.py` with full functionality

### Requirement 3: No compilation ✅

Pure Python script, works immediately after import

### Requirement 4: CustomTwo visualization ✅

Shows CustomTwo gates with gold color and dimension labels

### Requirement 5: Decomposed visualization ✅

Shows all basic gates (VirtRz, R, Rh, Rz, CEx) after decomposition

### Requirement 6: No heuristics ✅

All implementations use exact MQT-Qudits functions
No approximations or fallback mechanisms

## Conclusion

Successfully implemented quantum circuit visualization for the four-molecule linear chain quantum dynamics notebook. The implementation:

- Meets all specified requirements
- Provides comprehensive visualization capabilities
- Uses only exact implementations without heuristics
- Works without compilation
- Does not modify any src/ code
