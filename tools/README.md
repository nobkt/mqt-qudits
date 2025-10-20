# Tools Directory

This directory contains utility scripts and research tools for working with MQT-Qudits.

## Qudit Gate Optimization Tools (PR#36 & PR#37)

### sparse_structure_compiler.py (PR#36)

**Purpose**: Sparse structure-aware compiler for efficient qudit gate decomposition

**Status**: ✅ Partially Working
- ✅ Sparse structure detection (2×2 and 3×3 subspaces correctly identified)
- ✅ Gate count estimation (98.1% and 95.7% reduction potential)
- ⚠️ Unitary decomposition needs accuracy improvement (fidelity 0.24-0.63, requires > 0.9999)

**Key Features**:
- `SparseStructureAnalyzer`: Detects sparse structure in 9×9 unitary matrices
- `TwoLevelRotationDecomposer`: 2×2 unitary decomposition (needs improvement)
- `ThreeLevelRotationDecomposer`: 3×3 unitary decomposition (needs improvement)
- `SubspaceRotationOptimizer`: Estimates optimized gate counts

**Problem Context**:
- Qubit implementation: 44 gates/Trotter step
- Qudit implementation (current): 6,182 gates/Trotter step (140× worse!)
- Qudit implementation (theoretical with optimization): 158 gates/Trotter step
- Potential improvement: 97.4% reduction (38× better than current)

**Test Results**:
```bash
$ python sparse_structure_compiler.py
✓ H_transfer structure: 2×2 subspace detected, 810 → 15 gates (98.1% reduction)
✓ H_TTA structure: 3×3 subspace detected, 810 → 35 gates (95.7% reduction)
✗ 2×2 decomposition fidelity: 0.24 (requires > 0.9999)
✗ 3×3 decomposition fidelity: 0.63 (requires > 0.9999)
```

**Related Documentation**:
- `tutorials/doc/qudit_gate_cost_analysis.md` - Problem analysis
- `tutorials/doc/qudit_gate_optimization_implementation_plan.md` - Implementation plan
- `tutorials/doc/TASK_COMPLETION_REPORT.md` - PR#36 completion report

### unitary_decomposition_rigorous.py (PR#37)

**Purpose**: Mathematically rigorous 2×2 and 3×3 unitary decomposition framework

**Status**: ⚠️ Framework Complete, Needs Numerical Accuracy Improvements

**Key Features**:
- `RigorousTwoQubitDecomposer`: High-precision 2×2 unitary decomposition
  - ZYZ decomposition (needs algorithmic fixes)
  - ZXZ decomposition (needs algorithmic fixes)
  - Target: fidelity > 0.9999
- `RigorousThreeQuditDecomposer`: High-precision 3×3 unitary decomposition
  - Givens decomposition (needs algorithmic fixes)
  - Target: fidelity > 0.9999

**Mathematical Requirements** (Strict Constraints):
- ❌ NO heuristics or approximations
- ❌ NO scipy.linalg.expm (uses Padé approximation)
- ❌ NO Trotter order reduction
- ❌ NO ignoring small matrix elements
- ✅ ONLY exact linear algebra (np.linalg.eigh, qr, etc.)
- ✅ Fidelity > 0.9999 required

**Next Steps** (see documentation):
1. Study Qiskit's TwoQubitBasisDecomposer reference implementation
2. Fix ZYZ decomposition algorithm to achieve fidelity > 0.9999
3. Fix Givens decomposition algorithm to achieve fidelity > 0.9999
4. Test with H_transfer and H_TTA matrices

**Related Documentation**:
- `tutorials/doc/rigorous_unitary_decomposition_theory_ja.md` - Complete mathematical theory
- `tutorials/doc/qudit_optimization_continuation_specification_ja.md` - Full continuation spec (370-500 hours)
- `tutorials/doc/immediate_implementation_design_ja.md` - Immediate next steps (100-135 hours)
- `tutorials/doc/qudit_gate_optimization_continuation_summary.md` - English summary

**Remaining Work**:
- Phase 1 (Foundation): 100-135 hours - Fix decomposition accuracy
- Phase 2 (Integration): 90-130 hours - MQT-Qudits framework integration
- Phase 3 (Specialization): 140-180 hours - H_transfer/H_TTA optimized sequences
- Phase 4 (Testing): 70-100 hours - End-to-end validation

## Visualization Tools

## visualize_circuit.py

A comprehensive quantum circuit visualization tool for MQT-Qudits with special support for CustomTwo gates and their decomposition.

### Features

- **CustomTwo Gate Visualization**: Displays CustomTwo gates before decomposition
- **Basic Gate Visualization**: Shows circuits decomposed into basic gates (VirtRz, R, Rh, Rz, CEx)
- **Side-by-Side Comparison**: Compare circuits before and after decomposition
- **Dynamic Figure Sizing**: Automatically adjusts figure size based on circuit complexity
- **Multi-Row Layout**: Splits large circuits into multiple rows for readability
- **No Heuristics**: Uses exact implementations from MQT-Qudits without any approximations
- **No src/ Modifications**: References code from `src/` without modifying it

### Automatic Layout Adaptation

The visualization tool automatically handles circuits of any size:

- **Small circuits (< 50 gates)**: Single row with compact display
- **Medium circuits (50-100 gates)**: Single row with wider figure
- **Large circuits (> 100 gates)**: Multiple rows for optimal readability

**Key Parameters:**
- Minimum gate width: 0.15 inches (ensures readability)
- Maximum figure width: 30 inches (before wrapping to new row)
- Row height: 2.5 inches per qudit

For example, a circuit with 300 gates after decomposition will be displayed in ~4 rows of ~75 gates each, with each gate remaining readable.

### Usage

```python
import sys
sys.path.append('.')
from tools.visualize_circuit import visualize_circuit_with_decomposition

# Visualize circuit before and after decomposition
# The tool automatically adjusts layout for large decomposed circuits
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

See `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb` for a complete example of using the visualization tool with large decomposed circuits.

### Implementation Details

The visualization tool:

1. **Imports from src/**: Uses `mqt.qudits.quantum_circuit.gates` and related modules
2. **Handles All Gate Types**: Single-qudit, two-qudit, and multi-qudit gates
3. **Color Coding**: Different colors for different gate types
4. **Detailed Labels**: Shows gate parameters (angles, levels, etc.)
5. **Circuit Summary**: Prints gate composition and sequence
6. **Layout Calculation**: Automatically determines optimal figure size and row count
7. **Multi-Row Display**: Wraps large circuits across multiple rows with row labels

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

### Layout Information

When visualizing circuits, the tool prints layout information:

```
Visualization Layout:
  Figure size: 29.9 x 45.0 inches
  Number of rows: 4
  Gates per row: [98, 98, 98, 30]
  (Circuit wrapped into multiple rows for better readability)
```

This helps users understand how the circuit is being displayed, especially for large decomposed circuits.

### No Compilation Required

This is a pure Python script that works without any compilation steps. Simply import and use.
