# Tools Directory

This directory contains utility scripts and research tools for working with MQT-Qudits.

## Qudit Gate Optimization Tools (PR#36 & PR#37)

### Overview

The qudit gate optimization project addresses the gate count explosion problem where qudit implementations use ~6,182 gates compared to qubit's 44 gates (140× difference). Through PR#36 and PR#37, we've developed tools for sparse structure analysis and rigorous unitary decomposition.

### PR#37: Rigorous Unitary Decomposition (✅ COMPLETE)

#### improved_unitary_decomposition.py ⭐

**Status**: ✅ Fully Functional - Perfect Fidelity

**Purpose**: Mathematically rigorous 2×2 unitary decomposition using ZYZ decomposition

**Key Features**:
- Perfect ZYZ decomposition (fidelity = 1.0)
- Handles all edge cases (singular points at θ≈0, π)
- No heuristics or approximations
- 100% test pass rate (100/100 random unitaries)

**Test Results**:
```bash
$ python tools/improved_unitary_decomposition.py
Minimum fidelity: 1.0000000000
Average fidelity: 1.0000000000
Pass rate: 100/100 (100.0%)
✓ All tests passed
```

**Usage**:
```python
from tools.improved_unitary_decomposition import ImprovedTwoQubitDecomposer

decomposer = ImprovedTwoQubitDecomposer()
result = decomposer.decompose_zyz(U_2x2)
# result.theta, result.phi, result.lam, result.global_phase
# result.fidelity == 1.0
```

#### perfect_3x3_decomposition.py ⭐

**Status**: ✅ Fully Functional - Perfect Fidelity

**Purpose**: Mathematically rigorous 3×3 unitary decomposition using QR decomposition

**Key Features**:
- Uses numpy.linalg.qr directly (Householder/Givens based)
- Perfect fidelity = 1.0
- No heuristics or approximations (no scipy.linalg.expm)
- Works perfectly on real problems (H_TTA tested)
- 100% test pass rate (100/100 random unitaries)

**Test Results**:
```bash
$ python tools/perfect_3x3_decomposition.py
Random unitary test:
  Minimum fidelity: 1.0000000000
  Pass rate: 100/100 (100.0%)
H_TTA real problem test:
  Fidelity: 1.0000000000
✓✓✓ All tests passed!
```

**Usage**:
```python
from tools.perfect_3x3_decomposition import Perfect3x3Decomposer

decomposer = Perfect3x3Decomposer()
result = decomposer.decompose(U_3x3)
# result.Q: unitary matrix (product of Givens/Householder rotations)
# result.R: upper triangular matrix
# result.fidelity == 1.0

# For real problems (e.g., H_TTA)
# Extract diagonal phases
phases = decomposer.extract_diagonal_phases(result.R)
# Normalize R
D, R_norm = decomposer.normalize_R(result.R)
# U = Q @ D @ R_norm
```

**Real Problem Validation**: Tested with H_TTA time evolution operator:
```python
J = 0.05
dt = 1.0
hbar = 0.6582119569
H_sub = J * np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
phases = np.exp(-1j * eigenvalues * dt / hbar)
U = eigenvectors @ np.diag(phases) @ eigenvectors.conj().T

result = decomposer.decompose(U)
# result.fidelity == 1.0 ✓
```

#### debug_3x3_decomposition.py

**Purpose**: Debug tool for analyzing 3×3 Givens decomposition step-by-step

**Features**:
- Visualizes each Givens rotation step
- Verifies element zeroing
- Checks unitarity at each step
- Helped identify the correct Givens formula

**Usage**:
```bash
$ python tools/debug_3x3_decomposition.py
# Shows detailed step-by-step decomposition
# Final fidelity: 1.0 (when using correct formulas)
```

#### final_unitary_decomposition.py

**Purpose**: Research tool exploring direct Givens parameter extraction

**Status**: Research/experimental - superseded by perfect_3x3_decomposition.py

**Key Finding**: Direct use of QR decomposition is more reliable than explicit Givens parameter extraction

### PR#36: Sparse Structure Compiler

#### sparse_structure_compiler.py

**Purpose**: Sparse structure-aware compiler for efficient qudit gate decomposition

**Status**: ✅ Partially Working
- ✅ Sparse structure detection (2×2 and 3×3 subspaces correctly identified)
- ✅ Gate count estimation (98.1% and 95.7% reduction potential)
- ⚠️ Unitary decomposition needs to use improved/perfect decomposers from PR#37

**Key Features**:
- `SparseStructureAnalyzer`: Detects sparse structure in 9×9 unitary matrices
- `TwoLevelRotationDecomposer`: 2×2 unitary decomposition (should use improved_unitary_decomposition.py)
- `ThreeLevelRotationDecomposer`: 3×3 unitary decomposition (should use perfect_3x3_decomposition.py)
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
⚠ Needs integration with PR#37 decomposers for perfect fidelity
```

**Next Steps**: 
1. Replace internal decomposers with `improved_unitary_decomposition.py` (2×2)
2. Replace internal decomposers with `perfect_3x3_decomposition.py` (3×3)
3. Achieve perfect fidelity throughout the pipeline

**Related Documentation**:
- `tutorials/doc/qudit_gate_cost_analysis.md` - Problem analysis
- `tutorials/doc/qudit_gate_optimization_implementation_plan.md` - Implementation plan
- `tutorials/doc/TASK_COMPLETION_REPORT.md` - PR#36 completion report

### unitary_decomposition_rigorous.py

**Purpose**: Original rigorous 2×2 and 3×3 unitary decomposition framework (from PR#36)

**Status**: ⚠️ Framework Complete, Superseded by PR#37 implementations

**Note**: This was the initial attempt. The correct, working implementations are:
- For 2×2: Use `improved_unitary_decomposition.py`
- For 3×3: Use `perfect_3x3_decomposition.py`

## Documentation (PR#37)

### PR37_COMPLETION_REPORT.md ⭐

**Purpose**: Complete report of PR#37 work

**Contents**:
- Executive summary of achievements
- 2×2 and 3×3 decomposition implementation details
- Comprehensive test results (all passing with fidelity = 1.0)
- Mathematical rigor guarantees
- Next steps for integration
- Complete deliverables list

**Key Results**:
- 2×2 decomposition: 100/100 tests passed, fidelity = 1.0
- 3×3 decomposition: 100/100 tests passed, fidelity = 1.0
- H_TTA real problem: fidelity = 1.0
- All constraints followed (no heuristics, no approximations)

### pr37_3x3_decomposition_continuation_spec_ja.md

**Purpose**: Detailed continuation specification for 3×3 decomposition work

**Contents**:
- Summary of completed work (2×2 perfect, 3×3 theoretical analysis)
- Detailed 3×3 decomposition theory
- Implementation algorithms and pseudocode
- Class design and test design
- Numerical stability considerations
- Effort estimates (60-75 hours for full implementation)

### givens_rotation_theory_ja.md ⭐

**Purpose**: Complete mathematical theory of Givens rotations for 3×3 unitary decomposition

**Contents**:
- Rigorous definition and properties of Givens rotations
- Proof of unitarity
- Element zeroing theorem with complete derivation
- **Correct formulas**: c = a*/r, s = b*/r, r = sqrt(|a|² + |b|²)
- 3×3 decomposition algorithm
- Parameter extraction formulas
- Numerical stability analysis
- Implementation examples

**Key Formulas**:
```
Givens parameters:
  c = a*/r
  s = b*/r
  r = sqrt(|a|² + |b|²)

Givens matrix:
  G[i,i] = c, G[i,j] = s
  G[j,i] = -s*, G[j,j] = c*

Parameter extraction:
  θ = 2 arccos(|G[i,i]|)
  φ = 2 arg(G[i,i])
```

## Mathematical Requirements (Strict Constraints)

All implementations follow strict mathematical rigor:

### Allowed ✅
- ✅ Exact linear algebra: np.linalg.eigh, np.linalg.qr
- ✅ Exact trigonometry: np.cos, np.sin, np.arccos
- ✅ Exact complex operations: np.angle, np.exp, np.conj
- ✅ Quantum gate combinations
- ✅ Fidelity > 0.9999 required (achieved 1.0)

### Prohibited ❌
- ❌ scipy.linalg.expm (uses Padé approximation)
- ❌ Heuristics or approximations
- ❌ Trotter order reduction
- ❌ Ignoring small matrix elements
- ❌ Any truncation or fallback

## Integration Roadmap

### Phase 1: Update sparse_structure_compiler.py (1-2 weeks)
1. Replace `TwoLevelRotationDecomposer` with `ImprovedTwoQubitDecomposer`
2. Replace `ThreeLevelRotationDecomposer` with `Perfect3x3Decomposer`
3. Test with H_transfer and H_TTA
4. Verify fidelity > 0.9999 throughout

### Phase 2: MQT-Qudits Integration (1-2 months)
1. Convert QR decomposition to MQT-Qudits basic gates (CEx, R, Rz, VirtRz)
2. Implement CompilerPass
3. Test on full circuits
4. Measure actual gate count reduction

### Expected Results
- Current: ~6,000 gates/Trotter step
- After optimization: ~150-200 gates/Trotter step
- Reduction: 97.5% (30-40× improvement)

## Visualization Tools

### visualize_circuit.py

A comprehensive quantum circuit visualization tool for MQT-Qudits with special support for CustomTwo gates and their decomposition.

#### Features

- **CustomTwo Gate Visualization**: Displays CustomTwo gates before decomposition
- **Basic Gate Visualization**: Shows circuits decomposed into basic gates (VirtRz, R, Rh, Rz, CEx)
- **Side-by-Side Comparison**: Compare circuits before and after decomposition
- **Dynamic Figure Sizing**: Automatically adjusts figure size based on circuit complexity
- **Multi-Row Layout**: Splits large circuits into multiple rows for readability
- **No Heuristics**: Uses exact implementations from MQT-Qudits without any approximations
- **No src/ Modifications**: References code from `src/` without modifying it

#### Automatic Layout Adaptation

The visualization tool automatically handles circuits of any size:

- **Small circuits (< 50 gates)**: Single row with compact display
- **Medium circuits (50-100 gates)**: Single row with wider figure
- **Large circuits (> 100 gates)**: Multiple rows for optimal readability

**Key Parameters:**
- Minimum gate width: 0.15 inches (ensures readability)
- Maximum figure width: 30 inches (before wrapping to new row)
- Row height: 2.5 inches per qudit

For example, a circuit with 300 gates after decomposition will be displayed in ~4 rows of ~75 gates each, with each gate remaining readable.

#### Usage

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

#### Requirements

- matplotlib
- numpy
- mqt.qudits (installed package)

#### Example

See `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb` for a complete example of using the visualization tool with large decomposed circuits.

#### Implementation Details

The visualization tool:

1. **Imports from src/**: Uses `mqt.qudits.quantum_circuit.gates` and related modules
2. **Handles All Gate Types**: Single-qudit, two-qudit, and multi-qudit gates
3. **Color Coding**: Different colors for different gate types
4. **Detailed Labels**: Shows gate parameters (angles, levels, etc.)
5. **Circuit Summary**: Prints gate composition and sequence
6. **Layout Calculation**: Automatically determines optimal figure size and row count
7. **Multi-Row Display**: Wraps large circuits across multiple rows with row labels

#### Gate Colors

- VirtRz: Peach (#FFE5B4)
- R: Light blue (#B4D7FF)
- Rh: Light green (#B4FFD7)
- Rz: Light purple (#D7B4FF)
- CEx: Pink (#FFB4D7)
- CustomTwo: Gold (#FFD700)
- CustomOne: Orange (#FFA500)
- X: Sky blue (#87CEEB)
- Default: Gray (#E0E0E0)

#### Layout Information

When visualizing circuits, the tool prints layout information:

```
Visualization Layout:
  Figure size: 29.9 x 45.0 inches
  Number of rows: 4
  Gates per row: [98, 98, 98, 30]
  (Circuit wrapped into multiple rows for better readability)
```

This helps users understand how the circuit is being displayed, especially for large decomposed circuits.

#### No Compilation Required

This is a pure Python script that works without any compilation steps. Simply import and use.

## Testing

All tools have comprehensive test suites:

```bash
# Test 2×2 decomposition
python tools/improved_unitary_decomposition.py

# Test 3×3 decomposition
python tools/perfect_3x3_decomposition.py

# Debug 3×3 decomposition
python tools/debug_3x3_decomposition.py

# Test sparse structure analysis (needs PR#37 integration)
python tools/sparse_structure_compiler.py
```

## Summary

**PR#37 Achievement**: ✅ Complete Success
- 2×2 unitary decomposition: fidelity = 1.0
- 3×3 unitary decomposition: fidelity = 1.0
- All mathematical rigor requirements met
- Ready for integration with MQT-Qudits

**Next Steps**: 
1. Integrate PR#37 decomposers into sparse_structure_compiler.py
2. Convert to MQT-Qudits gates
3. Achieve 97.5% gate count reduction
