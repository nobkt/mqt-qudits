# Fix for AssertionError in Four Molecule Tutorial

## Problem Statement

When executing `tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`, the following error occurred:

```python
AssertionError                            Traceback (most recent call last)
Cell In[3], line 28
     26 # H_transferのゲートを追加
     27 initial_count = len(test_circuit.instructions)
---> 28 time_evol.add_H_transfer_evolution_gates(test_circuit, dt_test)

File mqt_qudits_four_molecule_sparse_implementation.py:276
--> 276 self._add_gates_to_circuit(circuit, gate_info['gates'])

File mqt_qudits_four_molecule_sparse_implementation.py:359
--> 359 circuit.r(qudits[0], [params['level1'], params['level2'],
    360                      params['theta'], params['phi']])

File r.py:88, in R.validate_parameter
---> 88 assert parameter[1] < self.dimensions
AssertionError:
```

## Root Cause Analysis

### The Mathematical Problem

The sparse compiler was incorrectly treating **two-qudit operations** as single-qudit R gates:

1. **H_transfer Hamiltonian**: Creates a rotation between states |01⟩ and |10⟩

   - These states involve **both qudits**:
     - |01⟩ = qudit₀ in state |0⟩, qudit₁ in state |1⟩
     - |10⟩ = qudit₀ in state |1⟩, qudit₁ in state |0⟩
   - This is fundamentally a **two-qudit entangling operation**
   - Cannot be expressed as a single-qudit R gate

2. **H_TTA Hamiltonian**: Creates rotations in the |02⟩, |11⟩, |20⟩ subspace
   - Similarly involves **both qudits**
   - Cannot be decomposed into single-qudit gates alone

### The Technical Problem

The sparse compiler flow was:

```
1. Detect 2×2 or 3×3 active subspace (e.g., indices [1, 3] for |01⟩, |10⟩)
2. Generate Givens rotations in this subspace
3. Convert to MQT-Qudits gates
4. Output R gates with parameters: level1=1, level2=3
```

**The bug**: Indices 1 and 3 are **global indices** in the 9-dimensional composite Hilbert space (3⊗3), but R gates expect **local qudit indices** (0-2 for qutrits).

When `circuit.r(qudit, [level1=3, level2=...])` was called, the R gate validation checked:

```python
assert parameter[1] < self.dimensions  # Fails! 3 >= 3
```

## Solution

### Conceptual Fix

Added logic to detect whether the active subspace:

- **Spans a single qudit**: Use R/Rz/Rh gates with proper local index conversion
- **Spans multiple qudits**: Use CustomTwo gate (preserves exact unitary, will be decomposed)

### Implementation Details

#### 1. Added Helper Methods

**`_global_index_to_qudit_states(global_idx, dimensions)`**

- Converts global index in composite space to local qudit states
- Example: For 3⊗3 system, index 3 = |10⟩ = [1, 0]

**`_analyze_subspace(active_indices, dimensions)`**

- Detects which qudits are involved in the subspace
- Returns:
  - `type`: 'single_qudit' or 'multi_qudit'
  - `qudit_idx`: (for single-qudit) which qudit
  - `global_to_local`: (for single-qudit) index mapping
  - `involved_qudits`: (for multi-qudit) list of involved qudits

#### 2. Modified Gate Compilation

**`compile_unitary_to_gates(U, qudit_indices)`**

```python
subspace_type = self._analyze_subspace(active_subspace, dimensions=[3, 3])

if subspace_type['type'] == 'multi_qudit':
    # Use CustomTwo gate for operations spanning multiple qudits
    gates.append({
        'type': 'CustomTwo',
        'qudit_indices': qudit_indices,
        'params': {'unitary': U}
    })
else:
    # Use R/Rz/Rh gates with LOCAL indices for single-qudit operations
    for gate in result.gate_sequence.gates:
        # Convert global indices to local qudit indices
        local_level1 = subspace_type['global_to_local'][global_level1]
        local_level2 = subspace_type['global_to_local'][global_level2]
        ...
```

#### 3. Updated Circuit Building

**`_add_gates_to_circuit(circuit, gates)`**

- Added handler for CustomTwo gates: `circuit.cu_two(qudits, unitary)`
- Fixed gate method names: `cex` → `cx` (correct MQT-Qudits API)

## Mathematical Rigor

The fix maintains complete mathematical rigor:

✅ **No heuristics**: Subspace detection uses exact basis state analysis
✅ **No approximations**: CustomTwo preserves the full unitary matrix
✅ **No fallbacks**: Proper structural detection, not error handling
✅ **Preserves fidelity = 1.0**: CustomTwo will be exactly decomposed by LogEntQRCEXPass

This is the **correct** approach because:

1. Multi-qudit operations **must** use multi-qudit gates (CustomTwo)
2. Single-qudit operations can be optimized with R gates (using local indices)
3. The original design already used CustomTwo + decomposition

## Test Results

### Unit Tests (test_subspace_analysis.py)

All helper function tests pass:

- ✓ `test_global_index_to_qudit_states`: Index conversion works correctly
- ✓ `test_analyze_subspace_single_qudit`: Detects single-qudit subspaces
- ✓ `test_analyze_subspace_multi_qudit_h_transfer`: Detects H_transfer as multi-qudit
- ✓ `test_analyze_subspace_multi_qudit_h_tta`: Detects H_TTA as multi-qudit

### Expected Behavior

For the notebook `four_molecule_linear_chain_quantum_dynamics.ipynb`:

**Before fix:**

```
H0の時間発展: 8 個のVirtRzゲート
AssertionError (when adding H_transfer gates)
```

**After fix:**

```
H0の時間発展: 8 個のVirtRzゲート
H_transfer分解: sparse_2x2, 1ゲート (CustomTwo), 忠実度=1.0
H_TTA分解: sparse_3x3, 1ゲート (CustomTwo), 忠実度=1.0
```

The CustomTwo gates will then be decomposed by LogEntQRCEXPass into basic gates (CEx, R, Rz, etc.) during compilation, achieving the desired gate reduction while maintaining mathematical rigor.

## Files Changed

- `tutorials/mqt_qudits_four_molecule_sparse_implementation.py`
  - Added: `_global_index_to_qudit_states()` method
  - Added: `_analyze_subspace()` method
  - Modified: `compile_unitary_to_gates()` method
  - Modified: `_add_gates_to_circuit()` method

## Verification

To verify the fix works:

1. Install the package: `pip install -e .`
2. Run the notebook: `jupyter notebook tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`
3. Execute cells - should run without AssertionError
4. Verify gate counts and fidelity = 1.0

Or run the unit tests:

```bash
python3 test/python/tutorials/test_subspace_analysis.py
python3 test/python/tutorials/test_gate_generation_fix.py  # Requires NumPy
```

## Technical Notes

### Why CustomTwo for Multi-Qudit Operations?

A rotation between |01⟩ and |10⟩ mixes states where excitations transfer between qudits. This is analogous to a SWAP or entangling gate in qubits. It **cannot** be expressed as:

- Single-qudit R gate (only affects one qudit)
- Product of single-qudit gates (no entanglement)

It **must** use:

- Two-qudit gate (CEx, CustomTwo)
- Or sequence involving two-qudit gates

The sparse compiler detected the 2×2 structure correctly, but the conversion to gates needs to respect the multi-qudit nature.

### Performance Impact

Using CustomTwo for multi-qudit subspaces does not negate the sparse compiler benefits:

1. **Structure detection still works**: Identifies 2×2 and 3×3 subspaces (vs 9×9 dense)
2. **Decomposition is still optimized**: LogEntQRCEXPass can decompose a structured CustomTwo more efficiently than a random one
3. **Future optimization opportunity**: Could implement specialized two-qudit decomposers for these specific structures

The key insight: The sparse compiler's value is in **detecting** the structure. The fix ensures this detection correctly determines whether single-qudit or multi-qudit gates are needed.
