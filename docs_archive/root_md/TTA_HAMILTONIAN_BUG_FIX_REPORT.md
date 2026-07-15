# TTA Hamiltonian Bug Fix Report

## Executive Summary

A **critical bug** was discovered and fixed in the TTA (Triplet-Triplet Annihilation) Hamiltonian implementation across multiple files in the quantum dynamics tutorial. The bug caused incorrect simulation results by missing half of the required coupling terms.

## Problem Description

### Physical Process

The TTA process describes the annihilation of two triplet excitations to produce one singlet excitation:

```
T₁ + T₁ → S₁ + S₀
```

This process has **two possible product states**:
1. `|T₁T₁⟩ → |S₀S₁⟩` (first molecule to ground, second to singlet)
2. `|T₁T₁⟩ → |S₁S₀⟩` (first molecule to singlet, second to ground)

### The Bug

The implementation was **missing the second coupling term**, only implementing:
- `|T₁T₁⟩ ↔ |S₀S₁⟩` ✓

But missing:
- `|T₁T₁⟩ ↔ |S₁S₀⟩` ✗

This caused the Hamiltonian to only represent half of the physical TTA process.

## Mathematical Details

### Correct TTA Hamiltonian

The complete TTA Hamiltonian for a pair of molecules (i, j) should be:

```
Ĥ_TTA = J [|S₀⟩ᵢ⟨T₁|ᵢ ⊗ |S₁⟩ⱼ⟨T₁|ⱼ 
         + |S₁⟩ᵢ⟨T₁|ᵢ ⊗ |S₀⟩ⱼ⟨T₁|ⱼ 
         + hermitian conjugates]
```

### Matrix Representation (9×9 for 2-qutrit system)

Basis ordering: `{|00⟩, |01⟩, |02⟩, |10⟩, |11⟩, |12⟩, |20⟩, |21⟩, |22⟩}`

Mapping: `|S₀⟩=|0⟩, |T₁⟩=|1⟩, |S₁⟩=|2⟩`

**CORRECT Matrix:**
```
H_TTA = J × ⎡0 0 0 0 0 0 0 0 0⎤
            ⎢0 0 0 0 0 0 0 0 0⎥
            ⎢0 0 0 0 1 0 0 0 0⎥  ← row 2: |S₀S₁⟩
            ⎢0 0 0 0 0 0 0 0 0⎥
            ⎢0 0 1 0 0 0 1 0 0⎥  ← row 4: |T₁T₁⟩
            ⎢0 0 0 0 0 0 0 0 0⎥
            ⎢0 0 0 0 1 0 0 0 0⎥  ← row 6: |S₁S₀⟩
            ⎢0 0 0 0 0 0 0 0 0⎥
            ⎣0 0 0 0 0 0 0 0 0⎦
                  ↑       ↑
              col 2   col 6
              |S₀S₁⟩  |S₁S₀⟩
```

Non-zero elements:
- `H[2,4] = J`: `|S₀S₁⟩ ↔ |T₁T₁⟩`
- `H[4,2] = J`: `|T₁T₁⟩ ↔ |S₀S₁⟩`
- `H[6,4] = J`: `|S₁S₀⟩ ↔ |T₁T₁⟩` ← **This was missing!**
- `H[4,6] = J`: `|T₁T₁⟩ ↔ |S₁S₀⟩` ← **This was missing!**

**INCORRECT Matrix (old implementation):**
```
H_TTA_old = J × ⎡0 0 0 0 0 0 0 0 0⎤
                ⎢0 0 0 0 0 0 0 0 0⎥
                ⎢0 0 0 0 1 0 0 0 0⎥
                ⎢0 0 0 0 0 0 0 0 0⎥
                ⎢0 0 1 0 0 0 0 0 0⎥  ← Missing elements at [4,6] and [6,4]
                ⎢0 0 0 0 0 0 0 0 0⎥
                ⎢0 0 0 0 0 0 0 0 0⎥  ← Row 6 was all zeros!
                ⎢0 0 0 0 0 0 0 0 0⎥
                ⎣0 0 0 0 0 0 0 0 0⎦
```

## Files Fixed

### 1. `tutorials/quantum_dynamics_complete_comparison.ipynb`

**Location:** Cell 5, `ClassicalSuzukiTrotterSimulator.build_H_TTA_pair()`

**Old Code:**
```python
def build_H_TTA_pair(self, mol_i: int, mol_j: int) -> np.ndarray:
    """1ペアのTTAハミルトニアンを構築（量子実装に合わせる）"""
    # J (|S0⟩_i|S1⟩_j⟨T1|_i⟨T1|_j + h.c.)
    op_i = self.S0_to_T1  # |S0⟩⟨T1|
    op_j = self.S1_to_T1  # |S1⟩⟨T1|
    term = self.build_two_site_operator(mol_i, mol_j, op_i, op_j)
    H_TTA = self.params.J * (term + term.conj().T)
    
    return H_TTA
```

**New Code:**
```python
def build_H_TTA_pair(self, mol_i: int, mol_j: int) -> np.ndarray:
    """1ペアのTTAハミルトニアンを構築（量子実装に合わせる）"""
    # FIXED: Complete TTA Hamiltonian with all 4 terms
    # Ĥ_TTA = J [|S0⟩_i⟨T1|_i ⊗ |S1⟩_j⟨T1|_j
    #          + |S1⟩_i⟨T1|_i ⊗ |S0⟩_j⟨T1|_j
    #          + hermitian conjugates]
    
    # Term 1: |S0⟩_i⟨T1|_i ⊗ |S1⟩_j⟨T1|_j
    term1 = self.build_two_site_operator(mol_i, mol_j, self.T1_to_S0, self.T1_to_S1)
    
    # Term 2: |S1⟩_i⟨T1|_i ⊗ |S0⟩_j⟨T1|_j
    term2 = self.build_two_site_operator(mol_i, mol_j, self.T1_to_S1, self.T1_to_S0)
    
    # Build Hamiltonian with all terms plus hermitian conjugates
    H_TTA = self.params.J * (term1 + term1.conj().T + term2 + term2.conj().T)
    
    return H_TTA
```

**Key Changes:**
- Added explicit `term2` for the second coupling
- Both terms now included with their hermitian conjugates
- Added detailed comments explaining the complete structure

### 2. `tutorials/exact_qubit_hamiltonians.py`

**Location:** `build_H_TTA_qubit_unitary()` function (lines 70-108)

**Old Code:**
```python
def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """..."""
    H = np.zeros((16, 16), dtype=complex)
    
    idx_01_01 = 0b0101  # = 5
    idx_00_10 = 0b0010  # = 2
    
    # Only one coupling term!
    H[idx_00_10, idx_01_01] = J
    H[idx_01_01, idx_00_10] = J
    
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    return U
```

**New Code:**
```python
def build_H_TTA_qubit_unitary(J: float, dt: float, hbar: float = 0.6582119569) -> np.ndarray:
    """..."""
    H = np.zeros((16, 16), dtype=complex)
    
    idx_T1_T1 = 0b0101  # = 5
    idx_S0_S1 = 0b0010  # = 2
    idx_S1_S0 = 0b1000  # = 8  ← NEW!
    
    # Term 1: |T1,T1⟩ ↔ |S0,S1⟩
    H[idx_S0_S1, idx_T1_T1] = J
    H[idx_T1_T1, idx_S0_S1] = J
    
    # Term 2: |T1,T1⟩ ↔ |S1,S0⟩ (THIS WAS MISSING!)
    H[idx_S1_S0, idx_T1_T1] = J
    H[idx_T1_T1, idx_S1_S0] = J
    
    U = scipy.linalg.expm(-1j * H * dt / hbar)
    return U
```

**Key Changes:**
- Added `idx_S1_S0` for the second product state
- Added the missing coupling terms at indices [5,8] and [8,5]
- Updated comments to clearly indicate both coupling terms

### 3. `tutorials/doc/theory_quantum_dynamics_complete_comparison.md`

**Fixed Sections:**

1. **TTA Matrix Representation (lines 750-805)**
   - Corrected the 9×9 matrix to show all 4 coupling terms
   - Fixed the non-zero element descriptions
   
2. **Eigenvalue Analysis (lines 775-805)**
   - Changed from "2 independent 2D subspaces" to "1 coherent 3D subspace"
   - Updated eigenvalues from `λ = ±J` to `λ = 0, ±√2 J`
   - Corrected eigenvectors to match the 3D subspace `{|02⟩, |11⟩, |20⟩}`

3. **Time Evolution Unitary (lines 807-836)**
   - Updated the 9×9 unitary matrix formula
   - Changed parameter from `φ = Jt/ℏ` to `θ = √2 Jt/ℏ`
   - Corrected matrix elements to match the eigenvalue structure

## Verification

### Test Script: `test_tta_fix.py`

Created comprehensive test script that verifies:
1. ✅ Correct 9×9 TTA Hamiltonian for 2-qutrit system
2. ✅ Correct 16×16 TTA Hamiltonian for 4-qubit system
3. ✅ Hermiticity of all Hamiltonians
4. ✅ Unitarity of time evolution operators
5. ✅ Physical coupling structure (|T₁T₁⟩ ↔ |S₀S₁⟩ + |S₁S₀⟩)

### Results

```
================================================================================
✅ ALL TESTS PASSED!
================================================================================

Summary:
- Qudit implementation: ✓ Correct (was already correct)
- Qubit implementation: ✓ Fixed
- Classical implementation: ✓ Fixed
- All three implementations are now mathematically consistent
- Physical TTA process |T₁T₁⟩ → |S₀S₁⟩ + |S₁S₀⟩ correctly modeled
```

## Impact

### Before Fix
- ❌ Only 50% of TTA process was modeled
- ❌ Simulation results were physically incorrect
- ❌ Qubit and classical implementations didn't match qudit implementation
- ❌ Theory document had wrong matrix representations

### After Fix
- ✅ Complete TTA process correctly modeled
- ✅ All three implementations (classical, qubit, qudit) are mathematically consistent
- ✅ Theory documentation is accurate
- ✅ Simulation results are physically correct

## Why This Bug Occurred

The bug likely originated from a common simplification error:

1. The physical process `T₁ + T₁ → S₁ + S₀` was correctly identified
2. One implementation (`|T₁T₁⟩ → |S₀S₁⟩`) was coded
3. The hermitian conjugate was added (`|S₀S₁⟩ → |T₁T₁⟩`)
4. **The second product state `|S₁S₀⟩` was forgotten**

This is because:
- `|S₀S₁⟩` and `|S₁S₀⟩` are **distinct states** (not symmetric)
- The full Hamiltonian needs **both** couplings to the same initial state `|T₁T₁⟩`
- Simply adding h.c. to term1 only gives 2 of the required 4 matrix elements

## Lessons Learned

1. **Always verify physical processes map to all required states**
   - Don't assume symmetry where there isn't any
   - Check all product states explicitly

2. **Theory documentation must match implementation**
   - Discrepancies between theory and code indicate bugs
   - Matrix representations should be explicitly verified

3. **Comprehensive testing is essential**
   - Test all matrix elements, not just overall behavior
   - Verify eigenvalue structure matches theory
   - Check consistency across different implementations

## References

- Theory document: `tutorials/doc/theory_quantum_dynamics_complete_comparison.md`
- Notebook: `tutorials/quantum_dynamics_complete_comparison.ipynb`
- Qubit implementation: `tutorials/exact_qubit_hamiltonians.py`
- Qudit implementation: `tutorials/exact_hamiltonian_builders.py` (was already correct)
- Test script: `test_tta_fix.py`

---

**Report Date:** 2025-11-13
**Fixed By:** GitHub Copilot Coding Agent
**Severity:** Critical - affects all quantum dynamics simulation results
**Status:** ✅ Fixed and Verified
