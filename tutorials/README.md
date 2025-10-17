# MQT Qudits Tutorials

This directory contains tutorial documentation for the MQT Qudits framework.

## Contents

### doc/

Theoretical documentation on quantum dynamics and Qudit-based quantum algorithms:

#### 1. **quantum_dynamics_molecular_triplet_states.md** (基礎理論)
Comprehensive theory of quantum dynamics for molecular systems with ground singlet state (S₀), excited triplet state (T₁), and excited singlet state (S₁). This foundational document covers:
- Triplet excitation energy transfer between adjacent molecules
- Triplet-triplet annihilation (TTA) processes  
- Fluorescence emission from excited singlet states
- Complete mathematical formulations with Hamiltonians, operators, and rate equations
- Master equations and rate equations for population dynamics

#### 2. **suzuki_trotter_decomposition_theory.md** (数値計算理論)
Detailed theory of Suzuki-Trotter decomposition for numerical simulation of quantum dynamics. This document provides:
- Baker-Campbell-Hausdorff formula and non-commutativity effects
- 1st, 2nd, and 4th order Suzuki-Trotter decomposition methods
- Error analysis and convergence properties
- Application to molecular triplet systems with explicit algorithms
- Implementation strategies including sparse matrices, symmetry utilization, and parallelization
- Complete pseudocode for programmatic implementation

#### 3. **qudit_quantum_algorithm_for_molecular_triplet_dynamics.md** (完全実装理論)
**Main comprehensive document** synthesizing the above theories into a complete, implementation-ready framework using MQT Qudits. This document provides:
- Qudit (Qutrit) representation of molecular 3-level systems
- Hamiltonian formulation in Qudit basis with explicit matrix representations
- Gate-level implementation mapping to MQT Qudits operations (VirtRz, R, cu_one, cu_two)
- Complete Suzuki-Trotter decomposition using Qudit gates
- Full mathematical derivations without omissions (eigenvalues, eigenvectors, time evolution operators)
- Multiple complete Python implementation examples using `mqt.qudits` API
- Observable calculations (populations, fluorescence intensity, correlation functions)
- Circuit optimization techniques (gate fusion, parallelization, adaptive time steps)
- Convergence analysis with benchmarks against exact diagonalization
- Complete simulation scripts ready to run (2213 lines, 67KB)
- Extensions for inhomogeneous systems and 2D lattices
- Comparison with qubit approaches demonstrating Qudit advantages
- Physical applications to organic photonics and TTA upconversion
- Future directions including VQA, tensor networks, and experimental implementation

#### 4. **mqt_qudits_gates_and_bases_reference.md** (ゲート完全リファレンス)
**Comprehensive reference for all quantum gates and bases** available in the MQT-Qudits framework. This document provides:
- Complete mathematical formulations with LaTeX for all 20 gate types
- Detailed explanations of computational basis and qudit representations
- Single-qudit gates: H, X, Z, S, R, Rz, Rh, VirtRz, Perm, NoiseX, NoiseY (11 gates)
- Two-qudit gates: CEx, CSum, LS, MS (4 gates)
- Multi-qudit gates: RandU (1 gate)
- Custom gates: CustomOne, CustomTwo, CustomMulti (3 gates)
- Auxiliary operations: Gell-Mann matrices (1 type)
- Matrix representations with concrete examples for various dimensions
- Properties, use cases, and Python implementation code
- Complete usage examples with `mqt.qudits` API (1125 lines, 28KB)

#### 5. **n_molecule_triplet_dynamics_basic_gates.md** (N-molecule system generalization)
**Complete generalization to N-molecule systems** using only basic quantum gates. This document provides:
- Generalization to arbitrary number of molecules N
- Implementation using only basic gates (VirtRz, R, RH, CEx, CSum) without CustomTwo
- Complete gate decomposition for all Hamiltonian terms
- 1D chain, 2D lattice, and arbitrary topology support
- Inhomogeneous systems with different molecular parameters
- Detailed gate sequences and implementation-ready code (2547 lines)

#### 6. **tutorial_completion_report.md** (実装完了報告)
**Comprehensive completion report** for the tutorial implementation. This document provides:
- Summary of all deliverables (notebook, implementation file, documentation)
- Verification of theoretical rigor (no heuristics, gates-only implementation)
- New features added (circuit visualization, Qudit count display)
- Integration of all reference documents
- Execution requirements and expected performance
- Future extension possibilities

#### 7. **exact_diagonalization_theory.md** (厳密対角化理論) ⭐ NEW
**Complete theory for analytical solution comparison** using exact diagonalization. This document provides:
- Mathematical formulation of exact diagonalization method
- Hamiltonian matrix construction (81×81 for 4-molecule system)
- Exact time evolution using eigenvalue decomposition
- Comparison methodology with Qudit quantum algorithm
- Fidelity calculation and error analysis
- Convergence theory and verification checklist
- Implementation algorithms with rigorous mathematical foundations
- No heuristic methods - only exact mathematical operations
- Complete Python implementation examples (350+ lines)

#### 8. **PR18_COMPLETION_VERIFICATION.md** (PR#18完了検証) ⭐ NEW
**Completion verification report for PR#18**. This document provides:
- Verification of notebook format fixes (proper Jupyter format with `\n` handling)
- Confirmation of tutorial completeness (all 12 cells verified)
- Verification that no heuristic methods are used
- Implementation file completeness check (843 lines, all required components)
- Documentation cross-reference verification
- Quality assurance and final approval for production use

**Reading order**:
- **For complete understanding**: Document 1 → Document 2 → Document 3
- **For implementation focus**: Start directly with Document 3
- **Gate reference**: Document 4 (comprehensive reference for all available gates)
- **N-molecule generalization**: Document 5
- **Implementation completion**: Document 6
- **Analytical solution verification**: Document 7
- **PR#18 completion verification**: Document 8

### Tutorial Notebook

#### **four_molecule_linear_chain_quantum_dynamics.ipynb** ✨ COMPLETE
**Complete tutorial implementation** of 4-molecule quantum dynamics using MQT-Qudits gates only (no heuristics). The notebook includes:

**Structure** (12 cells):
1. Title and Introduction with implementation policy
2. Theoretical Background (molecular states, Qudit representation, Hamiltonians)
3. Suzuki-Trotter Decomposition theory
4. Implementation Setup (imports and library loading)
5. Physical Parameters configuration
6. Quantum Circuit Construction with gate implementation
7. Circuit Visualization with Qudit count and gate statistics
8. Simulation Execution with progress tracking
9. Results Visualization and physical interpretation
10. **Exact Diagonalization for Analytical Solution** ⭐ NEW
11. **Comparison: Qudit Algorithm vs Analytical Solution** ⭐ NEW
12. Summary with references and future directions

**Key Features**:
- ✅ Uses only MQT-Qudits basic gates (VirtRz, CEx, R, Rh, Rz, X)
- ✅ CustomTwo gates automatically decomposed to basic gates via LogEntQRCEXPass
- ✅ No heuristic methods (no scipy.linalg.expm)
- ✅ Complete circuit visualization showing Qudit configuration and gate composition
- ✅ **Exact diagonalization for rigorous algorithm verification** ⭐ NEW
- ✅ **Fidelity calculation between Qudit and exact solutions** ⭐ NEW
- ✅ **Quantitative error analysis and convergence validation** ⭐ NEW
- ✅ Physical interpretation of quantum dynamics results
- ✅ Integration of all theoretical documentation

**Companion Implementation**: `mqt_qudits_four_molecule_implementation.py` (843 lines)

**Execution Time**: ~20-30 seconds for full notebook (N_steps=20)

**Status**: ✅ **COMPLETE WITH ANALYTICAL VERIFICATION**

### Supporting Files

- **NOTEBOOK_SUMMARY.md**: Overview of notebook structure and content
- **IMPLEMENTATION_VERIFICATION.md**: Detailed verification of implementation correctness
- **NOTEBOOK_MODIFICATION.md**: History of notebook modifications
- **REFACTORING_SUMMARY.md**: Summary of code refactoring

## Language

The documentation is written in Japanese (日本語) with complete LaTeX mathematical formulations and English-language Python code examples.
