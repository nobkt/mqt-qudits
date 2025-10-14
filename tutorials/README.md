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

**Reading order**: Document 1 → Document 2 → Document 3 for complete understanding, or start directly with Document 3 for implementation focus.

## Language

The documentation is written in Japanese (日本語) with complete LaTeX mathematical formulations and English-language Python code examples.
