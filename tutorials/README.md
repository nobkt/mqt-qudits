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
- ✅ Complete circuit visualization with dynamic sizing and multi-row layout for large circuits
- ✅ Qudit configuration and gate composition statistics
- ✅ **Exact diagonalization for rigorous algorithm verification** ⭐ NEW
- ✅ **Fidelity calculation between Qudit and exact solutions** ⭐ NEW
- ✅ **Quantitative error analysis and convergence validation** ⭐ NEW
- ✅ Physical interpretation of quantum dynamics results
- ✅ Integration of all theoretical documentation

**Companion Implementation**: `mqt_qudits_four_molecule_implementation.py` (843 lines)

**Execution Time**: ~20-30 seconds for full notebook (N_steps=20)

**Status**: ✅ **COMPLETE WITH ANALYTICAL VERIFICATION**

### Circuit Visualization Tool

The tutorial uses an enhanced circuit visualization tool (`tools/visualize_circuit.py`) that automatically adapts to circuit complexity:

**Features**:
- **Dynamic Figure Sizing**: Automatically adjusts figure dimensions based on number of gates
- **Multi-Row Layout**: Large circuits (>100 gates) are wrapped across multiple rows for readability
- **Minimum Readability**: Maintains minimum gate width of 0.15 inches per gate
- **No Compression**: All gates remain visible and readable, even for circuits with hundreds of gates

**Example**: After decomposing CustomTwo gates, a circuit with 300+ basic gates is automatically displayed across multiple rows (~4 rows × 75 gates each), with each gate remaining clearly visible and labeled.

This ensures that the decomposed circuits with hundreds of basic gates are properly displayed, addressing the issue where the previous fixed-size visualization made large circuits unreadable.

### Qubit Version Documentation

#### **QUBIT_TUTORIAL_INDEX.md** 📘 NEW - Complete Guide
**Comprehensive index to all Qubit (2-level system) implementation resources**. This master document provides:
- Complete navigation guide for all Qubit-related documentation (4,500+ lines)
- Three distinct reading paths (beginner, implementer, decision-maker)
- Qubit vs Qudit comparison analysis
- Implementation status and continuation planning

**Location**: `tutorials/QUBIT_TUTORIAL_INDEX.md`

#### **doc/qubit/** - Qubit Implementation Documentation (4,500+ lines)
Complete theoretical and technical foundation for implementing the same molecular triplet dynamics using **Qubits (2-level systems) and Qiskit framework**. This comprehensive documentation includes:

1. **qubit_quantum_dynamics_molecular_triplet_states_theory.md** (1,079 lines, 32KB)
   - 2-qubit encoding of 3-level molecular systems
   - Pauli operator representation of Hamiltonians
   - Physical subspace preservation theory
   - Suzuki-Trotter decomposition for qubits

2. **qubit_implementation_specification.md** (1,409 lines, 34KB)
   - Qiskit gate catalog (20+ gates with matrix representations)
   - State encoding: |S0⟩→|00⟩, |T1⟩→|01⟩, |S1⟩→|10⟩, |11⟩=unused
   - Complete gate decomposition for H0, H_transfer, H_TTA
   - Performance specifications and benchmarks

3. **qubit_detailed_design.md** (1,447 lines, 41KB)
   - Complete system architecture
   - 6 main classes (PhysicalParameters, StateEncoder, HamiltonianGates, etc.)
   - Implementation-ready Python code (~500 lines)
   - Algorithm details and convergence analysis

4. **README.md**, **COMPLETION_SUMMARY.md**, **CONTINUATION_PLAN.md**
   - Project overview and navigation
   - Completion report and statistics
   - Strategic continuation planning with 3 scenarios

#### **qubit/** - Implementation Resources
**Current Status**: Documentation complete, implementation requires Qiskit dependency

- **README.md**: Overview and usage guide
- **IMPLEMENTATION_STATUS.md**: Detailed status report and requirements
- **IMPLEMENTATION_GUIDE.md**: Step-by-step implementation guide with code examples

**Key Comparison** (Qubit vs Qudit):
- Qubit: ~430 gates/step (8x more), requires 2 qubits/molecule, 68% unused states
- Qudit: ~55 gates/step, 1 qutrit/molecule, 0% unused states
- Advantage: Qubit version runs on widely available hardware (IBMQ, Rigetti, IonQ)

**Next Steps**: Implementation requires adding Qiskit to dependencies and 7-9 days of development work.

### Supporting Files

- **NOTEBOOK_SUMMARY.md**: Overview of notebook structure and content
- **IMPLEMENTATION_VERIFICATION.md**: Detailed verification of implementation correctness
- **NOTEBOOK_MODIFICATION.md**: History of notebook modifications
- **REFACTORING_SUMMARY.md**: Summary of code refactoring

## 疎構造認識コンパイラを使用したQuditシミュレーション

### 概要

4分子線形鎖の量子ダイナミクスシミュレーションにおいて、疎構造認識コンパイラを使用することで、従来のLogEntQRCEXPass方式と比較して**99.6%のゲート数削減**を実現しました。

### 実装比較

| 実装方式 | Qubit数/Qudit数 | ゲート数/ステップ | 20ステップ総数 | 備考 |
|---------|----------------|-----------------|---------------|------|
| Qubit | 8 qubits | 112 | 2,240 | 標準的な実装 |
| Qudit（従来） | 4 qutrits | 6,182 | 123,640 | LogEntQRCEXPass使用 |
| **Qudit（改良）** | 4 qutrits | **29** | **580** | **疎構造認識使用** |

### 主な改善点

1. **疎構造の自動検出**
   - H_transfer: 2×2部分空間を自動検出
   - H_TTA: 3×3部分空間を自動検出

2. **最適化された分解**
   - 2×2部分空間: ~1,000ゲート → ~1ゲート (99.9%削減)
   - 3×3部分空間: ~1,000ゲート → ~6ゲート (99.4%削減)

3. **数学的厳密性の保証**
   - 忠実度 = 1.0 を保証
   - ヒューリスティック・近似を一切使用しない
   - 厳密な線形代数のみ使用

### 使用方法

#### 基本的な使用

```python
from tutorials.mqt_qudits_four_molecule_sparse_implementation import (
    PhysicalParameters,
    SparseAwareMQTQuditTimeEvolution
)

# パラメータ初期化
params = PhysicalParameters()

# 疎構造認識版時間発展演算子
time_evol = SparseAwareMQTQuditTimeEvolution(params)

# 回路構築（従来と同じインターフェース）
# ... 回路初期化 ...

# ゲート追加
dt = 10.0  # fs
time_evol.add_H0_evolution_gates(circuit, dt/2)
time_evol.add_H_transfer_evolution_gates(circuit, dt/2)
time_evol.add_H_TTA_evolution_gates(circuit, dt/2)

# 統計レポート
print(time_evol.get_compilation_report())
```

#### ノートブックでの使用

`tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`を参照してください。

### 理論的基盤

詳細な理論的基盤は以下の文書を参照:
- `tutorials/doc/SPARSE_COMPILER_THEORETICAL_FOUNDATION_JA.md`
- `tutorials/doc/PR46_FRAMEWORK_INTEGRATION_SPECIFICATION_JA.md`
- `tutorials/doc/PR46_IMPLEMENTATION_DESIGN.md`

### テストとベンチマーク

#### テストの実行

```bash
# 疎構造認識実装のテスト
python3 test/python/tutorials/test_sparse_aware_implementation.py

# 出力例:
# ✓ 忠実度保存テスト合格
#   忠実度: 1.0000000000
# ✓ 疎構造検出テスト合格
#   2×2部分空間: 検出成功
#   3×3部分空間: 検出成功
# ✓ ゲート数削減テスト合格
#   検出された2×2部分空間: 3
#   検出された3×3部分空間: 3
# ✓ 統計レポート生成テスト合格
```

#### ベンチマークの実行

```bash
# 性能ベンチマーク
python3 tools/benchmark_sparse_compiler.py

# 出力例:
# 【ゲート数】
#   H0: 8 ゲート
#   H_transfer: 3 ゲート
#   H_TTA: 18 ゲート
#   合計（疎構造認識）: 29 ゲート
#   合計（従来方式推定）: 6,008 ゲート
#   削減率: 99.5%
#
# 【性能】
#   コンパイル時間: 17.75 ms
#   ピークメモリ使用量: 0.02 MB
```

### 参考文献

- PR#42-46: 疎構造認識コンパイラの開発
- PR#47: チュートリアルへの統合

## トラブルシューティング

### Q: ゲート数が期待より多い

A: 以下を確認してください:
1. `IntegratedSparseCompilerV2`が正しくインポートされているか
2. `optimize_gates=True`が設定されているか
3. 統計レポートで疎構造が正しく検出されているか

### Q: 忠実度が1.0でない

A: これは通常発生すべきではありません。以下を確認:
1. 数値許容誤差の設定（デフォルト: 1e-10）
2. ユニタリ行列の構築が正しいか
3. Issue報告をお願いします

### Q: ImportErrorが発生する

A: 以下を確認:
1. `mqt.qudits`がインストールされているか
2. `tools/`ディレクトリへのパスが正しく設定されているか
3. `numpy`, `scipy`がインストールされているか

## Language

The documentation is written in Japanese (日本語) with complete LaTeX mathematical formulations and English-language Python code examples.
