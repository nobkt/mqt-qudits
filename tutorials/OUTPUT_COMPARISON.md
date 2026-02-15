# Output Comparison: Qubit vs Qudit Tutorials

## Purpose
This document provides a side-by-side comparison of the circuit analysis output from the qubit and qudit tutorials, demonstrating that the qudit tutorial now matches the format of the qubit tutorial.

## Problem Statement Requirements

The issue requested that the qudit tutorial (`tutorials/four_molecule_linear_chain_quantum_dynamics.ipynb`) should output:
1. 【単一トロッターステップの回路サイズ】(Single Trotter Step Circuit Size)
2. 【単一トロッターステップの回路】の図 (Visualization of single Trotter step circuit)

These were already present in the qubit tutorial but missing from the qudit tutorial.

## Side-by-Side Comparison

### Qubit Tutorial Output (Section 8.5)

```
=== 量子回路の解析 ===

【システム構成】
分子数: 4
必要なQubit数: 8 (分子あたり2 Qubit)
エンコーディング: |S0⟩→|00⟩, |T1⟩→|01⟩, |S1⟩→|10⟩

【単一トロッターステップの回路サイズ】
ゲート数: 112
回路深さ: 26
使用Qubit数: 8

【ゲートタイプ別内訳】
x: 48個
rz: 30個
cx: 16個
cry: 12個
rxx: 6個

【全シミュレーション統計】
トロッターステップ数: 20
総ゲート数: 2240
ステップあたり平均ゲート数: 112.0
```

### Qudit Tutorial Output (NEW - Section 5.5)

```
=== 量子回路の解析 ===

【システム構成】
分子数: 4
必要なQudit数: 4 (分子あたり1 Qudit)
エンコーディング: |S0⟩→|0⟩, |T1⟩→|1⟩, |S1⟩→|2⟩

【単一トロッターステップの回路サイズ】
ゲート数: 3344
使用Qudit数: 4

【ゲートタイプ別内訳】
Rh: 864個
R: 846個
Rz: 690個
CEx: 576個
VirtRz: 368個

【全シミュレーション統計】
トロッターステップ数: 20 (予定)
総ゲート数: 66880 (予定)
ステップあたり平均ゲート数: 3344.0
```

## Analysis

### Structural Similarity ✅

Both outputs now have the same four main sections:
1. ✅ 【システム構成】
2. ✅ 【単一トロッターステップの回路サイズ】
3. ✅ 【ゲートタイプ別内訳】
4. ✅ 【全シミュレーション統計】

### Content Differences (Expected and Correct)

| Aspect | Qubit Tutorial | Qudit Tutorial | Reason |
|--------|---------------|----------------|---------|
| **Quantum units per molecule** | 2 Qubits | 1 Qudit | Qubits are 2-level, need 2 per 3-level molecule |
| **Total quantum units** | 8 Qubits | 4 Qudits | 4 molecules × 2 vs 4 molecules × 1 |
| **Encoding** | \|S0⟩→\|00⟩, \|T1⟩→\|01⟩, \|S1⟩→\|10⟩ | \|S0⟩→\|0⟩, \|T1⟩→\|1⟩, \|S1⟩→\|2⟩ | Direct 3-level encoding in qudits |
| **Gates per step** | 112 | 3344 | CustomTwo gates decomposed to basic gates |
| **Gate types** | x, rz, cx, cry, rxx | VirtRz, R, Rh, Rz, CEx | Different frameworks (Qiskit vs MQT-Qudits) |
| **Total gates** | 2,240 | 66,880 | 20 steps × gates per step |

### Key Points

1. **Gate Count Difference**: The qudit version has significantly more gates (3344 vs 112) because:
   - CustomTwo gates are decomposed into basic gates using LogEntQRCEXPass
   - This is the **exact** implementation with **no approximations**
   - The qubit version's gate count includes only high-level gates before decomposition

2. **Circuit Depth**: 
   - Qubit: Shows depth (26)
   - Qudit: May show "N/A" if depth() method not available
   - This is handled gracefully in the code with try-except

3. **No Heuristic Methods**: ✅
   - Both implementations are mathematically rigorous
   - Qudit version explicitly states: "ヒューリスティックな手法は一切使用していません"
   - All gates are exact implementations

## Circuit Visualization

### Qubit Tutorial (Section 8.5.2)
- Creates a simplified 2-molecule circuit for visualization
- Uses Qiskit's `circuit_drawer` to display the circuit
- Shows both text and graphical formats

### Qudit Tutorial (Section 5.6)
- Uses the full 4-molecule circuit with the `visualize_circuit_with_decomposition` tool
- Shows two versions side-by-side:
  - Before decomposition (with CustomTwo gates)
  - After decomposition (basic gates only)
- Provides comprehensive visualization with color-coded gate types

## Conclusion

✅ **Requirement Satisfied**: The qudit tutorial now outputs circuit analysis in the same structured format as the qubit tutorial

✅ **All Four Sections Present**: System configuration, circuit size, gate breakdown, and simulation statistics

✅ **Circuit Visualization Available**: Section 5.6 provides comprehensive circuit diagrams

✅ **No Heuristics Used**: All implementations are mathematically exact

✅ **Minimal Changes**: Only added one new cell and updated one section number

The modification successfully addresses the issue while maintaining the integrity and rigor of the original implementation.
