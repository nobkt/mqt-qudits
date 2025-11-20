# Noise Model Integration Guide

This document describes how to integrate noise model simulations into the
`quantum_dynamics_complete_comparison.ipynb` notebook.

## Overview

We add noise simulations for all three methods (Classical, Qubit, Qudit) to
compare noisy vs noiseless quantum dynamics.

## Integration Points

### 1. After Cell 3 (Physical Parameters): Add Noise Parameters

**New Markdown Cell:**
```markdown
## 2.2 ノイズモデルパラメータの設定

ノイズ付きシミュレーションのために、以下の物理的なノイズパラメータを定義します：

### デコヒーレンスパラメータ

- **T1 (エネルギー緩和時間)**: 励起状態から基底状態への減衰  
- **T2 (位相緩和時間)**: 位相コヒーレンスの損失
- **減極性エラー確率**: 量子ゲート適用時のランダムエラー

これらのパラメータは実験的に観測される値に基づいています。
```

**New Code Cell:**
```python
# Import noise simulation modules
import sys
import site
user_site = site.getusersitepackages()
if user_site not in sys.path:
    sys.path.insert(0, user_site)

from noise_simulation_implementations import NoiseParameters
from extended_noise_simulators import ClassicalNoisySimulator

# Initialize noise parameters
noise_params = NoiseParameters()
noise_params.print_info()
```

### 2. After Cell 6 (Classical Visualization): Add Classical Noisy Simulation

**New Markdown Cell:**
```markdown
### 3.2 ノイズ付き古典的シミュレーション

現象論的なデコヒーレンス（T1, T2過程）を追加した古典的シミュレーションを実行します。
```

**New Code Cell:**
```python
# Run noisy classical simulation
classical_noisy_sim = ClassicalNoisySimulator(classical_sim, noise_params)
classical_noisy_results = classical_noisy_sim.simulate_with_noise(
    T_total=params.T_total,
    N_steps=params.N_steps,
    initial_state_type=params.initial_state_type
)

# Plot comparison: noiseless vs noisy
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 7))

# Noiseless
times_clean = classical_results['times']
pops_clean = classical_results['populations']
N_S0_clean = [p['N_S0'] for p in pops_clean]
N_T1_clean = [p['N_T1'] for p in pops_clean]
N_S1_clean = [p['N_S1'] for p in pops_clean]

ax1.plot(times_clean, N_S0_clean, 'b-', linewidth=2.5, label=r'$N_{S_0}$', marker='o', markersize=5, alpha=0.8)
ax1.plot(times_clean, N_T1_clean, 'r-', linewidth=2.5, label=r'$N_{T_1}$', marker='s', markersize=5, alpha=0.8)
ax1.plot(times_clean, N_S1_clean, 'g-', linewidth=2.5, label=r'$N_{S_1}$', marker='^', markersize=5, alpha=0.8)
ax1.set_xlabel('Time (fs)', fontsize=13, fontweight='bold')
ax1.set_ylabel('Population', fontsize=13, fontweight='bold')
ax1.set_title('Noiseless Classical Simulation', fontsize=15, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_xlim(0, max(times_clean))
ax1.set_ylim(0, 4.5)

# Noisy
times_noisy = classical_noisy_results['times']
pops_noisy = classical_noisy_results['populations']
N_S0_noisy = [p['N_S0'] for p in pops_noisy]
N_T1_noisy = [p['N_T1'] for p in pops_noisy]
N_S1_noisy = [p['N_S1'] for p in pops_noisy]

ax2.plot(times_noisy, N_S0_noisy, 'b-', linewidth=2.5, label=r'$N_{S_0}$', marker='o', markersize=5, alpha=0.8)
ax2.plot(times_noisy, N_T1_noisy, 'r-', linewidth=2.5, label=r'$N_{T_1}$', marker='s', markersize=5, alpha=0.8)
ax2.plot(times_noisy, N_S1_noisy, 'g-', linewidth=2.5, label=r'$N_{S_1}$', marker='^', markersize=5, alpha=0.8)
ax2.set_xlabel('Time (fs)', fontsize=13, fontweight='bold')
ax2.set_ylabel('Population', fontsize=13, fontweight='bold')
ax2.set_title(f'Noisy Classical Simulation (T1={noise_params.T1}fs, T2={noise_params.T2}fs)', fontsize=15, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_xlim(0, max(times_noisy))
ax2.set_ylim(0, 4.5)

plt.tight_layout()
plt.show()

print("\\n" + "="*70)
print("Noise Impact on Classical Simulation")
print("="*70)
print(f"Final N_T1 difference: {N_T1_clean[-1] - N_T1_noisy[-1]:.4f}")
print(f"Final N_S1 difference: {N_S1_clean[-1] - N_S1_noisy[-1]:.4f}")
print(f"Final N_S0 difference: {N_S0_clean[-1] - N_S0_noisy[-1]:.4f}")
print("="*70)
```

### 3. Similar additions for Qubit and Qudit simulations

For Qubit:
- Use Qiskit Aer with the noise model created by `create_qiskit_noise_model(noise_params)`
- Run the qubit simulation with noise_model parameter

For Qudit:
- Use MQT-Qudits MISim backend with the noise model created by `create_qudit_noise_model(noise_params)`
- Run the qudit simulation with noise_model parameter

### 4. Update Comparison Section

Add a new comparison table that includes noise results:

```python
# Comparison: Noiseless vs Noisy for all methods
comparison_noise = {
    'Method': [
        'Classical (Noiseless)',
        'Classical (Noisy)',
        'Qubit (Noiseless)',
        'Qubit (Noisy)',
        'Qudit (Noiseless)',
        'Qudit (Noisy)'
    ],
    'Final N_T1': [
        f"{classical_results['populations'][-1]['N_T1']:.4f}",
        f"{classical_noisy_results['populations'][-1]['N_T1']:.4f}",
        # ... similar for qubit and qudit
    ],
    # ... more columns
}
```

## Implementation Notes

1. **NO HEURISTICS**: All noise implementations use standard quantum noise models
   - Classical: Lindblad master equation for T1/T2
   - Qubit: Qiskit Aer standard noise model (depolarizing + damping)
   - Qudit: MQT-Qudits SubspaceNoise (exact noise channels)

2. **NO FALLBACKS**: If noise simulation fails, we report the error, not fall back to noiseless

3. **Preserve Existing**: All existing cells remain unchanged

4. **Side-by-side comparison**: Show noiseless and noisy results together for direct comparison
