"""Deep analysis of energy transfer dynamics.

This script investigates why the energy transfer term doesn't cause
significant population changes in the edge_triplet initial state.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_math_utils import (
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_populations_from_density_matrix,
)
from gksl_physical_parameters import GKSLPhysicalParameters


def analyze_transfer_hamiltonian_coupling():
    """Analyze which states are coupled by H_transfer."""
    print("="*80)
    print("エネルギー移動ハミルトニアンの結合構造分析")
    print("="*80)

    params = GKSLPhysicalParameters(
        E_T=1.5,
        E_S=3.0,
        V=0.1,
        gamma_TTA=0,
        Gamma_fl=0,
        Gamma_ph=0,
        k_IC=0,
        k_ISC_ST=0,
        k_ISC_TS=0,
        N_molecules=4,
    )

    d = params.d
    N = params.N_molecules
    dim = d**N

    # Build Hamiltonians
    H_0 = build_onsite_hamiltonian(params)
    H_transfer = build_transfer_hamiltonian(params)

    # Prepare initial state: edge_triplet |1,0,0,1>
    psi0 = np.zeros(dim, dtype=np.complex128)
    index_initial = 1 * (d ** (N - 1)) + 1  # |1,0,0,1>
    psi0[index_initial] = 1.0

    print(f"\n初期状態: |1,0,0,1> (index={index_initial})")
    print(f"状態表現: 分子0=T1, 分子1=S0, 分子2=S0, 分子3=T1")

    # Find which states are directly coupled to the initial state
    print(f"\nH_transfer により初期状態から直接結合される状態:")
    coupled_states = []
    for idx in range(dim):
        coupling_strength = abs(H_transfer[index_initial, idx])
        if coupling_strength > 1e-10 and idx != index_initial:
            # Decode the state
            remainder = idx
            state_config = []
            for mol in range(N - 1, -1, -1):
                local_state = remainder % d
                remainder //= d
                state_config.insert(0, local_state)

            coupled_states.append({
                "index": idx,
                "config": state_config,
                "coupling": coupling_strength,
            })

            state_str = ",".join(str(s) for s in state_config)
            print(f"  |{state_str}> (index={idx}): 結合強度 = {coupling_strength:.6f}")

    # Compute energy differences
    print(f"\n初期状態のエネルギー固有値:")
    E_initial = H_0[index_initial, index_initial].real
    print(f"  E_initial = {E_initial:.6f} eV")

    print(f"\n結合状態のエネルギー固有値と差:")
    for state in coupled_states:
        idx = state["index"]
        E_coupled = H_0[idx, idx].real
        delta_E = E_coupled - E_initial
        state_str = ",".join(str(s) for s in state["config"])
        print(f"  |{state_str}>: E = {E_coupled:.6f} eV, ΔE = {delta_E:.6f} eV")

    # Key insight: Check if energy transfer conserves total energy
    print(f"\n重要な観察:")
    print(f"  H_transfer の形式: V * (|0>_i<1| ⊗ |1>_j<0| + h.c.)")
    print(f"  これは隣接分子間で T1 ↔ S0 の状態交換を引き起こす")
    print(f"  初期状態 |1,0,0,1> から:")
    for state in coupled_states:
        state_str = ",".join(str(s) for s in state["config"])
        print(f"    → |{state_str}>")

    return {
        "initial_index": int(index_initial),
        "initial_energy": float(E_initial),
        "coupled_states": [
            {
                "index": int(s["index"]),
                "config": s["config"],
                "coupling": float(s["coupling"]),
                "energy": float(H_0[s["index"], s["index"]].real),
                "energy_diff": float(H_0[s["index"], s["index"]].real - E_initial),
            }
            for s in coupled_states
        ],
    }


def analyze_population_dynamics():
    """Analyze population dynamics in more detail."""
    print("\n" + "="*80)
    print("詳細なポピュレーションダイナミクス分析")
    print("="*80)

    params = GKSLPhysicalParameters(
        E_T=1.5,
        E_S=3.0,
        V=0.1,
        gamma_TTA=0,
        Gamma_fl=0,
        Gamma_ph=0,
        k_IC=0,
        k_ISC_ST=0,
        k_ISC_TS=0,
        N_molecules=4,
    )

    d = params.d
    N = params.N_molecules
    dim = d**N

    H_0 = build_onsite_hamiltonian(params)
    H_transfer = build_transfer_hamiltonian(params)
    H_total = H_0 + H_transfer

    # Prepare initial state
    psi0 = np.zeros(dim, dtype=np.complex128)
    index_initial = 1 * (d ** (N - 1)) + 1
    psi0[index_initial] = 1.0
    rho0 = np.outer(psi0, psi0.conj())

    # Time evolution with finer time steps
    t_max = 100.0
    n_steps = 1000
    dt = t_max / n_steps

    print(f"\n時間発展シミュレーション (t_max={t_max}, n_steps={n_steps})")

    times = []
    populations_total = []
    populations_per_molecule = []
    state_probabilities = []  # Track specific state probabilities

    for k in range(n_steps + 1):
        t = k * dt
        U = scipy.linalg.expm(-1j * H_total * t)
        rho_t = U @ rho0 @ U.conj().T

        pops = compute_populations_from_density_matrix(rho_t, params)

        times.append(t)
        populations_total.append({
            "N_S0": pops["N_S0"],
            "N_T1": pops["N_T1"],
            "N_S1": pops["N_S1"],
        })
        populations_per_molecule.append(pops["per_molecule_populations"])

        # Track probability of initial state
        prob_initial = abs(rho_t[index_initial, index_initial])
        state_probabilities.append(prob_initial)

    # Calculate max changes
    N_S0_values = [p["N_S0"] for p in populations_total]
    N_T1_values = [p["N_T1"] for p in populations_total]
    N_S1_values = [p["N_S1"] for p in populations_total]

    max_delta_S0 = max(N_S0_values) - min(N_S0_values)
    max_delta_T1 = max(N_T1_values) - min(N_T1_values)
    max_delta_S1 = max(N_S1_values) - min(N_S1_values)

    print(f"\n総ポピュレーションの変動:")
    print(f"  N_S0: 初期={N_S0_values[0]:.6f}, 最終={N_S0_values[-1]:.6f}, 最大変動={max_delta_S0:.6e}")
    print(f"  N_T1: 初期={N_T1_values[0]:.6f}, 最終={N_T1_values[-1]:.6f}, 最大変動={max_delta_T1:.6e}")
    print(f"  N_S1: 初期={N_S1_values[0]:.6f}, 最終={N_S1_values[-1]:.6f}, 最大変動={max_delta_S1:.6e}")

    print(f"\n分子ごとのポピュレーション (初期 vs 最終):")
    for mol_idx in range(N):
        initial_pops = populations_per_molecule[0][mol_idx]
        final_pops = populations_per_molecule[-1][mol_idx]
        print(f"  分子{mol_idx}: S0 {initial_pops[0]:.6f}→{final_pops[0]:.6f}, "
              f"T1 {initial_pops[1]:.6f}→{final_pops[1]:.6f}, "
              f"S1 {initial_pops[2]:.6f}→{final_pops[2]:.6f}")

    print(f"\n初期状態 |1,0,0,1> の占有確率:")
    print(f"  初期: {state_probabilities[0]:.6f}")
    print(f"  最終: {state_probabilities[-1]:.6f}")
    print(f"  最小: {min(state_probabilities):.6f}")

    # Create plots
    os.makedirs("developing/figures", exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Total populations
    ax = axes[0, 0]
    ax.plot(times, N_S0_values, label="N_S0", linewidth=2)
    ax.plot(times, N_T1_values, label="N_T1", linewidth=2)
    ax.plot(times, N_S1_values, label="N_S1", linewidth=2)
    ax.set_xlabel("Time (eV⁻¹·ħ)", fontsize=12)
    ax.set_ylabel("Population", fontsize=12)
    ax.set_title("Total Population Dynamics", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Plot 2: Per-molecule populations for molecule 0
    ax = axes[0, 1]
    mol0_S0 = [p[0][0] for p in populations_per_molecule]
    mol0_T1 = [p[0][1] for p in populations_per_molecule]
    mol0_S1 = [p[0][2] for p in populations_per_molecule]
    ax.plot(times, mol0_S0, label="S0", linewidth=2)
    ax.plot(times, mol0_T1, label="T1", linewidth=2)
    ax.plot(times, mol0_S1, label="S1", linewidth=2)
    ax.set_xlabel("Time (eV⁻¹·ħ)", fontsize=12)
    ax.set_ylabel("Population", fontsize=12)
    ax.set_title("Molecule 0 Population", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Plot 3: Per-molecule populations for molecule 1
    ax = axes[1, 0]
    mol1_S0 = [p[1][0] for p in populations_per_molecule]
    mol1_T1 = [p[1][1] for p in populations_per_molecule]
    mol1_S1 = [p[1][2] for p in populations_per_molecule]
    ax.plot(times, mol1_S0, label="S0", linewidth=2)
    ax.plot(times, mol1_T1, label="T1", linewidth=2)
    ax.plot(times, mol1_S1, label="S1", linewidth=2)
    ax.set_xlabel("Time (eV⁻¹·ħ)", fontsize=12)
    ax.set_ylabel("Population", fontsize=12)
    ax.set_title("Molecule 1 Population", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Plot 4: Initial state probability
    ax = axes[1, 1]
    ax.plot(times, state_probabilities, linewidth=2, color="purple")
    ax.set_xlabel("Time (eV⁻¹·ħ)", fontsize=12)
    ax.set_ylabel("Probability", fontsize=12)
    ax.set_title("Initial State |1,0,0,1> Occupation", fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("developing/figures/energy_transfer_dynamics_iteration29.png", dpi=150)
    print(f"\nプロットを保存しました: developing/figures/energy_transfer_dynamics_iteration29.png")

    return {
        "max_delta_S0": float(max_delta_S0),
        "max_delta_T1": float(max_delta_T1),
        "max_delta_S1": float(max_delta_S1),
        "initial_state_prob_final": float(state_probabilities[-1]),
    }


def main():
    """Run deep analysis of energy transfer."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    print("Energy Transfer Deep Analysis - Iteration 29")
    print(f"Timestamp: {timestamp}\n")

    results = {}

    # Analysis 1: Transfer Hamiltonian coupling structure
    results["coupling_analysis"] = analyze_transfer_hamiltonian_coupling()

    # Analysis 2: Population dynamics
    results["dynamics_analysis"] = analyze_population_dynamics()

    # Save results
    os.makedirs("developing/verification_results", exist_ok=True)
    output_file = f"developing/verification_results/energy_transfer_analysis_iteration29_{timestamp}.json"

    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n{'='*80}")
    print(f"結果を保存しました: {output_file}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
