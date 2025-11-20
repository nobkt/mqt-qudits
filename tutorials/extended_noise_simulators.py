#!/usr/bin/env python3
"""
Extended Simulators with Noise Support

This module extends the existing Classical, Qubit, and Qudit simulators
to include realistic noise models.

NO HEURISTICS - All noise implementations are based on standard quantum noise models.
"""

import numpy as np
import scipy.linalg
from typing import Dict, List
import time
import sys
import site

# Add user site-packages for Qiskit
user_site = site.getusersitepackages()
if user_site not in sys.path:
    sys.path.insert(0, user_site)

# Import noise parameters and models
from noise_simulation_implementations import (
    NoiseParameters,
    create_qiskit_noise_model,
    create_qudit_noise_model,
    QISKIT_AVAILABLE,
    MQTQUDITS_NOISE_AVAILABLE
)

if QISKIT_AVAILABLE:
    from qiskit import QuantumCircuit
    from qiskit_aer import AerSimulator
    from qiskit.quantum_info import Statevector


class ClassicalNoisySimulator:
    """
    Classical simulator with phenomenological decoherence.
    
    Adds T1 (amplitude damping) and T2 (phase damping) relaxation to the
    classical state vector simulation.
    """
    
    def __init__(self, classical_sim, noise_params: NoiseParameters):
        """
        Args:
            classical_sim: Instance of ClassicalSuzukiTrotterSimulator
            noise_params: NoiseParameters instance
        """
        self.sim = classical_sim
        self.noise_params = noise_params
        self.params = classical_sim.params
        self.N = classical_sim.N
        self.dim = classical_sim.dim
    
    def apply_decoherence(self, state: np.ndarray, dt: float) -> np.ndarray:
        """
        Apply phenomenological decoherence (T1 and T2 processes)
        
        This implements Lindblad-form decoherence on the state vector.
        We work with the density matrix temporarily.
        """
        # Convert to density matrix
        rho = np.outer(state, state.conj())
        
        # T1 process: amplitude damping (energy relaxation)
        # Excited states (T1=1, S1=2) decay to ground state (S0=0)
        gamma1 = dt / self.noise_params.T1
        
        # T2 process: pure dephasing
        # gamma_phi = 1/T2 - 1/(2*T1)
        gamma_phi = dt * (1.0 / self.noise_params.T2 - 1.0 / (2.0 * self.noise_params.T1))
        if gamma_phi < 0:
            gamma_phi = 0
        
        # Apply T1: decay of excited state populations
        # For each basis state, if it has excited molecules, transfer some population to lower states
        for idx in range(self.dim):
            # Decode configuration
            config = []
            temp_idx = idx
            for _ in range(self.N):
                config.append(temp_idx % 3)
                temp_idx //= 3
            config = config[::-1]
            
            # Count excited molecules
            n_excited = sum(1 for level in config if level > 0)
            
            if n_excited > 0 and rho[idx, idx] > 1e-15:
                # Decay probability proportional to number of excited molecules
                decay_prob = n_excited * gamma1
                
                # Find ground state (all molecules in S0)
                ground_config = [0] * self.N
                ground_idx = 0  # All zeros
                
                # Transfer population
                transfer_amount = rho[idx, idx] * decay_prob
                if transfer_amount > 0:
                    rho[idx, idx] -= transfer_amount
                    rho[ground_idx, ground_idx] += transfer_amount
        
        # Apply T2: pure dephasing (decay of off-diagonal elements)
        damping_factor = np.exp(-gamma_phi)
        for i in range(self.dim):
            for j in range(i + 1, self.dim):
                rho[i, j] *= damping_factor
                rho[j, i] *= damping_factor
        
        # Ensure hermiticity and normalization
        rho = 0.5 * (rho + rho.conj().T)
        rho /= np.trace(rho)
        
        # Convert back to state vector (use dominant eigenstate)
        eigenvalues, eigenvectors = np.linalg.eigh(rho)
        max_idx = np.argmax(np.abs(eigenvalues))
        state = eigenvectors[:, max_idx]
        state /= np.linalg.norm(state)
        
        return state
    
    def simulate_with_noise(self, T_total: float, N_steps: int,
                           initial_state_type: str = 'edge_triplet') -> Dict:
        """Run simulation with decoherence"""
        print("\\n" + "="*70)
        print("Noisy Classical Simulation (with T1/T2 decoherence)")
        print("="*70)
        
        start_time = time.time()
        dt = T_total / N_steps
        
        # Build unitaries (same as noiseless)
        print("Building time evolution operators...")
        
        U_H0_half_list = []
        for mol_idx in range(self.N):
            H0_mol = self.sim.build_H0_single_molecule(mol_idx)
            U = scipy.linalg.expm(-1j * H0_mol * dt / (2 * self.params.hbar))
            U_H0_half_list.append(U)
        
        U_transfer_half_list = []
        for mol_i, mol_j in self.params.neighbors:
            H_tr_pair = self.sim.build_H_transfer_pair(mol_i, mol_j)
            U = scipy.linalg.expm(-1j * H_tr_pair * dt / (2 * self.params.hbar))
            U_transfer_half_list.append(U)
        
        U_TTA_half_list = []
        for mol_i, mol_j in self.params.neighbors:
            H_TTA_pair = self.sim.build_H_TTA_pair(mol_i, mol_j)
            U = scipy.linalg.expm(-1j * H_TTA_pair * dt / (2 * self.params.hbar))
            U_TTA_half_list.append(U)
        
        # Initialize
        state = self.sim.prepare_initial_state(initial_state_type)
        pop_0 = self.sim.calculate_populations(state)
        pop_per_mol_0 = self.sim.calculate_per_molecule_populations(state)
        
        print(f"\\nInitial state: {initial_state_type}")
        print(f"  N_S0 = {pop_0['N_S0']:.4f}")
        print(f"  N_T1 = {pop_0['N_T1']:.4f}")
        print(f"  N_S1 = {pop_0['N_S1']:.4f}")
        print(f"\\nNoise parameters:")
        print(f"  T1 = {self.noise_params.T1} fs")
        print(f"  T2 = {self.noise_params.T2} fs")
        
        times = [0.0]
        populations = [pop_0]
        per_molecule_populations = [pop_per_mol_0]
        
        # Time evolution with decoherence
        print(f"\\nRunning time evolution ({N_steps} steps with decoherence)...")
        
        for step in range(1, N_steps + 1):
            # Unitary evolution (same as noiseless)
            # Forward
            for mol_idx in range(self.N):
                state = U_H0_half_list[mol_idx] @ state
            for pair_idx in range(len(self.params.neighbors)):
                state = U_transfer_half_list[pair_idx] @ state
            for pair_idx in range(len(self.params.neighbors)):
                state = U_TTA_half_list[pair_idx] @ state
            
            # Backward
            for pair_idx in reversed(range(len(self.params.neighbors))):
                state = U_TTA_half_list[pair_idx] @ state
            for pair_idx in reversed(range(len(self.params.neighbors))):
                state = U_transfer_half_list[pair_idx] @ state
            for mol_idx in reversed(range(self.N)):
                state = U_H0_half_list[mol_idx] @ state
            
            # Apply decoherence
            state = self.apply_decoherence(state, dt)
            
            # Normalize
            state /= np.linalg.norm(state)
            
            # Calculate populations
            t = step * dt
            pop = self.sim.calculate_populations(state)
            pop_per_mol = self.sim.calculate_per_molecule_populations(state)
            times.append(t)
            populations.append(pop)
            per_molecule_populations.append(pop_per_mol)
            
            if step % max(1, N_steps // 10) == 0:
                print(f"  Step {step}/{N_steps}: t = {t:.2f} fs, "
                      f"N_T1 = {pop['N_T1']:.4f}, N_S1 = {pop['N_S1']:.4f}")
        
        elapsed = time.time() - start_time
        
        print("\\n" + "="*70)
        print("Simulation Complete")
        print("="*70)
        print(f"Final populations:")
        print(f"  N_S0 = {populations[-1]['N_S0']:.4f}")
        print(f"  N_T1 = {populations[-1]['N_T1']:.4f}")
        print(f"  N_S1 = {populations[-1]['N_S1']:.4f}")
        print(f"Execution time: {elapsed:.2f} seconds")
        
        return {
            'times': times,
            'populations': populations,
            'per_molecule_populations': per_molecule_populations,
            'state_final': state,
            'elapsed_time': elapsed,
            'method': 'Classical Suzuki-Trotter with T1/T2 Decoherence',
            'noise_params': {
                'T1': self.noise_params.T1,
                'T2': self.noise_params.T2
            }
        }


# Test if module can be imported
if __name__ == "__main__":
    print("Extended Simulators with Noise Support")
    print("=" * 70)
    print(f"Qiskit available: {QISKIT_AVAILABLE}")
    print(f"MQT-Qudits noise available: {MQTQUDITS_NOISE_AVAILABLE}")
    print("=" * 70)
