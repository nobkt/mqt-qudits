#!/usr/bin/env python3
"""
Qudit Noise Simulator

This module provides noise simulation support for Qudit-based quantum dynamics.
Uses MQT-qudits NoiseModel and NoisyCircuitFactory - no heuristics or fallbacks.
"""

import numpy as np
from typing import Dict
from mqt.qudits.simulation.noise_tools import NoiseModel, SubspaceNoise, Noise, NoisyCircuitFactory


def create_qudit_noise_model():
    """
    Create a realistic noise model for Qutrit (3-level) systems.
    
    Applies level-specific noise to qutrit transitions:
    - (0,1): S₀ ↔ T₁ transition
    - (1,2): T₁ ↔ S₁ transition  
    - (0,2): S₀ ↔ S₁ transition (direct, typically weaker)
    
    Returns:
        NoiseModel: MQT-qudits NoiseModel instance
        dict: Noise parameters used
    """
    noise_model = NoiseModel()
    
    # Noise parameters (comparable to qubit noise levels)
    depol_01 = 0.001   # 0.1% for S₀↔T₁ transition
    depol_12 = 0.001   # 0.1% for T₁↔S₁ transition
    depol_02 = 0.002   # 0.2% for S₀↔S₁ (direct transition, higher error)
    dephase = 0.002    # 0.2% dephasing for all transitions
    
    # Define subspace noise for all level transitions in one object
    # We need to create a single SubspaceNoise with all transitions
    # to avoid "same level defined multiple times" error
    
    # Create combined subspace noise for local gates
    local_subspace = SubspaceNoise(0.0, 0.0, (0, 1))  # Initialize with one transition
    # Manually add all transitions with different parameters
    local_subspace.subspace_w_probs[(0, 1)] = Noise(depol_01, dephase)
    local_subspace.subspace_w_probs[(1, 2)] = Noise(depol_12, dephase)
    local_subspace.subspace_w_probs[(0, 2)] = Noise(depol_02, dephase)
    
    # Apply to local (single-qudit) gates
    local_gates = ['x', 'z', 'h', 'rz', 'r', 'virtrz', 's', 'rh']
    noise_model.add_quantum_error_locally(local_subspace, local_gates)
    
    # Apply to nonlocal (two-qudit) gates
    # These gates typically have higher error rates in real hardware
    nonlocal_gates = ['cx', 'cex', 'csum', 'ls', 'ms']
    
    # Target qudit noise (higher error for 2-qudit gates)
    target_subspace = SubspaceNoise(0.0, 0.0, (0, 1))
    target_subspace.subspace_w_probs[(0, 1)] = Noise(depol_01 * 10, dephase * 2)
    target_subspace.subspace_w_probs[(1, 2)] = Noise(depol_12 * 10, dephase * 2)
    noise_model.add_nonlocal_quantum_error_on_target(target_subspace, nonlocal_gates)
    
    # Control qudit noise
    control_subspace = SubspaceNoise(0.0, 0.0, (0, 1))
    control_subspace.subspace_w_probs[(0, 1)] = Noise(depol_01 * 5, dephase)
    control_subspace.subspace_w_probs[(1, 2)] = Noise(depol_12 * 5, dephase)
    noise_model.add_nonlocal_quantum_error_on_control(control_subspace, nonlocal_gates)
    
    params = {
        'depol_01': depol_01,
        'depol_12': depol_12,
        'depol_02': depol_02,
        'dephase': dephase,
        'two_qudit_factor': 10  # Multiplier for 2-qudit gate errors
    }
    
    return noise_model, params


class QuditNoiseSimulator:
    """
    Wrapper for Qudit simulation with noise support.
    
    This class applies MQT-qudits noise models to quantum circuits
    and simulates noisy evolution of molecular quantum dynamics.
    """
    
    def __init__(self, noise_model, noise_params):
        """
        Initialize Qudit noise simulator.
        
        Args:
            noise_model: MQT-qudits NoiseModel instance
            noise_params: Dict with noise parameters
        """
        self.noise_model = noise_model
        self.noise_params = noise_params
        
        print(f"\nQudit Noise Simulator initialized")
        print(f"  Noise model: SubspaceNoise for qutrits")
        print(f"  Basis gates: {noise_model.basis_gates}")
        print(f"  Subspace (0,1) depol: {noise_params['depol_01']}")
        print(f"  Subspace (1,2) depol: {noise_params['depol_12']}")
        print(f"  Dephasing rate: {noise_params['dephase']}")
    
    def apply_noise_to_circuit(self, clean_circuit):
        """
        Apply noise model to a clean quantum circuit.
        
        Args:
            clean_circuit: MQT-qudits QuantumCircuit (noise-free)
        
        Returns:
            Noisy circuit with noise gates inserted
        """
        factory = NoisyCircuitFactory(self.noise_model, clean_circuit)
        noisy_circuit = factory.generate_circuit()
        
        print(f"\nNoise application:")
        print(f"  Clean circuit gates: {clean_circuit.number_gates}")
        print(f"  Noisy circuit gates: {noisy_circuit.number_gates}")
        print(f"  Noise gates added: {noisy_circuit.number_gates - clean_circuit.number_gates}")
        
        return noisy_circuit
    
    def estimate_noise_impact(self, clean_results, noisy_results):
        """
        Estimate the impact of noise on simulation results.
        
        Args:
            clean_results: Results from noise-free simulation
            noisy_results: Results from noisy simulation
        
        Returns:
            Dict with noise impact metrics
        """
        times = clean_results['times']
        pops_clean = clean_results['populations']
        pops_noisy = noisy_results['populations']
        
        # Calculate population differences
        errors = []
        for pc, pn in zip(pops_clean, pops_noisy):
            error = np.sqrt(
                (pc['N_S0'] - pn['N_S0'])**2 +
                (pc['N_T1'] - pn['N_T1'])**2 +
                (pc['N_S1'] - pn['N_S1'])**2
            )
            errors.append(error)
        
        # Analyze error growth
        early_error = np.mean(errors[:len(errors)//4])  # First quarter
        late_error = np.mean(errors[3*len(errors)//4:])  # Last quarter
        error_growth = late_error / early_error if early_error > 1e-10 else 0
        
        return {
            'mean_error': np.mean(errors),
            'max_error': np.max(errors),
            'final_error': errors[-1],
            'error_growth_factor': error_growth,
            'early_error': early_error,
            'late_error': late_error,
            'errors_over_time': errors,
            'times': times
        }


def demonstrate_qudit_noise_model():
    """
    Demonstrate Qudit noise model setup and capabilities.
    
    This function shows how to:
    1. Create a noise model for qutrits
    2. Apply noise to circuits  
    3. Simulate with noise (conceptual)
    
    Note: Full integration requires existing Qudit simulation code
    """
    print("\n" + "="*70)
    print("Qudit Noise Model Demonstration")
    print("="*70)
    
    # Create noise model
    noise_model, params = create_qudit_noise_model()
    
    print("\n1. Noise Model Created:")
    print(f"   - Subspace-specific noise for qutrits")
    print(f"   - (0,1) transition: {params['depol_01']*100:.2f}% depol, {params['dephase']*100:.2f}% dephase")
    print(f"   - (1,2) transition: {params['depol_12']*100:.2f}% depol, {params['dephase']*100:.2f}% dephase")
    print(f"   - (0,2) transition: {params['depol_02']*100:.2f}% depol, {params['dephase']*100:.2f}% dephase")
    
    print("\n2. Noise Application Method:")
    print("   ```python")
    print("   from mqt.qudits.simulation.noise_tools import NoisyCircuitFactory")
    print("   factory = NoisyCircuitFactory(noise_model, clean_circuit)")
    print("   noisy_circuit = factory.generate_circuit()")
    print("   ```")
    
    print("\n3. Integration with Existing Qudit Simulation:")
    print("   - Build quantum circuit for dynamics (as in current implementation)")
    print("   - Apply NoisyCircuitFactory to insert noise gates")
    print("   - Simulate noisy circuit with existing backend")
    print("   - Compare results with noise-free simulation")
    
    print("\n4. Noise Model Details:")
    print(noise_model)
    
    print("\n" + "="*70)
    print("Qudit Noise Model Ready for Integration")
    print("="*70)
    
    return noise_model, params


# Placeholder for full Qudit noisy simulation
# (Requires integration with existing Qudit implementation)
def run_qudit_noisy_simulation_placeholder(params, noise_model):
    """
    Placeholder for full Qudit noisy simulation.
    
    Complete implementation requires:
    1. Access to Qudit circuit builder (from existing implementation)
    2. Application of NoisyCircuitFactory
    3. Simulation with MQT-qudits backend
    4. Population calculation from results
    
    Args:
        params: PhysicalParameters instance
        noise_model: MQT-qudits NoiseModel
    
    Returns:
        Dict with placeholder results
    """
    print("\n" + "="*70)
    print("Qudit Noisy Simulation (Placeholder)")
    print("="*70)
    print("\nNote: Full implementation requires integration with:")
    print("  - mqt_qudits_four_molecule_sparse_implementation.py")
    print("  - exact_hamiltonian_builders.py")
    print("  - NoisyCircuitFactory for each Trotter step")
    
    print("\nNoise model configured and ready for use:")
    print(f"  - Basis gates: {noise_model.basis_gates}")
    print(f"  - Subspace noise defined for (0,1), (1,2), (0,2) transitions")
    print(f"  - Applied to local and nonlocal gates")
    
    print("\n" + "="*70)
    
    return {
        'status': 'placeholder',
        'noise_model': noise_model,
        'message': 'Noise model ready - full simulation requires Qudit circuit integration'
    }
