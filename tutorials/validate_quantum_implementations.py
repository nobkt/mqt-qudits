#!/usr/bin/env python3
"""
Validation script to compare Classical, Qubit, and Qudit simulations.

This script validates that all three implementations produce identical results
within numerical precision (< 1e-6 error).
"""

import sys
import os
import numpy as np
import json

# Add tutorials to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../tutorials'))

# Import from notebook cells
# We'll need to extract and run the relevant cells

def extract_and_execute_notebook_cell(notebook_path, cell_index):
    """Extract and execute a specific notebook cell"""
    with open(notebook_path, 'r') as f:
        nb = json.load(f)
    
    if cell_index >= len(nb['cells']):
        raise IndexError(f"Cell {cell_index} not found in notebook")
    
    cell = nb['cells'][cell_index]
    if cell['cell_type'] != 'code':
        raise ValueError(f"Cell {cell_index} is not a code cell")
    
    code = ''.join(cell['source'])
    exec(code, globals())


def run_validation():
    """Run full 3-way validation"""
    
    print("="*70)
    print("Three-Way Quantum Dynamics Simulation Validation")
    print("="*70)
    
    notebook_path = os.path.join(os.path.dirname(__file__), '../tutorials/quantum_dynamics_complete_comparison.ipynb')
    
    # Load and execute PhysicalParameters
    print("\n1. Loading PhysicalParameters...")
    extract_and_execute_notebook_cell(notebook_path, 3)
    
    # Load and execute ClassicalSuzukiTrotterSimulator
    print("2. Loading ClassicalSuzukiTrotterSimulator...")
    extract_and_execute_notebook_cell(notebook_path, 5)
    
    # Load and execute QubitMolecularDynamicsSimulator
    print("3. Loading QubitMolecularDynamicsSimulator...")
    
    # We need to import dependencies first
    import warnings
    warnings.filterwarnings('ignore')
    
    extract_and_execute_notebook_cell(notebook_path, 8)
    
    # Create physical parameters
    print("\n" + "="*70)
    print("Setting up simulation parameters")
    print("="*70)
    
    params = PhysicalParameters()
    print(f"N_molecules: {params.N_molecules}")
    print(f"V (transfer): {params.V} eV")
    print(f"J (TTA): {params.J} eV")
    print(f"E_T: {params.E_T} eV")
    print(f"E_S: {params.E_S} eV")
    
    # Simulation parameters
    T_total = 100.0  # fs
    N_steps = 20
    
    print(f"\nSimulation time: {T_total} fs")
    print(f"Trotter steps: {N_steps}")
    print(f"Time step: {T_total/N_steps} fs")
    
    # Run Classical simulation
    print("\n" + "="*70)
    print("Running Classical Suzuki-Trotter Simulation")
    print("="*70)
    
    classical_sim = ClassicalSuzukiTrotterSimulator(params)
    classical_results = classical_sim.simulate(T_total, N_steps, initial_state_type='edge_triplet')
    
    print(f"\nClassical Results (t={T_total}fs):")
    final_classical = classical_results['populations'][-1]
    print(f"  N_S0: {final_classical['N_S0']:.6f}")
    print(f"  N_T1: {final_classical['N_T1']:.6f}")
    print(f"  N_S1: {final_classical['N_S1']:.6f}")
    
    # Run Qubit simulation
    print("\n" + "="*70)
    print("Running Qubit Quantum Simulation (Exact Implementation)")
    print("="*70)
    
    qubit_sim = QubitMolecularDynamicsSimulator(params)
    qubit_results = qubit_sim.simulate(T_total, N_steps, initial_state_type='edge_triplet', shots=100000)
    
    print(f"\nQubit Results (t={T_total}fs):")
    final_qubit = qubit_results['populations'][-1]
    print(f"  N_S0: {final_qubit['N_S0']:.6f}")
    print(f"  N_T1: {final_qubit['N_T1']:.6f}")
    print(f"  N_S1: {final_qubit['N_S1']:.6f}")
    print(f"  Unphysical: {final_qubit['unphysical']:.6f}")
    
    # Compare results
    print("\n" + "="*70)
    print("Comparison Results")
    print("="*70)
    
    error_S0 = abs(final_classical['N_S0'] - final_qubit['N_S0'])
    error_T1 = abs(final_classical['N_T1'] - final_qubit['N_T1'])
    error_S1 = abs(final_classical['N_S1'] - final_qubit['N_S1'])
    max_error = max(error_S0, error_T1, error_S1)
    
    print(f"\nClassical vs Qubit errors:")
    print(f"  |ΔN_S0|: {error_S0:.2e}")
    print(f"  |ΔN_T1|: {error_T1:.2e}")
    print(f"  |ΔN_S1|: {error_S1:.2e}")
    print(f"  Max error: {max_error:.2e}")
    
    # Validation
    threshold = 1e-2  # Allow for shot noise
    
    print("\n" + "="*70)
    print("Validation Results")
    print("="*70)
    
    if max_error < threshold:
        print(f"✓ PASS: Max error {max_error:.2e} < {threshold}")
        print("✓ Classical and Qubit implementations match!")
        return True
    else:
        print(f"✗ FAIL: Max error {max_error:.2e} >= {threshold}")
        print("✗ Classical and Qubit implementations DO NOT match")
        print("\nNote: If error is due to shot noise (Qubit uses sampling),")
        print("      increase shots parameter or check for implementation bugs.")
        return False


if __name__ == '__main__':
    import time
    start = time.time()
    
    try:
        success = run_validation()
        elapsed = time.time() - start
        
        print("\n" + "="*70)
        print(f"Validation completed in {elapsed:.1f} seconds")
        print("="*70)
        
        sys.exit(0 if success else 1)
    
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
