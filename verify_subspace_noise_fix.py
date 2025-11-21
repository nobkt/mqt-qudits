#!/usr/bin/env python3
"""
Verification script for SubspaceNoise fix in MQT-Qudits.

This script verifies that the fix for the AttributeError with SubspaceNoise
is working correctly by testing various scenarios.

Original error:
    AttributeError: 'SubspaceNoise' object has no attribute 'probability_depolarizing'

Fix:
    Updated src/python/bindings.cpp to handle SubspaceNoise objects by:
    1. Checking for 'subspace_w_probs' attribute
    2. Extracting and averaging noise probabilities from all subspaces
    3. Falling back to direct attribute access for regular Noise objects
"""

import sys
import numpy as np
from mqt.qudits.simulation.noise_tools import SubspaceNoise, Noise, NoiseModel
from mqt.qudits.quantum_circuit import QuantumCircuit, QuantumRegister
from mqt.qudits.simulation import MQTQuditProvider


def test_subspace_noise_creation():
    """Test 1: SubspaceNoise object creation."""
    print("Test 1: Creating SubspaceNoise object...")
    try:
        subspace_noise = SubspaceNoise(0.001, 0.001, [(0, 1), (0, 2), (1, 2)])
        print(f"  ✓ SubspaceNoise created: {subspace_noise}")
        return True
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def test_noise_model_with_subspace_noise():
    """Test 2: NoiseModel with SubspaceNoise."""
    print("\nTest 2: Creating NoiseModel with SubspaceNoise...")
    try:
        subspace_noise = SubspaceNoise(0.001, 0.001, [(0, 1), (0, 2), (1, 2)])
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(subspace_noise, ['h', 'x', 'cx'])
        print(f"  ✓ NoiseModel created with basis gates: {noise_model.basis_gates}")
        return True
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def test_circuit_execution_with_subspace_noise():
    """Test 3: Circuit execution with SubspaceNoise."""
    print("\nTest 3: Executing circuit with SubspaceNoise...")
    try:
        # Create SubspaceNoise
        subspace_noise = SubspaceNoise(0.001, 0.001, [(0, 1), (0, 2), (1, 2)])
        
        # Create NoiseModel
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(subspace_noise, ['h', 'x'])
        
        # Create circuit
        qreg = QuantumRegister("test", 1, [3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.x(0)
        
        # Execute
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        
        state_vector = result.get_state_vector()
        counts = result.get_counts()
        
        print(f"  ✓ Circuit executed successfully")
        print(f"    State vector shape: {state_vector.shape}")
        print(f"    Number of shots: {len(counts)}")
        return True
    except AttributeError as e:
        print(f"  ✗ Failed with AttributeError: {e}")
        print(f"    This is the original bug that should be fixed!")
        return False
    except Exception as e:
        print(f"  ✗ Failed with unexpected error: {e}")
        return False


def test_mixed_noise_types():
    """Test 4: Mixed Noise and SubspaceNoise."""
    print("\nTest 4: Using both Noise and SubspaceNoise in same model...")
    try:
        # Create both types
        simple_noise = Noise(0.01, 0.01)
        subspace_noise = SubspaceNoise(0.005, 0.005, [(0, 1), (0, 2), (1, 2)])
        
        # Create NoiseModel
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(simple_noise, ['h'])
        noise_model.add_quantum_error_locally(subspace_noise, ['x', 'cx'])
        
        # Create circuit
        qreg = QuantumRegister("test", 2, [3, 3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.x(1)
        circuit.cx([0, 1], [0, 1, 1, np.pi/2])
        
        # Execute
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        
        print(f"  ✓ Mixed noise types work correctly")
        return True
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def test_four_qutrit_system():
    """Test 5: 4-qutrit system (as in the tutorial notebook)."""
    print("\nTest 5: 4-qutrit system with SubspaceNoise (tutorial scenario)...")
    try:
        # Create SubspaceNoise as in tutorial
        subspace_noise = SubspaceNoise(0.001, 0.001, [(0, 1), (0, 2), (1, 2)])
        
        # Create NoiseModel
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(subspace_noise, 
                                             ['virtrz', 'r', 'rz', 'rh', 'cx', 'h', 'x', 'z', 's'])
        
        # Create 4-qutrit circuit
        qreg = QuantumRegister("molecules", 4, [3, 3, 3, 3])
        circuit = QuantumCircuit(qreg)
        circuit.x(0)  # Set molecule 0 to T1
        circuit.x(3)  # Set molecule 3 to T1
        
        # Execute
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        
        state_vector = result.get_state_vector()
        print(f"  ✓ 4-qutrit system executed successfully")
        print(f"    State vector shape: {state_vector.shape} (expected: (1, 81))")
        return True
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def main():
    """Run all verification tests."""
    print("="*70)
    print("SubspaceNoise Fix Verification")
    print("="*70)
    print()
    print("This script verifies the fix for:")
    print("  AttributeError: 'SubspaceNoise' object has no attribute 'probability_depolarizing'")
    print()
    
    tests = [
        test_subspace_noise_creation,
        test_noise_model_with_subspace_noise,
        test_circuit_execution_with_subspace_noise,
        test_mixed_noise_types,
        test_four_qutrit_system,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"\nUnexpected error in {test.__name__}: {e}")
            results.append(False)
    
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    passed = sum(results)
    total = len(results)
    print(f"\nTests passed: {passed}/{total}")
    
    if passed == total:
        print("\n✓ All tests passed! SubspaceNoise fix is working correctly.")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed. Please review the output above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
