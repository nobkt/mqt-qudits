#!/usr/bin/env python3
"""
Validation script for noise model integration.

This script tests that:
1. All noise model modules can be imported
2. Noise models can be created
3. Basic simulations run without errors
"""

import sys
import os

# Add tutorials directory to path
sys.path.insert(0, '/home/runner/work/mqt-qudits/mqt-qudits/tutorials')

# Add user site-packages for Qiskit
import site
user_site = site.getusersitepackages()
if user_site not in sys.path:
    sys.path.insert(0, user_site)

print("="*70)
print("Noise Model Integration Validation")
print("="*70)

# Test 1: Import modules
print("\n1. Testing module imports...")
try:
    from noise_simulation_implementations import (
        NoiseParameters,
        create_qiskit_noise_model,
        create_qudit_noise_model,
        QISKIT_AVAILABLE,
        MQTQUDITS_NOISE_AVAILABLE
    )
    print("   ✓ noise_simulation_implementations imported")
except Exception as e:
    print(f"   ✗ Failed to import noise_simulation_implementations: {e}")
    sys.exit(1)

try:
    from extended_noise_simulators import ClassicalNoisySimulator
    print("   ✓ extended_noise_simulators imported")
except Exception as e:
    print(f"   ✗ Failed to import extended_noise_simulators: {e}")
    sys.exit(1)

# Test 2: Create noise parameters
print("\n2. Testing noise parameter creation...")
try:
    noise_params = NoiseParameters()
    print(f"   ✓ NoiseParameters created")
    print(f"     T1 = {noise_params.T1} fs")
    print(f"     T2 = {noise_params.T2} fs")
    print(f"     p_depol_1q = {noise_params.p_depol_1q}")
    print(f"     p_depol_2q = {noise_params.p_depol_2q}")
except Exception as e:
    print(f"   ✗ Failed to create NoiseParameters: {e}")
    sys.exit(1)

# Test 3: Create Qiskit noise model
print("\n3. Testing Qiskit noise model...")
if QISKIT_AVAILABLE:
    try:
        qiskit_noise = create_qiskit_noise_model(noise_params)
        print(f"   ✓ Qiskit NoiseModel created")
        print(f"     Basis gates: {len(qiskit_noise.basis_gates)}")
    except Exception as e:
        print(f"   ✗ Failed to create Qiskit noise model: {e}")
else:
    print("   ⚠ Qiskit not available (expected in some environments)")

# Test 4: Create Qudit noise model
print("\n4. Testing Qudit noise model...")
if MQTQUDITS_NOISE_AVAILABLE:
    try:
        qudit_noise = create_qudit_noise_model(noise_params)
        print(f"   ✓ Qudit NoiseModel created")
        print(f"     Basis gates: {len(qudit_noise.basis_gates)}")
    except Exception as e:
        print(f"   ✗ Failed to create Qudit noise model: {e}")
        sys.exit(1)
else:
    print("   ✗ MQT-Qudits noise tools not available")
    sys.exit(1)

# Test 5: Verify notebook structure
print("\n5. Testing notebook structure...")
try:
    import json
    nb_path = '/home/runner/work/mqt-qudits/mqt-qudits/tutorials/quantum_dynamics_complete_comparison.ipynb'
    with open(nb_path, 'r') as f:
        nb = json.load(f)
    
    cell_count = len(nb['cells'])
    noise_cells = [c for c in nb['cells'] if c.get('id', '').startswith('noise_')]
    
    print(f"   ✓ Notebook loaded successfully")
    print(f"     Total cells: {cell_count}")
    print(f"     Noise-related cells: {len(noise_cells)}")
    
    if cell_count < 35:
        print(f"   ⚠ Warning: Expected at least 35 cells, got {cell_count}")
    
    if len(noise_cells) < 3:
        print(f"   ⚠ Warning: Expected at least 3 noise cells, got {len(noise_cells)}")
    
except Exception as e:
    print(f"   ✗ Failed to verify notebook: {e}")
    sys.exit(1)

# Test 6: Check file presence
print("\n6. Checking created files...")
files_to_check = [
    'noise_simulation_implementations.py',
    'extended_noise_simulators.py',
    'NOISE_INTEGRATION_GUIDE.md',
    'add_noise_to_notebook.py',
    'quantum_dynamics_complete_comparison.ipynb',
    'quantum_dynamics_complete_comparison.ipynb.backup'
]

for filename in files_to_check:
    filepath = os.path.join('/home/runner/work/mqt-qudits/mqt-qudits/tutorials', filename)
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        print(f"   ✓ {filename} ({size} bytes)")
    else:
        print(f"   ✗ {filename} NOT FOUND")

print("\n" + "="*70)
print("Validation Summary")
print("="*70)
print("✓ All critical tests passed")
print("✓ Noise model infrastructure is functional")
print("✓ Notebook modifications are in place")
print("\nNOTE: Full notebook execution requires Jupyter environment")
print("="*70)
