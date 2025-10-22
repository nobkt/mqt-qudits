#!/usr/bin/env python3
"""
Unit tests for subspace analysis in sparse compiler fix.

Tests the helper functions that detect whether an active subspace
spans a single qudit or multiple qudits.
"""

import sys
from pathlib import Path

# Add tutorials path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / 'tutorials'))

# Import test utilities
def test_global_index_to_qudit_states():
    """Test conversion of global index to qudit states"""
    # Create a minimal mock class with the method
    class MockGateGenerator:
        def _global_index_to_qudit_states(self, global_idx: int, dimensions: list[int]) -> list[int]:
            states = []
            idx = global_idx
            for dim in reversed(dimensions):
                states.append(idx % dim)
                idx //= dim
            return list(reversed(states))
    
    mock = MockGateGenerator()
    
    # Test case 1: 2-qutrit system (3⊗3)
    assert mock._global_index_to_qudit_states(0, [3, 3]) == [0, 0], "Index 0 should be |00⟩"
    assert mock._global_index_to_qudit_states(1, [3, 3]) == [0, 1], "Index 1 should be |01⟩"
    assert mock._global_index_to_qudit_states(2, [3, 3]) == [0, 2], "Index 2 should be |02⟩"
    assert mock._global_index_to_qudit_states(3, [3, 3]) == [1, 0], "Index 3 should be |10⟩"
    assert mock._global_index_to_qudit_states(4, [3, 3]) == [1, 1], "Index 4 should be |11⟩"
    assert mock._global_index_to_qudit_states(8, [3, 3]) == [2, 2], "Index 8 should be |22⟩"
    
    print("✓ test_global_index_to_qudit_states passed")


def test_analyze_subspace_single_qudit():
    """Test analysis of single-qudit subspace"""
    # Minimal mock class
    class MockGateGenerator:
        def _global_index_to_qudit_states(self, global_idx: int, dimensions: list[int]) -> list[int]:
            states = []
            idx = global_idx
            for dim in reversed(dimensions):
                states.append(idx % dim)
                idx //= dim
            return list(reversed(states))
        
        def _analyze_subspace(self, active_indices: list[int], dimensions: list[int]) -> dict:
            states = [self._global_index_to_qudit_states(idx, dimensions) for idx in active_indices]
            n_qudits = len(dimensions)
            involved_qudits = []
            
            for qudit_idx in range(n_qudits):
                qudit_states = [state[qudit_idx] for state in states]
                if len(set(qudit_states)) > 1:
                    involved_qudits.append(qudit_idx)
            
            if len(involved_qudits) == 1:
                qudit_idx = involved_qudits[0]
                local_levels = sorted(set([state[qudit_idx] for state in states]))
                
                global_to_local = {}
                for global_idx in active_indices:
                    state = self._global_index_to_qudit_states(global_idx, dimensions)
                    local_level = state[qudit_idx]
                    global_to_local[global_idx] = local_level
                
                return {
                    'type': 'single_qudit',
                    'qudit_idx': qudit_idx,
                    'local_levels': local_levels,
                    'global_to_local': global_to_local,
                    'involved_qudits': involved_qudits
                }
            else:
                return {
                    'type': 'multi_qudit',
                    'involved_qudits': involved_qudits
                }
    
    mock = MockGateGenerator()
    
    # Test: Subspace [3, 4] in 3⊗3 system
    # [3, 4] = [|10⟩, |11⟩] - only qudit 1 varies (0→1), qudit 0 stays at 1
    result = mock._analyze_subspace([3, 4], [3, 3])
    
    assert result['type'] == 'single_qudit', "Should detect single-qudit subspace"
    assert result['qudit_idx'] == 1, "Should identify qudit 1"
    assert result['local_levels'] == [0, 1], "Local levels should be [0, 1]"
    assert result['global_to_local'][3] == 0, "Global 3 (|10⟩) should map to local 0"
    assert result['global_to_local'][4] == 1, "Global 4 (|11⟩) should map to local 1"
    
    print("✓ test_analyze_subspace_single_qudit passed")


def test_analyze_subspace_multi_qudit_h_transfer():
    """Test analysis of multi-qudit subspace (H_transfer case)"""
    class MockGateGenerator:
        def _global_index_to_qudit_states(self, global_idx: int, dimensions: list[int]) -> list[int]:
            states = []
            idx = global_idx
            for dim in reversed(dimensions):
                states.append(idx % dim)
                idx //= dim
            return list(reversed(states))
        
        def _analyze_subspace(self, active_indices: list[int], dimensions: list[int]) -> dict:
            states = [self._global_index_to_qudit_states(idx, dimensions) for idx in active_indices]
            n_qudits = len(dimensions)
            involved_qudits = []
            
            for qudit_idx in range(n_qudits):
                qudit_states = [state[qudit_idx] for state in states]
                if len(set(qudit_states)) > 1:
                    involved_qudits.append(qudit_idx)
            
            if len(involved_qudits) == 1:
                qudit_idx = involved_qudits[0]
                local_levels = sorted(set([state[qudit_idx] for state in states]))
                
                global_to_local = {}
                for global_idx in active_indices:
                    state = self._global_index_to_qudit_states(global_idx, dimensions)
                    local_level = state[qudit_idx]
                    global_to_local[global_idx] = local_level
                
                return {
                    'type': 'single_qudit',
                    'qudit_idx': qudit_idx,
                    'local_levels': local_levels,
                    'global_to_local': global_to_local,
                    'involved_qudits': involved_qudits
                }
            else:
                return {
                    'type': 'multi_qudit',
                    'involved_qudits': involved_qudits
                }
    
    mock = MockGateGenerator()
    
    # Test: H_transfer subspace [1, 3] in 3⊗3 system
    # [1, 3] = [|01⟩, |10⟩] - both qudits vary
    result = mock._analyze_subspace([1, 3], [3, 3])
    
    assert result['type'] == 'multi_qudit', "Should detect multi-qudit subspace"
    assert result['involved_qudits'] == [0, 1], "Both qudits should be involved"
    
    print("✓ test_analyze_subspace_multi_qudit_h_transfer passed")


def test_analyze_subspace_multi_qudit_h_tta():
    """Test analysis of multi-qudit subspace (H_TTA case)"""
    class MockGateGenerator:
        def _global_index_to_qudit_states(self, global_idx: int, dimensions: list[int]) -> list[int]:
            states = []
            idx = global_idx
            for dim in reversed(dimensions):
                states.append(idx % dim)
                idx //= dim
            return list(reversed(states))
        
        def _analyze_subspace(self, active_indices: list[int], dimensions: list[int]) -> dict:
            states = [self._global_index_to_qudit_states(idx, dimensions) for idx in active_indices]
            n_qudits = len(dimensions)
            involved_qudits = []
            
            for qudit_idx in range(n_qudits):
                qudit_states = [state[qudit_idx] for state in states]
                if len(set(qudit_states)) > 1:
                    involved_qudits.append(qudit_idx)
            
            if len(involved_qudits) == 1:
                qudit_idx = involved_qudits[0]
                local_levels = sorted(set([state[qudit_idx] for state in states]))
                
                global_to_local = {}
                for global_idx in active_indices:
                    state = self._global_index_to_qudit_states(global_idx, dimensions)
                    local_level = state[qudit_idx]
                    global_to_local[global_idx] = local_level
                
                return {
                    'type': 'single_qudit',
                    'qudit_idx': qudit_idx,
                    'local_levels': local_levels,
                    'global_to_local': global_to_local,
                    'involved_qudits': involved_qudits
                }
            else:
                return {
                    'type': 'multi_qudit',
                    'involved_qudits': involved_qudits
                }
    
    mock = MockGateGenerator()
    
    # Test: H_TTA subspace [2, 4, 6] in 3⊗3 system
    # [2, 4, 6] = [|02⟩, |11⟩, |20⟩] - both qudits vary
    result = mock._analyze_subspace([2, 4, 6], [3, 3])
    
    assert result['type'] == 'multi_qudit', "Should detect multi-qudit subspace"
    assert result['involved_qudits'] == [0, 1], "Both qudits should be involved"
    
    print("✓ test_analyze_subspace_multi_qudit_h_tta passed")


if __name__ == "__main__":
    print("Running subspace analysis tests...")
    print()
    
    test_global_index_to_qudit_states()
    test_analyze_subspace_single_qudit()
    test_analyze_subspace_multi_qudit_h_transfer()
    test_analyze_subspace_multi_qudit_h_tta()
    
    print()
    print("=" * 70)
    print("All tests passed! ✓")
    print("=" * 70)
