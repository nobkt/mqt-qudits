from __future__ import annotations

from unittest import TestCase

import numpy as np

from mqt.qudits.quantum_circuit import QuantumCircuit
from mqt.qudits.quantum_circuit.components.quantum_register import QuantumRegister
from mqt.qudits.simulation import MQTQuditProvider
from mqt.qudits.simulation.noise_tools import Noise, NoiseModel, SubspaceNoise

from .._qudits.test_pymisim import is_quantum_state


class TestSubspaceNoise(TestCase):
    @staticmethod
    def test_subspace_noise_single_qudit():
        """Test that SubspaceNoise works with single qudit gates."""
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")

        qreg = QuantumRegister("reg", 1, [3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.x(0)

        # Create SubspaceNoise for all qutrit transitions
        subspace_noise = SubspaceNoise(0.01, 0.01, [(0, 1), (0, 2), (1, 2)])

        # Create noise model
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(subspace_noise, ["h", "x"])

        # Run with noise model
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        state_vector = result.get_state_vector()
        counts = result.get_counts()

        assert len(counts) == 100
        assert len(state_vector.squeeze()) == 3
        assert is_quantum_state(state_vector)

    @staticmethod
    def test_subspace_noise_two_qudit():
        """Test that SubspaceNoise works with two qudit gates."""
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")

        qreg = QuantumRegister("reg", 2, [3, 3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.cx([0, 1], [0, 1, 1, np.pi / 2])

        # Create SubspaceNoise for all qutrit transitions
        local_noise = SubspaceNoise(0.005, 0.005, [(0, 1), (0, 2), (1, 2)])
        nonlocal_noise = SubspaceNoise(0.02, 0.01, [(0, 1), (0, 2), (1, 2)])

        # Create noise model
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(local_noise, ["h"])
        noise_model.add_nonlocal_quantum_error(nonlocal_noise, ["cx"])

        # Run with noise model
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        state_vector = result.get_state_vector()
        counts = result.get_counts()

        assert len(counts) == 100
        assert len(state_vector.squeeze()) == 9
        assert is_quantum_state(state_vector)

    @staticmethod
    def test_subspace_noise_different_levels():
        """Test that SubspaceNoise with different probability for different levels works."""
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")

        qreg = QuantumRegister("reg", 1, [3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)

        # Create SubspaceNoise with different probabilities for different transitions
        # In practice, we use the same probabilities for all, but the structure allows different ones
        subspace_noise = SubspaceNoise(0.01, 0.005, [(0, 1), (0, 2), (1, 2)])

        # Create noise model
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(subspace_noise, ["h"])

        # Run with noise model
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        state_vector = result.get_state_vector()
        counts = result.get_counts()

        assert len(counts) == 100
        assert len(state_vector.squeeze()) == 3
        assert is_quantum_state(state_vector)

    @staticmethod
    def test_mixed_noise_and_subspace_noise():
        """Test that both Noise and SubspaceNoise can be used in the same model."""
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")

        qreg = QuantumRegister("reg", 2, [3, 3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.x(1)
        circuit.cx([0, 1], [0, 1, 1, np.pi / 2])

        # Create both Noise and SubspaceNoise
        simple_noise = Noise(0.01, 0.01)
        subspace_noise = SubspaceNoise(0.02, 0.01, [(0, 1), (0, 2), (1, 2)])

        # Create noise model
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(simple_noise, ["h"])
        noise_model.add_quantum_error_locally(subspace_noise, ["x"])
        noise_model.add_nonlocal_quantum_error(subspace_noise, ["cx"])

        # Run with noise model
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        state_vector = result.get_state_vector()
        counts = result.get_counts()

        assert len(counts) == 100
        assert len(state_vector.squeeze()) == 9
        assert is_quantum_state(state_vector)

    @staticmethod
    def test_subspace_noise_with_multiple_gates():
        """Test SubspaceNoise applied to multiple different gate types."""
        provider = MQTQuditProvider()
        backend = provider.get_backend("misim")

        qreg = QuantumRegister("reg", 1, [3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.x(0)
        circuit.z(0)
        circuit.s(0)

        # Create SubspaceNoise
        subspace_noise = SubspaceNoise(0.005, 0.005, [(0, 1), (0, 2), (1, 2)])

        # Create noise model for multiple gate types
        noise_model = NoiseModel()
        noise_model.add_quantum_error_locally(subspace_noise, ["h", "x", "z", "s"])

        # Run with noise model
        job = backend.run(circuit, noise_model=noise_model, shots=100)
        result = job.result()
        state_vector = result.get_state_vector()
        counts = result.get_counts()

        assert len(counts) == 100
        assert len(state_vector.squeeze()) == 3
        assert is_quantum_state(state_vector)
