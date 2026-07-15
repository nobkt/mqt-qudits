"""Tests for the density-matrix backend (DMSim) and KrausChannel instruction."""

from __future__ import annotations

from unittest import TestCase

import numpy as np
from scipy.linalg import expm

from mqt.qudits.quantum_circuit import QuantumCircuit
from mqt.qudits.quantum_circuit.components.quantum_register import QuantumRegister
from mqt.qudits.quantum_circuit.gates import KrausChannel
from mqt.qudits.simulation import MQTQuditProvider
from mqt.qudits.simulation.backends.dmsim import (
    apply_kraus_to_density,
    apply_unitary_to_density,
)


def _dm_from_ket(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


class TestDMSimUnitary(TestCase):
    """Unitary gates on DMSim must reproduce U |0…0⟩⟨0…0| U†."""

    @staticmethod
    def test_single_qutrit_h_on_dmsim() -> None:
        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        d = 3
        qreg = QuantumRegister("reg", 1, [d])
        circuit = QuantumCircuit(qreg)
        h = circuit.h(0)

        zero_state = np.zeros(d, dtype=np.complex128)
        zero_state[0] = 1.0
        u = h.to_matrix()
        rho_ref = u @ _dm_from_ket(zero_state) @ u.conj().T

        job = backend.run(circuit)
        rho = job.result().get_density_matrix()

        assert rho.shape == (d, d)
        assert np.allclose(rho, rho_ref, atol=1e-12)

    @staticmethod
    def test_two_qutrit_unitary_via_cu_two() -> None:
        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        d = 3
        qreg = QuantumRegister("reg", 2, [d, d])
        circuit = QuantumCircuit(qreg)
        rng = np.random.default_rng(42)
        # Build an arbitrary 9x9 unitary
        a = rng.standard_normal((d * d, d * d)) + 1j * rng.standard_normal((d * d, d * d))
        q, _ = np.linalg.qr(a)
        u9 = q.astype(np.complex128)
        circuit.cu_two([0, 1], u9)

        zero = np.zeros(d * d, dtype=np.complex128)
        zero[0] = 1.0
        rho_ref = u9 @ _dm_from_ket(zero) @ u9.conj().T

        job = backend.run(circuit)
        rho = job.result().get_density_matrix()

        assert rho.shape == (d * d, d * d)
        assert np.allclose(rho, rho_ref, atol=1e-12)

    @staticmethod
    def test_initial_density_matrix_is_used() -> None:
        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        d = 3
        qreg = QuantumRegister("reg", 1, [d])
        circuit = QuantumCircuit(qreg)
        # Empty circuit, with custom initial ρ
        rho_in = np.array(
            [[0.5, 0.1, 0.0], [0.1, 0.5, 0.0], [0.0, 0.0, 0.0]], dtype=np.complex128
        )
        job = backend.run(circuit, initial_density_matrix=rho_in)
        rho_out = job.result().get_density_matrix()
        assert np.allclose(rho_out, rho_in, atol=1e-14)


class TestKrausChannel(TestCase):
    """KrausChannel construction-time validation and DMSim execution."""

    @staticmethod
    def test_completeness_check_rejects_non_cptp() -> None:
        d = 2
        qreg = QuantumRegister("reg", 1, [d])
        circuit = QuantumCircuit(qreg)
        # K = 0.5 * I is NOT trace-preserving (Σ K† K = 0.25 I ≠ I).
        bad_k = 0.5 * np.eye(d, dtype=np.complex128)
        try:
            KrausChannel(circuit, "BadKraus", 0, [bad_k], d)
        except ValueError as exc:
            assert "trace-preserving" in str(exc)
        else:
            msg = "Expected ValueError for non-CPTP Kraus operator"
            raise AssertionError(msg)

    @staticmethod
    def test_to_matrix_raises() -> None:
        d = 3
        qreg = QuantumRegister("reg", 1, [d])
        circuit = QuantumCircuit(qreg)
        ident = [np.eye(d, dtype=np.complex128)]
        ch = KrausChannel(circuit, "Identity", 0, ident, d)
        try:
            ch.to_matrix(identities=0)
        except NotImplementedError:
            pass
        else:
            msg = "Expected NotImplementedError from KrausChannel.to_matrix"
            raise AssertionError(msg)

    @staticmethod
    def test_amplitude_damping_qubit_against_analytic() -> None:
        """Standard qubit amplitude-damping channel: compare DMSim with analytic ρ.

        K_0 = [[1, 0], [0, sqrt(1-γ)]]
        K_1 = [[0, sqrt(γ)], [0, 0]]
        Acts on initial ρ = |+⟩⟨+|.
        """
        d = 2
        gamma = 0.3
        k0 = np.array([[1.0, 0.0], [0.0, np.sqrt(1.0 - gamma)]], dtype=np.complex128)
        k1 = np.array([[0.0, np.sqrt(gamma)], [0.0, 0.0]], dtype=np.complex128)

        plus = np.array([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0)
        rho_in = _dm_from_ket(plus)
        rho_ref = k0 @ rho_in @ k0.conj().T + k1 @ rho_in @ k1.conj().T

        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        qreg = QuantumRegister("reg", 1, [d])
        circuit = QuantumCircuit(qreg)
        circuit.kraus_channel(0, [k0, k1])

        job = backend.run(circuit, initial_density_matrix=rho_in)
        rho_out = job.result().get_density_matrix()

        assert np.allclose(rho_out, rho_ref, atol=1e-12)
        # Trace preservation
        assert abs(np.trace(rho_out).real - 1.0) < 1e-12

    @staticmethod
    def test_kraus_channel_on_subset_of_qudits() -> None:
        """Apply qubit amplitude damping on qudit index 1 of a 2-qubit register."""
        d = 2
        gamma = 0.4
        k0 = np.array([[1.0, 0.0], [0.0, np.sqrt(1.0 - gamma)]], dtype=np.complex128)
        k1 = np.array([[0.0, np.sqrt(gamma)], [0.0, 0.0]], dtype=np.complex128)

        # Initial state: |1⟩ ⊗ |1⟩  (d=2)
        psi = np.zeros(4, dtype=np.complex128)
        psi[3] = 1.0  # |1⟩|1⟩ → index 1*2+1 = 3
        rho_in = _dm_from_ket(psi)

        # Reference: identity on qudit 0, amplitude damping on qudit 1
        i_d = np.eye(d, dtype=np.complex128)
        kraus_emb = [np.kron(i_d, k0), np.kron(i_d, k1)]
        rho_ref = sum(k @ rho_in @ k.conj().T for k in kraus_emb)

        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        qreg = QuantumRegister("reg", 2, [d, d])
        circuit = QuantumCircuit(qreg)
        circuit.kraus_channel(1, [k0, k1])

        job = backend.run(circuit, initial_density_matrix=rho_in)
        rho_out = job.result().get_density_matrix()

        assert np.allclose(rho_out, rho_ref, atol=1e-12)


class TestDMSimMixedCircuit(TestCase):
    """End-to-end: unitary + Kraus + unitary on a 2-qutrit circuit."""

    @staticmethod
    def test_unitary_then_kraus_then_unitary() -> None:
        d = 3
        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        qreg = QuantumRegister("reg", 2, [d, d])
        circuit = QuantumCircuit(qreg)

        rng = np.random.default_rng(7)
        # Random 9x9 unitary
        a = rng.standard_normal((9, 9)) + 1j * rng.standard_normal((9, 9))
        u9, _ = np.linalg.qr(a)
        u9 = u9.astype(np.complex128)
        # Random 3x3 unitary
        b = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
        u3, _ = np.linalg.qr(b)
        # Trivial CPTP set: identity (compose with unitaries to test the
        # plumbing through DMSim with a non-empty Kraus-channel instruction).
        kraus = [np.eye(d, dtype=np.complex128)]

        circuit.cu_two([0, 1], u9)
        circuit.kraus_channel(0, kraus)
        circuit.cu_one(1, u3)

        # Reference computation
        psi0 = np.zeros(9, dtype=np.complex128)
        psi0[0] = 1.0
        rho = _dm_from_ket(psi0)
        rho = u9 @ rho @ u9.conj().T
        rho = sum(np.kron(k, np.eye(d)) @ rho @ np.kron(k, np.eye(d)).conj().T for k in kraus)
        u3_emb = np.kron(np.eye(d), u3)
        rho = u3_emb @ rho @ u3_emb.conj().T

        job = backend.run(circuit)
        rho_out = job.result().get_density_matrix()

        assert np.allclose(rho_out, rho, atol=1e-12)
        assert abs(np.trace(rho_out).real - 1.0) < 1e-12


class TestApplyHelpers(TestCase):
    """Direct unit tests for apply_unitary_to_density / apply_kraus_to_density."""

    @staticmethod
    def test_apply_unitary_matches_kron_embedding() -> None:
        rng = np.random.default_rng(123)
        dims = [3, 3, 3]
        d_total = 27
        rho = rng.standard_normal((d_total, d_total)) + 1j * rng.standard_normal(
            (d_total, d_total)
        )
        rho = rho @ rho.conj().T  # PSD
        rho = rho / np.trace(rho)

        u3 = expm(1j * rng.standard_normal((3, 3)))  # arbitrary unitary on 1 qutrit
        u3 = u3 @ u3.conj().T  # not unitary, fix
        # Build a proper unitary
        a = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        u3, _ = np.linalg.qr(a)

        out = apply_unitary_to_density(rho, u3, [1], dims)

        # Reference: U embedded as I ⊗ U ⊗ I
        u_emb = np.kron(np.kron(np.eye(3), u3), np.eye(3))
        ref = u_emb @ rho @ u_emb.conj().T

        assert np.allclose(out, ref, atol=1e-12)
