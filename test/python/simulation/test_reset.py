"""Tests for the Reset instruction and mid-circuit stochastic channel
execution on the TNSim state-vector backend."""

from __future__ import annotations

from unittest import TestCase

import numpy as np

from mqt.qudits.quantum_circuit import QuantumCircuit
from mqt.qudits.quantum_circuit.components.quantum_register import QuantumRegister
from mqt.qudits.quantum_circuit.gates import KrausChannel, Reset
from mqt.qudits.simulation import MQTQuditProvider


def _reset_kraus(d: int) -> list[np.ndarray]:
    ops = []
    for k in range(d):
        m = np.zeros((d, d), dtype=np.complex128)
        m[0, k] = 1.0
        ops.append(m)
    return ops


class TestResetInstruction(TestCase):
    """Reset must be a KrausChannel with Kraus set {|0><k|}."""

    @staticmethod
    def test_reset_is_kraus_channel_with_reset_kraus_set() -> None:
        qreg = QuantumRegister("q", 1, [3])
        circuit = QuantumCircuit(qreg)
        gate = circuit.reset_qudit(0)

        assert isinstance(gate, Reset)
        assert isinstance(gate, KrausChannel)
        assert gate.qasm_tag == "reset"
        assert gate.kraus_rank == 3
        for k_op, ref in zip(gate.kraus_operators, _reset_kraus(3)):
            assert np.array_equal(k_op, ref)

    @staticmethod
    def test_reset_has_no_unitary_matrix() -> None:
        qreg = QuantumRegister("q", 1, [3])
        circuit = QuantumCircuit(qreg)
        gate = circuit.reset_qudit(0)
        try:
            gate.to_matrix()
        except NotImplementedError:
            pass
        else:  # pragma: no cover
            msg = "Reset.to_matrix() must raise NotImplementedError"
            raise AssertionError(msg)


class TestResetOnDMSim(TestCase):
    """DMSim executes Reset deterministically as the exact CPTP channel."""

    @staticmethod
    def test_reset_reinitialises_qudit_to_zero() -> None:
        provider = MQTQuditProvider()
        backend = provider.get_backend("dmsim")
        qreg = QuantumRegister("q", 2, [3, 3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.h(1)
        circuit.reset_qudit(0)

        rho = backend.run(circuit).result().get_density_matrix()
        # Marginal of qudit 0 must be |0><0|; qudit 1 stays H|0>.
        rho_q0 = np.einsum("ikjk->ij", rho.reshape(3, 3, 3, 3))
        assert np.allclose(rho_q0, np.diag([1.0, 0.0, 0.0]), atol=1e-12)
        rho_q1 = np.einsum("kikj->ij", rho.reshape(3, 3, 3, 3))
        h_matrix = QuantumCircuit(QuantumRegister("r", 1, [3])).h(0).to_matrix()
        psi1 = h_matrix[:, 0]
        assert np.allclose(rho_q1, np.outer(psi1, psi1.conj()), atol=1e-12)


class TestResetOnTNSim(TestCase):
    """TNSim executes Reset stochastically (single trajectory per run)."""

    @staticmethod
    def test_reset_collapses_target_to_zero() -> None:
        provider = MQTQuditProvider()
        backend = provider.get_backend("tnsim")
        qreg = QuantumRegister("q", 2, [3, 3])
        circuit = QuantumCircuit(qreg)
        circuit.h(0)
        circuit.h(1)
        circuit.reset_qudit(0)

        psi = backend.run(circuit, seed=11).result().get_state_vector().reshape(3, 3)
        assert np.isclose(np.linalg.norm(psi), 1.0, atol=1e-12)
        # Target qudit must be exactly |0>.
        assert np.allclose(psi[1:], 0.0)

    @staticmethod
    def test_seed_reproducibility() -> None:
        provider = MQTQuditProvider()
        backend = provider.get_backend("tnsim")

        def build() -> QuantumCircuit:
            circuit = QuantumCircuit(QuantumRegister("q", 1, [3]))
            circuit.h(0)
            circuit.reset_qudit(0)
            circuit.h(0)
            return circuit

        s1 = backend.run(build(), seed=7).result().get_state_vector()
        s2 = backend.run(build(), seed=7).result().get_state_vector()
        assert np.allclose(s1, s2)

    @staticmethod
    def test_trajectory_average_matches_dmsim() -> None:
        """Averaging TNSim reset trajectories reproduces the DMSim channel.

        Circuit: entangle two qutrits via H + CSum, then reset qudit 0.
        The trajectory average of |ψ⟩⟨ψ| over seeds must converge to the
        deterministic DMSim density matrix at the O(1/√n) statistical rate.
        """
        provider = MQTQuditProvider()
        tnsim = provider.get_backend("tnsim")
        dmsim = provider.get_backend("dmsim")

        def build() -> QuantumCircuit:
            circuit = QuantumCircuit(QuantumRegister("q", 2, [3, 3]))
            circuit.h(0)
            circuit.csum([0, 1])
            circuit.reset_qudit(0)
            return circuit

        rho_ref = dmsim.run(build()).result().get_density_matrix()

        n_traj = 600
        acc = np.zeros((9, 9), dtype=np.complex128)
        for seed in range(n_traj):
            psi = tnsim.run(build(), seed=seed).result().get_state_vector().reshape(-1)
            acc += np.outer(psi, psi.conj())
        acc /= n_traj

        err = float(np.linalg.norm(acc - rho_ref))
        assert err < 5.0 / np.sqrt(n_traj), f"trajectory average error {err:.3e}"

    @staticmethod
    def test_general_kraus_channel_stochastic_on_tnsim() -> None:
        """A non-reset CPTP channel (qutrit depolarising-like) also runs.

        The trajectory average must match DMSim's deterministic result.
        """
        provider = MQTQuditProvider()
        tnsim = provider.get_backend("tnsim")
        dmsim = provider.get_backend("dmsim")

        p = 0.3
        d = 3
        x_matrix = np.roll(np.eye(d), 1, axis=0).astype(np.complex128)
        kraus = [np.sqrt(1 - p) * np.eye(d, dtype=np.complex128), np.sqrt(p) * x_matrix]

        def build() -> QuantumCircuit:
            circuit = QuantumCircuit(QuantumRegister("q", 1, [d]))
            circuit.h(0)
            circuit.kraus_channel(0, kraus)
            return circuit

        rho_ref = dmsim.run(build()).result().get_density_matrix()

        n_traj = 600
        acc = np.zeros((d, d), dtype=np.complex128)
        for seed in range(n_traj):
            psi = tnsim.run(build(), seed=seed).result().get_state_vector().reshape(-1)
            acc += np.outer(psi, psi.conj())
        acc /= n_traj

        err = float(np.linalg.norm(acc - rho_ref))
        assert err < 5.0 / np.sqrt(n_traj), f"trajectory average error {err:.3e}"
