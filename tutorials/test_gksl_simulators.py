"""Comprehensive tests for GKSL-Lindblad simulation modules."""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pytest

from gksl_physical_parameters import GKSLPhysicalParameters
from gksl_math_utils import (
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    build_lindblad_operators,
    vectorize_density_matrix,
    unvectorize_density_matrix,
    compute_von_neumann_entropy,
    compute_purity,
    compute_populations_from_density_matrix,
)
from gksl_validation import (
    PhysicsViolationError,
    validate_density_matrix,
    validate_particle_conservation,
)
from stinespring_utils import (
    stinespring_unitary_from_lindblad,
    apply_stinespring_to_density_matrix,
)
from classical_gksl_simulator import ClassicalGKSLSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_simulator import QuditGKSLSimulator


# ---------------------------------------------------------------------------
# 1. TestGKSLPhysicalParameters
# ---------------------------------------------------------------------------
class TestGKSLPhysicalParameters:
    def test_default_initialization(self):
        p = GKSLPhysicalParameters()
        assert p.E_T == 1.5
        assert p.E_S == 3.0
        assert p.V == 0.1
        assert p.gamma_TTA == 0.05
        assert p.Gamma_fl == 0.01
        assert p.Gamma_ph == 1e-6
        assert p.k_IC == 0.005
        assert p.k_ISC_ST == 0.003
        assert p.k_ISC_TS == 1e-5
        assert p.N_molecules == 4
        assert p.d == 3

    def test_validate_valid(self):
        p = GKSLPhysicalParameters()
        errors = p.validate()
        assert errors == []

    def test_validate_energy_ratio(self):
        p = GKSLPhysicalParameters(E_T=5.0, E_S=3.0)
        errors = p.validate()
        assert any("energy" in e.lower() or "E_T" in e or "E_S" in e for e in errors)

    def test_validate_hierarchy(self):
        p = GKSLPhysicalParameters(Gamma_fl=1e-8, Gamma_ph=1.0)
        errors = p.validate()
        assert len(errors) > 0

    def test_hilbert_space_dim(self):
        p = GKSLPhysicalParameters()
        assert p.get_hilbert_space_dim() == 81  # 3^4

    def test_hilbert_space_dim_boson(self):
        p = GKSLPhysicalParameters(with_boson=True, n_max=2)
        dim = p.get_hilbert_space_dim()
        assert dim == 3**4 * 3**4  # 81 * 81 = 6561

    def test_to_dict(self):
        p = GKSLPhysicalParameters()
        d = p.to_dict()
        assert isinstance(d, dict)
        assert d["E_T"] == 1.5
        assert d["N_molecules"] == 4

    def test_neighbors(self):
        p = GKSLPhysicalParameters()
        nbrs = p.neighbors
        assert nbrs == [(0, 1), (1, 2), (2, 3)]


# ---------------------------------------------------------------------------
# 2. TestMathUtils
# ---------------------------------------------------------------------------
class TestMathUtils:
    @pytest.fixture()
    def params(self):
        return GKSLPhysicalParameters()

    def test_hamiltonian_hermiticity(self, params):
        H0 = build_onsite_hamiltonian(params)
        Ht = build_transfer_hamiltonian(params)
        assert np.allclose(H0, H0.conj().T, atol=1e-12)
        assert np.allclose(Ht, Ht.conj().T, atol=1e-12)

    def test_hamiltonian_eigenvalues(self, params):
        H0 = build_onsite_hamiltonian(params)
        eigs = np.linalg.eigvalsh(H0)
        assert abs(eigs.min()) < 1e-10
        assert abs(eigs.max() - 4 * params.E_S) < 1e-10

    def test_lindblad_count(self, params):
        ops = build_lindblad_operators(params)
        assert len(ops) == 26

    def test_vectorization_roundtrip(self):
        dim = 4
        rho = np.random.randn(dim, dim) + 1j * np.random.randn(dim, dim)
        vec = vectorize_density_matrix(rho)
        rho2 = unvectorize_density_matrix(vec, dim)
        assert np.allclose(rho, rho2)

    def test_entropy_pure_state(self):
        dim = 4
        psi = np.zeros(dim, dtype=complex)
        psi[0] = 1.0
        rho = np.outer(psi, psi.conj())
        S = compute_von_neumann_entropy(rho)
        assert abs(S) < 1e-10

    def test_entropy_mixed_state(self):
        dim = 81
        rho = np.eye(dim, dtype=complex) / dim
        S = compute_von_neumann_entropy(rho)
        assert abs(S - np.log(dim)) < 1e-10

    def test_purity_pure_state(self):
        dim = 81
        psi = np.zeros(dim, dtype=complex)
        psi[0] = 1.0
        rho = np.outer(psi, psi.conj())
        P = compute_purity(rho)
        assert abs(P - 1.0) < 1e-10

    def test_populations_conservation(self, params):
        dim = params.get_hilbert_space_dim()
        psi = np.zeros(dim, dtype=complex)
        psi[28] = 1.0  # |1,0,0,1> edge_triplet
        rho = np.outer(psi, psi.conj())
        pops = compute_populations_from_density_matrix(rho, params)
        total = pops["N_S0"] + pops["N_T1"] + pops["N_S1"]
        assert abs(total - params.N_molecules) < 1e-10


# ---------------------------------------------------------------------------
# 3. TestClassicalGKSLSimulator
# ---------------------------------------------------------------------------
class TestClassicalGKSLSimulator:
    @pytest.fixture()
    def sim(self):
        return ClassicalGKSLSimulator(GKSLPhysicalParameters())

    def test_initialization(self, sim):
        assert sim.H_0.shape == (81, 81)
        assert sim.H_transfer.shape == (81, 81)
        assert sim.H_total.shape == (81, 81)

    def test_trace_preservation(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-8

    def test_particle_conservation(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for pop in result["populations"]:
            total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
            assert abs(total - 4.0) < 1e-6

    def test_initial_populations(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        pop0 = result["populations"][0]
        assert abs(pop0["N_T1"] - 2.0) < 1e-10
        assert abs(pop0["N_S0"] - 2.0) < 1e-10
        assert abs(pop0["N_S1"] - 0.0) < 1e-10

    def test_entropy_non_negative(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for S in result["entropy"]:
            assert S >= -1e-10

    def test_physical_evolution(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        pops = result["populations"]
        # Triplets should decay
        assert pops[-1]["N_T1"] < pops[0]["N_T1"]
        # Ground state should increase
        assert pops[-1]["N_S0"] > pops[0]["N_S0"]


# ---------------------------------------------------------------------------
# 4. TestQubitGKSLSimulator
# ---------------------------------------------------------------------------
class TestQubitGKSLSimulator:
    @pytest.fixture()
    def sim(self):
        return QubitGKSLSimulator(GKSLPhysicalParameters())

    def test_initialization(self, sim):
        assert sim.n_sys_qubits == 8
        assert sim.n_ancilla == 26
        assert sim.n_total_qubits == 34

    def test_trace_preservation(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-4

    def test_particle_conservation(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for pop in result["populations"]:
            total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
            assert abs(total - 4.0) < 0.1


# ---------------------------------------------------------------------------
# 5. TestQuditGKSLSimulator
# ---------------------------------------------------------------------------
class TestQuditGKSLSimulator:
    @pytest.fixture()
    def sim(self):
        return QuditGKSLSimulator(GKSLPhysicalParameters())

    def test_initialization(self, sim):
        assert sim.n_system_qudits == 4
        assert sim.n_ancilla_qubits == 26

    def test_trace_preservation(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-4

    def test_particle_conservation(self, sim):
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for pop in result["populations"]:
            total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
            assert abs(total - 4.0) < 0.1

    def test_qubit_qudit_consistency(self):
        params = GKSLPhysicalParameters()
        q2 = QubitGKSLSimulator(params)
        qd = QuditGKSLSimulator(params)
        r2 = q2.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        rd = qd.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for p2, pd in zip(r2["populations"], rd["populations"]):
            assert abs(p2["N_T1"] - pd["N_T1"]) < 0.2
            assert abs(p2["N_S0"] - pd["N_S0"]) < 0.2


# ---------------------------------------------------------------------------
# 6. TestBosonSimulators
# ---------------------------------------------------------------------------
class TestBosonSimulators:
    @pytest.fixture()
    def boson_params(self):
        return GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1
        )

    def test_classical_boson_small(self, boson_params):
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator

        sim = ClassicalGKSLBosonSimulator(boson_params)
        result = sim.simulate(t_max=2.0, n_steps=3, initial_state="edge_triplet")
        assert len(result["times"]) >= 2

    def test_trace_preservation_boson(self, boson_params):
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator

        sim = ClassicalGKSLBosonSimulator(boson_params)
        result = sim.simulate(t_max=2.0, n_steps=3, initial_state="edge_triplet")
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-6

    def test_qubit_qudit_boson_match(self, boson_params):
        from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator

        q2 = QubitGKSLBosonSimulator(boson_params)
        qd = QuditGKSLBosonSimulator(boson_params)
        r2 = q2.simulate(t_max=2.0, n_steps=3, initial_state="edge_triplet")
        rd = qd.simulate(t_max=2.0, n_steps=3, initial_state="edge_triplet")
        for p2, pd in zip(r2["populations"], rd["populations"]):
            assert abs(p2["N_T1"] - pd["N_T1"]) < 0.3


# ---------------------------------------------------------------------------
# 7. TestValidation
# ---------------------------------------------------------------------------
class TestValidation:
    def _make_valid_rho(self, dim=4):
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        return np.outer(psi, psi.conj())

    def test_valid_density_matrix(self):
        rho = self._make_valid_rho(8)
        result = validate_density_matrix(rho)
        assert result["valid"]

    def test_invalid_trace(self):
        rho = np.eye(4, dtype=complex) * 0.5  # trace = 2
        with pytest.raises(PhysicsViolationError, match="[Tt]race"):
            validate_density_matrix(rho)

    def test_particle_conservation_valid(self):
        pops = {"N_S0": 2.0, "N_T1": 1.5, "N_S1": 0.5}
        assert validate_particle_conservation(pops, N_molecules=4)

    def test_particle_conservation_invalid(self):
        pops = {"N_S0": 2.0, "N_T1": 1.5, "N_S1": 1.5}
        with pytest.raises(PhysicsViolationError):
            validate_particle_conservation(pops, N_molecules=4)


# ---------------------------------------------------------------------------
# 8. TestStinespring
# ---------------------------------------------------------------------------
class TestStinespring:
    @pytest.fixture()
    def setup(self):
        params = GKSLPhysicalParameters()
        ops = build_lindblad_operators(params)
        L, _gamma = ops[0]
        dt = 0.1
        U = stinespring_unitary_from_lindblad(L, dt)
        return L, dt, U, params

    def test_unitarity(self, setup):
        _L, _dt, U, _params = setup
        eye = np.eye(U.shape[0], dtype=complex)
        assert np.allclose(U.conj().T @ U, eye, atol=1e-10)

    def test_trace_preservation(self, setup):
        _L, _dt, U, _params = setup
        dim = U.shape[0] // 2
        psi = np.zeros(dim, dtype=complex)
        psi[0] = 1.0
        rho = np.outer(psi, psi.conj())
        rho_out = apply_stinespring_to_density_matrix(rho, U)
        assert abs(np.trace(rho_out).real - 1.0) < 1e-10

    def test_hermiticity(self, setup):
        _L, _dt, U, _params = setup
        dim = U.shape[0] // 2
        psi = np.zeros(dim, dtype=complex)
        psi[0] = 1.0
        rho = np.outer(psi, psi.conj())
        rho_out = apply_stinespring_to_density_matrix(rho, U)
        assert np.allclose(rho_out, rho_out.conj().T, atol=1e-12)
