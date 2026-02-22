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


# ---------------------------------------------------------------------------
# 9. TestPhysicalLimits – Unitary limit, fluorescence analytical, steady state
# ---------------------------------------------------------------------------
class TestPhysicalLimits:
    def test_unitary_limit(self):
        """Unitary limit: all gamma=0, entropy must stay zero (pure state)."""
        params = GKSLPhysicalParameters(
            gamma_TTA=0, Gamma_fl=0, Gamma_ph=0,
            k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
        )
        sim = ClassicalGKSLSimulator(params)
        result = sim.simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")

        # ODE solver (RK45, rtol=1e-9) accumulates numerical errors over
        # the integration interval. For t_max=10, the accumulated error in
        # the density matrix eigenvalues leads to entropy ~1e-8.
        # Tolerance 1e-6 is still extremely stringent
        # (maximally mixed 81-dim entropy = ln(81) ≈ 4.4).
        for entropy in result["entropy"]:
            assert abs(entropy) < 1e-6

        # Purity must stay 1 (pure state)
        for purity in result["purity"]:
            assert abs(purity - 1.0) < 1e-6

        # Trace must stay 1
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-12

    def test_fluorescence_analytical(self):
        """Fluorescence-only: V=0, only Gamma_fl, compare with analytical exponential decay.

        The Stinespring+Trotter approach has O(dt^{3/2}) per-step error in
        the dissipator channel, which accumulates over 100 steps.  The
        tolerance reflects this algorithmic accuracy (vs the exact analytical
        solution) while still being a stringent physical check (<0.1%
        relative error).
        """
        Gamma_fl = 0.01
        params = GKSLPhysicalParameters(
            V=0, gamma_TTA=0, Gamma_fl=Gamma_fl, Gamma_ph=0,
            k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
        )
        sim = ClassicalGKSLSimulator(params)
        result = sim.simulate(t_max=50.0, n_steps=100, initial_state="all_singlet")

        for i, t in enumerate(result["times"]):
            # In natural units (hbar=1), the decay rate is simply Gamma_fl.
            # The division by hbar is included for generality but is 1.0.
            expected_N_S1 = 4.0 * np.exp(-Gamma_fl * t / params.hbar)
            actual_N_S1 = result["populations"][i]["N_S1"]
            assert abs(expected_N_S1 - actual_N_S1) < 2e-3, (
                f"t={t}: expected N_S1={expected_N_S1:.6f}, got {actual_N_S1:.6f}"
            )

    def test_steady_state(self):
        """Steady state: long-time evolution relaxes all molecules to ground state S0."""
        params = GKSLPhysicalParameters()
        sim = ClassicalGKSLSimulator(params)
        result = sim.simulate(t_max=1000.0, n_steps=100, initial_state="edge_triplet")

        pops_final = result["populations"][-1]
        assert pops_final["N_S0"] > 3.5, f"N_S0={pops_final['N_S0']:.4f} (expected > 3.5)"
        assert pops_final["N_T1"] < 0.5, f"N_T1={pops_final['N_T1']:.4f} (expected < 0.5)"
        assert pops_final["N_S1"] < 0.5, f"N_S1={pops_final['N_S1']:.4f} (expected < 0.5)"

    def test_validate_unitary_params(self):
        """Validate that all-zero dissipation passes parameter validation."""
        params = GKSLPhysicalParameters(
            gamma_TTA=0, Gamma_fl=0, Gamma_ph=0,
            k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
        )
        errors = params.validate()
        # No weak-coupling violation when all dissipation is zero
        assert not any("Weak coupling" in e for e in errors)

    def test_validate_zero_V_params(self):
        """Validate that V=0 with small dissipation passes (no spurious weak coupling error)."""
        params = GKSLPhysicalParameters(
            V=0, gamma_TTA=0, Gamma_fl=0.01, Gamma_ph=0,
            k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
        )
        errors = params.validate()
        # V=0 is valid; weak coupling should compare against E_T, not V
        assert not any("Weak coupling" in e for e in errors)


# ---------------------------------------------------------------------------
# 10. TestStinespringFidelity – Classical vs Qudit/Qubit fidelity comparison
# ---------------------------------------------------------------------------
class TestStinespringFidelity:
    @staticmethod
    def _quantum_fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
        """Compute quantum state fidelity F(rho, sigma) = (Tr[sqrt(sqrt(rho) sigma sqrt(rho))])^2.

        Uses eigendecomposition for numerical stability (avoids sqrtm on
        near-singular matrices).
        """
        # sqrt(rho) via eigendecomposition
        evals_rho, evecs_rho = np.linalg.eigh(rho)
        evals_rho = np.maximum(evals_rho, 0.0)
        sqrt_rho = evecs_rho @ np.diag(np.sqrt(evals_rho)) @ evecs_rho.conj().T

        M = sqrt_rho @ sigma @ sqrt_rho
        M = (M + M.conj().T) / 2  # enforce Hermiticity
        evals_M = np.linalg.eigvalsh(M)
        evals_M = np.maximum(evals_M, 0.0)
        return float(np.real(np.sum(np.sqrt(evals_M))) ** 2)

    def test_stinespring_fidelity_qudit(self):
        """Qudit Stinespring+Trotter fidelity vs classical ODE: F > 0.99 at small dt."""
        params = GKSLPhysicalParameters()

        sim_classical = ClassicalGKSLSimulator(params)
        result_classical = sim_classical.simulate(
            t_max=1.0, n_steps=20, initial_state="edge_triplet"
        )

        sim_qudit = QuditGKSLSimulator(params)
        result_qudit = sim_qudit.simulate(
            t_max=1.0, n_steps=20, initial_state="edge_triplet"
        )

        rho_classical = result_classical["rho_final"]
        rho_qudit = result_qudit["rho_final"]

        F = self._quantum_fidelity(rho_classical, rho_qudit)
        assert F > 0.99, f"Fidelity = {F:.6f} (expected > 0.99)"


# ---------------------------------------------------------------------------
# 11. TestEdgeTripletBoundary – N-variable edge_triplet and N<2 validation
# ---------------------------------------------------------------------------
class TestEdgeTripletBoundary:
    @staticmethod
    def _make_n1_params() -> GKSLPhysicalParameters:
        """Create a GKSLPhysicalParameters with N_molecules=1 (bypasses __init__ validation)."""
        params = GKSLPhysicalParameters.__new__(GKSLPhysicalParameters)
        params.d = 3
        params.N_molecules = 1
        params.E_T = 1.5
        params.E_S = 3.0
        params.V = 0.1
        params.gamma_TTA = 0.05
        params.Gamma_fl = 0.01
        params.Gamma_ph = 1e-6
        params.k_IC = 0.005
        params.k_ISC_ST = 0.003
        params.k_ISC_TS = 1e-5
        params.with_boson = False
        params.n_max = 2
        params.omega_ph = 0.15
        params.g_eph = 0.02
        return params

    def test_classical_edge_triplet_n2(self):
        """Classical simulator: edge_triplet with N=2 produces correct index."""
        params = GKSLPhysicalParameters(N_molecules=2)
        sim = ClassicalGKSLSimulator(params)
        rho = sim.prepare_initial_state("edge_triplet")
        # |1,1⟩ in base-3: index = 1*3 + 1 = 4
        assert abs(rho[4, 4] - 1.0) < 1e-12
        assert abs(np.trace(rho) - 1.0) < 1e-12

    def test_classical_edge_triplet_n4_consistency(self):
        """Classical simulator: N=4 edge_triplet must match other simulators' indexing."""
        params = GKSLPhysicalParameters(N_molecules=4)
        sim = ClassicalGKSLSimulator(params)
        rho = sim.prepare_initial_state("edge_triplet")
        # |1,0,0,1⟩ in base-3: index = 1*27 + 0*9 + 0*3 + 1 = 28
        assert abs(rho[28, 28] - 1.0) < 1e-12

    def test_classical_edge_triplet_rejects_n1(self):
        """Classical simulator: edge_triplet with N_molecules < 2 raises ValueError."""
        sim = ClassicalGKSLSimulator(self._make_n1_params())
        with pytest.raises(ValueError, match="edge_triplet requires N_molecules >= 2"):
            sim.prepare_initial_state("edge_triplet")

    def test_qubit_edge_triplet_rejects_n1(self):
        """Qubit simulator: edge_triplet with N < 2 raises ValueError."""
        sim = QubitGKSLSimulator(self._make_n1_params())
        with pytest.raises(ValueError, match="edge_triplet requires N_molecules >= 2"):
            sim.prepare_initial_state("edge_triplet")

    def test_qudit_edge_triplet_rejects_n1(self):
        """Qudit simulator: edge_triplet with N < 2 raises ValueError."""
        sim = QuditGKSLSimulator(self._make_n1_params())
        with pytest.raises(ValueError, match="edge_triplet requires N_molecules >= 2"):
            sim.prepare_initial_state("edge_triplet")

    def test_classical_boson_edge_triplet_n2_works(self):
        """Classical boson simulator: edge_triplet with N=2 produces valid state."""
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator

        params = GKSLPhysicalParameters(N_molecules=2, with_boson=True, n_max=1)
        sim = ClassicalGKSLBosonSimulator(params)
        rho = sim.prepare_initial_state("edge_triplet")
        assert abs(np.trace(rho) - 1.0) < 1e-12

    def test_qubit_boson_edge_triplet_n2_works(self):
        """Qubit boson simulator: edge_triplet with N=2 produces valid state."""
        from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator

        params = GKSLPhysicalParameters(N_molecules=2, with_boson=True, n_max=1)
        sim = QubitGKSLBosonSimulator(params)
        rho = sim.prepare_initial_state("edge_triplet")
        assert abs(np.trace(rho) - 1.0) < 1e-12

    def test_qudit_boson_edge_triplet_n2_works(self):
        """Qudit boson simulator: edge_triplet with N=2 produces valid state."""
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator

        params = GKSLPhysicalParameters(N_molecules=2, with_boson=True, n_max=1)
        sim = QuditGKSLBosonSimulator(params)
        rho = sim.prepare_initial_state("edge_triplet")
        assert abs(np.trace(rho) - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# 12. TestBosonGephZeroReduction – g_eph=0 exact reduction
# ---------------------------------------------------------------------------
class TestBosonGephZeroReduction:
    def test_boson_g_eph_zero_matches_non_boson(self):
        """g_eph=0 boson simulator must produce results identical to non-boson."""
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator

        params_nb = GKSLPhysicalParameters(N_molecules=4)
        params_b = GKSLPhysicalParameters(
            N_molecules=4, with_boson=True, n_max=1, g_eph=0.0
        )

        sim_nb = ClassicalGKSLSimulator(params_nb)
        sim_b = ClassicalGKSLBosonSimulator(params_b)

        result_nb = sim_nb.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        result_b = sim_b.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")

        for pop_nb, pop_b in zip(result_nb["populations"], result_b["populations"]):
            assert abs(pop_nb["N_S0"] - pop_b["N_S0"]) < 1e-6
            assert abs(pop_nb["N_T1"] - pop_b["N_T1"]) < 1e-6
            assert abs(pop_nb["N_S1"] - pop_b["N_S1"]) < 1e-6

    def test_boson_g_eph_zero_traces(self):
        """g_eph=0 boson simulator trace preservation matches non-boson exactly."""
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator

        params_b = GKSLPhysicalParameters(
            N_molecules=4, with_boson=True, n_max=2, g_eph=0.0
        )
        sim_b = ClassicalGKSLBosonSimulator(params_b)
        result = sim_b.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-8

    def test_boson_g_eph_zero_method_label(self):
        """g_eph=0 boson simulator reports its reduction in the method field."""
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator

        params_b = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, g_eph=0.0
        )
        sim_b = ClassicalGKSLBosonSimulator(params_b)
        result = sim_b.simulate(t_max=1.0, n_steps=2, initial_state="edge_triplet")
        assert "g_eph=0" in result["method"]


# ---------------------------------------------------------------------------
# 13. TestQuditGKSLNoisySimulator
# ---------------------------------------------------------------------------
class TestQuditGKSLNoisySimulator:
    @pytest.fixture()
    def params(self):
        return GKSLPhysicalParameters()

    def test_zero_noise_matches_ideal(self, params):
        """p_depol=0 must produce results identical to the ideal QuditGKSLSimulator."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim_ideal = QuditGKSLSimulator(params)
        sim_zero = QuditGKSLNoisySimulator(params, p_depol=0.0)
        r_ideal = sim_ideal.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        r_zero = sim_zero.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for pi, pz in zip(r_ideal["populations"], r_zero["populations"]):
            assert abs(pi["N_S0"] - pz["N_S0"]) < 1e-12
            assert abs(pi["N_T1"] - pz["N_T1"]) < 1e-12
            assert abs(pi["N_S1"] - pz["N_S1"]) < 1e-12

    def test_noise_reduces_purity(self, params):
        """Hardware noise must decrease purity (increase mixedness)."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim_ideal = QuditGKSLSimulator(params)
        sim_noisy = QuditGKSLNoisySimulator(params, p_depol=0.01)
        r_ideal = sim_ideal.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        r_noisy = sim_noisy.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        assert r_noisy["purity"][-1] < r_ideal["purity"][-1]

    def test_trace_preservation(self, params):
        """Noisy simulation must preserve trace within tolerance."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim = QuditGKSLNoisySimulator(params, p_depol=0.01)
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-10

    def test_method_label(self, params):
        """Method field must indicate noisy qudit simulation."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim = QuditGKSLNoisySimulator(params, p_depol=0.01)
        result = sim.simulate(t_max=2.0, n_steps=2, initial_state="edge_triplet")
        assert "noisy" in result["method"]

    def test_noise_params_in_result(self, params):
        """Result dict must contain noise_params."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim = QuditGKSLNoisySimulator(params, p_depol=0.02, p_dephasing=0.005)
        result = sim.simulate(t_max=2.0, n_steps=2, initial_state="edge_triplet")
        assert result["noise_params"]["p_depol"] == 0.02
        assert result["noise_params"]["p_dephasing"] == 0.005

    def test_invalid_p_depol(self, params):
        """Invalid depolarization probability must raise ValueError."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        with pytest.raises(ValueError, match="p_depol"):
            QuditGKSLNoisySimulator(params, p_depol=-0.1)
        with pytest.raises(ValueError, match="p_depol"):
            QuditGKSLNoisySimulator(params, p_depol=1.5)

    def test_dephasing_reduces_coherence(self, params):
        """Dephasing noise must reduce off-diagonal elements of density matrix."""
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim_no_deph = QuditGKSLNoisySimulator(params, p_depol=0.0, p_dephasing=0.0)
        sim_deph = QuditGKSLNoisySimulator(params, p_depol=0.0, p_dephasing=0.05)
        r_nd = sim_no_deph.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        r_d = sim_deph.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        # Dephasing reduces off-diagonal norm
        offdiag_nd = np.linalg.norm(r_nd["rho_final"] - np.diag(np.diag(r_nd["rho_final"])))
        offdiag_d = np.linalg.norm(r_d["rho_final"] - np.diag(np.diag(r_d["rho_final"])))
        assert offdiag_d < offdiag_nd


# ---------------------------------------------------------------------------
# 14. TestQubitGKSLNoisySimulator
# ---------------------------------------------------------------------------
class TestQubitGKSLNoisySimulator:
    @pytest.fixture()
    def params(self):
        return GKSLPhysicalParameters()

    def test_zero_noise_matches_ideal(self, params):
        """p_depol=0 must produce results identical to the ideal QubitGKSLSimulator."""
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

        sim_ideal = QubitGKSLSimulator(params)
        sim_zero = QubitGKSLNoisySimulator(params, p_depol=0.0)
        r_ideal = sim_ideal.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        r_zero = sim_zero.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for pi, pz in zip(r_ideal["populations"], r_zero["populations"]):
            assert abs(pi["N_S0"] - pz["N_S0"]) < 1e-12
            assert abs(pi["N_T1"] - pz["N_T1"]) < 1e-12
            assert abs(pi["N_S1"] - pz["N_S1"]) < 1e-12

    def test_noise_reduces_purity(self, params):
        """Hardware noise must decrease purity."""
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

        sim_ideal = QubitGKSLSimulator(params)
        sim_noisy = QubitGKSLNoisySimulator(params, p_depol=0.01)
        r_ideal = sim_ideal.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        r_noisy = sim_noisy.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        assert r_noisy["purity"][-1] < r_ideal["purity"][-1]

    def test_trace_preservation(self, params):
        """Noisy qubit trace + forbidden-state population must sum to 1."""
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

        sim = QubitGKSLNoisySimulator(params, p_depol=0.01)
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr, fp in zip(result["trace"], result["forbidden_state_population"]):
            # trace (qutrit subspace) + forbidden pop = 1 (total probability)
            assert abs(tr + fp - 1.0) < 1e-8
            # trace must not exceed 1 (no probability gain)
            assert tr <= 1.0 + 1e-10

    def test_method_label(self, params):
        """Method field must indicate noisy qubit simulation."""
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

        sim = QubitGKSLNoisySimulator(params, p_depol=0.01)
        result = sim.simulate(t_max=2.0, n_steps=2, initial_state="edge_triplet")
        assert "noisy" in result["method"]

    def test_noise_params_in_result(self, params):
        """Result dict must contain noise_params with T1, T2, p_reset."""
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

        sim = QubitGKSLNoisySimulator(params, p_depol=0.01, T1=5e10, T2=7e10)
        result = sim.simulate(t_max=2.0, n_steps=2, initial_state="edge_triplet")
        assert result["noise_params"]["p_depol"] == 0.01
        assert result["noise_params"]["T1"] == 5e10
        assert result["noise_params"]["p_reset"] > 0

    def test_thermal_relaxation_trace(self, params):
        """Thermal relaxation: trace + forbidden pop must sum to 1."""
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator

        sim = QubitGKSLNoisySimulator(params, p_depol=0.01, T1=5e10, T2=7e10, t_gate=300.0)
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for tr, fp in zip(result["trace"], result["forbidden_state_population"]):
            assert abs(tr + fp - 1.0) < 1e-8
            assert tr <= 1.0 + 1e-10

    def test_qubit_qudit_noisy_consistency(self, params):
        """Noisy qubit and qudit simulators should produce qualitatively similar results.

        Qubit uses d=4 Pauli noise (with leakage), qudit uses d=3 Weyl-Heisenberg
        noise (no leakage). Both use the same p_depol, so purity reduction is similar
        in magnitude but differs due to the noise model.
        """
        from qubit_gksl_noisy_simulator import QubitGKSLNoisySimulator
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        sim_qb = QubitGKSLNoisySimulator(params, p_depol=0.01)
        sim_qd = QuditGKSLNoisySimulator(params, p_depol=0.01)
        r_qb = sim_qb.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        r_qd = sim_qd.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")
        for p_qb, p_qd in zip(r_qb["populations"], r_qd["populations"]):
            # Wider tolerance due to different noise models and leakage
            assert abs(p_qb["N_T1"] - p_qd["N_T1"]) < 1.0
            assert abs(p_qb["N_S0"] - p_qd["N_S0"]) < 1.0


# ---------------------------------------------------------------------------
# 15. TestNoiseChannelProperties – mathematical rigor of noise channels
# ---------------------------------------------------------------------------
class TestNoiseChannelProperties:
    """Verify that each noise channel is a valid quantum channel (CPTP)."""

    def test_depolarization_single_trace(self):
        """Single-site depolarization preserves trace."""
        from qudit_gksl_noisy_simulator import _apply_local_depolarization_single

        d, N = 3, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_noisy = _apply_local_depolarization_single(rho, 1, d, N, 0.1)
        assert abs(np.trace(rho_noisy).real - 1.0) < 1e-12

    def test_depolarization_pair_trace(self):
        """Pair depolarization preserves trace."""
        from qudit_gksl_noisy_simulator import _apply_local_depolarization_pair

        d, N = 3, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_noisy = _apply_local_depolarization_pair(rho, 0, 2, d, N, 0.15)
        assert abs(np.trace(rho_noisy).real - 1.0) < 1e-12

    def test_depolarization_hermiticity(self):
        """Depolarization preserves Hermiticity."""
        from qudit_gksl_noisy_simulator import _apply_local_depolarization_single

        d, N = 3, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_noisy = _apply_local_depolarization_single(rho, 2, d, N, 0.05)
        assert np.allclose(rho_noisy, rho_noisy.conj().T, atol=1e-12)

    def test_depolarization_positivity(self):
        """Depolarization preserves positive semi-definiteness."""
        from qudit_gksl_noisy_simulator import _apply_local_depolarization_single

        d, N = 3, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_noisy = _apply_local_depolarization_single(rho, 0, d, N, 0.2)
        eigs = np.linalg.eigvalsh(rho_noisy)
        assert eigs.min() >= -1e-12

    def test_full_depolarization_gives_maximally_mixed(self):
        """p=1 depolarization on all sites gives maximally mixed state."""
        from qudit_gksl_noisy_simulator import _apply_local_depolarization_single

        d, N = 3, 4
        dim = d**N
        psi = np.zeros(dim, dtype=np.complex128)
        psi[0] = 1.0
        rho = np.outer(psi, psi.conj())
        for s in range(N):
            rho = _apply_local_depolarization_single(rho, s, d, N, 1.0)
        assert np.allclose(rho, np.eye(dim) / dim, atol=1e-10)

    def test_dephasing_trace(self):
        """Dephasing preserves trace."""
        from qudit_gksl_noisy_simulator import _apply_local_dephasing_single

        d, N = 3, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_deph = _apply_local_dephasing_single(rho, 0, d, N, 0.3)
        assert abs(np.trace(rho_deph).real - 1.0) < 1e-12

    def test_thermal_relaxation_trace(self):
        """Thermal relaxation preserves trace (d=4 qubit version)."""
        from qubit_gksl_noisy_simulator import _apply_thermal_relaxation_qubit

        d, N = 4, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_th = _apply_thermal_relaxation_qubit(rho, 1, N, 0.1)
        assert abs(np.trace(rho_th).real - 1.0) < 1e-12

    def test_thermal_relaxation_positivity(self):
        """Thermal relaxation preserves positive semi-definiteness (d=4 qubit version)."""
        from qubit_gksl_noisy_simulator import _apply_thermal_relaxation_qubit

        d, N = 4, 4
        dim = d**N
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        rho_th = _apply_thermal_relaxation_qubit(rho, 0, N, 0.5)
        eigs = np.linalg.eigvalsh(rho_th)
        assert eigs.min() >= -1e-12


# ================================================================
# Circuit simulator tests
# ================================================================


class TestQuditGKSLCircuitSimulator:
    """Tests for the MQT-Qudits circuit-based GKSL simulator."""

    def test_initialization(self):
        """Circuit simulator initializes with correct local operators."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        assert sim.N == 4
        assert sim.d == 3
        assert sim.dim == 81
        assert len(sim.lindblad_local_info) == 26

    def test_boson_rejected(self):
        """Circuit simulator rejects boson parameters."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters(with_boson=True)
        with pytest.raises(ValueError, match="non-boson"):
            QuditGKSLCircuitSimulator(params)

    def test_hamiltonian_circuit_unitarity(self):
        """Circuit-decomposed Hamiltonian is unitary."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        v = sim.verify_hamiltonian_circuit(0.5)
        assert v["circuit_is_unitary"]

    def test_hamiltonian_circuit_trotter_convergence(self):
        """Circuit Hamiltonian converges to exact as dt -> 0."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        d1 = sim.verify_hamiltonian_circuit(1.0)["frobenius_distance"]
        d2 = sim.verify_hamiltonian_circuit(0.1)["frobenius_distance"]
        # Distance should decrease with dt (Trotter error ~ dt^2)
        assert d2 < d1

    def test_stinespring_local_matches_full(self):
        """Local Stinespring unitaries match full-system computation exactly."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        v = sim.verify_stinespring_circuit(1.0)
        assert v["all_match"]
        assert v["max_frobenius_distance"] < 1e-10

    def test_mqt_circuit_execution(self):
        """MQT-Qudits circuit execution matches matrix computation."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        v = sim.verify_circuit_via_mqt(0.5)
        assert v["match"]
        assert v["statevector_distance"] < 1e-10
        expected_gates = params.N_molecules + len(params.neighbors)
        assert v["gate_count"] == expected_gates

    def test_trace_preservation(self):
        """Circuit-based simulation preserves trace."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5)
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-10

    def test_circuit_matches_matrix_simulator(self):
        """Circuit-based simulator matches matrix-level simulator closely."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_matrix = QuditGKSLSimulator(params)
        sim_circuit = QuditGKSLCircuitSimulator(params)

        r_m = sim_matrix.simulate(t_max=5.0, n_steps=5)
        r_c = sim_circuit.simulate(t_max=5.0, n_steps=5)

        # Stinespring channels match exactly; Hamiltonian has Trotter error
        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(r_m["populations"][-1][key] - r_c["populations"][-1][key])
            assert diff < 1e-3, f"{key} mismatch: {diff}"

    def test_gate_breakdown(self):
        """Gate breakdown is reported correctly."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5)
        gb = result["gate_breakdown"]
        N = params.N_molecules
        n_pairs = len(params.neighbors)
        n_single = 5 * N  # 5 single-site operator types × N molecules
        n_pair = 2 * n_pairs  # 2 TTA channels × number of pairs
        assert gb["cu_one_onsite"] == 2 * N
        assert gb["cu_two_transfer"] == 2 * n_pairs
        assert gb["cu_two_stinespring_single"] == n_single
        assert gb["cu_multi_stinespring_pair"] == n_pair
        assert result["gates_per_step"] == 2 * (N + n_pairs) + n_single + n_pair

    def test_method_label(self):
        """Result contains correct method label."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5)
        assert result["method"] == "qudit_gksl_circuit"

    def test_build_full_trotter_step_circuit(self):
        """Full Trotter step circuit construction reports correct counts."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        info = sim.build_full_trotter_step_circuit(dt=1.0)
        assert info["total_gates"] == 40
        assert info["n_stinespring_gates"] == 26


# ================================================================
# Convergence tests (Classical vs Qudit/Qubit with fine dt)
# ================================================================


class TestClassicalQuantumConvergence:
    """Tests verifying Classical-Qudit and Classical-Qubit convergence.

    The Stinespring+Trotter quantum simulators converge to the Classical ODE
    result as dt -> 0. With n_steps=20 (dt=0.25 for t_max=5), all three
    methods agree to within 1e-3 on population dynamics.
    """

    def test_classical_qudit_convergence(self):
        """Classical and Qudit agree within 1e-3 for fine dt."""
        from classical_gksl_simulator import ClassicalGKSLSimulator
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_cl = ClassicalGKSLSimulator(params)
        sim_qd = QuditGKSLSimulator(params)

        t_max = 5.0
        n_steps = 20  # dt = 0.25
        r_cl = sim_cl.simulate(t_max=t_max, n_steps=n_steps)
        r_qd = sim_qd.simulate(t_max=t_max, n_steps=n_steps)

        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(
                r_cl["populations"][-1][key] - r_qd["populations"][-1][key]
            )
            assert diff < 1e-3, f"{key} diff {diff} exceeds 1e-3"

    def test_classical_qubit_convergence(self):
        """Classical and Qubit agree within 1e-3 for fine dt."""
        from classical_gksl_simulator import ClassicalGKSLSimulator
        from qubit_gksl_simulator import QubitGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_cl = ClassicalGKSLSimulator(params)
        sim_qb = QubitGKSLSimulator(params)

        t_max = 5.0
        n_steps = 20  # dt = 0.25
        r_cl = sim_cl.simulate(t_max=t_max, n_steps=n_steps)
        r_qb = sim_qb.simulate(t_max=t_max, n_steps=n_steps)

        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(
                r_cl["populations"][-1][key] - r_qb["populations"][-1][key]
            )
            assert diff < 1e-3, f"{key} diff {diff} exceeds 1e-3"

    def test_convergence_improves_with_dt(self):
        """Convergence improves as dt decreases (Trotter error scaling)."""
        from classical_gksl_simulator import ClassicalGKSLSimulator
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_cl = ClassicalGKSLSimulator(params)
        sim_qd = QuditGKSLSimulator(params)

        t_max = 5.0
        diffs = []
        for n_steps in [5, 20]:
            r_cl = sim_cl.simulate(t_max=t_max, n_steps=n_steps)
            r_qd = sim_qd.simulate(t_max=t_max, n_steps=n_steps)
            max_diff = max(
                abs(r_cl["populations"][-1][k] - r_qd["populations"][-1][k])
                for k in ["N_S0", "N_T1", "N_S1"]
            )
            diffs.append(max_diff)

        # With more steps (smaller dt), difference should decrease
        assert diffs[1] < diffs[0]


# ================================================================
# Native gate compilation tests
# ================================================================


class TestNativeGateCompilation:
    """Tests for compileO0/compileO1 native gate decomposition."""

    def test_compile_to_native_gates_o0(self):
        """compileO0 produces native gate statistics."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.compile_to_native_gates(dt=0.5, optimization_level=0)

        assert result["hamiltonian"]["native_gates_total"] > 0
        assert result["stinespring_single"]["native_gates_total"] > 0
        assert result["stinespring_pair"]["uncompiled_cu_multi"] == 6
        assert result["per_step_summary"]["native_gates"] > 0
        assert result["per_step_summary"]["optimization_level"] == 0

    def test_compile_to_native_gates_o1(self):
        """compileO1 also produces native gate statistics."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.compile_to_native_gates(dt=0.5, optimization_level=1)

        assert result["hamiltonian"]["native_gates_total"] > 0
        assert result["per_step_summary"]["optimization_level"] == 1

    def test_cu_one_compiles_to_virtrz(self):
        """Diagonal on-site Hamiltonian compiles to VirtRz gates only."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.compile_to_native_gates(dt=0.5, optimization_level=0)

        # cu_one is diagonal -> should compile to VirtRz only
        breakdown = result["hamiltonian"]["cu_one_breakdown"]
        assert "VirtRz" in breakdown
        assert sum(v for k, v in breakdown.items() if k != "VirtRz") == 0

    def test_cu_two_compiles_to_native_set(self):
        """Transfer Hamiltonian compiles to native gate set."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.compile_to_native_gates(dt=0.5, optimization_level=0)

        # cu_two -> should have native gates from {R, Rz, Rh, VirtRz, CEx}
        native_types = {"R", "Rz", "Rh", "VirtRz", "CEx"}
        for pair_info in result["hamiltonian"]["cu_two_per_pair"]:
            for gate_name in pair_info["breakdown"]:
                assert gate_name in native_types, f"Unexpected gate: {gate_name}"

    def test_cu_multi_not_decomposed(self):
        """cu_multi (TTA pair Stinespring) gates are correctly reported as uncompiled."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.compile_to_native_gates(dt=0.5, optimization_level=0)

        n_pairs = len(params.neighbors)
        n_tta_channels = 2  # TTA channel 1 and 2
        expected_cu_multi = n_tta_channels * n_pairs
        assert result["stinespring_pair"]["uncompiled_cu_multi"] == expected_cu_multi
        for ch in result["stinespring_pair"]["per_channel"]:
            assert ch["status"] == "uncompiled_cu_multi"

    def test_verify_compiled_circuit(self):
        """Compiled and uncompiled circuits produce identical state vectors."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        result = sim.verify_compiled_circuit(dt=0.5, optimization_level=0)

        assert result["match"]
        assert result["statevector_distance"] < 1e-8
        assert result["n_compiled_gates"] >= result["n_original_gates"]

    def test_invalid_optimization_level(self):
        """Invalid optimization level raises ValueError."""
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLCircuitSimulator(params)
        with pytest.raises(ValueError, match="optimization_level"):
            sim.compile_to_native_gates(dt=0.5, optimization_level=2)


# ================================================================
# Boson circuit simulator tests
# ================================================================


class TestQuditGKSLCircuitBosonSimulator:
    """Tests for QuditGKSLCircuitBosonSimulator."""

    def test_requires_boson_params(self):
        """Non-boson params raise ValueError."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

        params = GKSLPhysicalParameters(with_boson=False)
        with pytest.raises(ValueError, match="with_boson"):
            QuditGKSLCircuitBosonSimulator(params)

    def test_g_eph_zero_exact_reduction(self):
        """g_eph=0 delegates to non-boson circuit simulator."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

        params = GKSLPhysicalParameters(N_molecules=2, with_boson=True, g_eph=0.0)
        sim = QuditGKSLCircuitBosonSimulator(params)
        result = sim.simulate(t_max=2.0, n_steps=5)

        assert "g_eph=0" in result["method"]
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-10

    def test_g_eph_zero_matches_non_boson(self):
        """g_eph=0 boson circuit matches non-boson circuit simulator."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator

        params_b = GKSLPhysicalParameters(N_molecules=2, with_boson=True, g_eph=0.0)
        params_nb = GKSLPhysicalParameters(N_molecules=2, with_boson=False)

        r_b = QuditGKSLCircuitBosonSimulator(params_b).simulate(t_max=2.0, n_steps=5)
        r_nb = QuditGKSLCircuitSimulator(params_nb).simulate(t_max=2.0, n_steps=5)

        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(r_b["populations"][-1][key] - r_nb["populations"][-1][key])
            assert diff < 1e-10, f"{key} mismatch: {diff}"

    def test_trace_preservation(self):
        """Boson circuit simulation preserves trace."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, g_eph=0.005
        )
        sim = QuditGKSLCircuitBosonSimulator(params)
        result = sim.simulate(t_max=2.0, n_steps=5)

        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-8, f"Trace violation: {tr}"

    def test_matches_matrix_boson_simulator(self):
        """Circuit boson simulator matches matrix-level boson simulator."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator

        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, g_eph=0.005
        )
        r_circuit = QuditGKSLCircuitBosonSimulator(params).simulate(
            t_max=2.0, n_steps=5
        )
        r_matrix = QuditGKSLBosonSimulator(params).simulate(
            t_max=2.0, n_steps=5
        )

        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(
                r_circuit["populations"][-1][key]
                - r_matrix["populations"][-1][key]
            )
            assert diff < 1e-10, f"{key} mismatch: {diff}"

    def test_gate_breakdown(self):
        """Gate breakdown is reported correctly for boson model."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, g_eph=0.005
        )
        sim = QuditGKSLCircuitBosonSimulator(params)
        result = sim.simulate(t_max=2.0, n_steps=5)

        gb = result["gate_breakdown"]
        N = params.N_molecules
        n_pairs = len(params.neighbors)
        assert gb["cu_one_el_onsite"] == 2 * N
        assert gb["cu_two_el_transfer"] == 2 * n_pairs
        assert gb["cu_one_ph_onsite"] == 2 * N
        assert gb["cu_two_eph_coupling"] == 2 * N  # g_eph > 0
        assert result["n_el_qutrits"] == N
        assert result["n_ph_qutrits"] == N

    def test_method_label(self):
        """Result contains correct method label."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, g_eph=0.005
        )
        sim = QuditGKSLCircuitBosonSimulator(params)
        result = sim.simulate(t_max=2.0, n_steps=5)
        assert result["method"] == "qudit_gksl_circuit_boson"

    def test_edge_triplet_requires_n2(self):
        """edge_triplet initial state requires N >= 2."""
        from qudit_gksl_circuit_boson_simulator import QuditGKSLCircuitBosonSimulator

        params = GKSLPhysicalParameters(
            N_molecules=1, with_boson=True, n_max=1, g_eph=0.0
        )
        sim = QuditGKSLCircuitBosonSimulator(params)
        with pytest.raises(ValueError, match="edge_triplet"):
            sim.prepare_initial_state("edge_triplet")


# ================================================================
# Qubit GKSL Circuit Simulator (Qiskit) tests
# ================================================================


class TestQubitGKSLCircuitSimulator:
    """Tests for the Qiskit circuit-based Qubit GKSL simulator."""

    def test_initialization(self):
        """Circuit simulator initializes with correct operators."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        assert sim.N == 4
        assert sim.d == 3
        assert sim.dim == 81
        assert len(sim.lindblad_local_info) == 26

    def test_boson_rejected(self):
        """Circuit simulator rejects boson parameters."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters(with_boson=True)
        with pytest.raises(ValueError, match="non-boson"):
            QubitGKSLCircuitSimulator(params)

    def test_embedding_unitarity(self):
        """All embedded qubit unitaries are unitary."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        v = sim.verify_embedding_unitarity(0.5)
        assert v["all_unitary"]
        assert v["onsite_unitarity_residual"] < 1e-10
        assert v["transfer_max_unitarity_residual"] < 1e-10
        assert v["stinespring_single_max_unitarity_residual"] < 1e-10
        assert v["stinespring_pair_max_unitarity_residual"] < 1e-10

    def test_qiskit_circuit_matches_operator(self):
        """Qiskit circuit operator matches direct computation."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        v = sim.verify_circuit_via_qiskit(0.5)
        assert v["match"]
        assert v["operator_distance"] < 1e-8

    def test_trace_preservation(self):
        """Circuit-based simulation preserves trace."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5)
        for tr in result["trace"]:
            assert abs(tr - 1.0) < 1e-10

    def test_circuit_matches_matrix_simulator(self):
        """Circuit simulator matches matrix-level QubitGKSLSimulator."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator
        from qubit_gksl_simulator import QubitGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_matrix = QubitGKSLSimulator(params)
        sim_circuit = QubitGKSLCircuitSimulator(params)

        r_m = sim_matrix.simulate(t_max=5.0, n_steps=5)
        r_c = sim_circuit.simulate(t_max=5.0, n_steps=5)

        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(r_m["populations"][-1][key] - r_c["populations"][-1][key])
            assert diff < 1e-3, f"{key} mismatch: {diff}"

    def test_gate_breakdown(self):
        """Gate breakdown is reported correctly."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5)
        gb = result["gate_breakdown"]
        N = params.N_molecules
        n_pairs = len(params.neighbors)
        n_single = 5 * N
        n_pair_ops = 2 * n_pairs
        assert gb["unitary_4x4_onsite"] == 2 * N
        assert gb["unitary_16x16_transfer"] == 2 * n_pairs
        assert gb["unitary_8x8_stinespring_single"] == n_single
        assert gb["unitary_32x32_stinespring_pair"] == n_pair_ops
        assert result["gates_per_step"] == 2 * (N + n_pairs) + n_single + n_pair_ops

    def test_method_label(self):
        """Result contains correct method label."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5)
        assert result["method"] == "qubit_gksl_circuit"

    def test_build_full_trotter_step_circuit(self):
        """Full Trotter step circuit construction reports correct counts."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        info = sim.build_full_trotter_step_circuit(dt=1.0)
        assert info["total_gates"] == 40
        assert info["n_stinespring_gates"] == 26

    def test_qubit_circuit_info(self):
        """Simulation result includes correct qubit resource info."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=1.0, n_steps=2)
        assert result["n_system_qubits"] == 2 * params.N_molecules
        assert result["n_ancilla_qubits"] == 26
        assert result["n_total_qubits"] == 2 * params.N_molecules + 26
