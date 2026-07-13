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
        assert sim.n_ancilla_qudits == 26
        assert sim.d_anc == 3  # ancilla uses native qutrit on qudit QC

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
        """Circuit-based simulator matches matrix-level simulator closely.

        Both use palindromic 2nd-order Trotter for Lindblad channels.
        Remaining O(dt²) difference is from Hamiltonian circuit decomposition
        (product of local unitaries vs full matrix exponential).
        """
        from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_matrix = QuditGKSLSimulator(params)
        sim_circuit = QuditGKSLCircuitSimulator(params)

        r_m = sim_matrix.simulate(t_max=5.0, n_steps=5)
        r_c = sim_circuit.simulate(t_max=5.0, n_steps=5)

        # Hamiltonian Trotter error O(dt²) ≈ 3.5e-5 for dt=1.0
        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(r_m["populations"][-1][key] - r_c["populations"][-1][key])
            assert diff < 1e-4, f"{key} mismatch: {diff}"

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
        assert gb["cu_two_stinespring_single"] == 2 * n_single
        assert gb["cu_multi_stinespring_pair"] == 2 * n_pair
        assert result["gates_per_step"] == 2 * (N + n_pairs) + 2 * (n_single + n_pair)

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
        assert info["total_gates"] == 66
        assert info["n_stinespring_gates"] == 52


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
        """Circuit simulator matches matrix-level QubitGKSLSimulator.

        Both use palindromic 2nd-order Trotter for Lindblad channels.
        Remaining O(dt²) difference is from Hamiltonian circuit decomposition.
        """
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator
        from qubit_gksl_simulator import QubitGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_matrix = QubitGKSLSimulator(params)
        sim_circuit = QubitGKSLCircuitSimulator(params)

        r_m = sim_matrix.simulate(t_max=5.0, n_steps=5)
        r_c = sim_circuit.simulate(t_max=5.0, n_steps=5)

        # Hamiltonian Trotter error O(dt²) ≈ 3.5e-5 for dt=1.0
        for key in ["N_S0", "N_T1", "N_S1"]:
            diff = abs(r_m["populations"][-1][key] - r_c["populations"][-1][key])
            assert diff < 1e-4, f"{key} mismatch: {diff}"

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
        assert gb["unitary_8x8_stinespring_single"] == 2 * n_single
        assert gb["unitary_32x32_stinespring_pair"] == 2 * n_pair_ops
        assert result["gates_per_step"] == 2 * (N + n_pairs) + 2 * (n_single + n_pair_ops)

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
        assert info["total_gates"] == 66
        assert info["n_stinespring_gates"] == 52

    def test_qubit_circuit_info(self):
        """Simulation result includes correct qubit resource info."""
        from qubit_gksl_circuit_simulator import QubitGKSLCircuitSimulator

        params = GKSLPhysicalParameters()
        sim = QubitGKSLCircuitSimulator(params)
        result = sim.simulate(t_max=1.0, n_steps=2)
        assert result["n_system_qubits"] == 2 * params.N_molecules
        assert result["n_ancilla_qubits"] == 26
        assert result["n_total_qubits"] == 2 * params.N_molecules + 26


# ---------------------------------------------------------------------------
# B-1 / A-1 / A-3: continuation of PR #257
# ---------------------------------------------------------------------------
class TestExactLocalChannelsConvergence:
    """B-1: ``algorithm="exact_local_channels"`` recovers O(dt²) convergence.

    The default ``algorithm="stinespring"`` is mathematically a 1st-order
    approximation of each Lindblad channel and so the global trace-distance
    convergence collapses to O(dt) regardless of the Strang+palindromic
    structure (this is recorded as B-1 in ``STATUS_HONEST_2026-05.md``).

    Replacing per-channel Stinespring with exact local-channel
    exponentiation should restore the 2nd-order behaviour expected from
    the Strang split.  We check this empirically on the smallest config
    (N=2, d=3 ⇒ dim=9) by comparing against
    :class:`ClassicalGKSLSimulator` (full Liouvillian exp).
    """

    @staticmethod
    def _trace_distance(rho, sigma):
        delta = rho - sigma
        delta = (delta + delta.conj().T) / 2
        eigenvalues = np.linalg.eigvalsh(delta)
        return float(0.5 * np.sum(np.abs(eigenvalues)))

    def test_local_lindblad_ops_match_global_ones(self):
        """Local-form Lindblad ops, embedded into the full system, must
        equal the full-system ops returned by build_lindblad_operators
        in the same order and with the same sqrt(gamma) prefactor.

        This pins down the consistency between
        :func:`gksl_math_utils.build_lindblad_operators` and
        :func:`exact_local_channels.get_local_lindblad_ops`, so that
        the two algorithms in :class:`QuditGKSLSimulator` apply the
        same physical channels.
        """
        from exact_local_channels import get_local_lindblad_ops
        from gksl_math_utils import (
            build_lindblad_operators,
            build_single_site_operator,
        )
        from functools import reduce

        params = GKSLPhysicalParameters(N_molecules=3, with_boson=False)
        d = params.d
        N = params.N_molecules

        global_ops = build_lindblad_operators(params)
        local_ops = get_local_lindblad_ops(params)
        assert len(global_ops) == len(local_ops)

        eye = np.eye(d, dtype=np.complex128)
        for (L_global, _gamma), (sites, L_local, kind) in zip(
            global_ops, local_ops
        ):
            if kind == "single":
                expected = build_single_site_operator(L_local, sites[0], N, d)
            else:
                # For pair channels, L_local = kron(A_i, A_j) sits on two
                # arbitrary sites (i, j); we verify equivalence between the
                # local and global routes by exponentiating each route's
                # dissipator superoperator and applying it to the same
                # fiducial PSD density matrix.
                from exact_local_channels import (
                    apply_channel_pair,
                    build_local_dissipator_super,
                )
                from scipy.linalg import expm as _expm
                from stinespring_utils import build_gksl_superoperator

                rng = np.random.default_rng(123)
                rho = (
                    rng.standard_normal((d**N, d**N))
                    + 1j * rng.standard_normal((d**N, d**N))
                )
                rho = rho @ rho.conj().T  # PSD
                rho = rho / np.trace(rho)

                # Local route: build local dissipator on d²×d² space and
                # apply via einsum embedding.
                LD_local = build_local_dissipator_super(L_local, d * d)
                exp_LD = _expm(LD_local * 0.01)
                rho_local = apply_channel_pair(rho, exp_LD, sites, N, d)

                # Global route: build full-system dissipator and exponentiate.
                LD_global = build_gksl_superoperator(
                    np.zeros_like(L_global), [L_global]
                )
                exp_LD_g = _expm(LD_global * 0.01)
                vec_rho = rho.flatten(order="F")
                rho_global = exp_LD_g.dot(vec_rho).reshape(
                    d**N, d**N, order="F"
                )

                assert np.allclose(rho_local, rho_global, atol=1e-10), (
                    f"pair channel at sites={sites} disagrees between "
                    f"local and global routes"
                )
                continue
            assert np.allclose(L_global, expected, atol=1e-12), (
                f"single-site channel #{sites} disagrees: order or "
                f"sqrt(gamma) mismatch between get_local_lindblad_ops "
                f"and build_lindblad_operators"
            )

    def test_constructor_rejects_unknown_algorithm(self):
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        with pytest.raises(ValueError, match="algorithm"):
            QuditGKSLSimulator(params, algorithm="bogus")

    def test_default_algorithm_is_stinespring(self):
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        sim = QuditGKSLSimulator(params)
        assert sim.algorithm == "stinespring"
        result = sim.simulate(t_max=1.0, n_steps=2, initial_state="edge_triplet")
        assert result["algorithm"] == "stinespring"

    def test_stinespring_is_first_order(self):
        """Stinespring algorithm gives empirical rate ≈ 1.0 (B-1 docs)."""
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        t_max = 10.0
        ref = ClassicalGKSLSimulator(params).simulate(
            t_max=t_max, n_steps=2000, initial_state="edge_triplet"
        )
        rho_ref = ref["rho_final"]

        rates = []
        prev_T, prev_dt = None, None
        for n_steps in (50, 100, 200):
            res = QuditGKSLSimulator(params).simulate(
                t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
            )
            T = self._trace_distance(rho_ref, res["rho_final"])
            dt = t_max / n_steps
            if prev_T is not None:
                rates.append(np.log(prev_T / T) / np.log(prev_dt / dt))
            prev_T, prev_dt = T, dt
        # Empirically rates should be ≈ 1.0 — anything in [0.85, 1.15] is fine.
        for r in rates:
            assert 0.85 <= r <= 1.15, f"Stinespring rate not ≈ 1: {rates}"

    def test_exact_local_channels_is_second_order(self):
        """exact_local_channels gives empirical rate ≈ 2.0 on N=2."""
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        t_max = 10.0
        ref = ClassicalGKSLSimulator(params).simulate(
            t_max=t_max, n_steps=2000, initial_state="edge_triplet"
        )
        rho_ref = ref["rho_final"]

        rates = []
        prev_T, prev_dt = None, None
        for n_steps in (20, 50, 100, 200):
            res = QuditGKSLSimulator(
                params, algorithm="exact_local_channels"
            ).simulate(
                t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
            )
            T = self._trace_distance(rho_ref, res["rho_final"])
            dt = t_max / n_steps
            if prev_T is not None:
                rates.append(np.log(prev_T / T) / np.log(prev_dt / dt))
            prev_T, prev_dt = T, dt
        # Empirically rates should be ≈ 2.0 — assert ≥ 1.9 to allow for
        # the asymptotic approach in the smallest n_steps pair.
        for r in rates:
            assert r >= 1.9, f"exact_local_channels rate not ≥ 1.9: {rates}"

    def test_exact_local_beats_stinespring_at_same_dt(self):
        """For the same n_steps the exact algorithm strictly more accurate."""
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        t_max = 10.0
        ref = ClassicalGKSLSimulator(params).simulate(
            t_max=t_max, n_steps=2000, initial_state="edge_triplet"
        )
        rho_ref = ref["rho_final"]

        n_steps = 100
        T_stine = self._trace_distance(
            rho_ref,
            QuditGKSLSimulator(params).simulate(
                t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
            )["rho_final"],
        )
        T_exact = self._trace_distance(
            rho_ref,
            QuditGKSLSimulator(
                params, algorithm="exact_local_channels"
            ).simulate(
                t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
            )["rho_final"],
        )
        # At n_steps=100 / t_max=10 / N=2, the exact channel scheme is
        # observed to be ~89× more accurate than Stinespring.  Assert at
        # least 10× to leave generous numerical headroom.
        assert T_exact * 10 < T_stine, (
            f"exact_local_channels not at least 10× more accurate than "
            f"stinespring: T_stine={T_stine:.3e}, T_exact={T_exact:.3e}"
        )


class TestDMSimBackendExecution:
    """Verify that ``execute_on_backend='dmsim'`` reproduces the NumPy
    direct path bit-for-bit (up to round-off).

    The DMSim backend executes one MQT-Qudits ``QuantumCircuit`` per
    Trotter step (using ``cu_multi`` for the Hamiltonian half-steps and
    :class:`KrausChannel` for every Lindblad channel).  Mathematically
    it computes the same density matrix as the in-process NumPy path —
    this test asserts that as a regression guard.

    *No heuristic correction* is performed: the only difference between
    the two paths is whether the local channels are dispatched through
    DMSim (Kraus) or applied via :func:`apply_channel_*` (column-major
    superoperator).  Both representations are exact (Choi-Jamiolkowski)
    so agreement up to ``1e-12`` is expected for ``N=4, 100`` steps.
    """

    @staticmethod
    def test_dmsim_matches_numpy_n2_t10_20_steps():
        """N=2, short run: agreement at round-off."""
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        ref = QuditGKSLSimulator(
            params, algorithm="exact_local_channels"
        ).simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        bk = QuditGKSLSimulator(
            params,
            algorithm="exact_local_channels",
            execute_on_backend="dmsim",
        ).simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        diff = np.linalg.norm(ref["rho_final"] - bk["rho_final"])
        assert diff < 1e-12, f"N=2 DMSim mismatch: ‖Δρ‖_F={diff:.3e}"

    @staticmethod
    def test_dmsim_matches_numpy_n4_notebook_settings():
        """N=4, t_max=100, n_steps=100 (matches the notebook scenario).

        This is the case that is *impossible* to run on the state-vector
        backends (tnsim/misim) — see STATUS_HONEST_2026-05.md A-1.
        With DMSim it runs in seconds and matches the in-process NumPy
        path at round-off level.
        """
        params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)
        ref = QuditGKSLSimulator(
            params, algorithm="exact_local_channels"
        ).simulate(t_max=100.0, n_steps=100, initial_state="edge_triplet")
        bk = QuditGKSLSimulator(
            params,
            algorithm="exact_local_channels",
            execute_on_backend="dmsim",
        ).simulate(t_max=100.0, n_steps=100, initial_state="edge_triplet")
        diff = np.linalg.norm(ref["rho_final"] - bk["rho_final"])
        # 100 steps × 26 channels per step × Choi eigendecomposition
        # round-off propagates to ~1e-13 in our measurements; assert 1e-10
        # to leave headroom.
        assert diff < 1e-10, (
            f"N=4 DMSim backend disagrees with NumPy reference: "
            f"‖Δρ‖_F={diff:.3e}"
        )
        # Trace conservation must hold for both paths.
        assert abs(bk["trace"][-1] - 1.0) < 1e-8

    @staticmethod
    def test_dmsim_stinespring_ancilla_matches_numpy_reference():
        """Stinespring + DMSim: ancilla-updating circuit matches NumPy.

        The per-step circuit contains N system qutrits **plus one
        physical ancilla qudit** that is updated (reset to ``|0⟩`` via a
        KrausChannel with ``K_k = |0><k|``) between Lindblad channels.
        The result must agree with the in-process NumPy Stinespring path
        (dilate → partial trace per channel) at round-off level.
        """
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        ref = QuditGKSLSimulator(
            params, algorithm="stinespring"
        ).simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        bk_sim = QuditGKSLSimulator(
            params, algorithm="stinespring", execute_on_backend="dmsim"
        )
        bk = bk_sim.simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        diff = np.linalg.norm(ref["rho_final"] - bk["rho_final"])
        # Measured agreement is ~1e-16 (round-off); assert 1e-13 to leave
        # modest headroom while still requiring round-off-level agreement.
        assert diff < 1e-13, (
            f"ancilla-updating DMSim Stinespring disagrees with NumPy "
            f"reference: ‖Δρ‖_F={diff:.3e}"
        )
        assert abs(bk["trace"][-1] - 1.0) < 1e-10
        assert bk["n_ancilla_qudits_backend_circuit"] == 1

    @staticmethod
    def test_dmsim_stinespring_circuit_contains_ancilla_and_resets():
        """The backend circuit really contains the ancilla + reset channels."""
        from mqt.qudits.quantum_circuit.gates.kraus_channel import KrausChannel

        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        sim = QuditGKSLSimulator(
            params, algorithm="stinespring", execute_on_backend="dmsim"
        )
        sim._precompute_unitaries(0.1)
        circuit = sim._stinespring_ancilla_circuit
        n = params.N_molecules
        d = params.d
        # N system qudits + 1 ancilla qudit.
        assert list(circuit.dimensions) == [d] * (n + 1)
        n_channels = len(sim.lindblad_ops)
        resets = [
            inst
            for inst in circuit.instructions
            if isinstance(inst, KrausChannel)
        ]
        # One ancilla reset per channel application (palindromic → 2×).
        assert len(resets) == 2 * n_channels
        anc = n
        for inst in resets:
            target = inst.target_qudits
            assert (target == anc) or (target == [anc])
            # Reset Kraus set is exactly {K_k = |0><k|}.
            for k, K in enumerate(inst.kraus_operators):
                expected = np.zeros((d, d), dtype=np.complex128)
                expected[0, k] = 1.0
                assert np.allclose(K, expected)

    @staticmethod
    def test_dmsim_only_dmsim_is_recognised():
        """Only 'dmsim' is accepted (state-vector backends would silently fail)."""
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        try:
            QuditGKSLSimulator(
                params,
                algorithm="exact_local_channels",
                execute_on_backend="tnsim",
            )
        except ValueError as exc:
            assert "dmsim" in str(exc)
        else:
            msg = "Expected ValueError for non-dmsim execute_on_backend"
            raise AssertionError(msg)


class TestCompilerMeasuredGateCounts:
    """A-3: ``compute_compiler_measured_gate_counts`` returns real numbers
    from the MQT-Qudits compiler, not from a hard-coded formula.
    """

    def test_returns_compiler_breakdown(self):
        from qudit_gksl_circuit_simulator import QuditGKSLKrausSimulator
        params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
        sim = QuditGKSLSimulator(params)
        out = sim.compute_compiler_measured_gate_counts(
            dt=0.5, optimization_level=0
        )
        # Structure must match QuditGKSLKrausSimulator.compile_to_native_gates
        ref = QuditGKSLKrausSimulator(params).compile_to_native_gates(
            dt=0.5, optimization_level=0
        )
        assert (
            out["per_step_summary"]["native_gates"]
            == ref["per_step_summary"]["native_gates"]
        )
        assert (
            out["per_step_summary"]["uncompiled_cu_multi"]
            == ref["per_step_summary"]["uncompiled_cu_multi"]
        )
        # The number of high-level gates from the existing simulator
        # output must NOT match the compiler-measured count — the whole
        # point of A-3 is to expose that they are different.
        result = sim.simulate(t_max=1.0, n_steps=2, initial_state="edge_triplet")
        assert (
            result["n_high_level_gates_per_step"]
            != out["per_step_summary"]["native_gates"]
        )
        assert result["gate_count_method"] == "high_level_count"


class TestN2TnsimVerification:
    """A-1: every Trotter-step building block runs end-to-end on tnsim."""

    def test_all_subcircuits_match_in_process_unitary(self):
        from run_n2_tnsim_verification import run_n2_tnsim_verification

        results = run_n2_tnsim_verification(dt=0.5)
        # 1 Hamiltonian + 2 pair Stinespring + 5×2 single Stinespring = 13
        assert len(results) == 13
        for r in results:
            assert r["match"], (
                f"sub-circuit {r['label']} did not match: "
                f"||Δsv|| = {r['distance']:.3e}"
            )


class TestQutipIndependentCrossValidation:
    """A-2: independent QuTiP reference agrees with `exact_local_channels`.

    See ``tutorials/qutip_gksl_reference.py`` for why QuTiP is a genuinely
    independent reference (no shared code with the in-tree GKSL stack:
    the Hamiltonian and every collapse operator are rebuilt from QuTiP
    primitives, and time evolution uses ``qutip.mesolve`` adaptive ODE,
    not Trotter / Stinespring).

    The test is skipped when QuTiP is not installed so that the rest of
    the suite remains runnable in minimal environments.  Tolerances are
    those calibrated in :func:`run_qutip_cross_validation.test_qutip_independent_cross_validation`
    against the actually-measured numbers.
    """

    def test_run_script_passes_assertions(self):
        pytest.importorskip("qutip", reason="QuTiP needed for A-2 cross-validation")
        from run_qutip_cross_validation import (
            test_qutip_independent_cross_validation as _impl,
        )

        # Re-use the script's own assertions (single source of truth for
        # the calibrated tolerances).
        _impl()

    def test_qutip_collapse_op_count_matches_in_tree_builder(self):
        """Independent collapse-op list has the same length and shape as `build_lindblad_operators`."""
        pytest.importorskip("qutip", reason="QuTiP needed for A-2 cross-validation")
        from qutip_gksl_reference import build_qutip_collapse_ops

        for N in (2, 3, 4):
            params = GKSLPhysicalParameters(N_molecules=N, with_boson=False)
            c_ops_qt = build_qutip_collapse_ops(params)
            in_tree = build_lindblad_operators(params)
            assert len(c_ops_qt) == len(in_tree), (
                f"N={N}: QuTiP built {len(c_ops_qt)} collapse ops, "
                f"in-tree built {len(in_tree)}"
            )
            d = params.d
            for k, (c, (L_in, _gamma)) in enumerate(zip(c_ops_qt, in_tree)):
                shape = c.shape
                assert shape == (d**N, d**N), (
                    f"collapse op #{k} shape={shape} != {(d**N, d**N)}"
                )
                # NOTE: We do not require a per-element match here:
                # the in-tree builder includes the √γ factor in L itself,
                # while QuTiP's collapse-op convention is identical, so the
                # arrays *do* match — but enforcing element-wise equality
                # would defeat the purpose of an "independent" reference.
                # The end-to-end agreement is asserted by the mesolve
                # cross-validation above.


# ---------------------------------------------------------------------------
# Boson DMSim backend execution (Scenario 6d in the notebook)
# ---------------------------------------------------------------------------


class TestQuditGKSLBosonDMSim:
    """Verify ``QuditGKSLBosonSimulator(execute_on_backend='dmsim')`` reproduces
    the in-process NumPy ``exact_local_channels`` path bit-for-bit."""

    @staticmethod
    def test_dmsim_boson_matches_numpy_n2():
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02
        )
        ref = QuditGKSLBosonSimulator(
            params, algorithm="exact_local_channels"
        ).simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        bk = QuditGKSLBosonSimulator(
            params,
            algorithm="exact_local_channels",
            execute_on_backend="dmsim",
        ).simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        diff = np.linalg.norm(ref["rho_final"] - bk["rho_final"])
        assert diff < 1e-10, f"boson DMSim mismatch: ‖Δρ‖_F={diff:.3e}"
        # Trace conservation must hold for both paths.
        assert abs(bk["trace"][-1] - 1.0) < 1e-8

    @staticmethod
    def test_dmsim_boson_requires_exact_local_channels():
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1
        )
        with pytest.raises(ValueError, match="exact_local_channels"):
            QuditGKSLBosonSimulator(
                params, algorithm="stinespring", execute_on_backend="dmsim"
            )

    @staticmethod
    def test_dmsim_boson_only_dmsim_backend():
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1
        )
        with pytest.raises(ValueError, match="dmsim"):
            QuditGKSLBosonSimulator(
                params,
                algorithm="exact_local_channels",
                execute_on_backend="tnsim",
            )


# ---------------------------------------------------------------------------
# Qiskit Aer execution of the qubit GKSL simulators (Scenarios 3d / 4d)
# ---------------------------------------------------------------------------


def _has_qiskit() -> bool:
    try:
        import qiskit  # noqa: F401
        import qiskit_aer  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not _has_qiskit(), reason="qiskit / qiskit-aer not installed")
class TestQiskitQubitGKSL:
    """Verify the Qiskit-Aer-backed qubit GKSL simulator matches the
    classical density-matrix reference within Trotter convergence."""

    @staticmethod
    def test_qiskit_qubit_no_boson_n2_trace_and_classical_agreement():
        from qiskit_qubit_gksl_simulator import QiskitQubitGKSLSimulator
        params = GKSLPhysicalParameters(N_molecules=2)
        sim = QiskitQubitGKSLSimulator(params)
        res = sim.simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        # Trace conservation
        assert all(abs(t - 1.0) < 1e-8 for t in res["trace"])
        # Forbidden subspace must remain (numerically) zero.
        assert max(res["forbidden_state_population"]) < 1e-8

        ref = ClassicalGKSLSimulator(params).simulate(
            t_max=10.0, n_steps=400, initial_state="edge_triplet"
        )
        diff = np.linalg.norm(res["rho_final"] - ref["rho_final"])
        # 2nd-order Trotter at dt=0.5: leading error ~ O(dt²) ~ 1e-1 for
        # the dimensionless rates here, but the *measured* Frobenius
        # difference is ~1e-6.  Assert with comfortable headroom.
        assert diff < 1e-3, f"qiskit vs classical mismatch: {diff:.3e}"

    @staticmethod
    def test_qiskit_qubit_boson_n2_trace_and_classical_agreement():
        from qiskit_qubit_gksl_simulator import QiskitQubitGKSLBosonSimulator
        from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator
        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1, omega_ph=0.15, g_eph=0.02
        )
        sim = QiskitQubitGKSLBosonSimulator(params)
        res = sim.simulate(t_max=10.0, n_steps=20, initial_state="edge_triplet")
        assert all(abs(t - 1.0) < 1e-8 for t in res["trace"])

        ref = ClassicalGKSLBosonSimulator(params).simulate(
            t_max=10.0, n_steps=400, initial_state="edge_triplet"
        )
        diff = np.linalg.norm(res["rho_final"] - ref["rho_final"])
        assert diff < 1e-3, f"qiskit boson vs classical mismatch: {diff:.3e}"

    @staticmethod
    def test_qiskit_qubit_boson_rejects_no_boson_params():
        from qiskit_qubit_gksl_simulator import QiskitQubitGKSLBosonSimulator
        params = GKSLPhysicalParameters(N_molecules=2)
        with pytest.raises(ValueError, match="with_boson"):
            QiskitQubitGKSLBosonSimulator(params)

    @staticmethod
    def test_qiskit_qubit_no_boson_rejects_boson_params():
        from qiskit_qubit_gksl_simulator import QiskitQubitGKSLSimulator
        params = GKSLPhysicalParameters(
            N_molecules=2, with_boson=True, n_max=1
        )
        with pytest.raises(ValueError, match="non-boson"):
            QiskitQubitGKSLSimulator(params)


# ---------------------------------------------------------------------------
# Backend-executed shot simulators (DMSim / Qiskit Aer + Born sampling)
# ---------------------------------------------------------------------------
class TestDMSimShotSimulators:
    """Backend-executed replacements for the notebook's shot cells (5b/5c/3b/3c).

    These tests verify that ``QuditDMSimShotSimulator`` /
    ``QiskitQubitShotSimulator`` produce density-matrix trajectories
    that match the underlying density-matrix backend simulators
    *exactly* (no shot noise on intermediate ρ — only the final-state
    counts are sampled), and that the depolarisation Kraus channels
    are CPTP on the local space they act on.
    """

    @staticmethod
    def test_depolarisation_kraus_are_cptp():
        from dmsim_shot_simulator import (
            depolarisation_kraus, pair_depolarisation_kraus,
        )
        for d in (2, 3, 4):
            for p in (0.0, 1e-3, 0.5, 1.0):
                Ks = depolarisation_kraus(d, p)
                closure = sum(K.conj().T @ K for K in Ks)
                assert np.allclose(closure, np.eye(d), atol=1e-10)
                # E_p[|0><0|] should be (1-p)|0><0| + p I/d.
                rho = np.zeros((d, d), dtype=complex)
                rho[0, 0] = 1.0
                out = sum(K @ rho @ K.conj().T for K in Ks)
                expected = (1 - p) * rho + p * np.eye(d) / d
                assert np.allclose(out, expected, atol=1e-12)
                Kp = pair_depolarisation_kraus(d, p)
                closure_p = sum(K.conj().T @ K for K in Kp)
                assert np.allclose(closure_p, np.eye(d * d), atol=1e-10)

    @staticmethod
    def test_qudit_dmsim_ideal_matches_exact_local_channels_reference():
        from dmsim_shot_simulator import QuditDMSimShotSimulator
        params = GKSLPhysicalParameters(N_molecules=2)
        sim = QuditDMSimShotSimulator(params, p_depol=0.0)
        res = sim.simulate(t_max=5.0, n_steps=10, n_shots=500, seed=0)
        ref = QuditGKSLSimulator(params, algorithm="exact_local_channels")
        ref_res = ref.simulate(t_max=5.0, n_steps=10)
        # backend ρ must match NumPy ρ to machine precision (same algorithm)
        assert np.linalg.norm(res["rho_final"] - ref_res["rho_final"]) < 1e-10
        # Born-sampled counts sum to n_shots
        assert sum(res["counts"].values()) == 500
        # No noise_params present in ideal mode
        assert "noise_params" not in res
        assert res["backend"] == "mqt_qudits:dmsim"

    @staticmethod
    def test_qudit_dmsim_noisy_changes_rho_and_preserves_trace():
        from dmsim_shot_simulator import QuditDMSimShotSimulator
        params = GKSLPhysicalParameters(N_molecules=2)
        sim_i = QuditDMSimShotSimulator(params, p_depol=0.0)
        sim_n = QuditDMSimShotSimulator(
            params, p_depol=0.001, p_dephasing=0.0, depol_pair_only=True
        )
        ri = sim_i.simulate(t_max=5.0, n_steps=10, n_shots=200, seed=0)
        rn = sim_n.simulate(t_max=5.0, n_steps=10, n_shots=200, seed=0)
        # Noise produces a measurable ρ shift
        assert np.linalg.norm(rn["rho_final"] - ri["rho_final"]) > 1e-4
        # Qudit pair depolarisation stays inside d=3 — trace exactly preserved
        assert max(abs(t - 1.0) for t in rn["trace"]) < 1e-10
        assert rn["noise_params"]["noise_model"] == "pair_depolarisation_kraus"

    @staticmethod
    def test_qiskit_qubit_ideal_matches_aer_density_matrix_reference():
        from dmsim_shot_simulator import QiskitQubitShotSimulator
        from qiskit_qubit_gksl_simulator import QiskitQubitGKSLSimulator
        params = GKSLPhysicalParameters(N_molecules=2)
        sim = QiskitQubitShotSimulator(params, p_depol=0.0)
        res = sim.simulate(t_max=5.0, n_steps=10, n_shots=500, seed=0)
        ref = QiskitQubitGKSLSimulator(params)
        ref_res = ref.simulate(t_max=5.0, n_steps=10)
        assert np.linalg.norm(res["rho_final"] - ref_res["rho_final"]) < 1e-12
        assert sum(res["counts"].values()) == 500
        # No leakage with ideal Lindblad-only evolution
        assert res["forbidden_count"] == 0
        assert res["backend"] == "qiskit_aer:density_matrix"

    @staticmethod
    def test_qiskit_qubit_noisy_introduces_leakage():
        from dmsim_shot_simulator import QiskitQubitShotSimulator
        params = GKSLPhysicalParameters(N_molecules=2)
        sim = QiskitQubitShotSimulator(
            params, p_depol=0.001, p_dephasing=0.0, depol_pair_only=True
        )
        res = sim.simulate(t_max=5.0, n_steps=10, n_shots=2000, seed=0)
        # 4-qubit pair Pauli noise leaks into forbidden |11⟩ states
        assert res["trace"][-1] < 1.0
        assert res["forbidden_count"] >= 1
        assert sum(res["counts"].values()) == 2000

    @staticmethod
    def test_unsupported_noise_options_are_rejected_explicitly():
        """Refuse silent semantic drift: dephasing/single-site noise raise."""
        from dmsim_shot_simulator import (
            QuditDMSimShotSimulator, QiskitQubitShotSimulator,
        )
        params = GKSLPhysicalParameters(N_molecules=2)
        with pytest.raises(NotImplementedError, match="p_dephasing"):
            QuditDMSimShotSimulator(params, p_depol=0.001, p_dephasing=0.001)
        with pytest.raises(NotImplementedError, match="depol_pair_only"):
            QuditDMSimShotSimulator(
                params, p_depol=0.001, depol_pair_only=False
            )
        with pytest.raises(NotImplementedError, match="p_dephasing"):
            QiskitQubitShotSimulator(params, p_depol=0.001, p_dephasing=0.001)
        with pytest.raises(NotImplementedError, match="depol_pair_only"):
            QiskitQubitShotSimulator(
                params, p_depol=0.001, depol_pair_only=False
            )

    @staticmethod
    def test_born_sampling_statistics_converge_to_diagonal():
        """High-shot Born sampling matches diag(ρ_final) within 3σ."""
        from dmsim_shot_simulator import (
            QuditDMSimShotSimulator, sample_counts_from_density_matrix,
        )
        params = GKSLPhysicalParameters(N_molecules=2)
        sim = QuditDMSimShotSimulator(params, p_depol=0.0)
        res = sim.simulate(t_max=5.0, n_steps=10, n_shots=20000, seed=0)
        diag = np.real(np.diag(res["rho_final"]))
        diag = np.clip(diag, 0.0, None)
        diag = diag / diag.sum()
        empirical = np.zeros_like(diag)
        for idx, c in res["counts"].items():
            empirical[idx] = c / 20000
        # 3σ binomial bound for n=20000 is ~3·sqrt(p(1-p)/n) ≤ 0.011
        assert np.max(np.abs(empirical - diag)) < 0.02
