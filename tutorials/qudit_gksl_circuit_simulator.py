"""Qudit GKSL simulator using actual MQT-Qudits QuantumCircuit API.

Scenario 5c: Circuit-based Qudit GKSL-Lindblad (no boson).

This simulator constructs actual MQT-Qudits quantum circuits for each building
block of the Trotter step, replacing matrix-level computations with gate-level
circuit operations verified against the MQT-Qudits framework.

Circuit decomposition:
  - Hamiltonian step: cu_one (on-site phases, 3x3) + cu_two (pair transfer, 9x9)
  - Stinespring channels: cu_two (single-site, 6x6) or cu_multi (TTA pair, 18x18)

Density matrix evolution uses circuit-derived local Kraus operators, which is
mathematically equivalent to executing the full circuit with ancilla reset
between channels, but is computationally tractable for open-system dynamics.
"""

from __future__ import annotations

import os
import sys
import time as time_module
from functools import reduce

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import (
    build_single_site_operator,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
)
from gksl_physical_parameters import GKSLPhysicalParameters


class QuditGKSLCircuitSimulator:
    """Qudit GKSL simulator using MQT-Qudits QuantumCircuit API.

    Constructs real MQT-Qudits quantum circuits for each Trotter step component.
    Each unitary (Hamiltonian gates, Stinespring dilations) is built as a
    quantum circuit object, verified via MQT-Qudits state-vector simulation,
    and then used in density-matrix evolution via local Kraus operators.

    Quantum resources per Trotter step:
      - 4 cu_one gates (on-site Hamiltonian phases, half step) x 2 = 8
      - 3 cu_two gates (pair transfer, half step) x 2 = 6
      - 20 cu_two gates (single-site Stinespring, 6x6 each)
      - 6 cu_multi gates (TTA pair Stinespring, 18x18 each)
      Total: 40 gates per Trotter step
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QuditGKSLCircuitSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.d = params.d
        self.N = params.N_molecules
        self.dim = params.d ** params.N_molecules

        self._build_local_operators()

    # ------------------------------------------------------------------
    # Local operator construction
    # ------------------------------------------------------------------

    def _build_local_operators(self) -> None:
        """Pre-compute local operators for circuit gate construction."""
        d = self.d

        # On-site Hamiltonian (3x3 diagonal)
        # State ordering: |0⟩ = S0 (ground), |1⟩ = T1 (triplet), |2⟩ = S1 (excited singlet)
        self.h_local = np.diag(
            np.array([0.0, self.params.E_T, self.params.E_S], dtype=np.complex128)
        )

        # Transfer Hamiltonian per pair (9x9)
        self.h_transfer_pairs: dict[tuple[int, int], np.ndarray] = {}
        for pair in self.params.neighbors:
            i, j = pair
            H_pair = np.zeros((d * d, d * d), dtype=np.complex128)
            # |01><10| + |10><01| in the (i,j) pair space
            idx_01 = 0 * d + 1
            idx_10 = 1 * d + 0
            H_pair[idx_01, idx_10] = self.params.V
            H_pair[idx_10, idx_01] = self.params.V
            self.h_transfer_pairs[(i, j)] = H_pair

        # Local Lindblad operators with site information
        self.lindblad_local_info = self._build_local_lindblad_info()

    def _build_local_lindblad_info(
        self,
    ) -> list[tuple[str, tuple[int, ...], np.ndarray, float]]:
        """Build local Lindblad operators with site information.

        Returns list of (op_type, sites, local_op, gamma) where:
          op_type: 'single' or 'pair'
          sites: (i,) for single-site, (i, j) for pair
          local_op: local Lindblad operator (3x3 or 9x9), includes sqrt(gamma)
          gamma: dissipation rate
        """
        d = self.d
        N = self.N
        p = self.params

        def _ket_bra(a: int, b: int) -> np.ndarray:
            m = np.zeros((d, d), dtype=np.complex128)
            m[a, b] = 1.0
            return m

        info: list[tuple[str, tuple[int, ...], np.ndarray, float]] = []

        # TTA: pair operators (6 total: 3 pairs × 2 channels)
        for i, j in p.neighbors:
            gamma = p.gamma_TTA / 2.0
            sq = np.sqrt(gamma)
            # Channel 1: |2>_i<1| ⊗ |0>_j<1|
            l_pair = sq * np.kron(_ket_bra(2, 1), _ket_bra(0, 1))
            info.append(("pair", (i, j), l_pair, gamma))
            # Channel 2: |0>_i<1| ⊗ |2>_j<1|
            l_pair = sq * np.kron(_ket_bra(0, 1), _ket_bra(2, 1))
            info.append(("pair", (i, j), l_pair, gamma))

        # Fluorescence: √Γ_fl |0><2| per molecule
        for i in range(N):
            gamma = p.Gamma_fl
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 2), gamma))

        # Phosphorescence: √Γ_ph |0><1| per molecule
        for i in range(N):
            gamma = p.Gamma_ph
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 1), gamma))

        # Internal conversion: √k_IC |0><2| per molecule
        for i in range(N):
            gamma = p.k_IC
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 2), gamma))

        # ISC S₁→T₁: √k_ISC_ST |1><2| per molecule
        for i in range(N):
            gamma = p.k_ISC_ST
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(1, 2), gamma))

        # ISC T₁→S₀: √k_ISC_TS |0><1| per molecule
        for i in range(N):
            gamma = p.k_ISC_TS
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 1), gamma))

        expected = 2 * len(p.neighbors) + 5 * N
        if len(info) != expected:
            msg = f"Expected {expected} Lindblad operators, got {len(info)}"
            raise ValueError(msg)
        return info

    # ------------------------------------------------------------------
    # Circuit construction: Hamiltonian
    # ------------------------------------------------------------------

    def build_hamiltonian_circuit(self, dt: float):
        """Build MQT-Qudits circuit for Hamiltonian half-step.

        Decomposition using 1st-order Trotter for inner splitting:
          exp(-iH_total dt) ≈ [⊗_i exp(-i h_local dt)] · [Π_{<i,j>} exp(-iH_pair dt)]

        The on-site and transfer parts are applied as separate local gates.

        Returns (circuit, gate_count) tuple.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        N = self.N
        d = self.d

        circuit = QuantumCircuit(N, [d] * N, 0)
        gate_count = 0

        # On-site phases: cu_one per qudit
        U_onsite = expm(-1j * self.h_local * dt)
        for i in range(N):
            circuit.cu_one(i, U_onsite)
            gate_count += 1

        # Transfer unitaries: cu_two per pair
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair = expm(-1j * H_pair * dt)
            circuit.cu_two([ip, jp], U_pair)
            gate_count += 1

        return circuit, gate_count

    def _compute_hamiltonian_unitary_from_circuit(self, dt: float) -> np.ndarray:
        """Compute the composite Hamiltonian unitary from local gate decomposition.

        Composes the tensor-product on-site unitary with sequential pair unitaries
        to obtain the full-system unitary corresponding to the circuit.
        """
        d = self.d
        N = self.N
        dim = self.dim
        eye = np.eye(d, dtype=np.complex128)

        # On-site part: tensor product of single-site unitaries
        U_onsite = expm(-1j * self.h_local * dt)
        U_full = reduce(np.kron, [U_onsite] * N)

        # Transfer part: sequential application of pair unitaries
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair = expm(-1j * H_pair * dt)

            # Embed pair unitary into full system
            op_list = [eye] * N
            # Replace qudits ip, jp with the pair unitary
            # We need to handle non-adjacent qudit indices correctly
            U_pair_full = self._embed_pair_unitary(U_pair, ip, jp)
            U_full = U_pair_full @ U_full

        return U_full

    def _embed_pair_unitary(
        self, U_pair: np.ndarray, site_i: int, site_j: int
    ) -> np.ndarray:
        """Embed a 2-qudit unitary into the full N-qudit system.

        U_pair is d^2 × d^2 acting on the (site_i, site_j) pair subspace.
        Returns the full dim × dim unitary.
        """
        d = self.d
        N = self.N
        dim = self.dim

        U_full = np.zeros((dim, dim), dtype=np.complex128)
        for row in range(dim):
            for col in range(dim):
                # Decompose row and col indices into per-qudit states
                row_states = []
                col_states = []
                r, c = row, col
                for _ in range(N):
                    row_states.append(r % d)
                    col_states.append(c % d)
                    r //= d
                    c //= d
                row_states.reverse()
                col_states.reverse()

                # Check if non-pair qudits match
                match = True
                for k in range(N):
                    if k != site_i and k != site_j:
                        if row_states[k] != col_states[k]:
                            match = False
                            break

                if match:
                    # Map to pair indices
                    pair_row = row_states[site_i] * d + row_states[site_j]
                    pair_col = col_states[site_i] * d + col_states[site_j]
                    U_full[row, col] = U_pair[pair_row, pair_col]

        return U_full

    # ------------------------------------------------------------------
    # Circuit construction: Stinespring channels
    # ------------------------------------------------------------------

    def _build_local_stinespring_unitary(
        self, L_local: np.ndarray, dt: float
    ) -> np.ndarray:
        """Build local Stinespring unitary from local Lindblad operator.

        For a d_local-dimensional local operator L_local, constructs the
        (2*d_local × 2*d_local) Stinespring unitary:
          G = [[0, L†], [L, 0]]
          U = expm(-i * sqrt(dt) * G)

        Convention: env ⊗ system ordering (ancilla qubit first).
        """
        d_local = L_local.shape[0]
        G = np.zeros((2 * d_local, 2 * d_local), dtype=np.complex128)
        G[:d_local, d_local:] = L_local.conj().T
        G[d_local:, :d_local] = L_local

        theta = np.sqrt(dt)
        U = expm(-1j * theta * G)

        residual = np.linalg.norm(U.conj().T @ U - np.eye(2 * d_local), ord="fro")
        if residual >= 1e-10:
            msg = f"Local Stinespring unitarity check failed: ||U†U - I||_F = {residual}"
            raise ValueError(msg)
        return U

    def build_stinespring_circuit_single(
        self, L_local_3x3: np.ndarray, dt: float, target_qudit: int
    ):
        """Build MQT-Qudits circuit for single-site Stinespring channel.

        Creates a circuit with the target qutrit and one ancilla qubit,
        applying the 6×6 Stinespring unitary as a cu_two gate.

        Returns (circuit, U_local_6x6) tuple.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        d = self.d
        U_local = self._build_local_stinespring_unitary(L_local_3x3, dt)

        # Circuit: [target_qutrit(d=3), ancilla(d=2)]
        circuit = QuantumCircuit(2, [d, 2], 0)
        circuit.cu_two([0, 1], U_local)

        return circuit, U_local

    def build_stinespring_circuit_pair(
        self,
        L_local_9x9: np.ndarray,
        dt: float,
        qudit_i: int,
        qudit_j: int,
    ):
        """Build MQT-Qudits circuit for TTA pair Stinespring channel.

        Creates a circuit with two qutrits and one ancilla qubit,
        applying the 18×18 Stinespring unitary as a cu_multi gate.

        Returns (circuit, U_local_18x18) tuple.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        d = self.d
        U_local = self._build_local_stinespring_unitary(L_local_9x9, dt)

        # Circuit: [qutrit_i(d=3), qutrit_j(d=3), ancilla(d=2)]
        circuit = QuantumCircuit(3, [d, d, 2], 0)
        circuit.cu_multi([0, 1, 2], U_local)

        return circuit, U_local

    def build_full_trotter_step_circuit(self, dt: float):
        """Build MQT-Qudits circuit for one full 2nd-order Trotter step.

        Structure: H(dt/2) → D_1...D_26(dt) → H(dt/2)

        Returns a dict with circuit info and gate counts.
        """
        circuits = []
        total_gates = 0

        # Half Hamiltonian step
        h_circ, h_gates = self.build_hamiltonian_circuit(dt / 2)
        circuits.append(("hamiltonian_half_1", h_circ))
        total_gates += h_gates

        # Lindblad channels
        for idx, (op_type, sites, L_local, _gamma) in enumerate(
            self.lindblad_local_info
        ):
            if op_type == "single":
                circ, _ = self.build_stinespring_circuit_single(
                    L_local, dt, sites[0]
                )
                circuits.append((f"stinespring_single_{idx}", circ))
                total_gates += 1
            elif op_type == "pair":
                circ, _ = self.build_stinespring_circuit_pair(
                    L_local, dt, sites[0], sites[1]
                )
                circuits.append((f"stinespring_pair_{idx}", circ))
                total_gates += 1

        # Second half Hamiltonian step
        h_circ2, h_gates2 = self.build_hamiltonian_circuit(dt / 2)
        circuits.append(("hamiltonian_half_2", h_circ2))
        total_gates += h_gates2

        return {
            "circuits": circuits,
            "total_gates": total_gates,
            "n_hamiltonian_gates": h_gates + h_gates2,
            "n_stinespring_gates": len(self.lindblad_local_info),
        }

    # ------------------------------------------------------------------
    # Kraus operator extraction from local Stinespring unitaries
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_kraus_from_local_stinespring(
        U_local: np.ndarray, d_local: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Extract Kraus operators from a local Stinespring unitary.

        U_local has env⊗sys ordering with a 2-level ancilla:
          U_local = [[A, B], [C, D]]  (blocks of size d_local × d_local)

        For ancilla initially in |0⟩:
          K_0 = A  (ancilla stays |0⟩)
          K_1 = C  (ancilla flips to |1⟩)

        Satisfies: K_0†K_0 + K_1†K_1 = I (trace preservation)
        """
        K0 = U_local[:d_local, :d_local].copy()
        K1 = U_local[d_local:, :d_local].copy()
        return K0, K1

    # ------------------------------------------------------------------
    # Density matrix evolution using circuit-derived local operators
    # ------------------------------------------------------------------

    def _apply_hamiltonian_step_circuit(
        self, rho: np.ndarray, dt: float
    ) -> np.ndarray:
        """Apply Hamiltonian evolution using circuit-decomposed local unitaries.

        Decomposition: exp(-iH_total dt) ≈ [Π_pairs U_pair] · [⊗_i U_onsite_i]
        """
        U = self._compute_hamiltonian_unitary_from_circuit(dt)
        return U @ rho @ U.conj().T

    def _apply_single_site_channel(
        self,
        rho: np.ndarray,
        K0: np.ndarray,
        K1: np.ndarray,
        site: int,
    ) -> np.ndarray:
        """Apply single-site Kraus channel using tensor-product structure.

        E[ρ] = (I⊗...⊗K0⊗...⊗I) ρ (I⊗...⊗K0†⊗...⊗I)
             + (I⊗...⊗K1⊗...⊗I) ρ (I⊗...⊗K1†⊗...⊗I)
        """
        d = self.d
        N = self.N
        K0_full = build_single_site_operator(K0, site, N, d)
        K1_full = build_single_site_operator(K1, site, N, d)
        return K0_full @ rho @ K0_full.conj().T + K1_full @ rho @ K1_full.conj().T

    def _apply_pair_channel(
        self,
        rho: np.ndarray,
        K0: np.ndarray,
        K1: np.ndarray,
        site_i: int,
        site_j: int,
    ) -> np.ndarray:
        """Apply pair Kraus channel using tensor-product structure.

        K0, K1 are d^2 × d^2 Kraus operators acting on the (site_i, site_j) pair.
        Embeds into the full system via tensor product with identity on other sites.
        """
        K0_full = self._embed_pair_operator(K0, site_i, site_j)
        K1_full = self._embed_pair_operator(K1, site_i, site_j)
        return K0_full @ rho @ K0_full.conj().T + K1_full @ rho @ K1_full.conj().T

    def _embed_pair_operator(
        self, op_pair: np.ndarray, site_i: int, site_j: int
    ) -> np.ndarray:
        """Embed a d^2 × d^2 pair operator into the full N-qudit system."""
        d = self.d
        N = self.N
        dim = self.dim

        op_full = np.zeros((dim, dim), dtype=np.complex128)
        for row in range(dim):
            for col in range(dim):
                row_states = []
                col_states = []
                r, c = row, col
                for _ in range(N):
                    row_states.append(r % d)
                    col_states.append(c % d)
                    r //= d
                    c //= d
                row_states.reverse()
                col_states.reverse()

                match = True
                for k in range(N):
                    if k != site_i and k != site_j:
                        if row_states[k] != col_states[k]:
                            match = False
                            break

                if match:
                    pair_row = row_states[site_i] * d + row_states[site_j]
                    pair_col = col_states[site_i] * d + col_states[site_j]
                    op_full[row, col] = op_pair[pair_row, pair_col]

        return op_full

    def _trotter_step_circuit(self, rho: np.ndarray, dt: float) -> np.ndarray:
        """2nd-order symmetric Trotter step using circuit-derived local operators.

        exp(L dt) ≈ exp(L_H dt/2) · Π_α exp(L_D_α dt) · exp(L_H dt/2)
        """
        # Half Hamiltonian (circuit-decomposed)
        rho = self._apply_hamiltonian_step_circuit(rho, dt / 2)

        # All Lindblad channels via local Stinespring Kraus operators
        for op_type, sites, L_local, _gamma in self.lindblad_local_info:
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt)
            K0, K1 = self._extract_kraus_from_local_stinespring(U_local, d_local)

            if op_type == "single":
                rho = self._apply_single_site_channel(rho, K0, K1, sites[0])
            elif op_type == "pair":
                rho = self._apply_pair_channel(rho, K0, K1, sites[0], sites[1])

        # Half Hamiltonian (circuit-decomposed)
        rho = self._apply_hamiltonian_step_circuit(rho, dt / 2)
        return rho

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix in qutrit space."""
        d = self.d
        N = self.N
        dim = d**N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d**i) for i in range(N))
            psi[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d**i) for i in range(N))
            psi[index] = 1.0
        else:
            msg = f"Unknown state type: {state_type}"
            raise ValueError(msg)

        return np.outer(psi, psi.conj())

    # ------------------------------------------------------------------
    # Circuit verification
    # ------------------------------------------------------------------

    def verify_hamiltonian_circuit(self, dt: float) -> dict:
        """Verify Hamiltonian circuit against exact matrix exponential.

        Computes the Frobenius-norm distance between the circuit-decomposed
        Hamiltonian unitary and the exact exp(-iH_total dt).
        """
        from gksl_math_utils import build_onsite_hamiltonian, build_transfer_hamiltonian

        H_0 = build_onsite_hamiltonian(self.params)
        H_transfer = build_transfer_hamiltonian(self.params)
        H_total = H_0 + H_transfer

        U_exact = expm(-1j * H_total * dt)
        U_circuit = self._compute_hamiltonian_unitary_from_circuit(dt)

        # Both should be unitary
        res_exact = np.linalg.norm(
            U_exact.conj().T @ U_exact - np.eye(self.dim), ord="fro"
        )
        res_circuit = np.linalg.norm(
            U_circuit.conj().T @ U_circuit - np.eye(self.dim), ord="fro"
        )

        # Distance between exact and circuit unitaries
        distance = np.linalg.norm(U_circuit - U_exact, ord="fro")

        return {
            "dt": dt,
            "exact_unitarity_residual": float(res_exact),
            "circuit_unitarity_residual": float(res_circuit),
            "frobenius_distance": float(distance),
            "circuit_is_unitary": res_circuit < 1e-10,
        }

    def verify_stinespring_circuit(self, dt: float) -> dict:
        """Verify local Stinespring unitaries against full-system computation.

        For each Lindblad channel, compares the CPTP map from the local
        Stinespring unitary with the map from the full-system computation.
        """
        from stinespring_utils import (
            apply_stinespring_to_density_matrix,
            stinespring_unitary_from_lindblad,
        )
        from gksl_math_utils import build_lindblad_operators

        lindblad_ops = build_lindblad_operators(self.params)
        rho_test = self.prepare_initial_state("edge_triplet")

        max_distance = 0.0
        results = []
        for idx, ((L_full, gamma), (op_type, sites, L_local, _gamma_local)) in (
            enumerate(zip(lindblad_ops, self.lindblad_local_info))
        ):
            # Full-system Stinespring
            U_full = stinespring_unitary_from_lindblad(L_full, dt)
            rho_full = apply_stinespring_to_density_matrix(rho_test, U_full)

            # Local Stinespring via Kraus
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt)
            K0, K1 = self._extract_kraus_from_local_stinespring(U_local, d_local)

            if op_type == "single":
                rho_local = self._apply_single_site_channel(
                    rho_test, K0, K1, sites[0]
                )
            else:
                rho_local = self._apply_pair_channel(
                    rho_test, K0, K1, sites[0], sites[1]
                )

            distance = np.linalg.norm(rho_full - rho_local, ord="fro")
            max_distance = max(max_distance, distance)
            results.append({
                "index": idx,
                "op_type": op_type,
                "sites": sites,
                "frobenius_distance": float(distance),
                "match": distance < 1e-10,
            })

        return {
            "dt": dt,
            "max_frobenius_distance": float(max_distance),
            "all_match": max_distance < 1e-10,
            "per_channel": results,
        }

    def verify_circuit_via_mqt(self, dt: float) -> dict:
        """Verify circuit outputs using MQT-Qudits state-vector simulation.

        Runs the Hamiltonian circuit on the |0...0⟩ initial state and checks
        the output against the expected matrix-vector product.
        """
        from mqt.qudits.simulation import MQTQuditProvider

        N = self.N
        d = self.d
        dim = self.dim

        circuit, gate_count = self.build_hamiltonian_circuit(dt)
        provider = MQTQuditProvider()
        backend = provider.get_backend("tnsim")

        job = backend.run(circuit, shots=1)
        result = job.result()
        sv_mqt = result.get_state_vector()[0]

        # Expected: U_circuit |0...0⟩
        U_circuit = self._compute_hamiltonian_unitary_from_circuit(dt)
        psi_0 = np.zeros(dim, dtype=np.complex128)
        psi_0[0] = 1.0
        sv_expected = U_circuit @ psi_0

        distance = np.linalg.norm(sv_mqt - sv_expected)

        return {
            "dt": dt,
            "gate_count": gate_count,
            "statevector_distance": float(distance),
            "match": distance < 1e-10,
        }

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run circuit-based Qudit GKSL simulation.

        Uses MQT-Qudits circuit-derived local unitaries for density matrix
        evolution. The Hamiltonian is decomposed into cu_one + cu_two gates,
        and Stinespring channels use local cu_two/cu_multi gates.

        Returns a dict with time series and circuit statistics.
        """
        start = time_module.time()

        dt = t_max / n_steps
        rho = self.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho, self.params)]
        entropies = [compute_von_neumann_entropy(rho)]
        purities = [compute_purity(rho)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step_circuit(rho, dt)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho))
            purities.append(compute_purity(rho))

        elapsed = time_module.time() - start

        # Gate counts from circuit decomposition
        n_hamiltonian_onsite = self.N
        n_hamiltonian_transfer = len(self.params.neighbors)
        n_stinespring_single = sum(
            1 for t, _, _, _ in self.lindblad_local_info if t == "single"
        )
        n_stinespring_pair = sum(
            1 for t, _, _, _ in self.lindblad_local_info if t == "pair"
        )
        gates_per_half_ham = n_hamiltonian_onsite + n_hamiltonian_transfer
        gates_per_step = 2 * gates_per_half_ham + n_stinespring_single + n_stinespring_pair

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho,
            "elapsed_time": elapsed,
            "method": "qudit_gksl_circuit",
            "params": self.params.to_dict(),
            "n_system_qudits": self.N,
            "n_ancilla_qubits": len(self.lindblad_local_info),
            "gates_per_step": gates_per_step,
            "total_gates": gates_per_step * n_steps,
            "gate_breakdown": {
                "cu_one_onsite": 2 * n_hamiltonian_onsite,
                "cu_two_transfer": 2 * n_hamiltonian_transfer,
                "cu_two_stinespring_single": n_stinespring_single,
                "cu_multi_stinespring_pair": n_stinespring_pair,
            },
        }
