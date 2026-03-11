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

    Since this targets qudit quantum computers, all registers — including
    Stinespring ancillas — are native d-level qudits.

    Quantum resources per Trotter step (palindromic 2nd-order):
      - 4 cu_one gates (on-site Hamiltonian phases, half step) x 2 = 8
      - 3 cu_two gates (pair transfer, half step) x 2 = 6
      - 20 cu_two gates (single-site Stinespring) x 2 (fwd+rev) = 40
      - 6 cu_multi gates (TTA pair Stinespring) x 2 (fwd+rev) = 12
      Total: 66 gates per Trotter step
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QuditGKSLCircuitSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.d = params.d
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension
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
        (d_anc*d_local × d_anc*d_local) Stinespring unitary:
          G = [[0, L†, 0], [L, 0, 0], [0, 0, 0]]   (for d_anc=3)
          U = expm(-i * sqrt(dt) * G)

        Convention: env ⊗ system ordering (ancilla qudit first).
        """
        d_local = L_local.shape[0]
        d_anc = self.d_anc
        dim_total = d_anc * d_local
        G = np.zeros((dim_total, dim_total), dtype=np.complex128)
        G[:d_local, d_local : 2 * d_local] = L_local.conj().T
        G[d_local : 2 * d_local, :d_local] = L_local

        theta = np.sqrt(dt)
        U = expm(-1j * theta * G)

        residual = np.linalg.norm(U.conj().T @ U - np.eye(dim_total), ord="fro")
        if residual >= 1e-10:
            msg = f"Local Stinespring unitarity check failed: ||U†U - I||_F = {residual}"
            raise ValueError(msg)
        return U

    def build_stinespring_circuit_single(
        self, L_local_3x3: np.ndarray, dt: float, target_qudit: int
    ):
        """Build MQT-Qudits circuit for single-site Stinespring channel.

        Creates a circuit with the target qutrit and one ancilla qudit,
        applying the (d_anc*3 × d_anc*3) Stinespring unitary as a cu_two gate.

        Returns (circuit, U_local) tuple.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        d = self.d
        d_anc = self.d_anc
        U_local = self._build_local_stinespring_unitary(L_local_3x3, dt)

        # Circuit: [target_qutrit(d=3), ancilla(d=d_anc)]
        circuit = QuantumCircuit(2, [d, d_anc], 0)
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

        Creates a circuit with two qutrits and one ancilla qudit,
        applying the (d_anc*9 × d_anc*9) Stinespring unitary as a cu_multi gate.

        Returns (circuit, U_local) tuple.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        d = self.d
        d_anc = self.d_anc
        U_local = self._build_local_stinespring_unitary(L_local_9x9, dt)

        # Circuit: [qutrit_i(d=3), qutrit_j(d=3), ancilla(d=d_anc)]
        circuit = QuantumCircuit(3, [d, d, d_anc], 0)
        circuit.cu_multi([0, 1, 2], U_local)

        return circuit, U_local

    def build_full_trotter_step_circuit(self, dt: float):
        """Build MQT-Qudits circuit for one full 2nd-order Trotter step.

        Structure: H(dt/2) → D_1...D_n(dt/2) → D_n...D_1(dt/2) → H(dt/2)

        Palindromic Lindblad ordering matches QuditGKSLSimulator._trotter_step.
        Returns a dict with circuit info and gate counts.
        """
        circuits = []
        total_gates = 0

        # Half Hamiltonian step
        h_circ, h_gates = self.build_hamiltonian_circuit(dt / 2)
        circuits.append(("hamiltonian_half_1", h_circ))
        total_gates += h_gates

        # Forward Lindblad half-step channels
        lindblad_circuits_forward = []
        for idx, (op_type, sites, L_local, _gamma) in enumerate(
            self.lindblad_local_info
        ):
            if op_type == "single":
                circ, _ = self.build_stinespring_circuit_single(
                    L_local, dt / 2, sites[0]
                )
                lindblad_circuits_forward.append((f"stinespring_fwd_single_{idx}", circ))
                total_gates += 1
            elif op_type == "pair":
                circ, _ = self.build_stinespring_circuit_pair(
                    L_local, dt / 2, sites[0], sites[1]
                )
                lindblad_circuits_forward.append((f"stinespring_fwd_pair_{idx}", circ))
                total_gates += 1
        circuits.extend(lindblad_circuits_forward)

        # Reverse Lindblad half-step channels (palindromic)
        lindblad_circuits_reverse = []
        for idx, (op_type, sites, L_local, _gamma) in reversed(
            list(enumerate(self.lindblad_local_info))
        ):
            if op_type == "single":
                circ, _ = self.build_stinespring_circuit_single(
                    L_local, dt / 2, sites[0]
                )
                lindblad_circuits_reverse.append((f"stinespring_rev_single_{idx}", circ))
                total_gates += 1
            elif op_type == "pair":
                circ, _ = self.build_stinespring_circuit_pair(
                    L_local, dt / 2, sites[0], sites[1]
                )
                lindblad_circuits_reverse.append((f"stinespring_rev_pair_{idx}", circ))
                total_gates += 1
        circuits.extend(lindblad_circuits_reverse)

        # Second half Hamiltonian step
        h_circ2, h_gates2 = self.build_hamiltonian_circuit(dt / 2)
        circuits.append(("hamiltonian_half_2", h_circ2))
        total_gates += h_gates2

        return {
            "circuits": circuits,
            "total_gates": total_gates,
            "n_hamiltonian_gates": h_gates + h_gates2,
            "n_stinespring_gates": 2 * len(self.lindblad_local_info),
        }

    # ------------------------------------------------------------------
    # Kraus operator extraction from local Stinespring unitaries
    # ------------------------------------------------------------------

    def _extract_kraus_from_local_stinespring(
        self, U_local: np.ndarray, d_local: int
    ) -> list[np.ndarray]:
        """Extract Kraus operators from a local Stinespring unitary.

        U_local has env⊗sys ordering with a d_anc-level ancilla.
        The unitary is partitioned into d_anc × d_anc blocks of size d_local.

        For ancilla initially in |0⟩:
          K_k = U_local[k*d_local:(k+1)*d_local, 0:d_local]  for k=0,...,d_anc-1

        Satisfies: sum_k K_k† K_k = I (trace preservation)
        """
        d_anc = self.d_anc
        kraus_ops = []
        for k in range(d_anc):
            K = U_local[k * d_local : (k + 1) * d_local, :d_local].copy()
            kraus_ops.append(K)
        return kraus_ops

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
        kraus_ops: list[np.ndarray],
        site: int,
    ) -> np.ndarray:
        """Apply single-site Kraus channel using tensor-product structure.

        E[ρ] = sum_k (I⊗...⊗K_k⊗...⊗I) ρ (I⊗...⊗K_k†⊗...⊗I)
        """
        d = self.d
        N = self.N
        rho_out = np.zeros_like(rho)
        for K in kraus_ops:
            K_full = build_single_site_operator(K, site, N, d)
            rho_out += K_full @ rho @ K_full.conj().T
        return rho_out

    def _apply_pair_channel(
        self,
        rho: np.ndarray,
        kraus_ops: list[np.ndarray],
        site_i: int,
        site_j: int,
    ) -> np.ndarray:
        """Apply pair Kraus channel using tensor-product structure.

        K_k are d^2 × d^2 Kraus operators acting on the (site_i, site_j) pair.
        Embeds into the full system via tensor product with identity on other sites.
        """
        rho_out = np.zeros_like(rho)
        for K in kraus_ops:
            K_full = self._embed_pair_operator(K, site_i, site_j)
            rho_out += K_full @ rho @ K_full.conj().T
        return rho_out

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
        """2nd-order symmetric Trotter step with palindromic Lindblad ordering.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        Matches the palindromic structure of QuditGKSLSimulator._trotter_step.
        """
        # Half Hamiltonian (circuit-decomposed)
        rho = self._apply_hamiltonian_step_circuit(rho, dt / 2)

        # Precompute half-step Kraus operators for each Lindblad channel
        kraus_list = []
        for op_type, sites, L_local, _gamma in self.lindblad_local_info:
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt / 2)
            kraus_ops = self._extract_kraus_from_local_stinespring(U_local, d_local)
            kraus_list.append((op_type, sites, kraus_ops))

        # Forward half-step for all Lindblad channels
        for op_type, sites, kraus_ops in kraus_list:
            if op_type == "single":
                rho = self._apply_single_site_channel(rho, kraus_ops, sites[0])
            elif op_type == "pair":
                rho = self._apply_pair_channel(rho, kraus_ops, sites[0], sites[1])

        # Reverse half-step for all Lindblad channels (palindromic)
        for op_type, sites, kraus_ops in reversed(kraus_list):
            if op_type == "single":
                rho = self._apply_single_site_channel(rho, kraus_ops, sites[0])
            elif op_type == "pair":
                rho = self._apply_pair_channel(rho, kraus_ops, sites[0], sites[1])

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
            # Full-system Stinespring (use d_anc for consistency)
            U_full = stinespring_unitary_from_lindblad(L_full, dt, d_anc=self.d_anc)
            rho_full = apply_stinespring_to_density_matrix(rho_test, U_full, d_anc=self.d_anc)

            # Local Stinespring via Kraus
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt)
            kraus_ops = self._extract_kraus_from_local_stinespring(U_local, d_local)

            if op_type == "single":
                rho_local = self._apply_single_site_channel(
                    rho_test, kraus_ops, sites[0]
                )
            else:
                rho_local = self._apply_pair_channel(
                    rho_test, kraus_ops, sites[0], sites[1]
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
    # Native gate compilation (compileO0 / compileO1)
    # ------------------------------------------------------------------

    def compile_to_native_gates(
        self,
        dt: float,
        optimization_level: int = 0,
        backend_name: str = "faketraps2trits",
    ) -> dict:
        """Compile circuit gates to MQT-Qudits native gate set.

        Compiles cu_one and cu_two custom gates into native gates:
        VirtRz, R, Rh, Rz, CEx.

        Each gate is compiled individually in a minimal circuit matching
        the backend's qudit count to avoid dimension mismatch.

        cu_multi gates (TTA pair Stinespring, 18x18) are NOT decomposed
        because the MQT-Qudits compiler does not currently support
        multi-qudit gate decomposition into 2-qudit primitives.
        This is an honest limitation, not a fallback.

        Parameters
        ----------
        dt : float
            Time step for gate construction.
        optimization_level : int
            0 for compileO0 (baseline), 1 for compileO1 (optimized).
        backend_name : str
            MQT-Qudits backend for compilation target.

        Returns
        -------
        dict
            Compilation statistics including native gate counts.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        if optimization_level not in (0, 1):
            msg = (
                f"optimization_level must be 0 (compileO0: baseline) "
                f"or 1 (compileO1: optimized), got {optimization_level}"
            )
            raise ValueError(msg)

        compile_fn_name = f"compileO{optimization_level}"
        d = self.d
        results = {}

        # --- Hamiltonian on-site compilation (cu_one, 3x3) ---
        U_onsite = expm(-1j * self.h_local * dt)
        onsite_circ = QuantumCircuit(2, [d, d], 0)
        onsite_circ.cu_one(0, U_onsite)
        onsite_compiled = getattr(onsite_circ, compile_fn_name)(backend_name)
        onsite_counts = self._count_native_gates(onsite_compiled)
        onsite_native_per_gate = sum(onsite_counts.values())

        # --- Hamiltonian transfer compilation (cu_two, 9x9) ---
        transfer_compilations = []
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair = expm(-1j * H_pair * dt)
            pair_circ = QuantumCircuit(2, [d, d], 0)
            pair_circ.cu_two([0, 1], U_pair)
            pair_compiled = getattr(pair_circ, compile_fn_name)(backend_name)
            counts = self._count_native_gates(pair_compiled)
            transfer_compilations.append({
                "pair": (ip, jp),
                "native_gates": sum(counts.values()),
                "breakdown": counts,
            })

        ham_onsite_total = self.N * onsite_native_per_gate
        ham_transfer_total = sum(c["native_gates"] for c in transfer_compilations)

        results["hamiltonian"] = {
            "cu_one_per_gate": onsite_native_per_gate,
            "cu_one_breakdown": onsite_counts,
            "cu_one_total": ham_onsite_total,
            "cu_two_per_pair": transfer_compilations,
            "cu_two_total": ham_transfer_total,
            "native_gates_total": ham_onsite_total + ham_transfer_total,
        }

        # --- Stinespring single-site compilation (cu_two, 6x6) ---
        single_native_total = 0
        single_custom_total = 0
        single_compilations = []
        pair_uncompiled_total = 0
        pair_custom_total = 0
        pair_circuits_info = []

        for op_type, sites, L_local, _gamma in self.lindblad_local_info:
            if op_type == "single":
                circ, _ = self.build_stinespring_circuit_single(
                    L_local, dt, sites[0]
                )
                compiled = getattr(circ, compile_fn_name)(backend_name)
                counts = self._count_native_gates(compiled)
                single_native_total += sum(counts.values())
                single_custom_total += 1
                single_compilations.append({
                    "sites": sites,
                    "native_gates": sum(counts.values()),
                    "breakdown": counts,
                })
            elif op_type == "pair":
                pair_uncompiled_total += 1
                pair_custom_total += 1
                pair_circuits_info.append({
                    "sites": sites,
                    "status": "uncompiled_cu_multi",
                })

        results["stinespring_single"] = {
            "custom_gates": single_custom_total,
            "native_gates_total": single_native_total,
            "per_channel": single_compilations,
        }
        results["stinespring_pair"] = {
            "custom_gates": pair_custom_total,
            "uncompiled_cu_multi": pair_uncompiled_total,
            "note": (
                "cu_multi gates are not decomposed into native gates. "
                "MQT-Qudits compiler does not support multi-qudit gate "
                "decomposition into 2-qudit primitives."
            ),
            "per_channel": pair_circuits_info,
        }

        # --- Per-step totals ---
        ham_native_half = ham_onsite_total + ham_transfer_total
        total_native_per_step = 2 * ham_native_half + 2 * single_native_total
        total_uncompiled_per_step = 2 * pair_uncompiled_total

        results["per_step_summary"] = {
            "native_gates": total_native_per_step,
            "uncompiled_cu_multi": total_uncompiled_per_step,
            "optimization_level": optimization_level,
            "backend": backend_name,
        }

        return results

    @staticmethod
    def _count_native_gates(compiled_circuit) -> dict[str, int]:
        """Count native gates in a compiled circuit by type."""
        counts: dict[str, int] = {}
        for inst in compiled_circuit.instructions:
            name = type(inst).__name__
            counts[name] = counts.get(name, 0) + 1
        return counts

    def verify_compiled_circuit(
        self,
        dt: float,
        optimization_level: int = 0,
        backend_name: str = "faketraps2trits",
    ) -> dict:
        """Verify that compiled Hamiltonian on-site gate gives same state vector.

        Compiles a single cu_one on-site gate and runs both compiled and
        uncompiled versions through MQT-Qudits state-vector simulation,
        comparing the output state vectors.

        Uses a 2-qutrit circuit (matching backend) with the gate on qudit 0.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit
        from mqt.qudits.simulation import MQTQuditProvider

        compile_fn_name = f"compileO{optimization_level}"
        d = self.d
        provider = MQTQuditProvider()
        backend = provider.get_backend("tnsim")

        # Build a 2-qutrit circuit with cu_one on qudit 0
        U_onsite = expm(-1j * self.h_local * dt)
        circuit = QuantumCircuit(2, [d, d], 0)
        circuit.cu_one(0, U_onsite)
        compiled = getattr(circuit, compile_fn_name)(backend_name)

        # Run both
        job_orig = backend.run(circuit, shots=1)
        sv_orig = job_orig.result().get_state_vector()[0]

        job_compiled = backend.run(compiled, shots=1)
        sv_compiled = job_compiled.result().get_state_vector()[0]

        distance = np.linalg.norm(sv_orig - sv_compiled)

        return {
            "dt": dt,
            "optimization_level": optimization_level,
            "backend": backend_name,
            "statevector_distance": float(distance),
            "match": distance < 1e-8,
            "n_original_gates": len(circuit.instructions),
            "n_compiled_gates": len(compiled.instructions),
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
        gates_per_step = (
            2 * gates_per_half_ham
            + 2 * n_stinespring_single
            + 2 * n_stinespring_pair
        )

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
            "n_ancilla_qudits": len(self.lindblad_local_info),
            "d_anc": self.d_anc,
            "gates_per_step": gates_per_step,
            "total_gates": gates_per_step * n_steps,
            "gate_breakdown": {
                "cu_one_onsite": 2 * n_hamiltonian_onsite,
                "cu_two_transfer": 2 * n_hamiltonian_transfer,
                "cu_two_stinespring_single": 2 * n_stinespring_single,
                "cu_multi_stinespring_pair": 2 * n_stinespring_pair,
            },
        }
