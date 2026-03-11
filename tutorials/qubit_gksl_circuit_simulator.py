"""Qubit GKSL simulator using Qiskit QuantumCircuit API.

Scenario 3c: Circuit-based Qubit GKSL-Lindblad (no boson).

This simulator constructs actual Qiskit quantum circuits for each building
block of the Trotter step, encoding each 3-level molecule using 2 qubits
(|S0>->|00>, |T1>->|01>, |S1>->|10>, |11> forbidden).

Circuit decomposition:
  - Hamiltonian step: UnitaryGate (4x4, on-site) + UnitaryGate (16x16, pair transfer)
  - Stinespring channels: UnitaryGate (8x8, single-site) or UnitaryGate (32x32, TTA pair)

Density matrix evolution uses circuit-derived local Kraus operators (in
81-dim qutrit space), which is mathematically equivalent to executing the
full qubit circuit with ancilla reset between channels, but computationally
tractable for open-system dynamics.

Qutrit-to-qubit embedding:
  |0> (S0) -> |00>  (index 0)
  |1> (T1) -> |01>  (index 1)
  |2> (S1) -> |10>  (index 2)
  forbidden  -> |11> (index 3)
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

# Qutrit level -> qubit pair state (2-bit)
_QT_TO_QB = {0: 0b00, 1: 0b01, 2: 0b10}
_QB_TO_QT = {v: k for k, v in _QT_TO_QB.items()}


class QubitGKSLCircuitSimulator:
    """Qubit GKSL simulator using Qiskit QuantumCircuit API.

    Constructs real Qiskit quantum circuits for each Trotter step component.
    Each unitary is embedded from qutrit space into qubit space and applied
    as a ``UnitaryGate``.  The density matrix evolution itself uses the
    same Kraus-operator machinery as the qutrit circuit simulator, operating
    in the 81-dim qutrit space.

    Quantum resources per Trotter step (palindromic 2nd-order):
      - 4 UnitaryGate (4x4 on-site Hamiltonian, per molecule) x 2 = 8
      - 3 UnitaryGate (16x16 pair transfer) x 2 = 6
      - 20 UnitaryGate (8x8 single-site Stinespring) x 2 (fwd+rev) = 40
      - 6 UnitaryGate (32x32 TTA pair Stinespring) x 2 (fwd+rev) = 12
      Total: 66 gates per Trotter step
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QubitGKSLCircuitSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.d = params.d  # 3 (qutrit)
        self.N = params.N_molecules
        self.dim = params.d ** params.N_molecules  # 81

        self._build_local_operators()

    # ------------------------------------------------------------------
    # Qutrit <-> qubit embedding
    # ------------------------------------------------------------------

    @staticmethod
    def _embed_single_site_unitary(U_qt: np.ndarray) -> np.ndarray:
        """Embed a 3x3 qutrit unitary into a 4x4 qubit-pair unitary.

        The forbidden state |11> (index 3) is mapped to itself (identity).
        """
        U_qb = np.eye(4, dtype=np.complex128)
        for i_qt in range(3):
            for j_qt in range(3):
                U_qb[_QT_TO_QB[i_qt], _QT_TO_QB[j_qt]] = U_qt[i_qt, j_qt]
        return U_qb

    @staticmethod
    def _embed_pair_unitary(U_qt_pair: np.ndarray) -> np.ndarray:
        """Embed a 9x9 qutrit-pair unitary into a 16x16 qubit-pair unitary.

        The 7 forbidden indices (any qubit pair in |11>) are mapped to themselves.
        """
        U_qb = np.eye(16, dtype=np.complex128)
        for i_a in range(3):
            for i_b in range(3):
                qt_row = i_a * 3 + i_b
                qb_row = _QT_TO_QB[i_a] * 4 + _QT_TO_QB[i_b]
                for j_a in range(3):
                    for j_b in range(3):
                        qt_col = j_a * 3 + j_b
                        qb_col = _QT_TO_QB[j_a] * 4 + _QT_TO_QB[j_b]
                        U_qb[qb_row, qb_col] = U_qt_pair[qt_row, qt_col]
        return U_qb

    @staticmethod
    def _embed_stinespring_single(U_st_6x6: np.ndarray) -> np.ndarray:
        """Embed a 6x6 Stinespring unitary (ancilla⊗qutrit) into 8x8 (ancilla⊗2-qubit).

        The Stinespring unitary has env(d=2) ⊗ system(d=3) ordering.
        In qubit space: env(d=2) ⊗ qubit-pair(d=4) = 8 states.
        Forbidden states: (anc, |11>) at indices 3 and 7.
        """
        U_qb = np.eye(8, dtype=np.complex128)
        for anc_r in range(2):
            for qt_r in range(3):
                st_row = anc_r * 3 + qt_r
                qb_row = anc_r * 4 + _QT_TO_QB[qt_r]
                for anc_c in range(2):
                    for qt_c in range(3):
                        st_col = anc_c * 3 + qt_c
                        qb_col = anc_c * 4 + _QT_TO_QB[qt_c]
                        U_qb[qb_row, qb_col] = U_st_6x6[st_row, st_col]
        return U_qb

    @staticmethod
    def _embed_stinespring_pair(U_st_18x18: np.ndarray) -> np.ndarray:
        """Embed an 18x18 Stinespring unitary (ancilla⊗qutrit-pair) into 32x32.

        In qubit space: env(d=2) ⊗ 4-qubit(d=16) = 32 states.
        14 forbidden states (any qubit pair in |11>).
        """
        U_qb = np.eye(32, dtype=np.complex128)
        for anc_r in range(2):
            for qa_r in range(3):
                for qb_r in range(3):
                    st_row = anc_r * 9 + qa_r * 3 + qb_r
                    qb_row = anc_r * 16 + _QT_TO_QB[qa_r] * 4 + _QT_TO_QB[qb_r]
                    for anc_c in range(2):
                        for qa_c in range(3):
                            for qb_c in range(3):
                                st_col = anc_c * 9 + qa_c * 3 + qb_c
                                qb_col = (
                                    anc_c * 16
                                    + _QT_TO_QB[qa_c] * 4
                                    + _QT_TO_QB[qb_c]
                                )
                                U_qb[qb_row, qb_col] = U_st_18x18[st_row, st_col]
        return U_qb

    # ------------------------------------------------------------------
    # Local operator construction (in qutrit space for Kraus evolution)
    # ------------------------------------------------------------------

    def _build_local_operators(self) -> None:
        """Pre-compute local operators for circuit gate construction."""
        d = self.d

        # On-site Hamiltonian (3x3 diagonal)
        self.h_local = np.diag(
            np.array([0.0, self.params.E_T, self.params.E_S], dtype=np.complex128)
        )

        # Transfer Hamiltonian per pair (9x9)
        self.h_transfer_pairs: dict[tuple[int, int], np.ndarray] = {}
        for pair in self.params.neighbors:
            i, j = pair
            H_pair = np.zeros((d * d, d * d), dtype=np.complex128)
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

        # TTA: pair operators (6 total: 3 pairs x 2 channels)
        for i, j in p.neighbors:
            gamma = p.gamma_TTA / 2.0
            sq = np.sqrt(gamma)
            # Channel 1: |2>_i<1| ⊗ |0>_j<1|
            l_pair = sq * np.kron(_ket_bra(2, 1), _ket_bra(0, 1))
            info.append(("pair", (i, j), l_pair, gamma))
            # Channel 2: |0>_i<1| ⊗ |2>_j<1|
            l_pair = sq * np.kron(_ket_bra(0, 1), _ket_bra(2, 1))
            info.append(("pair", (i, j), l_pair, gamma))

        # Fluorescence: sqrt(Gamma_fl) |0><2| per molecule
        for i in range(N):
            gamma = p.Gamma_fl
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 2), gamma))

        # Phosphorescence: sqrt(Gamma_ph) |0><1| per molecule
        for i in range(N):
            gamma = p.Gamma_ph
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 1), gamma))

        # Internal conversion: sqrt(k_IC) |0><2| per molecule
        for i in range(N):
            gamma = p.k_IC
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 2), gamma))

        # ISC S1->T1: sqrt(k_ISC_ST) |1><2| per molecule
        for i in range(N):
            gamma = p.k_ISC_ST
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(1, 2), gamma))

        # ISC T1->S0: sqrt(k_ISC_TS) |0><1| per molecule
        for i in range(N):
            gamma = p.k_ISC_TS
            info.append(("single", (i,), np.sqrt(gamma) * _ket_bra(0, 1), gamma))

        expected = 2 * len(p.neighbors) + 5 * N
        if len(info) != expected:
            msg = f"Expected {expected} Lindblad operators, got {len(info)}"
            raise ValueError(msg)
        return info

    # ------------------------------------------------------------------
    # Circuit construction: Hamiltonian (Qiskit)
    # ------------------------------------------------------------------

    def build_hamiltonian_circuit(self, dt: float):
        """Build Qiskit circuit for Hamiltonian half-step.

        Decomposition using 1st-order Trotter for inner splitting:
          exp(-iH_total dt) ~ [x_i exp(-i h_local dt)] . [Pi_{<i,j>} exp(-iH_pair dt)]

        Each gate is a UnitaryGate embedded in qubit space.

        Returns (circuit, gate_count) tuple.
        """
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate

        N = self.N
        n_qubits = 2 * N  # 8 system qubits

        circuit = QuantumCircuit(n_qubits)
        gate_count = 0

        # On-site phases: 4x4 UnitaryGate per molecule (2 qubits)
        U_onsite_qt = expm(-1j * self.h_local * dt)
        U_onsite_qb = self._embed_single_site_unitary(U_onsite_qt)
        onsite_gate = UnitaryGate(U_onsite_qb, label="H_onsite")
        for i in range(N):
            circuit.append(onsite_gate, [2 * i, 2 * i + 1])
            gate_count += 1

        # Transfer unitaries: 16x16 UnitaryGate per pair (4 qubits)
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair_qt = expm(-1j * H_pair * dt)
            U_pair_qb = self._embed_pair_unitary(U_pair_qt)
            pair_gate = UnitaryGate(U_pair_qb, label="H_transfer")
            circuit.append(pair_gate, [2 * ip, 2 * ip + 1, 2 * jp, 2 * jp + 1])
            gate_count += 1

        return circuit, gate_count

    def _compute_hamiltonian_unitary_from_circuit(self, dt: float) -> np.ndarray:
        """Compute the composite Hamiltonian unitary from local gate decomposition.

        Composes the tensor-product on-site unitary with sequential pair unitaries
        to obtain the full-system unitary in the 81-dim qutrit space.
        """
        d = self.d
        N = self.N
        eye = np.eye(d, dtype=np.complex128)

        # On-site: tensor product of single-site unitaries
        U_onsite = expm(-1j * self.h_local * dt)
        U_full = reduce(np.kron, [U_onsite] * N)

        # Transfer: sequential application of pair unitaries
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair = expm(-1j * H_pair * dt)
            U_pair_full = self._embed_pair_unitary_qutrit(U_pair, ip, jp)
            U_full = U_pair_full @ U_full

        return U_full

    def _embed_pair_unitary_qutrit(
        self, U_pair: np.ndarray, site_i: int, site_j: int
    ) -> np.ndarray:
        """Embed a 9x9 pair unitary into the full 81-dim qutrit system."""
        d = self.d
        N = self.N
        dim = self.dim

        U_full = np.zeros((dim, dim), dtype=np.complex128)
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
                    U_full[row, col] = U_pair[pair_row, pair_col]

        return U_full

    # ------------------------------------------------------------------
    # Circuit construction: Stinespring channels (Qiskit)
    # ------------------------------------------------------------------

    def _build_local_stinespring_unitary(
        self, L_local: np.ndarray, dt: float
    ) -> np.ndarray:
        """Build local Stinespring unitary from local Lindblad operator.

        For a d_local-dimensional operator, constructs the (2*d_local x 2*d_local)
        Stinespring unitary: G = [[0, L†], [L, 0]], U = expm(-i * sqrt(dt) * G).
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
        self, L_local_3x3: np.ndarray, dt: float, target_site: int
    ):
        """Build Qiskit circuit for single-site Stinespring channel.

        Creates a circuit with 2 system qubits (target molecule) and 1 ancilla qubit,
        applying the 8x8 embedded Stinespring unitary as a UnitaryGate.

        Returns (circuit, U_local_6x6_qutrit, U_embedded_8x8_qubit) tuple.
        """
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate

        U_local = self._build_local_stinespring_unitary(L_local_3x3, dt)
        U_embedded = self._embed_stinespring_single(U_local)

        # Circuit: [ancilla(1), sys_q0(1), sys_q1(1)] = 3 qubits
        circuit = QuantumCircuit(3)
        gate = UnitaryGate(U_embedded, label="Stinespring_single")
        circuit.append(gate, [0, 1, 2])

        return circuit, U_local, U_embedded

    def build_stinespring_circuit_pair(
        self,
        L_local_9x9: np.ndarray,
        dt: float,
        site_i: int,
        site_j: int,
    ):
        """Build Qiskit circuit for TTA pair Stinespring channel.

        Creates a circuit with 4 system qubits (2 molecules) and 1 ancilla qubit,
        applying the 32x32 embedded Stinespring unitary as a UnitaryGate.

        Returns (circuit, U_local_18x18_qutrit, U_embedded_32x32_qubit) tuple.
        """
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate

        U_local = self._build_local_stinespring_unitary(L_local_9x9, dt)
        U_embedded = self._embed_stinespring_pair(U_local)

        # Circuit: [ancilla(1), mol_i_q0, mol_i_q1, mol_j_q0, mol_j_q1] = 5 qubits
        circuit = QuantumCircuit(5)
        gate = UnitaryGate(U_embedded, label="Stinespring_pair")
        circuit.append(gate, [0, 1, 2, 3, 4])

        return circuit, U_local, U_embedded

    def build_full_trotter_step_circuit(self, dt: float):
        """Build Qiskit circuit for one full 2nd-order Trotter step.

        Structure: H(dt/2) -> D_1...D_n(dt/2) -> D_n...D_1(dt/2) -> H(dt/2)

        Palindromic Lindblad ordering matches QubitGKSLSimulator._trotter_step.
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
                circ, _, _ = self.build_stinespring_circuit_single(
                    L_local, dt / 2, sites[0]
                )
                lindblad_circuits_forward.append((f"stinespring_fwd_single_{idx}", circ))
                total_gates += 1
            elif op_type == "pair":
                circ, _, _ = self.build_stinespring_circuit_pair(
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
                circ, _, _ = self.build_stinespring_circuit_single(
                    L_local, dt / 2, sites[0]
                )
                lindblad_circuits_reverse.append((f"stinespring_rev_single_{idx}", circ))
                total_gates += 1
            elif op_type == "pair":
                circ, _, _ = self.build_stinespring_circuit_pair(
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

    @staticmethod
    def _extract_kraus_from_local_stinespring(
        U_local: np.ndarray, d_local: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """Extract Kraus operators from a local Stinespring unitary.

        U_local has env x system ordering with a 2-level ancilla:
          U_local = [[A, B], [C, D]]  (blocks of size d_local x d_local)

        For ancilla initially in |0>:
          K_0 = A  (ancilla stays |0>)
          K_1 = C  (ancilla flips to |1>)

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
        """Apply Hamiltonian evolution using circuit-decomposed local unitaries."""
        U = self._compute_hamiltonian_unitary_from_circuit(dt)
        return U @ rho @ U.conj().T

    def _apply_single_site_channel(
        self,
        rho: np.ndarray,
        K0: np.ndarray,
        K1: np.ndarray,
        site: int,
    ) -> np.ndarray:
        """Apply single-site Kraus channel using tensor-product structure."""
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
        """Apply pair Kraus channel using tensor-product structure."""
        K0_full = self._embed_pair_operator_qutrit(K0, site_i, site_j)
        K1_full = self._embed_pair_operator_qutrit(K1, site_i, site_j)
        return K0_full @ rho @ K0_full.conj().T + K1_full @ rho @ K1_full.conj().T

    def _embed_pair_operator_qutrit(
        self, op_pair: np.ndarray, site_i: int, site_j: int
    ) -> np.ndarray:
        """Embed a d^2 x d^2 pair operator into the full N-qutrit system."""
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

        Matches the palindromic structure of QubitGKSLSimulator._trotter_step.
        """
        # Half Hamiltonian (circuit-decomposed)
        rho = self._apply_hamiltonian_step_circuit(rho, dt / 2)

        # Precompute half-step Kraus operators for each Lindblad channel
        kraus_list = []
        for op_type, sites, L_local, _gamma in self.lindblad_local_info:
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt / 2)
            K0, K1 = self._extract_kraus_from_local_stinespring(U_local, d_local)
            kraus_list.append((op_type, sites, K0, K1))

        # Forward half-step for all Lindblad channels
        for op_type, sites, K0, K1 in kraus_list:
            if op_type == "single":
                rho = self._apply_single_site_channel(rho, K0, K1, sites[0])
            elif op_type == "pair":
                rho = self._apply_pair_channel(rho, K0, K1, sites[0], sites[1])

        # Reverse half-step for all Lindblad channels (palindromic)
        for op_type, sites, K0, K1 in reversed(kraus_list):
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
        """Prepare initial density matrix in the qutrit Hilbert space."""
        d = self.d
        N = self.N
        dim = d ** N
        psi = np.zeros(dim, dtype=np.complex128)

        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d ** i) for i in range(N))
            psi[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d ** i) for i in range(N))
            psi[index] = 1.0
        else:
            msg = f"Unknown state type: {state_type}"
            raise ValueError(msg)

        return np.outer(psi, psi.conj())

    # ------------------------------------------------------------------
    # Circuit verification
    # ------------------------------------------------------------------

    def verify_embedding_unitarity(self, dt: float) -> dict:
        """Verify that all embedded qubit unitaries are unitary.

        Checks unitarity of embedded on-site, transfer, and Stinespring gates.
        """
        results = {}

        # On-site
        U_onsite_qt = expm(-1j * self.h_local * dt)
        U_onsite_qb = self._embed_single_site_unitary(U_onsite_qt)
        res_onsite = np.linalg.norm(
            U_onsite_qb.conj().T @ U_onsite_qb - np.eye(4), ord="fro"
        )
        results["onsite_unitarity_residual"] = float(res_onsite)

        # Transfer
        max_transfer_res = 0.0
        for pair in self.params.neighbors:
            H_pair = self.h_transfer_pairs[pair]
            U_pair_qt = expm(-1j * H_pair * dt)
            U_pair_qb = self._embed_pair_unitary(U_pair_qt)
            res = np.linalg.norm(
                U_pair_qb.conj().T @ U_pair_qb - np.eye(16), ord="fro"
            )
            max_transfer_res = max(max_transfer_res, res)
        results["transfer_max_unitarity_residual"] = float(max_transfer_res)

        # Stinespring single-site
        max_ss_res = 0.0
        max_sp_res = 0.0
        for op_type, _sites, L_local, _gamma in self.lindblad_local_info:
            U_st = self._build_local_stinespring_unitary(L_local, dt)
            if op_type == "single":
                U_emb = self._embed_stinespring_single(U_st)
                res = np.linalg.norm(
                    U_emb.conj().T @ U_emb - np.eye(8), ord="fro"
                )
                max_ss_res = max(max_ss_res, res)
            elif op_type == "pair":
                U_emb = self._embed_stinespring_pair(U_st)
                res = np.linalg.norm(
                    U_emb.conj().T @ U_emb - np.eye(32), ord="fro"
                )
                max_sp_res = max(max_sp_res, res)

        results["stinespring_single_max_unitarity_residual"] = float(max_ss_res)
        results["stinespring_pair_max_unitarity_residual"] = float(max_sp_res)
        results["all_unitary"] = all(v < 1e-10 for v in results.values())

        return results

    def verify_circuit_via_qiskit(self, dt: float) -> dict:
        """Verify Hamiltonian circuit using Qiskit operator comparison.

        Computes the operator from the Qiskit circuit and compares with the
        Hamiltonian unitary computed directly in qubit space.
        """
        from qiskit.quantum_info import Operator

        circuit, gate_count = self.build_hamiltonian_circuit(dt)

        # Qiskit: get operator from circuit
        op_qiskit = Operator(circuit).data

        # Build expected operator directly in qubit space
        U_qb_expected = self._compute_hamiltonian_unitary_qubit_space(dt)

        distance = np.linalg.norm(op_qiskit - U_qb_expected, ord="fro")

        return {
            "dt": dt,
            "gate_count": gate_count,
            "operator_distance": float(distance),
            "match": distance < 1e-8,
            "n_qubits": 2 * self.N,
        }

    def _compute_hamiltonian_unitary_qubit_space(self, dt: float) -> np.ndarray:
        """Compute the Hamiltonian unitary directly in 256-dim qubit space.

        Matches the Qiskit circuit exactly (including identity action on
        forbidden |11> states).  The computation mirrors the Qiskit circuit:
        tensor-product on-site gates followed by sequential pair gates.
        """
        N = self.N
        dim_qb = 4 ** N

        # On-site: kron of embedded 4x4 gates
        U_onsite_qt = expm(-1j * self.h_local * dt)
        U_onsite_qb = self._embed_single_site_unitary(U_onsite_qt)
        U_full_qb = reduce(np.kron, [U_onsite_qb] * N)

        # Transfer: embed each 16x16 pair unitary into the full qubit space
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair_qt = expm(-1j * H_pair * dt)
            U_pair_qb_16 = self._embed_pair_unitary(U_pair_qt)
            U_pair_full = self._embed_pair_gate_in_full_qubit_space(
                U_pair_qb_16, ip, jp
            )
            U_full_qb = U_pair_full @ U_full_qb

        return U_full_qb

    def _embed_pair_gate_in_full_qubit_space(
        self, U_pair_qb: np.ndarray, mol_i: int, mol_j: int
    ) -> np.ndarray:
        """Embed a 16x16 qubit-pair gate into the full 4^N qubit space.

        U_pair_qb acts on molecule pair (mol_i, mol_j), each using 2 qubits.
        Other molecules are identity.

        In Qiskit ordering: molecule k contributes to 4^k position.
        The pair subspace spans 4^mol_i and 4^mol_j positions.
        """
        N = self.N
        dim_qb = 4 ** N

        U_full = np.zeros((dim_qb, dim_qb), dtype=np.complex128)
        for row in range(dim_qb):
            for col in range(dim_qb):
                # Decompose into per-molecule 4-ary digits
                row_mols = []
                col_mols = []
                r, c = row, col
                for _ in range(N):
                    row_mols.append(r % 4)
                    col_mols.append(c % 4)
                    r //= 4
                    c //= 4
                # row_mols[k] = 4-ary digit for molecule k (LSB-first = mol 0 first)

                # Non-pair molecules must match
                match = True
                for k in range(N):
                    if k != mol_i and k != mol_j:
                        if row_mols[k] != col_mols[k]:
                            match = False
                            break

                if match:
                    # Map to pair indices (mol_i = lower pair dimension, mol_j = upper)
                    pair_row = row_mols[mol_i] * 4 + row_mols[mol_j]
                    pair_col = col_mols[mol_i] * 4 + col_mols[mol_j]
                    U_full[row, col] = U_pair_qb[pair_row, pair_col]

        return U_full

    def verify_stinespring_embedding(self, dt: float) -> dict:
        """Verify embedded Stinespring unitaries preserve channel physics.

        For each Lindblad channel, compares the Kraus map from the qutrit-space
        Stinespring unitary with the map from the qubit-embedded version.
        """
        rho_test = self.prepare_initial_state("edge_triplet")

        max_distance = 0.0
        results = []
        for idx, (op_type, sites, L_local, _gamma) in enumerate(
            self.lindblad_local_info
        ):
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt)
            K0, K1 = self._extract_kraus_from_local_stinespring(U_local, d_local)

            # Qutrit-space Kraus evolution
            if op_type == "single":
                rho_kraus = self._apply_single_site_channel(
                    rho_test, K0, K1, sites[0]
                )
            else:
                rho_kraus = self._apply_pair_channel(
                    rho_test, K0, K1, sites[0], sites[1]
                )

            # Verify embedding is unitary
            if op_type == "single":
                U_emb = self._embed_stinespring_single(U_local)
                emb_size = 8
            else:
                U_emb = self._embed_stinespring_pair(U_local)
                emb_size = 32

            unitarity_res = np.linalg.norm(
                U_emb.conj().T @ U_emb - np.eye(emb_size), ord="fro"
            )

            results.append({
                "index": idx,
                "op_type": op_type,
                "sites": sites,
                "embedding_unitarity_residual": float(unitarity_res),
                "embedding_is_unitary": unitarity_res < 1e-10,
            })
            max_distance = max(max_distance, unitarity_res)

        return {
            "dt": dt,
            "max_unitarity_residual": float(max_distance),
            "all_unitary": max_distance < 1e-10,
            "per_channel": results,
        }

    # ------------------------------------------------------------------
    # Qiskit transpilation to basic gates
    # ------------------------------------------------------------------

    def transpile_to_basic_gates(
        self, dt: float, optimization_level: int = 1
    ) -> dict:
        """Transpile circuit gates to Qiskit's basic gate set.

        Each building-block circuit is transpiled individually using
        Qiskit's transpiler at the specified optimization level.

        Returns transpilation statistics including basic gate counts.
        """
        from qiskit import transpile

        results = {}

        # --- On-site Hamiltonian ---
        U_onsite_qt = expm(-1j * self.h_local * dt)
        U_onsite_qb = self._embed_single_site_unitary(U_onsite_qt)
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate

        onsite_circ = QuantumCircuit(2)
        onsite_circ.append(UnitaryGate(U_onsite_qb), [0, 1])
        onsite_transpiled = transpile(
            onsite_circ, optimization_level=optimization_level
        )
        onsite_counts = dict(onsite_transpiled.count_ops())
        onsite_total = sum(onsite_counts.values())

        results["onsite"] = {
            "basic_gates_per_gate": onsite_total,
            "breakdown": onsite_counts,
            "depth": onsite_transpiled.depth(),
        }

        # --- Transfer Hamiltonian ---
        transfer_results = []
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair_qt = expm(-1j * H_pair * dt)
            U_pair_qb = self._embed_pair_unitary(U_pair_qt)
            pair_circ = QuantumCircuit(4)
            pair_circ.append(UnitaryGate(U_pair_qb), [0, 1, 2, 3])
            pair_transpiled = transpile(
                pair_circ, optimization_level=optimization_level
            )
            pair_counts = dict(pair_transpiled.count_ops())
            transfer_results.append({
                "pair": (ip, jp),
                "basic_gates": sum(pair_counts.values()),
                "breakdown": pair_counts,
                "depth": pair_transpiled.depth(),
            })

        results["transfer"] = transfer_results

        # --- Stinespring ---
        single_results = []
        pair_results = []
        for op_type, sites, L_local, _gamma in self.lindblad_local_info:
            U_st = self._build_local_stinespring_unitary(L_local, dt)
            if op_type == "single":
                U_emb = self._embed_stinespring_single(U_st)
                ss_circ = QuantumCircuit(3)
                ss_circ.append(UnitaryGate(U_emb), [0, 1, 2])
                ss_transpiled = transpile(
                    ss_circ, optimization_level=optimization_level
                )
                ss_counts = dict(ss_transpiled.count_ops())
                single_results.append({
                    "sites": sites,
                    "basic_gates": sum(ss_counts.values()),
                    "breakdown": ss_counts,
                    "depth": ss_transpiled.depth(),
                })
            elif op_type == "pair":
                U_emb = self._embed_stinespring_pair(U_st)
                sp_circ = QuantumCircuit(5)
                sp_circ.append(UnitaryGate(U_emb), [0, 1, 2, 3, 4])
                sp_transpiled = transpile(
                    sp_circ, optimization_level=optimization_level
                )
                sp_counts = dict(sp_transpiled.count_ops())
                pair_results.append({
                    "sites": sites,
                    "basic_gates": sum(sp_counts.values()),
                    "breakdown": sp_counts,
                    "depth": sp_transpiled.depth(),
                })

        results["stinespring_single"] = single_results
        results["stinespring_pair"] = pair_results

        # --- Per-step totals ---
        total_onsite = self.N * onsite_total * 2  # 2 half-steps
        total_transfer = sum(r["basic_gates"] for r in transfer_results) * 2
        total_single = sum(r["basic_gates"] for r in single_results) * 2  # fwd+rev
        total_pair = sum(r["basic_gates"] for r in pair_results) * 2  # fwd+rev
        total_per_step = total_onsite + total_transfer + total_single + total_pair

        results["per_step_summary"] = {
            "basic_gates_total": total_per_step,
            "onsite_total": total_onsite,
            "transfer_total": total_transfer,
            "stinespring_single_total": total_single,
            "stinespring_pair_total": total_pair,
            "optimization_level": optimization_level,
        }

        return results

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run circuit-based Qubit GKSL simulation.

        Uses Qiskit circuit-derived local unitaries for density matrix
        evolution. The Hamiltonian is decomposed into per-molecule and
        per-pair UnitaryGates, and Stinespring channels use local
        embedded UnitaryGates.

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
            "method": "qubit_gksl_circuit",
            "params": self.params.to_dict(),
            "n_system_qubits": 2 * self.N,
            "n_ancilla_qubits": len(self.lindblad_local_info),
            "n_total_qubits": 2 * self.N + len(self.lindblad_local_info),
            "gates_per_step": gates_per_step,
            "total_gates": gates_per_step * n_steps,
            "gate_breakdown": {
                "unitary_4x4_onsite": 2 * n_hamiltonian_onsite,
                "unitary_16x16_transfer": 2 * n_hamiltonian_transfer,
                "unitary_8x8_stinespring_single": 2 * n_stinespring_single,
                "unitary_32x32_stinespring_pair": 2 * n_stinespring_pair,
            },
        }
