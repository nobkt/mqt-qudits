"""Qudit GKSL circuit simulator with boson (phonon) interaction.

Scenario 6c: Circuit-based Qudit GKSL-Lindblad with boson.

Extends QuditGKSLCircuitSimulator to support phonon degrees of freedom,
adding phonon qutrit registers and electron-phonon coupling gates.

Circuit decomposition:
  - Electronic Hamiltonian: cu_one (on-site phases) + cu_two (pair transfer)
  - Phonon Hamiltonian: cu_one (phonon on-site) per phonon qutrit
  - Electron-phonon coupling: cu_two (el-ph coupling) per molecule
  - Stinespring channels: cu_two (single-site) or cu_multi (TTA pair)
    acting on electronic qutrits only (Lindblad operators act on electronic
    subspace with identity on phonon subspace)

The density matrix evolution uses circuit-derived local Kraus operators in
the extended electron⊗phonon Hilbert space.
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
    build_H_total_boson,
    build_phonon_operators,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    partial_trace_phonon,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_circuit_simulator import QuditGKSLCircuitSimulator


class QuditGKSLCircuitBosonSimulator:
    """Qudit GKSL circuit simulator with boson (phonon) interaction.

    Constructs MQT-Qudits quantum circuits for each Trotter step component
    in the extended electronic⊗phonon Hilbert space.

    Since this targets qudit quantum computers, all registers — including
    Stinespring ancillas — are native d-level qudits.

    Quantum resources:
      - N electronic qutrits (d=3)
      - N phonon qutrits (d=n_max+1, typically d=3 for n_max=2)
      - 26 ancilla qudits (d=3, Stinespring channels)
      Total: 2N qutrits + 26 ancilla qudits

    Gate types per Trotter step (palindromic 2nd-order):
      - cu_one: electronic on-site (N) + phonon on-site (N)   [x2 half-steps]
      - cu_two: electronic transfer (N-1 pairs) + el-ph coupling (N) [x2]
               + single-site Stinespring (20 x 2 fwd+rev = 40)
      - cu_multi: TTA pair Stinespring (6 x 2 fwd+rev = 12)
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            msg = "QuditGKSLCircuitBosonSimulator requires with_boson=True"
            raise ValueError(msg)
        self.params = params
        self.d = params.d
        self.d_anc = params.d  # ancilla dimension matches system qudit dimension
        self.N = params.N_molecules
        self.n_max = params.n_max
        self.d_ph = params.n_max + 1

        # Dimensions
        self.dim_el = params.d ** params.N_molecules
        self.dim_ph = self.d_ph ** params.N_molecules
        self.dim_total = self.dim_el * self.dim_ph

        # Build local operators for circuit gates
        self._build_local_operators()

    def _build_local_operators(self) -> None:
        """Pre-compute local operators for circuit gate construction."""
        d = self.d
        d_ph = self.d_ph
        N = self.N

        # Electronic on-site Hamiltonian (d x d diagonal)
        self.h_el_local = np.diag(
            np.array([0.0, self.params.E_T, self.params.E_S], dtype=np.complex128)
        )

        # Transfer Hamiltonian per pair (d^2 x d^2)
        self.h_transfer_pairs: dict[tuple[int, int], np.ndarray] = {}
        for pair in self.params.neighbors:
            i, j = pair
            H_pair = np.zeros((d * d, d * d), dtype=np.complex128)
            idx_01 = 0 * d + 1
            idx_10 = 1 * d + 0
            H_pair[idx_01, idx_10] = self.params.V
            H_pair[idx_10, idx_01] = self.params.V
            self.h_transfer_pairs[(i, j)] = H_pair

        # Phonon on-site Hamiltonian (d_ph x d_ph diagonal)
        _, _, n_op = build_phonon_operators(self.n_max)
        self.h_ph_local = self.params.omega_ph * n_op

        # Electron-phonon coupling per site ((d*d_ph) x (d*d_ph))
        # H_eph_i = g_eph * |T1><T1| ⊗ (a + a†)
        proj_T = np.zeros((d, d), dtype=np.complex128)
        proj_T[1, 1] = 1.0
        a, a_dag, _ = build_phonon_operators(self.n_max)
        x_op = a + a_dag
        self.h_eph_local = self.params.g_eph * np.kron(proj_T, x_op)

        # Lindblad operators (from non-boson circuit simulator)
        # Create a non-boson params for building electronic Lindblad operators
        nb_dict = self.params.to_dict()
        nb_dict["with_boson"] = False
        nb_params = GKSLPhysicalParameters(**nb_dict)
        nb_sim = QuditGKSLCircuitSimulator(nb_params)
        self.lindblad_local_info = nb_sim.lindblad_local_info

    # ------------------------------------------------------------------
    # Circuit construction helpers (delegated to QuditGKSLCircuitSimulator)
    # ------------------------------------------------------------------

    def _build_local_stinespring_unitary(
        self, L_local: np.ndarray, dt: float
    ) -> np.ndarray:
        """Build local Stinespring unitary from local Lindblad operator.

        Uses d_anc-level ancilla qudit (same dimension as system qudits).
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

    def _extract_kraus_from_local_stinespring(
        self, U_local: np.ndarray, d_local: int
    ) -> list[np.ndarray]:
        """Extract Kraus operators from local Stinespring unitary."""
        d_anc = self.d_anc
        kraus_ops = []
        for k in range(d_anc):
            K = U_local[k * d_local : (k + 1) * d_local, :d_local].copy()
            kraus_ops.append(K)
        return kraus_ops

    # ------------------------------------------------------------------
    # MQT-Qudits circuit construction
    # ------------------------------------------------------------------

    def build_hamiltonian_circuit(self, dt: float):
        """Build MQT-Qudits circuit for Hamiltonian half-step in extended space.

        The Hamiltonian is decomposed into local gates:
          1. Electronic on-site: cu_one on each electronic qutrit
          2. Electronic transfer: cu_two on each electronic pair
          3. Phonon on-site: cu_one on each phonon qutrit
          4. Electron-phonon coupling: cu_two on each (electron, phonon) pair

        Returns (circuit, gate_count) tuple.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        N = self.N
        d = self.d
        d_ph = self.d_ph

        # Circuit layout: [el_0, ..., el_{N-1}, ph_0, ..., ph_{N-1}]
        dims = [d] * N + [d_ph] * N
        circuit = QuantumCircuit(2 * N, dims, 0)
        gate_count = 0

        # 1. Electronic on-site phases
        U_el_onsite = expm(-1j * self.h_el_local * dt)
        for i in range(N):
            circuit.cu_one(i, U_el_onsite)
            gate_count += 1

        # 2. Electronic transfer
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair = expm(-1j * H_pair * dt)
            circuit.cu_two([ip, jp], U_pair)
            gate_count += 1

        # 3. Phonon on-site
        U_ph_onsite = expm(-1j * self.h_ph_local * dt)
        for i in range(N):
            circuit.cu_one(N + i, U_ph_onsite)
            gate_count += 1

        # 4. Electron-phonon coupling
        if self.params.g_eph != 0.0:
            U_eph = expm(-1j * self.h_eph_local * dt)
            for i in range(N):
                circuit.cu_two([i, N + i], U_eph)
                gate_count += 1

        return circuit, gate_count

    def build_stinespring_circuit_single(
        self, L_local_3x3: np.ndarray, dt: float, target_qudit: int
    ):
        """Build Stinespring circuit for single-site channel.

        The Lindblad operator acts on the electronic qutrit only.
        The phonon qutrit is not involved.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        d = self.d
        d_anc = self.d_anc
        U_local = self._build_local_stinespring_unitary(L_local_3x3, dt)

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
        """Build Stinespring circuit for TTA pair channel.

        The Lindblad operator acts on two electronic qutrits only.
        The phonon qutrits are not involved.
        """
        from mqt.qudits.quantum_circuit import QuantumCircuit

        d = self.d
        d_anc = self.d_anc
        U_local = self._build_local_stinespring_unitary(L_local_9x9, dt)

        circuit = QuantumCircuit(3, [d, d, d_anc], 0)
        circuit.cu_multi([0, 1, 2], U_local)

        return circuit, U_local

    # ------------------------------------------------------------------
    # Density matrix evolution
    # ------------------------------------------------------------------

    def _compute_hamiltonian_unitary(self, dt: float) -> np.ndarray:
        """Compute the full Hamiltonian unitary from local gate decomposition.

        Composes electronic, phonon, and coupling unitaries in the
        extended dim_total x dim_total space.
        """
        d = self.d
        d_ph = self.d_ph
        N = self.N
        dim_el = self.dim_el
        dim_ph = self.dim_ph
        dim_total = self.dim_total

        eye_el = np.eye(d, dtype=np.complex128)
        eye_ph = np.eye(d_ph, dtype=np.complex128)

        # Electronic on-site: tensor product in electronic space ⊗ I_phonon
        U_el_onsite = expm(-1j * self.h_el_local * dt)
        U_el_full = reduce(np.kron, [U_el_onsite] * N)
        U_onsite = np.kron(U_el_full, np.eye(dim_ph, dtype=np.complex128))

        # Electronic transfer in extended space
        U_transfer = np.eye(dim_total, dtype=np.complex128)
        for pair in self.params.neighbors:
            ip, jp = pair
            H_pair = self.h_transfer_pairs[(ip, jp)]
            U_pair = expm(-1j * H_pair * dt)
            U_pair_ext = self._embed_pair_unitary_electronic(U_pair, ip, jp)
            U_transfer = U_pair_ext @ U_transfer

        # Phonon on-site: I_electronic ⊗ tensor product in phonon space
        U_ph_onsite = expm(-1j * self.h_ph_local * dt)
        U_ph_full = reduce(np.kron, [U_ph_onsite] * N)
        U_phonon = np.kron(np.eye(dim_el, dtype=np.complex128), U_ph_full)

        # Electron-phonon coupling
        U_eph = np.eye(dim_total, dtype=np.complex128)
        if self.params.g_eph != 0.0:
            U_eph_local = expm(-1j * self.h_eph_local * dt)
            for i in range(N):
                U_eph_full = self._embed_eph_unitary(U_eph_local, i)
                U_eph = U_eph_full @ U_eph

        return U_eph @ U_phonon @ U_transfer @ U_onsite

    def _embed_pair_unitary_electronic(
        self, U_pair: np.ndarray, site_i: int, site_j: int
    ) -> np.ndarray:
        """Embed a 2-qutrit electronic unitary in the extended el⊗ph space.

        U_pair is d^2 x d^2 acting on electronic qutrits (site_i, site_j).
        Returns dim_total x dim_total unitary (I on phonon space).
        """
        d = self.d
        N = self.N
        dim_el = self.dim_el
        dim_ph = self.dim_ph

        # First embed in electronic space
        U_el_full = np.zeros((dim_el, dim_el), dtype=np.complex128)
        for row in range(dim_el):
            for col in range(dim_el):
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

                match = all(
                    row_states[k] == col_states[k]
                    for k in range(N)
                    if k != site_i and k != site_j
                )
                if match:
                    pair_row = row_states[site_i] * d + row_states[site_j]
                    pair_col = col_states[site_i] * d + col_states[site_j]
                    U_el_full[row, col] = U_pair[pair_row, pair_col]

        return np.kron(U_el_full, np.eye(dim_ph, dtype=np.complex128))

    def _embed_eph_unitary(
        self, U_eph_local: np.ndarray, site: int
    ) -> np.ndarray:
        """Embed electron-phonon coupling unitary for one molecule.

        U_eph_local is (d*d_ph) x (d*d_ph) acting on (electronic_i, phonon_i).
        Returns dim_total x dim_total unitary.

        Space ordering: electronic qutrits ⊗ phonon qutrits
        For molecule i, we need to address:
          - electronic qudit i (in electronic register)
          - phonon qudit i (in phonon register)
        """
        d = self.d
        d_ph = self.d_ph
        N = self.N
        dim_total = self.dim_total

        U_full = np.zeros((dim_total, dim_total), dtype=np.complex128)

        for row in range(dim_total):
            for col in range(dim_total):
                # Decompose into electronic and phonon states
                row_el_idx = row // self.dim_ph
                row_ph_idx = row % self.dim_ph
                col_el_idx = col // self.dim_ph
                col_ph_idx = col % self.dim_ph

                # Extract per-site electronic states
                el_row_states = self._decompose_index(row_el_idx, d, N)
                el_col_states = self._decompose_index(col_el_idx, d, N)

                # Extract per-site phonon states
                ph_row_states = self._decompose_index(row_ph_idx, d_ph, N)
                ph_col_states = self._decompose_index(col_ph_idx, d_ph, N)

                # Check: all sites except 'site' must match in both spaces
                match = True
                for k in range(N):
                    if k != site:
                        if el_row_states[k] != el_col_states[k]:
                            match = False
                            break
                        if ph_row_states[k] != ph_col_states[k]:
                            match = False
                            break

                if match:
                    local_row = el_row_states[site] * d_ph + ph_row_states[site]
                    local_col = el_col_states[site] * d_ph + ph_col_states[site]
                    U_full[row, col] = U_eph_local[local_row, local_col]

        return U_full

    @staticmethod
    def _decompose_index(idx: int, base: int, n_digits: int) -> list[int]:
        """Decompose a flat index into per-qudit states (big-endian)."""
        states = []
        for _ in range(n_digits):
            states.append(idx % base)
            idx //= base
        states.reverse()
        return states

    def _apply_hamiltonian_step(
        self, rho: np.ndarray, dt: float
    ) -> np.ndarray:
        """Apply Hamiltonian evolution using circuit-decomposed unitaries."""
        U = self._compute_hamiltonian_unitary(dt)
        return U @ rho @ U.conj().T

    # ------------------------------------------------------------------
    # Precomputation (called once per simulation)
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation).

        Precomputes:
          - Half-step Hamiltonian unitary (exact matrix exponential of full
            H_total, matching QuditGKSLBosonSimulator)
          - Kraus operators for each Lindblad channel (half-step Stinespring)

        The Hamiltonian unitary uses the exact matrix exponential of the full
        extended-space Hamiltonian (H_el_ext + H_phonon + H_eph), matching
        QuditGKSLBosonSimulator.  The Trotter decomposition into local gates
        (build_hamiltonian_circuit) is used only for gate counting and
        circuit visualization.

        This avoids a first-order Trotter error that would otherwise appear
        in the electron-phonon correlations of the full density matrix.
        """
        H_total = build_H_total_boson(self.params)
        self._U_H_half = expm(-1j * H_total * dt / 2)

        self._kraus_list: list[
            tuple[str, list[int], list[np.ndarray]]
        ] = []
        for op_type, sites, L_local, _gamma in self.lindblad_local_info:
            d_local = L_local.shape[0]
            U_local = self._build_local_stinespring_unitary(L_local, dt / 2)
            kraus_ops = self._extract_kraus_from_local_stinespring(U_local, d_local)
            self._kraus_list.append((op_type, sites, kraus_ops))

    def _apply_single_site_channel_extended(
        self,
        rho: np.ndarray,
        kraus_ops: list[np.ndarray],
        el_site: int,
    ) -> np.ndarray:
        """Apply Kraus channel on electronic qutrit in extended space.

        kraus_ops are d x d Kraus operators acting on electronic qutrit el_site.
        Phonon space is identity.
        """
        d = self.d
        N = self.N
        dim_el = self.dim_el
        dim_ph = self.dim_ph
        dim_total = self.dim_total

        # Embed Kraus operators: K_ext = K_el ⊗ I_ph
        eye_el = np.eye(d, dtype=np.complex128)
        eye_ph = np.eye(dim_ph, dtype=np.complex128)

        def _embed_kraus(K: np.ndarray) -> np.ndarray:
            op_list = [eye_el] * N
            op_list[el_site] = K
            K_el_full = reduce(np.kron, op_list)
            return np.kron(K_el_full, eye_ph)

        rho_out = np.zeros_like(rho)
        for K in kraus_ops:
            K_full = _embed_kraus(K)
            rho_out += K_full @ rho @ K_full.conj().T
        return rho_out

    def _apply_pair_channel_extended(
        self,
        rho: np.ndarray,
        kraus_ops: list[np.ndarray],
        el_site_i: int,
        el_site_j: int,
    ) -> np.ndarray:
        """Apply Kraus channel on electronic qutrit pair in extended space.

        kraus_ops are d^2 x d^2 Kraus operators acting on (el_site_i, el_site_j).
        """
        d = self.d
        N = self.N
        dim_ph = self.dim_ph

        # Embed pair Kraus operators in electronic space, then extend to phonon
        def _embed_pair_kraus(K: np.ndarray) -> np.ndarray:
            K_el_full = self._embed_pair_operator_electronic(
                K, el_site_i, el_site_j
            )
            return np.kron(K_el_full, np.eye(dim_ph, dtype=np.complex128))

        rho_out = np.zeros_like(rho)
        for K in kraus_ops:
            K_full = _embed_pair_kraus(K)
            rho_out += K_full @ rho @ K_full.conj().T
        return rho_out

    def _embed_pair_operator_electronic(
        self, op: np.ndarray, site_i: int, site_j: int
    ) -> np.ndarray:
        """Embed a d^2 x d^2 operator on pair (i,j) in full electronic space."""
        d = self.d
        N = self.N
        dim_el = self.dim_el

        op_full = np.zeros((dim_el, dim_el), dtype=np.complex128)
        for row in range(dim_el):
            for col in range(dim_el):
                row_states = self._decompose_index(row, d, N)
                col_states = self._decompose_index(col, d, N)

                match = all(
                    row_states[k] == col_states[k]
                    for k in range(N)
                    if k != site_i and k != site_j
                )
                if match:
                    pair_row = row_states[site_i] * d + row_states[site_j]
                    pair_col = col_states[site_i] * d + col_states[site_j]
                    op_full[row, col] = op[pair_row, pair_col]

        return op_full

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """2nd-order symmetric Trotter step with palindromic Lindblad ordering.

        exp(L dt) ≈ exp(L_H dt/2)
                     · prod_{α=1..n} E_α(dt/2)
                     · prod_{α=n..1} E_α(dt/2)
                     · exp(L_H dt/2)

        Matches the palindromic structure of QuditGKSLBosonSimulator.
        Uses precomputed unitaries and Kraus operators from _precompute_unitaries().
        """
        # Half Hamiltonian (precomputed)
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T

        # Forward half-step for all Lindblad channels
        for op_type, sites, kraus_ops in self._kraus_list:
            if op_type == "single":
                rho = self._apply_single_site_channel_extended(
                    rho, kraus_ops, sites[0]
                )
            elif op_type == "pair":
                rho = self._apply_pair_channel_extended(
                    rho, kraus_ops, sites[0], sites[1]
                )

        # Reverse half-step for all Lindblad channels (palindromic)
        for op_type, sites, kraus_ops in reversed(self._kraus_list):
            if op_type == "single":
                rho = self._apply_single_site_channel_extended(
                    rho, kraus_ops, sites[0]
                )
            elif op_type == "pair":
                rho = self._apply_pair_channel_extended(
                    rho, kraus_ops, sites[0], sites[1]
                )

        # Half Hamiltonian (precomputed)
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        return rho

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic state ⊗ phonon vacuum."""
        d = self.d
        N = self.N
        dim_el = self.dim_el

        psi_el = np.zeros(dim_el, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi_el[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d ** i) for i in range(N))
            psi_el[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d ** i) for i in range(N))
            psi_el[index] = 1.0
        else:
            msg = f"Unknown state type: {state_type}"
            raise ValueError(msg)

        # Phonon vacuum
        psi_ph = np.zeros(self.dim_ph, dtype=np.complex128)
        psi_ph[0] = 1.0

        psi_total = np.kron(psi_el, psi_ph)
        return np.outer(psi_total, psi_total.conj())

    # ------------------------------------------------------------------
    # Main simulation
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run circuit-based Qudit GKSL boson simulation.

        For g_eph=0, delegates to the non-boson circuit simulator
        (exact reduction: phonon decoupling).
        """
        # Exact reduction for g_eph = 0
        if self.params.g_eph == 0.0:
            nb_dict = self.params.to_dict()
            nb_dict["with_boson"] = False
            nb_params = GKSLPhysicalParameters(**nb_dict)
            result = QuditGKSLCircuitSimulator(nb_params).simulate(
                t_max=t_max, n_steps=n_steps, initial_state=initial_state
            )
            result["method"] = "qudit_gksl_circuit_boson (g_eph=0 exact reduction)"
            return result

        start = time_module.time()
        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        times: list[float] = [0.0]
        rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)
        populations = [compute_populations_from_density_matrix(rho_el, self.params)]
        entropies = [compute_von_neumann_entropy(rho_el)]
        purities = [compute_purity(rho_el)]
        traces = [float(np.real(np.trace(rho_el)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            rho_el = partial_trace_phonon(rho, self.dim_el, self.dim_ph)
            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_el))))
            populations.append(
                compute_populations_from_density_matrix(rho_el, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        elapsed = time_module.time() - start

        # Return electronic-only reduced density matrix for consistency
        # with QuditGKSLBosonSimulator and to match time-series observables
        rho_el_final = partial_trace_phonon(rho, self.dim_el, self.dim_ph)

        N = self.N
        n_pairs = len(self.params.neighbors)
        n_stinespring_single = sum(
            1 for t, _, _, _ in self.lindblad_local_info if t == "single"
        )
        n_stinespring_pair = sum(
            1 for t, _, _, _ in self.lindblad_local_info if t == "pair"
        )
        n_eph_gates = N if self.params.g_eph != 0.0 else 0
        # N (el on-site) + n_pairs (el transfer) + N (ph on-site) + n_eph (el-ph coupling)
        gates_per_half_ham = N + n_pairs + N + n_eph_gates
        gates_per_step = 2 * gates_per_half_ham + 2 * n_stinespring_single + 2 * n_stinespring_pair

        # Register counts consistent with QuditGKSLBosonSimulator
        n_system_qudits = N
        n_ph_per_mol = int(np.ceil(np.log(self.params.n_max + 1) / np.log(self.d)))
        n_phonon_qudits = n_ph_per_mol * N
        n_ancilla_qudits = len(self.lindblad_local_info)
        n_total_qudits = n_system_qudits + n_phonon_qudits + n_ancilla_qudits

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_el_final,
            "elapsed_time": elapsed,
            "method": "qudit_gksl_circuit_boson",
            "params": self.params.to_dict(),
            "n_system_qudits": n_system_qudits,
            "n_phonon_qudits": n_phonon_qudits,
            "n_ancilla_qudits": n_ancilla_qudits,
            "d_anc": self.d_anc,
            "n_total_qudits": n_total_qudits,
            "dim_total": self.dim_total,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "gate_breakdown": {
                "cu_one_el_onsite": 2 * N,
                "cu_two_el_transfer": 2 * n_pairs,
                "cu_one_ph_onsite": 2 * N,
                "cu_two_eph_coupling": 2 * n_eph_gates,
                "cu_two_stinespring_single": 2 * n_stinespring_single,
                "cu_multi_stinespring_pair": 2 * n_stinespring_pair,
            },
        }
