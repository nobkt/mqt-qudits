"""Qubit GKSL simulator with boson (phonon) interaction using Stinespring + Trotter.

Scenario 4: Qubit-based GKSL-Lindblad with boson.

Encodes each 3-level molecule using 2 qubits and phonon modes using native
Fock space. Works in the extended qubit-electronic + phonon Hilbert space
dim_total = 4^N * (n_max+1)^N, using the Stinespring + 2nd-order Trotter approach.

The electronic degrees of freedom are embedded in the qubit-encoded space
(4^N instead of d^N), while phonon degrees of freedom remain in the native
Fock space. This faithfully represents what a qubit quantum computer would
compute for the boson-coupled TTA-UC system.
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import (
    build_H_total_boson,
    build_lindblad_operators,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    extend_lindblad_operators,
    partial_trace_phonon,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_simulator import build_qubit_qutrit_mapping
from stinespring_utils import (
    apply_stinespring_to_density_matrix,
    stinespring_unitary_from_lindblad,
)


def _build_extended_qubit_mapping(
    mapping_el: dict[int, int], dim_el_qubit: int, dim_ph: int
) -> dict[int, int]:
    """Build mapping from qutrit-extended index to qubit-extended index.

    For the tensor product structure el ⊗ ph:
      qutrit_ext_idx = el_qt * dim_ph + ph
      qubit_ext_idx  = el_qb * dim_ph + ph
    """
    ext_mapping: dict[int, int] = {}
    for el_qt, el_qb in mapping_el.items():
        for ph in range(dim_ph):
            qt_idx = el_qt * dim_ph + ph
            qb_idx = el_qb * dim_ph + ph
            ext_mapping[qt_idx] = qb_idx
    return ext_mapping


def _embed_extended_operator(
    op: np.ndarray, ext_mapping: dict[int, int], dim_qubit_ext: int
) -> np.ndarray:
    """Embed extended-space operator into qubit-encoded extended space."""
    op_q = np.zeros((dim_qubit_ext, dim_qubit_ext), dtype=np.complex128)
    for i_qt, i_qb in ext_mapping.items():
        for j_qt, j_qb in ext_mapping.items():
            op_q[i_qb, j_qb] = op[i_qt, j_qt]
    return op_q


def _embed_extended_density_matrix(
    rho: np.ndarray, ext_mapping: dict[int, int], dim_qubit_ext: int
) -> np.ndarray:
    """Embed extended-space density matrix into qubit-encoded extended space."""
    rho_q = np.zeros((dim_qubit_ext, dim_qubit_ext), dtype=np.complex128)
    for i_qt, i_qb in ext_mapping.items():
        for j_qt, j_qb in ext_mapping.items():
            rho_q[i_qb, j_qb] = rho[i_qt, j_qt]
    return rho_q


def _partial_trace_phonon_qubit(
    rho_total: np.ndarray, dim_el_qubit: int, dim_ph: int
) -> np.ndarray:
    """Trace out phonon DOFs from qubit-extended density matrix."""
    rho_reshaped = rho_total.reshape((dim_el_qubit, dim_ph, dim_el_qubit, dim_ph))
    return np.trace(rho_reshaped, axis1=1, axis2=3).astype(np.complex128)


class QubitGKSLBosonSimulator:
    """Qubit GKSL simulator with boson (phonon) interaction.

    Works in the extended qubit space dim_total = 4^N * (n_max+1)^N.
    Electronic DOFs are embedded in the 4^N qubit-encoded space.
    Uses Stinespring + 2nd-order Trotter in this extended space.
    After simulation, partial-traces over phonon, then extracts
    electronic observables from the qubit subspace.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            raise ValueError("QubitGKSLBosonSimulator requires with_boson=True")
        self.params = params
        self.N = params.N_molecules

        # Dimensions
        self.dim_el_qutrit = params.d ** params.N_molecules
        self.dim_el_qubit = (2**2) ** params.N_molecules
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_total_qutrit = self.dim_el_qutrit * self.dim_ph
        self.dim_total_qubit = self.dim_el_qubit * self.dim_ph

        # Qubit counts
        self.n_el_qubits = 2 * self.N
        self.n_ph_qubits = int(np.ceil(np.log2(params.n_max + 1))) * self.N
        self.n_ancilla = 2 * len(params.neighbors) + 5 * params.N_molecules
        self.n_total_qubits = self.n_el_qubits + self.n_ph_qubits + self.n_ancilla

        # Build qubit-qutrit mappings
        self._mapping_el = build_qubit_qutrit_mapping(self.N, params.d)
        self._mapping_ext = _build_extended_qubit_mapping(
            self._mapping_el, self.dim_el_qubit, self.dim_ph
        )

        # Build qutrit-space extended Hamiltonian, then embed
        H_total_qt = build_H_total_boson(params)
        self.H_total = _embed_extended_operator(
            H_total_qt, self._mapping_ext, self.dim_total_qubit
        )

        # Build extended Lindblad operators, then embed
        lindblad_ops_el = build_lindblad_operators(params)
        lindblad_ops_ext_qt = extend_lindblad_operators(lindblad_ops_el, self.dim_ph)
        self.lindblad_ops = [
            (
                _embed_extended_operator(L, self._mapping_ext, self.dim_total_qubit),
                gamma,
            )
            for L, gamma in lindblad_ops_ext_qt
        ]

    # ------------------------------------------------------------------
    # Trotter step primitives (in qubit-extended space)
    # ------------------------------------------------------------------

    def _precompute_unitaries(self, dt: float) -> None:
        """Pre-compute time-step-dependent unitaries (called once per simulation)."""
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines = [
            stinespring_unitary_from_lindblad(L_op, dt)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _trotter_step(self, rho: np.ndarray) -> np.ndarray:
        """2nd-order symmetric Trotter step in qubit-extended space."""
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        for U_stine in self._U_stines:
            rho = apply_stinespring_to_density_matrix(rho, U_stine)
        rho = self._U_H_half @ rho @ self._U_H_half.conj().T
        return rho

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Prepare initial density matrix: electronic ⊗ phonon vacuum, in qubit space."""
        d = self.params.d
        N = self.N
        dim_el_qt = d**N

        # Electronic initial state (qutrit space)
        psi_el = np.zeros(dim_el_qt, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            index = 1 * (d ** (N - 1)) + 1
            psi_el[index] = 1.0
        elif state_type == "all_triplet":
            index = sum(1 * (d**i) for i in range(N))
            psi_el[index] = 1.0
        elif state_type == "all_singlet":
            index = sum(2 * (d**i) for i in range(N))
            psi_el[index] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

        # Phonon vacuum
        psi_ph = np.zeros(self.dim_ph, dtype=np.complex128)
        psi_ph[0] = 1.0

        # Full state in qutrit-extended space
        psi_full_qt = np.kron(psi_el, psi_ph)
        rho_qt = np.outer(psi_full_qt, psi_full_qt.conj())

        # Embed in qubit-extended space
        return _embed_extended_density_matrix(
            rho_qt, self._mapping_ext, self.dim_total_qubit
        )

    # ------------------------------------------------------------------
    # Electronic observable extraction
    # ------------------------------------------------------------------

    def _extract_electronic_rho(self, rho_qubit_ext: np.ndarray) -> np.ndarray:
        """Extract electronic-only reduced density matrix in qutrit space.

        1. Partial trace over phonon DOFs (in qubit-extended space)
        2. Extract qutrit subspace from qubit-electronic space
        """
        # Partial trace over phonon
        rho_el_qubit = _partial_trace_phonon_qubit(
            rho_qubit_ext, self.dim_el_qubit, self.dim_ph
        )
        # Extract qutrit subspace
        from qubit_gksl_simulator import extract_density_matrix_from_qubit_space

        return extract_density_matrix_from_qubit_space(
            rho_el_qubit, self._mapping_el, self.dim_el_qutrit
        )

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run simulation in qubit-encoded extended space.

        Returns a dict with time series of populations, entropy, purity,
        trace, and qubit circuit statistics. Populations are computed from
        the qutrit electronic subspace after partial-tracing over phonons.
        """
        start = time_module.time()

        dt = t_max / n_steps
        self._precompute_unitaries(dt)
        rho = self.prepare_initial_state(initial_state)

        # Extract electronic observables
        rho_el = self._extract_electronic_rho(rho)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_el, self.params)]
        entropies = [compute_von_neumann_entropy(rho_el)]
        purities = [compute_purity(rho_el)]
        traces = [float(np.real(np.trace(rho)))]

        for step in range(n_steps):
            rho = self._trotter_step(rho)

            rho_el = self._extract_electronic_rho(rho)

            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho))))
            populations.append(
                compute_populations_from_density_matrix(rho_el, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el))
            purities.append(compute_purity(rho_el))

        elapsed = time_module.time() - start

        # Return electronic-only qutrit density matrix for cross-comparison
        rho_el_final = self._extract_electronic_rho(rho)

        # Gate count estimate for qubit circuit in extended space
        gates_per_step = (
            self.n_el_qubits
            + self.n_ph_qubits
            + len(self.params.neighbors) * 10
            + self.n_ancilla * 6
        )

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_el_final,
            "elapsed_time": elapsed,
            "method": "qubit_gksl_boson",
            "params": self.params.to_dict(),
            "n_el_qubits": self.n_el_qubits,
            "n_ph_qubits": self.n_ph_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "dim_total": self.dim_total_qubit,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
        }
