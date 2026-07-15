"""Deterministic statevector-ensemble GKSL simulators with explicit ancilla updates.

These simulators propagate the open-system GKSL dynamics **without ever
propagating a density matrix**.  The mixed state is represented as a weighted
ensemble of pure statevectors {(w_i, |psi_i>)}.  For every Lindblad channel a
fresh ancilla |0> is attached, the Stinespring dilation unitary is applied,
and the ancilla is *updated* by branching the statevector over all ancilla
computational-basis outcomes (deterministic branching instead of the Born
sampling used by the shot-based simulators).  Discarding the ancilla after
branching is mathematically the exact partial trace, so in exact arithmetic

    sum_i w_i |psi_i><psi_i|  ==  rho(t)  (density-matrix result, same Trotter)

To keep the ensemble size bounded, the ensemble is compressed after each
channel application using the Gram-matrix method: the eigen-decomposition of
the K x K Gram matrix G_ij = sqrt(w_i w_j) <psi_i|psi_j> yields an equivalent
ensemble of at most rank(rho) orthogonal pure states with the same mixture.
The Gram matrix has the same non-zero spectrum as rho, so entropy/purity are
also computed from it — a dim x dim density matrix is never formed during
propagation (rho_final is reconstructed once at the end, for comparison with
the density-matrix simulators only).

Honest notes:
  - This is exact (up to floating-point round-off) relative to the
    Stinespring + palindromic 2nd-order Trotter scheme of
    QuditGKSLSimulator / QubitGKSLSimulator; the O(dt) Trotter error is the
    same as for those simulators.
  - Branches whose probability is exactly below ``branch_tol`` (default
    1e-14) are dropped; the total dropped weight is tracked and reported in
    the result dict as ``discarded_weight`` so trace accounting stays honest.

Two classes are provided:
  - QuditGKSLStatevectorSimulator: native qutrit registers (d=3 ancillas)
  - QubitGKSLStatevectorSimulator: qubit-pair encoding (2-level ancillas)
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import eigh, expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_shot_simulator import _populations_from_diagonal_qubit
from qubit_gksl_simulator import (
    build_qubit_qutrit_mapping,
    embed_operator_in_qubit_space,
    embed_statevector_in_qubit_space,
    extract_density_matrix_from_qubit_space,
)
from qudit_gksl_shot_simulator import _populations_from_diagonal
from stinespring_utils import stinespring_unitary_from_lindblad

_DEFAULT_BRANCH_TOL = 1e-14
_EIG_KEEP_TOL = 1e-14


class _StatevectorEnsemble:
    """Weighted ensemble of pure statevectors representing a mixed state."""

    def __init__(self, weights: list[float], states: list[np.ndarray]) -> None:
        self.weights = weights
        self.states = states

    @property
    def size(self) -> int:
        return len(self.weights)

    def total_weight(self) -> float:
        return float(sum(self.weights))

    def gram_eigenvalues(self) -> np.ndarray:
        """Eigenvalues of the Gram matrix == non-zero eigenvalues of rho."""
        K = self.size
        sq = np.sqrt(np.asarray(self.weights))
        V = np.column_stack(self.states)  # dim x K
        G = (V.conj().T @ V) * np.outer(sq, sq)
        vals = eigh(G, eigvals_only=True)
        return np.clip(vals, 0.0, None)

    def compress(self, max_rank: int | None = None) -> None:
        """Replace the ensemble by an equivalent orthogonal one via the Gram matrix.

        Exact (round-off level): preserves sum_i w_i |psi_i><psi_i| and the
        total weight.  Never forms the dim x dim density matrix.
        """
        K = self.size
        if K <= 1:
            return
        sq = np.sqrt(np.asarray(self.weights))
        V = np.column_stack(self.states)  # dim x K
        G = (V.conj().T @ V) * np.outer(sq, sq)  # K x K Gram matrix
        vals, vecs = eigh(G)
        total_before = self.total_weight()

        order = np.argsort(vals)[::-1]
        vals = vals[order]
        vecs = vecs[:, order]

        keep = vals > _EIG_KEEP_TOL * max(vals[0], 1.0) if vals.size else vals > 0
        if max_rank is not None:
            keep_idx = np.nonzero(keep)[0][:max_rank]
        else:
            keep_idx = np.nonzero(keep)[0]

        new_weights: list[float] = []
        new_states: list[np.ndarray] = []
        for k in keep_idx:
            lam = float(vals[k])
            if lam <= 0.0:
                continue
            # |phi_k> = (1/sqrt(lam)) * sum_i sqrt(w_i) vecs[i,k] |psi_i>
            phi = V @ (sq * vecs[:, k]) / np.sqrt(lam)
            new_weights.append(lam)
            new_states.append(phi)

        # Renormalize total weight to be exactly preserved (round-off only)
        total_after = float(sum(new_weights))
        if total_after > 0.0 and total_before > 0.0:
            scale = total_before / total_after
            new_weights = [w * scale for w in new_weights]

        self.weights = new_weights
        self.states = new_states

    def weighted_diagonal(self, dim: int) -> np.ndarray:
        diag = np.zeros(dim)
        for w, psi in zip(self.weights, self.states):
            diag += w * np.abs(psi) ** 2
        return diag

    def density_matrix(self) -> np.ndarray:
        dim = len(self.states[0])
        rho = np.zeros((dim, dim), dtype=np.complex128)
        for w, psi in zip(self.weights, self.states):
            rho += w * np.outer(psi, psi.conj())
        return rho

    def entropy(self) -> float:
        vals = self.gram_eigenvalues()
        vals = vals[vals > 1e-15]
        total = vals.sum()
        if total <= 0.0:
            return 0.0
        p = vals / total
        return float(-np.sum(p * np.log(p)))

    def purity(self) -> float:
        vals = self.gram_eigenvalues()
        total = vals.sum()
        if total <= 0.0:
            return 0.0
        return float(np.sum((vals / total) ** 2))


def _apply_stinespring_with_branching(
    ensemble: _StatevectorEnsemble,
    U_stine: np.ndarray,
    d_anc: int,
    branch_tol: float,
) -> float:
    """Attach ancilla |0>, apply U_stine, branch on all ancilla outcomes.

    Updates ``ensemble`` in place with the (unnormalized-weight) branches and
    returns the total weight discarded due to ``branch_tol``.
    """
    new_weights: list[float] = []
    new_states: list[np.ndarray] = []
    discarded = 0.0
    for w, psi in zip(ensemble.weights, ensemble.states):
        d_sys = len(psi)
        psi_ext = np.zeros(d_anc * d_sys, dtype=np.complex128)
        psi_ext[:d_sys] = psi  # ancilla initialized to |0>
        psi_ext = U_stine @ psi_ext
        for k in range(d_anc):
            branch = psi_ext[k * d_sys : (k + 1) * d_sys]
            p_k = float(np.real(np.vdot(branch, branch)))
            w_branch = w * p_k
            if p_k <= branch_tol:
                discarded += w_branch
                continue
            new_weights.append(w_branch)
            new_states.append(branch / np.sqrt(p_k))
    ensemble.weights = new_weights
    ensemble.states = new_states
    return discarded


class QuditGKSLStatevectorSimulator:
    """Deterministic statevector-ensemble qudit GKSL simulator.

    Same Hamiltonian, Lindblad channels, Stinespring dilations and
    palindromic 2nd-order Trotter structure as QuditGKSLSimulator /
    QuditGKSLShotSimulator, but the state is a weighted ensemble of pure
    statevectors and the ancilla is updated by deterministic branching
    (no Born sampling, no density-matrix propagation).

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QuditGKSLStatevectorSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.dim = params.d ** params.N_molecules
        self.d_anc = params.d

        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total = self.H_0 + self.H_transfer
        self.lindblad_ops = build_lindblad_operators(params)

        self.n_system_qudits = params.N_molecules
        self.n_ancilla_qudits = len(self.lindblad_ops)
        self.n_total_qudits = self.n_system_qudits + self.n_ancilla_qudits

    def _precompute_unitaries(self, dt: float) -> None:
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines_half = [
            stinespring_unitary_from_lindblad(L_op, dt / 2, d_anc=self.d_anc)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _prepare_initial_statevector(self, state_type: str) -> np.ndarray:
        d = self.params.d
        N = self.params.N_molecules
        psi = np.zeros(self.dim, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            psi[1 * (d ** (N - 1)) + 1] = 1.0
        elif state_type == "all_triplet":
            psi[sum(1 * (d ** i) for i in range(N))] = 1.0
        elif state_type == "all_singlet":
            psi[sum(2 * (d ** i) for i in range(N))] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")
        return psi

    def _trotter_step_ensemble(
        self,
        ensemble: _StatevectorEnsemble,
        branch_tol: float,
        max_rank: int | None,
    ) -> float:
        """One palindromic 2nd-order Trotter step on the ensemble.

        Returns the total weight discarded by branch truncation in this step.
        """
        discarded = 0.0
        # Half Hamiltonian (unitary, applies to every ensemble member)
        ensemble.states = [self._U_H_half @ psi for psi in ensemble.states]
        # Forward half-step for all Lindblad channels
        for U_stine_half in self._U_stines_half:
            discarded += _apply_stinespring_with_branching(
                ensemble, U_stine_half, self.d_anc, branch_tol
            )
            if ensemble.size > (max_rank or self.dim):
                ensemble.compress(max_rank)
        # Reverse half-step (palindromic)
        for U_stine_half in reversed(self._U_stines_half):
            discarded += _apply_stinespring_with_branching(
                ensemble, U_stine_half, self.d_anc, branch_tol
            )
            if ensemble.size > (max_rank or self.dim):
                ensemble.compress(max_rank)
        # Half Hamiltonian
        ensemble.states = [self._U_H_half @ psi for psi in ensemble.states]
        return discarded

    def _populations(self, diag: np.ndarray) -> dict:
        return _populations_from_diagonal(diag, self.params)

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        branch_tol: float = _DEFAULT_BRANCH_TOL,
        max_rank: int | None = None,
    ) -> dict:
        """Run the deterministic statevector-ensemble GKSL simulation.

        Parameters:
            t_max: total simulation time
            n_steps: number of Trotter steps
            initial_state: initial state type
            branch_tol: ancilla branches with probability <= branch_tol are
                dropped (their weight is accumulated in discarded_weight)
            max_rank: maximum ensemble size after compression
                (default None = Hilbert space dimension, i.e. exact)

        Returns:
            dict compatible with the other GKSL simulators (times,
            populations, entropy, purity, trace, rho_final, ...), plus
            ensemble diagnostics (ensemble_size, discarded_weight).
        """
        start = time_module.time()
        dt = t_max / n_steps
        self._precompute_unitaries(dt)

        ensemble = _StatevectorEnsemble(
            [1.0], [self._prepare_initial_statevector(initial_state)]
        )

        times = [0.0]
        diag = ensemble.weighted_diagonal(self.dim)
        populations = [self._populations(diag)]
        entropies = [ensemble.entropy()]
        purities = [ensemble.purity()]
        traces = [ensemble.total_weight()]
        ensemble_sizes = [ensemble.size]
        discarded_total = 0.0

        for step in range(n_steps):
            discarded_total += self._trotter_step_ensemble(
                ensemble, branch_tol, max_rank
            )
            times.append((step + 1) * dt)
            diag = ensemble.weighted_diagonal(self.dim)
            populations.append(self._populations(diag))
            entropies.append(ensemble.entropy())
            purities.append(ensemble.purity())
            traces.append(ensemble.total_weight())
            ensemble_sizes.append(ensemble.size)

        rho_final = ensemble.density_matrix()
        elapsed = time_module.time() - start

        gates_per_step = (
            2 * (self.n_system_qudits + len(self.params.neighbors))
            + self.n_ancilla_qudits * 2
        )

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_final,
            "elapsed_time": elapsed,
            "method": "qudit_gksl_statevector_ensemble",
            "params": self.params.to_dict(),
            "n_system_qudits": self.n_system_qudits,
            "n_ancilla_qudits": self.n_ancilla_qudits,
            "n_total_qudits": self.n_total_qudits,
            "d_anc": self.d_anc,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "ensemble_size": ensemble_sizes,
            "final_ensemble_size": ensemble.size,
            "discarded_weight": discarded_total,
            "branch_tol": branch_tol,
        }


class QubitGKSLStatevectorSimulator:
    """Deterministic statevector-ensemble qubit GKSL simulator.

    Qubit-pair encoding (|00>=S0, |01>=T1, |10>=S1, |11>=forbidden) in the
    4^N-dimensional qubit space, with 2-level Stinespring ancillas updated by
    deterministic branching.  Mirrors QubitGKSLShotSimulator's operator
    construction; propagation is statevector-ensemble based (no
    density-matrix propagation).

    Parameters:
        params: GKSLPhysicalParameters (with_boson=False)
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QubitGKSLStatevectorSimulator is for non-boson model only"
            raise ValueError(msg)
        self.params = params
        self.N = params.N_molecules
        self.n_sys_qubits = 2 * self.N
        self.dim_qutrit = params.d ** params.N_molecules
        self.dim_qubit = (2 ** 2) ** params.N_molecules
        self.d_anc = 2

        self.H_0 = build_onsite_hamiltonian(params)
        self.H_transfer = build_transfer_hamiltonian(params)
        self.H_total_qutrit = self.H_0 + self.H_transfer
        self.lindblad_ops_qutrit = build_lindblad_operators(params)

        self._mapping = build_qubit_qutrit_mapping(self.N, params.d)
        self.H_total = embed_operator_in_qubit_space(
            self.H_total_qutrit, self._mapping, self.dim_qubit
        )
        self.lindblad_ops = [
            (embed_operator_in_qubit_space(L, self._mapping, self.dim_qubit), gamma)
            for L, gamma in self.lindblad_ops_qutrit
        ]

        self.n_ancilla = len(self.lindblad_ops)
        self.n_total_qubits = self.n_sys_qubits + self.n_ancilla

    def _precompute_unitaries(self, dt: float) -> None:
        self._U_H_half = expm(-1j * self.H_total * dt / 2)
        self._U_stines_half = [
            stinespring_unitary_from_lindblad(L_op, dt / 2)
            for L_op, _gamma in self.lindblad_ops
        ]

    def _prepare_initial_statevector(self, state_type: str) -> np.ndarray:
        d = self.params.d
        N = self.params.N_molecules
        psi_qt = np.zeros(self.dim_qutrit, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            psi_qt[1 * (d ** (N - 1)) + 1] = 1.0
        elif state_type == "all_triplet":
            psi_qt[sum(1 * (d ** i) for i in range(N))] = 1.0
        elif state_type == "all_singlet":
            psi_qt[sum(2 * (d ** i) for i in range(N))] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")
        return embed_statevector_in_qubit_space(psi_qt, self._mapping, self.dim_qubit)

    def _trotter_step_ensemble(
        self,
        ensemble: _StatevectorEnsemble,
        branch_tol: float,
        max_rank: int | None,
    ) -> float:
        discarded = 0.0
        rank_limit = max_rank or self.dim_qutrit  # rank(rho) <= 3^N always
        ensemble.states = [self._U_H_half @ psi for psi in ensemble.states]
        for U_stine_half in self._U_stines_half:
            discarded += _apply_stinespring_with_branching(
                ensemble, U_stine_half, self.d_anc, branch_tol
            )
            if ensemble.size > rank_limit:
                ensemble.compress(max_rank)
        for U_stine_half in reversed(self._U_stines_half):
            discarded += _apply_stinespring_with_branching(
                ensemble, U_stine_half, self.d_anc, branch_tol
            )
            if ensemble.size > rank_limit:
                ensemble.compress(max_rank)
        ensemble.states = [self._U_H_half @ psi for psi in ensemble.states]
        return discarded

    def _populations_and_forbidden(self, diag_qubit: np.ndarray) -> tuple[dict, float]:
        """Map qubit-space diagonal to qutrit populations + forbidden weight."""
        diag_qt = np.zeros(self.dim_qutrit)
        mapped_total = 0.0
        for qt_idx, qb_idx in self._mapping.items():
            diag_qt[qt_idx] = diag_qubit[qb_idx]
            mapped_total += diag_qubit[qb_idx]
        forbidden = float(diag_qubit.sum() - mapped_total)
        return _populations_from_diagonal_qubit(diag_qt, self.params), forbidden

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
        branch_tol: float = _DEFAULT_BRANCH_TOL,
        max_rank: int | None = None,
    ) -> dict:
        """Run the deterministic statevector-ensemble GKSL simulation (qubit)."""
        start = time_module.time()
        dt = t_max / n_steps
        self._precompute_unitaries(dt)

        ensemble = _StatevectorEnsemble(
            [1.0], [self._prepare_initial_statevector(initial_state)]
        )

        times = [0.0]
        diag = ensemble.weighted_diagonal(self.dim_qubit)
        pops0, forb0 = self._populations_and_forbidden(diag)
        populations = [pops0]
        forbidden_pops = [forb0]
        entropies = [ensemble.entropy()]
        purities = [ensemble.purity()]
        traces = [ensemble.total_weight()]
        ensemble_sizes = [ensemble.size]
        discarded_total = 0.0

        for step in range(n_steps):
            discarded_total += self._trotter_step_ensemble(
                ensemble, branch_tol, max_rank
            )
            times.append((step + 1) * dt)
            diag = ensemble.weighted_diagonal(self.dim_qubit)
            pops, forb = self._populations_and_forbidden(diag)
            populations.append(pops)
            forbidden_pops.append(forb)
            entropies.append(ensemble.entropy())
            purities.append(ensemble.purity())
            traces.append(ensemble.total_weight())
            ensemble_sizes.append(ensemble.size)

        rho_final_qubit = ensemble.density_matrix()
        # Extract to qutrit space for comparison with the other simulators
        # (same convention as QubitGKSLSimulator's rho_final)
        rho_final = extract_density_matrix_from_qubit_space(
            rho_final_qubit, self._mapping, self.dim_qutrit
        )
        elapsed = time_module.time() - start

        gates_per_step = (
            2 * (self.N + len(self.params.neighbors)) + self.n_ancilla * 2
        )

        return {
            "times": times,
            "populations": populations,
            "forbidden_populations": forbidden_pops,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_final,
            "rho_final_qubit": rho_final_qubit,
            "elapsed_time": elapsed,
            "method": "qubit_gksl_statevector_ensemble",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_sys_qubits,
            "n_ancilla": self.n_ancilla,
            "n_total_qubits": self.n_total_qubits,
            "estimated_gates_per_step": gates_per_step,
            "total_estimated_gates": gates_per_step * n_steps,
            "ensemble_size": ensemble_sizes,
            "final_ensemble_size": ensemble.size,
            "discarded_weight": discarded_total,
            "branch_tol": branch_tol,
        }
