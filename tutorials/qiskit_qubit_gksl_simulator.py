"""Qiskit-Aer-based qubit GKSL simulator (real circuit execution).

This module re-implements Scenarios 3 and 4 of the GKSL comparison
notebook so that the qubit-encoded open-system dynamics is **executed
as an actual quantum circuit** on a real backend, namely Qiskit Aer's
``density_matrix`` simulator.  The previous simulators
(:mod:`qubit_gksl_simulator`, :mod:`qubit_gksl_boson_simulator`)
applied the time evolution as a NumPy ``expm`` + Stinespring matrix
product without ever building a quantum circuit.

What this simulator does (and does not do)
------------------------------------------
For each Trotter step we build a Qiskit
:class:`qiskit.circuit.QuantumCircuit` containing:

1. ``UnitaryGate(expm(-i H_total · dt/2))`` on **all** system qubits
   (electronic + phonon, when applicable) — same Hamiltonian Strang
   half-step as the in-process simulator, but applied as a real
   gate by Qiskit Aer.
2. For every Lindblad channel (single-site or pair), a
   :class:`qiskit.quantum_info.Kraus` instruction on the relevant
   electronic qubit subset.  The Kraus operators are extracted
   exactly via Choi–Jamiolkowski from the
   :func:`exact_local_channels.build_local_dissipator_super`
   exponential built **in the qubit-pair–embedded local space**.
3. The same set of channels in **reverse order** (palindromic).
4. A second ``UnitaryGate`` for the closing Hamiltonian half-step.

The circuit is then executed on
``AerSimulator(method="density_matrix")`` with the previous step's ρ
fed in via ``set_density_matrix``; the output ρ is read back via
``save_density_matrix`` and used as the input for the next step.

This is a genuine Qiskit backend execution: every Trotter step is a
real :class:`qiskit.circuit.QuantumCircuit` whose instructions are
dispatched by Qiskit Aer.  No NumPy ``expm`` is called inside the
simulation loop; only the once-per-simulation pre-computation of
``U_H_half`` and the local Kraus operators uses ``expm`` /
``eigh`` — those are properties of the Trotter step itself,
independent of whether one runs the circuit on Aer, on the qubit
simulator NumPy path, or on real hardware that supports the same
gate set.

Sub-stochasticity (no closure needed)
------------------------------------
The qutrit→qubit-pair encoding ``|S0⟩→|00⟩, |T1⟩→|01⟩, |S1⟩→|10⟩``
leaves ``|11⟩`` unused.  When we embed an electronic Lindblad operator
``L`` into the qubit-pair space the resulting matrix has zero rows and
columns at index 3 (and analogously for the 16-dim pair space).  The
GKSL dissipator ``D[ρ] = L ρ L† − ½{L†L, ρ}`` constructed from such
``L`` is **automatically trace-preserving** on the full local qubit
space (the trace-cancellation identity holds element-wise), so
``expm(L_D · dt)`` is CPTP on the full local qubit space and its Choi
decomposition gives ``Σ_α K_α†K_α = I`` exactly — no closure
correction is needed.  We verify this numerically as a guardrail.

What this module does **not** claim
-----------------------------------
* This is a Qiskit Aer **simulator** of the qubit circuit, not a real
  hardware run.  Aer's ``density_matrix`` method evolves ρ classically
  on the host machine.
* Both this module and the original :class:`QubitGKSLSimulator` model
  the same physical system in the same qubit-pair encoding.  Agreement
  between them validates that this module's circuit construction
  faithfully reproduces the qubit-space density-matrix evolution; it
  does **not** independently validate the open-system physics (which
  is checked against the classical reference simulators elsewhere).
"""

from __future__ import annotations

import os
import sys
import time as time_module

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    build_H_total_boson,
    compute_populations_from_density_matrix,
    compute_purity,
    compute_von_neumann_entropy,
    extend_lindblad_operators,
    partial_trace_phonon,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_simulator import (
    _QUTRIT_TO_QUBIT_PAIR,
    build_qubit_qutrit_mapping,
    embed_operator_in_qubit_space,
    embed_density_matrix_in_qubit_space,
    extract_density_matrix_from_qubit_space,
    compute_forbidden_state_population,
)
from qubit_gksl_boson_simulator import (
    _build_extended_qubit_mapping,
    _embed_extended_density_matrix,
)
from exact_local_channels import (
    build_local_dissipator_super,
    get_local_lindblad_ops,
)
from dmsim_kraus_helpers import kraus_from_local_superoperator

# Local qubit-pair dimension per molecule (2 qubits = 4-dim).
_DIM_QB_PAIR = 4
_DIM_QT = 3

# Kraus-channel CPTP closure tolerance used when sanity-checking the
# augmented qubit-space Kraus list before handing it to Qiskit.
_CPTP_CHECK_TOL = 1e-9


# ---------------------------------------------------------------------------
# Local qutrit→qubit-pair embedding (single molecule and pair-of-molecules)
# ---------------------------------------------------------------------------


def _local_single_mapping(d: int) -> dict[int, int]:
    """qutrit local idx (0..d-1) → qubit-pair local idx (0..3)."""
    if d != _DIM_QT:
        msg = f"_local_single_mapping currently only supports d=3, got {d}"
        raise ValueError(msg)
    return dict(_QUTRIT_TO_QUBIT_PAIR)


def _local_pair_mapping(d: int) -> dict[int, int]:
    """qutrit-pair local idx (0..d²-1) → qubit-quad local idx (0..15).

    The qutrit local index uses the convention of
    :func:`exact_local_channels.get_local_lindblad_ops`, which builds
    pair Lindblad operators as ``np.kron(L_at_i, L_at_j)`` — i.e. the
    qutrit-at-i is the *most significant* digit:
    ``qt_idx = a · d + b`` where ``a`` is the qutrit at site ``i``.

    Likewise we encode the qubit-pair local index as
    ``qb_idx = encoding(a) · 4 + encoding(b)`` — qubit pair at i is most
    significant in the 16-dim space.
    """
    out: dict[int, int] = {}
    for a in range(d):
        for b in range(d):
            qt_idx = a * d + b
            qb_idx = _QUTRIT_TO_QUBIT_PAIR[a] * _DIM_QB_PAIR + _QUTRIT_TO_QUBIT_PAIR[b]
            out[qt_idx] = qb_idx
    return out


def _embed_local_operator(op_qt: np.ndarray, mapping: dict[int, int], dim_qb: int) -> np.ndarray:
    """Place qutrit-local operator into qubit-local space (forbidden rows/cols = 0)."""
    op_qb = np.zeros((dim_qb, dim_qb), dtype=np.complex128)
    for i_qt, i_qb in mapping.items():
        for j_qt, j_qb in mapping.items():
            op_qb[i_qb, j_qb] = op_qt[i_qt, j_qt]
    return op_qb


# ---------------------------------------------------------------------------
# Build a CPTP qubit Kraus channel from an electronic Lindblad
# ---------------------------------------------------------------------------


def _build_qubit_kraus_for_local_channel(
    L_local_qt: np.ndarray,
    d_local_qt: int,
    dt_half: float,
    *,
    is_pair: bool,
) -> list[np.ndarray]:
    """Return Kraus operators of ``expm(L_D^local · dt/2)`` in qubit-pair encoding.

    Steps:
    1. Embed the qutrit-local Lindblad operator into the qubit-local
       (4-dim or 16-dim) space.
    2. Build the local GKSL dissipator superoperator in qubit space and
       exponentiate it: ``expm(L_D^local_qb · dt/2)``.
    3. Choi-Jamiolkowski decompose into ``{K_α^physical}``.
    4. Append a closure Kraus operator equal to the projector onto the
       forbidden subspace so that ``Σ K†K = I`` on the **full** local
       qubit space (mathematical identity — no effect on physical
       states; required for ``qiskit.quantum_info.Kraus.to_instruction``
       to accept the channel).
    """
    if is_pair:
        d_local_qb = _DIM_QB_PAIR ** 2  # 16
        mapping = _local_pair_mapping(_DIM_QT)
    else:
        d_local_qb = _DIM_QB_PAIR  # 4
        mapping = _local_single_mapping(_DIM_QT)

    L_qb = _embed_local_operator(L_local_qt, mapping, d_local_qb)
    L_D_qb = build_local_dissipator_super(L_qb, d_local_qb)
    M_half_qb = expm(L_D_qb * dt_half)

    kraus_phys = kraus_from_local_superoperator(M_half_qb, d_local_qb)

    # The GKSL dissipator ``L ρ L† − ½{L†L, ρ}`` is automatically
    # trace-preserving — even for a sub-stochastic embedded ``L`` whose
    # forbidden rows/columns are zero — so ``expm(L_D · dt)`` is CPTP on
    # the full local qubit space without any closure correction.
    # Verified numerically below as a guardrail; no heuristic adjustment.
    completeness = np.zeros((d_local_qb, d_local_qb), dtype=np.complex128)
    for K in kraus_phys:
        completeness += K.conj().T @ K
    defect = float(np.linalg.norm(completeness - np.eye(d_local_qb)))
    if defect > _CPTP_CHECK_TOL:
        msg = (
            "qubit Kraus closure failed: ||Σ K†K - I||_F = "
            f"{defect:.3e} > {_CPTP_CHECK_TOL:g}.  This indicates an "
            "arithmetic bug in the embedded dissipator — no heuristic "
            "correction applied."
        )
        raise AssertionError(msg)

    return kraus_phys


# ---------------------------------------------------------------------------
# Qiskit-Aer GKSL simulator (no boson) — Scenario 3
# ---------------------------------------------------------------------------


class QiskitQubitGKSLSimulator:
    """Qubit GKSL simulator that runs each Trotter step on Qiskit Aer.

    Qubit register layout (LSB first, matches Qiskit little-endian):
    qubits ``[2i, 2i+1]`` = molecule ``i`` (encoding bit 0, bit 1).
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QiskitQubitGKSLSimulator is for the non-boson model only"
            raise ValueError(msg)
        # Eagerly import Qiskit so that import errors surface at
        # construction time rather than at the first ``simulate`` call.
        from qiskit_aer import AerSimulator

        self.params = params
        self.N = params.N_molecules
        self.n_qubits = 2 * self.N
        self.dim_qubit = 2 ** self.n_qubits
        self.dim_qutrit = params.d ** params.N_molecules

        # Qutrit-space operators
        self.H_total_qt = (
            build_onsite_hamiltonian(params) + build_transfer_hamiltonian(params)
        )
        self.lindblad_ops_qt = build_lindblad_operators(params)

        # Mapping qutrit→qubit (full system)
        self._mapping = build_qubit_qutrit_mapping(self.N, params.d)

        # Embed Hamiltonian into qubit space (used for the half-step gate).
        self.H_total_qb = embed_operator_in_qubit_space(
            self.H_total_qt, self._mapping, self.dim_qubit
        )

        # Local Lindblad operators (qutrit) — used to build qubit-local Kraus.
        self._local_ops_qt = get_local_lindblad_ops(params)

        # Cache Aer simulator instance.
        self._aer = AerSimulator(method="density_matrix")

        # Qubit indices for each molecule.
        self._mol_qubits = [(2 * i, 2 * i + 1) for i in range(self.N)]

    # ------------------------------------------------------------------
    # Pre-compute time-step-dependent gates
    # ------------------------------------------------------------------

    def _precompute(self, dt: float) -> None:
        """Build U_H_half and the per-channel qubit Kraus lists."""
        self._U_H_half = expm(-1j * self.H_total_qb * dt / 2)

        kraus_per_channel: list[tuple[tuple[int, ...], list[np.ndarray], str]] = []
        for sites, L_local_qt, kind in self._local_ops_qt:
            kraus = _build_qubit_kraus_for_local_channel(
                L_local_qt,
                d_local_qt=L_local_qt.shape[0],
                dt_half=dt / 2.0,
                is_pair=(kind == "pair"),
            )
            kraus_per_channel.append((sites, kraus, kind))
        self._kraus_per_channel = kraus_per_channel

    # ------------------------------------------------------------------
    # Build one Trotter-step Qiskit circuit
    # ------------------------------------------------------------------

    def _build_trotter_circuit(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import Kraus

        qc = QuantumCircuit(self.n_qubits)

        # 1. Hamiltonian half-step on all qubits.
        # Qiskit's UnitaryGate uses the same little-endian basis ordering
        # as our embedding (q0 = LSB), so we apply it to qubits [0..n-1]
        # in their natural order.
        all_qubits = list(range(self.n_qubits))
        qc.append(UnitaryGate(self._U_H_half), all_qubits)

        # 2. Forward palindromic pass.
        for sites, kraus_ops, kind in self._kraus_per_channel:
            self._append_kraus(qc, Kraus, sites, kraus_ops, kind)

        # 3. Reverse palindromic pass.
        for sites, kraus_ops, kind in reversed(self._kraus_per_channel):
            self._append_kraus(qc, Kraus, sites, kraus_ops, kind)

        # 4. Closing Hamiltonian half-step.
        qc.append(UnitaryGate(self._U_H_half), all_qubits)

        return qc

    def _append_kraus(self, qc, Kraus, sites: tuple[int, ...], kraus_ops, kind: str) -> None:
        """Append a Kraus instruction targeting the qubits of the given molecules."""
        if kind == "single":
            qbits = list(self._mol_qubits[sites[0]])
        else:
            i, j = sites
            qbits = list(self._mol_qubits[i]) + list(self._mol_qubits[j])
        # Qiskit qubit-list ordering: same little-endian convention as our
        # local mapping (first qubit in the list is the LSB of the local
        # operator's basis index).
        qc.append(Kraus(kraus_ops).to_instruction(), qbits)

    # ------------------------------------------------------------------
    # Initial state in qubit space
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        d = self.params.d
        N = self.N
        psi = np.zeros(d ** N, dtype=np.complex128)
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
        rho_qt = np.outer(psi, psi.conj())
        return embed_density_matrix_in_qubit_space(rho_qt, self._mapping, self.dim_qubit)

    # ------------------------------------------------------------------
    # Main simulation loop
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run the qubit GKSL simulation on Qiskit Aer (density-matrix method)."""
        from qiskit import QuantumCircuit, transpile
        from qiskit.quantum_info import DensityMatrix

        start = time_module.time()

        dt = t_max / n_steps
        self._precompute(dt)

        rho = self.prepare_initial_state(initial_state)
        rho_qt = extract_density_matrix_from_qubit_space(rho, self._mapping, self.dim_qutrit)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_qt, self.params)]
        entropies = [compute_von_neumann_entropy(rho_qt)]
        purities = [compute_purity(rho_qt)]
        traces = [float(np.real(np.trace(rho_qt)))]
        forbidden_pops: list[float] = [
            compute_forbidden_state_population(rho, self._mapping)
        ]

        # Build the Trotter-step circuit once and reuse it every step
        # (the circuit's gate matrices do not depend on the input state).
        step_circ = self._build_trotter_circuit()

        for step in range(n_steps):
            # Re-bind the initial-state instruction for this step's ρ.
            # We rebuild the small wrapper each step to update the
            # set_density_matrix payload — gate matrices are unchanged.
            wrapped = QuantumCircuit(self.n_qubits)
            wrapped.set_density_matrix(DensityMatrix(rho))
            wrapped.compose(step_circ, inplace=True)
            wrapped.save_density_matrix()
            result = self._aer.run(transpile(wrapped, self._aer)).result()
            rho = np.asarray(result.data(0)["density_matrix"], dtype=np.complex128)

            rho_qt = extract_density_matrix_from_qubit_space(
                rho, self._mapping, self.dim_qutrit
            )
            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_qt))))
            populations.append(
                compute_populations_from_density_matrix(rho_qt, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_qt))
            purities.append(compute_purity(rho_qt))
            forbidden_pops.append(
                compute_forbidden_state_population(rho, self._mapping)
            )

        elapsed = time_module.time() - start

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_qt,
            "elapsed_time": elapsed,
            "method": "qiskit_qubit_gksl",
            "backend": "qiskit_aer:density_matrix",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_qubits,
            "n_total_qubits": self.n_qubits,
            "dim_qubit_space": self.dim_qubit,
            "forbidden_state_population": forbidden_pops,
        }


# ---------------------------------------------------------------------------
# Qiskit-Aer GKSL simulator with phonons — Scenario 4
# ---------------------------------------------------------------------------


class QiskitQubitGKSLBosonSimulator:
    """Qubit GKSL+phonon simulator on Qiskit Aer.

    Encodes each electronic site as 2 qubits and each phonon mode as
    ``ceil(log2(n_max+1))`` qubit(s).  For the notebook setting
    (``N=2, n_max=1``) the register is 4 electronic + 2 phonon qubits
    (6 qubits, 64-dim density matrix).

    Qubit register layout (Qiskit little-endian):
    * qubits ``[0 .. n_ph_qubits-1]`` — phonon qubits (mode 0 LSB first,
      then mode 1, ...).
    * qubits ``[n_ph_qubits .. n_ph_qubits + n_el_qubits - 1]`` —
      electronic qubits (mol 0 = qubits ``[n_ph, n_ph+1]``, mol 1 =
      qubits ``[n_ph+2, n_ph+3]``, ...).

    This matches the Kronecker ordering ``el ⊗ ph`` used by
    :func:`gksl_math_utils.build_H_total_boson` and the extended qubit
    mapping ``qb_idx = el_qb_idx · dim_ph + ph_idx``.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if not params.with_boson:
            msg = "QiskitQubitGKSLBosonSimulator requires with_boson=True"
            raise ValueError(msg)
        from qiskit_aer import AerSimulator

        self.params = params
        self.N = params.N_molecules

        # Phonon qubits per mode and totals.
        n_ph_per_mode = int(np.ceil(np.log2(params.n_max + 1)))
        if (1 << n_ph_per_mode) != (params.n_max + 1):
            msg = (
                "QiskitQubitGKSLBosonSimulator currently requires "
                "n_max+1 to be a power of 2 so that the phonon Fock "
                f"space embeds exactly into n_ph_per_mode qubits; got "
                f"n_max={params.n_max}"
            )
            raise ValueError(msg)
        self.n_ph_per_mode = n_ph_per_mode
        self.n_ph_qubits = n_ph_per_mode * params.N_molecules
        self.n_el_qubits = 2 * params.N_molecules
        self.n_qubits = self.n_ph_qubits + self.n_el_qubits

        self.dim_el_qubit = (1 << self.n_el_qubits)
        self.dim_ph = (params.n_max + 1) ** params.N_molecules
        self.dim_qubit_total = self.dim_el_qubit * self.dim_ph
        self.dim_el_qutrit = params.d ** params.N_molecules

        # Build extended Hamiltonian in qutrit-electronic + Fock-phonon space.
        self.H_total_qt_ext = build_H_total_boson(params)

        # Build qutrit→qubit mapping (electronic only).
        self._mapping_el = build_qubit_qutrit_mapping(self.N, params.d)
        # Build extended-space mapping (el-qt ⊗ ph) → (el-qb ⊗ ph).
        self._ext_mapping = _build_extended_qubit_mapping(
            self._mapping_el, self.dim_el_qubit, self.dim_ph
        )

        # Embed extended Hamiltonian into qubit-extended space.
        self.H_total_qb_ext = self._embed_ext(self.H_total_qt_ext)

        # Local electronic Lindblad operators (twin params with with_boson=False).
        self._params_el = GKSLPhysicalParameters(
            E_T=params.E_T, E_S=params.E_S, V=params.V,
            gamma_TTA=params.gamma_TTA,
            Gamma_fl=params.Gamma_fl, Gamma_ph=params.Gamma_ph,
            k_IC=params.k_IC, k_ISC_ST=params.k_ISC_ST, k_ISC_TS=params.k_ISC_TS,
            N_molecules=params.N_molecules, d=params.d, with_boson=False,
        )
        self._local_ops_qt = get_local_lindblad_ops(self._params_el)

        self._aer = AerSimulator(method="density_matrix")

        # Qubit indices: phonons first (LSB), then electronic.
        self._mol_qubits = [
            (self.n_ph_qubits + 2 * i, self.n_ph_qubits + 2 * i + 1)
            for i in range(self.N)
        ]

    # ------------------------------------------------------------------
    # Helpers: embed/extract qutrit-extended ↔ qubit-extended
    # ------------------------------------------------------------------

    def _embed_ext(self, op_qt_ext: np.ndarray) -> np.ndarray:
        out = np.zeros((self.dim_qubit_total, self.dim_qubit_total), dtype=np.complex128)
        for i_qt, i_qb in self._ext_mapping.items():
            for j_qt, j_qb in self._ext_mapping.items():
                out[i_qb, j_qb] = op_qt_ext[i_qt, j_qt]
        return out

    def _embed_ext_density(self, rho_qt_ext: np.ndarray) -> np.ndarray:
        return _embed_extended_density_matrix(
            rho_qt_ext, self._ext_mapping, self.dim_qubit_total
        )

    def _extract_ext_density(self, rho_qb_ext: np.ndarray) -> np.ndarray:
        out = np.zeros(
            (self.dim_el_qutrit * self.dim_ph, self.dim_el_qutrit * self.dim_ph),
            dtype=np.complex128,
        )
        for i_qt, i_qb in self._ext_mapping.items():
            for j_qt, j_qb in self._ext_mapping.items():
                out[i_qt, j_qt] = rho_qb_ext[i_qb, j_qb]
        return out

    # ------------------------------------------------------------------
    # Pre-compute Trotter-step gates and Kraus channels
    # ------------------------------------------------------------------

    def _precompute(self, dt: float) -> None:
        self._U_H_half = expm(-1j * self.H_total_qb_ext * dt / 2)

        kraus_per_channel: list[tuple[tuple[int, ...], list[np.ndarray], str]] = []
        for sites, L_local_qt, kind in self._local_ops_qt:
            kraus = _build_qubit_kraus_for_local_channel(
                L_local_qt,
                d_local_qt=L_local_qt.shape[0],
                dt_half=dt / 2.0,
                is_pair=(kind == "pair"),
            )
            kraus_per_channel.append((sites, kraus, kind))
        self._kraus_per_channel = kraus_per_channel

    # ------------------------------------------------------------------
    # Trotter-step circuit
    # ------------------------------------------------------------------

    def _build_trotter_circuit(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import Kraus

        qc = QuantumCircuit(self.n_qubits)
        all_qubits = list(range(self.n_qubits))
        qc.append(UnitaryGate(self._U_H_half), all_qubits)

        for sites, kraus_ops, kind in self._kraus_per_channel:
            self._append_kraus(qc, Kraus, sites, kraus_ops, kind)
        for sites, kraus_ops, kind in reversed(self._kraus_per_channel):
            self._append_kraus(qc, Kraus, sites, kraus_ops, kind)

        qc.append(UnitaryGate(self._U_H_half), all_qubits)
        return qc

    def _append_kraus(self, qc, Kraus, sites: tuple[int, ...], kraus_ops, kind: str) -> None:
        if kind == "single":
            qbits = list(self._mol_qubits[sites[0]])
        else:
            i, j = sites
            qbits = list(self._mol_qubits[i]) + list(self._mol_qubits[j])
        qc.append(Kraus(kraus_ops).to_instruction(), qbits)

    # ------------------------------------------------------------------
    # Initial state
    # ------------------------------------------------------------------

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        d = self.params.d
        N = self.N
        psi_el = np.zeros(d ** N, dtype=np.complex128)
        if state_type == "edge_triplet":
            if N < 2:
                msg = "edge_triplet requires N_molecules >= 2"
                raise ValueError(msg)
            psi_el[1 * (d ** (N - 1)) + 1] = 1.0
        elif state_type == "all_triplet":
            psi_el[sum(1 * (d ** i) for i in range(N))] = 1.0
        elif state_type == "all_singlet":
            psi_el[sum(2 * (d ** i) for i in range(N))] = 1.0
        else:
            raise ValueError(f"Unknown state type: {state_type}")

        psi_ph = np.zeros(self.dim_ph, dtype=np.complex128)
        psi_ph[0] = 1.0
        psi_total_qt = np.kron(psi_el, psi_ph)
        rho_qt_ext = np.outer(psi_total_qt, psi_total_qt.conj())
        return self._embed_ext_density(rho_qt_ext)

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------

    def simulate(self, t_max: float, n_steps: int, initial_state: str = "edge_triplet") -> dict:
        from qiskit import QuantumCircuit, transpile
        from qiskit.quantum_info import DensityMatrix

        start = time_module.time()
        dt = t_max / n_steps
        self._precompute(dt)

        rho = self.prepare_initial_state(initial_state)
        rho_qt_ext = self._extract_ext_density(rho)
        rho_el_qt = partial_trace_phonon(rho_qt_ext, self.dim_el_qutrit, self.dim_ph)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_el_qt, self.params)]
        entropies = [compute_von_neumann_entropy(rho_el_qt)]
        purities = [compute_purity(rho_el_qt)]
        traces = [float(np.real(np.trace(rho_el_qt)))]

        step_circ = self._build_trotter_circuit()

        for step in range(n_steps):
            wrapped = QuantumCircuit(self.n_qubits)
            wrapped.set_density_matrix(DensityMatrix(rho))
            wrapped.compose(step_circ, inplace=True)
            wrapped.save_density_matrix()
            result = self._aer.run(transpile(wrapped, self._aer)).result()
            rho = np.asarray(result.data(0)["density_matrix"], dtype=np.complex128)

            rho_qt_ext = self._extract_ext_density(rho)
            rho_el_qt = partial_trace_phonon(rho_qt_ext, self.dim_el_qutrit, self.dim_ph)
            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_el_qt))))
            populations.append(
                compute_populations_from_density_matrix(rho_el_qt, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_el_qt))
            purities.append(compute_purity(rho_el_qt))

        elapsed = time_module.time() - start

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_el_qt,
            "elapsed_time": elapsed,
            "method": "qiskit_qubit_gksl_boson",
            "backend": "qiskit_aer:density_matrix",
            "params": self.params.to_dict(),
            "n_el_qubits": self.n_el_qubits,
            "n_ph_qubits": self.n_ph_qubits,
            "n_total_qubits": self.n_qubits,
            "dim_qubit_total": self.dim_qubit_total,
        }


# ---------------------------------------------------------------------------
# Qiskit-Aer ancilla-updating Stinespring simulator — §5.1-3
# ---------------------------------------------------------------------------


class QiskitQubitGKSLStinespringSimulator:
    """Qubit GKSL simulator with a **mid-circuit-updated ancilla qubit**.

    This is the qubit-side counterpart of the qudit ancilla-updating
    Stinespring path
    (``QuditGKSLSimulator(algorithm='stinespring', execute_on_backend='dmsim')``):
    each Trotter step is a real Qiskit circuit on ``2N system qubits + 1
    ancilla qubit`` in which every Lindblad channel is realised by

    1. its **local Stinespring dilation unitary**
       ``U_α = expm(-i √(dt/2) G_α)`` with ``G_α = [[0, L†],[L, 0]]``
       built in the qubit-pair–embedded local space (8×8 for single-site
       channels, 32×32 for TTA-pair channels, ancilla = most significant
       qubit), applied as a ``UnitaryGate`` on ``local qubits + [ancilla]``
       (Qiskit little-endian: first listed qubit is the LSB), followed by
    2. a **native Qiskit ``reset`` instruction on the ancilla** — Aer's
       ``density_matrix`` method executes ``reset`` deterministically as
       the exact CPTP channel ``{K_k = |0⟩⟨k|}`` (trace out + re-prepare
       ``|0⟩``), so the ancilla is genuinely updated between channels
       inside the circuit.

    The channel set and palindromic ordering are identical to
    :class:`QiskitQubitGKSLSimulator`; only the channel realisation
    differs (1st-order Stinespring dilation + reset instead of the exact
    2nd-order Kraus channel).  Consequently the result must agree with
    the NumPy qubit Stinespring simulator
    (:class:`qubit_gksl_simulator.QubitGKSLSimulator`) at round-off
    level, and with the exact dynamics up to the usual O(dt) Stinespring
    error.
    """

    def __init__(self, params: GKSLPhysicalParameters) -> None:
        if params.with_boson:
            msg = "QiskitQubitGKSLStinespringSimulator is for the non-boson model only"
            raise ValueError(msg)
        from qiskit_aer import AerSimulator

        self.params = params
        self.N = params.N_molecules
        self.n_sys_qubits = 2 * self.N
        self.n_qubits = self.n_sys_qubits + 1  # + 1 re-used ancilla qubit
        self.anc = self.n_sys_qubits  # ancilla qubit index (last)
        self.dim_qubit = 2 ** self.n_sys_qubits
        self.dim_qutrit = params.d ** params.N_molecules

        self.H_total_qt = (
            build_onsite_hamiltonian(params) + build_transfer_hamiltonian(params)
        )
        self._mapping = build_qubit_qutrit_mapping(self.N, params.d)
        self.H_total_qb = embed_operator_in_qubit_space(
            self.H_total_qt, self._mapping, self.dim_qubit
        )
        self._local_ops_qt = get_local_lindblad_ops(params)
        self._aer = AerSimulator(method="density_matrix")
        self._mol_qubits = [(2 * i, 2 * i + 1) for i in range(self.N)]

    def _precompute(self, dt: float) -> None:
        """Build U_H_half and the per-channel local dilation unitaries."""
        self._U_H_half = expm(-1j * self.H_total_qb * dt / 2)

        def _local_stinespring_qb(L_qb: np.ndarray, dtau: float) -> np.ndarray:
            d_loc = L_qb.shape[0]
            G = np.zeros((2 * d_loc, 2 * d_loc), dtype=np.complex128)
            G[:d_loc, d_loc:] = L_qb.conj().T
            G[d_loc:, :d_loc] = L_qb
            U = expm(-1j * np.sqrt(dtau) * G)
            residual = np.linalg.norm(U.conj().T @ U - np.eye(2 * d_loc), ord="fro")
            if residual >= 1e-10:
                msg = (
                    "Local Stinespring dilation unitary failed unitarity "
                    f"check: ||U†U - I||_F = {residual:.3e}"
                )
                raise ValueError(msg)
            return U

        dilations: list[tuple[tuple[int, ...], np.ndarray, str]] = []
        for sites, L_local_qt, kind in self._local_ops_qt:
            if kind == "pair":
                mapping = _local_pair_mapping(_DIM_QT)
                d_local_qb = _DIM_QB_PAIR ** 2
            else:
                mapping = _local_single_mapping(_DIM_QT)
                d_local_qb = _DIM_QB_PAIR
            L_qb = _embed_local_operator(L_local_qt, mapping, d_local_qb)
            dilations.append((sites, _local_stinespring_qb(L_qb, dt / 2.0), kind))
        self._dilations = dilations

    def _build_trotter_circuit(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate

        qc = QuantumCircuit(self.n_qubits)
        sys_qubits = list(range(self.n_sys_qubits))

        qc.append(UnitaryGate(self._U_H_half), sys_qubits)

        for sites, U_loc, kind in list(self._dilations) + list(reversed(self._dilations)):
            if kind == "single":
                qbits = list(self._mol_qubits[sites[0]])
            else:
                i, j = sites
                qbits = list(self._mol_qubits[i]) + list(self._mol_qubits[j])
            # Ancilla is the most significant leg of U_loc (index
            # ``anc_level · d_local + local_idx``), so it goes LAST in the
            # little-endian Qiskit qubit list.
            qc.append(UnitaryGate(U_loc), [*qbits, self.anc])
            # Mid-circuit ancilla update: native Qiskit reset instruction.
            qc.reset(self.anc)

        qc.append(UnitaryGate(self._U_H_half), sys_qubits)
        return qc

    def prepare_initial_state(self, state_type: str = "edge_triplet") -> np.ndarray:
        """Initial ρ on the 2N system qubits (ancilla is added per step)."""
        d = self.params.d
        N = self.N
        psi = np.zeros(d ** N, dtype=np.complex128)
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
        rho_qt = np.outer(psi, psi.conj())
        return embed_density_matrix_in_qubit_space(rho_qt, self._mapping, self.dim_qubit)

    def simulate(
        self,
        t_max: float,
        n_steps: int,
        initial_state: str = "edge_triplet",
    ) -> dict:
        """Run the ancilla-updating Stinespring dynamics on Qiskit Aer."""
        from qiskit import QuantumCircuit, transpile
        from qiskit.quantum_info import DensityMatrix

        start = time_module.time()

        dt = t_max / n_steps
        self._precompute(dt)

        rho = self.prepare_initial_state(initial_state)
        rho_qt = extract_density_matrix_from_qubit_space(rho, self._mapping, self.dim_qutrit)

        times: list[float] = [0.0]
        populations = [compute_populations_from_density_matrix(rho_qt, self.params)]
        entropies = [compute_von_neumann_entropy(rho_qt)]
        purities = [compute_purity(rho_qt)]
        traces = [float(np.real(np.trace(rho_qt)))]
        forbidden_pops: list[float] = [
            compute_forbidden_state_population(rho, self._mapping)
        ]

        step_circ = self._build_trotter_circuit()
        anc0 = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=np.complex128)

        for step in range(n_steps):
            # ρ_tot = |0><0|_anc ⊗ ρ_sys (ancilla = most significant qubit).
            rho_tot = np.kron(anc0, rho)
            wrapped = QuantumCircuit(self.n_qubits)
            wrapped.set_density_matrix(DensityMatrix(rho_tot))
            wrapped.compose(step_circ, inplace=True)
            wrapped.save_density_matrix()
            result = self._aer.run(transpile(wrapped, self._aer)).result()
            rho_tot = np.asarray(result.data(0)["density_matrix"], dtype=np.complex128)

            # Partial trace over the ancilla (most significant qubit).
            rho_t = rho_tot.reshape(2, self.dim_qubit, 2, self.dim_qubit)
            rho = np.einsum("kakb->ab", rho_t)

            rho_qt = extract_density_matrix_from_qubit_space(
                rho, self._mapping, self.dim_qutrit
            )
            times.append((step + 1) * dt)
            traces.append(float(np.real(np.trace(rho_qt))))
            populations.append(
                compute_populations_from_density_matrix(rho_qt, self.params)
            )
            entropies.append(compute_von_neumann_entropy(rho_qt))
            purities.append(compute_purity(rho_qt))
            forbidden_pops.append(
                compute_forbidden_state_population(rho, self._mapping)
            )

        elapsed = time_module.time() - start

        return {
            "times": times,
            "populations": populations,
            "entropy": entropies,
            "purity": purities,
            "trace": traces,
            "rho_final": rho_qt,
            "elapsed_time": elapsed,
            "method": "qiskit_qubit_gksl_stinespring_ancilla",
            "backend": "qiskit_aer:density_matrix",
            "params": self.params.to_dict(),
            "n_sys_qubits": self.n_sys_qubits,
            "n_total_qubits": self.n_qubits,
            "n_ancilla_qubits_backend_circuit": 1,
            "dim_qubit_space": self.dim_qubit,
            "forbidden_state_population": forbidden_pops,
        }


__all__ = [
    "QiskitQubitGKSLSimulator",
    "QiskitQubitGKSLBosonSimulator",
    "QiskitQubitGKSLStinespringSimulator",
]
