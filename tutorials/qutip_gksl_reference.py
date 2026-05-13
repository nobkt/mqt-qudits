"""Independent QuTiP-based GKSL reference solver (A-2).

Why this module exists
----------------------
``STATUS_HONEST_2026-05.md`` records under **A-2** that the existing
``QubitGKSLSimulator`` and ``QuditGKSLSimulator`` are **not** algorithmically
independent: they share

* :func:`gksl_math_utils.build_onsite_hamiltonian`,
* :func:`gksl_math_utils.build_transfer_hamiltonian`,
* :func:`gksl_math_utils.build_lindblad_operators`,
* :mod:`stinespring_utils` (Stinespring dilation),
* :mod:`exact_local_channels` (exact local channel exponentiation),
* the same Strang+palindromic Trotter outer loop.

So an agreement between Qubit and Qudit simulators only validates the
qutrit→qubit embedding, not the open-system dynamics.

This module provides a **completely independent** reference
implementation of the same physical model using only QuTiP_ primitives:

* The Hamiltonian and every collapse operator are built from
  ``qutip.basis``, ``qutip.qeye``, ``qutip.tensor`` — the same physics,
  but constructed without re-using a single line from
  :mod:`gksl_math_utils`.
* The time evolution is performed by :func:`qutip.mesolve`, an
  **adaptive ODE** integrator (Adams / BDF) that does not use Trotter
  splitting, does not use Stinespring dilation, and shares no code with
  this repository.

A measurable agreement between this QuTiP reference and
:class:`QuditGKSLSimulator` (or, after embedding,
:class:`QubitGKSLSimulator`) therefore is a genuine independent
cross-validation.

What this module does **not** claim
-----------------------------------
* It does **not** make the qubit and qudit simulators independent of
  each other — they still share the same internal core, and that fact
  is recorded as-is in ``STATUS_HONEST_2026-05.md``.
* The boson-extended model (``with_boson=True``) is **out of scope**
  here: the same independent treatment for the boson case can be added
  later, but is not implemented in this module.

.. _QuTiP: https://qutip.org/
"""

from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Sequence


__all__ = [
    "QutipGKSLReference",
    "build_qutip_hamiltonian",
    "build_qutip_collapse_ops",
]


def _require_qutip():
    """Import QuTiP lazily so the module can be imported without it.

    Returns the ``qutip`` module on success; raises a clear ``ImportError``
    otherwise.
    """
    try:
        import qutip
    except ImportError as exc:  # pragma: no cover - import-time guard
        msg = (
            "qutip is required for the independent GKSL reference. "
            "Install with `pip install qutip` (>=5.0)."
        )
        raise ImportError(msg) from exc
    return qutip


def _site_op(qt, op_local: np.ndarray, site: int, N: int, d: int):
    """Embed a local d×d operator on ``site`` into the N-qudit register.

    Site 0 is the most significant qudit, matching the
    ``np.kron``-based ordering used by
    :func:`gksl_math_utils.build_single_site_operator` and by
    :func:`QuditGKSLSimulator.prepare_initial_state` (which constructs
    ``index = i_0 · d^(N-1) + i_1 · d^(N-2) + … + i_{N-1}``).
    """
    if op_local.shape != (d, d):
        msg = f"local op must be ({d},{d}); got {op_local.shape}"
        raise ValueError(msg)
    factors = []
    for k in range(N):
        if k == site:
            factors.append(qt.Qobj(op_local))
        else:
            factors.append(qt.qeye(d))
    return qt.tensor(*factors)


def _pair_op(qt, op_i: np.ndarray, op_j: np.ndarray, i: int, j: int, N: int, d: int):
    """Embed a tensor product ``op_i ⊗ op_j`` on sites ``i < j``."""
    if not 0 <= i < j < N:
        msg = f"need 0 ≤ i < j < N; got i={i}, j={j}, N={N}"
        raise ValueError(msg)
    factors = []
    for k in range(N):
        if k == i:
            factors.append(qt.Qobj(op_i))
        elif k == j:
            factors.append(qt.Qobj(op_j))
        else:
            factors.append(qt.qeye(d))
    return qt.tensor(*factors)


def _ket_bra(d: int, a: int, b: int) -> np.ndarray:
    m = np.zeros((d, d), dtype=np.complex128)
    m[a, b] = 1.0
    return m


def build_qutip_hamiltonian(params: GKSLPhysicalParameters):
    """Return the system Hamiltonian as a QuTiP ``Qobj``.

    Independently constructed from QuTiP primitives only:

    * On-site:  ``H_0 = Σ_i (E_T |1⟩⟨1| + E_S |2⟩⟨2|)_i``
    * Transfer: ``H_t = V Σ_⟨i,j⟩ (|0⟩⟨1|_i ⊗ |1⟩⟨0|_j + h.c.)``

    No reference to :mod:`gksl_math_utils`.
    """
    if params.with_boson:
        msg = "QutipGKSLReference does not support with_boson=True"
        raise ValueError(msg)
    qt = _require_qutip()

    d = params.d
    N = params.N_molecules

    # On-site energies
    h_local = np.diag(np.array([0.0, params.E_T, params.E_S], dtype=np.complex128))
    H = sum(
        (_site_op(qt, h_local, i, N, d) for i in range(N)),
        start=qt.tensor(*([qt.qeye(d)] * N)) * 0,
    )

    # Transfer (V * (|0><1|_i ⊗ |1><0|_j + h.c.))
    ket0_bra1 = _ket_bra(d, 0, 1)
    ket1_bra0 = _ket_bra(d, 1, 0)
    for (i, j) in params.neighbors:
        fwd = _pair_op(qt, ket0_bra1, ket1_bra0, i, j, N, d)
        H = H + params.V * (fwd + fwd.dag())

    return H


def build_qutip_collapse_ops(params: GKSLPhysicalParameters) -> list:
    """Return all Lindblad collapse operators as a list of QuTiP ``Qobj``.

    The list ordering and the ``√γ`` prefactors match the convention of
    :func:`gksl_math_utils.build_lindblad_operators`, but the operators
    themselves are constructed entirely from QuTiP primitives:

    * 2 × |neighbors| TTA-pair collapse operators (γ = γ_TTA / 2 each)
        - ``L_1 = √γ · |2⟩⟨1|_i ⊗ |0⟩⟨1|_j``
        - ``L_2 = √γ · |0⟩⟨1|_i ⊗ |2⟩⟨1|_j``
    * 5 × N single-site collapse operators
        - fluorescence       ``√Γ_fl  · |0⟩⟨2|``
        - phosphorescence    ``√Γ_ph  · |0⟩⟨1|``
        - internal conversion``√k_IC  · |0⟩⟨2|``
        - ISC S→T            ``√k_ISC_ST · |1⟩⟨2|``
        - ISC T→S            ``√k_ISC_TS · |0⟩⟨1|``

    The outer order is "all 5 single-site channel kinds × all sites"
    (kind-major, site-minor), matching
    :func:`gksl_math_utils.build_lindblad_operators`.
    """
    if params.with_boson:
        msg = "QutipGKSLReference does not support with_boson=True"
        raise ValueError(msg)
    qt = _require_qutip()

    d = params.d
    N = params.N_molecules
    c_ops: list = []

    # --- TTA pair channels ---
    for (i, j) in params.neighbors:
        gamma = params.gamma_TTA / 2.0
        sq = np.sqrt(gamma)
        c_ops.append(sq * _pair_op(qt, _ket_bra(d, 2, 1), _ket_bra(d, 0, 1), i, j, N, d))
        c_ops.append(sq * _pair_op(qt, _ket_bra(d, 0, 1), _ket_bra(d, 2, 1), i, j, N, d))

    # --- Single-site channels (kind-major, site-minor) ---
    singles = (
        (params.Gamma_fl, _ket_bra(d, 0, 2)),
        (params.Gamma_ph, _ket_bra(d, 0, 1)),
        (params.k_IC,     _ket_bra(d, 0, 2)),
        (params.k_ISC_ST, _ket_bra(d, 1, 2)),
        (params.k_ISC_TS, _ket_bra(d, 0, 1)),
    )
    for (gamma, op_local) in singles:
        sq = np.sqrt(gamma)
        for i in range(N):
            c_ops.append(sq * _site_op(qt, op_local, i, N, d))

    expected = 2 * len(params.neighbors) + 5 * N
    if len(c_ops) != expected:
        msg = f"Expected {expected} collapse ops, built {len(c_ops)}"
        raise AssertionError(msg)
    return c_ops


def _initial_state_indices(state_type: str, N: int, d: int) -> list[int]:
    """Return the per-site computational-basis indices for ``state_type``.

    Mirrors :meth:`QuditGKSLSimulator.prepare_initial_state` exactly.
    """
    if state_type == "edge_triplet":
        if N < 2:
            msg = "edge_triplet requires N_molecules >= 2"
            raise ValueError(msg)
        idx = [0] * N
        idx[0] = 1
        idx[N - 1] = 1
        return idx
    if state_type == "all_triplet":
        return [1] * N
    if state_type == "all_singlet":
        return [2] * N
    msg = f"Unknown state type: {state_type!r}"
    raise ValueError(msg)


class QutipGKSLReference:
    """Independent GKSL reference solver using ``qutip.mesolve``.

    The class deliberately exposes a small surface — only what is needed
    to cross-validate :class:`QuditGKSLSimulator` (and the qubit-space
    simulator after embedding).  The outputs match the format of
    :meth:`QuditGKSLSimulator.simulate` *only for the keys we care
    about*; the QuTiP solver does not, of course, expose Trotter
    diagnostics like Stinespring gate counts.

    Parameters
    ----------
    params:
        Physical parameters.  ``with_boson=True`` is not supported.
    options:
        Optional dict forwarded to ``qutip.SolverOptions``.  Defaults to
        a tight-tolerance ODE: ``atol = rtol = 1e-12``, ``nsteps = 100000``.
    """

    def __init__(
        self,
        params: GKSLPhysicalParameters,
        options: dict | None = None,
    ) -> None:
        if params.with_boson:
            msg = "QutipGKSLReference does not support with_boson=True"
            raise ValueError(msg)
        self.params = params
        self.d = params.d
        self.N = params.N_molecules
        self.dim = self.d ** self.N

        qt = _require_qutip()
        self._qt = qt
        self.H = build_qutip_hamiltonian(params)
        self.c_ops = build_qutip_collapse_ops(params)

        if options is None:
            options = {"atol": 1e-12, "rtol": 1e-12, "nsteps": 100000}
        self.options = options

    def prepare_initial_state(self, state_type: str = "edge_triplet"):
        """Return the initial density matrix as a QuTiP ``Qobj``."""
        qt = self._qt
        idx = _initial_state_indices(state_type, self.N, self.d)
        psi = qt.basis([self.d] * self.N, idx)
        return psi * psi.dag()

    def simulate(
        self,
        t_max: float,
        n_times: int = 51,
        initial_state: str = "edge_triplet",
        times: Sequence[float] | None = None,
    ) -> dict:
        """Evolve under the Lindblad master equation with ``qutip.mesolve``.

        Parameters
        ----------
        t_max:
            Final time.  Ignored when ``times`` is given.
        n_times:
            Number of equally spaced output times in ``[0, t_max]``
            (including both endpoints).  Ignored when ``times`` is
            given.  This controls only the **output sampling** — the
            internal ODE step size is chosen adaptively by ``mesolve``.
        initial_state:
            One of ``"edge_triplet"``, ``"all_triplet"``, ``"all_singlet"``.
        times:
            Optional explicit array of output times.  When provided,
            ``t_max`` and ``n_times`` are ignored.

        Returns
        -------
        dict with keys
            ``"times"``: ``list[float]`` of output times.
            ``"rho"``: ``list[np.ndarray]`` of full ``(dim, dim)``
                density matrices at each output time, as plain NumPy
                arrays (so the caller does not need QuTiP at the
                comparison site).
            ``"rho_final"``: the final density matrix as ``np.ndarray``.
            ``"trace"``: ``list[float]`` traces of every output state
                (a built-in honesty check; should be ≈1 throughout).
            ``"method"``: ``"qutip_mesolve"``.
            ``"params"``: dict snapshot of the input parameters.
        """
        qt = self._qt

        if times is None:
            tlist = np.linspace(0.0, float(t_max), int(n_times))
        else:
            tlist = np.asarray(times, dtype=float)
        rho0 = self.prepare_initial_state(initial_state)

        # QuTiP 5.x deprecated the dedicated SolverOptions class in favour
        # of a plain dict.  Pass the dict directly to avoid the
        # FutureWarning and keep forward compatibility.
        result = qt.mesolve(
            self.H,
            rho0,
            tlist,
            c_ops=self.c_ops,
            options=dict(self.options),
        )

        rhos = [np.asarray(s.full(), dtype=np.complex128) for s in result.states]
        traces = [float(np.real(np.trace(r))) for r in rhos]

        return {
            "times": tlist.tolist(),
            "rho": rhos,
            "rho_final": rhos[-1],
            "trace": traces,
            "method": "qutip_mesolve",
            "params": self.params.to_dict(),
        }
