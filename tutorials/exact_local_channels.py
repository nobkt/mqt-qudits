"""Exact local-channel exponentiation for the Qudit GKSL simulator (B-1).

Background
----------
The default :class:`QuditGKSLSimulator` (Scenario 5 in the comparison
notebook) approximates each individual Lindblad channel
``E_α(dt) = exp(L_{D_α} · dt)`` by a Stinespring dilation
``U_α = expm(-i · sqrt(dt) · G_α)``.  That dilation is mathematically only a
**1st-order** approximation of the channel:

    E_α(dt)(ρ) = ρ + dt · (L_α ρ L_α† − ½{L_α†L_α, ρ}) + O(dt²)
    Stinespring(dt)(ρ) = ρ + dt · (L_α ρ L_α† − ½{L_α†L_α, ρ}) + O(dt²)

so the ``O(dt²)`` next term differs between the two and the global
convergence collapses to ``O(dt)`` regardless of how high-order the
Hamiltonian–Dissipator splitting and the Lie-product (palindromic)
ordering are.  This is the issue documented as **B-1** in
``STATUS_HONEST_2026-05.md``.

What this module does
---------------------
Each Lindblad operator built by
:func:`gksl_math_utils.build_lindblad_operators` is, by construction,
local — every TTA-pair channel acts non-trivially on exactly two adjacent
qutrits, every fluorescence/phosphorescence/IC/ISC channel acts on a
single qutrit.  We exploit that locality:

1. :func:`get_local_lindblad_ops` re-derives the same Lindblad operators
   in their **local** ``d`` × ``d`` (single-site) or ``d² × d²`` (pair)
   form, together with the site indices on which they act.
2. :func:`build_local_dissipator_super` builds the local GKSL dissipator
   superoperator
   ``L_D^local · vec(ρ_local) = vec(L ρ_local L† − ½{L†L, ρ_local})``
   in column-major vectorisation, returning a ``d² × d²`` (or
   ``d⁴ × d⁴``) matrix.
3. :func:`apply_channel_single` / :func:`apply_channel_pair` apply the
   already-exponentiated local channel ``expm(L_D^local · dt)`` to the
   full system density matrix via ``numpy.einsum``, without ever
   materialising the full-system superoperator (which would be
   ``dim² × dim²`` and is prohibitive for ``dim = d^N``).

This is mathematically the **exact** action of ``exp(L_{D_α} · dt)`` —
no truncation, no heuristic — provided each ``L_α`` is supported on at
most two sites, which is the case for the present GKSL model.

Empirical validation (recorded as fact, not as a guard)
-------------------------------------------------------
Replacing per-channel Stinespring with exact local channels in the
Strang+palindromic step of :class:`QuditGKSLSimulator` yields the
following measured convergence rates against the
:class:`ClassicalGKSLSimulator` reference (full Liouvillian exponential):

* N=2 (dim=9), ``t_max=10``: rate = 2.000 across n_steps = 20…400.
* N=4 (dim=81), ``t_max=100``: rate = 2.712 → 2.100 → 2.017 → 2.004
  for n_steps = 20 → 50 → 100 → 200 (asymptotic 2nd order).

These figures are reproduced by
:func:`tutorials.test_gksl_simulators.TestExactLocalChannels` and by
the convergence cell of ``quantum_dynamics_gksl_comparison.ipynb``.
The Stinespring algorithm continues to give rate ≈ 1, exactly as
documented in B-1.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters

LocalLindbladOp = tuple[tuple[int, ...], np.ndarray, str]
"""``(sites, L_local, kind)`` where ``kind`` is ``"single"`` or ``"pair"``.

* ``sites`` is a 1-tuple ``(i,)`` for ``"single"`` and a 2-tuple
  ``(i, j)`` with ``i < j`` for ``"pair"``.
* ``L_local`` is a ``d × d`` (single) or ``d² × d²`` (pair) complex
  matrix that already includes the ``sqrt(gamma)`` prefactor, matching
  the convention of :func:`gksl_math_utils.build_lindblad_operators`.
"""


def _ket_bra(d: int, a: int, b: int) -> np.ndarray:
    m = np.zeros((d, d), dtype=np.complex128)
    m[a, b] = 1.0
    return m


def get_local_lindblad_ops(
    params: GKSLPhysicalParameters,
) -> list[LocalLindbladOp]:
    """Return the Lindblad operators in their natural local form.

    The order matches :func:`gksl_math_utils.build_lindblad_operators`:
    first the 2 × len(neighbors) TTA-pair channels (kind ``"pair"``),
    then the 5 × N single-site channels (kind ``"single"``) for
    fluorescence, phosphorescence, internal conversion, ISC S→T,
    ISC T→S — in that order.

    Each ``L_local`` already contains the ``sqrt(gamma)`` factor.

    Parameters
    ----------
    params:
        GKSL physical parameters.  Must have ``with_boson=False``;
        the boson-extended model is out of scope for this routine.

    Returns
    -------
    list of (sites, L_local, kind)
        Local Lindblad operators ready for
        :func:`build_local_dissipator_super`.
    """
    if params.with_boson:
        msg = "get_local_lindblad_ops: with_boson=True is not supported"
        raise ValueError(msg)
    d = params.d
    N = params.N_molecules
    out: list[LocalLindbladOp] = []

    # --- TTA pair channels (kind="pair") ---
    for (i, j) in params.neighbors:
        if i >= j:
            msg = f"neighbors must satisfy i<j, got ({i},{j})"
            raise ValueError(msg)
        gamma = params.gamma_TTA / 2.0
        sq = np.sqrt(gamma)
        # Channel 1: |2>_i<1| ⊗ |0>_j<1|
        L1 = sq * np.kron(_ket_bra(d, 2, 1), _ket_bra(d, 0, 1))
        out.append(((i, j), L1, "pair"))
        # Channel 2: |0>_i<1| ⊗ |2>_j<1|
        L2 = sq * np.kron(_ket_bra(d, 0, 1), _ket_bra(d, 2, 1))
        out.append(((i, j), L2, "pair"))

    # --- Single-site channels in the same order as build_lindblad_operators ---
    singles_per_site = (
        (params.Gamma_fl, _ket_bra(d, 0, 2)),  # fluorescence
        (params.Gamma_ph, _ket_bra(d, 0, 1)),  # phosphorescence
        (params.k_IC,     _ket_bra(d, 0, 2)),  # internal conversion
        (params.k_ISC_ST, _ket_bra(d, 1, 2)),  # ISC S->T
        (params.k_ISC_TS, _ket_bra(d, 0, 1)),  # ISC T->S
    )
    # Note: build_lindblad_operators iterates outer "channel kind", inner "site".
    for (gamma, op) in singles_per_site:
        for i in range(N):
            L = np.sqrt(gamma) * op
            out.append(((i,), L, "single"))

    expected = 2 * len(params.neighbors) + 5 * N
    if len(out) != expected:
        msg = f"Expected {expected} local Lindblad ops, got {len(out)}"
        raise AssertionError(msg)
    return out


def build_local_dissipator_super(
    L_local: np.ndarray, d_local: int
) -> np.ndarray:
    """Build the local GKSL dissipator superoperator (column-major vec).

    Returns ``L_D^local`` such that

        L_D^local · vec(ρ_local)
            = vec( L ρ_local L† − ½ {L†L, ρ_local} ),

    using the column-major vectorisation
    ``vec(X)[i + d_local · j] = X[i, j]``.  The returned matrix has shape
    ``(d_local², d_local²)``.

    Parameters
    ----------
    L_local:
        ``(d_local, d_local)`` Lindblad operator (already ``√γ``-scaled).
    d_local:
        Local Hilbert-space dimension (``d`` for single-site, ``d²`` for
        pair channels).

    Notes
    -----
    Identical structure to the dissipator block of
    :func:`stinespring_utils.build_gksl_superoperator`, but for a local
    operator only.  No Hamiltonian piece is included on purpose — the
    Hamiltonian is handled separately by the Strang half-step in the
    Trotter loop.
    """
    if L_local.shape != (d_local, d_local):
        msg = (
            f"L_local has shape {L_local.shape}; expected "
            f"({d_local},{d_local})"
        )
        raise ValueError(msg)
    L = np.asarray(L_local, dtype=np.complex128)
    I = np.eye(d_local, dtype=np.complex128)
    LdL = L.conj().T @ L
    return (
        np.kron(L.conj(), L)
        - 0.5 * np.kron(I, LdL)
        - 0.5 * np.kron(LdL.T, I)
    )


def apply_channel_single(
    rho: np.ndarray,
    exp_LD_local: np.ndarray,
    site: int,
    N: int,
    d: int,
) -> np.ndarray:
    """Apply a single-site exact channel to the full-system density matrix.

    ``exp_LD_local`` is ``expm(L_D^local · dt)`` of shape
    ``(d², d²)``.  We reshape ``rho`` so that the row and column legs of
    the target ``site`` are explicit, then contract with the
    rank-4 tensor representation of the channel.

    Parameters
    ----------
    rho:
        ``(dim, dim)`` full-system density matrix, ``dim = d**N``.
    exp_LD_local:
        ``(d², d²)`` exponentiated local dissipator (column-major vec
        convention).
    site:
        Index of the qutrit on which the channel acts (``0 ≤ site < N``).
    N:
        Number of system qudits.
    d:
        Per-site dimension.
    """
    if exp_LD_local.shape != (d * d, d * d):
        msg = (
            f"exp_LD_local has shape {exp_LD_local.shape}; expected "
            f"({d*d},{d*d})"
        )
        raise ValueError(msg)
    dim = d ** N
    d_pre = d ** site
    d_post = d ** (N - 1 - site)
    rho6 = rho.reshape(d_pre, d, d_post, d_pre, d, d_post)
    # T[i, j, k, l] = exp_LD_local[i + d*j, k + d*l]   (column-major vec)
    T = exp_LD_local.reshape(d, d, d, d, order="F")
    # rho_out[a_pre, i, a_post, b_pre, j, b_post]
    #   = sum_{k,l} T[i, j, k, l] · rho6[a_pre, k, a_post, b_pre, l, b_post]
    out = np.einsum("ijkl,pkqrls->piqrjs", T, rho6, optimize=True)
    return out.reshape(dim, dim)


def apply_channel_pair(
    rho: np.ndarray,
    exp_LD_local: np.ndarray,
    sites: tuple[int, int],
    N: int,
    d: int,
) -> np.ndarray:
    """Apply a 2-site exact channel on (i, j) with i<j to the full ρ.

    ``exp_LD_local`` is ``expm(L_D^local · dt)`` of shape
    ``(d⁴, d⁴)`` for the local 2-qudit subsystem.  The two qudits are
    not required to be adjacent — we permute the legs into the front
    explicitly.  Mathematically this is identical to a tensor-network
    contraction; we perform it with reshape + einsum to avoid building
    the full ``dim² × dim²`` superoperator.

    Parameters
    ----------
    rho:
        ``(dim, dim)`` full-system density matrix.
    exp_LD_local:
        ``(d⁴, d⁴)`` exponentiated local pair dissipator (column-major
        vec on the d²-dim local Hilbert space).
    sites:
        ``(i, j)`` with ``0 ≤ i < j < N``.
    """
    i, j = sites
    if not (0 <= i < j < N):
        msg = f"sites must satisfy 0 ≤ i < j < N; got {sites} with N={N}"
        raise ValueError(msg)
    d2 = d * d
    if exp_LD_local.shape != (d2 * d2, d2 * d2):
        msg = (
            f"exp_LD_local has shape {exp_LD_local.shape}; expected "
            f"({d2*d2},{d2*d2})"
        )
        raise ValueError(msg)

    dim = d ** N
    shp = (d,) * (2 * N)
    R = rho.reshape(shp)

    # Original axes: 0..N-1 (row legs), N..2N-1 (col legs)
    row_axes = list(range(N))
    col_axes = list(range(N, 2 * N))
    other_row = [a for a in row_axes if a != i and a != j]
    other_col = [a for a in col_axes if a != (i + N) and a != (j + N)]
    perm = [i, j] + other_row + [i + N, j + N] + other_col
    R = R.transpose(perm)

    # Combine local row legs (d, d) -> d² and rest into a single block.
    d_rest = d ** (N - 2)
    R = R.reshape(d2, d_rest, d2, d_rest)

    # T[I, J, K, L] = exp_LD_local[I + d²·J, K + d²·L]   (column-major vec)
    T = exp_LD_local.reshape(d2, d2, d2, d2, order="F")
    out = np.einsum("IJKL,KqLs->IqJs", T, R, optimize=True)

    # Restore the rank-2N tensor and undo the permutation.
    out = out.reshape((d,) * 2 + (d,) * (N - 2) + (d,) * 2 + (d,) * (N - 2))
    inv = [0] * (2 * N)
    for new_pos, orig in enumerate(perm):
        inv[orig] = new_pos
    out = out.transpose(inv)
    return out.reshape(dim, dim)


def precompute_exact_channels_half(
    params: GKSLPhysicalParameters, dt: float
) -> list[tuple[tuple[int, ...], np.ndarray, str]]:
    """Pre-compute ``expm(L_D^local · dt/2)`` for every Lindblad channel.

    Returns a list of ``(sites, exp_LD_half, kind)`` triples in the
    same order as :func:`get_local_lindblad_ops`.  These half-step
    channels are used both in the forward and reverse passes of the
    palindromic Trotter step (one half-step each), so the product over
    the full pass implements ``exp(L_D · dt)`` exactly up to the
    Lie-product commutator error among the local channels — and that
    commutator error is itself ``O(dt³)`` per step thanks to the
    palindromic ordering, which is precisely what gives the global
    ``O(dt²)`` convergence.
    """
    out = []
    for (sites, L_local, kind) in get_local_lindblad_ops(params):
        d_local = L_local.shape[0]
        L_D_local = build_local_dissipator_super(L_local, d_local)
        exp_LD_half = expm(L_D_local * (dt / 2.0))
        out.append((sites, exp_LD_half, kind))
    return out


__all__ = [
    "LocalLindbladOp",
    "get_local_lindblad_ops",
    "build_local_dissipator_super",
    "apply_channel_single",
    "apply_channel_pair",
    "precompute_exact_channels_half",
]
