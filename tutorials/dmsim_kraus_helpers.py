"""Helpers for converting GKSL local-channel superoperators to Kraus form.

The :class:`tutorials.qudit_gksl_simulator.QuditGKSLSimulator` running with
``algorithm="exact_local_channels"`` already maintains a list of
exponentiated local dissipators
``M = expm(L_D^local · dt/2)`` for every Lindblad channel (single-site
``(d², d²)`` or pair ``(d⁴, d⁴)``) — see
:func:`tutorials.exact_local_channels.precompute_exact_channels_half`.

In order to feed those channels into a *backend execution* on the new
:class:`mqt.qudits.simulation.backends.DMSim` density-matrix backend
through MQT-Qudits :class:`~mqt.qudits.quantum_circuit.gates.KrausChannel`
instructions, we need to convert each local superoperator ``M`` into a
list of Kraus operators ``{K_α}`` such that
``M(ρ) = Σ_α K_α ρ K_α†``.

We do this exactly via the Choi-Jamiolkowski isomorphism: build the
Choi matrix of ``M``, eigen-decompose it (it is positive semidefinite
for a CPTP map), and read off ``K_α = sqrt(λ_α) · reshape(v_α, (d, d))``.
This is a textbook procedure — no heuristic, no truncation other than
discarding strictly non-positive eigenvalues whose magnitude is
indistinguishable from numerical noise (``< CHOI_EIG_TOL``).

Conventions
-----------
The local superoperators built by
:func:`tutorials.exact_local_channels.build_local_dissipator_super` use
**column-major vectorisation**: ``M[i + d·j, k + d·l]`` is the matrix
element such that ``M(ρ)[i, j] = Σ_{k,l} M[i + d·j, k + d·l] · ρ[k, l]``.

The Kraus operators returned by this module use the standard matrix
convention ``K[i, k]`` — row index = output, column index = input — so
they can be passed directly to
:meth:`QuantumCircuit.kraus_channel` and applied by DMSim with the
qudit ordering matching the original Lindblad operator's site ordering.
"""

from __future__ import annotations

import numpy as np

CHOI_EIG_TOL = 1e-12
"""Threshold below which Choi eigenvalues are treated as numerical noise.

Negative eigenvalues with magnitude above this tolerance indicate that
the input superoperator is not CPTP (e.g. an arithmetic bug); we raise
in that case rather than silently discarding them.
"""


def kraus_from_local_superoperator(
    m_local: np.ndarray,
    d_local: int,
    eig_tol: float = CHOI_EIG_TOL,
) -> list[np.ndarray]:
    """Return Kraus operators of the local CPTP map encoded by ``m_local``.

    Parameters
    ----------
    m_local:
        ``(d_local², d_local²)`` complex matrix in the column-major vec
        convention used by
        :func:`tutorials.exact_local_channels.build_local_dissipator_super`,
        already exponentiated (i.e. ``expm(L_D^local · dt)``).
    d_local:
        Local Hilbert-space dimension (``d`` for single-site channels,
        ``d²`` for pair channels).
    eig_tol:
        Tolerance for negative Choi eigenvalues.  Eigenvalues in
        ``(-eig_tol, eig_tol)`` are clamped to zero; eigenvalues
        ``< -eig_tol`` raise :class:`ValueError`.

    Returns
    -------
    list of (d_local, d_local) complex matrices
        Kraus operators ``{K_α}`` with ``Σ_α K_α ρ K_α† = M(ρ)`` for
        every density matrix ``ρ``.  Length equals the rank of the
        Choi matrix above ``eig_tol``.
    """
    d = int(d_local)
    if m_local.shape != (d * d, d * d):
        msg = (
            f"m_local has shape {m_local.shape}; expected "
            f"({d * d},{d * d}) for d_local={d}."
        )
        raise ValueError(msg)

    # Reshape M to T[i,j,k,l] = M[i + d·j, k + d·l]   (column-major vec).
    t_tensor = np.asarray(m_local, dtype=np.complex128).reshape(d, d, d, d, order="F")

    # Choi matrix: C[(i,k), (j,l)] = T[i,j,k,l].
    # Index mapping (row-major flattening): row = i·d + k, col = j·d + l.
    choi = t_tensor.transpose(0, 2, 1, 3).reshape(d * d, d * d)

    # Symmetrise to suppress numerical anti-Hermitian part (typically 1e-15).
    choi = 0.5 * (choi + choi.conj().T)

    eigvals, eigvecs = np.linalg.eigh(choi)

    if np.any(eigvals < -eig_tol):
        worst = float(np.min(eigvals.real))
        msg = (
            "Choi matrix has a strictly negative eigenvalue "
            f"{worst:.3e} below tol={eig_tol:.0e}; the input superoperator "
            "is not CPTP. No heuristic correction applied."
        )
        raise ValueError(msg)

    kraus_ops: list[np.ndarray] = []
    for alpha in range(eigvals.shape[0]):
        lam = eigvals[alpha].real
        if lam < eig_tol:
            continue
        v = eigvecs[:, alpha]
        # v[i·d + k] -> K[i, k] (row-major reshape).
        k_op = np.sqrt(lam) * v.reshape(d, d)
        kraus_ops.append(k_op.astype(np.complex128, copy=False))

    if not kraus_ops:
        # Degenerate case: zero map.  Return a single zero Kraus to keep
        # downstream code simple — but a CPTP map can never be the zero
        # map (it must preserve the trace), so flag this.
        msg = (
            "All Choi eigenvalues are below tol; the map appears to be "
            "the zero map, which is not CPTP."
        )
        raise ValueError(msg)

    return kraus_ops


def verify_kraus_channel_matches_superoperator(
    kraus_ops: list[np.ndarray],
    m_local: np.ndarray,
    d_local: int,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> float:
    """Return ``max |M(ρ) − Σ_α K_α ρ K_α†|`` over a fixed test basis.

    The test basis is the set of ``d_local²`` Hermitian matrices
    ``{|i⟩⟨j| + |j⟩⟨i|, i(|i⟩⟨j| − |j⟩⟨i|)}`` plus the diagonal projectors
    — sufficient to span the operator space.  Used as a numerical
    sanity check (not a heuristic correction).
    """
    d = int(d_local)
    t_tensor = m_local.reshape(d, d, d, d, order="F")

    max_err = 0.0
    for k in range(d):
        for ll in range(d):
            rho = np.zeros((d, d), dtype=np.complex128)
            rho[k, ll] = 1.0
            ref = t_tensor[:, :, k, ll].copy()
            kraus_action = sum(
                k_op @ rho @ k_op.conj().T for k_op in kraus_ops
            )
            err = float(np.max(np.abs(kraus_action - ref)))
            max_err = max(max_err, err)
    if max_err > atol + rtol * max(1.0, float(np.max(np.abs(m_local)))):
        msg = (
            f"Kraus reconstruction failed: max element error {max_err:.3e} "
            f"exceeds atol={atol:g} + rtol*||M||={rtol * float(np.max(np.abs(m_local))):.3e}"
        )
        raise AssertionError(msg)
    return max_err


__all__ = [
    "CHOI_EIG_TOL",
    "kraus_from_local_superoperator",
    "verify_kraus_channel_matches_superoperator",
]
