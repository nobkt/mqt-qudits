"""A-2: cross-validate QuditGKSL simulators against an independent QuTiP reference.

This script runs :class:`QutipGKSLReference` (which uses ``qutip.mesolve``,
sharing **no code** with the in-tree GKSL stack) against
:class:`QuditGKSLSimulator` for the two algorithms it supports
(``"exact_local_channels"`` and ``"stinespring"``), at ``N = 2`` and
``N = 4``, and prints the **measured** trace distances and
Frobenius-norm differences.

Why: ``STATUS_HONEST_2026-05.md`` records under A-2 that the existing
qubit/qudit simulators share their entire core (Hamiltonian / collapse-op
builders, Stinespring dilation, Strang+palindromic Trotter loop), so
their internal agreement only validates the qutrit→qubit embedding.
The QuTiP reference is genuinely independent.

Exit status: ``0`` on success (all final-state agreements better than the
asserted tolerances).  ``1`` on any tolerance violation.

The script is also picked up by ``pytest`` via
``test_qutip_independent_cross_validation`` below.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_simulator import QuditGKSLSimulator
from qutip_gksl_reference import QutipGKSLReference


# Tolerances are calibrated against the *measured* numbers reported in
# the script output below; they are NOT padded heuristics.  See the
# block comments next to each assertion in
# :func:`test_qutip_independent_cross_validation` for the actual
# observed margins.
TOL_EXACT_N2 = 1e-7   # measured ~2e-8 (n_steps=200, t_max=10)
TOL_EXACT_N4 = 5e-6   # measured ~8e-7 (n_steps=200, t_max=10)


def _trace_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Tr|a - b| / 2 via singular values (numerical-grade implementation)."""
    diff = a - b
    # 1/2 * sum of singular values = trace distance for Hermitian Δ
    s = np.linalg.svd(diff, compute_uv=False)
    return float(0.5 * np.sum(s))


def cross_validate(N: int, t_max: float, n_steps: int) -> dict:
    """Run QuTiP reference and both QuditGKSL algorithms; return diagnostics.

    Returns a dict with the measured Frobenius and trace distances and
    the final density matrices for visual inspection.
    """
    params = GKSLPhysicalParameters(N_molecules=N, with_boson=False)
    ref = QutipGKSLReference(params)
    out_qutip = ref.simulate(t_max=t_max, n_times=11)

    sim_exact = QuditGKSLSimulator(params, algorithm="exact_local_channels")
    res_exact = sim_exact.simulate(t_max=t_max, n_steps=n_steps)

    sim_stine = QuditGKSLSimulator(params, algorithm="stinespring")
    res_stine = sim_stine.simulate(t_max=t_max, n_steps=n_steps)

    rho_qutip = out_qutip["rho_final"]
    rho_exact = res_exact["rho_final"]
    rho_stine = res_stine["rho_final"]

    return {
        "N": N,
        "t_max": t_max,
        "n_steps": n_steps,
        "frob_qutip_vs_exact":  float(np.linalg.norm(rho_qutip - rho_exact, ord="fro")),
        "frob_qutip_vs_stine":  float(np.linalg.norm(rho_qutip - rho_stine, ord="fro")),
        "trd_qutip_vs_exact":   _trace_distance(rho_qutip, rho_exact),
        "trd_qutip_vs_stine":   _trace_distance(rho_qutip, rho_stine),
        "tr_qutip":  float(np.real(np.trace(rho_qutip))),
        "tr_exact":  float(np.real(np.trace(rho_exact))),
        "tr_stine":  float(np.real(np.trace(rho_stine))),
    }


def test_qutip_independent_cross_validation() -> None:
    """pytest entry point: independent QuTiP reference must agree with `exact_local_channels`.

    Tolerances reflect measured behaviour with a tight ODE (atol=rtol=1e-12)
    and ``n_steps=200`` for the in-tree Trotter loop:

      * N=2: trace distance ~1.8e-8 → asserted < 1e-7.
      * N=4: trace distance ~7.7e-7 → asserted < 5e-6.

    The Stinespring algorithm intentionally is **not** asserted to the
    same tolerance — it is a 1st-order channel approximation and its
    finite trace-distance to the QuTiP reference is precisely what
    B-1 documents.  We only print its diagnostic for context.
    """
    diag2 = cross_validate(N=2, t_max=10.0, n_steps=200)
    diag4 = cross_validate(N=4, t_max=10.0, n_steps=200)

    # Trace preservation in the QuTiP reference
    assert abs(diag2["tr_qutip"] - 1.0) < 1e-9, diag2
    assert abs(diag4["tr_qutip"] - 1.0) < 1e-9, diag4

    # Independent agreement with `exact_local_channels`
    assert diag2["trd_qutip_vs_exact"] < TOL_EXACT_N2, (
        f"N=2 trace distance {diag2['trd_qutip_vs_exact']:.3e} "
        f"exceeds tol {TOL_EXACT_N2}"
    )
    assert diag4["trd_qutip_vs_exact"] < TOL_EXACT_N4, (
        f"N=4 trace distance {diag4['trd_qutip_vs_exact']:.3e} "
        f"exceeds tol {TOL_EXACT_N4}"
    )


def _print_table(diag: dict) -> None:
    print(
        f"  N={diag['N']:>2d}  n_steps={diag['n_steps']:>3d}  t_max={diag['t_max']}"
    )
    print(
        f"     Tr|ρ_qutip - ρ_exact_local|/2 = {diag['trd_qutip_vs_exact']:.3e}"
        f"   ||·||_F = {diag['frob_qutip_vs_exact']:.3e}"
    )
    print(
        f"     Tr|ρ_qutip - ρ_stinespring|/2 = {diag['trd_qutip_vs_stine']:.3e}"
        f"   ||·||_F = {diag['frob_qutip_vs_stine']:.3e}"
    )
    print(
        f"     Tr(ρ_qutip)={diag['tr_qutip']:.12f}"
        f"   Tr(ρ_exact)={diag['tr_exact']:.12f}"
        f"   Tr(ρ_stine)={diag['tr_stine']:.12f}"
    )


def _main() -> int:
    print("=" * 76)
    print("A-2: independent QuTiP cross-validation of QuditGKSLSimulator")
    print("=" * 76)
    print()
    print("QuTiP reference uses qutip.mesolve (adaptive Adams/BDF ODE), shares")
    print("no code with the in-tree GKSL stack.  Same physical Hamiltonian and")
    print("the same 5×N + 2×|neighbors| Lindblad collapse operators are built")
    print("independently from QuTiP primitives.")
    print()

    overall_ok = True
    for N in (2, 4):
        diag = cross_validate(N=N, t_max=10.0, n_steps=200)
        _print_table(diag)
        tol = TOL_EXACT_N2 if N == 2 else TOL_EXACT_N4
        ok = diag["trd_qutip_vs_exact"] < tol
        overall_ok = overall_ok and ok
        status = "OK" if ok else "FAIL"
        print(f"     -> {status} (asserted Tr|Δρ|/2 < {tol:g})")
        print()

    if not overall_ok:
        print("FAILURE")
        return 1
    print("SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
