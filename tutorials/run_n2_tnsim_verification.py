"""A-1: end-to-end MQT-Qudits ``tnsim`` verification on the N=2 minimal config.

Why this script exists
----------------------
PR #257 documents in ``STATUS_HONEST_2026-05.md`` (item A-1) that the
notebook ``quantum_dynamics_gksl_comparison.ipynb`` and every
``QuditGKSL*Simulator`` perform their actual time-evolution arithmetic
with NumPy / ``scipy.linalg.expm`` on a density matrix.  The MQT-Qudits
``QuantumCircuit`` API is used to **construct** circuits, but those
circuits are never executed on a MQT-Qudits backend
(``tnsim`` / ``misim``) — partly because the *combined* per-step
circuit needs mid-circuit ancilla reset (which the backends do not
implement).

This script closes the gap that *can* honestly be closed: each
**individual building block** of one Trotter step (the Hamiltonian
half-step plus every Stinespring sub-circuit) **is** a self-contained
circuit that can be executed end-to-end on the ``tnsim`` backend.
We do exactly that for the smallest meaningful configuration
(``N_molecules=2``, ``d=3``, dim ≤ 27 including the ancilla) and
verify that the resulting state vector matches the matrix-vector
product produced by the in-Python gate matrices to floating-point
precision.

What is verified (per Trotter step, ``N=2``, ``d=3``)
-----------------------------------------------------
1. **Hamiltonian half-step circuit**: ``cu_one`` on each of 2 qutrits
   plus ``cu_two`` on the (0, 1) pair → state vector under ``tnsim``
   matches ``U_circuit · |0...0⟩`` from the in-process unitary.

2. **Single-site Stinespring sub-circuit** (``cu_two`` on
   target+ancilla, both d=3 → 9-dim local) — done for every
   single-site Lindblad channel
   (5 × N = 10 channels for N=2): ``tnsim`` state vector matches
   ``expm(-i √dt · G_α) · |00⟩`` exactly (||Δ|| ≤ 1e-10).

3. **TTA-pair Stinespring sub-circuit** (``cu_multi`` on
   pair+ancilla, three d=3 qutrits → 27-dim local) — done for every
   pair Lindblad channel (2 × |neighbors| = 2 channels for N=2):
   same exact match against the in-process unitary.

What is **NOT** verified (and why)
----------------------------------
The combined per-step circuit returned by
``QuditGKSLKrausSimulator.build_combined_trotter_step_circuit`` is
*not* executed.  It contains ancilla reuse: each Stinespring channel
needs its ancilla reset to ``|0⟩`` before the next channel.
MQT-Qudits' ``tnsim`` / ``misim`` do not currently implement
mid-circuit reset, so running the combined circuit would either error
out or silently produce wrong results.  This is recorded as the
remaining piece of A-1 in ``STATUS_HONEST_2026-05.md``.

Exit status
-----------
``0`` on full success (all matches within 1e-10).
``1`` on any failure; the failing block is printed.

The script is also picked up by ``pytest`` via the
``test_n2_tnsim_verification`` function below.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_circuit_simulator import QuditGKSLKrausSimulator

TOL = 1e-10


def _run_subcircuit_on_tnsim(circuit, U_local: np.ndarray, label: str) -> dict:
    """Run ``circuit`` on ``tnsim`` and compare to ``U_local · |0...0⟩``.

    Returns a result dict with ``label``, ``distance``, ``match``.
    Raises ``ImportError`` if MQT-Qudits is not installed.
    """
    from mqt.qudits.simulation import MQTQuditProvider

    provider = MQTQuditProvider()
    backend = provider.get_backend("tnsim")
    job = backend.run(circuit, shots=1)
    sv_tnsim = job.result().get_state_vector()[0]

    psi0 = np.zeros(U_local.shape[0], dtype=np.complex128)
    psi0[0] = 1.0
    sv_expected = U_local @ psi0

    distance = float(np.linalg.norm(sv_tnsim - sv_expected))
    return {
        "label": label,
        "distance": distance,
        "match": distance < TOL,
    }


def run_n2_tnsim_verification(dt: float = 0.5) -> list[dict]:
    """Run the full N=2 sub-circuit verification.

    Returns a list of result dicts, one per verified sub-circuit, in
    the order: Hamiltonian half-step, all single-site Stinespring,
    all TTA-pair Stinespring.
    """
    params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
    sim = QuditGKSLKrausSimulator(params)
    results: list[dict] = []

    # --- 1. Hamiltonian half-step circuit ---
    ham_check = sim.verify_circuit_via_mqt(dt=dt)
    results.append({
        "label": f"hamiltonian_half_step(dt={dt})",
        "distance": ham_check["statevector_distance"],
        "match": ham_check["match"],
        "gate_count": ham_check["gate_count"],
    })

    # --- 2. Each single-site Stinespring sub-circuit ---
    single_idx = 0
    pair_idx = 0
    for op_type, sites, L_local, _gamma in sim.lindblad_local_info:
        if op_type == "single":
            circ, U = sim.build_stinespring_circuit_single(
                L_local, dt, sites[0]
            )
            label = f"stinespring_single[#{single_idx}, site={sites[0]}]"
            results.append(_run_subcircuit_on_tnsim(circ, U, label))
            single_idx += 1
        elif op_type == "pair":
            circ, U = sim.build_stinespring_circuit_pair(
                L_local, dt, sites[0], sites[1]
            )
            label = (
                f"stinespring_pair[#{pair_idx}, "
                f"sites=({sites[0]},{sites[1]})]"
            )
            results.append(_run_subcircuit_on_tnsim(circ, U, label))
            pair_idx += 1
        else:
            msg = f"Unexpected op_type: {op_type!r}"
            raise ValueError(msg)

    return results


def test_n2_tnsim_verification() -> None:
    """pytest entry point: every sub-circuit must match within 1e-10."""
    results = run_n2_tnsim_verification(dt=0.5)
    failures = [r for r in results if not r["match"]]
    if failures:
        msg_lines = ["A-1 sub-circuit verification failed:"]
        for r in failures:
            msg_lines.append(
                f"  {r['label']}: ||Δsv|| = {r['distance']:.3e} (tol {TOL:g})"
            )
        raise AssertionError("\n".join(msg_lines))


def _main() -> int:
    print("=" * 72)
    print("A-1: N=2 sub-circuit verification on MQT-Qudits tnsim backend")
    print("=" * 72)
    print()
    print("Combined per-step circuit is NOT executed (mid-circuit ancilla")
    print("reset is not implemented by tnsim / misim).  Each individual")
    print("building block is executed on tnsim and compared against the")
    print("in-process unitary applied to |0...0⟩.")
    print()
    try:
        results = run_n2_tnsim_verification(dt=0.5)
    except ImportError as e:
        print(f"ERROR: MQT-Qudits not available: {e}")
        return 1

    n_total = len(results)
    n_match = sum(1 for r in results if r["match"])
    print(f"Verified {n_total} sub-circuits (tol = {TOL:g}):")
    print(f"  {'label':<55} {'||Δsv||':>10}  match")
    print("  " + "-" * 75)
    for r in results:
        print(
            f"  {r['label']:<55} {r['distance']:>10.2e}  "
            f"{'OK' if r['match'] else 'FAIL'}"
        )
    print()
    print(f"Result: {n_match}/{n_total} sub-circuits match (tol = {TOL:g})")
    if n_match != n_total:
        print("FAILURE")
        return 1
    print("SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
