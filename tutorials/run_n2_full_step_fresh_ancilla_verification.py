"""A-1 残課題: backend-executable per-step circuit via fresh ancillas.

This script targets the residual A-1 item recorded in
``STATUS_HONEST_2026-05.md``:

    The combined per-step circuit produced by
    ``QuditGKSLKrausSimulator.build_combined_trotter_step_circuit``
    re-uses one ancilla across every Stinespring channel and therefore
    requires mid-circuit ancilla reset, which ``tnsim`` / ``misim`` do
    not implement.  Hence the combined circuit is not executable on a
    MQT-Qudits backend in its current form.

Two independent honest contributions are reported below:

A. Mathematical equivalence (always runs, no backend involvement).
   We construct an alternative per-step circuit using
   :meth:`QuditGKSLKrausSimulator.build_executable_per_step_circuit_fresh_ancillas`,
   which allocates a **fresh ancilla per Stinespring channel** so that
   no reset is needed.  We then evolve the joint pure state in NumPy,
   trace out the ancillas, and compare the resulting reduced density
   matrix on the system register against the in-tree Kraus simulator
   applied to the same initial system state with the same gate
   sequence.  The two ρ_sys must agree to floating-point precision.

B. Backend feasibility probe (always runs, may report negative results
   honestly).  We attempt to execute progressively larger fresh-ancilla
   circuits on ``tnsim`` for ``N=2``:

     - 1 single-site Stinespring channel  ( 3 qudits)
     - 2 single-site Stinespring channels ( 4 qudits)
     - 5 single-site Stinespring channels ( 7 qudits)  -- "all single-site
       channels on site 0"
     - 10 single-site Stinespring channels (12 qudits) -- "all single-site
       channels on both sites"
     - the full forward-only per-step circuit (14 qudits, 12 channels)

   For each we record whether ``tnsim`` succeeds (and in that case
   compare the traced-out ρ_sys against the in-process Kraus result),
   or fails with a backend-side memory error and exactly how much
   memory the backend tried to allocate.  The reason the backend
   ultimately can't run the **full** N=2 per-step circuit is **not**
   "missing reset" — it is that ``tnsim`` expands every long-range or
   multi-qudit gate to a contiguous ``[min(qudits) … max(qudits)]``
   intermediate matrix (see ``src/mqt/qudits/simulation/backends/tnsim.py``,
   lines 110-118), and that intermediate matrix exceeds available
   memory once the qudit range crosses ~9 qutrits.  This is a
   different, separately documented obstacle.

Exit status: ``0`` on success of part (A) and on either success or
"honestly reported memory failure" of part (B).  ``1`` only if part
(A) fails (the mathematical equivalence is non-negotiable) or if part
(B) exhibits a *correctness* error (mismatch in traced ρ_sys).
"""

from __future__ import annotations

import gc
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_circuit_simulator import QuditGKSLKrausSimulator


TOL = 1e-10


# ---------------------------------------------------------------------------
# Part A: mathematical equivalence (in-process pure-state ↔ Kraus on system)
# ---------------------------------------------------------------------------


def _kraus_forward_step_on_system(
    sim: QuditGKSLKrausSimulator,
    rho_sys: np.ndarray,
    dt: float,
) -> np.ndarray:
    """Apply *one* forward-only step in pure NumPy density-matrix arithmetic.

    The exact analogue of running the circuit returned by
    ``build_executable_per_step_circuit_fresh_ancillas(dt, palindromic=False)``
    starting from ``ρ_sys ⊗ |0..0⟩⟨0..0|_anc`` and tracing out the
    ancillas:

        ρ_sys'  =  Π_α  E_α(dt)  ∘  U_H(dt)  [ ρ_sys ]

    where ``U_H(dt)`` is the same Trotter-decomposed Hamiltonian
    unitary the circuit applies (on-site tensor product, then sequential
    pair unitaries), and ``E_α`` is the Kraus channel obtained from the
    per-channel Stinespring unitary.
    """
    rho = sim._apply_hamiltonian_step_circuit(rho_sys, dt)
    for (op_type, sites, L_local, _gamma) in sim.lindblad_local_info:
        U_local = sim._build_local_stinespring_unitary(L_local, dt)
        d_local = L_local.shape[0]
        kraus = sim._extract_kraus_from_local_stinespring(U_local, d_local)
        if op_type == "single":
            rho = sim._apply_single_site_channel(rho, kraus, sites[0])
        elif op_type == "pair":
            rho = sim._apply_pair_channel(rho, kraus, sites[0], sites[1])
        else:
            msg = f"Unexpected op_type: {op_type!r}"
            raise ValueError(msg)
    return rho


def part_a_math_equivalence(dt: float = 0.5) -> dict:
    """Verify fresh-ancilla pure-state evolution + partial trace == Kraus map."""
    params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
    sim = QuditGKSLKrausSimulator(params)
    info = sim.build_executable_per_step_circuit_fresh_ancillas(dt=dt, palindromic=False)

    n_total = info["n_total_qudits"]
    d = sim.d
    N = sim.N

    # In-process pure-state evolution from |0...0⟩
    psi0 = np.zeros(d ** n_total, dtype=np.complex128)
    psi0[0] = 1.0
    psi = sim.evolve_pure_state_in_process(
        psi0, info["gate_sequence"], info["local_unitaries"], n_total
    )
    rho_sys_from_circuit = sim.reduced_density_matrix_on_system(psi, N, n_total)

    # Reference: Kraus channel composition on |0..0⟩⟨0..0|_sys
    rho_sys0 = np.zeros((d ** N, d ** N), dtype=np.complex128)
    rho_sys0[0, 0] = 1.0
    rho_sys_ref = _kraus_forward_step_on_system(sim, rho_sys0, dt)

    diff = np.linalg.norm(rho_sys_from_circuit - rho_sys_ref, ord="fro")
    tr_circuit = float(np.real(np.trace(rho_sys_from_circuit)))
    tr_ref = float(np.real(np.trace(rho_sys_ref)))

    return {
        "dt": dt,
        "n_total_qudits": n_total,
        "n_ancillas": info["n_ancilla_qudits"],
        "n_gates": len(info["gate_sequence"]),
        "frob_diff": float(diff),
        "tr_circuit": tr_circuit,
        "tr_ref": tr_ref,
        "match": diff < TOL,
    }


# ---------------------------------------------------------------------------
# Part B: backend feasibility probe on tnsim with subsets of channels
# ---------------------------------------------------------------------------


def _build_subset_circuit(
    sim: QuditGKSLKrausSimulator, dt: float, subset: list[int]
) -> dict:
    """Build a fresh-ancilla circuit using only the Lindblad channels at indices in ``subset``.

    Hamiltonian half-step is included.  ``subset`` selects from
    ``sim.lindblad_local_info`` (forward-only sweep, no reverse).

    The returned dict mirrors the keys of
    :meth:`QuditGKSLKrausSimulator.build_executable_per_step_circuit_fresh_ancillas`.
    """
    from mqt.qudits.quantum_circuit import QuantumCircuit
    from scipy.linalg import expm

    N = sim.N
    d = sim.d
    d_anc = sim.d_anc
    n_anc = len(subset)
    n_total = N + n_anc

    dims = [d] * N + [d_anc] * n_anc
    circuit = QuantumCircuit(n_total, dims, 0)

    gate_sequence: list[dict] = []
    local_unitaries: list[np.ndarray] = []

    # Hamiltonian (full dt because forward-only)
    U_onsite = expm(-1j * sim.h_local * dt)
    for i in range(N):
        circuit.cu_one(i, U_onsite)
        gate_sequence.append({"kind": "cu_one", "qudits": [i], "tag": "H_onsite"})
        local_unitaries.append(U_onsite)
    for (ip, jp) in sim.params.neighbors:
        U_pair = expm(-1j * sim.h_transfer_pairs[(ip, jp)] * dt)
        circuit.cu_two([ip, jp], U_pair)
        gate_sequence.append({"kind": "cu_two", "qudits": [ip, jp], "tag": "H_pair"})
        local_unitaries.append(U_pair)

    # Selected Stinespring channels
    next_anc = N
    for k in subset:
        op_type, sites, L_local, _gamma = sim.lindblad_local_info[k]
        U_local = sim._build_local_stinespring_unitary(L_local, dt)
        anc = next_anc
        next_anc += 1
        if op_type == "single":
            circuit.cu_two([sites[0], anc], U_local)
            gate_sequence.append({"kind": "cu_two", "qudits": [sites[0], anc], "tag": "stinespring_single_fwd"})
        elif op_type == "pair":
            circuit.cu_multi([sites[0], sites[1], anc], U_local)
            gate_sequence.append({"kind": "cu_multi", "qudits": [sites[0], sites[1], anc], "tag": "stinespring_pair_fwd"})
        else:
            msg = f"Unexpected op_type: {op_type!r}"
            raise ValueError(msg)
        local_unitaries.append(U_local)

    return {
        "circuit": circuit,
        "n_total_qudits": n_total,
        "n_ancilla_qudits": n_anc,
        "n_system_qudits": N,
        "gate_sequence": gate_sequence,
        "local_unitaries": local_unitaries,
    }


def _kraus_subset_step_on_system(
    sim: QuditGKSLKrausSimulator,
    rho_sys: np.ndarray,
    dt: float,
    subset: list[int],
) -> np.ndarray:
    """In-process Kraus forward-only step using only the channels in ``subset``."""
    rho = sim._apply_hamiltonian_step_circuit(rho_sys, dt)
    for k in subset:
        op_type, sites, L_local, _gamma = sim.lindblad_local_info[k]
        U_local = sim._build_local_stinespring_unitary(L_local, dt)
        d_local = L_local.shape[0]
        kraus = sim._extract_kraus_from_local_stinespring(U_local, d_local)
        if op_type == "single":
            rho = sim._apply_single_site_channel(rho, kraus, sites[0])
        elif op_type == "pair":
            rho = sim._apply_pair_channel(rho, kraus, sites[0], sites[1])
    return rho


def _try_tnsim(circuit) -> tuple[bool, str, np.ndarray | None]:
    """Run ``circuit`` on ``tnsim``; return (success, message, statevector_or_None)."""
    try:
        from mqt.qudits.simulation import MQTQuditProvider

        prov = MQTQuditProvider()
        be = prov.get_backend("tnsim")
        sv = np.asarray(be.run(circuit, shots=1).result().get_state_vector()[0])
    except MemoryError as e:
        return False, f"MemoryError: {e}", None
    except Exception as e:  # noqa: BLE001
        # Numpy raises numpy._core._exceptions._ArrayMemoryError, which is a
        # subclass of MemoryError, but some platforms surface it as plain
        # Exception.  We catch broadly here so that the probe always
        # produces a structured report rather than aborting.
        return False, f"{type(e).__name__}: {e}", None
    return True, "ok", sv


def part_b_backend_probe(dt: float = 0.5) -> list[dict]:
    """Probe what subsets of channels tnsim can execute end-to-end."""
    params = GKSLPhysicalParameters(N_molecules=2, with_boson=False)
    sim = QuditGKSLKrausSimulator(params)

    # Inspect the channel order so the chosen subsets are interpretable.
    # `lindblad_local_info` order:
    #   * pair channels first  (2 of them for N=2),
    #   * then single-site channels (5 kinds × 2 sites = 10), in
    #     kind-major / site-minor order.
    # We probe progressively: start from the lightest.
    n_lindblad = len(sim.lindblad_local_info)
    if n_lindblad != 12:
        msg = f"Expected 12 Lindblad channels for N=2, got {n_lindblad}"
        raise AssertionError(msg)

    # Indices of single-site channels (skip the 2 pair channels at the front)
    single_channel_indices = [
        k for k, info in enumerate(sim.lindblad_local_info) if info[0] == "single"
    ]

    probes: list[tuple[str, list[int]]] = [
        ("1× single-site channel",         single_channel_indices[:1]),
        ("2× single-site channels",        single_channel_indices[:2]),
        ("5× single-site channels (1 site)", single_channel_indices[:5]),
        ("10× single-site channels (both sites)", single_channel_indices),
        ("full forward-only per-step (12 channels)", list(range(n_lindblad))),
    ]

    results: list[dict] = []
    for label, subset in probes:
        info = _build_subset_circuit(sim, dt, subset)
        n_total = info["n_total_qudits"]
        # The intermediate-matrix size that tnsim allocates for the
        # *worst* gate in this circuit (worst = max over gates of
        # max(qudits)-min(qudits)+1).
        max_range = 0
        for g in info["gate_sequence"]:
            qd = g["qudits"]
            if len(qd) >= 2:
                rng = max(qd) - min(qd) + 1
                if rng > max_range:
                    max_range = rng
        worst_intermediate_size = sim.d ** max_range
        worst_intermediate_bytes = (worst_intermediate_size ** 2) * 16  # complex128

        ok, msg, sv = _try_tnsim(info["circuit"])
        gc.collect()

        record = {
            "label": label,
            "n_channels": len(subset),
            "n_total_qudits": n_total,
            "worst_qudit_range_in_gate": max_range,
            "worst_intermediate_matrix_bytes": int(worst_intermediate_bytes),
            "tnsim_ok": ok,
            "tnsim_msg": msg,
        }

        if ok:
            # Verify ρ_sys via partial trace == in-process Kraus result.
            rho_sys_circ = sim.reduced_density_matrix_on_system(
                sv, sim.N, n_total
            )
            rho_sys0 = np.zeros((sim.d ** sim.N, sim.d ** sim.N), dtype=np.complex128)
            rho_sys0[0, 0] = 1.0
            rho_sys_ref = _kraus_subset_step_on_system(sim, rho_sys0, dt, subset)
            diff = float(np.linalg.norm(rho_sys_circ - rho_sys_ref, ord="fro"))
            record["frob_diff_vs_kraus_ref"] = diff
            record["correctness_match"] = diff < TOL
        results.append(record)
    return results


def test_a1_fresh_ancilla_circuit() -> None:
    """pytest entry: math equivalence MUST pass; backend probe records facts."""
    a = part_a_math_equivalence(dt=0.5)
    assert a["match"], (
        f"Part A failed: ||Δρ_sys||_F = {a['frob_diff']:.3e} (tol {TOL:g})\n"
        f"  tr_circuit = {a['tr_circuit']}, tr_ref = {a['tr_ref']}"
    )

    probes = part_b_backend_probe(dt=0.5)
    # For every tnsim execution that *did* succeed, the traced ρ_sys
    # must match the Kraus reference (this is a correctness invariant —
    # tnsim is allowed to fail with a memory error, but it is not
    # allowed to produce a wrong density matrix).
    for r in probes:
        if r["tnsim_ok"]:
            assert r["correctness_match"], (
                f"tnsim ran on {r['label']} but ρ_sys mismatch: "
                f"||Δρ||_F = {r['frob_diff_vs_kraus_ref']:.3e}"
            )
    # We also require at least *some* non-trivial channel subset to
    # successfully execute on tnsim — otherwise the residual A-1 result
    # would be empty.  Measured (see _main below): the 1, 2, 5 channel
    # subsets all succeed.
    n_ok = sum(1 for r in probes if r["tnsim_ok"])
    assert n_ok >= 3, (
        f"Expected at least 3 tnsim successes among the probes, got {n_ok}.\n"
        + "\n".join(
            f"  {r['label']}: {'OK' if r['tnsim_ok'] else 'FAIL'} ({r['tnsim_msg']})"
            for r in probes
        )
    )


def _main() -> int:
    print("=" * 78)
    print("A-1 残課題: backend-executable per-step circuit via fresh ancillas")
    print("=" * 78)
    print()
    print("Part A — Mathematical equivalence")
    print("---------------------------------")
    a = part_a_math_equivalence(dt=0.5)
    print(
        f"  forward-only per-step circuit: {a['n_total_qudits']} qudits "
        f"({a['n_ancillas']} fresh ancillas), {a['n_gates']} gates"
    )
    print(
        f"  ||ρ_sys(circuit) - ρ_sys(Kraus)||_F = {a['frob_diff']:.3e}  "
        f"(tol {TOL:g})  -> {'OK' if a['match'] else 'FAIL'}"
    )
    print(
        f"  Tr(ρ_sys_circuit) = {a['tr_circuit']:.12f}, "
        f"Tr(ρ_sys_ref) = {a['tr_ref']:.12f}"
    )
    print()
    print("Part B — Backend feasibility probe on tnsim")
    print("-------------------------------------------")
    probes = part_b_backend_probe(dt=0.5)
    print(
        f"  {'subset':<48}{'qudits':>7}"
        f"{'gate_range':>12}{'inter_mat':>14}{'tnsim':>8}"
    )
    print("  " + "-" * 89)
    for r in probes:
        size_h = (
            f"{r['worst_intermediate_matrix_bytes']/1024/1024/1024:.1f}GB"
            if r["worst_intermediate_matrix_bytes"] >= 1024 ** 3
            else f"{r['worst_intermediate_matrix_bytes']/1024/1024:.1f}MB"
        )
        status = "OK" if r["tnsim_ok"] else "FAIL"
        print(
            f"  {r['label']:<48}{r['n_total_qudits']:>7d}"
            f"{r['worst_qudit_range_in_gate']:>12d}{size_h:>14}{status:>8}"
        )
        if r["tnsim_ok"]:
            print(
                f"    -> ||ρ_sys(circuit) - ρ_sys(Kraus)||_F = "
                f"{r['frob_diff_vs_kraus_ref']:.3e}"
            )
        else:
            print(f"    -> tnsim msg: {r['tnsim_msg'][:120]}")
    print()
    if not a["match"]:
        print("FAILURE (Part A correctness)")
        return 1
    n_ok = sum(1 for r in probes if r["tnsim_ok"])
    if n_ok < 3:
        print(f"FAILURE (only {n_ok} tnsim successes; expected ≥3)")
        return 1
    print(f"SUCCESS (Part A: math equivalence OK; Part B: {n_ok}/{len(probes)} probes ran on tnsim)")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
