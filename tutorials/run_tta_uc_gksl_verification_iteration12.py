#!/usr/bin/env python3
"""TTA-UC GKSL iteration 12 verification: noise channel mathematical verification.

This script provides independent numerical verification that the stochastic
noise implementations in the shot-based simulators correctly reproduce the
corresponding density-matrix-level noise channels.

Verification tests:
  1. Single-step noise channel: stochastic average vs analytical density matrix
  2. Per-step fidelity decay tracking
  3. Error budget analysis (noise gate count × parameters)
  4. Qubit vs Qudit fair comparison (normalized fidelity)

Results are written to developing/verification_results/ as JSON and Markdown.
"""

from __future__ import annotations

import json
import os
import sys
import time as time_module
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_math_utils import (
    build_lindblad_operators,
    build_onsite_hamiltonian,
    build_transfer_hamiltonian,
    compute_von_neumann_entropy,
)
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_shot_simulator import QubitGKSLNoisyShotSimulator, QubitGKSLShotSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator, QuditGKSLShotSimulator
from qudit_gksl_simulator import QuditGKSLSimulator
from stinespring_utils import apply_stinespring_to_density_matrix, stinespring_unitary_from_lindblad


def quantum_fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Compute quantum state fidelity F(rho, sigma) via eigendecomposition."""
    evals, evecs = np.linalg.eigh(rho)
    evals = np.maximum(evals, 0.0)
    sqrt_rho = evecs @ np.diag(np.sqrt(evals)) @ evecs.conj().T
    m = sqrt_rho @ sigma @ sqrt_rho
    m = (m + m.conj().T) / 2
    evals_m = np.linalg.eigvalsh(m)
    evals_m = np.maximum(evals_m, 0.0)
    return float(min(np.real(np.sum(np.sqrt(evals_m))) ** 2, 1.0))


def normalize_rho(rho: np.ndarray) -> np.ndarray:
    """Normalize density matrix to trace 1."""
    tr = float(np.real(np.trace(rho)))
    return rho / tr if tr > 1e-10 else rho


def apply_depolarization_channel_dm(
    rho: np.ndarray, site: int, d: int, N: int, p: float
) -> np.ndarray:
    """Apply single-site depolarization channel at density-matrix level.

    E(rho) = (1-p) rho + (p/d) I_site ⊗ Tr_site(rho)

    Uses tensor reshaping for correct Kronecker product ordering.
    """
    if p <= 0.0:
        return rho

    dim = d**N
    # Reshape to tensor: axes 0..N-1 are ket indices, N..2N-1 are bra indices
    # Axis k corresponds to site k in the Kronecker ordering
    rho_tensor = rho.reshape([d] * N + [d] * N)

    # Partial trace over site: contract ket[site] with bra[site]
    tr_site = np.trace(rho_tensor, axis1=site, axis2=site + N)

    # Build (I_site / d) ⊗ Tr_site(rho)
    mixed = np.zeros_like(rho_tensor)
    for a in range(d):
        idx = [slice(None)] * (2 * N)
        idx[site] = a
        idx[site + N] = a
        mixed[tuple(idx)] = (1.0 / d) * tr_site

    return ((1 - p) * rho_tensor + p * mixed).reshape(dim, dim)


def apply_dephasing_channel_dm(
    rho: np.ndarray, site: int, d: int, N: int, p: float
) -> np.ndarray:
    """Apply single-site dephasing channel at density-matrix level.

    E(rho) = (1-p) rho + p * sum_k P_k rho P_k
    where P_k = |k><k|_site ⊗ I_rest

    The dephased part zeroes out all off-diagonal elements in the site's subspace.
    """
    if p <= 0.0:
        return rho

    dim = d**N
    rho_tensor = rho.reshape([d] * N + [d] * N)

    # Dephased part: keep only elements where ket[site] == bra[site]
    rho_dephased = np.zeros_like(rho_tensor)
    for k in range(d):
        idx = [slice(None)] * (2 * N)
        idx[site] = k
        idx[site + N] = k
        rho_dephased[tuple(idx)] = rho_tensor[tuple(idx)]

    return ((1 - p) * rho_tensor + p * rho_dephased).reshape(dim, dim)


def apply_pair_depolarization_channel_dm(
    rho: np.ndarray, site_a: int, site_b: int, d: int, N: int, p: float
) -> np.ndarray:
    """Apply pair depolarization channel at density-matrix level.

    E(rho) = (1-p) rho + (p/d^2) I_{AB} ⊗ Tr_{AB}(rho)

    Uses tensor reshaping for correct Kronecker product ordering.
    """
    if p <= 0.0:
        return rho

    d_pair = d * d
    dim = d**N
    rho_tensor = rho.reshape([d] * N + [d] * N)

    # Partial trace over both sites a and b
    # Trace over the higher-indexed site first to keep index arithmetic simple
    s_hi = max(site_a, site_b)
    s_lo = min(site_a, site_b)

    tr_hi = np.trace(rho_tensor, axis1=s_hi, axis2=s_hi + N)
    # After removing s_hi, axis for s_lo is unchanged (it's below s_hi)
    # but the bra axis for s_lo has shifted down by 1
    tr_both = np.trace(tr_hi, axis1=s_lo, axis2=s_lo + (N - 1))

    # Build (I_AB / d^2) ⊗ Tr_AB(rho)
    mixed = np.zeros_like(rho_tensor)
    for a in range(d):
        for b in range(d):
            idx = [slice(None)] * (2 * N)
            idx[site_a] = a
            idx[site_a + N] = a
            idx[site_b] = b
            idx[site_b + N] = b
            mixed[tuple(idx)] = (1.0 / d_pair) * tr_both

    return ((1 - p) * rho_tensor + p * mixed).reshape(dim, dim)


def test_noise_channel_verification(params: GKSLPhysicalParameters) -> dict:
    """Test 1: Verify stochastic noise matches density-matrix channel.

    Runs 1 Trotter step with noise using many shots and compares with
    analytical density-matrix-level noise application.
    """
    print("\n" + "=" * 70)
    print("Test 1: Noise Channel Mathematical Verification")
    print("=" * 70)

    d = params.d
    N = params.N_molecules
    dim = d**N
    p_depol = 0.01
    p_dephasing = 0.005
    dt = 1.0  # Single step with dt=1.0
    n_shots_verify = 5000

    # Build operators
    H_0 = build_onsite_hamiltonian(params)
    H_transfer = build_transfer_hamiltonian(params)
    H_total = H_0 + H_transfer
    lindblad_ops = build_lindblad_operators(params)

    # Prepare initial state (edge_triplet)
    psi_init = np.zeros(dim, dtype=np.complex128)
    index = 1 * (d ** (N - 1)) + 1
    psi_init[index] = 1.0
    rho_init = np.outer(psi_init, psi_init.conj())

    # === Analytical: Apply 1 Trotter step then noise channels ===
    print("\n  Computing analytical (density-matrix) result...", flush=True)
    start = time_module.time()

    # Half Hamiltonian
    U_H_half = expm(-1j * H_total * dt / 2)
    U_stines = [stinespring_unitary_from_lindblad(L, dt) for L, _ in lindblad_ops]

    rho_dm = rho_init.copy()
    # Trotter step (noiseless)
    rho_dm = U_H_half @ rho_dm @ U_H_half.conj().T
    for U_s in U_stines:
        rho_dm = apply_stinespring_to_density_matrix(rho_dm, U_s)
    rho_dm = U_H_half @ rho_dm @ U_H_half.conj().T

    # Now apply noise channels analytically (same order as noisy simulator)
    # Phase 1: Hamiltonian noise (pair depol + single dephase for each neighbor)
    for i, j in params.neighbors:
        rho_dm = apply_pair_depolarization_channel_dm(rho_dm, i, j, d, N, p_depol)
        if p_dephasing > 0:
            rho_dm = apply_dephasing_channel_dm(rho_dm, i, d, N, p_dephasing)
            rho_dm = apply_dephasing_channel_dm(rho_dm, j, d, N, p_dephasing)

    # Phase 2: Lindblad channel noise
    lindblad_sites: list[list[int]] = []
    for i, j in params.neighbors:
        lindblad_sites.append([i, j])
        lindblad_sites.append([i, j])
    for _ in range(5):
        for mol in range(N):
            lindblad_sites.append([mol])

    for k in range(len(lindblad_ops)):
        sites = lindblad_sites[k]
        if len(sites) == 1:
            rho_dm = apply_depolarization_channel_dm(rho_dm, sites[0], d, N, p_depol)
            if p_dephasing > 0:
                rho_dm = apply_dephasing_channel_dm(rho_dm, sites[0], d, N, p_dephasing)
        else:
            rho_dm = apply_pair_depolarization_channel_dm(
                rho_dm, sites[0], sites[1], d, N, p_depol
            )
            if p_dephasing > 0:
                rho_dm = apply_dephasing_channel_dm(rho_dm, sites[0], d, N, p_dephasing)
                rho_dm = apply_dephasing_channel_dm(rho_dm, sites[1], d, N, p_dephasing)

    # Phase 3: Second Hamiltonian noise
    for i, j in params.neighbors:
        rho_dm = apply_pair_depolarization_channel_dm(rho_dm, i, j, d, N, p_depol)
        if p_dephasing > 0:
            rho_dm = apply_dephasing_channel_dm(rho_dm, i, d, N, p_dephasing)
            rho_dm = apply_dephasing_channel_dm(rho_dm, j, d, N, p_dephasing)

    t_dm = time_module.time() - start
    print(f"  Analytical computation: {t_dm:.2f}s")
    print(f"  Tr(rho_dm) = {np.real(np.trace(rho_dm)):.10f}")

    # === Stochastic: Run many shots ===
    print(f"  Computing stochastic ({n_shots_verify} shots) result...", flush=True)
    start = time_module.time()

    sim_noisy = QuditGKSLNoisyShotSimulator(params, p_depol=p_depol, p_dephasing=p_dephasing)
    result_noisy = sim_noisy.simulate(
        t_max=dt, n_steps=1, initial_state="edge_triplet",
        n_shots=n_shots_verify, seed=42,
    )
    rho_stoch = result_noisy["rho_final"]
    t_stoch = time_module.time() - start
    print(f"  Stochastic computation: {t_stoch:.2f}s")
    print(f"  Tr(rho_stoch) = {np.real(np.trace(rho_stoch)):.10f}")

    # === Compare ===
    fid = quantum_fidelity(rho_dm, rho_stoch)
    trace_diff = abs(np.real(np.trace(rho_dm)) - np.real(np.trace(rho_stoch)))
    frobenius_diff = np.linalg.norm(rho_dm - rho_stoch, "fro")

    print(f"\n  === Comparison ===")
    print(f"  Fidelity(analytical, stochastic): {fid:.6f}")
    print(f"  |Tr(analytical) - Tr(stochastic)|: {trace_diff:.2e}")
    print(f"  ||analytical - stochastic||_F: {frobenius_diff:.6f}")
    print(f"  Expected: F close to 1.0 (deviations from finite shot noise)")

    # Also run noiseless for reference
    sim_noiseless = QuditGKSLShotSimulator(params)
    result_noiseless = sim_noiseless.simulate(
        t_max=dt, n_steps=1, initial_state="edge_triplet",
        n_shots=n_shots_verify, seed=42,
    )
    rho_noiseless = result_noiseless["rho_final"]
    fid_noisy_vs_noiseless = quantum_fidelity(rho_noiseless, rho_stoch)
    print(f"  Fidelity(noiseless_shot, noisy_shot): {fid_noisy_vs_noiseless:.6f}")

    return {
        "test": "noise_channel_verification",
        "n_shots": n_shots_verify,
        "fidelity_analytical_vs_stochastic": fid,
        "trace_diff": trace_diff,
        "frobenius_diff": frobenius_diff,
        "fidelity_noiseless_vs_noisy": fid_noisy_vs_noiseless,
        "trace_analytical": float(np.real(np.trace(rho_dm))),
        "trace_stochastic": float(np.real(np.trace(rho_stoch))),
        "passed": fid > 0.95,
    }


def test_per_step_fidelity_decay(params: GKSLPhysicalParameters) -> dict:
    """Test 2: Track per-step fidelity decay for noisy simulations."""
    print("\n" + "=" * 70)
    print("Test 2: Per-Step Fidelity Decay Tracking")
    print("=" * 70)

    t_max = 10.0
    n_steps = 10
    n_shots = 1000
    seed = 42
    p_depol = 0.01
    p_dephasing = 0.005

    # Run noiseless density matrix simulations for reference
    print("\n  Running reference simulations...", flush=True)
    sim_qudit_dm = QuditGKSLSimulator(params)
    result_qudit_dm = sim_qudit_dm.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
    )

    sim_qubit_dm = QubitGKSLSimulator(params)
    result_qubit_dm = sim_qubit_dm.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
    )

    # Run noisy shot simulations
    print("  Running qudit noisy shot...", flush=True)
    sim_qudit_noisy = QuditGKSLNoisyShotSimulator(params, p_depol=p_depol, p_dephasing=p_dephasing)
    result_qudit_noisy = sim_qudit_noisy.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet",
        n_shots=n_shots, seed=seed,
    )

    print("  Running qubit noisy shot...", flush=True)
    sim_qubit_noisy = QubitGKSLNoisyShotSimulator(params, p_depol=p_depol, p_dephasing=p_dephasing)
    result_qubit_noisy = sim_qubit_noisy.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet",
        n_shots=n_shots, seed=seed,
    )

    # We need per-step density matrices, which are stored in rho_accum
    # Unfortunately the current simulators only return rho_final.
    # So we compute fidelity from the population and trace time series instead.

    # Track trace (leakage) per step
    print("\n  === Per-Step Trace ===")
    print(f"  {'Step':>4} | {'Qudit noisy Tr':>16} | {'Qubit noisy Tr':>16} | "
          f"{'Qubit deficit':>14}")
    print("  " + "-" * 70)

    step_data = []
    for s in range(n_steps + 1):
        tr_qd = result_qudit_noisy["trace"][s]
        tr_qb = result_qubit_noisy["trace"][s]
        deficit = 1.0 - tr_qb
        step_data.append({
            "step": s,
            "time": result_qudit_noisy["times"][s],
            "trace_qudit_noisy": float(tr_qd),
            "trace_qubit_noisy": float(tr_qb),
            "qubit_deficit": float(deficit),
        })
        print(f"  {s:4d} | {tr_qd:16.10f} | {tr_qb:16.10f} | {deficit:14.6f}")

    # Track population deviation per step
    print("\n  === Per-Step Population Deviation (vs DM reference) ===")
    print(f"  {'Step':>4} | {'Qudit |ΔN_S0|':>14} | {'Qudit |ΔN_T1|':>14} | "
          f"{'Qubit |ΔN_S0|':>14} | {'Qubit |ΔN_T1|':>14}")
    print("  " + "-" * 70)

    for s in range(n_steps + 1):
        pop_qd_noisy = result_qudit_noisy["populations"][s]
        pop_qd_dm = result_qudit_dm["populations"][s]
        pop_qb_noisy = result_qubit_noisy["populations"][s]
        pop_qb_dm = result_qubit_dm["populations"][s]

        d_s0_qd = abs(pop_qd_noisy["N_S0"] - pop_qd_dm["N_S0"])
        d_t1_qd = abs(pop_qd_noisy["N_T1"] - pop_qd_dm["N_T1"])
        d_s0_qb = abs(pop_qb_noisy["N_S0"] - pop_qb_dm["N_S0"])
        d_t1_qb = abs(pop_qb_noisy["N_T1"] - pop_qb_dm["N_T1"])

        step_data[s]["delta_N_S0_qudit"] = float(d_s0_qd)
        step_data[s]["delta_N_T1_qudit"] = float(d_t1_qd)
        step_data[s]["delta_N_S0_qubit"] = float(d_s0_qb)
        step_data[s]["delta_N_T1_qubit"] = float(d_t1_qb)

        print(f"  {s:4d} | {d_s0_qd:14.6f} | {d_t1_qd:14.6f} | "
              f"{d_s0_qb:14.6f} | {d_t1_qb:14.6f}")

    # Final state fidelity comparison
    print("\n  === Final State Fidelity ===")
    rho_qd_dm = result_qudit_dm["rho_final"]
    rho_qd_noisy = result_qudit_noisy["rho_final"]
    rho_qb_dm = result_qubit_dm["rho_final"]
    rho_qb_noisy = result_qubit_noisy["rho_final"]

    f_qd = quantum_fidelity(rho_qd_dm, rho_qd_noisy)
    f_qb_raw = quantum_fidelity(rho_qb_dm, rho_qb_noisy)
    f_qb_norm = quantum_fidelity(normalize_rho(rho_qb_dm), normalize_rho(rho_qb_noisy))
    f_cross_raw = quantum_fidelity(rho_qd_noisy, rho_qb_noisy)
    f_cross_norm = quantum_fidelity(normalize_rho(rho_qd_noisy), normalize_rho(rho_qb_noisy))

    tr_qd_final = float(np.real(np.trace(rho_qd_noisy)))
    tr_qb_final = float(np.real(np.trace(rho_qb_noisy)))

    print(f"  Qudit DM vs Qudit noisy:     F = {f_qd:.6f} (Tr={tr_qd_final:.6f})")
    print(f"  Qubit DM vs Qubit noisy raw: F = {f_qb_raw:.6f} (Tr={tr_qb_final:.6f})")
    print(f"  Qubit DM vs Qubit noisy norm:F = {f_qb_norm:.6f}")
    print(f"  Qudit noisy vs Qubit noisy:  F_raw = {f_cross_raw:.6f}, F_norm = {f_cross_norm:.6f}")

    fidelity_results = {
        "qudit_dm_vs_qudit_noisy": f_qd,
        "qubit_dm_vs_qubit_noisy_raw": f_qb_raw,
        "qubit_dm_vs_qubit_noisy_norm": f_qb_norm,
        "qudit_noisy_vs_qubit_noisy_raw": f_cross_raw,
        "qudit_noisy_vs_qubit_noisy_norm": f_cross_norm,
        "trace_qudit_noisy_final": tr_qd_final,
        "trace_qubit_noisy_final": tr_qb_final,
    }

    return {
        "test": "per_step_fidelity_decay",
        "config": {
            "t_max": t_max,
            "n_steps": n_steps,
            "n_shots": n_shots,
            "seed": seed,
            "p_depol": p_depol,
            "p_dephasing": p_dephasing,
        },
        "step_data": step_data,
        "fidelity_results": fidelity_results,
    }


def test_error_budget(params: GKSLPhysicalParameters) -> dict:
    """Test 3: Error budget analysis."""
    print("\n" + "=" * 70)
    print("Test 3: Error Budget Analysis")
    print("=" * 70)

    d = params.d
    N = params.N_molecules
    n_neighbors = len(params.neighbors)
    n_lindblad = 2 * n_neighbors + 5 * N

    # Count noise gates per Trotter step (qudit)
    # Phase 1 (Hamiltonian half): pair depol per neighbor + dephasing per site
    ham_pair_depol = n_neighbors  # 3
    ham_dephasing = 2 * n_neighbors  # 6

    # Phase 2 (Lindblad channels)
    pair_lindblad = 2 * n_neighbors  # 6 pair channels
    single_lindblad = 5 * N  # 20 single channels

    lindblad_pair_depol = pair_lindblad  # 6
    lindblad_pair_dephase = 2 * pair_lindblad  # 12
    lindblad_single_depol = single_lindblad  # 20
    lindblad_single_dephase = single_lindblad  # 20

    # Phase 3 (Hamiltonian half): same as Phase 1
    total_pair_depol = 2 * ham_pair_depol + lindblad_pair_depol  # 12
    total_single_depol = lindblad_single_depol  # 20
    total_dephasing = 2 * ham_dephasing + lindblad_pair_dephase + lindblad_single_dephase  # 44
    total_noise_gates = total_pair_depol + total_single_depol + total_dephasing  # 76

    print(f"\n  === Noise Gate Count Per Trotter Step (Qudit) ===")
    print(f"  Pair depolarization:   {total_pair_depol} applications")
    print(f"  Single depolarization: {total_single_depol} applications")
    print(f"  Dephasing:             {total_dephasing} applications")
    print(f"  Total noise gates:     {total_noise_gates} per step")

    p_depol = 0.01
    p_dephasing = 0.005

    # Expected errors per step
    # (1 - 1/d^n) = probability that a random Weyl-Heisenberg operator is non-identity
    # For pair (d^2 local dim): (d^4 - 1)/d^4 non-identity operators out of d^4 total
    # For single (d local dim): (d^2 - 1)/d^2 non-identity operators out of d^2 total
    exp_pair_errors = total_pair_depol * p_depol * (1 - 1 / (d**4))
    exp_single_errors = total_single_depol * p_depol * (1 - 1 / (d**2))
    exp_dephase_events = total_dephasing * p_dephasing

    print(f"\n  === Expected Errors Per Step (p_depol={p_depol}, p_dephasing={p_dephasing}) ===")
    print(f"  Pair depol errors:      {exp_pair_errors:.4f}")
    print(f"  Single depol errors:    {exp_single_errors:.4f}")
    print(f"  Dephasing events:       {exp_dephase_events:.4f}")
    print(f"  Total error events:     {exp_pair_errors + exp_single_errors + exp_dephase_events:.4f}")

    # For 100 steps (notebook settings)
    n_steps_notebook = 100
    total_errors_100 = (exp_pair_errors + exp_single_errors + exp_dephase_events) * n_steps_notebook
    print(f"\n  === Cumulative Errors for {n_steps_notebook} Steps ===")
    print(f"  Total expected error events: {total_errors_100:.1f}")
    print(f"  This is extremely high and explains the near-complete depolarization")
    print(f"  observed in the notebook (F ≈ 0.136 for qudit noisy)")

    # Qubit specific: leakage analysis
    # In 2-qubit encoding (|00>=S0, |01>=T1, |10>=S1, |11>=forbidden),
    # there are 4^2=16 two-qubit Pauli operators, of which 15 are non-identity.
    # Only 3 (I⊗Z, Z⊗I, Z⊗Z) preserve the physical subspace {|00>,|01>,|10>}.
    # The remaining 12 out of 15 non-identity Paulis cause leakage to |11>.
    n_paulis_total = 16  # 4^2 two-qubit Pauli operators
    n_paulis_nonidentity = n_paulis_total - 1  # 15 non-identity operators
    n_leakage_paulis = 12  # operators that map physical states to |11> (forbidden)
    p_leakage_per_single_depol = p_depol * (1 - 1 / n_paulis_total) * (n_leakage_paulis / n_paulis_nonidentity)
    p_leakage_per_step = 1 - (1 - p_leakage_per_single_depol) ** total_single_depol
    # This is approximate; pair depol also contributes

    print(f"\n  === Qubit Leakage Analysis ===")
    print(f"  Leakage-causing Paulis (per molecule): {n_leakage_paulis}/{n_paulis_nonidentity} (80%)")
    print(f"  p(leakage per single depol): {p_leakage_per_single_depol:.6f}")
    print(f"  p(any leakage per step, single only): ~{p_leakage_per_step:.4f}")

    return {
        "test": "error_budget",
        "noise_gates_per_step": {
            "pair_depol": total_pair_depol,
            "single_depol": total_single_depol,
            "dephasing": total_dephasing,
            "total": total_noise_gates,
        },
        "expected_errors_per_step": {
            "pair_depol": exp_pair_errors,
            "single_depol": exp_single_errors,
            "dephasing": exp_dephase_events,
            "total": exp_pair_errors + exp_single_errors + exp_dephase_events,
        },
        "cumulative_100_steps": total_errors_100,
        "qubit_leakage": {
            "leakage_paulis_fraction": n_leakage_paulis / 15,
            "p_leakage_per_single_depol": p_leakage_per_single_depol,
        },
    }


def main() -> int:
    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    output_dir = Path(__file__).resolve().parents[1] / "developing" / "verification_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    params = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
    )

    print("=" * 70)
    print("TTA-UC GKSL Iteration 12: Noise Channel Mathematical Verification")
    print(f"Timestamp: {timestamp}")
    print("=" * 70)

    results: dict = {
        "timestamp": timestamp,
        "iteration": 12,
        "purpose": "Noise channel mathematical verification + per-step fidelity tracking",
        "tests": [],
    }

    # Test 1: Noise channel verification
    t1_result = test_noise_channel_verification(params)
    results["tests"].append(t1_result)

    # Test 2: Per-step fidelity decay
    t2_result = test_per_step_fidelity_decay(params)
    results["tests"].append(t2_result)

    # Test 3: Error budget
    t3_result = test_error_budget(params)
    results["tests"].append(t3_result)

    # ===================================================================
    # Write JSON report
    # ===================================================================
    json_path = output_dir / f"iteration12_noise_verification_{timestamp}.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # ===================================================================
    # Write Markdown report
    # ===================================================================
    md_lines = [
        "# TTA-UC GKSL Iteration 12: ノイズチャネル数学的検証結果",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        "",
        "## Test 1: ノイズチャネルの確率的実装 vs 密度行列解析的チャネル",
        "",
        f"- ショット数: {t1_result['n_shots']}",
        f"- 忠実度(解析的 vs 確率的): **{t1_result['fidelity_analytical_vs_stochastic']:.6f}**",
        f"- トレース差: {t1_result['trace_diff']:.2e}",
        f"- Frobenius距離: {t1_result['frobenius_diff']:.6f}",
        f"- 忠実度(ノイズ無し vs ノイズ有り): {t1_result['fidelity_noiseless_vs_noisy']:.6f}",
        f"- 判定: **{'PASS' if t1_result['passed'] else 'FAIL'}**",
        "",
        "### 解釈",
        "",
        "忠実度が1.0に近ければ、確率的ノイズ実装が密度行列チャネルと一致している。",
        "有限ショット数による統計的揺らぎ分だけ1.0から乖離する。",
        "",
    ]

    # Test 2 results
    t2_fid = t2_result["fidelity_results"]
    t2_config = t2_result["config"]
    md_lines.extend([
        "## Test 2: ステップ毎忠実度劣化追跡",
        "",
        f"- t_max: {t2_config['t_max']}, n_steps: {t2_config['n_steps']}, n_shots: {t2_config['n_shots']}",
        f"- ノイズ: p_depol={t2_config['p_depol']}, p_dephasing={t2_config['p_dephasing']}",
        "",
        "### トレース時系列",
        "",
        "| Step | Time | Qudit noisy Tr | Qubit noisy Tr | Qubit deficit |",
        "|------|------|---------------|---------------|--------------|",
    ])
    for sd in t2_result["step_data"]:
        md_lines.append(
            f"| {sd['step']} | {sd['time']:.1f} | {sd['trace_qudit_noisy']:.10f} | "
            f"{sd['trace_qubit_noisy']:.10f} | {sd['qubit_deficit']:.6f} |"
        )
    md_lines.extend([
        "",
        "### 最終状態忠実度",
        "",
        f"- Qudit DM vs Qudit noisy: F = {t2_fid['qudit_dm_vs_qudit_noisy']:.6f}",
        f"- Qubit DM vs Qubit noisy (raw): F = {t2_fid['qubit_dm_vs_qubit_noisy_raw']:.6f}",
        f"- Qubit DM vs Qubit noisy (norm): F = {t2_fid['qubit_dm_vs_qubit_noisy_norm']:.6f}",
        f"- Qudit noisy vs Qubit noisy (raw): F = {t2_fid['qudit_noisy_vs_qubit_noisy_raw']:.6f}",
        f"- Qudit noisy vs Qubit noisy (norm): F = {t2_fid['qudit_noisy_vs_qubit_noisy_norm']:.6f}",
        "",
    ])

    # Test 3 results
    t3_gates = t3_result["noise_gates_per_step"]
    t3_errors = t3_result["expected_errors_per_step"]
    md_lines.extend([
        "## Test 3: エラーバジェット分析",
        "",
        "### 1ステップあたりのノイズゲート数（Qudit）",
        "",
        f"- ペア脱分極: {t3_gates['pair_depol']}回",
        f"- 単一サイト脱分極: {t3_gates['single_depol']}回",
        f"- 位相緩和: {t3_gates['dephasing']}回",
        f"- **合計: {t3_gates['total']}回/ステップ**",
        "",
        "### 1ステップあたりの期待エラー数",
        "",
        f"- ペア脱分極エラー: {t3_errors['pair_depol']:.4f}",
        f"- 単一サイト脱分極エラー: {t3_errors['single_depol']:.4f}",
        f"- 位相緩和イベント: {t3_errors['dephasing']:.4f}",
        f"- **合計: {t3_errors['total']:.4f}**",
        "",
        f"### 100ステップ累積エラー数: {t3_result['cumulative_100_steps']:.1f}",
        "",
        "これは極めて大量のノイズであり、量子状態がほぼ完全に",
        "最大混合状態に近づくことを意味する。",
        "",
    ])

    md_path = output_dir / f"iteration12_noise_verification_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Test 1 (Noise channel verification): {'PASS' if t1_result['passed'] else 'FAIL'}")
    print(f"    F(analytical, stochastic) = {t1_result['fidelity_analytical_vs_stochastic']:.6f}")
    print(f"  Test 2 (Per-step fidelity decay):")
    print(f"    Qudit F(DM, noisy) = {t2_fid['qudit_dm_vs_qudit_noisy']:.6f}")
    print(f"    Qubit F(DM, noisy) norm = {t2_fid['qubit_dm_vs_qubit_noisy_norm']:.6f}")
    print(f"    Cross F(qudit, qubit) norm = {t2_fid['qudit_noisy_vs_qubit_noisy_norm']:.6f}")
    print(f"  Test 3 (Error budget):")
    print(f"    Noise gates/step: {t3_gates['total']}")
    print(f"    Expected errors/step: {t3_errors['total']:.4f}")
    print(f"    Cumulative errors (100 steps): {t3_result['cumulative_100_steps']:.1f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
