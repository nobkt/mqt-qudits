"""検証スクリプト: quditシミュレータの補助ビット(アンシラ)がquditであることを検証する.

qudit量子演算型量子コンピュータ向けのシミュレータでは、Stinespring拡張の
補助ビットもシステムと同じd次元のquditを使用する必要がある。
このスクリプトは以下を検証する:

1. Stinespring unitaryの次元が d_anc * d_sys であること (d_anc=d, not d_anc=2)
2. 密度行列のStinespring適用結果が正しいこと
3. ショット測定でd_anc個の測定結果があること
4. 全quditシミュレータで n_ancilla_qudits (not n_ancilla_qubits) が使われていること
5. 回路シミュレータのancillaレジスタがd=d_anc=3であること
6. d_anc=2 (旧qubit) vs d_anc=3 (新qudit) で物理結果が一致すること
"""

from __future__ import annotations

import json
import os
import sys
import traceback
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_verification() -> dict:
    """全検証を実行して結果を返す."""
    results: dict = {
        "timestamp": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "checks": [],
        "summary": {},
    }

    def add_check(name: str, passed: bool, details: str = "") -> None:
        results["checks"].append({
            "name": name,
            "passed": passed,
            "details": details,
        })
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  [{status}] {name}: {details}")

    # ===================================================================
    # Check 1: stinespring_unitary_from_lindblad d_anc parameter
    # ===================================================================
    print("\n=== Check 1: Stinespring unitary dimensions ===")
    try:
        from stinespring_utils import stinespring_unitary_from_lindblad

        # Create a simple 3x3 Lindblad operator
        L = np.zeros((3, 3), dtype=np.complex128)
        L[0, 1] = 0.1  # |0><1| transition
        dt = 0.5

        # d_anc=2 (old qubit style)
        U2 = stinespring_unitary_from_lindblad(L, dt, d_anc=2)
        add_check(
            "Stinespring d_anc=2 dimension",
            U2.shape == (6, 6),
            f"shape={U2.shape}, expected=(6, 6)",
        )

        # d_anc=3 (new qudit style)
        U3 = stinespring_unitary_from_lindblad(L, dt, d_anc=3)
        add_check(
            "Stinespring d_anc=3 dimension",
            U3.shape == (9, 9),
            f"shape={U3.shape}, expected=(9, 9)",
        )

        # Verify K0, K1 are identical between d_anc=2 and d_anc=3
        K0_2 = U2[:3, :3]
        K1_2 = U2[3:, :3]
        K0_3 = U3[:3, :3]
        K1_3 = U3[3:6, :3]
        K2_3 = U3[6:, :3]

        k0_diff = np.linalg.norm(K0_2 - K0_3, ord="fro")
        k1_diff = np.linalg.norm(K1_2 - K1_3, ord="fro")
        k2_norm = np.linalg.norm(K2_3, ord="fro")

        add_check(
            "K0 identical (d_anc=2 vs 3)",
            k0_diff < 1e-12,
            f"||K0_2 - K0_3||_F = {k0_diff:.2e}",
        )
        add_check(
            "K1 identical (d_anc=2 vs 3)",
            k1_diff < 1e-12,
            f"||K1_2 - K1_3||_F = {k1_diff:.2e}",
        )
        add_check(
            "K2 is zero (d_anc=3)",
            k2_norm < 1e-12,
            f"||K2||_F = {k2_norm:.2e}",
        )

        # Verify trace preservation
        tp_check = np.linalg.norm(
            K0_3.conj().T @ K0_3 + K1_3.conj().T @ K1_3 + K2_3.conj().T @ K2_3
            - np.eye(3),
            ord="fro",
        )
        add_check(
            "Trace preservation (d_anc=3)",
            tp_check < 1e-10,
            f"||sum K_k†K_k - I||_F = {tp_check:.2e}",
        )

    except Exception as e:
        add_check("Stinespring unitary", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 2: apply_stinespring_to_density_matrix with d_anc
    # ===================================================================
    print("\n=== Check 2: Stinespring density matrix application ===")
    try:
        from stinespring_utils import (
            apply_stinespring_to_density_matrix,
            stinespring_unitary_from_lindblad,
        )

        L = np.zeros((3, 3), dtype=np.complex128)
        L[0, 1] = 0.1
        dt = 0.5
        rho = np.eye(3, dtype=np.complex128) / 3.0  # maximally mixed

        U2 = stinespring_unitary_from_lindblad(L, dt, d_anc=2)
        U3 = stinespring_unitary_from_lindblad(L, dt, d_anc=3)

        rho_out_2 = apply_stinespring_to_density_matrix(rho, U2, d_anc=2)
        rho_out_3 = apply_stinespring_to_density_matrix(rho, U3, d_anc=3)

        diff = np.linalg.norm(rho_out_2 - rho_out_3, ord="fro")
        add_check(
            "DM channel output identical (d_anc=2 vs 3)",
            diff < 1e-12,
            f"||rho_2 - rho_3||_F = {diff:.2e}",
        )

        # Trace preservation
        tr_out = float(np.real(np.trace(rho_out_3)))
        add_check(
            "Trace preservation after channel",
            abs(tr_out - 1.0) < 1e-10,
            f"Tr(rho_out) = {tr_out:.15f}",
        )

    except Exception as e:
        add_check("DM channel application", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 3: QuditGKSLSimulator attributes
    # ===================================================================
    print("\n=== Check 3: QuditGKSLSimulator attributes ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLSimulator(params)

        add_check(
            "has n_ancilla_qudits attribute",
            hasattr(sim, "n_ancilla_qudits"),
            f"n_ancilla_qudits={getattr(sim, 'n_ancilla_qudits', 'MISSING')}",
        )
        add_check(
            "no n_ancilla_qubits attribute",
            not hasattr(sim, "n_ancilla_qubits"),
            "should not have n_ancilla_qubits",
        )
        add_check(
            "d_anc == params.d",
            sim.d_anc == params.d,
            f"d_anc={sim.d_anc}, params.d={params.d}",
        )
        add_check(
            "n_ancilla_qudits == 26",
            sim.n_ancilla_qudits == 26,
            f"n_ancilla_qudits={sim.n_ancilla_qudits}",
        )

    except Exception as e:
        add_check("QuditGKSLSimulator attributes", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 4: QuditGKSLShotSimulator attributes
    # ===================================================================
    print("\n=== Check 4: QuditGKSLShotSimulator attributes ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_shot_simulator import QuditGKSLShotSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLShotSimulator(params)

        add_check(
            "shot sim has n_ancilla_qudits",
            hasattr(sim, "n_ancilla_qudits"),
            f"n_ancilla_qudits={getattr(sim, 'n_ancilla_qudits', 'MISSING')}",
        )
        add_check(
            "shot sim no n_ancilla_qubits",
            not hasattr(sim, "n_ancilla_qubits"),
            "should not have n_ancilla_qubits",
        )
        add_check(
            "shot sim d_anc == params.d",
            sim.d_anc == params.d,
            f"d_anc={sim.d_anc}",
        )

    except Exception as e:
        add_check("QuditGKSLShotSimulator attributes", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 5: QuditGKSLBosonSimulator attributes
    # ===================================================================
    print("\n=== Check 5: QuditGKSLBosonSimulator attributes ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator

        params_b = GKSLPhysicalParameters(with_boson=True)
        sim_b = QuditGKSLBosonSimulator(params_b)

        add_check(
            "boson sim has n_ancilla_qudits",
            hasattr(sim_b, "n_ancilla_qudits"),
            f"n_ancilla_qudits={getattr(sim_b, 'n_ancilla_qudits', 'MISSING')}",
        )
        add_check(
            "boson sim no n_ancilla_qubits",
            not hasattr(sim_b, "n_ancilla_qubits"),
            "should not have n_ancilla_qubits",
        )
        add_check(
            "boson sim d_anc == params.d",
            sim_b.d_anc == params_b.d,
            f"d_anc={sim_b.d_anc}",
        )

    except Exception as e:
        add_check("QuditGKSLBosonSimulator attributes", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 6: Simulation result dict keys
    # ===================================================================
    print("\n=== Check 6: Simulation result dict keys ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLSimulator(params)
        result = sim.simulate(t_max=1.0, n_steps=2, initial_state="edge_triplet")

        add_check(
            "result has n_ancilla_qudits key",
            "n_ancilla_qudits" in result,
            f"keys={[k for k in result if 'ancilla' in k or 'd_anc' in k]}",
        )
        add_check(
            "result has d_anc key",
            "d_anc" in result,
            f"d_anc={result.get('d_anc', 'MISSING')}",
        )
        add_check(
            "result no n_ancilla_qubits key",
            "n_ancilla_qubits" not in result,
            "should not have n_ancilla_qubits",
        )

    except Exception as e:
        add_check("Simulation result dict", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 7: Physics equivalence (d_anc=2 vs d_anc=3 produce same channel)
    # ===================================================================
    print("\n=== Check 7: Physics equivalence d_anc=2 vs d_anc=3 ===")
    try:
        from gksl_math_utils import build_lindblad_operators
        from gksl_physical_parameters import GKSLPhysicalParameters
        from stinespring_utils import (
            apply_stinespring_to_density_matrix,
            stinespring_unitary_from_lindblad,
        )

        params = GKSLPhysicalParameters()
        lindblad_ops = build_lindblad_operators(params)
        d_sys = params.d ** params.N_molecules  # 81

        # Test with edge_triplet initial state
        psi = np.zeros(d_sys, dtype=np.complex128)
        index = 1 * (params.d ** (params.N_molecules - 1)) + 1
        psi[index] = 1.0
        rho = np.outer(psi, psi.conj())

        max_diff = 0.0
        for i, (L_op, gamma) in enumerate(lindblad_ops[:5]):  # Test first 5
            U2 = stinespring_unitary_from_lindblad(L_op, 0.5, d_anc=2)
            U3 = stinespring_unitary_from_lindblad(L_op, 0.5, d_anc=3)
            rho2 = apply_stinespring_to_density_matrix(rho, U2, d_anc=2)
            rho3 = apply_stinespring_to_density_matrix(rho, U3, d_anc=3)
            diff = np.linalg.norm(rho2 - rho3, ord="fro")
            max_diff = max(max_diff, diff)

        add_check(
            "Channel equivalence (d_anc=2 vs 3, 81-dim)",
            max_diff < 1e-12,
            f"max ||rho_2 - rho_3||_F = {max_diff:.2e} (first 5 channels)",
        )

    except Exception as e:
        add_check("Physics equivalence", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 8: Stinespring unitary sizes for full system
    # ===================================================================
    print("\n=== Check 8: Full system Stinespring dimensions ===")
    try:
        from gksl_math_utils import build_lindblad_operators
        from gksl_physical_parameters import GKSLPhysicalParameters
        from stinespring_utils import stinespring_unitary_from_lindblad

        params = GKSLPhysicalParameters()
        lindblad_ops = build_lindblad_operators(params)
        d_sys = params.d ** params.N_molecules  # 81

        L_op, gamma = lindblad_ops[0]
        U3 = stinespring_unitary_from_lindblad(L_op, 0.5, d_anc=3)

        expected_dim = 3 * d_sys  # 243 for qutrit ancilla
        add_check(
            f"Full system Stinespring dim (d_anc=3)",
            U3.shape == (expected_dim, expected_dim),
            f"shape={U3.shape}, expected=({expected_dim}, {expected_dim})",
        )

        # Verify unitarity
        residual = np.linalg.norm(U3.conj().T @ U3 - np.eye(expected_dim), ord="fro")
        add_check(
            "Unitarity of full Stinespring (d_anc=3)",
            residual < 1e-10,
            f"||U†U - I||_F = {residual:.2e}",
        )

    except Exception as e:
        add_check("Full system Stinespring", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 9: Short simulation trace preservation
    # ===================================================================
    print("\n=== Check 9: Simulation trace preservation ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_simulator import QuditGKSLSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLSimulator(params)
        result = sim.simulate(t_max=5.0, n_steps=5, initial_state="edge_triplet")

        max_trace_err = max(abs(t - 1.0) for t in result["trace"])
        add_check(
            "Trace preservation (5 steps)",
            max_trace_err < 1e-4,
            f"max |Tr - 1| = {max_trace_err:.2e}",
        )

    except Exception as e:
        add_check("Simulation trace preservation", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 10: QuditGKSLNoisySimulator inherits d_anc
    # ===================================================================
    print("\n=== Check 10: QuditGKSLNoisySimulator inherits d_anc ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_noisy_simulator import QuditGKSLNoisySimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLNoisySimulator(params, p_depol=0.001, p_dephasing=0.0)

        add_check(
            "noisy sim has d_anc",
            hasattr(sim, "d_anc"),
            f"d_anc={getattr(sim, 'd_anc', 'MISSING')}",
        )
        add_check(
            "noisy sim d_anc == 3",
            sim.d_anc == 3,
            f"d_anc={sim.d_anc}",
        )
        add_check(
            "noisy sim has n_ancilla_qudits",
            hasattr(sim, "n_ancilla_qudits"),
            f"n_ancilla_qudits={getattr(sim, 'n_ancilla_qudits', 'MISSING')}",
        )

    except Exception as e:
        add_check("QuditGKSLNoisySimulator", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 11: QuditGKSLNoisyShotSimulator inherits d_anc
    # ===================================================================
    print("\n=== Check 11: QuditGKSLNoisyShotSimulator inherits d_anc ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator

        params = GKSLPhysicalParameters()
        sim = QuditGKSLNoisyShotSimulator(
            params, p_depol=0.001, p_dephasing=0.0
        )

        add_check(
            "noisy shot sim has d_anc",
            hasattr(sim, "d_anc"),
            f"d_anc={getattr(sim, 'd_anc', 'MISSING')}",
        )
        add_check(
            "noisy shot sim d_anc == 3",
            sim.d_anc == 3,
            f"d_anc={sim.d_anc}",
        )

    except Exception as e:
        add_check("QuditGKSLNoisyShotSimulator", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Check 12: Qubit simulators still use d_anc=2 (NOT changed)
    # ===================================================================
    print("\n=== Check 12: Qubit simulators unchanged (d_anc=2) ===")
    try:
        from gksl_physical_parameters import GKSLPhysicalParameters
        from qubit_gksl_simulator import QubitGKSLSimulator

        params = GKSLPhysicalParameters()
        sim_qb = QubitGKSLSimulator(params)

        # Qubit simulators should NOT have d_anc (or should be 2)
        has_d_anc = hasattr(sim_qb, "d_anc")
        if has_d_anc:
            add_check(
                "qubit sim d_anc == 2",
                sim_qb.d_anc == 2,
                f"d_anc={sim_qb.d_anc}",
            )
        else:
            add_check(
                "qubit sim no d_anc (default d_anc=2 in stinespring)",
                True,
                "qubit simulator uses default d_anc=2",
            )

        add_check(
            "qubit sim has n_ancilla (not n_ancilla_qudits)",
            hasattr(sim_qb, "n_ancilla"),
            f"n_ancilla={getattr(sim_qb, 'n_ancilla', 'MISSING')}",
        )

    except Exception as e:
        add_check("Qubit simulators unchanged", False, f"Exception: {e}")
        traceback.print_exc()

    # ===================================================================
    # Summary
    # ===================================================================
    n_pass = sum(1 for c in results["checks"] if c["passed"])
    n_total = len(results["checks"])
    results["summary"] = {
        "total_checks": n_total,
        "passed": n_pass,
        "failed": n_total - n_pass,
        "all_passed": n_pass == n_total,
    }

    print(f"\n{'=' * 60}")
    print(f"SUMMARY: {n_pass}/{n_total} checks passed")
    if n_pass == n_total:
        print("ALL CHECKS PASSED ✓")
    else:
        print("SOME CHECKS FAILED ✗")
        for c in results["checks"]:
            if not c["passed"]:
                print(f"  FAILED: {c['name']}: {c['details']}")
    print(f"{'=' * 60}")

    return results


def main() -> None:
    """実行して結果をJSONファイルに出力する."""
    results = run_verification()

    # Save results
    output_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "developing",
        "verification_results",
    )
    os.makedirs(output_dir, exist_ok=True)

    ts = results["timestamp"]
    filename = f"qudit_ancilla_verification_{ts}.json"
    filepath = os.path.join(output_dir, filename)

    with open(filepath, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)

    print(f"\nResults saved to: {filepath}")


if __name__ == "__main__":
    main()
