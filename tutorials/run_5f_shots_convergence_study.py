"""§5.1-4: n_shots convergence study for the quantum-trajectory method (5f).

Measures, at the notebook's production parameters (N=4, t_max=100,
n_steps=100, edge_triplet), the Frobenius distance between the
trajectory-averaged density matrix of :class:`QuditGKSLShotSimulator`
(the scenario-5f simulator) and the deterministic density-matrix
Stinespring solution of :class:`QuditGKSLSimulator` — which uses the
identical Trotter scheme, so the difference is purely the Monte-Carlo
statistical error and must scale as ``O(1/sqrt(n_shots))``.

Additionally cross-checks the *backend-executed* trajectory mode
(``QuditGKSLSimulator(execute_on_backend='tnsim')``) at N=4 with a
moderate trajectory count.

Run from the ``tutorials`` directory:

    python run_5f_shots_convergence_study.py

All printed numbers are actual measurements from this run.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_shot_simulator import QuditGKSLShotSimulator
from qudit_gksl_simulator import QuditGKSLSimulator

T_MAX = 100.0
N_STEPS = 100
N_SHOTS_LIST = [50, 200, 800, 3200]
SEED = 12345


def main() -> None:
    params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)

    print("=== reference: deterministic density-matrix Stinespring (NumPy) ===")
    t0 = time.time()
    ref = QuditGKSLSimulator(params, algorithm="stinespring").simulate(
        t_max=T_MAX, n_steps=N_STEPS, initial_state="edge_triplet"
    )
    rho_ref = ref["rho_final"]
    print(f"reference done in {time.time() - t0:.1f} s")

    print()
    print("=== n_shots convergence (QuditGKSLShotSimulator, scenario 5f) ===")
    print(f"{'n_shots':>8} {'‖Δρ‖_F':>12} {'1/√n':>10} {'ratio':>8} {'time [s]':>10}")
    errors: list[float] = []
    for n_shots in N_SHOTS_LIST:
        sim = QuditGKSLShotSimulator(params)
        t0 = time.time()
        res = sim.simulate(
            t_max=T_MAX,
            n_steps=N_STEPS,
            initial_state="edge_triplet",
            n_shots=n_shots,
            seed=SEED,
        )
        elapsed = time.time() - t0
        err = float(np.linalg.norm(res["rho_final"] - rho_ref))
        errors.append(err)
        inv_sqrt = 1.0 / np.sqrt(n_shots)
        print(
            f"{n_shots:>8d} {err:>12.4e} {inv_sqrt:>10.4e} "
            f"{err / inv_sqrt:>8.3f} {elapsed:>10.1f}"
        )

    # Log-log slope: err ~ C n^p  =>  p should be close to -1/2.
    log_n = np.log(np.asarray(N_SHOTS_LIST, dtype=float))
    log_e = np.log(np.asarray(errors))
    slope = float(np.polyfit(log_n, log_e, 1)[0])
    print(f"\nlog-log fitted slope p (err ~ n^p): {slope:.3f}  (theory: -0.5)")

    print()
    print("=== cross-check: backend-executed trajectories (tnsim), N=4 ===")
    n_traj = 50
    t0 = time.time()
    res_tn = QuditGKSLSimulator(
        params,
        algorithm="stinespring",
        execute_on_backend="tnsim",
        n_trajectories=n_traj,
        seed=SEED,
    ).simulate(t_max=T_MAX, n_steps=N_STEPS, initial_state="edge_triplet")
    elapsed = time.time() - t0
    err_tn = float(np.linalg.norm(res_tn["rho_final"] - rho_ref))
    bound = 5.0 / np.sqrt(n_traj)
    print(
        f"n_trajectories={n_traj}: ‖Δρ‖_F = {err_tn:.4e} "
        f"(5/√n bound {bound:.4e}), time {elapsed:.1f} s"
    )
    if err_tn >= bound:
        msg = "tnsim trajectory average outside statistical bound"
        raise AssertionError(msg)
    print("tnsim cross-check PASSED")


if __name__ == "__main__":
    main()
