"""§5.1-4: n_shots convergence study for the quantum-trajectory method (5f).

Measures the Frobenius distance between the trajectory-averaged density
matrix of :class:`QuditGKSLShotSimulator` (the scenario-5f simulator) and
the deterministic density-matrix Stinespring solution of
:class:`QuditGKSLSimulator` — which uses the identical Trotter scheme, so
the difference is purely the Monte-Carlo statistical error and must scale
as ``O(1/sqrt(n_shots))``.

Additionally cross-checks the *backend-executed* trajectory mode
(``QuditGKSLSimulator(execute_on_backend='tnsim')``) at N=4 with a
moderate trajectory count.

Run from the ``tutorials`` directory::

    python run_5f_shots_convergence_study.py                 # production params
    python run_5f_shots_convergence_study.py --n-steps 20 --processes 4

Notes on determinism: with ``--processes 1`` (default) the shot batches
are executed by a single :class:`QuditGKSLShotSimulator` run with the
given seed, exactly as before.  With ``--processes > 1`` each n_shots
value is split into per-process batches with independent seeds derived
from ``numpy.random.SeedSequence(seed)``, and the batch density matrices
are averaged with weights proportional to the batch sizes.  The result is
statistically equivalent but *not* bit-identical to the single-process
run; the actual configuration used is printed at the start of the run.

All printed numbers are actual measurements from this run.
"""

from __future__ import annotations

import argparse
import os
import sys

_BLAS_ENV_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def _requested_processes(argv: list[str]) -> int:
    """Parse only ``--processes`` from ``argv`` (before heavy imports)."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--processes", type=int, default=1)
    known, _ = pre.parse_known_args(argv)
    return known.processes


def _limit_blas_threads_if_parallel() -> None:
    """Force single-threaded BLAS when running with ``--processes > 1``.

    With ``--processes > 1`` the multi-threaded BLAS of each forked worker
    competes with the other workers for the same CPU cores, which was
    measured to slow down each shot by ~15x.  The thread-count environment
    variables must be set *before* NumPy is imported (fork children inherit
    the parent's already-initialised BLAS thread pool), hence this runs at
    module import time, ahead of ``import numpy``.
    """
    if _requested_processes(sys.argv[1:]) > 1:
        for var in _BLAS_ENV_VARS:
            os.environ.setdefault(var, "1")


_limit_blas_threads_if_parallel()

import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gksl_physical_parameters import GKSLPhysicalParameters
from qudit_gksl_shot_simulator import QuditGKSLShotSimulator
from qudit_gksl_simulator import QuditGKSLSimulator

DEFAULT_T_MAX = 100.0
DEFAULT_N_STEPS = 100
DEFAULT_N_SHOTS_LIST = [50, 200, 800, 3200]
DEFAULT_N_TRAJECTORIES = 50
DEFAULT_SEED = 12345


def _run_shot_batch(args: tuple[float, int, int, int]) -> np.ndarray:
    """Worker: run one batch of trajectories, return the averaged rho_final."""
    t_max, n_steps, n_shots, seed = args
    params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)
    sim = QuditGKSLShotSimulator(params)
    res = sim.simulate(
        t_max=t_max,
        n_steps=n_steps,
        initial_state="edge_triplet",
        n_shots=n_shots,
        seed=seed,
    )
    return np.asarray(res["rho_final"])


def _run_tnsim_batch(args: tuple[float, int, int, int]) -> np.ndarray:
    """Worker: run one batch of tnsim backend trajectories, return rho_final."""
    t_max, n_steps, n_traj, seed = args
    params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)
    res = QuditGKSLSimulator(
        params,
        algorithm="stinespring",
        execute_on_backend="tnsim",
        n_trajectories=n_traj,
        seed=seed,
    ).simulate(t_max=t_max, n_steps=n_steps, initial_state="edge_triplet")
    return np.asarray(res["rho_final"])


def _split_batches(total: int, n_batches: int) -> list[int]:
    """Split ``total`` trajectories into at most ``n_batches`` non-empty batches."""
    base, rem = divmod(total, n_batches)
    sizes = [base + (1 if i < rem else 0) for i in range(n_batches)]
    return [s for s in sizes if s > 0]


def _batch_seeds(seed: int, n_batches: int) -> list[int]:
    return [int(s) for s in np.random.SeedSequence(seed).generate_state(n_batches)]


def _averaged_rho(
    worker,  # noqa: ANN001
    t_max: float,
    n_steps: int,
    total: int,
    seed: int,
    processes: int,
) -> np.ndarray:
    """Run ``total`` trajectories (optionally in parallel) and average rho."""
    if processes <= 1:
        return worker((t_max, n_steps, total, seed))
    batches = _split_batches(total, processes)
    seeds = _batch_seeds(seed, len(batches))
    work = [(t_max, n_steps, b, s) for b, s in zip(batches, seeds)]
    with ProcessPoolExecutor(max_workers=processes) as pool:
        rhos = list(pool.map(worker, work))
    return sum(b * r for b, r in zip(batches, rhos)) / total


def main() -> None:
    parser = argparse.ArgumentParser(description="n_shots convergence study (scenario 5f)")
    parser.add_argument("--t-max", type=float, default=DEFAULT_T_MAX)
    parser.add_argument("--n-steps", type=int, default=DEFAULT_N_STEPS)
    parser.add_argument(
        "--n-shots",
        type=int,
        nargs="+",
        default=DEFAULT_N_SHOTS_LIST,
        help="list of n_shots values for the convergence table",
    )
    parser.add_argument(
        "--n-trajectories",
        type=int,
        default=DEFAULT_N_TRAJECTORIES,
        help="trajectory count for the tnsim backend cross-check (0 = skip)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="number of worker processes (1 = sequential, bit-identical to the original script)",
    )
    args = parser.parse_args()

    print("=== configuration (actual values used in this run) ===", flush=True)
    print(
        f"t_max={args.t_max}, n_steps={args.n_steps}, n_shots={args.n_shots}, "
        f"n_trajectories={args.n_trajectories}, seed={args.seed}, "
        f"processes={args.processes}",
        flush=True,
    )

    params = GKSLPhysicalParameters(N_molecules=4, with_boson=False)

    print("\n=== reference: deterministic density-matrix Stinespring (NumPy) ===", flush=True)
    t0 = time.time()
    ref = QuditGKSLSimulator(params, algorithm="stinespring").simulate(
        t_max=args.t_max, n_steps=args.n_steps, initial_state="edge_triplet"
    )
    rho_ref = ref["rho_final"]
    print(f"reference done in {time.time() - t0:.1f} s", flush=True)

    print()
    print("=== n_shots convergence (QuditGKSLShotSimulator, scenario 5f) ===", flush=True)
    print(f"{'n_shots':>8} {'‖Δρ‖_F':>12} {'1/√n':>10} {'ratio':>8} {'time [s]':>10}", flush=True)
    errors: list[float] = []
    for n_shots in args.n_shots:
        t0 = time.time()
        rho = _averaged_rho(_run_shot_batch, args.t_max, args.n_steps, n_shots, args.seed, args.processes)
        elapsed = time.time() - t0
        err = float(np.linalg.norm(rho - rho_ref))
        errors.append(err)
        inv_sqrt = 1.0 / np.sqrt(n_shots)
        print(
            f"{n_shots:>8d} {err:>12.4e} {inv_sqrt:>10.4e} "
            f"{err / inv_sqrt:>8.3f} {elapsed:>10.1f}",
            flush=True,
        )

    # Log-log slope: err ~ C n^p  =>  p should be close to -1/2.
    if len(args.n_shots) >= 2:
        log_n = np.log(np.asarray(args.n_shots, dtype=float))
        log_e = np.log(np.asarray(errors))
        slope = float(np.polyfit(log_n, log_e, 1)[0])
        print(f"\nlog-log fitted slope p (err ~ n^p): {slope:.3f}  (theory: -0.5)", flush=True)
    else:
        print("\nlog-log slope fit skipped (needs >= 2 n_shots values)", flush=True)

    if args.n_trajectories <= 0:
        print("\ntnsim cross-check skipped (--n-trajectories 0)", flush=True)
        return

    print()
    print("=== cross-check: backend-executed trajectories (tnsim), N=4 ===", flush=True)
    n_traj = args.n_trajectories
    t0 = time.time()
    rho_tn = _averaged_rho(_run_tnsim_batch, args.t_max, args.n_steps, n_traj, args.seed + 1, args.processes)
    elapsed = time.time() - t0
    err_tn = float(np.linalg.norm(rho_tn - rho_ref))
    bound = 5.0 / np.sqrt(n_traj)
    print(
        f"n_trajectories={n_traj}: ‖Δρ‖_F = {err_tn:.4e} "
        f"(5/√n bound {bound:.4e}), time {elapsed:.1f} s",
        flush=True,
    )
    if err_tn >= bound:
        msg = "tnsim trajectory average outside statistical bound"
        raise AssertionError(msg)
    print("tnsim cross-check PASSED", flush=True)


if __name__ == "__main__":
    main()
