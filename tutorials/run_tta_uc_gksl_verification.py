#!/usr/bin/env python3
"""Run TTA-UC GKSL-Lindblad verification and write file-based reports.

This script covers ALL scenarios from quantum_dynamics_gksl_comparison.ipynb:
  - Scenario 1: Classical GKSL (no boson, N=4)
  - Scenario 2: Classical GKSL (with boson, N=2, n_max=1)
  - Scenario 3: Qubit GKSL (no boson, N=4)
  - Scenario 4: Qubit GKSL (with boson, N=2, n_max=1)
  - Scenario 5: Qudit GKSL (no boson, N=4)
  - Scenario 6: Qudit GKSL (with boson, N=2, n_max=1)
  - Scenario 3b: Qubit shot-based (no noise)
  - Scenario 3c: Qubit shot-based (noisy)
  - Scenario 5b: Qudit shot-based (no noise)
  - Scenario 5c: Qudit shot-based (noisy)
  - Unitary limit comparison
  - Cross-scenario fidelity comparison

This script intentionally performs only physics-based checks and does not use
heuristic shortcuts or fallback logic.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_boson_simulator import ClassicalGKSLBosonSimulator
from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_physical_parameters import GKSLPhysicalParameters
from gksl_validation import validate_density_matrix, validate_particle_conservation
from qubit_gksl_boson_simulator import QubitGKSLBosonSimulator
from qubit_gksl_shot_simulator import QubitGKSLNoisyShotSimulator, QubitGKSLShotSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_boson_simulator import QuditGKSLBosonSimulator
from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator, QuditGKSLShotSimulator
from qudit_gksl_simulator import QuditGKSLSimulator

# Non-boson simulators (N=4, d=3)
_NONBOSON_SIMULATOR_MAP: dict[str, type] = {
    "classical": ClassicalGKSLSimulator,
    "qubit": QubitGKSLSimulator,
    "qudit": QuditGKSLSimulator,
}

# Boson simulators (N=2, n_max=1)
_BOSON_SIMULATOR_MAP: dict[str, type] = {
    "classical_boson": ClassicalGKSLBosonSimulator,
    "qubit_boson": QubitGKSLBosonSimulator,
    "qudit_boson": QuditGKSLBosonSimulator,
}

# Shot-based simulators (non-boson, N=4)
_SHOT_SIMULATOR_MAP: dict[str, type] = {
    "qudit_shot": QuditGKSLShotSimulator,
    "qubit_shot": QubitGKSLShotSimulator,
}

_NOISY_SHOT_SIMULATOR_MAP: dict[str, type] = {
    "qudit_noisy_shot": QuditGKSLNoisyShotSimulator,
    "qubit_noisy_shot": QubitGKSLNoisyShotSimulator,
}

ALL_SCENARIO_NAMES = (
    list(_NONBOSON_SIMULATOR_MAP)
    + list(_BOSON_SIMULATOR_MAP)
    + list(_SHOT_SIMULATOR_MAP)
    + list(_NOISY_SHOT_SIMULATOR_MAP)
    + ["unitary"]
)


def _to_float_dict(population: dict[str, Any]) -> dict[str, float]:
    return {
        "N_S0": float(population["N_S0"]),
        "N_T1": float(population["N_T1"]),
        "N_S1": float(population["N_S1"]),
    }


def quantum_fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Compute quantum state fidelity F(rho, sigma) via eigendecomposition."""
    evals, evecs = np.linalg.eigh(rho)
    evals = np.maximum(evals, 0.0)
    sqrt_rho = evecs @ np.diag(np.sqrt(evals)) @ evecs.conj().T
    m = sqrt_rho @ sigma @ sqrt_rho
    m = (m + m.conj().T) / 2
    evals_m = np.linalg.eigvalsh(m)
    evals_m = np.maximum(evals_m, 0.0)
    return float(np.real(np.sum(np.sqrt(evals_m))) ** 2)


def _validate_result(
    result: dict[str, Any],
    params: GKSLPhysicalParameters,
    n_steps: int,
) -> dict[str, Any]:
    """Validate a simulation result and return metrics."""
    max_trace_error = 0.0
    max_particle_error = 0.0

    for pop, trace_val in zip(result["populations"], result["trace"]):
        max_trace_error = max(max_trace_error, abs(float(trace_val) - 1.0))

        pop_float = _to_float_dict(pop)
        validate_particle_conservation(
            pop_float, N_molecules=float(params.N_molecules), tolerance=1e-6
        )
        particle_total = pop_float["N_S0"] + pop_float["N_T1"] + pop_float["N_S1"]
        max_particle_error = max(
            max_particle_error, abs(particle_total - float(params.N_molecules))
        )

    density_validation = validate_density_matrix(result["rho_final"], step=n_steps)

    return {
        "final_density_validation": {
            "trace": float(density_validation["trace"]),
            "hermiticity_error": float(density_validation["hermiticity_error"]),
            "min_eigenvalue": float(density_validation["min_eigenvalue"]),
            "entropy": float(density_validation["entropy"]),
            "purity": float(density_validation["purity"]),
            "valid": bool(density_validation["valid"]),
        },
        "max_trace_error": max_trace_error,
        "max_particle_error": max_particle_error,
    }


def _build_scenario_report(
    scenario: str,
    result: dict[str, Any],
    validation: dict[str, Any],
) -> dict[str, Any]:
    """Build a scenario report dict from simulation result and validation."""
    report: dict[str, Any] = {
        "scenario": scenario,
        "method": result["method"],
        "elapsed_time": float(result["elapsed_time"]),
        "times": [float(t) for t in result["times"]],
        "trace": [float(v) for v in result["trace"]],
        "entropy": [float(v) for v in result["entropy"]],
        "purity": [float(v) for v in result["purity"]],
        "populations": [_to_float_dict(p) for p in result["populations"]],
    }
    report.update(validation)

    # Include shot-based specific fields
    if "n_shots" in result:
        report["n_shots"] = int(result["n_shots"])
    if "seed" in result:
        report["seed"] = result["seed"]
    if "counts" in result:
        report["counts_top10"] = dict(
            sorted(result["counts"].items(), key=lambda x: -x[1])[:10]
        )
    if "noise_params" in result:
        report["noise_params"] = {
            k: float(v) if v is not None else None
            for k, v in result["noise_params"].items()
        }

    return report


def _run_single_scenario(
    scenario: str,
    params: GKSLPhysicalParameters,
    t_max: float,
    n_steps: int,
    initial_state: str,
    n_shots: int = 1000,
    seed: int = 42,
) -> tuple[dict[str, Any], np.ndarray]:
    """Run a single scenario and return (report_dict, rho_final).

    Returns the raw rho_final separately for cross-scenario fidelity computation.
    """
    # Instantiate simulator
    if scenario in _NONBOSON_SIMULATOR_MAP:
        simulator = _NONBOSON_SIMULATOR_MAP[scenario](params)
        result = simulator.simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state
        )
    elif scenario in _BOSON_SIMULATOR_MAP:
        simulator = _BOSON_SIMULATOR_MAP[scenario](params)
        result = simulator.simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state
        )
    elif scenario in _SHOT_SIMULATOR_MAP:
        simulator = _SHOT_SIMULATOR_MAP[scenario](params)
        result = simulator.simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state,
            n_shots=n_shots, seed=seed,
        )
    elif scenario in _NOISY_SHOT_SIMULATOR_MAP:
        if scenario == "qudit_noisy_shot":
            simulator = QuditGKSLNoisyShotSimulator(
                params, p_depol=0.01, p_dephasing=0.005
            )
        elif scenario == "qubit_noisy_shot":
            simulator = QubitGKSLNoisyShotSimulator(
                params, p_depol=0.01, T1=50000.0, t_gate=300.0
            )
        else:
            msg = f"Unknown noisy shot scenario: {scenario}"
            raise ValueError(msg)
        result = simulator.simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state,
            n_shots=n_shots, seed=seed,
        )
    elif scenario == "unitary":
        params_unitary = GKSLPhysicalParameters(
            E_T=params.E_T, E_S=params.E_S, V=params.V,
            N_molecules=params.N_molecules, d=params.d,
            gamma_TTA=0, Gamma_fl=0, Gamma_ph=0,
            k_IC=0, k_ISC_ST=0, k_ISC_TS=0,
        )
        simulator = ClassicalGKSLSimulator(params_unitary)
        result = simulator.simulate(
            t_max=t_max, n_steps=n_steps, initial_state=initial_state
        )
    else:
        msg = f"Unknown scenario: {scenario}"
        raise ValueError(msg)

    rho_final = result["rho_final"]
    validation = _validate_result(result, params, n_steps)
    report = _build_scenario_report(scenario, result, validation)

    return report, rho_final


def _compute_fidelities(
    rho_map: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    """Compute pairwise fidelities for scenarios with matching Hilbert space dims."""
    fidelities: list[dict[str, Any]] = []
    names = list(rho_map.keys())
    for i, name_a in enumerate(names):
        for name_b in names[i + 1 :]:
            rho_a = rho_map[name_a]
            rho_b = rho_map[name_b]
            if rho_a.shape != rho_b.shape:
                continue
            error_msg = None
            try:
                f_val = quantum_fidelity(rho_a, rho_b)
            except Exception as exc:  # noqa: BLE001
                f_val = None
                error_msg = f"{type(exc).__name__}: {exc}"
                print(f"  WARNING: fidelity({name_a}, {name_b}) failed: {error_msg}")
            fidelities.append({
                "scenario_a": name_a,
                "scenario_b": name_b,
                "dim": int(rho_a.shape[0]),
                "fidelity": float(f_val) if f_val is not None else None,
                "error": error_msg,
            })
    return fidelities


def run_verification_workflow(
    output_dir: Path,
    scenarios: tuple[str, ...] = ("classical", "qubit", "qudit"),
    t_max: float = 5.0,
    n_steps: int = 5,
    initial_state: str = "edge_triplet",
    n_shots: int = 1000,
    seed: int = 42,
) -> tuple[Path, Path, dict[str, Any]]:
    """Run all requested scenarios and write JSON/Markdown reports."""
    # Parameters matching the notebook
    params_nonboson = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
    )
    params_boson = GKSLPhysicalParameters(
        N_molecules=2, with_boson=True, n_max=1,
        omega_ph=0.15, g_eph=0.02,
    )

    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    output_dir.mkdir(parents=True, exist_ok=True)

    scenario_reports: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    rho_map: dict[str, np.ndarray] = {}

    for scenario in scenarios:
        print(f"[{scenario}] running...", flush=True)
        try:
            # Select appropriate parameters
            if scenario in _BOSON_SIMULATOR_MAP:
                params = params_boson
            elif scenario == "unitary":
                params = params_nonboson
            else:
                params = params_nonboson

            report, rho_final = _run_single_scenario(
                scenario=scenario,
                params=params,
                t_max=t_max,
                n_steps=n_steps,
                initial_state=initial_state,
                n_shots=n_shots,
                seed=seed,
            )
            scenario_reports.append(report)
            rho_map[scenario] = rho_final
            status = "PASS" if report["final_density_validation"]["valid"] else "FAIL"
            print(f"[{scenario}] done ({report['elapsed_time']:.2f}s) -> {status}")
        except Exception as exc:  # noqa: BLE001
            tb = traceback.format_exc()
            errors.append({
                "scenario": scenario,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": tb,
            })
            print(f"[{scenario}] ERROR: {exc}")

    # Cross-scenario fidelity comparison
    fidelities = _compute_fidelities(rho_map)

    # Unitary-specific analysis
    unitary_analysis: dict[str, Any] = {}
    if "unitary" in rho_map and "classical" in rho_map:
        unitary_report = next(
            (r for r in scenario_reports if r["scenario"] == "unitary"), None
        )
        classical_report = next(
            (r for r in scenario_reports if r["scenario"] == "classical"), None
        )
        if unitary_report and classical_report:
            unitary_analysis = {
                "unitary_max_entropy": max(unitary_report["entropy"]),
                "gksl_final_entropy": classical_report["entropy"][-1],
                "comment": (
                    "Unitary limit (all dissipation=0) should have near-zero entropy. "
                    "GKSL with dissipation should show entropy increase over time."
                ),
            }

    all_valid = all(
        r["final_density_validation"]["valid"] for r in scenario_reports
    )
    overall_success = len(errors) == 0 and all_valid

    report_out = {
        "timestamp": timestamp,
        "success": overall_success,
        "config": {
            "scenarios": list(scenarios),
            "t_max": t_max,
            "n_steps": n_steps,
            "initial_state": initial_state,
            "n_shots": n_shots,
            "seed": seed,
            "params_nonboson": params_nonboson.to_dict(),
            "params_boson": params_boson.to_dict(),
        },
        "results": scenario_reports,
        "fidelities": fidelities,
        "unitary_analysis": unitary_analysis,
        "errors": errors,
    }

    json_path = output_dir / f"tta_uc_gksl_verification_{timestamp}.json"
    markdown_path = output_dir / f"tta_uc_gksl_verification_{timestamp}.md"

    json_path.write_text(
        json.dumps(report_out, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    md_lines = _build_markdown_report(report_out, scenarios)
    markdown_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"\nJSON report: {json_path}")
    print(f"Markdown report: {markdown_path}")
    print(f"Overall: {'PASS' if overall_success else 'FAIL'}")

    return json_path, markdown_path, report_out


def _build_markdown_report(
    report: dict[str, Any],
    scenarios: tuple[str, ...],
) -> list[str]:
    """Build Markdown report lines."""
    md = [
        "# TTA-UC GKSL-Lindblad 検証結果",
        "",
        f"- 実行時刻(UTC): {report['timestamp']}",
        f"- 総合判定: **{'PASS' if report['success'] else 'FAIL'}**",
        f"- シナリオ: {', '.join(scenarios)}",
        f"- t_max: {report['config']['t_max']}",
        f"- n_steps: {report['config']['n_steps']}",
        f"- initial_state: {report['config']['initial_state']}",
        f"- n_shots (shot系): {report['config']['n_shots']}",
        f"- seed: {report['config']['seed']}",
        "",
        "## 評価基準",
        "",
        "- max_trace_error: Tr(ρ)=1からの最大偏差",
        "- max_particle_error: N_S0+N_T1+N_S1=N_moleculesからの最大偏差",
        "- density_validation: 密度行列の物理的妥当性(エルミート性、正値性、トレース)",
        "- fidelity: 異シナリオ間の量子状態忠実度",
        "",
        "## シナリオ別結果",
        "",
    ]

    for sr in report["results"]:
        v = sr["final_density_validation"]
        status = "PASS ✅" if v["valid"] else "FAIL ❌"
        md.extend([
            f"### {sr['scenario']}",
            "",
            f"| 項目 | 値 |",
            f"|------|-----|",
            f"| method | {sr['method']} |",
            f"| elapsed_time | {sr['elapsed_time']:.6f} s |",
            f"| max_trace_error | {sr['max_trace_error']:.3e} |",
            f"| max_particle_error | {sr['max_particle_error']:.3e} |",
            f"| final_trace | {v['trace']:.12f} |",
            f"| final_hermiticity_error | {v['hermiticity_error']:.3e} |",
            f"| final_min_eigenvalue | {v['min_eigenvalue']:.3e} |",
            f"| final_entropy | {v['entropy']:.12f} |",
            f"| final_purity | {v['purity']:.12f} |",
            f"| density_validation | {status} |",
        ])
        if "n_shots" in sr:
            md.append(f"| n_shots | {sr['n_shots']} |")
        if "noise_params" in sr:
            for k, val in sr["noise_params"].items():
                md.append(f"| noise_{k} | {val} |")
        if "counts_top10" in sr:
            top = ", ".join(f"{k}:{v}" for k, v in sr["counts_top10"].items())
            md.append(f"| counts_top10 | {top} |")

        # Population at initial, mid, and final times
        pops = sr["populations"]
        if len(pops) >= 3:
            mid_idx = len(pops) // 2
            for label, idx in [("initial", 0), ("mid", mid_idx), ("final", -1)]:
                p = pops[idx]
                md.append(
                    f"| pop_{label} | S0={p['N_S0']:.4f}, T1={p['N_T1']:.4f}, S1={p['N_S1']:.4f} |"
                )
        md.append("")

    # Fidelity comparison
    if report.get("fidelities"):
        md.extend(["## シナリオ間忠実度 (Fidelity)", ""])
        md.extend([
            "| Scenario A | Scenario B | dim | Fidelity |",
            "|------------|------------|-----|----------|",
        ])
        for f in report["fidelities"]:
            f_str = f"{f['fidelity']:.6f}" if f["fidelity"] is not None else "N/A"
            err_note = f" ({f['error']})" if f.get("error") else ""
            md.append(f"| {f['scenario_a']} | {f['scenario_b']} | {f['dim']} | {f_str}{err_note} |")
        has_na = any(f["fidelity"] is None for f in report["fidelities"])
        if has_na:
            md.extend([
                "",
                "> N/A: 忠実度計算が失敗した組合せ。原因は数値的不安定性や計算エラー。",
            ])
        md.append("")

    # Unitary analysis
    if report.get("unitary_analysis"):
        ua = report["unitary_analysis"]
        md.extend([
            "## Unitary vs GKSL 比較",
            "",
            f"- Unitary最大エントロピー: {ua.get('unitary_max_entropy', 'N/A')}",
            f"- GKSL最終エントロピー: {ua.get('gksl_final_entropy', 'N/A')}",
            f"- 備考: {ua.get('comment', '')}",
            "",
        ])

    # Errors
    if report.get("errors"):
        md.extend(["## エラー", ""])
        for err in report["errors"]:
            md.extend([
                f"### {err['scenario']}",
                "",
                f"```",
                err["error"],
                f"```",
                "",
            ])
            if "traceback" in err:
                md.extend([
                    "<details><summary>Traceback</summary>",
                    "",
                    "```",
                    err["traceback"],
                    "```",
                    "",
                    "</details>",
                    "",
                ])

    return md


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="TTA-UC GKSL-Lindblad検証スクリプト")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "developing" / "verification_results",
        help="検証結果(JSON/Markdown)の出力先ディレクトリ",
    )
    parser.add_argument("--t-max", type=float, default=5.0)
    parser.add_argument("--n-steps", type=int, default=5)
    parser.add_argument("--n-shots", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--initial-state",
        choices=("edge_triplet", "all_triplet", "all_singlet"),
        default="edge_triplet",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=ALL_SCENARIO_NAMES,
        default=ALL_SCENARIO_NAMES,
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    _, _, report = run_verification_workflow(
        output_dir=args.output_dir,
        scenarios=tuple(args.scenarios),
        t_max=args.t_max,
        n_steps=args.n_steps,
        initial_state=args.initial_state,
        n_shots=args.n_shots,
        seed=args.seed,
    )
    return 0 if report["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
