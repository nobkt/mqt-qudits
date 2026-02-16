#!/usr/bin/env python3
"""Run TTA-UC GKSL-Lindblad verification and write file-based reports.

This script intentionally performs only physics-based checks and does not use
heuristic shortcuts or fallback logic.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_physical_parameters import GKSLPhysicalParameters
from gksl_validation import validate_density_matrix, validate_particle_conservation
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_simulator import QuditGKSLSimulator

_SIMULATOR_MAP = {
    "classical": ClassicalGKSLSimulator,
    "qubit": QubitGKSLSimulator,
    "qudit": QuditGKSLSimulator,
}


def _to_float_dict(population: dict[str, Any]) -> dict[str, float]:
    return {
        "N_S0": float(population["N_S0"]),
        "N_T1": float(population["N_T1"]),
        "N_S1": float(population["N_S1"]),
    }


def _run_single_scenario(
    scenario: str,
    params: GKSLPhysicalParameters,
    t_max: float,
    n_steps: int,
    initial_state: str,
) -> dict[str, Any]:
    simulator = _SIMULATOR_MAP[scenario](params)
    result = simulator.simulate(t_max=t_max, n_steps=n_steps, initial_state=initial_state)

    max_trace_error = 0.0
    max_particle_error = 0.0

    for pop, trace in zip(result["populations"], result["trace"]):
        max_trace_error = max(max_trace_error, abs(float(trace) - 1.0))

        pop_float = _to_float_dict(pop)
        validate_particle_conservation(pop_float, N_molecules=float(params.N_molecules), tolerance=1e-6)
        particle_total = pop_float["N_S0"] + pop_float["N_T1"] + pop_float["N_S1"]
        max_particle_error = max(max_particle_error, abs(particle_total - float(params.N_molecules)))

    density_validation = validate_density_matrix(result["rho_final"], step=n_steps)

    return {
        "scenario": scenario,
        "method": result["method"],
        "elapsed_time": float(result["elapsed_time"]),
        "times": [float(t) for t in result["times"]],
        "trace": [float(v) for v in result["trace"]],
        "entropy": [float(v) for v in result["entropy"]],
        "purity": [float(v) for v in result["purity"]],
        "populations": [_to_float_dict(p) for p in result["populations"]],
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


def run_verification_workflow(
    output_dir: Path,
    scenarios: tuple[str, ...] = ("classical", "qubit", "qudit"),
    t_max: float = 5.0,
    n_steps: int = 5,
    initial_state: str = "edge_triplet",
) -> tuple[Path, Path, dict[str, Any]]:
    """Run all requested scenarios and write JSON/Markdown reports."""
    params = GKSLPhysicalParameters()

    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    output_dir.mkdir(parents=True, exist_ok=True)

    scenario_reports: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []

    for scenario in scenarios:
        try:
            scenario_reports.append(
                _run_single_scenario(
                    scenario=scenario,
                    params=params,
                    t_max=t_max,
                    n_steps=n_steps,
                    initial_state=initial_state,
                )
            )
        except Exception as exc:  # noqa: BLE001
            errors.append({"scenario": scenario, "error": str(exc)})

    report = {
        "timestamp": timestamp,
        "success": len(errors) == 0,
        "config": {
            "scenarios": list(scenarios),
            "t_max": t_max,
            "n_steps": n_steps,
            "initial_state": initial_state,
            "params": params.to_dict(),
        },
        "results": scenario_reports,
        "errors": errors,
    }

    json_path = output_dir / f"tta_uc_gksl_verification_{timestamp}.json"
    markdown_path = output_dir / f"tta_uc_gksl_verification_{timestamp}.md"

    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# TTA-UC GKSL-Lindblad 検証結果",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        f"- 成功判定: {'PASS' if report['success'] else 'FAIL'}",
        f"- シナリオ: {', '.join(scenarios)}",
        f"- t_max: {t_max}",
        f"- n_steps: {n_steps}",
        f"- initial_state: {initial_state}",
        "",
        "## シナリオ別結果",
        "",
    ]

    for scenario_report in scenario_reports:
        validation = scenario_report["final_density_validation"]
        md_lines.extend(
            [
                f"### {scenario_report['scenario']}",
                f"- method: {scenario_report['method']}",
                f"- elapsed_time: {scenario_report['elapsed_time']:.6f} s",
                f"- max_trace_error: {scenario_report['max_trace_error']:.3e}",
                f"- max_particle_error: {scenario_report['max_particle_error']:.3e}",
                f"- final_trace: {validation['trace']:.12f}",
                f"- final_hermiticity_error: {validation['hermiticity_error']:.3e}",
                f"- final_min_eigenvalue: {validation['min_eigenvalue']:.3e}",
                f"- final_entropy: {validation['entropy']:.12f}",
                f"- final_purity: {validation['purity']:.12f}",
                f"- density_validation: {'PASS' if validation['valid'] else 'FAIL'}",
                "",
            ]
        )

    if errors:
        md_lines.extend(["## エラー", ""])
        for error in errors:
            md_lines.append(f"- {error['scenario']}: {error['error']}")

    markdown_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    return json_path, markdown_path, report


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
    parser.add_argument(
        "--initial-state",
        choices=("edge_triplet", "all_triplet", "all_singlet"),
        default="edge_triplet",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=tuple(_SIMULATOR_MAP.keys()),
        default=("classical", "qubit", "qudit"),
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
    )
    return 0 if report["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
