from __future__ import annotations

from tutorials.gksl_visualization import (
    plot_6scenario_comparison,
    plot_entropy_dynamics,
    plot_per_molecule_populations,
)


def _create_sample_gksl_result() -> dict:
    return {
        "times": [0.0, 1.0],
        "populations": [
            {"N_S0": 2.0, "N_T1": 1.0, "N_S1": 1.0},
            {"N_S0": 2.2, "N_T1": 0.9, "N_S1": 0.9},
        ],
        "entropy": [0.0, 0.1],
    }


def test_plot_per_molecule_populations_returns_figure() -> None:
    result = _create_sample_gksl_result()
    result["per_molecule_populations"] = [
        {"S0_per_mol": [1.0, 1.0], "T1_per_mol": [0.5, 0.5], "S1_per_mol": [0.5, 0.5]},
        {"S0_per_mol": [1.1, 1.1], "T1_per_mol": [0.45, 0.45], "S1_per_mol": [0.45, 0.45]},
    ]
    fig = plot_per_molecule_populations(result, title="per-molecule")
    assert fig is not None


def test_plot_entropy_dynamics_returns_figure() -> None:
    result_a = _create_sample_gksl_result()
    result_b = _create_sample_gksl_result()
    result_b["entropy"] = [0.0, 0.2]
    fig = plot_entropy_dynamics([result_a, result_b], ["A", "B"], title="entropy")
    assert fig is not None


def test_plot_6scenario_comparison_returns_figure() -> None:
    result = _create_sample_gksl_result()
    fig = plot_6scenario_comparison({"Classical NB": result, "Qudit NB": result}, title="6 scenarios")
    assert fig is not None
