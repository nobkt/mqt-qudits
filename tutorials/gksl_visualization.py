from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


def plot_population_dynamics(result, title=None, save_path=None):
    """Plot N_S0, N_T1, N_S1 vs time from a single simulation result."""
    times = result["times"]
    N_S0 = [p["N_S0"] for p in result["populations"]]
    N_T1 = [p["N_T1"] for p in result["populations"]]
    N_S1 = [p["N_S1"] for p in result["populations"]]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, N_S0, "o-", label="$N_{S_0}$ (ground)", linewidth=2, markersize=4)
    ax.plot(times, N_T1, "s-", label="$N_{T_1}$ (triplet)", linewidth=2, markersize=4)
    ax.plot(times, N_S1, "^-", label="$N_{S_1}$ (singlet)", linewidth=2, markersize=4)
    ax.set_xlabel("Time ($\\hbar$/eV)", fontsize=14)
    ax.set_ylabel("Population", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    if title:
        ax.set_title(title, fontsize=16)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_entropy_and_purity(result, title=None, save_path=None):
    """Plot entropy and purity evolution side by side."""
    times = result["times"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(times, result["entropy"], "o-", linewidth=2, markersize=4)
    ax1.set_xlabel("Time ($\\hbar$/eV)")
    ax1.set_ylabel("von Neumann Entropy")
    ax1.grid(True, alpha=0.3)
    ax1.set_title("Entropy Evolution")

    ax2.plot(times, result["purity"], "s-", linewidth=2, markersize=4, color="orange")
    ax2.set_xlabel("Time ($\\hbar$/eV)")
    ax2.set_ylabel("Purity $\\mathrm{Tr}[\\rho^2]$")
    ax2.grid(True, alpha=0.3)
    ax2.set_title("Purity Evolution")

    if title:
        fig.suptitle(title, fontsize=16)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def compare_multiple_scenarios(results_dict, title=None, save_path=None):
    """Compare N_S0, N_T1, N_S1 across multiple scenarios.

    results_dict: {scenario_name: result_dict}
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    labels = ["$N_{S_0}$", "$N_{T_1}$", "$N_{S_1}$"]
    keys = ["N_S0", "N_T1", "N_S1"]

    for name, result in results_dict.items():
        times = result["times"]
        for ax, key in zip(axes, keys):
            vals = [p[key] for p in result["populations"]]
            ax.plot(times, vals, label=name, linewidth=2)

    for ax, label in zip(axes, labels):
        ax.set_xlabel("Time ($\\hbar$/eV)")
        ax.set_ylabel(label)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    if title:
        fig.suptitle(title, fontsize=16)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_gksl_comparison(result_unitary, result_gksl, title=None, save_path=None):
    """Compare unitary (closed system) vs GKSL (open system) dynamics."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    for result, style, prefix in [(result_gksl, "-", "GKSL"), (result_unitary, "--", "Unitary")]:
        times = result["times"]
        ax1.plot(
            times,
            [p["N_S0"] for p in result["populations"]],
            f"b{style}",
            label=f"{prefix} $N_{{S_0}}$",
            linewidth=2,
        )
        ax1.plot(
            times,
            [p["N_T1"] for p in result["populations"]],
            f"r{style}",
            label=f"{prefix} $N_{{T_1}}$",
            linewidth=2,
        )
        ax1.plot(
            times,
            [p["N_S1"] for p in result["populations"]],
            f"g{style}",
            label=f"{prefix} $N_{{S_1}}$",
            linewidth=2,
        )

    ax1.set_xlabel("Time ($\\hbar$/eV)")
    ax1.set_ylabel("Population")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_title("Unitary vs GKSL Dynamics")

    times_g = result_gksl["times"]
    ax2.plot(times_g, result_gksl["entropy"], "o-", label="Entropy", linewidth=2, markersize=3)
    ax2_twin = ax2.twinx()
    ax2_twin.plot(
        times_g,
        result_gksl["purity"],
        "s-",
        color="orange",
        label="Purity",
        linewidth=2,
        markersize=3,
    )
    ax2.set_xlabel("Time ($\\hbar$/eV)")
    ax2.set_ylabel("Entropy")
    ax2_twin.set_ylabel("Purity")
    ax2.legend(loc="upper left")
    ax2_twin.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)
    ax2.set_title("GKSL: Entropy & Purity")

    if title:
        fig.suptitle(title, fontsize=16)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_trace_conservation(result, title=None, save_path=None):
    """Plot trace deviation from 1.0 over time."""
    fig, ax = plt.subplots(figsize=(10, 4))
    times = result["times"]
    trace_dev = [abs(t - 1.0) for t in result["trace"]]
    ax.semilogy(times, trace_dev, "o-", linewidth=2, markersize=4)
    ax.set_xlabel("Time ($\\hbar$/eV)")
    ax.set_ylabel("$|\\mathrm{Tr}[\\rho] - 1|$")
    ax.grid(True, alpha=0.3)
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Trace Conservation")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig
