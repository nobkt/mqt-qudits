#!/usr/bin/env python3
"""TTA-UC GKSL iteration 11 verification: forbidden-state leakage diagnostics.

This script provides detailed numerical analysis of the qubit noisy shot
simulator's forbidden-state leakage issue identified in iteration 11.

Outputs:
  1. Trace deficit time series for qubit_noisy_shot
  2. Normalized vs raw fidelity comparison
  3. Forbidden-state leakage rate quantification
  4. Fair comparison: qubit_noisy vs qudit_noisy after normalization

Results are written to developing/verification_results/ as JSON and Markdown.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "."))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from classical_gksl_simulator import ClassicalGKSLSimulator
from gksl_physical_parameters import GKSLPhysicalParameters
from qubit_gksl_shot_simulator import QubitGKSLNoisyShotSimulator, QubitGKSLShotSimulator
from qubit_gksl_simulator import QubitGKSLSimulator
from qudit_gksl_shot_simulator import QuditGKSLNoisyShotSimulator, QuditGKSLShotSimulator
from qudit_gksl_simulator import QuditGKSLSimulator


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
    """Normalize density matrix to trace 1 for fidelity computation."""
    tr = float(np.real(np.trace(rho)))
    return rho / tr if tr > 1e-10 else rho


def main() -> int:
    timestamp = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ")
    output_dir = Path(__file__).resolve().parents[1] / "developing" / "verification_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parameters matching the notebook
    params = GKSLPhysicalParameters(
        E_T=1.5, E_S=3.0, V=0.1,
        gamma_TTA=0.05, Gamma_fl=0.01, Gamma_ph=1e-6,
        k_IC=0.005, k_ISC_ST=0.003, k_ISC_TS=1e-5,
    )

    t_max = 5.0
    n_steps = 5
    n_shots = 1000
    seed = 42

    print("=" * 70)
    print("TTA-UC GKSL Iteration 11: Forbidden-State Leakage Diagnostics")
    print("=" * 70)

    # --- Run classical reference ---
    print("\n[1/6] Classical reference...", flush=True)
    sim_classical = ClassicalGKSLSimulator(params)
    result_classical = sim_classical.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
    )
    print(f"  Done ({result_classical['elapsed_time']:.2f}s)")

    # --- Run qubit DM ---
    print("[2/6] Qubit DM (statevector)...", flush=True)
    sim_qubit = QubitGKSLSimulator(params)
    result_qubit = sim_qubit.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
    )
    print(f"  Done ({result_qubit['elapsed_time']:.2f}s)")

    # --- Run qudit DM ---
    print("[3/6] Qudit DM (statevector)...", flush=True)
    sim_qudit = QuditGKSLSimulator(params)
    result_qudit = sim_qudit.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet"
    )
    print(f"  Done ({result_qudit['elapsed_time']:.2f}s)")

    # --- Run qudit noisy shot ---
    print("[4/6] Qudit noisy shot...", flush=True)
    sim_qudit_noisy = QuditGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.005)
    result_qudit_noisy = sim_qudit_noisy.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet",
        n_shots=n_shots, seed=seed,
    )
    print(f"  Done ({result_qudit_noisy['elapsed_time']:.2f}s)")

    # --- Run qubit noisy shot ---
    print("[5/6] Qubit noisy shot...", flush=True)
    sim_qubit_noisy = QubitGKSLNoisyShotSimulator(params, p_depol=0.01, p_dephasing=0.005)
    result_qubit_noisy = sim_qubit_noisy.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet",
        n_shots=n_shots, seed=seed,
    )
    print(f"  Done ({result_qubit_noisy['elapsed_time']:.2f}s)")

    # --- Run qubit shot (no noise) for baseline ---
    print("[6/6] Qubit shot (no noise)...", flush=True)
    sim_qubit_shot = QubitGKSLShotSimulator(params)
    result_qubit_shot = sim_qubit_shot.simulate(
        t_max=t_max, n_steps=n_steps, initial_state="edge_triplet",
        n_shots=n_shots, seed=seed,
    )
    print(f"  Done ({result_qubit_shot['elapsed_time']:.2f}s)")

    # =====================================================================
    # Analysis 1: Trace deficit time series
    # =====================================================================
    print("\n" + "=" * 70)
    print("Analysis 1: Trace Deficit Time Series (Qubit Noisy Shot)")
    print("=" * 70)

    trace_qubit_noisy = result_qubit_noisy["trace"]
    trace_qubit_shot = result_qubit_shot["trace"]
    trace_qudit_noisy = result_qudit_noisy["trace"]
    times = result_qubit_noisy["times"]

    print(f"\n{'Time':>8} | {'Tr(qubit noisy)':>16} | {'Deficit':>10} | {'Tr(qubit ideal)':>16} | {'Tr(qudit noisy)':>16}")
    print("-" * 80)
    trace_deficit_series = []
    for i in range(len(times)):
        tr_qn = trace_qubit_noisy[i]
        tr_qi = trace_qubit_shot[i]
        tr_dn = trace_qudit_noisy[i]
        deficit = 1.0 - tr_qn
        trace_deficit_series.append({
            "time": float(times[i]),
            "trace_qubit_noisy": float(tr_qn),
            "trace_deficit": float(deficit),
            "trace_qubit_ideal": float(tr_qi),
            "trace_qudit_noisy": float(tr_dn),
        })
        print(f"{times[i]:8.2f} | {tr_qn:16.10f} | {deficit:10.4f} | {tr_qi:16.10f} | {tr_dn:16.10f}")

    # =====================================================================
    # Analysis 2: Forbidden state measurement counts
    # =====================================================================
    print("\n" + "=" * 70)
    print("Analysis 2: Forbidden State Measurement Counts")
    print("=" * 70)

    forbidden_qubit_noisy = result_qubit_noisy.get("forbidden_count", 0)
    forbidden_qubit_ideal = result_qubit_shot.get("forbidden_count", 0)

    print(f"  Qubit noisy shot: {forbidden_qubit_noisy}/{n_shots} ({100*forbidden_qubit_noisy/n_shots:.1f}%)")
    print(f"  Qubit ideal shot: {forbidden_qubit_ideal}/{n_shots} ({100*forbidden_qubit_ideal/n_shots:.1f}%)")

    # =====================================================================
    # Analysis 3: Population conservation
    # =====================================================================
    print("\n" + "=" * 70)
    print("Analysis 3: Population Conservation (Final Step)")
    print("=" * 70)

    for name, result in [
        ("Classical", result_classical),
        ("Qubit DM", result_qubit),
        ("Qudit DM", result_qudit),
        ("Qudit noisy shot", result_qudit_noisy),
        ("Qubit noisy shot", result_qubit_noisy),
        ("Qubit ideal shot", result_qubit_shot),
    ]:
        pop = result["populations"][-1]
        total = pop["N_S0"] + pop["N_T1"] + pop["N_S1"]
        tr = result["trace"][-1]
        print(f"  {name:20s}: S0={pop['N_S0']:.4f}, T1={pop['N_T1']:.4f}, "
              f"S1={pop['N_S1']:.4f}, sum={total:.4f}, Tr={tr:.6f}")

    # =====================================================================
    # Analysis 4: Fidelity comparison (raw vs normalized)
    # =====================================================================
    print("\n" + "=" * 70)
    print("Analysis 4: Fidelity Comparison (Raw vs Normalized)")
    print("=" * 70)

    rho_classical = result_classical["rho_final"]
    rho_qubit = result_qubit["rho_final"]
    rho_qudit = result_qudit["rho_final"]
    rho_qudit_noisy = result_qudit_noisy["rho_final"]
    rho_qubit_noisy = result_qubit_noisy["rho_final"]

    fidelity_results = []

    pairs = [
        ("classical", "qubit_noisy", rho_classical, rho_qubit_noisy),
        ("classical", "qudit_noisy", rho_classical, rho_qudit_noisy),
        ("qubit_DM", "qubit_noisy", rho_qubit, rho_qubit_noisy),
        ("qudit_DM", "qudit_noisy", rho_qudit, rho_qudit_noisy),
        ("qubit_noisy", "qudit_noisy", rho_qubit_noisy, rho_qudit_noisy),
    ]

    print(f"\n{'Pair':>30s} | {'Tr_A':>8} | {'Tr_B':>8} | {'F_raw':>10} | {'F_norm':>10} | {'Delta':>10}")
    print("-" * 90)

    for name_a, name_b, rho_a, rho_b in pairs:
        if rho_a.shape != rho_b.shape:
            continue
        tr_a = float(np.real(np.trace(rho_a)))
        tr_b = float(np.real(np.trace(rho_b)))

        f_raw = quantum_fidelity(rho_a, rho_b)
        f_norm = quantum_fidelity(normalize_rho(rho_a), normalize_rho(rho_b))
        delta = f_norm - f_raw

        fidelity_results.append({
            "pair": f"{name_a} vs {name_b}",
            "trace_a": tr_a,
            "trace_b": tr_b,
            "fidelity_raw": f_raw,
            "fidelity_normalized": f_norm,
            "delta": delta,
        })
        print(f"  {name_a} vs {name_b:>15s} | {tr_a:8.4f} | {tr_b:8.4f} | "
              f"{f_raw:10.6f} | {f_norm:10.6f} | {delta:+10.6f}")

    # =====================================================================
    # Analysis 5: Leakage rate per Trotter step
    # =====================================================================
    print("\n" + "=" * 70)
    print("Analysis 5: Leakage Rate Analysis")
    print("=" * 70)

    trace_vals = trace_qubit_noisy
    leakage_per_step = []
    for i in range(1, len(trace_vals)):
        if trace_vals[i - 1] > 0:
            rate = (trace_vals[i - 1] - trace_vals[i]) / trace_vals[i - 1]
        else:
            rate = 0.0
        leakage_per_step.append(rate)
        print(f"  Step {i}: Tr={trace_vals[i]:.6f}, "
              f"delta_Tr={trace_vals[i]-trace_vals[i-1]:.6f}, "
              f"rate={rate:.4f} ({rate*100:.2f}%)")

    avg_leakage_rate = np.mean(leakage_per_step) if leakage_per_step else 0.0
    print(f"\n  Average leakage rate per step: {avg_leakage_rate:.4f} ({avg_leakage_rate*100:.2f}%)")
    print(f"  Total leakage after {n_steps} steps: {1.0-trace_vals[-1]:.4f} ({(1.0-trace_vals[-1])*100:.2f}%)")

    # =====================================================================
    # Write results
    # =====================================================================
    report = {
        "timestamp": timestamp,
        "iteration": 11,
        "purpose": "Forbidden-state leakage diagnostics",
        "config": {
            "t_max": t_max,
            "n_steps": n_steps,
            "n_shots": n_shots,
            "seed": seed,
            "p_depol": 0.01,
            "p_dephasing": 0.005,
        },
        "trace_deficit_series": trace_deficit_series,
        "forbidden_counts": {
            "qubit_noisy": forbidden_qubit_noisy,
            "qubit_ideal": forbidden_qubit_ideal,
            "total_shots": n_shots,
        },
        "fidelity_comparison": fidelity_results,
        "leakage_analysis": {
            "per_step_rates": leakage_per_step,
            "average_rate": float(avg_leakage_rate),
            "total_leakage": float(1.0 - trace_vals[-1]),
        },
    }

    json_path = output_dir / f"iteration11_leakage_diagnostics_{timestamp}.json"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nJSON report: {json_path}")

    # Markdown report
    md_lines = [
        "# TTA-UC GKSL Iteration 11: 禁止状態リーケージ診断結果",
        "",
        f"- 実行時刻(UTC): {timestamp}",
        f"- t_max: {t_max}, n_steps: {n_steps}, n_shots: {n_shots}, seed: {seed}",
        f"- ノイズ: p_depol={0.01}, p_dephasing={0.005}",
        "",
        "## 1. トレース欠損時系列（Qubit Noisy Shot）",
        "",
        "| Time | Tr(qubit noisy) | Deficit | Tr(qubit ideal) | Tr(qudit noisy) |",
        "|------|----------------|---------|----------------|----------------|",
    ]
    for td in trace_deficit_series:
        md_lines.append(
            f"| {td['time']:.2f} | {td['trace_qubit_noisy']:.10f} | "
            f"{td['trace_deficit']:.4f} | {td['trace_qubit_ideal']:.10f} | "
            f"{td['trace_qudit_noisy']:.10f} |"
        )
    md_lines.extend([
        "",
        "## 2. 禁止状態測定カウント",
        "",
        f"- Qubit noisy shot: {forbidden_qubit_noisy}/{n_shots} ({100*forbidden_qubit_noisy/n_shots:.1f}%)",
        f"- Qubit ideal shot: {forbidden_qubit_ideal}/{n_shots} ({100*forbidden_qubit_ideal/n_shots:.1f}%)",
        "",
        "## 3. 忠実度比較（Raw vs Normalized）",
        "",
        "| Pair | Tr_A | Tr_B | F_raw | F_norm | Delta |",
        "|------|------|------|-------|--------|-------|",
    ])
    for fr in fidelity_results:
        md_lines.append(
            f"| {fr['pair']} | {fr['trace_a']:.4f} | {fr['trace_b']:.4f} | "
            f"{fr['fidelity_raw']:.6f} | {fr['fidelity_normalized']:.6f} | "
            f"{fr['delta']:+.6f} |"
        )
    md_lines.extend([
        "",
        "## 4. リーケージ率分析",
        "",
        f"- 平均リーケージ率/ステップ: {avg_leakage_rate:.4f} ({avg_leakage_rate*100:.2f}%)",
        f"- 総リーケージ（{n_steps}ステップ後）: {1.0-trace_vals[-1]:.4f} ({(1.0-trace_vals[-1])*100:.2f}%)",
        "",
        "## 5. 考察",
        "",
        "### 問題の本質",
        "",
        "Qubitエンコーディング（2-qubit/分子: |00⟩=S0, |01⟩=T1, |10⟩=S1, |11⟩=禁止）では、",
        "2-qubitパウリノイズが物理状態→禁止状態 |11⟩ へのリーケージを引き起こす。",
        "15個の非恒等2-qubitパウリ中12個（80%）がリーケージを発生させる。",
        "",
        "これはqubit量子コンピュータの回路レベルノイズを正確にモデル化したものであり、",
        "バグではなく物理的に正しい挙動である。",
        "",
        "### Qudit noisy との比較",
        "",
        "- **Quditノイズ**: Weyl-Heisenberg演算子はd=3空間内で完結。",
        "  禁止状態が存在しないため、Tr(ρ)=1が厳密保存",
        "- **Qubitノイズ**: 2-qubitパウリ演算子はd=4空間で作用。",
        "  リーケージによりTr(ρ)<1",
        "",
        "F_rawはリーケージの影響で誤解を招く低い値を示すため、",
        "サブノーマライズ行列をノーマライズしてからのF_normが正しい比較指標。",
        "",
    ])

    md_path = output_dir / f"iteration11_leakage_diagnostics_{timestamp}.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Markdown report: {md_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
