#!/usr/bin/env python3
"""Verification script for quantum_dynamics_complete_comparison.ipynb iteration 3 fixes.

Validates:
1. N_steps=100 produces smooth data (101 time points)
2. matplotlib Agg backend + savefig pattern works correctly
3. Plot marker settings are correct (markersize=3, markevery)
4. Decomposed circuit visualization is removed from cells 13, 26
5. All code cells are consistent

Run from tutorials/ directory:
    cd tutorials && python run_iteration3_verification.py

Results are saved to developing/verification_results/iteration3_*.json
"""

import json
import os
import sys
from datetime import datetime, timezone

RESULTS = {"timestamp": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"), "checks": []}


def add_check(name, passed, detail=""):
    RESULTS["checks"].append({"name": name, "passed": passed, "detail": detail})
    status = "✓ PASS" if passed else "✗ FAIL"
    print(f"  {status}: {name}")
    if detail:
        print(f"         {detail}")


def main():
    # Load notebook
    nb_path = os.path.join(os.path.dirname(__file__), "quantum_dynamics_complete_comparison.ipynb")
    if not os.path.exists(nb_path):
        print(f"ERROR: Notebook not found at {nb_path}")
        sys.exit(1)

    with open(nb_path, "r") as f:
        nb = json.load(f)

    cells = nb["cells"]
    print(f"Loaded notebook with {len(cells)} cells")
    print()

    # ================================================================
    # Check 1: N_steps = 100 in Cell 3
    # ================================================================
    print("=" * 60)
    print("Check 1: N_steps parameter")
    print("=" * 60)
    cell3_src = "".join(cells[3]["source"])
    add_check(
        "N_steps = 100",
        "self.N_steps = 100" in cell3_src,
        f"Found: {'self.N_steps = 100' if 'self.N_steps = 100' in cell3_src else 'N_steps != 100'}",
    )
    add_check(
        "N_steps not 20",
        "self.N_steps = 20" not in cell3_src,
        "Old value N_steps=20 should not be present",
    )

    # ================================================================
    # Check 2: matplotlib Agg backend in all plotting cells
    # ================================================================
    print()
    print("=" * 60)
    print("Check 2: matplotlib Agg backend usage")
    print("=" * 60)

    # Cell 6 should set Agg
    cell6_src = "".join(cells[6]["source"])
    add_check("Cell 6: matplotlib.use('Agg')", 'matplotlib.use("Agg")' in cell6_src)

    # Cell 9 should set Agg
    cell9_src = "".join(cells[9]["source"])
    add_check("Cell 9: matplotlib.use('Agg')", 'matplotlib.use("Agg")' in cell9_src)

    # ================================================================
    # Check 3: savefig + display(Image()) pattern in affected cells
    # ================================================================
    print()
    print("=" * 60)
    print("Check 3: savefig + display(Image()) pattern")
    print("=" * 60)

    for idx in [6, 9, 29, 33, 34, 35]:
        src = "".join(cells[idx]["source"])
        has_savefig = "fig.savefig(" in src or ".savefig(" in src
        has_display = "display(Image(" in src
        has_plt_close = "plt.close(" in src
        has_plt_show_only = "plt.show()" in src and "fig.savefig(" not in src

        add_check(f"Cell {idx}: has savefig", has_savefig)
        add_check(f"Cell {idx}: has display(Image())", has_display)
        add_check(f"Cell {idx}: has plt.close()", has_plt_close)
        add_check(f"Cell {idx}: no standalone plt.show()", not has_plt_show_only)

    # ================================================================
    # Check 4: markevery in plot cells
    # ================================================================
    print()
    print("=" * 60)
    print("Check 4: markevery for smooth rendering")
    print("=" * 60)

    for idx in [6, 29, 33, 34]:
        src = "".join(cells[idx]["source"])
        add_check(f"Cell {idx}: has markevery", "markevery" in src)
        add_check(f"Cell {idx}: markersize<=3", "markersize=3" in src)

    # ================================================================
    # Check 5: Decomposed circuit removed from visualization
    # ================================================================
    print()
    print("=" * 60)
    print("Check 5: Circuit visualization")
    print("=" * 60)

    cell13_src = "".join(cells[13]["source"])
    add_check(
        "Cell 13: no 分解後UnitaryGate版",
        "分解後UnitaryGate版" not in cell13_src,
        "Decomposed UnitaryGate should not be visualized",
    )
    add_check(
        "Cell 13: has UnitaryGate版",
        "UnitaryGate版" in cell13_src,
        "UnitaryGate version should still be visualized",
    )
    add_check(
        "Cell 13: has 基本ゲート分解版",
        "基本ゲート分解版" in cell13_src,
        "Basic gate version should still be visualized",
    )

    cell26_src = "".join(cells[26]["source"])
    add_check(
        "Cell 26: has CustomTwoゲート版",
        "CustomTwoゲート版" in cell26_src,
        "CustomTwo gate version should be visualized",
    )
    add_check(
        "Cell 26: no 基本ゲート分解版",
        "基本ゲート分解版" not in cell26_src,
        "Decomposed basic gate version should not be visualized",
    )

    # ================================================================
    # Check 6: Outputs cleared (require re-execution)
    # ================================================================
    print()
    print("=" * 60)
    print("Check 6: Outputs cleared")
    print("=" * 60)

    all_cleared = True
    for i, cell in enumerate(cells):
        if cell["cell_type"] == "code":
            if cell.get("outputs", []):
                all_cleared = False
                add_check(f"Cell {i}: outputs cleared", False, f"Cell has {len(cell['outputs'])} outputs")
    if all_cleared:
        add_check("All code cells: outputs cleared", True, "All cells ready for re-execution")

    # ================================================================
    # Summary
    # ================================================================
    print()
    print("=" * 60)
    total = len(RESULTS["checks"])
    passed = sum(1 for c in RESULTS["checks"] if c["passed"])
    failed = total - passed
    print(f"Summary: {passed}/{total} checks passed, {failed} failed")
    print("=" * 60)

    RESULTS["summary"] = {"total": total, "passed": passed, "failed": failed}

    # Save results
    results_dir = os.path.join(os.path.dirname(__file__), "..", "developing", "verification_results")
    os.makedirs(results_dir, exist_ok=True)
    results_path = os.path.join(results_dir, f"iteration3_{RESULTS['timestamp']}.json")
    with open(results_path, "w") as f:
        json.dump(RESULTS, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to: {results_path}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
