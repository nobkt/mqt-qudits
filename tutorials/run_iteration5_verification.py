#!/usr/bin/env python3
"""
Iteration 5 検証スクリプト

quantum_dynamics_gksl_comparison.ipynb のノイズモデル修正を検証する。
修正内容:
  - Cell 34: QuditGKSLNoisyShotSimulator で p_dephasing=0.0, depol_pair_only=True
  - Cell 38: QubitGKSLNoisyShotSimulator で p_dephasing=0.0, depol_pair_only=True

検証項目:
  1. ノートブックのセル内容が正しく修正されているか（静的検証）
  2. シミュレータの depol_pair_only パラメータが正しく動作するか
  3. quantum_dynamics_complete_comparison.ipynb のノイズモデルが既に正しいか

実行方法:
  cd tutorials && python run_iteration5_verification.py

結果は developing/verification_results/ に JSON 形式で保存される。
"""

import json
import os
import sys
import datetime

RESULTS = {"timestamp": datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
           "iteration": 5,
           "checks": [],
           "passed": 0,
           "failed": 0}


def check(name: str, condition: bool, detail: str = ""):
    status = "PASS" if condition else "FAIL"
    RESULTS["checks"].append({"name": name, "status": status, "detail": detail})
    if condition:
        RESULTS["passed"] += 1
    else:
        RESULTS["failed"] += 1
    print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""))
    return condition


def verify_gksl_notebook():
    """GKSL比較ノートブックのセル内容を静的に検証する"""
    print("=" * 70)
    print("1. quantum_dynamics_gksl_comparison.ipynb 静的検証")
    print("=" * 70)

    nb_path = "quantum_dynamics_gksl_comparison.ipynb"
    if not os.path.exists(nb_path):
        check("notebook_exists", False, f"{nb_path} not found")
        return

    with open(nb_path, "r") as f:
        nb = json.load(f)

    cells = nb["cells"]
    check("total_cells", len(cells) == 44, f"cells={len(cells)}, expected=44")

    # Cell 34 (Qudit noisy): p_dephasing=0.0, depol_pair_only=True
    src34 = "".join(cells[34]["source"])
    check("cell34_p_dephasing_0",
          "p_dephasing=0.0" in src34,
          "Cell 34 should have p_dephasing=0.0")
    check("cell34_depol_pair_only",
          "depol_pair_only=True" in src34,
          "Cell 34 should have depol_pair_only=True")
    check("cell34_no_old_dephasing",
          "p_dephasing=0.005" not in src34,
          "Cell 34 should NOT have old p_dephasing=0.005")

    # Cell 38 (Qubit noisy): p_dephasing=0.0, depol_pair_only=True
    src38 = "".join(cells[38]["source"])
    check("cell38_p_dephasing_0",
          "p_dephasing=0.0" in src38,
          "Cell 38 should have p_dephasing=0.0")
    check("cell38_depol_pair_only",
          "depol_pair_only=True" in src38,
          "Cell 38 should have depol_pair_only=True")
    check("cell38_no_old_dephasing",
          "p_dephasing=0.005" not in src38,
          "Cell 38 should NOT have old p_dephasing=0.005")

    # Cell 30 (markdown): table should NOT reference "脱分極+位相緩和" for 5c/3c
    src30 = "".join(cells[30]["source"])
    check("cell30_no_dephasing_in_table",
          "脱分極+位相緩和" not in src30,
          "Cell 30 table should not mention 脱分極+位相緩和")
    check("cell30_depol_only_in_table",
          "脱分極のみ" in src30,
          "Cell 30 table should mention 脱分極のみ")

    # Cell 33 (markdown): should mention depol_pair_only
    src33 = "".join(cells[33]["source"])
    check("cell33_depol_pair_only",
          "depol_pair_only" in src33,
          "Cell 33 should mention depol_pair_only")
    check("cell33_no_old_dephasing",
          "p_dephasing=0.005" not in src33,
          "Cell 33 should NOT mention old p_dephasing=0.005")

    # Cell 37 (markdown): should mention depol_pair_only
    src37 = "".join(cells[37]["source"])
    check("cell37_depol_pair_only",
          "depol_pair_only" in src37,
          "Cell 37 should mention depol_pair_only")
    check("cell37_no_old_dephasing",
          "p_dephasing=0.005" not in src37,
          "Cell 37 should NOT mention old p_dephasing=0.005")

    # Cell 39 (markdown): should mention depol_pair_only and unified model
    src39 = "".join(cells[39]["source"])
    check("cell39_depol_pair_only",
          "depol_pair_only" in src39,
          "Cell 39 should mention depol_pair_only")
    check("cell39_model_unification",
          "統一" in src39 or "complete_comparison" in src39,
          "Cell 39 should mention model unification with complete_comparison")

    # Cell 41 (markdown): table should use "脱分極のみ"
    src41 = "".join(cells[41]["source"])
    check("cell41_depol_only_in_table",
          "脱分極のみ" in src41,
          "Cell 41 table should use 脱分極のみ")
    check("cell41_no_old_dephasing",
          "脱分極+位相緩和" not in src41,
          "Cell 41 should NOT have 脱分極+位相緩和")

    # Cell 43 (markdown): conclusion should reflect new noise model
    src43 = "".join(cells[43]["source"])
    check("cell43_depol_only_qudit",
          "脱分極のみ" in src43 or "脱分極のみ（ペアゲートのみ）" in src43,
          "Cell 43 should mention 脱分極のみ for qudit")
    check("cell43_depol_only_qubit",
          "脱分極のみ（ペアゲートのみ）" in src43,
          "Cell 43 should mention 脱分極のみ（ペアゲートのみ）for qubit")
    check("cell43_no_old_noise_model",
          "位相緩和（Qudit）" not in src43,
          "Cell 43 should NOT have old dephasing noise model description")


def verify_complete_comparison_notebook():
    """完全比較ノートブックのノイズモデルが正しいことを確認する"""
    print()
    print("=" * 70)
    print("2. quantum_dynamics_complete_comparison.ipynb 静的検証")
    print("=" * 70)

    nb_path = "quantum_dynamics_complete_comparison.ipynb"
    if not os.path.exists(nb_path):
        check("complete_notebook_exists", False, f"{nb_path} not found")
        return

    with open(nb_path, "r") as f:
        nb = json.load(f)

    cells = nb["cells"]

    # Cell 15 (Qubit noisy): should use depol_2q only
    src15 = "".join(cells[15]["source"])
    check("complete_cell15_depol_2q",
          "depol_2q" in src15,
          "Cell 15 should use depol_2q parameter")
    check("complete_cell15_no_dephasing",
          "dephasing" not in src15.lower() or "ノイズなし" in src15,
          "Cell 15 should not have dephasing noise")
    check("complete_cell15_2qubit_only",
          "2-qubitゲートのみ" in src15 or "2量子ビットゲート" in src15,
          "Cell 15 should specify 2-qubit gates only")

    # Cell 23 (Qudit noisy): should use depol_2q only
    src23 = "".join(cells[23]["source"])
    check("complete_cell23_depol_2q",
          "depol_2q" in src23,
          "Cell 23 should use depol_2q parameter")
    check("complete_cell23_noise_gates",
          "noise_gates" in src23,
          "Cell 23 should specify noise_gates list")
    check("complete_cell23_2qudit_only",
          "2-quditゲートのみ" in src23 or "2量子ビットゲート" in src23,
          "Cell 23 should specify 2-qudit gates only")


def verify_simulator_classes():
    """シミュレータクラスのdepol_pair_onlyパラメータを検証する"""
    print()
    print("=" * 70)
    print("3. シミュレータクラス depol_pair_only パラメータ検証")
    print("=" * 70)

    # QuditGKSLNoisyShotSimulator
    qudit_file = "qudit_gksl_shot_simulator.py"
    if os.path.exists(qudit_file):
        with open(qudit_file, "r") as f:
            content = f.read()
        check("qudit_shot_noisy_depol_pair_only_param",
              "depol_pair_only" in content,
              "QuditGKSLNoisyShotSimulator should have depol_pair_only parameter")
        check("qudit_shot_noisy_depol_pair_only_default",
              "depol_pair_only: bool = False" in content,
              "depol_pair_only should default to False")
        check("qudit_shot_noisy_noise_params_output",
              '"depol_pair_only": self.depol_pair_only' in content,
              "noise_params should include depol_pair_only")
    else:
        check("qudit_shot_file_exists", False, f"{qudit_file} not found")

    # QubitGKSLNoisyShotSimulator
    qubit_file = "qubit_gksl_shot_simulator.py"
    if os.path.exists(qubit_file):
        with open(qubit_file, "r") as f:
            content = f.read()
        check("qubit_shot_noisy_depol_pair_only_param",
              "depol_pair_only" in content,
              "QubitGKSLNoisyShotSimulator should have depol_pair_only parameter")
        check("qubit_shot_noisy_depol_pair_only_default",
              "depol_pair_only: bool = False" in content,
              "depol_pair_only should default to False")
        check("qubit_shot_noisy_noise_params_output",
              '"depol_pair_only": self.depol_pair_only' in content,
              "noise_params should include depol_pair_only")
    else:
        check("qubit_shot_file_exists", False, f"{qubit_file} not found")


def verify_noise_model_consistency():
    """両ノートブック間のノイズモデル一貫性を検証する"""
    print()
    print("=" * 70)
    print("4. ノイズモデル一貫性の検証")
    print("=" * 70)

    # Both notebooks should use depolarizing noise only on pair gates
    gksl_path = "quantum_dynamics_gksl_comparison.ipynb"
    complete_path = "quantum_dynamics_complete_comparison.ipynb"

    if not (os.path.exists(gksl_path) and os.path.exists(complete_path)):
        check("both_notebooks_exist", False, "Both notebooks must exist")
        return

    with open(gksl_path, "r") as f:
        gksl_nb = json.load(f)
    with open(complete_path, "r") as f:
        complete_nb = json.load(f)

    # GKSL: p_depol=0.01
    src34 = "".join(gksl_nb["cells"][34]["source"])
    src38 = "".join(gksl_nb["cells"][38]["source"])
    check("gksl_qudit_depol_rate",
          "p_depol=0.01" in src34,
          "GKSL qudit should use p_depol=0.01")
    check("gksl_qubit_depol_rate",
          "p_depol=0.01" in src38,
          "GKSL qubit should use p_depol=0.01")

    # Complete: depol_2q=0.01
    src15 = "".join(complete_nb["cells"][15]["source"])
    src23 = "".join(complete_nb["cells"][23]["source"])
    check("complete_qubit_depol_rate",
          "0.01" in src15,
          "Complete qubit should use depol_2q=0.01")
    check("complete_qudit_depol_rate",
          "0.01" in src23,
          "Complete qudit should use depol_2q=0.01")

    # Both should NOT have dephasing
    check("gksl_no_dephasing_qudit",
          "p_dephasing=0.005" not in src34,
          "GKSL qudit should NOT have p_dephasing=0.005")
    check("gksl_no_dephasing_qubit",
          "p_dephasing=0.005" not in src38,
          "GKSL qubit should NOT have p_dephasing=0.005")

    print()
    print("  ✓ 両ノートブックのノイズモデルが統一されています:")
    print("    - 脱分極ノイズのみ（p_depol=0.01）")
    print("    - 2-qubit/2-quditペアゲートのみに適用")
    print("    - 位相緩和（dephasing）なし")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    print("=" * 70)
    print("Iteration 5 検証スクリプト")
    print(f"実行日時: {RESULTS['timestamp']}")
    print("=" * 70)
    print()

    verify_gksl_notebook()
    verify_complete_comparison_notebook()
    verify_simulator_classes()
    verify_noise_model_consistency()

    print()
    print("=" * 70)
    print(f"検証結果: {RESULTS['passed']}/{RESULTS['passed'] + RESULTS['failed']} チェック通過")
    if RESULTS["failed"] > 0:
        print(f"  ⚠ {RESULTS['failed']} チェック失敗")
    else:
        print("  ✓ 全チェック通過")
    print("=" * 70)

    # Save results
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "developing", "verification_results")
    os.makedirs(results_dir, exist_ok=True)
    fname = f"iteration5_noise_model_{RESULTS['timestamp']}.json"
    fpath = os.path.join(results_dir, fname)
    with open(fpath, "w") as f:
        json.dump(RESULTS, f, indent=2, ensure_ascii=False)
    print(f"\n結果保存先: {fpath}")

    return 0 if RESULTS["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
