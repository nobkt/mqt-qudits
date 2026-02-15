#!/usr/bin/env python3
"""Comprehensive Molecular Hamiltonian Test Suite - PR#44.

このツールは、実際の4分子鎖シミュレーションパラメータを使用して、
様々な条件下でのツールの性能を包括的にテストします。

テスト項目:
1. 複数の時間刻み幅でのテスト (dt = 0.1, 0.5, 1.0, 2.0 fs)
2. 複数のパラメータセット (V, Jの異なる値)
3. 統計的分析（平均ゲート削減率、忠実度分布）
4. 実際の4分子シミュレーションでの推定ゲート数
5. 計算時間の測定とスケーラビリティ分析

制約:
- 既存ソースコード（src/）の修正なし
- ヒューリスティック・Fallback絶対なし
- 数学的に完全に厳密な実装のみ
"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Add tools directory to path
sys.path.insert(0, str(Path(__file__).parent))

from real_hamiltonian_analyzer import DecompositionTestResult, PhysicalParameters, RealHamiltonianAnalyzer


@dataclass
class ComprehensiveTestResult:
    """包括的テスト結果."""

    test_name: str
    parameter_set: str
    time_steps: list[float]
    h_transfer_results: list[DecompositionTestResult]
    h_tta_results: list[DecompositionTestResult]

    def get_average_fidelity(self, matrix_type: str) -> float:
        """平均忠実度を計算."""
        results = self.h_transfer_results if matrix_type == "H_transfer" else self.h_tta_results
        fidelities = [r.fidelity for r in results if r.success]
        return np.mean(fidelities) if fidelities else 0.0

    def get_average_reduction_rate(self, matrix_type: str) -> float:
        """平均ゲート削減率を計算."""
        results = self.h_transfer_results if matrix_type == "H_transfer" else self.h_tta_results
        rates = [r.reduction_rate for r in results if r.success]
        return np.mean(rates) if rates else 0.0

    def get_min_fidelity(self, matrix_type: str) -> float:
        """最小忠実度を計算."""
        results = self.h_transfer_results if matrix_type == "H_transfer" else self.h_tta_results
        fidelities = [r.fidelity for r in results if r.success]
        return np.min(fidelities) if fidelities else 0.0

    def get_success_rate(self, matrix_type: str) -> float:
        """成功率を計算."""
        results = self.h_transfer_results if matrix_type == "H_transfer" else self.h_tta_results
        total = len(results)
        success = sum(1 for r in results if r.success)
        return success / total if total > 0 else 0.0


class ComprehensiveMolecularTester:
    """包括的な分子ハミルトニアンテスター."""

    def __init__(self, tolerance: float = 1e-10) -> None:
        """Args:
        tolerance: 数値許容誤差.
        """
        self.tolerance = tolerance
        self.analyzer = None  # Will be initialized for each parameter set

    def run_time_step_analysis(
        self, params: PhysicalParameters, time_steps: list[float], test_name: str = "Default"
    ) -> ComprehensiveTestResult:
        """複数の時間刻み幅でテストを実行.

        Args:
            params: 物理パラメータ
            time_steps: テストする時間刻み幅のリスト (fs)
            test_name: テスト名

        Returns:
            ComprehensiveTestResult: 包括的テスト結果
        """
        # アナライザー初期化
        self.analyzer = RealHamiltonianAnalyzer(params=params, tolerance=self.tolerance)

        h_transfer_results = []
        h_tta_results = []

        for dt in time_steps:
            # H_transferテスト
            h_transfer = self.analyzer.generate_H_transfer_unitary(dt=dt)
            decomp_transfer = self.analyzer.test_decomposition(h_transfer)
            h_transfer_results.append(decomp_transfer)

            # H_TTAテスト
            h_tta = self.analyzer.generate_H_TTA_unitary(dt=dt)
            decomp_tta = self.analyzer.test_decomposition(h_tta)
            h_tta_results.append(decomp_tta)

        return ComprehensiveTestResult(
            test_name=test_name,
            parameter_set=f"V={params.V[0]:.2f}, J={params.J[0]:.2f}",
            time_steps=time_steps,
            h_transfer_results=h_transfer_results,
            h_tta_results=h_tta_results,
        )

    def estimate_4_molecule_simulation_gates(self, result: ComprehensiveTestResult, n_steps: int = 100) -> dict:
        """4分子鎖シミュレーションの推定ゲート数を計算.

        4分子鎖シミュレーション構成:
        - 3つの隣接ペア: (0,1), (1,2), (2,3)
        - 各ペアに対してH_transferとH_TTAを適用
        - 鈴木トロッター分解: 各ステップで 2 × (H_0 + H_transfer + H_TTA)

        Args:
            result: 包括的テスト結果
            n_steps: シミュレーションステップ数

        Returns:
            Dict: 推定ゲート数の詳細
        """
        # 最後の時間刻み幅の結果を使用（典型的な値）
        h_transfer_last = result.h_transfer_results[-1]
        h_tta_last = result.h_tta_results[-1]

        # ペアあたりのゲート数
        gates_per_pair_transfer_v1 = h_transfer_last.gate_count_v1
        gates_per_pair_transfer_v2 = h_transfer_last.gate_count_v2
        gates_per_pair_tta_v1 = h_tta_last.gate_count_v1
        gates_per_pair_tta_v2 = h_tta_last.gate_count_v2

        # 3ペアでのゲート数
        n_pairs = 3
        gates_transfer_v1 = gates_per_pair_transfer_v1 * n_pairs
        gates_transfer_v2 = gates_per_pair_transfer_v2 * n_pairs
        gates_tta_v1 = gates_per_pair_tta_v1 * n_pairs
        gates_tta_v2 = gates_per_pair_tta_v2 * n_pairs

        # トロッターステップあたりのゲート数
        # H_0 は最適化の対象外（単純なVirtRzゲート8個程度）
        gates_h0 = 8

        gates_per_step_v1 = gates_h0 + gates_transfer_v1 + gates_tta_v1
        gates_per_step_v2 = gates_h0 + gates_transfer_v2 + gates_tta_v2

        # シミュレーション全体のゲート数
        # 鈴木トロッター: 2回の適用
        gates_total_v1 = 2 * gates_per_step_v1 * n_steps
        gates_total_v2 = 2 * gates_per_step_v2 * n_steps

        reduction_rate = (gates_total_v1 - gates_total_v2) / gates_total_v1 if gates_total_v1 > 0 else 0.0

        return {
            "n_steps": n_steps,
            "n_pairs": n_pairs,
            "gates_per_pair": {
                "H_transfer_v1": gates_per_pair_transfer_v1,
                "H_transfer_v2": gates_per_pair_transfer_v2,
                "H_TTA_v1": gates_per_pair_tta_v1,
                "H_TTA_v2": gates_per_pair_tta_v2,
            },
            "gates_per_step": {
                "v1": gates_per_step_v1,
                "v2": gates_per_step_v2,
            },
            "gates_total": {
                "v1": gates_total_v1,
                "v2": gates_total_v2,
            },
            "reduction_rate": reduction_rate,
        }

    def print_comprehensive_result(self, result: ComprehensiveTestResult) -> None:
        """包括的テスト結果を表示."""
        # H_transfer結果

        # H_TTA結果

        # 詳細な結果表（時間刻み幅ごと）
        for i, _dt in enumerate(result.time_steps):
            result.h_transfer_results[i].fidelity
            result.h_tta_results[i].fidelity

    def print_simulation_estimate(self, estimate: dict) -> None:
        """シミュレーション推定結果を表示."""
        estimate["gates_per_pair"]

        estimate["gates_per_step"]

        estimate["gates_total"]


def run_comprehensive_tests():
    """包括的テストを実行."""
    tester = ComprehensiveMolecularTester()

    # テスト1: デフォルトパラメータ、複数の時間刻み幅

    params_default = PhysicalParameters()
    time_steps_1 = [0.1, 0.5, 1.0, 2.0]

    start_time = time.time()
    result_1 = tester.run_time_step_analysis(
        params=params_default, time_steps=time_steps_1, test_name="Default Parameters"
    )
    time.time() - start_time

    tester.print_comprehensive_result(result_1)

    # 4分子シミュレーション推定
    estimate_1 = tester.estimate_4_molecule_simulation_gates(result_1, n_steps=100)
    tester.print_simulation_estimate(estimate_1)

    # テスト2: 強い相互作用 (V=0.20, J=0.10)

    params_strong = PhysicalParameters(V=np.array([0.20, 0.20, 0.20]), J=np.array([0.10, 0.10, 0.10]))
    time_steps_2 = [0.5, 1.0]

    start_time = time.time()
    result_2 = tester.run_time_step_analysis(
        params=params_strong, time_steps=time_steps_2, test_name="Strong Interaction"
    )
    time.time() - start_time

    tester.print_comprehensive_result(result_2)

    estimate_2 = tester.estimate_4_molecule_simulation_gates(result_2, n_steps=100)
    tester.print_simulation_estimate(estimate_2)

    # テスト3: 弱い相互作用 (V=0.05, J=0.025)

    params_weak = PhysicalParameters(V=np.array([0.05, 0.05, 0.05]), J=np.array([0.025, 0.025, 0.025]))
    time_steps_3 = [0.5, 1.0]

    start_time = time.time()
    result_3 = tester.run_time_step_analysis(params=params_weak, time_steps=time_steps_3, test_name="Weak Interaction")
    time.time() - start_time

    tester.print_comprehensive_result(result_3)

    estimate_3 = tester.estimate_4_molecule_simulation_gates(result_3, n_steps=100)
    tester.print_simulation_estimate(estimate_3)

    # 総合サマリー

    all_results = [result_1, result_2, result_3]

    # すべてのテストで成功したか
    all_success = all(
        result.get_success_rate("H_transfer") == 1.0 and result.get_success_rate("H_TTA") == 1.0
        for result in all_results
    )

    # 最小忠実度
    min(result.get_min_fidelity("H_transfer") for result in all_results)
    min(result.get_min_fidelity("H_TTA") for result in all_results)

    # 平均削減率
    np.mean([result.get_average_reduction_rate("H_transfer") for result in all_results])
    np.mean([result.get_average_reduction_rate("H_TTA") for result in all_results])

    # 4分子シミュレーション比較
    for _i, (_result, estimate) in enumerate(
        [(result_1, estimate_1), (result_2, estimate_2), (result_3, estimate_3)], 1
    ):
        estimate["gates_total"]

    if all_success:
        pass
    else:
        pass

    return all_success


if __name__ == "__main__":
    success = run_comprehensive_tests()
    sys.exit(0 if success else 1)
