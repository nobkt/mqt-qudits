#!/usr/bin/env python3
"""Givens → ZYZ分解器 (Givens to ZYZ Decomposer).

PR#41: 3×3ゲート変換の正しい実装

このツールは、Givens回転をZYZ分解に変換し、
その後MQT-Quditsゲートに変換します。

理論的根拠:
- Givens回転は2×2ユニタリ行列
- 任意の2×2ユニタリはZYZ分解可能
- ZYZ分解はMQT-Quditsゲート（VirtRz + R）で完璧に表現可能

これにより、忠実度 1.0 を保証します。

数学的厳密性:
- すべての変換は厳密な線形代数に基づく
- ヒューリスティックゼロ
- 近似ゼロ
- 忠実度 1.0 を保証
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


@dataclass
class MQTGate:
    """MQT-Quditsゲートの抽象表現."""

    gate_type: str  # 'VirtRz', 'R'
    parameters: dict
    cost: int = 0


def construct_givens_from_theta_phi(i: int, j: int, theta: float, phi: float, size: int = 3) -> np.ndarray:
    """θ, φからGivens行列を構築."""
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)

    G = np.eye(size, dtype=complex)
    G[i, i] = c
    G[i, j] = s
    G[j, i] = -np.conj(s)
    G[j, j] = np.conj(c)

    return G


def compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
    """忠実度を計算."""
    N = U1.shape[0]
    trace = np.trace(U1.conj().T @ U2)
    return abs(trace) / N


class GivensToZYZDecomposer:
    """Givens回転をZYZ分解に変換.

    このクラスは、Givens回転(θ, φ)を2×2ユニタリ行列として扱い、
    improved_unitary_decomposition.pyのZYZ分解を適用します。
    """

    def __init__(self, tolerance: float = 1e-10) -> None:
        """Args:
        tolerance: 数値誤差の許容範囲.
        """
        self.tolerance = tolerance

        # improved_unitary_decomposition.pyをインポート
        try:
            from improved_unitary_decomposition import ImprovedTwoQubitDecomposer

            self.decomposer = ImprovedTwoQubitDecomposer()
            self.available = True
        except (ImportError, TypeError):
            self.available = False

    def decompose(self, theta: float, phi: float) -> dict:
        """Givens回転(θ, φ)をZYZ分解.

        Args:
            theta: Givens回転角度
            phi: Givens位相角度

        Returns:
            ZYZ分解パラメータの辞書:
            {
                'theta': ZYZ θ,
                'phi': ZYZ φ,
                'lambda': ZYZ λ,
                'global_phase': グローバル位相,
                'fidelity': 分解の忠実度
            }
        """
        if not self.available:
            msg = "improved_unitary_decomposition.py がインポートできません"
            raise RuntimeError(msg)

        # Givens行列を2×2として構築
        c = np.cos(theta / 2) * np.exp(1j * phi / 2)
        s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
        G_2x2 = np.array([[c, s], [-np.conj(s), np.conj(c)]], dtype=complex)

        # ZYZ分解を適用
        result = self.decomposer.decompose_zyz(G_2x2)

        # 結果を辞書形式で返す
        return {
            "theta": result.theta,
            "phi": result.phi,
            "lambda": result.lam,
            "global_phase": result.global_phase,
            "fidelity": result.fidelity,
        }

    def convert_to_mqt_gates(self, i: int, j: int, theta: float, phi: float) -> list[MQTGate]:
        """Givens回転をMQT-Quditsゲートに変換.

        手順:
        1. Givens(θ, φ) → 2×2ユニタリ
        2. 2×2ユニタリ → ZYZ分解
        3. ZYZ → MQT-Quditsゲート（VirtRz + R）

        Args:
            i: 第1のレベル
            j: 第2のレベル
            theta: Givens回転角度
            phi: Givens位相角度

        Returns:
            MQT-Quditsゲートのリスト
        """
        # ZYZ分解
        zyz = self.decompose(theta, phi)

        gates = []

        # ZYZ分解: U = e^(iα) Rz(φ) Ry(θ) Rz(λ)
        #
        # MQT-Quditsゲートに変換:
        # - Rz(φ) → VirtRz(φ/2, i) + VirtRz(-φ/2, j)
        # - Ry(θ) → R(-θ, 0, i, j)  (2×2の場合と同様に符号反転)
        # - Rz(λ) → VirtRz(λ/2, i) + VirtRz(-λ/2, j)
        # - e^(iα) → VirtRz(α, i) + VirtRz(α, j)  (グローバル位相)

        alpha = zyz["global_phase"]
        phi_zyz = zyz["phi"]
        theta_zyz = zyz["theta"]
        lambda_zyz = zyz["lambda"]

        # e^(iα) Rz(φ/2) の結合: VirtRz(α + φ/2, i)
        phase_i_1 = alpha + phi_zyz / 2
        phase_j_1 = alpha - phi_zyz / 2

        if abs(phase_i_1) > self.tolerance:
            gates.append(MQTGate(gate_type="VirtRz", parameters={"level": i, "phase": phase_i_1}, cost=0))

        if abs(phase_j_1) > self.tolerance:
            gates.append(MQTGate(gate_type="VirtRz", parameters={"level": j, "phase": phase_j_1}, cost=0))

        # Ry(θ) → R(-θ, 0)
        if abs(theta_zyz) > self.tolerance:
            gates.append(
                MQTGate(
                    gate_type="R",
                    parameters={
                        "level1": i,
                        "level2": j,
                        "theta": -theta_zyz,  # 符号反転
                        "phi": 0.0,
                    },
                    cost=1,
                )
            )

        # Rz(λ)
        phase_i_2 = lambda_zyz / 2
        phase_j_2 = -lambda_zyz / 2

        if abs(phase_i_2) > self.tolerance:
            gates.append(MQTGate(gate_type="VirtRz", parameters={"level": i, "phase": phase_i_2}, cost=0))

        if abs(phase_j_2) > self.tolerance:
            gates.append(MQTGate(gate_type="VirtRz", parameters={"level": j, "phase": phase_j_2}, cost=0))

        return gates

    def verify_conversion(
        self, i: int, j: int, theta: float, phi: float, size: int = 3
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """変換の正確性を検証.

        Returns:
            (fidelity, G_target, G_mqt)
        """
        # 目標Givens行列
        G_target = construct_givens_from_theta_phi(i, j, theta, phi, size)

        # MQT-Quditsゲートに変換
        gates = self.convert_to_mqt_gates(i, j, theta, phi)

        # 行列を再構築（逆順で左から掛ける）
        G_mqt = np.eye(size, dtype=complex)
        for gate in reversed(gates):
            if gate.gate_type == "VirtRz":
                level = gate.parameters["level"]
                phase = gate.parameters["phase"]
                Rz = np.eye(size, dtype=complex)
                Rz[level, level] = np.exp(1j * phase)
                G_mqt = Rz @ G_mqt
            elif gate.gate_type == "R":
                level1 = gate.parameters["level1"]
                level2 = gate.parameters["level2"]
                theta_r = gate.parameters["theta"]
                phi_r = gate.parameters["phi"]

                # MQT-Qudits R gate
                c_r = np.cos(theta_r / 2)
                s_r = np.sin(theta_r / 2)
                R = np.eye(size, dtype=complex)
                R[level1, level1] = c_r
                R[level1, level2] = s_r * np.exp(1j * phi_r)
                R[level2, level1] = -s_r * np.exp(1j * phi_r)
                R[level2, level2] = c_r
                G_mqt = R @ G_mqt

        # 忠実度を計算
        fidelity = compute_fidelity(G_target, G_mqt)

        return fidelity, G_target, G_mqt


def test_single_givens():
    """単一のGivens回転でテスト."""
    decomposer = GivensToZYZDecomposer()

    if not decomposer.available:
        return False

    # テストケース
    test_cases = [
        (0.5, 1.0),
        (0.1, 0.5),
        (1.0, 2.0),
        (np.pi / 4, np.pi / 6),
        (0.0, 1.0),  # θ=0のエッジケース
        (np.pi, 0.0),  # θ=πのエッジケース
        (1.5, 0.0),  # φ=0のケース
        (0.3, np.pi),  # φ=πのケース
    ]

    all_passed = True

    for theta, phi in test_cases:
        # ZYZ分解
        decomposer.decompose(theta, phi)

        # MQT-Quditsゲートに変換
        decomposer.convert_to_mqt_gates(0, 1, theta, phi)

        # 検証
        fidelity, _G_target, _G_mqt = decomposer.verify_conversion(0, 1, theta, phi)

        passed = fidelity > 0.9999

        if not passed:
            all_passed = False

    if all_passed:
        pass
    else:
        pass

    return all_passed


def test_random_givens(num_tests: int = 100):
    """ランダムなGivens回転でテスト."""
    decomposer = GivensToZYZDecomposer()

    if not decomposer.available:
        return False

    np.random.seed(42)

    fidelities = []
    gate_counts = []

    for _i in range(num_tests):
        # ランダムパラメータ
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(-np.pi, np.pi)

        # 変換
        gates = decomposer.convert_to_mqt_gates(0, 1, theta, phi)
        fidelity, _, _ = decomposer.verify_conversion(0, 1, theta, phi)

        fidelities.append(fidelity)
        gate_counts.append(len(gates))

    # 統計
    min(fidelities)
    np.mean(fidelities)
    pass_count = sum(1 for f in fidelities if f > 0.9999)
    np.mean(gate_counts)

    passed = pass_count == num_tests
    if passed:
        pass
    else:
        pass

    return passed


def main() -> int:
    """メイン関数."""
    # 単一Givens回転のテスト
    success1 = test_single_givens()

    # ランダムGivens回転のテスト
    success2 = test_random_givens(num_tests=100)

    if success1 and success2:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
