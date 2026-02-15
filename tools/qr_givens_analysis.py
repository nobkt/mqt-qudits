#!/usr/bin/env python3
"""QR分解とGivensパラメータの関係分析.

PR#41: integrated_sparse_compiler.pyが抽出するGivensパラメータと
理論的なパラメータ化の関係を明らかにする

重要な問題:
- integrated_sparse_compiler.pyはQR分解からGivensパラメータを抽出
- しかし、抽出されたパラメータがgivens_rotation_theory_ja.mdの
  パラメータ化とどう対応するか不明確

このツールは、QR分解から正しいGivensパラメータを抽出する方法を検証する
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


def compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
    """忠実度を計算."""
    N = U1.shape[0]
    trace = np.trace(U1.conj().T @ U2)
    return abs(trace) / N


def construct_givens_from_cs(i: int, j: int, c: complex, s: complex, size: int = 3) -> np.ndarray:
    """c, sパラメータからGivens行列を構築.

    G[i,i] = c, G[i,j] = s
    G[j,i] = -s*, G[j,j] = c*
    """
    G = np.eye(size, dtype=complex)
    G[i, i] = c
    G[i, j] = s
    G[j, i] = -np.conj(s)
    G[j, j] = np.conj(c)
    return G


def extract_theta_phi_from_cs(c: complex, s: complex):
    """c, sからθ, φを抽出.

    理論:
    c = cos(θ/2) e^(iφ/2)
    s = sin(θ/2) e^(-iφ/2)

    したがって:
    |c| = cos(θ/2) => θ = 2 arccos(|c|)
    arg(c) = φ/2 => φ = 2 arg(c)

    検証:
    |s| = sin(θ/2) = sin(arccos(|c|)) = sqrt(1 - |c|²)
    arg(s) = -φ/2 = -arg(c)
    """
    # θの抽出
    theta = 2 * np.arccos(np.clip(abs(c), 0, 1))

    # φの抽出
    phi = 2 * np.angle(c)

    return theta, phi


def test_givens_extraction() -> None:
    """Givensパラメータ抽出のテスト."""
    # テストケース1: 単純なGivens回転
    theta_input = 0.5
    phi_input = 1.0

    c_input = np.cos(theta_input / 2) * np.exp(1j * phi_input / 2)
    s_input = np.sin(theta_input / 2) * np.exp(-1j * phi_input / 2)

    # Givens行列を構築
    G_input = construct_givens_from_cs(0, 1, c_input, s_input, size=3)

    # c, sから逆算
    theta_extract, phi_extract = extract_theta_phi_from_cs(c_input, s_input)

    # 再構築して確認
    c_reconstruct = np.cos(theta_extract / 2) * np.exp(1j * phi_extract / 2)
    s_reconstruct = np.sin(theta_extract / 2) * np.exp(-1j * phi_extract / 2)
    G_reconstruct = construct_givens_from_cs(0, 1, c_reconstruct, s_reconstruct, size=3)

    compute_fidelity(G_input, G_reconstruct)

    # テストケース2: QR分解から抽出

    # ランダムユニタリ
    np.random.seed(42)
    U = np.array([
        [0.8 + 0.2j, 0.5 - 0.1j, 0.1 + 0.2j],
        [0.3 + 0.4j, -0.6 + 0.3j, 0.4 - 0.2j],
        [-0.2 + 0.3j, 0.3 + 0.4j, 0.7 + 0.1j],
    ])
    # 正規化してユニタリにする
    U, _ = np.linalg.qr(U)

    # QR分解
    Q, _R = np.linalg.qr(U)

    # integrated_sparse_compiler.pyの方法を模倣
    # Q = G(0,1) @ G(0,2) @ G(1,2)
    # Q†を右から掛けていく

    Q_work = Q.copy()
    givens_list = []

    # Step 1: G(1,2)を抽出してQ[2,1]をゼロ化
    a = Q_work[1, 1]
    b = Q_work[2, 1]
    r = np.sqrt(abs(a) ** 2 + abs(b) ** 2)
    if r > 1e-10:
        c1 = np.conj(a) / r
        s1 = np.conj(b) / r

        theta1, phi1 = extract_theta_phi_from_cs(c1, s1)
        givens_list.append((1, 2, theta1, phi1))

        # G† を適用
        G1 = construct_givens_from_cs(1, 2, c1, s1, size=3)
        Q_work = G1.conj().T @ Q_work

    # Step 2: G(0,2)を抽出してQ[2,0]をゼロ化
    a = Q_work[0, 0]
    b = Q_work[2, 0]
    r = np.sqrt(abs(a) ** 2 + abs(b) ** 2)
    if r > 1e-10:
        c2 = np.conj(a) / r
        s2 = np.conj(b) / r

        theta2, phi2 = extract_theta_phi_from_cs(c2, s2)
        givens_list.append((0, 2, theta2, phi2))

        # G† を適用
        G2 = construct_givens_from_cs(0, 2, c2, s2, size=3)
        Q_work = G2.conj().T @ Q_work

    # Step 3: G(0,1)を抽出してQ[1,0]をゼロ化
    a = Q_work[0, 0]
    b = Q_work[1, 0]
    r = np.sqrt(abs(a) ** 2 + abs(b) ** 2)
    if r > 1e-10:
        c3 = np.conj(a) / r
        s3 = np.conj(b) / r

        theta3, phi3 = extract_theta_phi_from_cs(c3, s3)
        givens_list.append((0, 1, theta3, phi3))

        # G† を適用
        G3 = construct_givens_from_cs(0, 1, c3, s3, size=3)
        Q_work = G3.conj().T @ Q_work

    # 再構築

    Q_reconstruct = np.eye(3, dtype=complex)
    # 逆順で適用: G(0,1) @ G(0,2) @ G(1,2)
    for i, j, theta, phi in reversed(givens_list):
        c = np.cos(theta / 2) * np.exp(1j * phi / 2)
        s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
        G = construct_givens_from_cs(i, j, c, s, size=3)
        Q_reconstruct = G @ Q_reconstruct

    compute_fidelity(Q, Q_reconstruct)

    # MQT-Quditsゲートへの変換をテスト

    # 単一のGivens回転でテスト
    i, j, theta, phi = 0, 1, 0.5, 1.0

    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    G_target = construct_givens_from_cs(i, j, c, s, size=3)

    # MQT-Qudits Rゲートの定義を確認
    def construct_mqt_r(level1, level2, theta_r, phi_r, size=3):
        R = np.eye(size, dtype=complex)
        c_r = np.cos(theta_r / 2)
        s_r = np.sin(theta_r / 2)
        R[level1, level1] = c_r
        R[level1, level2] = s_r * np.exp(1j * phi_r)
        R[level2, level1] = -s_r * np.exp(1j * phi_r)
        R[level2, level2] = c_r
        return R

    # R(θ, -φ)を試す
    R_mqt = construct_mqt_r(i, j, theta, -phi, size=3)
    compute_fidelity(G_target, R_mqt)

    # 実際、MQT R gate の定義を考慮すると:
    # R(θ, φ): R[0,1] = sin(θ/2) e^(iφ)
    # Givens: G[0,1] = sin(θ/2) e^(-iφ/2)
    #
    # マッチさせるには: e^(iφ_R) = e^(-iφ_G/2)
    # => φ_R = -φ_G/2

    R_mqt2 = construct_mqt_r(i, j, theta, -phi / 2, size=3)
    compute_fidelity(G_target, R_mqt2)

    # さらに、対角要素も考慮する必要がある
    # Givens: G[0,0] = cos(θ/2) e^(iφ/2)
    # MQT R: R[0,0] = cos(θ/2)
    #
    # 位相差 e^(iφ/2) を補正するには、VirtRzが必要

    def construct_virtrz(level, phase, size=3):
        Rz = np.eye(size, dtype=complex)
        Rz[level, level] = np.exp(1j * phase)
        return Rz

    gates = [construct_virtrz(i, phi / 2), construct_virtrz(j, phi / 2), construct_mqt_r(i, j, theta, -phi)]
    G_mqt = np.eye(3, dtype=complex)
    for gate in reversed(gates):
        G_mqt = gate @ G_mqt
    compute_fidelity(G_target, G_mqt)


def main() -> int:
    """メイン関数."""
    test_givens_extraction()
    return 0


if __name__ == "__main__":
    sys.exit(main())
