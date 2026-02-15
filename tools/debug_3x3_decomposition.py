#!/usr/bin/env python3
"""3×3ユニタリ分解のデバッグツール.

Givens分解の各ステップを詳細に検証します。
"""

from __future__ import annotations

import numpy as np
from scipy.stats import unitary_group


def compute_givens_standard(a, b):
    """標準的なGivens回転パラメータ（QR分解の教科書的方法）.

    目標: G† [a, b]^T = [r, 0]^T

    G = [[c,  s],
         [-s*, c*]]

    c = a*/r, s = b*/r
    """
    # 複素数に変換
    a = complex(a)
    b = complex(b)

    r = np.sqrt(abs(a) ** 2 + abs(b) ** 2)

    if r < 1e-15:
        c, s = 1.0, 0.0
    else:
        c = np.conj(a) / r
        s = np.conj(b) / r

    # 検証
    G = np.array([[c, s], [-s.conj(), c.conj()]], dtype=complex)
    v = np.array([a, b], dtype=complex)
    G.conj().T @ v

    return G, r


def test_givens_basic() -> None:
    """基本的なGivens回転のテスト."""
    # テストケース1: 実数
    a, b = 3.0, 4.0
    _G, _r = compute_givens_standard(a, b)

    # テストケース2: 複素数
    a, b = 3.0 + 4.0j, 1.0 + 2.0j
    _G, _r = compute_givens_standard(a, b)

    # テストケース3: 位相が大きい
    a, b = np.exp(1j * 2.5), np.exp(1j * 1.3)
    _G, _r = compute_givens_standard(a, b)


def test_3x3_decomposition_step_by_step():
    """3×3分解をステップごとに検証."""
    # ランダムな3×3ユニタリ
    U_original = unitary_group.rvs(3)

    U_work = U_original.copy()
    rotations = []

    # ステップ1: U[1,0]をゼロにする
    a, b = U_work[0, 0], U_work[1, 0]

    G1_2x2, _r1 = compute_givens_standard(a, b)

    # 3×3に拡張
    G1 = np.eye(3, dtype=complex)
    G1[0:2, 0:2] = G1_2x2

    U_work = G1.conj().T @ U_work

    rotations.append(("G1", G1))

    # ステップ2: U[2,0]をゼロにする
    a, b = U_work[0, 0], U_work[2, 0]

    G2_2x2, _r2 = compute_givens_standard(a, b)

    # 3×3に拡張（準位0と2）
    G2 = np.eye(3, dtype=complex)
    G2[0, 0] = G2_2x2[0, 0]
    G2[0, 2] = G2_2x2[0, 1]
    G2[2, 0] = G2_2x2[1, 0]
    G2[2, 2] = G2_2x2[1, 1]

    U_work = G2.conj().T @ U_work

    rotations.append(("G2", G2))

    # ステップ3: U[2,1]をゼロにする
    a, b = U_work[1, 1], U_work[2, 1]

    G3_2x2, _r3 = compute_givens_standard(a, b)

    # 3×3に拡張（準位1と2）
    G3 = np.eye(3, dtype=complex)
    G3[1:3, 1:3] = G3_2x2

    U_work = G3.conj().T @ U_work

    rotations.append(("G3", G3))

    # 再構築

    # U_work は R（上三角）
    R = U_work.copy()

    # U_original = G1 @ G2 @ G3 @ R
    U_reconstructed = G1 @ G2 @ G3 @ R

    np.linalg.norm(U_reconstructed - U_original)

    # 忠実度
    return abs(np.trace(U_original.conj().T @ U_reconstructed)) / 3


if __name__ == "__main__":
    # 基本テスト
    test_givens_basic()

    # 3×3分解の詳細テスト
    fidelity = test_3x3_decomposition_step_by_step()

    if fidelity > 0.9999:
        pass
    else:
        pass
