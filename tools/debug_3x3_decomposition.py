#!/usr/bin/env python3
"""
3×3ユニタリ分解のデバッグツール

Givens分解の各ステップを詳細に検証します。
"""

import numpy as np
from scipy.stats import unitary_group


def compute_givens_standard(a, b):
    """
    標準的なGivens回転パラメータ（QR分解の教科書的方法）
    
    目標: G† [a, b]^T = [r, 0]^T
    
    G = [[c,  s],
         [-s*, c*]]
    
    c = a*/r, s = b*/r
    """
    # 複素数に変換
    a = complex(a)
    b = complex(b)
    
    r = np.sqrt(abs(a)**2 + abs(b)**2)
    
    if r < 1e-15:
        c, s = 1.0, 0.0
    else:
        c = np.conj(a) / r
        s = np.conj(b) / r
    
    # 検証
    G = np.array([[c, s], [-s.conj(), c.conj()]], dtype=complex)
    v = np.array([a, b], dtype=complex)
    result = G.conj().T @ v
    
    print(f"  入力: a={a:.4f}, b={b:.4f}")
    print(f"  r={r:.4f}, c={c:.4f}, s={s:.4f}")
    print(f"  G†v = [{result[0]:.4f}, {result[1]:.4f}]")
    print(f"  期待: [{r:.4f}, 0.0000]")
    print(f"  誤差: |result[1]| = {abs(result[1]):.10f}")
    print(f"  ユニタリ性: ||G†G - I|| = {np.linalg.norm(G.conj().T @ G - np.eye(2)):.10f}")
    
    return G, r


def test_givens_basic():
    """基本的なGivens回転のテスト"""
    print("="*70)
    print("基本的なGivens回転のテスト")
    print("="*70)
    
    # テストケース1: 実数
    print("\nテスト1: 実数の場合")
    a, b = 3.0, 4.0
    G, r = compute_givens_standard(a, b)
    
    # テストケース2: 複素数
    print("\nテスト2: 複素数の場合")
    a, b = 3.0 + 4.0j, 1.0 + 2.0j
    G, r = compute_givens_standard(a, b)
    
    # テストケース3: 位相が大きい
    print("\nテスト3: 位相が大きい場合")
    a, b = np.exp(1j * 2.5), np.exp(1j * 1.3)
    G, r = compute_givens_standard(a, b)


def test_3x3_decomposition_step_by_step():
    """3×3分解をステップごとに検証"""
    print("\n" + "="*70)
    print("3×3ユニタリ分解のステップバイステップ検証")
    print("="*70)
    
    # ランダムな3×3ユニタリ
    U_original = unitary_group.rvs(3)
    print(f"\n元のユニタリ行列 U:")
    print(U_original)
    print(f"ユニタリ性チェック: ||U†U - I|| = {np.linalg.norm(U_original.conj().T @ U_original - np.eye(3)):.10f}")
    
    U_work = U_original.copy()
    rotations = []
    
    # ステップ1: U[1,0]をゼロにする
    print("\n" + "-"*70)
    print("ステップ1: U[1,0]をゼロにする（準位0と1の回転）")
    print("-"*70)
    a, b = U_work[0, 0], U_work[1, 0]
    print(f"U[0,0] = {a:.4f}, U[1,0] = {b:.4f}")
    
    G1_2x2, r1 = compute_givens_standard(a, b)
    
    # 3×3に拡張
    G1 = np.eye(3, dtype=complex)
    G1[0:2, 0:2] = G1_2x2
    
    print(f"\nG1 (3×3):")
    print(G1)
    print(f"G1のユニタリ性: ||G1†G1 - I|| = {np.linalg.norm(G1.conj().T @ G1 - np.eye(3)):.10f}")
    
    U_work = G1.conj().T @ U_work
    print(f"\nU after G1†:")
    print(U_work)
    print(f"U[1,0] = {U_work[1,0]:.10f} (期待: ≈0)")
    
    rotations.append(('G1', G1))
    
    # ステップ2: U[2,0]をゼロにする
    print("\n" + "-"*70)
    print("ステップ2: U[2,0]をゼロにする（準位0と2の回転）")
    print("-"*70)
    a, b = U_work[0, 0], U_work[2, 0]
    print(f"U[0,0] = {a:.4f}, U[2,0] = {b:.4f}")
    
    G2_2x2, r2 = compute_givens_standard(a, b)
    
    # 3×3に拡張（準位0と2）
    G2 = np.eye(3, dtype=complex)
    G2[0, 0] = G2_2x2[0, 0]
    G2[0, 2] = G2_2x2[0, 1]
    G2[2, 0] = G2_2x2[1, 0]
    G2[2, 2] = G2_2x2[1, 1]
    
    print(f"\nG2 (3×3):")
    print(G2)
    print(f"G2のユニタリ性: ||G2†G2 - I|| = {np.linalg.norm(G2.conj().T @ G2 - np.eye(3)):.10f}")
    
    U_work = G2.conj().T @ U_work
    print(f"\nU after G2†:")
    print(U_work)
    print(f"U[2,0] = {U_work[2,0]:.10f} (期待: ≈0)")
    
    rotations.append(('G2', G2))
    
    # ステップ3: U[2,1]をゼロにする
    print("\n" + "-"*70)
    print("ステップ3: U[2,1]をゼロにする（準位1と2の回転）")
    print("-"*70)
    a, b = U_work[1, 1], U_work[2, 1]
    print(f"U[1,1] = {a:.4f}, U[2,1] = {b:.4f}")
    
    G3_2x2, r3 = compute_givens_standard(a, b)
    
    # 3×3に拡張（準位1と2）
    G3 = np.eye(3, dtype=complex)
    G3[1:3, 1:3] = G3_2x2
    
    print(f"\nG3 (3×3):")
    print(G3)
    print(f"G3のユニタリ性: ||G3†G3 - I|| = {np.linalg.norm(G3.conj().T @ G3 - np.eye(3)):.10f}")
    
    U_work = G3.conj().T @ U_work
    print(f"\nU after G3†:")
    print(U_work)
    print(f"U[2,1] = {U_work[2,1]:.10f} (期待: ≈0)")
    print(f"上三角化誤差: {np.linalg.norm(np.tril(U_work, -1)):.10f}")
    
    rotations.append(('G3', G3))
    
    # 再構築
    print("\n" + "="*70)
    print("再構築の検証")
    print("="*70)
    
    # U_work は R（上三角）
    R = U_work.copy()
    
    # U_original = G1 @ G2 @ G3 @ R
    U_reconstructed = G1 @ G2 @ G3 @ R
    
    print(f"\n再構築されたユニタリ U_reconstructed:")
    print(U_reconstructed)
    
    print(f"\n元のユニタリ U_original:")
    print(U_original)
    
    print(f"\n差分 ||U_reconstructed - U_original||:")
    diff = np.linalg.norm(U_reconstructed - U_original)
    print(f"  {diff:.10f}")
    
    # 忠実度
    fidelity = abs(np.trace(U_original.conj().T @ U_reconstructed)) / 3
    print(f"\n忠実度: {fidelity:.10f}")
    
    return fidelity


if __name__ == '__main__':
    # 基本テスト
    test_givens_basic()
    
    # 3×3分解の詳細テスト
    fidelity = test_3x3_decomposition_step_by_step()
    
    print("\n" + "="*70)
    print(f"最終結果: 忠実度 = {fidelity:.10f}")
    if fidelity > 0.9999:
        print("✓ 合格")
    else:
        print("✗ 不合格 - さらなる調査が必要")
    print("="*70)
