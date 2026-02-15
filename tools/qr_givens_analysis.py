#!/usr/bin/env python3
"""
QR分解とGivensパラメータの関係分析

PR#41: integrated_sparse_compiler.pyが抽出するGivensパラメータと
理論的なパラメータ化の関係を明らかにする

重要な問題:
- integrated_sparse_compiler.pyはQR分解からGivensパラメータを抽出
- しかし、抽出されたパラメータがgivens_rotation_theory_ja.mdの
  パラメータ化とどう対応するか不明確

このツールは、QR分解から正しいGivensパラメータを抽出する方法を検証する
"""

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


def compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
    """忠実度を計算"""
    N = U1.shape[0]
    trace = np.trace(U1.conj().T @ U2)
    return abs(trace) / N


def construct_givens_from_cs(i: int, j: int, c: complex, s: complex, size: int = 3) -> np.ndarray:
    """
    c, sパラメータからGivens行列を構築
    
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
    """
    c, sからθ, φを抽出
    
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


def test_givens_extraction():
    """
    Givensパラメータ抽出のテスト
    """
    print("="*70)
    print("QR分解とGivensパラメータの関係")
    print("="*70)
    
    # テストケース1: 単純なGivens回転
    print("\nテストケース1: 既知のGivens回転")
    theta_input = 0.5
    phi_input = 1.0
    
    c_input = np.cos(theta_input / 2) * np.exp(1j * phi_input / 2)
    s_input = np.sin(theta_input / 2) * np.exp(-1j * phi_input / 2)
    
    print(f"入力パラメータ: θ={theta_input:.6f}, φ={phi_input:.6f}")
    print(f"c = {c_input:.6f}, s = {s_input:.6f}")
    
    # Givens行列を構築
    G_input = construct_givens_from_cs(0, 1, c_input, s_input, size=3)
    print(f"\nGivens行列:")
    print(G_input)
    
    # c, sから逆算
    theta_extract, phi_extract = extract_theta_phi_from_cs(c_input, s_input)
    print(f"\n抽出されたパラメータ: θ={theta_extract:.6f}, φ={phi_extract:.6f}")
    print(f"誤差: Δθ={abs(theta_extract - theta_input):.2e}, Δφ={abs(phi_extract - phi_input):.2e}")
    
    # 再構築して確認
    c_reconstruct = np.cos(theta_extract / 2) * np.exp(1j * phi_extract / 2)
    s_reconstruct = np.sin(theta_extract / 2) * np.exp(-1j * phi_extract / 2)
    G_reconstruct = construct_givens_from_cs(0, 1, c_reconstruct, s_reconstruct, size=3)
    
    fid = compute_fidelity(G_input, G_reconstruct)
    print(f"再構築忠実度: {fid:.10f}")
    
    # テストケース2: QR分解から抽出
    print("\n" + "="*70)
    print("テストケース2: QR分解からの抽出")
    print("="*70)
    
    # ランダムユニタリ
    np.random.seed(42)
    U = np.array([
        [0.8+0.2j, 0.5-0.1j, 0.1+0.2j],
        [0.3+0.4j, -0.6+0.3j, 0.4-0.2j],
        [-0.2+0.3j, 0.3+0.4j, 0.7+0.1j]
    ])
    # 正規化してユニタリにする
    U, _ = np.linalg.qr(U)
    
    print(f"元のユニタリ行列:")
    print(U)
    
    # QR分解
    Q, R = np.linalg.qr(U)
    print(f"\nQR分解:")
    print(f"Q:")
    print(Q)
    print(f"\nR:")
    print(R)
    
    # integrated_sparse_compiler.pyの方法を模倣
    # Q = G(0,1) @ G(0,2) @ G(1,2)
    # Q†を右から掛けていく
    
    print("\n" + "="*70)
    print("Givens回転の抽出（QR分解の逆過程）")
    print("="*70)
    
    Q_work = Q.copy()
    givens_list = []
    
    # Step 1: G(1,2)を抽出してQ[2,1]をゼロ化
    print("\nステップ1: G(1,2)を抽出")
    a = Q_work[1, 1]
    b = Q_work[2, 1]
    r = np.sqrt(abs(a)**2 + abs(b)**2)
    if r > 1e-10:
        c1 = np.conj(a) / r
        s1 = np.conj(b) / r
        print(f"  a={a:.6f}, b={b:.6f}, r={r:.6f}")
        print(f"  c={c1:.6f}, s={s1:.6f}")
        
        theta1, phi1 = extract_theta_phi_from_cs(c1, s1)
        print(f"  θ={theta1:.6f}, φ={phi1:.6f}")
        givens_list.append((1, 2, theta1, phi1))
        
        # G† を適用
        G1 = construct_givens_from_cs(1, 2, c1, s1, size=3)
        Q_work = G1.conj().T @ Q_work
        print(f"  適用後のQ_work[2,1]={Q_work[2,1]:.6f} (ゼロに近いはず)")
    
    # Step 2: G(0,2)を抽出してQ[2,0]をゼロ化
    print("\nステップ2: G(0,2)を抽出")
    a = Q_work[0, 0]
    b = Q_work[2, 0]
    r = np.sqrt(abs(a)**2 + abs(b)**2)
    if r > 1e-10:
        c2 = np.conj(a) / r
        s2 = np.conj(b) / r
        print(f"  a={a:.6f}, b={b:.6f}, r={r:.6f}")
        print(f"  c={c2:.6f}, s={s2:.6f}")
        
        theta2, phi2 = extract_theta_phi_from_cs(c2, s2)
        print(f"  θ={theta2:.6f}, φ={phi2:.6f}")
        givens_list.append((0, 2, theta2, phi2))
        
        # G† を適用
        G2 = construct_givens_from_cs(0, 2, c2, s2, size=3)
        Q_work = G2.conj().T @ Q_work
        print(f"  適用後のQ_work[2,0]={Q_work[2,0]:.6f} (ゼロに近いはず)")
    
    # Step 3: G(0,1)を抽出してQ[1,0]をゼロ化
    print("\nステップ3: G(0,1)を抽出")
    a = Q_work[0, 0]
    b = Q_work[1, 0]
    r = np.sqrt(abs(a)**2 + abs(b)**2)
    if r > 1e-10:
        c3 = np.conj(a) / r
        s3 = np.conj(b) / r
        print(f"  a={a:.6f}, b={b:.6f}, r={r:.6f}")
        print(f"  c={c3:.6f}, s={s3:.6f}")
        
        theta3, phi3 = extract_theta_phi_from_cs(c3, s3)
        print(f"  θ={theta3:.6f}, φ={phi3:.6f}")
        givens_list.append((0, 1, theta3, phi3))
        
        # G† を適用
        G3 = construct_givens_from_cs(0, 1, c3, s3, size=3)
        Q_work = G3.conj().T @ Q_work
        print(f"  適用後のQ_work[1,0]={Q_work[1,0]:.6f} (ゼロに近いはず)")
    
    print(f"\n最終的なQ_work (対角行列に近いはず):")
    print(Q_work)
    
    # 再構築
    print("\n" + "="*70)
    print("Givens回転から再構築")
    print("="*70)
    
    Q_reconstruct = np.eye(3, dtype=complex)
    # 逆順で適用: G(0,1) @ G(0,2) @ G(1,2)
    for i, j, theta, phi in reversed(givens_list):
        c = np.cos(theta / 2) * np.exp(1j * phi / 2)
        s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
        G = construct_givens_from_cs(i, j, c, s, size=3)
        Q_reconstruct = G @ Q_reconstruct
        print(f"G({i},{j}; θ={theta:.6f}, φ={phi:.6f}) を適用")
    
    print(f"\n再構築されたQ:")
    print(Q_reconstruct)
    
    print(f"\n元のQ:")
    print(Q)
    
    fid = compute_fidelity(Q, Q_reconstruct)
    print(f"\n再構築忠実度: {fid:.10f}")
    
    # MQT-Quditsゲートへの変換をテスト
    print("\n" + "="*70)
    print("MQT-Quditsゲートへの変換テスト")
    print("="*70)
    
    # 単一のGivens回転でテスト
    i, j, theta, phi = 0, 1, 0.5, 1.0
    print(f"\n単一Givens: G({i},{j}; θ={theta:.6f}, φ={phi:.6f})")
    
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    G_target = construct_givens_from_cs(i, j, c, s, size=3)
    
    print(f"\n目標Givens行列:")
    print(G_target)
    
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
    print(f"\n試行: R(θ, -φ)")
    R_mqt = construct_mqt_r(i, j, theta, -phi, size=3)
    print(R_mqt)
    fid_r = compute_fidelity(G_target, R_mqt)
    print(f"忠実度: {fid_r:.10f}")
    
    # 実際、MQT R gate の定義を考慮すると:
    # R(θ, φ): R[0,1] = sin(θ/2) e^(iφ)
    # Givens: G[0,1] = sin(θ/2) e^(-iφ/2)
    # 
    # マッチさせるには: e^(iφ_R) = e^(-iφ_G/2)
    # => φ_R = -φ_G/2
    
    print(f"\n試行: R(θ, -φ/2)")
    R_mqt2 = construct_mqt_r(i, j, theta, -phi/2, size=3)
    print(R_mqt2)
    fid_r2 = compute_fidelity(G_target, R_mqt2)
    print(f"忠実度: {fid_r2:.10f}")
    
    # さらに、対角要素も考慮する必要がある
    # Givens: G[0,0] = cos(θ/2) e^(iφ/2)
    # MQT R: R[0,0] = cos(θ/2)
    #
    # 位相差 e^(iφ/2) を補正するには、VirtRzが必要
    
    print(f"\n試行: Rz(φ/2)_0 @ Rz(φ/2)_1 @ R(θ, -φ)")
    def construct_virtrz(level, phase, size=3):
        Rz = np.eye(size, dtype=complex)
        Rz[level, level] = np.exp(1j * phase)
        return Rz
    
    gates = [
        construct_virtrz(i, phi/2),
        construct_virtrz(j, phi/2),
        construct_mqt_r(i, j, theta, -phi)
    ]
    G_mqt = np.eye(3, dtype=complex)
    for gate in reversed(gates):
        G_mqt = gate @ G_mqt
    print(G_mqt)
    fid_mqt = compute_fidelity(G_target, G_mqt)
    print(f"忠実度: {fid_mqt:.10f}")


def main():
    """メイン関数"""
    test_givens_extraction()
    return 0


if __name__ == '__main__':
    sys.exit(main())
