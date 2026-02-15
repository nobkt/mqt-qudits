#!/usr/bin/env python3
"""
正しいGivens回転分解 (Correct Givens Rotation Decomposition)

PR#41: 3×3ゲート変換の正しい実装

このツールは、Givens回転の正しい分解方法を実装します。

重要な発見:
givens_rotation_theory_ja.mdのパラメータ化された定義:
  c = cos(θ/2) e^(iφ/2)
  s = sin(θ/2) e^(-iφ/2)

Givens行列:
  G[i,i] = c
  G[i,j] = s  
  G[j,i] = -s*
  G[j,j] = c*

この形式では、G ≠ Rz(φ/2)_i @ Ry(θ) @ Rz(-φ/2)_j

正しい分解:
  G(i,j; θ, φ) = Rz(α)_i @ Rz(β)_j @ R(θ, γ) @ Rz(δ)_i @ Rz(ε)_j

または、より直接的に:
  G行列を直接MQT-Quditsゲートに変換
"""

import sys
from pathlib import Path
import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent))


def construct_givens_theory(i: int, j: int, theta: float, phi: float, size: int = 3) -> np.ndarray:
    """
    givens_rotation_theory_ja.md の定義に基づくGivens行列
    
    c = cos(θ/2) e^(iφ/2)
    s = sin(θ/2) e^(-iφ/2)
    
    G[i,i] = c, G[i,j] = s
    G[j,i] = -s*, G[j,j] = c*
    """
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    
    G = np.eye(size, dtype=complex)
    G[i, i] = c
    G[i, j] = s
    G[j, i] = -np.conj(s)
    G[j, j] = np.conj(c)
    
    return G


def compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
    """忠実度を計算"""
    N = U1.shape[0]
    trace = np.trace(U1.conj().T @ U2)
    return abs(trace) / N


def analyze_givens_structure(theta: float, phi: float):
    """Givens回転の構造を分析"""
    
    print("="*70)
    print("Givens回転の構造分析")
    print("="*70)
    print(f"\nパラメータ: θ={theta:.6f}, φ={phi:.6f}")
    
    # パラメータ化されたGivens
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    
    print(f"\nパラメータ化:")
    print(f"  c = cos(θ/2) e^(iφ/2) = {c:.6f}")
    print(f"  s = sin(θ/2) e^(-iφ/2) = {s:.6f}")
    print(f"  -s* = {-np.conj(s):.6f}")
    print(f"  c* = {np.conj(c):.6f}")
    
    # Givens行列
    G = construct_givens_theory(0, 1, theta, phi, size=3)
    print(f"\nGivens行列 G(0,1; θ, φ):")
    print(G)
    
    # 位相と振幅を分離
    print(f"\n各要素の振幅と位相:")
    print(f"  G[0,0]: |{abs(G[0,0]):.6f}| × e^(i×{np.angle(G[0,0]):.6f})")
    print(f"  G[0,1]: |{abs(G[0,1]):.6f}| × e^(i×{np.angle(G[0,1]):.6f})")
    print(f"  G[1,0]: |{abs(G[1,0]):.6f}| × e^(i×{np.angle(G[1,0]):.6f})")
    print(f"  G[1,1]: |{abs(G[1,1]):.6f}| × e^(i×{np.angle(G[1,1]):.6f})")
    
    # ZYZ形式との比較
    print(f"\n観察:")
    print(f"  - G[0,0]とG[1,1]は共役複素数")
    print(f"  - G[0,1]とG[1,0]は -共役 の関係")
    print(f"  - 対角要素の位相: ±φ/2")
    print(f"  - 非対角要素の位相: -φ/2 (G[0,1]), -φ/2+π (G[1,0])")


def find_correct_mqt_decomposition(theta: float, phi: float):
    """
    正しいMQT-Quditsゲート分解を探す
    
    試行:
    1. 対角と非対角の位相を分離
    2. 全体の位相回転を考慮
    3. 複数のVirtRzとRゲートの組み合わせ
    """
    
    print("\n" + "="*70)
    print("MQT-Quditsゲート分解の探索")
    print("="*70)
    
    G_target = construct_givens_theory(0, 1, theta, phi, size=3)
    
    # 試行1: 位相を分離した分解
    print("\n試行1: 位相分離分解")
    print("  G = Rz(φ/2)_0 @ Rz(-φ/2)_1 @ R(θ, φ) @ Rz(0)_0 @ Rz(0)_1")
    
    def construct_r_gate(theta_r, phi_r):
        """MQT-Qudits R ゲート"""
        c_r = np.cos(theta_r / 2)
        s_r = np.sin(theta_r / 2)
        R = np.array([
            [c_r, s_r * np.exp(1j * phi_r), 0],
            [-s_r * np.exp(1j * phi_r), c_r, 0],
            [0, 0, 1]
        ], dtype=complex)
        return R
    
    def construct_virtrz(level, phase, size=3):
        """VirtRz ゲート"""
        Rz = np.eye(size, dtype=complex)
        Rz[level, level] = np.exp(1j * phase)
        return Rz
    
    # 試行1a: R(θ, φ)を直接使用
    gates_1a = [
        construct_virtrz(0, 0.0),
        construct_virtrz(1, 0.0),
        construct_r_gate(theta, phi),
        construct_virtrz(0, 0.0),
        construct_virtrz(1, 0.0)
    ]
    G_1a = np.eye(3, dtype=complex)
    for gate in reversed(gates_1a):
        G_1a = gate @ G_1a
    fid_1a = compute_fidelity(G_target, G_1a)
    print(f"  1a. R(θ, φ): 忠実度 = {fid_1a:.10f}")
    
    # 試行1b: 前後に位相回転を追加
    print("\n試行2: 前後に位相回転")
    
    # 理論的な分析:
    # G[0,0] = cos(θ/2) e^(iφ/2)
    # G[0,1] = sin(θ/2) e^(-iφ/2)
    # 
    # MQT R(θ, φ):
    # R[0,0] = cos(θ/2)
    # R[0,1] = sin(θ/2) e^(iφ)
    #
    # 位相差:
    # G[0,0] = e^(iφ/2) × cos(θ/2)
    # G[0,1] = e^(-iφ/2) × sin(θ/2)
    #
    # R[0,1] = sin(θ/2) e^(iφ) ≠ sin(θ/2) e^(-iφ/2)
    
    # 試行2a: Rz(φ/2)_0 @ R(θ, -φ) @ Rz(0)_1
    gates_2a = [
        construct_r_gate(theta, -phi),
        construct_virtrz(0, phi / 2),
        construct_virtrz(1, 0.0)
    ]
    G_2a = np.eye(3, dtype=complex)
    for gate in reversed(gates_2a):
        G_2a = gate @ G_2a
    fid_2a = compute_fidelity(G_target, G_2a)
    print(f"  2a. Rz(φ/2)_0 @ R(θ, -φ) @ Rz(0)_1: 忠実度 = {fid_2a:.10f}")
    
    # 試行2b: 両側に位相
    gates_2b = [
        construct_virtrz(0, phi / 2),
        construct_virtrz(1, -phi / 2),
        construct_r_gate(theta, 0.0)
    ]
    G_2b = np.eye(3, dtype=complex)
    for gate in reversed(gates_2b):
        G_2b = gate @ G_2b
    fid_2b = compute_fidelity(G_target, G_2b)
    print(f"  2b. Rz(φ/2)_0 @ Rz(-φ/2)_1 @ R(θ, 0): 忠実度 = {fid_2b:.10f}")
    
    # 試行3: 数値的に最適解を探索
    print("\n試行3: 数値探索")
    
    best_fid = 0.0
    best_params = None
    
    # パラメータ探索: Rz(α)_0 @ Rz(β)_1 @ R(θ_r, φ_r)
    for alpha in np.linspace(-np.pi, np.pi, 21):
        for beta in np.linspace(-np.pi, np.pi, 21):
            for theta_r in np.linspace(-np.pi, np.pi, 21):
                for phi_r in np.linspace(-np.pi, np.pi, 21):
                    gates = [
                        construct_virtrz(0, alpha),
                        construct_virtrz(1, beta),
                        construct_r_gate(theta_r, phi_r)
                    ]
                    G_test = np.eye(3, dtype=complex)
                    for gate in reversed(gates):
                        G_test = gate @ G_test
                    fid = compute_fidelity(G_target, G_test)
                    
                    if fid > best_fid:
                        best_fid = fid
                        best_params = (alpha, beta, theta_r, phi_r)
    
    if best_params:
        alpha, beta, theta_r, phi_r = best_params
        print(f"  最良の結果:")
        print(f"    α = {alpha:.6f}, β = {beta:.6f}")
        print(f"    θ_r = {theta_r:.6f}, φ_r = {phi_r:.6f}")
        print(f"    忠実度 = {best_fid:.10f}")
        
        # 検証
        gates_best = [
            construct_virtrz(0, alpha),
            construct_virtrz(1, beta),
            construct_r_gate(theta_r, phi_r)
        ]
        G_best = np.eye(3, dtype=complex)
        for gate in reversed(gates_best):
            G_best = gate @ G_best
        
        print(f"\n  再構築された行列:")
        print(G_best)
        print(f"\n  目標行列:")
        print(G_target)
        print(f"\n  差分:")
        print(G_target - G_best)


def main():
    """メイン関数"""
    
    # テストパラメータ
    theta = 0.5
    phi = 1.0
    
    # 構造分析
    analyze_givens_structure(theta, phi)
    
    # 正しい分解を探索
    find_correct_mqt_decomposition(theta, phi)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
