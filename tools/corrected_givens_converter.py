#!/usr/bin/env python3
"""
修正されたGivens変換器 (Corrected Givens Converter)

PR#41: 3×3ゲート変換の正しい実装

このツールは、診断結果に基づいて正しいGivens → MQT-Quditsゲート変換を実装します。

重要な発見:
Givens回転の定義:
  G[i,i] = c = cos(θ/2) e^(iφ/2)
  G[i,j] = s = sin(θ/2) e^(-iφ/2)
  G[j,i] = -s* = -sin(θ/2) e^(iφ/2)
  G[j,j] = c* = cos(θ/2) e^(-iφ/2)

MQT-Qudits Rゲートの定義:
  R[i,i] = cos(θ_R/2)
  R[i,j] = sin(θ_R/2) e^(iφ_R)
  R[j,i] = -sin(θ_R/2) e^(iφ_R)
  R[j,j] = cos(θ_R/2)

正しい変換:
  G(i,j; θ, φ) = VirtRz(φ/2, i) @ VirtRz(-φ/2, j) @ R(θ, -φ, i, j)

検証:
  VirtRz(φ/2, i) @ VirtRz(-φ/2, j) @ R(θ, -φ):
  - 対角要素: e^(iφ/2) × cos(θ/2) ✓
  - 非対角要素: e^(-iφ/2) × sin(θ/2) ✓
"""

import sys
from pathlib import Path
import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent))


@dataclass
class MQTGate:
    """MQT-Quditsゲートの抽象表現"""
    gate_type: str  # 'VirtRz', 'R'
    parameters: Dict
    cost: int = 0


def construct_givens_from_theta_phi(i: int, j: int, theta: float, phi: float, size: int = 3) -> np.ndarray:
    """
    θ, φからGivens行列を構築
    
    c = cos(θ/2) e^(iφ/2)
    s = sin(θ/2) e^(-iφ/2)
    """
    c = np.cos(theta / 2) * np.exp(1j * phi / 2)
    s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
    
    G = np.eye(size, dtype=complex)
    G[i, i] = c
    G[i, j] = s
    G[j, i] = -np.conj(s)
    G[j, j] = np.conj(c)
    
    return G


def construct_mqt_r_gate(level1: int, level2: int, theta: float, phi: float, size: int = 3) -> np.ndarray:
    """
    MQT-Qudits R ゲートの行列表現
    """
    R = np.eye(size, dtype=complex)
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    R[level1, level1] = c
    R[level1, level2] = s * np.exp(1j * phi)
    R[level2, level1] = -s * np.exp(1j * phi)
    R[level2, level2] = c
    return R


def construct_virtrz(level: int, phase: float, size: int = 3) -> np.ndarray:
    """VirtRz ゲートの行列表現"""
    Rz = np.eye(size, dtype=complex)
    Rz[level, level] = np.exp(1j * phase)
    return Rz


def compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
    """忠実度を計算"""
    N = U1.shape[0]
    trace = np.trace(U1.conj().T @ U2)
    return abs(trace) / N


def convert_givens_to_mqt_gates(i: int, j: int, theta: float, phi: float) -> List[MQTGate]:
    """
    Givens回転をMQT-Quditsゲートに変換
    
    Args:
        i: 第1のレベル
        j: 第2のレベル
        theta: 回転角度
        phi: 位相角度
        
    Returns:
        MQT-Quditsゲートのリスト
    """
    gates = []
    
    # VirtRz(φ/2, i)
    gates.append(MQTGate(
        gate_type='VirtRz',
        parameters={'level': i, 'phase': phi / 2},
        cost=0
    ))
    
    # VirtRz(-φ/2, j)
    gates.append(MQTGate(
        gate_type='VirtRz',
        parameters={'level': j, 'phase': -phi / 2},
        cost=0
    ))
    
    # R(θ, -φ, i, j)
    gates.append(MQTGate(
        gate_type='R',
        parameters={'level1': i, 'level2': j, 'theta': theta, 'phi': -phi},
        cost=1
    ))
    
    return gates


def verify_conversion(theta: float, phi: float, size: int = 3) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    変換の正確性を検証
    
    Returns:
        (fidelity, G_target, G_mqt)
    """
    # 目標Givens行列
    G_target = construct_givens_from_theta_phi(0, 1, theta, phi, size)
    
    # MQT-Quditsゲートに変換
    gates = convert_givens_to_mqt_gates(0, 1, theta, phi)
    
    # 行列を再構築（逆順で左から掛ける）
    G_mqt = np.eye(size, dtype=complex)
    for gate in reversed(gates):
        if gate.gate_type == 'VirtRz':
            level = gate.parameters['level']
            phase = gate.parameters['phase']
            Rz = construct_virtrz(level, phase, size)
            G_mqt = Rz @ G_mqt
        elif gate.gate_type == 'R':
            level1 = gate.parameters['level1']
            level2 = gate.parameters['level2']
            theta_r = gate.parameters['theta']
            phi_r = gate.parameters['phi']
            R = construct_mqt_r_gate(level1, level2, theta_r, phi_r, size)
            G_mqt = R @ G_mqt
    
    # 忠実度を計算
    fidelity = compute_fidelity(G_target, G_mqt)
    
    return fidelity, G_target, G_mqt


def test_single_givens():
    """単一のGivens回転でテスト"""
    print("="*70)
    print("単一Givens回転のテスト")
    print("="*70)
    
    # テストケース
    test_cases = [
        (0.5, 1.0),
        (0.1, 0.5),
        (1.0, 2.0),
        (np.pi/4, np.pi/6),
        (0.0, 1.0),  # θ=0のエッジケース
        (np.pi, 0.0),  # θ=πのエッジケース
    ]
    
    all_passed = True
    
    for idx, (theta, phi) in enumerate(test_cases, 1):
        print(f"\nテストケース {idx}: θ={theta:.6f}, φ={phi:.6f}")
        
        fidelity, G_target, G_mqt = verify_conversion(theta, phi)
        
        passed = fidelity > 0.9999
        status = "✓" if passed else "✗"
        
        print(f"  忠実度: {fidelity:.10f} {status}")
        
        if not passed:
            all_passed = False
            print(f"\n  目標行列:")
            print(f"  {G_target}")
            print(f"\n  再構築行列:")
            print(f"  {G_mqt}")
            print(f"\n  差分:")
            print(f"  {G_target - G_mqt}")
    
    print("\n" + "="*70)
    if all_passed:
        print("✓✓✓ すべてのテストに合格")
    else:
        print("⚠⚠⚠ 一部のテストが不合格")
    print("="*70)
    
    return all_passed


def test_with_integrated_sparse_compiler():
    """integrated_sparse_compiler.pyとの統合テスト"""
    print("\n" + "="*70)
    print("integrated_sparse_compiler.pyとの統合テスト")
    print("="*70)
    
    try:
        from integrated_sparse_compiler import IntegratedThreeLevelDecomposer
    except ImportError as e:
        print(f"ERROR: integrated_sparse_compiler.pyをインポートできません: {e}")
        return False
    
    # テストユニタリ（H_TTA風）
    # 簡単な3×3ユニタリ
    U = np.array([
        [0.8+0.1j, 0.5-0.2j, 0.1+0.1j],
        [-0.3+0.2j, 0.6+0.3j, 0.5-0.1j],
        [0.2-0.3j, -0.4+0.2j, 0.7+0.2j]
    ])
    # 正規化してユニタリにする
    U, _ = np.linalg.qr(U)
    
    print(f"\nテストユニタリ行列:")
    print(U)
    
    # integrated_sparse_compiler.pyで分解
    decomposer = IntegratedThreeLevelDecomposer()
    result = decomposer.decompose(U)
    
    print(f"\n分解結果:")
    print(f"  QR忠実度: {result['fidelity']:.10f}")
    print(f"  Givens回転数: {len(result['rotations'])}")
    
    # 各Givens回転をMQT-Quditsゲートに変換
    all_gates = []
    Q_reconstruct = np.eye(3, dtype=complex)
    
    for idx, (i, j, theta, phi) in enumerate(result['rotations']):
        print(f"\n  Givens {idx+1}: G({i},{j}; θ={theta:.6f}, φ={phi:.6f})")
        
        # MQT-Quditsゲートに変換
        gates = convert_givens_to_mqt_gates(i, j, theta, phi)
        all_gates.extend(gates)
        
        # 検証: 単一のGivens回転
        fid_single, G_target, G_mqt = verify_conversion(theta, phi)
        print(f"    単一変換忠実度: {fid_single:.10f}")
        
        # Qを再構築
        G = construct_givens_from_theta_phi(i, j, theta, phi, size=3)
        Q_reconstruct = G @ Q_reconstruct
    
    # 対角位相を追加
    if 'diagonal_phases' in result:
        for level, phase in enumerate(result['diagonal_phases']):
            if abs(phase) > 1e-10:
                all_gates.append(MQTGate(
                    gate_type='VirtRz',
                    parameters={'level': level, 'phase': phase},
                    cost=0
                ))
    
    print(f"\n総MQT-Quditsゲート数: {len(all_gates)}")
    print(f"  物理ゲート数: {sum(1 for g in all_gates if g.cost > 0)}")
    
    # 全体の忠実度を検証
    # U = Q @ D @ R_norm (integrated_sparse_compiler.pyの結果)
    # ここではQの部分のみ検証
    Q = result['Q']
    fid_q = compute_fidelity(Q, Q_reconstruct)
    print(f"\nQ行列の再構築忠実度: {fid_q:.10f}")
    
    passed = fid_q > 0.9999
    if passed:
        print("✓ 統合テスト合格")
    else:
        print("✗ 統合テスト不合格")
        print(f"\n元のQ:")
        print(Q)
        print(f"\n再構築されたQ:")
        print(Q_reconstruct)
    
    return passed


def main():
    """メイン関数"""
    print("="*70)
    print("修正されたGivens変換器のテスト")
    print("="*70)
    print("\nPR#41: 3×3ゲート変換の正しい実装")
    print("\n変換式: G(i,j; θ, φ) = VirtRz(φ/2, i) @ VirtRz(-φ/2, j) @ R(θ, -φ, i, j)")
    
    # 単一Givens回転のテスト
    success1 = test_single_givens()
    
    # integrated_sparse_compiler.pyとの統合テスト
    success2 = test_with_integrated_sparse_compiler()
    
    print("\n" + "="*70)
    print("最終結果")
    print("="*70)
    print(f"単一Givens回転: {'✓ 合格' if success1 else '✗ 不合格'}")
    print(f"統合テスト: {'✓ 合格' if success2 else '✗ 不合格'}")
    
    if success1 and success2:
        print("\n✓✓✓ すべてのテストに合格")
        print("\n次のステップ:")
        print("1. gate_converter.pyのThreeLevelGateConverterを更新")
        print("2. H_TTAでの検証")
        print("3. ランダムユニタリテスト（100個）")
        return 0
    else:
        print("\n⚠⚠⚠ 一部のテストが不合格")
        print("\n修正が必要な箇所:")
        if not success1:
            print("- 単一Givens回転の変換式を再確認")
        if not success2:
            print("- integrated_sparse_compiler.pyからの抽出方法を再確認")
        return 1


if __name__ == '__main__':
    sys.exit(main())
