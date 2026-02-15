#!/usr/bin/env python3
"""
Givens回転診断ツール (Givens Rotation Diagnostic Tool)

PR#41: 3×3ゲート変換の問題を診断するツール

このツールは、3×3 Givens回転のMQT-Quditsゲートへの変換が
なぜ忠実度 0.68 しか達成できないのかを診断します。

診断項目:
1. Givens回転の定義の確認
2. ZYZ分解との対応の確認
3. MQT-Quditsゲートへの変換の確認
4. integrated_sparse_compiler.pyの抽出方法の検証

理論的根拠:
- givens_rotation_theory_ja.md の定義
- integrated_sparse_compiler.py の実装
- gate_converter.py の変換ロジック

数学的厳密性:
- すべてのテストは厳密な線形代数に基づく
- ヒューリスティックゼロ
- 近似ゼロ
- 忠実度 > 0.9999 を要求
"""

import sys
from pathlib import Path
import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass

# tools/ディレクトリからのインポートを有効化
sys.path.insert(0, str(Path(__file__).parent))


def compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
    """
    2つのユニタリ行列の忠実度を計算
    
    F = |Tr(U1† @ U2)| / N
    
    Args:
        U1: 第1のユニタリ行列
        U2: 第2のユニタリ行列
        
    Returns:
        忠実度 (0.0 ~ 1.0)
    """
    N = U1.shape[0]
    trace = np.trace(U1.conj().T @ U2)
    return abs(trace) / N


def construct_givens_theory(i: int, j: int, theta: float, phi: float, size: int = 3) -> np.ndarray:
    """
    givens_rotation_theory_ja.md の定義に基づくGivens行列
    
    G(i,j; θ, φ) = I + (c-1)(|i⟩⟨i| + |j⟩⟨j|) - s*|i⟩⟨j| + s|j⟩⟨i|
    
    ここで:
    c = cos(θ/2) e^(iφ/2)
    s = sin(θ/2) e^(-iφ/2)
    
    Args:
        i: 第1のレベル
        j: 第2のレベル
        theta: 回転角度
        phi: 位相角度
        size: 行列サイズ（デフォルト: 3）
        
    Returns:
        Givens行列 (size × size)
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
    
    R(θ, φ) = [[cos(θ/2), sin(θ/2)e^(iφ)],
               [-sin(θ/2)e^(iφ), cos(θ/2)]]
    
    Args:
        level1: 第1のレベル
        level2: 第2のレベル
        theta: 回転角度
        phi: 位相角度
        size: 行列サイズ（デフォルト: 3）
        
    Returns:
        R ゲートの行列表現 (size × size)
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
    """
    VirtRz ゲートの行列表現
    
    VirtRz(φ) = I の level 番目の対角要素に e^(iφ) を適用
    
    Args:
        level: 対象レベル
        phase: 位相
        size: 行列サイズ（デフォルト: 3）
        
    Returns:
        VirtRz ゲートの行列表現 (size × size)
    """
    Rz = np.eye(size, dtype=complex)
    Rz[level, level] = np.exp(1j * phase)
    return Rz


@dataclass
class DiagnosticResult:
    """診断結果"""
    test_name: str
    passed: bool
    fidelity: float
    details: str


class GivensDiagnostic:
    """Givens回転診断器"""
    
    def __init__(self, tolerance: float = 1e-10):
        """
        Args:
            tolerance: 数値誤差の許容範囲
        """
        self.tolerance = tolerance
        self.results: List[DiagnosticResult] = []
    
    def test_givens_definition(self) -> DiagnosticResult:
        """
        テスト1: Givens回転の定義を確認
        
        integrated_sparse_compiler.pyの抽出結果と
        givens_rotation_theory_ja.mdの定義を比較
        """
        print("\n" + "="*70)
        print("テスト1: Givens回転の定義を確認")
        print("="*70)
        
        # テストパラメータ
        theta = 0.5
        phi = 1.0
        
        # 理論的なGivens行列
        G_theory = construct_givens_theory(0, 1, theta, phi, size=3)
        
        print(f"\nパラメータ:")
        print(f"  θ = {theta:.6f}")
        print(f"  φ = {phi:.6f}")
        print(f"\n理論的なGivens行列 G(0,1; θ, φ):")
        print(G_theory)
        
        # ユニタリ性の確認
        G_H_G = G_theory.conj().T @ G_theory
        unitarity_error = np.max(np.abs(G_H_G - np.eye(3)))
        print(f"\nユニタリ性検証: ||G†G - I||_max = {unitarity_error:.2e}")
        
        # integrated_sparse_compiler.pyを使って確認
        try:
            from integrated_sparse_compiler import IntegratedThreeLevelDecomposer
            
            # 単純な3×3ユニタリでテスト
            U = G_theory.copy()
            
            decomposer = IntegratedThreeLevelDecomposer()
            result = decomposer.decompose(U)
            
            print(f"\nintegrated_sparse_compiler.pyによる分解:")
            print(f"  QR忠実度: {result['fidelity']:.10f}")
            print(f"  Givens回転数: {len(result['rotations'])}")
            
            if len(result['rotations']) > 0:
                print(f"\n抽出されたGivensパラメータ:")
                for idx, (i, j, t, p) in enumerate(result['rotations']):
                    print(f"  回転{idx+1}: G({i},{j}; θ={t:.6f}, φ={p:.6f})")
            
            fidelity = result['fidelity']
            passed = fidelity > 0.9999
            
        except ImportError as e:
            print(f"\nERROR: integrated_sparse_compiler.pyをインポートできません: {e}")
            fidelity = 0.0
            passed = False
        
        result = DiagnosticResult(
            test_name="Givens定義検証",
            passed=passed,
            fidelity=fidelity,
            details=f"ユニタリ性誤差: {unitarity_error:.2e}"
        )
        
        self.results.append(result)
        return result
    
    def test_givens_zyz_decomposition(self) -> DiagnosticResult:
        """
        テスト2: Givens回転のZYZ分解を確認
        
        理論的な関係:
        G(i,j; θ, φ) = Rz(φ/2)_i @ Ry(θ)_{i,j} @ Rz(-φ/2)_j
        """
        print("\n" + "="*70)
        print("テスト2: Givens回転のZYZ分解を確認")
        print("="*70)
        
        # テストパラメータ
        theta = 0.5
        phi = 1.0
        
        # 理論的なGivens
        G_theory = construct_givens_theory(0, 1, theta, phi, size=3)
        
        print(f"\nパラメータ:")
        print(f"  θ = {theta:.6f}")
        print(f"  φ = {phi:.6f}")
        
        # ZYZ分解: G = Rz(φ/2)_0 @ Ry(θ)_{0,1} @ Rz(-φ/2)_1
        
        # Rz(φ/2) on level 0
        Rz_phi_2 = np.diag([np.exp(1j * phi / 2), 1.0, 1.0])
        
        # Ry(θ) on levels (0,1) - 標準Ryゲート
        cos_t2 = np.cos(theta / 2)
        sin_t2 = np.sin(theta / 2)
        Ry = np.array([
            [cos_t2, -sin_t2, 0],
            [sin_t2, cos_t2, 0],
            [0, 0, 1]
        ], dtype=complex)
        
        # Rz(-φ/2) on level 1
        Rz_minus_phi_2 = np.diag([1.0, np.exp(-1j * phi / 2), 1.0])
        
        # 積
        G_zyz = Rz_phi_2 @ Ry @ Rz_minus_phi_2
        
        print(f"\nZYZ分解:")
        print(f"  Rz(φ/2)_0:")
        print(Rz_phi_2)
        print(f"\n  Ry(θ)_(0,1):")
        print(Ry)
        print(f"\n  Rz(-φ/2)_1:")
        print(Rz_minus_phi_2)
        print(f"\n  積 G = Rz(φ/2)_0 @ Ry(θ) @ Rz(-φ/2)_1:")
        print(G_zyz)
        
        # 忠実度
        fidelity = compute_fidelity(G_theory, G_zyz)
        print(f"\n忠実度: {fidelity:.10f}")
        
        passed = fidelity > 0.9999
        
        if not passed:
            print("\nERROR: ZYZ分解が理論値と一致しません")
            print(f"G_theory:\n{G_theory}")
            print(f"G_zyz:\n{G_zyz}")
            print(f"差分:\n{G_theory - G_zyz}")
        
        result = DiagnosticResult(
            test_name="Givens ZYZ分解検証",
            passed=passed,
            fidelity=fidelity,
            details=f"Rz(φ/2)_i @ Ry(θ) @ Rz(-φ/2)_j"
        )
        
        self.results.append(result)
        return result
    
    def test_givens_to_mqt_gates_standard_ry(self) -> DiagnosticResult:
        """
        テスト3a: Givens → MQT-Quditsゲート変換（標準Ry使用）
        
        標準Ry(θ)を使った場合の変換をテスト
        """
        print("\n" + "="*70)
        print("テスト3a: Givens → MQT-Quditsゲート（標準Ry）")
        print("="*70)
        
        # テストパラメータ
        theta = 0.5
        phi = 1.0
        
        # 理論的なGivens
        G_theory = construct_givens_theory(0, 1, theta, phi, size=3)
        
        print(f"\nパラメータ:")
        print(f"  θ = {theta:.6f}")
        print(f"  φ = {phi:.6f}")
        
        # MQT-Quditsゲート変換（標準Ryを使用）
        # G = Rz(φ/2)_0 @ Ry(θ)_{0,1} @ Rz(-φ/2)_1
        # 
        # 標準 Ry(θ) = [[cos(θ/2), -sin(θ/2)],
        #               [sin(θ/2),  cos(θ/2)]]
        
        # ゲートシーケンス
        gates = [
            ('VirtRz', 0, phi / 2),
            ('Ry_standard', 0, 1, theta),  # 標準Ry
            ('VirtRz', 1, -phi / 2)
        ]
        
        # 再構築（逆順で左から掛ける）
        G_mqt = np.eye(3, dtype=complex)
        for gate in reversed(gates):
            if gate[0] == 'VirtRz':
                level, phase = gate[1], gate[2]
                Rz = construct_virtrz(level, phase, size=3)
                G_mqt = Rz @ G_mqt
            elif gate[0] == 'Ry_standard':
                level1, level2, theta_g = gate[1], gate[2], gate[3]
                # 標準Ry
                cos_t2 = np.cos(theta_g / 2)
                sin_t2 = np.sin(theta_g / 2)
                Ry = np.eye(3, dtype=complex)
                Ry[level1, level1] = cos_t2
                Ry[level1, level2] = -sin_t2
                Ry[level2, level1] = sin_t2
                Ry[level2, level2] = cos_t2
                G_mqt = Ry @ G_mqt
        
        print(f"\n再構築されたGivens行列 G_mqt:")
        print(G_mqt)
        
        # 忠実度
        fidelity = compute_fidelity(G_theory, G_mqt)
        print(f"\n忠実度: {fidelity:.10f}")
        
        passed = fidelity > 0.9999
        
        result = DiagnosticResult(
            test_name="Givens → MQT (標準Ry)",
            passed=passed,
            fidelity=fidelity,
            details="VirtRz(φ/2)_i + Ry(θ) + VirtRz(-φ/2)_j"
        )
        
        self.results.append(result)
        return result
    
    def test_givens_to_mqt_gates_r_negative_theta(self) -> DiagnosticResult:
        """
        テスト3b: Givens → MQT-Quditsゲート変換（R(-θ, 0)使用）
        
        MQT-Qudits R(-θ, 0)を使った場合の変換をテスト
        （gate_converter.pyの現在の実装）
        """
        print("\n" + "="*70)
        print("テスト3b: Givens → MQT-Quditsゲート（R(-θ, 0)）")
        print("="*70)
        
        # テストパラメータ
        theta = 0.5
        phi = 1.0
        
        # 理論的なGivens
        G_theory = construct_givens_theory(0, 1, theta, phi, size=3)
        
        print(f"\nパラメータ:")
        print(f"  θ = {theta:.6f}")
        print(f"  φ = {phi:.6f}")
        
        # MQT-Quditsゲート変換（R(-θ, 0)を使用）
        # Ry(θ) = R(-θ, 0)
        
        # ゲートシーケンス
        gates = [
            ('VirtRz', 0, phi / 2),
            ('R', 0, 1, -theta, 0.0),  # R(-θ, 0)
            ('VirtRz', 1, -phi / 2)
        ]
        
        # 再構築（逆順で左から掛ける）
        G_mqt = np.eye(3, dtype=complex)
        for gate in reversed(gates):
            if gate[0] == 'VirtRz':
                level, phase = gate[1], gate[2]
                Rz = construct_virtrz(level, phase, size=3)
                G_mqt = Rz @ G_mqt
            elif gate[0] == 'R':
                level1, level2, theta_r, phi_r = gate[1], gate[2], gate[3], gate[4]
                R = construct_mqt_r_gate(level1, level2, theta_r, phi_r, size=3)
                G_mqt = R @ G_mqt
        
        print(f"\n再構築されたGivens行列 G_mqt:")
        print(G_mqt)
        
        # 忠実度
        fidelity = compute_fidelity(G_theory, G_mqt)
        print(f"\n忠実度: {fidelity:.10f}")
        
        passed = fidelity > 0.9999
        
        result = DiagnosticResult(
            test_name="Givens → MQT (R(-θ, 0))",
            passed=passed,
            fidelity=fidelity,
            details="VirtRz(φ/2)_i + R(-θ, 0) + VirtRz(-φ/2)_j"
        )
        
        self.results.append(result)
        return result
    
    def test_givens_to_mqt_gates_r_positive_theta(self) -> DiagnosticResult:
        """
        テスト3c: Givens → MQT-Quditsゲート変換（R(θ, 0)使用）
        
        MQT-Qudits R(θ, 0)を使った場合の変換をテスト
        """
        print("\n" + "="*70)
        print("テスト3c: Givens → MQT-Quditsゲート（R(θ, 0)）")
        print("="*70)
        
        # テストパラメータ
        theta = 0.5
        phi = 1.0
        
        # 理論的なGivens
        G_theory = construct_givens_theory(0, 1, theta, phi, size=3)
        
        print(f"\nパラメータ:")
        print(f"  θ = {theta:.6f}")
        print(f"  φ = {phi:.6f}")
        
        # MQT-Quditsゲート変換（R(θ, 0)を使用）
        
        # ゲートシーケンス
        gates = [
            ('VirtRz', 0, phi / 2),
            ('R', 0, 1, theta, 0.0),  # R(θ, 0) - 符号そのまま
            ('VirtRz', 1, -phi / 2)
        ]
        
        # 再構築（逆順で左から掛ける）
        G_mqt = np.eye(3, dtype=complex)
        for gate in reversed(gates):
            if gate[0] == 'VirtRz':
                level, phase = gate[1], gate[2]
                Rz = construct_virtrz(level, phase, size=3)
                G_mqt = Rz @ G_mqt
            elif gate[0] == 'R':
                level1, level2, theta_r, phi_r = gate[1], gate[2], gate[3], gate[4]
                R = construct_mqt_r_gate(level1, level2, theta_r, phi_r, size=3)
                G_mqt = R @ G_mqt
        
        print(f"\n再構築されたGivens行列 G_mqt:")
        print(G_mqt)
        
        # 忠実度
        fidelity = compute_fidelity(G_theory, G_mqt)
        print(f"\n忠実度: {fidelity:.10f}")
        
        passed = fidelity > 0.9999
        
        result = DiagnosticResult(
            test_name="Givens → MQT (R(θ, 0))",
            passed=passed,
            fidelity=fidelity,
            details="VirtRz(φ/2)_i + R(θ, 0) + VirtRz(-φ/2)_j"
        )
        
        self.results.append(result)
        return result
    
    def run_all_tests(self) -> bool:
        """すべての診断テストを実行"""
        print("\n" + "="*70)
        print("Givens回転診断ツール")
        print("="*70)
        print("\nPR#41: 3×3ゲート変換の問題を診断")
        print("\n目的: 忠実度 0.68 → 1.0 を達成するための問題特定")
        
        # テスト実行
        self.test_givens_definition()
        self.test_givens_zyz_decomposition()
        self.test_givens_to_mqt_gates_standard_ry()
        self.test_givens_to_mqt_gates_r_negative_theta()
        self.test_givens_to_mqt_gates_r_positive_theta()
        
        # サマリー
        print("\n" + "="*70)
        print("診断結果サマリー")
        print("="*70)
        
        all_passed = True
        for result in self.results:
            status = "✓" if result.passed else "✗"
            print(f"\n{status} {result.test_name}:")
            print(f"  忠実度: {result.fidelity:.10f}")
            print(f"  詳細: {result.details}")
            if not result.passed:
                all_passed = False
        
        print("\n" + "="*70)
        if all_passed:
            print("✓✓✓ すべての診断テストに合格")
        else:
            print("⚠⚠⚠ 一部の診断テストが不合格")
            print("\n推奨アクション:")
            for result in self.results:
                if not result.passed:
                    print(f"  - {result.test_name}の問題を修正")
        print("="*70)
        
        return all_passed


def main():
    """メイン関数"""
    diagnostic = GivensDiagnostic()
    success = diagnostic.run_all_tests()
    
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
