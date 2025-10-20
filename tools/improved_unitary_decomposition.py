#!/usr/bin/env python3
"""
改良版厳密ユニタリ分解ツール (Improved Rigorous Unitary Decomposition)

このモジュールは、既存のunitary_decomposition_rigorous.pyを分析し、
数学的に正確な2×2および3×3ユニタリ分解を実装します。

主な改良点:
1. 2×2 ZYZ分解の数学的に正しいパラメータ抽出
2. 3×3 Givens分解の数値安定性の向上
3. 忠実度 > 0.9999 の達成

理論的根拠:
- tutorials/doc/rigorous_unitary_decomposition_theory_ja.md に完全に記載
- すべての実装はヒューリスティックや近似を使用しない
- scipy.linalg.expmは使用しない（Padé近似を含むため）
"""

import numpy as np
from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass


@dataclass
class ImprovedDecomposition2x2:
    """改良版2×2ユニタリ分解の結果"""
    theta: float  # Y軸回転角
    phi: float    # 第1のZ軸回転角
    lam: float    # 第2のZ軸回転角
    global_phase: float  # グローバル位相
    fidelity: float  # 元の行列との忠実度
    method: str = 'ZYZ_improved'
    
    def __post_init__(self):
        """検証: 忠実度が要求を満たしているか"""
        if self.fidelity < 0.9999:
            import warnings
            warnings.warn(
                f"忠実度 {self.fidelity:.6f} が要求値 0.9999 を下回っています",
                UserWarning
            )


@dataclass
class ImprovedDecomposition3x3:
    """改良版3×3ユニタリ分解の結果"""
    rotations: List[Tuple[int, int, float, float]]  # (level1, level2, theta, phi)
    diagonal_phases: np.ndarray  # 対角位相
    fidelity: float  # 元の行列との忠実度
    method: str = 'Givens_improved'
    
    def __post_init__(self):
        """検証: 忠実度が要求を満たしているか"""
        if self.fidelity < 0.9999:
            import warnings
            warnings.warn(
                f"忠実度 {self.fidelity:.6f} が要求値 0.9999 を下回っています",
                UserWarning
            )


class ImprovedTwoQubitDecomposer:
    """
    改良版2×2ユニタリ分解器
    
    数学的に正確なZYZ分解を実装します。
    参考: Qiskitの実装とrigorous_unitary_decomposition_theory_ja.md
    """
    
    @staticmethod
    def decompose_zyz(U: np.ndarray, tolerance: float = 1e-10) -> ImprovedDecomposition2x2:
        """
        改良されたZYZ分解
        
        U = e^(iα) Rz(φ) Ry(θ) Rz(λ)
        
        数学的根拠:
        - SU(2)の一般形: U_SU2 = [[e^(i(φ+λ)/2)cos(θ/2), -e^(i(φ-λ)/2)sin(θ/2)],
                                    [e^(i(λ-φ)/2)sin(θ/2),  e^(-i(φ+λ)/2)cos(θ/2)]]
        
        Args:
            U: 2×2ユニタリ行列
            tolerance: 数値誤差の許容範囲
            
        Returns:
            ImprovedDecomposition2x2
        """
        if U.shape != (2, 2):
            raise ValueError(f"2×2行列である必要があります。入力形状: {U.shape}")
        
        # ユニタリ性の検証
        if not ImprovedTwoQubitDecomposer._is_unitary(U, tolerance):
            raise ValueError("入力行列はユニタリではありません")
        
        # Step 1: グローバル位相を抽出
        # det(U) = e^(2iα) なので α = arg(det(U)) / 2
        det_U = np.linalg.det(U)
        global_phase = np.angle(det_U) / 2.0
        
        # Step 2: SU(2)に正規化
        U_su2 = U * np.exp(-1j * global_phase)
        
        # Step 3: 要素を取得
        a = U_su2[0, 0]  # e^(i(φ+λ)/2) cos(θ/2)
        b = U_su2[0, 1]  # -e^(i(φ-λ)/2) sin(θ/2)
        c = U_su2[1, 0]  # e^(i(λ-φ)/2) sin(θ/2)
        d = U_su2[1, 1]  # e^(-i(φ+λ)/2) cos(θ/2)
        
        # Step 4: θを計算
        # |a|^2 = cos^2(θ/2) より cos(θ/2) = |a|
        cos_theta_2 = abs(a)
        cos_theta_2 = np.clip(cos_theta_2, 0.0, 1.0)  # 数値誤差対策
        theta = 2.0 * np.arccos(cos_theta_2)
        
        # Step 5: φとλを計算
        sin_theta_2 = np.sin(theta / 2.0)
        
        if sin_theta_2 > tolerance:
            # 一般的なケース
            # a = e^(i(φ+λ)/2) cos(θ/2)
            # d = e^(-i(φ+λ)/2) cos(θ/2)
            # a / d = e^(i(φ+λ))
            # したがって φ+λ = arg(a / d) = arg(a) - arg(d)
            
            # c = e^(i(λ-φ)/2) sin(θ/2)
            # -b = e^(i(φ-λ)/2) sin(θ/2)
            # c / (-b) = e^(i(λ-φ))
            # したがって λ-φ = arg(c / (-b)) = arg(c) - arg(-b)
            
            phi_plus_lambda = np.angle(a) - np.angle(d)
            lambda_minus_phi = np.angle(c) - np.angle(-b)
            
            # 連立方程式を解く
            # φ + λ = phi_plus_lambda
            # λ - φ = lambda_minus_phi
            # => 2φ = phi_plus_lambda - lambda_minus_phi
            # => 2λ = phi_plus_lambda + lambda_minus_phi
            
            phi = (phi_plus_lambda - lambda_minus_phi) / 2.0
            lam = (phi_plus_lambda + lambda_minus_phi) / 2.0
            
        else:
            # 特異点: θ ≈ 0 または π
            if abs(cos_theta_2 - 1.0) < tolerance:
                # θ ≈ 0: U ≈ e^(i(φ+λ)/2) I
                # この場合、φとλは個別に定義できず、和のみが意味を持つ
                # 慣例として φ = 0 と設定
                phi = 0.0
                lam = np.angle(a) - np.angle(d)
            else:
                # θ ≈ π: U ≈ e^(i(φ+λ)/2) [[0, -e^(i(φ-λ)/2)],
                #                              [e^(i(λ-φ)/2), 0]]
                # この場合も φとλの和のみが意味を持つ
                phi = 0.0
                lam = np.angle(-b) - np.angle(c)
        
        # 正規化（-π ~ π）
        phi = np.angle(np.exp(1j * phi))
        lam = np.angle(np.exp(1j * lam))
        
        # Step 6: 再構築して忠実度を計算
        U_reconstructed = ImprovedTwoQubitDecomposer._construct_zyz(
            theta, phi, lam, global_phase
        )
        fidelity = ImprovedTwoQubitDecomposer._compute_fidelity(U, U_reconstructed)
        
        return ImprovedDecomposition2x2(
            theta=theta,
            phi=phi,
            lam=lam,
            global_phase=global_phase,
            fidelity=fidelity,
            method='ZYZ_improved'
        )
    
    @staticmethod
    def _construct_zyz(theta: float, phi: float, lam: float, 
                       global_phase: float = 0.0) -> np.ndarray:
        """
        ZYZ分解からユニタリ行列を構築
        
        U = e^(iα) Rz(φ) Ry(θ) Rz(λ)
        
        Args:
            theta: Y軸回転角
            phi: 第1のZ軸回転角
            lam: 第2のZ軸回転角
            global_phase: グローバル位相
            
        Returns:
            2×2ユニタリ行列
        """
        # Rz(φ) = [[e^(iφ/2), 0], [0, e^(-iφ/2)]]
        Rz_phi = np.array([
            [np.exp(1j * phi / 2), 0],
            [0, np.exp(-1j * phi / 2)]
        ], dtype=complex)
        
        # Ry(θ) = [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
        Ry_theta = np.array([
            [np.cos(theta / 2), -np.sin(theta / 2)],
            [np.sin(theta / 2), np.cos(theta / 2)]
        ], dtype=complex)
        
        # Rz(λ) = [[e^(iλ/2), 0], [0, e^(-iλ/2)]]
        Rz_lam = np.array([
            [np.exp(1j * lam / 2), 0],
            [0, np.exp(-1j * lam / 2)]
        ], dtype=complex)
        
        # U = e^(iα) Rz(φ) @ Ry(θ) @ Rz(λ)
        U = Rz_phi @ Ry_theta @ Rz_lam
        U = U * np.exp(1j * global_phase)
        
        return U
    
    @staticmethod
    def _is_unitary(U: np.ndarray, tolerance: float = 1e-10) -> bool:
        """
        行列がユニタリであるかを検証
        
        U† U = I を確認
        """
        d = U.shape[0]
        identity = np.eye(d, dtype=complex)
        product = U.conj().T @ U
        return np.allclose(product, identity, atol=tolerance)
    
    @staticmethod
    def _compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
        """
        2つのユニタリ行列の忠実度を計算
        
        F = |Tr(U1† U2)| / d
        
        グローバル位相を無視した忠実度
        """
        d = U1.shape[0]
        trace = np.trace(U1.conj().T @ U2)
        fidelity = abs(trace) / d
        return float(fidelity)


class ImprovedThreeQuditDecomposer:
    """
    改良版3×3ユニタリ分解器
    
    数値的に安定したGivens分解を実装します。
    """
    
    @staticmethod
    def decompose_givens(U: np.ndarray, tolerance: float = 1e-10) -> ImprovedDecomposition3x3:
        """
        改良されたGivens分解
        
        3×3ユニタリを2準位回転の列に分解します。
        
        Args:
            U: 3×3ユニタリ行列
            tolerance: 数値誤差の許容範囲
            
        Returns:
            ImprovedDecomposition3x3
        """
        if U.shape != (3, 3):
            raise ValueError(f"3×3行列である必要があります。入力形状: {U.shape}")
        
        # ユニタリ性の検証
        if not ImprovedThreeQuditDecomposer._is_unitary(U, tolerance):
            raise ValueError("入力行列はユニタリではありません")
        
        rotations = []
        U_work = U.copy()
        
        # Givens分解: QR分解ベース
        # 目標: U_work を上三角行列にする
        
        # Step 1: U[1,0]をゼロにする（準位0と1の回転）
        if abs(U_work[1, 0]) > tolerance:
            theta, phi = ImprovedThreeQuditDecomposer._compute_givens_params(
                U_work[0, 0], U_work[1, 0]
            )
            rotations.append((0, 1, theta, phi))
            
            # Givens回転を適用
            G = ImprovedThreeQuditDecomposer._construct_givens_matrix(
                3, 0, 1, theta, phi
            )
            U_work = G.conj().T @ U_work
        
        # Step 2: U[2,0]をゼロにする（準位0と2の回転）
        if abs(U_work[2, 0]) > tolerance:
            theta, phi = ImprovedThreeQuditDecomposer._compute_givens_params(
                U_work[0, 0], U_work[2, 0]
            )
            rotations.append((0, 2, theta, phi))
            
            G = ImprovedThreeQuditDecomposer._construct_givens_matrix(
                3, 0, 2, theta, phi
            )
            U_work = G.conj().T @ U_work
        
        # Step 3: U[2,1]をゼロにする（準位1と2の回転）
        if abs(U_work[2, 1]) > tolerance:
            theta, phi = ImprovedThreeQuditDecomposer._compute_givens_params(
                U_work[1, 1], U_work[2, 1]
            )
            rotations.append((1, 2, theta, phi))
            
            G = ImprovedThreeQuditDecomposer._construct_givens_matrix(
                3, 1, 2, theta, phi
            )
            U_work = G.conj().T @ U_work
        
        # Step 4: 対角位相を抽出
        # U_workは今や上三角（実際には対角+上三角の小さな要素）
        diagonal_phases = np.angle(np.diag(U_work))
        
        # Step 5: 再構築して忠実度を計算
        U_reconstructed = ImprovedThreeQuditDecomposer._reconstruct(
            rotations, diagonal_phases
        )
        fidelity = ImprovedThreeQuditDecomposer._compute_fidelity(U, U_reconstructed)
        
        return ImprovedDecomposition3x3(
            rotations=rotations,
            diagonal_phases=diagonal_phases,
            fidelity=fidelity,
            method='Givens_improved'
        )
    
    @staticmethod
    def _compute_givens_params(a: complex, b: complex) -> Tuple[float, float]:
        """
        Givens回転のパラメータを計算
        
        目標: G @ [a, b]^T = [r, 0]^T
        
        G = [[c, -s*],
             [s,  c ]]
        
        ここで:
        c = cos(θ/2) e^(iφ/2)
        s = sin(θ/2) e^(-iφ/2)
        
        Args:
            a, b: 複素数要素
            
        Returns:
            (theta, phi): 回転角と位相
        """
        r = np.sqrt(abs(a)**2 + abs(b)**2)
        
        if r < 1e-15:
            return 0.0, 0.0
        
        # 正規化
        a_norm = a / r
        b_norm = b / r
        
        # θの計算
        # |c|^2 = cos^2(θ/2)
        # |s|^2 = sin^2(θ/2)
        # c a_norm - s* b_norm = 1
        # より tan(θ/2) = |b_norm| / |a_norm|
        
        theta = 2.0 * np.arctan2(abs(b_norm), abs(a_norm))
        
        # φの計算
        if abs(b_norm) > 1e-10:
            # c / s* = (a_norm) / (b_norm*) / tan(θ/2)
            # arg(c) - arg(s*) = arg(a_norm) - arg(-b_norm)
            # 2 arg(c) = φ (c = cos(θ/2) e^(iφ/2) より)
            # 2 arg(s*) = -φ (s = sin(θ/2) e^(-iφ/2) より)
            # したがって arg(c) - arg(s*) = φ
            
            # より直接的な方法:
            # c = a_norm / cos(θ/2) の位相部分
            # s = -b_norm / sin(θ/2) の位相部分
            
            cos_theta_2 = np.cos(theta / 2)
            sin_theta_2 = np.sin(theta / 2)
            
            if abs(cos_theta_2) > 1e-10:
                c_phase = np.angle(a_norm) - np.angle(cos_theta_2)
            else:
                c_phase = 0.0
            
            if abs(sin_theta_2) > 1e-10:
                s_phase = np.angle(-b_norm) - np.angle(sin_theta_2)
            else:
                s_phase = 0.0
            
            # c = cos(θ/2) e^(iφ/2) => arg(c) = φ/2
            # s = sin(θ/2) e^(-iφ/2) => arg(s) = -φ/2
            # したがって φ = 2 arg(c) = -2 arg(s)
            
            phi = c_phase * 2
            phi = np.angle(np.exp(1j * phi))  # 正規化
        else:
            phi = 0.0
        
        return theta, phi
    
    @staticmethod
    def _construct_givens_matrix(d: int, level1: int, level2: int,
                                 theta: float, phi: float) -> np.ndarray:
        """
        Givens回転行列を構築
        
        G = I + (c-1)(|level1⟩⟨level1| + |level2⟩⟨level2|)
            - s* |level1⟩⟨level2|
            + s |level2⟩⟨level1|
        
        Args:
            d: 行列の次元
            level1, level2: 回転する準位
            theta: 回転角
            phi: 位相
            
        Returns:
            Givens回転行列
        """
        G = np.eye(d, dtype=complex)
        
        c = np.cos(theta / 2) * np.exp(1j * phi / 2)
        s = np.sin(theta / 2) * np.exp(-1j * phi / 2)
        
        G[level1, level1] = c
        G[level1, level2] = -s.conj()
        G[level2, level1] = s
        G[level2, level2] = c.conj()
        
        return G
    
    @staticmethod
    def _reconstruct(rotations: List[Tuple[int, int, float, float]],
                    diagonal_phases: np.ndarray) -> np.ndarray:
        """
        回転列と対角位相からユニタリを再構築
        
        Args:
            rotations: [(level1, level2, theta, phi), ...]
            diagonal_phases: [φ₀, φ₁, φ₂]
            
        Returns:
            3×3ユニタリ行列
        """
        # 対角位相行列
        D = np.diag(np.exp(1j * diagonal_phases))
        
        # Givens回転を順に適用
        U = D.copy()
        for level1, level2, theta, phi in reversed(rotations):
            G = ImprovedThreeQuditDecomposer._construct_givens_matrix(
                3, level1, level2, theta, phi
            )
            U = G @ U
        
        return U
    
    @staticmethod
    def _is_unitary(U: np.ndarray, tolerance: float = 1e-10) -> bool:
        """行列がユニタリであるかを検証"""
        d = U.shape[0]
        identity = np.eye(d, dtype=complex)
        product = U.conj().T @ U
        return np.allclose(product, identity, atol=tolerance)
    
    @staticmethod
    def _compute_fidelity(U1: np.ndarray, U2: np.ndarray) -> float:
        """2つのユニタリ行列の忠実度を計算"""
        d = U1.shape[0]
        trace = np.trace(U1.conj().T @ U2)
        fidelity = abs(trace) / d
        return float(fidelity)


def test_improved_2x2_decomposition():
    """改良版2×2分解のテスト"""
    print("="*70)
    print("改良版2×2ユニタリ分解テスト")
    print("="*70)
    
    # ランダムなユニタリでテスト
    try:
        from scipy.stats import unitary_group
        num_tests = 100
    except ImportError:
        print("警告: scipy.stats.unitary_groupが利用できません")
        print("手動生成のユニタリでテスト")
        num_tests = 10
    
    decomposer = ImprovedTwoQubitDecomposer()
    fidelities = []
    failed_tests = []
    
    for i in range(num_tests):
        # ランダムな2×2ユニタリを生成
        if 'unitary_group' in locals():
            U = unitary_group.rvs(2)
        else:
            # 手動生成: ランダムなパラメータからユニタリを構築
            theta = np.random.uniform(0, np.pi)
            phi = np.random.uniform(-np.pi, np.pi)
            lam = np.random.uniform(-np.pi, np.pi)
            alpha = np.random.uniform(-np.pi, np.pi)
            U = decomposer._construct_zyz(theta, phi, lam, alpha)
        
        # 分解
        result = decomposer.decompose_zyz(U)
        fidelities.append(result.fidelity)
        
        if result.fidelity < 0.9999:
            failed_tests.append((i, result.fidelity))
    
    # 結果表示
    min_fidelity = min(fidelities)
    avg_fidelity = sum(fidelities) / len(fidelities)
    passed = sum(1 for f in fidelities if f >= 0.9999)
    
    print(f"\nテスト結果:")
    print(f"  - テスト数: {num_tests}")
    print(f"  - 最小忠実度: {min_fidelity:.10f}")
    print(f"  - 平均忠実度: {avg_fidelity:.10f}")
    print(f"  - 合格数: {passed}/{num_tests} ({passed/num_tests*100:.1f}%)")
    
    if failed_tests:
        print(f"\n失敗したテスト:")
        for idx, fid in failed_tests[:5]:  # 最初の5個のみ表示
            print(f"    テスト#{idx}: 忠実度 = {fid:.10f}")
    
    if min_fidelity >= 0.9999:
        print("\n✓ すべてのテストに合格（忠実度 > 0.9999）")
        return True
    else:
        print(f"\n✗ 一部のテストが不合格（最小忠実度 {min_fidelity:.10f} < 0.9999）")
        return False


def test_improved_3x3_decomposition():
    """改良版3×3分解のテスト"""
    print("\n" + "="*70)
    print("改良版3×3ユニタリ分解テスト")
    print("="*70)
    
    # ランダムなユニタリでテスト
    try:
        from scipy.stats import unitary_group
        num_tests = 100
    except ImportError:
        print("警告: scipy.stats.unitary_groupが利用できません")
        print("手動生成のユニタリでテスト")
        num_tests = 10
    
    decomposer = ImprovedThreeQuditDecomposer()
    fidelities = []
    failed_tests = []
    
    for i in range(num_tests):
        # ランダムな3×3ユニタリを生成
        if 'unitary_group' in locals():
            U = unitary_group.rvs(3)
        else:
            # 手動生成: ランダムなGivens回転から構築
            theta1 = np.random.uniform(0, np.pi)
            phi1 = np.random.uniform(-np.pi, np.pi)
            theta2 = np.random.uniform(0, np.pi)
            phi2 = np.random.uniform(-np.pi, np.pi)
            theta3 = np.random.uniform(0, np.pi)
            phi3 = np.random.uniform(-np.pi, np.pi)
            
            G1 = decomposer._construct_givens_matrix(3, 0, 1, theta1, phi1)
            G2 = decomposer._construct_givens_matrix(3, 0, 2, theta2, phi2)
            G3 = decomposer._construct_givens_matrix(3, 1, 2, theta3, phi3)
            
            phases = np.random.uniform(-np.pi, np.pi, 3)
            D = np.diag(np.exp(1j * phases))
            
            U = G1 @ G2 @ G3 @ D
        
        # 分解
        result = decomposer.decompose_givens(U)
        fidelities.append(result.fidelity)
        
        if result.fidelity < 0.9999:
            failed_tests.append((i, result.fidelity))
    
    # 結果表示
    min_fidelity = min(fidelities)
    avg_fidelity = sum(fidelities) / len(fidelities)
    passed = sum(1 for f in fidelities if f >= 0.9999)
    
    print(f"\nテスト結果:")
    print(f"  - テスト数: {num_tests}")
    print(f"  - 最小忠実度: {min_fidelity:.10f}")
    print(f"  - 平均忠実度: {avg_fidelity:.10f}")
    print(f"  - 合格数: {passed}/{num_tests} ({passed/num_tests*100:.1f}%)")
    
    if failed_tests:
        print(f"\n失敗したテスト:")
        for idx, fid in failed_tests[:5]:
            print(f"    テスト#{idx}: 忠実度 = {fid:.10f}")
    
    if min_fidelity >= 0.9999:
        print("\n✓ すべてのテストに合格（忠実度 > 0.9999）")
        return True
    else:
        print(f"\n✗ 一部のテストが不合格（最小忠実度 {min_fidelity:.10f} < 0.9999）")
        return False


if __name__ == '__main__':
    print("改良版厳密ユニタリ分解ツール - テスト実行")
    print("="*70)
    
    # 2×2分解テスト
    test_2x2_passed = test_improved_2x2_decomposition()
    
    # 3×3分解テスト
    test_3x3_passed = test_improved_3x3_decomposition()
    
    # 総合結果
    print("\n" + "="*70)
    print("総合テスト結果")
    print("="*70)
    print(f"2×2分解: {'✓ 合格' if test_2x2_passed else '✗ 不合格'}")
    print(f"3×3分解: {'✓ 合格' if test_3x3_passed else '✗ 不合格'}")
    
    if test_2x2_passed and test_3x3_passed:
        print("\n✓ すべてのテストに合格しました！")
        print("数学的厳密性が保証されています（忠実度 > 0.9999）")
    else:
        print("\n⚠ 一部のテストが不合格です")
        print("さらなる改良が必要です")
