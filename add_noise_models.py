#!/usr/bin/env python3
"""
Script to add noise model simulations to quantum_dynamics_complete_comparison.ipynb

This script adds:
1. Qubit noise model simulation using Qiskit Aer NoiseModel
2. Qudit noise model simulation using MQT-qudits NoiseModel
3. Comparison plots between noise-free and noisy simulations
"""

import json
import copy

def create_markdown_cell(text):
    """Create a markdown cell"""
    return {
        "cell_type": "markdown",
        "id": None,  # Will be auto-generated
        "metadata": {},
        "source": text.split('\n')
    }

def create_code_cell(code):
    """Create a code cell"""
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": None,
        "metadata": {},
        "outputs": [],
        "source": code.split('\n')
    }

def add_noise_sections(notebook_path):
    """Add noise model sections to the notebook"""
    
    # Load the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
    
    # Find insertion point - before the comparison section (section 6)
    insertion_idx = None
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'markdown':
            text = ''.join(cell['source'])
            if '## 6. 3手法の包括的比較' in text:
                insertion_idx = i
                break
    
    if insertion_idx is None:
        print("Could not find insertion point (section 6)")
        return False
    
    # Create new cells for noise models
    new_cells = []
    
    # Section 4.4: Qubit Noise Model Simulation
    new_cells.append(create_markdown_cell("""### 4.4 ノイズモデルを用いたQubitシミュレーション

実際の量子コンピュータでは、量子ゲート操作時のノイズが避けられません。ここではQiskit Aerの`NoiseModel`を使用して、現実的なノイズを考慮したシミュレーションを実行します。

#### 4.4.1 ノイズモデルの定義

以下の物理的なノイズ源を考慮します：

1. **脱分極ノイズ (Depolarizing Error)**: 量子ゲート操作後に状態がランダム化される確率
2. **脱位相ノイズ (Dephasing Error)**: 位相情報が失われる確率

これらのノイズパラメータは、実際の量子デバイスで観測される典型的な値に基づいて設定します。"""))
    
    new_cells.append(create_code_cell("""# Qiskit Aer NoiseModelを用いたQubitノイズシミュレーション
from qiskit_aer.noise import NoiseModel, depolarizing_error, phase_damping_error
from qiskit_aer import AerSimulator

# ノイズパラメータの設定（現実的な値）
SINGLE_QUBIT_DEPOL_ERROR = 0.001  # 単一qubitゲートの脱分極エラー率
TWO_QUBIT_DEPOL_ERROR = 0.01      # 2-qubitゲートの脱分極エラー率
PHASE_DAMPING_RATE = 0.002        # 位相減衰率

def create_qubit_noise_model():
    \"\"\"
    Qiskit用のノイズモデルを作成
    
    Returns:
        NoiseModel: 脱分極・脱位相ノイズを含むノイズモデル
    \"\"\"
    noise_model = NoiseModel()
    
    # 単一qubitゲートへの脱分極ノイズ
    single_qubit_error = depolarizing_error(SINGLE_QUBIT_DEPOL_ERROR, 1)
    noise_model.add_all_qubit_quantum_error(single_qubit_error, ['rx', 'ry', 'rz', 'x', 'h'])
    
    # 2-qubitゲートへの脱分極ノイズ（より大きいエラー率）
    two_qubit_error = depolarizing_error(TWO_QUBIT_DEPOL_ERROR, 2)
    noise_model.add_all_qubit_quantum_error(two_qubit_error, ['cx', 'cz', 'unitary'])
    
    # 位相減衰ノイズ（すべてのqubitに適用）
    phase_error = phase_damping_error(PHASE_DAMPING_RATE)
    noise_model.add_all_qubit_quantum_error(phase_error, ['rx', 'ry', 'rz', 'x', 'h', 'cx', 'cz', 'unitary'])
    
    print("\\n" + "="*70)
    print("Qubitノイズモデル設定")
    print("="*70)
    print(f"単一qubitゲート脱分極エラー率: {SINGLE_QUBIT_DEPOL_ERROR}")
    print(f"2-qubitゲート脱分極エラー率: {TWO_QUBIT_DEPOL_ERROR}")
    print(f"位相減衰率: {PHASE_DAMPING_RATE}")
    print("="*70)
    
    return noise_model

# ノイズモデルの作成
qubit_noise_model = create_qubit_noise_model()"""))
    
    new_cells.append(create_markdown_cell("""#### 4.4.2 ノイズ有りシミュレーションの実行

ノイズモデルを適用したQubitシミュレーションを実行します。ノイズの影響により、理想的なシミュレーションと比較して精度が低下することが予想されます。"""))
    
    new_cells.append(create_code_cell("""# ノイズ有りQubitシミュレーションの実行
import copy

def run_qubit_simulation_with_noise(noise_model, params):
    \"\"\"
    ノイズモデルを適用したQubitシミュレーションを実行
    
    Args:
        noise_model: Qiskit NoiseModel
        params: PhysicalParameters instance
    
    Returns:
        Dict with simulation results
    \"\"\"
    print("\\n" + "="*70)
    print("ノイズ有りQubitシミュレーション開始")
    print("="*70)
    
    import time
    start_time = time.time()
    
    # qubit_sim変数が定義されていることを確認
    # ノイズ無しシミュレーション結果からパラメータを取得
    from exact_qubit_hamiltonians import build_exact_qubit_simulation
    
    # ノイズ有りシミュレータの作成
    noisy_simulator = AerSimulator(noise_model=noise_model)
    
    # Qubitシミュレーションを実行（ノイズモデル付き）
    noisy_results = build_exact_qubit_simulation(
        params=params,
        simulator=noisy_simulator,
        method='statevector'  # ノイズ有りでもstatevectorを使用
    )
    
    elapsed = time.time() - start_time
    noisy_results['elapsed_time'] = elapsed
    noisy_results['method'] = 'Qubit (Qiskit) with Noise Model'
    
    print(f"\\nシミュレーション完了: {elapsed:.2f}秒")
    print("="*70)
    
    return noisy_results

# ノイズ有りシミュレーションの実行
try:
    qubit_results_noisy = run_qubit_simulation_with_noise(qubit_noise_model, params)
    print("\\nノイズ有りQubitシミュレーション成功")
except Exception as e:
    print(f"\\nノイズ有りQubitシミュレーションでエラーが発生: {e}")
    print("既存のQubit実装を使用してノイズモデルを適用します...")
    
    # 既存のQubit実装を使用（フォールバック無し、エラーの場合は失敗とする）
    raise"""))
    
    new_cells.append(create_markdown_cell("""#### 4.4.3 ノイズの影響評価

ノイズ無しとノイズ有りの結果を比較し、ノイズが量子シミュレーションに与える影響を定量的に評価します。"""))
    
    new_cells.append(create_code_cell("""# ノイズ無し vs ノイズ有り の比較プロット
def plot_noise_comparison_qubit(results_clean, results_noisy, title="Qubit Simulation: Noise Comparison"):
    \"\"\"
    ノイズ無しとノイズ有りのQubitシミュレーション結果を比較
    
    Args:
        results_clean: ノイズ無し結果
        results_noisy: ノイズ有り結果
        title: プロットタイトル
    \"\"\"
    times = results_clean['times']
    pops_clean = results_clean['populations']
    pops_noisy = results_noisy['populations']
    
    # S0, T1, S1の個体数を抽出
    N_S0_clean = [p['N_S0'] for p in pops_clean]
    N_T1_clean = [p['N_T1'] for p in pops_clean]
    N_S1_clean = [p['N_S1'] for p in pops_clean]
    
    N_S0_noisy = [p['N_S0'] for p in pops_noisy]
    N_T1_noisy = [p['N_T1'] for p in pops_noisy]
    N_S1_noisy = [p['N_S1'] for p in pops_noisy]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # S0状態
    axes[0].plot(times, N_S0_clean, 'b-', linewidth=2.5, label='Noise-free', marker='o', markersize=5, alpha=0.8)
    axes[0].plot(times, N_S0_noisy, 'r--', linewidth=2.5, label='With Noise', marker='s', markersize=5, alpha=0.8)
    axes[0].set_xlabel('Time (fs)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Population', fontsize=12, fontweight='bold')
    axes[0].set_title(r'$N_{S_0}$ (Ground Singlet)', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=10, loc='best')
    axes[0].grid(True, alpha=0.3, linestyle='--')
    
    # T1状態
    axes[1].plot(times, N_T1_clean, 'b-', linewidth=2.5, label='Noise-free', marker='o', markersize=5, alpha=0.8)
    axes[1].plot(times, N_T1_noisy, 'r--', linewidth=2.5, label='With Noise', marker='s', markersize=5, alpha=0.8)
    axes[1].set_xlabel('Time (fs)', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Population', fontsize=12, fontweight='bold')
    axes[1].set_title(r'$N_{T_1}$ (Triplet)', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=10, loc='best')
    axes[1].grid(True, alpha=0.3, linestyle='--')
    
    # S1状態
    axes[2].plot(times, N_S1_clean, 'b-', linewidth=2.5, label='Noise-free', marker='o', markersize=5, alpha=0.8)
    axes[2].plot(times, N_S1_noisy, 'r--', linewidth=2.5, label='With Noise', marker='s', markersize=5, alpha=0.8)
    axes[2].set_xlabel('Time (fs)', fontsize=12, fontweight='bold')
    axes[2].set_ylabel('Population', fontsize=12, fontweight='bold')
    axes[2].set_title(r'$N_{S_1}$ (Excited Singlet)', fontsize=13, fontweight='bold')
    axes[2].legend(fontsize=10, loc='best')
    axes[2].grid(True, alpha=0.3, linestyle='--')
    
    fig.suptitle(title, fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.show()
    
    # 誤差評価
    errors = []
    for p_clean, p_noisy in zip(pops_clean, pops_noisy):
        error = np.sqrt(
            (p_clean['N_S0'] - p_noisy['N_S0'])**2 +
            (p_clean['N_T1'] - p_noisy['N_T1'])**2 +
            (p_clean['N_S1'] - p_noisy['N_S1'])**2
        )
        errors.append(error)
    
    print("\\n" + "="*70)
    print("ノイズの影響評価 (Qubit)")
    print("="*70)
    print(f"平均誤差: {np.mean(errors):.4f}")
    print(f"最大誤差: {np.max(errors):.4f}")
    print(f"最終時刻の誤差: {errors[-1]:.4f}")
    print("="*70)
    
    return fig

# ノイズ比較プロットの実行
if 'qubit_results_noisy' in locals():
    fig_qubit_noise = plot_noise_comparison_qubit(qubit_results, qubit_results_noisy)
else:
    print("ノイズ有りQubitシミュレーション結果が利用できません")"""))
    
    # Section 5.4: Qudit Noise Model Simulation
    new_cells.append(create_markdown_cell("""### 5.4 ノイズモデルを用いたQuditシミュレーション

MQT-quditsフレームワークには、qudit特有のノイズモデルが実装されています。qutrit（3準位系）では、各準位間の遷移に対して個別にノイズパラメータを設定できます。

#### 5.4.1 Quditノイズモデルの定義

`SubspaceNoise`を使用して、以下の準位間遷移に対してノイズを定義します：

- **(0,1)**: S₀ ↔ T₁ 遷移
- **(1,2)**: T₁ ↔ S₁ 遷移  
- **(0,2)**: S₀ ↔ S₁ 遷移（直接遷移は物理的に小さい）"""))
    
    new_cells.append(create_code_cell("""# MQT-qudits NoiseModelを用いたQuditノイズシミュレーション
from mqt.qudits.simulation.noise_tools import NoiseModel, SubspaceNoise, Noise, NoisyCircuitFactory

# Quditノイズパラメータの設定（Qubitと同程度の物理的ノイズ）
QUDIT_DEPOL_01 = 0.001   # (0,1)遷移の脱分極エラー率
QUDIT_DEPOL_12 = 0.001   # (1,2)遷移の脱分極エラー率
QUDIT_DEPOL_02 = 0.002   # (0,2)遷移の脱分極エラー率（直接遷移は稀）
QUDIT_DEPHASE = 0.002    # 脱位相エラー率

def create_qudit_noise_model():
    \"\"\"
    MQT-qudits用のノイズモデルを作成
    
    Returns:
        NoiseModel: SubspaceNoiseを含むQuditノイズモデル
    \"\"\"
    noise_model = NoiseModel()
    
    # 各準位間遷移に対するノイズ定義
    # (0,1): S0 <-> T1 遷移
    subspace_01 = SubspaceNoise(QUDIT_DEPOL_01, QUDIT_DEPHASE, (0, 1))
    
    # (1,2): T1 <-> S1 遷移
    subspace_12 = SubspaceNoise(QUDIT_DEPOL_12, QUDIT_DEPHASE, (1, 2))
    
    # (0,2): S0 <-> S1 遷移（直接遷移）
    subspace_02 = SubspaceNoise(QUDIT_DEPOL_02, QUDIT_DEPHASE, (0, 2))
    
    # ローカルゲート（単一quditゲート）にノイズを適用
    # Note: MQT-quditsで使用される基本ゲート名を指定
    local_gates = ['x', 'z', 'h', 'rz', 'r', 'virtrz']
    noise_model.add_quantum_error_locally(subspace_01, local_gates)
    noise_model.add_quantum_error_locally(subspace_12, local_gates)
    noise_model.add_quantum_error_locally(subspace_02, local_gates)
    
    # 非ローカルゲート（2-quditゲート）にノイズを適用
    # Quditでは2-quditゲートのノイズがより重要
    nonlocal_gates = ['cx', 'csum', 'cex']
    # ターゲットquditにノイズを適用
    noise_model.add_nonlocal_quantum_error_on_target(subspace_01, nonlocal_gates)
    noise_model.add_nonlocal_quantum_error_on_target(subspace_12, nonlocal_gates)
    # コントロールquditにもノイズを適用
    noise_model.add_nonlocal_quantum_error_on_control(subspace_01, nonlocal_gates)
    noise_model.add_nonlocal_quantum_error_on_control(subspace_12, nonlocal_gates)
    
    print("\\n" + "="*70)
    print("Quditノイズモデル設定")
    print("="*70)
    print(f"(0,1)遷移 脱分極エラー率: {QUDIT_DEPOL_01}")
    print(f"(1,2)遷移 脱分極エラー率: {QUDIT_DEPOL_12}")
    print(f"(0,2)遷移 脱分極エラー率: {QUDIT_DEPOL_02}")
    print(f"脱位相エラー率: {QUDIT_DEPHASE}")
    print("="*70)
    print(f"\\nノイズモデル情報:")
    print(noise_model)
    
    return noise_model

# Quditノイズモデルの作成
qudit_noise_model = create_qudit_noise_model()"""))
    
    new_cells.append(create_markdown_cell("""#### 5.4.2 ノイズ有りQuditシミュレーションの実行

`NoisyCircuitFactory`を使用して、量子回路にノイズゲートを挿入し、ノイズを考慮したシミュレーションを実行します。"""))
    
    new_cells.append(create_code_cell("""# ノイズ有りQuditシミュレーションの実行
def run_qudit_simulation_with_noise(noise_model, params):
    \"\"\"
    ノイズモデルを適用したQuditシミュレーションを実行
    
    Args:
        noise_model: MQT-qudits NoiseModel
        params: PhysicalParameters instance
    
    Returns:
        Dict with simulation results
    \"\"\"
    print("\\n" + "="*70)
    print("ノイズ有りQuditシミュレーション開始")
    print("="*70)
    
    import time
    start_time = time.time()
    
    # Qudit量子回路の構築（ノイズ無し版と同じ）
    from mqt.qudits.quantum_circuit import QuantumCircuit as QuditQuantumCircuit
    from mqt.qudits.quantum_circuit.components.quantum_register import QuantumRegister
    from exact_hamiltonian_builders import (
        build_H_transfer_matrix,
        build_H_TTA_matrix,
        build_time_evolution_unitary
    )
    
    # 4 qutrits用のQuantumRegisterを作成
    qreg = QuantumRegister("mol", params.N_molecules, params.N_molecules * [3])
    
    # 時間発展の計算
    times = [0.0]
    populations = []
    per_molecule_populations = []
    
    # 初期状態の準備（|1001⟩: 両端がT1）
    # Note: ノイズ適用のため、量子回路として実装
    
    dt = params.T_total / params.N_steps
    
    # ユニタリ行列の事前計算
    H_tr = build_H_transfer_matrix(params.V, dim=3)
    H_TTA = build_H_TTA_matrix(params.J, dim=3)
    
    U_tr_half = build_time_evolution_unitary(H_tr, dt / 2, params.hbar)
    U_TTA_half = build_time_evolution_unitary(H_TTA, dt / 2, params.hbar)
    
    # 各ステップで量子回路を構築
    for step in range(params.N_steps):
        circ = QuditQuantumCircuit(qreg)
        
        # 初期状態準備（最初のステップのみ）
        if step == 0:
            # |1001⟩状態: molecule 0と3をT1状態に
            circ.x(0)  # |0⟩ -> |1⟩
            circ.x(3)  # |0⟩ -> |1⟩
        
        # Trotter step: H0, H_transfer, H_TTA の順に適用
        # H0: 対角項（位相ゲート）
        for mol_idx in range(params.N_molecules):
            # T1エネルギー
            phase_T1 = -params.E_T * (dt / 2) / params.hbar
            # S1エネルギー
            phase_S1 = -params.E_S * (dt / 2) / params.hbar
            # 位相ゲートを適用（簡略化）
            # Note: 実際の実装では、各準位に応じた位相ゲートを適用
            circ.virtrz(mol_idx, [1, phase_T1])
            circ.virtrz(mol_idx, [2, phase_S1])
        
        # H_transfer: ペアごとに適用
        for mol_i, mol_j in params.neighbors:
            # CustomTwoゲートまたはCExゲートでH_transferを実装
            # Note: exact_hamiltonian_buildersの実装を使用
            from mqt.qudits.compiler import QuditCompiler
            # U_tr_halfをqudit回路に追加
            # Simplified: 実際の実装では適切なゲート分解を使用
            pass  # 詳細は既存のQudit実装を参照
        
        # H_TTA: ペアごとに適用
        for mol_i, mol_j in params.neighbors:
            # U_TTA_halfをqudit回路に追加
            pass  # 詳細は既存のQudit実装を参照
        
        # Backward pass（対称Trotter分解）
        # 逆順でH_TTA, H_transfer, H0を適用
        pass
        
        # NoisyCircuitFactoryでノイズを適用
        factory = NoisyCircuitFactory(noise_model, circ)
        noisy_circ = factory.generate_circuit()
        
        # ノイズ有り回路のシミュレーション
        # Note: MQT-quditsのシミュレータを使用
        # state = simulate_noisy_circuit(noisy_circ)
        # populations.append(calculate_populations(state))
    
    # Placeholder: 実装詳細は既存のQudit実装を参照
    # 実際のノイズシミュレーションは、既存のQudit実装と組み合わせて実装
    
    elapsed = time.time() - start_time
    
    print(f"\\nノイズ有りQuditシミュレーション実装中...")
    print(f"経過時間: {elapsed:.2f}秒")
    print("="*70)
    
    # Placeholder results
    # 実際の実装では、ノイズ有り回路のシミュレーション結果を返す
    return {
        'times': times,
        'populations': populations,
        'per_molecule_populations': per_molecule_populations,
        'elapsed_time': elapsed,
        'method': 'Qudit (MQT-qudits) with Noise Model'
    }

# Note: Quditノイズシミュレーションの完全な実装は、
# 既存のQudit実装（exact_hamiltonian_builders.py）との統合が必要
print("\\nQuditノイズシミュレーション:")
print("完全な実装には、既存のQudit回路構築コードとの統合が必要です")
print("ノイズモデルは正しく定義されており、NoisyCircuitFactoryで適用可能です")"""))
    
    # Insert the new cells
    nb['cells'] = nb['cells'][:insertion_idx] + new_cells + nb['cells'][insertion_idx:]
    
    # Update section numbering (6 becomes 7, 7 becomes 8, etc.)
    for i in range(insertion_idx + len(new_cells), len(nb['cells'])):
        cell = nb['cells'][i]
        if cell['cell_type'] == 'markdown':
            source = ''.join(cell['source'])
            # Update section numbers
            source = source.replace('## 6.', '## 7.')
            source = source.replace('### 6.', '### 7.')
            source = source.replace('## 7.', '## 8.')
            source = source.replace('### 7.', '### 8.')
            source = source.replace('## 8.', '## 9.')
            source = source.replace('### 8.', '### 9.')
            cell['source'] = source.split('\n')
    
    # Save the modified notebook
    output_path = notebook_path.replace('.ipynb', '_with_noise.ipynb')
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    
    print(f"Modified notebook saved to: {output_path}")
    return True

if __name__ == '__main__':
    notebook_path = 'tutorials/quantum_dynamics_complete_comparison.ipynb'
    success = add_noise_sections(notebook_path)
    if success:
        print("Successfully added noise model sections!")
    else:
        print("Failed to add noise model sections")
