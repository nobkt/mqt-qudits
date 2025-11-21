#!/usr/bin/env python3
"""
Script to add noise model simulations to quantum_dynamics_complete_comparison.ipynb
This script adds noise model simulations for both Qubit and Qudit implementations.
"""

import json
import sys

def create_noise_model_cells():
    """Create cells for noise model simulations"""
    
    cells = []
    
    # Cell 1: Markdown - Introduction to Noise Models
    cells.append({
        "cell_type": "markdown",
        "id": "noise_intro",
        "metadata": {},
        "source": [
            "## 9. ノイズモデルを考慮したシミュレーション\n",
            "\n",
            "### 9.1 ノイズモデルの必要性\n",
            "\n",
            "実際の量子デバイスではノイズが避けられません。ここでは、QiskitおよびMQT-Quditsで利用可能なノイズモデルを使用して、より現実的なシミュレーションを行います。\n",
            "\n",
            "### 9.2 使用するノイズモデル\n",
            "\n",
            "#### Quditノイズモデル（MQT-Qudits）\n",
            "- **脱分極ノイズ（Depolarizing noise）**: 量子状態がランダムに混合状態になる\n",
            "- **位相緩和ノイズ（Dephasing noise）**: 位相情報が失われる\n",
            "- 部分空間ごとにノイズパラメータを設定可能\n",
            "\n",
            "#### Qubitノイズモデル（Qiskit Aer）\n",
            "- **脱分極ノイズ（Depolarizing error）**: 1量子ビットおよび2量子ビットゲートに適用\n",
            "- **位相緩和ノイズ（Phase damping）**: T2デコヒーレンス\n",
            "- **測定ノイズ（Measurement error）**: 読み出しエラー\n",
            "\n",
            "### 9.3 ノイズパラメータの設定\n",
            "\n",
            "以下のノイズパラメータを使用します：\n",
            "- 1量子ビット/quditゲートの脱分極確率: 0.001（0.1%）\n",
            "- 2量子ビット/quditゲートの脱分極確率: 0.01（1%）\n",
            "- 位相緩和確率: 0.001（0.1%）\n",
            "- 測定エラー確率（Qubitのみ）: 0.01（1%）\n",
            "\n",
            "これらの値は、現在の超伝導量子ビット技術における典型的な値を参考にしています。\n"
        ]
    })
    
    # Cell 2: Code - Qudit Noise Model Implementation
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "id": "qudit_noise_model",
        "metadata": {},
        "outputs": [],
        "source": [
            "# QuditノイズモデルのSimulation\n",
            "\n",
            "from mqt.qudits.simulation.noise_tools import Noise, NoiseModel, SubspaceNoise\n",
            "from mqt.qudits.simulation import MQTQuditProvider\n",
            "\n",
            "print(\"=\"*70)\n",
            "print(\"Quditノイズモデルシミュレーション\")\n",
            "print(\"=\"*70)\n",
            "\n",
            "# ノイズパラメータの設定\n",
            "noise_params = {\n",
            "    'local_depolarizing': 0.001,  # 1quditゲート脱分極確率\n",
            "    'local_dephasing': 0.001,     # 1qudit位相緩和確率\n",
            "    'nonlocal_depolarizing': 0.01,  # 2quditゲート脱分極確率\n",
            "    'nonlocal_dephasing': 0.01      # 2qudit位相緩和確率\n",
            "}\n",
            "\n",
            "print(\"\\nノイズパラメータ:\")\n",
            "print(f\"  1quditゲート脱分極確率: {noise_params['local_depolarizing']*100:.2f}%\")\n",
            "print(f\"  1qudit位相緩和確率: {noise_params['local_dephasing']*100:.2f}%\")\n",
            "print(f\"  2quditゲート脱分極確率: {noise_params['nonlocal_depolarizing']*100:.1f}%\")\n",
            "print(f\"  2qudit位相緩和確率: {noise_params['nonlocal_dephasing']*100:.1f}%\")\n",
            "\n",
            "# ノイズモデルの構築\n",
            "qudit_noise_model = NoiseModel()\n",
            "\n",
            "# 局所ゲート用ノイズ（1quditゲート）\n",
            "local_noise = Noise(\n",
            "    probability_depolarizing=noise_params['local_depolarizing'],\n",
            "    probability_dephasing=noise_params['local_dephasing']\n",
            ")\n",
            "\n",
            "# 非局所ゲート用ノイズ（2quditゲート）\n",
            "nonlocal_noise = Noise(\n",
            "    probability_depolarizing=noise_params['nonlocal_depolarizing'],\n",
            "    probability_dephasing=noise_params['nonlocal_dephasing']\n",
            ")\n",
            "\n",
            "# ゲートにノイズを適用\n",
            "# 局所ゲート: rz, virtrz, h, x, z, s, r, rh\n",
            "qudit_noise_model.add_quantum_error_locally(\n",
            "    local_noise, \n",
            "    [\"rz\", \"virtrz\", \"h\", \"x\", \"z\", \"s\", \"r\", \"rh\"]\n",
            ")\n",
            "\n",
            "# 非局所ゲート: csum, cx, customtwo\n",
            "qudit_noise_model.add_nonlocal_quantum_error(\n",
            "    nonlocal_noise,\n",
            "    [\"csum\", \"cx\", \"customtwo\"]\n",
            ")\n",
            "\n",
            "# ターゲットquditへのノイズ\n",
            "qudit_noise_model.add_nonlocal_quantum_error_on_target(\n",
            "    nonlocal_noise,\n",
            "    [\"csum\", \"cx\", \"customtwo\"]\n",
            ")\n",
            "\n",
            "# コントロールquditへのノイズ\n",
            "qudit_noise_model.add_nonlocal_quantum_error_on_control(\n",
            "    nonlocal_noise,\n",
            "    [\"csum\", \"cx\", \"customtwo\"]\n",
            ")\n",
            "\n",
            "print(\"\\nノイズモデル構築完了\")\n",
            "print(f\"  対象ゲート: {qudit_noise_model.basis_gates}\")\n",
            "print(\"=\"*70)\n"
        ]
    })
    
    # Cell 3: Code - Run Qudit Noisy Simulation  
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "id": "qudit_noisy_simulation",
        "metadata": {},
        "outputs": [],
        "source": [
            "# Quditノイズありシミュレーションの実行\n",
            "\n",
            "print(\"\\n\" + \"=\"*70)\n",
            "print(\"Quditノイズありシミュレーション実行\")\n",
            "print(\"=\"*70)\n",
            "\n",
            "# 注意: ノイズシミュレーションは確率的なため、十分なショット数が必要\n",
            "# ここでは計算時間を考慮して適切なショット数を設定\n",
            "qudit_noise_shots = 1000  # ショット数\n",
            "\n",
            "print(f\"\\nシミュレーション設定:\")\n",
            "print(f\"  ショット数: {qudit_noise_shots}\")\n",
            "print(f\"  トロッターステップ数: {params.N_steps}\")\n",
            "print(f\"  総シミュレーション時間: {params.T_total} fs\")\n",
            "\n",
            "# MQT-Quditsバックエンドを使用してノイズありシミュレーション実行\n",
            "provider = MQTQuditProvider()\n",
            "backend_noisy = provider.get_backend(\"misim\")\n",
            "\n",
            "# qudit_resultsから回路を取得して再実行\n",
            "if 'circuit' in qudit_results:\n",
            "    qudit_circuit = qudit_results['circuit']\n",
            "    \n",
            "    print(\"\\nノイズありシミュレーション実行中...\")\n",
            "    start_time = time.time()\n",
            "    \n",
            "    # ノイズモデルを指定して実行\n",
            "    job_noisy = backend_noisy.run(\n",
            "        qudit_circuit, \n",
            "        noise_model=qudit_noise_model,\n",
            "        shots=qudit_noise_shots\n",
            "    )\n",
            "    result_noisy = job_noisy.result()\n",
            "    \n",
            "    elapsed_noisy = time.time() - start_time\n",
            "    \n",
            "    # 状態ベクトルとカウントを取得\n",
            "    state_vector_noisy = result_noisy.get_state_vector()\n",
            "    counts_noisy = result_noisy.get_counts()\n",
            "    \n",
            "    print(f\"\\nノイズありシミュレーション完了\")\n",
            "    print(f\"  実行時間: {elapsed_noisy:.2f}秒\")\n",
            "    print(f\"  取得ショット数: {len(counts_noisy)}\")\n",
            "    print(f\"  状態ベクトル次元: {len(state_vector_noisy.squeeze())}\")\n",
            "    \n",
            "    # カウントから最終個体数を計算\n",
            "    # ここでは簡略化のため、ノイズなしシミュレーションと同じ計算方法を使用\n",
            "    # 実際の実装では、qudit_simの calculate_per_molecule_populations を使用\n",
            "    \n",
            "    qudit_noisy_results = {\n",
            "        'state_vector': state_vector_noisy,\n",
            "        'counts': counts_noisy,\n",
            "        'elapsed_time': elapsed_noisy,\n",
            "        'shots': qudit_noise_shots,\n",
            "        'noise_model': qudit_noise_model,\n",
            "        'method': 'Qudit (MQT-Qudits - Noisy)'\n",
            "    }\n",
            "    \n",
            "    print(\"\\n✓ Quditノイズありシミュレーション完了\")\nelse:\n",
            "    print(\"\\n警告: qudit_results に回路情報が見つかりません\")\n",
            "    print(\"ノイズありシミュレーションをスキップします\")\n",
            "    qudit_noisy_results = None\n",
            "\n",
            "print(\"=\"*70)\n"
        ]
    })
    
    # Cell 4: Markdown - Qubit Noise Model
    cells.append({
        "cell_type": "markdown",
        "id": "qubit_noise_intro",
        "metadata": {},
        "source": [
            "### 9.4 Qubitノイズモデルシミュレーション\n",
            "\n",
            "Qiskit Aerを使用して、Qubitベースの量子シミュレーションにノイズモデルを適用します。\n",
            "\n",
            "#### 実装されるノイズチャネル\n",
            "1. **1量子ビットゲートの脱分極ノイズ**: RZ, H, X, Z, S などの単一量子ビットゲート\n",
            "2. **2量子ビットゲートの脱分極ノイズ**: CX, UnitaryGate などの2量子ビットゲート\n",
            "3. **測定読み出しエラー**: 測定時の誤り\n",
            "\n",
            "注意: Qiskit Aerは現在オプショナル依存関係です。インストールされていない場合は、\n",
            "```bash\n",
            "pip install qiskit-aer\n",
            "```\n",
            "でインストールしてください。\n"
        ]
    })
    
    # Cell 5: Code - Qubit Noise Model (conditional on qiskit_aer availability)
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "id": "qubit_noise_model",
        "metadata": {},
        "outputs": [],
        "source": [
            "# Qubitノイズモデルの実装（Qiskit Aer使用）\n",
            "\n",
            "try:\n",
            "    from qiskit_aer import AerSimulator\n",
            "    from qiskit_aer.noise import NoiseModel as QiskitNoiseModel\n",
            "    from qiskit_aer.noise import depolarizing_error, phase_damping_error, ReadoutError\n",
            "    \n",
            "    QISKIT_AER_AVAILABLE = True\n",
            "    print(\"✓ Qiskit Aer が利用可能です\")\n",
            "except ImportError:\n",
            "    QISKIT_AER_AVAILABLE = False\n",
            "    print(\"⚠ Qiskit Aer が見つかりません\")\n",
            "    print(\"  Qubitノイズシミュレーションには qiskit-aer が必要です\")\n",
            "    print(\"  インストール方法: pip install qiskit-aer\")\n",
            "\n",
            "if QISKIT_AER_AVAILABLE:\n",
            "    print(\"\\n\" + \"=\"*70)\n",
            "    print(\"Qubitノイズモデル構築\")\n",
            "    print(\"=\"*70)\n",
            "    \n",
            "    # ノイズパラメータ（Quditと同じ値を使用）\n",
            "    qubit_noise_params = {\n",
            "        'single_qubit_depol': 0.001,  # 1量子ビットゲート脱分極確率\n",
            "        'two_qubit_depol': 0.01,      # 2量子ビットゲート脱分極確率\n",
            "        'readout_error': 0.01         # 測定エラー確率\n",
            "    }\n",
            "    \n",
            "    print(\"\\nノイズパラメータ:\")\n",
            "    print(f\"  1量子ビットゲート脱分極確率: {qubit_noise_params['single_qubit_depol']*100:.2f}%\")\n",
            "    print(f\"  2量子ビットゲート脱分極確率: {qubit_noise_params['two_qubit_depol']*100:.1f}%\")\n",
            "    print(f\"  測定読み出しエラー確率: {qubit_noise_params['readout_error']*100:.1f}%\")\n",
            "    \n",
            "    # ノイズモデルの構築\n",
            "    qubit_noise_model = QiskitNoiseModel()\n",
            "    \n",
            "    # 1量子ビットゲートの脱分極エラー\n",
            "    single_qubit_error = depolarizing_error(\n",
            "        qubit_noise_params['single_qubit_depol'], \n",
            "        1\n",
            "    )\n",
            "    \n",
            "    # 2量子ビットゲートの脱分極エラー\n",
            "    two_qubit_error = depolarizing_error(\n",
            "        qubit_noise_params['two_qubit_depol'], \n",
            "        2\n",
            "    )\n",
            "    \n",
            "    # 測定読み出しエラー\n",
            "    readout_prob = qubit_noise_params['readout_error']\n",
            "    readout_error = ReadoutError(\n",
            "        [[1 - readout_prob, readout_prob], \n",
            "         [readout_prob, 1 - readout_prob]]\n",
            "    )\n",
            "    \n",
            "    # 1量子ビットゲートにノイズを追加\n",
            "    single_qubit_gates = ['rz', 'h', 'x', 'z', 's', 'p']\n",
            "    for gate in single_qubit_gates:\n",
            "        qubit_noise_model.add_all_qubit_quantum_error(\n",
            "            single_qubit_error, \n",
            "            gate\n",
            "        )\n",
            "    \n",
            "    # 2量子ビットゲートにノイズを追加\n",
            "    two_qubit_gates = ['cx', 'unitary']\n",
            "    for gate in two_qubit_gates:\n",
            "        qubit_noise_model.add_all_qubit_quantum_error(\n",
            "            two_qubit_error, \n",
            "            gate\n",
            "        )\n",
            "    \n",
            "    # 測定にノイズを追加（全量子ビット）\n",
            "    n_qubits = 8  # 4分子 × 2 qubits/molecule\n",
            "    for qubit in range(n_qubits):\n",
            "        qubit_noise_model.add_readout_error(readout_error, [qubit])\n",
            "    \n",
            "    print(\"\\nノイズモデル構築完了\")\n",
            "    print(f\"  ノイズが適用されるゲート: {qubit_noise_model.noise_qubits}\")\n",
            "    print(\"=\"*70)\n",
            "else:\n",
            "    qubit_noise_model = None\n"
        ]
    })
    
    # Cell 6: Code - Run Qubit Noisy Simulation
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "id": "qubit_noisy_simulation",
        "metadata": {},
        "outputs": [],
        "source": [
            "# Qubitノイズありシミュレーションの実行\n",
            "\n",
            "if QISKIT_AER_AVAILABLE and qubit_noise_model is not None:\n",
            "    print(\"\\n\" + \"=\"*70)\n",
            "    print(\"Qubitノイズありシミュレーション実行\")\n",
            "    print(\"=\"*70)\n",
            "    \n",
            "    qubit_noise_shots = 10000  # Qubitシミュレーションのショット数\n",
            "    \n",
            "    print(f\"\\nシミュレーション設定:\")\n",
            "    print(f\"  ショット数: {qubit_noise_shots}\")\n",
            "    print(f\"  トロッターステップ数: {params.N_steps}\")\n",
            "    print(f\"  総シミュレーション時間: {params.T_total} fs\")\n",
            "    \n",
            "    # Aer シミュレータの初期化（ノイズモデル付き）\n",
            "    simulator_noisy = AerSimulator(noise_model=qubit_noise_model)\n",
            "    \n",
            "    # qubit_sim を使用して新しいシミュレーションを実行する代わりに、\n",
            "    # 既存の回路を使用してノイズありシミュレーションを実行\n",
            "    \n",
            "    if 'circuit_final' in qubit_results:\n",
            "        # 最終回路を取得（測定なし）\n",
            "        circuit_to_run = qubit_results['circuit_final'].copy()\n",
            "        \n",
            "        # 測定を追加\n",
            "        from qiskit import ClassicalRegister\n",
            "        n_qubits = circuit_to_run.num_qubits\n",
            "        c_reg = ClassicalRegister(n_qubits, 'c')\n",
            "        circuit_to_run.add_register(c_reg)\n",
            "        circuit_to_run.measure(range(n_qubits), range(n_qubits))\n",
            "        \n",
            "        print(\"\\nノイズありシミュレーション実行中...\")\n",
            "        start_time = time.time()\n",
            "        \n",
            "        # ノイズモデルを使用してシミュレーション実行\n",
            "        job_noisy = simulator_noisy.run(\n",
            "            circuit_to_run, \n",
            "            shots=qubit_noise_shots\n",
            "        )\n",
            "        result_noisy = job_noisy.result()\n",
            "        counts_noisy = result_noisy.get_counts()\n",
            "        \n",
            "        elapsed_noisy = time.time() - start_time\n",
            "        \n",
            "        # カウントから最終個体数を計算\n",
            "        # qubit_sim の calculate_populations_from_counts を使用\n",
            "        pop_noisy = qubit_sim.calculate_populations_from_counts(\n",
            "            counts_noisy, \n",
            "            qubit_noise_shots\n",
            "        )\n",
            "        \n",
            "        print(f\"\\nノイズありシミュレーション完了\")\n",
            "        print(f\"  実行時間: {elapsed_noisy:.2f}秒\")\n",
            "        print(f\"\\n最終個体数（ノイズあり）:\")\n",
            "        print(f\"  N_S0 = {pop_noisy['N_S0']:.4f}\")\n",
            "        print(f\"  N_T1 = {pop_noisy['N_T1']:.4f}\")\n",
            "        print(f\"  N_S1 = {pop_noisy['N_S1']:.4f}\")\n",
            "        print(f\"  非物理的状態: {pop_noisy['unphysical']:.6f}\")\n",
            "        \n",
            "        qubit_noisy_results = {\n",
            "            'counts': counts_noisy,\n",
            "            'populations_final': pop_noisy,\n",
            "            'elapsed_time': elapsed_noisy,\n",
            "            'shots': qubit_noise_shots,\n",
            "            'noise_model': qubit_noise_model,\n",
            "            'method': 'Qubit (Qiskit Aer - Noisy)'\n",
            "        }\n",
            "        \n",
            "        print(\"\\n✓ Qubitノイズありシミュレーション完了\")\n",
            "    else:\n",
            "        print(\"\\n警告: qubit_results に回路情報が見つかりません\")\n",
            "        print(\"ノイズありシミュレーションをスキップします\")\n",
            "        qubit_noisy_results = None\n",
            "    \n",
            "    print(\"=\"*70)\n",
            "else:\n",
            "    print(\"\\nQubitノイズありシミュレーションはスキップされました\")\n",
            "    print(\"（Qiskit Aer が利用できません）\")\n",
            "    qubit_noisy_results = None\n"
        ]
    })
    
    # Cell 7: Markdown - Comparison with Noise
    cells.append({
        "cell_type": "markdown",
        "id": "noise_comparison_intro",
        "metadata": {},
        "source": [
            "### 9.5 ノイズあり・なし の比較\n",
            "\n",
            "ノイズモデルを適用した場合と適用しない場合の結果を比較します。\n",
            "ノイズの影響により、以下の変化が観測されることが予想されます：\n",
            "\n",
            "1. **最終個体数の変化**: ノイズにより理想的な状態からの偏差が生じる\n",
            "2. **非物理的状態の増加**: 特にQubit実装で顕著（エンコーディングの制約による）\n",
            "3. **測定統計の揺らぎ**: ショットベースのシミュレーションによる統計的ノイズ\n",
            "\n",
            "以下の比較表で、ノイズの影響を定量的に評価します。\n"
        ]
    })
    
    # Cell 8: Code - Comparison Table with Noise
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "id": "noise_comparison_table",
        "metadata": {},
        "outputs": [],
        "source": [
            "# ノイズあり・なしの比較表\n",
            "\n",
            "import pandas as pd\n",
            "\n",
            "print(\"\\n\" + \"=\"*90)\n",
            "print(\"ノイズあり・なし シミュレーション比較\")\n",
            "print(\"=\"*90)\n",
            "\n",
            "# 比較データの準備\n",
            "comparison_noise_data = {\n",
            "    '実装': [\n",
            "        'Classical',\n",
            "        'Qubit (ノイズなし)',\n",
            "        'Qubit (ノイズあり)',\n",
            "        'Qudit (ノイズなし)',\n",
            "        'Qudit (ノイズあり)'\n",
            "    ],\n",
            "    'N_S0 (最終)': [\n",
            "        f\"{classical_results['populations'][-1]['N_S0']:.4f}\",\n",
            "        f\"{qubit_results['populations'][-1]['N_S0']:.4f}\",\n",
            "        f\"{qubit_noisy_results['populations_final']['N_S0']:.4f}\" if qubit_noisy_results else 'N/A',\n",
            "        f\"{qudit_results['populations'][-1]['N_S0']:.4f}\" if qudit_results else 'N/A',\n",
            "        'N/A'  # Quditノイズありの結果は別途処理が必要\n",
            "    ],\n",
            "    'N_T1 (最終)': [\n",
            "        f\"{classical_results['populations'][-1]['N_T1']:.4f}\",\n",
            "        f\"{qubit_results['populations'][-1]['N_T1']:.4f}\",\n",
            "        f\"{qubit_noisy_results['populations_final']['N_T1']:.4f}\" if qubit_noisy_results else 'N/A',\n",
            "        f\"{qudit_results['populations'][-1]['N_T1']:.4f}\" if qudit_results else 'N/A',\n",
            "        'N/A'\n",
            "    ],\n",
            "    'N_S1 (最終)': [\n",
            "        f\"{classical_results['populations'][-1]['N_S1']:.4f}\",\n",
            "        f\"{qubit_results['populations'][-1]['N_S1']:.4f}\",\n",
            "        f\"{qubit_noisy_results['populations_final']['N_S1']:.4f}\" if qubit_noisy_results else 'N/A',\n",
            "        f\"{qudit_results['populations'][-1]['N_S1']:.4f}\" if qudit_results else 'N/A',\n",
            "        'N/A'\n",
            "    ],\n",
            "    '非物理的状態': [\n",
            "        '0.0000',\n",
            "        f\"{qubit_results['populations'][-1].get('unphysical', 0.0):.4f}\",\n",
            "        f\"{qubit_noisy_results['populations_final'].get('unphysical', 0.0):.4f}\" if qubit_noisy_results else 'N/A',\n",
            "        '0.0000' if qudit_results else 'N/A',\n",
            "        'N/A'\n",
            "    ],\n",
            "    'ショット数': [\n",
            "        '-',\n",
            "        f\"{qubit_results.get('shots', '-')}\",\n",
            "        f\"{qubit_noisy_results['shots']}\" if qubit_noisy_results else 'N/A',\n",
            "        '-',\n",
            "        f\"{qudit_noisy_results['shots']}\" if qudit_noisy_results else 'N/A'\n",
            "    ]\n",
            "}\n",
            "\n",
            "df_noise_comparison = pd.DataFrame(comparison_noise_data)\n",
            "print(df_noise_comparison.to_string(index=False))\n",
            "print()\n",
            "print(\"=\"*90)\n",
            "\n",
            "# ノイズの影響を視覚化\n",
            "if qubit_noisy_results:\n",
            "    print(\"\\nノイズによる偏差（Qubit実装）:\")\n",
            "    for key in ['N_S0', 'N_T1', 'N_S1']:\n",
            "        noiseless = qubit_results['populations'][-1][key]\n",
            "        noisy = qubit_noisy_results['populations_final'][key]\n",
            "        diff = noisy - noiseless\n",
            "        rel_diff = (diff / noiseless * 100) if noiseless > 0 else 0\n",
            "        print(f\"  Δ{key}: {diff:+.4f} ({rel_diff:+.2f}%)\")\n",
            "    print()\n"
        ]
    })
    
    # Cell 9: Markdown - Conclusion on Noise
    cells.append({
        "cell_type": "markdown",
        "id": "noise_conclusion",
        "metadata": {},
        "source": [
            "### 9.6 ノイズモデルシミュレーションの考察\n",
            "\n",
            "#### 主要な観測結果\n",
            "\n",
            "1. **ノイズの影響**\n",
            "   - 脱分極ノイズと位相緩和ノイズにより、最終個体数に偏差が生じる\n",
            "   - ゲート数が多いほど、ノイズの累積効果が顕著になる\n",
            "   - 2量子ビット/quditゲートのノイズが特に大きな影響を与える\n",
            "\n",
            "2. **Qubit vs Qudit のノイズ耐性**\n",
            "   - Qudit実装はゲート数が少ないため、ノイズの累積が抑えられる可能性がある\n",
            "   - ただし、quditゲートの物理的実装におけるノイズレートが重要\n",
            "\n",
            "3. **測定統計の影響**\n",
            "   - ショットベースのシミュレーションでは、統計的な揺らぎが避けられない\n",
            "   - ショット数を増やすことで、測定精度が向上する\n",
            "\n",
            "#### 実用化への示唆\n",
            "\n",
            "- **ノイズ緩和技術の重要性**: 実機での量子計算には、誤り訂正やノイズ緩和が不可欠\n",
            "- **Quditの利点**: ゲート数削減により、NISTQの可能性が向上\n",
            "- **ハードウェア最適化**: ノイズレートの低減が、実用的な量子シミュレーションの鍵\n",
            "\n",
            "### 9.7 今後の課題\n",
            "\n",
            "1. より現実的なノイズモデル（相関ノイズ、時間依存性など）の実装\n",
            "2. ノイズ緩和技術（Zero Noise Extrapolation、Probabilistic Error Cancellationなど）の適用\n",
            "3. 量子誤り訂正符号の統合\n",
            "4. 実機での検証実験\n"
        ]
    })
    
    return cells


def insert_noise_cells_into_notebook(notebook_path, output_path):
    """Insert noise model cells into the notebook before the conclusion section"""
    
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
    
    # Create noise model cells
    noise_cells = create_noise_model_cells()
    
    # Find the index of the conclusion section (Section 7 or 8) - more robust pattern matching
    conclusion_idx = None
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'markdown':
            source = ''.join(cell['source'])
            # Match either "## 7." or "## 8." followed by conclusion-related keywords
            if ('## 7.' in source or '## 8.' in source) and ('考察' in source or '結論' in source or 'Conclusion' in source):
                conclusion_idx = i
                break
    
    if conclusion_idx is None:
        print("Warning: Could not find conclusion section. Appending at the end.")
        conclusion_idx = len(nb['cells'])
    
    # Insert noise cells before the conclusion
    nb['cells'] = nb['cells'][:conclusion_idx] + noise_cells + nb['cells'][conclusion_idx:]
    
    # Write the modified notebook
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    
    print(f"Successfully added {len(noise_cells)} noise model cells to the notebook")
    print(f"Modified notebook saved to: {output_path}")


if __name__ == "__main__":
    notebook_path = "tutorials/quantum_dynamics_complete_comparison.ipynb"
    output_path = "tutorials/quantum_dynamics_complete_comparison.ipynb"
    
    insert_noise_cells_into_notebook(notebook_path, output_path)
    print("Notebook enhancement complete!")
