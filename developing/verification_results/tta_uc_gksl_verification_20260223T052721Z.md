# TTA-UC GKSL-Lindblad 検証結果

- 実行時刻(UTC): 20260223T052721Z
- 総合判定: **PASS**
- シナリオ: classical, qubit, qudit, classical_boson, qubit_boson, qudit_boson, qudit_shot, qubit_shot, qudit_noisy_shot, qubit_noisy_shot, unitary
- t_max: 5.0
- n_steps: 5
- initial_state: edge_triplet
- n_shots (shot系): 1000
- seed: 42

## 評価基準

- max_trace_error: Tr(ρ)=1からの最大偏差
- max_particle_error: N_S0+N_T1+N_S1=N_moleculesからの最大偏差
- density_validation: 密度行列の物理的妥当性(エルミート性、正値性、トレース)
- fidelity: 異シナリオ間の量子状態忠実度

## シナリオ別結果

### classical

| 項目 | 値 |
|------|-----|
| method | classical_gksl |
| elapsed_time | 1.728863 s |
| max_trace_error | 8.882e-16 |
| max_particle_error | 3.553e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 8.107e-17 |
| final_min_eigenvalue | -2.596e-16 |
| final_entropy | 0.030569978936 |
| final_purity | 0.992258104547 |
| density_validation | PASS ✅ |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=2.0004, T1=1.9993, S1=0.0003 |
| pop_final | S0=2.0039, T1=1.9924, S1=0.0037 |

### qubit

| 項目 | 値 |
|------|-----|
| method | qubit_gksl |
| elapsed_time | 7.806767 s |
| max_trace_error | 1.332e-15 |
| max_particle_error | 4.441e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 2.944e-16 |
| final_min_eigenvalue | -1.767e-16 |
| final_entropy | 0.030046553479 |
| final_purity | 0.992445599524 |
| density_validation | PASS ✅ |
| dim_qubit_space | 256 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=2.0004, T1=1.9993, S1=0.0003 |
| pop_final | S0=2.0039, T1=1.9926, S1=0.0036 |

### qudit

| 項目 | 値 |
|------|-----|
| method | qudit_gksl |
| elapsed_time | 6.517806 s |
| max_trace_error | 1.332e-15 |
| max_particle_error | 4.441e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 2.848e-16 |
| final_min_eigenvalue | -2.905e-16 |
| final_entropy | 0.030046553479 |
| final_purity | 0.992445599524 |
| density_validation | PASS ✅ |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=2.0004, T1=1.9993, S1=0.0003 |
| pop_final | S0=2.0039, T1=1.9926, S1=0.0036 |

### classical_boson

| 項目 | 値 |
|------|-----|
| method | classical_gksl_boson |
| elapsed_time | 0.026429 s |
| max_trace_error | 4.441e-16 |
| max_particle_error | 8.882e-16 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 6.648e-20 |
| final_min_eigenvalue | 0.000e+00 |
| final_entropy | 0.722017046577 |
| final_purity | 0.628756985844 |
| density_validation | PASS ✅ |
| pop_initial | S0=0.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=0.1425, T1=1.7220, S1=0.1355 |
| pop_final | S0=0.2297, T1=1.5592, S1=0.2111 |

### qubit_boson

| 項目 | 値 |
|------|-----|
| method | qubit_gksl_boson |
| elapsed_time | 2.271600 s |
| max_trace_error | 1.110e-15 |
| max_particle_error | 2.220e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 1.769e-17 |
| final_min_eigenvalue | 0.000e+00 |
| final_entropy | 0.729042696056 |
| final_purity | 0.627284692103 |
| density_validation | PASS ✅ |
| pop_initial | S0=0.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=0.1441, T1=1.7211, S1=0.1348 |
| pop_final | S0=0.2321, T1=1.5579, S1=0.2100 |

### qudit_boson

| 項目 | 値 |
|------|-----|
| method | qudit_gksl_boson |
| elapsed_time | 0.564178 s |
| max_trace_error | 1.110e-15 |
| max_particle_error | 2.220e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 1.769e-17 |
| final_min_eigenvalue | 0.000e+00 |
| final_entropy | 0.729042696056 |
| final_purity | 0.627284692103 |
| density_validation | PASS ✅ |
| pop_initial | S0=0.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=0.1441, T1=1.7211, S1=0.1348 |
| pop_final | S0=0.2321, T1=1.5579, S1=0.2100 |

### qudit_shot

| 項目 | 値 |
|------|-----|
| method | qudit_gksl_shot |
| elapsed_time | 8.861834 s |
| max_trace_error | 2.065e-14 |
| max_particle_error | 8.171e-14 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 1.158e-16 |
| final_min_eigenvalue | -5.013e-15 |
| final_entropy | 0.014427214862 |
| final_purity | 0.996008000000 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| counts_top10 | 28:594, 30:173, 10:139, 12:66, 4:13, 36:13, 18:2 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=2.0010, T1=1.9980, S1=0.0010 |
| pop_final | S0=2.0020, T1=1.9960, S1=0.0020 |

### qubit_shot

| 項目 | 値 |
|------|-----|
| method | qubit_gksl_shot |
| elapsed_time | 13.799886 s |
| max_trace_error | 2.065e-14 |
| max_particle_error | 8.171e-14 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 2.911e-17 |
| final_min_eigenvalue | -5.052e-15 |
| final_entropy | 0.014427214862 |
| final_purity | 0.996008000000 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| counts_top10 | 28:594, 30:173, 10:139, 12:66, 4:13, 36:13, 18:2 |
| forbidden_count | 0 |
| dim_qubit_space | 256 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=2.0010, T1=1.9980, S1=0.0010 |
| pop_final | S0=2.0020, T1=1.9960, S1=0.0020 |

### qudit_noisy_shot

| 項目 | 値 |
|------|-----|
| method | qudit_gksl_noisy_shot |
| elapsed_time | 9.178898 s |
| max_trace_error | 8.660e-15 |
| max_particle_error | 3.464e-14 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 2.999e-17 |
| final_min_eigenvalue | 3.283e-06 |
| final_entropy | 3.241609225898 |
| final_purity | 0.103154143759 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| noise_p_depol | 0.01 |
| noise_p_dephasing | 0.005 |
| counts_top10 | 28:232, 10:59, 30:54, 46:39, 27:39, 29:37, 37:35, 1:34, 34:30, 31:30 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=1.7597, T1=1.8854, S1=0.3549 |
| pop_final | S0=1.6964, T1=1.7529, S1=0.5506 |

### qubit_noisy_shot

| 項目 | 値 |
|------|-----|
| method | qubit_gksl_noisy_shot |
| elapsed_time | 14.615161 s |
| max_trace_error | 0.000e+00 |
| max_particle_error | 1.776e-15 |
| max_trace_deficit | 3.441e-01 |
| final_trace | 0.655942855770 |
| final_hermiticity_error | 2.002e-17 |
| final_min_eigenvalue | 0.000e+00 |
| final_entropy | 2.130089624521 |
| final_purity | 0.073634831730 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| noise_p_depol | 0.01 |
| noise_p_dephasing | 0.005 |
| noise_T1 | None |
| noise_t_gate | 300.0 |
| noise_p_reset | 0.0 |
| counts_top10 | 28:202, 10:52, 30:40, 34:26, 46:25, 1:24, 31:22, 27:20, 16:15, 55:15 |
| forbidden_count | 347 |
| dim_qubit_space | 256 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=1.3859, T1=1.4490, S1=0.2025 |
| pop_final | S0=1.1389, T1=1.1972, S1=0.2877 |

### unitary

| 項目 | 値 |
|------|-----|
| method | classical_gksl |
| elapsed_time | 1.570344 s |
| max_trace_error | 4.441e-16 |
| max_particle_error | 1.776e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 9.061e-17 |
| final_min_eigenvalue | -2.709e-16 |
| final_entropy | 0.000000000000 |
| final_purity | 1.000000000000 |
| density_validation | PASS ✅ |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_final | S0=2.0000, T1=2.0000, S1=0.0000 |

## シナリオ間忠実度 (Fidelity)

| Scenario A | Scenario B | dim | Fidelity |
|------------|------------|-----|----------|
| classical | qubit | 81 | 0.999995 |
| classical | qudit | 81 | 0.999995 |
| classical | qudit_shot | 81 | 0.997638 |
| classical | qubit_shot | 81 | 0.997638 |
| classical | qudit_noisy_shot | 81 | 0.295712 |
| classical | qubit_noisy_shot | 81 | 0.393445 |
| classical | unitary | 81 | 0.995992 |
| qubit | qudit | 81 | 1.000000 |
| qubit | qudit_shot | 81 | 0.997683 |
| qubit | qubit_shot | 81 | 0.997683 |
| qubit | qudit_noisy_shot | 81 | 0.295836 |
| qubit | qubit_noisy_shot | 81 | 0.393550 |
| qubit | unitary | 81 | 0.996087 |
| qudit | qudit_shot | 81 | 0.997683 |
| qudit | qubit_shot | 81 | 0.997683 |
| qudit | qudit_noisy_shot | 81 | 0.295836 |
| qudit | qubit_noisy_shot | 81 | 0.393550 |
| qudit | unitary | 81 | 0.996087 |
| classical_boson | qubit_boson | 9 | 0.999905 |
| classical_boson | qudit_boson | 9 | 0.999905 |
| qubit_boson | qudit_boson | 9 | 1.000000 |
| qudit_shot | qubit_shot | 81 | 1.000000 |
| qudit_shot | qudit_noisy_shot | 81 | 0.287033 |
| qudit_shot | qubit_noisy_shot | 81 | 0.385865 |
| qudit_shot | unitary | 81 | 0.997874 |
| qubit_shot | qudit_noisy_shot | 81 | 0.287033 |
| qubit_shot | qubit_noisy_shot | 81 | 0.385865 |
| qubit_shot | unitary | 81 | 0.997874 |
| qudit_noisy_shot | qubit_noisy_shot | 81 | 0.938379 |
| qudit_noisy_shot | unitary | 81 | 0.283639 |
| qubit_noisy_shot | unitary | 81 | 0.382658 |

## Unitary vs GKSL 比較

- Unitary最大エントロピー: 2.2204460492503126e-16
- GKSL最終エントロピー: 0.030569978935907162
- 備考: Unitary limit (all dissipation=0) should have near-zero entropy. GKSL with dissipation should show entropy increase over time.

