# TTA-UC GKSL-Lindblad 検証結果

- 実行時刻(UTC): 20260222T135417Z
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
| elapsed_time | 1.572462 s |
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
| elapsed_time | 39.928484 s |
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
| elapsed_time | 37.244368 s |
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
| elapsed_time | 0.026983 s |
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
| elapsed_time | 16.317347 s |
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
| elapsed_time | 2.579332 s |
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
| elapsed_time | 8.754039 s |
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
| elapsed_time | 12.736956 s |
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
| elapsed_time | 9.288660 s |
| max_trace_error | 8.660e-15 |
| max_particle_error | 3.375e-14 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 2.570e-17 |
| final_min_eigenvalue | 7.976e-05 |
| final_entropy | 3.243625725917 |
| final_purity | 0.103592684591 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| noise_p_depol | 0.01 |
| noise_p_dephasing | 0.005 |
| counts_top10 | 28:225, 30:72, 10:50, 34:39, 31:39, 29:39, 37:38, 1:35, 46:27, 27:25 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=1.7586, T1=1.8882, S1=0.3533 |
| pop_final | S0=1.6724, T1=1.7983, S1=0.5293 |

### qubit_noisy_shot

| 項目 | 値 |
|------|-----|
| method | qubit_gksl_noisy_shot |
| elapsed_time | 15.409039 s |
| max_trace_error | 0.000e+00 |
| max_particle_error | 1.332e-15 |
| final_trace | 0.672713441501 |
| final_hermiticity_error | 4.681e-18 |
| final_min_eigenvalue | 2.586e-06 |
| final_entropy | 2.110852813148 |
| final_purity | 0.084269601954 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| noise_p_depol | 0.01 |
| noise_T1 | 50000.0 |
| noise_t_gate | 300.0 |
| noise_p_reset | 0.005982035946064723 |
| counts_top10 | 28:186, 30:68, 10:35, 1:28, 37:26, 27:25, 31:24, 34:23, 46:21, 29:17 |
| forbidden_count | 332 |
| dim_qubit_space | 256 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=1.3845, T1=1.4497, S1=0.1951 |
| pop_final | S0=1.1947, T1=1.2360, S1=0.2602 |

### unitary

| 項目 | 値 |
|------|-----|
| method | classical_gksl |
| elapsed_time | 1.559972 s |
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
| classical | qudit_noisy_shot | 81 | 0.300322 |
| classical | qubit_noisy_shot | 81 | 0.423945 |
| classical | unitary | 81 | 0.995992 |
| qubit | qudit | 81 | 1.000000 |
| qubit | qudit_shot | 81 | 0.997683 |
| qubit | qubit_shot | 81 | 0.997683 |
| qubit | qudit_noisy_shot | 81 | 0.300467 |
| qubit | qubit_noisy_shot | 81 | 0.424092 |
| qubit | unitary | 81 | 0.996087 |
| qudit | qudit_shot | 81 | 0.997683 |
| qudit | qubit_shot | 81 | 0.997683 |
| qudit | qudit_noisy_shot | 81 | 0.300467 |
| qudit | qubit_noisy_shot | 81 | 0.424092 |
| qudit | unitary | 81 | 0.996087 |
| classical_boson | qubit_boson | 9 | 0.999905 |
| classical_boson | qudit_boson | 9 | 0.999905 |
| qubit_boson | qudit_boson | 9 | 1.000000 |
| qudit_shot | qubit_shot | 81 | 1.000000 |
| qudit_shot | qudit_noisy_shot | 81 | 0.291688 |
| qudit_shot | qubit_noisy_shot | 81 | 0.414682 |
| qudit_shot | unitary | 81 | 0.997874 |
| qubit_shot | qudit_noisy_shot | 81 | 0.291688 |
| qubit_shot | qubit_noisy_shot | 81 | 0.414682 |
| qubit_shot | unitary | 81 | 0.997874 |
| qudit_noisy_shot | qubit_noisy_shot | 81 | 0.943264 |
| qudit_noisy_shot | unitary | 81 | 0.290089 |
| qubit_noisy_shot | unitary | 81 | 0.412671 |

## Unitary vs GKSL 比較

- Unitary最大エントロピー: 2.2204460492503126e-16
- GKSL最終エントロピー: 0.030569978935907162
- 備考: Unitary limit (all dissipation=0) should have near-zero entropy. GKSL with dissipation should show entropy increase over time.

