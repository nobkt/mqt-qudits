# TTA-UC GKSL-Lindblad 検証結果

- 実行時刻(UTC): 20260222T032847Z
- 総合判定: **FAIL**
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

### qubit

| 項目 | 値 |
|------|-----|
| method | qubit_gksl |
| elapsed_time | 36.040261 s |
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

### qudit

| 項目 | 値 |
|------|-----|
| method | qudit_gksl |
| elapsed_time | 36.940134 s |
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
| elapsed_time | 2.471817 s |
| max_trace_error | 3.273e-13 |
| max_particle_error | 6.546e-13 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 9.199e-22 |
| final_min_eigenvalue | 0.000e+00 |
| final_entropy | 0.722017046655 |
| final_purity | 0.628756985525 |
| density_validation | PASS ✅ |
| pop_initial | S0=0.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=0.1425, T1=1.7220, S1=0.1355 |
| pop_final | S0=0.2297, T1=1.5592, S1=0.2111 |

### qubit_boson

| 項目 | 値 |
|------|-----|
| method | qubit_gksl_boson |
| elapsed_time | 2.536081 s |
| max_trace_error | 1.110e-15 |
| max_particle_error | 2.220e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 5.034e-17 |
| final_min_eigenvalue | -1.140e-18 |
| final_entropy | 0.731610751487 |
| final_purity | 0.627215433153 |
| density_validation | PASS ✅ |
| pop_initial | S0=0.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=0.1441, T1=1.7211, S1=0.1348 |
| pop_final | S0=0.2321, T1=1.5579, S1=0.2100 |

### qudit_boson

| 項目 | 値 |
|------|-----|
| method | qudit_gksl_boson |
| elapsed_time | 2.076101 s |
| max_trace_error | 1.110e-15 |
| max_particle_error | 2.220e-15 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 5.034e-17 |
| final_min_eigenvalue | -1.140e-18 |
| final_entropy | 0.731610751487 |
| final_purity | 0.627215433153 |
| density_validation | PASS ✅ |
| pop_initial | S0=0.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=0.1441, T1=1.7211, S1=0.1348 |
| pop_final | S0=0.2321, T1=1.5579, S1=0.2100 |

### qudit_shot

| 項目 | 値 |
|------|-----|
| method | qudit_gksl_shot |
| elapsed_time | 8.413427 s |
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
| elapsed_time | 8.599035 s |
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

### qudit_noisy_shot

| 項目 | 値 |
|------|-----|
| method | qudit_gksl_noisy_shot |
| elapsed_time | 9.261956 s |
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
| elapsed_time | 9.786740 s |
| max_trace_error | 8.660e-15 |
| max_particle_error | 3.197e-14 |
| final_trace | 1.000000000000 |
| final_hermiticity_error | 1.535e-17 |
| final_min_eigenvalue | 9.917e-05 |
| final_entropy | 3.162288251635 |
| final_purity | 0.112598853086 |
| density_validation | PASS ✅ |
| n_shots | 1000 |
| noise_p_depol | 0.01 |
| noise_T1 | 50000.0 |
| noise_t_gate | 300.0 |
| noise_p_reset | 0.005982035946064723 |
| counts_top10 | 28:198, 30:68, 10:64, 1:41, 34:39, 27:39, 46:38, 37:35, 29:35, 55:30 |
| pop_initial | S0=2.0000, T1=2.0000, S1=0.0000 |
| pop_mid | S0=1.7761, T1=1.8262, S1=0.3977 |
| pop_final | S0=1.7220, T1=1.7487, S1=0.5293 |

## シナリオ間忠実度 (Fidelity)

| Scenario A | Scenario B | dim | Fidelity |
|------------|------------|-----|----------|
| qubit | qudit | 81 | 1.000000 |
| qubit | qudit_shot | 81 | 0.997683 |
| qubit | qubit_shot | 81 | 0.997683 |
| qubit | qudit_noisy_shot | 81 | 0.300467 |
| qubit | qubit_noisy_shot | 81 | 0.319660 |
| qudit | qudit_shot | 81 | 0.997683 |
| qudit | qubit_shot | 81 | 0.997683 |
| qudit | qudit_noisy_shot | 81 | 0.300467 |
| qudit | qubit_noisy_shot | 81 | 0.319660 |
| qubit_boson | qudit_boson | 36 | 1.000000 |
| qudit_shot | qubit_shot | 81 | 1.000000 |
| qudit_shot | qudit_noisy_shot | 81 | 0.291688 |
| qudit_shot | qubit_noisy_shot | 81 | 0.311509 |
| qubit_shot | qudit_noisy_shot | 81 | 0.291688 |
| qubit_shot | qubit_noisy_shot | 81 | 0.311509 |
| qudit_noisy_shot | qubit_noisy_shot | 81 | 0.954653 |

## エラー

### classical

```
PhysicsViolationError: Step 5: Min eigenvalue = -1.83443704434125e-10 (negative)
```

<details><summary>Traceback</summary>

```
Traceback (most recent call last):
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py", line 311, in run_verification_workflow
    report, rho_final = _run_single_scenario(
                        ^^^^^^^^^^^^^^^^^^^^^
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py", line 237, in _run_single_scenario
    validation = _validate_result(result, params, n_steps)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py", line 123, in _validate_result
    density_validation = validate_density_matrix(result["rho_final"], step=n_steps)
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/gksl_validation.py", line 82, in validate_density_matrix
    raise PhysicsViolationError("; ".join(errors))
gksl_validation.PhysicsViolationError: Step 5: Min eigenvalue = -1.83443704434125e-10 (negative)

```

</details>

### unitary

```
PhysicsViolationError: Step 5: Min eigenvalue = -2.2472774401048736e-10 (negative)
```

<details><summary>Traceback</summary>

```
Traceback (most recent call last):
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py", line 311, in run_verification_workflow
    report, rho_final = _run_single_scenario(
                        ^^^^^^^^^^^^^^^^^^^^^
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py", line 237, in _run_single_scenario
    validation = _validate_result(result, params, n_steps)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/run_tta_uc_gksl_verification.py", line 123, in _validate_result
    density_validation = validate_density_matrix(result["rho_final"], step=n_steps)
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/A23321P/work/myQudit/mqt-qudits/tutorials/gksl_validation.py", line 82, in validate_density_matrix
    raise PhysicsViolationError("; ".join(errors))
gksl_validation.PhysicsViolationError: Step 5: Min eigenvalue = -2.2472774401048736e-10 (negative)

```

</details>

