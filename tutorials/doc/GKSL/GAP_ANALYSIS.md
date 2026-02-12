# Gap Analysis: Implementation Spec vs Theory Document vs Existing Notebook

**Date**: Gap analysis performed against:
- **Implementation Spec** (実装詳細仕様書): `TTA-UC現象GKSL-Lindblad量子ダイナミクス実装詳細仕様書.md`
- **Theory Document** (完全理論書): `TTA-UC現象のGKSL-Lindblad量子ダイナミクス完全理論書.md`
- **Existing Notebook**: `tutorials/quantum_dynamics_complete_comparison.ipynb`

---

## Part A: Theory Document Requirements NOT Covered in the Implementation Spec

### A1. Section 2 (TTA-UC Physical Foundation) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| TTA-UC process overview (4 steps: absorption → ISC → TTA → fluorescence) | §2.1.1 | **NOT COVERED** | The spec jumps directly to GKSL formulation without explaining the full TTA-UC photophysical cycle. §1.1 mentions the system briefly but omits the sensitizer/annihilator distinction and the TTET step. |
| Energy relationship $2E_{T_1} \geq E_{S_1}$ | §2.1.2 | Covered in §4.3 (parameter consistency check) | ✅ Covered |
| Molecular electronic state definitions with physical details (spin multiplicity, electron configuration, typical lifetimes) | §2.2.1–2.2.3 | **PARTIALLY** — Spec §3.2.2 defines matrix representations but omits spin multiplicity, electron configuration, and typical lifetime values | Gap: missing spin physics context (S=0, S=1), triplet sublevel discussion, typical lifetime ranges |
| Spin statistics factor for TTA (1/9 or 1/5) | §2.3.1 | **NOT COVERED** | The spec defines $\gamma_{\text{TTA}}$ but never discusses the spin statistics factor that modifies TTA rates |
| Dexter exchange mechanism for TTA and energy transfer | §2.3.1, §2.3.6 | Spec §4.2.1 mentions Dexter briefly for V, but not for TTA mechanism | Gap: no discussion of Dexter mechanism for TTA |
| Distance dependence of energy transfer $V_{ij} \propto \exp(-2r_{ij}/L)$ | §2.3.6 | **NOT COVERED** | Spec uses constant $V=0.1$ eV without discussing distance dependence |
| Energy gap law for IC $k_{\text{IC}} \propto \exp(-\gamma \Delta E / \hbar\omega_{\text{vib}})$ | §2.3.4 | **NOT COVERED** | Spec §4.2.4 states IC value but not the underlying energy gap law |

### A2. Section 3 (Mathematical Foundations of Open Quantum Systems) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Pure vs mixed state definitions with properties ($\hat{\rho}^2 = \hat{\rho}$, $\text{Tr}[\hat{\rho}^2]$, entropy) | §3.1.1 | **NOT COVERED explicitly** | The spec uses density matrices but never defines the pure/mixed state distinction or the axioms |
| Density operator axioms (Hermiticity, positive semi-definiteness, trace=1) | §3.1.2 | Covered indirectly in §13.1 (verification) but not as foundational axioms | Gap: verification checks exist but foundational axioms are not stated |
| System-environment separation ($\mathcal{H}_{\text{total}} = \mathcal{H}_S \otimes \mathcal{H}_E$) | §3.2.1 | **NOT COVERED** | The spec goes directly to GKSL without establishing the system-environment framework |
| Reduced density operator via partial trace | §3.2.2 | Covered in §7.8 (Statevector method) for ancilla but not as a general foundation | Gap: missing as foundational concept |
| CPTP map definition (linearity, trace preservation, complete positivity) | §3.3.1 | **NOT COVERED** as a standalone definition | Spec mentions CPTP in context of Stinespring but never formally defines the three conditions |
| Kraus representation theorem | §3.3.2 | **NOT COVERED** | Spec uses Stinespring but never mentions Kraus operators or the Kraus theorem |
| Markov approximation conditions (weak coupling, timescale separation) | §3.4.1 | **NOT COVERED** | The spec assumes Markov validity without stating conditions |
| Semigroup property $\mathcal{E}_{t+s} = \mathcal{E}_t \circ \mathcal{E}_s$ | §3.4.2 | **NOT COVERED** | This mathematical property is absent from the spec |

### A3. Section 4.2 (GKSL Mathematical Properties) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Trace preservation proof | §4.2.1 | Covered as verification check (§13.1.1) but **proof not included** | Gap: the theory provides a formal mathematical proof; the spec only has a tolerance check |
| Complete positivity guarantee | §4.2.2 | **NOT COVERED** | The spec never discusses that GKSL form automatically guarantees CPTP |
| Entropy increase / Second law of thermodynamics | §4.2.3 | Covered in §5.6.1 (von Neumann entropy) and §13.1.4 (entropy non-decrease check) | ✅ Covered (as verification) |
| GKSL theorem statement (Gorini-Kossakowski-Sudarshan-Lindblad, 1976) | §4.1.1 | **NOT COVERED** — Spec uses GKSL equation directly without stating the theorem | Gap: missing the formal theorem statement and attribution |
| Anti-commutator notation $\{A,B\} = AB + BA$ | §4.1.1 | **NOT COVERED** — Spec uses expanded form $\frac{1}{2}L^\dagger L \rho + \frac{1}{2}\rho L^\dagger L$ instead | Minor notational difference |

### A4. Section 5.3 (Parameter Physical Ranges) — SIGNIFICANT GAP

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Typical parameter value table with physical origins | §5.3.1 | **PARTIALLY** — Spec §4.2 gives rationale for chosen values but does not provide the full range table | **Key Gap**: Theory gives ranges like $\gamma_{\text{TTA}} = 10^{-3}\text{-}10^{-1}$ eV/ℏ, $\Gamma_{\text{fl}} = 10^{-7}$ eV/ℏ, while spec uses $\Gamma_{\text{fl}} = 0.01$ eV/ℏ (6 orders of magnitude larger than typical). The spec acknowledges this discrepancy in §4.2.2 ("idealized model") but does **not** provide the complete table of typical experimental ranges |
| Unit system conversion table (natural ↔ SI) | §5.3.2 | Covered in Appendix C | ✅ Covered |
| **Parameter discrepancy**: Theory §5.3.1 gives $V = 0.01$ eV as typical, but spec uses $V = 0.1$ eV | §5.3.1 vs Spec §4.1 | **NOT FLAGGED** | The spec uses 10× the theory's typical value without noting this discrepancy |

### A5. Section 6 (Boson Interaction Model) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Physical motivation for boson modes (4 reasons: non-Markov, temperature, vibronic, coherence) | §6.1.1 | Spec §6.1 gives brief overview | ✅ Covered (briefly) |
| Phonon Hilbert space with creation/annihilation operators and commutation relations | §6.1.2 | Spec §6.3.2 provides phonon annihilation operator matrix | **PARTIALLY**: commutation relations $[\hat{a}_k, \hat{a}_{k'}^\dagger] = \delta_{kk'}$ not explicitly stated |
| **Photon mode quantization** ($\hat{\mathbf{E}}(\mathbf{r})$ expansion) | §6.1.2 | **NOT COVERED** | Spec only discusses phonons, completely omits photon field quantization |
| Holstein-type coupling | §6.2.1 | Covered in Spec §6.3.3 | ✅ Covered |
| Huang-Rhys parameter definition | §6.2.1 | Covered in Spec §6.3.3 | ✅ Covered |
| **Peierls-type coupling** ($V_{ij}(1 + \hat{u}_{ij})$) | §6.2.2 | **NOT COVERED** | The spec only covers Holstein coupling, not Peierls coupling for hopping modulation |
| **Electron-photon coupling** (electric dipole interaction, transition dipole operator) | §6.3.1 | **NOT COVERED** | Entirely absent from the spec |
| **Rotating Wave Approximation (RWA)** for electron-photon coupling | §6.3.2 | **NOT COVERED** | Absent from spec |
| **Complete Hamiltonian** (5 terms: $H_{\text{el}} + H_{\text{phonon}} + H_{\text{photon}} + H_{e\text{-ph}} + H_{e\text{-photon}}$) | §6.4.1 | **PARTIALLY** — Spec §6.3.4 has 3 terms only ($H_{\text{el}} + H_{\text{phonon}} + H_{e\text{-ph}}$), missing $H_{\text{photon}}$ and $H_{e\text{-photon}}$ | Gap: photon terms missing |
| **Spectral density function** (Drude-Lorentz, Ohmic models) | §6.5.1 | **NOT COVERED** | Spec does not discuss spectral density models at all |
| **Finite temperature effects** (Bose-Einstein distribution, detailed balance) | §6.5.2 | **NOT COVERED** | No temperature-dependent rates in the spec |
| **HEOM** (basic form, hierarchical index, computational cost, pros/cons) | §6.6 | **NOT COVERED** | Spec §6.6 mentions ODE solver but does not discuss HEOM at all |

### A6. Section 7 (Classical Numerical Methods) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Superoperator formalism (vectorization, superoperator matrix) | §7.1 | Covered in Spec §5.2 | ✅ Covered |
| Matrix exponential method | §7.2.1 | Covered in Spec §5.3.1 | ✅ Covered |
| ODE solver (RK45, BDF) | §7.2.2 | Covered in Spec §5.3.2 | ✅ Covered |
| Sparse matrix methods | §7.2.3 | Covered in Spec §5.3.3 | ✅ Covered |
| **Tensor network methods** (MPS, MPO, TEBD) | §7.3.1 | **NOT COVERED** | Spec does not mention tensor network methods at all |
| **HEOM numerical implementation** (hierarchy enumeration, auxiliary density matrices, convergence) | §7.3.2 | **NOT COVERED** | Not in the spec |
| **Quantum Monte Carlo method** (stochastic wavefunction, quantum jumps, statistical error) | §7.3.3 | **NOT COVERED** | Not in the spec |

### A7. Section 8 (Qubit Implementation) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| 2-qubit encoding | §8.1.1 | Covered in Spec §7.2 | ✅ Covered |
| **Physical subspace projection operator** $\hat{P}_{\text{phys}}^{(i)} = \|00\rangle\langle 00\| + \|01\rangle\langle 01\| + \|10\rangle\langle 10\|$ | §8.1.2 | **NOT COVERED explicitly** | Spec §7.2 mentions "forbidden state |11⟩" but does not define the projection operator |
| Stinespring representation theorem | §8.2.1 | Covered in Spec §7.3 | ✅ Covered |
| Stinespring Lindblad implementation | §8.2.2 | Covered in Spec §7.3.2 | ✅ Covered |
| Unitary circuit (H0 + H_transfer) | §8.3.1 | Covered in Spec §7.5 | ✅ Covered |
| Lindblad circuit | §8.3.2 | Covered in Spec §7.4 | ✅ Covered |
| Trotter decomposition for GKSL | §8.3.3 | Covered in Spec §7.5 | ✅ Covered |
| **Unary encoding** as alternative for boson modes | §8.4.1 | **NOT COVERED** | Spec §8.2.1 only discusses binary encoding |
| **Displacement operator** $\hat{D}(\alpha)$ for phonon modes | §8.4.2 | **NOT COVERED** | Spec §8.4 uses simpler C-RX approach |
| **Forbidden state verification** $P_{\text{forbidden}} < 10^{-8}$ | §11.4.1 | **NOT COVERED in spec** (theory §11.4.1 has it, but spec §13 doesn't include it as a verification step) | Gap: spec has general verification but doesn't include qubit-specific forbidden state check |

### A8. Section 9 (Qudit Implementation) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Qutrit basis | §9.1.1 | Covered in Spec §9.2 | ✅ Covered |
| Qudit advantages (4 items: dimension match, qudit reduction, no forbidden states, natural operators) | §9.1.2 | Covered in Spec §9.8 | ✅ Covered |
| **MQT-Qudits gate set: Generalized Pauli-X** $\hat{X}_d$ | §9.2.1 | **NOT COVERED** | Spec mentions VirtRz, CEx, CustomTwo but does not define generalized Pauli operators |
| **MQT-Qudits gate set: Generalized Pauli-Z** $\hat{Z}_d$ | §9.2.1 | **NOT COVERED** | Not defined in spec |
| **Subspace rotation** $\hat{R}_{mn}(\theta, \phi)$ with Pauli-like generators $\hat{\sigma}_{mn}^x$, $\hat{\sigma}_{mn}^y$ | §9.2.1 | **NOT COVERED** in full generality | Spec uses rotations implicitly but doesn't define the general formalism |
| **2-Qudit gates**: Control gate, SWAP gate | §9.2.2 | **PARTIALLY** — CEx is mentioned but formal definitions of control gate and SWAP are missing | |
| **Givens rotation** definition $G_{mn}(\theta, \phi)$ | §9.3.3 | **NOT COVERED** | Theory defines Givens rotation explicitly; spec mentions "Givens rotation" in §9.8 but provides no definition |
| Stinespring for TTA on qutrits (18×18 unitary) | §9.3.2 | Covered in Spec §9.4.4 | ✅ Covered |
| Boson-mode qudit encoding | §9.4.1 | Covered in Spec §10.2 | ✅ Covered |
| Boson electron-phonon qudit circuit | §9.4.2 | Covered in Spec §10.4 | ✅ Covered |
| **Resource estimation table** (theory gives different numbers from spec) | §9.4.3 | Covered but **values differ**: Theory gives ~7 auxiliary qudit for N=4, Spec gives 26 qubit ancilla | Values are internally consistent within each document but represent different accounting |

### A9. Section 10 (6 Scenario Comparison) — COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Scenario list (6 scenarios) | §10.1 | Covered in Spec §12.1 (expanded to 9 with unitary baselines) | ✅ Covered (and expanded) |
| Classical-NB detailed formulation + cost | §10.2 | Covered in Spec §5 | ✅ Covered |
| Classical-B detailed formulation + cost | §10.3 | Covered in Spec §6 | ✅ Covered |
| Qubit-NB detailed formulation + cost | §10.4 | Covered in Spec §7 | ✅ Covered |
| Qubit-B detailed formulation + cost | §10.5 | Covered in Spec §8 | ✅ Covered |
| Qudit-NB detailed formulation + cost | §10.6 | Covered in Spec §9 | ✅ Covered |
| Qudit-B detailed formulation + cost | §10.7 | Covered in Spec §10 | ✅ Covered |
| **Comprehensive comparison table** (resources, runtime, accuracy, scalability) | §10.8 | Covered in Spec §12.3 | ✅ Covered |
| **Scalability comparison** (N-dependence, practical limits) | §10.8.2 | **NOT COVERED** | Spec doesn't discuss scalability limits per scenario |

### A10. Section 11 (Verification) — MOSTLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Trace preservation | §11.1.1 | Spec §13.1.1 | ✅ Covered |
| Positive semi-definiteness | §11.1.2 | Spec §13.1.2 | ✅ Covered |
| Hermiticity | §11.1.3 | Spec §13.1.3 | ✅ Covered |
| Entropy non-decrease | §11.2.1 | Spec §13.1.4 | ✅ Covered |
| Particle number conservation | §11.2.2 | Spec §13.2.1 | ✅ Covered |
| Trotter error bounds (1st and 2nd order) | §11.3 | **NOT COVERED in spec** | Theory provides explicit commutator-based error bounds; spec Appendix B covers Stinespring error but not Trotter error bounds |
| **Practical Trotter error estimate** with specific parameters | §11.3.2 | **NOT COVERED** | Theory gives concrete numerical estimate ($\epsilon \sim 0.03$); spec doesn't |
| **Forbidden state transition probability** | §11.4.1 | **NOT COVERED** | Missing from spec verification |
| **Stinespring fidelity** $F > 0.99$ | §11.4.2 | **NOT COVERED** | Missing from spec verification |

### A11. Section 12 (Conclusion/Outlook) — PARTIALLY COVERED

| Theory Item | Theory Section | Spec Coverage | Gap |
|---|---|---|---|
| Summary of contributions | §12.1 | Spec has implementation roadmap (§15) instead | Different focus (theory summary vs implementation plan) |
| Heuristic exclusion statement | §12.2 | Covered in Spec §1.4 | ✅ Covered |
| **Future extensions: Non-Markov quantum computation** | §12.3.1 | **NOT COVERED** | |
| **Future extensions: Larger molecular systems** | §12.3.1 | **NOT COVERED** | |
| **Future extensions: Spatial inhomogeneity** | §12.3.1 | **NOT COVERED** | |
| **Experimental verification on quantum hardware** | §12.3.2 | **NOT COVERED** | |
| **Algorithm optimization (VQE, QPE integration)** | §12.3.3 | **NOT COVERED** | |

---

## Part B: Existing Notebook Features NOT Reflected in the Implementation Spec

### B1. Classical Simulator Features

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| `ClassicalSuzukiTrotterSimulator` class with `scipy.linalg.expm` | Cell 5 | Spec §2.4 describes the algorithm but as **unitary baseline only**. The GKSL spec defines a new `ClassicalGKSLSimulator` in §5.5 | ✅ Adequately covered (unitary version documented as baseline, GKSL version as extension) |
| Per-pair 2nd-order symmetric Suzuki-Trotter decomposition (forward/backward H0 → transfer → TTA) | Cell 5, Spec §2.4.1–2.4.2 | Covered in Spec §2.4 | ✅ Covered |
| `expm`-based exact unitary calculation for each pair | Cell 5 | Covered in Spec §2.4.1 | ✅ Covered |

### B2. Qubit Simulator Features

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **Statevector simulation mode** (exact, no shots) | Cell 8 (uses `Statevector`) | Spec §7.8 mentions `Statevector` but the GKSL spec doesn't explicitly specify statevector vs shot-based as separate modes | **Minor gap**: Spec should explicitly distinguish statevector mode (for density matrix reconstruction) vs shot-based mode |
| **Shot-based simulation** with configurable `shots` parameter | Cell 8, Cell 11 | Covered in Spec §7.7 (`simulate` method has `shots` parameter) | ✅ Covered |
| **UnitaryGate version** (`QubitMolecularDynamicsSimulatorUnitary`) | Cell 11 | **NOT COVERED in GKSL spec** | Gap: The spec does not describe a UnitaryGate-based GKSL simulator. The existing notebook's UnitaryGate approach is documented in Spec §2 as baseline but the GKSL extension doesn't address whether UnitaryGate wrapping is used for the Stinespring unitaries |
| **Gate decomposition comparison** (UnitaryGate vs basic gates via KAK decomposition) | Cell 12 | **NOT COVERED in GKSL spec** | Gap: Spec doesn't specify how Stinespring unitaries should be decomposed to basic gates (KAK, etc.) or whether a UnitaryGate wrapper option should exist |

### B3. Qudit Simulator Features

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **Shot-based simulation** (`simulate_shot_based` method) | Cell 17 | Covered in Spec §9.7 | ✅ Covered |
| **`SuzukiTrotterMQTQuditSimulator`** with sparse-aware implementation | Cell 17 | Spec §9.7 defines new `QuditGKSLSimulator` | ✅ Covered (as extension) |
| **`SparseAwareMQTQuditTimeEvolution`** for circuit construction | Cell 17 | Not explicitly mentioned but implied | Minor gap |

### B4. Noise Models

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **Qubit depolarizing noise** (2-qubit gates only, 1.0%) | Cell 14–15 | **NOT COVERED** | **Significant gap**: The GKSL spec has no section on hardware noise models. It defines Lindblad dissipation (physical noise) but not gate-level hardware noise (depolarizing, thermal relaxation). These are distinct: Lindblad = physical TTA/fluorescence dissipation; hardware noise = gate imperfections |
| **Qubit thermal relaxation** (T₁=50μs, T₂=70μs) | Cell 14 | **NOT COVERED** | Same gap as above |
| **`QubitMolecularDynamicsSimulatorNoisy`** class | Cell 15 | **NOT COVERED** | No noisy GKSL simulator class defined |
| **Qudit depolarizing noise** (2-qudit gates only, 1.0%) | Cell 22–23 | **NOT COVERED** | Same gap: no hardware noise model for GKSL qudit circuits |
| **Qudit dephasing noise** | Cell 22 | **NOT COVERED** | Same gap |
| **`NoisyQuditMolecularDynamicsSimulator`** class | Cell 23 | **NOT COVERED** | No noisy GKSL qudit simulator class defined |
| **Noise impact comparison** (noiseless vs noisy) | Cell 28–30 | **NOT COVERED** | Spec doesn't include noise impact evaluation framework |

### B5. Circuit Visualization

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **Qubit circuit visualization** via `circuit_drawer` (matplotlib output) | Cell 13 | **NOT COVERED** | Spec §12 covers population/entropy plots but not circuit diagram visualization |
| **Qudit circuit visualization** via `visualize_circuit` | Cell 20, 26 | **NOT COVERED** | Same gap |
| **Side-by-side circuit comparison** (UnitaryGate vs decomposed, CustomTwo vs decomposed) | Cell 13, 26 | **NOT COVERED** | Same gap |

### B6. Gate Decomposition Details

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **`comparison_helpers` module** with `count_gates_by_type`, `decompose_qiskit_unitary_gates`, `decompose_qudit_customtwo_gates_to_circuit` | Cells 12, 24, 25 | **NOT COVERED** | Gap: The spec doesn't define helper functions for gate counting or decomposition analysis |
| **Qubit: UnitaryGate vs basic gate decomposition** (KAK decomposition) | Cell 10, 12 | **NOT COVERED** | Spec doesn't specify decomposition strategy for GKSL circuits |
| **Qudit: CustomTwo vs basic gate decomposition** (IntegratedSparseCompilerV2) | Cell 21, 24 | Spec §9.4 mentions "CustomTwo gate" and decomposition to basic gates but **doesn't specify the compiler** | Gap: no reference to `IntegratedSparseCompilerV2` or ~6 gates/CustomTwo metric |
| **Gate count comparison tables** | Cells 12, 24 | **NOT COVERED** | No gate counting specification for GKSL circuits |

### B7. Comparison Framework

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **`calculate_population_error`** function (pairwise error between methods) | Cell 32 | **NOT COVERED explicitly** — Spec §13.3.1 defines tolerance but no specific function | Minor gap |
| **Error time-evolution plots** (semilogy) | Cell 33 | **NOT COVERED** | Spec §12.2 defines GKSL-specific plots but not error time-evolution |
| **Side-by-side population dynamics** (1×3 subplot: Classical/Qubit/Qudit) | Cell 34 | **PARTIALLY** — Spec §12.2.3 defines 2×3 grid for 6 scenarios | Covered with different layout |
| **Comprehensive comparison table** (pandas DataFrame with method, resources, dimensions, gates, depth, accuracy, time) | Cell 35 | Covered in Spec §12.3 | ✅ Covered |
| **Accuracy comparison** (Classical as ground truth) | Cell 32 | Covered in Spec §13.3 | ✅ Covered |
| **Resource comparison** | Cell 35 | Covered in Spec §12.3 | ✅ Covered |

### B8. Per-molecule Population Tracking

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **`per_molecule_populations`** dict with `S0_per_mol`, `T1_per_mol`, `S1_per_mol` arrays | Cell 3 (output format), Spec §2.2.1 | Covered in Spec §2.2.1 (data structure) and §5.4 (GKSL calculation) | ✅ Covered |
| **`plot_per_molecule_populations`** (2×2 subplot per molecule) | Cell 6 | Covered in Spec §2.6.2 | ✅ Covered |

### B9. Per-pair Trotter Decomposition

| Notebook Feature | Notebook Location | Spec Coverage | Gap |
|---|---|---|---|
| **Per-pair Trotter decomposition** (H_transfer and H_TTA applied per neighbor pair, not as full-system operators) | Cell 5 | Covered in Spec §2.4 | ✅ Covered |
| **Forward/backward ordering** (symmetric 2nd-order) | Cell 5 | Covered in Spec §2.4.2 | ✅ Covered |

---

## Summary of Critical Gaps

### HIGH PRIORITY (Must be added to the implementation spec to satisfy the theory document):

1. **Section 3 foundations** (Theory §3): Missing density operator axioms, CPTP maps, Kraus theorem, Markov approximation conditions, semigroup property
2. **Boson interaction gaps** (Theory §6): Missing photon quantization, Peierls coupling, electron-photon coupling, spectral density functions, finite temperature effects, HEOM
3. **Classical numerical methods** (Theory §7): Missing tensor network, HEOM implementation, quantum Monte Carlo descriptions
4. **Parameter physical ranges** (Theory §5.3): Missing comprehensive table of typical experimental ranges and discrepancy notes
5. **Hardware noise models** (Notebook): Completely missing gate-level noise (depolarizing, thermal relaxation, dephasing) for both qubit and qudit GKSL implementations
6. **GKSL theorem formal statement** (Theory §4.1.1): Missing theorem attribution and statement

### MEDIUM PRIORITY:

7. **Gate decomposition strategy** for GKSL circuits (UnitaryGate vs basic gates for qubit; CustomTwo vs basic gates for qudit)
8. **Physical subspace projection** and forbidden state verification for qubit implementation
9. **MQT-Qudits formal gate set** (generalized Pauli-X, Pauli-Z, subspace rotation formalism, Givens rotation definition)
10. **Trotter error bounds** (commutator-based formal bounds)
11. **Stinespring fidelity** target ($F > 0.99$)
12. **Scalability comparison** across scenarios
13. **Future extensions** (non-Markov, larger systems, VQE/QPE integration)
14. **Circuit visualization** specifications

### LOW PRIORITY:

15. **TTA-UC photophysical cycle overview** (sensitizer/annihilator, TTET step)
16. **Spin statistics factor** for TTA
17. **Distance dependence** of coupling
18. **Energy gap law** for IC
19. **Unary encoding** as alternative for boson modes
20. **Displacement operator** for phonon modes
