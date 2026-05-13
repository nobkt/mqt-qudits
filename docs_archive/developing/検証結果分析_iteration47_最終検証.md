# Iteration 47 最終検証分析: 全コード・全出力の詳細検証

## 分析日時
2026-03-06

## 1. 結論（先頭に記載）

**quantum_dynamics_gksl_comparison.ipynb および呼び出す全ての機能に対して、これ以上の修正は不要です。**

以下の詳細検証に基づき、GKSLシミュレーションコード全体、ノートブックの全セル、回帰テストの全31チェックが数学的に正確であることを独立に確認しました。ヒューリスティックな処理やfallbackは一切使用されていません。

## 2. 検証方法

### 2.1 独立実行による数値再現
- `run_tta_uc_gksl_regression_test.py` を独立実行し、全31チェックがPASSすることを確認
- 前回（20260306T004704Z）と今回（20260306T020159Z）の数値が完全一致（決定論的再現性確認）

### 2.2 数学的コード精査
以下の全モジュールのソースコードを行単位で精査:
- `stinespring_utils.py` (104行)
- `gksl_math_utils.py` (293行)
- `gksl_physical_parameters.py` (175行)
- `classical_gksl_simulator.py` (178行)
- `qudit_gksl_simulator.py` (212行)
- `qubit_gksl_simulator.py` (389行)
- `gksl_validation.py` (150行)
- `run_tta_uc_gksl_regression_test.py` (1182行)
- `quantum_dynamics_gksl_comparison.ipynb` (全40セル)

### 2.3 ノートブック出力の検証
- 全20コードセル（Cell 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 30, 32, 34, 36, 38）の出力を確認
- 全20マークダウンセル（Cell 0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 28, 29, 31, 33, 35, 37, 39）の内容を確認

## 3. 数学的正確性の詳細検証

### 3.1 Stinespring dilation (`stinespring_utils.py`)

**構成**: $G = \begin{pmatrix} 0 & L^\dagger \\ L & 0 \end{pmatrix}$, $U = \exp(-i\sqrt{\Delta t} \cdot G)$

**検証**:
$$E(\rho) = \text{Tr}_{\text{env}}[U(|0\rangle\langle 0| \otimes \rho)U^\dagger] = U_{00}\rho U_{00}^\dagger + U_{10}\rho U_{10}^\dagger$$

展開（$\theta = \sqrt{\Delta t}$）:
- $U_{00} = \cos(\theta\sqrt{L^\dagger L}) \approx I - \frac{\Delta t}{2}L^\dagger L + O(\Delta t^2)$
- $U_{10} = -iL(L^\dagger L)^{-1/2}\sin(\theta\sqrt{L^\dagger L}) \approx -i\sqrt{\Delta t} \cdot L + O(\Delta t^{3/2})$

$$E(\rho) \approx \rho - \frac{\Delta t}{2}(L^\dagger L \rho + \rho L^\dagger L) + \Delta t \cdot L\rho L^\dagger = \rho + \Delta t(L\rho L^\dagger - \frac{1}{2}\{L^\dagger L, \rho\})$$

✅ **GKSL dissipator のリーディングオーダー再現を確認**

**コード実装確認**:
- `env⊗sys` 順序の kron が `G` のブロック構造と一致 ✅
- 部分トレース `rho_prime[:d,:d] + rho_prime[d:,d:]` が $\sum_k \langle k|_\text{env} \rho' |k\rangle_\text{env}$ と一致 ✅
- ユニタリ性チェック `||U†U - I||_F < 1e-10` あり ✅

### 3.2 GKSL超演算子 (`stinespring_utils.py: build_gksl_superoperator`)

**ベクトル化**: 列優先 (Fortran-order) $\text{vec}(X) = X.\text{flatten}(\text{order='F'})$

恒等式: $\text{vec}(AXB) = (B^T \otimes A)\text{vec}(X)$

**ハミルトニアン部分**:
$$\mathcal{L}_H = -\frac{i}{\hbar}(I \otimes H - H^T \otimes I)$$

検証: $\text{vec}(-i[H,\rho]/\hbar) = -\frac{i}{\hbar}(\text{vec}(H\rho I) - \text{vec}(I\rho H)) = -\frac{i}{\hbar}(I^T \otimes H - H^T \otimes I)\text{vec}(\rho)$

コード: `(-1j / hbar) * (np.kron(I, H_total) - np.kron(H_total.T, I))` ✅

**散逸部分**:
$$\mathcal{L}_D = \sum_\alpha [L_\alpha^* \otimes L_\alpha - \frac{1}{2}(I \otimes L_\alpha^\dagger L_\alpha + (L_\alpha^\dagger L_\alpha)^T \otimes I)]$$

検証:
- $\text{vec}(L\rho L^\dagger) = (L^\dagger)^T \otimes L = L^* \otimes L$ → `np.kron(L_op.conj(), L_op)` ✅
- $\text{vec}(-\frac{1}{2}L^\dagger L\rho) = -\frac{1}{2}I \otimes L^\dagger L$ → `-0.5 * np.kron(I, LdL)` ✅
- $\text{vec}(-\frac{1}{2}\rho L^\dagger L) = -\frac{1}{2}(L^\dagger L)^T \otimes I$ → `-0.5 * np.kron(LdL.T, I)` ✅

### 3.3 Lindblad演算子 (`gksl_math_utils.py: build_lindblad_operators`)

| チャネル | 演算子 | 数 | レート |
|----------|--------|-----|--------|
| TTA | $\|2\rangle_i\langle 1\| \otimes \|0\rangle_j\langle 1\|$ (+ 逆) | $2 \times (N-1)$ | $\gamma_\text{TTA}/2$ |
| 蛍光 | $\|0\rangle_i\langle 2\|$ | $N$ | $\Gamma_\text{fl}$ |
| 燐光 | $\|0\rangle_i\langle 1\|$ | $N$ | $\Gamma_\text{ph}$ |
| 内部変換 | $\|0\rangle_i\langle 2\|$ | $N$ | $k_\text{IC}$ |
| ISC S→T | $\|1\rangle_i\langle 2\|$ | $N$ | $k_\text{ISC\_ST}$ |
| ISC T→S | $\|0\rangle_i\langle 1\|$ | $N$ | $k_\text{ISC\_TS}$ |

合計: $2(N-1) + 5N = 26$ (N=4の場合) ✅

各 $L_\alpha = \sqrt{\gamma_\alpha} \cdot l_\alpha$ で $\sqrt{\gamma}$ 因子を含む ✅

### 3.4 Trotter分割 (`qudit_gksl_simulator.py: _trotter_step`)

```
exp(L dt) ≈ exp(L_H dt/2) · ∏_{α=1→n} E_α(dt/2) · ∏_{α=n→1} E_α(dt/2) · exp(L_H dt/2)
```

**Strang分割** (H-D): $O(\Delta t^3)$/ステップ → $O(\Delta t^2)$ 全体 ✅
**回文順序積** (Lindblad): Lie-Trotter交換子誤差消去 → $O(\Delta t^3)$/ステップ ✅
**Stinespring近似** (各チャネル): $O(\Delta t^2)$/チャネル/ステップ → **$O(\Delta t)$ 全体** ← 支配的 ✅

### 3.5 初期状態 (`prepare_initial_state`)

"edge_triplet" → $|1,0,0,1\rangle$: 分子0,3がT₁、分子1,2がS₀

インデックス: $1 \times 3^3 + 1 = 28$ ✅

### 3.6 人口計算 (`compute_populations_from_density_matrix`)

各基底状態 $|m_0, m_1, ..., m_{N-1}\rangle$ について:
- インデックス $\text{idx} = \sum_i m_i \cdot d^i$
- ループ `for mol in range(N-1, -1, -1): remainder % d` で正しく分解 ✅
- 粒子数保存: $N_{S_0} + N_{T_1} + N_{S_1} = N \cdot \text{Tr}[\rho] = N$ ✅

### 3.7 古典シミュレータ (`classical_gksl_simulator.py`)

- `scipy.sparse.linalg.expm_multiply` による厳密な $\exp(\mathcal{L}t)\text{vec}(\rho_0)$ 計算 ✅
- CPTP性はLindblad-GKS定理により構造的に保証 ✅
- 独立参照解として正しく機能 ✅

### 3.8 ベクトル化の一貫性

- `vectorize_density_matrix`: `rho.flatten(order='F')` (列優先) ✅
- `unvectorize_density_matrix`: `vec.reshape((dim, dim), order='F')` ✅
- `build_gksl_superoperator`: $B^T \otimes A$ 恒等式に基づく構成 ✅
- 全コンポーネントで同一のベクトル化規約を使用 ✅

## 4. 回帰テスト全31チェックの検証

### 独立実行結果（2回目: 20260306T020159Z）

| # | チェック項目 | 結果 | 数値 |
|---|---|---|---|
| 1 | トレース距離単調減少 | ✅ PASS | — |
| 2 | 収束次数 ≈ 1.0 | ✅ PASS | 平均: 1.2010 |
| 3 | 忠実度単調増加 | ✅ PASS | — |
| 4 | 密度行列品質 | ✅ PASS | tr_err < 1e-14 |
| 5 | n_steps=100 で T < 1e-3 | ✅ PASS | T = 7.32e-04 |
| 6 | Cell 27 テーブル値整合性 | ✅ PASS | 全値比率 0.5-2.0 以内 |
| 9 | Cell 27 n_steps=200 値 | ✅ PASS | ~3.6e-04 確認 |
| 10 | 粒子数保存 | ✅ PASS | Cl: 2.22e-14, Qd: 1.24e-14 |
| 11 | 中間時刻人口動態 | ✅ PASS | max diff: 1.60e-03 |
| 12 | 2分子系収束次数 | ✅ PASS | Rate: 1.0022 |
| 13 | ユニタリ限界精度 | ✅ PASS | T: 2.22e-14 |
| 14 | 固有値スペクトル一致 | ✅ PASS | max diff: 6.42e-04 |
| 15 | Qudit Shot-DM一致 | ✅ PASS | T: 2.81e-03 |
| 16 | Qubit Shot-DM一致 | ✅ PASS | T: 2.81e-03 |
| 17 | Qudit Noisy(p=0)=DM | ✅ PASS | T: 0.00e+00 |
| 18 | Qubit Noisy(p=0)=DM | ✅ PASS | T: 0.00e+00 |
| 19 | Qudit DM アンシラ数 (4分子) | ✅ PASS | 26=26 |
| 20 | Qudit DM アンシラ数 (2分子) | ✅ PASS | 12=12 |
| 21 | Qudit DM ゲート数 (4分子) | ✅ PASS | 59=59 |
| 22 | Qudit DM ゲート数 (2分子) | ✅ PASS | 27=27 |
| 23 | Qudit Shot ゲート数 (4分子) | ✅ PASS | 59=59 |
| 24 | Qudit Shot ゲート数 (2分子) | ✅ PASS | 27=27 |
| 25 | Qubit Shot ゲート数推定あり | ✅ PASS | 158=158 |
| 26 | Qudit-Qubit Shot 同一結果 | ✅ PASS | T: 2.96e-20 |
| 27 | Boson Qudit 収束次数 | ✅ PASS | Rate: 1.0009 |
| 28 | Boson Qubit 収束次数 | ✅ PASS | Rate: 1.0009 |
| 29 | Boson Qudit-Qubit一致 | ✅ PASS | T: 1.69e-21 |
| 30 | Boson 密度行列品質 | ✅ PASS | tr_err: 4.44e-16 |
| 32 | 非Bosonに_U_stinesなし | ✅ PASS | — |
| 33 | 誤差係数 T/dt の安定性 | ✅ PASS | 偏差: 9.02e-03 |
| 34 | Stinespring チャネル誤差 O(dt²) | ✅ PASS | 偏差: 1.65e-04 |

### 数値再現性の確認

前回実行（20260306T004704Z）と今回実行（20260306T020159Z）の主要数値を比較:

| 項目 | 前回 | 今回 | 差異 |
|------|------|------|------|
| T(n=100) | 7.319238e-04 | 7.319238e-04 | 0 |
| T(n=200) | 3.626913e-04 | 3.626913e-04 | 0 |
| Rate avg | 1.2010 | 1.2010 | 0 |
| Qudit-Qubit Shot T | 3.34e-19 | 2.96e-20 | 乱数シード依存の微小差 |

決定論的な数値は完全一致。ショットベースの微小差はシード依存の正常な挙動。

## 5. ノートブック全セル出力の検証

### Cell 2（パラメータ設定）
- `params`: 4分子、d=3、dim=81 ✅
- `params_boson`: 2分子、ボソンn_max=1、dim=36 ✅
- `validation: []`（パラメータ検証エラーなし）✅

### Cell 4（古典GKSL）
- N_S0=3.2762, N_T1=0.4184, N_S1=0.3054 ✅
- Tr保存: max|Tr-1| = 5.33e-15（機械精度）✅

### Cell 6（Qudit GKSL）
- System qutrits: 4, Ancilla qubits: 26 ✅
- Gates/step: 59 = 4(VirtRz) + 3(CustomTwo) + 26×2(Stinespring) ✅
- Tr保存: max|Tr-1| = 2.02e-14 ✅

### Cell 8（Qubit GKSL）
- System qubits: 8, Ancilla: 26, Total: 34 ✅
- Tr保存: max|Tr-1| = 2.04e-14 ✅

### Cell 10-14（ボソンシミュレータ）
- 全てTr保存が機械精度 ✅

### Cell 16（回路可視化）
- 全ゲートの量子回路図が正しく生成 ✅

### Cell 20（ユニタリ比較）
- Unitary max entropy = 2.22e-16 ≈ 0（純粋状態維持）✅
- GKSL final entropy = 1.7048（散逸による混合状態への遷移）✅

### Cell 22（保存則）
- 全シミュレータでTr保存 ✅
- 粒子数保存（Classical: 2.22e-14, Qubit: 8.35e-14, Qudit: 8.17e-14）✅

### Cell 24（忠実度）
- F(Classical vs Qudit) = 0.999999 ✅
- F(Classical vs Qubit) = 0.999999 ✅
- F(Qubit vs Qudit) = 1.000000（機械精度で一致）✅

### Cell 26（収束解析）
- Rate(T) → 1.0 確認（O(dt)収束）✅
- 回帰テスト結果と完全一致 ✅

### Cell 30-36（ショットベース）
- Qudit Shot: N_S0=3.24, N_T1=0.46, N_S1=0.30（1000 shots, DM結果に近い）✅
- Qubit Shot: DM結果と一致（同一シード、同一物理部分空間）✅
- Noisy: 重度デポーラリゼーションによる品質低下は物理的に正しい ✅
- Qubit Noisy: Tr(rho)=0.3006（禁止状態リーケージ 70%）は p_depol=0.01 × 5900ゲートで期待される ✅

### Cell 38（忠実度比較）
- F(DM vs Shot_no_noise) ≈ 0.998（1000ショットの統計的揺らぎ）✅
- F(DM vs Shot_noisy) ≈ 0.11（重度ノイズによる忠実度低下）✅
- 正規化忠実度 F_norm の計算が正しい ✅

## 6. 追加の確認事項

### 6.1 ヒューリスティック処理の不在
以下を確認:
- 密度行列の強制正規化は使用されていない ✅
- 固有値クリッピングは使用されていない ✅
- トレースの強制補正は使用されていない ✅
- 人工的なfallback処理は存在しない ✅

### 6.2 エラー処理
- `PhysicsViolationError`: 物理制約違反時に例外送出 ✅
- `NumericalInstabilityError`: 数値不安定時に例外送出 ✅
- `ValueError`: パラメータ不正時に例外送出 ✅

### 6.3 自己参照ループの解消状態
- ノートブック Cell 25/39: 「複数回の検証イテレーション」（固定表記、ハードコードなし）✅
- 回帰テスト: Check 7/8/31 除外済み ✅
- 再実行しても新しい「問題」を生成しない ✅

## 7. 最終判定

### 修正不要の根拠

1. **数学的正確性**: Stinespring dilation、GKSL超演算子、Trotter分割、Lindblad演算子構成、ベクトル化規約の全てが数学的に正確
2. **数値的正確性**: 回帰テスト31チェック全PASS、2回の独立実行で数値再現
3. **物理的整合性**: トレース保存、粒子数保存、正定値性、収束次数が全て理論的期待値と一致
4. **ノートブック品質**: 全セルの出力が正確、ドキュメントが正確、ハードコードのイテレーション回数なし
5. **コード品質**: ヒューリスティック処理なし、適切なエラー処理、一貫したベクトル化規約

### 今後の対応

- `run_tta_uc_gksl_regression_test.py` はイテレーション回数に依存しないため、コード変更時の回帰テストとして永続的に使用可能
- これ以上のイテレーションは不要
