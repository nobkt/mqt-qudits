# TTA-UC現象のGKSL-Lindblad量子ダイナミクス完全理論書

## 文書情報

**作成日**: 2026年1月25日
**バージョン**: 1.0.0
**対象フレームワーク**: MQT-Qudits / Qiskit
**理論的基礎**: Gorini-Kossakowski-Sudarshan-Lindblad (GKSL) 方程式
**適用系**: 三重項-三重項消滅アップコンバージョン (TTA-UC) 現象

---

## 目次

1. [はじめに](#1-はじめに)
2. [TTA-UC現象の物理的基礎](#2-tta-uc現象の物理的基礎)
3. [開放量子系の数学的基礎](#3-開放量子系の数学的基礎)
4. [GKSL-Lindblad方程式の完全定式化](#4-gksl-lindblad方程式の完全定式化)
5. [ボソン相互作用を含まないモデル](#5-ボソン相互作用を含まないモデル)
6. [ボソン相互作用を含むモデル](#6-ボソン相互作用を含むモデル)
7. [古典計算による数値シミュレーション手法](#7-古典計算による数値シミュレーション手法)
8. [Qubit量子アルゴリズムによる実装](#8-qubit量子アルゴリズムによる実装)
9. [Qudit量子アルゴリズムによる実装](#9-qudit量子アルゴリズムによる実装)
10. [6つの実装シナリオの詳細比較](#10-6つの実装シナリオの詳細比較)
11. [数値精度と物理的整合性の検証](#11-数値精度と物理的整合性の検証)
12. [結論](#12-結論)
13. [参考文献](#13-参考文献)

---

## 1. はじめに

### 1.1 本文書の目的

本文書は、**三重項-三重項消滅アップコンバージョン（Triplet-Triplet Annihilation Upconversion, TTA-UC）現象**のGKSL-Lindblad量子ダイナミクスを、以下の6つのシナリオにおいて理論的に完全かつ厳密に定式化することを目的とする：

1. **古典計算、ボソン相互作用無し**
2. **古典計算、ボソン相互作用有り**
3. **Qubit量子アルゴリズム、ボソン相互作用無し**
4. **Qubit量子アルゴリズム、ボソン相互作用有り**
5. **Qudit量子アルゴリズム、ボソン相互作用無し**
6. **Qudit量子アルゴリズム、ボソン相互作用有り**

### 1.2 厳密性の原則

本文書は以下の原則に厳格に従う：

**✅ 許容される手法**：

- GKSL定理に基づく数学的に厳密な導出
- Stinespring dilationによる非ユニタリ演算の量子回路表現
- 実験的に測定可能なパラメータのみの使用
- 制御可能な近似誤差（Trotter分解など）の明示的評価

**❌ 禁止される手法**：

- ヒューリスティックな近似や経験的フィッティング
- Fallback処理（計算失敗時の「適当な値」への置き換え）
- 物理的根拠のない簡略化
- ごまかしや真実を隠蔽する記述
- ユーザーへの迎合（真実ベースの厳密性を優先）

### 1.3 理論的位置づけ

TTA-UC現象は本質的に**開放量子系**（Open Quantum System）のダイナミクスである。分子系（系）は周囲の環境（フォノンバス、光子場など）と相互作用しており、この相互作用により系のダイナミクスは非ユニタリとなる。

**閉じた系（Closed System）の限界**：

- Schrödinger方程式 $i\hbar \frac{d|\psi\rangle}{dt} = \hat{H}|\psi\rangle$ はユニタリ時間発展のみを記述
- エネルギー散逸、デコヒーレンス、不可逆過程を表現できない
- TTA過程の不可逆性を正しく捉えられない

**開放量子系（Open Quantum System）の必要性**：

- 密度演算子 $\hat{\rho}$ による混合状態の記述
- GKSL-Lindblad方程式による非ユニタリ時間発展
- エントロピー増大（熱力学第二法則）との整合性
- 実験的速度定数との直接対応

---

## 2. TTA-UC現象の物理的基礎

### 2.1 TTA-UC過程の概要

#### 2.1.1 アップコンバージョンの基本原理

TTA-UC（三重項-三重項消滅アップコンバージョン）は、低エネルギーの光子を吸収して高エネルギーの光子を放出する光物理過程である。この過程は以下のステップで進行する：

1. **増感剤による光吸収と項間交差**:

   $$
   S_0^{\text{sensitizer}} + h\nu_{\text{low}} \xrightarrow{\text{absorption}} S_1^{\text{sensitizer}} \xrightarrow{\text{ISC}} T_1^{\text{sensitizer}}
   $$

2. **三重項エネルギー移動（TTET）**:

   $$
   T_1^{\text{sensitizer}} + S_0^{\text{annihilator}} \xrightarrow{\text{TTET}} S_0^{\text{sensitizer}} + T_1^{\text{annihilator}}
   $$

3. **三重項-三重項消滅（TTA）**:

   $$
   T_1^{\text{ann}} + T_1^{\text{ann}} \xrightarrow{\text{TTA}} S_1^{\text{ann}} + S_0^{\text{ann}}
   $$

4. **アップコンバージョン蛍光発光**:
   $$
   S_1^{\text{annihilator}} \xrightarrow{\text{fluorescence}} S_0^{\text{annihilator}} + h\nu_{\text{high}}
   $$

#### 2.1.2 エネルギー関係

TTA過程がエネルギー的に許容される条件：

$$
2E_{T_1} \geq E_{S_1}
$$

本モデルでの具体的数値：

- 三重項エネルギー：$E_{T_1} = 1.5$ eV
- 一重項エネルギー：$E_{S_1} = 3.0$ eV
- 基底状態エネルギー：$E_{S_0} = 0$ eV（基準）

したがって：

$$
2 \times 1.5 \text{ eV} = 3.0 \text{ eV} = E_{S_1}
$$

等号が成立し、エネルギー的に厳密に許容される。

### 2.2 分子電子状態の定義

各分子 $i$ は以下の3つの電子状態を持つ：

#### 2.2.1 基底一重項状態 $|S_0\rangle_i$

$$
|S_0\rangle_i
$$

- **エネルギー**: $E_{S_0} = 0$ eV
- **スピン多重度**: 1（singlet, $S = 0$）
- **電子配置**: 全電子がスピン対を形成
- **Qutrit表現**: $|0\rangle_i$
- **寿命**: 無限大（安定）

#### 2.2.2 励起三重項状態 $|T_1\rangle_i$

$$
|T_1\rangle_i
$$

- **エネルギー**: $E_{T_1} = 1.5$ eV
- **スピン多重度**: 3（triplet, $S = 1$）
- **電子配置**: 不対電子2個（同方向スピン）
- **Qutrit表現**: $|1\rangle_i$
- **寿命**: μs〜ms（スピン禁制遷移のため長寿命）

三重項状態の3つの副準位（$M_S = -1, 0, +1$）は、ゼロ磁場分裂により縮退が解けるが、本モデルでは単一の有効準位として扱う。

#### 2.2.3 励起一重項状態 $|S_1\rangle_i$

$$
|S_1\rangle_i
$$

- **エネルギー**: $E_{S_1} = 3.0$ eV
- **スピン多重度**: 1（singlet, $S = 0$）
- **電子配置**: 不対電子2個（逆方向スピン）
- **Qutrit表現**: $|2\rangle_i$
- **寿命**: ns（許容遷移のため短寿命）

### 2.3 関連する物理過程

TTA-UC現象に関連するすべての物理過程を列挙する：

#### 2.3.1 三重項-三重項消滅（TTA）

$$
|T_1\rangle_i |T_1\rangle_j \to |S_1\rangle_i |S_0\rangle_j \quad \text{or} \quad |S_0\rangle_i |S_1\rangle_j
$$

- **速度定数**: $\gamma_{\text{TTA}}$ または $k_{\text{TTA}}$
- **機構**: 電子交換相互作用（Dexter型）
- **スピン統計因子**: 1/9（単純モデル）または 1/5（より精密なモデル）

**重要**: この過程は不可逆であり、Lindblad演算子で記述される。

#### 2.3.2 蛍光発光（Fluorescence）

$$
|S_1\rangle_i \to |S_0\rangle_i + h\nu_{\text{fl}}
$$

- **速度定数**: $\Gamma_{\text{fl}}$
- **光子エネルギー**: $h\nu_{\text{fl}} = E_{S_1} - E_{S_0} = 3.0$ eV
- **波長**: $\lambda_{\text{fl}} = 1240/3.0 \approx 413$ nm（紫色）
- **典型的寿命**: 1-10 ns

#### 2.3.3 燐光発光（Phosphorescence）

$$
|T_1\rangle_i \to |S_0\rangle_i + h\nu_{\text{ph}}
$$

- **速度定数**: $\Gamma_{\text{ph}}$
- **光子エネルギー**: $h\nu_{\text{ph}} = E_{T_1} - E_{S_0} = 1.5$ eV
- **波長**: $\lambda_{\text{ph}} = 1240/1.5 \approx 827$ nm（近赤外）
- **典型的寿命**: μs〜s（スピン禁制）

#### 2.3.4 内部転換（Internal Conversion, IC）

$$
|S_1\rangle_i \to |S_0\rangle_i + \text{phonons}
$$

- **速度定数**: $k_{\text{IC}}$
- **エネルギーギャップ則**: $k_{\text{IC}} \propto \exp(-\gamma \Delta E / \hbar \omega_{\text{vib}})$

#### 2.3.5 項間交差（Intersystem Crossing, ISC）

**S₁ → T₁**:

$$
|S_1\rangle_i \to |T_1\rangle_i
$$

- **速度定数**: $k_{\text{ISC}}^{S \to T}$

**T₁ → S₀**:

$$
|T_1\rangle_i \to |S_0\rangle_i + \text{phonons}
$$

- **速度定数**: $k_{\text{ISC}}^{T \to S}$

#### 2.3.6 三重項エネルギー移動（Triplet Energy Transfer）

$$
|T_1\rangle_i |S_0\rangle_j \leftrightarrow |S_0\rangle_i |T_1\rangle_j
$$

- **結合定数**: $V_{ij}$
- **機構**: Dexter交換機構
- **距離依存性**: $V_{ij} \propto \exp(-2r_{ij}/L)$

---

## 3. 開放量子系の数学的基礎

### 3.1 密度演算子

#### 3.1.1 純粋状態と混合状態

**純粋状態**：

$$
\hat{\rho} = |\psi\rangle\langle\psi|
$$

性質：

- $\hat{\rho}^2 = \hat{\rho}$（べき等）
- $\text{Tr}[\hat{\rho}^2] = 1$
- $S(\hat{\rho}) = -\text{Tr}[\hat{\rho} \ln \hat{\rho}] = 0$（エントロピー）

**混合状態**：

$$
\hat{\rho} = \sum_k p_k |\psi_k\rangle\langle\psi_k|
$$

ここで、$p_k \geq 0$ かつ $\sum_k p_k = 1$。

性質：

- $\hat{\rho}^2 \neq \hat{\rho}$
- $\text{Tr}[\hat{\rho}^2] < 1$
- $S(\hat{\rho}) > 0$

#### 3.1.2 密度演算子の公理

任意の物理的に許容される密度演算子 $\hat{\rho}$ は以下を満たす：

1. **エルミート性**：$\hat{\rho} = \hat{\rho}^\dagger$
2. **正定値性**：$\langle\phi|\hat{\rho}|\phi\rangle \geq 0$ for all $|\phi\rangle$
3. **トレース1**：$\text{Tr}[\hat{\rho}] = 1$

### 3.2 系と環境の分離

#### 3.2.1 複合系の記述

系（System）$\mathcal{S}$ と環境（Environment）$\mathcal{E}$ の複合系：

$$
\mathcal{H}_{\text{total}} = \mathcal{H}_{\mathcal{S}} \otimes \mathcal{H}_{\mathcal{E}}
$$

全ハミルトニアン：

$$
\hat{H}_{\text{total}} = \hat{H}_{\mathcal{S}} \otimes \hat{I}_{\mathcal{E}} + \hat{I}_{\mathcal{S}} \otimes \hat{H}_{\mathcal{E}} + \hat{H}_{\text{int}}
$$

ここで：

- $\hat{H}_{\mathcal{S}}$：系のハミルトニアン
- $\hat{H}_{\mathcal{E}}$：環境のハミルトニアン
- $\hat{H}_{\text{int}}$：系-環境相互作用ハミルトニアン

#### 3.2.2 縮約密度演算子

系のみの情報を抽出するため、環境を部分トレースで消去：

$$
\hat{\rho}_{\mathcal{S}} = \text{Tr}_{\mathcal{E}}[\hat{\rho}_{\mathcal{SE}}]
$$

**重要**：$\hat{\rho}_{\mathcal{S}}$ の時間発展は一般に**非ユニタリ**である。

### 3.3 動力学写像

#### 3.3.1 CPTP写像の定義

系の時間発展を記述する写像 $\mathcal{E}_t$：

$$
\hat{\rho}_{\mathcal{S}}(t) = \mathcal{E}_t[\hat{\rho}_{\mathcal{S}}(0)]
$$

物理的に許容される写像は**完全正値トレース保存（CPTP）写像**でなければならない。

**CPTP条件**：

1. **線形性**：

   $$
   \mathcal{E}_t[a\hat{\rho}_1 + b\hat{\rho}_2] = a\mathcal{E}_t[\hat{\rho}_1] + b\mathcal{E}_t[\hat{\rho}_2]
   $$

2. **トレース保存**：

   $$
   \text{Tr}[\mathcal{E}_t[\hat{\rho}]] = \text{Tr}[\hat{\rho}] = 1
   $$

3. **完全正値性**：任意の補助系 $\mathcal{A}$ に対して
   $$
   (\mathcal{E}_t \otimes \mathcal{I}_{\mathcal{A}})[\hat{\rho}_{\mathcal{SA}}] \geq 0
   $$

#### 3.3.2 Kraus表現定理

**定理（Kraus, 1983）**：

任意のCPTP写像 $\mathcal{E}$ は以下の形式で表現できる：

$$
\mathcal{E}[\hat{\rho}] = \sum_{\alpha} \hat{K}_\alpha \hat{\rho} \hat{K}_\alpha^\dagger
$$

ここで、Kraus演算子 $\{\hat{K}_\alpha\}$ は完全性条件を満たす：

$$
\sum_{\alpha} \hat{K}_\alpha^\dagger \hat{K}_\alpha = \hat{I}
$$

### 3.4 マルコフ近似

#### 3.4.1 マルコフ過程の条件

以下の条件下でマルコフ近似が正当化される：

1. **弱結合極限**：系-環境結合が十分弱い
2. **時間スケール分離**：環境の相関時間 $\tau_{\mathcal{E}}$ が系の特徴的時間 $\tau_{\mathcal{S}}$ より十分短い
   $$
   \tau_{\mathcal{E}} \ll \tau_{\mathcal{S}}
   $$

#### 3.4.2 半群性

マルコフ近似が成立する場合、動力学写像は半群性を満たす：

$$
\mathcal{E}_{t+s} = \mathcal{E}_t \circ \mathcal{E}_s
$$

これにより、微分形式のマスター方程式が導出される。

---

## 4. GKSL-Lindblad方程式の完全定式化

### 4.1 GKSL定理

#### 4.1.1 定理の主張

**Gorini-Kossakowski-Sudarshan-Lindblad (GKSL) 定理（1976）**：

マルコフ的CPTP半群の最も一般的な生成子は、以下のLindblad方程式で与えられる：

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}, \hat{\rho}] + \sum_{\alpha} \gamma_\alpha \left( \hat{L}_\alpha \hat{\rho} \hat{L}_\alpha^\dagger - \frac{1}{2}\{\hat{L}_\alpha^\dagger \hat{L}_\alpha, \hat{\rho}\} \right)
$$

ここで：

- $\hat{H}$：系のハミルトニアン（エルミート）
- $\hat{L}_\alpha$：Lindblad演算子（ジャンプ演算子）
- $\gamma_\alpha > 0$：散逸速度定数
- $\{A, B\} = AB + BA$：反交換子

#### 4.1.2 Lindblad超演算子

散逸項を超演算子形式で書く：

$$
\mathcal{D}[\hat{L}][\hat{\rho}] = \hat{L}\hat{\rho}\hat{L}^\dagger - \frac{1}{2}\hat{L}^\dagger\hat{L}\hat{\rho} - \frac{1}{2}\hat{\rho}\hat{L}^\dagger\hat{L}
$$

これを用いると：

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}, \hat{\rho}] + \sum_{\alpha} \gamma_\alpha \mathcal{D}[\hat{L}_\alpha][\hat{\rho}]
$$

### 4.2 GKSL方程式の数学的性質

#### 4.2.1 トレース保存

$$
\frac{d}{dt}\text{Tr}[\hat{\rho}] = 0
$$

**証明**：

$$
\begin{align}
\frac{d}{dt}\text{Tr}[\hat{\rho}] &= -\frac{i}{\hbar}\text{Tr}[[\hat{H}, \hat{\rho}]] + \sum_{\alpha}\gamma_\alpha \text{Tr}\left[\hat{L}_\alpha\hat{\rho}\hat{L}_\alpha^\dagger - \frac{1}{2}\{\hat{L}_\alpha^\dagger\hat{L}_\alpha, \hat{\rho}\}\right] \\
&= 0 + \sum_{\alpha}\gamma_\alpha \left(\text{Tr}[\hat{L}_\alpha^\dagger\hat{L}_\alpha\hat{\rho}] - \text{Tr}[\hat{L}_\alpha^\dagger\hat{L}_\alpha\hat{\rho}]\right) \\
&= 0 \quad \square
\end{align}
$$

#### 4.2.2 完全正値性

GKSL形式は自動的にCPTP性を保証する。これは量子もつれ状態を含む任意の状態に対して物理的に正しい確率解釈が可能であることを意味する。

#### 4.2.3 エントロピー増大（熱力学第二法則）

von Neumannエントロピー：

$$
S(\hat{\rho}) = -\text{Tr}[\hat{\rho}\ln\hat{\rho}]
$$

GKSL方程式の下で：

$$
\frac{dS}{dt} \geq 0
$$

等号は平衡状態でのみ成立。

### 4.3 TTA-UC系のGKSL方程式

#### 4.3.1 完全なGKSL方程式

N分子系の完全なGKSL-Lindblad方程式：

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}_{\text{system}}, \hat{\rho}] + \mathcal{L}_{\text{total}}[\hat{\rho}]
$$

ここで：

$$
\hat{H}_{\text{system}} = \hat{H}_0 + \hat{H}_{\text{transfer}}
$$

$$
\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{TTA}} + \mathcal{L}_{\text{fl}} + \mathcal{L}_{\text{ph}} + \mathcal{L}_{\text{IC}} + \mathcal{L}_{\text{ISC}}^{S\to T} + \mathcal{L}_{\text{ISC}}^{T\to S}
$$

#### 4.3.2 ハミルトニアンの定義

**オンサイトエネルギー項**：

$$
\hat{H}_0 = \sum_{i=0}^{N-1} \left( E_{T_1} |T_1\rangle_i\langle T_1|_i + E_{S_1} |S_1\rangle_i\langle S_1|_i \right)
$$

Qutrit基底 $\{|0\rangle, |1\rangle, |2\rangle\} = \{|S_0\rangle, |T_1\rangle, |S_1\rangle\}$ では：

$$
\hat{H}_0 = \sum_{i=0}^{N-1} \left( E_T |1\rangle_i\langle 1|_i + E_S |2\rangle_i\langle 2|_i \right)
$$

**三重項エネルギー移動項**：

$$
\hat{H}_{\text{transfer}} = \sum_{\langle i,j \rangle} V_{ij} \left( |0\rangle_i\langle 1|_i \otimes |1\rangle_j\langle 0|_j + |1\rangle_i\langle 0|_i \otimes |0\rangle_j\langle 1|_j \right)
$$

物理的意味：

$$
|T_1\rangle_i |S_0\rangle_j \leftrightarrow |S_0\rangle_i |T_1\rangle_j
$$

#### 4.3.3 Lindblad演算子の完全リスト

**TTA過程**（隣接ペア $(i,j)$ に対して）：

$$
\hat{L}_{\text{TTA},1}^{(ij)} = \sqrt{\frac{\gamma_{\text{TTA}}}{2}} |2\rangle_i\langle 1|_i \otimes |0\rangle_j\langle 1|_j
$$

$$
\hat{L}_{\text{TTA},2}^{(ij)} = \sqrt{\frac{\gamma_{\text{TTA}}}{2}} |0\rangle_i\langle 1|_i \otimes |2\rangle_j\langle 1|_j
$$

**蛍光発光**（分子 $i$ に対して）：

$$
\hat{L}_{\text{fl}}^{(i)} = \sqrt{\Gamma_{\text{fl}}} |0\rangle_i\langle 2|_i
$$

**燐光発光**（分子 $i$ に対して）：

$$
\hat{L}_{\text{ph}}^{(i)} = \sqrt{\Gamma_{\text{ph}}} |0\rangle_i\langle 1|_i
$$

**内部転換**（分子 $i$ に対して）：

$$
\hat{L}_{\text{IC}}^{(i)} = \sqrt{k_{\text{IC}}} |0\rangle_i\langle 2|_i
$$

**項間交差 S₁→T₁**（分子 $i$ に対して）：

$$
\hat{L}_{\text{ISC}}^{S\to T,(i)} = \sqrt{k_{\text{ISC}}^{S\to T}} |1\rangle_i\langle 2|_i
$$

**項間交差 T₁→S₀**（分子 $i$ に対して）：

$$
\hat{L}_{\text{ISC}}^{T\to S,(i)} = \sqrt{k_{\text{ISC}}^{T\to S}} |0\rangle_i\langle 1|_i
$$

---

## 5. ボソン相互作用を含まないモデル

### 5.1 モデルの定義

ボソン相互作用を含まないモデルでは、フォノン（格子振動）や光子場との明示的な結合を省略し、それらの効果をLindblad散逸項の速度定数に繰り込む。

#### 5.1.1 有効モデルの妥当性

以下の条件下で有効モデルが正当化される：

1. **Born-Markov近似**：環境（フォノン、光子）との結合が弱く、環境の相関時間が系の動力学時間より十分短い

2. **時間スケール分離**：

   - フォノン緩和時間：$\tau_{\text{phonon}} \sim 10^{-13}$ s
   - 系の動力学時間：$\tau_{\text{system}} \sim 10^{-9} - 10^{-6}$ s
   - $\tau_{\text{phonon}} \ll \tau_{\text{system}}$ が成立

3. **詳細釣り合い条件**：有限温度効果を必要に応じて速度定数に含める

#### 5.1.2 系のヒルベルト空間

N分子系の状態空間：

$$
\mathcal{H}_{\text{system}} = \bigotimes_{i=0}^{N-1} \mathcal{H}^{(3)}_i
$$

ここで、$\mathcal{H}^{(3)}_i$ は分子 $i$ の3次元ヒルベルト空間。

状態空間の次元：

$$
\dim(\mathcal{H}_{\text{system}}) = 3^N
$$

### 5.2 完全なGKSL方程式（ボソン無し）

#### 5.2.1 マスター方程式

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}_0 + \hat{H}_{\text{transfer}}, \hat{\rho}] + \sum_{\alpha} \gamma_\alpha \mathcal{D}[\hat{L}_\alpha][\hat{\rho}]
$$

#### 5.2.2 ハミルトニアン項の行列表現

**単一分子のオンサイトハミルトニアン**：

Qutrit基底 $\{|0\rangle, |1\rangle, |2\rangle\}$ での行列表現：

$$
\hat{H}_0^{(i)} = \begin{pmatrix}
0 & 0 & 0 \\
0 & E_T & 0 \\
0 & 0 & E_S
\end{pmatrix} = E_T |1\rangle\langle 1| + E_S |2\rangle\langle 2|
$$

**2分子間のエネルギー移動ハミルトニアン**：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V_{ij} \left( |01\rangle\langle 10| + |10\rangle\langle 01| \right)_{ij}
$$

9×9行列（基底：$\{|00\rangle, |01\rangle, |02\rangle, |10\rangle, |11\rangle, |12\rangle, |20\rangle, |21\rangle, |22\rangle\}$）：

$$
\hat{H}_{\text{transfer}}^{(ij)} = V_{ij} \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

#### 5.2.3 Lindblad項の行列表現

**TTA Lindblad演算子**：

分子ペア $(i,j)$ に対する $\hat{L}_{\text{TTA},1}^{(ij)}$：

$$
\hat{L}_{\text{TTA},1}^{(ij)} = \sqrt{\frac{\gamma_{\text{TTA}}}{2}} |20\rangle\langle 11|_{ij}
$$

9×9行列表現：

$$
\hat{L}_{\text{TTA},1}^{(ij)} = \sqrt{\frac{\gamma_{\text{TTA}}}{2}} \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

**蛍光Lindblad演算子**：

$$
\hat{L}_{\text{fl}}^{(i)} = \sqrt{\Gamma_{\text{fl}}} |0\rangle\langle 2|_i = \sqrt{\Gamma_{\text{fl}}} \begin{pmatrix}
0 & 0 & 1 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}
$$

### 5.3 パラメータの物理的範囲

#### 5.3.1 典型的パラメータ値

| パラメータ         | 記号                      | 典型値              | 単位 | 物理的起源           |
| ------------------ | ------------------------- | ------------------- | ---- | -------------------- |
| 三重項エネルギー   | $E_T$                     | 1.5                 | eV   | 実験値               |
| 一重項エネルギー   | $E_S$                     | 3.0                 | eV   | 実験値               |
| エネルギー移動積分 | $V$                       | 0.01                | eV   | Dexter機構           |
| TTA速度定数        | $\gamma_{\text{TTA}}$     | $10^{-3}-10^{-1}$   | eV/ℏ | 拡散制御             |
| 蛍光速度           | $\Gamma_{\text{fl}}$      | $10^{-7}$           | eV/ℏ | 許容遷移             |
| 燐光速度           | $\Gamma_{\text{ph}}$      | $10^{-12}-10^{-9}$  | eV/ℏ | 禁制遷移             |
| IC速度             | $k_{\text{IC}}$           | $10^{-8}-10^{-7}$   | eV/ℏ | エネルギーギャップ則 |
| ISC S→T            | $k_{\text{ISC}}^{S\to T}$ | $10^{-8}-10^{-7}$   | eV/ℏ | スピン-軌道結合      |
| ISC T→S            | $k_{\text{ISC}}^{T\to S}$ | $10^{-13}-10^{-10}$ | eV/ℏ | 大ギャップ           |

#### 5.3.2 単位系

自然単位系（$\hbar = 1$）を使用する場合：

- エネルギー：eV
- 時間：ℏ/eV ≈ 0.658 fs
- 速度定数：eV/ℏ ≈ 1.52 fs⁻¹

SI単位系との対応：

- 1 eV/ℏ = 1.52 × 10¹⁵ s⁻¹
- 1 ns⁻¹ = 6.58 × 10⁻¹⁰ eV/ℏ

---

## 6. ボソン相互作用を含むモデル

### 6.1 ボソンモードの導入

#### 6.1.1 物理的動機

実際の分子系では、電子励起状態はフォノン（格子振動）および光子場と強く結合している。この相互作用を明示的に取り扱うことで：

1. 非マルコフ効果の記述が可能
2. 温度依存性の微視的理解
3. 振電結合の効果
4. 量子コヒーレンス効果

が得られる。

#### 6.1.2 ボソンヒルベルト空間

**フォノンモード**：

各分子に局所振動モード $\omega_k$ を導入：

$$
\mathcal{H}_{\text{phonon}} = \bigotimes_{k} \mathcal{H}_{\text{Fock},k}
$$

フォノンの生成・消滅演算子：

$$
[\hat{a}_k, \hat{a}_{k'}^\dagger] = \delta_{kk'}, \quad [\hat{a}_k, \hat{a}_{k'}] = 0
$$

数状態基底：

$$
|n_k\rangle, \quad n_k = 0, 1, 2, \ldots
$$

**光子モード**：

電磁場の量子化：

$$
\hat{\mathbf{E}}(\mathbf{r}) = \sum_{\mathbf{k},\lambda} \sqrt{\frac{\hbar\omega_k}{2\epsilon_0 V}} \left( \hat{b}_{\mathbf{k}\lambda} \mathbf{e}_{\mathbf{k}\lambda} e^{i\mathbf{k}\cdot\mathbf{r}} + \text{h.c.} \right)
$$

### 6.2 電子-フォノン結合

#### 6.2.1 Holstein型結合

電子状態とローカルフォノンモードの結合：

$$
\hat{H}_{e\text{-ph}} = \sum_{i,n} g_n^{(i)} (\hat{a}_n^{(i)} + \hat{a}_n^{(i)\dagger}) |n\rangle_i\langle n|_i
$$

ここで：

- $g_n^{(i)}$：状態 $|n\rangle$ のフォノン結合定数
- $\hat{a}_n^{(i)}$：分子 $i$ のフォノン消滅演算子

**Huang-Rhysパラメータ**：

$$
S_n = \frac{(g_n)^2}{(\hbar\omega)^2}
$$

これは各振動モードへの励起確率の平均値を表す。

#### 6.2.2 Peierls型結合

分子間ホッピング積分の格子振動への依存：

$$
\hat{H}_{\text{Peierls}} = \sum_{\langle i,j \rangle} V_{ij}(1 + \hat{u}_{ij}) \left( |01\rangle\langle 10| + |10\rangle\langle 01| \right)_{ij}
$$

ここで、$\hat{u}_{ij}$ は格子変位演算子：

$$
\hat{u}_{ij} = \sum_q \alpha_q (\hat{a}_q + \hat{a}_q^\dagger)
$$

### 6.3 電子-光子結合

#### 6.3.1 電気双極子相互作用

$$
\hat{H}_{e\text{-photon}} = -\sum_i \hat{\mathbf{d}}_i \cdot \hat{\mathbf{E}}(\mathbf{r}_i)
$$

ここで、$\hat{\mathbf{d}}_i$ は分子 $i$ の遷移双極子演算子：

$$
\hat{\mathbf{d}}_i = \sum_{m \neq n} \mathbf{d}_{mn}^{(i)} |m\rangle_i\langle n|_i
$$

#### 6.3.2 回転波近似（RWA）

高周波振動項を無視する回転波近似の下で：

$$
\hat{H}_{e\text{-photon}}^{\text{RWA}} = \sum_i \sum_{\mathbf{k},\lambda} \left( g_{\mathbf{k}\lambda}^{(i)} \hat{b}_{\mathbf{k}\lambda}^\dagger |0\rangle_i\langle 2|_i + \text{h.c.} \right)
$$

### 6.4 完全なハミルトニアン（ボソン有り）

#### 6.4.1 全ハミルトニアン

$$
\hat{H}_{\text{total}} = \hat{H}_{\text{el}} + \hat{H}_{\text{phonon}} + \hat{H}_{\text{photon}} + \hat{H}_{e\text{-ph}} + \hat{H}_{e\text{-photon}}
$$

各項の定義：

**電子系ハミルトニアン**：

$$
\hat{H}_{\text{el}} = \hat{H}_0 + \hat{H}_{\text{transfer}}
$$

**フォノンハミルトニアン**：

$$
\hat{H}_{\text{phonon}} = \sum_k \hbar\omega_k \hat{a}_k^\dagger \hat{a}_k
$$

**光子ハミルトニアン**：

$$
\hat{H}_{\text{photon}} = \sum_{\mathbf{k},\lambda} \hbar\omega_k \hat{b}_{\mathbf{k}\lambda}^\dagger \hat{b}_{\mathbf{k}\lambda}
$$

#### 6.4.2 状態空間の次元

完全な状態空間：

$$
\mathcal{H}_{\text{total}} = \mathcal{H}_{\text{el}} \otimes \mathcal{H}_{\text{phonon}} \otimes \mathcal{H}_{\text{photon}}
$$

フォノンと光子のFock空間を $n_{\max}$ で切断した場合：

$$
\dim(\mathcal{H}_{\text{total}}) = 3^N \times (n_{\max}+1)^{N_{\text{phonon}}} \times (n_{\max}+1)^{N_{\text{photon}}}
$$

**計算コストの爆発**：

- $N = 4$ 分子
- $N_{\text{phonon}} = 4$ モード
- $n_{\max} = 5$

の場合：$\dim = 81 \times 6^4 \times (\text{photon}) = 10^5$ 以上

### 6.5 ボソン浴の取り扱い

#### 6.5.1 スペクトル密度関数

フォノン浴のスペクトル密度：

$$
J(\omega) = \sum_k |g_k|^2 \delta(\omega - \omega_k)
$$

**Drude-Lorentzモデル**（多く使用される形式）：

$$
J(\omega) = \frac{2\lambda\omega_c\omega}{\omega^2 + \omega_c^2}
$$

ここで：

- $\lambda$：再配置エネルギー
- $\omega_c$：カットオフ周波数

**Ohmicモデル**：

$$
J(\omega) = \eta\omega e^{-\omega/\omega_c}
$$

#### 6.5.2 有限温度効果

温度 $T$ でのボソン占有数：

$$
\bar{n}(\omega) = \frac{1}{e^{\hbar\omega/(k_B T)} - 1}
$$

**詳細釣り合い条件**：

遷移速度に温度依存性が入る：

$$
\frac{\gamma_{n \to m}}{\gamma_{m \to n}} = e^{-(E_m - E_n)/(k_B T)}
$$

### 6.6 階層的運動方程式（HEOM）

#### 6.6.1 HEOMの基本形式

ボソン浴との相互作用を厳密に取り扱う非摂動的手法：

$$
\frac{\partial}{\partial t} \hat{\rho}_{\mathbf{n}}(t) = -\left(\frac{i}{\hbar}\hat{H}_{\text{el}}^{\times} + \sum_k n_k \nu_k\right) \hat{\rho}_{\mathbf{n}}(t)
$$

$$
+ \sum_k \hat{V}_k^{\times} \hat{\rho}_{\mathbf{n}+\mathbf{e}_k}(t) + \sum_k n_k \hat{\Phi}_k \hat{\rho}_{\mathbf{n}-\mathbf{e}_k}(t)
$$

ここで、$\mathbf{n} = (n_1, n_2, \ldots)$ は階層インデックス。

#### 6.6.2 計算コスト

階層数 $L$、各階層での次元数 $K$ とすると：

$$
\text{補助密度行列の数} = \binom{L + K - 1}{K - 1}
$$

$L = 10$, $K = 4$ の場合：286個の補助密度行列が必要。

**HEOMの利点**：

- 非マルコフ効果を厳密に取り扱い
- 任意の温度で適用可能
- 系統的な収束確認が可能

**HEOMの欠点**：

- 計算コストが高い
- 低温で収束が遅い
- 量子計算への直接適用が困難

---

## 7. 古典計算による数値シミュレーション手法

### 7.1 超演算子形式

#### 7.1.1 ベクトル化

密度行列 $\hat{\rho}$ を列ベクトルに変換：

$$
|\hat{\rho}\rangle\rangle = \text{vec}(\hat{\rho})
$$

$d \times d$ 行列を $d^2$ 次元ベクトルに変換する。

**変換規則**：

$$
(\hat{\rho})_{ij} \mapsto |\hat{\rho}\rangle\rangle_{i + j \cdot d}
$$

#### 7.1.2 超演算子行列

GKSL方程式を超演算子形式で書く：

$$
\frac{d|\hat{\rho}\rangle\rangle}{dt} = \mathcal{L}_{\text{super}} |\hat{\rho}\rangle\rangle
$$

超演算子 $\mathcal{L}_{\text{super}}$ は $d^2 \times d^2$ 行列。

**ハミルトニアン項**：

$$
\mathcal{L}_{\text{H}} = -\frac{i}{\hbar} (\hat{H} \otimes \hat{I} - \hat{I} \otimes \hat{H}^T)
$$

**Lindblad項**：

$$
\mathcal{L}_{\text{Lindblad},\alpha} = \gamma_\alpha \left( \hat{L}_\alpha \otimes \bar{\hat{L}}_\alpha - \frac{1}{2} \hat{L}_\alpha^\dagger \hat{L}_\alpha \otimes \hat{I} - \frac{1}{2} \hat{I} \otimes (\hat{L}_\alpha^\dagger \hat{L}_\alpha)^T \right)
$$

**全超演算子**：

$$
\mathcal{L}_{\text{super}} = \mathcal{L}_{\text{H}} + \sum_{\alpha} \mathcal{L}_{\text{Lindblad},\alpha}
$$

### 7.2 ボソン相互作用無しの古典計算

#### 7.2.1 行列指数関数法

形式解：

$$
|\hat{\rho}(t)\rangle\rangle = e^{\mathcal{L}_{\text{super}} t} |\hat{\rho}(0)\rangle\rangle
$$

数値計算：

```python
import scipy.linalg as la

rho_t = la.expm(L_super * t) @ rho_0_vec
```

計算量：$O(d^6)$（行列指数関数の計算）

**適用範囲**：$d = 3^N \leq 27$（N ≤ 3分子）

#### 7.2.2 常微分方程式ソルバー

より大きな系では、ODEソルバーを使用：

$$
\frac{d|\hat{\rho}\rangle\rangle}{dt} = f(t, |\hat{\rho}\rangle\rangle) = \mathcal{L}_{\text{super}} |\hat{\rho}\rangle\rangle
$$

**Runge-Kutta法（RK45）**：

```python
from scipy.integrate import solve_ivp

sol = solve_ivp(lambda t, y: L_super @ y, [0, t_final], rho_0_vec, method="RK45")
```

**BDF法（スティッフ系向け）**：

```python
sol = solve_ivp(lambda t, y: L_super @ y, [0, t_final], rho_0_vec, method="BDF")
```

#### 7.2.3 疎行列の活用

超演算子は高度に疎である：

$$
\text{非ゼロ要素数} \sim O(d^2)
$$

$$
\text{全要素数} = d^4
$$

スパース率：

$$
\text{sparsity} = 1 - \frac{O(d^2)}{d^4} \approx 1 - O(d^{-2})
$$

**疎行列演算の実装**：

```python
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

L_super_sparse = csr_matrix(L_super)
rho_t = expm_multiply(L_super_sparse * t, rho_0_vec)
```

計算量：$O(d^2 \cdot N_{\text{Krylov}})$

### 7.3 ボソン相互作用有りの古典計算

#### 7.3.1 テンソルネットワーク法

状態空間の指数的成長を抑制するため、テンソルネットワーク表現を使用：

**行列積状態（MPS）**：

$$
|\psi\rangle = \sum_{s_1, \ldots, s_N} A^{[1]}_{s_1} A^{[2]}_{s_2} \cdots A^{[N]}_{s_N} |s_1 s_2 \cdots s_N\rangle
$$

**行列積演算子（MPO）**：

$$
\hat{O} = \sum_{\{s\}, \{s'\}} W^{[1]}_{s_1,s'_1} W^{[2]}_{s_2,s'_2} \cdots W^{[N]}_{s_N,s'_N} |s_1 \cdots s_N\rangle\langle s'_1 \cdots s'_N|
$$

**TEBD (Time Evolving Block Decimation)**：

時間発展演算子をTrotter分解し、2サイトゲートとして適用：

$$
e^{-i\hat{H}\tau} \approx \prod_{\text{odd } j} e^{-i\hat{h}_{j,j+1}\tau} \prod_{\text{even } j} e^{-i\hat{h}_{j,j+1}\tau}
$$

#### 7.3.2 HEOM法の数値実装

階層的運動方程式の数値解法：

1. 階層インデックス $\mathbf{n}$ の列挙
2. 補助密度行列 $\hat{\rho}_{\mathbf{n}}$ の初期化
3. 時間発展の計算（RK法など）
4. 物理量の抽出（$\hat{\rho}_{\mathbf{0}}$から）

**実装上の注意**：

- 階層の打ち切り：$|\mathbf{n}| \leq L_{\max}$
- 収束確認：$L_{\max}$ を増やして結果が変化しないことを確認

#### 7.3.3 量子モンテカルロ法

確率的波動関数法：

1. 初期純粋状態 $|\psi(0)\rangle$ を準備
2. 確率的にジャンプを適用
3. 多数の軌跡を平均

**量子ジャンプ**：

- 確率 $dp_\alpha = \gamma_\alpha \langle\hat{L}_\alpha^\dagger \hat{L}_\alpha\rangle dt$ でジャンプ
- ジャンプ時：$|\psi\rangle \to \hat{L}_\alpha |\psi\rangle / \|\hat{L}_\alpha |\psi\rangle\|$
- ジャンプなし：非ユニタリ時間発展

**統計誤差**：$O(1/\sqrt{N_{\text{trajectories}}})$

---

## 8. Qubit量子アルゴリズムによる実装

### 8.1 Qubit表現による3準位系のエンコーディング

#### 8.1.1 2-Qubitエンコーディング

各分子 $i$ の3準位状態を2個のqubit $(q_{2i}, q_{2i+1})$ で表現：

$$
\begin{align}
|S_0\rangle_i &\longleftrightarrow |00\rangle_{2i,2i+1} \\
|T_1\rangle_i &\longleftrightarrow |01\rangle_{2i,2i+1} \\
|S_1\rangle_i &\longleftrightarrow |10\rangle_{2i,2i+1}
\end{align}
$$

**禁止状態**：

$$
|11\rangle_{2i,2i+1} \quad \text{（物理的意味なし）}
$$

N分子系のqubit数：$2N$

状態空間次元：

- 完全Qubit空間：$4^N$
- 物理的部分空間：$3^N$
- 非物理状態数：$4^N - 3^N$

#### 8.1.2 物理的部分空間への射影

物理的部分空間への射影演算子：

$$
\hat{P}_{\text{phys}}^{(i)} = |00\rangle\langle 00| + |01\rangle\langle 01| + |10\rangle\langle 10|
$$

全射影演算子：

$$
\hat{P}_{\text{phys}} = \bigotimes_{i=0}^{N-1} \hat{P}_{\text{phys}}^{(i)}
$$

### 8.2 Stinespring Dilationによる非ユニタリ演算の実装

#### 8.2.1 Stinespring表現定理

**定理**：任意のCPTP写像 $\mathcal{E}$ は、補助系を用いたユニタリ演算と部分トレースで表現できる：

$$
\mathcal{E}[\hat{\rho}_S] = \text{Tr}_E \left[ \hat{U}_{SE} (\hat{\rho}_S \otimes |0\rangle_E\langle 0|) \hat{U}_{SE}^\dagger \right]
$$

#### 8.2.2 Lindblad演算子のStinespring実装

単一Lindblad演算子 $\hat{L}$ に対するStinespring構成：

補助qubit $E$ を導入し、以下のユニタリを構築：

$$
\hat{U}_{\text{Lindblad}}(\theta) = \exp\left(-i\theta \hat{G}_{\text{Lindblad}}\right)
$$

ここで：

$$
\hat{G}_{\text{Lindblad}} = \frac{1}{\sqrt{2}} \left( \hat{L} \otimes \hat{\sigma}^-_E + \hat{L}^\dagger \otimes \hat{\sigma}^+_E \right)
$$

$\theta = \sqrt{\gamma\tau}$ で時間発展パラメータを設定。

### 8.3 ボソン相互作用無しのQubit実装

#### 8.3.1 ユニタリ部分の量子回路

**オンサイトエネルギー項**：

$\hat{H}_0$ は対角なので、$Z$回転ゲートで実装：

$$
e^{-i\hat{H}_0\tau/\hbar} = \prod_i e^{-i(E_T |01\rangle\langle 01| + E_S |10\rangle\langle 10|)_i \tau/\hbar}
$$

各分子に対して：

$$
\text{RZ}_{q_{2i}}(-E_S\tau/\hbar) \cdot \text{RZ}_{q_{2i+1}}(-E_T\tau/\hbar) \cdot \text{制御位相ゲート}
$$

**エネルギー移動項**：

$\hat{H}_{\text{transfer}}$ は4-qubit演算（2分子 = 4qubit）：

$$
e^{-iV\tau(\hat{\sigma}^+_i\hat{\sigma}^-_j + \text{h.c.})/\hbar}
$$

これをCNOTと単一qubit回転に分解。

#### 8.3.2 Lindblad項の量子回路

**蛍光Lindblad演算子のStinespring実装**：

$$
\hat{L}_{\text{fl}} = \sqrt{\Gamma_{\text{fl}}} |00\rangle\langle 10|
$$

Qubit表現で：

- 系qubit：$(q_0, q_1)$
- 補助qubit：$q_E$

制御回転ゲートで実装：

$$
\text{C-RY}(2\arcsin(\sqrt{\Gamma_{\text{fl}}\tau}))
$$

制御：$q_0 = 1, q_1 = 0$（状態 $|10\rangle$）
ターゲット：$q_E$

その後、$q_E$ を測定または破棄。

#### 8.3.3 Trotter分解

完全な時間発展を時間刻み $\tau$ で分割：

**2次Trotter分解**：

$$
e^{\mathcal{L}_{\text{GKSL}}\tau} \approx e^{\mathcal{L}_H \tau/2} \cdot e^{\mathcal{L}_{\text{diss}}\tau} \cdot e^{\mathcal{L}_H \tau/2}
$$

誤差：$O(\tau^3)$ per step、全体で $O(t^3/N^2)$

### 8.4 ボソン相互作用有りのQubit実装

#### 8.4.1 ボソンモードのQubit表現

フォノン/光子モードをqubitでエンコード：

**バイナリエンコーディング**：

$$
|n\rangle_{\text{Fock}} \longleftrightarrow |n_1 n_2 \cdots n_k\rangle_{\text{qubit}}
$$

$n_{\max}$ まで表現するのに必要なqubit数：$\lceil \log_2(n_{\max}+1) \rceil$

**ユナリエンコーディング**：

$$
|n\rangle_{\text{Fock}} \longleftrightarrow |\underbrace{1\cdots1}_{n}\underbrace{0\cdots0}_{n_{\max}-n}\rangle
$$

必要qubit数：$n_{\max}$

#### 8.4.2 電子-フォノン結合の量子回路

Holstein型結合：

$$
\hat{H}_{e\text{-ph}} = g(\hat{a} + \hat{a}^\dagger)|1\rangle\langle 1|
$$

量子回路実装：

1. 電子状態 $|1\rangle$ を制御として使用
2. フォノンモードに変位演算子を適用

$$
\hat{D}(\alpha) = e^{\alpha\hat{a}^\dagger - \alpha^*\hat{a}}
$$

ただし、変位演算子のqubit分解は非自明であり、近似的実装が必要。

#### 8.4.3 リソース見積もり

N分子系でボソンモードを含む場合：

| リソース                 | ボソン無し            | ボソン有り                               |
| ------------------------ | --------------------- | ---------------------------------------- |
| 系qubit                  | $2N$                  | $2N$                                     |
| フォノンqubit            | 0                     | $N \cdot \lceil\log_2(n_{\max}+1)\rceil$ |
| 補助qubit（Stinespring） | $\sim 5N$             | $\sim 5N$                                |
| 全qubit数                | $\sim 7N$             | $\sim 7N + N\log_2 n_{\max}$             |
| ゲート深さ               | $O(N^2 \cdot t/\tau)$ | $O(N^2 n_{\max} \cdot t/\tau)$           |

---

## 9. Qudit量子アルゴリズムによる実装

### 9.1 Qutrit（3準位系）による自然な表現

#### 9.1.1 Qutrit基底

各分子の状態を直接Qutritで表現：

$$
|S_0\rangle_i \longleftrightarrow |0\rangle_i, \quad |T_1\rangle_i \longleftrightarrow |1\rangle_i, \quad |S_1\rangle_i \longleftrightarrow |2\rangle_i
$$

N分子系のqutrit数：$N$（qubitの半分）

状態空間次元：$3^N$（禁止状態なし）

#### 9.1.2 Quditの利点

1. **次元の一致**：物理的部分空間と計算空間が完全一致
2. **qudit数削減**：$N$ qutrit vs $2N$ qubit
3. **禁止状態回避**：非物理状態への遷移リスクがない
4. **自然な演算子表現**：Lindblad演算子が直接的に表現可能

### 9.2 MQT-Quditsの基本ゲートセット

#### 9.2.1 単一Quditゲート

**一般化Pauli-X演算子**：

$$
\hat{X}_d = \sum_{j=0}^{d-1} |j+1 \mod d\rangle\langle j|
$$

**一般化Pauli-Z演算子**：

$$
\hat{Z}_d = \sum_{j=0}^{d-1} \omega^j |j\rangle\langle j|, \quad \omega = e^{2\pi i/d}
$$

**部分空間回転**：

$$
\hat{R}_{mn}(\theta, \phi) = e^{-i\theta(\cos\phi \hat{\sigma}_{mn}^x + \sin\phi \hat{\sigma}_{mn}^y)/2}
$$

ここで、$\hat{\sigma}_{mn}^x = |m\rangle\langle n| + |n\rangle\langle m|$、$\hat{\sigma}_{mn}^y = -i|m\rangle\langle n| + i|n\rangle\langle m|$

#### 9.2.2 2-Quditゲート

**制御ゲート**：

$$
\text{C-}U = |0\rangle\langle 0| \otimes \hat{I} + |1\rangle\langle 1| \otimes \hat{U}_{01} + |2\rangle\langle 2| \otimes \hat{U}_{02}
$$

**SWAPゲート**：

$$
\text{SWAP} = \sum_{j,k} |jk\rangle\langle kj|
$$

### 9.3 ボソン相互作用無しのQudit実装

#### 9.3.1 ハミルトニアンの量子回路化

**オンサイトエネルギー項**：

$$
e^{-i\hat{H}_0^{(i)}\tau/\hbar} = \text{diag}(1, e^{-iE_T\tau/\hbar}, e^{-iE_S\tau/\hbar})
$$

単一qutritの対角ゲートとして直接実装。

**エネルギー移動項**：

2-qutrit演算子：

$$
e^{-iV\tau(|01\rangle\langle 10| + |10\rangle\langle 01|)/\hbar}
$$

$9 \times 9$ ユニタリ行列を基本ゲートに分解。

#### 9.3.2 Lindblad演算子の量子回路化

**TTA Lindblad演算子のStinespring実装**：

$$
\hat{L}_{\text{TTA},1}^{(ij)} = \sqrt{\frac{\gamma_{\text{TTA}}}{2}} |2\rangle_i\langle 1| \otimes |0\rangle_j\langle 1|
$$

Stinespring構成：

- 系：2 qutrit（分子 $i, j$）
- 環境：1 qubit（または1 qutrit）

ユニタリ演算子：

$$
\hat{U} = \exp\left(-i\sqrt{\gamma_{\text{TTA}}\tau/2} \cdot \hat{G}\right)
$$

$$
\hat{G} = |20\rangle\langle 11| \otimes |1\rangle_E\langle 0| + |11\rangle\langle 20| \otimes |0\rangle_E\langle 1|
$$

**蛍光Lindblad演算子**：

$$
\hat{L}_{\text{fl}}^{(i)} = \sqrt{\Gamma_{\text{fl}}} |0\rangle\langle 2|_i
$$

単一qutritと補助qubitの2体ゲート：

$$
\hat{U}_{\text{fl}} = \exp\left(-i\sqrt{\Gamma_{\text{fl}}\tau} (|0\rangle\langle 2| \otimes |1\rangle_E\langle 0| + \text{h.c.})\right)
$$

#### 9.3.3 ゲート分解

MQT-Quditsフレームワークでの分解：

1. **Givens回転**による部分空間回転の連鎖
2. **対角ゲート**による位相の付加
3. **制御ゲート**による条件付き演算

Givens回転 $G_{mn}(\theta, \phi)$：

$$
G_{mn}(\theta, \phi) = I + (\cos\theta - 1)(|m\rangle\langle m| + |n\rangle\langle n|) + \sin\theta(e^{i\phi}|m\rangle\langle n| - e^{-i\phi}|n\rangle\langle m|)
$$

### 9.4 ボソン相互作用有りのQudit実装

#### 9.4.1 ボソンモードのQudit表現

フォノン/光子Fock空間を高次元quditで表現：

$$
|n\rangle_{\text{Fock}} \longleftrightarrow |n\rangle_{d_{\text{phonon}}}
$$

$d_{\text{phonon}} = n_{\max} + 1$ 次元qudit

**利点**：バイナリエンコーディング不要、自然なFock空間表現

#### 9.4.2 電子-フォノン結合の量子回路

Holstein型結合：

$$
\hat{H}_{e\text{-ph}} = g(\hat{a} + \hat{a}^\dagger) |1\rangle\langle 1|_{\text{el}}
$$

ここで：

$$
\hat{a} = \sum_{n=0}^{d-2} \sqrt{n+1} |n\rangle\langle n+1|
$$

制御付き昇降演算子として実装：

$$
\text{C-}a = |0\rangle\langle 0|_{\text{el}} \otimes \hat{I}_{\text{ph}} + |1\rangle\langle 1|_{\text{el}} \otimes \hat{a}_{\text{ph}} + |2\rangle\langle 2|_{\text{el}} \otimes \hat{I}_{\text{ph}}
$$

#### 9.4.3 リソース見積もり

| リソース      | ボソン無し            | ボソン有り                     |
| ------------- | --------------------- | ------------------------------ |
| 電子qutrit    | $N$                   | $N$                            |
| フォノンqudit | 0                     | $N$ (次元 $n_{\max}+1$)        |
| 補助qubit     | $\sim 3N$             | $\sim 3N$                      |
| 全qudit等価数 | $N + 2N = 3N$         | $2N + 2N = 4N$                 |
| ゲート深さ    | $O(N^2 \cdot t/\tau)$ | $O(N^2 n_{\max} \cdot t/\tau)$ |

---

## 10. 6つの実装シナリオの詳細比較

本章では、問題文で要求された6つのシナリオを詳細に比較する。

### 10.1 シナリオ一覧

| #   | 計算方式      | ボソン相互作用 | 略称         |
| --- | ------------- | -------------- | ------------ |
| 1   | 古典計算      | 無し           | Classical-NB |
| 2   | 古典計算      | 有り           | Classical-B  |
| 3   | Qubit量子計算 | 無し           | Qubit-NB     |
| 4   | Qubit量子計算 | 有り           | Qubit-B      |
| 5   | Qudit量子計算 | 無し           | Qudit-NB     |
| 6   | Qudit量子計算 | 有り           | Qudit-B      |

### 10.2 シナリオ1：古典計算、ボソン相互作用無し（Classical-NB）

#### 10.2.1 理論的定式化

**GKSL方程式**：

$$
\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}_0 + \hat{H}_{\text{transfer}}, \hat{\rho}] + \sum_{\alpha} \gamma_\alpha \mathcal{D}[\hat{L}_\alpha][\hat{\rho}]
$$

**状態空間**：$3^N$ 次元

**超演算子形式**：

$$
\frac{d|\hat{\rho}\rangle\rangle}{dt} = \mathcal{L}_{\text{super}} |\hat{\rho}\rangle\rangle
$$

$\mathcal{L}_{\text{super}}$ は $(3^N)^2 \times (3^N)^2$ 行列

#### 10.2.2 数値手法

1. **行列指数関数**（小規模系 $N \leq 3$）：

   $$
   |\hat{\rho}(t)\rangle\rangle = e^{\mathcal{L}_{\text{super}} t} |\hat{\rho}(0)\rangle\rangle
   $$

2. **ODEソルバー**（中規模系 $N \leq 5$）：

   - RK45法
   - BDF法（スティッフ系）

3. **疎行列＋Krylov法**（大規模系 $N \geq 6$）：
   $$
   e^{\mathcal{L}t} |\rho\rangle \approx \sum_{k=0}^{K} \frac{t^k}{k!} \mathcal{L}^k |\rho\rangle
   $$

#### 10.2.3 計算コスト

| 分子数 $N$ | 状態次元 $3^N$ | 密度行列サイズ | メモリ（double complex） |
| ---------- | -------------- | -------------- | ------------------------ |
| 2          | 9              | 81             | 1.3 KB                   |
| 3          | 27             | 729            | 12 KB                    |
| 4          | 81             | 6,561          | 105 KB                   |
| 5          | 243            | 59,049         | 0.9 MB                   |
| 6          | 729            | 531,441        | 8.5 MB                   |

計算時間（目安）：

- $N = 4$：秒オーダー
- $N = 5$：分オーダー
- $N = 6$：時間オーダー

### 10.3 シナリオ2：古典計算、ボソン相互作用有り（Classical-B）

#### 10.3.1 理論的定式化

**拡張ハミルトニアン**：

$$
\hat{H}_{\text{total}} = \hat{H}_{\text{el}} + \hat{H}_{\text{phonon}} + \hat{H}_{e\text{-ph}}
$$

**完全な状態空間**：

$$
\mathcal{H} = \mathcal{H}_{\text{el}} \otimes \mathcal{H}_{\text{phonon}}
$$

次元：$3^N \times (n_{\max}+1)^{N_{\text{ph}}}$

#### 10.3.2 数値手法

1. **HEOM（階層的運動方程式）**：

   - 非マルコフ効果を厳密に取り扱い
   - 階層数 $L$ で精度を制御
   - 計算コスト：$O(\binom{L+K}{K} \cdot d^4)$

2. **テンソルネットワーク法**：

   - MPS/MPO表現
   - 結合次元 $\chi$ で精度を制御
   - 計算コスト：$O(N \cdot d \cdot \chi^3)$

3. **量子モンテカルロ法**：
   - 確率的波動関数法
   - 軌跡数 $N_{\text{traj}}$ で精度を制御
   - 計算コスト：$O(N_{\text{traj}} \cdot d \cdot t/\tau)$

#### 10.3.3 計算コスト比較

| 手法                 | 精度               | メモリ                        | 計算時間                        |
| -------------------- | ------------------ | ----------------------------- | ------------------------------- |
| HEOM                 | 厳密               | $O(n_{\text{aux}} \cdot d^2)$ | $O(n_{\text{aux}}^2 \cdot d^2)$ |
| テンソルネットワーク | 近似（$\chi$制御） | $O(N \cdot d \cdot \chi^2)$   | $O(N \cdot d \cdot \chi^3)$     |
| 量子モンテカルロ     | 統計誤差           | $O(d)$                        | $O(N_{\text{traj}} \cdot d)$    |

### 10.4 シナリオ3：Qubit量子計算、ボソン相互作用無し（Qubit-NB）

#### 10.4.1 理論的定式化

**Qubitエンコーディング**：

- 分子状態：2 qubit/分子
- 状態空間：$4^N$（物理部分空間：$3^N$）

**Stinespring実装**：

- 各Lindblad演算子に補助qubitを導入
- ユニタリ演算後に部分トレース

#### 10.4.2 量子回路構成

**時間ステップあたりの回路**：

1. **ユニタリ部分**（$\hat{H}_0 + \hat{H}_{\text{transfer}}$）：

   - 対角ゲート：$O(N)$ 深さ
   - エネルギー移動：$O(N)$ 深さ

2. **散逸部分**（Stinespring）：
   - 各Lindblad項：$O(1)$ 深さ
   - 全体：$O(N + N(N-1)/2)$ 深さ（TTAペア数）

**総ゲート深さ**：

$$
D = O\left(\frac{t}{\tau} \cdot (N + N^2)\right) = O\left(\frac{t \cdot N^2}{\tau}\right)
$$

#### 10.4.3 リソース見積もり

| 分子数 $N$ | 系qubit | 補助qubit | 全qubit | ゲート数/ステップ |
| ---------- | ------- | --------- | ------- | ----------------- |
| 2          | 4       | ~6        | ~10     | ~50               |
| 3          | 6       | ~10       | ~16     | ~100              |
| 4          | 8       | ~15       | ~23     | ~200              |

### 10.5 シナリオ4：Qubit量子計算、ボソン相互作用有り（Qubit-B）

#### 10.5.1 理論的定式化

**拡張エンコーディング**：

- 電子状態：2 qubit/分子
- フォノンモード：$\lceil\log_2(n_{\max}+1)\rceil$ qubit/モード

**状態空間**：

$$
2^{2N} \times 2^{N \cdot \lceil\log_2(n_{\max}+1)\rceil}
$$

#### 10.5.2 量子回路構成

**電子-フォノン結合ゲート**：

$$
e^{-ig\tau(\hat{a}+\hat{a}^\dagger)|T_1\rangle\langle T_1|}
$$

qubitでの実装：

1. フォノンモードの加算/減算回路
2. 電子状態による制御

ゲート深さ：$O(n_{\max})$ per coupling

#### 10.5.3 リソース見積もり

| 分子数 $N$ | $n_{\max}$ | 電子qubit | フォノンqubit | 全qubit |
| ---------- | ---------- | --------- | ------------- | ------- |
| 2          | 4          | 4         | 6             | ~16     |
| 3          | 4          | 6         | 9             | ~25     |
| 4          | 4          | 8         | 12            | ~35     |

### 10.6 シナリオ5：Qudit量子計算、ボソン相互作用無し（Qudit-NB）

#### 10.6.1 理論的定式化

**Qutritエンコーディング**：

- 分子状態：1 qutrit/分子
- 状態空間：$3^N$（完全一致）

**Stinespring実装**：

- 補助qubit（または補助qutrit）を導入
- より自然なLindblad演算子表現

#### 10.6.2 量子回路構成

**単一qutritゲート**：

- Givens回転 $G_{01}(\theta)$, $G_{12}(\theta)$, $G_{02}(\theta)$
- 対角位相ゲート

**2-qutritゲート**：

- 制御Givens回転
- SWAP風ゲート

#### 10.6.3 リソース見積もり

| 分子数 $N$ | 系qutrit | 補助qudit | 等価qubit数 |
| ---------- | -------- | --------- | ----------- |
| 2          | 2        | ~3        | ~8          |
| 3          | 3        | ~5        | ~13         |
| 4          | 4        | ~7        | ~18         |

**Qubit比較**：

- Qubit-NB（$N=4$）：~23 qubit
- Qudit-NB（$N=4$）：~18 qubit等価（約22%削減）

### 10.7 シナリオ6：Qudit量子計算、ボソン相互作用有り（Qudit-B）

#### 10.7.1 理論的定式化

**拡張エンコーディング**：

- 電子状態：1 qutrit/分子（$d=3$）
- フォノンモード：1 qudit/モード（$d=n_{\max}+1$）

**状態空間**：$3^N \times (n_{\max}+1)^{N_{\text{ph}}}$

#### 10.7.2 量子回路構成

**電子-フォノン結合**：

より自然な実装が可能：

$$
\hat{a} = \sum_{n=0}^{d-2} \sqrt{n+1}|n\rangle\langle n+1|
$$

はqudit上で直接的に表現される。

#### 10.7.3 リソース見積もり

| 分子数 $N$ | $n_{\max}$ | 電子qutrit | フォノンqudit | 等価qubit数 |
| ---------- | ---------- | ---------- | ------------- | ----------- |
| 2          | 4          | 2          | 2 (d=5)       | ~11         |
| 3          | 4          | 3          | 3 (d=5)       | ~17         |
| 4          | 4          | 4          | 4 (d=5)       | ~23         |

### 10.8 6シナリオの総合比較

#### 10.8.1 計算資源比較（N=4分子、t=100 fs）

| シナリオ     | 計算資源        | 実行時間目安 | 精度        |
| ------------ | --------------- | ------------ | ----------- |
| Classical-NB | RAM 100KB       | 秒           | 機械精度    |
| Classical-B  | RAM 10MB (HEOM) | 分〜時       | 収束次第    |
| Qubit-NB     | 23 qubit        | 量子時間     | Trotter誤差 |
| Qubit-B      | 35 qubit        | 量子時間     | Trotter誤差 |
| Qudit-NB     | 18 qubit等価    | 量子時間     | Trotter誤差 |
| Qudit-B      | 23 qubit等価    | 量子時間     | Trotter誤差 |

#### 10.8.2 スケーラビリティ比較

| シナリオ     | $N$ 依存性                         | 実用限界             |
| ------------ | ---------------------------------- | -------------------- |
| Classical-NB | $O(9^N)$                           | $N \sim 6$           |
| Classical-B  | $O(9^N \cdot d_{\text{ph}}^N)$     | $N \sim 4$           |
| Qubit-NB     | $O(N^2)$ qubit                     | $N \sim 100+$ (NISQ) |
| Qubit-B      | $O(N^2 \cdot \log n_{\max})$ qubit | $N \sim 50+$ (NISQ)  |
| Qudit-NB     | $O(N)$ qutrit                      | $N \sim 100+$        |
| Qudit-B      | $O(N)$ qudit                       | $N \sim 50+$         |

---

## 11. 数値精度と物理的整合性の検証

### 11.1 保存量の検証

#### 11.1.1 トレース保存

$$
\text{Tr}[\hat{\rho}(t)] = 1 \pm \epsilon_{\text{trace}}
$$

許容誤差：$\epsilon_{\text{trace}} < 10^{-10}$

#### 11.1.2 正定値性

密度行列の全固有値：

$$
\lambda_k \geq -\epsilon_{\text{pos}}, \quad \forall k
$$

許容誤差：$\epsilon_{\text{pos}} < 10^{-10}$

#### 11.1.3 エルミート性

$$
\|\hat{\rho} - \hat{\rho}^\dagger\|_F < \epsilon_{\text{herm}}
$$

許容誤差：$\epsilon_{\text{herm}} < 10^{-10}$

### 11.2 物理法則の検証

#### 11.2.1 熱力学第二法則

von Neumannエントロピー：

$$
S(t) = -\text{Tr}[\hat{\rho}(t)\ln\hat{\rho}(t)]
$$

検証条件：

$$
\frac{dS}{dt} \geq -\epsilon_{\text{entropy}}
$$

（数値誤差による微小な負値は許容）

#### 11.2.2 粒子数保存

全分子数：

$$
N_{\text{total}}(t) = \sum_{i,n} n_i(t) = N
$$

ここで、$n_i(t) = \text{Tr}[\hat{P}_i \hat{\rho}(t)]$

### 11.3 Trotter誤差の評価

#### 11.3.1 理論的誤差限界

**1次Trotter**：

$$
\|e^{(\hat{A}+\hat{B})t} - (e^{\hat{A}t/N}e^{\hat{B}t/N})^N\| \leq \frac{t^2}{2N}\|[\hat{A},\hat{B}]\|
$$

**2次Trotter**：

$$
\|e^{(\hat{A}+\hat{B})t} - (e^{\hat{A}t/2N}e^{\hat{B}t/N}e^{\hat{A}t/2N})^N\| \leq \frac{t^3}{12N^2}\|[[\hat{A},\hat{B}],\hat{A}+\hat{B}]\|
$$

#### 11.3.2 実用的誤差見積もり

典型的パラメータ（$E_T = 1.5$ eV, $V = 0.01$ eV）：

$$
\|[\hat{H}_0, \hat{H}_{\text{transfer}}]\| \sim V \cdot (E_T - E_{S_0}) = 0.015 \text{ eV}^2
$$

時間 $t = 100$ fs、ステップ数 $N = 100$ の場合：

$$
\epsilon_{\text{Trotter}} \sim \frac{(100 \times 0.658)^2}{2 \times 100} \times 0.015 \sim 0.03
$$

より高精度が必要な場合は $N$ を増やす。

### 11.4 量子計算特有の検証

#### 11.4.1 禁止状態への遷移（Qubit実装）

$$
P_{\text{forbidden}}(t) = \text{Tr}[\hat{P}_{\text{non-phys}}\hat{\rho}(t)] < \epsilon_{\text{forbidden}}
$$

許容誤差：$\epsilon_{\text{forbidden}} < 10^{-8}$

#### 11.4.2 Stinespring実装の忠実度

理想的Lindblad時間発展と比較：

$$
F = \text{Tr}\sqrt{\sqrt{\hat{\rho}_{\text{ideal}}}\hat{\rho}_{\text{Stinespring}}\sqrt{\hat{\rho}_{\text{ideal}}}}^2
$$

目標：$F > 0.99$

---

## 12. 結論

### 12.1 本文書の成果

本文書では、**TTA-UC（三重項-三重項消滅アップコンバージョン）現象のGKSL-Lindblad量子ダイナミクス**を、6つの異なる実装シナリオにおいて完全かつ厳密に定式化した。

#### 12.1.1 理論的貢献

1. **GKSL-Lindblad方程式の完全定式化**

   - 開放量子系理論に基づく厳密な数学的基礎
   - TTA過程のLindblad演算子表現
   - 放射減衰・無放射遷移の完全な組み込み

2. **ボソン相互作用の系統的取り扱い**

   - 電子-フォノン結合（Holstein型、Peierls型）
   - 電子-光子結合（電気双極子相互作用）
   - 有限温度効果と詳細釣り合い

3. **量子計算への完全な理論的橋渡し**
   - Stinespring dilationによる非ユニタリ演算のユニタリ化
   - Qubit/Quditエンコーディングの詳細
   - 量子回路構成の具体的手順

#### 12.1.2 実装シナリオの比較

| シナリオ     | 適用領域       | 主な利点         |
| ------------ | -------------- | ---------------- |
| Classical-NB | 小規模系       | 高精度、低コスト |
| Classical-B  | 非マルコフ効果 | 物理的厳密性     |
| Qubit-NB     | NISQ量子計算   | 既存ハードウェア |
| Qubit-B      | 将来の量子計算 | 完全な物理モデル |
| Qudit-NB     | 次世代量子計算 | 資源効率         |
| Qudit-B      | 将来の量子計算 | 最も自然な表現   |

### 12.2 ヒューリスティック手法の完全排除

本文書は以下の原則を厳格に遵守した：

**❌ 排除された手法**：

- 物理的根拠のない近似パラメータ
- 計算失敗時のFallback処理
- 非物理的な状態への遷移
- ごまかしや隠蔽

**✅ 採用された手法**：

- GKSL定理に基づく厳密な理論
- Stinespring表現による厳密なユニタリ化
- 制御可能な近似誤差（Trotter分解）
- 物理法則（熱力学第二法則）との整合性

### 12.3 今後の展望

1. **理論的拡張**

   - 非マルコフ効果の量子計算実装
   - より大規模な分子系への適用
   - 空間的不均一性の導入

2. **実験的検証**

   - 実際の量子ハードウェアでの実装
   - 古典計算結果との比較
   - 新材料設計への応用

3. **アルゴリズム最適化**
   - 変分量子固有値法（VQE）との統合
   - 量子位相推定の適用
   - ノイズ耐性の向上

---

## 13. 参考文献

### 開放量子系理論

1. Breuer, H.-P., & Petruccione, F. (2002). _The Theory of Open Quantum Systems_. Oxford University Press.
2. Gorini, V., Kossakowski, A., & Sudarshan, E. C. G. (1976). "Completely positive dynamical semigroups of N-level systems." _J. Math. Phys._, 17, 821.
3. Lindblad, G. (1976). "On the generators of quantum dynamical semigroups." _Commun. Math. Phys._, 48, 119.

### TTA-UCとフォトフィジクス

4. Singh-Rachford, T. N., & Castellano, F. N. (2010). "Photon upconversion based on sensitized triplet–triplet annihilation." _Coord. Chem. Rev._, 254, 2560.
5. Monguzzi, A., et al. (2012). "Upconversion-induced fluorescence in multicomponent systems." _Phys. Chem. Chem. Phys._, 14, 4322.
6. Turro, N. J., Ramamurthy, V., & Scaiano, J. C. (2010). _Modern Molecular Photochemistry of Organic Molecules_. University Science Books.

### 量子計算

7. Nielsen, M. A., & Chuang, I. L. (2010). _Quantum Computation and Quantum Information_. Cambridge University Press.
8. Stinespring, W. F. (1955). "Positive functions on C*-algebras." *Proc. Amer. Math. Soc.\*, 6, 211.
9. Lloyd, S. (1996). "Universal quantum simulators." _Science_, 273, 1073.

### ボソン相互作用と数値手法

10. Tanimura, Y. (2020). "Numerically 'exact' approach to open quantum dynamics." _J. Chem. Phys._, 153, 020901.
11. May, V., & Kühn, O. (2011). _Charge and Energy Transfer Dynamics in Molecular Systems_. Wiley-VCH.

### MQT-Qudits関連

12. MQT-Qudits Documentation: https://github.com/cda-tum/mqt-qudits
13. 本リポジトリ内の関連文書：
    - `tutorials/doc/GKSL/量子ダイナミクスGKSL-Lindblad理論完全定式化.md`
    - `tutorials/doc/GKSL/GKSL-Lindblad量子ダイナミクスQiskit-Qubit完全実装理論.md`
    - `tutorials/doc/GKSL/GKSL-Lindblad量子ダイナミクスQudit完全実装理論.md`

---

## 付録A：記号一覧

| 記号                       | 意味                  | 単位/値 |
| -------------------------- | --------------------- | ------- |
| $\|S_0\rangle$             | 基底一重項状態        | -       |
| $\|T_1\rangle$             | 励起三重項状態        | -       |
| $\|S_1\rangle$             | 励起一重項状態        | -       |
| $E_T$                      | 三重項エネルギー      | 1.5 eV  |
| $E_S$                      | 一重項エネルギー      | 3.0 eV  |
| $\hat{\rho}$               | 密度演算子            | -       |
| $\hat{H}$                  | ハミルトニアン        | eV      |
| $\hat{L}_\alpha$           | Lindblad演算子        | -       |
| $\gamma_\alpha$            | 散逸速度定数          | eV/ℏ    |
| $\mathcal{D}[\hat{L}]$     | Lindblad超演算子      | -       |
| $\Gamma_{\text{fl}}$       | 蛍光発光速度          | eV/ℏ    |
| $\Gamma_{\text{ph}}$       | 燐光発光速度          | eV/ℏ    |
| $\gamma_{\text{TTA}}$      | TTA速度定数           | eV/ℏ    |
| $V_{ij}$                   | エネルギー移動積分    | eV      |
| $\hat{a}, \hat{a}^\dagger$ | ボソン消滅/生成演算子 | -       |
| $\omega$                   | 振動周波数            | eV/ℏ    |
| $g$                        | 電子-フォノン結合定数 | eV      |
| $S$                        | Huang-Rhysパラメータ  | 無次元  |
| $\tau$                     | 時間刻み（Trotter）   | fs      |

---

## 付録B：単位系変換

### B.1 エネルギー・時間

| 量         | 自然単位 (ℏ=1) | SI単位                    |
| ---------- | -------------- | ------------------------- |
| エネルギー | 1 eV           | 1.602 × 10⁻¹⁹ J           |
| 時間       | 1 ℏ/eV         | 6.58 × 10⁻¹⁶ s ≈ 0.658 fs |
| 速度定数   | 1 eV/ℏ         | 1.52 × 10¹⁵ s⁻¹           |

### B.2 よく使う変換

| 変換            | 関係式         |
| --------------- | -------------- |
| ns⁻¹ → eV/ℏ     | × 6.58 × 10⁻¹⁰ |
| eV → nm（波長） | λ = 1240/E(eV) |
| K → eV          | × 8.62 × 10⁻⁵  |

---

**文書終了**

作成日: 2026年1月25日
バージョン: 1.0.0
著者: MQT-Qudits研究グループ
ライセンス: MIT License
