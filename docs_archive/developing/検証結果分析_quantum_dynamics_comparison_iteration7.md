# 量子ダイナミクス完全比較ノートブック 検証結果分析レポート（Iteration 7）

## 作成日: 2026-03-10
## 対象ファイル:
- `tutorials/qubit_noisy_simulator.py`
- `tutorials/mqt_qudits_noisy_simulator.py`
- `tutorials/run_iteration7_verification.py`
## 前回参照: `developing/検証結果分析_quantum_dynamics_comparison_iteration6.md`

---

## 1. 問題の再整理

### 1.1 PR#222 Iteration 6での修正内容

Iteration 6では `qubit_noisy_simulator.py` を以下のように書き換えた:
- **変更前**: Qiskit Aerゲートレベルノイズ（~814 CXゲート/ステップ → 81,400ノイズイベント/100ステップ）
- **変更後**: 密度行列ペアレベルノイズ（12ペアゲート/ステップ → 1,200ノイズイベント/100ステップ）

この書き換えの目的は、qubitとquditのノイズシミュレータの挙動を「平等にする」ことだった。

### 1.2 問題の本質

**ユーザーの指摘**: qubit表現とqudit表現でのノイズシミュレータの挙動を平等にする必要はない。

**ノイズシミュレータの目的**: 実機での挙動をシミュレーションすること。

**物理的事実**:
1. qudit（qutrit、d=3）表現は、この量子ダイナミクス（TTA-UC 4分子系）に対して物理的に自然な表現
2. 各分子の3つの状態（S0, T1, S1）が直接1つのqutritに対応
3. 2分子間相互作用は単一の2-qutritネイティブゲートで実装可能

一方、qubit（d=4）表現では:
1. 各分子に2量子ビットが必要（|00⟩=S0, |01⟩=T1, |10⟩=S1, |11⟩=forbidden）
2. 2分子間相互作用は複数のCXゲートに分解が必要
3. 脱分極ノイズにより禁止状態 |11⟩ へのリーケージが発生

### 1.3 結論: qudit表現の方が精度が向上する理由

qudit表現がqubit表現より優れる理由は以下の2つの独立した物理的効果による:

**効果1: ゲートオーバーヘッド（Gate Overhead）**
- qudit: 1ネイティブ2-qutritゲート/ペア相互作用
- qubit: N_CX個のCXゲート/ペア相互作用
- 同一の物理ゲートエラー率 p_phys に対して:
  - qudit有効ペアエラー = p_phys
  - qubit有効ペアエラー = 1 - (1-p_phys)^N_CX >> p_phys

**効果2: 禁止状態リーケージ（Forbidden State Leakage）**
- qudit (d=3): 全9状態が物理的 → リーケージ = 0
- qubit (d=4): 16状態中7状態が非物理的 → リーケージ = 7/16 = 43.75%/ノイズイベント

---

## 2. 修正内容

### 2.1 `qubit_noisy_simulator.py` の修正

**変更の要約**: `simulate()` メソッドに `cx_per_pair_gate` パラメータを追加し、ゲートオーバーヘッドを正確にモデル化。

**詳細**:

1. **`cx_per_pair_gate` パラメータ**: ペア相互作用あたりのCXゲート数
   - デフォルト値: 1（quditと同一、後方互換性を維持）
   - 実際の値: Qiskit transpilationにより計算（`estimate_cx_per_pair_gate()` メソッド）

2. **有効脱分極率の計算**:
   ```python
   depol_2q_eff = 1 - (1 - depol_2q) ** cx_per_pair_gate
   ```

3. **`estimate_cx_per_pair_gate()` 静的メソッド**:
   - Qiskitの transpile (optimization_level=3) を使用
   - H_transfer と H_TTA の16×16ユニタリを分解
   - CXゲート数をカウント

4. **出力の拡張**: 
   - `depol_2q_physical`: 物理CXゲートエラー率
   - `depol_2q_effective`: 有効ペア脱分極率
   - `cx_per_pair_gate`: ペアあたりのCXゲート数

### 2.2 `mqt_qudits_noisy_simulator.py` は変更なし

quditシミュレータは既に物理的に正しいモデルを使用:
- 各ペア相互作用に1つのネイティブ2-qutritゲート
- `depol_2q` はそのまま有効ペアエラー率として使用
- 変更不要

---

## 3. ヒューリスティック・フォールバックの不使用

**有効脱分極率の計算 `p_eff = 1-(1-p)^N` は標準的な量子情報理論の公式であり、ヒューリスティックではない。**

根拠:
- N個の独立な脱分極チャネルの合成に対する正確な公式（Nielsen & Chuang, 2000）
- 各CXゲートが独立に脱分極エラーを起こす場合の正確な有効確率
- 注: 厳密には、異なるqubitペアに作用するCXゲートの合成は単純な脱分極チャネルにはならない。しかし、ペアレベルの脱分極モデルとして有効率を使用することは、物理的に妥当な近似である

---

## 4. 検証スクリプト

```bash
cd tutorials && python run_iteration7_verification.py
```

結果は `developing/verification_results/iteration7_*.json` に保存される。

### 検証項目

1. **コード構造検証** (Section 1-3):
   - `cx_per_pair_gate` パラメータの存在
   - 有効脱分極率の計算式
   - quditシミュレータの非変更
   - CPTP特性

2. **CXゲート数計測** (Section 5):
   - Qiskitが利用可能な場合のみ実行
   - H_transfer と H_TTA のCX数を計算

3. **有効率分析** (Section 6):
   - 各CX数に対する有効率の計算
   - 単調増加性の検証

4. **リーケージ分析** (Section 7):
   - qubit: 7/16 = 43.75% 理論値との整合
   - qudit: 0% リーケージ

---

## 5. 次回の検証で確認すべき事項

1. **CXゲート数の実測値**: ユーザー環境（Qiskit installed）で `estimate_cx_per_pair_gate()` を実行し、H_transfer と H_TTA のCXゲート数を確認
2. **ノイズ有りシミュレーション結果**: 
   - `cx_per_pair_gate=1` (現行): qubit/quditの差異はリーケージのみ
   - `cx_per_pair_gate=N_CX` (実機): qubitはさらに劣化
3. **ノートブックの更新**: CXゲート数に基づいてノートブックのパラメータを更新
4. **全セルの正常実行**: 更新されたパラメータでノートブックが正常に動作すること

---

## 6. Iteration 6からの変更点

| 項目 | Iteration 6 | Iteration 7 |
|------|------------|------------|
| qubit_noisy_simulator | 密度行列（ペアレベル） | 密度行列（ペアレベル + **CXオーバーヘッド**） |
| ノイズ粒度の「平等化」 | 両方12イベント/ステップ | **不平等は物理的事実として反映** |
| depol_2qの意味 | ペアあたりエラー率 | **物理CXゲートエラー率** |
| 有効ペアエラー率 | depol_2q (同一) | qubit: 1-(1-p)^N_CX, qudit: p |
| CXゲート数推定 | なし | **Qiskit transpileによる計測** |
| quditの変更 | Iteration 5から変更なし | **変更なし**（正しいため） |
