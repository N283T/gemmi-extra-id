# Molecule-ID 付加プログラム 計画

## 目的
- gemmi を使って mmCIF を読み込み、glossary の命名規則に基づく `chain_*` / `pn_unit_*` / `molecule_*` の ID を付与する。
- 付加結果は mmCIF への追記、またはサイドカー(JSON/CSV)として出力できるようにする。

## 前提整理 (glossary 準拠)
- `entity` = 共有の化学的同一性 (座標は問わない)
- `instance` = 3D 空間上の実体 (transform で増える)
- `chain`: `(asym_id, transformation_id)` の組
- `pn_unit`: 同種 chain の共有結合で連結した集合
- `molecule`: 共有結合グラフの連結成分

## 参考リポジトリ確認 (最小限)
1. gemmi:
   - mmCIF 読み込み API (Python: `gemmi.cif`, `gemmi.read_structure`)
   - `struct_conn`, `chem_comp_bond` 等の結合情報の取り扱い
   - assembly/transform の展開方法
2. atomworks:
   - molecule-id の付与ロジックと ID 命名
   - 実装上のデータ構造 (chain/instance の識別方法)

## 実装方針 (段階的)
1. **入力/出力仕様の確定**
   - 入力: 1 mmCIF (ASU/assembly どちらを対象にするか)
   - 出力: mmCIF への追記 or JSON/CSV サイドカー
   - ID の具体形式 (数値/文字列) と安定性方針
2. **mmCIF パース層**
   - gemmi で構造体/カテゴリを読み込み
   - chain 単位の原子集合と `entity_id`, `asym_id` を取得
   - assembly 対応が必要なら transformation を展開
3. **結合グラフ構築**
   - 共有結合の情報源を整理 (`struct_conn`/`chem_comp_bond` 等)
   - atom レベルの共有結合グラフを構築
4. **ID 付与ロジック**
   - `chain_iid`: `(asym_id, transform_id)` を基準に採番
   - `pn_unit_iid`: 同種 chain の共有結合で連結した集合
   - `molecule_iid`: 共有結合グラフの連結成分
   - `*_entity` は `entity_id` ベースで付与
5. **出力層**
   - mmCIF のカスタムカテゴリ or サイドカーへ出力
   - 追記する場合は `pdbx_` など衝突回避プレフィックスを検討
6. **検証・テスト**
   - 最小ケース (単一 chain) / 複数 chain / 共有結合のある複合体
   - assembly あり/なしの比較
   - atomworks の期待結果と比較できるサンプルの整備

## 追加で確認したい点
- 対象は ASU のみか、assembly まで展開するか
- 出力形式の希望 (mmCIF 追記 vs JSON/CSV)
- `struct_conn` が無い場合の扱い (推定結合を行うか)

## ざっくりスケジュール (目安)
- 1日目: 参考リポジトリ・仕様整理
- 2日目: パース層 + 結合グラフ
- 3日目: ID 付与 + 出力 + テスト
