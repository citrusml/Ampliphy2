# Ampliphy 変更点まとめ (2026-10-26)

比較: **旧** `~/amplify/{ampliphy, ampliphy-analysis}` → **新** `~/Ampliphy2/{ampliphy, ampliphy-analysis}`

> [!NOTE]
> 指定の `~/ampliphy` は存在せず、同名ディレクトリが `~/amplify/` 配下にあったため、これを旧版として比較した。
> `.git`, `__pycache__`, `.DS_Store` は比較対象外。

---

## 1. `ampliphy` (Nextflow パイプライン本体)

### 1.1 概要

| 観点         | 旧 (amplify)                                         | 新 (Ampliphy2)                                                                    |
| ------------ | ---------------------------------------------------- | --------------------------------------------------------------------------------- |
| 主眼         | ホモログ選択モード・MaxAlign・再現性オプションが豊富 | **分類 (taxonomy) レポートと TCS 算出**を追加。選択オプション類は削除しシンプル化 |
| ホモログ選択 | `enrich_mode` (default/random/choose), 各種フィルタ  | `max_depth` / `max_seqs` のみ                                                     |
| MaxAlign     | あり (before/after)                                  | **削除**                                                                          |
| Taxonomy     | なし                                                 | MMseqs2 の taxonomy から taxid 取得、レポート出力                                 |
| TCS          | なし                                                 | 元 / 増幅ツリーの Taxonomy Overlap Score を計算                                   |
| 最終ツリー   | `{id}.nwk`                                           | `{id}.amp.nwk` (+ 任意で `{id}.ori` 系)                                           |

### 1.2 追加されたファイル

| ファイル                                     | 役割                                                                                      |
| -------------------------------------------- | ----------------------------------------------------------------------------------------- |
| `bin/ampliphy_homolog_taxonomy.py` (292行)   | ホモログの taxid から NCBI taxdump を用いて系統情報を集計 (`homolog_taxonomy_report.tsv`) |
| `bin/ampliphy_amplified_taxonomy.py` (357行) | 入力 `.tax` + 追加ホモログから、増幅後の種/属/科/目/綱の distinct 数を集計                |
| `bin/ampliphy_input_lineage.py` (203行)      | 入力配列の taxonomy から TCS 用 lineage ファイルを作成                                    |
| `bin/tcs.py` (352行)                         | `ampliphy-analysis/scripts/tcs.py` 由来の TCS 計算 (ツリー + lineage)                     |
| `modules/homolog_taxonomy.nf`                | 上記ホモログ taxonomy レポート process                                                    |
| `modules/amplified_taxonomy.nf`              | 増幅 taxonomy レポート process                                                            |
| `modules/input_lineage.nf`                   | lineage 生成 process                                                                      |
| `modules/amplified_tcs.nf`                   | ori / amp ツリーの TCS を算出 (`{id}.tcs.tsv`)                                            |
| `sample/input/*.fa`, `sample/taxonomy/*.tax` | サンプル入力 (E0754546, E0762353, E0793441, E0802500, E0803542)                           |

### 1.3 削除されたファイル

| ファイル                                   | 内容                                                    |
| ------------------------------------------ | ------------------------------------------------------- |
| `bin/maxalign.pl`, `modules/maxalign.nf`   | MaxAlign によるアライメント前/後フィルタ                |
| `modules/degap.nf`                         | MaxAlign 出力のデギャップ (MMseqs クエリ用)             |
| `bin/filter_homologs.py`                   | ホモログのフィルタ (e-value / pident / random / choose) |
| `.gitattributes`, `.github/`, `.gitignore` | Ampliphy2 側にはコピーされていない                      |

### 1.4 `ampliphy.nf` (ワークフロー)

- MaxAlign / degap 関連の include と分岐を全削除 (`maxalign_before/after`, `maxalign_keep_original`)。
- DB 準備 channel を `multiMap` で `mafft / mmseqs / taxonomy` に分岐し、`taxonomy_status` を検索 process に渡す。
- `!params.no_taxonomy` のとき `homolog_taxonomy` を実行。
- `--input_taxonomy` 指定時:
  - `--no_taxonomy` との併用はエラー (増幅 taxonomy に taxid が必須のため)。
  - `<input_taxonomy>/*.tax` を読み、`amplified_taxonomy` を実行。ファイルが無ければエラー。
  - 元 MSA (`ori`) も IQ-TREE に流し、`input_lineage` → `amplified_tcs` で TCS を計算。
  - `amplified_tcs_report.tsv` (`family_id, tcs_ori, tcs_amp`) を集約出力。
- IQ-TREE への入力を `tuple(id, tree_type, msa, nseq)` に変更 (`tree_type` = `amp.unpruned` または `ori`、`nseq` は配列数)。
- 増幅ツリーはいったん **unpruned** で推定し、`root_and_prune` で rooting + pruning。

### 1.5 `nextflow.config`

| 項目                                                                                                           | 旧                            | 新                                           |
| -------------------------------------------------------------------------------------------------------------- | ----------------------------- | -------------------------------------------- |
| `mafft_threadit`                                                                                               | 0                             | **削除** (`--threadit 0` を付けなくなった)   |
| `max_seqs_per_query`, `enrich_mode`, `mmseqs_evalue`, `filter_evalue`, `min_pw_id`, `max_pw_id`, `random_seed` | あり                          | **削除**                                     |
| `maxalign_*` (before/after/script/options/keep_original)                                                       | あり                          | **削除**                                     |
| `no_taxonomy`, `input_taxonomy`, `ncbi_dir` (`./ncbi_taxonomy`)                                                | なし                          | **追加**                                     |
| `tcs_script` (`bin/tcs.py`), `keep_unpruned_tree` (false)                                                      | なし                          | **追加**                                     |
| `iqtree_options`                                                                                               | `-m JTT+I+G4 -B 1000`         | `-m JTT+I+G4 -B 1000 -keep-ident`            |
| ラベル `tree`                                                                                                  | なし                          | 追加 (16 GB / 24 h / threads)                |
| ラベル `mmseqs` memory                                                                                         | `params.max_memory` / `'12h'` | `16.GB` / `12.h` (固定値、Nextflow 単位表記) |
| ラベル `short`/`minimal`                                                                                       | `'4 GB'`, `'2h'` (文字列)     | `4.GB`, `2.h` (単位表記)                     |

### 1.6 modules / bin の変更

| ファイル                                     | 変更点                                                                                                                                                                                                                                                                                                        |
| -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `modules/mmseqs_prepare_db.nf`               | label を `params.minimal ? 'minimal' : 'mmseqs'` に。`no_taxonomy` / `input_taxonomy` を渡し、出力に `mmseqs_taxonomy_status.txt` を追加                                                                                                                                                                      |
| `bin/ampliphy_mmseqs_prepare_db.sh`          | `--no-taxonomy`, `--require-taxonomy`, `--taxonomy-status` 追加。DB の `_mapping` / `_taxonomy` の有無で `enabled` / `disabled` を判定。`--input_taxonomy` 指定なのに taxonomy 無し、または `--no_taxonomy` 併用は **エラー終了**。それ以外で taxonomy 無しは警告のみ                                         |
| `modules/mmseqs_search.nf`                   | 入力に `taxonomy_status_file` 追加。出力に `{id}.homolog_search_stats.tsv`、`{id}.homolog_taxids.tsv` (optional) 追加。渡す引数を `--taxonomy-status` 追加、`--mmseqs-evalue`, `--max-seqs-per-query`, `--mode`, `--seed`, filter 系を削除                                                                    |
| `bin/ampliphy_mmseqs_search.sh`              | 引数削減 (mode/filter/seed など削除)、`--taxonomy-status` 追加、taxonomy 有効時に taxid と検索統計を出力 (約 124 行追加 / 82 行削除)                                                                                                                                                                          |
| `modules/mafft_align.nf`, `mafft_amplify.nf` | `--threadit` 削除。`mafft_amplify` は **ホモログが空の場合は元 MSA をそのままコピー**して警告 (旧はエラーの可能性)                                                                                                                                                                                            |
| `modules/iqtree_inference.nf`                | 入出力を `(id, tree_type, msa, nseq)` → `(id, tree_type, {id}.{tree_type}.nwk)` に変更。**配列数 < 4 ならブートストラップ系オプション (`-B`, `-b`, `--boot`, `--ufboot`, `-bnni`) を自動除去**して警告。`-T AUTO --threads-max` → `-T ${threads}`。`keep_unpruned_tree` 時のみ `.amp.unpruned.nwk` を publish |
| `modules/root_and_prune.nf`                  | label `minimal`→`short`。入力/出力を tuple 化 (`emit: amplified_trees`)。MAD に `-n` を付与、`chmod +x`。ホモログが空のときは pruning せず rooting 結果をそのまま出力。出力名 `{id}.nwk`→`{id}.amp.nwk`                                                                                                       |
| `envs/ampliphy.yml`                          | `python>=3.9` を追加                                                                                                                                                                                                                                                                                          |
| `README.md`                                  | 入力形式の例、`--input_taxonomy` の `.tax` 形式 (2列 TSV) などを追記。出力レイアウト節・Main options 節を書き換え (旧の maxalign / enrich_mode 説明は削除)                                                                                                                                                    |

> [!WARNING]
> 新版では旧版の `enrich_mode`, `filter_evalue`, `min/max_pw_id`, `mmseqs_evalue`, `random_seed`, MaxAlign, `mafft_threadit` (再現性のための `--threadit 0`) が使えない。
> 旧オプション前提の実験を再実行する場合は旧版 (`~/amplify/ampliphy`) を使うこと。`mafft --threadit 0` が無くなったため、マルチスレッド時の MSA 再現性にも注意。

---

## 2. `ampliphy-analysis`

| ファイル               | 状態                        | 内容                                                                                                                                                                  |
| ---------------------- | --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `scripts/tcs_batch.py` | 差分あり (**実質変更なし**) | AST 比較 (docstring 除く) で **ロジック同一**。差分は整形 (black 風 `"` クォート・空行・折返し) のみ。旧側が整形済み (git 未コミット変更あり)、Ampliphy2 側は元の書式 |
| `scripts/mad`          | 旧のみに存在 (git 未追跡)   | LSF 環境で `foldtree` conda/micromamba 環境の Python から `~/tools/mad2.2/mad/mad.py` を起動するラッパー                                                              |
| `scripts/__pycache__/` | 旧のみ                      | キャッシュ (無視可)                                                                                                                                                   |
| `.git`                 | 旧のみ                      | Ampliphy2 側は git 管理外                                                                                                                                             |
| `.DS_Store`            | 差異あり                    | 無視可                                                                                                                                                                |

上記以外の `scripts/` (aligner_*, builder_*, run_*, plot_*, get_*, write_hogs* など) と `dataset/`, `LICENSE` は **両者で同一**。

---

## 3. 全体所見

- 新版 `ampliphy` は「ホモログ選択の多様なオプション」から「**taxonomy 情報の付与と TCS による評価**」へ方針転換している。
- `ampliphy-analysis` に実質的な変更はない (`tcs_batch.py` の整形差のみ)。
- 新版 `bin/tcs.py` は `ampliphy-analysis/scripts/tcs.py` と同系統のスクリプトをパイプライン内に取り込んだもの。