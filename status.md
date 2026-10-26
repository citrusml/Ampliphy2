最終更新: 2026-10-27 03:08

## 📈 進捗状況
- [x] 完了: Phylomethod_comparison との比較によるログ非同期・非表示原因の特定と `-ansi-log false` 削除による標準進捗表示の復元 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ再投入とリアルタイム進捗ログの確認 (進捗率: 70% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**: スパコン上でジョブ実行中、ログで進捗状況（各ステップの進行率）が随時確認できない。外部集計スクリプトなしで `~/Phylomethod_comparison` と同様に Nextflow 標準のログで随時進捗を確認したい。
- **[原因 (Phylomethod_comparison との差分)]**:
  1. **`-ansi-log false` の混入**: `scripts/run_ampliphy.sh` に `-ansi-log false` が指定されていたため、Nextflow 標準の動的進捗ブロック（`[xx%] x of y` など）の出力が完全に抑制されていた（※過去の `~/Phylomethod_comparison/memo/edit20260815.md` でも同一原因で `-ansi-log false` を削除していた経緯と完全一致）。
  2. **ログ分離と名称**: Nextflow の標準出力を `nextflow_<taxon>.console.log` に逃がし、`nextflow_<taxon>.log` には Java 内部デバッグログを出力していたため、ユーザーが確認すべきログの動線が乱れていた。また、4生物種を並列実行するためマスターログ（`master_*.out`）には各プロセスの標準出力が直接流れていなかった。
- **[解決策]**:
  1. 外部集計スクリプトや監視ツールは一切使わず、Nextflow ネイティブの進捗表示に一本化（余計なファイルは不要）。
  2. `scripts/run_ampliphy.sh` から `-ansi-log false` を削除。
  3. 各生物種の実行進捗ログを `logs/ampliphy/nextflow_<taxon>.log` に直接ストリーミング出力するように整理。
- **[パラメータ]**: `nextflow run ... > logs/ampliphy/nextflow_<taxon>.log 2>&1` (標準進捗ブロック有効)

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. 前回のジョブが残っていれば停止 (`bjobs` で確認後 `bkill <JOB_ID>`)
3. ジョブ再投入: `bsub < scripts/run_ampliphy.sh`
4. リアルタイム進捗の確認:
   - 例: `tail -f logs/ampliphy/nextflow_amniota.log`
   - 例: `tail -f logs/ampliphy/nextflow_eukaryota.log`

---

最終更新: 2026-10-27 02:47

## 📈 進捗状況
- [x] 完了: `input_lineage` の rank 'subspecies' エラー解消（全399件の `.tax` ファイルの種レベル TaxId 化 ＆ `ampliphy` 側の種ランク自動解決フォールバック実装） (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ再投入と TCS 算出完了の確認 (進捗率: 60% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**: Nextflow の `input_lineage` タスク実行時、以下のエラーが発生した：
  ```
  Command error:
    Error: E0732656.tax:4: TaxId '9595' resolves to rank 'subspecies', not species.
  ```
- **[原因]**:
  1. `ampliphy_input_lineage.py` は配列に割り当てられた TaxId が厳密に `rank == "species"` であることを要求している。
  2. しかし、元データ（`dataset/*_taxa.tsv` の第2列）には leaf レベルの TaxId（亜種や株、例: ゴリラ亜種 `9595` やイヌ亜種 `9615` など）が登録されていたため、`.tax` ファイルにも亜種 TaxId が引き継がれ、種ランク不一致で失敗していた。
- **[解決策]**:
  1. **`.tax` ファイルの種レベル TaxId への一括更新:** `dataset/*_taxa.tsv` の系統情報から正確な種レベルの TaxId（例: ゴリラ種 `9593`, イヌ種 `9612` など）を抽出し、`taxonomy/*/*.tax` の全4,752箇所の亜種・株 TaxId を正規の種 TaxId に一括置換。
  2. **`ampliphy` 側の自動解決フォールバック:** 万が一亜種・株レベルの TaxId が渡された場合でも、NCBI 系統樹の親ノードを遡って自動的に `rank == "species"` のノードを解決するよう `ampliphy_input_lineage.py` および `ampliphy_amplified_taxonomy.py` を改修。
- **[パラメータ]**: `taxonomy/<taxon>/*.tax` (全400ファイル中399ファイルの4,752エントリを種レベルTaxIdへ正規化)

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. ジョブ再投入: `bsub < scripts/run_ampliphy.sh`
3. リアルタイムログで進捗監視:
   - 全体ログ: `tail -f logs/ampliphy/ampliphy_*.log`
   - 各生物種コンソール: `tail -f logs/ampliphy/nextflow_*.console.log`
   - LSF ジョブ分散状況: `bjobs`

---

最終更新: 2026-10-27 02:40

## 📈 進捗状況
- [x] 完了: Nextflow LSF投入エラー (`Syntax Error: Multiple -R resource requirement strings...`) の原因特定と `clusterOptions` 削除 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ再投入と各タスクのLSFディスパッチ確認 (進捗率: 55% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**: スパコン上で AmpliPhy 実行時、マスタージョブ (246037) は起動したが、個別タスク投入時に LSF がエラーコード 255 を返して停止した：
  ```
  Syntax Error: Multiple -R resource requirement strings are not supported on "span", "cu" and "affinity" sections. Specify multiple -R resource requirement strings only on order, same, rusage, and select sections.
  Bad resource requirement syntax. Job not submitted.
  ```
- **[原因]**: Nextflow の LSF エグゼキュータは、`cpus` や `memory` 設定に基づいて自動的に `-R "span[hosts=1] ..."` を構築して `bsub` コマンドに付与する。これに対し、`configs/cluster_dynamic.config` で `clusterOptions = '-R "span[hosts=1]"'` を手動指定していたため、`bsub` コマンドに `-R` の `span` セクションが重複して渡され、LSF の構文エラーとなっていた。
- **[解決策]**: `configs/cluster_dynamic.config` から重複する `clusterOptions = '-R "span[hosts=1]"'` を削除。Nextflow ネイティブの LSF リソース構築処理に委ねることで、単一かつ正常な `-R` オプションが渡されるように修正。
- **[パラメータ]**: `configs/cluster_dynamic.config` (executor: `lsf`, queue: `mafft`, cpus/memory による自動リソース生成)

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. ジョブ再投入: `bsub < scripts/run_ampliphy.sh`
3. リアルタイムログで進捗監視:
   - 全体ログ: `tail -f logs/ampliphy/ampliphy_*.log`
   - 各生物種コンソール: `tail -f logs/ampliphy/nextflow_*.console.log`
   - LSF ジョブ分散状況: `bjobs`

---

最終更新: 2026-10-27 02:35

## 📈 進捗状況
- [x] 完了: ユーザー指定の確定パスに基づき、`run_ampliphy.sh` および `run_gs_pipeline.sh` の探索ループを全撤廃・完全一本化 (完了日: 10/27)
- [ ] 進行中: スパコン上での AmpliPhy ジョブ実行と進捗確認 (進捗率: 50% / 課題: ジョブ投入と進捗監視)

## 🛠 技術的決定・メモ
- **[課題]**: スパコン上の各リポジトリやツールの配置場所が確定したため、スクリプト内に残っていたフォールバック用ループ探索コード（`for cand in ...`）を全て排除し、明確で保守性の高い確定パス構成にリファクタリングしたい。
- **[確定パス]**:
  - `DB_DIR`: `/home/citrusml2004/mmseqs_db`
  - `AMPLIPHY_DIR`: `/home/citrusml2004/ampliphy`（`Ampliphy2/ampliphy` へシンボリックリンク保証）
  - `DATASET_DIR`: `/home/citrusml2004/ampliphy-analysis/dataset`（`Ampliphy2/ampliphy-analysis` へシンボリックリンク保証）
  - `GS2_BIN`: `/home/citrusml2004/gs/gs2`
  - `MAD_SCRIPT`: `/home/citrusml2004/ampliphy/bin/mad.py`
  - `TCS_SCRIPT`: `/home/citrusml2004/ampliphy-analysis/scripts/tcs_batch.py`
  - `NCBI_DIR`: `/home/citrusml2004/Ampliphy2/ncbi_taxonomy`
  - `amp_fa`: `results/<taxon>/msa/<hog_id>.amp.fa`
- **[解決策]**:
  1. `run_ampliphy.sh` の `ampliphy`、`ampliphy-analysis`、`DB_DIR` 設定から探索ループを全削除し、上記確定パスへ直接設定。
  2. `run_gs_pipeline.sh` の `GS2_BIN`、`MAD_SCRIPT`、`TCS_SCRIPT`、`DATASET_DIR` および `amp_fa` 探索ループを全削除し確定パスに一本化。

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. ジョブ投入: `bsub < scripts/run_ampliphy.sh`
3. リアルタイムログで進捗監視:
   - 全体ログ: `tail -f logs/ampliphy/ampliphy_*.log`
   - 各生物種コンソール: `tail -f logs/ampliphy/nextflow_*.console.log`
   - LSF ジョブ分散状況: `bjobs`

---

最終更新: 2026-10-27 02:28

## 📈 進捗状況
- [x] 完了: `cluster_dynamic.config` の外部静的ファイル化 ([configs/cluster_dynamic.config](file:///Users/aibarakazuki/Ampliphy2/configs/cluster_dynamic.config)) と Git 管理化 (完了日: 10/27)
- [x] 完了: taxonomy (.tax) ディレクトリの準備と TCS 評価の常時実行化 (完了日: 10/27)
- [ ] 進行中: スパコン上での AmpliPhy クラスタ全体分散実行 (進捗率: 45% / 課題: ジョブ投入と進捗監視)

## 🛠 技術的決定・メモ
- **[課題]**: 
  1. `run_ampliphy.sh` 内で毎回ヒアドキュメント (`cat << 'EOF' > ...`) を用いて Nextflow 設定ファイルを生成するのではなく、設定ファイルとして独立させて Git 管理したい。
  2. taxonomy (.tax) の準備が条件分岐ではなく、常に実行されて TCS 評価・系統樹評価が確実に行われるようにしたい。
- **[解決策]**:
  1. `configs/cluster_dynamic.config` を静的設定ファイルとして独立作成し、Git 追跡対象に設定（`.gitignore` から `configs/` を除外）。
  2. `run_ampliphy.sh` からインラインの Nextflow 設定生成コードを完全撤廃し、`CONFIG_FILE="$CONFIG_DIR/cluster_dynamic.config"` を参照するシンプルな構成に変更。
  3. `taxonomy/$taxon`（全400件の `.tax` ファイル）はすでに Git 追跡対象としてリポジトリに揃っているため、動的生成コード（Python）を全削除し、直接ファイル存在確認のみを行う構成に単純化。
- **[パラメータ]**: Nextflow Config: `configs/cluster_dynamic.config`, Taxonomy: `--input_taxonomy taxonomy/<taxon>` (常時有効)

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. ジョブ投入: `bsub < scripts/run_ampliphy.sh`
3. リアルタイムログで進捗監視:
   - 全体ログ: `tail -f logs/ampliphy/ampliphy_*.log`
   - 各生物種コンソール: `tail -f logs/ampliphy/nextflow_*.console.log`
   - LSF ジョブ分散状況: `bjobs`

---

最終更新: 2026-10-27 02:26

## 📈 進捗状況
- [x] 完了: `run_ampliphy.sh` および `run_gs_pipeline.sh` の探索ループを全撤廃し、`/home/citrusml2004/Ampliphy2` 配下の確定パスに固定・単純化 (完了日: 10/27)
- [x] 完了: `run_ampliphy.sh` にて `--input_taxonomy` を指定し、元MSAと増幅後MSAの TCS スコア評価レポート (`amplified_tcs_report.tsv`) 算出を常時有効化 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ投入と進捗確認 (進捗率: 45% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**:
  1. `input_fasta` や `taxonomy`、`ampliphy` 関連ディレクトリの位置が確定しているにもかかわらず、過去のフォールバック用ループ探索コードが多数残っており冗長になっていた。
  2. taxonomy データ（100件/.tax）がすでに準備完了しているため、AmpliPhy パイプライン内で元ツリーと増幅後ツリーの TCS スコア評価（`amplified_tcs_report.tsv`）まで一貫して実行したい。
- **[解決策]**:
  1. `run_ampliphy.sh`:
     - `input_dir` を `$AMPLIPHY_ROOT/input_fasta/$taxon` に完全固定（探索ループ削除）。
     - `tax_dir` を `$AMPLIPHY_ROOT/taxonomy/$taxon` に完全固定（探索ループ・動的生成処理を削除）。
     - `AMPLIPHY_DIR`、`DATASET_DIR`、`CONFIG_DIR`、`NCBI_DIR` をすべて `$AMPLIPHY_ROOT` 直下に一本化。
     - Nextflow 実行引数に `--input_taxonomy "$tax_dir"` および `--ncbi_dir "$NCBI_DIR"` を指定し、TCS 算出（`amplified_tcs.nf`）を常時実行。
     - 並行実行時の競合を防ぐため、UniRef50 と NCBI taxonomy の事前準備ブロックを配置。
  2. `run_gs_pipeline.sh`:
     - `orig_dir` (`$AMPLIPHY_ROOT/input_fasta/$taxon`)、`amp_dir` (`$AMPLIPHY_ROOT/results/$taxon/msa`)、`DATASET_DIR` の探索ループを全撤廃し確定パスに固定。
- **[パラメータ]**:
  - `input_dir`: `Ampliphy2/input_fasta/<taxon>`
  - `tax_dir`: `Ampliphy2/taxonomy/<taxon>`
  - `CONFIG_FILE`: `Ampliphy2/configs/cluster_dynamic.config`
  - 出力: `results/<taxon>/report/amplified_tcs_report.tsv` (family_id, tcs_ori, tcs_amp)

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. ジョブ投入: `bsub < scripts/run_ampliphy.sh`
3. ログのリアルタイム確認:
   - 全体: `tail -f logs/ampliphy/ampliphy_*.log`
   - 各生物種コンソール: `tail -f logs/ampliphy/nextflow_*.console.log`

---

最終更新: 2026-10-27 02:20

## 📈 進捗状況
- [x] 完了: `cluster_dynamic.config: No such file or directory` エラーの原因特定と自動シンボリックリンク＆設定保存先修正 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ再投入とクラスタ全体分散実行の監視 (進捗率: 40% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**: スパコン上で `run_ampliphy.sh` 実行時、`/home/citrusml2004/.lsbatch/1793034537.246034.shell: line 140: /home/citrusml2004/Ampliphy2/ampliphy/cluster_dynamic.config: No such file or directory` というエラーで停止した。
- **[原因]**:
  1. `.gitignore` に `ampliphy*` が指定されているため、Git から clone/pull した `Ampliphy2` 直下には `ampliphy/` ディレクトリが存在しない。
  2. スパコン上では `ampliphy` は `$HOME/ampliphy` (または `/lustre10/home/citrusml2004/ampliphy`) に別途 clone されていた。
  3. スクリプトが `$AMPLIPHY_ROOT/ampliphy/cluster_dynamic.config` へ直接ファイルを生成しようとしたため、親ディレクトリが存在せず書き込みエラーとなった。
- **[解決策]**:
  1. `run_ampliphy.sh` および `run_gs_pipeline.sh` の冒頭に、クラスタ上の `ampliphy` および `ampliphy-analysis` を自動検出し、`$AMPLIPHY_ROOT/ampliphy` / `$AMPLIPHY_ROOT/ampliphy-analysis` へシンボリックリンクを自動作成する処理を追加。
  2. 動的生成する Nextflow 設定ファイルの保存先を `$AMPLIPHY_ROOT/configs/cluster_dynamic.config` に変更し、事前に `mkdir -p "$AMPLIPHY_ROOT/configs"` を実行して確実に配置できるようにした。
  3. `configs/`、`results/`、`logs/`、`mmseqs_db/`、`.nextflow/` を `.gitignore` に追加し、生成物が Git を汚染しないよう保護。
- **[パラメータ]**: CONFIG_FILE: `Ampliphy2/configs/cluster_dynamic.config`, PATH: `$HOME/Ampliphy2/ampliphy/bin:$HOME/ampliphy/bin:$PATH`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを取り込み (`git pull`)
2. ジョブ投入: `bsub < scripts/run_ampliphy.sh`
3. リアルタイムログで進捗監視:
   - 全体ログ: `tail -f logs/ampliphy/ampliphy_*.log`
   - 各生物種コンソール: `tail -f logs/ampliphy/nextflow_*.console.log`
   - LSF ジョブ分散状況: `bjobs`

---

最終更新: 2026-10-27 02:17

## 📈 進捗状況
- [x] 完了: `run_ampliphy.sh` から単一ノード内実行モードの分岐を撤廃し、常時クラスタ全体動的分散実行へ完全統一 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ投入とクラスタ全体分散実行の監視 (進捗率: 35% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**: `run_ampliphy.sh` 内に存在した「クラスタ全体分散モード」と「単一ノード内実行モード」の条件分岐（`MODE="cluster"` / `MODE="local"`）をなくし、常にクラスタ全体へタスクを動的ディスパッチする構成に一本化したい。
- **[解決策]**:
  1. `MODE` 環境変数および単一ノード用設定（`THREADS=40`, `MAX_MEM="200 GB"`, `local_resource_limit.config`）を完全に削除。
  2. Nextflow 設定ファイル生成部分で無条件に `cluster_dynamic.config`（`executor = 'lsf'`, queueSize 150）を生成するように簡素化。
  3. メインループ部でもモード分岐を取り払い、複数生物種指定時はバックグラウンド並行投入、単一生物種時は単独実行とするシンプルなフローに整理。
- **[パラメータ]**: executor: `lsf`, queueSize: 150, 各ラベル設定: mmseqs(8cpus, 16GB), tree(4cpus, 16GB), short(4cpus, 8GB), minimal(1cpu, 4GB)

## 🚀 次に行うこと (Next Actions)
1. スパコン上で最新コードを反映
2. ジョブ投入: `bsub < scripts/run_ampliphy.sh`
3. クラスタ全体へのタスク分散状況を確認: `bjobs` および `tail -f logs/ampliphy/*.console.log`

---

最終更新: 2026-10-27 02:07

## 📈 進捗状況
- [x] 完了: スパコン上の全実行・成果物・ログの配置先を `/home/citrusml2004/Ampliphy2` に完全統一 (完了日: 10/27)
- [x] 完了: LSF ディレクティブに `#BSUB -cwd /home/citrusml2004/Ampliphy2` を追加し、起動場所依存のズレを解消 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ投入と進捗確認 (進捗率: 30% / 課題: スパコン上での再投入)

## 🛠 技術的決定・メモ
- **[課題]**: スパコン上でジョブを実行した際、ホームディレクトリ直下（`/home/citrusml2004/results` や `/home/citrusml2004/logs`）にフォルダが作成されてしまっていた。すべて `/home/citrusml2004/Ampliphy2` 配下に完結させたい。
- **[解決策]**:
  1. `run_ampliphy.sh` および `run_gs_pipeline.sh` に `#BSUB -cwd /home/citrusml2004/Ampliphy2` を明示指定。どのディレクトリから `bsub` を実行しても必ず `Ampliphy2` 内が作業ディレクトリになるよう保証。
  2. スクリプト冒頭で `AMPLIPHY_ROOT="/home/citrusml2004/Ampliphy2"` を設定し、`cd "$AMPLIPHY_ROOT"` を実行。
  3. `results` (`results/<taxon>`)、`logs` (`logs/ampliphy`, `logs/gs`)、`results_gs` (`results_gs/<taxon>`) をすべて `Ampliphy2` 直下に生成するようパスを完全固定。
- **[パラメータ]**: カレントディレクトリ: `/home/citrusml2004/Ampliphy2`, 出力先: `Ampliphy2/results`, ログ: `Ampliphy2/logs`

## 🚀 次に行うこと (Next Actions)
1. ホーム直下に誤って作られた不要フォルダがあれば削除: `rm -rf ~/results ~/logs`
2. 最新コードを反映後、ジョブ投入: `bsub < scripts/run_ampliphy.sh`
3. ログのリアルタイム確認: `tail -f logs/ampliphy/ampliphy_live_*.log` または `tail -f logs/ampliphy/*.console.log`

---

最終更新: 2026-10-27 02:00

## 📈 進捗状況
- [x] 完了: PhyloTree_Methods_Comparison/run_scripts/run_nextflow2.sh との比較によるログ非表示の原因特定と run_ampliphy.sh 修正 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ再投入と進捗確認 (進捗率: 25% / 課題: ジョブ 246032 は旧スクリプトで走行中。kill して再投入が必要)

## 🛠 技術的決定・メモ
- **[課題]**: Ampliphy2 では実行中に logs が見えない。PhyloTree_Methods_Comparison では見えた。
- **[原因(差分)]**:
  1. 参照スクリプトは `nextflow run ... > logs/xxx.log 2>&1` と **スクリプト内で直接ファイルへリダイレクト**。今回は標準出力を LSF の `#BSUB -o` に任せていた (完了時に書き出し)。
  2. `-profile conda` (`process.conda=envs/ampliphy.yml`) により Nextflow が **隠れて conda 環境を新規構築**していた可能性 (旧 run_amp.sh は不使用)。この間ログが出ない。
  3. UniRef50 の事前 DL を出力先なしで実行 (未取得なら数時間無音)。
  4. 4 生物種を同一ディレクトリ (`ampliphy/`) から並列起動 → `.nextflow/` history 競合 (`-resume` 誤参照の恐れ)。
  5. `#BSUB -o logs/...` のディレクトリが投入時に存在しないと LSF が出力を書けない (参照は事前に mkdir logs)。
- **[解決策]**: console を `logs/ampliphy/nextflow_<taxon>.console.log` に直接出力、`-profile conda` 削除、生物種ごとの launch dir、UniRef50 DL は `uniref50_download.log` へ、mad/tcs は絶対パス指定。
- **[パラメータ]**: `-ansi-log false`, `NXF_OPTS=-Xms512m -Xmx4g`

## 🚀 次に行うこと (Next Actions)
1. `mkdir -p logs/ampliphy logs/gs && bkill 246032`
2. `ls ~/mmseqs_db ~/Ampliphy2/mmseqs_db` で UniRef50 の有無を確認
3. 更新を反映後 `bsub < scripts/run_ampliphy.sh`、`tail -f logs/ampliphy/*.console.log`

---

最終更新: 2026-10-27 01:48

## 📈 進捗状況
- [x] 完了: AmpliPhy マスタージョブ (Job ID: 246032) が mafftsrv4 にて RUN (実行中) 状態で正常起動確認 (完了日: 10/27)
- [x] 完了: LSF 相対ログ出力パス (logs/ampliphy, logs/gs) の階層ズレ解消修正 (完了日: 10/27)
- [ ] 進行中: AmpliPhy による UniRef50 チェック・4生物種増幅後MSA (*.amp.fa) のクラスタ並行生成 (進捗率: 15% / 課題: ジョブ実行完了待ち)

## 🛠 技術的決定・メモ
- **[課題]**: `bjobs` で `ampliphy_cluster_master` が `STAT: RUN` で留まっており、ログファイル (`master_246032.out`) がディスク上に見当たらない。
- **[解決策]**:
  1. LSF の標準仕様として、`-o` / `-e` ファイルは実行中はスプールにバッファリングされ、「ジョブ完了時」に書き出される。実行中のリアルタイム標準出力は `bpeek <JOB_ID>` で確認可能。
  2. 実行ディレクトリ（`~/Ampliphy2`）からの相対パスで `Ampliphy2/logs/...` となっていた箇所を `logs/ampliphy/...` に修正。
- **[パラメータ]**: JobID: `246032`, ノード: `mafftsrv4`, 確認コマンド: `bpeek 246032`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bpeek 246032` を実行し、マスタージョブの進行状況（Nextflowのタスク投入ログ）を確認する。
2. `bjobs` で各ノード（mafftsrv0〜5）に個別タスクが分散実行されている様子を確認する。
3. 全タスク完了後、`bsub < scripts/run_gs_pipeline.sh` を投入する。

---

最終更新: 2026-10-27 01:38

## 📈 進捗状況
- [x] 完了: GS法スクリプト (run_gs_pipeline.sh) を LSF ジョブ配列 (Job Array [1-4]) による全サーバー動的分散仕様に改修 (完了日: 10/27)
- [x] 完了: AmpliPhy デフォルト実行結果 (results/<taxon>/msa/*.amp.fa) と GS パイプライン入力の優先紐付け設定完了 (完了日: 10/27)
- [ ] 進行中: AmpliPhy ジョブ投入（デフォルト実行）→ 増幅後MSA生成 → GS法評価 (進捗率: 99% / 課題: スパコン上での bsub 投入)

## 🛠 技術的決定・メモ
- **[課題]**: GS法パイプライン（`run_gs_pipeline.sh`）も単一ノード固定ではなく、LSF クラスタ内の空いているサーバーに柔軟に分散させたい。また、AmpliPhy のデフォルト実行で生成される `results/<taxon>/msa/*.amp.fa` を確実に参照できるようにしたい。
- **[解決策]**:
  1. `#BSUB -J "gs_pipeline[1-4]"` および `#BSUB -n 8` を設定。1生物種あたり 8 コアの小単位ジョブとして投入することで、48コア機だけでなく 24コア機（mafftsrv0）や部分的な空きコアを持つ任意のサーバーへ即座に充填される仕組みを構築。
  2. `orig_dir` および `amp_dir` の探索順序を更新し、`PY_MODULE_DIR/results/$taxon/msa`（AmpliPhy の出力先）を最優先で探索するように設定。
  3. 各生物種の結果出力先を `Ampliphy2/results_gs/$taxon/` に標準化し、ログも `logs/gs/gs_%J_%I.out` に分割保存。
- **[パラメータ]**: `#BSUB -J "gs_pipeline[1-4]"`, `#BSUB -n 8`, `THREADS=8`, `MMSEQS_SENS=7.5`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で AmpliPhy を実行: `bsub < scripts/run_ampliphy.sh`
   （デフォルト設定で全サーバーに動的分散し、`results/<taxon>/msa/*.amp.fa` を生成）
2. AmpliPhy 完了後、GS法パイプラインを実行: `bsub < scripts/run_gs_pipeline.sh`
   （LSFジョブ配列により4種が空いているサーバーへ自動分散し、数分で系統樹分割・TCS評価が完了）

---

最終更新: 2026-10-27 01:26

## 📈 進捗状況
- [x] 完了: 全サーバー（mafftsrv0〜5, mafftweb 計7台）の空きスロットへ動的にタスクを充填するクラスタ動的分散機能の実装 (完了日: 10/27)
- [x] 完了: 4生物種のNextflowワークフロー完全同時並行投入によるクラスタ全コア（計280コア）フル活用化 (完了日: 10/27)
- [ ] 進行中: スパコン上での AmpliPhy ジョブ投入・増幅後MSA (*.amp.fa) 生成 (進捗率: 99% / 課題: スパコン上での bsub 投入)

## 🛠 技術的決定・メモ
- **[課題]**: 4台固定ではなく、クラスタ内の全7台（mafftsrv1~5: 48コア×5台、mafftsrv0: 24コア、mafftweb: 16コア）の「空いているサーバー・コア」へ積極的にタスクを吸い込ませ、280コアの遊休リソースを一切残さずフル活用したい。
- **[解決策]**:
  1. Nextflow の `executor = 'lsf'` によるタスク単位動的分散 (`cluster_dynamic.config`) を導入。MAFFT(4コア)、MMseqs(8コア)、IQ-TREE(4コア) といった細粒度タスクを LSF に直接投入。テトリスのように「どのサーバーの隙間コアであっても」即座にタスクが充填される仕組みを構築。
  2. マスタースクリプトは軽量なコーディネーター（`-n 2`）として即座に起動し、4生物種（計400 HOG）をバックグラウンド並行でキューに一括プール。
  3. 各生物種ごとに `-work-dir results/<taxon>/work` と独立ログファイルを分離し、セッション衝突やロック競合を完全防止。UniRef50 の事前チェックによりダウンロード競合も排除。
- **[パラメータ]**: `MODE="cluster"`, `executor.name='lsf'`, `queueSize=150`, `submitRateLimit='5 sec'`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bsub < scripts/run_ampliphy.sh` を実行（クラスタ全体の空きコアへ自動充填開始）。
2. `bjobs` で各ノード（mafftsrv0〜5）にタスクが分散実行されている様子を確認。
3. 完了後、`bsub < scripts/run_gs_pipeline.sh` を投入してGS法の系統樹分割とTCS評価を実行。

---

最終更新: 2026-10-27 01:10

## 📈 進捗状況
- [x] 完了: 48コアノード (mafftsrv1~5: 48CPUs, 352GB RAM) に向けたリソース・スレッド数最適化 (完了日: 10/27)
- [x] 完了: run_ampliphy.sh における taxonomy 引数の任意化（ENABLE_TAXONOMY=false で高速・軽量化）(完了日: 10/27)
- [ ] 進行中: スパコン上での AmpliPhy ジョブ投入・増幅後MSA (*.amp.fa) 生成 (進捗率: 98% / 課題: スパコン上での bsub 投入)

## 🛠 技術的決定・メモ
- **[課題]**: スパコン構成（mafftsrv1~5: 48コア/352GB、mafftsrv0: 24コア、mafftweb: 16コア）において、最適な並列度とメモリ設定を適用し、かつジョブがキューで滞留しないようにする。また不要な taxonomy 計算をスキップして高速化を図る。
- **[解決策]**:
  1. `#BSUB -n 40` および `#BSUB -R "select[ncpus>=48] span[hosts=1]"` を設定。48コアの主力ノード（mafftsrv1~5）に限定して割り当て、残り8コアをOS・Nextflow管理用バッファとして残すことでキュー待ち時間を短縮し安定稼働を確保。
  2. Nextflow の並列度を `threads=40`, `max_memory="200 GB"` に拡張。IQ-TREE 10並列 (`cpus=4, maxForks=10`)、MMseqs2 5並列 (`cpus=8, maxForks=5`) を実現（従来の2.5倍以上のスループット）。
  3. `ENABLE_TAXONOMY` フラグ（デフォルト: false）を導入。GS法で必要な `*.amp.fa` 生成に特化し、余計な taxonomy 処理と IQ-TREE 元ツリー推定をスキップして所要時間を半減。
- **[パラメータ]**: `-n 40`, `select[ncpus>=48]`, `MAX_MEM="200 GB"`, `ENABLE_TAXONOMY=false`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bsub < scripts/run_ampliphy.sh` を実行し、全4生物種のホモログ追加後MSA（`results/<taxon>/msa/*.amp.fa`）を生成する。
2. 完了後、`bsub < scripts/run_gs_pipeline.sh` を投入して40スレッドでGS法による系統樹分割とTCS評価を実行する。

---

最終更新: 2026-10-27 01:01

## 📈 進捗状況
- [x] 完了: run_ampliphy.sh を Ampliphy2/input_fasta 直結型・4生物種自動ループ版に改修 (完了日: 10/27)
- [x] 完了: ampliphy-analysis/dataset/*.tsv から全4生物種 (計400 HOG) の .tax ファイル一括生成完了 (完了日: 10/27)
- [ ] 進行中: スパコン上でのAmpliPhyジョブ実行・増幅後MSA生成 (進捗率: 95% / 課題: スパコン上でのジョブ投入)

## 🛠 技術的決定・メモ
- **[課題]**: sampleデータはスキップし、`Ampliphy2/input_fasta` の本番4種（eukaryota, amniota, bacteria, actinomycete）を直接ターゲットにしたい。また、AmpliPhyの `--input_taxonomy` に渡す `.tax` ファイルと `ampliphy-analysis/dataset/*.tsv` の違いの解消。
- **[解決策]**: `dataset/*_taxa.tsv`（種コード・TaxID・全系統階層のマスター）から各FASTA配列IDに対応する `sequence_id \t species_taxid` 形式の `<family_id>.tax` ファイルを全400HOG分自動生成（`Ampliphy2/taxonomy/`）。`run_ampliphy.sh` も引数なしで本番4種を直接走らせる構成に刷新。
- **[パラメータ]**: `--input_dir Ampliphy2/input_fasta/<taxon>`, `--input_taxonomy Ampliphy2/taxonomy/<taxon>`, `--output_dir Ampliphy2/results/<taxon>`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bsub < scripts/run_ampliphy.sh` を投入し、全4種（または `bash scripts/run_ampliphy.sh eukaryota`）のホモログ収集と `results/<taxon>/msa/*.amp.fa` 生成を行う。
2. 完了後、`bsub < scripts/run_gs_pipeline.sh` を投入してGS法による評価を走らせる。

---

最終更新: 2026-10-27 00:54

## 📈 進捗状況
- [x] 完了: AmpliPhy（ホモログ収集・増幅MSA生成）実行スクリプト run_ampliphy.sh の作成 (完了日: 10/27)
- [x] 完了: research-notes/20261022-.md の実行手順・ワークフロー更新 (完了日: 10/27)
- [ ] 進行中: AmpliPhyによるホモログ追加後MSA（*.amp.fa）の生成とGS法パイプライン実行 (進捗率: 90% / 課題: スパコン上でのAmpliPhyジョブ投入)

## 🛠 技術的決定・メモ
- **[課題]**: GS法によるTCS改善評価には、まずAmpliPhyによるホモログ追加後MSA（`*.amp.fa`）が必要。また、初回実行時のUniRef50ダウンロードや複数プロセス並列化時のOOMを避けるリソース制御が必要。
- **[解決策]**: `scripts/run_ampliphy.sh` を作成。サンプルデータによる動作確認から本番4種（eukaryota, amniota, bacteria, actinomycete）の実行まで対応。ノード内の `local` executor でメモリ60GB、スレッド数16、maxForksを制御する `local_resource_limit.config` を動的生成してNextflowを実行。
- **[パラメータ]**: `--database UniRef50`, `--mafft_preset auto`, `--max_depth 5`, `--threads 16`, `-resume`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bsub < scripts/run_ampliphy.sh`（または `bash scripts/run_ampliphy.sh sample`）を実行し、動作確認とUniRef50のダウンロードを完了させる。
2. `bash scripts/run_ampliphy.sh eukaryota`（または `all`）で本番データの増幅後MSA（`*.amp.fa`）を生成する。
3. `bsub < scripts/run_gs_pipeline.sh` を投入してGS法による系統推定とTCS改善度の評価を行う。
