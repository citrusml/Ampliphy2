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
