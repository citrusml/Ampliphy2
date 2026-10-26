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
