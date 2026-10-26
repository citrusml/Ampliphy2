最終更新: 2026-10-27 00:25

## 📈 進捗状況
- [x] 完了: GS法パイプラインの4分類群（Eukaryota, Amniota, Bacteria, Actinomycete）一括実行対応 (完了日: 10/27)
- [x] 完了: 各種メタデータからのlineage.tsv自動生成モジュール実装・単体テスト (完了日: 10/27)
- [ ] 進行中: スパコン上での4分類群一括実行と改善度集計 (進捗率: 85% / 課題: スパコン上でのジョブ投入)

## 🛠 技術的決定・メモ
- **[課題]**: Eukaryotaだけでなく、Amniota, Bacteria, Actinomycete の全4種を一括で評価したい。また、分類群によって lineage.tsv が未生成の場合がある。
- **[解決策]**: `run_gs_pipeline.sh` をマルチ生物種対応に拡張。引数無しの場合は4種を順次自動実行し、引数で特定種のみの実行も可能にした。さらに `dataset/{short}_taxa.tsv` から各配列の `lineage.tsv` を動的自動生成する処理を `src/gs_pipeline.py` に実装し、`pytest` で検証完了。
- **[パラメータ]**: 対象生物種: `eukaryota`, `amniota`, `bacteria`, `actinomycete`。各100HOG。

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bsub < scripts/run_gs_pipeline.sh` を投入し、4分類群すべて（計400 HOG）の系統推定とTCS改善度を算出する。
2. 出力された各分類群の `gs_tcs_summary.txt` および `gs_tcs_comparison.tsv` を確認する。

---

最終更新: 2026-10-27 00:18

## 📈 進捗状況
- [x] 完了: LSF非対話型実行時における仮想環境有効化エラー（activate_env）の修正 (完了日: 10/27)
- [ ] 進行中: スパコン上でのGS法パイプライン実行・動作確認 (進捗率: 80% / 課題: エラーログの確認と再投入)

## 🛠 技術的決定・メモ
- **[課題]**: LSFバッチジョブ（非対話型シェル）では `~/.bashrc` が自動読み込みされず、`micromamba` コマンドが未検出になる、もしくは `set -u` の影響で conda/micromamba のシェル関数内部の未定義変数チェック（`PS1` 等）に引っかかりジョブ開始直後に Exit 1 で即終了する問題。
- **[解決策]**: `activate_env` 関数内で一時的に `set +u` を適用し、`~/.bashrc` の明示的 source および `~/.local/bin` / `~/micromamba/bin` の PATH 追加を安全に行うよう改修。

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `cat Ampliphy2/logs/gs/gs_pipeline_246023.err` を確認してエラー内容を照合する。
2. 修正した `scripts/run_gs_pipeline.sh` をスパコン側に反映し、再度 `bsub < scripts/run_gs_pipeline.sh` を投入する。

---

最終更新: 2026-10-27 00:13

## 📈 進捗状況
- [x] 完了: AmpliPhy新旧差分の調査・ChangeAmpliphy.md作成 (完了日: 10/26)
- [x] 完了: TCSスコア計算アルゴリズムの数理・実装差異の整理 (完了日: 10/26)
- [x] 完了: GS法実行パイプラインのパス解決・Pythonモジュール参照強化 (完了日: 10/27)
- [ ] 進行中: スパコン上でのGS法パイプライン実行・動作確認 (進捗率: 70% / 課題: 入力データの配置確認)

## 🛠 技術的決定・メモ
- **[課題]**: `run_gs_pipeline.sh` をホームディレクトリや任意ディレクトリから `bsub` 投入した際、`PY_MODULE_DIR="$PWD"` だと `from src.gs_pipeline import ...` で ModuleNotFoundError になるリスク。
- **[解決策]**: `SCRIPT_DIR` から親ディレクトリ（`Ampliphy2/`）の絶対パスを動的に取得して `PY_MODULE_DIR` に設定し、`PYTHONPATH` を確実に通すように修正。また、LSFログ出力用の `Ampliphy2/logs/gs/` ディレクトリを作成。
- **[パラメータ]**: GS2感度 `MMSEQS_SENS=7.5`, スレッド数 `THREADS=16`, ジョブ要求 `span[hosts=1]`, 時間制限 `24:00`

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `mkdir -p Ampliphy2/logs/gs` があることを確認し、`bsub < Ampliphy2/scripts/run_gs_pipeline.sh`（または `Ampliphy2` ディレクトリ内で実行）を投入する。
2. 入力ディレクトリ (`input_fasta/eukaryota` および `results/eukaryota_new/msa`) に必要な FASTA ファイルが配置されていることを確認する。
