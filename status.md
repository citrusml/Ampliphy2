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
