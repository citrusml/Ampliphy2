#!/bin/bash
# ==============================================================================
# AmpliPhy (Nextflow) 実行ジョブスクリプト
# Ampliphy2/input_fasta 配下の4生物種グループ (eukaryota, amniota, bacteria, actinomycete)
# を対象とし、ホモログ収集・増幅MSA作成・系統樹推定・TCSレポート算出を実行
#
# 使い方:
#   1) 全4種を一括実行する場合 (デフォルト):
#      bsub < scripts/run_ampliphy.sh
#      または bash scripts/run_ampliphy.sh
#
#   2) 特定の生物種のみ実行する場合:
#      bash scripts/run_ampliphy.sh eukaryota
#      bash scripts/run_ampliphy.sh bacteria actinomycete
# ==============================================================================

#BSUB -J ampliphy_cluster_master
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -W 72:00
#BSUB -cwd /home/citrusml2004/Ampliphy2
#BSUB -o logs/ampliphy/master_%J.out
#BSUB -e logs/ampliphy/master_%J.err

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. パスとディレクトリの設定 (すべて /home/citrusml2004/Ampliphy2 配下に完全統一)
# ------------------------------------------------------------------------------
if [ -d "/home/citrusml2004/Ampliphy2" ]; then
    AMPLIPHY_ROOT="/home/citrusml2004/Ampliphy2"
else
    AMPLIPHY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

cd "$AMPLIPHY_ROOT"

# ampliphy リポジトリの配置 (/home/citrusml2004/ampliphy へのシンボリックリンク)
AMPLIPHY_SOURCE="/home/citrusml2004/ampliphy"
if [ ! -d "$AMPLIPHY_SOURCE" ] && [ -d "$HOME/ampliphy" ]; then
    AMPLIPHY_SOURCE="$HOME/ampliphy"
fi
if [ ! -d "$AMPLIPHY_ROOT/ampliphy" ] && [ -d "$AMPLIPHY_SOURCE" ]; then
    ln -sfn "$AMPLIPHY_SOURCE" "$AMPLIPHY_ROOT/ampliphy"
fi
AMPLIPHY_DIR="$AMPLIPHY_ROOT/ampliphy"

# ampliphy-analysis リポジトリの配置 (/home/citrusml2004/ampliphy-analysis へのシンボリックリンク)
ANALYSIS_SOURCE="/home/citrusml2004/ampliphy-analysis"
if [ ! -d "$ANALYSIS_SOURCE" ] && [ -d "$HOME/ampliphy-analysis" ]; then
    ANALYSIS_SOURCE="$HOME/ampliphy-analysis"
fi
if [ ! -d "$AMPLIPHY_ROOT/ampliphy-analysis" ] && [ -d "$ANALYSIS_SOURCE" ]; then
    ln -sfn "$ANALYSIS_SOURCE" "$AMPLIPHY_ROOT/ampliphy-analysis"
fi
DATASET_DIR="$AMPLIPHY_ROOT/ampliphy-analysis/dataset"

CONFIG_DIR="$AMPLIPHY_ROOT/configs"
NCBI_DIR="$AMPLIPHY_ROOT/ncbi_taxonomy"

# MMseqs2 データベースの保存先 (/home/citrusml2004 直下の mmseqs_db)
if [ -d "/home/citrusml2004/mmseqs_db" ]; then
    DB_DIR="/home/citrusml2004/mmseqs_db"
else
    DB_DIR="$HOME/mmseqs_db"
fi

# 必要なディレクトリ群を作成
mkdir -p "$NCBI_DIR" "$CONFIG_DIR" "$AMPLIPHY_ROOT/logs/ampliphy" "$AMPLIPHY_ROOT/results"

# ------------------------------------------------------------------------------
# リアルタイムログ出力設定 (LSFの完了バッファリングを回避し、即座にディスクへ出力)
# ------------------------------------------------------------------------------
REALTIME_LOG="$AMPLIPHY_ROOT/logs/ampliphy/ampliphy_${LSB_JOBID:-master}.log"
exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

# 実行対象生物種 (引数がなければ input_fasta 配下の4種すべて)
if [ "$#" -gt 0 ]; then
    TARGET_TAXA=("$@")
else
    TARGET_TAXA=("eukaryota" "amniota" "bacteria" "actinomycete")
fi

echo "=============================================================================="
echo "AmpliPhy Pipeline Execution (Cluster-wide Dynamic Scheduling + TCS Evaluation)"
echo "Started at:       $(date)"
echo "Host:             $(hostname)"
echo "Cluster Mode:     LSF クラスタ全体分散 (mafftsrv0〜5, mafftweb 空きコア動的充填)"
echo "AMPLIPHY_ROOT:    $AMPLIPHY_ROOT"
echo "AMPLIPHY_DIR:     $AMPLIPHY_DIR"
echo "DB_DIR:           $DB_DIR"
echo "NCBI_DIR:         $NCBI_DIR"
echo "Target Taxa:      ${TARGET_TAXA[*]}"
echo "TCS Evaluation:   有効 (taxonomy/<taxon> を使用して元/増幅後ツリーのTCSを自動算出)"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 2. 環境有効化 (biotools: Nextflow + OpenJDK 17 + MMseqs2 + MAFFT + IQ-TREE 2)
# ------------------------------------------------------------------------------
activate_env() {
    local env_name="$1"
    set +u

    # ~/.bashrc の読み込み
    if [ -f "$HOME/.bashrc" ]; then
        source "$HOME/.bashrc" || true
    fi

    # PATH に ~/.local/bin, ~/micromamba/bin を追加
    if [ -d "$HOME/.local/bin" ]; then
        export PATH="$HOME/.local/bin:$PATH"
    fi
    if [ -d "$HOME/micromamba/bin" ]; then
        export PATH="$HOME/micromamba/bin:$PATH"
    fi

    if command -v micromamba &>/dev/null; then
        eval "$(micromamba shell hook --shell bash 2>/dev/null)" || true
        micromamba activate "$env_name"
    elif command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook 2>/dev/null)" || true
        conda activate "$env_name"
    else
        echo "Error: Neither micromamba nor conda found in PATH." >&2
        exit 1
    fi

    set -u
}

activate_env "biotools"

# ------------------------------------------------------------------------------
# 3. リソース設定ファイル (Nextflow Config: クラスタ全体分散)
# ------------------------------------------------------------------------------
CONFIG_FILE="$CONFIG_DIR/cluster_dynamic.config"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "[ERROR] 設定ファイルが見つかりません: $CONFIG_FILE" >&2
    exit 1
fi
echo "Nextflow Config:  $CONFIG_FILE"

# ------------------------------------------------------------------------------
# 4. 生物種ごとの実行関数
# ------------------------------------------------------------------------------
run_taxon() {
    local taxon="$1"
    local short_name=""
    case "$taxon" in
        "eukaryota")   short_name="euk" ;;
        "amniota")     short_name="amn" ;;
        "bacteria")    short_name="bac" ;;
        "actinomycete") short_name="act" ;;
        *)             short_name="${taxon:0:3}" ;;
    esac

    echo ""
    echo "##############################################################################"
    echo "# AmpliPhy 実行開始: $taxon ($short_name)"
    echo "##############################################################################"

    # 入力FASTAディレクトリ (探索なし・Ampliphy2/input_fasta/<taxon> 固定)
    local input_dir="$AMPLIPHY_ROOT/input_fasta/$taxon"
    if [ ! -d "$input_dir" ]; then
        echo "[ERROR] $taxon の入力FASTAフォルダ ($input_dir) が見つかりません。" >&2
        return 1
    fi

    # taxonomy (.tax) ディレクトリ (常に有効・Ampliphy2/taxonomy/<taxon> 固定)
    local tax_dir="$AMPLIPHY_ROOT/taxonomy/$taxon"
    if [ ! -d "$tax_dir" ] || ! compgen -G "$tax_dir/*.tax" > /dev/null; then
        echo "[ERROR] $taxon の taxonomy ファイル群 ($tax_dir/*.tax) が見つかりません。" >&2
        return 1
    fi

    # 出力先ディレクトリ
    local output_dir="$AMPLIPHY_ROOT/results/$taxon"
    local work_dir="$output_dir/work"
    local launch_dir="$output_dir/launch"   # 生物種ごとに独立した起動ディレクトリ (.nextflow/ の競合防止)
    local nxf_log="$AMPLIPHY_ROOT/logs/ampliphy/nextflow_${taxon}.log" # 進捗ログ (tail -f で確認可)
    mkdir -p "$output_dir" "$work_dir" "$launch_dir"

    echo "  INPUT_DIR:   $input_dir"
    echo "  TAX_DIR:     $tax_dir (TCS評価有効)"
    echo "  OUTPUT_DIR:  $output_dir"
    echo "  WORK_DIR:    $work_dir"
    echo "  DB_DIR:      $DB_DIR"
    echo "  NCBI_DIR:    $NCBI_DIR"
    echo "  LOG_FILE:    $nxf_log  (tail -f で進捗確認可)"

    cd "$launch_dir"
    export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx4g}"

    # Nextflow 実行:
    # -ansi-log false は指定せず、Nextflow 標準のリアルタイム進捗ブロックを直接ログへストリーミング出力
    # --input_taxonomy を指定することで、homolog_taxonomy, amplified_taxonomy,
    # input_lineage, amplified_tcs が自動実行され amplified_tcs_report.tsv が生成されます。
    # shellcheck disable=SC2086
    nextflow run "$AMPLIPHY_DIR/ampliphy.nf" \
        -work-dir "$work_dir" \
        -c "$CONFIG_FILE" \
        --input_dir "$input_dir" \
        --input_taxonomy "$tax_dir" \
        --ncbi_dir "$NCBI_DIR" \
        --output_dir "$output_dir" \
        --database UniRef50 \
        --database_dir "$DB_DIR" \
        --threads 4 \
        --max_depth 5 \
        --mafft_preset auto \
        --mad_script "$AMPLIPHY_DIR/bin/mad.py" \
        --tcs_script "$AMPLIPHY_DIR/bin/tcs.py" \
        -resume > "$nxf_log" 2>&1

    echo "AmpliPhy 完了: $taxon"
    echo "  -> 増幅後MSA: $output_dir/msa/*.amp.fa"
    echo "  -> 系統樹:   $output_dir/tree/*.amp.nwk"
    echo "  -> TCSレポ:  $output_dir/report/amplified_tcs_report.tsv"
}

# ------------------------------------------------------------------------------
# 5. メインループ実行 (クラスタ全体分散投入)
# ------------------------------------------------------------------------------
# UniRef50 データベースの事前準備 (並行実行時の重複ダウンロード・競合防止)
if [ ! -e "$DB_DIR/UniRef50.dbtype" ] && [ ! -e "$DB_DIR/UniRef50" ]; then
    DB_LOG="$AMPLIPHY_ROOT/logs/ampliphy/uniref50_download.log"
    echo "UniRef50 データベースが未準備のため、初回ダウンロードを実行します (数時間かかる場合があります)"
    echo "  進捗: tail -f $DB_LOG"
    mkdir -p "$AMPLIPHY_ROOT/tmp"
    mmseqs databases UniRef50 "$DB_DIR/UniRef50" "$AMPLIPHY_ROOT/tmp" > "$DB_LOG" 2>&1
    echo "UniRef50 データベースの準備が完了しました。"
else
    echo "UniRef50 データベース検出済み: $DB_DIR/UniRef50"
fi

# NCBI Taxonomy dump (nodes.dmp, names.dmp) の事前準備 (並行タスク間の重複DL・競合防止)
if [ ! -f "$NCBI_DIR/nodes.dmp" ] || [ ! -f "$NCBI_DIR/names.dmp" ]; then
    echo "NCBI Taxonomy データを準備中 ($NCBI_DIR)..."
    python3 -c "import sys; sys.path.insert(0, '$AMPLIPHY_DIR/bin'); from ampliphy_homolog_taxonomy import ensure_taxdump; from pathlib import Path; ensure_taxdump(Path('$NCBI_DIR'))" || true
    echo "NCBI Taxonomy データの準備が完了しました。"
else
    echo "NCBI Taxonomy 検出済み: $NCBI_DIR"
fi

echo ""
echo "=============================================================================="
echo "クラスタ全体分散投入: 対象生物種 (${TARGET_TAXA[*]})"
echo "各生物種のNextflowタスク (MAFFT, MMseqs2, IQ-TREE, TCS) が全サーバー (mafftsrv0〜5, mafftweb) へ自動分散されます。"
echo "=============================================================================="

if [ "${#TARGET_TAXA[@]}" -gt 1 ]; then
    pids=()
    for taxon in "${TARGET_TAXA[@]}"; do
        run_taxon "$taxon" &
        pids+=($!)
    done

    failed=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed=$((failed + 1))
        fi
    done

    if [ "$failed" -gt 0 ]; then
        echo "[ERROR] $failed 件の生物種タスクでエラーが発生しました。ログを確認してください。" >&2
        exit 1
    fi
else
    run_taxon "${TARGET_TAXA[0]}"
fi

echo ""
echo "=============================================================================="
echo "すべての AmpliPhy パイプライン実行が完了しました: $(date)"
echo "各 results/<taxon>/ 配下に増幅後MSA (*.amp.fa)、系統樹 (*.amp.nwk)、TCSレポート (report/amplified_tcs_report.tsv) が揃いました。"
echo "続いて GS法評価 (bsub < scripts/run_gs_pipeline.sh) を実行できます。"
echo "=============================================================================="
