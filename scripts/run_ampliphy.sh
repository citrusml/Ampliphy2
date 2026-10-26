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
elif [ -d "$HOME/Ampliphy2" ]; then
    AMPLIPHY_ROOT="$HOME/Ampliphy2"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    AMPLIPHY_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

# スパコン上での実行カレントディレクトリを確実に Ampliphy2 に設定
cd "$AMPLIPHY_ROOT"
WORK_DIR="$AMPLIPHY_ROOT"

AMPLIPHY_DIR="$AMPLIPHY_ROOT/ampliphy"
DATASET_DIR="$AMPLIPHY_ROOT/ampliphy-analysis/dataset"

# MMseqs2 データベースの保存先 (Ampliphy2/mmseqs_db または既存の mmseqs_db)
DB_DIR="$AMPLIPHY_ROOT/mmseqs_db"
if [ ! -d "$DB_DIR" ] && [ -d "$HOME/mmseqs_db" ]; then
    DB_DIR="$HOME/mmseqs_db"
fi
mkdir -p "$DB_DIR" "$AMPLIPHY_ROOT/logs/ampliphy" "$AMPLIPHY_ROOT/results"

# ------------------------------------------------------------------------------
# リアルタイムログ出力設定 (LSFの完了バッファリングを回避し、即座にディスクへ出力)
# ------------------------------------------------------------------------------
REALTIME_LOG="$AMPLIPHY_ROOT/logs/ampliphy/ampliphy_live_${LSB_JOBID:-master}.log"
exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

# ------------------------------------------------------------------------------
# 実行モード設定:
#   MODE="cluster" (デフォルト・推奨):
#     Nextflow が個別の計算タスク (MAFFT, MMseqs2, IQ-TREE) を LSF に直接投入。
#     空いているすべてのサーバー (mafftsrv0〜5, mafftweb) の空きコアへ自動的に充填されます。
#   MODE="local":
#     単一ノード内 (40コア) でローカル並列実行。
# ------------------------------------------------------------------------------
MODE="${MODE:-cluster}"

# 各プロセスの並列スレッド・メモリ割当設定
THREADS=40
MAX_MEM="200 GB"

# Taxonomy 評価を行うかどうか (デフォルト: false = 旧版同様にオフにして高速化)
# GS法パイプライン (run_gs_pipeline.sh) で必要なのは *.amp.fa のみなため通常は false で十分です。
ENABLE_TAXONOMY="${ENABLE_TAXONOMY:-false}"

# 実行対象生物種 (引数がなければ input_fasta 配下の4種すべて)
if [ "$#" -gt 0 ]; then
    TARGET_TAXA=("$@")
else
    TARGET_TAXA=("eukaryota" "amniota" "bacteria" "actinomycete")
fi

echo "=============================================================================="
echo "AmpliPhy Pipeline Execution (Cluster-wide Dynamic Scheduling)"
echo "Started at:       $(date)"
echo "Host:             $(hostname)"
echo "Execution Mode:   $MODE (全サーバー空きスロット動的充填)"
echo "ENABLE_TAXONOMY:  $ENABLE_TAXONOMY"
echo "AMPLIPHY_ROOT:    $AMPLIPHY_ROOT"
echo "DB_DIR:           $DB_DIR"
echo "Target Taxa:      ${TARGET_TAXA[*]}"
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
# 3. リソース設定ファイル生成 (Nextflow Config)
# ------------------------------------------------------------------------------
if [ "$MODE" = "cluster" ]; then
    # 【クラスタ全体分散モード】
    # 各タスク (MAFFT=4, MMseqs=8, IQ-TREE=4) を個別ジョブとして LSF に投入。
    # 空きコアのある任意のサーバー (mafftsrv0〜5, mafftweb) へ自動ディスパッチ。
    CONFIG_FILE="$AMPLIPHY_DIR/cluster_dynamic.config"
    cat << 'EOF' > "$CONFIG_FILE"
executor {
    name            = 'lsf'
    queueSize       = 150
    pollInterval    = '10 sec'
    submitRateLimit = '5 sec'
}

process {
    executor       = 'lsf'
    clusterOptions = '-R "span[hosts=1]"'
    beforeScript   = '''
        set +u
        source ~/.bashrc 2>/dev/null || true
        export PATH="$HOME/.local/bin:$HOME/micromamba/bin:$PATH"
        if command -v micromamba &>/dev/null; then
            eval "$(micromamba shell hook --shell bash 2>/dev/null)" || true
            micromamba activate biotools 2>/dev/null || true
        elif command -v conda &>/dev/null; then
            eval "$(conda shell.bash hook 2>/dev/null)" || true
            conda activate biotools 2>/dev/null || true
        fi
        set -u
    '''

    withLabel: 'mmseqs' {
        cpus   = 8
        memory = 16.GB
        time   = 6.h
    }
    withLabel: 'tree' {
        cpus   = 4
        memory = 16.GB
        time   = 12.h
    }
    withLabel: 'short' {
        cpus   = 4
        memory = 8.GB
        time   = 2.h
    }
    withLabel: 'minimal' {
        cpus   = 1
        memory = 4.GB
        time   = 1.h
    }
}
EOF
else
    # 【単一ノード内実行モード】
    CONFIG_FILE="$AMPLIPHY_DIR/local_resource_limit.config"
    cat << EOF > "$CONFIG_FILE"
executor {
    name   = 'local'
    cpus   = ${THREADS}
    memory = '${MAX_MEM}'
}

process {
    withLabel: 'mmseqs' {
        cpus     = 8
        maxForks = 5
    }
    withLabel: 'tree' {
        cpus     = 4
        maxForks = 10
    }
    withLabel: 'short' {
        cpus     = 4
        maxForks = 10
    }
}
EOF
fi

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

    # input_fasta ディレクトリ探索
    local input_dir=""
    for cand in \
        "$AMPLIPHY_ROOT/input_fasta/$taxon" \
        "$WORK_DIR/Ampliphy2/input_fasta/$taxon" \
        "$WORK_DIR/input_fasta/$taxon"; do
        if [ -d "$cand" ] && compgen -G "$cand/*.fa" > /dev/null; then
            input_dir="$cand"
            break
        fi
    done

    if [ -z "$input_dir" ]; then
        echo "[ERROR] $taxon の入力FASTAフォルダ (input_fasta/$taxon) が見つかりません。" >&2
        return 1
    fi

    # taxonomy (.tax) ディレクトリの準備 (ENABLE_TAXONOMY=true の場合のみ)
    local tax_arg=""
    if [ "$ENABLE_TAXONOMY" = "true" ]; then
        local tax_dir="$AMPLIPHY_ROOT/taxonomy/$taxon"
        mkdir -p "$tax_dir"

        if ! compgen -G "$tax_dir/*.tax" > /dev/null; then
            local meta_tsv="$DATASET_DIR/${short_name}_taxa.tsv"
            if [ -s "$meta_tsv" ]; then
                echo "  [Taxonomy] ${short_name}_taxa.tsv から各配列の .tax ファイルを生成中..."
                python3 -c "
import glob
from pathlib import Path

taxa_map = dict(line.strip().split('\t')[:2] for line in open('$meta_tsv') if line.strip())
fa_files = sorted(glob.glob('$input_dir/*.fa'))
out_dir = Path('$tax_dir')
for fa in fa_files:
    tax_path = out_dir / (Path(fa).stem + '.tax')
    with open(tax_path, 'w', encoding='utf-8') as out_f, open(fa, 'r', encoding='utf-8') as in_f:
        for line in in_f:
            if line.startswith('>'):
                qid = line[1:].strip().split()[0]
                code = qid[:5]
                if code in taxa_map:
                    out_f.write(f'{qid}\t{taxa_map[code]}\n')
print(f'  [Taxonomy] {len(fa_files)} 件の .tax ファイルを $tax_dir に生成しました。')
"
            fi
        fi
        tax_arg="--input_taxonomy $tax_dir"
        echo "  TAX_DIR:    $tax_dir (有効)"
    else
        echo "  TAXONOMY:   無効 (旧版同様にホモログ追加MSA作成に専念して高速化)"
    fi

    # 出力先ディレクトリ
    local output_dir="$AMPLIPHY_ROOT/results/$taxon"
    local work_dir="$output_dir/work"
    local launch_dir="$output_dir/launch"   # 生物種ごとに独立した起動ディレクトリ (.nextflow/ の競合防止)
    local log_file="$AMPLIPHY_ROOT/logs/ampliphy/nextflow_${taxon}.log"          # .nextflow.log 相当
    local console_log="$AMPLIPHY_ROOT/logs/ampliphy/nextflow_${taxon}.console.log" # 進捗表示 (tail -f 用)
    mkdir -p "$output_dir" "$work_dir" "$launch_dir"

    echo "  INPUT_DIR:   $input_dir"
    echo "  OUTPUT_DIR:  $output_dir"
    echo "  WORK_DIR:    $work_dir"
    echo "  DB_DIR:      $DB_DIR"
    echo "  CONSOLE_LOG: $console_log  (tail -f で確認可)"

    cd "$launch_dir"
    export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx4g}"

    # Nextflow の標準出力を直接ファイルへリダイレクトする (LSF の -o バッファを経由しない)。
    # -profile conda は使わない: biotools 環境 (および LSF タスクの beforeScript) を使用する。
    # shellcheck disable=SC2086
    nextflow -log "$log_file" run "$AMPLIPHY_DIR/ampliphy.nf" \
        -ansi-log false \
        -work-dir "$work_dir" \
        -c "$CONFIG_FILE" \
        --input_dir "$input_dir" \
        $tax_arg \
        --output_dir "$output_dir" \
        --database UniRef50 \
        --database_dir "$DB_DIR" \
        --threads 4 \
        --max_depth 5 \
        --mafft_preset auto \
        --mad_script "$AMPLIPHY_DIR/bin/mad.py" \
        --tcs_script "$AMPLIPHY_DIR/bin/tcs.py" \
        -resume > "$console_log" 2>&1

    echo "AmpliPhy 完了: $taxon"
    echo "  -> 増幅後MSA: $output_dir/msa/*.amp.fa"
    echo "  -> 系統樹:   $output_dir/tree/*.amp.nwk"
    if [ "$ENABLE_TAXONOMY" = "true" ]; then
        echo "  -> TCSレポ:  $output_dir/report/amplified_tcs_report.tsv"
    fi
}

# ------------------------------------------------------------------------------
# 5. メインループ実行 (全生物種を同時投入してクラスタ全体をフル稼働)
# ------------------------------------------------------------------------------
# UniRef50 データベースの事前準備 (並行実行時の重複ダウンロード・競合防止)
if [ ! -e "$DB_DIR/UniRef50.dbtype" ] && [ ! -e "$DB_DIR/UniRef50" ]; then
    DB_LOG="$AMPLIPHY_ROOT/logs/ampliphy/uniref50_download.log"
    echo "UniRef50 データベースが未準備のため、初回ダウンロードを実行します (数時間かかる場合があります)"
    echo "  進捗: tail -f $DB_LOG"
    mkdir -p "$WORK_DIR/tmp"
    mmseqs databases UniRef50 "$DB_DIR/UniRef50" "$WORK_DIR/tmp" > "$DB_LOG" 2>&1
    echo "UniRef50 データベースの準備が完了しました。"
else
    echo "UniRef50 データベース検出済み: $DB_DIR/UniRef50"
fi

if [ "$MODE" = "cluster" ] && [ "${#TARGET_TAXA[@]}" -gt 1 ]; then
    echo ""
    echo "=============================================================================="
    echo "【クラスタ動的分散モード】全生物種 (${TARGET_TAXA[*]}) を同時並行投入します。"
    echo "空いている全サーバー (mafftsrv0〜5, mafftweb) の空きコアへ自動的にタスクが分散されます。"
    echo "=============================================================================="
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
    for taxon in "${TARGET_TAXA[@]}"; do
        run_taxon "$taxon"
    done
fi

echo ""
echo "=============================================================================="
echo "すべての AmpliPhy パイプライン実行が完了しました: $(date)"
echo "これで各 results/<taxon>/msa/ 配下にホモログ追加後MSA (*.amp.fa) が揃いました。"
echo "続いて GS法評価 (bsub < scripts/run_gs_pipeline.sh) を実行できます。"
echo "=============================================================================="
