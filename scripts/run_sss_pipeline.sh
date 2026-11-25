#!/bin/bash
# ==============================================================================
# SSS (Symmetrical Sequence Similarity score) 系統樹作成＆配列増幅精度評価 パイプライン
# (LSF クラスタ全体分散 Nextflow 版)
#
# 対象生物種: 4種 (eukaryota, amniota, bacteria, actinomycete)
# 手法:
#   1) sss_linear : d_ij = 1.0 - w_ij
#   2) sss_log    : d_ij = -ln(max(w_ij, 1e-4))
#
# 使い方:
#   1) 全4種 × 2手法 を一括実行する場合 (デフォルト):
#      bsub < scripts/run_sss_pipeline.sh
#      または bash scripts/run_sss_pipeline.sh
#
#   2) 特定の生物種のみ実行する場合:
#      bash scripts/run_sss_pipeline.sh eukaryota
#      bash scripts/run_sss_pipeline.sh eukaryota amniota
#
#   3) 特定の手法のみ実行する場合:
#      MODELS="linear" bash scripts/run_sss_pipeline.sh eukaryota
#      MODELS="log" bsub < scripts/run_sss_pipeline.sh
# ==============================================================================

#BSUB -J sss_pipeline
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -W 48:00
#BSUB -cwd /home/citrusml2004/Ampliphy2
#BSUB -o logs/sss/master_%J.out
#BSUB -e logs/sss/master_%J.err

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. パスとディレクトリの設定
# ------------------------------------------------------------------------------
if [ -d "/home/citrusml2004/Ampliphy2" ]; then
    AMPLIPHY_ROOT="/home/citrusml2004/Ampliphy2"
else
    AMPLIPHY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

cd "$AMPLIPHY_ROOT"

# ampliphy リポジトリの配置
AMPLIPHY_DIR="$AMPLIPHY_ROOT/ampliphy"
if [ ! -d "$AMPLIPHY_DIR" ]; then
    if [ -d "/home/citrusml2004/ampliphy" ]; then
        ln -sfn "/home/citrusml2004/ampliphy" "$AMPLIPHY_DIR"
    elif [ -d "$HOME/ampliphy" ]; then
        ln -sfn "$HOME/ampliphy" "$AMPLIPHY_DIR"
    fi
fi

# ampliphy-analysis リポジトリの配置
DATASET_DIR="$AMPLIPHY_ROOT/ampliphy-analysis/dataset"
if [ ! -d "$DATASET_DIR" ]; then
    if [ -d "/home/citrusml2004/ampliphy-analysis/dataset" ]; then
        ln -sfn "/home/citrusml2004/ampliphy-analysis" "$AMPLIPHY_ROOT/ampliphy-analysis"
        DATASET_DIR="$AMPLIPHY_ROOT/ampliphy-analysis/dataset"
    elif [ -d "$HOME/ampliphy-analysis/dataset" ]; then
        ln -sfn "$HOME/ampliphy-analysis" "$AMPLIPHY_ROOT/ampliphy-analysis"
        DATASET_DIR="$AMPLIPHY_ROOT/ampliphy-analysis/dataset"
    fi
fi

RESULTS_BASE_DIR="${RESULTS_DIR:-$AMPLIPHY_ROOT/results_sss}"
CONFIG_FILE="$AMPLIPHY_ROOT/configs/cluster_dynamic_sss.config"
PIPELINE_NF="$AMPLIPHY_ROOT/configs/sss_pipeline/main.nf"

mkdir -p "$AMPLIPHY_ROOT/logs/sss" "$RESULTS_BASE_DIR"

# ------------------------------------------------------------------------------
# リアルタイムログ出力設定
# ------------------------------------------------------------------------------
REALTIME_LOG="$AMPLIPHY_ROOT/logs/sss/sss_pipeline_${LSB_JOBID:-master}.log"
exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

# 実行対象生物種 (引数がなければ input_fasta 配下の4種すべて)
if [ "$#" -gt 0 ]; then
    TARGET_TAXA=("$@")
else
    TARGET_TAXA=("eukaryota" "amniota" "bacteria" "actinomycete")
fi

# 実行対象モデル ('linear', 'log', または 'linear log')
TARGET_MODELS=(${MODELS:-"linear" "log"})

echo "=============================================================================="
echo "AmpliPhy SSS (Symmetrical Sequence Similarity) Nextflow Pipeline"
echo "Started at:       $(date)"
echo "Host:             $(hostname)"
echo "AMPLIPHY_ROOT:    $AMPLIPHY_ROOT"
echo "Target Taxa:      ${TARGET_TAXA[*]}"
echo "Target Models:    ${TARGET_MODELS[*]}"
echo "Output Base Dir:  $RESULTS_BASE_DIR"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 2. 環境有効化 (biotools)
# ------------------------------------------------------------------------------
activate_env() {
    local env_name="$1"
    set +u

    if [ -f "$HOME/.bashrc" ]; then
        source "$HOME/.bashrc" 2>/dev/null || true
    fi

    if [ -d "$HOME/.local/bin" ]; then
        export PATH="$HOME/.local/bin:$PATH"
    fi
    if [ -d "$HOME/micromamba/bin" ]; then
        export PATH="$HOME/micromamba/bin:$PATH"
    fi

    if command -v micromamba &>/dev/null; then
        eval "$(micromamba shell hook --shell bash 2>/dev/null)" || true
        micromamba activate "$env_name" 2>/dev/null || true
    elif command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook 2>/dev/null)" || true
        conda activate "$env_name" 2>/dev/null || true
    fi

    set -u
}

activate_env "biotools"
export PYTHONUNBUFFERED=1
export PYTHONPATH="$AMPLIPHY_ROOT:${PYTHONPATH:-}"

# スクリプト類の検出
MAD_PY=""
for cand in \
    "$AMPLIPHY_DIR/bin/mad.py" \
    "/home/citrusml2004/ampliphy/bin/mad.py" \
    "$HOME/ampliphy/bin/mad.py" \
    "$HOME/tools/mad2.2/mad/mad.py"; do
    if [ -f "$cand" ]; then
        MAD_PY="$cand"
        break
    fi
done

TCS_PY=""
for cand in \
    "$AMPLIPHY_DIR/bin/tcs.py" \
    "/home/citrusml2004/ampliphy/bin/tcs.py" \
    "$HOME/ampliphy/bin/tcs.py"; do
    if [ -f "$cand" ]; then
        TCS_PY="$cand"
        break
    fi
done

if [ -z "$MAD_PY" ]; then
    echo "[ERROR] mad.py スクリプトが見つかりません。" >&2
    exit 1
fi
if [ -z "$TCS_PY" ]; then
    echo "[ERROR] tcs.py スクリプトが見つかりません。" >&2
    exit 1
fi

echo "MAD スクリプト:       $MAD_PY"
echo "TCS スクリプト:       $TCS_PY"

# ------------------------------------------------------------------------------
# 3. 実行関数
# ------------------------------------------------------------------------------
run_taxon_model() {
    local taxon="$1"
    local model="$2"

    local short_name=""
    case "$taxon" in
        "eukaryota")    short_name="euk" ;;
        "amniota")      short_name="amn" ;;
        "bacteria")     short_name="bac" ;;
        "actinomycete") short_name="act" ;;
        *)              short_name="${taxon:0:3}" ;;
    esac

    echo ""
    echo "=============================================================================="
    echo "# 実行開始: $taxon ($short_name) | 手法: sss_${model}"
    echo "=============================================================================="

    local orig_dir="$AMPLIPHY_ROOT/input_fasta/$taxon"
    if [ ! -d "$orig_dir" ]; then
        echo "[ERROR] $taxon の入力FASTAフォルダ ($orig_dir) が見つかりません。" >&2
        return 1
    fi

    local hom_dir="$AMPLIPHY_ROOT/results/$taxon/homologs"
    if [ ! -d "$hom_dir" ] && [ -d "$AMPLIPHY_ROOT/results_max/$taxon/homologs" ]; then
        hom_dir="$AMPLIPHY_ROOT/results_max/$taxon/homologs"
    elif [ ! -d "$hom_dir" ] && [ -d "$AMPLIPHY_ROOT/results_denovo/$taxon/homologs" ]; then
        hom_dir="$AMPLIPHY_ROOT/results_denovo/$taxon/homologs"
    fi

    if [ ! -d "$hom_dir" ]; then
        echo "[ERROR] $taxon のホモログフォルダ ($hom_dir) が見つかりません。" >&2
        return 1
    fi

    local out_dir="$RESULTS_BASE_DIR/$taxon/$model"
    local work_dir="$out_dir/work"
    local launch_dir="$out_dir/launch"
    local nxf_log="$AMPLIPHY_ROOT/logs/sss/nextflow_${taxon}_${model}.log"
    mkdir -p "$out_dir" "$work_dir" "$launch_dir"

    # Lineage TSV
    local lineage_tsv="$RESULTS_BASE_DIR/$taxon/lineage.tsv"
    if [ ! -s "$lineage_tsv" ]; then
        local taxa_tsv="$DATASET_DIR/${short_name}_taxa.tsv"
        if [ -s "$taxa_tsv" ]; then
            echo "  [Lineage] $taxa_tsv から lineage.tsv を生成中..."
            python3 -c "
from src.pwa_pipeline import generate_lineage_tsv
import glob
fastas = sorted(glob.glob('$orig_dir/*.fa'))
generate_lineage_tsv('$taxa_tsv', fastas, '$lineage_tsv')
"
        elif [ -s "$DATASET_DIR/${short_name}_lineage.tsv" ]; then
            cp -f "$DATASET_DIR/${short_name}_lineage.tsv" "$lineage_tsv"
        fi
    fi

    cd "$launch_dir"
    export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx4g}"

    nextflow run "$PIPELINE_NF" \
        -work-dir "$work_dir" \
        -c "$CONFIG_FILE" \
        --input_dir "$orig_dir" \
        --homolog_dir "$hom_dir" \
        --lineage_file "$lineage_tsv" \
        --output_dir "$out_dir" \
        --model "$model" \
        --mad_script "$MAD_PY" \
        --tcs_script "$TCS_PY" \
        -resume > "$nxf_log" 2>&1

    echo "完了: $taxon (sss_${model})"
    echo "  -> 増幅後系統樹 (amp):      $out_dir/tree_amp/*.amp.nwk"
    echo "  -> 根付け後樹 (amp.rooted): $out_dir/tree_rooted/*.amp.rooted.nwk"
    echo "  -> 枝刈り後樹 (enr):        $out_dir/tree_enr/*.enr.nwk"
    echo "  -> TCSレポート:             $out_dir/report/tcs_report.tsv"
}

# ------------------------------------------------------------------------------
# 4. メインループ実行
# ------------------------------------------------------------------------------
for taxon in "${TARGET_TAXA[@]}"; do
    for model in "${TARGET_MODELS[@]}"; do
        run_taxon_model "$taxon" "$model"
    done
done

echo ""
echo "=============================================================================="
echo "すべての SSS (linear / log) パイプライン実行が完了しました: $(date)"
echo "出力先: $RESULTS_BASE_DIR/<taxon>/<model>/"
echo "=============================================================================="
