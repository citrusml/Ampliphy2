#!/bin/bash
# ==============================================================================
# E-value × addfragments (keeplength有無) × MaxAlign (有無) × 4系統樹推定手法
# ベンチマークパイプライン (LSF クラスタ全体分散 Nextflow 版)
#
# 対象生物種: eukaryota, bacteria
# E-value 閾値: 1e-100, 1e-10, 1e-3, 1e-1, 1e-0 (固定シード 42 でランダムサンプリング)
# MSA バリアント:
#   - maxa_on_keep_on   (MaxAlignあり + addfragments + keeplengthあり)
#   - maxa_on_keep_off  (MaxAlignあり + addfragments + keeplengthなし)
#   - maxa_off_keep_on  (MaxAlignなし + addfragments + keeplengthあり)
#   - maxa_off_keep_off (MaxAlignなし + addfragments + keeplengthなし)
# 系統樹推定手法:
#   - msa_ml         : IQ-TREE 2 (-m JTT+I+G4 -B 1000 -keep-ident) + MAD + GoTree
#   - fastme_poisson : MSA Poisson距離行列 + FastME (-m N -n -s) + MAD + GoTree
#   - msa_nj         : MSA Poisson距離行列 + RapidNJ (-i pd) + MAD + GoTree
#   - gs_poisson     : MSA Poisson距離行列 + GS2 (-d -s -l) + GoTree
#
# 使い方:
#   1) 全条件一括実行 (eukaryota & bacteria):
#      bsub < scripts/run_evalue_benchmark.sh
#   2) 特定の生物種・条件のみ実行:
#      bash scripts/run_evalue_benchmark.sh eukaryota
#      EVALUES="1e-3" VARIANTS="maxa_on_keep_on" bash scripts/run_evalue_benchmark.sh eukaryota
# ==============================================================================

#BSUB -J eval_bench
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -W 72:00
#BSUB -cwd /home/citrusml2004/Ampliphy2
#BSUB -o logs/evalue_bench/master_%J.out
#BSUB -e logs/evalue_bench/master_%J.err

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

# ampliphy リポジトリの配置 & 最新スクリプト (ランダムサンプリング版等) の自動同期
AMPLIPHY_DIR="$AMPLIPHY_ROOT/ampliphy"
if [ ! -d "$AMPLIPHY_DIR" ]; then
    if [ -d "/home/citrusml2004/ampliphy" ]; then
        ln -sfn "/home/citrusml2004/ampliphy" "$AMPLIPHY_DIR"
    elif [ -d "$HOME/ampliphy" ]; then
        ln -sfn "$HOME/ampliphy" "$AMPLIPHY_DIR"
    fi
fi

if [ -d "$AMPLIPHY_ROOT/configs/ampliphy_nf/bin" ] && [ -d "$AMPLIPHY_DIR/bin" ]; then
    cp -f "$AMPLIPHY_ROOT/configs/ampliphy_nf/bin/"* "$AMPLIPHY_DIR/bin/" 2>/dev/null || true
    chmod +x "$AMPLIPHY_DIR/bin/"* 2>/dev/null || true
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

if [ -d "/home/citrusml2004/mmseqs_db" ]; then
    DB_DIR="/home/citrusml2004/mmseqs_db"
else
    DB_DIR="$HOME/mmseqs_db"
fi

RESULTS_BASE_DIR="${RESULTS_DIR:-$AMPLIPHY_ROOT/results_evalue_bench}"
CONFIG_FILE="$AMPLIPHY_ROOT/configs/cluster_dynamic_evalue.config"
PIPELINE_NF="$AMPLIPHY_ROOT/configs/evalue_benchmark/main.nf"

mkdir -p "$AMPLIPHY_ROOT/logs/evalue_bench" "$RESULTS_BASE_DIR"

# ------------------------------------------------------------------------------
# リアルタイムログ出力設定
# ------------------------------------------------------------------------------
REALTIME_LOG="$AMPLIPHY_ROOT/logs/evalue_bench/evalue_bench_${LSB_JOBID:-master}.log"
exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

if [ "$#" -gt 0 ]; then
    TARGET_TAXA=("$@")
else
    TARGET_TAXA=("eukaryota" "bacteria")
fi

EVALUES="${EVALUES:-1e-100,1e-10,1e-3,1e-1,1e-0}"
VARIANTS="${VARIANTS:-maxa_on_keep_on,maxa_on_keep_off,maxa_off_keep_on,maxa_off_keep_off}"
RANDOM_SEED="${RANDOM_SEED:-42}"

echo "=============================================================================="
echo "AmpliPhy E-value × addfragments × MaxAlign Benchmark Pipeline"
echo "Started at:       $(date)"
echo "Host:             $(hostname)"
echo "AMPLIPHY_ROOT:    $AMPLIPHY_ROOT"
echo "DB_DIR:           $DB_DIR"
echo "Target Taxa:      ${TARGET_TAXA[*]}"
echo "E-values:         $EVALUES (Random Seed: $RANDOM_SEED)"
echo "Variants:         $VARIANTS"
echo "Methods:          msa_ml, fastme_poisson (-n -s), msa_nj (RapidNJ), gs_poisson (gs2 -d), psa_nj_nw, psa_fastme_nw, psa_nj_sss_linear, psa_fastme_sss_linear"
echo "Output Base Dir:  $RESULTS_BASE_DIR"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 2. 環境有効化 (biotools) & 必要ツール確認
# ------------------------------------------------------------------------------
activate_env() {
    local env_name="$1"
    set +u
    [ -f "$HOME/.bashrc" ] && source "$HOME/.bashrc" 2>/dev/null || true
    [ -d "$HOME/.local/bin" ] && export PATH="$HOME/.local/bin:$PATH"
    [ -d "$HOME/micromamba/bin" ] && export PATH="$HOME/micromamba/bin:$PATH"
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
export PATH="$AMPLIPHY_DIR/bin:$HOME/gs:$HOME/.local/bin:$PATH"

# fastme / rapidnj / toytree / biopython が未導入の場合は自動インストール
MISSING_TOOLS=()
command -v fastme &>/dev/null || MISSING_TOOLS+=("fastme")
command -v rapidnj &>/dev/null || MISSING_TOOLS+=("rapidnj")
python3 -c "import toytree" &>/dev/null || MISSING_TOOLS+=("toytree")
python3 -c "import Bio" &>/dev/null || MISSING_TOOLS+=("biopython")
if [ "${#MISSING_TOOLS[@]}" -gt 0 ]; then
    echo "[INFO] biotools 環境に不足ツール/パッケージ (${MISSING_TOOLS[*]}) をインストールします..."
    if command -v micromamba &>/dev/null; then
        micromamba install -n biotools "${MISSING_TOOLS[@]}" -c bioconda -c conda-forge -y 2>&1 || true
    elif command -v conda &>/dev/null; then
        conda install -n biotools "${MISSING_TOOLS[@]}" -c bioconda -c conda-forge -y 2>&1 || true
    fi
fi

# 必要ツール・スクリプトの検出
GS2_BIN=""
for cand in "/home/citrusml2004/gs/gs2" "$HOME/gs/gs2" "$(which gs2 2>/dev/null || true)"; do
    if [ -n "$cand" ] && [ -x "$cand" ]; then GS2_BIN="$cand"; break; fi
done

MAD_PY="${AMPLIPHY_DIR}/bin/mad.py"
TCS_PY="${AMPLIPHY_DIR}/bin/tcs.py"
MAXALIGN_PL="${AMPLIPHY_DIR}/bin/maxalign.pl"

for tool_check in fastme rapidnj iqtree mafft mmseqs gotree nextflow; do
    if ! command -v "$tool_check" &>/dev/null; then
        echo "[ERROR] 必要コマンドが見つかりません: $tool_check" >&2
        exit 1
    fi
done

if [ -z "$GS2_BIN" ]; then
    echo "[ERROR] gs2 バイナリが見つかりません。" >&2
    exit 1
fi
for f in "$MAD_PY" "$TCS_PY" "$MAXALIGN_PL"; do
    if [ ! -f "$f" ]; then
        echo "[ERROR] 必要スクリプトが見つかりません: $f" >&2
        exit 1
    fi
done

echo "GS2 Binary:       $GS2_BIN"
echo "MAD Script:       $MAD_PY"
echo "TCS Script:       $TCS_PY"
echo "MaxAlign Script:  $MAXALIGN_PL"

# ------------------------------------------------------------------------------
# 3. 生物種ごとの実行関数
# ------------------------------------------------------------------------------
run_taxon() {
    local taxon="$1"
    local short_name=""
    case "$taxon" in
        "eukaryota") short_name="euk" ;;
        "bacteria")  short_name="bac" ;;
        *)           short_name="${taxon:0:3}" ;;
    esac

    echo ""
    echo "##############################################################################"
    echo "# E-value ベンチマーク実行開始: $taxon ($short_name)"
    echo "##############################################################################"

    local input_dir="$AMPLIPHY_ROOT/input_fasta/$taxon"
    if [ ! -d "$input_dir" ]; then
        echo "[ERROR] $taxon の入力FASTAフォルダ ($input_dir) が見つかりません。" >&2
        return 1
    fi

    local output_dir="$RESULTS_BASE_DIR/$taxon"
    local work_dir="$output_dir/work"
    local launch_dir="$output_dir/launch"
    local nxf_log="$AMPLIPHY_ROOT/logs/evalue_bench/nextflow_${taxon}.log"
    mkdir -p "$output_dir" "$work_dir" "$launch_dir"

    # Lineage TSV の生成 (TCS 評価用)
    local lineage_tsv="$output_dir/lineage.tsv"
    if [ ! -s "$lineage_tsv" ]; then
        local taxa_tsv="$DATASET_DIR/${short_name}_taxa.tsv"
        if [ -s "$taxa_tsv" ]; then
            echo "  [Lineage] $taxa_tsv から lineage.tsv を生成中..."
            python3 -c "
from src.pwa_pipeline import generate_lineage_tsv
import glob
fastas = sorted(glob.glob('$input_dir/*.fa'))
generate_lineage_tsv('$taxa_tsv', fastas, '$lineage_tsv')
"
        elif [ -s "$DATASET_DIR/${short_name}_lineage.tsv" ]; then
            cp -f "$DATASET_DIR/${short_name}_lineage.tsv" "$lineage_tsv"
        fi
    fi

    echo "  INPUT_DIR:   $input_dir"
    echo "  OUTPUT_DIR:  $output_dir"
    echo "  LINEAGE_TSV: $lineage_tsv"
    echo "  LOG_FILE:    $nxf_log (tail -f で進捗確認可)"

    cd "$launch_dir"
    export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx4g}"

    nextflow run "$PIPELINE_NF" \
        -work-dir "$work_dir" \
        -c "$CONFIG_FILE" \
        --input_dir "$input_dir" \
        --output_dir "$output_dir" \
        --database UniRef50 \
        --database_dir "$DB_DIR" \
        --tmp_dir "/tmp/ampliphy_ebench_${USER:-citrusml2004}_${taxon}" \
        --max_depth 5 \
        --random_seed "$RANDOM_SEED" \
        --evalues "$EVALUES" \
        --variants "$VARIANTS" \
        --lineage_file "$lineage_tsv" \
        --gs2_bin "$GS2_BIN" \
        --mad_script "$MAD_PY" \
        --tcs_script "$TCS_PY" \
        --maxalign_script "$MAXALIGN_PL" \
        -resume > "$nxf_log" 2>&1

    echo "完了: $taxon"
    echo "  -> TCSサマリー: $output_dir/report/tcs_summary.tsv"
}

# ------------------------------------------------------------------------------
# 4. メイン実行 (対象生物種をクラスタ並列投入)
# ------------------------------------------------------------------------------
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
echo "すべての E-value ベンチマーク実行が完了しました: $(date)"
echo "出力先: $RESULTS_BASE_DIR/<taxon>/report/tcs_summary.tsv"
echo "=============================================================================="
