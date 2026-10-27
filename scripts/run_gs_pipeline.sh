#!/bin/bash
# ==============================================================================
# Graph Splitting (GS) 法による配列追加前後の系統樹推定・TCSスコア向上検証パイプライン
# 4分類群一括実行対応版 (Eukaryota, Amniota, Bacteria, Actinomycete)
#
# LSFバッチジョブ投入例 (全4種を一括実行):
#   bsub < scripts/run_gs_pipeline.sh
# 特定の生物種のみ実行する例:
#   bash scripts/run_gs_pipeline.sh bacteria
#   bash scripts/run_gs_pipeline.sh eukaryota amniota
# ==============================================================================

#BSUB -J "gs_pipeline[1-4]"
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -W 12:00
#BSUB -cwd /home/citrusml2004/Ampliphy2
#BSUB -o logs/gs/gs_%J_%I.out
#BSUB -e logs/gs/gs_%J_%I.err

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. 基本パスとツールの設定 (すべて /home/citrusml2004/Ampliphy2 配下に完全統一)
# ------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AMPLIPHY_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
if [ ! -d "$AMPLIPHY_ROOT/input_fasta" ]; then
    if [ -d "/home/citrusml2004/Ampliphy2" ]; then
        AMPLIPHY_ROOT="/home/citrusml2004/Ampliphy2"
    elif [ -d "$HOME/Ampliphy2" ]; then
        AMPLIPHY_ROOT="$HOME/Ampliphy2"
    fi
fi

cd "$AMPLIPHY_ROOT"
WORK_DIR="$AMPLIPHY_ROOT"
PY_MODULE_DIR="$AMPLIPHY_ROOT"

mkdir -p "$AMPLIPHY_ROOT/logs/gs" "$AMPLIPHY_ROOT/results_gs"

# 外部ツール・スクリプトのパス (/home/citrusml2004 直下の確定パス)
GS2_BIN="/home/citrusml2004/gs/gs2"
if [ ! -f "$GS2_BIN" ]; then
    GS2_BIN="$HOME/gs/gs2"
fi

AMPLIPHY_DIR="/home/citrusml2004/ampliphy"
if [ ! -d "$AMPLIPHY_DIR" ]; then
    AMPLIPHY_DIR="$HOME/ampliphy"
fi
# GS法では MAD は使用せず、PhyloMethod_Comparison 準拠で gs2 の第一分割（トップダウン2分割）を根としてそのまま使用
# MAD_SCRIPT="$AMPLIPHY_DIR/bin/mad.py"

DATASET_DIR="/home/citrusml2004/ampliphy-analysis/dataset"
if [ ! -d "$DATASET_DIR" ]; then
    DATASET_DIR="$HOME/ampliphy-analysis/dataset"
fi

TCS_SCRIPT="/home/citrusml2004/ampliphy-analysis/scripts/tcs_batch.py"
if [ ! -f "$TCS_SCRIPT" ]; then
    TCS_SCRIPT="$HOME/ampliphy-analysis/scripts/tcs_batch.py"
fi

# スレッド数 (LSF -n 8 割当にあわせてデフォルト 8 スレッド)
# 8コア要求のため、クラスタ内 (mafftsrv0〜5, mafftweb) の空きコアへ即座に充填されます
THREADS="${THREADS:-8}"
MMSEQS_SENS=7.5

# 全4生物種の定義
ALL_TAXA=("eukaryota" "amniota" "bacteria" "actinomycete")

# 実行対象生物種の決定 (LSF Job Array / CLI引数 / 全種一括に対応)
if [ "$#" -gt 0 ]; then
    # 手動指定時 (例: bash scripts/run_gs_pipeline.sh bacteria)
    TARGET_TAXA=("$@")
elif [ -n "${LSB_JOBINDEX:-}" ] && [ "$LSB_JOBINDEX" -ge 1 ] && [ "$LSB_JOBINDEX" -le "${#ALL_TAXA[@]}" ]; then
    # LSF ジョブ配列投入時 (bsub < scripts/run_gs_pipeline.sh)
    # 4台の空きサーバーへ1生物種ずつ自動分配
    TARGET_TAXA=("${ALL_TAXA[$((LSB_JOBINDEX - 1))]}")
else
    # 非LSF対話的実行時 (全種順番に実行)
    TARGET_TAXA=("${ALL_TAXA[@]}")
fi

# ------------------------------------------------------------------------------
# リアルタイムログ出力設定 (LSFバッファリング回避 + tail -f で即時進捗確認可能)
# ------------------------------------------------------------------------------
mkdir -p "$AMPLIPHY_ROOT/logs/gs"
if [ -n "${LSB_JOBID:-}" ] && [ -n "${LSB_JOBINDEX:-}" ]; then
    REALTIME_LOG="$AMPLIPHY_ROOT/logs/gs/gs_${LSB_JOBID}_${LSB_JOBINDEX}.log"
elif [ -n "${LSB_JOBID:-}" ]; then
    REALTIME_LOG="$AMPLIPHY_ROOT/logs/gs/gs_${LSB_JOBID}.log"
elif [ "${#TARGET_TAXA[@]}" -eq 1 ]; then
    REALTIME_LOG="$AMPLIPHY_ROOT/logs/gs/gs_${TARGET_TAXA[0]}.log"
else
    REALTIME_LOG="$AMPLIPHY_ROOT/logs/gs/gs_master.log"
fi

# Pythonの標準出力バッファリングを無効化 (リアルタイムストリーミング保証)
export PYTHONUNBUFFERED=1

exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

echo "=============================================================================="
echo "GS Method (gs2) Multi-Taxa TCS Evaluation Pipeline"
echo "Started at:    $(date)"
echo "Host:          $(hostname)"
echo "Allocated CPUs:$THREADS"
echo "WORK_DIR:      $WORK_DIR"
echo "PY_MODULE_DIR: $PY_MODULE_DIR"
echo "GS2 Binary:    $GS2_BIN"
echo "Target Taxa:   ${TARGET_TAXA[*]}"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 2. 環境有効化の準備関数 (LSF非対話型ジョブ対応)
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

# ------------------------------------------------------------------------------
# 3. 生物種ごとの処理関数
# ------------------------------------------------------------------------------
process_taxon() {
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
    echo "# 生物種グループ開始: $taxon ($short_name)"
    echo "# リアルタイム監視: tail -f logs/gs/gs_${taxon}.log"
    echo "##############################################################################"

    # 対象生物種のリアルタイムログシンボリックリンクを作成
    ln -sfn "$REALTIME_LOG" "$AMPLIPHY_ROOT/logs/gs/gs_${taxon}.log"

    # 入力FASTAディレクトリ (探索なし・Ampliphy2/input_fasta/<taxon> 固定)
    local orig_dir="$AMPLIPHY_ROOT/input_fasta/$taxon"
    if [ ! -d "$orig_dir" ]; then
        echo "[ERROR] $taxon のオリジナルFASTAディレクトリ ($orig_dir) が見つかりません。" >&2
        return 1
    fi

    # AmpliPhy 実行後の増幅後MSA (Ampliphy2/results/<taxon>/msa/*.amp.fa)
    local amp_dir="$AMPLIPHY_ROOT/results/$taxon/msa"
    if [ ! -d "$amp_dir" ] && [ -d "results/$taxon/msa" ]; then
        amp_dir="$(pwd)/results/$taxon/msa"
    elif [ ! -d "$amp_dir" ] && [ -d "$HOME/Ampliphy2/results/$taxon/msa" ]; then
        amp_dir="$HOME/Ampliphy2/results/$taxon/msa"
    fi

    shopt -s nullglob
    local amp_check_files=("$amp_dir"/*.amp.fa)
    if [ ! -d "$amp_dir" ] || [ ${#amp_check_files[@]} -eq 0 ]; then
        echo "[ERROR] $taxon のホモログ追加後FASTAディレクトリ ($amp_dir/*.amp.fa) にファイルが見つかりません。" >&2
        echo "       AmpliPhy がまだ実行されていないか、results/$taxon/msa/*.amp.fa が未生成です。" >&2
        echo "       先に 'bsub < scripts/run_ampliphy.sh' を実行してください。" >&2
        return 1
    fi

    # 出力ディレクトリ定義
    local results_dir="$AMPLIPHY_ROOT/results_gs/$taxon"

    local log_dir="${results_dir}/logs"
    local tree_ori_dir="${results_dir}/tree_ori"
    local tree_amp_dir="${results_dir}/tree_amp"
    local eval_dir="${results_dir}/eval"

    mkdir -p "$log_dir" "$tree_ori_dir" "$tree_amp_dir" "${eval_dir}/aug" "${eval_dir}/enr"

    echo "  ORIG_INPUT_DIR: $orig_dir"
    echo "  AMP_INPUT_DIR:  $amp_dir (${#amp_check_files[@]} files found)"
    echo "  RESULTS_DIR:    $results_dir"

    # Lineage TSV (データセット配下、または自動生成)
    local lineage_tsv=""
    if [ -s "$DATASET_DIR/${short_name}_lineage.tsv" ]; then
        lineage_tsv="$DATASET_DIR/${short_name}_lineage.tsv"
    fi

    if [ -z "$lineage_tsv" ]; then
        local taxa_meta="$DATASET_DIR/${short_name}_taxa.tsv"
        if [ -s "$taxa_meta" ]; then
            echo "  [Lineage] ${short_name}_taxa.tsv から lineage.tsv を自動生成します..."
            lineage_tsv="${results_dir}/lineage.tsv"
            python3 -c "
import glob
from src.gs_pipeline import generate_lineage_tsv
fastas = sorted(glob.glob('$orig_dir/*.fa'))
generate_lineage_tsv('$taxa_meta', fastas, '$lineage_tsv')
"
        fi
    fi

    # 入力ファイルの走査
    shopt -s nullglob
    local orig_files=("$orig_dir"/*.fa)
    local total_hogs=${#orig_files[@]}

    echo "  検出されたHOG総数: $total_hogs 件"

    # biotools 環境（gs2, mmseqs2, gotree）を有効化
    activate_env "biotools"
    export PATH="$(dirname "$GS2_BIN"):$PATH"
    export PYTHONPATH="${PY_MODULE_DIR}:${PYTHONPATH:-}"

    local success_count=0
    local count=0

    for orig_fa in "${orig_files[@]}"; do
        ((count++))
        local fname=$(basename "$orig_fa")
        local hog_id="${fname%.fa}"

        printf "  [%3d/%3d] %s: 系統樹推定 (gs2) & 枝刈り..." "$count" "$total_hogs" "$hog_id"

        # ホモログ追加後FASTA (AmpliPhy出力: results/<taxon>/msa/<hog_id>.amp.fa)
        local amp_fa="${amp_dir}/${hog_id}.amp.fa"
        if [ ! -f "$amp_fa" ]; then
            echo " [スキップ: $amp_fa 未検出]"
            continue
        fi

        local clean_orig_fa="${tree_ori_dir}/${hog_id}.unaligned.fa"
        local clean_amp_fa="${tree_amp_dir}/${hog_id}.unaligned.fa"
        local ori_nwk="${tree_ori_dir}/${hog_id}.ori.nwk"
        local amp_nwk="${tree_amp_dir}/${hog_id}.amp.nwk"
        local prune_list="${tree_amp_dir}/${hog_id}.homolog_ids.txt"
        local amp_pruned_nwk="${tree_amp_dir}/${hog_id}.amp.pruned.nwk"

        # Step 1: ギャップ除去（Unalign）
        python3 -c "
from src.gs_pipeline import prepare_unaligned_fasta
prepare_unaligned_fasta('$orig_fa', '$clean_orig_fa')
prepare_unaligned_fasta('$amp_fa', '$clean_amp_fa')
"

        # Step 2: GS法 (gs2) 実行
        if [ ! -s "$ori_nwk" ]; then
            "$GS2_BIN" -l -t "$THREADS" -m "$MMSEQS_SENS" -s "$clean_orig_fa" > "$ori_nwk" 2>/dev/null || true
        fi

        if [ ! -s "$amp_nwk" ]; then
            "$GS2_BIN" -l -t "$THREADS" -m "$MMSEQS_SENS" -s "$clean_amp_fa" > "$amp_nwk" 2>/dev/null || true
        fi

        if [ ! -s "$ori_nwk" ] || [ ! -s "$amp_nwk" ]; then
            echo " [エラー: gs2系統樹出力が空]"
            continue
        fi

        # Step 3: クリーンアップ & 第一分割ルートの採用（MAD不使用・PhyloMethod_Comparison準拠）
        python3 -c "
from src.gs_pipeline import clean_newick_tree
clean_newick_tree('$ori_nwk')
clean_newick_tree('$amp_nwk')
"

        # Step 4: 追加ホモログの枝刈り (GoTree Pruning)
        python3 -c "
from src.gs_pipeline import extract_homolog_ids
homologs = extract_homolog_ids('$clean_orig_fa', '$clean_amp_fa')
with open('$prune_list', 'w', encoding='utf-8') as f:
    for hid in homologs:
        f.write(f'{hid}\n')
"

        if [ -s "$prune_list" ]; then
            tr -d "'" < "$amp_nwk" | gotree prune \
                -f "$prune_list" \
                -o "$amp_pruned_nwk" 2> /dev/null || true
        else
            cp "$amp_nwk" "$amp_pruned_nwk"
        fi

        if [ ! -s "$amp_pruned_nwk" ]; then
            echo " [エラー: GoTree枝刈りに失敗]"
            continue
        fi

        # 枝刈り後Newickの構文クリーンアップ（セミコロン保証等）
        python3 -c "
from src.gs_pipeline import clean_newick_tree
clean_newick_tree('$amp_pruned_nwk')
"

        # 評価用ディレクトリへのシンボリックリンク（第一分割根付き樹）
        ln -sf "$ori_nwk" "${eval_dir}/aug/${hog_id}_gs_aug.nwk"
        ln -sf "$amp_pruned_nwk" "${eval_dir}/enr/${hog_id}_gs_enr.nwk"

        echo " 完了"
        ((success_count++))
    done

    echo "  系統樹構築完了: $success_count / $total_hogs HOGs ($taxon)"

    # Step 5: TCS スコア計算 (foldtree 環境)
    local out_scores="${results_dir}/tcs_scores_gs.tsv"
    activate_env "foldtree"
    cat /dev/null > "$out_scores"

    if [ -s "$lineage_tsv" ]; then
        echo "  [Step 5/6] TCS スコア計算中 (eval/aug & eval/enr)..."
        python "$TCS_SCRIPT" "${eval_dir}/enr" "$lineage_tsv" >> "$out_scores"
        python "$TCS_SCRIPT" "${eval_dir}/aug" "$lineage_tsv" >> "$out_scores"
        echo "  [Step 5/6] TCS スコア計算完了"
    else
        echo "  [WARN] Lineageファイルが見つかりません。TCS計算をスキップ。"
    fi

    # Step 6: 結果集計とサマリーレポート
    local report_txt="${results_dir}/gs_tcs_summary.txt"
    if [ -s "$out_scores" ]; then
        echo "  [Step 6/6] TCS 結果集計およびサマリー出力中..."
        python3 -c "
import pandas as pd
from src.gs_pipeline import parse_gs_tcs_scores, summarize_gs_results

df = parse_gs_tcs_scores('$out_scores')
summary = summarize_gs_results(df)

report = f'''
============================================================
GS法 (gs2) TCS改善度 サマリー: $taxon
============================================================
評価HOG数            : {summary['n_hogs']}
Original TCS 平均    : {summary['mean_tcs_ori']:.4f}
Enriched TCS 平均    : {summary['mean_tcs_enr']:.4f}
Delta TCS 比率 (平均): {summary['mean_delta_ratio']:.4f}
Delta TCS 差分 (平均): {summary['mean_delta_diff']:.4f}
スコア改善HOG割合    : {summary['pct_improved']:.2f} %
Wilcoxon 検定 p-value: {summary['wilcoxon_pvalue']:.4e}
============================================================
'''
print(report)
with open('$report_txt', 'w', encoding='utf-8') as f:
    f.write(report)
df.to_csv('${results_dir}/gs_tcs_comparison.tsv', sep='\t', index=False)
"
    fi
}

# ------------------------------------------------------------------------------
# 4. メインループ実行
# ------------------------------------------------------------------------------
for taxon in "${TARGET_TAXA[@]}"; do
    process_taxon "$taxon" || echo "[WARN] $taxon の処理中にエラーが発生しました。"
done

echo ""
echo "=============================================================================="
echo "全対象生物種の GS Pipeline 実行が完了しました: $(date)"
echo "=============================================================================="
