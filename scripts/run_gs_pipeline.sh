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
#BSUB -o logs/gs/gs_%J_%I.out
#BSUB -e logs/gs/gs_%J_%I.err

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. 基本パスとツールの設定
# ------------------------------------------------------------------------------
WORK_DIR="${WORK_DIR:-$HOME}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_MODULE_DIR="${PY_MODULE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"

mkdir -p "$PY_MODULE_DIR/logs/gs"

# 外部ツール・スクリプトのパス
GS2_BIN="${GS2_BIN:-$WORK_DIR/gs/gs2}"
MAD_SCRIPT="${MAD_SCRIPT:-$WORK_DIR/tools/mad2.2/mad/mad.py}"
TCS_SCRIPT="${TCS_SCRIPT:-$WORK_DIR/ampliphy-analysis/scripts/tcs_batch.py}"
DATASET_DIR="${DATASET_DIR:-$PY_MODULE_DIR/ampliphy-analysis/dataset}"

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
    echo "##############################################################################"

    # 入力ディレクトリの自動探索 (Ampliphy2/input_fasta/<taxon> を最優先)
    local orig_dir=""
    for cand in \
        "$PY_MODULE_DIR/input_fasta/$taxon" \
        "$WORK_DIR/Ampliphy2/input_fasta/$taxon" \
        "$WORK_DIR/input_fasta/$taxon" \
        "$WORK_DIR/archive/input_fasta/$taxon"; do
        if [ -d "$cand" ] && compgen -G "$cand/*.fa" > /dev/null; then
            orig_dir="$cand"
            break
        fi
    done

    if [ -z "$orig_dir" ]; then
        echo "[ERROR] $taxon のオリジナルFASTAディレクトリが見つかりません。" >&2
        return 1
    fi

    # AmpliPhy 実行後の増幅後MSA (results/<taxon>/msa/*.amp.fa) を最優先探索
    local amp_dir=""
    for cand in \
        "$PY_MODULE_DIR/results/$taxon/msa" \
        "$WORK_DIR/Ampliphy2/results/$taxon/msa" \
        "$WORK_DIR/results/$taxon/msa" \
        "$WORK_DIR/Ampliphy2/results/${taxon}_new/msa" \
        "$WORK_DIR/results/${taxon}_new/msa" \
        "$WORK_DIR/results/${taxon}_new/${short_name}3/msa" \
        "$WORK_DIR/Ampliphy2/results/${taxon}/${short_name}3/msa" \
        "$WORK_DIR/results/${taxon}/${short_name}3/msa" \
        "$WORK_DIR/results/${taxon}/1e-3/msa"; do
        if [ -d "$cand" ] && compgen -G "$cand/*.amp.fa" > /dev/null; then
            amp_dir="$cand"
            break
        elif [ -d "$cand" ] && compgen -G "$cand/*.fa" > /dev/null; then
            amp_dir="$cand"
            break
        fi
    done

    if [ -z "$amp_dir" ]; then
        echo "[ERROR] $taxon のホモログ追加後FASTAディレクトリ (AMP_INPUT_DIR) が見つかりません。" >&2
        echo "       AmpliPhy がまだ実行されていないか、results/$taxon/msa/*.amp.fa が未生成です。" >&2
        echo "       先に 'bsub < scripts/run_ampliphy.sh' を実行してください。" >&2
        return 1
    fi

    # 出力ディレクトリ定義
    local results_dir="$PY_MODULE_DIR/results_gs/$taxon"

    local log_dir="${results_dir}/logs"
    local tree_ori_dir="${results_dir}/tree_ori"
    local tree_amp_dir="${results_dir}/tree_amp"
    local eval_dir="${results_dir}/eval"

    mkdir -p "$log_dir" "$tree_ori_dir" "$tree_amp_dir" "${eval_dir}/aug" "${eval_dir}/enr"

    echo "  ORIG_INPUT_DIR: $orig_dir"
    echo "  AMP_INPUT_DIR:  $amp_dir"
    echo "  RESULTS_DIR:    $results_dir"

    # Lineage TSV の探索・自動生成
    local lineage_tsv=""
    for cand in \
        "$DATASET_DIR/${short_name}_lineage.tsv" \
        "$WORK_DIR/results/${taxon}/${short_name}3/lineage.tsv" \
        "$WORK_DIR/results/${taxon}_new/${short_name}3/lineage.tsv" \
        "$WORK_DIR/results/${taxon}/lineage.tsv"; do
        if [ -s "$cand" ]; then
            lineage_tsv="$cand"
            break
        fi
    done

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

        # ホモログ追加後FASTAの探索
        local amp_fa=""
        for candidate in \
            "$amp_dir/${hog_id}.amp.fa" \
            "$amp_dir/${hog_id}_amp.fa" \
            "$amp_dir/HOG_${hog_id}.amp.fa" \
            "$amp_dir/${hog_id}.fa"; do
            if [ -f "$candidate" ]; then
                amp_fa="$candidate"
                break
            fi
        done

        if [ -z "$amp_fa" ]; then
            echo "    [WARN] 対応するホモログ追加後FASTAが見つかりません。スキップ: $hog_id"
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
            echo "    [ERROR] gs2 系統樹出力が空です: $hog_id"
            continue
        fi

        # Step 3: クリーンアップ & MAD ルート化
        python3 -c "
from src.gs_pipeline import clean_newick_tree
clean_newick_tree('$ori_nwk')
clean_newick_tree('$amp_nwk')
"

        local ori_rooted="${ori_nwk}.rooted"
        if [ ! -s "$ori_rooted" ]; then
            python "$MAD_SCRIPT" "$ori_nwk" &> /dev/null || true
        fi

        local amp_rooted="${amp_nwk}.rooted"
        if [ ! -s "$amp_rooted" ]; then
            python "$MAD_SCRIPT" "$amp_nwk" &> /dev/null || true
        fi

        if [ ! -s "$ori_rooted" ] || [ ! -s "$amp_rooted" ]; then
            echo "    [ERROR] MADルート化に失敗しました: $hog_id"
            continue
        fi

        # Step 4: 追加ホモログの枝刈り (GoTree Pruning)
        python3 -c "
from src.gs_pipeline import extract_homolog_ids
homologs = extract_homolog_ids('$clean_orig_fa', '$clean_amp_fa')
with open('$prune_list', 'w', encoding='utf-8') as f:
    for hid in homologs:
        f.write(f'{hid}\n')
"

        if [ -s "$prune_list" ]; then
            tr -d "'" < "$amp_rooted" | gotree prune \
                -f "$prune_list" \
                -o "$amp_pruned_nwk" 2> /dev/null || true
        else
            cp "$amp_rooted" "$amp_pruned_nwk"
        fi

        if [ ! -s "$amp_pruned_nwk" ]; then
            echo "    [ERROR] GoTree枝刈りに失敗しました: $hog_id"
            continue
        fi

        # 評価用ディレクトリへのシンボリックリンク
        ln -sf "$ori_rooted" "${eval_dir}/aug/${hog_id}_gs_aug.nwk"
        ln -sf "$amp_pruned_nwk" "${eval_dir}/enr/${hog_id}_gs_enr.nwk"

        ((success_count++))
    done

    echo "  系統樹構築完了: $success_count / $total_hogs HOGs ($taxon)"

    # Step 5: TCS スコア計算 (foldtree 環境)
    local out_scores="${results_dir}/tcs_scores_gs.tsv"
    activate_env "foldtree"
    cat /dev/null > "$out_scores"

    if [ -s "$lineage_tsv" ]; then
        echo "  [TCS] スコア計算中..."
        python "$TCS_SCRIPT" "${eval_dir}/enr" "$lineage_tsv" >> "$out_scores"
        python "$TCS_SCRIPT" "${eval_dir}/aug" "$lineage_tsv" >> "$out_scores"
    else
        echo "  [WARN] Lineageファイルが見つかりません。TCS計算をスキップ。"
    fi

    # Step 6: 結果集計とサマリーレポート
    local report_txt="${results_dir}/gs_tcs_summary.txt"
    if [ -s "$out_scores" ]; then
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
