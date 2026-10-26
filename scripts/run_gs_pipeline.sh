#!/bin/bash
# ==============================================================================
# Graph Splitting (GS) 法による配列追加前後の系統樹推定・TCSスコア向上検証パイプライン
#
# LSFバッチジョブ投入例:
#   bsub < scripts/run_gs_pipeline.sh
# 対話型シェル実行例:
#   bash scripts/run_gs_pipeline.sh
# ==============================================================================

#BSUB -J gs_tcs_pipeline
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -W 24:00
#BSUB -o Ampliphy2/logs/gs/gs_pipeline_%J.out
#BSUB -e Ampliphy2/logs/gs/gs_pipeline_%J.err

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. パスとディレクトリの設定
# ------------------------------------------------------------------------------
WORK_DIR="${WORK_DIR:-$HOME/Ampliphy2}"

# 入力ディレクトリ
# 1) オリジナルHOG配列 (FASTA)
ORIG_INPUT_DIR="${ORIG_INPUT_DIR:-$WORK_DIR/input_fasta/eukaryota}"
# 2) ホモログ追加後配列 (Amplified FASTA / MSA)
AMP_INPUT_DIR="${AMP_INPUT_DIR:-$WORK_DIR/results/eukaryota_new/msa}"

# 出力ディレクトリ
RESULTS_DIR="${RESULTS_DIR:-$WORK_DIR/results_gs/eukaryota}"
LOG_DIR="${RESULTS_DIR}/logs"
TREE_ORI_DIR="${RESULTS_DIR}/tree_ori"
TREE_AMP_DIR="${RESULTS_DIR}/tree_amp"
EVAL_DIR="${RESULTS_DIR}/eval"

mkdir -p "$LOG_DIR" "$TREE_ORI_DIR" "$TREE_AMP_DIR" "${EVAL_DIR}/aug" "${EVAL_DIR}/enr"

# 外部ツール・スクリプトのパス
GS2_BIN="${GS2_BIN:-$WORK_DIR/gs/gs2}"
MAD_SCRIPT="${MAD_SCRIPT:-$WORK_DIR/tools/mad2.2/mad/mad.py}"
TCS_SCRIPT="${TCS_SCRIPT:-$WORK_DIR/ampliphy-analysis/scripts/tcs_batch.py}"
LINEAGE_TSV="${LINEAGE_TSV:-$WORK_DIR/ampliphy-analysis/dataset/euk_lineage.tsv}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_MODULE_DIR="${PY_MODULE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"

THREADS=8
MMSEQS_SENS=7.5

echo "=============================================================================="
echo "GS Method (gs2) TCS Improvement Evaluation Pipeline"
echo "Started at: $(date)"
echo "Host: $(hostname)"
echo "WORK_DIR: $WORK_DIR"
echo "ORIG_INPUT_DIR: $ORIG_INPUT_DIR"
echo "AMP_INPUT_DIR:  $AMP_INPUT_DIR"
echo "RESULTS_DIR:    $RESULTS_DIR"
echo "GS2 Binary:     $GS2_BIN"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 2. 環境有効化の準備関数
# ------------------------------------------------------------------------------
activate_env() {
    local env_name="$1"
    if command -v micromamba &>/dev/null; then
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate "$env_name"
    elif command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "$env_name"
    else
        echo "Error: Neither micromamba nor conda found." >&2
        exit 1
    fi
}

# ------------------------------------------------------------------------------
# 3. 入力FASTAファイルの走査
# ------------------------------------------------------------------------------
shopt -s nullglob
orig_files=("$ORIG_INPUT_DIR"/*.fa)
total_hogs=${#orig_files[@]}

if [ "$total_hogs" -eq 0 ]; then
    echo "エラー: オリジナル配列ファイルが見つかりません: $ORIG_INPUT_DIR" >&2
    exit 1
fi

echo "検出されたHOG総数: $total_hogs 件"

# biotools 環境（gs2, mmseqs2, gotree）を有効化
activate_env "biotools"

# PATHに gs2 のディレクトリを追加
export PATH="$(dirname "$GS2_BIN"):$PATH"
export PYTHONPATH="${PY_MODULE_DIR}:${PYTHONPATH:-}"

success_count=0
count=0

# ------------------------------------------------------------------------------
# 4. 各HOGに対するGS法実行、MADルート化、GoTree枝刈り
# ------------------------------------------------------------------------------
for orig_fa in "${orig_files[@]}"; do
    ((count++))
    fname=$(basename "$orig_fa")
    hog_id="${fname%.fa}"

    echo ""
    echo "------------------------------------------------------------------------------"
    echo "[$count/$total_hogs] 処理中: $hog_id"
    echo "------------------------------------------------------------------------------"

    # ホモログ追加後FASTAの探索
    amp_fa=""
    for candidate in \
        "$AMP_INPUT_DIR/${hog_id}.amp.fa" \
        "$AMP_INPUT_DIR/${hog_id}_amp.fa" \
        "$AMP_INPUT_DIR/HOG_${hog_id}.amp.fa" \
        "$AMP_INPUT_DIR/${hog_id}.fa"; do
        if [ -f "$candidate" ]; then
            amp_fa="$candidate"
            break
        fi
    done

    if [ -z "$amp_fa" ]; then
        echo "  [WARN] 対応するホモログ追加後FASTAが見つかりません。スキップします: $hog_id"
        continue
    fi

    # 中間ファイル定義
    clean_orig_fa="${TREE_ORI_DIR}/${hog_id}.unaligned.fa"
    clean_amp_fa="${TREE_AMP_DIR}/${hog_id}.unaligned.fa"
    ori_nwk="${TREE_ORI_DIR}/${hog_id}.ori.nwk"
    amp_nwk="${TREE_AMP_DIR}/${hog_id}.amp.nwk"
    prune_list="${TREE_AMP_DIR}/${hog_id}.homolog_ids.txt"
    amp_pruned_nwk="${TREE_AMP_DIR}/${hog_id}.amp.pruned.nwk"

    # Step 1: ギャップ除去（Unalign）とヘッダー整形
    echo "  [Step 1] FASTAのギャップ除去とヘッダー整形..."
    python3 -c "
from src.gs_pipeline import prepare_unaligned_fasta
prepare_unaligned_fasta('$orig_fa', '$clean_orig_fa')
prepare_unaligned_fasta('$amp_fa', '$clean_amp_fa')
"

    # Step 2: GS法 (gs2) による系統樹推定
    echo "  [Step 2] GS法 (gs2) による系統樹推定..."
    # 1) オリジナル樹の推定
    if [ ! -s "$ori_nwk" ]; then
        echo "    - Original樹推定中..."
        "$GS2_BIN" -l -t "$THREADS" -m "$MMSEQS_SENS" -s "$clean_orig_fa" > "$ori_nwk"
    else
        echo "    - Original樹は既に存在します: $ori_nwk"
    fi

    # 2) ホモログ追加後樹の推定
    if [ ! -s "$amp_nwk" ]; then
        echo "    - Amplified樹推定中..."
        "$GS2_BIN" -l -t "$THREADS" -m "$MMSEQS_SENS" -s "$clean_amp_fa" > "$amp_nwk"
    else
        echo "    - Amplified樹は既に存在します: $amp_nwk"
    fi

    # Step 3: 系統樹のクリーンアップ & MADルート決定
    echo "  [Step 3] MADによるルート決定 (Rooting)..."
    python3 -c "
from src.gs_pipeline import clean_newick_tree
clean_newick_tree('$ori_nwk')
clean_newick_tree('$amp_nwk')
"

    # MAD実行
    ori_rooted="${ori_nwk}.rooted"
    if [ ! -s "$ori_rooted" ]; then
        python "$MAD_SCRIPT" "$ori_nwk" &> /dev/null || true
    fi

    amp_rooted="${amp_nwk}.rooted"
    if [ ! -s "$amp_rooted" ]; then
        python "$MAD_SCRIPT" "$amp_nwk" &> /dev/null || true
    fi

    if [ ! -s "$ori_rooted" ] || [ ! -s "$amp_rooted" ]; then
        echo "  [ERROR] MADルート化に失敗しました: $hog_id"
        continue
    fi

    # Step 4: 追加ホモログの枝刈り (GoTree Pruning)
    echo "  [Step 4] GoTreeによるホモログ枝刈り (Pruning)..."
    python3 -c "
from src.gs_pipeline import extract_homolog_ids
homologs = extract_homolog_ids('$clean_orig_fa', '$clean_amp_fa')
with open('$prune_list', 'w', encoding='utf-8') as f:
    for hid in homologs:
        f.write(f'{hid}\n')
"

    if [ -s "$prune_list" ]; then
        # 追加ホモログが存在する場合、枝刈りを実行
        tr -d "'" < "$amp_rooted" | gotree prune \
            -f "$prune_list" \
            -o "$amp_pruned_nwk" 2> /dev/null || true
    else
        # ホモログがない場合はそのまま使用
        cp "$amp_rooted" "$amp_pruned_nwk"
    fi

    if [ ! -s "$amp_pruned_nwk" ]; then
        echo "  [ERROR] GoTree枝刈りに失敗しました: $hog_id"
        continue
    fi

    # Step 5: 評価用ディレクトリへのシンボリックリンク配置
    ln -sf "$ori_rooted" "${EVAL_DIR}/aug/${hog_id}_gs_aug.nwk"
    ln -sf "$amp_pruned_nwk" "${EVAL_DIR}/enr/${hog_id}_gs_enr.nwk"

    ((success_count++))
    echo "  [SUCCESS] $hog_id の系統樹構築・ルート化・枝刈りが完了しました。"
done

echo ""
echo "=============================================================================="
echo "系統樹構築完了: $success_count / $total_hogs HOGs"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 5. TCS スコアの計算 (foldtree 環境)
# ------------------------------------------------------------------------------
OUT_SCORES="${RESULTS_DIR}/tcs_scores_gs.tsv"
echo "=== TCSスコアの計算を開始します (foldtree 環境) ==="
activate_env "foldtree"

cat /dev/null > "$OUT_SCORES"

if [ -s "$LINEAGE_TSV" ]; then
    echo "  [TCS] Enriched (enr) 木のスコア計算中..."
    python "$TCS_SCRIPT" "${EVAL_DIR}/enr" "$LINEAGE_TSV" >> "$OUT_SCORES"

    echo "  [TCS] Original (aug) 木のスコア計算中..."
    python "$TCS_SCRIPT" "${EVAL_DIR}/aug" "$LINEAGE_TSV" >> "$OUT_SCORES"
else
    echo "  [WARN] Lineageファイルが見つかりません: $LINEAGE_TSV"
    echo "  TCSスコアの算出をスキップします。必要に応じて後から手動実行してください。"
fi

# ------------------------------------------------------------------------------
# 6. 結果集計とサマリーレポート出力
# ------------------------------------------------------------------------------
REPORT_TXT="${RESULTS_DIR}/gs_tcs_summary.txt"
if [ -s "$OUT_SCORES" ]; then
    echo ""
    echo "=== TCSスコアの集計と改善度 (Delta TCS) の算出 ==="
    python3 -c "
import pandas as pd
from src.gs_pipeline import parse_gs_tcs_scores, summarize_gs_results

df = parse_gs_tcs_scores('$OUT_SCORES')
summary = summarize_gs_results(df)

report = f'''
============================================================
GS法 (gs2) TCS改善度 評価結果サマリー
============================================================
評価HOG数            : {summary['n_hogs']}
Original TCS 平均    : {summary['mean_tcs_ori']:.4f}
Enriched TCS 平均    : {summary['mean_tcs_enr']:.4f}
Delta TCS 比率 (平均): {summary['mean_delta_ratio']:.4f}
Delta TCS 比率 (中央): {summary['median_delta_ratio']:.4f}
Delta TCS 差分 (平均): {summary['mean_delta_diff']:.4f}
スコア改善HOG割合    : {summary['pct_improved']:.2f} %
Wilcoxon 検定 p-value: {summary['wilcoxon_pvalue']:.4e}
============================================================
'''
print(report)
with open('$REPORT_TXT', 'w', encoding='utf-8') as f:
    f.write(report)
df.to_csv('${RESULTS_DIR}/gs_tcs_comparison.tsv', sep='\t', index=False)
"
    echo "サマリーレポート保存先: $REPORT_TXT"
fi

echo "=============================================================================="
echo "GS Pipeline 完了: $(date)"
echo "結果ディレクトリ: $RESULTS_DIR"
echo "=============================================================================="
