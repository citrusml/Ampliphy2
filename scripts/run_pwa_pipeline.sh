#!/bin/bash
# ==============================================================================
# PWA (Pairwise Sequence Alignment) 系統樹作成＆配列増幅精度評価 パイプライン
#
# 対象生物種: 4種 (eukaryota, amniota, bacteria, actinomycete)
# 評価手法:
#   1) global : Needleman-Wunsch 全長アラインメント + 補正ポアソン距離 (ver2) + panjep NJ
#   2) local  : MMseqs2 局所アラインメント + 補正ポアソン距離 (ver3) + panjep NJ
#   3) panjep : panjep -d poisson 直接推定 (C++ 高速NJ)
#
# 使い方:
#   1) 全4種 × 全3手法 を一括実行する場合 (デフォルト):
#      bsub < scripts/run_pwa_pipeline.sh
#      または bash scripts/run_pwa_pipeline.sh
#
#   2) 特定の生物種のみ実行する場合:
#      bash scripts/run_pwa_pipeline.sh eukaryota
#      bash scripts/run_pwa_pipeline.sh eukaryota amniota
#
#   3) 特定の手法のみ実行する場合:
#      METHODS="panjep" bash scripts/run_pwa_pipeline.sh eukaryota
#      METHODS="global local" bsub < scripts/run_pwa_pipeline.sh
# ==============================================================================

#BSUB -J pwa_pipeline
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -W 48:00
#BSUB -cwd /home/citrusml2004/Ampliphy2
#BSUB -o logs/pwa/master_%J.out
#BSUB -e logs/pwa/master_%J.err

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
ANALYSIS_DIR="$AMPLIPHY_ROOT/ampliphy-analysis"
if [ ! -d "$ANALYSIS_DIR" ]; then
    if [ -d "/home/citrusml2004/ampliphy-analysis" ]; then
        ln -sfn "/home/citrusml2004/ampliphy-analysis" "$ANALYSIS_DIR"
    elif [ -d "$HOME/ampliphy-analysis" ]; then
        ln -sfn "$HOME/ampliphy-analysis" "$ANALYSIS_DIR"
    fi
fi
DATASET_DIR="$ANALYSIS_DIR/dataset"

# pairwise_alignment ディレクトリの配置
PAIRWISE_DIR="$AMPLIPHY_ROOT/pairwise_alignment"
if [ ! -d "$PAIRWISE_DIR" ]; then
    if [ -d "/home/citrusml2004/pairwise_alignment" ]; then
        ln -sfn "/home/citrusml2004/pairwise_alignment" "$PAIRWISE_DIR"
    elif [ -d "$HOME/pairwise_alignment" ]; then
        ln -sfn "$HOME/pairwise_alignment" "$PAIRWISE_DIR"
    fi
fi

# 必要なディレクトリ群を作成
mkdir -p "$AMPLIPHY_ROOT/logs/pwa" "$AMPLIPHY_ROOT/results_pwa" "$AMPLIPHY_ROOT/data/amplified_fasta"

# ------------------------------------------------------------------------------
# リアルタイムログ出力設定
# ------------------------------------------------------------------------------
REALTIME_LOG="$AMPLIPHY_ROOT/logs/pwa/pwa_pipeline_${LSB_JOBID:-master}.log"
exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

# ------------------------------------------------------------------------------
# 2. 環境有効化 (biotools: python, mmseqs, gotree, panjep)
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

    # 必要なPythonパッケージ (biopython, pandas, scipy, numpy) の自動確認・補完
    python3 -c "import Bio, pandas, scipy, numpy" 2>/dev/null || {
        echo "[INFO] 必要なPythonパッケージ (biopython, pandas, scipy, numpy) をインストール中..."
        pip install --quiet biopython pandas scipy numpy 2>/dev/null || true
    }
}

activate_env "biotools"
export PYTHONUNBUFFERED=1
export PYTHONPATH="$AMPLIPHY_ROOT:${PYTHONPATH:-}"

# ------------------------------------------------------------------------------
# 3. 必要ツールの検出と準備
# ------------------------------------------------------------------------------
# panjep バイナリの探索
PANJEP_BIN=""
for cand in \
    "$PAIRWISE_DIR/panjep/panjep" \
    "/home/citrusml2004/pairwise_alignment/panjep/panjep" \
    "$HOME/pairwise_alignment/panjep/panjep" \
    "$HOME/.local/bin/panjep" \
    "$(which panjep 2>/dev/null || true)"; do
    if [ -n "$cand" ] && [ -x "$cand" ]; then
        PANJEP_BIN="$cand"
        break
    fi
done

# panjep がビルドされていない場合は自動ビルド
if [ -z "$PANJEP_BIN" ] && [ -d "$PAIRWISE_DIR/panjep" ]; then
    echo "panjep バイナリをビルド中 ($PAIRWISE_DIR/panjep)..."
    (cd "$PAIRWISE_DIR/panjep" && make -j4) 2>&1 || true
    if [ -x "$PAIRWISE_DIR/panjep/panjep" ]; then
        PANJEP_BIN="$PAIRWISE_DIR/panjep/panjep"
    fi
fi

if [ -z "$PANJEP_BIN" ] || [ ! -x "$PANJEP_BIN" ]; then
    echo "[ERROR] panjep バイナリが見つかりません。ビルドまたはインストールを確認してください。" >&2
    exit 1
fi
echo "panjep バイナリ:    $PANJEP_BIN"

# MAD (Minimal Ancestor Deviation) スクリプトの探索
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

if [ -z "$MAD_PY" ]; then
    echo "[ERROR] mad.py スクリプトが見つかりません。" >&2
    exit 1
fi
echo "MAD スクリプト:       $MAD_PY"

# TCS スクリプトの探索
TCS_PY=""
for cand in \
    "$AMPLIPHY_DIR/bin/tcs.py" \
    "$ANALYSIS_DIR/scripts/tcs.py" \
    "$ANALYSIS_DIR/scripts/tcs_batch.py"; do
    if [ -f "$cand" ]; then
        TCS_PY="$cand"
        break
    fi
done

if [ -z "$TCS_PY" ]; then
    echo "[ERROR] tcs.py スクリプトが見つかりません。" >&2
    exit 1
fi
echo "TCS スクリプト:       $TCS_PY"

# 実行対象生物種 (引数がなければ4種すべて)
if [ "$#" -gt 0 ]; then
    TARGET_TAXA=("$@")
else
    TARGET_TAXA=("eukaryota" "amniota" "bacteria" "actinomycete")
fi

# 実行対象手法 (環境変数 METHODS で指定可能、未指定時は全3手法)
TARGET_METHODS=(${METHODS:-"global" "local" "panjep"})

echo "=============================================================================="
echo "AmpliPhy PWA (Pairwise Alignment) Pipeline Execution"
echo "Started at:       $(date)"
echo "Host:             $(hostname)"
echo "AMPLIPHY_ROOT:    $AMPLIPHY_ROOT"
echo "Target Taxa:      ${TARGET_TAXA[*]}"
echo "Target Methods:   ${TARGET_METHODS[*]}"
echo "panjep Bin:       $PANJEP_BIN"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 4. 単一生物種・手法の実行関数
# ------------------------------------------------------------------------------
process_taxon_method() {
    local taxon="$1"
    local method="$2"
    local short_name=""
    case "$taxon" in
        "eukaryota")   short_name="euk" ;;
        "amniota")     short_name="amn" ;;
        "bacteria")    short_name="bac" ;;
        "actinomycete") short_name="act" ;;
        *)             short_name="${taxon:0:3}" ;;
    esac

    echo ""
    echo "=============================================================================="
    echo "# 実行開始: $taxon ($short_name) | 手法: $method"
    echo "=============================================================================="

    local orig_dir="$AMPLIPHY_ROOT/input_fasta/$taxon"
    if [ ! -d "$orig_dir" ]; then
        echo "[ERROR] $taxon のオリジナルFASTAフォルダ ($orig_dir) が見つかりません。" >&2
        return 1
    fi

    # ホモログ取得元 (results/<taxon>/homologs 優先、なければ results_max/<taxon>/homologs)
    local hom_dir="$AMPLIPHY_ROOT/results/$taxon/homologs"
    if [ ! -d "$hom_dir" ] && [ -d "$AMPLIPHY_ROOT/results_max/$taxon/homologs" ]; then
        hom_dir="$AMPLIPHY_ROOT/results_max/$taxon/homologs"
    fi

    # 増幅配列出力ディレクトリ (data/amplified_fasta/<taxon>)
    local amp_dir="$AMPLIPHY_ROOT/data/amplified_fasta/$taxon"
    mkdir -p "$amp_dir"

    # 手法別出力ディレクトリ (results_pwa/<taxon>/<method>)
    local out_dir="$AMPLIPHY_ROOT/results_pwa/$taxon/$method"
    local matrix_dir="$out_dir/matrix"
    local tree_ori_dir="$out_dir/tree_ori"
    local tree_amp_dir="$out_dir/tree_amp"
    local tree_pruned_dir="$out_dir/tree_amp_pruned"
    local report_dir="$out_dir/report"
    mkdir -p "$matrix_dir" "$tree_ori_dir" "$tree_amp_dir" "$tree_pruned_dir" "$report_dir"

    # Lineage TSV (TCS評価用)
    local lineage_tsv="$AMPLIPHY_ROOT/results_pwa/$taxon/lineage.tsv"
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

    # TCS 出力レポートTSV
    local tcs_report="$report_dir/tcs_report.tsv"
    printf "family_id\ttcs_ori\ttcs_amp\tdelta_tcs\tratio_tcs\n" > "$tcs_report"

    shopt -s nullglob
    local orig_files=("$orig_dir"/*.fa)
    local total_hogs=${#orig_files[@]}

    if [ "$total_hogs" -eq 0 ]; then
        echo "[ERROR] $taxon のオリジナルFASTAファイルが見つかりません: $orig_dir" >&2
        return 1
    fi

    echo "  総HOG数: $total_hogs 件"
    echo "  出力先:  $out_dir"

    local count=0
    local success_count=0

    for ori_fa in "${orig_files[@]}"; do
        count=$((count + 1))
        local filename
        filename=$(basename "$ori_fa")
        local hog_id="${filename%.*}"

        printf "  [%3d/%3d] %s: 系統樹推定 & TCS評価..." "$count" "$total_hogs" "$hog_id"

        # 1. ホモログ配列とマージして amp.fa を作成
        local hom_fa="$hom_dir/${hog_id}.homologs.fa"
        local amp_fa="$amp_dir/${hog_id}.amp.fa"
        python3 -m src.pwa_pipeline merge "$ori_fa" "$hom_fa" "$amp_fa" >/dev/null 2>&1 || true

        # 2. 系統樹推定 (ori & amp)
        local ori_nwk="$tree_ori_dir/${hog_id}.ori.nwk"
        local amp_nwk="$tree_amp_dir/${hog_id}.amp.nwk"

        if [ "$method" = "panjep" ]; then
            # panjep -d poisson による直接推定
            "$PANJEP_BIN" -t 4 -d poisson -e 0 "$ori_fa" > "$ori_nwk" 2>/dev/null || true
            "$PANJEP_BIN" -t 4 -d poisson -e 0 "$amp_fa" > "$amp_nwk" 2>/dev/null || true

        elif [ "$method" = "global" ]; then
            # Needleman-Wunsch グローバルアラインメント + 補正ポアソン距離 (ver2)
            local ori_phy="$matrix_dir/${hog_id}.ori.phy"
            local amp_phy="$matrix_dir/${hog_id}.amp.phy"

            python3 -m src.pwa_pipeline matrix "$ori_fa" "$ori_phy" --method ver2 -t 4 >/dev/null 2>&1 || true
            python3 -m src.pwa_pipeline matrix "$amp_fa" "$amp_phy" --method ver2 -t 4 >/dev/null 2>&1 || true

            if [ -s "$ori_phy" ]; then
                "$PANJEP_BIN" -t 4 "$ori_phy" > "$ori_nwk" 2>/dev/null || true
            fi
            if [ -s "$amp_phy" ]; then
                "$PANJEP_BIN" -t 4 "$amp_phy" > "$amp_nwk" 2>/dev/null || true
            fi

        elif [ "$method" = "local" ]; then
            # MMseqs2 局所アラインメント + 補正ポアソン距離 (ver3)
            local ori_phy="$matrix_dir/${hog_id}.ori.phy"
            local amp_phy="$matrix_dir/${hog_id}.amp.phy"

            python3 -m src.pwa_pipeline matrix "$ori_fa" "$ori_phy" --method ver3 -t 4 >/dev/null 2>&1 || true
            python3 -m src.pwa_pipeline matrix "$amp_fa" "$amp_phy" --method ver3 -t 4 >/dev/null 2>&1 || true

            if [ -s "$ori_phy" ]; then
                "$PANJEP_BIN" -t 4 "$ori_phy" > "$ori_nwk" 2>/dev/null || true
            fi
            if [ -s "$amp_phy" ]; then
                "$PANJEP_BIN" -t 4 "$amp_phy" > "$amp_nwk" 2>/dev/null || true
            fi
        fi

        if [ ! -s "$ori_nwk" ] || [ ! -s "$amp_nwk" ]; then
            echo " [スキップ: 系統樹出力が空]"
            continue
        fi

        # 3. 負の枝長のクリーニング (NJ法で稀に生じる負の枝長を :0.0 に置換)
        python3 -m src.pwa_pipeline clean-tree "$ori_nwk" >/dev/null 2>&1 || true
        python3 -m src.pwa_pipeline clean-tree "$amp_nwk" >/dev/null 2>&1 || true

        # 4. MAD による根付け (Rooting)
        local ori_rooted="${ori_nwk}.rooted"
        local amp_rooted="${amp_nwk}.rooted"

        python3 "$MAD_PY" "$ori_nwk" -n >/dev/null 2>&1 || true
        python3 "$MAD_PY" "$amp_nwk" -n >/dev/null 2>&1 || true

        if [ ! -s "$ori_rooted" ] || [ ! -s "$amp_rooted" ]; then
            echo " [スキップ: MAD根付けに失敗]"
            continue
        fi

        # 同点根で複数行出力された場合は最良の1行目のみを採用
        head -n 1 "$ori_rooted" > "${ori_rooted}.tmp" && mv "${ori_rooted}.tmp" "$ori_rooted"
        head -n 1 "$amp_rooted" > "${amp_rooted}.tmp" && mv "${amp_rooted}.tmp" "$amp_rooted"

        # 5. GoTree によるホモログ枝刈り (Pruning)
        local amp_pruned="$tree_pruned_dir/${hog_id}.amp.pruned.nwk"
        if [ -s "$hom_fa" ]; then
            local tips_file="/tmp/${hog_id}_hom_tips_${USER:-citrusml2004}.txt"
            awk '/^>/{print substr($1,2)}' "$hom_fa" > "$tips_file"

            if command -v gotree &>/dev/null; then
                gotree prune -i "$amp_rooted" -f "$tips_file" -o "$amp_pruned" >/dev/null 2>&1 || true
            fi
            rm -f "$tips_file"
        fi

        # ホモログがない場合、またはgotree未導入時はrootedツリーをそのままコピー
        if [ ! -s "$amp_pruned" ]; then
            cp -f "$amp_rooted" "$amp_pruned"
        fi

        # 枝刈り後ツリーの Newick 整形
        python3 -m src.pwa_pipeline clean-tree "$amp_pruned" >/dev/null 2>&1 || true

        # 6. TCS スコアの計算
        local tcs_ori="0.0"
        local tcs_amp="0.0"
        if [ -s "$lineage_tsv" ]; then
            tcs_ori=$(python3 "$TCS_PY" "$ori_rooted" "$lineage_tsv" 2>/dev/null | cut -f2) || tcs_ori="0.0"
            tcs_amp=$(python3 "$TCS_PY" "$amp_pruned" "$lineage_tsv" 2>/dev/null | cut -f2) || tcs_amp="0.0"
        fi

        # 数値妥当性チェックと delta / ratio の計算 (awk で安全・高速に計算)
        local delta="0.000000"
        local ratio="nan"
        read -r delta ratio <<< $(awk -v o="${tcs_ori:-0.0}" -v a="${tcs_amp:-0.0}" 'BEGIN {
            vo = (o ~ /^[0-9]+(\.[0-9]+)?$/) ? o + 0.0 : 0.0;
            va = (a ~ /^[0-9]+(\.[0-9]+)?$/) ? a + 0.0 : 0.0;
            d = va - vo;
            r = (vo > 0) ? sprintf("%.6f", va / vo) : "nan";
            printf "%.6f %s\n", d, r;
        }')

        printf "%s\t%s\t%s\t%s\t%s\n" "$hog_id" "${tcs_ori:-0.0}" "${tcs_amp:-0.0}" "$delta" "$ratio" >> "$tcs_report"
        echo " 完了 (TCS ori: ${tcs_ori:-0.0} -> amp: ${tcs_amp:-0.0} | ratio: $ratio)"
        success_count=$((success_count + 1))
    done

    # 7. 統計サマリーの集計と表示
    local summary_file="$report_dir/tcs_summary.txt"
    python3 -c "
from src.pwa_pipeline import parse_pwa_tcs_scores, summarize_pwa_results
import pandas as pd

df = parse_pwa_tcs_scores('$tcs_report')
stats = summarize_pwa_results(df)

with open('$summary_file', 'w', encoding='utf-8') as f:
    f.write('=======================================================\\n')
    f.write(f'PWA TCS Summary: $taxon | Method: $method\\n')
    f.write('=======================================================\\n')
    f.write(f'Total Analyzed HOGs : {stats[\"total_hogs\"]}\\n')
    f.write(f'Improved HOGs (Amp > Ori) : {stats[\"improved_count\"]}\\n')
    f.write(f'Declined HOGs (Amp < Ori) : {stats[\"declined_count\"]}\\n')
    f.write(f'Unchanged HOGs            : {stats[\"unchanged_count\"]}\\n')
    f.write(f'Mean Delta TCS (Amp - Ori): {stats[\"mean_delta_tcs\"]:+.6f} ± {stats[\"std_delta_tcs\"]:.6f}\\n')
    f.write(f'Mean TCS Ratio (Amp / Ori): {stats[\"mean_ratio_tcs\"]:.6f} ± {stats[\"std_ratio_tcs\"]:.6f}\\n')
    f.write(f'Paired t-test p-value     : {stats[\"p_value\"]:.4e}\\n')
    f.write('=======================================================\\n')
" 2>/dev/null || true

    if [ -f "$summary_file" ]; then
        cat "$summary_file"
    fi
}

# ------------------------------------------------------------------------------
# 5. メインループ実行 (生物種 × 手法)
# ------------------------------------------------------------------------------
for taxon in "${TARGET_TAXA[@]}"; do
    for method in "${TARGET_METHODS[@]}"; do
        process_taxon_method "$taxon" "$method"
    done
done

# ------------------------------------------------------------------------------
# 6. 全体比較サマリーテーブルの出力
# ------------------------------------------------------------------------------
echo ""
echo "=============================================================================="
echo "すべての PWA パイプライン実行が完了しました: $(date)"
echo "成果物ディレクトリ: results_pwa/<taxon>/<method>/"
echo "=============================================================================="

COMPARISON_SUMMARY="$AMPLIPHY_ROOT/results_pwa/overall_comparison_summary.tsv"
printf "Taxon\tMethod\tTotal_HOGs\tImproved\tDeclined\tMean_Delta_TCS\tMean_TCS_Ratio\tp_value\n" > "$COMPARISON_SUMMARY"

python3 -c "
import glob, os
import pandas as pd
from src.pwa_pipeline import parse_pwa_tcs_scores, summarize_pwa_results

reports = sorted(glob.glob('$AMPLIPHY_ROOT/results_pwa/*/*/report/tcs_report.tsv'))
if reports:
    print('\n[全体比較サマリーテーブル]')
    print(f'{\"Taxon\":<15} {\"Method\":<10} {\"Total\":<6} {\"Imp\":<5} {\"Dec\":<5} {\"Mean ΔTCS\":<14} {\"Mean Ratio\":<14} {\"p-value\":<10}')
    print('-' * 85)
    for rep in reports:
        parts = rep.split(os.sep)
        taxon = parts[-4]
        method = parts[-3]
        try:
            df = parse_pwa_tcs_scores(rep)
            st = summarize_pwa_results(df)
            print(f'{taxon:<15} {method:<10} {st[\"total_hogs\"]:>5} {st[\"improved_count\"]:>5} {st[\"declined_count\"]:>5} {st[\"mean_delta_tcs\"]:>+13.4f} {st[\"mean_ratio_tcs\"]:>13.4f} {st[\"p_value\"]:>10.3e}')
            with open('$COMPARISON_SUMMARY', 'a', encoding='utf-8') as out:
                out.write(f'{taxon}\\t{method}\\t{st[\"total_hogs\"]}\\t{st[\"improved_count\"]}\\t{st[\"declined_count\"]}\\t{st[\"mean_delta_tcs\"]:.6f}\\t{st[\"mean_ratio_tcs\"]:.6f}\\t{st[\"p_value\"]:.4e}\\n')
        except Exception:
            continue
    print('-' * 85)
" 2>/dev/null || true

echo "サマリーTSV保存完了: $COMPARISON_SUMMARY"

