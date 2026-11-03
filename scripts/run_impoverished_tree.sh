#!/bin/bash
# ==============================================================================
# Impoverished tree 作成パイプライン (NIGSC / LSF)
#
# Impoverished tree = 増幅MSAから元配列(HOGメンバー)のみを抜き出し (全ギャップ列は除去)、
#                     ori / enrich 木と同条件の IQ-TREE で再推定した未ルート系統樹。
#
# 入力 (既存成果物):
#   results[_max]/<taxon>/msa/<hog>.{maxaligned,after.maxaligned,amp}.fa  (この優先順で自動選択)
#   input_fasta/<taxon>/<hog>.fa                                          (元配列ID)
#   results_pwa/<taxon>/<method>/matrix/<hog>.amp.phy                     (global / local)
# 出力:
#   results[_max]/<taxon>/msa_imp/<hog>.imp.fa
#   results[_max]/<taxon>/tree_imp/<hog>.imp.nwk
#   results_pwa/<taxon>/<method>/tree_imp/<hog>.imp.nwk   (+ matrix/<hog>.imp.phy)
#   ※ results_pwa の panjep 手法は距離行列が無いため対象外 (global/local のみ)
#   ※ 既存の出力ファイルは上書きしません (再実行時はスキップ)
#
# 使い方:
#   mkdir -p logs/imp && bsub < scripts/run_impoverished_tree.sh        # 全て
#   RESULT_ROOTS="results_max" bash scripts/run_impoverished_tree.sh eukaryota
#   RESULT_ROOTS="results_pwa" METHODS="global local" bash scripts/run_impoverished_tree.sh bacteria
#   環境変数: RESULT_ROOTS (既定 "results results_max results_pwa")
#             METHODS (PWA用, 既定 "global local")
#             IQTREE_THREADS (既定 4) / PARALLEL_JOBS (既定 4)
# ==============================================================================

#BSUB -J imp_tree
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -W 72:00
#BSUB -cwd /home/citrusml2004/Ampliphy2
#BSUB -o logs/imp/master_%J.out
#BSUB -e logs/imp/master_%J.err

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. パス設定
# ------------------------------------------------------------------------------
if [ -d "/home/citrusml2004/Ampliphy2" ]; then
    AMPLIPHY_ROOT="/home/citrusml2004/Ampliphy2"
else
    AMPLIPHY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
cd "$AMPLIPHY_ROOT"

PAIRWISE_DIR="$AMPLIPHY_ROOT/pairwise_alignment"
if [ ! -d "$PAIRWISE_DIR" ]; then
    if [ -d "/home/citrusml2004/pairwise_alignment" ]; then
        ln -sfn "/home/citrusml2004/pairwise_alignment" "$PAIRWISE_DIR"
    elif [ -d "$HOME/pairwise_alignment" ]; then
        ln -sfn "$HOME/pairwise_alignment" "$PAIRWISE_DIR"
    fi
fi

LOG_DIR="$AMPLIPHY_ROOT/logs/imp"
mkdir -p "$LOG_DIR"

REALTIME_LOG="$LOG_DIR/imp_tree_${LSB_JOBID:-master}.log"
exec > >(tee -a "$REALTIME_LOG") 2>&1
echo "リアルタイムログ開始: $(date) -> $REALTIME_LOG"

# ------------------------------------------------------------------------------
# 2. 環境有効化 (biotools: iqtree, python)
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

IQTREE_BIN=""
for cand in iqtree2 iqtree; do
    if command -v "$cand" &>/dev/null; then IQTREE_BIN="$(command -v "$cand")"; break; fi
done

PANJEP_BIN=""
for cand in \
    "$PAIRWISE_DIR/panjep/panjep" \
    "/home/citrusml2004/pairwise_alignment/panjep/panjep" \
    "$HOME/pairwise_alignment/panjep/panjep" \
    "$HOME/.local/bin/panjep" \
    "$(which panjep 2>/dev/null || true)"; do
    if [ -n "$cand" ] && [ -x "$cand" ]; then PANJEP_BIN="$cand"; break; fi
done

# ------------------------------------------------------------------------------
# 3. 設定
# ------------------------------------------------------------------------------
if [ "$#" -gt 0 ]; then
    TARGET_TAXA=("$@")
else
    TARGET_TAXA=("eukaryota" "amniota" "bacteria" "actinomycete")
fi
RESULT_ROOTS="${RESULT_ROOTS:-results results_max results_pwa}"
METHODS="${METHODS:-global local}"
IQTREE_THREADS="${IQTREE_THREADS:-4}"
PARALLEL_JOBS="${PARALLEL_JOBS:-4}"
# ori / enrich 木と同一の IQ-TREE オプション (ampliphy/nextflow.config)
IQTREE_OPTIONS="${IQTREE_OPTIONS:--m JTT+I+G4 -B 1000 -keep-ident}"

if [[ " $RESULT_ROOTS " == *" results "* || " $RESULT_ROOTS " == *" results_max "* ]]; then
    if [ -z "$IQTREE_BIN" ]; then
        echo "[ERROR] iqtree が見つかりません (biotools 環境を確認してください)" >&2
        exit 1
    fi
fi
if [[ " $RESULT_ROOTS " == *" results_pwa "* ]] && [ -z "$PANJEP_BIN" ]; then
    echo "[ERROR] panjep が見つかりません" >&2
    exit 1
fi

echo "=============================================================================="
echo "Impoverished tree 作成"
echo "Started at:     $(date)"
echo "Host:           $(hostname)"
echo "AMPLIPHY_ROOT:  $AMPLIPHY_ROOT"
echo "Taxa:           ${TARGET_TAXA[*]}"
echo "Result roots:   $RESULT_ROOTS"
echo "PWA methods:    $METHODS (panjep 手法は対象外)"
echo "IQ-TREE:        ${IQTREE_BIN:-N/A} ($IQTREE_OPTIONS, -T $IQTREE_THREADS)"
echo "panjep:         ${PANJEP_BIN:-N/A}"
echo "Parallel jobs:  $PARALLEL_JOBS"
echo "=============================================================================="

# ------------------------------------------------------------------------------
# 4. HOG 単位のワーカー (xargs -P から呼ばれる)
#    引数: kind(msa|pwa) root taxon method hog
#    戻り値: 0=成功/スキップ, 1=失敗
# ------------------------------------------------------------------------------
process_hog() {
    local kind="$1" root="$2" taxon="$3" method="$4" hog="$5"
    local ori_fa="$AMPLIPHY_ROOT/input_fasta/$taxon/${hog}.fa"
    local errlog="$LOG_DIR/${root}_${taxon}_${method}_${hog}.log"

    if [ ! -s "$ori_fa" ]; then
        echo "[SKIP] $root/$taxon/$hog: 元FASTAなし" ; return 0
    fi

    if [ "$kind" = "msa" ]; then
        local base="$AMPLIPHY_ROOT/$root/$taxon"
        local out_tree="$base/tree_imp/${hog}.imp.nwk"
        local out_msa="$base/msa_imp/${hog}.imp.fa"
        mkdir -p "$base/tree_imp" "$base/msa_imp"
        if [ -s "$out_tree" ]; then echo "[SKIP] $root/$taxon/$hog: 既存"; return 0; fi

        local src_msa
        if ! src_msa=$(python3 -m src.impoverished_tree select "$base/msa" "$hog" 2>"$errlog"); then
            echo "[SKIP] $root/$taxon/$hog: 増幅MSAなし"; return 0
        fi
        if ! python3 -m src.impoverished_tree msa "$src_msa" "$ori_fa" "$out_msa" >"$errlog" 2>&1; then
            echo "[FAIL] $root/$taxon/$hog: 元配列抽出失敗 ($(tail -n 1 "$errlog"))"; return 1
        fi

        local nseq
        nseq=$(grep -c '^>' "$out_msa" || true)
        if [ "$nseq" -lt 3 ]; then
            echo "[SKIP] $root/$taxon/$hog: 配列数 $nseq (<3) のため系統樹推定不可"; return 0
        fi

        local tmpd opts="$IQTREE_OPTIONS"
        tmpd=$(mktemp -d "${TMPDIR:-/tmp}/imp_${hog}_XXXXXX")
        # 配列数が少ない場合は ampliphy 本体と同様にブートストラップを無効化
        if [ "$nseq" -lt 4 ]; then
            opts=$(echo "$opts" | sed -E 's/(^|[[:space:]])(-B|--ufboot|-b|--boot)[[:space:]]+[^[:space:]]+//g; s/(^|[[:space:]])-bnni([[:space:]]|$)/ /g')
        fi
        # shellcheck disable=SC2086
        if "$IQTREE_BIN" -s "$out_msa" -pre "$tmpd/${hog}.imp" -T "$IQTREE_THREADS" $opts >"$errlog" 2>&1 \
           && [ -s "$tmpd/${hog}.imp.treefile" ]; then
            cp "$tmpd/${hog}.imp.treefile" "$out_tree"
            rm -rf "$tmpd"; rm -f "$errlog"
            echo "[OK]   $root/$taxon/$hog"; return 0
        fi
        rm -rf "$tmpd"
        echo "[FAIL] $root/$taxon/$hog: IQ-TREE 失敗 (ログ: $errlog)"; return 1

    else
        local base="$AMPLIPHY_ROOT/$root/$taxon/$method"
        local amp_phy="$base/matrix/${hog}.amp.phy"
        local imp_phy="$base/matrix/${hog}.imp.phy"
        local out_tree="$base/tree_imp/${hog}.imp.nwk"
        mkdir -p "$base/tree_imp"
        if [ -s "$out_tree" ]; then echo "[SKIP] $root/$taxon/$method/$hog: 既存"; return 0; fi
        if [ ! -s "$amp_phy" ]; then echo "[SKIP] $root/$taxon/$method/$hog: amp.phy なし"; return 0; fi
        if ! python3 -m src.impoverished_tree matrix "$amp_phy" "$ori_fa" "$imp_phy" >"$errlog" 2>&1; then
            echo "[FAIL] $root/$taxon/$method/$hog: 部分行列失敗 ($(tail -n 1 "$errlog"))"; return 1
        fi
        if "$PANJEP_BIN" -t 1 "$imp_phy" >"$out_tree" 2>>"$errlog" && [ -s "$out_tree" ]; then
            rm -f "$errlog"
            echo "[OK]   $root/$taxon/$method/$hog"; return 0
        fi
        rm -f "$out_tree"
        echo "[FAIL] $root/$taxon/$method/$hog: panjep 失敗 (ログ: $errlog)"; return 1
    fi
}

export -f process_hog
export AMPLIPHY_ROOT LOG_DIR IQTREE_BIN PANJEP_BIN IQTREE_THREADS IQTREE_OPTIONS

# ------------------------------------------------------------------------------
# 5. ジョブリスト作成と並列実行
# ------------------------------------------------------------------------------
JOB_LIST="$(mktemp "${TMPDIR:-/tmp}/imp_jobs_XXXXXX")"
trap 'rm -f "$JOB_LIST"' EXIT

for root in $RESULT_ROOTS; do
    for taxon in "${TARGET_TAXA[@]}"; do
        ori_dir="$AMPLIPHY_ROOT/input_fasta/$taxon"
        [ -d "$ori_dir" ] || { echo "[WARN] $ori_dir がありません。スキップ"; continue; }

        if [ "$root" = "results_pwa" ]; then
            for method in $METHODS; do
                [ "$method" = "panjep" ] && { echo "[INFO] panjep 手法は対象外"; continue; }
                [ -d "$AMPLIPHY_ROOT/$root/$taxon/$method/matrix" ] || { echo "[WARN] $root/$taxon/$method/matrix なし"; continue; }
                for f in "$ori_dir"/*.fa; do
                    printf '%s %s %s %s %s\n' pwa "$root" "$taxon" "$method" "$(basename "${f%.fa}")" >> "$JOB_LIST"
                done
            done
        else
            [ -d "$AMPLIPHY_ROOT/$root/$taxon/msa" ] || { echo "[WARN] $root/$taxon/msa なし"; continue; }
            for f in "$ori_dir"/*.fa; do
                printf '%s %s %s %s %s\n' msa "$root" "$taxon" "-" "$(basename "${f%.fa}")" >> "$JOB_LIST"
            done
        fi
    done
done

TOTAL=$(wc -l < "$JOB_LIST")
echo "総ジョブ数: $TOTAL"

RESULT_LOG="$LOG_DIR/imp_tree_results_${LSB_JOBID:-master}.txt"
: > "$RESULT_LOG"
# xargs は 1 つでも失敗すると非0を返すため、集計は結果ログから行う
xargs -P "$PARALLEL_JOBS" -L 1 bash -c 'process_hog "$@"' _ < "$JOB_LIST" | tee -a "$RESULT_LOG" || true

N_OK=$(grep -c '^\[OK\]' "$RESULT_LOG" || true)
N_SKIP=$(grep -c '^\[SKIP\]' "$RESULT_LOG" || true)
N_FAIL=$(grep -c '^\[FAIL\]' "$RESULT_LOG" || true)

echo ""
echo "=============================================================================="
echo "完了: $(date)"
echo "  成功: $N_OK / スキップ: $N_SKIP / 失敗: $N_FAIL  (総数: $TOTAL)"
echo "  失敗一覧: grep '^\[FAIL\]' $RESULT_LOG"
echo "=============================================================================="
[ "$N_FAIL" -eq 0 ]

