from __future__ import annotations

# -*- coding: utf-8 -*-
"""Graph Splitting (GS) Method Evaluation Pipeline Module.

This module provides core utilities to prepare sequences, execute the GS method
(gs2), perform tree cleaning, handle leaf pruning, and compute TCS (Tree Certainty
Score) improvement metrics comparing original trees and homolog-added trees.
"""

import math
import os
import re
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple, Union

if TYPE_CHECKING:
    import pandas as pd


def prepare_unaligned_fasta(
    input_fasta: Union[str, Path],
    output_fasta: Union[str, Path],
) -> int:
    """Read a FASTA file, remove gap characters, and write unaligned sequences.

    Parameters
    ----------
    input_fasta : Union[str, Path]
        Path to the input FASTA file (can be aligned or unaligned). (unitless)
    output_fasta : Union[str, Path]
        Path to the output unaligned FASTA file. (unitless)

    Returns
    -------
    int
        The total count of valid sequences processed and written. (count)

    Raises
    ------
    FileNotFoundError
        If input_fasta does not exist.
    ValueError
        If no valid FASTA sequences are found in input_fasta.
    """
    in_path = Path(input_fasta)
    out_path = Path(output_fasta)

    if not in_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {in_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    seq_count: int = 0
    current_header: Optional[str] = None
    current_seq_parts: List[str] = []

    with open(in_path, "r", encoding="utf-8") as in_f, open(
        out_path, "w", encoding="utf-8"
    ) as out_f:
        for line in in_f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    full_seq = "".join(current_seq_parts)
                    clean_seq = re.sub(r"[\s\-.]", "", full_seq)
                    if clean_seq:
                        out_f.write(f"{current_header}\n{clean_seq}\n")
                        seq_count += 1
                # Standardize header: take only first token up to whitespace to avoid parser issues
                header_token = line[1:].split()[0]
                current_header = f">{header_token}"
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

        # Flush final sequence
        if current_header is not None:
            full_seq = "".join(current_seq_parts)
            clean_seq = re.sub(r"[\s\-.]", "", full_seq)
            if clean_seq:
                out_f.write(f"{current_header}\n{clean_seq}\n")
                seq_count += 1

    if seq_count == 0:
        raise ValueError(f"No valid sequences found in {in_path}")

    return seq_count


def extract_fasta_sequence_ids(fasta_path: Union[str, Path]) -> List[str]:
    """Extract sequence IDs (first token after '>') from a FASTA file.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to the FASTA file. (unitless)

    Returns
    -------
    List[str]
        List of unique sequence identifier strings in order of appearance. (unitless)
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    ids: List[str] = []
    seen: Set[str] = set()

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id and seq_id not in seen:
                    seen.add(seq_id)
                    ids.append(seq_id)
    return ids


def extract_homolog_ids(
    orig_fasta: Union[str, Path],
    amp_fasta: Union[str, Path],
) -> List[str]:
    """Identify added homolog IDs by comparing amplified and original FASTA files.

    Parameters
    ----------
    orig_fasta : Union[str, Path]
        Path to the original FASTA file containing baseline sequences. (unitless)
    amp_fasta : Union[str, Path]
        Path to the amplified FASTA file containing original and added homologs. (unitless)

    Returns
    -------
    List[str]
        List of sequence IDs present in amp_fasta but not in orig_fasta. (unitless)
    """
    orig_ids = set(extract_fasta_sequence_ids(orig_fasta))
    amp_ids = extract_fasta_sequence_ids(amp_fasta)

    added_homologs: List[str] = [sid for sid in amp_ids if sid not in orig_ids]
    return added_homologs


def build_gs_command(
    fasta_path: Union[str, Path],
    output_nwk: Union[str, Path],
    gs2_bin: str = "gs2",
    threads: int = 4,
    mmseqs_sensitivity: float = 7.5,
    name_output: bool = True,
    ep_replicates: int = 0,
) -> str:
    """Build the command line invocation string for gs2.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to the input FASTA file. (unitless)
    output_nwk : Union[str, Path]
        Path to the output Newick tree file. (unitless)
    gs2_bin : str, optional
        Path or command name for the gs2 executable, by default "gs2". (unitless)
    threads : int, optional
        Number of CPU threads for MMseqs2 pairwise alignment, by default 4. (count)
    mmseqs_sensitivity : float, optional
        Sensitivity parameter for MMseqs2 (range 1.0 to 7.5), by default 7.5. (unitless)
    name_output : bool, optional
        If True, use '-l' flag to output sequence names instead of IDs, by default True.
    ep_replicates : int, optional
        Number of Edge Perturbation (EP) replicates for branch supports, by default 0. (count)

    Returns
    -------
    str
        The shell command string to execute gs2. (unitless)
    """
    flags: List[str] = []
    if name_output:
        flags.append("-l")
    flags.append(f"-t {threads}")
    flags.append(f"-m {mmseqs_sensitivity}")
    if ep_replicates > 0:
        flags.append(f"-e {ep_replicates}")

    flags_str = " ".join(flags)
    cmd = f"{gs2_bin} {flags_str} {fasta_path} > {output_nwk}"
    return cmd


def clean_newick_tree(
    tree_nwk_path: Union[str, Path],
    output_clean_path: Optional[Union[str, Path]] = None,
) -> str:
    """Clean Newick tree format to prevent parsing errors in MAD, GoTree, and fold_tree.

    Strips single and double quotation marks from node labels and ensures
    a valid trailing semicolon.

    Parameters
    ----------
    tree_nwk_path : Union[str, Path]
        Path to the raw Newick tree file. (unitless)
    output_clean_path : Optional[Union[str, Path]], optional
        Path to write cleaned Newick tree. If None, overwrites tree_nwk_path. (unitless)

    Returns
    -------
    str
        Cleaned Newick tree string. (unitless)
    """
    in_path = Path(tree_nwk_path)
    if not in_path.exists():
        raise FileNotFoundError(f"Tree file not found: {in_path}")

    with open(in_path, "r", encoding="utf-8") as f:
        content = f.read().strip()

    # Extract first line if multi-line (e.g. from MAD rooted outputs with multiple tied roots)
    first_line = content.splitlines()[0].strip() if content else ""

    # Remove single and double quotes
    cleaned = first_line.replace("'", "").replace('"', "")

    # Ensure semicolon at the end
    if cleaned and not cleaned.endswith(";"):
        cleaned += ";"

    target_path = Path(output_clean_path) if output_clean_path else in_path
    target_path.parent.mkdir(parents=True, exist_ok=True)
    with open(target_path, "w", encoding="utf-8") as f:
        f.write(f"{cleaned}\n")

    return cleaned


def calculate_delta_tcs(ori_score: float, enr_score: float) -> Dict[str, float]:
    """Calculate TCS score ratio and difference metrics between original and enriched trees.

    Parameters
    ----------
    ori_score : float
        TCS score of original (augmented/unmodified) tree. (unitless)
    enr_score : float
        TCS score of enriched (homolog-amplified & pruned) tree. (unitless)

    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'delta_ratio': enr_score / ori_score (unitless)
        - 'delta_diff': enr_score - ori_score (unitless)
        - 'is_improved': 1.0 if enr_score > ori_score else 0.0 (unitless)
    """
    if ori_score > 0.0:
        ratio = enr_score / ori_score
    else:
        ratio = float("nan")

    diff = enr_score - ori_score
    is_improved = 1.0 if enr_score > ori_score else 0.0

    return {
        "delta_ratio": ratio,
        "delta_diff": diff,
        "is_improved": is_improved,
    }


def parse_gs_tcs_scores(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Parse TCS score TSV file generated by tcs_batch.py for GS evaluation.

    The input TSV contains lines of:
    <filepath>\t<score>

    Trees are classified into:
    - 'ori' (or 'aug'): original tree
    - 'enr': enriched tree

    Parameters
    ----------
    tsv_path : Union[str, Path]
        Path to the TCS score TSV file. (unitless)

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by hog_id with columns:
        - 'tcs_ori': float (unitless)
        - 'tcs_enr': float (unitless)
        - 'delta_ratio': float, tcs_enr / tcs_ori (unitless)
        - 'delta_diff': float, tcs_enr - tcs_ori (unitless)
        - 'is_improved': bool (unitless)
    """
    import pandas as pd

    path = Path(tsv_path)
    if not path.exists():
        raise FileNotFoundError(f"TCS TSV file not found: {path}")

    scores: Dict[str, Dict[str, float]] = {}

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            file_str, score_str = parts[0], parts[1]
            try:
                score_val = float(score_str)
            except ValueError:
                continue

            fname = os.path.basename(file_str)
            # Identify HOG ID
            # Matches formats like HOG_E0736661 or E0736661
            match = re.search(r"(HOG_[A-Za-z0-9]+|[A-Za-z]\d{6,})", fname)
            if not match:
                continue
            hog_id = match.group(1)

            if hog_id not in scores:
                scores[hog_id] = {}

            # Classify tree type
            lower_str = file_str.lower()
            if "enr" in lower_str or "amp" in lower_str:
                scores[hog_id]["enr"] = score_val
            elif "aug" in lower_str or "ori" in lower_str:
                scores[hog_id]["ori"] = score_val

    rows: List[Dict[str, Union[str, float, bool]]] = []
    for hog_id, d in scores.items():
        if "ori" in d and "enr" in d:
            ori_s = d["ori"]
            enr_s = d["enr"]
            metrics = calculate_delta_tcs(ori_s, enr_s)
            rows.append(
                {
                    "hog_id": hog_id,
                    "tcs_ori": ori_s,
                    "tcs_enr": enr_s,
                    "delta_ratio": metrics["delta_ratio"],
                    "delta_diff": metrics["delta_diff"],
                    "is_improved": bool(metrics["is_improved"]),
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values(by="hog_id", inplace=True)
        df.reset_index(drop=True, inplace=True)
    return df


def summarize_gs_results(df: pd.DataFrame) -> Dict[str, Union[float, int]]:
    """Compute overall summary statistics of GS method evaluation across HOGs.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame returned by parse_gs_tcs_scores. (unitless)

    Returns
    -------
    Dict[str, Union[float, int]]
        Statistical summary dictionary containing:
        - 'n_hogs': int, number of evaluated HOGs (count)
        - 'mean_tcs_ori': float, mean TCS of original trees (unitless)
        - 'mean_tcs_enr': float, mean TCS of enriched trees (unitless)
        - 'mean_delta_ratio': float, mean ratio of TCS improvement (unitless)
        - 'median_delta_ratio': float, median ratio of TCS improvement (unitless)
        - 'mean_delta_diff': float, mean absolute difference in TCS (unitless)
        - 'pct_improved': float, percentage of HOGs showing improved TCS (%)
        - 'wilcoxon_pvalue': float, p-value from Wilcoxon signed-rank test (unitless)
    """
    if df.empty:
        return {
            "n_hogs": 0,
            "mean_tcs_ori": 0.0,
            "mean_tcs_enr": 0.0,
            "mean_delta_ratio": 0.0,
            "median_delta_ratio": 0.0,
            "mean_delta_diff": 0.0,
            "pct_improved": 0.0,
            "wilcoxon_pvalue": 1.0,
        }

    n = len(df)
    mean_ori = float(df["tcs_ori"].mean())
    mean_enr = float(df["tcs_enr"].mean())
    mean_ratio = float(df["delta_ratio"].mean())
    median_ratio = float(df["delta_ratio"].median())
    mean_diff = float(df["delta_diff"].mean())
    pct_imp = float((df["is_improved"].sum() / n) * 100.0)

    # Wilcoxon signed-rank test between enriched and original scores
    try:
        from scipy import stats
        diffs = df["tcs_enr"] - df["tcs_ori"]
        if (diffs != 0).any():
            stat_res = stats.wilcoxon(df["tcs_enr"], df["tcs_ori"])
            pval = float(stat_res.pvalue)
        else:
            pval = 1.0
    except Exception:
        pval = 1.0

    return {
        "n_hogs": n,
        "mean_tcs_ori": mean_ori,
        "mean_tcs_enr": mean_enr,
        "mean_delta_ratio": mean_ratio,
        "median_delta_ratio": median_ratio,
        "mean_delta_diff": mean_diff,
        "pct_improved": pct_imp,
        "wilcoxon_pvalue": pval,
    }


def generate_lineage_tsv(
    taxa_tsv: Union[str, Path],
    fasta_files: List[Union[str, Path]],
    output_tsv: Union[str, Path],
    species_code_len: int = 5,
) -> int:
    """Generate a lineage TSV mapping query ID to lineage using species taxonomy metadata.

    Parameters
    ----------
    taxa_tsv : Union[str, Path]
        Path to the species taxonomy metadata TSV file (e.g. euk_taxa.tsv). (unitless)
    fasta_files : List[Union[str, Path]]
        List of paths to FASTA files to extract sequence query IDs from. (unitless)
    output_tsv : Union[str, Path]
        Path to the output lineage TSV file. (unitless)
    species_code_len : int, optional
        Length of the UniProt species mnemonic code prefix, by default 5. (count)

    Returns
    -------
    int
        Total number of unique query-to-lineage mappings written. (count)

    Raises
    ------
    FileNotFoundError
        If taxa_tsv does not exist.
    """
    taxa_path = Path(taxa_tsv)
    out_path = Path(output_tsv)

    if not taxa_path.exists():
        raise FileNotFoundError(f"Taxa metadata file not found: {taxa_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Read mapping: 5-letter species code -> lineage string
    taxa_map: Dict[str, str] = {}
    with open(taxa_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                taxa_map[parts[0].strip()] = parts[2].strip()

    written: int = 0
    seen_queries: Set[str] = set()

    with open(out_path, "w", encoding="utf-8") as out_f:
        out_f.write("query\tlineage\n")
        for fasta_file in fasta_files:
            p = Path(fasta_file)
            if not p.exists():
                continue
            with open(p, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith(">"):
                        query = line[1:].strip().split()[0]
                        if query in seen_queries:
                            continue
                        seen_queries.add(query)
                        sp_code = query[:species_code_len]
                        if sp_code in taxa_map:
                            out_f.write(f"{query}\t{taxa_map[sp_code]}\n")
                            written += 1

    return written
