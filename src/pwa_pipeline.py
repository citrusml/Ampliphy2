#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pwa_pipeline.py
===============
Core module for Pairwise Sequence Alignment (PSA/PWA) phylogenetic pipeline.

Features:
- Pairwise alignment (Global via Needleman-Wunsch, Local via MMseqs2)
- Evolutionary Poisson distance matrix calculation (Zuckerkandl & Pauling, 1965)
- PHYLIP distance matrix export
- Negative branch length cleaning for NJ trees
- Sequence preparation and homolog merging
- Lineage extraction for Taxonomy Overlap Score (TCS)
- Results summary and statistics calculation
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Union

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd


# ------------------------------------------------------------------------------
# 1. Sequence & lookup tables
# ------------------------------------------------------------------------------

def get_amino_acid_table() -> List[int]:
    """
    Lookup table for 20 standard amino acids.

    Returns
    -------
    table : List[int]
        Lookup table of length 256 mapping ASCII values to 0..19 (-1 for invalid/gap).
        Physical unit: unitless.
    """
    table = [-1] * 256
    aa_str = "ARNDCQEGHILKMFPSTWYV"
    for i, aa in enumerate(aa_str):
        table[ord(aa)] = i
        table[ord(aa.lower())] = i
    return table


def get_nucleotide_table() -> List[int]:
    """
    Lookup table for 4 nucleotides (ACGT/U).

    Returns
    -------
    table : List[int]
        Lookup table of length 256 mapping ASCII values to 0..3 (-1 for invalid/gap).
        Physical unit: unitless.
    """
    table = [-1] * 256
    table[ord("A")] = table[ord("a")] = 0
    table[ord("C")] = table[ord("c")] = 1
    table[ord("G")] = table[ord("g")] = 2
    table[ord("T")] = table[ord("t")] = 3
    table[ord("U")] = table[ord("u")] = 3
    return table


def detect_is_nucleotide(fasta_path: Union[str, Path]) -> bool:
    """
    Detect whether the input FASTA consists of nucleotide or protein sequences.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.

    Returns
    -------
    is_nucl : bool
        True if >= 90% of characters are nucleotide characters.
        Physical unit: unitless.
    """
    total = 0
    nucl = 0
    nucl_chars = set("acgtnuwsmkrybdhv")
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                continue
            for char in line.strip():
                if char in ("-", ".", " "):
                    continue
                total += 1
                if char.lower() in nucl_chars:
                    nucl += 1
            if total > 50000:
                break
    return total > 0 and (nucl * 100 // total) >= 90


# ------------------------------------------------------------------------------
# 2. Distance kernels
# ------------------------------------------------------------------------------

def calculate_poisson_distance(
    aligned_seq1: str,
    aligned_seq2: str,
    is_nucleotide: bool = False,
    method: str = "ver2",
) -> float:
    """
    Calculate evolutionary distance (substitutions per site) between two aligned sequences.

    Parameters
    ----------
    aligned_seq1 : str
        First aligned sequence string (containing gap characters '-').
    aligned_seq2 : str
        Second aligned sequence string (containing gap characters '-').
    is_nucleotide : bool, default=False
        True if sequence is DNA/RNA, False if protein.
    method : str, default='ver2'
        Distance method ('ver1' for raw Poisson, 'ver2' for multi-hit corrected Poisson).

    Returns
    -------
    distance : float
        Evolutionary distance clamped to max 5.0.
        Physical unit: substitutions per site.
    """
    if method == "ver1":
        # Raw Poisson correction: d = -ln(1 - p)
        mismatches = 0
        aligned_len = 0
        min_len = min(len(aligned_seq1), len(aligned_seq2))
        for idx in range(min_len):
            c1 = aligned_seq1[idx]
            c2 = aligned_seq2[idx]
            if c1 in ("-", ".") or c2 in ("-", "."):
                continue
            aligned_len += 1
            if c1.upper() != c2.upper():
                mismatches += 1
        if aligned_len == 0:
            return 5.0
        p = mismatches / aligned_len
        if p >= 0.99:
            return 5.0
        return -math.log(1.0 - p)

    # ver2 / panjep multi-hit corrected Poisson (Zuckerkandl & Pauling, 1965)
    nc = 4 if is_nucleotide else 20
    tab = get_nucleotide_table() if is_nucleotide else get_amino_acid_table()
    total = 0
    matches = 0
    min_len = min(len(aligned_seq1), len(aligned_seq2))
    for idx in range(min_len):
        c1 = aligned_seq1[idx]
        c2 = aligned_seq2[idx]
        o1 = ord(c1)
        o2 = ord(c2)
        if o1 >= 256 or o2 >= 256:
            continue
        pa = tab[o1]
        pb = tab[o2]
        if pa < 0 or pb < 0:
            continue
        total += 1
        if pa == pb:
            matches += 1

    if total == 0:
        return 5.0
    s = 1.0 - float(matches) / float(total)
    if s <= 0.0:
        return 0.0
    lim = float(nc - 1) / float(nc)
    if s >= lim:
        return 5.0
    dist = -lim * math.log(1.0 - s / lim)
    return min(dist, 5.0)


# ------------------------------------------------------------------------------
# 3. Pairwise alignment and distance matrix generation
# ------------------------------------------------------------------------------

def parse_fasta(fasta_path: Union[str, Path]) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file into a list of (id, sequence) tuples.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.

    Returns
    -------
    records : List[Tuple[str, str]]
        List of (seq_id, sequence).
    """
    records: List[Tuple[str, str]] = []
    current_id: Optional[str] = None
    current_seq: List[str] = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            records.append((current_id, "".join(current_seq)))
    return records


def run_mmseqs_local(
    fasta_path: Union[str, Path],
    temp_dir: str,
    threads: int = 4,
    sensitivity: float = 7.5,
    is_nucleotide: bool = False,
) -> Dict[Tuple[str, str], Tuple[str, str, float]]:
    """
    Execute all-vs-all local alignment with MMseqs2.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to FASTA file.
    temp_dir : str
        Temporary directory path.
    threads : int, default=4
        Number of threads.
    sensitivity : float, default=7.5
        MMseqs2 search sensitivity.
    is_nucleotide : bool, default=False
        True if nucleotide.

    Returns
    -------
    best_alns : Dict[Tuple[str, str], Tuple[str, str, float]]
        Mapping of (query_id, target_id) to (qaln, taln, bitscore).
    """
    if not shutil.which("mmseqs"):
        raise RuntimeError("mmseqs command not found in PATH.")

    db = os.path.join(temp_dir, "db")
    res = os.path.join(temp_dir, "result")
    tmp = os.path.join(temp_dir, "tmp")
    m8 = os.path.join(temp_dir, "result.m8")
    os.makedirs(tmp, exist_ok=True)

    cmd_createdb = ["mmseqs", "createdb", str(fasta_path), db]
    if is_nucleotide:
        cmd_createdb += ["--dbtype", "2"]
    subprocess.run(cmd_createdb, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd_search = [
        "mmseqs", "search", db, db, res, tmp,
        "--threads", str(threads), "-e", "10.0", "-s", f"{sensitivity:.1f}", "-a", "1"
    ]
    if is_nucleotide:
        cmd_search += ["--search-type", "3"]
    subprocess.run(cmd_search, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd_convert = [
        "mmseqs", "convertalis", db, db, res, m8,
        "--format-output", "query,target,qaln,taln,bits"
    ]
    subprocess.run(cmd_convert, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    best_alns: Dict[Tuple[str, str], Tuple[str, str, float]] = {}
    if os.path.exists(m8) and os.path.getsize(m8) > 0:
        with open(m8, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                q, t, qaln, taln = parts[0], parts[1], parts[2], parts[3]
                try:
                    bitscore = float(parts[4])
                except ValueError:
                    bitscore = 0.0
                pair = (q, t)
                if pair not in best_alns or bitscore > best_alns[pair][2]:
                    best_alns[pair] = (qaln, taln, bitscore)

    return best_alns


# BLOSUM62 matrix for fallback standalone Needleman-Wunsch alignment
BLOSUM62_AA = "ARNDCQEGHILKMFPSTWYVBZX*"
BLOSUM62_SCORES = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
    [-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
    [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1],
]
BLOSUM62_MAP = {c: i for i, c in enumerate(BLOSUM62_AA)}


def needleman_wunsch_align(
    s1: str, s2: str, is_nucleotide: bool = False, gap_open: int = -11, gap_extend: int = -1
) -> Tuple[str, str]:
    """
    Standalone fallback Needleman-Wunsch global alignment.

    Parameters
    ----------
    s1 : str
        First sequence string.
    s2 : str
        Second sequence string.
    is_nucleotide : bool, default=False
        True if nucleotide.
    gap_open : int, default=-11
        Gap open penalty.
    gap_extend : int, default=-1
        Gap extend penalty.

    Returns
    -------
    aligned_s1, aligned_s2 : Tuple[str, str]
        Pair of aligned sequence strings.
    """
    n, m = len(s1), len(s2)
    if n == 0 and m == 0:
        return "", ""
    if n == 0:
        return "-" * m, s2
    if m == 0:
        return s1, "-" * n

    gap_penalty = gap_open  # linear approximation for robust fast fallback
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty

    if is_nucleotide:
        def score_fn(a: str, b: str) -> int:
            return 1 if a.upper() == b.upper() else -3
    else:
        def score_fn(a: str, b: str) -> int:
            ia = BLOSUM62_MAP.get(a.upper(), 22)
            ib = BLOSUM62_MAP.get(b.upper(), 22)
            return BLOSUM62_SCORES[ia][ib]

    for i in range(1, n + 1):
        c1 = s1[i - 1]
        for j in range(1, m + 1):
            c2 = s2[j - 1]
            s = score_fn(c1, c2)
            dp[i][j] = max(dp[i - 1][j - 1] + s, dp[i - 1][j] + gap_penalty, dp[i][j - 1] + gap_penalty)

    i, j = n, m
    a1, a2 = [], []
    while i > 0 or j > 0:
        s = score_fn(s1[i - 1], s2[j - 1]) if (i > 0 and j > 0) else -9999
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + s:
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
        else:
            a1.append("-")
            a2.append(s2[j - 1])
            j -= 1

    return "".join(reversed(a1)), "".join(reversed(a2))


def compute_pwa_matrix(
    fasta_path: Union[str, Path],
    method: str = "ver2",
    threads: int = 4,
    sensitivity: float = 7.5,
) -> Tuple[List[str], List[List[float]]]:
    """
    Compute pairwise evolutionary distance matrix from FASTA.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to input FASTA file.
    method : str, default='ver2'
        'ver2' for global alignment (Needleman-Wunsch + Poisson),
        'ver3' for local alignment (MMseqs2 + Poisson).
    threads : int, default=4
        Number of threads for parallel alignment.
    sensitivity : float, default=7.5
        Sensitivity for MMseqs2 search.

    Returns
    -------
    seq_ids : List[str]
        List of sequence identifiers.
    matrix : List[List[float]]
        Full symmetric n x n distance matrix.
        Physical unit: substitutions per site.
    """
    records = parse_fasta(fasta_path)
    n = len(records)
    if n == 0:
        raise ValueError(f"No sequences found in {fasta_path}")

    seq_ids = [r[0] for r in records]
    seq_dict = {r[0]: r[1] for r in records}
    is_nucl = detect_is_nucleotide(fasta_path)

    # Initialize n x n matrix with 5.0 (saturation max), diagonal = 0.0
    matrix = [[0.0 if i == j else 5.0 for j in range(n)] for i in range(n)]

    if method in ("ver1", "ver2"):
        # Global pairwise alignment
        use_biopython = False
        aligner = None
        try:
            from Bio.Align import PairwiseAligner, substitution_matrices
            aligner = PairwiseAligner()
            aligner.mode = "global"
            if is_nucl:
                aligner.match_score = 1.0
                aligner.mismatch_score = -3.0
                aligner.open_gap_score = -5.0
                aligner.extend_gap_score = -2.0
            else:
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -11.0
                aligner.extend_gap_score = -1.0
            use_biopython = True
        except ImportError:
            use_biopython = False

        for i in range(n):
            for j in range(i + 1, n):
                s1 = seq_dict[seq_ids[i]]
                s2 = seq_dict[seq_ids[j]]
                if use_biopython and aligner is not None:
                    alns = aligner.align(s1, s2)
                    best_aln = alns[0]
                    aln1, aln2 = best_aln[0], best_aln[1]
                else:
                    aln1, aln2 = needleman_wunsch_align(s1, s2, is_nucleotide=is_nucl)
                d = calculate_poisson_distance(aln1, aln2, is_nucleotide=is_nucl, method=method)
                matrix[i][j] = d
                matrix[j][i] = d
                matrix[i][j] = d
                matrix[j][i] = d

    elif method in ("ver3", "panjep", "local"):
        # Local pairwise alignment via MMseqs2
        with tempfile.TemporaryDirectory() as temp_dir:
            best_alns = run_mmseqs_local(
                fasta_path, temp_dir, threads=threads, sensitivity=sensitivity, is_nucleotide=is_nucl
            )

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                pair = (seq_ids[i], seq_ids[j])
                if pair in best_alns:
                    qaln, taln, _ = best_alns[pair]
                    d = calculate_poisson_distance(qaln, taln, is_nucleotide=is_nucl, method="ver2")
                    matrix[i][j] = d

    return seq_ids, matrix


def save_as_phylip(
    seq_ids: List[str],
    matrix: List[List[float]],
    output_path: Union[str, Path],
) -> None:
    """
    Save distance matrix in standard PHYLIP format.

    Parameters
    ----------
    seq_ids : List[str]
        List of sequence identifiers.
    matrix : List[List[float]]
        n x n distance matrix.
    output_path : Union[str, Path]
        Path to output file.
    """
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"{len(seq_ids)}\n")
        for i, name in enumerate(seq_ids):
            row_str = " ".join(f"{x:.6f}" for x in matrix[i])
            f.write(f"{name} {row_str}\n")


# ------------------------------------------------------------------------------
# 4. Sequence merging and tree post-processing
# ------------------------------------------------------------------------------

def prepare_amplified_fasta(
    ori_fasta: Union[str, Path],
    homolog_fasta: Union[str, Path],
    output_fasta: Union[str, Path],
) -> int:
    """
    Merge original sequence FASTA and homolog FASTA into a single amplified FASTA.

    If homolog FASTA does not exist or is empty, the original FASTA is copied.

    Parameters
    ----------
    ori_fasta : Union[str, Path]
        Path to original FASTA.
    homolog_fasta : Union[str, Path]
        Path to homolog FASTA.
    output_fasta : Union[str, Path]
        Path to output merged FASTA.

    Returns
    -------
    num_seqs : int
        Number of sequences in the output FASTA.
        Physical unit: count.
    """
    ori_records = parse_fasta(ori_fasta)
    hom_records: List[Tuple[str, str]] = []
    if os.path.exists(homolog_fasta) and os.path.getsize(homolog_fasta) > 0:
        hom_records = parse_fasta(homolog_fasta)

    all_records = ori_records + hom_records
    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)
    with open(output_fasta, "w", encoding="utf-8") as f:
        for seq_id, seq in all_records:
            f.write(f">{seq_id}\n{seq}\n")

    return len(all_records)


def clean_negative_branches(
    tree_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
) -> str:
    """
    Replace negative branch lengths in Newick tree with :0.0 and strip single quotes.

    Parameters
    ----------
    tree_path : Union[str, Path]
        Path to input Newick tree.
    output_path : Optional[Union[str, Path]], default=None
        Path to output Newick tree. If None, overwrites input file.

    Returns
    -------
    cleaned_tree : str
        Cleaned Newick tree string.
    """
    with open(tree_path, "r", encoding="utf-8") as f:
        tree_str = f.read().strip()

    # Strip single quotes and double quotes around taxon labels
    tree_str = tree_str.replace("'", "").replace('"', "")

    # Replace negative branch lengths with 0.0 (e.g. :-0.005 -> :0.0)
    tree_str = re.sub(r":-[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?", ":0.0", tree_str)

    # Ensure ending semicolon
    if not tree_str.endswith(";"):
        tree_str += ";"

    dest = output_path or tree_path
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    with open(dest, "w", encoding="utf-8") as f:
        f.write(tree_str + "\n")

    return tree_str


# ------------------------------------------------------------------------------
# 5. Taxonomy lineage & TCS score evaluation
# ------------------------------------------------------------------------------

def generate_lineage_tsv(
    taxa_tsv: Union[str, Path],
    fasta_files: List[Union[str, Path]],
    output_tsv: Union[str, Path],
) -> int:
    """
    Generate lineage TSV file mapping query ID to taxonomic lineage for TCS evaluation.

    Parameters
    ----------
    taxa_tsv : Union[str, Path]
        Path to taxa.tsv containing organism code and taxonomy lineage.
    fasta_files : List[Union[str, Path]]
        List of input FASTA files.
    output_tsv : Union[str, Path]
        Path to output lineage.tsv.

    Returns
    -------
    matched_count : int
        Number of matched query sequence records.
        Physical unit: count.
    """
    org_to_lineage: Dict[str, str] = {}
    with open(taxa_tsv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                org_to_lineage[parts[0].strip()] = parts[2].strip()

    seen_ids: Set[str] = set()
    matched_lines: List[str] = []

    for fasta in fasta_files:
        for seq_id, _ in parse_fasta(fasta):
            if seq_id in seen_ids:
                continue
            seen_ids.add(seq_id)
            # Match 5-letter or alphanumeric organism prefix
            matched = False
            for org, lin in org_to_lineage.items():
                if seq_id.startswith(org):
                    matched_lines.append(f"{seq_id}\t{lin}")
                    matched = True
                    break
            if not matched:
                prefix = seq_id[:5]
                if prefix in org_to_lineage:
                    matched_lines.append(f"{seq_id}\t{org_to_lineage[prefix]}")

    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)
    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("query\tlineage\n")
        for line in matched_lines:
            f.write(line + "\n")

    return len(matched_lines)


def calculate_delta_tcs(ori_score: float, enr_score: float) -> Dict[str, float]:
    """
    Calculate Delta TCS and TCS ratio.

    Parameters
    ----------
    ori_score : float
        TCS score of original tree. Physical unit: unitless [0.0 - 1.0].
    enr_score : float
        TCS score of amplified/enriched tree. Physical unit: unitless [0.0 - 1.0].

    Returns
    -------
    metrics : Dict[str, float]
        Dictionary with 'delta_tcs' (enr - ori) and 'ratio_tcs' (enr / ori).
    """
    delta = enr_score - ori_score
    ratio = enr_score / ori_score if ori_score > 0.0 else float("nan")
    return {"delta_tcs": delta, "ratio_tcs": ratio}


def parse_pwa_tcs_scores(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """
    Parse TCS score TSV file into a pandas DataFrame.

    Parameters
    ----------
    tsv_path : Union[str, Path]
        Path to TCS report TSV file.

    Returns
    -------
    df : pd.DataFrame
        DataFrame with columns ['family_id', 'tcs_ori', 'tcs_amp', 'delta_tcs', 'ratio_tcs'].
    """
    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t")
    if "delta_tcs" not in df.columns:
        df["delta_tcs"] = df["tcs_amp"] - df["tcs_ori"]
    if "ratio_tcs" not in df.columns:
        df["ratio_tcs"] = df.apply(
            lambda r: r["tcs_amp"] / r["tcs_ori"] if r["tcs_ori"] > 0 else float("nan"),
            axis=1,
        )
    return df


def summarize_pwa_results(df: pd.DataFrame) -> Dict[str, Union[float, int]]:
    """
    Compute statistical summary of TCS improvement metrics.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'delta_tcs' and 'ratio_tcs'.

    Returns
    -------
    summary : Dict[str, Union[float, int]]
        Summary dictionary with mean, std, median, improved count, p-value.
    """
    import numpy as np
    from scipy import stats

    valid_delta = df["delta_tcs"].dropna()
    valid_ratio = df["ratio_tcs"].dropna()
    n_total = len(df)
    n_improved = int((df["delta_tcs"] > 0).sum())
    n_declined = int((df["delta_tcs"] < 0).sum())
    n_unchanged = int((df["delta_tcs"] == 0).sum())

    mean_delta = float(valid_delta.mean()) if len(valid_delta) > 0 else float("nan")
    std_delta = float(valid_delta.std()) if len(valid_delta) > 1 else 0.0
    mean_ratio = float(valid_ratio.mean()) if len(valid_ratio) > 0 else float("nan")
    std_ratio = float(valid_ratio.std()) if len(valid_ratio) > 1 else 0.0

    p_val = float("nan")
    if len(valid_delta) >= 5 and np.any(valid_delta != 0):
        try:
            _, p_val = stats.ttest_1samp(valid_delta, 0.0)
            p_val = float(p_val)
        except Exception:
            p_val = float("nan")

    return {
        "total_hogs": n_total,
        "improved_count": n_improved,
        "declined_count": n_declined,
        "unchanged_count": n_unchanged,
        "mean_delta_tcs": mean_delta,
        "std_delta_tcs": std_delta,
        "mean_ratio_tcs": mean_ratio,
        "std_ratio_tcs": std_ratio,
        "p_value": p_val,
    }


# ------------------------------------------------------------------------------
# 6. CLI entrypoint
# ------------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="AmpliPhy PWA Pipeline Utility")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: matrix
    p_mat = subparsers.add_parser("matrix", help="Compute PWA distance matrix")
    p_mat.add_argument("input_fasta", help="Input FASTA")
    p_mat.add_argument("output_matrix", help="Output matrix (.phy)")
    p_mat.add_argument("--method", choices=["ver2", "ver3", "local", "global"], default="ver2")
    p_mat.add_argument("-t", "--threads", type=int, default=4)
    p_mat.add_argument("-s", "--sensitivity", type=float, default=7.5)

    # Subcommand: merge
    p_merge = subparsers.add_parser("merge", help="Merge original and homolog FASTA")
    p_merge.add_argument("ori_fasta", help="Original FASTA")
    p_merge.add_argument("hom_fasta", help="Homolog FASTA")
    p_merge.add_argument("out_fasta", help="Merged output FASTA")

    # Subcommand: clean-tree
    p_clean = subparsers.add_parser("clean-tree", help="Clean negative branch lengths in Newick tree")
    p_clean.add_argument("input_tree", help="Input Newick tree")
    p_clean.add_argument("output_tree", nargs="?", default=None, help="Output Newick tree")

    args = parser.parse_args()

    if args.command == "matrix":
        seq_ids, matrix = compute_pwa_matrix(
            args.input_fasta, method=args.method, threads=args.threads, sensitivity=args.sensitivity
        )
        save_as_phylip(seq_ids, matrix, args.output_matrix)
        print(f"Matrix saved to {args.output_matrix} ({len(seq_ids)} sequences)")

    elif args.command == "merge":
        count = prepare_amplified_fasta(args.ori_fasta, args.hom_fasta, args.out_fasta)
        print(f"Merged {count} sequences into {args.out_fasta}")

    elif args.command == "clean-tree":
        clean_negative_branches(args.input_tree, args.output_tree)
        print(f"Cleaned tree saved to {args.output_tree or args.input_tree}")


if __name__ == "__main__":
    main()

