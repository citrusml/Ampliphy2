#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
msa_distance.py
===============
Core module for computing pairwise Poisson evolutionary distance matrices from
a Multiple Sequence Alignment (MSA) and exporting them in PHYLIP format,
plus Newick branch-length cleaning utilities.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np


def parse_msa_fasta(fasta_path: Union[str, Path]) -> Tuple[List[str], List[str]]:
    """
    Parse an aligned FASTA file and return sequence IDs and aligned sequences.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to the input MSA FASTA file.

    Returns
    -------
    ids : List[str]
        Ordered list of sequence identifiers (without leading '>').
    seqs : List[str]
        Ordered list of uppercase aligned sequences preserving gap characters.
    """
    ids: List[str] = []
    seqs: List[str] = []
    curr_id: Optional[str] = None
    curr_seq: List[str] = []

    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_id is not None:
                    ids.append(curr_id)
                    seqs.append("".join(curr_seq).upper())
                curr_id = line[1:].split()[0]
                curr_seq = []
            else:
                curr_seq.append(line)

        if curr_id is not None:
            ids.append(curr_id)
            seqs.append("".join(curr_seq).upper())

    return ids, seqs


def compute_msa_poisson_matrix(
    fasta_path: Union[str, Path],
    epsilon: float = 1e-4,
) -> Tuple[List[str], np.ndarray]:
    """
    Compute a pairwise Poisson evolutionary distance matrix from an MSA using
    pairwise gap deletion.

    For each pair of sequences (i, j), non-gap positions in both sequences are
    counted to obtain proportion of differing sites p_ij = diff / valid_len.
    The Poisson distance is d_ij = -ln(max(1.0 - p_ij, epsilon)).

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to input aligned FASTA file.
    epsilon : float, default=1e-4
        Minimum floor for (1 - p_ij) to prevent log(0). Physical unit: unitless.

    Returns
    -------
    taxa_names : List[str]
        Ordered list of taxon identifiers.
    d_mat : np.ndarray
        (N, N) symmetric pairwise Poisson distance matrix with zero diagonal.
        Physical unit: substitutions per site (unitless).
    """
    taxa_names, seqs = parse_msa_fasta(fasta_path)
    if not taxa_names:
        raise ValueError(f"No sequences found in {fasta_path}")

    n = len(taxa_names)
    d_mat = np.zeros((n, n), dtype=np.float64)
    max_dist = -math.log(epsilon)
    gap_chars = {"-", ".", "?", "*"}

    # Encode sequences as byte arrays for fast vectorized comparison
    max_len = max(len(s) for s in seqs)
    padded = [s.ljust(max_len, "-") for s in seqs]
    arr = np.frombuffer("".join(padded).encode("ascii", errors="replace"), dtype=np.uint8).reshape(n, max_len)

    is_gap = np.zeros((n, max_len), dtype=bool)
    for g in gap_chars:
        is_gap |= arr == ord(g)

    for i in range(n):
        valid_i = ~is_gap[i]
        row_i = arr[i]
        for j in range(i + 1, n):
            both_valid = valid_i & (~is_gap[j])
            valid_len = int(np.count_nonzero(both_valid))
            if valid_len == 0:
                d_val = max_dist
            else:
                matches = int(np.count_nonzero((row_i == arr[j]) & both_valid))
                identity = matches / valid_len
                d_val = max(0.0, -math.log(max(identity, epsilon)))
            d_mat[i, j] = d_val
            d_mat[j, i] = d_val

    return taxa_names, d_mat


def save_as_phylip(
    taxa_names: List[str],
    dist_matrix: np.ndarray,
    out_path: Union[str, Path],
) -> None:
    """
    Save a distance matrix in square PHYLIP format.

    Parameters
    ----------
    taxa_names : List[str]
        Ordered list of taxon identifiers.
    dist_matrix : np.ndarray
        (N, N) symmetric distance matrix. Physical unit: substitutions per site.
    out_path : Union[str, Path]
        Destination PHYLIP file path.
    """
    n = len(taxa_names)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"    {n}\n")
        for i, name in enumerate(taxa_names):
            vals = " ".join(f"{dist_matrix[i, j]:.6f}" for j in range(n))
            f.write(f"{name:<10} {vals}\n")


def clean_negative_branches(tree_newick: str) -> str:
    """
    Replace negative branch lengths in a Newick string with ':0.0'.

    Parameters
    ----------
    tree_newick : str
        Input Newick tree string.

    Returns
    -------
    cleaned_tree : str
        Newick string with all negative branch lengths replaced by ':0.0'.
    """
    pattern = r":-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?"

    def repl(m: re.Match[str]) -> str:
        s = m.group(0)
        try:
            if float(s[1:]) < 0.0:
                return ":0.0"
        except ValueError:
            pass
        return s

    return re.sub(pattern, repl, tree_newick)


def clean_tree_file(tree_path: Union[str, Path]) -> None:
    """
    Clean negative branch lengths and quotation marks in a Newick file in-place.

    Parameters
    ----------
    tree_path : Union[str, Path]
        Path to the Newick tree file.
    """
    with open(tree_path, "r", encoding="utf-8") as f:
        content = f.read().strip().replace('"', "").replace("'", "")
    cleaned = clean_negative_branches(content)
    with open(tree_path, "w", encoding="utf-8") as f:
        f.write(cleaned + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="MSA Poisson distance matrix utilities")
    sub = parser.add_subparsers(dest="command", required=True)

    p_mat = sub.add_parser("matrix", help="Compute Poisson distance matrix from MSA")
    p_mat.add_argument("msa", help="Input aligned FASTA file")
    p_mat.add_argument("output", help="Output PHYLIP distance matrix file")
    p_mat.add_argument("--epsilon", type=float, default=1e-4, help="Floor for identity (default: 1e-4)")

    p_cln = sub.add_parser("clean-tree", help="Clean quotes and negative branch lengths in Newick file")
    p_cln.add_argument("tree", help="Path to Newick tree file")

    args = parser.parse_args()
    if args.command == "matrix":
        taxa, d_mat = compute_msa_poisson_matrix(args.msa, epsilon=args.epsilon)
        save_as_phylip(taxa, d_mat, args.output)
        print(f"Saved Poisson distance matrix ({len(taxa)} taxa) -> {args.output}")
    elif args.command == "clean-tree":
        clean_tree_file(args.tree)


if __name__ == "__main__":
    main()
