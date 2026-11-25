#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sss_distance.py
===============
Core module for SSS (Symmetrical Sequence Similarity score) phylogenetic distance
calculation and phylogenetic tree processing.

Implements the symmetrical relative sequence similarity score w_ij
(Matsui & Iwasaki 2020, Systematic Biology, syz049; Dufour et al. 2010):
    w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))

Distance models:
- sss_linear: d_ij = 1.0 - w_ij
- sss_log:    d_ij = -ln(max(w_ij, epsilon)), epsilon = 1e-4
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np


def find_mmseqs_binary() -> str:
    """
    Find MMseqs2 executable in PATH or standard installation locations.

    Returns
    -------
    bin_path : str
        Path to the mmseqs executable.
    """
    for cand in [
        shutil.which("mmseqs"),
        os.path.expanduser("~/.local/bin/mmseqs"),
        os.path.expanduser("~/micromamba/bin/mmseqs"),
        "/usr/local/bin/mmseqs",
    ]:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return str(cand)
    return "mmseqs"


def parse_fasta(fasta_path: Union[str, Path]) -> Tuple[List[str], List[str]]:
    """
    Parses a FASTA file and returns identifiers and uppercase gapless sequences.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to input FASTA file.

    Returns
    -------
    ids : List[str]
        List of sequence identifiers (without leading '>').
    seqs : List[str]
        List of uppercase sequences with gaps removed.
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
                    seqs.append("".join(curr_seq).replace("-", "").upper())
                curr_id = line[1:].split()[0]
                curr_seq = []
            else:
                curr_seq.append(line)

        if curr_id is not None:
            ids.append(curr_id)
            seqs.append("".join(curr_seq).replace("-", "").upper())

    return ids, seqs


def compute_mmseqs_sss_matrix(
    fasta_path: Union[str, Path],
    sensitivity: float = 7.5,
    evalue: float = 10.0,
    threads: int = 4,
    tmp_dir: Optional[Union[str, Path]] = None,
) -> Tuple[List[str], np.ndarray]:
    """
    Computes all-to-all symmetrical sequence similarity matrix w_ij.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to input FASTA file.
    sensitivity : float, default=7.5
        MMseqs2 sensitivity parameter (-s). Physical unit: unitless.
    evalue : float, default=10.0
        MMseqs2 E-value threshold (-e). Physical unit: unitless.
    threads : int, default=4
        Number of CPU threads. Physical unit: count.
    tmp_dir : Optional[Union[str, Path]], default=None
        Directory for intermediate MMseqs2 databases.

    Returns
    -------
    taxa_names : List[str]
        Ordered list of sequence headers.
    w_mat : np.ndarray
        (N, N) symmetric similarity matrix with values in [0.0, 1.0].
        Physical unit: unitless.
    """
    taxa_names, clean_seqs = parse_fasta(fasta_path)
    if not taxa_names:
        raise ValueError(f"No sequences found in {fasta_path}")

    n = len(taxa_names)
    id_map = {name: i for i, name in enumerate(taxa_names)}
    mmseqs_bin = find_mmseqs_binary()
    s_bits = np.zeros((n, n), dtype=np.float64)

    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="mmseqs_sss_") as td:
        clean_fa = os.path.join(td, "clean.fa")
        with open(clean_fa, "w", encoding="utf-8") as wf:
            for name, seq in zip(taxa_names, clean_seqs):
                wf.write(f">{name}\n{seq}\n")

        db = os.path.join(td, "db")
        res = os.path.join(td, "res")
        sub_tmp = os.path.join(td, "tmp")
        out_m8 = os.path.join(td, "out.m8")
        os.makedirs(sub_tmp, exist_ok=True)

        # 1. MMseqs2 createdb
        subprocess.run(
            [mmseqs_bin, "createdb", clean_fa, db],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )

        # 2. MMseqs2 all-to-all search
        subprocess.run(
            [
                mmseqs_bin,
                "search",
                db,
                db,
                res,
                sub_tmp,
                "--threads",
                str(threads),
                "-s",
                str(sensitivity),
                "-e",
                str(evalue),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )

        # 3. MMseqs2 convertalis
        subprocess.run(
            [mmseqs_bin, "convertalis", db, db, res, out_m8, "--format-mode", "0"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )

        if os.path.exists(out_m8):
            with open(out_m8, "r", encoding="utf-8") as rf:
                for line in rf:
                    parts = line.strip().split("\t")
                    if len(parts) >= 12:
                        q, t = parts[0], parts[1]
                        if q in id_map and t in id_map:
                            try:
                                bits = float(parts[11])
                            except ValueError:
                                bits = 0.0
                            qi, ti = id_map[q], id_map[t]
                            if bits > s_bits[qi, ti]:
                                s_bits[qi, ti] = bits

    # Compute w_ij
    w_mat = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        w_mat[i, i] = 1.0
        for j in range(i + 1, n):
            self_score = s_bits[i, i] + s_bits[j, j]
            comp_score = s_bits[i, j] + s_bits[j, i]
            val = (comp_score / self_score) if self_score > 0.0 else 0.0
            val = min(1.0, max(0.0, float(val)))
            w_mat[i, j] = val
            w_mat[j, i] = val

    return taxa_names, w_mat


def similarity_to_distance_matrix(
    w_mat: np.ndarray,
    model: str = "linear",
    epsilon: float = 1e-4,
) -> np.ndarray:
    """
    Transforms similarity matrix w_ij into distance matrix d_ij.

    Parameters
    ----------
    w_mat : np.ndarray
        (N, N) symmetric similarity matrix with values in [0.0, 1.0].
        Physical unit: unitless.
    model : str, default='linear'
        Distance model: 'linear' (d_ij = 1 - w_ij) or 'log' (d_ij = -ln(max(w_ij, eps))).
    epsilon : float, default=1e-4
        Floor threshold for log transformation. Physical unit: unitless.

    Returns
    -------
    d_mat : np.ndarray
        (N, N) distance matrix with 0.0 on diagonal.
        Physical unit: unitless.
    """
    n = w_mat.shape[0]
    d_mat = np.zeros((n, n), dtype=np.float64)
    is_linear = model.lower() in ("linear", "sss_linear", "lin")

    for i in range(n):
        for j in range(i + 1, n):
            w = float(w_mat[i, j])
            if is_linear:
                d = max(0.0, 1.0 - w)
            else:
                eff_w = max(w, epsilon)
                d = max(0.0, -math.log(eff_w))
            d_mat[i, j] = d
            d_mat[j, i] = d

    return d_mat


def save_as_phylip(
    taxa_names: List[str],
    dist_matrix: np.ndarray,
    out_path: Union[str, Path],
) -> None:
    """
    Saves pairwise distance matrix in square PHYLIP format.

    Parameters
    ----------
    taxa_names : List[str]
        List of taxon names.
    dist_matrix : np.ndarray
        (N, N) distance matrix.
    out_path : Union[str, Path]
        Path to output PHYLIP file.
    """
    n = len(taxa_names)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"    {n}\n")
        for i, name in enumerate(taxa_names):
            row_vals = " ".join(f"{dist_matrix[i, j]:.6f}" for j in range(n))
            f.write(f"{name:<10} {row_vals}\n")


def clean_negative_branches(tree_newick: str) -> str:
    """
    Replaces negative branch lengths in Newick string with :0.0.

    Parameters
    ----------
    tree_newick : str
        Newick tree string.

    Returns
    -------
    cleaned_tree : str
        Newick string with non-negative branch lengths.
    """
    pattern = r":-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?"

    def repl(m: re.Match[str]) -> str:
        s = m.group(0)
        try:
            val = float(s[1:])
            if val < 0:
                return ":0.0"
            return s
        except ValueError:
            return s

    return re.sub(pattern, repl, tree_newick)


def clean_tree_file(tree_path: Union[str, Path]) -> None:
    """
    Reads a tree file, replaces negative branch lengths, and writes it back.

    Parameters
    ----------
    tree_path : Union[str, Path]
        Path to Newick tree file.
    """
    with open(tree_path, "r", encoding="utf-8") as f:
        content = f.read().strip()
    cleaned = clean_negative_branches(content)
    with open(tree_path, "w", encoding="utf-8") as f:
        f.write(cleaned + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="SSS distance matrix utilities")
    sub = parser.add_subparsers(dest="command", required=True)

    # Subcommand: matrix
    p_mat = sub.add_parser("matrix", help="Compute SSS distance matrix")
    p_mat.add_argument("fasta", help="Input FASTA file")
    p_mat.add_argument("output", help="Output PHYLIP distance matrix file")
    p_mat.add_argument("--model", choices=["linear", "log"], default="linear", help="SSS distance model")
    p_mat.add_argument("--threads", type=int, default=4, help="MMseqs2 threads")
    p_mat.add_argument("--sensitivity", type=float, default=7.5, help="MMseqs2 sensitivity (-s)")
    p_mat.add_argument("--evalue", type=float, default=10.0, help="MMseqs2 evalue (-e)")

    # Subcommand: clean-tree
    p_cln = sub.add_parser("clean-tree", help="Replace negative branch lengths with :0.0")
    p_cln.add_argument("tree", help="Path to Newick tree file")

    args = parser.parse_args()

    if args.command == "matrix":
        taxa, w_mat = compute_mmseqs_sss_matrix(
            args.fasta,
            sensitivity=args.sensitivity,
            evalue=args.evalue,
            threads=args.threads,
        )
        d_mat = similarity_to_distance_matrix(w_mat, model=args.model)
        save_as_phylip(taxa, d_mat, args.output)
        print(f"Saved {args.model} SSS distance matrix ({len(taxa)} taxa) -> {args.output}")

    elif args.command == "clean-tree":
        clean_tree_file(args.tree)


if __name__ == "__main__":
    main()
