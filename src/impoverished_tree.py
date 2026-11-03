#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
impoverished_tree.py
====================
Utilities for building Impoverished trees.

Impoverished tree = tree inferred from the amplified MSA restricted to the
original (HOG member) sequences, keeping the alignment column structure
(all-gap columns removed).

Features:
- Select the MSA actually used for the enrich tree (maxaligned > amp)
- Extract original sequences from the amplified MSA
- Subset a PHYLIP distance matrix to original sequences (for results_pwa)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple, Union

PathLike = Union[str, Path]

# Candidate suffixes in priority order (MaxAlign output is what the enrich tree used)
MSA_CANDIDATE_SUFFIXES: Tuple[str, ...] = (
    ".maxaligned.fa",
    ".after.maxaligned.fa",
    ".amp.fa",
)


def parse_fasta(fasta_path: PathLike) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file into (id, sequence) tuples.

    Parameters
    ----------
    fasta_path : PathLike
        Path to FASTA file.

    Returns
    -------
    records : List[Tuple[str, str]]
        List of (seq_id, sequence). The id is the first whitespace-delimited
        token of the header. Sequence length unit: residues (or alignment columns).
    """
    records: List[Tuple[str, str]] = []
    current_id: Optional[str] = None
    chunks: List[str] = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks)))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, "".join(chunks)))
    return records


def select_enrich_msa(msa_dir: PathLike, hog_id: str) -> Optional[Path]:
    """
    Select the MSA that was used to infer the enrich tree.

    Parameters
    ----------
    msa_dir : PathLike
        Directory containing ``<hog>.amp.fa`` / ``<hog>.maxaligned.fa``.
    hog_id : str
        HOG identifier (unitless), e.g. ``HOG_E0736661``.

    Returns
    -------
    path : Optional[Path]
        Selected MSA path (priority: maxaligned > after.maxaligned > amp), or None.
    """
    base = Path(msa_dir)
    for suffix in MSA_CANDIDATE_SUFFIXES:
        cand = base / f"{hog_id}{suffix}"
        if cand.is_file() and cand.stat().st_size > 0:
            return cand
    return None


def make_impoverished_msa(
    msa_path: PathLike,
    ori_fasta: PathLike,
    out_fasta: PathLike,
) -> Tuple[int, int]:
    """
    Extract original sequences from an amplified MSA and drop all-gap columns.

    Parameters
    ----------
    msa_path : PathLike
        Amplified MSA (aligned FASTA).
    ori_fasta : PathLike
        Original HOG FASTA providing the IDs to keep.
    out_fasta : PathLike
        Output impoverished MSA path.

    Returns
    -------
    n_seqs, aln_len : Tuple[int, int]
        Number of kept sequences (count) and alignment length (columns).
    """
    keep_ids: Set[str] = {r[0] for r in parse_fasta(ori_fasta)}
    kept = [(i, s) for i, s in parse_fasta(msa_path) if i in keep_ids]
    if not kept:
        return 0, 0

    aln_len = max(len(s) for _, s in kept)
    padded = [(i, s.ljust(aln_len, "-")) for i, s in kept]
    keep_cols = [
        c for c in range(aln_len) if any(s[c] not in ("-", ".") for _, s in padded)
    ]
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w", encoding="utf-8") as f:
        for seq_id, seq in padded:
            f.write(f">{seq_id}\n{''.join(seq[c] for c in keep_cols)}\n")
    return len(padded), len(keep_cols)


def subset_phylip_matrix(
    phy_path: PathLike,
    keep_ids: Set[str],
    out_phy: PathLike,
) -> int:
    """
    Subset a square PHYLIP distance matrix to the given sequence IDs.

    Parameters
    ----------
    phy_path : PathLike
        Input square PHYLIP matrix (distance unit: substitutions per site).
    keep_ids : Set[str]
        IDs to keep.
    out_phy : PathLike
        Output PHYLIP matrix path.

    Returns
    -------
    n_kept : int
        Number of sequences kept (count).
    """
    with open(phy_path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    names: List[str] = []
    rows: List[List[str]] = []
    for ln in lines[1:]:
        parts = ln.split()
        names.append(parts[0])
        rows.append(parts[1:])

    idx = [k for k, n in enumerate(names) if n in keep_ids]
    with open(out_phy, "w", encoding="utf-8") as f:
        f.write(f"{len(idx)}\n")
        for i in idx:
            f.write(f"{names[i]} {' '.join(rows[i][j] for j in idx)}\n")
    return len(idx)


def main() -> None:
    parser = argparse.ArgumentParser(description="Impoverished tree utilities")
    sub = parser.add_subparsers(dest="command", required=True)

    p_sel = sub.add_parser("select", help="Print the enrich MSA path for a HOG")
    p_sel.add_argument("msa_dir")
    p_sel.add_argument("hog_id")

    p_msa = sub.add_parser("msa", help="Build impoverished MSA")
    p_msa.add_argument("msa_pos", nargs="?", default=None, help="Amplified MSA path")
    p_msa.add_argument("ori_pos", nargs="?", default=None, help="Original MSA path")
    p_msa.add_argument("out_pos", nargs="?", default=None, help="Output path")
    p_msa.add_argument("--msa", default=None, help="Amplified MSA path (flag)")
    p_msa.add_argument("--ori", default=None, help="Original MSA path (flag)")
    p_msa.add_argument("--out", default=None, help="Output path (flag)")

    p_mat = sub.add_parser("matrix", help="Subset PHYLIP matrix to original IDs")
    p_mat.add_argument("amp_phy")
    p_mat.add_argument("ori_fasta")
    p_mat.add_argument("out_phy")

    args = parser.parse_args()
    if args.command == "select":
        path = select_enrich_msa(args.msa_dir, args.hog_id)
        if path is None:
            sys.exit(1)
        print(path)
    elif args.command == "msa":
        msa_path = args.msa or args.msa_pos
        ori_path = args.ori or args.ori_pos
        out_path = args.out or args.out_pos
        if not (msa_path and ori_path and out_path):
            p_msa.error("msa requires msa, ori, and out paths (either positional or as flags)")
        n, length = make_impoverished_msa(msa_path, ori_path, out_path)
        if n == 0:
            print("No original sequences found in MSA", file=sys.stderr)
            sys.exit(1)
        print(f"{n} seqs x {length} cols -> {out_path}")
    elif args.command == "matrix":
        ids = {r[0] for r in parse_fasta(args.ori_fasta)}
        n = subset_phylip_matrix(args.amp_phy, ids, args.out_phy)
        if n < 2:
            print("Fewer than 2 original sequences in matrix", file=sys.stderr)
            sys.exit(1)
        print(f"{n} seqs -> {args.out_phy}")


if __name__ == "__main__":
    main()

