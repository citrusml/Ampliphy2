#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for src/impoverished_tree.py."""

from pathlib import Path

from src.impoverished_tree import (
    make_impoverished_msa,
    parse_fasta,
    select_enrich_msa,
    subset_phylip_matrix,
)


def test_select_enrich_msa_priority(tmp_path: Path) -> None:
    (tmp_path / "H1.amp.fa").write_text(">a\nAC\n")
    assert select_enrich_msa(tmp_path, "H1") == tmp_path / "H1.amp.fa"
    (tmp_path / "H1.maxaligned.fa").write_text(">a\nAC\n")
    assert select_enrich_msa(tmp_path, "H1") == tmp_path / "H1.maxaligned.fa"
    assert select_enrich_msa(tmp_path, "NONE") is None


def test_make_impoverished_msa_drops_gap_columns(tmp_path: Path) -> None:
    msa = tmp_path / "amp.fa"
    msa.write_text(">s1\nAC-D\n>s2\nA--D\n>h1\nACWD\n")
    ori = tmp_path / "ori.fa"
    ori.write_text(">s1\nACD\n>s2\nAD\n")
    out = tmp_path / "imp.fa"
    n, length = make_impoverished_msa(msa, ori, out)
    assert (n, length) == (2, 3)  # column 3 is all-gap in s1/s2 -> removed
    assert dict(parse_fasta(out)) == {"s1": "ACD", "s2": "A-D"}


def test_subset_phylip_matrix(tmp_path: Path) -> None:
    phy = tmp_path / "amp.phy"
    phy.write_text("3\ns1 0.0 0.1 0.5\ns2 0.1 0.0 0.6\nh1 0.5 0.6 0.0\n")
    out = tmp_path / "imp.phy"
    assert subset_phylip_matrix(phy, {"s1", "s2"}, out) == 2
    assert out.read_text().splitlines() == ["2", "s1 0.0 0.1", "s2 0.1 0.0"]

