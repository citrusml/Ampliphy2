#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_msa_distance.py
====================
Unit tests for src/msa_distance.py and src/tcsr2.py.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from src.msa_distance import (
    clean_negative_branches,
    compute_msa_poisson_matrix,
    parse_msa_fasta,
    save_as_phylip,
)
from src.tcsr2 import compute_tree_tcsr2


def test_parse_msa_fasta(tmp_path: Path) -> None:
    fa = tmp_path / "test.msa.fa"
    fa.write_text(">seq1 desc\nA-CD\nEF\n>seq2\nABCD\nEF\n")
    ids, seqs = parse_msa_fasta(fa)
    assert ids == ["seq1", "seq2"]
    assert seqs == ["A-CDEF", "ABCDEF"]


def test_compute_msa_poisson_matrix_pairwise_deletion(tmp_path: Path) -> None:
    fa = tmp_path / "test.msa.fa"
    # seq1 vs seq2: pos 0 (A,A), pos 2 (C,C), pos 3 (D,X) -> 3 valid sites, 2 matches -> identity=2/3
    # seq1 vs seq3: pos 0..3 -> 3 valid sites, 3 matches -> d=0
    fa.write_text(">seq1\nA-CD\n>seq2\nABCX\n>seq3\nA-CD\n")
    taxa, d_mat = compute_msa_poisson_matrix(fa)

    assert taxa == ["seq1", "seq2", "seq3"]
    assert d_mat.shape == (3, 3)
    assert np.allclose(np.diag(d_mat), 0.0)
    assert np.isclose(d_mat[0, 1], -math.log(2.0 / 3.0))
    assert np.isclose(d_mat[0, 2], 0.0)
    assert np.allclose(d_mat, d_mat.T)


def test_save_as_phylip(tmp_path: Path) -> None:
    taxa = ["T1", "T2"]
    d_mat = np.array([[0.0, 0.25], [0.25, 0.0]], dtype=np.float64)
    out = tmp_path / "dist.phy"
    save_as_phylip(taxa, d_mat, out)
    lines = out.read_text().strip().splitlines()
    assert lines[0].strip() == "2"
    assert "0.250000" in lines[1]


def test_clean_negative_branches() -> None:
    assert clean_negative_branches("((A:0.1,B:-0.05):0.2,C:0.3);") == "((A:0.1,B:0.0):0.2,C:0.3);"


def test_compute_tree_tcsr2() -> None:
    taxa_map = {
        "SPAAA": {"Mammalia", "Vertebrata"},
        "SPBBB": {"Mammalia", "Vertebrata"},
        "SPCCC": {"Actinopterygii", "Vertebrata"},
        "SPDDD": {"Actinopterygii", "Vertebrata"},
    }
    # Perfect split: (SPAAA, SPBBB) vs (SPCCC, SPDDD)
    nwk_perfect = "((SPAAA_1:0.1,SPBBB_1:0.1):0.2,(SPCCC_1:0.1,SPDDD_1:0.1):0.2);"
    assert np.isclose(compute_tree_tcsr2(nwk_perfect, taxa_map), 1.0)

    # Split mammals: (SPAAA, SPCCC) vs (SPBBB, SPDDD) -> penalty > 0
    nwk_bad = "((SPAAA_1:0.1,SPCCC_1:0.1):0.2,(SPBBB_1:0.1,SPDDD_1:0.1):0.2);"
    score_bad = compute_tree_tcsr2(nwk_bad, taxa_map)
    assert 0.0 <= score_bad < 1.0
