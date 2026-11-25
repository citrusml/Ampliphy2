#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_sss_distance.py
====================
Unit tests for SSS distance calculation module (src/sss_distance.py).
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path
from typing import List

import numpy as np
import pytest

from src.sss_distance import (
    clean_negative_branches,
    parse_fasta,
    save_as_phylip,
    similarity_to_distance_matrix,
)


def test_parse_fasta(tmp_path: Path) -> None:
    fa_file = tmp_path / "test.fa"
    fa_file.write_text(
        ">seq1 description\n"
        "M-K-L-V\n"
        ">seq2\n"
        "ACDEFGHIK\n"
        ">seq3\n"
        "-\n"
    )

    ids, seqs = parse_fasta(fa_file)
    assert ids == ["seq1", "seq2", "seq3"]
    assert seqs == ["MKLV", "ACDEFGHIK", ""]


def test_similarity_to_distance_matrix_linear() -> None:
    # 3x3 identity matrix with similarities 1.0, 0.8, 0.4
    w_mat = np.array([
        [1.0, 0.8, 0.4],
        [0.8, 1.0, 0.1],
        [0.4, 0.1, 1.0],
    ], dtype=np.float64)

    d_mat = similarity_to_distance_matrix(w_mat, model="linear")

    assert d_mat.shape == (3, 3)
    assert np.allclose(np.diag(d_mat), 0.0)
    assert np.isclose(d_mat[0, 1], 1.0 - 0.8)
    assert np.isclose(d_mat[0, 2], 1.0 - 0.4)
    assert np.isclose(d_mat[1, 2], 1.0 - 0.1)
    assert np.allclose(d_mat, d_mat.T)


def test_similarity_to_distance_matrix_log() -> None:
    w_mat = np.array([
        [1.0, 0.5, 0.0],
        [0.5, 1.0, 0.2],
        [0.0, 0.2, 1.0],
    ], dtype=np.float64)

    eps = 1e-4
    d_mat = similarity_to_distance_matrix(w_mat, model="log", epsilon=eps)

    assert d_mat.shape == (3, 3)
    assert np.allclose(np.diag(d_mat), 0.0)
    assert np.isclose(d_mat[0, 1], -math.log(0.5))
    assert np.isclose(d_mat[0, 2], -math.log(eps))
    assert np.isclose(d_mat[1, 2], -math.log(0.2))
    assert np.allclose(d_mat, d_mat.T)


def test_save_as_phylip(tmp_path: Path) -> None:
    taxa = ["A", "B", "C"]
    d_mat = np.array([
        [0.0, 0.2, 0.5],
        [0.2, 0.0, 0.3],
        [0.5, 0.3, 0.0],
    ], dtype=np.float64)

    phy_path = tmp_path / "dist.phy"
    save_as_phylip(taxa, d_mat, phy_path)

    assert phy_path.exists()
    lines = phy_path.read_text().strip().split("\n")
    assert lines[0].strip() == "3"
    assert lines[1].startswith("A")
    assert "0.200000" in lines[1]


def test_clean_negative_branches() -> None:
    raw_nwk = "((A:0.05,B:-0.0034):0.12,C:0.45);"
    cleaned = clean_negative_branches(raw_nwk)
    assert ":-0.0034" not in cleaned
    assert ":0.0" in cleaned
    assert "A:0.05" in cleaned
