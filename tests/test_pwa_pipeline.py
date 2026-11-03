#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_pwa_pipeline.py
====================
Unit tests for src/pwa_pipeline.py.
"""

import math
import tempfile
from pathlib import Path
from src.pwa_pipeline import (
    calculate_poisson_distance,
    clean_negative_branches,
    compute_pwa_matrix,
    detect_is_nucleotide,
    generate_lineage_tsv,
    parse_fasta,
    prepare_amplified_fasta,
    save_as_phylip,
    calculate_delta_tcs,
    summarize_pwa_results,
)
import pandas as pd


def test_calculate_poisson_distance():
    # Identical sequences -> 0.0
    seq1 = "ACDEFGHIKLMNPQRSTVWY"
    seq2 = "ACDEFGHIKLMNPQRSTVWY"
    assert calculate_poisson_distance(seq1, seq2, is_nucleotide=False, method="ver2") == 0.0

    # 1 mismatch out of 20
    seq3 = "ACDEFGHIKLMNPQRSTVWX" # X is invalid, effectively 19 matches out of 19
    seq4 = "ACDEFGHIKLMNPQRSTVWA" # 1 mismatch (Y -> A)
    d = calculate_poisson_distance(seq1, seq4, is_nucleotide=False, method="ver2")
    assert 0.0 < d < 0.2

    # Gapped positions are skipped
    seq5 = "ACDEFGHIKL-NPQRSTVWY"
    seq6 = "ACDEFGHIKL-NPQRSTVWY"
    assert calculate_poisson_distance(seq5, seq6, is_nucleotide=False, method="ver2") == 0.0


def test_prepare_amplified_fasta():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        ori = tmp / "ori.fa"
        hom = tmp / "hom.fa"
        out = tmp / "amp.fa"

        ori.write_text(">seq1\nACDEF\n>seq2\nGHIKL\n", encoding="utf-8")
        hom.write_text(">hom1_enr\nMNPQR\n", encoding="utf-8")

        count = prepare_amplified_fasta(ori, hom, out)
        assert count == 3
        records = parse_fasta(out)
        assert len(records) == 3
        assert records[0][0] == "seq1"
        assert records[2][0] == "hom1_enr"

        # Missing homolog fallback
        out2 = tmp / "amp2.fa"
        count2 = prepare_amplified_fasta(ori, tmp / "nonexistent.fa", out2)
        assert count2 == 2
        assert len(parse_fasta(out2)) == 2


def test_clean_negative_branches():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        tree = tmp / "tree.nwk"
        out = tmp / "cleaned.nwk"

        # Newick with single quotes and negative branch length
        tree.write_text("('seq1':-0.0053, 'seq2':0.12):-1.2e-4", encoding="utf-8")
        cleaned = clean_negative_branches(tree, out)

        assert cleaned == "(seq1:0.0, seq2:0.12):0.0;"
        assert out.read_text(encoding="utf-8").strip() == "(seq1:0.0, seq2:0.12):0.0;"


def test_generate_lineage_tsv():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        taxa = tmp / "taxa.tsv"
        taxa.write_text("HOMSA\t9606\tEukaryota, Mammalia, Homo sapiens\n", encoding="utf-8")

        fasta = tmp / "sample.fa"
        fasta.write_text(">HOMSA_P01234\nACDEF\n", encoding="utf-8")

        lineage_out = tmp / "lineage.tsv"
        matched = generate_lineage_tsv(taxa, [fasta], lineage_out)
        assert matched == 1

        lines = lineage_out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "query\tlineage"
        assert "HOMSA_P01234\tEukaryota, Mammalia, Homo sapiens" in lines[1]


def test_calculate_delta_tcs():
    res = calculate_delta_tcs(0.5, 0.75)
    assert math.isclose(res["delta_tcs"], 0.25)
    assert math.isclose(res["ratio_tcs"], 1.5)


def test_summarize_pwa_results():
    df = pd.DataFrame({
        "family_id": ["HOG1", "HOG2", "HOG3"],
        "tcs_ori": [0.5, 0.6, 0.4],
        "tcs_amp": [0.6, 0.6, 0.5],
        "delta_tcs": [0.1, 0.0, 0.1],
        "ratio_tcs": [1.2, 1.0, 1.25],
    })
    summary = summarize_pwa_results(df)
    assert summary["total_hogs"] == 3
    assert summary["improved_count"] == 2
    assert summary["unchanged_count"] == 1
    assert summary["declined_count"] == 0
    assert summary["mean_delta_tcs"] > 0


def test_compute_pwa_matrix_global():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        fa = tmp / "test.fa"
        fa.write_text(
            ">seqA\nACDEFGHIKLMNPQRSTVWY\n"
            ">seqB\nACDEFGHIKLMNPQRSTVWY\n"
            ">seqC\nACDEFGHIKLMNPQRAAAAA\n",
            encoding="utf-8"
        )
        seq_ids, matrix = compute_pwa_matrix(fa, method="ver2")
        assert len(seq_ids) == 3
        assert matrix[0][0] == 0.0
        assert matrix[0][1] == 0.0  # Identical
        assert matrix[0][2] > 0.0   # Diverged

        phy_out = tmp / "test.phy"
        save_as_phylip(seq_ids, matrix, phy_out)
        assert phy_out.exists()
        lines = phy_out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "3"
        assert lines[1].startswith("seqA ")

