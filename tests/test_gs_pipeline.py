import tempfile
from pathlib import Path
from src.gs_pipeline import generate_lineage_tsv


def test_generate_lineage_tsv():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        taxa_file = tmp_path / "taxa.tsv"
        taxa_file.write_text(
            "DANRE\t7955\t131567 (no rank), 2759 (domain), 7955 (leaf)\n"
            "MOUSE\t10090\t131567 (no rank), 2759 (domain), 10090 (leaf)\n",
            encoding="utf-8",
        )

        fasta1 = tmp_path / "f1.fa"
        fasta1.write_text(
            ">DANRE001 foo\nATGC\n>MOUSE002 bar\nATGC\n", encoding="utf-8"
        )

        out_tsv = tmp_path / "lineage.tsv"
        count = generate_lineage_tsv(taxa_file, [fasta1], out_tsv)

        assert count == 2
        content = out_tsv.read_text(encoding="utf-8").splitlines()
        assert content[0] == "query\tlineage"
        assert "DANRE001\t131567 (no rank), 2759 (domain), 7955 (leaf)" in content
        assert "MOUSE002\t131567 (no rank), 2759 (domain), 10090 (leaf)" in content


def test_clean_newick_tree():
    """Verify that clean_newick_tree removes single and double quotes and ensures a trailing semicolon."""
    from src.gs_pipeline import clean_newick_tree

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        tree_file = tmp_path / "tree.nwk"
        tree_file.write_text("('A':0.1, \"B\":0.2):0.05", encoding="utf-8")

        cleaned = clean_newick_tree(tree_file)
        assert cleaned == "(A:0.1, B:0.2):0.05;"
        assert tree_file.read_text(encoding="utf-8").strip() == "(A:0.1, B:0.2):0.05;"


def test_extract_homolog_ids():
    """Verify that extract_homolog_ids correctly identifies sequences added to the amplified FASTA."""
    from src.gs_pipeline import extract_homolog_ids

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        orig_fa = tmp_path / "orig.fa"
        orig_fa.write_text(">seq1 description1\nATGC\n>seq2 description2\nCGTA\n", encoding="utf-8")

        amp_fa = tmp_path / "amp.fa"
        amp_fa.write_text(
            ">seq1 description1\nATGC\n>homolog1 desc\nATGG\n>seq2 description2\nCGTA\n>homolog2\nAAAA\n",
            encoding="utf-8",
        )

        homologs = extract_homolog_ids(orig_fa, amp_fa)
        assert homologs == ["homolog1", "homolog2"]


def test_parse_and_summarize_gs_tcs_scores():
    """Verify parsing of TCS score output and computation of summary metrics."""
    from src.gs_pipeline import parse_gs_tcs_scores, summarize_gs_results

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        score_file = tmp_path / "scores.tsv"
        score_file.write_text(
            "/path/to/eval/aug/HOG_E001_gs_aug.nwk\t0.500\n"
            "/path/to/eval/enr/HOG_E001_gs_enr.nwk\t0.600\n"
            "/path/to/eval/aug/HOG_E002_gs_aug.nwk\t0.400\n"
            "/path/to/eval/enr/HOG_E002_gs_enr.nwk\t0.400\n",
            encoding="utf-8",
        )

        df = parse_gs_tcs_scores(score_file)
        assert len(df) == 2
        assert "HOG_E001" in df["hog_id"].values
        assert "HOG_E002" in df["hog_id"].values

        summary = summarize_gs_results(df)
        assert summary["n_hogs"] == 2
        assert summary["mean_tcs_ori"] == 0.45
        assert summary["mean_tcs_enr"] == 0.50
        assert summary["pct_improved"] == 50.0


