import json
import math
from pathlib import Path

import pytest

from ecomp.cli import main as ecomp_main
from ecomp.io import read_alignment
from ecomp.storage import read_metadata

FASTA_BODY = ">s1\nACGTACGT\n>s2\nACGTTCGT\n"
TREE_BODY = "(s1:0.1,s2:0.2);\n"


@pytest.mark.integration
def test_cli_compress_and_decompress_alignment(tmp_path: Path) -> None:
    alignment = tmp_path / "example.fasta"
    alignment.write_text(FASTA_BODY)

    archive = tmp_path / "example.ecomp"
    restored = tmp_path / "restored.fasta"
    metadata = tmp_path / "example.json"

    assert ecomp_main([
        "zip",
        str(alignment),
        "-o",
        str(archive),
        "-m",
        str(metadata),
        "--stats",
    ]) == 0

    assert archive.exists()
    assert metadata.exists()

    assert ecomp_main([
        "unzip",
        str(archive),
        "-m",
        str(metadata),
        "-o",
        str(restored),
    ]) == 0

    original = read_alignment(alignment)
    round_tripped = read_alignment(restored)
    assert original.sequences == round_tripped.sequences
    assert original.ids == round_tripped.ids


@pytest.mark.integration
def test_cli_compress_with_tree_for_ordering(tmp_path: Path) -> None:
    alignment = tmp_path / "with_tree.fasta"
    alignment.write_text(FASTA_BODY)
    tree = tmp_path / "with_tree.tree"
    tree.write_text(TREE_BODY)

    archive = tmp_path / "with_tree.ecomp"
    metadata = tmp_path / "with_tree.json"

    assert ecomp_main([
        "zip",
        str(alignment),
        "--tree",
        str(tree),
        "-o",
        str(archive),
        "-m",
        str(metadata),
    ]) == 0

    restored = tmp_path / "with_tree.restored.fasta"
    assert ecomp_main([
        "unzip",
        str(archive),
        "-m",
        str(metadata),
        "-o",
        str(restored),
    ]) == 0

    original = read_alignment(alignment)
    round_tripped = read_alignment(restored)
    assert original.sequences == round_tripped.sequences
    assert original.ids == round_tripped.ids

    stored_metadata = read_metadata(metadata)
    assert "tree_newick" not in stored_metadata


@pytest.mark.integration
def test_cli_inspect_summary(tmp_path: Path, capsys) -> None:
    alignment = tmp_path / "inspect.fasta"
    alignment.write_text(FASTA_BODY)
    archive = tmp_path / "inspect.ecomp"

    assert ecomp_main([
        "zip",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0

    capsys.readouterr()

    capsys.readouterr()  # clear zip output before metrics

    assert ecomp_main([
        "inspect",
        str(archive),
        "--summary",
    ]) == 0

    summary = capsys.readouterr().out.strip().splitlines()
    assert any(line.startswith("Codec: ecomp") for line in summary)
    assert any("Sequences" in line for line in summary)
    assert any("Alignment columns" in line for line in summary)


@pytest.mark.integration
def test_cli_inspect_raw_json(tmp_path: Path, capsys) -> None:
    alignment = tmp_path / "inspect_raw.fasta"
    alignment.write_text(FASTA_BODY)
    archive = tmp_path / "inspect_raw.ecomp"

    assert ecomp_main([
        "zip",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0

    capsys.readouterr()  # clear zip output

    assert ecomp_main([
        "inspect",
        str(archive),
    ]) == 0

    metadata = json.loads(capsys.readouterr().out)
    assert metadata["codec"] == "ecomp"
    assert metadata["num_sequences"] == 2


@pytest.mark.integration
def test_cli_metrics_commands(tmp_path: Path, capsys) -> None:
    alignment = tmp_path / "metrics.fasta"
    alignment.write_text(
        ">a\nACGTT\n>b\nA-GTT\n>c\nATGCT\n>d\nATGCT\n"
    )
    archive = tmp_path / "metrics.ecomp"

    assert ecomp_main([
        "zip",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0

    capsys.readouterr()

    assert ecomp_main([
        "consensus_sequence",
        str(archive),
    ]) == 0
    output = capsys.readouterr().out.strip().splitlines()
    assert output == [">consensus", "ATGYT"]

    assert ecomp_main([
        "con_seq",
        str(archive),
        "--header",
        "alias",
        "--no-checksum",
    ]) == 0
    alias_output = capsys.readouterr().out.strip().splitlines()
    assert alias_output == [">alias", "ATGYT"]

    assert ecomp_main([
        "column_base_counts",
        str(archive),
        "--indent",
        "0",
    ]) == 0
    counts = json.loads(capsys.readouterr().out)
    assert counts[0]["counts"] == {"A": 4}
    assert counts[1]["counts"] == {"C": 1, "T": 2}

    assert ecomp_main([
        "col_counts",
        str(archive),
        "--include-gaps",
        "--indent",
        "0",
    ]) == 0
    counts_with_gaps = json.loads(capsys.readouterr().out)
    assert counts_with_gaps[1]["counts"] == {"C": 1, "T": 2, "-": 1}

    assert ecomp_main([
        "gap_fraction",
        str(archive),
    ]) == 0
    gap_lines = capsys.readouterr().out.strip().splitlines()
    frac_values = {int(idx): float(value) for idx, value in (line.split("\t") for line in gap_lines)}
    assert frac_values[1] == 0.0
    assert math.isclose(frac_values[2], 0.25, rel_tol=1e-9)

    assert ecomp_main([
        "shannon_entropy",
        str(archive),
    ]) == 0
    entropy_lines = capsys.readouterr().out.strip().splitlines()
    entropy_values = {int(idx): float(value) for idx, value in (line.split("\t") for line in entropy_lines)}
    assert entropy_values[1] == 0.0
    assert math.isclose(entropy_values[2], 0.918295, rel_tol=1e-5)
    assert math.isclose(entropy_values[4], 1.0, rel_tol=1e-6)

    assert ecomp_main([
        "parsimony_informative_sites",
        str(archive),
    ]) == 0
    pis_lines = capsys.readouterr().out.strip().splitlines()
    assert pis_lines[0] == "total\t1"
    assert pis_lines[1] == "indices\t4"

    assert ecomp_main([
        "constant_columns",
        str(archive),
    ]) == 0
    const_lines = capsys.readouterr().out.strip().splitlines()
    assert const_lines[0] == "total\t3"
    assert const_lines[1] == "indices\t1 3 5"

    assert ecomp_main([
        "pairwise_identity",
        str(archive),
    ]) == 0
    pid_output = capsys.readouterr().out.strip().splitlines()
    assert pid_output[0].startswith("id\t")
    identity_header = pid_output[0].split("\t")
    assert identity_header[1:] == ["a", "b", "c", "d"]
    first_row = pid_output[1].split("\t")
    assert first_row[0] == "a"
    assert first_row[1] == "1.000000"
    assert any(line.startswith("# coverage") for line in pid_output)

    assert ecomp_main([
        "pid",
        str(archive),
        "--no-checksum",
    ]) == 0
    pid_alias_output = capsys.readouterr().out.strip().splitlines()
    assert pid_alias_output[0] == pid_output[0]

    variable_alignment = tmp_path / "variable.fasta"
    variable_alignment.write_text(">x\nAC\n>y\nCA\n")
    variable_archive = variable_alignment.with_suffix(".ecomp")
    assert ecomp_main([
        "zip",
        str(variable_alignment),
    ]) == 0
    capsys.readouterr()

    assert ecomp_main([
        "parsimony_informative_sites",
        str(variable_archive),
    ]) == 0
    pis_zero = capsys.readouterr().out.strip().splitlines()
    assert pis_zero == ["total\t0"]

    assert ecomp_main([
        "constant_columns",
        str(variable_archive),
    ]) == 0
    const_zero = capsys.readouterr().out.strip().splitlines()
    assert const_zero == ["total\t0"]

    assert ecomp_main([
        "alignment_length_excluding_gaps",
        str(archive),
    ]) == 0
    len_no_gaps = capsys.readouterr().out.strip()
    assert len_no_gaps == "5"

    assert ecomp_main([
        "alignment_compressed_length",
        str(archive),
    ]) == 0
    compressed_len = capsys.readouterr().out.strip()
    assert compressed_len == "2"

    assert ecomp_main([
        "variable_sites",
        str(archive),
    ]) == 0
    var_sites = capsys.readouterr().out.strip()
    assert var_sites == "2"

    assert ecomp_main([
        "percentage_identity",
        str(archive),
    ]) == 0
    pct_identity = float(capsys.readouterr().out.strip())
    assert 50.0 < pct_identity < 100.0

    assert ecomp_main([
        "relative_composition_variability",
        str(archive),
    ]) == 0
    rcv_value = float(capsys.readouterr().out.strip())
    assert rcv_value >= 0.0

    single_alignment = tmp_path / "single.fasta"
    single_alignment.write_text(">solo\nAAAA\n")
    single_archive = single_alignment.with_suffix(".ecomp")
    assert ecomp_main([
        "zip",
        str(single_alignment),
    ]) == 0
    capsys.readouterr()

    assert ecomp_main([
        "percentage_identity",
        str(single_archive),
    ]) == 0
    assert capsys.readouterr().out.strip() == "nan"

    assert ecomp_main([
        "relative_composition_variability",
        str(single_archive),
    ]) == 0
    assert float(capsys.readouterr().out.strip()) == 0.0
