import json
import math
from pathlib import Path

import pytest

from ecomp.cli import main as ecomp_main
from ecomp.diagnostics.metrics import (
    column_base_counts as metrics_column_base_counts,
    column_gap_fraction as metrics_column_gap_fraction,
    column_shannon_entropy as metrics_column_shannon_entropy,
    constant_columns as metrics_constant_columns,
    majority_rule_consensus,
    pairwise_identity_matrix,
    parsimony_informative_columns,
    percentage_identity as metrics_percentage_identity,
    relative_composition_variability,
    variable_site_count as metrics_variable_site_count,
)
from ecomp.io import read_alignment
from ecomp.storage import read_archive, write_archive

FASTA_CONTENT = ">tax1\nACGTACGT\n>tax2\nACGTTCGT\n"
TREE_CONTENT = "(tax1:0.1,tax2:0.2);\n"


def _write_alignment(path: Path) -> None:
    path.write_text(FASTA_CONTENT)


def test_cli_happy_path_round_trip(tmp_path, capsys):
    alignment = tmp_path / "toy.fasta"
    _write_alignment(alignment)
    archive = tmp_path / "toy.ecomp"
    metadata = tmp_path / "toy.json"
    restored = tmp_path / "toy_restored.fasta"

    assert ecomp_main([
        "compress",
        str(alignment),
        "-o",
        str(archive),
        "-m",
        str(metadata),
        "--stats",
    ]) == 0
    compress_out = capsys.readouterr().out
    assert "Created" in compress_out
    assert "Stats:" in compress_out

    assert "Metadata copy" in compress_out

    assert archive.exists() and metadata.exists()

    assert ecomp_main([
        "decompress",
        str(archive),
        "-m",
        str(metadata),
        "-o",
        str(restored),
    ]) == 0
    decompress_out = capsys.readouterr().out
    assert "Wrote alignment" in decompress_out

    original = read_alignment(alignment)
    round_trip = read_alignment(restored)
    assert original.sequences == round_trip.sequences
    assert original.ids == round_trip.ids

    assert ecomp_main([
        "inspect",
        str(archive),
        "--summary",
    ]) == 0
    summary = capsys.readouterr().out
    assert "Codec:" in summary
    assert "Sequences:" in summary

    capsys.readouterr()
    assert ecomp_main(["inspect", str(archive)]) == 0
    raw_metadata = json.loads(capsys.readouterr().out)
    assert raw_metadata["num_sequences"] == len(original.ids)
    assert raw_metadata["alignment_length"] == len(original.sequences[0])


def test_cli_decompress_no_checksum_allows_mismatch(tmp_path, capsys):
    alignment = tmp_path / "with_tree.fasta"
    _write_alignment(alignment)
    (tmp_path / "with_tree.tree").write_text(TREE_CONTENT)

    # Use defaults for output paths to exercise derive_metadata_path
    assert ecomp_main([
        "compress",
        str(alignment),
        "--tree",
        str(tmp_path / "with_tree.tree"),
        "--stats",
    ]) == 0
    archive = alignment.with_suffix(".ecomp")
    payload, metadata, _ = read_archive(archive)
    metadata["checksum_sha256"] = "deadbeef"
    write_archive(archive, payload, metadata)

    capsys.readouterr()
    assert ecomp_main([
        "decompress",
        str(archive),
        "--no-checksum",
    ]) == 0
    output = alignment.with_suffix(".fasta")
    assert output.exists()
    round_trip = read_alignment(output)
    assert len(round_trip.ids) == 2

    # Ensure inspect raw JSON still works with derived metadata path
    capsys.readouterr()
    assert ecomp_main(["inspect", str(archive)]) == 0
    metadata_json = json.loads(capsys.readouterr().out)
    assert metadata_json["ordering_strategy"]


def test_cli_zip_raises_when_alignment_missing(tmp_path):
    missing = tmp_path / "missing.fasta"
    with pytest.raises(SystemExit, match="Alignment not found"):
        ecomp_main(["zip", str(missing)])


def test_cli_zip_raises_when_tree_missing(tmp_path):
    alignment = tmp_path / "toy.fasta"
    _write_alignment(alignment)
    with pytest.raises(SystemExit, match="Tree file not found"):
        ecomp_main([
            "zip",
            str(alignment),
            "--tree",
            str(tmp_path / "nope.tree"),
        ])


def test_cli_zip_raises_when_tree_unreadable(tmp_path):
    alignment = tmp_path / "toy.fasta"
    _write_alignment(alignment)
    bad_tree = tmp_path / "bad.tree"
    bad_tree.mkdir()
    with pytest.raises(SystemExit, match="Failed to read tree file"):
        ecomp_main([
            "zip",
            str(alignment),
            "--tree",
            str(bad_tree),
        ])


def test_cli_inspect_raises_when_archive_missing(tmp_path):
    with pytest.raises(SystemExit, match="Archive not found"):
        ecomp_main(["inspect", str(tmp_path / "ghost.ecomp")])


def test_cli_alignment_length_alias(tmp_path, capsys):
    alignment = tmp_path / "toy.fasta"
    alignment.write_text(">s1\nAC--\n>s2\nACGT\n")
    archive = tmp_path / "toy.ecomp"
    assert ecomp_main([
        "zip",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0
    capsys.readouterr()

    assert ecomp_main([
        "alignment_length",
        str(archive),
    ]) == 0
    assert capsys.readouterr().out.strip() == "4"

    assert ecomp_main([
        "alignment_length_excluding_gaps",
        str(archive),
    ]) == 0
    assert capsys.readouterr().out.strip() == "4"

    assert ecomp_main([
        "len_total",
        str(archive),
        "--no-checksum",
    ]) == 0
    assert capsys.readouterr().out.strip() == "4"

    assert ecomp_main([
        "len_no_gaps",
        str(archive),
        "--no-checksum",
    ]) == 0
    assert capsys.readouterr().out.strip() == "4"


def test_cli_metric_subcommands(tmp_path, capsys):
    alignment = tmp_path / "metrics.fasta"
    alignment.write_text(">s1\nACGT\n>s2\nA-GT\n>s3\nATGT\n")
    archive = tmp_path / "metrics.ecomp"

    assert ecomp_main([
        "zip",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0
    capsys.readouterr()

    frame = read_alignment(alignment)

    assert ecomp_main([
        "consensus_sequence",
        str(archive),
        "--header",
        "metric_consensus",
    ]) == 0
    consensus_output = capsys.readouterr().out.strip().splitlines()
    assert consensus_output == [">metric_consensus", majority_rule_consensus(frame)]

    assert ecomp_main([
        "column_base_counts",
        str(archive),
        "--include-gaps",
        "--indent",
        "0",
    ]) == 0
    counts_output = json.loads(capsys.readouterr().out)
    expected_counts = [
        {"column": idx + 1, "counts": dict(counter)}
        for idx, counter in enumerate(
            metrics_column_base_counts(frame, include_gaps=True)
        )
    ]
    assert counts_output == expected_counts

    assert ecomp_main([
        "gap_fraction",
        str(archive),
    ]) == 0
    gap_lines = capsys.readouterr().out.strip().splitlines()
    expected_gap = metrics_column_gap_fraction(frame)
    for idx, line in enumerate(gap_lines, start=1):
        column, value = line.split("\t")
        assert int(column) == idx
        assert value == f"{expected_gap[idx - 1]:.6f}"

    assert ecomp_main([
        "shannon_entropy",
        str(archive),
    ]) == 0
    entropy_lines = capsys.readouterr().out.strip().splitlines()
    expected_entropy = metrics_column_shannon_entropy(frame)
    for idx, line in enumerate(entropy_lines, start=1):
        column, value = line.split("\t")
        assert int(column) == idx
        assert value == f"{expected_entropy[idx - 1]:.6f}"

    assert ecomp_main([
        "parsimony_informative_sites",
        str(archive),
    ]) == 0
    parsimony_lines = capsys.readouterr().out.strip().splitlines()
    parsimony_mask = parsimony_informative_columns(frame)
    total_parsimony = sum(parsimony_mask)
    assert parsimony_lines[0] == f"total\t{total_parsimony}"
    if total_parsimony:
        expected_indices = " ".join(
            str(idx + 1) for idx, value in enumerate(parsimony_mask) if value
        )
        assert parsimony_lines[1] == f"indices\t{expected_indices}"
    else:
        assert len(parsimony_lines) == 1

    assert ecomp_main([
        "constant_columns",
        str(archive),
    ]) == 0
    constant_lines = capsys.readouterr().out.strip().splitlines()
    constant_mask = metrics_constant_columns(frame)
    total_constant = sum(1 for value in constant_mask if value)
    assert constant_lines[0] == f"total\t{total_constant}"
    if total_constant:
        expected_indices = " ".join(
            str(idx + 1) for idx, value in enumerate(constant_mask) if value
        )
        assert constant_lines[1] == f"indices\t{expected_indices}"

    assert ecomp_main([
        "pairwise_identity",
        str(archive),
    ]) == 0
    pairwise_lines = capsys.readouterr().out.strip().splitlines()
    result = pairwise_identity_matrix(frame)
    ids = frame.ids
    assert pairwise_lines[0].split("\t") == ["id", *ids]
    for row_index, line in enumerate(pairwise_lines[1 : 1 + len(ids)]):
        parts = line.split("\t")
        assert parts[0] == ids[row_index]
        expected = [
            "nan" if math.isnan(value) else f"{value:.6f}"
            for value in result.matrix[row_index]
        ]
        assert parts[1:] == expected

    coverage_offset = 1 + len(ids)
    assert pairwise_lines[coverage_offset] == "# coverage"
    assert pairwise_lines[coverage_offset + 1].split("\t") == ["id", *ids]
    for row_index, line in enumerate(
        pairwise_lines[coverage_offset + 2 : coverage_offset + 2 + len(ids)]
    ):
        parts = line.split("\t")
        assert parts[0] == ids[row_index]
        assert parts[1:] == [str(value) for value in result.coverage[row_index]]

    assert ecomp_main([
        "variable_sites",
        str(archive),
    ]) == 0
    variable_output = capsys.readouterr().out.strip()
    assert variable_output == str(metrics_variable_site_count(frame))

    assert ecomp_main([
        "percentage_identity",
        str(archive),
    ]) == 0
    percentage_output = capsys.readouterr().out.strip()
    expected_percentage = metrics_percentage_identity(frame)
    expected_percentage_str = (
        "nan"
        if math.isnan(expected_percentage)
        else f"{expected_percentage:.6f}"
    )
    assert percentage_output == expected_percentage_str

    assert ecomp_main([
        "rcv",
        str(archive),
    ]) == 0
    rcv_output = capsys.readouterr().out.strip()
    expected_rcv = relative_composition_variability(frame)
    assert rcv_output == f"{expected_rcv:.6f}"
