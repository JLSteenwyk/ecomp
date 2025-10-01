import json
from pathlib import Path

import pytest

from ecomp.cli import main as ecomp_main
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
        "alignment_length_excluding_gaps",
        str(archive),
    ]) == 0
    assert capsys.readouterr().out.strip() == "4"

    assert ecomp_main([
        "alignment_compressed_length",
        str(archive),
    ]) == 0
    assert capsys.readouterr().out.strip().isdigit()

    assert ecomp_main([
        "compressed_len",
        str(archive),
        "--no-checksum",
    ]) == 0
    assert capsys.readouterr().out.strip().isdigit()

    assert ecomp_main([
        "len_no_gaps",
        str(archive),
        "--no-checksum",
    ]) == 0
    assert capsys.readouterr().out.strip() == "4"
