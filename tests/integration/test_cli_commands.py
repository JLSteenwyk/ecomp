import json
from pathlib import Path

import pytest

from evolutionary_compression.cli import main as ecomp_main
from evolutionary_compression.io import read_alignment
from evolutionary_compression.storage import read_metadata

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
        "compress",
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
        "decompress",
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
        "compress",
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
        "decompress",
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
        "compress",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0

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
        "compress",
        str(alignment),
        "-o",
        str(archive),
    ]) == 0

    capsys.readouterr()  # clear compress output

    assert ecomp_main([
        "inspect",
        str(archive),
    ]) == 0

    metadata = json.loads(capsys.readouterr().out)
    assert metadata["codec"] == "ecomp"
    assert metadata["num_sequences"] == 2
