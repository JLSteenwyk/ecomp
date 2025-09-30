import json
from pathlib import Path

from evolutionary_compression.cli import main as ecomp_main
from evolutionary_compression.io import read_alignment
from evolutionary_compression.storage import read_metadata, write_metadata

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
    metadata = alignment.with_suffix(".json")

    # Corrupt checksum to ensure --no-checksum takes the skip branch
    meta = read_metadata(metadata)
    meta["checksum_sha256"] = "deadbeef"
    write_metadata(metadata, meta)

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
