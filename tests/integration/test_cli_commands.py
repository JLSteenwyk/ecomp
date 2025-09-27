import json
from pathlib import Path

import pytest

from evolutionary_compression.cli import main as ecomp_main
from evolutionary_compression.io import read_alignment
from evolutionary_compression.storage import derive_metadata_path, read_metadata
from evolutionary_compression.compression.phylo_bundle import _parse_newick  # type: ignore

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
def test_cli_phylo_bundle_round_trip(tmp_path: Path) -> None:
    alignment = tmp_path / "bundle.fasta"
    alignment.write_text(FASTA_BODY)
    tree = tmp_path / "bundle.tree"
    tree.write_text(TREE_BODY)

    archive = tmp_path / "bundle.ecbt"
    metadata = derive_metadata_path(archive)

    assert ecomp_main([
        "compress",
        str(alignment),
        str(tree),
        "-o",
        str(archive),
        "--stats",
    ]) == 0

    assert archive.exists()
    assert metadata.exists()

    restored_alignment = tmp_path / "roundtrip.fasta"
    restored_tree = tmp_path / "roundtrip.tree"

    assert ecomp_main([
        "decompress",
        str(archive),
        "-o",
        str(restored_alignment),
        "-t",
        str(restored_tree),
    ]) == 0

    original = read_alignment(alignment)
    round_tripped = read_alignment(restored_alignment)
    assert original.sequences == round_tripped.sequences
    assert original.ids == round_tripped.ids

    original_tree = _parse_newick(TREE_BODY)
    restored_tree = _parse_newick((tmp_path / "roundtrip.tree").read_text())
    assert _tree_signature(original_tree) == _tree_signature(restored_tree)


@pytest.mark.integration
def test_cli_alignment_codec_ignores_tree_when_forced(tmp_path: Path, capsys) -> None:
    alignment = tmp_path / "forced.fasta"
    alignment.write_text(FASTA_BODY)
    tree = tmp_path / "forced.tree"
    tree.write_text(TREE_BODY)

    archive = tmp_path / "forced.ecomp"

    assert ecomp_main([
        "compress",
        str(alignment),
        str(tree),
        "--codec",
        "alignment",
        "-o",
        str(archive),
    ]) == 0

    captured = capsys.readouterr()
    assert "Tree provided" in captured.err

    metadata = read_metadata(derive_metadata_path(archive))
    assert metadata["codec"] == "ecomp"


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


@pytest.mark.integration
def test_cli_requires_tree_when_phylo_codec_requested(tmp_path: Path) -> None:
    alignment = tmp_path / "phylo_required.fasta"
    alignment.write_text(FASTA_BODY)

    with pytest.raises(SystemExit):
        ecomp_main([
            "compress",
            str(alignment),
            "--codec",
            "phylo",
        ])


@pytest.mark.integration
def test_cli_errors_when_tree_missing(tmp_path: Path) -> None:
    alignment = tmp_path / "missing_tree.fasta"
    alignment.write_text(FASTA_BODY)
    missing_tree = tmp_path / "missing.tree"

    with pytest.raises(SystemExit):
        ecomp_main([
            "compress",
            str(alignment),
            str(missing_tree),
        ])


@pytest.mark.integration
def test_cli_respects_custom_bundle_suffix(tmp_path: Path) -> None:
    alignment = tmp_path / "custom.fasta"
    alignment.write_text(FASTA_BODY)
    tree = tmp_path / "custom.tree"
    tree.write_text(TREE_BODY)

    assert ecomp_main([
        "compress",
        str(alignment),
        str(tree),
        "--bundle-suffix",
        "ecbundle",
    ]) == 0

    archive = alignment.with_suffix(".ecbundle")
    assert archive.exists()
    metadata = derive_metadata_path(archive)
    assert metadata.exists()
    assert read_metadata(metadata)["codec"] == "phylo-bundle"


@pytest.mark.integration
def test_cli_decompress_phylo_uses_default_outputs(tmp_path: Path) -> None:
    alignment = tmp_path / "default_outputs.fasta"
    alignment.write_text(FASTA_BODY)
    tree = tmp_path / "default_outputs.tree"
    tree.write_text(TREE_BODY)

    assert ecomp_main([
        "compress",
        str(alignment),
        str(tree),
    ]) == 0

    archive = alignment.with_suffix(".ecbt")
    assert ecomp_main([
        "decompress",
        str(archive),
        "--no-checksum",
    ]) == 0

    assert archive.with_suffix(".fasta").exists()
    assert archive.with_suffix(".tree").exists()


@pytest.mark.integration
def test_cli_requires_existing_alignment_path(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.fasta"
    with pytest.raises(SystemExit):
        ecomp_main([
            "compress",
            str(missing),
        ])


@pytest.mark.integration
def test_cli_decompress_requires_existing_archive(tmp_path: Path) -> None:
    missing = tmp_path / "missing.ecomp"
    with pytest.raises(SystemExit):
        ecomp_main([
            "decompress",
            str(missing),
        ])


@pytest.mark.integration
def test_cli_inspect_with_explicit_metadata_path(tmp_path: Path, capsys) -> None:
    alignment = tmp_path / "inspect_meta.fasta"
    alignment.write_text(FASTA_BODY)
    archive = tmp_path / "inspect_meta.ecomp"
    metadata = tmp_path / "inspect_meta.json"

    assert ecomp_main([
        "compress",
        str(alignment),
        "-o",
        str(archive),
        "-m",
        str(metadata),
    ]) == 0

    capsys.readouterr()

    assert ecomp_main([
        "inspect",
        str(archive),
        "-m",
        str(metadata),
    ]) == 0
    metadata_out = json.loads(capsys.readouterr().out)
    assert metadata_out["sequence_ids"] == ["s1", "s2"]


@pytest.mark.integration
def test_cli_inspect_summary_reports_bundle_leaf_count(tmp_path: Path, capsys) -> None:
    alignment = tmp_path / "inspect_bundle.fasta"
    alignment.write_text(FASTA_BODY)
    tree = tmp_path / "inspect_bundle.tree"
    tree.write_text(TREE_BODY)

    assert ecomp_main([
        "compress",
        str(alignment),
        str(tree),
    ]) == 0

    capsys.readouterr()

    archive = alignment.with_suffix(".ecbt")
    assert ecomp_main([
        "inspect",
        str(archive),
        "--summary",
    ]) == 0

    summary = capsys.readouterr().out.splitlines()
    assert any("Bundle leaf count" in line for line in summary)

def _tree_signature(node) -> list[tuple[str, float]]:
    edges: list[tuple[str, float]] = []

    def visit(current):
        label = current.label or ""
        edges.append((label, round(current.length, 4)))
        for child in current.children:
            visit(child)

    visit(node)
    return sorted(edges)
