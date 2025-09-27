from pathlib import Path

import pytest

from evolutionary_compression.cli import (
    _normalize_suffix,
    _total_input_size,
    _verify_checksum,
)
from evolutionary_compression.diagnostics.checksums import alignment_checksum


def test_normalize_suffix_defaults_to_ecbt():
    assert _normalize_suffix("") == ".ecbt"


def test_normalize_suffix_adds_leading_dot():
    assert _normalize_suffix("bundle") == ".bundle"


def test_total_input_size_includes_tree(tmp_path: Path):
    alignment = tmp_path / "aln.fasta"
    alignment.write_bytes(b"ACGT")
    tree = tmp_path / "tree.nwk"
    tree.write_bytes(b"(a:0.1,b:0.2);")
    total = _total_input_size(alignment, tree)
    assert total == alignment.stat().st_size + tree.stat().st_size


def test_total_input_size_without_tree(tmp_path: Path):
    alignment = tmp_path / "aln.fasta"
    alignment.write_bytes(b"ACGT")
    assert _total_input_size(alignment, None) == alignment.stat().st_size


def test_verify_checksum_passes():
    sequences = ["ACGT", "ACGA"]
    checksum = alignment_checksum(sequences)
    metadata = {"checksum_sha256": checksum}
    _verify_checksum(sequences, metadata)


def test_verify_checksum_raises_on_mismatch():
    sequences = ["ACGT"]
    metadata = {"checksum_sha256": alignment_checksum(["AAAA"])}
    with pytest.raises(SystemExit):
        _verify_checksum(sequences, metadata)


def test_verify_checksum_ignored_when_missing():
    _verify_checksum(["ACGT"], {})
