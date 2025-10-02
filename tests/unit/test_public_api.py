from __future__ import annotations

from pathlib import Path

import pytest

from ecomp import __version__, compress_file, decompress_file, eunzip, ezip
from ecomp.storage import read_archive, read_metadata, write_archive


def _write_alignment(path: Path, body: str) -> None:
    path.write_text(body)


def test_compress_file_writes_default_metadata(tmp_path: Path) -> None:
    alignment = tmp_path / "example.fasta"
    _write_alignment(alignment, ">s1\nAC\n>s2\nAC\n")

    metadata_copy = tmp_path / "example.json"
    archive_path, metadata_path = compress_file(alignment, metadata_path=metadata_copy)
    assert archive_path.exists()
    assert metadata_path.exists()

    metadata = read_metadata(metadata_path)
    assert metadata["num_sequences"] == 2
    assert metadata["alignment_length"] == 2


def test_decompress_file_detects_checksum_mismatch(tmp_path: Path) -> None:
    alignment = tmp_path / "example.fasta"
    _write_alignment(alignment, ">s1\nAA\n>s2\nAA\n")

    metadata_copy = tmp_path / "example.json"
    archive_path, _ = compress_file(alignment, metadata_path=metadata_copy)
    payload, metadata, _ = read_archive(archive_path)
    metadata["checksum_sha256"] = "deadbeef"
    write_archive(archive_path, payload, metadata)

    with pytest.raises(ValueError, match="Checksum mismatch"):
        decompress_file(archive_path)


def test_version_attribute_present() -> None:
    assert isinstance(__version__, str)
    assert __version__


def test_ezip_eunzip_round_trip(tmp_path: Path) -> None:
    alignment = tmp_path / "alias.fasta"
    alignment.write_text(">s1\nACGT\n>s2\nACGT\n")

    metadata_copy = tmp_path / "alias.json"
    archive_path, metadata_path = ezip(alignment, metadata_path=metadata_copy)
    assert archive_path.exists()
    assert metadata_path.exists()

    output_path = eunzip(archive_path)
    assert output_path.exists()
