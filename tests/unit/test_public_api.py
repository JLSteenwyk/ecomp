from __future__ import annotations

from pathlib import Path

import pytest

from ecomp import __version__, compress_file, decompress_file
from ecomp.storage import read_archive, read_metadata, write_archive


def _write_alignment(path: Path, body: str) -> None:
    path.write_text(body)


def test_compress_file_writes_default_metadata(tmp_path: Path) -> None:
    alignment = tmp_path / "example.fasta"
    _write_alignment(alignment, ">s1\nAC\n>s2\nAC\n")

    archive_path, metadata_path = compress_file(alignment)
    assert archive_path.exists()
    assert metadata_path.exists()

    metadata = read_metadata(metadata_path)
    assert metadata["num_sequences"] == 2
    assert metadata["alignment_length"] == 2


def test_decompress_file_detects_checksum_mismatch(tmp_path: Path) -> None:
    alignment = tmp_path / "example.fasta"
    _write_alignment(alignment, ">s1\nAA\n>s2\nAA\n")

    archive_path, _ = compress_file(alignment)
    payload, metadata, _ = read_archive(archive_path)
    metadata["checksum_sha256"] = "deadbeef"
    write_archive(archive_path, payload, metadata)

    with pytest.raises(ValueError, match="Checksum mismatch"):
        decompress_file(archive_path)


def test_version_attribute_present() -> None:
    assert isinstance(__version__, str)
    assert __version__
