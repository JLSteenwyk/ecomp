import struct
from pathlib import Path

import pytest

from evolutionary_compression.storage import (
    derive_metadata_path,
    read_metadata,
    read_payload,
    write_metadata,
    write_payload,
)
from evolutionary_compression.config import HEADER_MAGIC, HEADER_STRUCT


def test_write_and_read_payload(tmp_path: Path):
    payload = b"hello-world"
    path = tmp_path / "example.ecomp"
    result_path = write_payload(path, payload)
    assert result_path == path
    assert read_payload(path) == payload


def test_read_payload_rejects_short_files(tmp_path: Path):
    path = tmp_path / "broken.ecomp"
    path.write_bytes(b"ECOMP")
    with pytest.raises(ValueError):
        read_payload(path)


def test_read_payload_rejects_bad_magic(tmp_path: Path):
    # craft valid-sized header but with different magic
    path = tmp_path / "bad_magic.ecomp"
    bad_magic = HEADER_MAGIC[:-1] + b"Z"
    header = struct.pack(HEADER_STRUCT, bad_magic, 0, 1, 0, 0)
    path.write_bytes(header)
    with pytest.raises(ValueError):
        read_payload(path)


def test_read_payload_length_mismatch(tmp_path: Path):
    path = tmp_path / "bad_length.ecomp"
    header = struct.pack(HEADER_STRUCT, HEADER_MAGIC, 0, 1, 0, 10)
    path.write_bytes(header + b"abc")
    with pytest.raises(ValueError):
        read_payload(path)


def test_metadata_round_trip(tmp_path: Path):
    metadata = {"codec": "ecomp", "alignment_length": 100}
    path = tmp_path / "example.json"
    result = write_metadata(path, metadata)
    assert result == path
    assert read_metadata(path) == metadata


def test_derive_metadata_path_uses_json_suffix(tmp_path: Path):
    payload_path = tmp_path / "example.ecomp"
    assert derive_metadata_path(payload_path) == payload_path.with_suffix(".json")
