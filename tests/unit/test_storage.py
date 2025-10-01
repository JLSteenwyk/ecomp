import struct
from pathlib import Path

import pytest

from ecomp.config import (
    FORMAT_VERSION_TUPLE,
    HEADER_MAGIC,
    HEADER_STRUCT,
    LEGACY_HEADER_STRUCT,
)
from ecomp.storage import (
    derive_metadata_path,
    read_archive,
    read_metadata,
    read_payload,
    write_archive,
    write_metadata,
    write_payload,
)


def test_write_and_read_archive_round_trip(tmp_path: Path):
    payload = b"hello-world"
    metadata = {"codec": "ecomp", "alignment_length": 123}
    path = tmp_path / "example.ecomp"

    result_path = write_archive(path, payload, metadata)
    assert result_path == path

    restored_payload, restored_metadata, version = read_archive(path)
    assert restored_payload == payload
    assert restored_metadata == metadata
    assert version == FORMAT_VERSION_TUPLE
    assert read_payload(path) == payload


def test_write_payload_alias(tmp_path: Path):
    payload = b"alias"
    metadata = {"foo": "bar"}
    path = tmp_path / "alias.ecomp"

    write_payload(path, payload, metadata)
    _, restored_metadata, _ = read_archive(path)
    assert restored_metadata == metadata


def test_read_archive_rejects_short_files(tmp_path: Path):
    path = tmp_path / "broken.ecomp"
    path.write_bytes(b"ECOMP")
    with pytest.raises(ValueError):
        read_archive(path)


def test_read_archive_rejects_bad_magic(tmp_path: Path):
    path = tmp_path / "bad_magic.ecomp"
    bad_magic = HEADER_MAGIC[:-1] + b"Z"
    header = struct.pack(
        HEADER_STRUCT,
        bad_magic,
        *FORMAT_VERSION_TUPLE,
        0,
        0,
    )
    path.write_bytes(header)
    with pytest.raises(ValueError):
        read_archive(path)


def test_read_archive_length_mismatch(tmp_path: Path):
    path = tmp_path / "bad_length.ecomp"
    header = struct.pack(
        HEADER_STRUCT,
        HEADER_MAGIC,
        *FORMAT_VERSION_TUPLE,
        10,
        0,
    )
    path.write_bytes(header + b"abc")
    with pytest.raises(ValueError):
        read_archive(path)


def test_metadata_round_trip(tmp_path: Path):
    metadata = {"codec": "ecomp", "alignment_length": 100}
    path = tmp_path / "example.json"
    result = write_metadata(path, metadata)
    assert result == path
    assert read_metadata(path) == metadata


def test_legacy_archive_requires_metadata(tmp_path: Path):
    payload = b"legacy"
    archive = tmp_path / "legacy.ecomp"
    header = struct.pack(
        LEGACY_HEADER_STRUCT,
        HEADER_MAGIC,
        0,
        1,
        0,
        len(payload),
    )
    archive.write_bytes(header + payload)

    with pytest.raises(FileNotFoundError):
        read_archive(archive)


def test_legacy_archive_uses_sidecar_metadata(tmp_path: Path):
    payload = b"legacy"
    metadata = {"codec": "legacy", "alignment_length": 42}
    archive = tmp_path / "legacy.ecomp"
    metadata_path = archive.with_suffix(".json")

    header = struct.pack(
        LEGACY_HEADER_STRUCT,
        HEADER_MAGIC,
        0,
        1,
        0,
        len(payload),
    )
    archive.write_bytes(header + payload)
    write_metadata(metadata_path, metadata)

    restored_payload, restored_metadata, version = read_archive(archive)
    assert restored_payload == payload
    assert restored_metadata == metadata
    assert version == (0, 1, 0)


def test_derive_metadata_path_uses_json_suffix(tmp_path: Path):
    payload_path = tmp_path / "example.ecomp"
    assert derive_metadata_path(payload_path) == payload_path.with_suffix(".json")
