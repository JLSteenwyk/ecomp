"""Read/write helpers for the binary `.ecomp` payload."""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path
from typing import Any, Tuple

from .config import (
    FORMAT_VERSION_TUPLE,
    HEADER_MAGIC,
    HEADER_STRUCT,
    LEGACY_HEADER_STRUCT,
    METADATA_SUFFIX,
)

_HEADER_LEGACY_STRUCT = LEGACY_HEADER_STRUCT
_HEADER_LEGACY_SIZE = struct.calcsize(_HEADER_LEGACY_STRUCT)
_HEADER_CURRENT_STRUCT = HEADER_STRUCT
_HEADER_CURRENT_SIZE = struct.calcsize(_HEADER_CURRENT_STRUCT)
_INLINE_METADATA_VERSION: Tuple[int, int, int] = (0, 2, 0)

_METADATA_COMPRESSED_MAGIC = b"ECMZ"
_METADATA_CODEC_VERSION = 1


def write_payload(path: str | Path, payload: bytes, metadata: dict[str, Any]) -> Path:
    """Write payload and metadata to *path* using the current `.ecomp` format."""

    return write_archive(path, payload, metadata)


def read_payload(path: str | Path) -> bytes:
    """Load payload bytes and validate the `.ecomp` header."""

    payload, _, _ = read_archive(path)
    return payload


def write_archive(path: str | Path, payload: bytes, metadata: dict[str, Any]) -> Path:
    """Persist *payload* and *metadata* in a single `.ecomp` file."""

    path = Path(path)
    metadata_bytes = _encode_metadata(metadata, add_trailing_newline=False)
    header = struct.pack(
        _HEADER_CURRENT_STRUCT,
        HEADER_MAGIC,
        *FORMAT_VERSION_TUPLE,
        len(payload),
        len(metadata_bytes),
    )

    with path.open("wb") as handle:
        handle.write(header)
        handle.write(payload)
        handle.write(metadata_bytes)
    return path


def read_archive(
    path: str | Path,
    *,
    metadata_path: str | Path | None = None,
) -> tuple[bytes, dict[str, Any], Tuple[int, int, int]]:
    """Return ``(payload, metadata, version)`` for *path*.

    For legacy archives (format < 0.2.0) this will read metadata from *metadata_path*
    or the default JSON sidecar if present.
    """

    path = Path(path)
    data = path.read_bytes()
    if len(data) < _HEADER_LEGACY_SIZE:
        raise ValueError("File is too short to be a valid .ecomp payload")

    magic, major, minor, patch, payload_length = struct.unpack(
        _HEADER_LEGACY_STRUCT, data[:_HEADER_LEGACY_SIZE]
    )
    if magic != HEADER_MAGIC:
        raise ValueError("Invalid .ecomp magic header")

    version = (major, minor, patch)
    offset = _HEADER_LEGACY_SIZE

    if version >= _INLINE_METADATA_VERSION:
        if len(data) < _HEADER_CURRENT_SIZE:
            raise ValueError("File is truncated; missing metadata header")
        metadata_length = struct.unpack(">Q", data[offset:_HEADER_CURRENT_SIZE])[0]
        offset = _HEADER_CURRENT_SIZE
    else:
        metadata_length = None

    payload_end = offset + payload_length
    if len(data) < payload_end:
        raise ValueError("Payload length does not match header metadata")
    payload = data[offset:payload_end]

    if metadata_length is not None:
        metadata_end = payload_end + metadata_length
        if len(data) < metadata_end:
            raise ValueError("Metadata length does not match header metadata")
        metadata_bytes = data[payload_end:metadata_end]
        metadata_dict = _decode_metadata(metadata_bytes)
    else:
        metadata_location = Path(metadata_path) if metadata_path else derive_metadata_path(path)
        if not metadata_location.exists():
            raise FileNotFoundError(
                "Metadata sidecar is required for legacy archive formats"
            )
        metadata_dict = read_metadata(metadata_location)

    return payload, metadata_dict, version


def write_metadata(path: str | Path, metadata: dict[str, Any]) -> Path:
    """Persist metadata, optionally applying compression to shrink overhead."""

    path = Path(path)
    data = _encode_metadata(metadata, add_trailing_newline=True)
    path.write_bytes(data)
    return path


def read_metadata(path: str | Path) -> dict[str, Any]:
    """Load metadata JSON from disk."""

    path = Path(path)
    return _decode_metadata(path.read_bytes())


def derive_metadata_path(ecomp_path: Path) -> Path:
    """Return the default metadata path derived from *ecomp_path*."""

    return ecomp_path.with_suffix(METADATA_SUFFIX)


def _encode_metadata(metadata: dict[str, Any], *, add_trailing_newline: bool) -> bytes:
    json_bytes = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode("utf-8")
    compressed = zlib.compress(json_bytes, level=9)
    use_compressed = len(compressed) + len(_METADATA_COMPRESSED_MAGIC) + 1 < len(json_bytes)
    if use_compressed:
        return _METADATA_COMPRESSED_MAGIC + bytes([_METADATA_CODEC_VERSION]) + compressed
    if add_trailing_newline:
        return json_bytes + b"\n"
    return json_bytes


def _decode_metadata(data: bytes) -> dict[str, Any]:
    if data.startswith(_METADATA_COMPRESSED_MAGIC):
        if len(data) < len(_METADATA_COMPRESSED_MAGIC) + 1:
            raise ValueError("Compressed metadata header truncated")
        codec_version = data[len(_METADATA_COMPRESSED_MAGIC)]
        if codec_version != _METADATA_CODEC_VERSION:
            raise ValueError(f"Unsupported compressed metadata version: {codec_version}")
        json_bytes = zlib.decompress(data[len(_METADATA_COMPRESSED_MAGIC) + 1 :])
    else:
        json_bytes = data
    return json.loads(json_bytes.decode("utf-8"))
