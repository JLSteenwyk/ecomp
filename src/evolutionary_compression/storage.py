"""Read/write helpers for the binary `.ecomp` payload."""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path
from typing import Any, Tuple

from .config import FORMAT_VERSION_TUPLE, HEADER_MAGIC, HEADER_STRUCT, METADATA_SUFFIX

HEADER_SIZE = struct.calcsize(HEADER_STRUCT)
_METADATA_COMPRESSED_MAGIC = b"ECMZ"
_METADATA_CODEC_VERSION = 1


def write_payload(path: str | Path, payload: bytes) -> Path:
    """Write payload bytes to *path* with an `.ecomp` header."""

    path = Path(path)
    header = struct.pack(HEADER_STRUCT, HEADER_MAGIC, *FORMAT_VERSION_TUPLE, len(payload))
    with path.open("wb") as handle:
        handle.write(header)
        handle.write(payload)
    return path


def read_payload(path: str | Path) -> bytes:
    """Load payload bytes and validate the `.ecomp` header."""

    path = Path(path)
    data = path.read_bytes()
    if len(data) < HEADER_SIZE:
        raise ValueError("File is too short to be a valid .ecomp payload")
    magic, major, minor, patch, length = struct.unpack(HEADER_STRUCT, data[:HEADER_SIZE])
    if magic != HEADER_MAGIC:
        raise ValueError("Invalid .ecomp magic header")
    payload = data[HEADER_SIZE:]
    if len(payload) != length:
        raise ValueError("Payload length does not match header metadata")
    return payload


def write_metadata(path: str | Path, metadata: dict[str, Any]) -> Path:
    """Persist metadata, optionally applying compression to shrink overhead."""

    path = Path(path)
    json_bytes = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode("utf-8")

    # Attempt to compress; fall back to plain JSON if it does not help.
    compressed = zlib.compress(json_bytes, level=9)
    use_compressed = len(compressed) + len(_METADATA_COMPRESSED_MAGIC) + 1 < len(json_bytes)

    if use_compressed:
        payload = _METADATA_COMPRESSED_MAGIC + bytes([_METADATA_CODEC_VERSION]) + compressed
        path.write_bytes(payload)
    else:
        path.write_bytes(json_bytes + b"\n")
    return path


def read_metadata(path: str | Path) -> dict[str, Any]:
    """Load metadata JSON from disk."""

    path = Path(path)
    data = path.read_bytes()

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


def derive_metadata_path(ecomp_path: Path) -> Path:
    """Return the default metadata path derived from *ecomp_path*."""

    return ecomp_path.with_suffix(METADATA_SUFFIX)
