"""Read/write helpers for the binary `.ecomp` payload."""

from __future__ import annotations

import json
import struct
from pathlib import Path
from typing import Any, Tuple

from .config import FORMAT_VERSION_TUPLE, HEADER_MAGIC, HEADER_STRUCT, METADATA_SUFFIX

HEADER_SIZE = struct.calcsize(HEADER_STRUCT)


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
    """Persist metadata as canonical JSON."""

    path = Path(path)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


def read_metadata(path: str | Path) -> dict[str, Any]:
    """Load metadata JSON from disk."""

    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def derive_metadata_path(ecomp_path: Path) -> Path:
    """Return the default metadata path derived from *ecomp_path*."""

    return ecomp_path.with_suffix(METADATA_SUFFIX)
