"""Public API for the eComp evolutionary compression toolkit."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

from .compression.pipeline import CompressedAlignment, compress_alignment, decompress_alignment
from .config import DEFAULT_OUTPUT_FORMAT, METADATA_SUFFIX
from .io import AlignmentFrame, read_alignment, write_alignment
from .storage import derive_metadata_path, read_metadata, read_payload, write_metadata, write_payload
from .diagnostics.checksums import alignment_checksum


__all__ = [
    "AlignmentFrame",
    "CompressedAlignment",
    "compress_file",
    "decompress_file",
    "compress_alignment",
    "decompress_alignment",
    "read_alignment",
    "write_alignment",
    "alignment_checksum",
]


def compress_file(
    input_path: str | Path,
    output_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
    input_format: str | None = None,
) -> Tuple[Path, Path]:
    """Compress *input_path* producing `.ecomp` and metadata files."""

    input_path = Path(input_path)
    frame = read_alignment(input_path, fmt=input_format)
    compressed = compress_alignment(frame)

    target_path = Path(output_path) if output_path else input_path.with_suffix(".ecomp")
    metadata_file = (
        Path(metadata_path)
        if metadata_path
        else Path(target_path).with_suffix(METADATA_SUFFIX)
    )

    write_payload(target_path, compressed.payload)
    write_metadata(metadata_file, compressed.metadata)
    return target_path, metadata_file


def decompress_file(
    ecomp_path: str | Path,
    output_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
    output_format: str | None = None,
    validate_checksum: bool = True,
) -> Path:
    """Decompress an `.ecomp` payload back into an alignment file."""

    ecomp_path = Path(ecomp_path)
    payload = read_payload(ecomp_path)
    metadata_file = Path(metadata_path) if metadata_path else derive_metadata_path(ecomp_path)
    metadata = read_metadata(metadata_file)

    frame = decompress_alignment(payload, metadata)

    if validate_checksum:
        checksum = alignment_checksum(frame.sequences)
        expected = metadata.get("checksum_sha256")
        if expected and checksum != expected:
            raise ValueError(
                "Checksum mismatch after decompression: "
                f"expected {expected}, observed {checksum}"
            )

    destination = Path(output_path) if output_path else ecomp_path.with_suffix(f".{output_format or DEFAULT_OUTPUT_FORMAT}")
    write_alignment(frame, destination, fmt=output_format)
    return destination
