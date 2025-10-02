"""Public API for the eComp evolutionary compression toolkit."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

from .compression.pipeline import CompressedAlignment, compress_alignment, decompress_alignment
from .config import DEFAULT_OUTPUT_FORMAT
from .io import AlignmentFrame, alignment_from_sequences, read_alignment, write_alignment
from .storage import derive_metadata_path, read_archive, write_archive, write_metadata
from .diagnostics.checksums import alignment_checksum
from .diagnostics.metrics import (
    PairwiseIdentityResult,
    alignment_length,
    alignment_length_excluding_gaps,
    column_base_counts,
    column_gap_fraction,
    column_shannon_entropy,
    constant_columns,
    majority_rule_consensus,
    pairwise_identity_matrix,
    parsimony_informative_columns,
    parsimony_informative_site_count,
    percentage_identity,
    relative_composition_variability,
    variable_site_count,
)
from ._version import __version__


__all__ = [
    "AlignmentFrame",
    "CompressedAlignment",
    "compress_file",
    "decompress_file",
    "ezip",
    "eunzip",
    "compress_alignment",
    "decompress_alignment",
    "read_alignment",
    "write_alignment",
    "alignment_from_sequences",
    "alignment_checksum",
    "column_base_counts",
    "column_gap_fraction",
    "column_shannon_entropy",
    "parsimony_informative_columns",
    "parsimony_informative_site_count",
    "constant_columns",
    "majority_rule_consensus",
    "alignment_length",
    "alignment_length_excluding_gaps",
    "variable_site_count",
    "percentage_identity",
    "relative_composition_variability",
    "pairwise_identity_matrix",
    "PairwiseIdentityResult",
    "__version__",
]


def compress_file(
    input_path: str | Path,
    output_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
    input_format: str | None = None,
) -> Tuple[Path, Path | None]:
    """Compress *input_path* producing an `.ecomp` archive.

    Returns ``(archive_path, metadata_copy_path)``. Metadata is embedded inside the
    archive by default; pass *metadata_path* to also persist a JSON sidecar.
    """

    input_path = Path(input_path)
    frame = read_alignment(input_path, fmt=input_format)
    compressed = compress_alignment(frame)

    target_path = Path(output_path) if output_path else input_path.with_suffix(".ecomp")
    metadata_file = Path(metadata_path) if metadata_path else None

    write_archive(target_path, compressed.payload, compressed.metadata)

    if metadata_file is not None:
        metadata_file.parent.mkdir(parents=True, exist_ok=True)
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
    metadata_file = Path(metadata_path) if metadata_path else None

    try:
        payload, metadata, _ = read_archive(ecomp_path, metadata_path=metadata_file)
    except FileNotFoundError as exc:
        fallback = metadata_file or derive_metadata_path(ecomp_path)
        raise FileNotFoundError(
            "Metadata sidecar missing for legacy archive; specify it explicitly or "
            f"restore {fallback}"
        ) from exc

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


def ezip(
    input_path: str | Path,
    output_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
    input_format: str | None = None,
) -> Tuple[Path, Path | None]:
    """Alias for :func:`compress_file` mirroring the CLI verb."""

    return compress_file(
        input_path=input_path,
        output_path=output_path,
        metadata_path=metadata_path,
        input_format=input_format,
    )


def eunzip(
    ecomp_path: str | Path,
    output_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
    output_format: str | None = None,
    validate_checksum: bool = True,
) -> Path:
    """Alias for :func:`decompress_file` mirroring the CLI verb."""

    return decompress_file(
        ecomp_path=ecomp_path,
        output_path=output_path,
        metadata_path=metadata_path,
        output_format=output_format,
        validate_checksum=validate_checksum,
    )
