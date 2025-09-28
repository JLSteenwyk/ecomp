"""Unified command-line interface for evolutionary compression workflows."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

from .compression.pipeline import compress_alignment, decompress_alignment
from .config import DEFAULT_OUTPUT_FORMAT
from .diagnostics.checksums import alignment_checksum
from .io import read_alignment, write_alignment
from .storage import (
    derive_metadata_path,
    read_metadata,
    read_payload,
    write_metadata,
    write_payload,
)

ALIGNMENT_SUFFIX = ".ecomp"


# ---------------------------------------------------------------------------
# Parser construction (ClipKIT/PhyKIT style: subcommands with handler binding)
# ---------------------------------------------------------------------------


def _add_compress_arguments(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "compress",
        help="Compress an alignment (optionally bundling a companion Newick tree)",
    )
    parser.add_argument(
        "alignment",
        help="Input alignment in FASTA/PHYLIP format",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Destination archive path (default: alignment stem + .ecomp)",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_path",
        help="Metadata JSON path (default: alongside the archive)",
    )
    parser.add_argument(
        "-f",
        "--input-format",
        dest="alignment_format",
        help="Alignment format hint passed to the parser",
    )
    parser.add_argument(
        "--tree",
        dest="tree_path",
        help="Optional Newick tree used only to guide sequence ordering",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print compression statistics (sizes and ratio)",
    )
    parser.set_defaults(handler=_cmd_compress)


def _add_decompress_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    parser = subparsers.add_parser(
        "decompress",
        help="Restore data from an evolutionary compression archive",
    )
    parser.add_argument("archive", help="Compressed archive produced by `ecomp compress`")
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_path",
        help="Metadata JSON path (default: alongside archive)",
    )
    parser.add_argument(
        "-o",
        "--alignment-output",
        dest="alignment_output",
        help="Alignment output path (default: archive stem + .fasta)",
    )
    parser.add_argument(
        "-F",
        "--format",
        dest="alignment_format",
        default=DEFAULT_OUTPUT_FORMAT,
        help=f"Alignment output format (default: {DEFAULT_OUTPUT_FORMAT})",
    )
    parser.add_argument(
        "--no-checksum",
        action="store_true",
        help="Skip checksum validation during decompression",
    )
    parser.set_defaults(handler=_cmd_decompress)


def _add_inspect_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    parser = subparsers.add_parser(
        "inspect",
        help="Display metadata for an evolutionary compression archive",
    )
    parser.add_argument("archive", help="Archive to inspect")
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_path",
        help="Optional metadata path (default: alongside archive)",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print a human-readable summary instead of raw JSON",
    )
    parser.set_defaults(handler=_cmd_inspect)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ecomp",
        description="Evolutionary compression toolkit for multiple sequence alignments",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_compress_arguments(subparsers)
    _add_decompress_arguments(subparsers)
    _add_inspect_arguments(subparsers)
    return parser


# ---------------------------------------------------------------------------
# Command handlers
# ---------------------------------------------------------------------------


def _cmd_compress(args: argparse.Namespace) -> int:
    alignment_path = Path(args.alignment).expanduser().resolve()
    if not alignment_path.exists():
        raise SystemExit(f"Alignment not found: {alignment_path}")

    tree_path = Path(args.tree_path).expanduser().resolve() if args.tree_path else None
    if tree_path and not tree_path.exists():
        raise SystemExit(f"Tree file not found: {tree_path}")

    frame = read_alignment(alignment_path, fmt=args.alignment_format)
    if tree_path is not None:
        try:
            frame.metadata["tree_newick"] = tree_path.read_text()
        except OSError as exc:
            raise SystemExit(f"Failed to read tree file: {tree_path}") from exc

    compressed = compress_alignment(frame)
    payload = compressed.payload
    metadata = compressed.metadata
    suffix = ALIGNMENT_SUFFIX

    output_path = Path(args.output).expanduser().resolve() if args.output else alignment_path.with_suffix(suffix)
    metadata_path = (
        Path(args.metadata_path).expanduser().resolve()
        if args.metadata_path
        else derive_metadata_path(output_path)
    )

    write_payload(output_path, payload)
    write_metadata(metadata_path, metadata)

    print(f"Created {output_path}")
    print(f"Metadata recorded at {metadata_path}")

    if args.stats:
        original_size = alignment_path.stat().st_size
        compressed_size = output_path.stat().st_size
        ratio = (original_size / compressed_size) if compressed_size else float("inf")
        print(
            f"Stats: codec={metadata.get('codec', 'ecomp')} original={original_size}B "
            f"compressed={compressed_size}B ratio={ratio:.3f}x"
        )

    return 0


def _cmd_decompress(args: argparse.Namespace) -> int:
    archive_path = Path(args.archive).expanduser().resolve()
    if not archive_path.exists():
        raise SystemExit(f"Archive not found: {archive_path}")

    metadata_path = (
        Path(args.metadata_path).expanduser().resolve()
        if args.metadata_path
        else derive_metadata_path(archive_path)
    )
    metadata = read_metadata(metadata_path)
    payload = read_payload(archive_path)

    frame = decompress_alignment(payload, metadata)
    if not args.no_checksum:
        _verify_checksum(frame.sequences, metadata)

    alignment_output = (
        Path(args.alignment_output).expanduser().resolve()
        if args.alignment_output
        else archive_path.with_suffix(f".{args.alignment_format}")
    )
    alignment_output.parent.mkdir(parents=True, exist_ok=True)
    write_alignment(frame, alignment_output, fmt=args.alignment_format)
    print(f"Wrote alignment to {alignment_output}")
    return 0


def _cmd_inspect(args: argparse.Namespace) -> int:
    archive_path = Path(args.archive).expanduser().resolve()
    metadata_path = (
        Path(args.metadata_path).expanduser().resolve()
        if args.metadata_path
        else derive_metadata_path(archive_path)
    )
    metadata = read_metadata(metadata_path)

    if args.summary:
        codec = metadata.get("codec", "ecomp")
        num_sequences = metadata.get("num_sequences")
        alignment_length = metadata.get("alignment_length")
        payload_encoding = metadata.get("payload_encoding", "raw")
        print(f"Codec: {codec}")
        print(f"Sequences: {num_sequences}")
        print(f"Alignment columns: {alignment_length}")
        print(f"Payload encoding: {payload_encoding}")
        return 0

    json.dump(metadata, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


def _verify_checksum(sequences: Iterable[str], metadata: dict[str, object]) -> None:
    expected = metadata.get("checksum_sha256")
    if not expected:
        return
    observed = alignment_checksum(sequences)
    if observed != expected:
        raise SystemExit(
            "Checksum mismatch after decompression: "
            f"expected {expected}, observed {observed}"
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    handler = getattr(args, "handler", None)
    if handler is None:
        parser.error("No command handler registered")
    return handler(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
