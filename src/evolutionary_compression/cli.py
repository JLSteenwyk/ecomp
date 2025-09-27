"""Unified command-line interface for evolutionary compression workflows."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

from .compression.pipeline import compress_alignment, decompress_alignment
from .compression.phylo_bundle import (
    compress_alignment_with_tree,
    decompress_alignment_with_tree,
)
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
BUNDLE_SUFFIX = ".ecbt"
CODEC_CHOICES = ("auto", "alignment", "phylo")


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
        "tree",
        nargs="?",
        help="Optional Newick tree; enables phylogenetic bundle compression",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Destination archive path (default: alignment stem + .ecomp or .ecbt)",
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
        "--codec",
        choices=CODEC_CHOICES,
        default="auto",
        help="Compression strategy: auto-select, alignment-only, or phylo bundle",
    )
    parser.add_argument(
        "--bundle-suffix",
        default=BUNDLE_SUFFIX,
        help=f"Suffix when writing phylo bundles (default: {BUNDLE_SUFFIX})",
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
        "-t",
        "--tree-output",
        dest="tree_output",
        help="Tree output path when decompressing a phylo bundle (default: archive stem + .tree)",
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
        description="Evolutionary compression toolkit for alignments and phylogenies",
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

    tree_path = Path(args.tree).expanduser().resolve() if args.tree else None
    codec_mode = args.codec
    if codec_mode == "auto":
        codec_mode = "phylo" if tree_path is not None else "alignment"

    if codec_mode == "phylo" and tree_path is None:
        raise SystemExit("Phylo codec requested but no tree file provided")

    if tree_path and not tree_path.exists():
        raise SystemExit(f"Tree file not found: {tree_path}")

    frame = read_alignment(alignment_path, fmt=args.alignment_format)

    if codec_mode == "phylo":
        newick = tree_path.read_text()
        payload, metadata = compress_alignment_with_tree(frame, newick)
        suffix = _normalize_suffix(args.bundle_suffix)
    else:
        compressed = compress_alignment(frame)
        payload = compressed.payload
        metadata = compressed.metadata
        suffix = ALIGNMENT_SUFFIX
        if tree_path is not None:
            print(
                "[warning] Tree provided but alignment codec was requested; tree will be ignored",
                file=sys.stderr,
            )

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
        original_size = _total_input_size(alignment_path, tree_path)
        compressed_size = output_path.stat().st_size
        ratio = (original_size / compressed_size) if compressed_size else float("inf")
        label = "alignment+tree" if codec_mode == "phylo" else "alignment"
        print(
            f"Stats: codec={metadata.get('codec', 'ecomp')} dataset={label} "
            f"original={original_size}B compressed={compressed_size}B ratio={ratio:.3f}x"
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

    codec = metadata.get("codec", "ecomp")
    if codec == "phylo-bundle":
        frame, newick = decompress_alignment_with_tree(payload, metadata)
        if not args.no_checksum:
            _verify_checksum(frame.sequences, metadata)

        alignment_output = (
            Path(args.alignment_output).expanduser().resolve()
            if args.alignment_output
            else archive_path.with_suffix(f".{args.alignment_format}")
        )
        alignment_output.parent.mkdir(parents=True, exist_ok=True)
        write_alignment(frame, alignment_output, fmt=args.alignment_format)

        tree_output = (
            Path(args.tree_output).expanduser().resolve()
            if args.tree_output
            else archive_path.with_suffix(".tree")
        )
        tree_output.parent.mkdir(parents=True, exist_ok=True)
        tree_output.write_text(newick)

        print(f"Wrote alignment to {alignment_output}")
        print(f"Wrote tree to {tree_output}")
        return 0

    # default: standard alignment archive handled by helper
    from . import decompress_file as _legacy_decompress_file  # local import to avoid cycle

    output_path = _legacy_decompress_file(
        archive_path,
        output_path=args.alignment_output,
        metadata_path=metadata_path,
        output_format=args.alignment_format,
        validate_checksum=not args.no_checksum,
    )
    print(f"Wrote alignment to {output_path}")
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
        if codec == "phylo-bundle":
            leaf_count = len(metadata.get("sequence_ids", []))
            print(f"Bundle leaf count: {leaf_count}")
        return 0

    json.dump(metadata, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _normalize_suffix(suffix: str) -> str:
    suffix = suffix.strip()
    if not suffix:
        return BUNDLE_SUFFIX
    if not suffix.startswith("."):
        suffix = f".{suffix}"
    return suffix


def _total_input_size(alignment_path: Path, tree_path: Path | None) -> int:
    total = alignment_path.stat().st_size
    if tree_path is not None:
        total += tree_path.stat().st_size
    return total


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
