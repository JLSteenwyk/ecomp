"""Unified command-line interface for evolutionary compression workflows."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

from ._compat import zip_strict
from .compression.pipeline import compress_alignment, decompress_alignment
from .config import DEFAULT_OUTPUT_FORMAT
from .diagnostics.checksums import alignment_checksum
from .diagnostics.metrics import (
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
from .phylo import infer_distance_tree_from_frame, tree_to_newick
from .io import read_alignment, write_alignment
from .storage import derive_metadata_path, read_archive, write_archive, write_metadata

ALIGNMENT_SUFFIX = ".ecomp"

ASCII_BANNER = (
    "        _____                      "\
    "\n       / ____|                     "\
    "\n   ___| |     ___  _ __ ___  _ __  "\
    "\n  / _ \\ |    / _ \\| '_ ` _ \\| '_ \\ "\
    "\n |  __/ |___| (_) | | | | | | |_) |"\
    "\n  \\___|\\_____\\___/|_| |_| |_| .__/ "\
    "\n                            | |    "\
    "\n                            |_|    "
)

COMMAND_CATEGORIES = [
    "Compression",
    "Diagnostics",
    "Phylogenetics",
    "Utilities",
]

COMMAND_REGISTRY: list[dict[str, Any]] = []
_COMMAND_COUNTER = 0


def _attach_help_banner(parser: argparse.ArgumentParser) -> None:
    original = parser.format_help

    def _wrapped() -> str:
        return f"{ASCII_BANNER}\n\n{original()}"

    parser.format_help = _wrapped


def _register_command(
    *,
    category: str,
    name: str,
    aliases: list[str] | None,
    help_text: str,
) -> None:
    global _COMMAND_COUNTER
    COMMAND_REGISTRY.append(
        {
            "category": category,
            "name": name,
            "aliases": aliases or [],
            "help": help_text,
            "order": _COMMAND_COUNTER,
        }
    )
    _COMMAND_COUNTER += 1


# ---------------------------------------------------------------------------
# Parser construction (ClipKIT/PhyKIT style: subcommands with handler binding)
# ---------------------------------------------------------------------------


def _resolve_archive_args(
    archive: str,
    metadata: str | None,
) -> tuple[Path, Path | None]:
    archive_path = Path(archive).expanduser().resolve()
    if not archive_path.exists():
        raise SystemExit(f"Archive not found: {archive_path}")

    metadata_path = (
        Path(metadata).expanduser().resolve()
        if metadata
        else None
    )
    return archive_path, metadata_path


def _load_alignment_from_archive(
    archive_path: Path,
    metadata_path: Path | None,
    *,
    validate_checksum: bool = True,
):
    try:
        payload, metadata, _ = read_archive(archive_path, metadata_path=metadata_path)
    except FileNotFoundError as exc:
        fallback = metadata_path or derive_metadata_path(archive_path)
        raise SystemExit(
            "Metadata sidecar missing for legacy archive; specify it with --metadata "
            f"or restore {fallback}"
        ) from exc

    frame = decompress_alignment(payload, metadata)
    if validate_checksum:
        _verify_checksum(frame.sequences, metadata)
    return frame, metadata


def _add_archive_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("archive", help="Compressed archive produced by `ecomp zip`")
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_path",
        help="Metadata JSON path for legacy archives (default: alongside archive)",
    )
    parser.add_argument(
        "--no-checksum",
        action="store_true",
        help="Skip checksum validation during analysis",
    )


def _add_consensus_sequence_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Emit a majority-rule consensus sequence in FASTA format"
    parser = subparsers.add_parser(
        "consensus_sequence",
        aliases=["con_seq"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="consensus_sequence",
        aliases=["con_seq"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.add_argument(
        "--header",
        default="consensus",
        help="FASTA header to use for the emitted consensus sequence",
    )
    parser.set_defaults(handler=_cmd_consensus_sequence)


def _add_column_base_counts_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Report per-column residue counts as JSON"
    parser = subparsers.add_parser(
        "column_base_counts",
        aliases=["col_counts"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="column_base_counts",
        aliases=["col_counts"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.add_argument(
        "--include-gaps",
        action="store_true",
        help="Include gap symbols in the reported tallies",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="Indent width for the JSON output (default: 2)",
    )
    parser.set_defaults(handler=_cmd_column_base_counts)


def _add_gap_fraction_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Print per-column gap fractions"
    parser = subparsers.add_parser(
        "gap_fraction",
        aliases=["gap_frac"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="gap_fraction",
        aliases=["gap_frac"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_gap_fraction)


def _add_shannon_entropy_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Print per-column Shannon entropy"
    parser = subparsers.add_parser(
        "shannon_entropy",
        aliases=["entropy"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="shannon_entropy",
        aliases=["entropy"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_shannon_entropy)


def _add_parsimony_informative_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Identify parsimony-informative columns"
    parser = subparsers.add_parser(
        "parsimony_informative_sites",
        aliases=["parsimony"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="parsimony_informative_sites",
        aliases=["parsimony"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_parsimony_informative)


def _add_constant_columns_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Flag columns that are constant after removing gaps"
    parser = subparsers.add_parser(
        "constant_columns",
        aliases=["const_cols"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="constant_columns",
        aliases=["const_cols"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_constant_columns)


def _add_pairwise_identity_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Report the pairwise identity matrix"
    parser = subparsers.add_parser(
        "pairwise_identity",
        aliases=["pid"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="pairwise_identity",
        aliases=["pid"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_pairwise_identity)


def _add_alignment_length_no_gaps_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Count alignment columns that contain at least one non-gap"
    parser = subparsers.add_parser(
        "alignment_length_excluding_gaps",
        aliases=["len_no_gaps"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="alignment_length_excluding_gaps",
        aliases=["len_no_gaps"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_alignment_length_no_gaps)


def _add_alignment_length_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Report total alignment columns including gaps"
    parser = subparsers.add_parser(
        "alignment_length",
        aliases=["len_total"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="alignment_length",
        aliases=["len_total"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_alignment_length)


def _add_variable_sites_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Count sites with more than one residue ignoring gaps"
    parser = subparsers.add_parser(
        "variable_sites",
        aliases=["var_sites"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="variable_sites",
        aliases=["var_sites"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_variable_sites)


def _add_percentage_identity_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Compute the mean pairwise identity across sequences"
    parser = subparsers.add_parser(
        "percentage_identity",
        aliases=["pct_id"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="percentage_identity",
        aliases=["pct_id"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_percentage_identity)


def _add_rcv_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Compute the relative composition variability (RCV)"
    parser = subparsers.add_parser(
        "relative_composition_variability",
        aliases=["rcv"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Diagnostics",
        name="relative_composition_variability",
        aliases=["rcv"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.set_defaults(handler=_cmd_rcv)


def _add_distance_tree_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Infer a distance-based phylogenetic tree and emit Newick"
    parser = subparsers.add_parser(
        "distance_tree",
        aliases=["dist_tree"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Phylogenetics",
        name="distance_tree",
        aliases=["dist_tree"],
        help_text=help_text,
    )
    _add_archive_options(parser)
    parser.add_argument(
        "--method",
        choices=["nj", "upgma"],
        default="nj",
        help="Tree construction method to use (default: nj)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Write the Newick tree to this path instead of stdout",
    )
    parser.set_defaults(handler=_cmd_distance_tree)


def _add_zip_arguments(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    help_text = "Zip an alignment (optionally bundling a companion Newick tree)"
    parser = subparsers.add_parser(
        "zip",
        aliases=["compress"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Compression",
        name="zip",
        aliases=["compress"],
        help_text=help_text,
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
        help="Optional path to write a JSON metadata copy",
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
    parser.set_defaults(handler=_cmd_zip)


def _add_unzip_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Restore data from an evolutionary compression archive"
    parser = subparsers.add_parser(
        "unzip",
        aliases=["decompress"],
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Compression",
        name="unzip",
        aliases=["decompress"],
        help_text=help_text,
    )
    parser.add_argument("archive", help="Compressed archive produced by `ecomp zip`")
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_path",
        help="Metadata JSON path for legacy archives (default: alongside archive)",
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
    parser.set_defaults(handler=_cmd_unzip)


def _add_inspect_arguments(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    help_text = "Display metadata for an evolutionary compression archive"
    parser = subparsers.add_parser(
        "inspect",
        help=help_text,
        description=help_text,
    )
    _attach_help_banner(parser)
    _register_command(
        category="Utilities",
        name="inspect",
        aliases=[],
        help_text=help_text,
    )
    parser.add_argument("archive", help="Archive to inspect")
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_path",
        help="Metadata JSON path for legacy archives (default: alongside archive)",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print a human-readable summary instead of raw JSON",
    )
    parser.set_defaults(handler=_cmd_inspect)


def _format_main_help(parser: argparse.ArgumentParser) -> str:
    usage = f"usage: {parser.prog} <command> [options]"
    description = parser.description or ""

    lines = [ASCII_BANNER, "", usage, "", description.rstrip(), ""]

    width = max(
        (len(_command_display(entry)) for entry in COMMAND_REGISTRY),
        default=0,
    )
    width = max(width, 20)

    for category in COMMAND_CATEGORIES:
        entries = [entry for entry in COMMAND_REGISTRY if entry["category"] == category]
        if not entries:
            continue
        lines.append(f"{category}:")
        for entry in sorted(entries, key=lambda item: item["order"]):
            display = _command_display(entry)
            help_text = entry["help"]
            lines.append(f"  {display.ljust(width)}  {help_text}")
        lines.append("")

    option_label = "-h, --help"
    lines.extend(
        [
            "General options:",
            f"  {option_label.ljust(width)}  Show this message and exit",
            "",
            "Run 'ecomp <command> --help' for details on a specific command.",
        ]
    )

    return "\n".join(line for line in lines if line is not None)


def _command_display(entry: dict[str, Any]) -> str:
    if entry["aliases"]:
        alias_text = ", ".join(entry["aliases"])
        return f"{entry['name']} ({alias_text})"
    return entry["name"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ecomp",
        description="Evolutionary compression toolkit for multiple sequence alignments",
        add_help=False,
    )
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS)
    COMMAND_REGISTRY.clear()
    global _COMMAND_COUNTER
    _COMMAND_COUNTER = 0
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_consensus_sequence_arguments(subparsers)
    _add_column_base_counts_arguments(subparsers)
    _add_gap_fraction_arguments(subparsers)
    _add_shannon_entropy_arguments(subparsers)
    _add_parsimony_informative_arguments(subparsers)
    _add_constant_columns_arguments(subparsers)
    _add_pairwise_identity_arguments(subparsers)
    _add_alignment_length_no_gaps_arguments(subparsers)
    _add_alignment_length_arguments(subparsers)
    _add_variable_sites_arguments(subparsers)
    _add_percentage_identity_arguments(subparsers)
    _add_rcv_arguments(subparsers)
    _add_distance_tree_arguments(subparsers)
    _add_zip_arguments(subparsers)
    _add_unzip_arguments(subparsers)
    _add_inspect_arguments(subparsers)
    parser.format_help = lambda: _format_main_help(parser)
    return parser


# ---------------------------------------------------------------------------
# Command handlers
# ---------------------------------------------------------------------------


def _cmd_zip(args: argparse.Namespace) -> int:
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

    output_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else alignment_path.with_suffix(suffix)
    )
    metadata_path = (
        Path(args.metadata_path).expanduser().resolve()
        if args.metadata_path
        else None
    )

    write_archive(output_path, payload, metadata)
    print(f"Created {output_path} (metadata embedded)")

    if metadata_path is not None:
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        write_metadata(metadata_path, metadata)
        print(f"Metadata copy written to {metadata_path}")

    if args.stats:
        original_size = alignment_path.stat().st_size
        compressed_size = output_path.stat().st_size
        ratio = (original_size / compressed_size) if compressed_size else float("inf")
        print(
            f"Stats: codec={metadata.get('codec', 'ecomp')} original={original_size}B "
            f"compressed={compressed_size}B ratio={ratio:.3f}x"
        )

    return 0


def _cmd_unzip(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, metadata = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )

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
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    try:
        _, metadata, _ = read_archive(archive_path, metadata_path=metadata_path)
    except FileNotFoundError as exc:
        fallback = metadata_path or derive_metadata_path(archive_path)
        raise SystemExit(
            "Metadata sidecar missing for legacy archive; specify it with --metadata "
            f"or restore {fallback}"
        ) from exc

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


def _cmd_consensus_sequence(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    consensus = majority_rule_consensus(frame)
    print(f">{args.header}")
    print(consensus)
    return 0


def _cmd_column_base_counts(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    counts = column_base_counts(frame, include_gaps=args.include_gaps)
    payload = [
        {"column": index + 1, "counts": dict(counter)}
        for index, counter in enumerate(counts)
    ]
    json.dump(payload, sys.stdout, indent=args.indent)
    sys.stdout.write("\n")
    return 0


def _cmd_gap_fraction(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    fractions = column_gap_fraction(frame)
    for index, fraction in enumerate(fractions, start=1):
        print(f"{index}\t{fraction:.6f}")
    return 0


def _cmd_shannon_entropy(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    entropies = column_shannon_entropy(frame)
    for index, entropy in enumerate(entropies, start=1):
        print(f"{index}\t{entropy:.6f}")
    return 0


def _cmd_parsimony_informative(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    mask = parsimony_informative_columns(frame)
    indices = [index + 1 for index, value in enumerate(mask) if value]
    print(f"total\t{len(indices)}")
    if indices:
        print("indices\t" + " ".join(str(i) for i in indices))
    return 0


def _cmd_constant_columns(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    mask = constant_columns(frame)
    indices = [index + 1 for index, value in enumerate(mask) if value]
    print(f"total\t{len(indices)}")
    if indices:
        print("indices\t" + " ".join(str(i) for i in indices))
    return 0


def _cmd_pairwise_identity(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    result = pairwise_identity_matrix(frame)
    ids = frame.ids
    header = "\t".join(["id", *ids])
    print(header)
    for seq_id, row in zip_strict(ids, result.matrix):
        formatted = [
            "nan" if math.isnan(value) else f"{value:.6f}"
            for value in row
        ]
        print("\t".join([seq_id, *formatted]))

    print("# coverage")
    print(header)
    for seq_id, row in zip_strict(ids, result.coverage):
        print("\t".join([seq_id, *[str(value) for value in row]]))
    return 0


def _cmd_alignment_length_no_gaps(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    length = alignment_length_excluding_gaps(frame)
    print(length)
    return 0


def _cmd_alignment_length(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    length = alignment_length(frame)
    print(length)
    return 0


def _cmd_variable_sites(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    count = variable_site_count(frame)
    print(count)
    return 0


def _cmd_percentage_identity(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    value = percentage_identity(frame)
    if math.isnan(value):
        print("nan")
    else:
        print(f"{value:.6f}")
    return 0


def _cmd_rcv(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    value = relative_composition_variability(frame)
    print(f"{value:.6f}")
    return 0


def _cmd_distance_tree(args: argparse.Namespace) -> int:
    archive_path, metadata_path = _resolve_archive_args(args.archive, args.metadata_path)
    frame, _ = _load_alignment_from_archive(
        archive_path, metadata_path, validate_checksum=not args.no_checksum
    )
    tree = infer_distance_tree_from_frame(frame, method=args.method)
    newick = tree_to_newick(tree)
    if args.output:
        output_path = Path(args.output).expanduser().resolve()
        output_path.write_text(newick + ("\n" if not newick.endswith("\n") else ""))
        print(f"Wrote Newick tree to {output_path}")
    else:
        print(newick)
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
