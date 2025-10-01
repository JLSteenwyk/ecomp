#!/usr/bin/env python3
"""Regenerate benchmark ratios for the EVOCOMP manuscript corpus."""

from __future__ import annotations

import argparse
import bz2
import contextlib
import csv
import gzip
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from evolutionary_compression.compression.pipeline import (
    _tree_guided_order,
    compress_alignment,
)
from evolutionary_compression.io import alignment_from_sequences, read_alignment
from evolutionary_compression.storage import write_archive


@dataclass
class Measurement:
    ratio_ecomp: float
    ratio_tree: float | None
    ratio_gzip: float
    ratio_bzip2: float
    ordering_label: str
    tree_order_label: str | None


@contextlib.contextmanager
def override_env(name: str, value: str) -> Iterable[None]:
    original = os.environ.get(name)
    os.environ[name] = value
    try:
        yield
    finally:
        if original is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = original


def _format_float(value: float) -> str:
    return f"{value:.15f}".rstrip("0").rstrip(".")


def _gzip_ratio(data: bytes) -> float:
    compressed = gzip.compress(data, compresslevel=6)
    return len(data) / len(compressed)


def _bzip2_ratio(data: bytes) -> float:
    compressed = bz2.compress(data, compresslevel=9)
    return len(data) / len(compressed)


def _write_archive_size(payload: bytes, metadata: dict[str, object], path: Path) -> int:
    write_archive(path, payload, metadata)
    size = path.stat().st_size
    path.unlink()
    return size


def measure_alignment(alignment_path: Path, tree_path: Path | None, temp_dir: Path) -> Measurement:
    frame = read_alignment(alignment_path)
    alignment_bytes = alignment_path.read_bytes()
    metadata = dict(frame.metadata)
    tree_text = None
    if tree_path and tree_path.exists():
        tree_text = tree_path.read_text()
        metadata["tree_newick"] = tree_text

    alphabet = list(frame.alphabet)
    ids = list(frame.ids)
    sequences = list(frame.sequences)

    auto_frame = alignment_from_sequences(ids, sequences, alphabet=alphabet, metadata=dict(metadata))
    auto_result = compress_alignment(auto_frame)
    auto_size = _write_archive_size(
        auto_result.payload, auto_result.metadata, temp_dir / f"{alignment_path.stem}.auto.ecomp"
    )
    auto_ratio = len(alignment_bytes) / auto_size
    ordering = str(auto_result.metadata.get("ordering_strategy", "unknown"))

    tree_ratio = None
    tree_label = None
    if tree_text is not None:
        guide_frame = alignment_from_sequences(ids, sequences, alphabet=alphabet, metadata=dict(metadata))
        order = _tree_guided_order(guide_frame)
        if order:
            ordered_ids = [ids[idx] for idx in order]
            ordered_sequences = [sequences[idx] for idx in order]
            tree_frame = alignment_from_sequences(
                ordered_ids,
                ordered_sequences,
                alphabet=alphabet,
                metadata=dict(metadata),
            )
            with override_env("ECOMP_SEQUENCE_ORDER", "baseline"):
                tree_result = compress_alignment(tree_frame)
            tree_size = _write_archive_size(
                tree_result.payload,
                tree_result.metadata,
                temp_dir / f"{alignment_path.stem}.tree.ecomp",
            )
            tree_ratio = len(alignment_bytes) / tree_size
            tree_label = "tree"

    gzip_ratio = _gzip_ratio(alignment_bytes)
    bzip2_ratio = _bzip2_ratio(alignment_bytes)

    return Measurement(
        ratio_ecomp=auto_ratio,
        ratio_tree=tree_ratio,
        ratio_gzip=gzip_ratio,
        ratio_bzip2=bzip2_ratio,
        ordering_label=ordering,
        tree_order_label=tree_label,
    )


def update_csv(data_root: Path, csv_path: Path) -> None:
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise SystemExit("Benchmark CSV is empty; nothing to update")

    required_columns = {
        "dataset",
        "file",
        "ratio_ecomp",
        "ratio_gzip",
        "gain_vs_gzip",
        "ratio_bzip2",
        "gain_vs_bzip2",
        "ratio_tree",
        "gain_tree_vs_ecomp",
        "tree_used",
        "ordering",
        "tree_order",
    }
    missing = required_columns - set(rows[0].keys())
    if missing:
        raise SystemExit(f"CSV missing required columns: {sorted(missing)}")

    with tempfile.TemporaryDirectory() as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        for row in rows:
            dataset = row["dataset"].strip()
            filename = row["file"].strip()
            alignment_path = data_root / dataset / filename
            if not alignment_path.exists():
                raise FileNotFoundError(f"Alignment not found: {alignment_path}")
            tree_path = alignment_path.parent / f"{alignment_path.name}.tre"
            measurement = measure_alignment(alignment_path, tree_path, tmpdir)

            row["ratio_ecomp"] = _format_float(measurement.ratio_ecomp)
            row["ratio_gzip"] = _format_float(measurement.ratio_gzip)
            row["gain_vs_gzip"] = _format_float(
                measurement.ratio_ecomp - measurement.ratio_gzip
            )
            row["ratio_bzip2"] = _format_float(measurement.ratio_bzip2)
            row["gain_vs_bzip2"] = _format_float(
                measurement.ratio_ecomp - measurement.ratio_bzip2
            )

            if measurement.ratio_tree is not None:
                row["ratio_tree"] = _format_float(measurement.ratio_tree)
                row["gain_tree_vs_ecomp"] = _format_float(
                    measurement.ratio_tree - measurement.ratio_ecomp
                )
                row["tree_used"] = "True"
                row["tree_order"] = measurement.tree_order_label or "tree"
            else:
                row["ratio_tree"] = row["ratio_ecomp"]
                row["gain_tree_vs_ecomp"] = _format_float(0.0)
                row["tree_used"] = "False"
                row["tree_order"] = "n/a"

            row["ordering"] = measurement.ordering_label

    fieldnames = rows[0].keys()
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Recompute benchmark ratios for data_large")
    parser.add_argument(
        "--data-root",
        type=Path,
        required=True,
        help="Path to the data_large directory from the manuscript archive",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        required=True,
        help="Path to data_large_benchmark_analysis.csv to update",
    )
    args = parser.parse_args()
    if args.data_root.is_file():
        parser.error("--data-root must be a directory, not a file")
    if args.csv.is_dir():
        args.csv = args.csv / "data_large_benchmark_analysis.csv"
    if not args.csv.exists():
        parser.error(f"CSV not found: {args.csv}")
    return args


def main() -> None:
    args = parse_args()
    update_csv(args.data_root, args.csv)


if __name__ == "__main__":
    main()
