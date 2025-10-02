#!/usr/bin/env python3
"""Benchmark eComp metrics versus equivalent PhyKIT commands.

This script measures wall-clock time for a set of alignment diagnostics when
executed directly on a compressed ``.ecomp`` archive (via the ``ecomp`` CLI) and
when executed on a decompressed FASTA with ``phykit``.  It reports average and
best timings so you can quantify the speed-up gained by staying in the compressed
representation.

Examples
========

Benchmark the default operation subset (5 repeats, 1 warm-up) on an archive::

    python scripts/benchmark_metrics.py data/fixtures/small_phylo.ecomp

Benchmark a custom subset and export JSON::

    python scripts/benchmark_metrics.py archive.ecomp \
        --operations consensus shannon_entropy percentage_identity \
        --repeat 10 --json

By default the script decompresses the archive once into a temporary directory
for the PhyKIT runs.  You can provide an explicit FASTA via ``--alignment`` if
you already have one on disk.

Requirements
============

* ``ecomp`` CLI available on the PATH (provided by this project).
* ``phykit`` CLI installed and available on the PATH for comparison runs.
* Python 3.11+.

"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import shutil

# Mapping of benchmark operations: name -> (ecomp subcommand, phykit command)
DEFAULT_OPERATIONS: dict[str, tuple[str, str]] = {
    "consensus": ("consensus_sequence", "consensus"),
    "gap_fraction": ("gap_fraction", "gap_fraction"),
    "shannon_entropy": ("shannon_entropy", "alignment_entropy"),
    "parsimony_informative_sites": ("parsimony_informative_sites", "parsimony_informative_sites"),
    "constant_columns": ("constant_columns", "constant_sites"),
    "alignment_length_excluding_gaps": (
        "alignment_length_excluding_gaps",
        "alignment_length_excluding_gaps",
    ),
    "percentage_identity": ("percentage_identity", "matrix_identity"),
    "relative_composition_variability": (
        "relative_composition_variability",
        "relative_composition_variability",
    ),
}


@dataclass
class TimingResult:
    operation: str
    tool: str
    mean: float
    best: float
    repetitions: int
    success: bool
    message: str | None = None


def _run_command(cmd: list[str]) -> None:
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _time_command(cmd: list[str], warmup: int, repeat: int) -> tuple[list[float], str | None]:
    durations: list[float] = []
    try:
        for _ in range(max(0, warmup)):
            _run_command(cmd)
        for _ in range(repeat):
            start = time.perf_counter()
            _run_command(cmd)
            durations.append(time.perf_counter() - start)
    except subprocess.CalledProcessError as exc:  # pragma: no cover - passthrough
        return [], f"command failed with exit code {exc.returncode}"
    return durations, None


def _benchmark_operation(
    operation: str,
    ecomp_cmd: list[str],
    phykit_cmd: list[str],
    warmup: int,
    repeat: int,
    phykit_available: bool,
) -> list[TimingResult]:
    results: list[TimingResult] = []

    ecomp_times, err = _time_command(ecomp_cmd, warmup, repeat)
    if ecomp_times:
        results.append(
            TimingResult(
                operation=operation,
                tool="ecomp",
                mean=statistics.mean(ecomp_times),
                best=min(ecomp_times),
                repetitions=len(ecomp_times),
                success=True,
            )
        )
    else:
        results.append(
            TimingResult(
                operation=operation,
                tool="ecomp",
                mean=float("nan"),
                best=float("nan"),
                repetitions=0,
                success=False,
                message=err or "command failed",
            )
        )

    if phykit_available:
        phykit_times, err = _time_command(phykit_cmd, warmup, repeat)
        if phykit_times:
            results.append(
                TimingResult(
                    operation=operation,
                    tool="phykit",
                    mean=statistics.mean(phykit_times),
                    best=min(phykit_times),
                    repetitions=len(phykit_times),
                    success=True,
                )
            )
        else:
            results.append(
                TimingResult(
                    operation=operation,
                    tool="phykit",
                    mean=float("nan"),
                    best=float("nan"),
                    repetitions=0,
                    success=False,
                    message=err or "command failed",
                )
            )
    return results


def _format_seconds(value: float) -> str:
    if value != value:  # NaN
        return "n/a"
    if value < 1e-3:
        return f"{value * 1e6:6.1f} µs"
    if value < 1:
        return f"{value * 1e3:6.1f} ms"
    return f"{value:6.3f} s"


def _print_table(results: Iterable[TimingResult]) -> None:
    rows: list[TimingResult] = list(results)
    if not rows:
        return
    print("\nBenchmark results (mean ± best):")
    header = f"{'operation':30} {'tool':10} {'mean':>12} {'best':>12}"
    print(header)
    print("-" * len(header))
    for row in rows:
        mean_str = _format_seconds(row.mean)
        best_str = _format_seconds(row.best)
        label = row.message if row.message else ""
        print(f"{row.operation:30} {row.tool:10} {mean_str:>12} {best_str:>12} {label}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Path to a .ecomp archive")
    parser.add_argument(
        "--alignment",
        type=Path,
        help="Existing alignment file to reuse for PhyKIT runs (FASTA/PHYLIP).",
    )
    parser.add_argument(
        "--operations",
        nargs="+",
        choices=sorted(DEFAULT_OPERATIONS.keys()),
        help="Subset of operations to benchmark (default: all).",
    )
    parser.add_argument("--repeat", type=int, default=5, help="Number of timed repetitions")
    parser.add_argument("--warmup", type=int, default=1, help="Number of warm-up runs before timing")
    parser.add_argument("--json", type=Path, help="Optional path to emit benchmark results as JSON")
    parser.add_argument("--csv", type=Path, help="Optional path to emit benchmark results as CSV")
    args = parser.parse_args()

    archive_path: Path = args.archive.expanduser().resolve()
    if not archive_path.exists():
        raise SystemExit(f"Archive not found: {archive_path}")

    ecomp_exe = shutil.which("ecomp")
    if not ecomp_exe:
        raise SystemExit("`ecomp` CLI not found on PATH")

    phykit_exe = shutil.which("phykit")
    if not phykit_exe:
        print("[warn] `phykit` CLI not found on PATH; skipping PhyKIT comparisons")

    operations = args.operations or sorted(DEFAULT_OPERATIONS.keys())

    temp_dir: tempfile.TemporaryDirectory[str] | None = None
    alignment_path = args.alignment
    if phykit_exe and alignment_path is None:
        temp_dir = tempfile.TemporaryDirectory(prefix="ecomp-benchmark-")
        alignment_path = Path(temp_dir.name) / "alignment.fasta"
        unzip_cmd = [
            ecomp_exe,
            "unzip",
            str(archive_path),
            "--alignment-output",
            str(alignment_path),
            "--no-checksum",
        ]
        print(f"[info] Decompressing archive once for PhyKIT: {' '.join(unzip_cmd)}")
        try:
            _run_command(unzip_cmd)
        except subprocess.CalledProcessError as exc:  # pragma: no cover - passthrough
            raise SystemExit(f"Failed to decompress archive for PhyKIT: {exc}") from exc

    results: list[TimingResult] = []
    for op in operations:
        ecomp_cmd, phykit_cmd = DEFAULT_OPERATIONS[op]
        ecomp_args = [ecomp_exe, ecomp_cmd, str(archive_path)]
        phykit_args: list[str] = []
        if phykit_exe and alignment_path is not None:
            phykit_args = [phykit_exe, phykit_cmd, str(alignment_path)]

        results.extend(
            _benchmark_operation(
                operation=op,
                ecomp_cmd=ecomp_args,
                phykit_cmd=phykit_args,
                warmup=args.warmup,
                repeat=args.repeat,
                phykit_available=bool(phykit_exe and phykit_args),
            )
        )

    _print_table(results)

    if args.json:
        payload = [
            {
                "operation": row.operation,
                "tool": row.tool,
                "mean_seconds": row.mean,
                "best_seconds": row.best,
                "repetitions": row.repetitions,
                "success": row.success,
                "message": row.message,
            }
            for row in results
        ]
        args.json.write_text(json.dumps(payload, indent=2))
        print(f"[info] JSON results written to {args.json}")

    if args.csv:
        header = [
            "operation",
            "tool",
            "mean_seconds",
            "best_seconds",
            "repetitions",
            "success",
            "message",
        ]
        with args.csv.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            for row in results:
                mean_val = "" if math.isnan(row.mean) else f"{row.mean:.6f}"
                best_val = "" if math.isnan(row.best) else f"{row.best:.6f}"
                writer.writerow(
                    [
                        row.operation,
                        row.tool,
                        mean_val,
                        best_val,
                        row.repetitions,
                        row.success,
                        row.message or "",
                    ]
                )
        print(f"[info] CSV results written to {args.csv}")

    if temp_dir is not None:
        temp_dir.cleanup()


if __name__ == "__main__":  # pragma: no cover - script entry point
    main()
