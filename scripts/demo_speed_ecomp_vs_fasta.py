#!/usr/bin/env python3
"""Demonstrate speed differences between eComp archives and raw FASTA files.

For a given alignment this script measures how long it takes to compute a
selection of diagnostics

  1. Directly on the ``.ecomp`` archive (via the ``ecomp`` CLI).
  2. On the decompressed FASTA using ``phykit``.

It prints a summary table and writes a CSV file (`speed_comparison.csv`) for
further analysis.
"""
from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

# Mapping operation -> (ecomp command, phykit command)
OPERATIONS = {
    "consensus": ("consensus_sequence", "consensus"),
    "shannon_entropy": ("shannon_entropy", "alignment_entropy"),
    "variable_sites": ("variable_sites", "variable_sites"),
    "percentage_identity": ("percentage_identity", "matrix_identity"),
}


def _detect(exe: str) -> str:
    path = shutil.which(exe)
    if not path:
        raise SystemExit(f"Required executable `{exe}` not found on PATH")
    return path


def _time_command(cmd: list[str], repeat: int = 5) -> float:
    durations = []
    for _ in range(repeat):
        start = time.perf_counter()
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        durations.append(time.perf_counter() - start)
    return sum(durations) / len(durations)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Path to the `.ecomp` archive")
    parser.add_argument("--alignment", type=Path, help="Optional FASTA to reuse for PhyKIT")
    parser.add_argument("--repeat", type=int, default=5, help="Number of repetitions per command")
    parser.add_argument("--csv", type=Path, default=Path("speed_comparison.csv"))
    args = parser.parse_args()

    archive = args.archive.expanduser().resolve()
    if not archive.exists():
        raise SystemExit(f"Archive not found: {archive}")

    ecomp_exe = _detect("ecomp")
    phykit_exe = _detect("phykit")

    temp_dir = None
    fasta_path = args.alignment
    if fasta_path is None:
        temp_dir = tempfile.TemporaryDirectory(prefix="ecomp-demo-")
        fasta_path = Path(temp_dir.name) / "alignment.fasta"
        decompress_cmd = [
            ecomp_exe,
            "unzip",
            str(archive),
            "--alignment-output",
            str(fasta_path),
        ]
        print(f"[info] Decompressing once for phykit: {' '.join(decompress_cmd)}")
        subprocess.run(decompress_cmd, check=True)

    results = []
    for op, (ecomp_cmd, phykit_cmd) in OPERATIONS.items():
        ecomp_time = _time_command([ecomp_exe, ecomp_cmd, str(archive)], repeat=args.repeat)
        phykit_time = _time_command([phykit_exe, phykit_cmd, str(fasta_path)], repeat=args.repeat)
        speedup = math.inf if phykit_time == 0 else phykit_time / ecomp_time
        results.append((op, ecomp_time, phykit_time, speedup))

    print("\nOperation                eComp (s)   FASTA (s)   Speedup")
    print("-------------------------------------------------------")
    for op, e_time, f_time, speedup in results:
        speed_str = "inf" if math.isinf(speedup) else f"{speedup:7.2f}x"
        print(f"{op:24} {e_time:10.4f} {f_time:10.4f} {speed_str}")

    if args.csv:
        with args.csv.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["operation", "ecomp_seconds", "fasta_seconds", "speedup"])
            for row in results:
                writer.writerow(row)
        print(f"[info] CSV written to {args.csv}")

    if temp_dir is not None:
        temp_dir.cleanup()


if __name__ == "__main__":
    main()
