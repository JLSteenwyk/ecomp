#!/usr/bin/env python3
"""Benchmark helper comparing eComp to standard compressors."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

from evolutionary_compression import read_alignment
from evolutionary_compression.compression.phylo_bundle import (
    compress_alignment_with_tree,
    decompress_alignment_with_tree,
)
from evolutionary_compression.storage import (
    read_metadata,
    read_payload,
    write_metadata,
    write_payload,
)

TREE_SUFFIXES = (".tree", ".nwk", ".newick")

REQUIRED_EXTERNAL = {
    "gzip": {"compress": ["gzip", "-k"], "decompress": ["gzip", "-dk"]},
    "bzip2": {"compress": ["bzip2", "-k"], "decompress": ["bzip2", "-dk"]},
    "xz": {"compress": ["xz", "-k"], "decompress": ["xz", "-dk"]},
}

@dataclass
class BenchmarkResult:
    codec: str
    original_size: int
    compressed_size: int
    compression_ratio: float
    compress_seconds: float
    decompress_seconds: float
    notes: str = ""

def run_command(cmd: list[str], cwd: Path | None = None) -> float:
    start = time.perf_counter()
    env = os.environ.copy()
    src_path = Path(__file__).resolve().parents[1] / "src"
    existing = env.get("PYTHONPATH")
    paths = [str(src_path)]
    venv_site = Path(__file__).resolve().parents[1] / "venv" / "lib"
    if venv_site.exists():
        for site_candidate in venv_site.iterdir():
            candidate = site_candidate / "site-packages"
            if candidate.exists():
                paths.append(str(candidate))
    if existing:
        paths.append(existing)
    env["PYTHONPATH"] = ":".join(paths)
    subprocess.run(
        cmd,
        cwd=cwd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )
    return time.perf_counter() - start

def benchmark_ecomp(
    input_path: Path, workdir: Path
) -> BenchmarkResult:
    ecomp_path = workdir / (input_path.name + ".ecomp")
    metadata_path = ecomp_path.with_suffix(".json")
    restored_path = workdir / (input_path.stem + ".restored.fasta")

    compress_seconds = run_command(
        [
            sys.executable,
            "-m",
            "evolutionary_compression.cli",
            "compress",
            str(input_path),
            "-o",
            str(ecomp_path),
            "-m",
            str(metadata_path),
        ]
    )
    decompress_seconds = run_command(
        [
            sys.executable,
            "-m",
            "evolutionary_compression.cli",
            "decompress",
            str(ecomp_path),
            "-m",
            str(metadata_path),
            "-o",
            str(restored_path),
        ]
    )

    original_size = input_path.stat().st_size
    compressed_size = ecomp_path.stat().st_size
    ratio = original_size / compressed_size if compressed_size else float("inf")

    return BenchmarkResult(
        codec="ecomp",
        original_size=original_size,
        compressed_size=compressed_size,
        compression_ratio=ratio,
        compress_seconds=compress_seconds,
        decompress_seconds=decompress_seconds,
    )

def _find_tree_path(input_path: Path) -> Path | None:
    for suffix in TREE_SUFFIXES:
        direct = input_path.with_suffix(suffix)
        if direct.exists():
            return direct
        stacked = input_path.with_suffix(input_path.suffix + suffix)
        if stacked.exists():
            return stacked

    # fall back to glob search in same directory
    stem = input_path.stem
    for suffix in TREE_SUFFIXES:
        for candidate in input_path.parent.glob(f"{stem}*{suffix}"):
            if candidate.exists():
                return candidate
    return None

def benchmark_phylo_bundle(input_path: Path, workdir: Path) -> BenchmarkResult:
    tree_path = _find_tree_path(input_path)
    if tree_path is None:
        original_size = input_path.stat().st_size
        return BenchmarkResult(
            codec="phylo-bundle",
            original_size=original_size,
            compressed_size=0,
            compression_ratio=0.0,
            compress_seconds=0.0,
            decompress_seconds=0.0,
            notes="tree_not_found",
        )

    bundle_path = workdir / (input_path.stem + ".ecbt")
    metadata_path = bundle_path.with_suffix(".json")

    frame = read_alignment(input_path)
    tree_text = tree_path.read_text()

    start = time.perf_counter()
    payload, metadata = compress_alignment_with_tree(frame, tree_text)
    compress_seconds = time.perf_counter() - start

    write_payload(bundle_path, payload)
    write_metadata(metadata_path, metadata)

    start = time.perf_counter()
    restored_frame, _restored_newick = decompress_alignment_with_tree(payload, metadata)
    decompress_seconds = time.perf_counter() - start

    original_size = input_path.stat().st_size + tree_path.stat().st_size
    compressed_size = bundle_path.stat().st_size
    ratio = original_size / compressed_size if compressed_size else 0.0
    note = f"tree={tree_path.name}"

    if restored_frame.sequences != frame.sequences:
        note = "mismatch"

    return BenchmarkResult(
        codec="phylo-bundle",
        original_size=original_size,
        compressed_size=compressed_size,
        compression_ratio=ratio,
        compress_seconds=compress_seconds,
        decompress_seconds=decompress_seconds,
        notes=note,
    )

def benchmark_external(codec: str, input_path: Path, workdir: Path) -> BenchmarkResult:
    mapping = REQUIRED_EXTERNAL[codec]
    binary = mapping["compress"][0]
    if shutil.which(binary) is None:
        return BenchmarkResult(
            codec=codec,
            original_size=input_path.stat().st_size,
            compressed_size=0,
            compression_ratio=0.0,
            compress_seconds=0.0,
            decompress_seconds=0.0,
            notes=f"{binary} not available on PATH",
        )

    target_path = workdir / (input_path.name + f".{codec}")
    temp_input = workdir / input_path.name
    shutil.copy2(input_path, temp_input)

    compress_cmd = mapping["compress"] + [str(temp_input)]
    compress_seconds = run_command(compress_cmd, cwd=workdir)

    if not target_path.exists():
        # some tools append their own suffix
        fallback = {
            "gzip": temp_input.with_suffix(temp_input.suffix + ".gz"),
            "bzip2": temp_input.with_suffix(temp_input.suffix + ".bz2"),
            "xz": temp_input.with_suffix(temp_input.suffix + ".xz"),
        }
        target_path = fallback.get(codec, target_path)

    if temp_input.exists():
        temp_input.unlink()

    restored_path = workdir / (input_path.name + ".restored")
    decompress_cmd = mapping["decompress"] + [str(target_path)]
    decompress_seconds = run_command(decompress_cmd, cwd=workdir)
    if not restored_path.exists():
        possible = target_path.with_suffix("")
        if possible.exists():
            restored_path = possible

    original_size = input_path.stat().st_size
    compressed_size = target_path.stat().st_size if target_path.exists() else 0
    ratio = original_size / compressed_size if compressed_size else 0.0

    return BenchmarkResult(
        codec=codec,
        original_size=original_size,
        compressed_size=compressed_size,
        compression_ratio=ratio,
        compress_seconds=compress_seconds,
        decompress_seconds=decompress_seconds,
    )

def run_benchmarks(input_paths: Iterable[Path], codecs: list[str]) -> list[BenchmarkResult]:
    results: list[BenchmarkResult] = []
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)
        for input_path in input_paths:
            for codec in codecs:
                if codec == "ecomp":
                    results.append(benchmark_ecomp(input_path, workdir))
                elif codec == "phylo-bundle":
                    results.append(benchmark_phylo_bundle(input_path, workdir))
                elif codec in REQUIRED_EXTERNAL:
                    results.append(benchmark_external(codec, input_path, workdir))
                else:
                    results.append(
                        BenchmarkResult(
                            codec=codec,
                            original_size=input_path.stat().st_size,
                            compressed_size=0,
                            compression_ratio=0.0,
                            compress_seconds=0.0,
                            decompress_seconds=0.0,
                            notes="Codec not recognized",
                        )
                    )
    return results

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare eComp against other compressors")
    parser.add_argument("inputs", nargs="+", help="Alignment files to benchmark")
    parser.add_argument(
        "--codecs",
        nargs="+",
        default=["ecomp", "phylo-bundle", "gzip", "bzip2", "xz"],
        help="Codec names to benchmark",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional JSON output path to persist benchmark results",
    )
    return parser.parse_args(argv)

def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_paths = [Path(path).resolve() for path in args.inputs]
    for path in input_paths:
        if not path.exists():
            raise SystemExit(f"Input path not found: {path}")

    results = run_benchmarks(input_paths, args.codecs)
    payload = [asdict(result) for result in results]

    if args.output:
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    else:
        json.dump(payload, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
    return 0

if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
