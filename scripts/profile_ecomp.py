#!/usr/bin/env python3
"""Profile eComp compression stages for one or more alignments."""

from __future__ import annotations

import argparse
import json
import lzma
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

from evolutionary_compression import read_alignment
from evolutionary_compression.compression import pipeline
from evolutionary_compression.compression.consensus import collect_column_profiles
from evolutionary_compression.compression.encoding import encode_blocks
from evolutionary_compression.compression.rle import collect_run_length_blocks
from evolutionary_compression.storage import write_metadata


try:  # pragma: no cover - optional dependency
    import zstandard as zstd  # noqa: F401

    HAS_ZSTD = True
except ModuleNotFoundError:  # pragma: no cover
    HAS_ZSTD = False


@dataclass
class ProfileResult:
    path: str
    num_sequences: int
    alignment_length: int
    reorder_seconds: float
    column_seconds: float
    rle_seconds: float
    encode_seconds: float
    compress_seconds: float
    decompress_seconds: float
    seq_id_bytes: int
    residue_payload_bytes: int
    raw_payload_bytes: int
    zstd_bytes: int | None
    xz_bytes: int
    zlib_bytes: int
    final_payload_bytes: int
    metadata_bytes: int
    payload_encoding: str
    permutation_applied: bool
    ordering_strategy: str


def _metadata_encoded_size(metadata: dict[str, Any]) -> int:
    with NamedTemporaryFile(prefix="ecomp_profile_", suffix=".json", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        write_metadata(tmp_path, metadata)
        return tmp_path.stat().st_size
    finally:
        try:
            tmp_path.unlink()
        except FileNotFoundError:  # pragma: no cover - best effort cleanup
            pass


def profile_alignment(path: Path) -> ProfileResult:
    frame = read_alignment(path)
    original_frame = frame

    start = time.perf_counter()
    frame, permutation, order_label = pipeline._compute_similarity_order(frame)  # type: ignore[attr-defined]
    reorder_seconds = time.perf_counter() - start
    permutation_changed = permutation != list(range(frame.num_sequences))

    start = time.perf_counter()
    column_profiles = collect_column_profiles(frame)
    column_seconds = time.perf_counter() - start

    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, (len(alphabet) - 1).bit_length() or 1)

    start = time.perf_counter()
    run_length_blocks = collect_run_length_blocks(
        column_profiles, frame.num_sequences, symbol_lookup, bits_per_symbol
    )
    rle_seconds = time.perf_counter() - start

    bitmask_bytes = (frame.num_sequences + 7) // 8

    start = time.perf_counter()
    run_length_payload = encode_blocks(
        run_length_blocks, bitmask_bytes, bits_per_symbol, alphabet
    )
    encode_seconds = time.perf_counter() - start

    seq_id_block = pipeline._encode_sequence_ids(frame.ids)  # type: ignore[attr-defined]
    raw_payload = seq_id_block + run_length_payload

    payload_candidates: list[tuple[str, bytes, float]] = [("raw", raw_payload, 0.0)]

    start = time.perf_counter()
    zlib_payload = pipeline.zlib.compress(raw_payload, level=9)  # type: ignore[attr-defined]
    zlib_seconds = time.perf_counter() - start
    payload_candidates.append(("zlib", zlib_payload, zlib_seconds))

    zstd_bytes: int | None = None
    if pipeline._ZSTD_COMPRESSOR is not None:  # type: ignore[attr-defined]
        start = time.perf_counter()
        zstd_payload = pipeline._ZSTD_COMPRESSOR.compress(raw_payload)  # type: ignore[attr-defined]
        zstd_seconds = time.perf_counter() - start
        payload_candidates.append(("zstd", zstd_payload, zstd_seconds))
        zstd_bytes = len(zstd_payload)
    else:
        zstd_seconds = 0.0

    start = time.perf_counter()
    xz_payload = lzma.compress(raw_payload, preset=6)
    xz_seconds = time.perf_counter() - start
    payload_candidates.append(("xz", xz_payload, xz_seconds))

    payload_encoding, payload_bytes, compress_seconds = min(
        payload_candidates, key=lambda item: len(item[1])
    )

    metadata = {
        "format_version": pipeline.FORMAT_VERSION,  # type: ignore[attr-defined]
        "codec": "ecomp",
        "num_sequences": frame.num_sequences,
        "alignment_length": frame.alignment_length,
        "alphabet": alphabet,
        "source_format": frame.metadata.get("source_format", "unknown"),
        "checksum_sha256": pipeline.alignment_checksum(original_frame.sequences),  # type: ignore[attr-defined]
        "run_length_blocks": len(run_length_blocks),
        "max_run_length": max((blk.run_length for blk in run_length_blocks), default=0),
        "columns_with_deviations": sum(1 for blk in column_profiles if blk.deviations),
        "bitmask_bytes": bitmask_bytes,
        "bits_per_symbol": bits_per_symbol,
        "payload_encoding": payload_encoding,
        "payload_encoded_bytes": len(payload_bytes),
        "payload_raw_bytes": len(raw_payload),
        "sequence_id_codec": "inline",
    }
    if permutation_changed:
        metadata["sequence_permutation"] = permutation

    final_payload, metadata = pipeline._maybe_use_gzip_fallback(  # type: ignore[attr-defined]
        original_frame, payload_bytes, metadata
    )

    metadata_bytes = _metadata_encoded_size(metadata)

    start = time.perf_counter()
    restored = pipeline.decompress_alignment(final_payload, metadata)
    decompress_seconds = time.perf_counter() - start
    assert restored.sequences == original_frame.sequences

    return ProfileResult(
        path=str(path),
        num_sequences=frame.num_sequences,
        alignment_length=frame.alignment_length,
        reorder_seconds=reorder_seconds,
        column_seconds=column_seconds,
        rle_seconds=rle_seconds,
        encode_seconds=encode_seconds,
        compress_seconds=compress_seconds,
        decompress_seconds=decompress_seconds,
        seq_id_bytes=len(seq_id_block),
        residue_payload_bytes=len(run_length_payload),
        raw_payload_bytes=len(raw_payload),
        zstd_bytes=zstd_bytes,
        zlib_bytes=len(zlib_payload),
        xz_bytes=len(xz_payload),
        final_payload_bytes=len(final_payload),
        metadata_bytes=metadata_bytes,
        payload_encoding=metadata.get("payload_encoding", "raw"),
        permutation_applied=permutation_changed,
        ordering_strategy=metadata.get("ordering_strategy", order_label),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Profile eComp compression stages")
    parser.add_argument("inputs", nargs="+", help="Alignment files to profile")
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional JSON output path; defaults to stdout",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results = [profile_alignment(Path(p)) for p in args.inputs]
    payload = [asdict(result) for result in results]
    if args.output:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    else:
        print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
