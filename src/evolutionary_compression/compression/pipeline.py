"""High-level compression and decompression orchestration."""

from __future__ import annotations

import gzip
import io
import math
import zlib
from dataclasses import dataclass
from typing import Any, Iterable, List

from ..diagnostics.checksums import alignment_checksum
from ..io import AlignmentFrame, alignment_from_sequences
from ..config import FORMAT_VERSION
from .consensus import collect_column_profiles
from .encoding import decode_blocks, encode_blocks
from .rle import RunLengthBlock, collect_run_length_blocks


@dataclass(slots=True)
class CompressedAlignment:
    """Payload plus metadata produced by the compression pipeline."""

    payload: bytes
    metadata: dict[str, Any]


def compress_alignment(frame: AlignmentFrame) -> CompressedAlignment:
    """Compress an alignment into a binary payload and structured metadata."""

    column_profiles = collect_column_profiles(frame)
    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, math.ceil(math.log2(max(len(alphabet), 1))))

    run_length_blocks = collect_run_length_blocks(
        column_profiles, frame.num_sequences, symbol_lookup, bits_per_symbol
    )
    bitmask_bytes = (frame.num_sequences + 7) // 8
    raw_payload = encode_blocks(run_length_blocks, bitmask_bytes, bits_per_symbol)
    compressed_payload = zlib.compress(raw_payload, level=9)
    if len(compressed_payload) < len(raw_payload):
        payload_bytes = compressed_payload
        payload_encoding = "zlib"
    else:
        payload_bytes = raw_payload
        payload_encoding = "raw"
    max_run_length = max((block.run_length for block in run_length_blocks), default=0)
    deviation_columns = sum(1 for block in column_profiles if block.deviations)
    metadata = {
        "format_version": FORMAT_VERSION,
        "codec": "ecomp",
        "num_sequences": frame.num_sequences,
        "alignment_length": frame.alignment_length,
        "sequence_ids": frame.ids,
        "alphabet": frame.alphabet,
        "source_format": frame.metadata.get("source_format", "unknown"),
        "checksum_sha256": alignment_checksum(frame.sequences),
        "run_length_blocks": len(run_length_blocks),
        "max_run_length": max_run_length,
        "columns_with_deviations": deviation_columns,
        "bitmask_bytes": bitmask_bytes,
        "bits_per_symbol": bits_per_symbol,
        "payload_encoding": payload_encoding,
        "payload_encoded_bytes": len(payload_bytes),
        "payload_raw_bytes": len(raw_payload),
    }
    payload_bytes, metadata = _maybe_use_gzip_fallback(frame, payload_bytes, metadata)
    return CompressedAlignment(payload=payload_bytes, metadata=metadata)


def decompress_alignment(payload: bytes, metadata: dict[str, Any]) -> AlignmentFrame:
    """Reconstruct an :class:`AlignmentFrame` from payload and metadata."""

    fallback_info = metadata.get("fallback")
    if fallback_info:
        return _decompress_fallback(payload, metadata)

    expected_columns = metadata["alignment_length"]
    num_sequences = metadata["num_sequences"]
    sequence_ids = metadata["sequence_ids"]
    alphabet = metadata.get("alphabet", [])
    bitmask_bytes = metadata.get("bitmask_bytes")
    bits_per_symbol = metadata.get("bits_per_symbol")

    if len(sequence_ids) != num_sequences:
        raise ValueError(
            "Metadata sequence count does not match sequence IDs provided"
        )
    if not isinstance(bitmask_bytes, int) or bitmask_bytes <= 0:
        raise ValueError("Metadata missing valid 'bitmask_bytes' entry")
    if not isinstance(bits_per_symbol, int) or bits_per_symbol <= 0:
        if alphabet:
            bits_per_symbol = max(1, math.ceil(math.log2(len(alphabet))))
        else:
            raise ValueError("Metadata missing valid 'bits_per_symbol' entry")

    payload_encoding = metadata.get("payload_encoding", "raw")
    if payload_encoding == "zlib":
        decoded_payload = zlib.decompress(payload)
    elif payload_encoding in {"raw", None}:
        decoded_payload = payload
    else:
        raise ValueError(f"Unsupported payload encoding: {payload_encoding}")

    blocks = decode_blocks(decoded_payload, bitmask_bytes, bits_per_symbol)
    sequences = [["" for _ in range(expected_columns)] for _ in range(num_sequences)]

    try:
        symbol_table = list(alphabet)
    except TypeError as exc:  # pragma: no cover - guard against malformed metadata
        raise ValueError("Alphabet metadata is not iterable") from exc

    column_index = 0
    for block in blocks:
        consensus = block.consensus
        for _ in range(block.run_length):
            if column_index >= expected_columns:
                raise ValueError(
                    "Decoded columns exceed expected alignment length"
                )
            residue_indices = _iter_deviation_indices(block.bitmask, num_sequences)
            residues = _decode_residues(
                block.residues,
                len(residue_indices),
                bits_per_symbol,
                symbol_table,
            )
            for seq_list in sequences:
                seq_list[column_index] = consensus
            for seq_index, residue in zip(residue_indices, residues, strict=True):
                sequences[seq_index][column_index] = residue
            column_index += 1

    if column_index != expected_columns:
        raise ValueError(
            f"Decoded columns ({column_index}) do not match expected length {expected_columns}"
        )

    reconstructed = ["".join(row) for row in sequences]
    return alignment_from_sequences(
        ids=sequence_ids,
        sequences=reconstructed,
        alphabet=alphabet,
        metadata={"source_format": metadata.get("source_format", "unknown")},
    )


def _iter_deviation_indices(bitmask: bytes, num_sequences: int) -> List[int]:
    indices: list[int] = []
    for seq_index in range(num_sequences):
        byte_index = seq_index // 8
        bit_index = seq_index % 8
        if bitmask[byte_index] & (1 << bit_index):
            indices.append(seq_index)
    return indices


def _decode_residues(
    data: bytes,
    count: int,
    bits_per_symbol: int,
    alphabet: Iterable[str],
) -> list[str]:
    if count == 0:
        return []
    mask = (1 << bits_per_symbol) - 1
    values: list[int] = []
    buffer = 0
    bits_in_buffer = 0
    data_iter = iter(data)
    while len(values) < count:
        while bits_in_buffer < bits_per_symbol:
            try:
                byte = next(data_iter)
            except StopIteration as exc:  # pragma: no cover - corruption guard
                raise ValueError("Insufficient residue data during decode") from exc
            buffer = (buffer << 8) | byte
            bits_in_buffer += 8
        shift = bits_in_buffer - bits_per_symbol
        value = (buffer >> shift) & mask
        buffer &= (1 << shift) - 1
        bits_in_buffer -= bits_per_symbol
        values.append(value)
    alphabet_list = list(alphabet)
    try:
        return [alphabet_list[value] for value in values]
    except IndexError as exc:  # pragma: no cover - corruption guard
        raise ValueError("Residue code exceeds alphabet size") from exc


def _alignment_to_fasta_bytes(frame: AlignmentFrame) -> bytes:
    buffer = io.StringIO()
    for seq_id, sequence in zip(frame.ids, frame.sequences, strict=True):
        buffer.write(f">{seq_id}\n{sequence}\n")
    return buffer.getvalue().encode("utf-8")


def _parse_fasta_bytes(data: bytes) -> AlignmentFrame:
    ids: list[str] = []
    sequences: list[str] = []
    current_id: str | None = None
    current_seq: list[str] = []
    for raw_line in data.decode("utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequences.append("".join(current_seq))
            current_id = line[1:]
            ids.append(current_id)
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        sequences.append("".join(current_seq))
    return alignment_from_sequences(ids=ids, sequences=sequences)


def _maybe_use_gzip_fallback(
    frame: AlignmentFrame, payload: bytes, metadata: dict[str, Any]
) -> tuple[bytes, dict[str, Any]]:
    fasta_bytes = _alignment_to_fasta_bytes(frame)
    gzip_payload = gzip.compress(fasta_bytes)
    if len(gzip_payload) + 1 < len(payload):
        updated = dict(metadata)
        updated["fallback"] = {
            "type": "gzip",
            "format": metadata.get("source_format", "fasta"),
        }
        updated["payload_encoding"] = "gzip"
        updated["payload_encoded_bytes"] = len(gzip_payload)
        updated["payload_raw_bytes"] = len(fasta_bytes)
        return gzip_payload, updated
    return payload, metadata


def _decompress_fallback(payload: bytes, metadata: dict[str, Any]) -> AlignmentFrame:
    info = metadata.get("fallback", {})
    fallback_type = info.get("type")
    if fallback_type == "gzip":
        data = gzip.decompress(payload)
        frame = _parse_fasta_bytes(data)
        frame.metadata["source_format"] = info.get("format", "fasta")
        return frame
    raise ValueError(f"Unsupported fallback type: {fallback_type}")
