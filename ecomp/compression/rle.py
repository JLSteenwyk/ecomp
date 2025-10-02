"""Run-length encoding of consensus column profiles."""

from __future__ import annotations

from .._compat import dataclass
from typing import Dict, Iterable, Iterator, Tuple

from .consensus import ColumnProfile


@dataclass(slots=True)
class RunLengthBlock:
    """A run of identical consensus columns encoded with bitmasks."""

    consensus: str
    bitmask: bytes
    residues: bytes
    run_length: int


MAX_RUN_LENGTH = 255


def _column_signature(
    column: ColumnProfile,
    num_sequences: int,
    symbol_lookup: Dict[str, int],
    bits_per_symbol: int,
) -> Tuple[bytes, bytes]:
    bitmask_length = (num_sequences + 7) // 8
    bitmask = bytearray(bitmask_length)
    residues: list[int] = []
    for seq_idx, residue in column.deviations:
        byte_index = seq_idx // 8
        bit_index = seq_idx % 8
        bitmask[byte_index] |= 1 << bit_index
        try:
            code = symbol_lookup[residue]
        except KeyError as exc:
            raise ValueError(f"Residue {residue!r} not present in alphabet mapping") from exc
        residues.append(code)
    return bytes(bitmask), _pack_codes(residues, bits_per_symbol)


def _pack_codes(codes: Iterable[int], bits_per_symbol: int) -> bytes:
    buffer = 0
    bits_in_buffer = 0
    output = bytearray()
    for code in codes:
        buffer = (buffer << bits_per_symbol) | code
        bits_in_buffer += bits_per_symbol
        while bits_in_buffer >= 8:
            bits_in_buffer -= 8
            output.append((buffer >> bits_in_buffer) & 0xFF)
            buffer &= (1 << bits_in_buffer) - 1
    if bits_in_buffer:
        output.append(buffer << (8 - bits_in_buffer))
    return bytes(output)


def iter_run_length_blocks(
    columns: Iterable[ColumnProfile],
    num_sequences: int,
    symbol_lookup: Dict[str, int],
    bits_per_symbol: int,
) -> Iterator[RunLengthBlock]:
    """Group identical columns into run-length blocks using bitmask signatures."""

    previous_key: tuple[str, bytes, bytes] | None = None
    run_length = 0
    previous_payload: tuple[str, bytes, bytes] | None = None
    for column in columns:
        bitmask, residues = _column_signature(
            column, num_sequences, symbol_lookup, bits_per_symbol
        )
        key = (column.consensus, bitmask, residues)
        if key == previous_key and run_length < MAX_RUN_LENGTH:
            run_length += 1
        else:
            if previous_payload is not None:
                consensus, prev_mask, prev_residues = previous_payload
                yield RunLengthBlock(
                    consensus=consensus,
                    bitmask=prev_mask,
                    residues=prev_residues,
                    run_length=run_length,
                )
            previous_key = key
            previous_payload = key
            run_length = 1
    if previous_payload is not None:
        consensus, prev_mask, prev_residues = previous_payload
        yield RunLengthBlock(
            consensus=consensus,
            bitmask=prev_mask,
            residues=prev_residues,
            run_length=run_length,
        )


def collect_run_length_blocks(
    columns: Iterable[ColumnProfile],
    num_sequences: int,
    symbol_lookup: Dict[str, int],
    bits_per_symbol: int,
) -> list[RunLengthBlock]:
    """Return run-length encoded blocks as a list."""

    return list(
        iter_run_length_blocks(columns, num_sequences, symbol_lookup, bits_per_symbol)
    )
