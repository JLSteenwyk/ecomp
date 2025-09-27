"""Binary encoding helpers for run-length blocks."""

from __future__ import annotations

import struct
from typing import Sequence

from collections import Counter

from .rle import RunLengthBlock

BLOCK_HEADER_STRUCT = ">BBB"
BLOCK_HEADER_SIZE = struct.calcsize(BLOCK_HEADER_STRUCT)


class EncodingError(RuntimeError):
    """Raised when binary encoding encounters an invalid condition."""


class DecodingError(RuntimeError):
    """Raised when binary decoding encounters an invalid condition."""


def encode_blocks(
    blocks: Sequence[RunLengthBlock], bitmask_bytes: int, bits_per_symbol: int
) -> bytes:
    """Serialize run-length blocks into a binary payload using dictionary coding."""

    dictionary, dictionary_map = _build_dictionary(
        blocks, bitmask_bytes, bits_per_symbol
    )

    payload = bytearray()
    payload.append(len(dictionary))
    for consensus, trimmed_mask, residues in dictionary:
        payload.append(ord(consensus))
        mask_len = len(trimmed_mask)
        payload.append(mask_len)
        payload.extend(trimmed_mask)
        payload.extend(struct.pack(">H", len(residues)))
        payload.extend(residues)

    payload.extend(struct.pack(">I", len(blocks)))
    for block in blocks:
        if not (1 <= block.run_length <= 255):
            raise EncodingError("Run length must be within 1..255 for encoding")
        key = (block.consensus, block.bitmask, block.residues)
        dict_id = dictionary_map.get(key)
        if dict_id is not None:
            payload.append(1)  # dictionary reference marker
            payload.append(dict_id)
            payload.append(block.run_length)
        else:
            payload.append(0)  # literal
            payload.append(block.run_length)
            payload.append(ord(block.consensus))
            trimmed_mask, mask_len = _trim_bitmask(block.bitmask)
            payload.append(mask_len)
            payload.extend(trimmed_mask)
            payload.extend(struct.pack(">H", len(block.residues)))
            payload.extend(block.residues)
    return bytes(payload)


def decode_blocks(
    payload: bytes, bitmask_bytes: int, bits_per_symbol: int
) -> list[RunLengthBlock]:
    """Parse binary payload into run-length blocks."""

    blocks: list[RunLengthBlock] = []
    cursor = 0
    payload_length = len(payload)
    if cursor >= payload_length:
        return blocks

    dict_count = payload[cursor]
    cursor += 1
    dictionary: list[tuple[str, bytes, bytes]] = []
    for _ in range(dict_count):
        if cursor + 2 > payload_length:
            raise DecodingError("Dictionary entry truncated")
        consensus_value = payload[cursor]
        cursor += 1
        mask_len = payload[cursor]
        cursor += 1
        if cursor + mask_len > payload_length:
            raise DecodingError("Dictionary mask truncated")
        mask_slice = payload[cursor : cursor + mask_len]
        cursor += mask_len
        if cursor + 2 > payload_length:
            raise DecodingError("Dictionary residue length truncated")
        residues_len = struct.unpack(">H", payload[cursor : cursor + 2])[0]
        cursor += 2
        if cursor + residues_len > payload_length:
            raise DecodingError("Dictionary residues truncated")
        residues = payload[cursor : cursor + residues_len]
        cursor += residues_len
        bitmask = mask_slice + bytes(bitmask_bytes - mask_len)
        dictionary.append((bytes([consensus_value]).decode("ascii"), bitmask, residues))

    if cursor + 4 > payload_length:
        raise DecodingError("Missing block count")
    (block_count,) = struct.unpack(">I", payload[cursor : cursor + 4])
    cursor += 4

    for _ in range(block_count):
        if cursor >= payload_length:
            raise DecodingError("Block data truncated")
        marker = payload[cursor]
        cursor += 1
        if marker == 1:
            if cursor + 2 > payload_length:
                raise DecodingError("Dictionary block truncated")
            dict_id = payload[cursor]
            cursor += 1
            run_length = payload[cursor]
            cursor += 1
            try:
                consensus, bitmask, residues = dictionary[dict_id]
            except IndexError as exc:  # pragma: no cover - guard
                raise DecodingError("Dictionary index out of range") from exc
            blocks.append(
                RunLengthBlock(
                    consensus=consensus,
                    bitmask=bitmask,
                    residues=residues,
                    run_length=run_length,
                )
            )
        elif marker == 0:
            if cursor + 3 > payload_length:
                raise DecodingError("Literal block truncated")
            run_length = payload[cursor]
            cursor += 1
            consensus_value = payload[cursor]
            cursor += 1
            mask_len = payload[cursor]
            cursor += 1
            if cursor + mask_len > payload_length:
                raise DecodingError("Literal mask truncated")
            mask_slice = payload[cursor : cursor + mask_len]
            cursor += mask_len
            if cursor + 2 > payload_length:
                raise DecodingError("Literal residue length truncated")
            residues_len = struct.unpack(">H", payload[cursor : cursor + 2])[0]
            cursor += 2
            if cursor + residues_len > payload_length:
                raise DecodingError("Literal residues truncated")
            residues = payload[cursor : cursor + residues_len]
            cursor += residues_len
            bitmask = mask_slice + bytes(bitmask_bytes - mask_len)
            blocks.append(
                RunLengthBlock(
                    consensus=bytes([consensus_value]).decode("ascii"),
                    bitmask=bitmask,
                    residues=residues,
                    run_length=run_length,
                )
            )
        else:
            raise DecodingError(f"Unknown block marker {marker}")
    return blocks


def _encode_char(char: str) -> bytes:
    if len(char) != 1:
        raise EncodingError(f"Expected single-character residue, received {char!r}")
    try:
        return char.encode("ascii")
    except UnicodeEncodeError as exc:  # pragma: no cover - guard for unexpected alphabets
        raise EncodingError(f"Non-ASCII residue encountered: {char!r}") from exc


def _popcount(data: bytes) -> int:
    return sum(bin(byte).count("1") for byte in data)


def _trim_bitmask(bitmask: bytes) -> tuple[bytes, int]:
    length = len(bitmask)
    while length > 0 and bitmask[length - 1] == 0:
        length -= 1
    return bitmask[:length], length


def _build_dictionary(
    blocks: Sequence[RunLengthBlock], bitmask_bytes: int, bits_per_symbol: int
) -> tuple[list[tuple[str, bytes, bytes]], dict[tuple[str, bytes, bytes], int]]:
    counter: Counter[tuple[str, bytes, bytes]] = Counter()
    for block in blocks:
        key = (block.consensus, block.bitmask, block.residues)
        counter[key] += block.run_length

    dictionary: list[tuple[str, bytes, bytes]] = []
    dictionary_map: dict[tuple[str, bytes, bytes], int] = {}
    for key, freq in counter.most_common(255):
        consensus, bitmask, residues = key
        trimmed_mask, mask_len = _trim_bitmask(bitmask)
        literal_size = 1 + 1 + 1 + 1 + mask_len + 2 + len(residues)
        entry_size = 1 + 1 + mask_len + 2 + len(residues)
        reference_size = 1 + 1 + 1
        if freq * literal_size <= entry_size + freq * reference_size:
            continue
        dictionary_map[key] = len(dictionary)
        dictionary.append((consensus, trimmed_mask, residues))
    return dictionary, dictionary_map
