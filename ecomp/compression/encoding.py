"""Binary encoding helpers for run-length blocks."""

from __future__ import annotations

import heapq
import math
import struct
from collections import Counter
from typing import Any, Sequence

from .rle import RunLengthBlock
from .._compat import zip_strict

BLOCK_HEADER_STRUCT = ">BBB"
BLOCK_HEADER_SIZE = struct.calcsize(BLOCK_HEADER_STRUCT)


class EncodingError(RuntimeError):
    """Raised when binary encoding encounters an invalid condition."""


class DecodingError(RuntimeError):
    """Raised when binary decoding encounters an invalid condition."""


def encode_blocks(
    blocks: Sequence[RunLengthBlock],
    bitmask_bytes: int,
    bits_per_symbol: int,
    alphabet: Sequence[str],
) -> bytes:
    """Serialize run-length blocks using dictionary coding and consensus-aware residue packing."""

    alphabet_list = list(alphabet)
    consensus_sets: dict[str, set[str]] = {}
    consensus_freqs: dict[str, Counter[str]] = {}
    block_residue_chars: list[list[str]] = []
    encoded_bitmasks: list[tuple[int, int, bytes]] = []

    for block in blocks:
        deviation_count = _popcount(block.bitmask)
        if deviation_count:
            codes = _unpack_codes(block.residues, deviation_count, bits_per_symbol)
            chars = [alphabet_list[code] for code in codes]
        else:
            chars = []
        consensus_sets.setdefault(block.consensus, set()).update(chars)
        if chars:
            freq_counter = consensus_freqs.setdefault(block.consensus, Counter())
            freq_counter.update(chars)
        block_residue_chars.append(chars)

        encoded_bitmasks.append(_encode_bitmask(block.bitmask, bitmask_bytes))

    consensus_models: dict[str, dict[str, Any]] = {}
    for consensus, residues_set in consensus_sets.items():
        if not residues_set:
            continue
        residues = sorted(residues_set)
        freq_counter = consensus_freqs.get(consensus, Counter())
        total = sum(freq_counter[res] for res in residues)
        unique = len(residues)

        mode = 0
        width = max(1, math.ceil(math.log2(unique))) if unique > 0 else 1
        encode_map = {residue: idx for idx, residue in enumerate(residues)}
        lengths_list = None
        decode_map = None
        max_len = None

        if total > 0 and unique > 1:
            huffman_lengths, huffman_bits = _build_huffman_lengths(residues, freq_counter)
            if huffman_lengths is not None:
                fixed_bits = width * total + 8  # width stored as one byte
                huffman_total_bits = huffman_bits + 8 * len(residues)
                if huffman_total_bits < fixed_bits:
                    mode = 1
                    encode_map, decode_map, max_len = _canonical_code_maps(residues, huffman_lengths)
                    lengths_list = huffman_lengths

        consensus_models[consensus] = {
            "mode": mode,
            "residues": residues,
            "width": max(1, width),
            "encode_map": encode_map,
            "lengths": lengths_list,
            "decode_map": decode_map,
            "max_code_len": max_len,
        }

    prepared_blocks: list[RunLengthBlock] = []
    for block, chars in zip_strict(blocks, block_residue_chars):
        if not chars:
            encoded_residues = b""
        else:
            model = consensus_models[block.consensus]
            if model["mode"] == 0:
                local_codes = [model["encode_map"][char] for char in chars]
                encoded_residues = _pack_codes(local_codes, model["width"])
            else:
                encoded_residues = _pack_huffman_codes(model["encode_map"], chars)
        prepared_blocks.append(
            RunLengthBlock(
                consensus=block.consensus,
                bitmask=block.bitmask,
                residues=encoded_residues,
                run_length=block.run_length,
            )
        )

    dictionary_entries, dictionary_map = _build_dictionary(
        prepared_blocks, encoded_bitmasks
    )

    payload = bytearray()
    payload.append(len(consensus_models))
    for consensus in sorted(consensus_models):
        model = consensus_models[consensus]
        residues = model["residues"]
        payload.append(ord(consensus))
        payload.append(model["mode"])
        payload.append(len(residues))
        for residue in residues:
            payload.extend(residue.encode("ascii"))
        if model["mode"] == 0:
            payload.append(model["width"])
        else:
            for length in model["lengths"]:
                payload.append(length)

    payload.append(len(dictionary_entries))
    for entry in dictionary_entries:
        payload.append(ord(entry["consensus"]))
        payload.append(entry["mode"])
        mask_payload = entry["mask_payload"]
        payload.extend(_write_varint(entry["deviation_count"]))
        payload.extend(_write_varint(len(mask_payload)))
        payload.extend(mask_payload)
        residues = entry["residues"]
        payload.extend(struct.pack(">H", len(residues)))
        payload.extend(residues)

    payload.extend(struct.pack(">I", len(prepared_blocks)))
    for block, encoded in zip_strict(prepared_blocks, encoded_bitmasks):
        if not (1 <= block.run_length <= 255):
            raise EncodingError("Run length must be within 1..255 for encoding")
        key = (block.consensus, block.bitmask, block.residues)
        dict_id = dictionary_map.get(key)
        if dict_id is not None:
            payload.append(1)
            payload.append(dict_id)
            payload.append(block.run_length)
        else:
            payload.append(0)
            payload.append(block.run_length)
            payload.append(ord(block.consensus))
            mode, deviation_count, mask_payload = encoded
            payload.append(mode)
            payload.extend(_write_varint(deviation_count))
            payload.extend(_write_varint(len(mask_payload)))
            payload.extend(mask_payload)
            payload.extend(struct.pack(">H", len(block.residues)))
            payload.extend(block.residues)
    return bytes(payload)


def decode_blocks(
    payload: bytes,
    bitmask_bytes: int,
    bits_per_symbol: int,
    alphabet: Sequence[str],
) -> list[RunLengthBlock]:
    """Parse binary payload into run-length blocks."""

    blocks: list[RunLengthBlock] = []
    cursor = 0
    payload_length = len(payload)
    if cursor >= payload_length:
        return blocks

    table_count = payload[cursor]
    cursor += 1
    consensus_models: dict[str, dict[str, Any]] = {}
    for _ in range(table_count):
        if cursor + 3 > payload_length:
            raise DecodingError("Consensus table truncated")
        consensus_value = payload[cursor]
        cursor += 1
        mode = payload[cursor]
        cursor += 1
        entry_count = payload[cursor]
        cursor += 1
        residues: list[str] = []
        for _ in range(entry_count):
            if cursor >= payload_length:
                raise DecodingError("Consensus table residues truncated")
            residues.append(bytes([payload[cursor]]).decode("ascii"))
            cursor += 1
        if mode == 0:
            if cursor >= payload_length:
                raise DecodingError("Consensus width missing")
            width = payload[cursor]
            cursor += 1
            model = {
                "mode": 0,
                "residues": residues,
                "width": max(1, width),
                "decode_map": {idx: res for idx, res in enumerate(residues)},
            }
        elif mode == 1:
            lengths = []
            for _ in range(entry_count):
                if cursor >= payload_length:
                    raise DecodingError("Consensus code lengths truncated")
                lengths.append(payload[cursor])
                cursor += 1
            encode_map, decode_map, max_len = _canonical_code_maps(residues, lengths)
            model = {
                "mode": 1,
                "residues": residues,
                "lengths": lengths,
                "decode_map": decode_map,
                "max_code_len": max_len,
            }
        else:
            raise DecodingError(f"Unknown residue encoding mode {mode}")
        consensus_models[bytes([consensus_value]).decode("ascii")] = model

    if cursor >= payload_length:
        raise DecodingError("Missing dictionary section")

    dict_count = payload[cursor]
    cursor += 1
    dictionary: list[tuple[str, bytes, bytes]] = []
    for _ in range(dict_count):
        if cursor + 2 > payload_length:
            raise DecodingError("Dictionary entry truncated")
        consensus_value = payload[cursor]
        cursor += 1
        mode = payload[cursor]
        cursor += 1
        deviation_count, cursor = _read_varint(payload, cursor)
        mask_len, cursor = _read_varint(payload, cursor)
        if cursor + mask_len > payload_length:
            raise DecodingError("Dictionary mask payload truncated")
        mask_payload = payload[cursor : cursor + mask_len]
        cursor += mask_len
        if cursor + 2 > payload_length:
            raise DecodingError("Dictionary residue length truncated")
        residues_len = struct.unpack(">H", payload[cursor : cursor + 2])[0]
        cursor += 2
        if cursor + residues_len > payload_length:
            raise DecodingError("Dictionary residues truncated")
        residues = payload[cursor : cursor + residues_len]
        cursor += residues_len
        bitmask = _decode_bitmask(mode, mask_payload, deviation_count, bitmask_bytes)
        dictionary.append((bytes([consensus_value]).decode("ascii"), bitmask, residues))

    if cursor + 4 > payload_length:
        raise DecodingError("Missing block count")
    (block_count,) = struct.unpack(">I", payload[cursor : cursor + 4])
    cursor += 4

    alphabet_lookup = {char: index for index, char in enumerate(alphabet)}

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
            residues = _decode_residue_stream(
                consensus,
                bitmask,
                residues,
                consensus_models,
                bits_per_symbol,
                alphabet_lookup,
            )
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
            mode = payload[cursor]
            cursor += 1
            deviation_count, cursor = _read_varint(payload, cursor)
            mask_len, cursor = _read_varint(payload, cursor)
            if cursor + mask_len > payload_length:
                raise DecodingError("Literal mask truncated")
            mask_payload = payload[cursor : cursor + mask_len]
            cursor += mask_len
            if cursor + 2 > payload_length:
                raise DecodingError("Literal residue length truncated")
            residues_len = struct.unpack(">H", payload[cursor : cursor + 2])[0]
            cursor += 2
            if cursor + residues_len > payload_length:
                raise DecodingError("Literal residues truncated")
            residues = payload[cursor : cursor + residues_len]
            cursor += residues_len
            bitmask = _decode_bitmask(mode, mask_payload, deviation_count, bitmask_bytes)
            consensus = bytes([consensus_value]).decode("ascii")
            residues = _decode_residue_stream(
                consensus,
                bitmask,
                residues,
                consensus_models,
                bits_per_symbol,
                alphabet_lookup,
            )
            blocks.append(
                RunLengthBlock(
                    consensus=consensus,
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


def _varint_size(value: int) -> int:
    size = 1
    while value >= 0x80:
        value >>= 7
        size += 1
    return size


def _write_varint(value: int) -> bytes:
    if value < 0:
        raise EncodingError("varint cannot encode negative values")
    out = bytearray()
    while True:
        byte = value & 0x7F
        value >>= 7
        if value:
            out.append(byte | 0x80)
        else:
            out.append(byte)
            break
    return bytes(out)


def _read_varint(data: bytes, cursor: int) -> tuple[int, int]:
    shift = 0
    result = 0
    while True:
        if cursor >= len(data):
            raise DecodingError("Unexpected end of data while decoding varint")
        byte = data[cursor]
        cursor += 1
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return result, cursor
        shift += 7
        if shift > 56:
            raise DecodingError("Varint too long")


def _pack_codes(codes: Sequence[int], bits_per_symbol: int) -> bytes:
    if not codes:
        return b""
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
        output.append((buffer << (8 - bits_in_buffer)) & 0xFF)
    return bytes(output)


def _unpack_codes(data: bytes, count: int, bits_per_symbol: int) -> list[int]:
    if count == 0:
        return []
    values: list[int] = []
    buffer = 0
    bits_in_buffer = 0
    mask = (1 << bits_per_symbol) - 1
    for byte in data:
        buffer = (buffer << 8) | byte
        bits_in_buffer += 8
        while bits_in_buffer >= bits_per_symbol and len(values) < count:
            bits_in_buffer -= bits_per_symbol
            values.append((buffer >> bits_in_buffer) & mask)
            buffer &= (1 << bits_in_buffer) - 1
        if len(values) == count:
            break
    if len(values) != count:
        raise EncodingError("Residue stream truncated while unpacking")
    return values


def _decode_residue_stream(
    consensus: str,
    bitmask: bytes,
    encoded: bytes,
    models: dict[str, dict[str, Any]],
    bits_per_symbol: int,
    alphabet_lookup: dict[str, int],
) -> bytes:
    deviation_count = _popcount(bitmask)
    if deviation_count == 0:
        return b""
    model = models.get(consensus)
    if model is None:
        raise DecodingError("Missing residue model for consensus")

    if model["mode"] == 0:
        width = model["width"]
        residues = model["residues"]
        local_codes = _unpack_codes(encoded, deviation_count, width)
        try:
            chars = [residues[code] for code in local_codes]
        except IndexError as exc:
            raise DecodingError("Residue code exceeds table size") from exc
    else:
        decode_map = model["decode_map"]
        max_len = model["max_code_len"]
        chars: list[str] = []
        total_bits = len(encoded) * 8
        bit_pos = 0
        for _ in range(deviation_count):
            current = 0
            length = 0
            while True:
                if bit_pos >= total_bits:
                    raise DecodingError("Insufficient Huffman bits for residue stream")
                byte = encoded[bit_pos // 8]
                bit = (byte >> (7 - (bit_pos % 8))) & 1
                bit_pos += 1
                current = (current << 1) | bit
                length += 1
                residue = decode_map.get((length, current))
                if residue is not None:
                    chars.append(residue)
                    break
                if length > max_len:
                    raise DecodingError("Invalid Huffman code in residue stream")

    try:
        global_codes = [alphabet_lookup[char] for char in chars]
    except KeyError as exc:
        raise DecodingError("Residue not present in alphabet lookup") from exc
    return _pack_codes(global_codes, bits_per_symbol)


def _bit_positions(bitmask: bytes) -> list[int]:
    positions: list[int] = []
    bit_index = 0
    for byte in bitmask:
        for shift in range(8):
            if byte & (1 << shift):
                positions.append(bit_index)
            bit_index += 1
    return positions


def _run_length_encode(bitmask: bytes) -> bytes:
    out = bytearray()
    prev = bitmask[0]
    count = 1
    for byte in bitmask[1:]:
        if byte == prev and count < 255:
            count += 1
        else:
            out.append(prev)
            out.append(count)
            prev = byte
            count = 1
    out.append(prev)
    out.append(count)
    return bytes(out)


def _encode_bitmask(bitmask: bytes, bitmask_bytes: int) -> tuple[int, int, bytes]:
    deviation_count = _popcount(bitmask)
    if deviation_count == 0:
        return 0, 0, b""

    trimmed, _ = _trim_bitmask(bitmask)
    raw_payload = trimmed

    positions = _bit_positions(bitmask)
    sparse_buffer = bytearray()
    prev = -1
    for pos in positions:
        delta = pos - prev
        sparse_buffer.extend(_write_varint(delta))
        prev = pos
    sparse_payload = bytes(sparse_buffer)

    rle_payload = _run_length_encode(trimmed)

    candidates = [
        (0, raw_payload),
        (1, sparse_payload),
        (2, rle_payload),
    ]
    mode, payload_bytes = min(candidates, key=lambda item: len(item[1]))
    if mode == 1:
        wrapper = _write_varint(len(payload_bytes)) + payload_bytes
        return 1, deviation_count, wrapper
    return mode, deviation_count, payload_bytes


def _decode_positions(encoded: bytes, count: int) -> list[int]:
    positions: list[int] = []
    cursor = 0
    prev = -1
    for _ in range(count):
        value, cursor = _read_varint(encoded, cursor)
        pos = prev + value
        positions.append(pos)
        prev = pos
    return positions


def _decode_bitmask(
    mode: int,
    payload: bytes,
    deviation_count: int,
    bitmask_bytes: int,
) -> bytes:
    if deviation_count == 0:
        return bytes(bitmask_bytes)
    if mode == 0:
        bitmask = bytearray(bitmask_bytes)
        length = min(len(payload), bitmask_bytes)
        bitmask[:length] = payload[:length]
        return bytes(bitmask)
    if mode == 1:
        length, cursor = _read_varint(payload, 0)
        encoded = payload[cursor: cursor + length]
        positions = _decode_positions(encoded, deviation_count)
        bitmask = bytearray(bitmask_bytes)
        max_bits = bitmask_bytes * 8
        for pos in positions:
            if pos >= max_bits:
                raise DecodingError("Bitmask position exceeds sequence count")
            byte_index = pos // 8
            bit_index = pos % 8
            bitmask[byte_index] |= 1 << bit_index
        return bytes(bitmask)
    if mode == 2:
        bitmask = bytearray()
        cursor = 0
        while cursor < len(payload) and len(bitmask) < bitmask_bytes:
            if cursor + 2 > len(payload):
                raise DecodingError("RLE bitmask truncated")
            byte = payload[cursor]
            count = payload[cursor + 1]
            cursor += 2
            bitmask.extend([byte] * count)
        if len(bitmask) < bitmask_bytes:
            bitmask.extend([0] * (bitmask_bytes - len(bitmask)))
        return bytes(bitmask[:bitmask_bytes])
    raise DecodingError(f"Unknown bitmask encoding mode {mode}")


def _literal_block_size(
    mode: int, deviation_count: int, mask_payload_len: int, residue_len: int
) -> int:
    return (
        1  # marker
        + 1  # run length
        + 1  # consensus
        + 1  # mode
        + _varint_size(deviation_count)
        + _varint_size(mask_payload_len)
        + mask_payload_len
        + 2  # residue length header
        + residue_len
    )


def _dictionary_entry_size(
    mode: int, deviation_count: int, mask_payload_len: int, residue_len: int
) -> int:
    return (
        1  # consensus
        + 1  # mode
        + _varint_size(deviation_count)
        + _varint_size(mask_payload_len)
        + mask_payload_len
        + 2
        + residue_len
    )


def _build_dictionary(
    blocks: Sequence[RunLengthBlock],
    encoded_bitmasks: Sequence[tuple[int, int, bytes]],
) -> tuple[list[dict[str, Any]], dict[tuple[str, bytes, bytes], int]]:
    counter: Counter[tuple[str, bytes, bytes]] = Counter()
    info_lookup: dict[tuple[str, bytes, bytes], tuple[int, int, bytes]] = {}
    residue_len_lookup: dict[tuple[str, bytes, bytes], int] = {}
    bitmask_lookup: dict[tuple[str, bytes, bytes], bytes] = {}

    for block, encoded in zip_strict(blocks, encoded_bitmasks):
        key = (block.consensus, block.bitmask, block.residues)
        counter[key] += block.run_length
        if key not in info_lookup:
            info_lookup[key] = encoded
            residue_len_lookup[key] = len(block.residues)
            bitmask_lookup[key] = block.bitmask

    dictionary_entries: list[dict[str, Any]] = []
    dictionary_map: dict[tuple[str, bytes, bytes], int] = {}
    reference_size = 3  # marker + dict id + run length

    for key, freq in counter.most_common(255):
        mode, deviation_count, mask_payload = info_lookup[key]
        mask_len = len(mask_payload)
        residue_len = residue_len_lookup[key]
        literal_size = _literal_block_size(mode, deviation_count, mask_len, residue_len)
        entry_size = _dictionary_entry_size(mode, deviation_count, mask_len, residue_len)
        if freq * literal_size <= entry_size + freq * reference_size:
            continue
        dictionary_map[key] = len(dictionary_entries)
        dictionary_entries.append(
            {
                "consensus": key[0],
                "mode": mode,
                "deviation_count": deviation_count,
                "mask_payload": mask_payload,
                "bitmask": bitmask_lookup[key],
                "residues": key[2],
            }
        )

    return dictionary_entries, dictionary_map


def _build_huffman_lengths(
    residues: Sequence[str], freqs: Counter[str]
) -> tuple[list[int] | None, int | None]:
    if len(residues) <= 1:
        return None, None

    heap: list[tuple[int, int, Any]] = []
    counter = 0
    for residue in residues:
        freq = freqs.get(residue, 0)
        if freq <= 0:
            return None, None
        heapq.heappush(heap, (freq, counter, {"symbol": residue}))
        counter += 1

    while len(heap) > 1:
        freq1, _, left = heapq.heappop(heap)
        freq2, _, right = heapq.heappop(heap)
        node = {"left": left, "right": right}
        heapq.heappush(heap, (freq1 + freq2, counter, node))
        counter += 1

    root = heap[0][2]
    lengths: dict[str, int] = {}

    def assign(node: dict[str, Any], depth: int) -> None:
        if "symbol" in node:
            lengths[node["symbol"]] = max(depth, 1)
            return
        assign(node["left"], depth + 1)
        assign(node["right"], depth + 1)

    assign(root, 0)

    total_bits = sum(freqs[res] * lengths[res] for res in residues)
    length_list = [lengths[res] for res in residues]
    return length_list, total_bits


def _canonical_code_maps(
    residues: Sequence[str], lengths: Sequence[int]
) -> tuple[dict[str, tuple[int, int]], dict[tuple[int, int], str], int]:
    items = sorted(zip(residues, lengths), key=lambda item: (item[1], item[0]))
    encode_map: dict[str, tuple[int, int]] = {}
    decode_map: dict[tuple[int, int], str] = {}
    code = 0
    prev_len = 0
    max_len = 0
    for residue, length in items:
        if length <= 0:
            continue
        code <<= length - prev_len
        encode_map[residue] = (code, length)
        decode_map[(length, code)] = residue
        prev_len = length
        code += 1
        if length > max_len:
            max_len = length
    return encode_map, decode_map, max_len


def _pack_huffman_codes(
    encode_map: dict[str, tuple[int, int]], chars: Sequence[str]
) -> bytes:
    buffer = 0
    bits_in_buffer = 0
    output = bytearray()
    for char in chars:
        code, length = encode_map[char]
        buffer = (buffer << length) | code
        bits_in_buffer += length
        while bits_in_buffer >= 8:
            bits_in_buffer -= 8
            output.append((buffer >> bits_in_buffer) & 0xFF)
            buffer &= (1 << bits_in_buffer) - 1
    if bits_in_buffer:
        output.append((buffer << (8 - bits_in_buffer)) & 0xFF)
    return bytes(output)
