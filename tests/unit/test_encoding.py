import math

from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.compression.consensus import collect_column_profiles
import struct

import pytest

from evolutionary_compression.compression.encoding import (
    DecodingError,
    EncodingError,
    _build_dictionary,
    _decode_bitmask,
    _decode_residue_stream,
    _canonical_code_maps,
    _pack_huffman_codes,
    _pack_codes,
    _unpack_codes,
    _encode_char,
    _encode_bitmask,
    _popcount,
    _trim_bitmask,
    _write_varint,
    decode_blocks,
    encode_blocks,
)
from evolutionary_compression.compression.rle import RunLengthBlock, collect_run_length_blocks
from evolutionary_compression.compression import pipeline


def test_encode_decode_round_trip():
    frame = alignment_from_sequences(
        ids=["seq1", "seq2", "seq3"],
        sequences=[
            "ACGT",
            "ACGA",
            "ACGG",
        ],
    )
    columns = collect_column_profiles(frame)
    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, math.ceil(math.log2(max(len(alphabet), 1))))
    blocks = collect_run_length_blocks(
        columns, frame.num_sequences, symbol_lookup, bits_per_symbol
    )
    bitmask_bytes = (frame.num_sequences + 7) // 8
    payload = encode_blocks(
        blocks,
        bitmask_bytes=bitmask_bytes,
        bits_per_symbol=bits_per_symbol,
        alphabet=alphabet,
    )
    decoded = decode_blocks(
        payload,
        bitmask_bytes=bitmask_bytes,
        bits_per_symbol=bits_per_symbol,
        alphabet=alphabet,
    )
    assert decoded == blocks


def test_encode_blocks_raises_for_invalid_run_length():
    block = RunLengthBlock(consensus="A", bitmask=b"\x00", residues=b"", run_length=0)
    with pytest.raises(EncodingError):
        encode_blocks([block], bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_encode_blocks_emits_dictionary_entries():
    block = RunLengthBlock(consensus="A", bitmask=b"\x00", residues=b"", run_length=3)
    blocks = [block] * 6
    payload = encode_blocks(blocks, bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])
    decoded = decode_blocks(
        payload, bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"]
    )
    assert decoded == blocks


def test_decode_blocks_detects_unknown_marker():
    payload = bytearray()
    payload.append(0)  # consensus tables
    payload.append(0)  # dictionary size
    payload.extend(struct.pack(">I", 1))  # block count
    payload.append(99)  # invalid marker
    with pytest.raises(DecodingError):
        decode_blocks(bytes(payload), bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_detects_truncated_literal():
    # dictionary count zero, block count one, literal marker but truncated mask
    payload = bytearray()
    payload.append(0)  # consensus tables
    payload.append(0)  # dictionary size
    payload.extend((0, 0, 0, 1))  # block count = 1
    payload.append(0)  # literal marker
    payload.append(1)  # run length
    payload.append(ord("A"))
    payload.append(0)  # mode
    payload.append(0)  # deviation count varint
    payload.append(1)  # mask length varint (encoded as single byte 1)
    with pytest.raises(DecodingError):
        decode_blocks(bytes(payload), bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_errors_on_truncated_consensus_table():
    payload = bytes([1, ord("A")])
    with pytest.raises(DecodingError):
        decode_blocks(payload, bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_errors_on_missing_block_count():
    payload = bytes([0, 0])
    with pytest.raises(DecodingError):
        decode_blocks(payload, bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_errors_on_dictionary_truncation():
    payload = bytearray()
    payload.append(0)
    payload.append(1)
    payload.append(ord("A"))
    payload.append(0)
    payload.append(1)
    payload.append(1)
    payload.extend(b"\x01")
    with pytest.raises(DecodingError):
        decode_blocks(bytes(payload), bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_errors_on_truncated_consensus_table():
    payload = bytes([1, ord("A")])  # table count followed by incomplete entry
    with pytest.raises(DecodingError):
        decode_blocks(payload, bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_errors_on_missing_block_count():
    payload = bytes([0, 0])  # no consensus tables, no dictionary, but missing block count
    with pytest.raises(DecodingError):
        decode_blocks(payload, bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_decode_blocks_errors_on_dictionary_truncation():
    payload = bytearray()
    payload.append(0)  # consensus tables
    payload.append(1)  # dictionary count
    payload.append(ord("A"))  # consensus
    payload.append(0)  # mode
    payload.append(1)  # deviation count (varint)
    payload.append(1)  # mask length (varint)
    payload.extend(b"\x01")  # mask payload but missing residues length
    with pytest.raises(DecodingError):
        decode_blocks(bytes(payload), bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"])


def test_trim_bitmask_and_popcount_helpers():
    mask = b"\xff\x00\x10\x00\x00"
    trimmed, length = _trim_bitmask(mask)
    assert trimmed == b"\xff\x00\x10"
    assert length == 3
    assert _popcount(trimmed) == 9


def test_encode_char_validates_single_ascii_character():
    assert _encode_char("A") == b"A"
    with pytest.raises(EncodingError):
        _encode_char("BC")
    with pytest.raises(EncodingError):
        _encode_char("Ω")


def test_build_dictionary_skips_unprofitable_entries():
    block = RunLengthBlock(consensus="A", bitmask=b"\x00", residues=b"", run_length=1)
    dictionary, mapping = _build_dictionary([block], [(0, 0, b"")])
    assert dictionary == []
    assert mapping == {}


def test_sequence_id_block_roundtrip():
    ids = [f"taxon_{i:03d}" for i in range(50)]
    block = pipeline._encode_sequence_ids(ids)
    decoded, remaining = pipeline._decode_sequence_ids(block)
    assert decoded == ids
    assert remaining == b""


def test_huffman_residue_encoding_header():
    frame = alignment_from_sequences(
        ids=[f"s{i}" for i in range(6)],
        sequences=[
            "A" * 12,
            "ACACACACACAC",
            "AGAAAGAAAGAA",
            "ATAAAATAAAAT",
            "ACAAAAACAAAA",
            "AGAAAAAGAAAA",
        ],
    )
    columns = collect_column_profiles(frame)
    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, math.ceil(math.log2(max(len(alphabet), 1))))
    blocks = collect_run_length_blocks(
        columns, frame.num_sequences, symbol_lookup, bits_per_symbol
    )
    bitmask_bytes = (frame.num_sequences + 7) // 8
    payload = encode_blocks(blocks, bitmask_bytes=bitmask_bytes, bits_per_symbol=bits_per_symbol, alphabet=alphabet)
    assert payload[0] >= 1
    mode = payload[2]
    assert mode in {0, 1}
    decoded = decode_blocks(payload, bitmask_bytes=bitmask_bytes, bits_per_symbol=bits_per_symbol, alphabet=alphabet)
    assert decoded == blocks


def test_decode_residue_stream_missing_model_raises():
    bitmask = bytes([0b00000001])
    alphabet_lookup = {"A": 0}
    with pytest.raises(DecodingError):
        _decode_residue_stream("A", bitmask, b"\x00", {}, bits_per_symbol=1, alphabet_lookup=alphabet_lookup)


def test_decode_residue_stream_invalid_huffman_code():
    residues = ["A", "C"]
    lengths = [1, 3]
    encode_map, decode_map, max_len = _canonical_code_maps(residues, lengths)
    encoded = bytes([0b11111111])
    model = {"mode": 1, "decode_map": decode_map, "max_code_len": max_len}
    alphabet_lookup = {"A": 0, "C": 1}
    bitmask = bytes([0b00000001])
    with pytest.raises(DecodingError):
        _decode_residue_stream("A", bitmask, encoded, {"A": model}, bits_per_symbol=1, alphabet_lookup=alphabet_lookup)


def test_decode_residue_stream_huffman_path():
    residues = ["A", "C"]
    lengths = [1, 3]
    encode_map, decode_map, max_len = _canonical_code_maps(residues, lengths)
    encoded = _pack_huffman_codes(encode_map, ["A", "C", "C"])
    model = {"mode": 1, "decode_map": decode_map, "max_code_len": max_len}
    alphabet_lookup = {"A": 0, "C": 1}
    bitmask = bytes([0b00000111])
    result = _decode_residue_stream("A", bitmask, encoded, {"A": model}, bits_per_symbol=1, alphabet_lookup=alphabet_lookup)
    unpacked = _unpack_codes(result, 3, 1)
    assert unpacked == [0, 1, 1]


def test_decode_blocks_handles_empty_payload():
    assert decode_blocks(b"", bitmask_bytes=1, bits_per_symbol=1, alphabet=["A"]) == []


def test_decode_bitmask_round_trip_modes():
    # mode 0 (raw)
    bitmask_raw = bytes([0b10101010, 0b01010101])
    mode, count, payload = _encode_bitmask(bitmask_raw, bitmask_bytes=len(bitmask_raw))
    assert mode == 0
    decoded = _decode_bitmask(mode, payload, count, len(bitmask_raw))
    assert decoded[: len(bitmask_raw)] == bitmask_raw

    # mode 1 (sparse varint encoding)
    deltas = _write_varint(1) + _write_varint(8)
    wrapper = _write_varint(len(deltas)) + deltas
    decoded_sparse = _decode_bitmask(1, wrapper, deviation_count=2, bitmask_bytes=2)
    assert decoded_sparse[:2] == bytes([0b00000001, 0b00000001])

    # mode 2 (run-length encoded)
    rle_payload = bytes([0b00000001, 4])
    decoded_rle = _decode_bitmask(2, rle_payload, deviation_count=4, bitmask_bytes=4)
    assert decoded_rle == bytes([0b00000001] * 4)

    # zero deviations short-circuits to zeros
    zero = _decode_bitmask(0, b"", deviation_count=0, bitmask_bytes=3)
    assert zero == bytes(3)
