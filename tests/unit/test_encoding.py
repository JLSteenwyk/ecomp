import math

from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.compression.consensus import collect_column_profiles
import struct

import pytest

from evolutionary_compression.compression.encoding import (
    DecodingError,
    EncodingError,
    _build_dictionary,
    _encode_char,
    _popcount,
    _trim_bitmask,
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
