from pathlib import Path

import base64
import zlib
import pytest

from ecomp.io import alignment_from_sequences

from ecomp.compression.pipeline import (
    _SEQ_ID_MAGIC,
    _ZSTD_DECOMPRESSOR,
    _WIDTH_TO_CODE,
    _alignment_to_fasta_bytes,
    _build_distance_matrix,
    _build_permutation_chunk,
    _compute_alignment_stats,
    _decode_permutation,
    _decode_sequence_ids,
    _decode_residues,
    _encode_sequence_ids,
    _encode_varint,
    _extract_permutation_chunk,
    _iter_deviation_indices,
    _maybe_use_gzip_fallback,
    _parse_fasta_bytes,
    _select_sample_indices,
)
from ecomp.compression.encoding import (
    DecodingError,
    EncodingError,
    RunLengthBlock,
    _decode_bitmask,
    _decode_residue_stream,
    _encode_bitmask,
    _pack_codes,
)


def test_iter_deviation_indices_returns_expected_positions():
    bitmask = bytes([0b00010101])
    indices = _iter_deviation_indices(bitmask, num_sequences=8)
    assert indices == [0, 2, 4]


def test_decode_residues_round_trip():
    data = bytes([0b00011011])  # packed values 0,1,2,3 for bits_per_symbol=2
    alphabet = ["A", "C", "G", "T"]
    residues = _decode_residues(data, count=4, bits_per_symbol=2, alphabet=alphabet)
    assert residues == ["A", "C", "G", "T"]


def test_decode_residues_raises_on_insufficient_data():
    with pytest.raises(ValueError):
        _decode_residues(b"\x00", count=5, bits_per_symbol=3, alphabet=["A", "C", "G", "T"])


def test_decode_residue_stream_mode_zero_round_trip():
    bitmask = bytes([0b00000001])
    models = {
        "A": {
            "mode": 0,
            "width": 1,
            "residues": ["G", "T"],
        }
    }
    encoded = _pack_codes([1], bits_per_symbol=1)
    alphabet_lookup = {char: idx for idx, char in enumerate(["A", "C", "G", "T"])}
    packed = _decode_residue_stream(
        "A",
        bitmask,
        encoded,
        models,
        bits_per_symbol=2,
        alphabet_lookup=alphabet_lookup,
    )
    decoded = _decode_residues(packed, count=1, bits_per_symbol=2, alphabet=["A", "C", "G", "T"])
    assert decoded == ["T"]


def test_decode_residue_stream_huffman_round_trip():
    bitmask = bytes([0b00000011])
    models = {
        "A": {
            "mode": 1,
            "decode_map": {
                (1, 0): "G",
                (2, 2): "T",
            },
            "max_code_len": 2,
        }
    }
    # Bits: 0 (for G), 1 0 (for T) -> 0b01000000
    encoded = bytes([0b01000000])
    alphabet_lookup = {char: idx for idx, char in enumerate(["A", "C", "G", "T"])}
    packed = _decode_residue_stream(
        "A",
        bitmask,
        encoded,
        models,
        bits_per_symbol=2,
        alphabet_lookup=alphabet_lookup,
    )
    decoded = _decode_residues(packed, count=2, bits_per_symbol=2, alphabet=["A", "C", "G", "T"])
    assert decoded == ["G", "T"]


def test_decode_residue_stream_errors_when_model_missing():
    bitmask = bytes([0b00000001])
    with pytest.raises(DecodingError):
        _decode_residue_stream(
            "A",
            bitmask,
            b"\x00",
            {},
            bits_per_symbol=2,
            alphabet_lookup={"A": 0},
        )


def test_decode_bitmask_sparse_and_rle_modes():
    bitmask = bytes([0b00010001])
    mode, deviation_count, payload = _encode_bitmask(bitmask, bitmask_bytes=len(bitmask))
    assert mode in {0, 1, 2}
    decoded = _decode_bitmask(mode, payload, deviation_count, len(bitmask))
    assert decoded == bitmask


def test_decode_bitmask_raises_on_unknown_mode():
    with pytest.raises(DecodingError):
        _decode_bitmask(5, b"payload", deviation_count=1, bitmask_bytes=1)


def test_alignment_to_fasta_bytes_round_trip(tmp_path: Path):
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["AC", "AG"])
    fasta_bytes = _alignment_to_fasta_bytes(frame)
    restored = _parse_fasta_bytes(fasta_bytes)
    assert restored.ids == frame.ids
    assert restored.sequences == frame.sequences


def test_maybe_use_gzip_fallback_replaces_payload(tmp_path: Path):
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["ACGT" * 10, "ACGT" * 10])
    payload = b"X" * 100
    metadata = {"codec": "ecomp", "source_format": "fasta"}
    new_payload, new_metadata = _maybe_use_gzip_fallback(frame, payload, metadata)
    assert new_payload != payload
    assert new_metadata["fallback"]["type"] == "gzip"


def test_permutation_chunk_round_trip_various_widths():
    permutations = [
        [2, 0, 1],
        list(range(400, -1, -1)),
        [0, 1, 70000, 2],
    ]
    for permutation in permutations:
        chunk, meta = _build_permutation_chunk(permutation)
        assert meta is not None
        payload = chunk + b"payload"
        remainder, decoded = _extract_permutation_chunk(payload, meta)
        assert remainder == b"payload"
        assert decoded == permutation


def test_build_permutation_chunk_handles_empty_permutation():
    chunk, meta = _build_permutation_chunk([])
    assert chunk is None
    assert meta is None


def test_decode_permutation_dict_format_round_trip():
    permutation = [3, 1, 4, 2]
    raw = b"".join(int(value).to_bytes(1, "little") for value in permutation)
    meta = {
        "encoding": "map",
        "version": 1,
        "dtype": "uint8",
        "compression": "none",
        "size": len(permutation),
        "data": base64.b64encode(raw).decode("ascii"),
    }
    assert _decode_permutation(permutation) == permutation
    assert _decode_permutation(meta) == permutation


def _decode_mode(encoded: bytes) -> int:
    idx = len(_SEQ_ID_MAGIC)
    idx += 1  # version byte
    while True:
        byte = encoded[idx]
        idx += 1
        if not byte & 0x80:
            break
    return encoded[idx]


def test_encode_sequence_ids_round_trip_zlib_mode():
    ids = [f"id{i}" for i in range(8)]
    encoded = _encode_sequence_ids(ids)
    mode = _decode_mode(encoded)
    assert mode in {0, 2}
    decoded, remainder = _decode_sequence_ids(encoded)
    assert remainder == b""
    assert decoded == ids


@pytest.mark.skipif(_ZSTD_DECOMPRESSOR is None, reason="zstd optional dependency not available")
def test_encode_sequence_ids_round_trip_zstd_mode():
    ids = [f"tax{str(i).zfill(4)}" for i in range(512)]
    encoded = _encode_sequence_ids(ids)
    mode = _decode_mode(encoded)
    assert mode in {0, 1}  # expect zstd when beneficial
    decoded, remainder = _decode_sequence_ids(encoded)
    assert decoded == ids
    assert remainder == b""


def test_extract_permutation_chunk_validates_metadata():
    with pytest.raises(ValueError):
        _extract_permutation_chunk(b"payload", {"length": 0})
    with pytest.raises(ValueError):
        _extract_permutation_chunk(b"short", {"length": 10})


def test_extract_permutation_chunk_checks_magic():
    payload = b"\x00" * 10
    with pytest.raises(ValueError):
        _extract_permutation_chunk(payload, {"length": 6})


def test_extract_permutation_chunk_rejects_unknown_width():
    chunk, meta = _build_permutation_chunk([0, 1, 2])
    mutated = bytearray(chunk)
    original_flags = mutated[5]
    mutated[5] = (original_flags & 0x01) | (0b11 << 1)
    with pytest.raises(ValueError):
        _extract_permutation_chunk(bytes(mutated) + b"suffix", {"length": len(mutated)})


def test_extract_permutation_chunk_rejects_bad_version():
    chunk, meta = _build_permutation_chunk([0, 1, 2])
    mutated = bytearray(chunk)
    mutated[4] = 0xFF
    with pytest.raises(ValueError):
        _extract_permutation_chunk(bytes(mutated) + b"rest", {"length": len(mutated)})


def test_extract_permutation_chunk_handles_compressed_payload():
    permutation = list(range(64))
    raw = bytearray()
    width = 2  # force uint16 width for better compression
    for value in permutation:
        raw.extend(value.to_bytes(width, "little"))
    compressed = zlib.compress(bytes(raw), level=9)
    header = bytearray()
    header.extend(b"ECPE")
    header.append(1)  # version
    flags = (_WIDTH_TO_CODE[width] << 1) | 1
    header.append(flags)
    header.extend(_encode_varint(len(permutation)))
    header.extend(_encode_varint(len(compressed)))
    chunk = bytes(header + compressed)
    remainder, decoded = _extract_permutation_chunk(chunk + b"tail", {"length": len(chunk)})
    assert remainder == b"tail"
    assert decoded == permutation


def test_extract_permutation_chunk_detects_length_mismatch():
    permutation = [0, 1, 2]
    chunk, meta = _build_permutation_chunk(permutation)
    truncated = chunk[:-1]
    with pytest.raises(ValueError):
        _extract_permutation_chunk(truncated, {"length": len(chunk)})


def test_decode_permutation_handles_error_conditions():
    with pytest.raises(ValueError):
        _decode_permutation({"encoding": "payload"})
    with pytest.raises(ValueError):
        _decode_permutation({"version": 99})
    payload = base64.b64encode(b"\x00\x01").decode("ascii")
    with pytest.raises(ValueError):
        _decode_permutation({"version": 1, "dtype": "uint64", "compression": "none", "size": 2, "data": payload})
    with pytest.raises(ValueError):
        _decode_permutation({"version": 1, "dtype": "uint8", "compression": "brotli", "size": 2, "data": payload})


def test_decode_sequence_ids_error_conditions():
    with pytest.raises(ValueError):
        _decode_sequence_ids(b"short")
    with pytest.raises(ValueError):
        _decode_sequence_ids(b"ECID\x01\x01")  # truncated block
    block = bytearray(b"ECID\x02")
    block.extend(_encode_varint(1))
    block.extend(b"\x03")  # unsupported mode byte without payload
    with pytest.raises(ValueError):
        _decode_sequence_ids(bytes(block))


def test_select_sample_indices_and_distance_matrix():
    indices = _select_sample_indices(1000, cap=10)
    assert indices[0] == 0 and indices[-1] == 999
    matrix = _build_distance_matrix(["AAAA", "AAAT", "AATT"], sample_indices=[2, 3])
    assert matrix[1][2] == 1


def test_compute_alignment_stats_handles_empty_input():
    empty = alignment_from_sequences(ids=["s1"], sequences=[""])
    assert _compute_alignment_stats(empty) is None
