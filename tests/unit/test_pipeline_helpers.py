from pathlib import Path

import pytest

from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.compression.pipeline import (
    _alignment_to_fasta_bytes,
    _decode_residues,
    _iter_deviation_indices,
    _maybe_use_gzip_fallback,
    _parse_fasta_bytes,
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
