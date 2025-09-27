from pathlib import Path

import pytest

from evolutionary_compression.io import AlignmentFrame, alignment_from_sequences, read_alignment, write_alignment


def test_alignment_frame_alphabet_string_returns_sorted_unique():
    frame = alignment_from_sequences(ids=["a", "b"], sequences=["AC", "AG"], alphabet=["G", "A", "C"])
    assert frame.alphabet_string() == "ACG"


def test_write_and_read_alignment_round_trip(tmp_path: Path):
    frame = alignment_from_sequences(
        ids=["seq1", "seq2"],
        sequences=["ACGT", "ACGA"],
        alphabet=["A", "C", "G", "T"],
    )
    output = tmp_path / "roundtrip.fasta"
    write_alignment(frame, output)
    restored = read_alignment(output, fmt="fasta")
    assert restored.ids == frame.ids
    assert restored.sequences == frame.sequences


def test_alignment_frame_validation_errors():
    frame = AlignmentFrame(ids=["s1"], sequences=[""], alphabet=["A"])
    assert frame.num_sequences == 1
    assert frame.alignment_length == 0

    with pytest.raises(ValueError):
        AlignmentFrame(ids=[], sequences=[], alphabet=[])

    with pytest.raises(ValueError):
        AlignmentFrame(ids=["s1"], sequences=[], alphabet=["A"])

    with pytest.raises(ValueError):
        AlignmentFrame(ids=["s1", "s2"], sequences=["A", "AA"], alphabet=["A"])
