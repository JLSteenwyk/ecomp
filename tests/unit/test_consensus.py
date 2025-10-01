from ecomp.io import alignment_from_sequences
from ecomp.compression.consensus import collect_column_profiles


def test_collect_column_profiles_identifies_consensus_and_deviations():
    frame = alignment_from_sequences(
        ids=["seq1", "seq2", "seq3"],
        sequences=["AAAA", "AAAT", "AATA"],
    )
    columns = collect_column_profiles(frame)
    assert len(columns) == 4
    assert columns[0].consensus == "A"
    assert columns[0].deviations == ()
    assert columns[2].consensus == "A"
    assert columns[2].deviations == ((2, "T"),)
    assert columns[3].consensus == "A"
    assert columns[3].deviations == ((1, "T"),)
    assert columns[2].equivalent_key() == ("A", ((2, "T"),))


def test_consensus_prefers_lexicographically_smaller_char_on_tie():
    frame = alignment_from_sequences(
        ids=["seq1", "seq2"],
        sequences=["AT", "TA"],
    )
    columns = collect_column_profiles(frame)
    assert columns[0].consensus == "A"
    assert columns[1].consensus == "A"


def test_collect_column_profiles_handles_zero_length_sequences():
    frame = alignment_from_sequences(ids=["seq1"], sequences=[""])
    assert collect_column_profiles(frame) == []
