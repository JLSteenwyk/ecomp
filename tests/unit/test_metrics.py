from __future__ import annotations

import math
from collections import Counter

from ecomp import (
    PairwiseIdentityResult,
    alignment_compressed_length,
    alignment_length_excluding_gaps,
    alignment_from_sequences,
    column_base_counts,
    column_gap_fraction,
    column_shannon_entropy,
    constant_columns,
    majority_rule_consensus,
    pairwise_identity_matrix,
    parsimony_informative_columns,
    parsimony_informative_site_count,
    percentage_identity,
    relative_composition_variability,
    variable_site_count,
)


def _example_alignment():
    ids = ["s1", "s2", "s3", "s4"]
    sequences = [
        "ACGTT",
        "A-GTT",
        "ATGCT",
        "ATGCT",
    ]
    return alignment_from_sequences(ids, sequences)


def test_column_base_counts_supports_gaps_toggle():
    frame = _example_alignment()
    counts = column_base_counts(frame, include_gaps=False)
    counts_with_gaps = column_base_counts(frame, include_gaps=True)

    assert counts[0] == Counter({"A": 4})
    assert counts[1] == Counter({"T": 2, "C": 1})
    assert counts_with_gaps[1] == Counter({"T": 2, "C": 1, "-": 1})


def test_majority_rule_consensus_uses_frequency_then_ordering():
    frame = _example_alignment()
    consensus = majority_rule_consensus(frame)
    assert consensus == "ATGYT"


def test_majority_rule_consensus_uses_iupac_for_ties():
    frame = alignment_from_sequences(
        ["s1", "s2", "s3", "s4"],
        [
            "AC",
            "AT",
            "GC",
            "GT",
        ],
    )
    consensus = majority_rule_consensus(frame)
    assert consensus == "RY"


def test_majority_rule_consensus_emits_n_for_full_mixture():
    frame = alignment_from_sequences(
        ["s1", "s2", "s3", "s4"],
        [
            "AC",
            "GT",
            "CA",
            "TG",
        ],
    )
    consensus = majority_rule_consensus(frame)
    assert consensus == "NN"


def test_majority_rule_consensus_returns_x_for_protein_ties():
    frame = alignment_from_sequences(
        ["s1", "s2", "s3", "s4"],
        [
            "AA",
            "AL",
            "LA",
            "LL",
        ],
    )
    consensus = majority_rule_consensus(frame)
    assert consensus == "XX"


def test_column_gap_fraction_values():
    frame = _example_alignment()
    fractions = column_gap_fraction(frame)
    assert fractions[0] == 0.0
    assert math.isclose(fractions[1], 0.25)
    assert fractions[2] == 0.0


def test_column_gap_fraction_treats_dot_as_residue():
    frame = alignment_from_sequences(
        ["a", "b"],
        ["A.", "A-"],
    )
    fractions = column_gap_fraction(frame)
    assert fractions == [0.0, 0.5]


def test_majority_rule_consensus_adds_gap_placeholder():
    frame = alignment_from_sequences(
        ["g1", "g2"],
        ["-", "-"],
    )
    consensus = majority_rule_consensus(frame)
    assert consensus == "-"


def test_alignment_length_and_variable_sites():
    frame = _example_alignment()
    assert alignment_length_excluding_gaps(frame) == 5
    assert variable_site_count(frame) == 2
    assert alignment_compressed_length(frame) == 2


def test_column_shannon_entropy_ignores_gaps():
    frame = _example_alignment()
    entropies = column_shannon_entropy(frame)
    assert entropies[0] == 0.0
    assert math.isclose(entropies[1], 0.918295, rel_tol=1e-6)
    assert math.isclose(entropies[3], 1.0, rel_tol=1e-6)


def test_parsimony_and_constant_masks():
    frame = _example_alignment()
    informative = parsimony_informative_columns(frame)
    constant = constant_columns(frame)

    assert informative == [False, False, False, True, False]
    assert constant == [True, False, True, False, True]
    assert parsimony_informative_site_count(frame) == 1


def test_pairwise_identity_matrix_matches_expected():
    frame = _example_alignment()
    result = pairwise_identity_matrix(frame)
    assert isinstance(result, PairwiseIdentityResult)
    matrix = result.matrix
    coverage = result.coverage

    expected = {
        (0, 1): (1.0, 4),
        (0, 2): (0.6, 5),
        (0, 3): (0.6, 5),
        (1, 2): (0.75, 4),
        (1, 3): (0.75, 4),
        (2, 3): (1.0, 5),
    }

    for i in range(len(frame.ids)):
        assert matrix[i][i] == 1.0
        for j in range(i + 1, len(frame.ids)):
            identity, cov = expected[(i, j)]
            assert math.isclose(matrix[i][j], identity, rel_tol=1e-9)
            assert math.isclose(matrix[j][i], identity, rel_tol=1e-9)
            assert coverage[i][j] == cov
            assert coverage[j][i] == cov


def test_pairwise_identity_nan_when_no_overlap():
    frame = alignment_from_sequences(
        ["s1", "s2"],
        ["A-", "-T"],
    )
    result = pairwise_identity_matrix(frame)
    assert math.isnan(result.matrix[0][1])
    assert result.coverage[0][1] == 0


def test_constant_columns_handles_all_gap_column():
    frame = alignment_from_sequences(
        ["s1", "s2"],
        ["-A", "-T"],
    )
    mask = constant_columns(frame)
    assert mask == [False, False]


def test_percentage_identity_and_rcv():
    frame = _example_alignment()
    pct = percentage_identity(frame)
    assert 50.0 < pct < 100.0
    rcv = relative_composition_variability(frame)
    assert rcv >= 0.0


def test_percentage_identity_nan_for_single_sequence():
    frame = alignment_from_sequences(ids=["only"], sequences=["AAAA"], alphabet=["A"], metadata={})
    value = percentage_identity(frame)
    assert math.isnan(value)
