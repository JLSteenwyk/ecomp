"""Alignment metrics that operate on :class:`AlignmentFrame` objects."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

from .._compat import zip_strict
from ..io import AlignmentFrame

DEFAULT_GAP_CHARACTERS = {"-"}

_IUPAC_NUCLEOTIDE_CODES = {
    frozenset({"A"}): "A",
    frozenset({"C"}): "C",
    frozenset({"G"}): "G",
    frozenset({"T"}): "T",
    frozenset({"U"}): "T",
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"C", "U"}): "Y",
    frozenset({"A", "C"}): "M",
    frozenset({"G", "T"}): "K",
    frozenset({"G", "U"}): "K",
    frozenset({"A", "T"}): "W",
    frozenset({"A", "U"}): "W",
    frozenset({"C", "G"}): "S",
    frozenset({"A", "C", "G"}): "V",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "U"}): "H",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "G", "U"}): "D",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"C", "G", "U"}): "B",
    frozenset({"A", "C", "G", "T"}): "N",
    frozenset({"A", "C", "G", "U"}): "N",
}


def column_base_counts(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
    include_gaps: bool = False,
) -> list[Counter[str]]:
    """Return per-column residue counts.

    ``include_gaps`` controls whether the returned tallies contain gap symbols.
    When ``False`` (default) the counts only reflect non-gap residues.
    """

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    counts: list[Counter[str]] = []
    for column in zip_strict(*frame.sequences):
        if include_gaps:
            counts.append(Counter(column))
        else:
            counts.append(Counter(char for char in column if char not in gap_chars))
    return counts


def majority_rule_consensus(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
    tie_breaker: Sequence[str] | None = None,
    gap_placeholder: str = "-",
) -> str:
    """Return a majority-rule consensus string for *frame*.

    ``tie_breaker`` provides a deterministic ordering of residues when multiple
    symbols share the top frequency.  It defaults to the frame alphabet.
    ``gap_placeholder`` is emitted when a column contains only gap characters.
    """

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    ordering = list(tie_breaker or frame.alphabet)
    ordering_index = {char: idx for idx, char in enumerate(ordering)}
    consensus: list[str] = []

    for column in zip_strict(*frame.sequences):
        counter = Counter(char for char in column if char not in gap_chars)
        if not counter:
            consensus.append(gap_placeholder)
            continue
        # Determine residues with top frequency
        top_count = max(counter.values())
        top_residues = [res for res, count in counter.items() if count == top_count]
        if len(top_residues) == 1:
            consensus.append(top_residues[0])
            continue

        residue_set = frozenset(top_residues)
        # Normalise uracil representation for mapping
        if "U" in residue_set and "T" not in residue_set:
            residue_set = frozenset({("T" if res == "U" else res) for res in residue_set})

        if residue_set in _IUPAC_NUCLEOTIDE_CODES:
            consensus.append(_IUPAC_NUCLEOTIDE_CODES[residue_set])
            continue

        if "N" in residue_set:
            consensus.append("N")
            continue

        # Fall back to alphabetical ordering for protein or unrecognised mixtures
        consensus.append("X")
    return "".join(consensus)


def column_gap_fraction(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> list[float]:
    """Compute the fraction of gap symbols per column."""

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    total = frame.num_sequences
    fractions: list[float] = []
    for column in zip_strict(*frame.sequences):
        gap_count = sum(1 for char in column if char in gap_chars)
        fractions.append(gap_count / total)
    return fractions


def column_shannon_entropy(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> list[float]:
    """Return per-column Shannon entropy (base 2) ignoring gaps."""

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    entropies: list[float] = []
    for column_counts in column_base_counts(frame, gap_characters=gap_chars, include_gaps=False):
        total = sum(column_counts.values())
        if total == 0:
            entropies.append(0.0)
            continue
        entropy = 0.0
        for count in column_counts.values():
            if count == 0:
                continue
            p = count / total
            entropy -= p * math.log2(p)
        entropies.append(entropy)
    return entropies


def parsimony_informative_columns(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> list[bool]:
    """Flag columns that are parsimony-informative.

    A column is parsimony-informative when at least two residues (ignoring gaps)
    occur with a frequency of two or more.
    """

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    informative: list[bool] = []
    for column_counts in column_base_counts(frame, gap_characters=gap_chars, include_gaps=False):
        qualifying = sum(1 for count in column_counts.values() if count >= 2)
        informative.append(qualifying >= 2)
    return informative


def constant_columns(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> list[bool]:
    """Return a mask indicating which columns are constant (ignoring gaps)."""

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    constant_mask: list[bool] = []
    for column_counts in column_base_counts(frame, gap_characters=gap_chars, include_gaps=False):
        non_zero = [count for count in column_counts.values() if count > 0]
        constant_mask.append(len(non_zero) == 1)
    return constant_mask


@dataclass(frozen=True)
class PairwiseIdentityResult:
    matrix: list[list[float]]
    coverage: list[list[int]]


def pairwise_identity_matrix(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> PairwiseIdentityResult:
    """Compute pairwise identity between sequences.

    The matrix contains values in ``[0, 1]``; entries are ``NaN`` when two
    sequences share no overlapping non-gap positions.  Diagonal entries are ``1``.
    ``coverage`` reports the number of positions compared for each pair.
    """

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    sequences = frame.sequences
    n = len(sequences)
    matches = [[0] * n for _ in range(n)]
    coverage = [[0] * n for _ in range(n)]

    columns = zip_strict(*sequences)
    for column in columns:
        for i in range(n):
            residue_i = column[i]
            if residue_i in gap_chars:
                continue
            for j in range(i + 1, n):
                residue_j = column[j]
                if residue_j in gap_chars:
                    continue
                coverage[i][j] += 1
                if residue_i == residue_j:
                    matches[i][j] += 1

    matrix = [[1.0 if i == j else float("nan") for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if coverage[i][j]:
                identity = matches[i][j] / coverage[i][j]
                matrix[i][j] = matrix[j][i] = identity
            else:
                matrix[i][j] = matrix[j][i] = float("nan")

    # Mirror lower-triangular counts into the symmetric position.
    for i in range(n):
        for j in range(i):
            coverage[i][j] = coverage[j][i]
            matches[i][j] = matches[j][i]

    return PairwiseIdentityResult(matrix=matrix, coverage=coverage)


def alignment_length_excluding_gaps(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> int:
    """Number of columns that contain at least one non-gap character."""

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    return sum(
        1
        for column in zip_strict(*frame.sequences)
        if any(char not in gap_chars for char in column)
    )


def variable_site_count(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> int:
    """Count columns that exhibit more than one residue ignoring gaps."""

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    variable = 0
    for column in zip_strict(*frame.sequences):
        residues = {char for char in column if char not in gap_chars}
        if len(residues) >= 2:
            variable += 1
    return variable


def alignment_length(frame: AlignmentFrame) -> int:
    """Total alignment length in columns (gaps retained)."""

    return frame.alignment_length


def parsimony_informative_site_count(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> int:
    """Total number of parsimony-informative columns."""

    return sum(parsimony_informative_columns(frame, gap_characters=gap_characters))


def percentage_identity(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> float:
    """Mean pairwise identity across all sequence pairs (percentage)."""

    result = pairwise_identity_matrix(frame, gap_characters=gap_characters)
    n = len(frame.ids)
    if n < 2:
        return float("nan")
    identities: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            value = result.matrix[i][j]
            if math.isnan(value):
                continue
            identities.append(value)
    if not identities:
        return float("nan")
    return sum(identities) / len(identities) * 100.0


def relative_composition_variability(
    frame: AlignmentFrame,
    *,
    gap_characters: Iterable[str] | None = None,
) -> float:
    """Relative composition variability (RCV) expressed as a percentage."""

    gap_chars = set(gap_characters or DEFAULT_GAP_CHARACTERS)
    alphabet = sorted({char for seq in frame.sequences for char in seq if char not in gap_chars})
    if not alphabet:
        return 0.0

    compositions: list[dict[str, float]] = []
    totals: list[int] = []
    for seq in frame.sequences:
        counts: dict[str, float] = {char: 0.0 for char in alphabet}
        total = 0
        for char in seq:
            if char in gap_chars:
                continue
            counts[char] += 1.0
            total += 1
        totals.append(total)
        if total:
            for char in alphabet:
                counts[char] /= total
        compositions.append(counts)

    if not any(totals):
        return 0.0

    mean_composition = {char: sum(comp[char] for comp in compositions) / len(compositions) for char in alphabet}
    rcv = 0.0
    for comp in compositions:
        for char in alphabet:
            rcv += abs(comp[char] - mean_composition[char])
    return (rcv / len(compositions)) * 100.0
