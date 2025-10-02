"""Consensus computation helpers for the compression pipeline."""

from __future__ import annotations

from collections import Counter
from typing import Iterator

from .._compat import dataclass, zip_strict

from ..io import AlignmentFrame


@dataclass(slots=True)
class ColumnProfile:
    """Consensus profile for a single alignment column."""

    consensus: str
    deviations: tuple[tuple[int, str], ...]

    def equivalent_key(self) -> tuple[str, tuple[tuple[int, str], ...]]:
        return self.consensus, self.deviations


def iter_column_profiles(frame: AlignmentFrame) -> Iterator[ColumnProfile]:
    """Yield :class:`ColumnProfile` objects for each column in *frame*."""

    sequences = frame.sequences
    if not sequences:
        return

    for column in zip_strict(*sequences):
        counts = Counter(column)
        consensus_char = max(counts.items(), key=lambda item: (item[1], -ord(item[0])))[0]
        deviations = tuple(
            (seq_idx, residue)
            for seq_idx, residue in enumerate(column)
            if residue != consensus_char
        )
        yield ColumnProfile(consensus=consensus_char, deviations=deviations)


def collect_column_profiles(frame: AlignmentFrame) -> list[ColumnProfile]:
    """Compute column profiles eagerly and return them as a list."""

    return list(iter_column_profiles(frame))
