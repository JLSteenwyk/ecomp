"""Utilities for reading and writing multiple sequence alignments."""

from __future__ import annotations

from dataclasses import field
from pathlib import Path
from typing import Iterable, Sequence

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .config import DEFAULT_OUTPUT_FORMAT, detect_format_from_suffix
from ._compat import dataclass, zip_strict


@dataclass(slots=True)
class AlignmentFrame:
    """In-memory representation of an alignment used throughout the pipeline."""

    ids: list[str]
    sequences: list[str]
    alphabet: list[str]
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.ids:
            raise ValueError("Alignment must contain at least one sequence")
        if len(self.ids) != len(self.sequences):
            raise ValueError("Sequence IDs and sequences length mismatch")
        lengths = {len(seq) for seq in self.sequences}
        if len(lengths) != 1:
            raise ValueError("All sequences must be the same length")

    @property
    def num_sequences(self) -> int:
        return len(self.sequences)

    @property
    def alignment_length(self) -> int:
        return len(self.sequences[0]) if self.sequences else 0

    def alphabet_string(self) -> str:
        """Return the alphabet as a deterministic string representation."""

        return "".join(sorted(set(self.alphabet)))


def _to_alignment(frame: AlignmentFrame) -> MultipleSeqAlignment:
    records = [
        SeqRecord(Seq(seq), id=seq_id, description=str(frame.metadata.get("description", "")))
        for seq_id, seq in zip_strict(frame.ids, frame.sequences)
    ]
    return MultipleSeqAlignment(records)


def read_alignment(path: str | Path, fmt: str | None = None) -> AlignmentFrame:
    """Load an alignment file into an :class:`AlignmentFrame`."""

    path = Path(path)
    if fmt is None:
        fmt = detect_format_from_suffix(str(path)) or DEFAULT_OUTPUT_FORMAT
    alignment = AlignIO.read(str(path), fmt)
    ids = [record.id for record in alignment]
    sequences = [str(record.seq) for record in alignment]
    alphabet = sorted({char for seq in sequences for char in seq})
    metadata = {"source_path": str(path), "source_format": fmt}
    return AlignmentFrame(ids=ids, sequences=sequences, alphabet=alphabet, metadata=metadata)


def write_alignment(
    frame: AlignmentFrame,
    path: str | Path,
    fmt: str | None = None,
    wrap: int | None = None,
) -> Path:
    """Persist an alignment to disk in the requested *fmt*."""

    path = Path(path)
    target_fmt = fmt or detect_format_from_suffix(str(path)) or DEFAULT_OUTPUT_FORMAT
    alignment = _to_alignment(frame)
    if wrap is not None:
        AlignIO.write(alignment, str(path), target_fmt, wrap=wrap)
    else:
        AlignIO.write(alignment, str(path), target_fmt)
    return path


def alignment_from_sequences(
    ids: Sequence[str],
    sequences: Sequence[str],
    alphabet: Iterable[str] | None = None,
    metadata: dict[str, object] | None = None,
) -> AlignmentFrame:
    """Construct an :class:`AlignmentFrame` directly from sequences."""

    alphabet_values = alphabet or {char for seq in sequences for char in seq}
    return AlignmentFrame(
        ids=list(ids),
        sequences=list(sequences),
        alphabet=sorted(set(alphabet_values)),
        metadata=dict(metadata or {}),
    )
