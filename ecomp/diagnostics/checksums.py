"""Checksum utilities used to validate round-trip fidelity."""

from __future__ import annotations

import hashlib
from typing import Iterable


def alignment_checksum(sequences: Iterable[str]) -> str:
    """Return a SHA256 checksum for the provided sequences."""

    digest = hashlib.sha256()
    for seq in sequences:
        digest.update(seq.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()
