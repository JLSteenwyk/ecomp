import pytest

from evolutionary_compression.cli import _verify_checksum
from evolutionary_compression.diagnostics.checksums import alignment_checksum


def test_verify_checksum_passes():
    sequences = ["ACGT", "ACGA"]
    checksum = alignment_checksum(sequences)
    metadata = {"checksum_sha256": checksum}
    _verify_checksum(sequences, metadata)


def test_verify_checksum_raises_on_mismatch():
    sequences = ["ACGT"]
    metadata = {"checksum_sha256": alignment_checksum(["AAAA"])}
    with pytest.raises(SystemExit):
        _verify_checksum(sequences, metadata)


def test_verify_checksum_ignored_when_missing():
    _verify_checksum(["ACGT"], {})
