import pytest

from evolutionary_compression.cli import _verify_checksum, build_parser, main
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


def test_build_parser_registers_subcommands():
    parser = build_parser()
    # The parser stores subcommands under the first subparsers action.
    subparsers_actions = [action for action in parser._subparsers._group_actions]
    assert subparsers_actions, "Expected at least one subparsers action"
    choices = subparsers_actions[0].choices
    assert {"compress", "decompress", "inspect"} <= set(choices)


def test_main_errors_on_missing_alignment(tmp_path):
    missing = tmp_path / "nope.fasta"
    with pytest.raises(SystemExit) as exc:
        main(["compress", str(missing)])
    assert "Alignment not found" in str(exc.value)


def test_main_errors_on_missing_archive(tmp_path):
    missing = tmp_path / "missing.ecomp"
    with pytest.raises(SystemExit) as exc:
        main(["decompress", str(missing)])
    assert "Archive not found" in str(exc.value)
