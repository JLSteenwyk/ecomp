from evolutionary_compression.config import (
    DEFAULT_OUTPUT_FORMAT,
    METADATA_SUFFIX,
    SUPPORTED_INPUT_FORMATS,
    detect_format_from_suffix,
)


def test_detect_format_from_suffix_matches_known_extensions():
    assert detect_format_from_suffix("example.fasta") == "fasta"
    assert detect_format_from_suffix("example.PHY") == "phylip"


def test_detect_format_from_suffix_returns_none_for_unknown_extension():
    assert detect_format_from_suffix("example.unknown") is None
    assert detect_format_from_suffix("example") is None


def test_config_constants_align_with_expectations():
    assert METADATA_SUFFIX == ".json"
    assert DEFAULT_OUTPUT_FORMAT in SUPPORTED_INPUT_FORMATS
