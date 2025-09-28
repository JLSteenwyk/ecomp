import random

import pytest

from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.compression import pipeline
from evolutionary_compression.compression.pipeline import compress_alignment, decompress_alignment


def test_compress_decompress_preserves_sequences():
    frame = alignment_from_sequences(
        ids=["seq1", "seq2", "seq3"],
        sequences=[
            "ACGTACGT",
            "ACGTACGA",
            "ACGTTCGT",
        ],
    )

    compressed = compress_alignment(frame)
    reconstructed = decompress_alignment(compressed.payload, compressed.metadata)

    assert reconstructed.ids == frame.ids
    assert reconstructed.sequences == frame.sequences
    assert compressed.metadata["max_run_length"] >= 1
    assert compressed.metadata["columns_with_deviations"] >= 0
    assert compressed.metadata["bitmask_bytes"] == 1
    assert compressed.metadata["bits_per_symbol"] == 2
    assert compressed.metadata["payload_encoded_bytes"] <= compressed.metadata["payload_raw_bytes"]
    assert compressed.metadata["payload_encoding"] in {"raw", "zlib"}
    assert compressed.metadata["codec"] == "ecomp"


def test_decompress_raises_when_metadata_sequence_ids_mismatch():
    frame = alignment_from_sequences(ids=["seq1"], sequences=["AAAA"])
    compressed = compress_alignment(frame)
    bad_metadata = dict(compressed.metadata)
    bad_metadata["sequence_ids"] = []

    try:
        decompress_alignment(compressed.payload, bad_metadata)
    except ValueError as exc:
        assert "Metadata sequence count" in str(exc)
    else:  # pragma: no cover - defensive
        raise AssertionError("Expected ValueError for mismatched metadata")


def test_decompress_raises_when_bitmask_bytes_missing():
    frame = alignment_from_sequences(ids=["seq1"], sequences=["AAAA"])
    compressed = compress_alignment(frame)
    bad_metadata = dict(compressed.metadata)
    bad_metadata.pop("bitmask_bytes", None)

    try:
        decompress_alignment(compressed.payload, bad_metadata)
    except ValueError as exc:
        assert "bitmask_bytes" in str(exc)
    else:  # pragma: no cover - defensive
        raise AssertionError("Expected ValueError for missing bitmask metadata")


def test_decompress_allows_missing_bits_per_symbol():
    frame = alignment_from_sequences(ids=["seq1"], sequences=["AAAA"])
    compressed = compress_alignment(frame)
    metadata = dict(compressed.metadata)
    metadata.pop("bits_per_symbol", None)

    reconstructed = decompress_alignment(compressed.payload, metadata)
    assert reconstructed.sequences == frame.sequences


def test_fallback_to_gzip_when_smaller():
    random.seed(0)
    alphabet = "ACGT"
    ids = [f"seq{i}" for i in range(6)]
    sequences = [
        "".join(random.choice(alphabet) for _ in range(200))
        for _ in ids
    ]
    frame = alignment_from_sequences(ids=ids, sequences=sequences)
    compressed = compress_alignment(frame)
    assert compressed.metadata.get("fallback", {}).get("type") == "gzip"
    restored = decompress_alignment(compressed.payload, compressed.metadata)
    assert restored.sequences == frame.sequences


def test_similarity_order_prefers_nonbaseline():
    ids = ["s1", "s3", "s2", "s4"]
    sequences = [
        "AAAAAAAA",
        "TTTTTTTT",
        "AAAATTTT",
        "TTTTAAAA",
    ]
    frame = alignment_from_sequences(ids=ids, sequences=sequences)

    reordered, permutation, label = pipeline._compute_similarity_order(frame)  # type: ignore[attr-defined]
    assert label.startswith("auto-")
    assert permutation != list(range(len(ids)))

    compressed = compress_alignment(frame)
    if compressed.metadata.get("sequence_permutation") is not None:
        assert compressed.metadata["sequence_permutation"] == permutation
    assert compressed.metadata.get("ordering_strategy", "").startswith("auto-")

    restored = decompress_alignment(compressed.payload, compressed.metadata)
    assert restored.sequences == frame.sequences


def test_tree_guided_sequence_ordering():
    ids = ["taxC", "taxA", "taxD", "taxB"]
    sequences = [
        "ACGTACGT",
        "ACGTACGA",
        "ACGTACGG",
        "ACGTACGC",
    ]
    newick = "((taxA:0.1,taxB:0.1):0.2,(taxC:0.1,taxD:0.1):0.2);"
    frame = alignment_from_sequences(
        ids=ids,
        sequences=sequences,
        metadata={"tree_newick": newick},
    )

    reordered, permutation, label = pipeline._compute_similarity_order(frame)  # type: ignore[attr-defined]
    assert label == "tree"
    assert permutation == [1, 3, 0, 2]
    assert reordered.ids == ["taxA", "taxB", "taxC", "taxD"]

    compressed = compress_alignment(frame)
    assert compressed.metadata.get("ordering_strategy", "").startswith("auto") or compressed.metadata.get("ordering_strategy") == "tree"

    restored = decompress_alignment(compressed.payload, compressed.metadata)
    assert restored.ids == ids
    assert restored.sequences == sequences


def test_decompress_alignment_raises_on_column_length_mismatch():
    frame = alignment_from_sequences(
        ids=["s1", "s2"],
        sequences=["ACGT", "ACGA"],
    )
    compressed = compress_alignment(frame)
    bad_metadata = dict(compressed.metadata)
    bad_metadata["alignment_length"] = bad_metadata["alignment_length"] + 1

    with pytest.raises(ValueError):
        decompress_alignment(compressed.payload, bad_metadata)

