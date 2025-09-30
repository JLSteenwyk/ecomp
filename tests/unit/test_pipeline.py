import os
import random
import lzma

import math

import pytest

from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.compression import pipeline
from evolutionary_compression.compression.pipeline import compress_alignment, decompress_alignment
from evolutionary_compression.diagnostics.checksums import alignment_checksum
from evolutionary_compression.compression.consensus import collect_column_profiles
from evolutionary_compression.compression.rle import collect_run_length_blocks
from evolutionary_compression.compression.encoding import encode_blocks


def _build_raw_payload_and_metadata(frame):
    column_profiles = collect_column_profiles(frame)
    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, math.ceil(math.log2(max(len(alphabet), 1))))
    run_length_blocks = collect_run_length_blocks(
        column_profiles, frame.num_sequences, symbol_lookup, bits_per_symbol
    )
    bitmask_bytes = (frame.num_sequences + 7) // 8
    run_length_payload = encode_blocks(
        run_length_blocks, bitmask_bytes, bits_per_symbol, alphabet
    )
    seq_block = pipeline._encode_sequence_ids(frame.ids)
    raw_payload = seq_block + run_length_payload
    deviation_columns = sum(1 for block in column_profiles if block.deviations)
    metadata = {
        "format_version": pipeline.FORMAT_VERSION,
        "codec": "ecomp",
        "num_sequences": frame.num_sequences,
        "alignment_length": frame.alignment_length,
        "alphabet": frame.alphabet,
        "source_format": frame.metadata.get("source_format", "fasta"),
        "checksum_sha256": alignment_checksum(frame.sequences),
        "run_length_blocks": len(run_length_blocks),
        "max_run_length": max((block.run_length for block in run_length_blocks), default=0),
        "columns_with_deviations": deviation_columns,
        "bitmask_bytes": bitmask_bytes,
        "bits_per_symbol": bits_per_symbol,
        "payload_encoding": "raw",
        "payload_encoded_bytes": len(raw_payload),
        "payload_raw_bytes": len(raw_payload),
        "sequence_id_codec": "inline",
        "ordering_strategy": "baseline",
    }
    return raw_payload, metadata


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
    perm_meta = compressed.metadata.get("sequence_permutation")
    if perm_meta is not None:
        assert perm_meta.get("encoding") == "payload"
        assert perm_meta.get("length", 0) > 0
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


def test_tree_guided_order_skipped_when_gap_heavy():
    ids = ["tax1", "tax2", "tax3", "tax4"]
    sequences = [
        "-A-A-A",
        "A-A-A-",
        "-A-A-A",
        "A-A-A-",
    ]
    newick = "((tax1:0.1,tax2:0.1):0.2,(tax3:0.1,tax4:0.1):0.2);"
    frame = alignment_from_sequences(
        ids=ids,
        sequences=sequences,
        metadata={"tree_newick": newick},
    )

    compressed = compress_alignment(frame)
    assert compressed.metadata.get("ordering_strategy") != "tree"
    restored = decompress_alignment(compressed.payload, compressed.metadata)
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


def test_decompress_alignment_requires_bits_per_symbol_without_alphabet():
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["AAAA", "AAAT"])
    raw_payload, metadata = _build_raw_payload_and_metadata(frame)
    bad_meta = dict(metadata)
    bad_meta.pop("bits_per_symbol", None)
    bad_meta["alphabet"] = []
    with pytest.raises(ValueError):
        decompress_alignment(raw_payload, bad_meta)


def test_decompress_alignment_rejects_unknown_payload_encoding():
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["AAAA", "AAAT"])
    raw_payload, metadata = _build_raw_payload_and_metadata(frame)
    bad_meta = dict(metadata)
    bad_meta["payload_encoding"] = "unknown"
    with pytest.raises(ValueError):
        decompress_alignment(raw_payload, bad_meta)


def test_decompress_alignment_supports_xz_encoding():
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["AAAA", "AAAT"])
    raw_payload, metadata = _build_raw_payload_and_metadata(frame)
    metadata_xz = dict(metadata)
    payload_xz = lzma.compress(raw_payload)
    metadata_xz["payload_encoding"] = "xz"
    metadata_xz["payload_encoded_bytes"] = len(payload_xz)
    restored = decompress_alignment(payload_xz, metadata_xz)
    assert restored.sequences == frame.sequences
    assert restored.ids == frame.ids


def test_decompress_alignment_handles_unknown_fallback_type():
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["AAAA", "AAAT"])
    raw_payload, metadata = _build_raw_payload_and_metadata(frame)
    bad_meta = dict(metadata)
    bad_meta["fallback"] = {"type": "zip"}
    with pytest.raises(ValueError):
        decompress_alignment(raw_payload, bad_meta)


def test_decompress_alignment_raises_when_columns_exceed_expected():
    frame = alignment_from_sequences(ids=["s1", "s2"], sequences=["AAAA", "AAAT"])
    raw_payload, metadata = _build_raw_payload_and_metadata(frame)
    bad_meta = dict(metadata)
    bad_meta["alignment_length"] = metadata["alignment_length"] - 1
    with pytest.raises(ValueError):
        decompress_alignment(raw_payload, bad_meta)


def test_choose_order_respects_env_override(monkeypatch):
    monkeypatch.setenv("ECOMP_SEQUENCE_ORDER", "mst")
    ids = ["s1", "s2", "s3"]
    sequences = [
        "AAAA",
        "TTTT",
        "AAAT",
    ]
    frame = alignment_from_sequences(ids=ids, sequences=sequences)
    _, permutation, label = pipeline._compute_similarity_order(frame)  # type: ignore[attr-defined]
    assert label == "mst"
    assert len(permutation) == len(ids)
