"""High-level compression and decompression orchestration."""

from __future__ import annotations

import base64
import gzip
import io
import lzma
import math
import os
import time
import zlib
from typing import Any, Iterable, List, Sequence, Tuple

from .._compat import dataclass, zip_strict

from ..diagnostics.checksums import alignment_checksum
from ..io import AlignmentFrame, alignment_from_sequences
from ..config import FORMAT_VERSION
from .consensus import collect_column_profiles
from .encoding import decode_blocks, encode_blocks
from .rle import RunLengthBlock, collect_run_length_blocks

try:  # pragma: no cover - optional dependency
    import zstandard as zstd

    _ZSTD_COMPRESSOR = zstd.ZstdCompressor(level=5)
    _ZSTD_DECOMPRESSOR = zstd.ZstdDecompressor()
except ModuleNotFoundError:  # pragma: no cover - environment without zstd
    zstd = None
    _ZSTD_COMPRESSOR = None
    _ZSTD_DECOMPRESSOR = None

_SEQ_ID_MAGIC = b"ECID"
_SEQ_ID_VERSION = 2
_SAMPLE_CAP = 256


_PERM_MAGIC = b"ECPE"
_PERM_VERSION = 1
_WIDTH_TO_CODE = {1: 0, 2: 1, 4: 2}
_CODE_TO_WIDTH = {value: key for key, value in _WIDTH_TO_CODE.items()}


def _build_permutation_chunk(
    permutation: Sequence[int],
) -> tuple[bytes | None, dict[str, Any] | None]:
    """Encode *permutation* as a binary chunk stored alongside the payload."""

    size = len(permutation)
    if size == 0:
        return None, None

    max_value = max(permutation)
    if max_value < 256:
        width = 1
    elif max_value < 65536:
        width = 2
    else:
        width = 4
    width_code = _WIDTH_TO_CODE[width]

    raw = bytearray()
    for value in permutation:
        raw.extend(int(value).to_bytes(width, "little", signed=False))

    compressed = zlib.compress(bytes(raw), level=9)
    if len(compressed) + 8 < len(raw):
        payload = compressed
        compression_flag = 1
    else:
        payload = bytes(raw)
        compression_flag = 0

    chunk = bytearray()
    chunk.extend(_PERM_MAGIC)
    chunk.append(_PERM_VERSION)
    chunk.append((width_code << 1) | compression_flag)
    chunk.extend(_encode_varint(size))
    chunk.extend(_encode_varint(len(payload)))
    chunk.extend(payload)

    chunk_bytes = bytes(chunk)
    metadata = {"encoding": "payload", "length": len(chunk_bytes)}
    return chunk_bytes, metadata


def _extract_permutation_chunk(
    payload_data: bytes, metadata: dict[str, Any]
) -> tuple[bytes, list[int]]:
    """Remove and decode the permutation chunk from *payload_data*."""

    length = metadata.get("length")
    if not isinstance(length, int) or length <= 0:
        raise ValueError("Invalid permutation metadata length")
    if len(payload_data) < length:
        raise ValueError("Permutation chunk exceeds payload size")

    chunk = memoryview(payload_data[:length])
    cursor = 0
    if chunk[cursor : cursor + 4].tobytes() != _PERM_MAGIC:
        raise ValueError("Permutation chunk missing magic header")
    cursor += 4

    version = chunk[cursor]
    cursor += 1
    if version != _PERM_VERSION:
        raise ValueError(f"Unsupported permutation chunk version: {version}")

    flags = chunk[cursor]
    cursor += 1
    compression_flag = flags & 0x01
    width_code = (flags >> 1) & 0x03
    try:
        width = _CODE_TO_WIDTH[width_code]
    except KeyError as exc:
        raise ValueError(f"Unsupported permutation width code: {width_code}") from exc

    size, cursor = _decode_varint(chunk, cursor)
    payload_len, cursor = _decode_varint(chunk, cursor)
    end = cursor + payload_len
    if end != length:
        raise ValueError("Permutation chunk length mismatch")

    payload = bytes(chunk[cursor:end])
    if compression_flag:
        payload = zlib.decompress(payload)

    if len(payload) != size * width:
        raise ValueError("Permutation payload size mismatch")

    permutation = [
        int.from_bytes(payload[i : i + width], "little", signed=False)
        for i in range(0, len(payload), width)
    ]

    return payload_data[length:], permutation


def _decode_permutation(value: Any) -> list[int]:
    """Decode legacy (pre-payload) permutation metadata formats."""

    if value is None:
        return []
    if isinstance(value, list):
        return [int(v) for v in value]
    if not isinstance(value, dict):
        raise ValueError("Unsupported permutation metadata format")

    if value.get("encoding") == "payload":  # handled elsewhere
        raise ValueError("Payload-encoded permutation should be extracted before decoding")

    version = value.get("version", 0)
    if version != _PERM_VERSION:
        raise ValueError(f"Unsupported permutation metadata version: {version}")

    size = int(value.get("size", 0))
    dtype = value.get("dtype")
    compression = value.get("compression", "none")
    data = value.get("data", "")

    payload = base64.b64decode(data.encode("ascii")) if data else b""
    if compression == "zlib" and payload:
        payload = zlib.decompress(payload)
    elif compression not in {"none", "zlib"}:
        raise ValueError(f"Unsupported permutation compression: {compression}")

    width_lookup = {"uint8": 1, "uint16": 2, "uint32": 4}
    try:
        width = width_lookup[dtype]
    except KeyError as exc:
        raise ValueError(f"Unsupported permutation dtype: {dtype}") from exc

    expected = size * width
    if len(payload) != expected:
        raise ValueError("Permutation payload length mismatch")

    permutation = [
        int.from_bytes(payload[i : i + width], "little", signed=False)
        for i in range(0, expected, width)
    ]
    return permutation


_GAP_CHARACTERS = {"-", ".", "?", "*", "X", "x"}


def _compute_alignment_stats(frame: AlignmentFrame) -> dict[str, float] | None:
    num_sequences = frame.num_sequences
    alignment_length = frame.alignment_length
    if num_sequences == 0 or alignment_length == 0:
        return None

    total_pairs = num_sequences * (num_sequences - 1) / 2
    if total_pairs == 0:
        total_pairs = 1.0

    gap_free_columns = 0
    variable_columns = 0
    identity_sum = 0.0
    identity_sumsq = 0.0

    sequences = frame.sequences

    for column_index in range(alignment_length):
        column_counts: dict[str, int] = {}
        has_gap = False
        for seq in sequences:
            char = seq[column_index]
            if char in _GAP_CHARACTERS:
                has_gap = True
            column_counts[char] = column_counts.get(char, 0) + 1

        non_gap_counts = [count for char, count in column_counts.items() if char not in _GAP_CHARACTERS]
        if not has_gap:
            gap_free_columns += 1
        if len([count for count in non_gap_counts if count > 0]) > 1:
            variable_columns += 1

        identical_pairs = sum(count * (count - 1) / 2 for count in non_gap_counts)
        identity = identical_pairs / total_pairs
        identity_sum += identity
        identity_sumsq += identity * identity

    pair_identity_mean = identity_sum / alignment_length
    variance = max(identity_sumsq / alignment_length - pair_identity_mean * pair_identity_mean, 0.0)
    pair_identity_sd = math.sqrt(variance)

    alphabet = sorted({char for seq in sequences for char in seq})
    alphabet_index = {char: idx for idx, char in enumerate(alphabet)}
    totals = [0.0 for _ in alphabet]
    per_sequence_counts: list[list[int]] = []
    for seq in sequences:
        counts = [0] * len(alphabet)
        for char in seq:
            counts[alphabet_index[char]] += 1
        per_sequence_counts.append(counts)
        for idx, value in enumerate(counts):
            totals[idx] += value

    averages = [value / num_sequences for value in totals]
    rcv_accumulator = 0.0
    for counts in per_sequence_counts:
        deviation = sum(abs(count - avg) for count, avg in zip(counts, averages))
        rcv_accumulator += deviation
    rcv_total = rcv_accumulator / (num_sequences * alignment_length)

    return {
        "gap_free_pct": 100.0 * gap_free_columns / alignment_length,
        "variable_pct": 100.0 * variable_columns / alignment_length,
        "pair_identity_mean": pair_identity_mean,
        "pair_identity_sd": pair_identity_sd,
        "rcv_total": rcv_total,
    }


def _should_use_tree(stats: dict[str, float] | None) -> bool:
    if not stats:
        return False
    if stats.get("gap_free_pct", 0.0) < 5.0:
        return False
    if stats.get("rcv_total", 1.0) > 0.25:
        return False
    if stats.get("pair_identity_mean", 0.0) < 0.45:
        return False
    if stats.get("pair_identity_sd", 0.0) > 0.35:
        return False
    return True


def _encode_varint(value: int) -> bytes:
    if value < 0:
        raise ValueError("varint cannot encode negative values")
    out = bytearray()
    while True:
        byte = value & 0x7F
        value >>= 7
        if value:
            out.append(byte | 0x80)
        else:
            out.append(byte)
            break
    return bytes(out)


def _decode_varint(buffer: memoryview, cursor: int) -> tuple[int, int]:
    shift = 0
    result = 0
    while True:
        if cursor >= len(buffer):
            raise ValueError("Unexpected end of data while decoding varint")
        byte = buffer[cursor]
        cursor += 1
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return result, cursor
        shift += 7
        if shift > 56:
            raise ValueError("Varint too long")


def _encode_sequence_ids(ids: Iterable[str]) -> bytes:
    names = [name.encode("utf-8") for name in ids]
    plain = bytearray()
    plain.extend(_encode_varint(len(names)))
    for name_bytes in names:
        plain.extend(_encode_varint(len(name_bytes)))
        plain.extend(name_bytes)

    mode = 0
    payload_bytes = bytes(plain)

    if _ZSTD_COMPRESSOR is not None:
        compressed = _ZSTD_COMPRESSOR.compress(payload_bytes)
        if len(compressed) + 1 < len(payload_bytes):
            mode = 1
            payload_bytes = compressed
    else:
        compressed = None

    if mode == 0:
        zlib_compressed = zlib.compress(payload_bytes, level=9)
        if len(zlib_compressed) + 1 < len(payload_bytes):
            mode = 2
            payload_bytes = zlib_compressed

    encoded = bytearray()
    encoded.append(mode)
    encoded.extend(payload_bytes)

    payload = bytearray(_SEQ_ID_MAGIC)
    payload.append(_SEQ_ID_VERSION)
    payload.extend(_encode_varint(len(encoded)))
    payload.extend(encoded)
    return bytes(payload)


def _decode_sequence_ids(buffer: bytes) -> Tuple[list[str], bytes]:
    view = memoryview(buffer)
    if len(view) < len(_SEQ_ID_MAGIC) + 1:
        raise ValueError("Sequence ID block truncated")
    if not view.tobytes().startswith(_SEQ_ID_MAGIC):
        raise ValueError("Sequence ID magic missing")
    cursor = len(_SEQ_ID_MAGIC)
    version = view[cursor]
    cursor += 1
    block_length, cursor = _decode_varint(view, cursor)
    end = cursor + block_length
    if end > len(view):
        raise ValueError("Sequence ID block length exceeds payload size")

    if version == 1:
        count, cursor = _decode_varint(view, cursor)
        ids: list[str] = []
        for _ in range(count):
            name_len, cursor = _decode_varint(view, cursor)
            if cursor + name_len > end:
                raise ValueError("Sequence ID entry exceeds declared block length")
            name_bytes = view[cursor : cursor + name_len].tobytes()
            cursor += name_len
            ids.append(name_bytes.decode("utf-8"))
    elif version == 2:
        if cursor >= end:
            raise ValueError("Sequence ID block missing mode byte")
        mode = view[cursor]
        cursor += 1
        encoded = view[cursor:end].tobytes()
        if mode == 0:
            output = encoded
        elif mode == 1:
            if _ZSTD_DECOMPRESSOR is None:
                raise RuntimeError(
                    "Cannot decompress sequence IDs: zstandard module not available"
                )
            output = _ZSTD_DECOMPRESSOR.decompress(encoded)
        elif mode == 2:
            output = zlib.decompress(encoded)
        else:
            raise ValueError(f"Unsupported sequence ID compression mode: {mode}")

        data_view = memoryview(output)
        cursor_out = 0
        count, cursor_out = _decode_varint(data_view, cursor_out)
        ids: list[str] = []
        for _ in range(count):
            name_len, cursor_out = _decode_varint(data_view, cursor_out)
            if cursor_out + name_len > len(data_view):
                raise ValueError("Sequence ID entry exceeds decoded block length")
            name_bytes = data_view[cursor_out : cursor_out + name_len].tobytes()
            cursor_out += name_len
            ids.append(name_bytes.decode("utf-8"))
        cursor = end
    else:
        raise ValueError(f"Unsupported sequence ID block version: {version}")

    if cursor != end:
        raise ValueError("Sequence ID block contains trailing data")

    remaining = view[end:].tobytes()
    return ids, remaining


def _select_sample_indices(length: int, cap: int = _SAMPLE_CAP) -> list[int]:
    if length <= 0:
        return []
    if length <= cap:
        return list(range(length))
    step = max(1, length // cap)
    indices = list(range(0, length, step))
    if indices[-1] != length - 1:
        indices.append(length - 1)
    return indices


def _build_distance_matrix(sequences: list[str], sample_indices: list[int]) -> list[list[int]]:
    num_sequences = len(sequences)
    matrix = [[0] * num_sequences for _ in range(num_sequences)]
    if not sample_indices:
        return matrix
    for i in range(num_sequences):
        seq_i = sequences[i]
        for j in range(i + 1, num_sequences):
            seq_j = sequences[j]
            mismatches = 0
            for idx in sample_indices:
                if seq_i[idx] != seq_j[idx]:
                    mismatches += 1
            matrix[i][j] = mismatches
            matrix[j][i] = mismatches
    return matrix


def _mst_sequence_order(dist_matrix: list[list[int]]) -> list[int]:
    num_sequences = len(dist_matrix)
    if num_sequences == 0:
        return []
    visited = [False] * num_sequences
    key = [float("inf")] * num_sequences
    parent = [-1] * num_sequences
    key[0] = 0
    for _ in range(num_sequences):
        candidates = [idx for idx in range(num_sequences) if not visited[idx]]
        u = min(candidates, key=lambda idx: key[idx])
        visited[u] = True
        for v in range(num_sequences):
            if visited[v] or v == u:
                continue
            dist = dist_matrix[u][v]
            if dist < key[v]:
                key[v] = dist
                parent[v] = u

    adjacency: list[list[int]] = [[] for _ in range(num_sequences)]
    for child, par in enumerate(parent):
        if par >= 0:
            adjacency[par].append(child)

    order: list[int] = []
    stack = [0]
    while stack:
        node = stack.pop()
        order.append(node)
        children = sorted(adjacency[node], key=lambda idx: dist_matrix[node][idx], reverse=True)
        stack.extend(children)
    return order


def _greedy_sequence_order(dist_matrix: list[list[int]]) -> list[int]:
    num_sequences = len(dist_matrix)
    if num_sequences == 0:
        return []
    remaining = set(range(num_sequences))
    start = min(remaining, key=lambda idx: sum(dist_matrix[idx]))
    order = [start]
    remaining.remove(start)
    current = start
    while remaining:
        next_node = min(remaining, key=lambda idx: dist_matrix[current][idx])
        order.append(next_node)
        remaining.remove(next_node)
        current = next_node
    return order


def _order_cost(order: list[int], dist_matrix: list[list[int]]) -> int:
    if len(order) <= 1:
        return 0
    cost = 0
    for i in range(1, len(order)):
        cost += dist_matrix[order[i - 1]][order[i]]
    return cost


def _choose_order(dist_matrix: list[list[int]], candidates: list[tuple[str, list[int]]]) -> tuple[list[int], str]:
    env = os.environ.get("ECOMP_SEQUENCE_ORDER", "auto").lower()
    unique: dict[tuple[int, ...], tuple[str, list[int]]] = {}
    ordered_keys: list[tuple[int, ...]] = []
    for label, order in candidates:
        key = tuple(order)
        if key not in unique:
            unique[key] = (label, order)
            ordered_keys.append(key)

    explicit_labels = {"baseline", "mst", "greedy"}
    if env in explicit_labels:
        for key in ordered_keys:
            label, order = unique[key]
            if label == env:
                return order, label
        # fall back to baseline if explicit request is missing
        base_label, base_order = unique[ordered_keys[0]]
        return base_order, base_label

    best_label = None
    best_order: list[int] | None = None
    best_cost = None
    for key in ordered_keys:
        label, order = unique[key]
        cost = _order_cost(order, dist_matrix)
        if best_cost is None or cost < best_cost or (cost == best_cost and label < best_label):
            best_cost = cost
            best_label = label
            best_order = order

    assert best_order is not None and best_label is not None  # for mypy
    return best_order, f"auto-{best_label}"

def _approximate_distance(
    seq_a: str, seq_b: str, sample_indices: list[int]
) -> int:
    mismatches = 0
    for idx in sample_indices:
        if seq_a[idx] != seq_b[idx]:
            mismatches += 1
    return mismatches


def _tree_guided_order(frame: AlignmentFrame) -> list[int] | None:
    """Return a leaf order derived from an associated phylogenetic tree if present."""

    metadata = frame.metadata or {}
    tree_newick = metadata.get("tree_newick")
    if not isinstance(tree_newick, str) or not tree_newick.strip():
        return None

    try:
        root = _parse_newick(tree_newick)
    except Exception:  # pragma: no cover - fall back on parse failures
        return None

    id_to_index = {seq_id: idx for idx, seq_id in enumerate(frame.ids)}
    order: list[int] = []
    seen: set[int] = set()

    def traverse(node) -> bool:
        children = getattr(node, "children", None)
        if children:
            for child in children:
                if not traverse(child):
                    return False
            return True

        label = getattr(node, "label", None)
        if not isinstance(label, str):
            return False
        try:
            index = id_to_index[label]
        except KeyError:
            return False
        if index in seen:
            return False
        seen.add(index)
        order.append(index)
        return True

    if not traverse(root):
        return None

    if len(order) != frame.num_sequences:
        return None

    if sorted(order) != list(range(frame.num_sequences)):
        return None

    return order


def _parse_newick(newick: str):
    text = newick.strip()
    if not text.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    idx = 0

    class _Node:
        __slots__ = ("label", "children")

        def __init__(self, label=None, children=None):
            self.label = label
            self.children = children or []

    def parse_subtree():
        nonlocal idx
        if idx >= len(text):
            raise ValueError("Unexpected end of Newick string")
        if text[idx] == "(":
            idx += 1
            children: list[_Node] = []
            while True:
                children.append(parse_subtree())
                if text[idx] == ",":
                    idx += 1
                    continue
                if text[idx] == ")":
                    idx += 1
                    break
                raise ValueError("Malformed Newick: expected ',' or ')' ")
            label, _ = parse_label_length()
            return _Node(label=label, children=children)
        label, _ = parse_label_length()
        return _Node(label=label, children=[])

    def parse_label_length():
        nonlocal idx
        label_chars: list[str] = []
        while idx < len(text) and text[idx] not in ":,);":
            label_chars.append(text[idx])
            idx += 1
        label = "".join(label_chars) or None
        if idx < len(text) and text[idx] == ":":
            idx += 1
            while idx < len(text) and text[idx] not in ",);":
                idx += 1
        return label, None

    root = parse_subtree()
    if idx >= len(text) or text[idx] != ";":
        raise ValueError("Unexpected content after Newick tree")
    return root


def _compute_similarity_order(frame: AlignmentFrame) -> tuple[AlignmentFrame, list[int], str]:
    num_sequences = frame.num_sequences
    baseline_order = list(range(num_sequences))
    if num_sequences <= 2:
        return frame, baseline_order, "baseline"

    stats = _compute_alignment_stats(frame)

    tree_order = None
    if stats is None or _should_use_tree(stats):
        tree_order = _tree_guided_order(frame)

    if tree_order is not None:
        if tree_order == baseline_order:
            return frame, tree_order, "tree"
        reordered = alignment_from_sequences(
            ids=[frame.ids[idx] for idx in tree_order],
            sequences=[frame.sequences[idx] for idx in tree_order],
            alphabet=frame.alphabet,
            metadata=dict(frame.metadata),
        )
        return reordered, tree_order, "tree"

    length = frame.alignment_length
    if length == 0:
        return frame, baseline_order, "baseline"

    sample_indices = _select_sample_indices(length)
    dist_matrix = _build_distance_matrix(frame.sequences, sample_indices)

    candidates = [
        ("baseline", baseline_order),
        ("mst", _mst_sequence_order(dist_matrix)),
        ("greedy", _greedy_sequence_order(dist_matrix)),
    ]

    best_order, label = _choose_order(dist_matrix, candidates)
    if best_order == baseline_order:
        return frame, best_order, label

    reordered = alignment_from_sequences(
        ids=[frame.ids[idx] for idx in best_order],
        sequences=[frame.sequences[idx] for idx in best_order],
        alphabet=frame.alphabet,
        metadata=dict(frame.metadata),
    )
    return reordered, best_order, label


@dataclass(slots=True)
class CompressedAlignment:
    """Payload plus metadata produced by the compression pipeline."""

    payload: bytes
    metadata: dict[str, Any]


def compress_alignment(frame: AlignmentFrame) -> CompressedAlignment:
    """Compress an alignment into a binary payload and structured metadata."""

    original_frame = frame
    frame, permutation, order_label = _compute_similarity_order(frame)
    permutation_changed = permutation != list(range(frame.num_sequences))

    checksum_value = alignment_checksum(original_frame.sequences)

    column_profiles = collect_column_profiles(frame)
    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, math.ceil(math.log2(max(len(alphabet), 1))))

    run_length_blocks = collect_run_length_blocks(
        column_profiles, frame.num_sequences, symbol_lookup, bits_per_symbol
    )
    bitmask_bytes = (frame.num_sequences + 7) // 8
    run_length_payload = encode_blocks(
        run_length_blocks, bitmask_bytes, bits_per_symbol, frame.alphabet
    )
    seq_id_block = _encode_sequence_ids(frame.ids)
    raw_payload = seq_id_block + run_length_payload
    perm_chunk: bytes | None = None
    perm_meta: dict[str, Any] | None = None
    if permutation_changed:
        perm_chunk, perm_meta = _build_permutation_chunk(permutation)
        if perm_chunk:
            raw_payload = perm_chunk + raw_payload

    payload_candidates: list[tuple[str, bytes, float]] = [("raw", raw_payload, 0.0)]
    # Try zstandard first; its ratio is typically better than zlib for genomic data.
    if _ZSTD_COMPRESSOR is not None:
        zstd_start = time.perf_counter()
        zstd_payload = _ZSTD_COMPRESSOR.compress(raw_payload)
        zstd_seconds = time.perf_counter() - zstd_start
        payload_candidates.append(("zstd", zstd_payload, zstd_seconds))

    zlib_start = time.perf_counter()
    zlib_payload = zlib.compress(raw_payload, level=9)
    zlib_seconds = time.perf_counter() - zlib_start
    payload_candidates.append(("zlib", zlib_payload, zlib_seconds))

    lzma_start = time.perf_counter()
    lzma_payload = lzma.compress(raw_payload, preset=6)
    lzma_seconds = time.perf_counter() - lzma_start
    payload_candidates.append(("xz", lzma_payload, lzma_seconds))

    payload_encoding, payload_bytes, _ = min(payload_candidates, key=lambda item: len(item[1]))
    max_run_length = max((block.run_length for block in run_length_blocks), default=0)
    deviation_columns = sum(1 for block in column_profiles if block.deviations)
    metadata = {
        "format_version": FORMAT_VERSION,
        "codec": "ecomp",
        "num_sequences": frame.num_sequences,
        "alignment_length": frame.alignment_length,
        "alphabet": frame.alphabet,
        "source_format": frame.metadata.get("source_format", "unknown"),
        "checksum_sha256": checksum_value,
        "run_length_blocks": len(run_length_blocks),
        "max_run_length": max_run_length,
        "columns_with_deviations": deviation_columns,
        "bitmask_bytes": bitmask_bytes,
        "bits_per_symbol": bits_per_symbol,
        "payload_encoding": payload_encoding,
        "payload_encoded_bytes": len(payload_bytes),
        "payload_raw_bytes": len(raw_payload),
        "sequence_id_codec": "inline",
        "ordering_strategy": order_label,
    }
    if perm_meta:
        metadata["sequence_permutation"] = perm_meta
    metadata.pop("sequence_ids", None)
    payload_bytes, metadata = _maybe_use_gzip_fallback(
        original_frame, payload_bytes, metadata
    )
    return CompressedAlignment(payload=payload_bytes, metadata=metadata)


def decompress_alignment(payload: bytes, metadata: dict[str, Any]) -> AlignmentFrame:
    """Reconstruct an :class:`AlignmentFrame` from payload and metadata."""

    fallback_info = metadata.get("fallback")
    if fallback_info:
        return _decompress_fallback(payload, metadata)

    expected_columns = metadata["alignment_length"]
    num_sequences = metadata["num_sequences"]
    sequence_ids = metadata.get("sequence_ids")
    if sequence_ids is not None and len(sequence_ids) != num_sequences:
        raise ValueError("Metadata sequence count does not match sequence IDs provided")
    alphabet = metadata.get("alphabet", [])
    bitmask_bytes = metadata.get("bitmask_bytes")
    bits_per_symbol = metadata.get("bits_per_symbol")

    if not isinstance(bitmask_bytes, int) or bitmask_bytes <= 0:
        raise ValueError("Metadata missing valid 'bitmask_bytes' entry")
    if not isinstance(bits_per_symbol, int) or bits_per_symbol <= 0:
        if alphabet:
            bits_per_symbol = max(1, math.ceil(math.log2(len(alphabet))))
        else:
            raise ValueError("Metadata missing valid 'bits_per_symbol' entry")

    payload_encoding = metadata.get("payload_encoding", "raw")
    if payload_encoding == "zstd":
        if _ZSTD_DECOMPRESSOR is None:
            raise RuntimeError(
                "Cannot decompress zstd payload: zstandard module not available"
            )
        decoded_payload = _ZSTD_DECOMPRESSOR.decompress(payload)
    elif payload_encoding == "xz":
        decoded_payload = lzma.decompress(payload)
    elif payload_encoding == "zlib":
        decoded_payload = zlib.decompress(payload)
    elif payload_encoding in {"raw", None}:
        decoded_payload = payload
    else:
        raise ValueError(f"Unsupported payload encoding: {payload_encoding}")

    payload_data = decoded_payload
    permutation: list[int] | None = None
    perm_meta = metadata.get("sequence_permutation")
    if isinstance(perm_meta, dict) and perm_meta.get("encoding") == "payload":
        payload_data, permutation = _extract_permutation_chunk(payload_data, perm_meta)
        metadata["sequence_permutation"] = permutation
    elif perm_meta is not None:
        permutation = _decode_permutation(perm_meta)
        metadata["sequence_permutation"] = permutation

    if payload_data.startswith(_SEQ_ID_MAGIC):
        seq_ids, payload_data = _decode_sequence_ids(payload_data)
        if sequence_ids is None:
            sequence_ids = seq_ids
            metadata["sequence_ids"] = seq_ids
        elif list(sequence_ids) != seq_ids:
            raise ValueError("Sequence IDs mismatch between metadata and payload")
    if sequence_ids is None:
        raise ValueError("Metadata missing sequence identifiers")
    if len(sequence_ids) != num_sequences:
        raise ValueError(
            "Metadata sequence count does not match sequence IDs provided"
        )

    blocks = decode_blocks(payload_data, bitmask_bytes, bits_per_symbol, alphabet)
    sequences = [["" for _ in range(expected_columns)] for _ in range(num_sequences)]

    try:
        symbol_table = list(alphabet)
    except TypeError as exc:  # pragma: no cover - guard against malformed metadata
        raise ValueError("Alphabet metadata is not iterable") from exc

    column_index = 0
    for block in blocks:
        consensus = block.consensus
        for _ in range(block.run_length):
            if column_index >= expected_columns:
                raise ValueError(
                    "Decoded columns exceed expected alignment length"
                )
            residue_indices = _iter_deviation_indices(block.bitmask, num_sequences)
            residues = _decode_residues(
                block.residues,
                len(residue_indices),
                bits_per_symbol,
                symbol_table,
            )
            for seq_list in sequences:
                seq_list[column_index] = consensus
            for seq_index, residue in zip_strict(residue_indices, residues):
                sequences[seq_index][column_index] = residue
            column_index += 1

    if column_index != expected_columns:
        raise ValueError(
            f"Decoded columns ({column_index}) do not match expected length {expected_columns}"
        )

    reconstructed = ["".join(row) for row in sequences]

    permutation = metadata.get("sequence_permutation")
    if permutation:
        inverse = [0] * len(permutation)
        for new_pos, original_index in enumerate(permutation):
            inverse[original_index] = new_pos
        reconstructed = [reconstructed[inverse[idx]] for idx in range(len(inverse))]
        sequence_ids = [sequence_ids[inverse[idx]] for idx in range(len(inverse))]

    return alignment_from_sequences(
        ids=sequence_ids,
        sequences=reconstructed,
        alphabet=alphabet,
        metadata={"source_format": metadata.get("source_format", "unknown")},
    )


def _iter_deviation_indices(bitmask: bytes, num_sequences: int) -> List[int]:
    indices: list[int] = []
    for seq_index in range(num_sequences):
        byte_index = seq_index // 8
        bit_index = seq_index % 8
        if bitmask[byte_index] & (1 << bit_index):
            indices.append(seq_index)
    return indices


def _decode_residues(
    data: bytes,
    count: int,
    bits_per_symbol: int,
    alphabet: Iterable[str],
) -> list[str]:
    if count == 0:
        return []
    mask = (1 << bits_per_symbol) - 1
    values: list[int] = []
    buffer = 0
    bits_in_buffer = 0
    data_iter = iter(data)
    while len(values) < count:
        while bits_in_buffer < bits_per_symbol:
            try:
                byte = next(data_iter)
            except StopIteration as exc:  # pragma: no cover - corruption guard
                raise ValueError("Insufficient residue data during decode") from exc
            buffer = (buffer << 8) | byte
            bits_in_buffer += 8
        shift = bits_in_buffer - bits_per_symbol
        value = (buffer >> shift) & mask
        buffer &= (1 << shift) - 1
        bits_in_buffer -= bits_per_symbol
        values.append(value)
    alphabet_list = list(alphabet)
    try:
        return [alphabet_list[value] for value in values]
    except IndexError as exc:  # pragma: no cover - corruption guard
        raise ValueError("Residue code exceeds alphabet size") from exc


def _alignment_to_fasta_bytes(frame: AlignmentFrame) -> bytes:
    buffer = io.StringIO()
    for seq_id, sequence in zip_strict(frame.ids, frame.sequences):
        buffer.write(f">{seq_id}\n{sequence}\n")
    return buffer.getvalue().encode("utf-8")


def _parse_fasta_bytes(data: bytes) -> AlignmentFrame:
    ids: list[str] = []
    sequences: list[str] = []
    current_id: str | None = None
    current_seq: list[str] = []
    for raw_line in data.decode("utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequences.append("".join(current_seq))
            current_id = line[1:]
            ids.append(current_id)
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        sequences.append("".join(current_seq))
    return alignment_from_sequences(ids=ids, sequences=sequences)


def _maybe_use_gzip_fallback(
    original_frame: AlignmentFrame, payload: bytes, metadata: dict[str, Any]
) -> tuple[bytes, dict[str, Any]]:
    fasta_bytes = _alignment_to_fasta_bytes(original_frame)
    gzip_payload = gzip.compress(fasta_bytes)
    if len(gzip_payload) + 1 < len(payload) and len(gzip_payload) < len(fasta_bytes):
        updated = dict(metadata)
        updated.pop("sequence_permutation", None)
        updated["fallback"] = {
            "type": "gzip",
            "format": metadata.get("source_format", "fasta"),
        }
        updated["payload_encoding"] = "gzip"
        updated["payload_encoded_bytes"] = len(gzip_payload)
        updated["payload_raw_bytes"] = len(fasta_bytes)
        return gzip_payload, updated
    return payload, metadata


def _decompress_fallback(payload: bytes, metadata: dict[str, Any]) -> AlignmentFrame:
    info = metadata.get("fallback", {})
    fallback_type = info.get("type")
    if fallback_type == "gzip":
        data = gzip.decompress(payload)
        frame = _parse_fasta_bytes(data)
        frame.metadata["source_format"] = info.get("format", "fasta")
        return frame
    raise ValueError(f"Unsupported fallback type: {fallback_type}")
