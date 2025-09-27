"""Co-compress an alignment and its associated phylogenetic tree."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..diagnostics.checksums import alignment_checksum
from ..io import AlignmentFrame, alignment_from_sequences


BRANCH_SCALE = 10000  # store branch lengths as integer millimeters


# ---------------------------------------------------------------------------
# Newick parsing
# ---------------------------------------------------------------------------


@dataclass
class _ParserNode:
    label: Optional[str]
    length: float
    children: List["_ParserNode"]


def _parse_newick(newick: str) -> _ParserNode:
    text = newick.strip()
    if not text.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    idx = 0

    def parse_subtree() -> _ParserNode:
        nonlocal idx
        if text[idx] == "(":
            idx += 1
            children: List[_ParserNode] = []
            while True:
                children.append(parse_subtree())
                if text[idx] == ",":
                    idx += 1
                    continue
                if text[idx] == ")":
                    idx += 1
                    break
                raise ValueError("Malformed Newick: expected ',' or ')' ")
            label, length = parse_label_length()
            return _ParserNode(label=label, length=length, children=children)
        else:
            label, length = parse_label_length()
            return _ParserNode(label=label, length=length, children=[])

    def parse_label_length() -> Tuple[Optional[str], float]:
        nonlocal idx
        label_chars: List[str] = []
        while idx < len(text) and text[idx] not in ":,);":
            label_chars.append(text[idx])
            idx += 1
        label = "".join(label_chars) or None
        length = 0.0
        if idx < len(text) and text[idx] == ":":
            idx += 1
            length_chars: List[str] = []
            while idx < len(text) and text[idx] not in ",);":
                length_chars.append(text[idx])
                idx += 1
            try:
                length = float("".join(length_chars))
            except ValueError as exc:  # pragma: no cover - guard
                raise ValueError("Invalid branch length in Newick") from exc
        return label, length

    root = parse_subtree()
    if text[idx] != ";":
        raise ValueError("Unexpected content after Newick tree")
    return root


# ---------------------------------------------------------------------------
# Tree structure helpers
# ---------------------------------------------------------------------------


@dataclass
class _TreeNode:
    index: int
    label: Optional[str]
    parent: Optional[int]
    children: List[int]
    length: float


def _collect_nodes(root: _ParserNode) -> Tuple[List[_TreeNode], List[int]]:
    nodes: List[_TreeNode] = []
    leaf_indices: List[int] = []

    def visit(parser_node: _ParserNode, parent: Optional[int]) -> int:
        index = len(nodes)
        node = _TreeNode(
            index=index,
            label=parser_node.label,
            parent=parent,
            children=[],
            length=parser_node.length,
        )
        nodes.append(node)
        for child in parser_node.children:
            child_index = visit(child, index)
            node.children.append(child_index)
        if not parser_node.children:
            leaf_indices.append(index)
        return index

    visit(root, None)
    return nodes, leaf_indices


# ---------------------------------------------------------------------------
# Fitch parsimony to assign internal sequences
# ---------------------------------------------------------------------------


def _assign_sequences(
    nodes: List[_TreeNode],
    leaf_indices: Iterable[int],
    leaf_sequences: Dict[int, str],
) -> List[str]:
    length = len(next(iter(leaf_sequences.values())))
    seq_buffers = [list("" for _ in range(length)) for _ in nodes]

    for pos in range(length):
        sets: List[set] = [set() for _ in nodes]

        def postorder(idx: int) -> set:
            node = nodes[idx]
            if not node.children:
                sets[idx] = {leaf_sequences[idx][pos]}
                return sets[idx]
            child_sets = [postorder(child) for child in node.children]
            intersect = set.intersection(*child_sets)
            sets[idx] = intersect if intersect else set.union(*child_sets)
            return sets[idx]

        postorder(0)

        def preorder(idx: int, parent_char: Optional[str]) -> None:
            node_set = sets[idx]
            char = parent_char if parent_char in node_set else next(iter(node_set))
            seq_buffers[idx][pos] = char
            for child in nodes[idx].children:
                preorder(child, char)

        preorder(0, next(iter(sets[0])))

    return ["".join(seq_buffers[idx]) for idx in range(len(nodes))]


# ---------------------------------------------------------------------------
# Varint helpers
# ---------------------------------------------------------------------------


def _encode_varint(value: int) -> bytes:
    if value < 0:
        raise ValueError("varint cannot encode negative values")
    out = bytearray()
    while True:
        to_write = value & 0x7F
        value >>= 7
        if value:
            out.append(to_write | 0x80)
        else:
            out.append(to_write)
            break
    return bytes(out)


def _decode_varint(data: memoryview, cursor: int) -> Tuple[int, int]:
    shift = 0
    result = 0
    while True:
        byte = data[cursor]
        cursor += 1
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return result, cursor
        shift += 7
        if shift > 56:  # pragma: no cover - guard
            raise ValueError("Varint too long")


# ---------------------------------------------------------------------------
# Main compression / decompression
# ---------------------------------------------------------------------------


def compress_alignment_with_tree(
    frame: AlignmentFrame, newick: str
) -> Tuple[bytes, dict[str, object]]:
    parser_root = _parse_newick(newick)
    nodes, leaf_indices = _collect_nodes(parser_root)

    leaf_label_to_index = {nodes[idx].label: idx for idx in leaf_indices}
    if None in leaf_label_to_index:
        raise ValueError("All leaves must be labelled in the Newick tree")

    if set(frame.ids) != set(leaf_label_to_index.keys()):
        raise ValueError("Alignment sequence IDs must match Newick leaf labels")

    leaf_sequences = {
        leaf_label_to_index[id_]: seq for id_, seq in zip(frame.ids, frame.sequences)
    }

    node_sequences = _assign_sequences(nodes, leaf_indices, leaf_sequences)

    payload = bytearray()
    payload.extend(b"PB01")
    payload.extend(_encode_varint(len(nodes)))
    for node in nodes:
        parent_encoded = 0 if node.parent is None else node.parent + 1
        payload.extend(_encode_varint(parent_encoded))
        length_q = int(round(node.length * BRANCH_SCALE))
        payload.extend(_encode_varint(length_q))

    payload.extend(_encode_varint(len(leaf_indices)))
    for idx in leaf_indices:
        label = nodes[idx].label or ""
        payload.extend(_encode_varint(idx))
        payload.extend(_encode_varint(len(label)))
        payload.extend(label.encode("utf-8"))

    root_sequence = node_sequences[0]
    payload.extend(_encode_varint(len(root_sequence)))
    payload.extend(root_sequence.encode("ascii"))

    for node_index in range(1, len(nodes)):
        parent_index = nodes[node_index].parent
        parent_seq = node_sequences[parent_index]
        node_seq = node_sequences[node_index]
        diffs: List[Tuple[int, str]] = []
        for pos, (p_char, n_char) in enumerate(zip(parent_seq, node_seq)):
            if p_char != n_char:
                diffs.append((pos, n_char))
        payload.extend(_encode_varint(len(diffs)))
        prev = 0
        for pos, char in diffs:
            delta = pos - prev
            payload.extend(_encode_varint(delta))
            payload.append(ord(char))
            prev = pos

    metadata = {
        "codec": "phylo-bundle",
        "num_sequences": frame.num_sequences,
        "alignment_length": frame.alignment_length,
        "sequence_ids": frame.ids,
        "alphabet": frame.alphabet,
        "source_format": frame.metadata.get("source_format", "unknown"),
        "checksum_sha256": alignment_checksum(frame.sequences),
        "branch_scale": BRANCH_SCALE,
    }
    return bytes(payload), metadata


def decompress_alignment_with_tree(
    payload: bytes, metadata: dict[str, object]
) -> Tuple[AlignmentFrame, str]:
    data = memoryview(payload)
    cursor = 0
    if data[:4].tobytes() != b"PB01":
        raise ValueError("Invalid phylo bundle header")
    cursor += 4

    node_count, cursor = _decode_varint(data, cursor)
    parents: List[Optional[int]] = []
    lengths: List[float] = []
    for _ in range(node_count):
        parent_encoded, cursor = _decode_varint(data, cursor)
        parent = None if parent_encoded == 0 else parent_encoded - 1
        parent_val = min(parent, node_count - 1) if parent is not None else None
        parents.append(parent_val)
        length_q, cursor = _decode_varint(data, cursor)
        lengths.append(length_q / BRANCH_SCALE)

    children: List[List[int]] = [[] for _ in range(node_count)]
    for idx, parent in enumerate(parents):
        if parent is not None:
            children[parent].append(idx)

    leaf_count, cursor = _decode_varint(data, cursor)
    leaf_labels: Dict[int, str] = {}
    for _ in range(leaf_count):
        leaf_idx, cursor = _decode_varint(data, cursor)
        label_len, cursor = _decode_varint(data, cursor)
        label_bytes = data[cursor:cursor + label_len].tobytes()
        cursor += label_len
        leaf_labels[leaf_idx] = label_bytes.decode("utf-8")

    alignment_length, cursor = _decode_varint(data, cursor)
    root_sequence = data[cursor:cursor + alignment_length].tobytes().decode("ascii")
    cursor += alignment_length

    node_sequences: List[List[str]] = [[ch for ch in root_sequence]] + [None] * (node_count - 1)

    for node_index in range(1, node_count):
        diff_count, cursor = _decode_varint(data, cursor)
        parent_index = parents[node_index]
        parent_seq = node_sequences[parent_index].copy()
        pos = 0
        for _ in range(diff_count):
            delta, cursor = _decode_varint(data, cursor)
            pos += delta
            residue = chr(data[cursor])
            cursor += 1
            parent_seq[pos] = residue
        node_sequences[node_index] = parent_seq

    leaf_sequences = {
        leaf_labels[idx]: "".join(node_sequences[idx])
        for idx in leaf_labels
    }

    ordered_sequences = [leaf_sequences[id_] for id_ in metadata["sequence_ids"]]
    frame = alignment_from_sequences(
        ids=metadata["sequence_ids"],
        sequences=ordered_sequences,
        alphabet=metadata.get("alphabet", []),
        metadata={"source_format": metadata.get("source_format", "unknown")},
    )

    newick = _build_newick(children, parents, leaf_labels, lengths)
    return frame, newick


# ---------------------------------------------------------------------------
# Newick reconstruction
# ---------------------------------------------------------------------------


def _build_newick(
    children: Sequence[Sequence[int]],
    parents: Sequence[Optional[int]],
    leaf_labels: Dict[int, str],
    lengths: Sequence[float],
) -> str:
    def recurse(node: int) -> str:
        child_list = children[node]
        if child_list:
            inner = ",".join(recurse(ch) for ch in child_list)
            label = '' if node not in leaf_labels else leaf_labels[node]
            part = f"({inner}){label if label else ''}:{lengths[node]:.6f}"
        else:
            label = leaf_labels.get(node, '')
            part = f"{label}:{lengths[node]:.6f}"
        return part

    return recurse(0) + ";"

