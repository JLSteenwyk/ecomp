#!/usr/bin/env python3
"""Benchmark eComp compression with different sequence orderings."""

from __future__ import annotations

import argparse
import json
import os
import random
import tempfile
import time
from pathlib import Path
from typing import Iterable, List

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.BaseTree import Clade

from evolutionary_compression import read_alignment
from evolutionary_compression.compression.pipeline import (
    compress_alignment,
    decompress_alignment,
)
from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.storage import write_metadata, write_payload


def _metadata_size(metadata: dict) -> int:
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        temp_path = Path(tmp.name)
    try:
        write_metadata(temp_path, metadata)
        return temp_path.stat().st_size
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _compress_case(frame, original_size: int) -> dict[str, float | int]:
    start = time.perf_counter()
    compressed = compress_alignment(frame)
    compress_seconds = time.perf_counter() - start

    start = time.perf_counter()
    restored = decompress_alignment(compressed.payload, compressed.metadata)
    decompress_seconds = time.perf_counter() - start

    if restored.sequences != frame.sequences:
        raise RuntimeError("Round-trip mismatch detected during benchmarking")

    fd, payload_name = tempfile.mkstemp(suffix=".ecomp")
    os.close(fd)
    payload_path = Path(payload_name)
    try:
        write_payload(payload_path, compressed.payload)
        compressed_size = payload_path.stat().st_size
    finally:
        payload_path.unlink(missing_ok=True)

    metadata_bytes = _metadata_size(compressed.metadata)
    ratio = original_size / compressed_size if compressed_size else float("inf")

    return {
        "compressed_size": compressed_size,
        "metadata_bytes": metadata_bytes,
        "compression_ratio": ratio,
        "compress_seconds": compress_seconds,
        "decompress_seconds": decompress_seconds,
        "ordering_strategy": compressed.metadata.get("ordering_strategy"),
    }


def _compress_with_strategy(frame, original_size: int, strategy: str | None) -> dict[str, float | int | None]:
    previous = os.environ.get("ECOMP_SEQUENCE_ORDER")
    try:
        if strategy in {None, "auto"}:
            os.environ.pop("ECOMP_SEQUENCE_ORDER", None)
        else:
            os.environ["ECOMP_SEQUENCE_ORDER"] = strategy
        result = _compress_case(frame, original_size)
    finally:
        if previous is None:
            os.environ.pop("ECOMP_SEQUENCE_ORDER", None)
        else:
            os.environ["ECOMP_SEQUENCE_ORDER"] = previous
    result.setdefault("ordering_strategy", strategy)
    return result


def _sample_columns(sequences: List[str], max_columns: int = 200) -> tuple[List[str], List[int]]:
    if not sequences:
        return sequences, []
    length = len(sequences[0])
    if length <= max_columns:
        indices = list(range(length))
    else:
        rng = random.Random(length)
        indices = sorted(rng.sample(range(length), max_columns))
    sampled = ["".join(seq[idx] for idx in indices) for seq in sequences]
    return sampled, indices


def _find_parent(root: Clade, target_name: str) -> tuple[Clade | None, Clade | None]:
    for idx, child in enumerate(root.clades):
        if child.name == target_name:
            return root, child
        parent, node = _find_parent(child, target_name)
        if node is not None:
            return parent, node
    return None, None


def _build_nj_tree(frame, max_taxa: int = 200) -> str:
    num_sequences = frame.num_sequences
    sampled_sequences, _ = _sample_columns(frame.sequences)

    if num_sequences > max_taxa:
        rng = random.Random(num_sequences)
        subset_indices = sorted(rng.sample(range(num_sequences), max_taxa))
    else:
        subset_indices = list(range(num_sequences))

    subset_records = [
        SeqRecord(Seq(sampled_sequences[idx]), id=frame.ids[idx], description="")
        for idx in subset_indices
    ]
    alignment = MultipleSeqAlignment(subset_records)
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, "nj")
    tree = constructor.build_tree(alignment)

    subset_set = set(subset_indices)
    leftovers = [idx for idx in range(num_sequences) if idx not in subset_set]

    # Precompute sampled sequences for leftover nearest-neighbour mapping
    subset_samples = {idx: sampled_sequences[idx] for idx in subset_indices}
    leftover_assignments: dict[int, list[int]] = {idx: [] for idx in subset_indices}

    for lf_idx in leftovers:
        seq = sampled_sequences[lf_idx]
        best_idx = None
        best_score = float("inf")
        for anchor_idx, anchor_seq in subset_samples.items():
            mismatches = sum(ch1 != ch2 for ch1, ch2 in zip(seq, anchor_seq))
            if mismatches < best_score:
                best_score = mismatches
                best_idx = anchor_idx
        if best_idx is not None:
            leftover_assignments[best_idx].append(lf_idx)

    for anchor_idx, assigned in leftover_assignments.items():
        if not assigned:
            continue
        anchor_id = frame.ids[anchor_idx]
        parent, anchor_clade = _find_parent(tree.root, anchor_id)
        if parent is None or anchor_clade is None:
            continue
        new_children = [Clade(branch_length=0.1, name=anchor_id)]
        for lf_idx in sorted(assigned):
            new_children.append(Clade(branch_length=0.1, name=frame.ids[lf_idx]))
        replacement = Clade(branch_length=anchor_clade.branch_length or 0.1, clades=new_children)
        target_idx = parent.clades.index(anchor_clade)
        parent.clades[target_idx] = replacement

    return tree.format("newick")


def run(paths: Iterable[Path]) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for path in paths:
        alignment_path = path.expanduser().resolve()
        if not alignment_path.exists():
            raise FileNotFoundError(alignment_path)

        base_frame = read_alignment(alignment_path)
        original_size = alignment_path.stat().st_size
        base_ids = list(base_frame.ids)
        base_sequences = list(base_frame.sequences)
        base_alphabet = base_frame.alphabet
        base_metadata = dict(base_frame.metadata)

        def make_frame(extra_metadata: dict | None = None):
            metadata = dict(base_metadata)
            if extra_metadata:
                metadata.update(extra_metadata)
            return alignment_from_sequences(
                ids=base_ids,
                sequences=base_sequences,
                alphabet=base_alphabet,
                metadata=metadata,
            )

        for strategy in ("baseline", "mst", "greedy", "auto"):
            frame_variant = make_frame()
            result = _compress_with_strategy(frame_variant, original_size, strategy if strategy != "auto" else None)
            result.update(
                {
                    "alignment": str(alignment_path),
                    "mode": f"strategy_{strategy}",
                    "tree_source": None,
                }
            )
            results.append(result)

        provided_tree_path = alignment_path.with_suffix(alignment_path.suffix + ".tre")
        if provided_tree_path.exists():
            tree_text = provided_tree_path.read_text()
            frame_with_tree = make_frame({"tree_newick": tree_text})
            provided = _compress_with_strategy(frame_with_tree, original_size, None)
            provided.update(
                {
                    "alignment": str(alignment_path),
                    "mode": "tree_provided_auto",
                    "tree_source": str(provided_tree_path),
                }
            )
            results.append(provided)

        nj_newick = _build_nj_tree(make_frame())
        frame_with_nj = make_frame({"tree_newick": nj_newick})
        nj_result = _compress_with_strategy(frame_with_nj, original_size, None)
        nj_result.update(
            {
                "alignment": str(alignment_path),
                "mode": "tree_nj_auto",
                "tree_source": "nj",
            }
        )
        results.append(nj_result)

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark tree-guided ordering strategies")
    parser.add_argument("inputs", nargs="+", help="Alignment paths to benchmark")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/tree_ordering_nj_comparison.json"),
        help="JSON file to store benchmark results",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paths = [Path(p) for p in args.inputs]
    results = run(paths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(f"Wrote {len(results)} records to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
