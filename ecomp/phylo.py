"""Phylogenetic helpers that operate directly on eComp archives."""

from __future__ import annotations

import math
from io import StringIO
from pathlib import Path
from typing import Iterable, Literal

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from .compression.pipeline import decompress_alignment
from .diagnostics.metrics import pairwise_identity_matrix
from .io import AlignmentFrame
from .storage import read_archive

TreeMethod = Literal["nj", "upgma"]


def _distance_matrix_from_frame(frame: AlignmentFrame) -> DistanceMatrix:
    pid = pairwise_identity_matrix(frame)
    names = frame.ids
    matrix: list[list[float]] = []
    for i, name_i in enumerate(names):
        row: list[float] = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
                continue
            identity = pid.matrix[i][j]
            if math.isnan(identity):
                distance = 1.0
            else:
                distance = max(0.0, min(1.0, 1.0 - identity))
            row.append(distance)
        matrix.append(row)
    return DistanceMatrix(names=names, matrix=matrix)


def infer_distance_tree_from_frame(
    frame: AlignmentFrame,
    *,
    method: TreeMethod = "nj",
) -> Tree:
    """Infer a distance-based tree from an in-memory alignment."""

    dm = _distance_matrix_from_frame(frame)
    constructor = DistanceTreeConstructor()
    if method == "upgma":
        tree = constructor.upgma(dm)
    else:
        tree = constructor.nj(dm)
    return tree


def infer_distance_tree(
    archive_path: str | Path,
    *,
    metadata_path: str | Path | None = None,
    method: TreeMethod = "nj",
) -> Tree:
    """Infer a distance tree directly from an `.ecomp` archive."""

    payload, metadata, _ = read_archive(Path(archive_path), metadata_path=metadata_path)
    frame = decompress_alignment(payload, metadata)
    return infer_distance_tree_from_frame(frame, method=method)


def tree_to_newick(tree: Tree) -> str:
    """Serialise *tree* to a Newick string."""

    buffer = StringIO()
    Phylo.write(tree, buffer, "newick")
    return buffer.getvalue().strip()


__all__ = [
    "infer_distance_tree",
    "infer_distance_tree_from_frame",
    "tree_to_newick",
]
