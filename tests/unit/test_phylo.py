from __future__ import annotations

from io import StringIO

from Bio import Phylo

from ecomp.phylo import infer_distance_tree_from_frame, tree_to_newick
from ecomp.io import alignment_from_sequences


def test_distance_tree_groups_similar_sequences() -> None:
    frame = alignment_from_sequences(
        ["s1", "s2", "s3"],
        [
            "AAAA",
            "AAAT",
            "TTTT",
        ],
    )

    tree = infer_distance_tree_from_frame(frame)
    newick = tree_to_newick(tree)
    assert "s1" in newick and "s2" in newick and "s3" in newick

    parsed = Phylo.read(StringIO(newick), "newick")
    assert {term.name for term in parsed.get_terminals()} == {"s1", "s2", "s3"}

    # Ensure alternate method also works
    upgma_tree = infer_distance_tree_from_frame(frame, method="upgma")
    upgma_newick = tree_to_newick(upgma_tree)
    assert upgma_newick.endswith(";")
