from __future__ import annotations

from io import StringIO

from Bio import Phylo

from ecomp import ezip
from ecomp.io import alignment_from_sequences
from ecomp.phylo import infer_distance_tree, infer_distance_tree_from_frame, tree_to_newick


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


def test_infer_distance_tree_from_archive(tmp_path) -> None:
    alignment = tmp_path / "example.fasta"
    alignment.write_text(">a\nAAAA\n>b\nAAAT\n>c\nTTTT\n")
    metadata_path = tmp_path / "example.json"
    archive, returned_metadata = ezip(alignment, metadata_path=metadata_path)
    tree = infer_distance_tree(archive)
    newick = tree_to_newick(tree)
    assert newick.endswith(";")
    assert {term.name for term in tree.get_terminals()} == {"a", "b", "c"}

    # Ensure metadata path works when provided explicitly
    tree_with_metadata = infer_distance_tree(archive, metadata_path=returned_metadata)
    assert {term.name for term in tree_with_metadata.get_terminals()} == {"a", "b", "c"}
