Related Software
================

PhyKIT
  Provides a wide catalogue of alignment and tree diagnostics.  Many of the
  alignment statistics (entropy, constant sites, gap fractions) can be computed
  directly on eComp payloads by streaming over the run-length blocks.

IQ-TREE / RAxML-NG
  Use these to infer phylogenies after decompressing an alignment, or to supply
  trees that improve compression (pass ``--tree path/to/tree`` to ``ecomp compress``
  or set ``frame.metadata['tree_newick']`` in the Python API).

MAFFT / ClipKIT / AMAS
  eComp plays nicely with alignment filtering and concatenation tools.  Run your
  preferred pre-processing pipeline first, then archive the resulting alignment
  for reproducibility.

Standard compressors (gzip, bzip2, xz)
  Included in the benchmark harness as baselines.  Tree-guided eComp typically
  outperforms gzip and bzip2 on evolutionary datasets; xz remains competitive for
  some large protein families.
