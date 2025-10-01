.. role:: raw-html(raw)
    :format: html

About
=====

eComp is a lossless compressor for multiple sequence alignments (MSAs) that leans
on evolutionary structure:

- A **consensus-first pipeline** stores only the deviations from each column's
  majority residue.
- Optional **tree-guided ordering** places related taxa next to each other so the
  deviation sets shrink even further (the tree is used transiently and is not
  stored in the archive).
- Every archive carries rich metadata (alphabet, checksums, permutation, optional
  tree) so downstream tooling can validate and analyse without guesswork.

:raw-html:`<br />`

Suggested reading order:

- :doc:`../usage/index` – quickstart guide, benchmarking tips, API examples, and
  advanced metrics.
- :doc:`../change_log/index` – release highlights and compatibility notes.
- :doc:`../other_software/index` – tools that pair well with eComp in
  phylogenomic pipelines.
- :doc:`../faq/index` – troubleshooting and best practices.

New users should begin with the quickstart section under :doc:`../usage/index`;
returning users can jump straight to the FAQ or change log as needed.
