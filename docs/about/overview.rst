Overview
========

eComp is a lossless compressor for multiple sequence alignments (MSAs) that leans
on evolutionary structure:

- A **consensus-first pipeline** stores only the deviations from each column's
  majority residue.
- Optional **tree-guided ordering** places related taxa next to each other so the
  deviation sets shrink even further (the tree is used transiently and is not
  stored in the archive).
- Every archive carries rich metadata (alphabet, checksums, permutation, optional
  tree) so downstream tooling can validate and analyse without guesswork.

The rest of the documentation is organised as follows:

- :doc:`../usage/quickstart` – installation, CLI basics, and Python API entry points.
- :doc:`../usage/benchmarking` – how to reproduce the comparisons against gzip,
  bzip2, xz, and tree-guided experiments.
- :doc:`../tutorials/workflows` – end-to-end examples (automating tree ordering,
  integrating with existing phylogenomic workflows).
- :doc:`../frequently_asked_questions/index` – troubleshooting and best practices.
- :doc:`../change_log/index` – release highlights and compatibility notes.
- :doc:`../other_software/index` – companion tools used alongside eComp.
- :doc:`../data/index` – datasets bundled with the repository and how to obtain
  larger benchmark panels.

New users should start with the quickstart, while contributors can dive straight
into the workflow and API sections.
