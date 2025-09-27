Command Line Interface
======================

The ``ecomp`` CLI follows the ClipKIT convention of subcommands with concise
flags. Each invocation accepts an alignment in FASTA/PHYLIP format and optionally
an accompanying Newick tree.

Quick Reference
---------------

.. code-block:: bash

   # Compress alignment only (produces example.ecomp + metadata JSON)
   ecomp compress example.fasta

   # Compress alignment + tree bundle (produces example.ecbt + metadata JSON)
   ecomp compress example.fasta example.tree

   # Decompress archives (alignment codec auto-detected via metadata)
   ecomp decompress example.ecbt --alignment-output restored.fasta --tree-output restored.tree

   # Inspect metadata in either JSON or summary form
   ecomp inspect example.ecomp --summary

Subcommand Summary
------------------

``compress``
   - Auto-selects the codec (alignment vs. phylo bundle) unless ``--codec`` is
     provided.
   - Supports explicit ``--bundle-suffix`` to override the ``.ecbt`` default.
   - ``--stats`` prints original vs. compressed byte counts for quick feedback.

``decompress``
   - Validates checksums by default; disable with ``--no-checksum`` when
     benchmarking.
   - ``--alignment-output`` and ``--tree-output`` control output paths.

``inspect``
   - ``--summary`` emits human readable key fields (codec, sequence count,
     payload encoding).
   - Without flags the command prints the full JSON metadata.

Codec Selection
---------------

- ``alignment`` – standard evolutionary codec producing ``.ecomp`` payloads.
- ``phylo`` – bundle codec storing alignment + tree; outputs ``.ecbt`` by
  default.
- ``auto`` (default) – alignment-only unless a tree path is supplied.

Exit Codes
----------

- ``0`` – success.
- ``>0`` – fatal errors (missing files, checksum mismatch, invalid metadata).
  Diagnostics are printed to stderr with actionable guidance.
