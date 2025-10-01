Python API
==========

The CLI is a thin wrapper around a small Python surface that lives in
``ecomp.__init__``.  Import these helpers when you want to plug
compression into a larger workflow without spawning subprocesses.

.. code-block:: python

    from ecomp import (
        read_alignment,
        compress_alignment,
        decompress_alignment,
        compress_file,
        decompress_file,
    )

Common entry points
-------------------

``compress_file(input_path, output_path=None, metadata_path=None, input_format=None)``
    Reads an alignment from ``input_path`` and writes an ``.ecomp`` payload plus a
    JSON metadata file.  If ``output_path`` or ``metadata_path`` are omitted they are
    derived from the input name (exactly like the CLI).

``decompress_file(ecomp_path, output_path=None, metadata_path=None, output_format=None)``
    Reconstructs an alignment file from a payload + metadata pair.  When
    ``output_path`` is omitted a ``.fasta`` file is written next to the archive.

``compress_alignment(frame)`` / ``decompress_alignment(payload, metadata)``
    Lower-level helpers that operate on :class:`AlignmentFrame` objects.  Use these
    when the alignment is already in memory—e.g., you obtained it from Biopython or
    from a custom generator.

AlignmentFrame
--------------

``read_alignment(path, fmt=None)`` returns an :class:`AlignmentFrame` with four
fields:

- ``ids`` – list of sequence identifiers in read order.
- ``sequences`` – list of equally long strings.
- ``alphabet`` – sorted unique symbol list (useful for diagnostics).
- ``metadata`` – dictionary with ``source_path`` and ``source_format``; you can
  attach extra keys (for example ``tree_newick`` to enable tree-guided ordering).

The same structure can be constructed manually via
``alignment_from_sequences(ids, sequences, alphabet=None, metadata=None)`` when you
want to reorder sequences or inject metadata before a compression run.

Tree-guided ordering hook
-------------------------

If ``frame.metadata`` contains ``tree_newick`` (a Newick string whose leaf labels
match the sequence IDs) the compressor will try to derive a depth-first ordering
guided by that tree before emitting blocks.  For noisy alignments—lots of gaps,
high composition heterogeneity, or very low pairwise identity—the optimizer now
skips the tree and falls back to the similarity-based order.  Set the environment
variable ``ECOMP_SEQUENCE_ORDER=tree`` if you want to force tree usage regardless
of those heuristics.  The tree itself is never persisted; ordering hints do not
inflate archive sizes.

Error handling
--------------

All helpers raise ``ValueError`` for malformed metadata, alignment length
mismatches, or corrupted payloads.  Wrap calls in ``try/except`` if you need
custom recovery logic.  The low-level API never writes to disk until all checks pass,
so a failure leaves the original archive untouched.
