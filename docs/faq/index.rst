FAQ
===

How do I install eComp offline?
-------------------------------

Install the Python dependencies from a local wheelhouse or mirror, then run
``pip install .`` inside the repository.  eComp only depends on Biopython, NumPy,
and bitarray at runtime.

Why does compression fail with "sequence IDs mismatch"?
-------------------------------------------------------

This indicates the metadata (usually the JSON file) was edited or mismatched with
the payload.  Run ``ecomp inspect`` to confirm the checksum, or re-run
``ecomp zip`` to regenerate both files together.

Can I supply a tree to improve compression?
-------------------------------------------

Yes.  Pass ``--tree path/to/tree.nwk`` when running ``ecomp zip``.  The tree is
used only to reorder the sequences; it is *not* stored in the archive.  This often
reduces the number of deviations per column and improves the compression ratio.

Does eComp ever fall back to gzip?
---------------------------------

If the custom payload is larger than ``gzip`` would produce, the compressor stores
the gzip result instead.  The metadata still records the original checksum so you
can verify the round trip.

How can I compute statistics without fully decompressing?
---------------------------------------------------------

Use the Python API.  ``decompress_alignment`` returns an ``AlignmentFrame``; you can
modify that helper to stream over the run-length blocks and compute per-column
metrics (entropy, gap fraction, etc.) without writing intermediate FASTA files.
