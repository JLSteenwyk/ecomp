Benchmarks
==========

The benchmarking helper mirrors ClipKIT's philosophy of reproducible comparisons
against widely used compressors.

Quick Start
-----------

.. code-block:: bash

   ecomp compress example.fasta --stats
   python scripts/compare_compressors.py example_data/*.fasta \
       --codecs ecomp phylo-bundle gzip bzip2 xz \
       --output benchmark_results.json

Command Reference
-----------------

``scripts/compare_compressors.py`` accepts one or more input alignments and a
set of codec names.

- ``ecomp`` – standard alignment codec.
- ``phylo-bundle`` – alignment + tree bundle (requires ``<alignment>.tree`` or
  ``.nwk`` in the same directory).
- ``gzip``/``bzip2``/``xz`` – external compressors referenced via the system
  ``PATH``.

The script reports compression ratios, timings, and notes (e.g., fallback to the
legacy gzip payload). Results can be written to JSON for dashboarding.

Interpreting Results
--------------------

- Ratios > 1.0 indicate that the codec outperforms the original size.
- When comparing against gzip, sum the compressed alignment and tree sizes for a
  fair comparison with the phylo bundle.
- Metadata now records codec names (``ecomp`` vs ``phylo-bundle``) to simplify
  downstream aggregation.

Extending Benchmarks
--------------------

1. Collect MSAs of varying depth/length in ``example_data`` or ``synthetic_data``.
2. Attach matching tree files where available to exercise the bundle codec.
3. Regenerate ``benchmark_*.json`` artefacts with the helper above.

These reports feed into the documentation and release notes much like ClipKIT's
performance assessment notebooks.
