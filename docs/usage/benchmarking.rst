Benchmarking
============

The standalone benchmarking helpers previously shipped under ``scripts/`` have
been retired to keep the release footprint lean. You can still capture the same
signal by timing the CLI directly and parsing its JSON metadata.

Quick comparisons
-----------------

Use your shell's timing utility to compare eComp to other codecs:

.. code-block:: bash

    /usr/bin/time -p ecomp compress data/fixtures/small_phylo.fasta \
        --output out.ecomp
    /usr/bin/time -p gzip -k data/fixtures/small_phylo.fasta

The ``ecomp`` command prints archive size and compression ratio; record the wall
clock figures reported by ``time`` for performance tracking.

Decompression checks
--------------------

Always validate the round trip after timing a run:

.. code-block:: bash

    ecomp decompress out.ecomp --output restored.fasta
    diff -u data/fixtures/small_phylo.fasta restored.fasta

Supplementary data
-------------------

The manuscript archive (``../EVOCOMP_MANUSCRIPT`` or the published supplement)
contains large-scale benchmarking tables and example alignments. Point
``EVOCOMP_DATA_ROOT`` at that directory to reuse the workflows below.

- ``data_large_benchmark_analysis.csv`` — runtime and ratio data for the
  ``data_large`` corpus.
- ``hmmer_bwt_comparison.json`` — reference results contrasting HMMER and eComp.

Reuse these files for plots or include new measurements alongside them to retain
consistency between releases.