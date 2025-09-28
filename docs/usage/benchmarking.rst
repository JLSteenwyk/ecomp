Benchmarking
============

This project provides two helper scripts for reproducible measurements.

compare_compressors.py
----------------------

``scripts/compare_compressors.py`` measures compression ratio, archive size, and
run time for any combination of codecs:

.. code-block:: bash

    python scripts/compare_compressors.py alignment1.fasta alignment2.fasta \
        --codecs ecomp gzip bzip2 xz --runner internal --output docs/bench.json

Key options:

- ``--runner cli`` forces the script to call the ``ecomp`` CLI instead of the
  Python API (useful when profiling end-to-end behaviour).
- When ``--runner internal`` is used, the script automatically feeds companion
  ``*.tre`` files (if present) into eComp to guide sequence ordering.
- ``--output`` writes a JSON summary suitable for plotting or archiving.

You can fix the reordering heuristic by setting the ``ECOMP_SEQUENCE_ORDER``
environment variable to ``baseline``, ``mst``, or ``greedy`` before running the
script.  Leave it unset (the default) to let eComp choose the lowest-cost ordering
automatically.

Tree-ordering experiments
-------------------------

``scripts/benchmark_tree_ordering.py`` compares three modes for each alignment:

1. Baseline (no tree metadata).
2. Tree-guided ordering using a provided Newick file.
3. Neighbor-joining tree inferred on the fly (the default behaviour when no tree
   is present).

.. code-block:: bash

    python scripts/benchmark_tree_ordering.py --output docs/tree_ordering.json \
        data_large/ascomycota_data/*.fasta

The output JSON records compressed size, metadata bytes, and timing for each
scenario, making it easy to quantify how much the phylogeny helps.  Supplied trees
are used only to guide ordering—they are not stored in the resulting archives.

Reporting tips
--------------

When publishing results, capture the command lines, dataset description (taxa,
columns, alphabet), the measured ratios for eComp and the baseline codecs, and any
notes about tree usage.  Store the JSON artifacts under ``docs/`` so they can be
referenced from tutorials or the change log.

For full end-to-end automation examples see :doc:`../tutorials/workflows`.
