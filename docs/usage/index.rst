.. role:: raw-html(raw)
    :format: html

Usage
=====

This section collects the most common workflows for working with eComp. Start with
the installation steps below, then explore CLI examples and the Python API. The
pages in the toctree dive deeper into quickstarts, metrics, benchmarking, and
programmatic usage.

:raw-html:`<br />`

Installation
============

.. code-block:: bash

   # 1) Create and activate a virtual environment (Python 3.11+)
   python3 -m venv venv
   source venv/bin/activate

   # 2) Install dependencies and eComp (runtime + developer extras)
   pip install -r requirements.txt
   pip install .[dev]

   # 3) (Optional) run the fast test suite to confirm setup
   make test.fast

:raw-html:`<br />`

CLI Usage Overview
==================

The ``ecomp`` command exposes compression, decompression, inspection, and
PhyKIT-style diagnostics. Short aliases (mirroring the CLI verbs) are shown in
parentheses.

.. code-block:: bash

   # Archive an alignment (produces example.ecomp + metadata sidecar if requested)
   ecomp zip example.fasta --metadata example.json

   # Restore FASTA (automatically infers the output extension)
   ecomp unzip example.ecomp --alignment-output restored.fasta

   # Display metadata (summary or full JSON)
   ecomp inspect example.ecomp --summary
   ecomp inspect example.ecomp

   # Diagnostics and metrics
   ecomp consensus_sequence example.ecomp             # (con_seq)
   ecomp column_base_counts example.ecomp             # (col_counts)
   ecomp gap_fraction example.ecomp                   # (gap_frac)
   ecomp shannon_entropy example.ecomp                # (entropy)
   ecomp parsimony_informative_sites example.ecomp    # (parsimony)
   ecomp constant_columns example.ecomp               # (const_cols)
   ecomp pairwise_identity example.ecomp              # (pid)
   ecomp alignment_length_excluding_gaps example.ecomp    # (len_no_gaps)
   ecomp alignment_compressed_length example.ecomp        # (compressed_len)
   ecomp variable_sites example.ecomp                     # (var_sites)
   ecomp percentage_identity example.ecomp                # (pct_id)
   ecomp relative_composition_variability example.ecomp   # (rcv)

:raw-html:`<br />`

Python API Highlight
====================

All CLI behaviour is accessible via ``ecomp``. Core helpers are grouped below for
quick reference (see :doc:`api` for full details).

File-based workflow
-------------------

.. code-block:: python

    from ecomp import zip, unzip, read_alignment, percentage_identity, column_base_counts

    # Compress an alignment file (optionally mirroring the CLI flags)
    archive_path, metadata_path = zip(
        "data/example.fasta",
        metadata_path="data/example.json",  # optional JSON copy
    )

    # Decompress, writing FASTA to disk
    restored_path = unzip(
        archive_path,
        output_path="data/restored.fasta",
    )

    # Load into an AlignmentFrame to run diagnostics
    frame = read_alignment("data/example.fasta")
    pct_identity = percentage_identity(frame)
    base_counts = column_base_counts(frame)

    print(f"Mean pairwise identity: {pct_identity:.2f}%")
    print(f"Column counts for column 1: {base_counts[0]}")

In-memory usage
---------------

.. code-block:: python

    from ecomp import AlignmentFrame, compress_alignment, decompress_alignment

    frame = AlignmentFrame(
        ids=["s1", "s2"],
        sequences=["ACGT", "ACGA"],
        alphabet=["A", "C", "G", "T"],
    )
    compressed = compress_alignment(frame)
    restored = decompress_alignment(compressed.payload, compressed.metadata)
    assert restored.sequences == frame.sequences

.. toctree::
   :maxdepth: 1

   quickstart
   api
   benchmarking
   ecomp_metrics
