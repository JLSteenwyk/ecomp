Quickstart
==========

Installation
------------

Use a virtual environment (recommended) and install the dependencies listed in
``requirements.txt``.  The project supports Python 3.11–3.13.

.. code-block:: bash

    python -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    pip install .

The optional development extras (linting, tests) can be installed via

.. code-block:: bash

    pip install .[dev]

Command-line usage
------------------

Three sub-commands cover the common tasks (legacy ``ecomp`` / ``ec`` aliases
remain for convenience):

.. code-block:: bash

    # Zip an alignment (auto-detect output names)
    codex zip data/alignment.fasta

    # Supply a tree only for sequence ordering (tree is NOT stored)
    codex zip data/alignment.fasta --tree data/alignment.fasta.tre

eComp now checks a few alignment heuristics (gap density, relative composition
variability, and pairwise identity spread) before using the tree ordering.  If
the data look noisy, it falls back to the similarity-based reorder automatically;
set ``ECOMP_SEQUENCE_ORDER=tree`` to force tree usage when you know it is safe.

    # Unzip back to FASTA (auto-picks .fasta extension)
    codex unzip data/alignment.ecomp

    # Inspect metadata without full decompression
    codex inspect data/alignment.ecomp --summary

    # Alignment diagnostics (aliases mirror PhyKIT-style shorthand)
    codex consensus_sequence data/alignment.ecomp
    codex column_base_counts data/alignment.ecomp
    codex gap_fraction data/alignment.ecomp
    codex shannon_entropy data/alignment.ecomp
    codex parsimony_informative_sites data/alignment.ecomp
    codex constant_columns data/alignment.ecomp
    codex pairwise_identity data/alignment.ecomp
    codex alignment_length_excluding_gaps data/alignment.ecomp
    codex alignment_compressed_length data/alignment.ecomp
    codex variable_sites data/alignment.ecomp
    codex percentage_identity data/alignment.ecomp
    codex relative_composition_variability data/alignment.ecomp

To forward arguments explicitly (for example, to choose an output name):

.. code-block:: bash

    codex zip alignment.fa --tree alignment.fa.tre \
        -o results/alignment.ecomp -m results/alignment.json
    codex unzip results/alignment.ecomp -m results/alignment.json -o restored.fa

Python API
----------

If you prefer to stay in Python, import the helpers described in
:doc:`api`. They wrap the same logic as the CLI and interoperate with Biopython's
``MultipleSeqAlignment`` objects via ``read_alignment`` / ``AlignmentFrame``.
