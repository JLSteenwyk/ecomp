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

Three sub-commands cover the common tasks (the ``ec`` alias is provided for
convenience):

.. code-block:: bash

    # Compress an alignment (auto-detect output names)
    ec compress data/alignment.fasta

    # Supply a tree only for sequence ordering (tree is NOT stored)
    ec compress data/alignment.fasta --tree data/alignment.fasta.tre

    # Decompress back to FASTA (auto-picks .fasta extension)
    ec decompress data/alignment.ecomp

    # Inspect metadata without full decompression
    ec inspect data/alignment.ecomp --summary

To forward arguments explicitly (for example, to choose an output name):

.. code-block:: bash

    ec compress alignment.fa --tree alignment.fa.tre \
        -o results/alignment.ecomp -m results/alignment.json
    ec decompress results/alignment.ecomp -m results/alignment.json -o restored.fa

Python API
----------

If you prefer to stay in Python, import the helpers described in
:doc:`api`. They wrap the same logic as the CLI and interoperate with Biopython's
``MultipleSeqAlignment`` objects via ``read_alignment`` / ``AlignmentFrame``.
