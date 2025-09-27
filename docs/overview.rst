Overview
========

|project| provides lossless compression for multiple sequence alignments (MSAs)
with evolutionary context in mind. The core pipeline discovers consensus
columns, tracks deviations, and encodes the result with a compact binary format.
A companion codec bundles alignments with their phylogenetic trees so that both
artefacts can be archived together.

Highlights
----------

- **Alignment codec** – column-wise consensus discovery with run-length blocks
  and opportunistic entropy coding via zlib.
- **Phylo bundle** – joint storage of an alignment and its Newick tree using a
  Fitch-style reconstruction of ancestral states.
- **CLI-first workflow** – the ``ecomp`` command mirrors ClipKIT's subcommand
  layout (``compress``, ``decompress``, ``inspect``) to keep muscle memory
  intact.
- **Python API** – ``evolutionary_compression`` exposes ``compress_file``,
  ``decompress_file``, and tree-aware helpers for scripted workflows.

Getting Started
---------------

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   pip install .[dev]

Once installed, run ``make help`` to see frequently used developer commands.
