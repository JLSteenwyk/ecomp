.. role:: raw-html(raw)
    :format: html

.. image:: _static/img/logo_top_only.png 
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/ecomp

^^^^^

eComp uses evolution-based principles to conduct lossless compression of multiple sequence alignments.

The resulting files can be efficiently parsed for key information, such as the consensus sequence, number of parsimony-informative sites, and other metrics.

:raw-html:`<br />`

Quick start
===========

**Installation**

1. Install eComp:

   .. code-block:: bash

      # Create and activate a virtual environment (Python 3.11+)
      python3 -m venv venv
      source venv/bin/activate

      # Install dependencies and eComp (runtime + developer extras)
      pip install ecomp

:raw-html:`<br />`

**CLI Usage**

The ``ecomp`` command exposes all compression, decompression, and metrics:

.. code-block:: bash

   # Compress an alignment (produces example.ecomp)
   ecomp zip example.fasta

   # Decompress (writes FASTA by default)
   ecomp unzip example.ecomp --alignment-output restored.fasta

   # Inspect metadata (summary or JSON)
   ecomp inspect example.ecomp --summary

   # PhyKIT-style diagnostics (alias shown in parentheses)
   ecomp consensus_sequence example.ecomp         # (con_seq)
   ecomp pairwise_identity example.ecomp          # (pid)
   ecomp variable_sites example.ecomp             # (var_sites)

:raw-html:`<br />`

**API Usage**

All CLI functionality is re-exported via ``ecomp``:

.. code-block:: python

   from ecomp import read_alignment, zip, unzip

   frame = read_alignment("example.fasta")
   archive_path, metadata_path = zip(
       "example.fasta",
       metadata_path="example.json",  # optional JSON sidecar
   )
   restored_path = unzip(archive_path)

   from ecomp import percentage_identity
   pct_identity = percentage_identity(frame)

:raw-html:`<br />`:raw-html:`<br />`

Quick Start
===========

.. toctree::
   :maxdepth: 1

   about/index
   usage/index
   change_log/index
   other_software/index
   faq/index
