Change Log
==========

.. note::
   Keep this file brief and focused on user-facing changes.  The project has not
   tagged a formal release yet; the entry below captures the current milestone.

Unreleased
----------

- Implemented consensus + run-length pipeline with tree-guided ordering.
- Added CLI commands (``compress``, ``decompress``, ``inspect``) and matching Python API.
- Introduced ``ec`` command alias for quicker typing.
- Added in-house sequence reordering heuristics (MST + greedy) with automatic selection; configurable via the ``ECOMP_SEQUENCE_ORDER`` environment variable.
- Published benchmarking harness and documentation site scaffold (mirrors PhyKIT layout).
- Tree metadata now accepted via ``frame.metadata['tree_newick']`` to improve compression ratios.
