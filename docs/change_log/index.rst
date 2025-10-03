Change Log
==========

v0.1.0
------

* Renamed the public compression helpers to ``ezip`` / ``eunzip`` to avoid
  shadowing Python's built-in ``zip``.
* Added ``alignment_length`` as the canonical total column metric and removed
  the previous ``alignment_compressed_length`` helper.
* Expanded the CLI metrics suite with matching coverage tests and documentation
  updates.
* Hardened compatibility shims and added unit coverage for Python 3.9 support.
* Introduced ``distance_tree`` helpers (API + CLI) to infer distance-based
  trees directly from ``.ecomp`` archives and emit Newick output.


v0.1.2
------

* Enhanced the CLI UX with an ASCII banner, grouped command listings, and
  ``--list-commands`` / ``--version`` flags.
* Added richer unit coverage for the phylogenetics helpers and CLI utilities.
* General documentation and README updates to reflect the new interface.
