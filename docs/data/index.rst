Data Resources
==============

``data_large/``
  Real-world benchmark panels (e.g., Ascomycota, Bees, Birds).  Use these to
  reproduce the figures in the benchmarking section.  Large files are tracked in
  Git LFS or supplied separately—check the repository README for download
  instructions.

``example_data/``
  Small toy alignments used in unit tests and quick tutorials.  Safe to commit
  and share.

``scripts/prepare_fixtures.py`` *(planned)*
  Placeholder for automated dataset downloads.  Until it lands, manage external
  datasets manually and document provenance in your benchmark reports.

When adding new datasets, include a short README describing the source,
licensing, preprocessing steps, and any accompanying trees or partitions.
