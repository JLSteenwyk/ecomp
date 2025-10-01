Per-Column Metrics on ``.ecomp`` Archives
========================================

This note sketches a first batch of analytics routines that can operate directly
on the evolutionary compression (``.ecomp``) payloads.  The guiding idea is to
stream over the compressed column blocks to avoid materialising full FASTA
alignments when you only need summary statistics.

Shared Inputs
-------------

* ``payload``: byte stream returned by :func:`ecomp.storage.read_archive`.
* ``metadata``: dictionary with column and sequence counts, codec hints, and any
  embedded tree information.
* ``AlignmentFrame`` helpers (from :mod:`ecomp.io`), which
  expose decoded columns and metadata when a higher-level iterator is more
  convenient than working with raw blocks.

Command-Line Access
-------------------

Each metric is exposed directly via the CLI with descriptive names (and short
aliases for parity with PhyKIT): ``consensus_sequence`` (``con_seq``),
``column_base_counts`` (``col_counts``), ``gap_fraction`` (``gap_frac``),
``shannon_entropy`` (``entropy``), ``parsimony_informative_sites``
(``parsimony``), ``constant_columns`` (``const_cols``), ``pairwise_identity``
(``pid``), ``alignment_length_excluding_gaps`` (``len_no_gaps``),
``alignment_compressed_length`` (``compressed_len``), ``variable_sites``
(``var_sites``), ``percentage_identity`` (``pct_id``), and
``relative_composition_variability`` (``rcv``).  All accept the path to a
``.ecomp`` archive (plus optional ``--metadata`` for legacy sidecars) and print
tabular or JSON summaries to stdout.

Proposed Functions
------------------

Per-column base counts
    Traverse each column block, tally residue characters (A/C/G/T or amino-acid
    alphabet), and return a list of ``Counter`` objects or a ``numpy`` array of
    shape ``(columns, alphabet_size)``.  These counts underpin several of the
    downstream metrics, so design the routine to emit both raw tallies and a
    normalised frequency view.

Majority-rule consensus
    Build on the base-count iterator to select the most frequent non-gap symbol
    per column.  Ties mirror PhyKIT's behaviour: nucleotide mixtures are
    represented with IUPAC ambiguity codes (e.g., ``R``, ``Y``), while protein
    mixtures fall back to ``X``.  Emit a FASTA record or a simple string paired
    with a per-column support score.

Gap fraction
    Compute the ratio of gap characters to total sequences for each column.  The
    result can be used to mask gappy regions, feed trimming heuristics, or
    populate summary tables similar to PhyKIT's ``gap_fraction`` output.

Shannon entropy
    Using the frequency table from the base-count pass, evaluate
    ``entropy = -Σ p_i log2(p_i)`` per column, skipping zero frequencies.  High
    entropy points to variable sites, while low entropy highlights conserved
    columns.

Parsimony-informative sites
    Classify columns where at least two residues occur with minimum frequency of
    two sequences each.  Consume the base counts and flag boolean indicators or
    return the absolute count of qualifying columns.

Constant columns
    Flag columns where exactly one non-gap residue is present (after ignoring
    gaps).  This can either reuse the parsimony-informative scan or be exposed as
    a dedicated helper that also reports which residue dominates.

Pairwise identity
    Derive an ``n × n`` matrix (or condensed vector) capturing the fraction of
    identical positions between every pair of sequences.  To keep the pass
    scalable, stream column by column, increment matches for pairs that share the
    same non-gap residue, and maintain a coverage counter for each comparison to
    handle gaps gracefully.  Combine with metadata ordering information to emit a
    matrix aligned with the original sequence order.

Alignment summaries
    Mirror PhyKIT's aggregate metrics by computing: ``alignment_length_excluding_gaps``
    (count of columns with at least one non-gap residue),
    ``alignment_compressed_length`` (columns remaining after removing constant
    and all-gap sites), ``variable_sites`` (columns with more than one residue
    once gaps are ignored), ``percentage_identity`` (mean of all pairwise
    identities expressed as a percentage), and ``relative_composition_variability``
    (RCV) which averages per-sequence compositional deviations.

Implementation Notes
--------------------

* Cache the alphabet derived from ``metadata['alphabet']`` when available to
  avoid repeated detection of residue categories.
* Prefer streaming generators that yield one column at a time; most metrics can
  be layered without storing the entire table in memory.
* Consider returning ``xarray`` / ``pandas`` objects for downstream plotting, but
  keep the core functions dependency-light so they can run in minimal
  environments.
* Expose a small CLI shim (for example ``codex metrics``) once these helpers are
  implemented to mirror PhyKIT-style subcommands.

This document will expand as additional analytics—such as alignment-effective
length, homoplasy scores, or tree-linked diagnostics—are scoped.
