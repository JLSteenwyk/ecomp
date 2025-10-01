<p align="center">
  <a href="https://github.com/jlsteenwyk/ecomp">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/ecomp/master/docs/_static/img/logo.jpg" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/ecomp/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/ecomp/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/ecomp/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://github.com/JLSteenwyk/ecomp/actions" alt="Build">
            <img src="https://img.shields.io/github/actions/workflow/status/JLSteenwyk/ecomp/ci.yml?branch=master">
        </a>
        <a href="https://codecov.io/gh/jlsteenwyk/ecomp" alt="Coverage">
          <img src="https://codecov.io/gh/jlsteenwyk/ecomp/branch/master/graph/badge.svg?token=0J49I6441V">
        </a>
        <a href="https://github.com/jlsteenwyk/ecomp/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/ecomp">
        </a>
        <a href="https://bsky.app/profile/jlsteenwyk.bsky.social" target="_blank" rel="noopener noreferrer">
          <img src="https://img.shields.io/badge/Bluesky-0285FF?logo=bluesky&logoColor=fff">
        </a>
        <br />
        <a href="https://pepy.tech/badge/ecomp">
          <img src="https://static.pepy.tech/personalized-badge/cliecompkit?period=total&units=international_system&left_color=grey&right_color=blue&left_text=PyPi%20Downloads">
        </a>
        <!-- <a href="https://anaconda.org/bioconda/ecomp">
          <img src="https://img.shields.io/conda/dn/bioconda/ecomp?label=bioconda%20downloads" alt="Bioconda Downloads">
        </a> -->
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <br />
        <a href="https://pypi.org/project/ecomp/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/ecomp">
        </a>
        <!-- <a href="https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007">
          <img src="https://zenodo.org/badge/DOI/10.1371/journal.pbio.3001007.svg">  
        </a>    -->
    </p>
</p>


Evolution-informed lossless compression of multiple-sequence alignments (MSAs).

---

## Installation

From PyPI (recommended for users):

```bash
pip install ecomp
```

From source:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install .[dev]
```

> Offline? Pre-install `biopython`, `numpy`, `bitarray`, and the dev tools
> (`pytest`, `ruff`, `black`, `mypy`, …) inside your environment.

---

## CLI Quickstart

All commands are exposed through the `ecomp` entry point.

```bash
# Compress an alignment (produces example.ecomp, optional JSON sidecar)
ecomp zip example.fasta --metadata example.json

# Decompress (writes FASTA by default)
ecomp unzip example.ecomp --alignment-output restored.fasta

# Inspect metadata (summary or JSON)
ecomp inspect example.ecomp --summary

# Diagnostics (Phykit-style aliases in parentheses)
ecomp consensus_sequence example.ecomp             # con_seq
ecomp column_base_counts example.ecomp             # col_counts
ecomp gap_fraction example.ecomp                   # gap_frac
ecomp shannon_entropy example.ecomp                # entropy
ecomp parsimony_informative_sites example.ecomp    # parsimony
ecomp constant_columns example.ecomp               # const_cols
ecomp pairwise_identity example.ecomp              # pid
ecomp alignment_length_excluding_gaps example.ecomp    # len_no_gaps
ecomp alignment_compressed_length example.ecomp        # compressed_len
ecomp variable_sites example.ecomp                     # var_sites
ecomp percentage_identity example.ecomp                # pct_id
ecomp relative_composition_variability example.ecomp   # rcv
```

Benchmarks mirror standard codec comparisons:

```bash
/usr/bin/time -p ecomp zip data/fixtures/small_phylo.fasta --output out.ecomp
/usr/bin/time -p gzip  -k data/fixtures/small_phylo.fasta
/usr/bin/time -p bzip2 -k data/fixtures/small_phylo.fasta
```

---

## Python API

Everything the CLI does is re-exported in `ecomp`.

```python
from ecomp import zip, unzip, read_alignment, percentage_identity, column_base_counts

# File-based workflow
archive_path, metadata_path = zip(
    "data/example.fasta",
    metadata_path="data/example.json",  # optional JSON copy
)
restored_path = unzip(archive_path, output_path="data/restored.fasta")

# Diagnostics on an AlignmentFrame
frame = read_alignment("data/example.fasta")
pct_identity = percentage_identity(frame)
base_counts = column_base_counts(frame)

print(f"Mean pairwise identity: {pct_identity:.2f}%")
print("Column 1 counts:", base_counts[0])
```

In-memory usage (no intermediate files):

```python
from ecomp import AlignmentFrame, compress_alignment, decompress_alignment

frame = AlignmentFrame(
    ids=["s1", "s2"],
    sequences=["ACGT", "ACGA"],
    alphabet=["A", "C", "G", "T"],
)
compressed = compress_alignment(frame)
restored = decompress_alignment(compressed.payload, compressed.metadata)
assert restored.sequences == frame.sequences
```
