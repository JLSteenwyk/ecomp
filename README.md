<p align="center">
  <a href="https://github.com/jlsteenwyk/ecomp">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/ecomp/master/docs/_static/img/logo_transparent_background.png" alt="Logo" width="400">
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
            <img src="https://img.shields.io/github/actions/workflow/status/JLSteenwyk/ecomp/ci.yml?branch=main">
        </a>
        <a href="https://codecov.io/gh/jlsteenwyk/ecomp" alt="Coverage">
          <img src="https://codecov.io/gh/jlsteenwyk/ecomp/branch/main/graph/badge.svg?token=0J49I6441V">
        </a>
        <a href="https://github.com/jlsteenwyk/ecomp/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/ecomp">
        </a>
        <a href="https://bsky.app/profile/jlsteenwyk.bsky.social" target="_blank" rel="noopener noreferrer">
          <img src="https://img.shields.io/badge/Bluesky-0285FF?logo=bluesky&logoColor=fff">
        </a>
        <br />
        <a href="https://pepy.tech/badge/ecomp">
          <img src="https://static.pepy.tech/personalized-badge/ecomp?period=total&units=international_system&left_color=grey&right_color=blue&left_text=PyPi%20Downloads">
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
# create and activate a virtual environment
python -m venv venv
source venv/bin/activate

# install ecomp
pip install ecomp
```

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
ecomp alignment_length example.ecomp                   # len_total
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
from ecomp import ezip, eunzip, read_alignment, percentage_identity, column_base_counts

# File-based workflow
archive_path, metadata_path = ezip(
    "data/example.fasta",
    metadata_path="data/example.json",  # optional JSON copy
)
restored_path = eunzip(archive_path, output_path="data/restored.fasta")

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

### Available functions

**Compression & I/O** — `ezip`, `eunzip`, `compress_file`, `decompress_file`,
`compress_alignment`, `decompress_alignment`, `read_alignment`,
`write_alignment`, `alignment_from_sequences`, `alignment_checksum`

**Diagnostics & metrics** — `column_base_counts`, `column_gap_fraction`,
`column_shannon_entropy`, `parsimony_informative_columns`,
`parsimony_informative_site_count`, `constant_columns`,
`majority_rule_consensus`, `alignment_length`,
`alignment_length_excluding_gaps`, `variable_site_count`, `percentage_identity`,
`relative_composition_variability`, `pairwise_identity_matrix`

**Supporting types** — `AlignmentFrame`, `CompressedAlignment`,
`PairwiseIdentityResult`, `__version__`

---

## Development

```bash
make test.fast        # unit + non-slow integration tests
make test             # full test matrix
make lint             # lint checks (ruff, black, isort)
make format           # auto-formatting
mypy ecomp            # optional type checking
```

Build docs locally:

```bash
make docs
open docs/_build/html/index.html
```

Build and publish distributions:

```bash
pip install build twine
python -m build
python -m twine check dist/*
python -m twine upload dist/*
```

### Benchmarking eComp vs. PhyKIT

```bash
python scripts/benchmark_metrics.py data/example.ecomp \
    --operations consensus shannon_entropy variable_sites \
    --repeat 5 --warmup 1 --json results.json --csv results.csv
```

The script runs each metric via the `ecomp` CLI (on the compressed archive) and
the corresponding `phykit` command on a decompressed alignment, then reports
average and best runtimes. Add `--json`/`--csv` to emit machine-readable output.

---

## License

eComp is released under the MIT License. See [`LICENSE`](LICENSE).
