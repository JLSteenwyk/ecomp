# Evolutionary Compression (`ecomp`)

`ecomp` is a Python toolkit for lossless compression and decompression of multiple
sequence alignments (MSAs) tuned for evolutionary genomics workflows. The
pipeline combines column-wise consensus discovery, deviation tracking, and
run-length encoding to produce compact `.ecomp` payloads plus JSON metadata.

## Quickstart
```bash
python3 -m venv venv --system-site-packages  # or omit flag if you can install deps
source venv/bin/activate
pip install -r requirements.txt
pip install .[dev]
```
> NOTE: Installing dependencies requires outbound network access to PyPI.
> If the environment is offline, ensure `biopython`, `numpy`, `bitarray`, and
> dev tools (`pytest`, `ruff`, `black`, `mypy`, etc.) are provisioned manually.

## CLI Usage
```bash
# Compress an alignment (writes example.ecomp + metadata JSON)
ecomp compress example.fasta

# Compress alignment + tree bundle (writes example.ecbt + metadata)
ecomp compress example.fasta example.tree

# Decompress (auto-detects codec from metadata)
ecomp decompress example.ecbt --alignment-output restored.fasta --tree-output restored.tree

# Inspect metadata (JSON or short summary)
ecomp inspect example.ecomp --summary
```
The `ecomp` entry point mirrors the public Python API (`compress_file`, `decompress_file`, `compress_alignment_with_tree`, `decompress_alignment_with_tree`).

## Development Workflow
- Run the fast test suite (unit + non-slow integration):
  ```bash
  make test.fast
  ```
- Execute the full test matrix:
  ```bash
  make test
  ```
- Generate coverage reports for Codecov uploads:
  ```bash
  make test.coverage
  ```
- Lint and format:
  ```bash
  make lint
  make format
  ```
- Type-check:
  ```bash
  mypy src
  ```
- Pre-commit:
  ```bash
  pre-commit install
  pre-commit run --all-files
  ```

## Benchmarking
Use `scripts/compare_compressors.py` to compare eComp against `gzip`, `bzip2`,
`xz`, and the experimental `phylo-bundle` codec (alignment + Newick tree when
available):
```bash
python scripts/compare_compressors.py data/fixtures/small_orthogroup.fasta \
  --codecs ecomp phylo-bundle gzip bzip2 xz --output benchmark.json
```
If a companion tree file (e.g., `small_orthogroup.tree`) exists alongside an
alignment, the script will automatically include the bundled codec in the
results and compute ratios using the combined alignment + tree size.
The script emits compression ratios and runtime metrics in JSON format.

Additional roadmap milestones and contribution practices are documented in
`ECOMP_CODEBASE_PLAN.md` and `AGENTS.md`.
