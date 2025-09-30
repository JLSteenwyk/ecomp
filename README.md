# Evolutionary Compression (`ecomp`)

`ecomp` is a Python toolkit for lossless compression and decompression of multiple
sequence alignments (MSAs) tuned for evolutionary genomics workflows. The
pipeline combines column-wise consensus discovery, deviation tracking, and
run-length encoding to produce compact `.ecomp` payloads plus JSON metadata.

## How eComp Differs From gzip and bzip2

- **gzip** streams bytes through a dictionary-based (DEFLATE) coder. Repeated
  byte sequences shorten, but the algorithm has no awareness of alignment
  columns or alphabets. Each FASTA character is treated in isolation, so long
  stretches of gaps and the repeating headers dominate the output size.
- **bzip2** blocks data, applies Burrows–Wheeler and Huffman coding, and then
  run-length encodes the transformed bytes. It improves on gzip for highly
  repetitive text, yet again operates purely on byte patterns without knowing
  which symbols belong to which sequence.
- **eComp** exploits MSA structure. Every column is summarised by a consensus
  residue plus a sparse list of deviations. The consensus stream and deviation
  bitmasks are run-length encoded, and the residual payload is pushed through a
  generic compressor (zstd/zlib/xz) only after the biological redundancy has
  been stripped away.

### Toy Example

Given a tiny alignment:

```
>seq1
ACGTACGT
>seq2
ACGTACGA
>seq3
ACGTACGG
```

- **gzip / bzip2**: serialize the FASTA exactly as written and compress the raw
  bytes. Repeated `ACGT` substrings help, but headers and per-sequence gaps are
  still stored explicitly. Deviations in the final column (A/G) appear as
  independent characters, so both compressors must encode them in full.
- **eComp**:
  1. Reads column 1 (`A/A/A`), stores consensus `A` and marks “no deviations”.
  2. Repeats until the final column where the consensus is `A` with two
     deviations (`seq2= A`, `seq3= G`). Only the consensus symbol and the two
     deviating residues are emitted; all other positions inherit the consensus
     for free.
  3. Packs the deviation bitmask (which sequences differ), run-length encodes
     blocks of identical columns, and finally compresses the already-compact
     payload with zstd/zlib/xz.

Because most columns are perfectly conserved, eComp emits a handful of bits per
column, whereas gzip/bzip2 still allocate bytes for every residue in every
sequence. On real MSAs (thousands of taxa × columns) this structural awareness
translates into 2–3× better compression versus canonical codecs.

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
# The CLI is available as both `ecomp` and the shorter alias `ec`.
```bash
# Compress an alignment (writes example.ecomp + metadata JSON)
ec compress example.fasta

# Optionally supply a tree to guide ordering (tree is not stored)
ec compress example.fasta --tree example.tree

# Decompress (auto-detects codec from metadata)
ec decompress example.ecomp --alignment-output restored.fasta

# Inspect metadata (JSON or short summary)
ec inspect example.ecomp --summary
```
The `ecomp` entry point mirrors the public Python API (`compress_file`, `decompress_file`, `compress_alignment`, `decompress_alignment`).

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
Use the CLI together with standard timing tools to compare eComp to other
codecs. A quick local comparison works with the bundled fixture:
```bash
/usr/bin/time -p ecomp compress data/fixtures/small_phylo.fasta \
  --output out.ecomp
/usr/bin/time -p gzip -k data/fixtures/small_phylo.fasta
```
Round-trip the archive to confirm correctness before discarding temporary files.

Larger alignments and manuscript figures are published alongside the paper in
`../EVOCOMP_MANUSCRIPT/` (or the companion data archive). Set an environment
variable such as `EVOCOMP_DATA_ROOT` to point at that directory when running the
workflows under `docs/tutorials/`.

Additional roadmap milestones and contributor practices are documented in
`ECOMP_CODEBASE_PLAN.md` and `AGENTS.md`.
