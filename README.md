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

### How the Toy Alignment Is Stored

To see the storage format end to end, consider this 3x8 alignment:

```
>seq1
ACGTACGA
>seq2
ACGTTCGA
>seq3
ACGTACGG
```

- **Sequence ordering**: eComp may permute sequences to cluster similar rows.
  Here the input order is already good, so no permutation needs to be written.

- **Column modelling**: each column becomes a consensus residue plus a mask of
  which sequences deviate. The mask uses the order `(seq1, seq2, seq3)` where a
  `1` marks a deviation and the payload column lists those explicit residues.

| Column | Residues (seq1/seq2/seq3) | Consensus | Deviation mask | Payload |
|--------|---------------------------|-----------|----------------|---------|
| 1      | A/A/A                     | A         | 000            | -       |
| 2      | C/C/C                     | C         | 000            | -       |
| 3      | G/G/G                     | G         | 000            | -       |
| 4      | T/T/T                     | T         | 000            | -       |
| 5      | A/T/A                     | A         | 010            | T       |
| 6      | C/C/C                     | C         | 000            | -       |
| 7      | G/G/G                     | G         | 000            | -       |
| 8      | A/A/G                     | A         | 001            | G       |

- **Run-length bundling**: contiguous columns sharing the same consensus and
  mask are merged. Columns 1–4 become a single run with length 4, columns 6–7
  form another run of length 2, and the columns with deviations (5 and 8) are
  stored as single-column runs.

- **What gzip sees**: gzip streams the FASTA bytes exactly as they appear. Every
  header line, newline, and residue flows into the DEFLATE encoder with no
  column awareness. The same alignment therefore expands to 24 residue bytes
  plus six header/newline blocks before compression kicks in.

| Stream chunk | Literal bytes | Description |
|--------------|---------------|-------------|
| `>seq1\n`    | 6             | Header for `seq1` (leading `>` and trailing newline). |
| `ACGTACGA\n` | 9             | Row for `seq1`; every nucleotide stays literal. |
| `>seq2\n`    | 6             | Header for `seq2`. |
| `ACGTTCGA\n` | 9             | Row for `seq2`; mismatch `T` is just another byte. |
| `>seq3\n`    | 6             | Header for `seq3`. |
| `ACGTACGG\n` | 9             | Row for `seq3`; terminal `G` stored independently. |

  Adding those chunks yields 45 literal bytes before DEFLATE coding.

- **Binary layout**: the `.ecomp` payload records the sequence dictionary, the
  consensus stream (`ACGTACGA`), the deviation bitmasks (`0000 010 00 001`), the
  payload residues (`T` and `G`), and run-length counters (4,1,2,1). Only after
  these structures are assembled does eComp hand the compact stream to zstd
  (or another backend) for generic compression. The companion JSON stores the
  alignment dimensions, alphabet, checksum, backend codec, and the no-op
  permutation so the decompressor can rebuild the exact FASTA.

  Those structures total eight consensus bytes, two deviation bytes, and four
  run-length integers before the backend codec—around one third of the literal
  payload gzip must stream.

Compared with gzip, which must repeat every residue and header verbatim, eComp
stores the consensus once, notes the sparse deviations, and reuses run lengths
to avoid re-emitting conserved columns. Even this toy alignment shrinks from 24
characters down to a handful of bytes before the final backend compression is
applied.

### Comparing Archive Layouts (gzip vs. eComp)

```
alignment.fasta          # source file (FASTA or PHYLIP)

# gzip
alignment.fasta.gz       # single compressed stream, no extra metadata

# eComp
alignment.ecomp          # binary payload with consensus + run-length stream
alignment.json           # metadata sidecar (alignment stats, checksums, ordering)
```

gzip serializes the raw FASTA bytes and applies DEFLATE; decompression restores
exactly the original text but requires no added structure awareness. eComp first
reorders sequences (if beneficial), factors columns into consensus models, packs
deviations, and only then applies a general-purpose backend (zstd/zlib/xz). The
payload lives in `alignment.ecomp`, while `alignment.json` records everything
needed to reconstruct the alignment (dimensions, alphabet, permutation, checksum
and payload codec). Together these two files guarantee a lossless round trip
even when the compressor chooses different ordering or encoding strategies per
dataset.

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
# The CLI is available as `ecomp`.
```bash
# Zip an alignment (writes example.ecomp + metadata JSON)
ecomp zip example.fasta --metadata example.json  # metadata flag optional

# Optionally supply a tree to guide ordering (tree is not stored)
ecomp zip example.fasta --tree example.tree

# Unzip (auto-detects codec from metadata)
ecomp unzip example.ecomp --alignment-output restored.fasta

# Inspect metadata (JSON or short summary)
ecomp inspect example.ecomp --summary

# Alignment diagnostics (Phykit-style names with short aliases)
ecomp consensus_sequence example.ecomp             # alias: con_seq
ecomp column_base_counts example.ecomp             # alias: col_counts
ecomp gap_fraction example.ecomp                   # alias: gap_frac
ecomp shannon_entropy example.ecomp                # alias: entropy
ecomp parsimony_informative_sites example.ecomp    # alias: parsimony
ecomp constant_columns example.ecomp               # alias: const_cols
ecomp pairwise_identity example.ecomp              # alias: pid
ecomp alignment_length_excluding_gaps example.ecomp    # alias: len_no_gaps
ecomp alignment_compressed_length example.ecomp        # alias: compressed_len
ecomp variable_sites example.ecomp                     # alias: var_sites
ecomp percentage_identity example.ecomp                # alias: pct_id
ecomp relative_composition_variability example.ecomp   # alias: rcv
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
  mypy ecomp
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
/usr/bin/time -p ecomp zip data/fixtures/small_phylo.fasta \
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
