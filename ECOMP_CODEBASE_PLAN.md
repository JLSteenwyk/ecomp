# eComp Codebase Build Plan

## Guiding Objectives
- Deliver a lossless MSA compressor/decompressor optimized for evolutionary datasets, producing `.ecomp` bundles plus JSON metadata.
- Maintain transparency and reproducibility through deterministic algorithms, explicit metadata, and round-trip validation.
- Provide an ergonomic CLI and Python API suitable for workflow automation and benchmarking against standard compressors.

## Repository Layout Blueprint
```
/
├── src/
│   └── evolutionary_compression/
│       ├── __init__.py
│       ├── cli.py
│       ├── config.py
│       ├── io.py
│       ├── compression/
│       │   ├── __init__.py
│       │   ├── consensus.py
│       │   ├── rle.py
│       │   ├── encoding.py
│       │   └── pipeline.py
│       └── diagnostics/
│           ├── __init__.py
│           ├── checksums.py
│           └── benchmarks.py
├── tests/
│   ├── unit/
│   ├── integration/
│   └── regression/
├── data/
│   ├── fixtures/
│   └── benchmarks/
├── docs/
│   ├── architecture.md
│   └── api_reference.md
├── scripts/
│   ├── prepare_fixtures.py
│   └── compare_compressors.py
├── pyproject.toml
├── requirements.txt
└── README.md
```

## Phase Roadmap & Deliverables
1. **Phase 0 – Environment & Scaffolding (Week 1)**
   - Initialize Python package with `pyproject.toml`, `src/` layout, and tooling (`black`, `ruff`, `isort`, `mypy`).
   - Configure `pytest`, coverage (`coverage.py`), and pre-commit hooks.
   - Acceptance: `pytest` sample test passes; `pre-commit run --all-files` clean; CI placeholder workflow drafted.
2. **Phase 1 – I/O Foundations (Weeks 1–2)**
   - Implement `io.py` to load FASTA/PHYLIP via Biopython with schema validation and canonical in-memory representation (`AlignmentFrame`).
   - Define serialization contracts for `.ecomp` payload and metadata (JSON schema in `docs/` and `src/evolutionary_compression/config.py`).
   - Acceptance: Round-trip tests for input parsing and JSON schema validation.
3. **Phase 2 – Compression Core (Weeks 2–4)**
   - Build consensus extraction (`consensus.py`), deviation encoding, and columnar run-length encoding (`rle.py`).
   - Implement binary packing utilities (`encoding.py`) leveraging `bitarray` or `numpy`.
   - Establish pipeline orchestrator (`pipeline.py`) that composes steps and yields structured blocks.
   - Acceptance: `tests/unit/test_consensus.py` and `tests/unit/test_rle.py` achieve ≥95% coverage; property tests ensure idempotent encode/decode cycles for synthetic alignments.
4. **Phase 3 – Decompression Symmetry (Weeks 4–5)**
   - Implement inverse operations mirroring compression pipeline with explicit invariants.
   - Add checksum verification (`diagnostics/checksums.py`) using SHA256 and optional per-record checks.
   - Acceptance: Integration tests (`tests/integration/test_round_trip.py`) verify byte-for-byte equality on fixtures; mismatched checksum raises descriptive error.
5. **Phase 4 – CLI & API Surface (Week 5)**
   - Build Click/Typer-based CLI in `cli.py` exposing `compress`, `decompress`, and `inspect` commands.
   - Provide Python API wrappers in `__init__.py` (e.g., `compress_file`, `decompress_file`, `compress_alignment`).
   - Acceptance: CLI smoke tests (`pytest --pyargs evolutionary_compression`) run in CI; documentation includes command examples.
6. **Phase 5 – Diagnostics & Benchmarking (Weeks 6–7)**
   - Create benchmarking harness (`diagnostics/benchmarks.py` and `scripts/compare_compressors.py`) to evaluate speed/ratio vs. gzip, bzip2, xz.
   - Record benchmark protocol in `docs/architecture.md`; support CSV/Markdown reporting.
   - Acceptance: Benchmarks run against sample datasets with results captured in `docs/benchmark_report.md`.
7. **Phase 6 – Hardening & Release Preparation (Week 8)**
   - Add regression fixtures (edge cases, large alignments, non-standard alphabets).
   - Implement configuration overrides (chunk size, parallelism) and environment variable support.
   - Finalize documentation (`README`, `docs/api_reference.md`), license, and versioning policy.
   - Acceptance: `pytest --cov` ≥90%, lint passes, CLI help text tested, release checklist completed.

## Detailed Workstreams
### Data Model & Serialization
- Define `AlignmentFrame` dataclass capturing sequences, taxon order, alphabet, and annotations.
- Specify `.ecomp` container: header magic string, version, alignment dimensions, consensus block index, compressed payload, checksum.
- Use `struct` for binary packing; store deviations as variable-length records with prefix codes (sequence index, residue encoded in 5 bits for DNA/protein).
- Document schema in `docs/architecture.md` with byte offsets and upgrade strategy.

### Compression Logic Enhancements
- Optimize consensus detection with NumPy counts; fallback to Python for ambiguous characters (`-`, `N`).
- Implement block-level RLE to combine consecutive columns sharing consensus and identical deviation sets.
- Add optional entropy coding extension point (stub) for future Huffman/ANS modules.
- Introduce streaming interface for large alignments using generators and memory-mapped buffers.

### Testing Strategy
- Unit: Focus on deterministic outputs for consensus selection, RLE packing/unpacking, metadata validation.
- Property: Hypothesis-based tests generating random alignments with controlled mutation rates.
- Integration: Round-trip tests covering FASTA/PHYLIP, varying taxon counts, gaps, and ambiguous residues.
- Regression: Reproduce known bug scenarios; lock fixtures using Git LFS if files exceed 100 MB (request user approval before enabling).
- Performance: Automated benchmark suite gated behind `PYTEST_ADDOPTS="-m 'not benchmark'"` with optional `pytest-benchmark` marks.

### Tooling & Automation
- Configure GitHub Actions workflow with matrix over Python 3.11/3.12, running lint, mypy, tests, and packaging build.
- Provide `make` or `tox` entry points (`make test`, `make lint`, `tox -e pytest`) for consistent execution.
- Add `pre-commit` config covering formatting, linting, YAML checks, and docstring style (e.g., `pydocstyle` optional).
- Set up `dependabot.yml` for dependency updates and `codespell` to catch typos in metadata/docstrings.

### Documentation & Knowledge Base
- Maintain `docs/architecture.md` (high-level diagrams, data flow), `docs/storage_format.md` (binary layout), and `docs/benchmark_report.md`.
- Update `README` with quickstart, CLI usage, and performance summary once benchmarks exist.
- Embed docstrings and type hints to enable auto-generated API docs (future Sphinx or MkDocs integration).

### Data Management
- Curate small canonical fixtures in `data/fixtures/` (e.g., `small_orthogroup.fasta`, `protein_family.phylip`).
- Store larger benchmark datasets in `data/benchmarks/` with metadata README covering source, license, preprocessing steps.
- Provide `scripts/prepare_fixtures.py` to download/validate external datasets when permissible; ensure offline operation by caching checksums.

### Risk Mitigation
- **Performance bottlenecks**: Profile with `cProfile` and `line_profiler`; schedule optimization sprint if throughput < gzip baseline.
- **Binary compatibility**: Version header and JSON schema enforce forward/backward compatibility; add migration utilities before breaking changes.
- **External dependencies**: Pin critical libraries (`biopython`, `numpy`, `bitarray`) and test across OS/Python versions.
- **Large file handling**: Implement chunked streaming to avoid memory blow-ups; add stress tests with >10k sequences × 100k columns.

### Release & Distribution
- Package project for PyPI via `pyproject.toml` (set `project.scripts` for CLI).
- Draft CHANGELOG starting at v0.1.0-alpha aligned with Conventional Commits.
- Sign release artifacts, publish wheels/sdist, and attach benchmark summary for each release candidate.

## Milestone Acceptance Checklist
- ✅ Phase exit criteria met (unit/integration coverage, documentation updated).
- ✅ Issue tracker updated with tasks, dependencies, and owners.
- ✅ Benchmark results archived and compared against previous runs.
- ✅ Release artifacts generated and validated where applicable.
- ✅ Stakeholder review (code + docs) complete before advancing to next phase.

## Next Immediate Actions
1. Create initial `pyproject.toml`, `src/`, and `tests/` scaffolding (Phase 0 kick-off).
2. Add `pre-commit` configuration and seed smoke test to ensure CI viability.
3. Populate issue tracker with Phase 0/1 tasks referencing this plan for traceability.
