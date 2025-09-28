# eComp Documentation Notes

## Current Live Docs
- `docs/index.rst`: Table of contents (matches PhyKIT layout).
- `docs/usage/quickstart.rst`: basic commands and API pointer.
- `docs/usage/api.rst`: Python API summary.
- `docs/usage/benchmarking.rst`: instructions + tree-ordering script mention and `ECOMP_SEQUENCE_ORDER` env var.
- `docs/tutorials/workflows.rst`: CLI examples (compress dir, benchmarking, tree sweep, Python snippet).
- `docs/frequently_asked_questions/index.rst`: covers installation, tree option, fallback.
- `docs/change_log/index.rst`: “Unreleased” entry summarising removal of bundle, new heuristics, CLI alias.
- `docs/other_software/index.rst`: pointers to PhyKIT, IQ-TREE, etc.
- `docs/data/index.rst`: summary of packaged datasets.

## Open Follow-ups
- Fill the placeholder sections (`Change Log` with future releases, `Workflows` with richer tutorials).
- Document benchmark results in a dedicated page (e.g., summarize findings from synthetic and bee/bird JSON outputs).
- Confirm final architecture/API docs (currently removed old markdown `architecture.rst` etc.).

## Command alias
- `ecomp` and `ec` both point to the CLI entry point (pyproject updated).

