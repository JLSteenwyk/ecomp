# Alignment Fixtures

Place lightweight FASTA/PHYLIP alignments here for unit and integration tests.

Suggested pattern:
- `small_orthogroup.fasta` – tens of sequences, short length for smoke tests.
- `small_phylo.fasta` + `small_phylo.tree` – paired alignment + Newick tree for phylo bundle integration tests.
- `protein_family.phylip` – moderate size (hundreds of columns) to exercise gap handling.
- `ambiguous_cases.fasta` – include gaps (`-`) and ambiguous symbols (`N`, `X`).

Document provenance and licensing for any externally sourced data in this folder.
Checksum files with `shasum -a 256 <file>` and record them in `FIXTURES.md` when added.
