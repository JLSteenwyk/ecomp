Architecture
============

Data Flow
---------

#. **Parsing** – ``io.read_alignment`` loads FASTA or PHYLIP files into an
   :class:`~evolutionary_compression.io.AlignmentFrame` while normalising
   metadata such as sequence identifiers and alphabet.
#. **Profiling** – ``compression.consensus.iter_column_profiles`` produces
   deterministic consensus calls and deviation bitmasks for each alignment
   column.
#. **Run-length encoding** – ``compression.rle.collect_run_length_blocks``
   groups consecutive columns that share consensus/deviation signatures.
#. **Binary encoding** – ``compression.encoding.encode_blocks`` serialises each
   block into a compact payload (run length, consensus byte, deviation mask,
   bit-packed residues).
#. **Entropy stage** – ``zlib`` compresses the raw payload opportunistically;
   whichever representation is smaller is persisted and the choice recorded in
   metadata.
#. **Metadata capture** – ``compression.pipeline.compress_alignment`` stores
   structural metadata (sequence count, alignment length, checksum) alongside
   the payload, written by ``storage.write_payload`` / ``write_metadata``.
#. **Decompression** – ``compression.pipeline.decompress_alignment`` reverses
   the process, reconstructs sequences, and verifies optional checksums.

Extensibility Points
--------------------

- **Entropy coding** – add alternative compressors by extending
  ``compression.encoding`` with pluggable post-RLE coders.
- **Streaming** – adapt ``iter_column_profiles`` to stream columns lazily for
  very large alignments.
- **Metadata evolution** – JSON metadata is versioned; additive fields are
  backwards-compatible while breaking changes require a migration path.

Performance Considerations
--------------------------

- Alignment columns are iterated with ``zip(*sequences, strict=True)`` to avoid
  per-character slicing.
- Run-length blocks are capped at 65,535 columns to keep storage predictable.
- The current implementation is pure Python for clarity; hotspots are isolated
  so that future NumPy or Cython accelerations can be dropped in.

Validation Strategy
-------------------

- Unit tests assert consensus selection, RLE packing, frequency metadata, and
  binary encode/decode symmetry.
- Integration tests exercises CLI round-trips and checksum validation on
  synthetic alignments.
- ``scripts/compare_compressors.py`` provides repeatable benchmarks against
  ``gzip``, ``bzip2``, and ``xz``.
