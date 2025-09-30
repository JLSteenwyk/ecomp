# eComp Format README

This document explains the `.ecomp` container so another project can parse or
emit archives compatible with the Evolutionary Compression (eComp) tool.

## Artifact Overview
- Each archive ships as two files:
  - `<name>.ecomp`: binary payload prefixed with an 8-byte header.
  - `<name>.json`: metadata in JSON (optionally zlib-compressed).
- The metadata file stores alignment dimensions, codec choices, and other hints
  required to rebuild the original multiple sequence alignment.

## `.ecomp` Container Header
The file begins with a fixed-size header defined by the struct format `>8sBBBQ`:
```
+------------+------+-------+
| Field      | Size | Notes |
+------------+------+-------+
| magic      | 8    | ASCII `ECOMP001`
| version    | 3×1  | major, minor, patch bytes (e.g., 0,1,0)
| payloadLen | 8    | big-endian uint64 length of the payload
+------------+------+-------+
```
Immediately following the header is `payloadLen` bytes of encoded data.

## Metadata Sidecar
- Stored alongside the payload (`<name>.json`).
- Default suffix: `.json` (`METADATA_SUFFIX`).
- Content is UTF-8 JSON with sorted keys. When compression helps, the file is
  wrapped as `b"ECMZ" + version byte + zlib(data)`.
- Core keys emitted by the compressor:
  - `format_version`: semantic version string.
  - `codec`: always `"ecomp"` for native payloads.
  - `num_sequences`, `alignment_length`.
  - `alphabet`: list of characters used in the alignment.
  - `source_format`: original alignment format (e.g., `fasta`).
  - `checksum_sha256`: digest of the original sequences (order-sensitive).
  - `run_length_blocks`, `max_run_length`, `columns_with_deviations`.
  - `bitmask_bytes`: bytes needed for per-column deviation bitmasks.
  - `bits_per_symbol`: width used when packing residues.
  - `payload_encoding`: one of `raw`, `zlib`, `zstd`, `xz`, or `gzip` (fallback).
  - `payload_encoded_bytes`, `payload_raw_bytes`.
  - `sequence_id_codec`: currently `"inline"` to signal the ID block lives in
    the payload.
  - `ordering_strategy`: `baseline`, `tree`, `auto-*`, etc.
  - Optional `sequence_permutation` and `sequence_ids` entries (see below).
  - Optional `fallback` map when gzip was chosen (`{"type": "gzip", ...}`).

## Decoding Workflow
1. Load metadata (decompress if it begins with `ECMZ`).
2. Read the `.ecomp` header and slice the payload bytes.
3. If metadata contains `fallback.type == "gzip"`, treat the payload as a gzip
   stream of FASTA bytes, ignoring the native structures.
4. Otherwise decode according to `payload_encoding`:
   - `raw`: leave as-is.
   - `zlib`: `zlib.decompress`.
   - `zstd`: requires the `zstandard` Python module if you follow our code.
   - `xz`: `lzma.decompress`.
5. Inspect the decoded payload for optional chunks and the run-length stream in
   this order:
   1. **Permutation chunk** (optional) — only present when metadata contains
      `sequence_permutation` with `{"encoding": "payload"}`. Structure:
      - Magic `ECPE` (4 bytes) + version byte (`0x01`).
      - Flag byte: bit0 indicates zlib compression of the following payload,
        bits1-2 encode index width (`0:uint8`, `1:uint16`, `2:uint32`).
      - Varint number of indices followed by varint payload length.
      - Payload: contiguous index bytes (possibly zlib-compressed).
      After removal, update the metadata permutation to the decoded list and
      continue with the remainder of the payload.
   2. **Sequence ID chunk** — always present and identified by magic `ECID`.
      - Version 2 layout: `ECID` + `0x02` + varint block length + mode byte.
      - Mode `0`: raw bytes; `1`: zstd-compressed; `2`: zlib-compressed.
      - Decoded block: varint sequence count followed by length-prefixed UTF-8
        identifiers (each length is a varint). Metadata `sequence_ids` may be
        absent; if so, populate it with the decoded list.
   3. **Run-length stream** — the remaining bytes encode consensus blocks.

## Run-Length Stream Structure
The stream is decoded according to `encoding.py`:
1. **Consensus models**
   - Leading byte: number of distinct consensus residues with deviations.
   - For each entry:
     - Consensus character (1 byte ASCII).
     - Mode (1 byte): `0` = fixed-width packing, `1` = Huffman.
     - Residue count (1 byte) + that many residue characters (ASCII bytes).
     - Mode-specific data:
       - Mode 0: single byte width (bits per local symbol).
       - Mode 1: residue count bytes of Huffman code lengths (canonical order).
2. **Dictionary section**
   - Byte: dictionary size (`<=255`).
   - Each entry stores a reusable column pattern:
     - Consensus (1 byte) and bitmask mode (1 byte).
     - Deviation count (varint) and bitmask payload length (varint).
     - Bitmask payload bytes (encoding variant described below).
     - Residue payload length (uint16 big-endian) and payload bytes.
3. **Block records**
   - 4-byte big-endian uint32: number of blocks.
   - For each block:
     - Marker byte: `1` = reference dictionary entry, `0` = literal.
     - If marker `1`: dictionary id (1 byte) + run length (1 byte).
     - If marker `0`:
       - Run length (1 byte) and consensus (1 byte).
       - Bitmask mode (1 byte), deviation count (varint), mask length (varint),
         mask payload bytes.
       - Residue payload length (uint16) and residue bytes.

### Bitmask Encoding Modes
Bitmasks flag which sequences deviate from the consensus in each column. Three
encodings are used; the compressor always chooses the smallest:
- Mode `0`: raw bytes (tail trimmed to the last non-zero byte).
- Mode `1`: varint-coded delta positions. The payload begins with a varint
  length followed by that many varints representing positional deltas.
- Mode `2`: byte-wise run-length encoding (pairs of byte value + repeat count).

### Residue Payloads
Residual symbols are stored relative to the consensus:
- For mode `0` consensus models, residues map to a small alphabet with fixed
  width (1–8 bits). The payload packs those codes tightly.
- For mode `1`, residues are encoded with canonical Huffman codes whose lengths
  appear in the consensus table.
- The decoder expands the local symbols and then maps them back to the global
  alphabet indices provided in metadata.

### Varint Convention
All variable-length integers use a little-endian base-128 scheme (7 data bits,
continuation bit set on intermediate bytes). This matches the helpers in
`encoding.py` (`_write_varint`, `_read_varint`).

## Reconstructing the Alignment
1. Iterate blocks in order. Each block represents `run_length` consecutive
   columns sharing one consensus character and an identical deviation pattern.
2. For every column:
   - Start with an array filled with the consensus character.
   - Use the decoded bitmask to locate sequences with deviations.
   - Inject the decoded residues (already in global alphabet order) at those
     positions.
3. Append each column to per-sequence buffers. After all blocks, verify that the
   number of generated columns equals `alignment_length`.
4. If `sequence_permutation` is a list, invert it to restore the original
   sequence order. The permutation lists the source index for each reordered
   sequence; the inverse permutation puts sequences back to their original
   positions.
5. Validate the SHA-256 checksum against `checksum_sha256` for safety.

## Fallback Mode
The compressor may prefer gzip for highly random alignments. When metadata
contains `{"fallback": {"type": "gzip", "format": "fasta"}}`, the payload
bytes are simply `gzip.compress(fasta_bytes)`. Metadata fields such as
`sequence_permutation` and consensus statistics are omitted in this branch.

## Practical Tips
- Treat all reserved strings (`ECOMP001`, `ECMZ`, `ECPE`, `ECID`) as
  case-sensitive ASCII.
- The format currently assumes single-byte residue characters (ASCII). If you
  encounter non-standard alphabets, the compressor will raise an error.
- Zstandard support is optional; if your environment lacks it and you see
  `payload_encoding == "zstd"`, you must link against a zstd decoder.
- `sequence_ids` may appear both in metadata (explicit list) and in the payload
  (inline block). The post-decompression metadata should reflect the definitive
  ID order.
- New metadata keys may appear in future versions. Preserve unknown fields when
  forwarding metadata.

Following this README, another agent can implement a bit-for-bit compatible
encoder or decoder for `.ecomp` archives.
