# Evolutionary MSA Compression – Developer Instructions

## Overview
We want to build a **Python-based compressor/decompressor** for multiple sequence alignments (MSAs).  
The focus is on **lossless compression** that safely reduces space usage, even across **deep evolutionary timescales** where signals are weak.  

We will start with a **Column-wise Run-Length Encoding (RLE) + Consensus Anchoring** approach.  
This balances compression gains, speed, and safety.

---

## Goals
1. **Input**: FASTA/PHYLIP MSA.  
2. **Output**: Custom `.ecomp` compressed binary format + metadata (JSON).  
3. **Compression Strategy**:
   - Identify consensus per column.  
   - Store consensus + deltas (differences from consensus).  
   - Apply run-length encoding (RLE) for long stretches of identical states.  
4. **Decompression**: Must fully restore original alignment (lossless).  
5. **Validation**: Ensure compressed → decompressed alignment is byte-for-byte identical.  

---

## Step 1. Input/Output Handling
- Use **Biopython** (`Bio.AlignIO`) to parse MSAs.  
- Output:
  - `.ecomp` binary file containing compressed payload.  
  - `.json` sidecar file with metadata (number of sequences, alignment length, alphabet, method version).  

---

## Step 2. Compression Algorithm
**a. Consensus Anchoring**
- For each column, determine the most frequent residue (consensus).  
- Store the consensus as a single character.  
- Encode deviations (taxon index + residue).  

**b. Column-Wise Run-Length Encoding**
- If consecutive columns share the same consensus character, run-length encode them.  
- Store as: `(consensus_char, run_length, deviations_list[])`.  

**c. Binary Encoding**
- Use `bitarray` or `numpy.packbits` to efficiently store deltas.  

---

## Step 3. Decompression Algorithm
- Read metadata + compressed blocks.  
- Expand RLE runs into consensus columns.  
- Apply stored deviations to reconstruct the full alignment.  
- Output FASTA or PHYLIP identical to input.  

---

## Step 4. Validation
- Compute SHA256 checksum of input alignment.  
- After round-trip (compress → decompress), recompute checksum.  
- They must match exactly.  

---

## Step 5. Benchmarking
- Compare against `gzip`, `bzip2`, `xz`.  
- Metrics:
  - Compression ratio = original_size / compressed_size.  
  - Compression/decompression time.  
- Test datasets:
  - Small orthogroup (20 sequences × 1k bp).  
  - Large phylogenomic MSA (10k sequences × 100k bp).  

---

## Step 6. Extensions (Optional / Future Work)
- Add entropy-based lossy filtering (optional mode).  
- Add streaming decompression (yield sites without full decompression).  
- Standardize `.ecomp` file format for archival.  

---

## Suggested Python Libraries
- **Biopython** → MSA parsing.  
- **NumPy** → efficient array ops.  
- **bitarray** or `numpy.packbits` → binary packing.  
- **hashlib** → SHA256 validation.  
- **multiprocessing** → optional parallel column processing.  

---

## Minimal Prototype Pseudocode

```python
from Bio import AlignIO
import numpy as np

def compress_msa(input_fasta, output_file):
    alignment = AlignIO.read(input_fasta, "fasta")
    n_seq = len(alignment)
    aln_len = alignment.get_alignment_length()
    
    compressed = []
    for col_idx in range(aln_len):
        column = [record.seq[col_idx] for record in alignment]
        consensus = max(set(column), key=column.count)
        deviations = [(i, c) for i, c in enumerate(column) if c != consensus]
        compressed.append((consensus, deviations))
    
    # TODO: Apply RLE here
    # TODO: Write to binary + metadata
    
    return compressed

def decompress_msa(compressed, n_seq, aln_len):
    # TODO: Expand RLE, apply deviations
    pass
