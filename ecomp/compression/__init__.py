"""Compression-related helpers exposed for advanced integrations."""

from .consensus import ColumnProfile, collect_column_profiles, iter_column_profiles
from .encoding import decode_blocks, encode_blocks
from .pipeline import CompressedAlignment, compress_alignment, decompress_alignment
from .rle import RunLengthBlock, collect_run_length_blocks, iter_run_length_blocks

__all__ = [
    "ColumnProfile",
    "collect_column_profiles",
    "iter_column_profiles",
    "decode_blocks",
    "encode_blocks",
    "CompressedAlignment",
    "compress_alignment",
    "decompress_alignment",
    "RunLengthBlock",
    "collect_run_length_blocks",
    "iter_run_length_blocks",
]
