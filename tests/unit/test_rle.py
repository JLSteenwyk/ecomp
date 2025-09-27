import math

from evolutionary_compression.io import alignment_from_sequences
from evolutionary_compression.compression.consensus import collect_column_profiles
from evolutionary_compression.compression.rle import collect_run_length_blocks


def test_run_length_blocks_group_identical_columns():
    frame = alignment_from_sequences(
        ids=["seq1", "seq2"],
        sequences=["AAAAA", "AAAAT"],
    )
    columns = collect_column_profiles(frame)
    alphabet = frame.alphabet
    symbol_lookup = {symbol: index for index, symbol in enumerate(alphabet)}
    bits_per_symbol = max(1, math.ceil(math.log2(max(len(alphabet), 1))))
    blocks = collect_run_length_blocks(
        columns, frame.num_sequences, symbol_lookup, bits_per_symbol
    )

    assert len(blocks) == 2
    first_block, second_block = blocks
    assert first_block.run_length == 4
    assert int.from_bytes(first_block.bitmask, "little") == 0
    assert first_block.residues == b""
    assert second_block.run_length == 1
    assert int.from_bytes(second_block.bitmask, "little") == 0b10
    assert second_block.residues == b"\x80"
