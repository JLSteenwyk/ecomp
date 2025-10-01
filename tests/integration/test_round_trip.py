from pathlib import Path

from ecomp import compress_file, decompress_file, read_alignment


FASTA_CONTENT = ">seq1\nACGTACGT\n>seq2\nACGTTCGT\n>seq3\nACGTACGA\n"


def test_compress_and_decompress_round_trip(tmp_path: Path):
    input_path = tmp_path / "example.fasta"
    input_path.write_text(FASTA_CONTENT)

    ecomp_path, metadata_path = compress_file(input_path)
    assert ecomp_path.exists()
    assert metadata_path.exists()

    output_path = decompress_file(ecomp_path, metadata_path=metadata_path)
    assert output_path.exists()

    original = read_alignment(input_path)
    restored = read_alignment(output_path)

    assert original.sequences == restored.sequences
    assert original.ids == restored.ids
