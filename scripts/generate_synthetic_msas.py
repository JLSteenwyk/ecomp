#!/usr/bin/env python3
"""Generate synthetic multiple sequence alignments for benchmarking."""

from __future__ import annotations

import argparse
import random
from pathlib import Path
from typing import Iterable, Sequence

DNA_ALPHABET = "ACGT"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def random_consensus(length: int, alphabet: Sequence[str]) -> list[str]:
    return [random.choice(alphabet) for _ in range(length)]


def mutate_sequence(
    consensus: Sequence[str],
    mutation_rate: float,
    alphabet: Sequence[str],
) -> str:
    seq_chars: list[str] = []
    for char in consensus:
        if random.random() < mutation_rate:
            alternatives = [a for a in alphabet if a != char]
            seq_chars.append(random.choice(alternatives))
        else:
            seq_chars.append(char)
    return "".join(seq_chars)


def inject_motifs(
    consensus: list[str],
    motif: Sequence[str],
    repetitions: int,
) -> None:
    length = len(consensus)
    motif_length = len(motif)
    for _ in range(repetitions):
        if motif_length >= length:
            start = 0
        else:
            start = random.randint(0, length - motif_length)
        consensus[start : start + motif_length] = motif


def generate_alignment(
    num_taxa: int,
    length: int,
    alphabet: Sequence[str],
    mutation_rate: float,
    motif_probability: float,
) -> list[tuple[str, str]]:
    consensus = random_consensus(length, alphabet)

    if random.random() < motif_probability:
        motif_length = random.randint(8, 24)
        motif = [random.choice(alphabet) for _ in range(motif_length)]
        repetitions = max(1, length // (motif_length * 4))
        inject_motifs(consensus, motif, repetitions)

    sequences = []
    for idx in range(num_taxa):
        rate = mutation_rate * random.uniform(0.6, 1.4)
        seq = mutate_sequence(consensus, min(rate, 0.45), alphabet)
        sequences.append((f"taxon_{idx+1}", seq))
    return sequences


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate synthetic MSAs.")
    parser.add_argument(
        "output",
        type=Path,
        help="Directory to write FASTA alignments into",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
        help="Random seed (default: 1234)",
    )
    parser.add_argument(
        "--alphabet",
        choices=["dna", "protein"],
        default="protein",
        help="Alphabet for sequences (default: protein)",
    )
    parser.add_argument(
        "--mutation",
        type=float,
        default=0.08,
        help="Base mutation rate relative to consensus (default: 0.08)",
    )
    parser.add_argument(
        "--motif-probability",
        type=float,
        default=0.6,
        help="Probability of injecting conserved motifs (default: 0.6)",
    )
    parser.add_argument(
        "--sizes",
        nargs="*",
        default=[
            "10x200",
            "20x500",
            "20x1500",
            "50x500",
            "50x2000",
            "100x500",
            "100x2000",
            "200x1000",
            "200x5000",
            "500x1000",
        ],
        help="Alignment sizes formatted as '<taxa>x<length>'",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    random.seed(args.seed)

    alphabet = PROTEIN_ALPHABET if args.alphabet == "protein" else DNA_ALPHABET

    for size in args.sizes:
        try:
            taxa_str, length_str = size.lower().split("x")
            num_taxa = int(taxa_str)
            length = int(length_str)
        except ValueError as exc:  # pragma: no cover - CLI guard
            raise SystemExit(f"Invalid size specification: {size}") from exc

        records = generate_alignment(
            num_taxa=num_taxa,
            length=length,
            alphabet=alphabet,
            mutation_rate=args.mutation,
            motif_probability=args.motif_probability,
        )
        filename = (
            f"synthetic_{args.alphabet}_taxa{num_taxa}_len{length}_seed{args.seed}.fasta"
        )
        write_fasta(args.output / filename, records)

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
