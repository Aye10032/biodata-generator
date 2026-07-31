import random
from dataclasses import dataclass

from .random_values import weighted_index

DNA_ALPHABET = frozenset('ACGT')


@dataclass(frozen=True)
class GenomicTemplate:
    source_id: str
    start: int
    end: int
    strand: str
    sequence: str


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def sample_genomic_template(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    length: int,
) -> GenomicTemplate:
    eligible = [
        (identifier, sequence)
        for identifier, sequence in sequences
        if len(sequence) >= length
    ]
    if not eligible:
        raise ValueError(f'no reference sequence is long enough for a {length} bp template')
    index = weighted_index(rng, [len(sequence) - length + 1 for _, sequence in eligible])
    identifier, source = eligible[index]
    start = rng.randrange(0, len(source) - length + 1)
    template = source[start : start + length]
    if rng.random() < 0.5:
        template = reverse_complement(template)
        strand = '-'
    else:
        strand = '+'
    return GenomicTemplate(identifier, start, start + length, strand, template)
