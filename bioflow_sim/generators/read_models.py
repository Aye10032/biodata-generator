import random
from dataclasses import dataclass

from .sequences import DNA_ALPHABET


@dataclass(frozen=True)
class Technology:
    name: str
    paired: bool
    mean_length: int
    length_sd: int
    substitution_rate: float
    insertion_rate: float
    deletion_rate: float
    quality_char: str


TECHNOLOGIES: dict[str, Technology] = {
    'illumina-se': Technology('illumina-se', False, 150, 0, 0.001, 0.0, 0.0, 'I'),
    'illumina-pe': Technology('illumina-pe', True, 150, 0, 0.001, 0.0, 0.0, 'I'),
    'pacbio-clr': Technology('pacbio-clr', False, 12_000, 6_000, 0.035, 0.055, 0.055, '5'),
    'pacbio-hifi': Technology('pacbio-hifi', False, 15_000, 3_000, 0.0007, 0.0006, 0.0007, '?'),
    'ont': Technology('ont', False, 15_000, 8_000, 0.025, 0.012, 0.018, '8'),
}


def introduce_errors(
    sequence: str,
    rng: random.Random,
    technology: Technology,
) -> tuple[str, int, int, int]:
    output: list[str] = []
    substitutions = insertions = deletions = 0
    bases = 'ACGT'

    for base in sequence:
        if base not in DNA_ALPHABET:
            base = rng.choice(bases)
        if rng.random() < technology.deletion_rate:
            deletions += 1
            continue
        if rng.random() < technology.insertion_rate:
            output.append(rng.choice(bases))
            insertions += 1
        if rng.random() < technology.substitution_rate:
            output.append(rng.choice(tuple(b for b in bases if b != base)))
            substitutions += 1
        else:
            output.append(base)
    return ''.join(output), substitutions, insertions, deletions
