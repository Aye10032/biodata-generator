import random
from dataclasses import dataclass

from bioflow_sim.generators.random_values import weighted_index
from bioflow_sim.generators.sequences import DNA_ALPHABET


@dataclass(frozen=True)
class Snv:
    contig: str
    position: int
    reference: str
    alternate: str


def generate_snvs(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    count: int,
) -> list[Snv]:
    if count < 0:
        raise ValueError('SNV count cannot be negative')
    usable = [(name, sequence) for name, sequence in sequences if sequence]
    available = sum(sum(base in DNA_ALPHABET for base in sequence) for _, sequence in usable)
    if count > available:
        raise ValueError(f'requested {count} SNVs but only {available} reference bases are usable')

    variants: dict[tuple[str, int], Snv] = {}
    attempts = 0
    while len(variants) < count:
        attempts += 1
        if attempts > max(1000, count * 100):
            raise ValueError('could not place the requested number of unique SNVs')
        index = weighted_index(rng, [len(sequence) for _, sequence in usable])
        contig, sequence = usable[index]
        position = rng.randrange(len(sequence))
        reference = sequence[position]
        if reference not in DNA_ALPHABET or (contig, position) in variants:
            continue
        alternate = rng.choice(tuple(base for base in 'ACGT' if base != reference))
        variants[(contig, position)] = Snv(contig, position, reference, alternate)
    return sorted(variants.values(), key=lambda variant: (variant.contig, variant.position))


def apply_snvs(
    sequences: list[tuple[str, str]],
    variants: list[Snv],
) -> list[tuple[str, str]]:
    by_contig: dict[str, list[Snv]] = {}
    for variant in variants:
        by_contig.setdefault(variant.contig, []).append(variant)

    mutated = []
    for contig, sequence in sequences:
        bases = list(sequence)
        for variant in by_contig.get(contig, []):
            if bases[variant.position] != variant.reference:
                raise ValueError(f'{contig}:{variant.position + 1} reference allele mismatch')
            bases[variant.position] = variant.alternate
        mutated.append((contig, ''.join(bases)))
    return mutated
