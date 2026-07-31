import random
from dataclasses import dataclass

from bioflow_sim.generators.sequences import GenomicTemplate, reverse_complement


@dataclass(frozen=True)
class MethylationSite:
    contig: str
    position: int
    methylated: bool


def generate_cpg_methylation(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    count: int,
    methylation_rate: float,
) -> list[MethylationSite]:
    if count < 1 or not 0 <= methylation_rate <= 1:
        raise ValueError('site count must be positive and methylation rate must be between 0 and 1')
    candidates = [
        (contig, position)
        for contig, sequence in sequences
        for position in range(len(sequence) - 1)
        if sequence[position : position + 2] == 'CG'
    ]
    if count > len(candidates):
        raise ValueError(f'requested {count} CpG sites but reference contains {len(candidates)}')
    selected = rng.sample(candidates, count)
    sites = [MethylationSite(contig, position, rng.random() < methylation_rate) for contig, position in selected]
    return sorted(sites, key=lambda site: (site.contig, site.position))


def sample_template_around_cpg(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    sites: list[MethylationSite],
    length: int,
) -> tuple[GenomicTemplate, MethylationSite]:
    site = rng.choice(sites)
    source = dict(sequences)[site.contig]
    lower = max(0, site.position - length + 1)
    upper = min(site.position, len(source) - length)
    if upper < lower:
        raise ValueError(f'{site.contig}:{site.position + 1} cannot support a {length} bp template')
    start = rng.randint(lower, upper)
    sequence = source[start : start + length]
    if rng.random() < 0.5:
        sequence = reverse_complement(sequence)
        strand = '-'
    else:
        strand = '+'
    return GenomicTemplate(site.contig, start, start + length, strand, sequence), site


def convert_unmethylated_cytosines(
    rng: random.Random,
    template: GenomicTemplate,
    sites: dict[tuple[str, int], MethylationSite],
    conversion_rate: float,
) -> tuple[str, int]:
    if not 0 <= conversion_rate <= 1:
        raise ValueError('conversion rate must be between 0 and 1')
    converted = list(template.sequence)
    conversions = 0
    for index, base in enumerate(template.sequence):
        if template.strand == '+':
            position = template.start + index
            site = sites.get((template.source_id, position))
            convertible = base == 'C'
            replacement = 'T'
        else:
            position = template.end - index - 1
            site = sites.get((template.source_id, position))
            convertible = base == 'G'
            replacement = 'A'
        if site is not None and not site.methylated and convertible and rng.random() < conversion_rate:
            converted[index] = replacement
            conversions += 1
    return ''.join(converted), conversions
