import random
from dataclasses import dataclass

from bioflow_sim.generators.random_values import weighted_index
from bioflow_sim.generators.sequences import DNA_ALPHABET


@dataclass(frozen=True)
class TargetRegion:
    contig: str
    start: int
    end: int
    name: str


@dataclass(frozen=True)
class TumorNormalVariant:
    identifier: str
    contig: str
    position: int
    reference: str
    alternate: str
    category: str
    haplotype: int
    cancer_cell_fraction: float

    def normal_vaf(self) -> float:
        return 0.5 if self.category == 'germline' else 0.0

    def tumor_vaf(self, purity: float) -> float:
        if self.category == 'germline':
            return 0.5
        return purity * self.cancer_cell_fraction / 2


@dataclass(frozen=True)
class TargetFragment:
    target: TargetRegion
    start: int
    end: int
    sequence: str


def select_target_regions(
    rng: random.Random,
    bed_records: list[tuple[str, int, int, str]],
    sequences: list[tuple[str, str]],
    count: int,
    width: int,
) -> list[TargetRegion]:
    if count < 1 or width < 1:
        raise ValueError('target count and width must be positive')
    sequence_lengths = {contig: len(sequence) for contig, sequence in sequences}
    candidates = list(bed_records)
    rng.shuffle(candidates)
    selected: list[TargetRegion] = []
    for contig, bed_start, bed_end, name in candidates:
        contig_length = sequence_lengths.get(contig)
        if contig_length is None or bed_end - bed_start < width:
            continue
        midpoint = (bed_start + bed_end) // 2
        start = min(max(0, midpoint - width // 2), contig_length - width)
        end = start + width
        if any(region.contig == contig and start < region.end and end > region.start for region in selected):
            continue
        selected.append(TargetRegion(contig, start, end, f'TARGET_{len(selected) + 1:03d}_{name}'))
        if len(selected) == count:
            return sorted(selected, key=lambda region: (region.contig, region.start))
    raise ValueError(f'could only select {len(selected)} non-overlapping target regions, requested {count}')


def generate_tumor_normal_variants(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    targets: list[TargetRegion],
    germline_count: int,
    clonal_count: int,
    subclonal_count: int,
    subclone_fraction: float,
) -> list[TumorNormalVariant]:
    if min(germline_count, clonal_count, subclonal_count) < 0:
        raise ValueError('variant counts cannot be negative')
    if not 0 <= subclone_fraction <= 1:
        raise ValueError('subclone fraction must be between 0 and 1')
    references = dict(sequences)
    positions = [
        (target.contig, position)
        for target in targets
        for position in range(target.start, target.end)
        if references[target.contig][position] in DNA_ALPHABET
    ]
    total = germline_count + clonal_count + subclonal_count
    if total > len(positions):
        raise ValueError(f'requested {total} variants but targets contain only {len(positions)} usable bases')
    chosen = rng.sample(positions, total)
    specifications = (
        [('germline', 1.0)] * germline_count
        + [('somatic_clonal', 1.0)] * clonal_count
        + [('somatic_subclonal', subclone_fraction)] * subclonal_count
    )
    variants = []
    for index, ((contig, position), (category, fraction)) in enumerate(
        zip(chosen, specifications, strict=True),
        1,
    ):
        reference = references[contig][position]
        alternate = rng.choice(tuple(base for base in 'ACGT' if base != reference))
        variants.append(
            TumorNormalVariant(
                f'VAR{index:05d}',
                contig,
                position,
                reference,
                alternate,
                category,
                rng.randrange(2),
                fraction,
            )
        )
    return sorted(variants, key=lambda variant: (variant.contig, variant.position))


def sample_target_fragment(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    targets: list[TargetRegion],
    length: int,
) -> TargetFragment:
    eligible = [target for target in targets if target.end - target.start >= length]
    if not eligible:
        raise ValueError(f'no target region can support a {length} bp fragment')
    index = weighted_index(rng, [target.end - target.start - length + 1 for target in eligible])
    target = eligible[index]
    start = rng.randrange(target.start, target.end - length + 1)
    end = start + length
    sequence = dict(sequences)[target.contig][start:end]
    return TargetFragment(target, start, end, sequence)


def apply_variants_to_fragment(
    fragment: TargetFragment,
    variants: list[TumorNormalVariant],
    *,
    haplotype: int,
    cellular_source: str,
    subclone: bool,
) -> tuple[str, list[str]]:
    sequence = list(fragment.sequence)
    applied: list[str] = []
    for variant in variants:
        if (
            variant.contig != fragment.target.contig
            or not fragment.start <= variant.position < fragment.end
            or variant.haplotype != haplotype
        ):
            continue
        present = variant.category == 'germline' or (
            cellular_source == 'tumor'
            and (variant.category == 'somatic_clonal' or (variant.category == 'somatic_subclonal' and subclone))
        )
        if present:
            sequence[variant.position - fragment.start] = variant.alternate
            applied.append(variant.identifier)
    return ''.join(sequence), applied
