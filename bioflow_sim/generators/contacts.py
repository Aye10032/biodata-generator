import bisect
import random
from dataclasses import dataclass

from bioflow_sim.generators.sequences import reverse_complement


@dataclass(frozen=True)
class RestrictionEnzyme:
    name: str
    motif: str
    cut_offset: int


@dataclass(frozen=True)
class CutSite:
    contig: str
    position: int


@dataclass(frozen=True)
class Contact:
    left: CutSite
    right: CutSite
    contact_type: str


ENZYMES: dict[str, RestrictionEnzyme] = {
    'mboi': RestrictionEnzyme('MboI', 'GATC', 0),
    'dpnii': RestrictionEnzyme('DpnII', 'GATC', 0),
    'hindiii': RestrictionEnzyme('HindIII', 'AAGCTT', 1),
}


def find_restriction_sites(
    sequences: list[tuple[str, str]],
    enzyme: RestrictionEnzyme,
) -> list[CutSite]:
    sites: list[CutSite] = []
    for contig, sequence in sequences:
        start = 0
        while True:
            match = sequence.find(enzyme.motif, start)
            if match < 0:
                break
            sites.append(CutSite(contig, match + enzyme.cut_offset))
            start = match + 1
    if len(sites) < 2:
        raise ValueError(f'reference contains fewer than two {enzyme.name} cut sites')
    return sites


def sample_contact(
    rng: random.Random,
    sites: list[CutSite],
    intra_chromosomal_rate: float,
    mean_distance: int,
) -> Contact:
    if not 0 <= intra_chromosomal_rate <= 1:
        raise ValueError('intra-chromosomal rate must be between 0 and 1')
    by_contig: dict[str, list[CutSite]] = {}
    for site in sites:
        by_contig.setdefault(site.contig, []).append(site)
    left = rng.choice(sites)

    if rng.random() < intra_chromosomal_rate or len(by_contig) == 1:
        same = by_contig[left.contig]
        positions = [site.position for site in same]
        direction = -1 if rng.random() < 0.5 else 1
        target = max(0, left.position + direction * max(1, round(rng.expovariate(1 / mean_distance))))
        index = bisect.bisect_left(positions, target)
        index = min(index, len(same) - 1)
        right = same[index]
        if right == left and len(same) > 1:
            right = same[(index + 1) % len(same)]
        contact_type = 'intra'
    else:
        other_contigs = [contig for contig in by_contig if contig != left.contig]
        right = rng.choice(by_contig[rng.choice(other_contigs)])
        contact_type = 'inter'
    return Contact(left, right, contact_type)


def contact_read_sequences(
    sequences: list[tuple[str, str]],
    contact: Contact,
    read_length: int,
) -> tuple[str, str]:
    references = dict(sequences)
    left_source = references[contact.left.contig]
    right_source = references[contact.right.contig]
    left_start = min(max(0, contact.left.position), len(left_source) - read_length)
    right_end = min(max(read_length, contact.right.position), len(right_source))
    read1 = left_source[left_start : left_start + read_length]
    read2 = reverse_complement(right_source[right_end - read_length : right_end])
    return read1, read2
