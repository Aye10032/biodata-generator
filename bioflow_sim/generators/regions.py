import random
from dataclasses import dataclass

from .random_values import weighted_index
from .sequences import GenomicTemplate, reverse_complement


@dataclass(frozen=True)
class GenomicInterval:
    contig: str
    start: int
    end: int
    name: str
    score: int


def generate_intervals(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    count: int,
    width: int,
    prefix: str = 'peak',
) -> list[GenomicInterval]:
    if count < 1 or width < 1:
        raise ValueError('interval count and width must be positive')
    eligible = [(name, sequence) for name, sequence in sequences if len(sequence) >= width]
    intervals: list[GenomicInterval] = []
    occupied: dict[str, list[tuple[int, int]]] = {}
    attempts = 0
    while len(intervals) < count:
        attempts += 1
        if attempts > count * 1000:
            raise ValueError('could not place non-overlapping intervals')
        index = weighted_index(rng, [len(sequence) - width + 1 for _, sequence in eligible])
        contig, sequence = eligible[index]
        start = rng.randrange(0, len(sequence) - width + 1)
        end = start + width
        if any(start < old_end and end > old_start for old_start, old_end in occupied.get(contig, [])):
            continue
        occupied.setdefault(contig, []).append((start, end))
        intervals.append(
            GenomicInterval(
                contig,
                start,
                end,
                f'{prefix}_{len(intervals) + 1:05d}',
                rng.randint(100, 1000),
            )
        )
    return sorted(intervals, key=lambda interval: (interval.contig, interval.start))


def sample_template_near_interval(
    rng: random.Random,
    sequences: list[tuple[str, str]],
    intervals: list[GenomicInterval],
    length: int,
) -> tuple[GenomicTemplate, GenomicInterval]:
    interval = rng.choice(intervals)
    source = dict(sequences)[interval.contig]
    lower = max(0, interval.start - length + 1)
    upper = min(len(source) - length, interval.end - 1)
    if upper < lower:
        raise ValueError(f'{interval.name}: interval cannot support a {length} bp fragment')
    start = rng.randint(lower, upper)
    sequence = source[start : start + length]
    if rng.random() < 0.5:
        sequence = reverse_complement(sequence)
        strand = '-'
    else:
        strand = '+'
    return (
        GenomicTemplate(interval.contig, start, start + length, strand, sequence),
        interval,
    )
