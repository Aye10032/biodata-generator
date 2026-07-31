import random
import re
from dataclasses import dataclass

from bioflow_sim.generators.random_values import random_dna, weighted_index

GENE_PATTERNS = (
    re.compile(r'(?:gene_id|gene)[=:]([^\s\]]+)'),
    re.compile(r'\[gene=([^\]]+)\]'),
)


@dataclass(frozen=True)
class Transcript:
    identifier: str
    feature_id: str
    sequence: str


@dataclass(frozen=True)
class TranscriptFragment:
    transcript: Transcript
    start: int
    sequence: str


def feature_id_from_header(identifier: str, header: str) -> str:
    for pattern in GENE_PATTERNS:
        match = pattern.search(header)
        if match:
            return match.group(1).strip('"')
    return identifier


def build_transcripts(
    fasta_records: list[tuple[str, str, str]],
    annotation_mapping: dict[str, str],
) -> list[Transcript]:
    transcripts = []
    for identifier, header, sequence in fasta_records:
        feature_id = annotation_mapping.get(identifier, feature_id_from_header(identifier, header))
        transcripts.append(Transcript(identifier, feature_id, sequence))
    return sorted(transcripts, key=lambda transcript: transcript.identifier)


def unique_barcodes(
    rng: random.Random,
    count: int,
    length: int = 16,
) -> list[str]:
    barcodes: set[str] = set()
    while len(barcodes) < count:
        barcodes.add(random_dna(rng, length))
    return sorted(barcodes)


def choose_transcript(
    rng: random.Random,
    transcripts: list[Transcript],
    weights: list[float],
    minimum_length: int,
) -> Transcript:
    eligible = [
        (transcript, weight)
        for transcript, weight in zip(transcripts, weights, strict=True)
        if len(transcript.sequence) >= minimum_length
    ]
    if not eligible:
        raise ValueError(f'no transcript is at least {minimum_length} bp')
    index = weighted_index(rng, [weight for _, weight in eligible])
    return eligible[index][0]


def cell_expression_weights(base: list[float], cell_type: int) -> list[float]:
    """Create one of two deterministic synthetic expression states."""
    state = base.copy()
    block = max(1, len(state) // 10)
    start = 0 if cell_type == 1 else block
    for index in range(start, min(start + block, len(state))):
        state[index] *= 5.0
    return state


def differential_expression_weights(
    base: list[float],
    group: str,
    fold_change: float,
) -> list[float]:
    """Boost disjoint transcript blocks in control and treatment groups."""
    if fold_change <= 0:
        raise ValueError('fold change must be positive')
    state = base.copy()
    block = max(1, len(state) // 10)
    start = 0 if group == 'control' else block
    for index in range(start, min(start + block, len(state))):
        state[index] *= fold_change
    return state


def sample_three_prime_fragment(
    rng: random.Random,
    transcript: Transcript,
    read_length: int,
    window: int = 500,
) -> TranscriptFragment:
    start_min = max(0, len(transcript.sequence) - window)
    start_max = len(transcript.sequence) - read_length
    start = rng.randint(start_min, start_max)
    sequence = transcript.sequence[start : start + read_length]
    return TranscriptFragment(transcript, start, sequence)


def sample_transcript_fragment(
    rng: random.Random,
    transcript: Transcript,
    minimum_length: int,
    mean: int,
    sd: int,
) -> TranscriptFragment:
    fragment_length = min(
        len(transcript.sequence),
        max(minimum_length, round(rng.gauss(mean, sd))),
    )
    start = rng.randrange(0, len(transcript.sequence) - fragment_length + 1)
    sequence = transcript.sequence[start : start + fragment_length]
    return TranscriptFragment(transcript, start, sequence)
