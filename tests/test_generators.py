import random

from bioflow_sim.generators.expression import (
    Transcript,
    build_transcripts,
    cell_expression_weights,
    choose_transcript,
    sample_three_prime_fragment,
    unique_barcodes,
)
from bioflow_sim.generators.random_values import sample_positive_normal
from bioflow_sim.generators.read_models import TECHNOLOGIES, introduce_errors
from bioflow_sim.generators.sequences import sample_genomic_template


def test_expression_generators_return_values_without_io() -> None:
    records = [
        ('tx2', 'tx2 gene=GENE2', 'ACGT' * 30),
        ('tx1', 'tx1 gene=GENE1', 'TGCA' * 40),
    ]
    transcripts = build_transcripts(records, {})
    assert transcripts == [
        Transcript('tx1', 'GENE1', 'TGCA' * 40),
        Transcript('tx2', 'GENE2', 'ACGT' * 30),
    ]

    rng = random.Random(4)
    barcodes = unique_barcodes(rng, 5)
    assert len(barcodes) == len(set(barcodes)) == 5
    assert {len(barcode) for barcode in barcodes} == {16}

    weights = cell_expression_weights([1.0, 1.0], cell_type=1)
    selected = choose_transcript(rng, transcripts, weights, minimum_length=100)
    assert isinstance(selected, Transcript)
    fragment = sample_three_prime_fragment(rng, selected, read_length=50)
    assert len(fragment.sequence) == 50
    assert fragment.transcript == selected


def test_read_generators_are_seeded_and_bounded() -> None:
    first = random.Random(9)
    second = random.Random(9)
    assert sample_positive_normal(first, 500, 50, 100, 1000) == sample_positive_normal(second, 500, 50, 100, 1000)

    technology = TECHNOLOGIES['pacbio-hifi']
    observed1 = introduce_errors('ACGT' * 100, first, technology)
    observed2 = introduce_errors('ACGT' * 100, second, technology)
    assert observed1 == observed2

    template = sample_genomic_template(random.Random(2), [('chr1', 'ACGT' * 100)], 80)
    assert template.source_id == 'chr1'
    assert template.end - template.start == 80
