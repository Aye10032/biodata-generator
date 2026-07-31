import random
from contextlib import ExitStack
from pathlib import Path

from loguru import logger

from ..core.io import (
    fastq_writer,
    prepare_output_directory,
    read_fasta,
    validate_sample_name,
    write_fastq_record,
    write_json,
    write_rows,
    write_tsv,
)
from ..generators.random_values import sample_positive_normal
from ..generators.read_models import TECHNOLOGIES, introduce_errors
from ..generators.regions import generate_intervals, sample_template_near_interval
from ..generators.sequences import reverse_complement, sample_genomic_template

CHROMATIN_ASSAYS = ('atac', 'chip', 'cuttag')


def simulate_chromatin(
    *,
    reference: Path,
    output_dir: Path,
    sample: str,
    assay: str,
    reads: int,
    peaks: int,
    peak_width: int,
    enrichment: float,
    read_length: int,
    fragment_mean: int,
    fragment_sd: int,
    seed: int,
) -> dict[str, object]:
    if reads < 1 or not 0 <= enrichment <= 1:
        raise ValueError('reads must be positive and enrichment must be between 0 and 1')
    validate_sample_name(sample)
    sequences = [(name, sequence) for name, _, sequence in read_fasta(reference)]
    rng = random.Random(seed)
    truth_peaks = generate_intervals(rng, sequences, peaks, peak_width, assay)
    prepare_output_directory(output_dir)
    raw_dir = output_dir / 'raw'
    illumina = TECHNOLOGIES['illumina-pe']
    truth_rows: list[list[object]] = []
    r1_path = raw_dir / f'{sample}_R1.fastq.gz'
    r2_path = raw_dir / f'{sample}_R2.fastq.gz'

    with ExitStack() as stack:
        r1 = stack.enter_context(fastq_writer(r1_path))
        r2 = stack.enter_context(fastq_writer(r2_path))
        for ordinal in range(1, reads + 1):
            max_length = max(len(sequence) for _, sequence in sequences)
            fragment_length = sample_positive_normal(
                rng,
                fragment_mean,
                fragment_sd,
                read_length * 2,
                max_length,
            )
            if rng.random() < enrichment:
                template, peak = sample_template_near_interval(rng, sequences, truth_peaks, fragment_length)
                peak_name = peak.name
                enriched = 1
            else:
                template = sample_genomic_template(rng, sequences, fragment_length)
                peak_name = '.'
                enriched = 0
            read1, s1, i1, d1 = introduce_errors(template.sequence[:read_length], rng, illumina)
            read2, s2, i2, d2 = introduce_errors(reverse_complement(template.sequence[-read_length:]), rng, illumina)
            read_id = f'{sample}:{ordinal:09d}'
            write_fastq_record(r1, f'{read_id}/1', read1, 'I' * len(read1))
            write_fastq_record(r2, f'{read_id}/2', read2, 'I' * len(read2))
            truth_rows.append(
                [
                    read_id,
                    template.source_id,
                    template.start,
                    template.end,
                    template.strand,
                    enriched,
                    peak_name,
                    s1 + s2,
                    i1 + i2,
                    d1 + d2,
                ]
            )

    write_tsv(
        raw_dir / 'library.tsv',
        ['sample_id', 'assay', 'layout', 'read_length'],
        [[sample, assay, 'paired-end', read_length]],
    )
    write_rows(
        output_dir / 'truth' / 'peaks.bed',
        ([peak.contig, peak.start, peak.end, peak.name, peak.score] for peak in truth_peaks),
    )
    write_tsv(
        output_dir / 'truth' / 'reads.tsv',
        [
            'read_id',
            'contig',
            'start_0based',
            'end_0based_exclusive',
            'strand',
            'enriched',
            'peak_id',
            'substitutions',
            'insertions',
            'deletions',
        ],
        truth_rows,
    )
    metadata: dict[str, object] = {
        'schema_version': 1,
        'assay': f'{assay}-seq' if assay != 'cuttag' else 'cut-and-tag',
        'sample': sample,
        'reference': str(reference.resolve()),
        'seed': seed,
        'read_pairs': reads,
        'output_files': [
            str(r1_path.relative_to(output_dir)),
            str(r2_path.relative_to(output_dir)),
            'raw/library.tsv',
        ],
        'truth_files': ['truth/peaks.bed', 'truth/reads.tsv'],
        'parameters': {
            'peaks': peaks,
            'peak_width': peak_width,
            'enrichment': enrichment,
            'read_length': read_length,
            'fragment_mean': fragment_mean,
            'fragment_sd': fragment_sd,
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success('Generated {} {} read pairs in {}', reads, assay, output_dir)
    return metadata
