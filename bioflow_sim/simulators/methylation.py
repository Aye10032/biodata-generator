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
    write_tsv,
)
from ..generators.methylation import (
    convert_unmethylated_cytosines,
    generate_cpg_methylation,
    sample_template_around_cpg,
)
from ..generators.random_values import sample_positive_normal
from ..generators.read_models import TECHNOLOGIES, introduce_errors
from ..generators.sequences import reverse_complement

METHYLATION_PROTOCOLS = ('wgbs', 'emseq')


def simulate_methylation(
    *,
    reference: Path,
    output_dir: Path,
    sample: str,
    protocol: str,
    reads: int,
    sites: int,
    methylation_rate: float,
    conversion_rate: float,
    read_length: int,
    fragment_mean: int,
    fragment_sd: int,
    seed: int,
) -> dict[str, object]:
    if reads < 1:
        raise ValueError('reads must be positive')
    validate_sample_name(sample)
    sequences = [(name, sequence) for name, _, sequence in read_fasta(reference)]
    rng = random.Random(seed)
    truth_sites = generate_cpg_methylation(rng, sequences, sites, methylation_rate)
    site_lookup = {(site.contig, site.position): site for site in truth_sites}
    prepare_output_directory(output_dir)
    raw_dir = output_dir / 'raw'
    r1_path = raw_dir / f'{sample}_R1.fastq.gz'
    r2_path = raw_dir / f'{sample}_R2.fastq.gz'
    illumina = TECHNOLOGIES['illumina-pe']
    read_truth: list[list[object]] = []

    with ExitStack() as stack:
        r1 = stack.enter_context(fastq_writer(r1_path))
        r2 = stack.enter_context(fastq_writer(r2_path))
        for ordinal in range(1, reads + 1):
            fragment_length = sample_positive_normal(
                rng,
                fragment_mean,
                fragment_sd,
                read_length * 2,
                max(len(sequence) for _, sequence in sequences),
            )
            template, focal_site = sample_template_around_cpg(rng, sequences, truth_sites, fragment_length)
            converted, conversions = convert_unmethylated_cytosines(rng, template, site_lookup, conversion_rate)
            read1, s1, i1, d1 = introduce_errors(converted[:read_length], rng, illumina)
            read2, s2, i2, d2 = introduce_errors(reverse_complement(converted[-read_length:]), rng, illumina)
            read_id = f'{sample}:{ordinal:09d}'
            write_fastq_record(r1, f'{read_id}/1', read1, 'I' * len(read1))
            write_fastq_record(r2, f'{read_id}/2', read2, 'I' * len(read2))
            read_truth.append(
                [
                    read_id,
                    template.source_id,
                    template.start,
                    template.end,
                    template.strand,
                    focal_site.position,
                    conversions,
                    s1 + s2,
                    i1 + i2,
                    d1 + d2,
                ]
            )

    write_tsv(
        raw_dir / 'library.tsv',
        ['sample_id', 'protocol', 'layout', 'conversion_rate'],
        [[sample, protocol, 'paired-end', conversion_rate]],
    )
    write_tsv(
        output_dir / 'truth' / 'methylation.tsv',
        ['contig', 'position_0based', 'context', 'methylated'],
        ([site.contig, site.position, 'CG', int(site.methylated)] for site in truth_sites),
    )
    write_tsv(
        output_dir / 'truth' / 'reads.tsv',
        [
            'read_id',
            'contig',
            'start_0based',
            'end_0based_exclusive',
            'strand',
            'focal_cpg_0based',
            'conversions',
            'substitutions',
            'insertions',
            'deletions',
        ],
        read_truth,
    )
    metadata: dict[str, object] = {
        'schema_version': 1,
        'assay': protocol,
        'sample': sample,
        'reference': str(reference.resolve()),
        'seed': seed,
        'read_pairs': reads,
        'output_files': [
            str(r1_path.relative_to(output_dir)),
            str(r2_path.relative_to(output_dir)),
            'raw/library.tsv',
        ],
        'truth_files': ['truth/methylation.tsv', 'truth/reads.tsv'],
        'parameters': {
            'sites': sites,
            'methylation_rate': methylation_rate,
            'conversion_rate': conversion_rate,
            'read_length': read_length,
            'fragment_mean': fragment_mean,
            'fragment_sd': fragment_sd,
            'unlisted_cytosines': 'treated as methylated',
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success('Generated {} {} read pairs in {}', reads, protocol, output_dir)
    return metadata
