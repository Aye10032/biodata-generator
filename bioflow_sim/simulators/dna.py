import random
from contextlib import ExitStack
from pathlib import Path

from loguru import logger

from bioflow_sim.core.io import (
    fastq_writer,
    prepare_output_directory,
    read_fasta,
    validate_sample_name,
    write_fastq_record,
    write_json,
    write_rows,
    write_tsv,
)
from bioflow_sim.generators.random_values import sample_positive_normal
from bioflow_sim.generators.read_models import TECHNOLOGIES, introduce_errors
from bioflow_sim.generators.sequences import (
    DNA_ALPHABET,
    reverse_complement,
    sample_genomic_template,
)
from bioflow_sim.generators.variants import apply_snvs, generate_snvs

DNA_TECHNOLOGY_NAMES = tuple(sorted(TECHNOLOGIES))


def simulate_dna(
    *,
    reference: Path,
    output_dir: Path,
    sample: str,
    technology_name: str,
    reads: int,
    seed: int,
    read_length: int | None,
    fragment_mean: int,
    fragment_sd: int,
    long_read_mean: int | None,
    long_read_sd: int | None,
    snvs: int,
) -> dict[str, object]:
    if reads < 1:
        raise ValueError('reads must be at least 1')
    validate_sample_name(sample)
    technology = TECHNOLOGIES[technology_name]
    records = read_fasta(reference)
    contigs = [(name, sequence) for name, _, sequence in records]
    usable = [(name, seq) for name, seq in contigs if any(b in DNA_ALPHABET for b in seq)]
    if not usable:
        raise ValueError('reference has no usable DNA sequences')

    rng = random.Random(seed)
    variants = generate_snvs(rng, usable, snvs)
    simulation_sequences = apply_snvs(usable, variants)
    prepare_output_directory(output_dir)
    raw_dir = output_dir / 'raw'
    truth_path = output_dir / 'truth' / 'reads.tsv'
    truth_rows: list[list[object]] = []
    logger.info('Loaded {} reference sequences from {}', len(usable), reference)

    with ExitStack() as stack:
        if technology.paired:
            r1_path = raw_dir / f'{sample}_R1.fastq.gz'
            r2_path = raw_dir / f'{sample}_R2.fastq.gz'
            r1 = stack.enter_context(fastq_writer(r1_path))
            r2 = stack.enter_context(fastq_writer(r2_path))
            output_files = [str(r1_path.relative_to(output_dir)), str(r2_path.relative_to(output_dir))]
        else:
            read_path = raw_dir / f'{sample}.fastq.gz'
            single = stack.enter_context(fastq_writer(read_path))
            output_files = [str(read_path.relative_to(output_dir))]

        for number in range(1, reads + 1):
            if technology.paired:
                length = read_length or technology.mean_length
                eligible = [(name, seq) for name, seq in simulation_sequences if len(seq) >= length * 2]
                if not eligible:
                    raise ValueError(f'no contig is long enough for paired reads of {length} bp')
                max_fragment = max(len(seq) for _, seq in eligible)
                fragment_length = sample_positive_normal(rng, fragment_mean, fragment_sd, length * 2, max_fragment)
                template = sample_genomic_template(rng, eligible, fragment_length)
                fragment = template.sequence
                read1, s1, i1, d1 = introduce_errors(fragment[:length], rng, technology)
                read2, s2, i2, d2 = introduce_errors(reverse_complement(fragment[-length:]), rng, technology)
                read_id = f'{sample}:{number:08d}'
                write_fastq_record(r1, f'{read_id}/1', read1, technology.quality_char * len(read1))
                write_fastq_record(r2, f'{read_id}/2', read2, technology.quality_char * len(read2))
                truth_rows.append(
                    [
                        read_id,
                        template.source_id,
                        template.start,
                        template.end,
                        template.strand,
                        fragment_length,
                        len(read1) + len(read2),
                        s1 + s2,
                        i1 + i2,
                        d1 + d2,
                    ]
                )
            else:
                if technology_name == 'illumina-se':
                    mean = read_length or technology.mean_length
                    sd = 0
                    minimum_length = mean
                else:
                    mean = long_read_mean or read_length or technology.mean_length
                    sd = long_read_sd if long_read_sd is not None else technology.length_sd
                    minimum_length = 100
                max_length = max(len(seq) for _, seq in simulation_sequences)
                template_length = sample_positive_normal(rng, mean, sd, minimum_length, max_length)
                template = sample_genomic_template(rng, simulation_sequences, template_length)
                observed, substitutions, insertions, deletions = introduce_errors(template.sequence, rng, technology)
                read_id = f'{sample}:{number:08d}'
                write_fastq_record(
                    single,
                    read_id,
                    observed,
                    technology.quality_char * len(observed),
                )
                truth_rows.append(
                    [
                        read_id,
                        template.source_id,
                        template.start,
                        template.end,
                        template.strand,
                        template_length,
                        len(observed),
                        substitutions,
                        insertions,
                        deletions,
                    ]
                )

    write_tsv(
        truth_path,
        [
            'read_id',
            'contig',
            'start_0based',
            'end_0based_exclusive',
            'strand',
            'template_length',
            'observed_bases',
            'substitutions',
            'insertions',
            'deletions',
        ],
        truth_rows,
    )
    if variants:
        variant_path = output_dir / 'truth' / 'variants.vcf'
        write_rows(
            variant_path,
            [
                ['##fileformat=VCFv4.3'],
                ['##source=bioflow-sim'],
                ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'],
                *[
                    [
                        variant.contig,
                        variant.position + 1,
                        f'SNV{index:06d}',
                        variant.reference,
                        variant.alternate,
                        '.',
                        'PASS',
                        'SIMULATED=1',
                    ]
                    for index, variant in enumerate(variants, 1)
                ],
            ],
        )
    metadata: dict[str, object] = {
        'schema_version': 1,
        'assay': 'dna-sequencing',
        'sample': sample,
        'technology': technology_name,
        'reference': str(reference.resolve()),
        'seed': seed,
        'read_units': reads,
        'read_unit_definition': 'read pairs' if technology.paired else 'reads',
        'output_files': output_files,
        'truth_file': str(truth_path.relative_to(output_dir)),
        'variant_truth': 'truth/variants.vcf' if variants else None,
        'coordinate_system': '0-based half-open',
        'parameters': {
            'read_length': read_length,
            'fragment_mean': fragment_mean,
            'fragment_sd': fragment_sd,
            'long_read_mean': long_read_mean,
            'long_read_sd': long_read_sd,
            'snvs': snvs,
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success('Generated {} {} in {}', reads, metadata['read_unit_definition'], output_dir)
    return metadata
