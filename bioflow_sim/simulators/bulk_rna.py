import random
from collections import Counter
from contextlib import ExitStack
from pathlib import Path

from loguru import logger

from bioflow_sim.core.io import (
    fastq_writer,
    prepare_output_directory,
    read_fasta,
    read_gtf_transcript_gene_map,
    validate_sample_name,
    write_fastq_record,
    write_json,
    write_tsv,
)
from bioflow_sim.generators.expression import (
    build_transcripts,
    choose_transcript,
    differential_expression_weights,
    sample_transcript_fragment,
)
from bioflow_sim.generators.random_values import zipf_weights
from bioflow_sim.generators.read_models import TECHNOLOGIES, introduce_errors
from bioflow_sim.generators.sequences import reverse_complement

BULK_RNA_LAYOUTS = ('se', 'pe')
BULK_RNA_STRANDEDNESS = ('unstranded', 'forward', 'reverse')


def simulate_bulk_rna(
    *,
    transcripts_path: Path,
    annotation_path: Path | None,
    output_dir: Path,
    sample_prefix: str,
    samples_per_group: int,
    reads_per_sample: int,
    layout: str,
    strandedness: str,
    read_length: int,
    fragment_mean: int,
    fragment_sd: int,
    fold_change: float,
    seed: int,
) -> dict[str, object]:
    if samples_per_group < 1 or reads_per_sample < 1:
        raise ValueError('samples-per-group and reads-per-sample must be at least 1')
    validate_sample_name(sample_prefix)
    fasta_records = read_fasta(transcripts_path)
    annotation_mapping = read_gtf_transcript_gene_map(annotation_path) if annotation_path is not None else {}
    transcripts = build_transcripts(fasta_records, annotation_mapping)
    minimum_length = read_length * 2 if layout == 'pe' else read_length
    if not any(len(transcript.sequence) >= minimum_length for transcript in transcripts):
        raise ValueError(f'no transcript is at least {minimum_length} bp')

    prepare_output_directory(output_dir)
    raw_dir = output_dir / 'raw'
    rng = random.Random(seed)
    illumina = TECHNOLOGIES['illumina-pe' if layout == 'pe' else 'illumina-se']
    base_weights = zipf_weights(len(transcripts))
    samples = [
        (f'{sample_prefix}_{group}_{replicate}', group, replicate)
        for group in ('control', 'treatment')
        for replicate in range(1, samples_per_group + 1)
    ]
    output_files: list[str] = []
    read_truth: list[list[object]] = []
    counts: Counter[tuple[str, str]] = Counter()
    expression_rows: list[list[object]] = []

    for sample_id, group, _ in samples:
        weights = differential_expression_weights(base_weights, group, fold_change)
        with ExitStack() as stack:
            if layout == 'pe':
                r1_path = raw_dir / f'{sample_id}_R1.fastq.gz'
                r2_path = raw_dir / f'{sample_id}_R2.fastq.gz'
                r1 = stack.enter_context(fastq_writer(r1_path))
                r2 = stack.enter_context(fastq_writer(r2_path))
                output_files.extend(
                    [
                        str(r1_path.relative_to(output_dir)),
                        str(r2_path.relative_to(output_dir)),
                    ]
                )
            else:
                read_path = raw_dir / f'{sample_id}.fastq.gz'
                single = stack.enter_context(fastq_writer(read_path))
                output_files.append(str(read_path.relative_to(output_dir)))

            for ordinal in range(1, reads_per_sample + 1):
                transcript = choose_transcript(rng, transcripts, weights, minimum_length)
                fragment = sample_transcript_fragment(
                    rng,
                    transcript,
                    minimum_length,
                    fragment_mean,
                    fragment_sd,
                )
                sequence = fragment.sequence
                reverse = strandedness == 'reverse' or (strandedness == 'unstranded' and rng.random() < 0.5)
                if reverse:
                    sequence = reverse_complement(sequence)
                    strand = '-'
                else:
                    strand = '+'
                read_id = f'{sample_id}:{ordinal:09d}'
                if layout == 'pe':
                    read1, s1, i1, d1 = introduce_errors(sequence[:read_length], rng, illumina)
                    read2, s2, i2, d2 = introduce_errors(reverse_complement(sequence[-read_length:]), rng, illumina)
                    write_fastq_record(r1, f'{read_id}/1', read1, 'I' * len(read1))
                    write_fastq_record(r2, f'{read_id}/2', read2, 'I' * len(read2))
                    errors = (s1 + s2, i1 + i2, d1 + d2)
                else:
                    read, substitutions, insertions, deletions = introduce_errors(sequence[:read_length], rng, illumina)
                    write_fastq_record(single, read_id, read, 'I' * len(read))
                    errors = (substitutions, insertions, deletions)
                counts[(sample_id, transcript.identifier)] += 1
                read_truth.append(
                    [
                        read_id,
                        sample_id,
                        group,
                        transcript.identifier,
                        transcript.feature_id,
                        fragment.start,
                        strand,
                        *errors,
                    ]
                )

        for transcript, weight in zip(transcripts, weights, strict=True):
            expression_rows.append(
                [
                    sample_id,
                    group,
                    transcript.identifier,
                    transcript.feature_id,
                    f'{weight:.8g}',
                    counts[(sample_id, transcript.identifier)],
                ]
            )

    sample_sheet = raw_dir / 'samples.tsv'
    write_tsv(
        sample_sheet,
        ['sample_id', 'group', 'replicate', 'layout', 'strandedness'],
        ([sample_id, group, replicate, layout, strandedness] for sample_id, group, replicate in samples),
    )
    output_files.append(str(sample_sheet.relative_to(output_dir)))
    write_tsv(
        output_dir / 'truth' / 'expression.tsv',
        ['sample_id', 'group', 'transcript_id', 'gene_id', 'weight', 'read_count'],
        expression_rows,
    )
    write_tsv(
        output_dir / 'truth' / 'reads.tsv',
        [
            'read_id',
            'sample_id',
            'group',
            'transcript_id',
            'gene_id',
            'transcript_start_0based',
            'strand',
            'substitutions',
            'insertions',
            'deletions',
        ],
        read_truth,
    )
    metadata: dict[str, object] = {
        'schema_version': 1,
        'assay': 'bulk-rna-sequencing',
        'transcript_reference': str(transcripts_path.resolve()),
        'annotation': str(annotation_path.resolve()) if annotation_path else None,
        'sample_prefix': sample_prefix,
        'samples_per_group': samples_per_group,
        'reads_per_sample': reads_per_sample,
        'groups': ['control', 'treatment'],
        'layout': layout,
        'strandedness': strandedness,
        'seed': seed,
        'output_files': output_files,
        'truth_files': ['truth/expression.tsv', 'truth/reads.tsv'],
        'parameters': {
            'read_length': read_length,
            'fragment_mean': fragment_mean,
            'fragment_sd': fragment_sd,
            'fold_change': fold_change,
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success('Generated {} bulk RNA samples in {}', len(samples), output_dir)
    return metadata
