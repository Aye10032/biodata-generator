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
    write_tsv,
)
from bioflow_sim.generators.contacts import (
    ENZYMES,
    contact_read_sequences,
    find_restriction_sites,
    sample_contact,
)
from bioflow_sim.generators.read_models import TECHNOLOGIES, introduce_errors

HIC_ENZYMES = tuple(sorted(ENZYMES))


def simulate_hic(
    *,
    reference: Path,
    output_dir: Path,
    sample: str,
    enzyme_name: str,
    reads: int,
    read_length: int,
    intra_rate: float,
    mean_distance: int,
    seed: int,
) -> dict[str, object]:
    if reads < 1 or mean_distance < 1:
        raise ValueError('reads and mean distance must be positive')
    validate_sample_name(sample)
    sequences = [(name, sequence) for name, _, sequence in read_fasta(reference)]
    enzyme = ENZYMES[enzyme_name]
    cut_sites = find_restriction_sites(sequences, enzyme)
    prepare_output_directory(output_dir)
    raw_dir = output_dir / 'raw'
    rng = random.Random(seed)
    illumina = TECHNOLOGIES['illumina-pe']
    r1_path = raw_dir / f'{sample}_R1.fastq.gz'
    r2_path = raw_dir / f'{sample}_R2.fastq.gz'
    truth_rows: list[list[object]] = []

    with ExitStack() as stack:
        r1 = stack.enter_context(fastq_writer(r1_path))
        r2 = stack.enter_context(fastq_writer(r2_path))
        for ordinal in range(1, reads + 1):
            contact = sample_contact(rng, cut_sites, intra_rate, mean_distance)
            template1, template2 = contact_read_sequences(sequences, contact, read_length)
            read1, s1, i1, d1 = introduce_errors(template1, rng, illumina)
            read2, s2, i2, d2 = introduce_errors(template2, rng, illumina)
            read_id = f'{sample}:{ordinal:09d}'
            write_fastq_record(r1, f'{read_id}/1', read1, 'I' * len(read1))
            write_fastq_record(r2, f'{read_id}/2', read2, 'I' * len(read2))
            distance = abs(contact.right.position - contact.left.position) if contact.contact_type == 'intra' else '.'
            truth_rows.append(
                [
                    read_id,
                    contact.left.contig,
                    contact.left.position,
                    contact.right.contig,
                    contact.right.position,
                    contact.contact_type,
                    distance,
                    s1 + s2,
                    i1 + i2,
                    d1 + d2,
                ]
            )

    write_tsv(
        raw_dir / 'restriction_sites.tsv',
        ['contig', 'cut_position_0based', 'enzyme', 'motif'],
        ([site.contig, site.position, enzyme.name, enzyme.motif] for site in cut_sites),
    )
    write_tsv(
        raw_dir / 'library.tsv',
        ['sample_id', 'assay', 'enzyme', 'layout', 'read_length'],
        [[sample, 'Hi-C', enzyme.name, 'paired-end', read_length]],
    )
    write_tsv(
        output_dir / 'truth' / 'contacts.tsv',
        [
            'read_id',
            'left_contig',
            'left_cut_0based',
            'right_contig',
            'right_cut_0based',
            'contact_type',
            'distance',
            'substitutions',
            'insertions',
            'deletions',
        ],
        truth_rows,
    )
    metadata: dict[str, object] = {
        'schema_version': 1,
        'assay': 'hi-c',
        'sample': sample,
        'reference': str(reference.resolve()),
        'seed': seed,
        'read_pairs': reads,
        'output_files': [
            str(r1_path.relative_to(output_dir)),
            str(r2_path.relative_to(output_dir)),
            'raw/restriction_sites.tsv',
            'raw/library.tsv',
        ],
        'truth_files': ['truth/contacts.tsv'],
        'parameters': {
            'enzyme': enzyme.name,
            'motif': enzyme.motif,
            'read_length': read_length,
            'intra_rate': intra_rate,
            'mean_distance': mean_distance,
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success('Generated {} Hi-C read pairs in {}', reads, output_dir)
    return metadata
