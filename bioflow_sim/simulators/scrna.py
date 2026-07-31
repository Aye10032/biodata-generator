import random
from collections import Counter
from contextlib import ExitStack
from pathlib import Path

from loguru import logger

from ..core.io import (
    fastq_writer,
    prepare_output_directory,
    read_fasta,
    read_gtf_transcript_gene_map,
    validate_sample_name,
    write_fastq_record,
    write_json,
    write_tsv,
)
from ..generators.expression import (
    build_transcripts,
    cell_expression_weights,
    choose_transcript,
    sample_three_prime_fragment,
    sample_transcript_fragment,
    unique_barcodes,
)
from ..generators.random_values import random_dna, zipf_weights
from ..generators.read_models import TECHNOLOGIES, introduce_errors
from ..generators.sequences import reverse_complement

SCRNA_PROTOCOLS = ('10x-3prime', 'smartseq2')


def _write_matrix(
    output_dir: Path,
    cell_ids: list[str],
    features: list[tuple[str, str]],
    counts: Counter[tuple[int, int]],
    feature_type: str,
) -> None:
    matrix_dir = output_dir / 'truth' / 'expression_matrix'
    matrix_dir.mkdir(parents=True, exist_ok=True)
    with (matrix_dir / 'features.tsv').open('w', encoding='utf-8') as handle:
        for feature_id, gene_id in features:
            handle.write(f'{feature_id}\t{gene_id}\t{feature_type}\n')
    nonzero = [(feature, cell, value) for (feature, cell), value in counts.items() if value]
    nonzero.sort()
    with (matrix_dir / 'matrix.mtx').open('w', encoding='ascii') as handle:
        handle.write('%%MatrixMarket matrix coordinate integer general\n%\n')
        handle.write(f'{len(features)} {len(cell_ids)} {len(nonzero)}\n')
        for feature, cell, value in nonzero:
            handle.write(f'{feature + 1} {cell + 1} {value}\n')


def simulate_scrna(
    *,
    transcripts_path: Path,
    annotation_path: Path | None,
    output_dir: Path,
    sample: str,
    protocol: str,
    cells: int,
    reads_per_cell: int,
    seed: int,
    cdna_read_length: int,
    smartseq_read_length: int,
    fragment_mean: int,
    fragment_sd: int,
) -> dict[str, object]:
    if cells < 1 or reads_per_cell < 1:
        raise ValueError('cells and reads-per-cell must be at least 1')
    validate_sample_name(sample)
    fasta_records = read_fasta(transcripts_path)
    annotation_mapping = read_gtf_transcript_gene_map(annotation_path) if annotation_path is not None else {}
    if annotation_path is not None:
        mapped = sum(identifier in annotation_mapping for identifier, _, _ in fasta_records)
        mapping_rate = mapped / len(fasta_records)
        if mapping_rate < 0.95:
            raise ValueError(
                f'only {mapped}/{len(fasta_records)} transcript FASTA identifiers map to GTF '
                'transcript_id values; check that reference versions match'
            )
        feature_type = 'Gene Expression'
    else:
        mapping_rate = 0.0
        feature_type = 'Transcript Expression'
        logger.warning('No GTF annotation supplied; truth matrix features are transcripts, not genes')
    transcripts = build_transcripts(fasta_records, annotation_mapping)
    weights = zipf_weights(len(transcripts))
    rng = random.Random(seed)
    barcodes = unique_barcodes(rng, cells)
    cell_ids = [f'CELL{index:04d}-{barcode}' for index, barcode in enumerate(barcodes, 1)]
    cell_types = [index % 2 + 1 for index in range(cells)]
    counts: Counter[tuple[int, int]] = Counter()
    genes = sorted({record.feature_id for record in transcripts})
    gene_index = {gene: index for index, gene in enumerate(genes)}
    prepare_output_directory(output_dir)
    logger.info('Loaded {} transcript sequences from {}', len(transcripts), transcripts_path)

    read_truth: list[list[object]] = []
    illumina = TECHNOLOGIES['illumina-pe']
    if protocol == '10x-3prime':
        raw_dir = output_dir / 'raw'
        with ExitStack() as stack:
            r1 = stack.enter_context(fastq_writer(raw_dir / f'{sample}_R1.fastq.gz'))
            r2 = stack.enter_context(fastq_writer(raw_dir / f'{sample}_R2.fastq.gz'))
            ordinal = 0
            for cell_index, (cell_id, barcode) in enumerate(zip(cell_ids, barcodes, strict=True)):
                current_weights = cell_expression_weights(weights, cell_types[cell_index])
                for molecule in range(1, reads_per_cell + 1):
                    ordinal += 1
                    transcript = choose_transcript(
                        rng, transcripts, current_weights, cdna_read_length
                    )
                    transcript_id = transcript.identifier
                    gene_id = transcript.feature_id
                    umi = random_dna(rng, 12)
                    fragment = sample_three_prime_fragment(
                        rng, transcript, cdna_read_length
                    )
                    observed, substitutions, insertions, deletions = introduce_errors(
                        fragment.sequence, rng, illumina
                    )
                    read_id = f'{sample}:{ordinal:09d}'
                    barcode_read = barcode + umi
                    write_fastq_record(r1, f'{read_id}/1', barcode_read, 'I' * len(barcode_read))
                    write_fastq_record(r2, f'{read_id}/2', observed, 'I' * len(observed))
                    counts[(gene_index[gene_id], cell_index)] += 1
                    read_truth.append(
                        [
                            read_id,
                            cell_id,
                            barcode,
                            umi,
                            transcript_id,
                            gene_id,
                            fragment.start,
                            substitutions,
                            insertions,
                            deletions,
                        ]
                    )
        output_files = [f'raw/{sample}_R1.fastq.gz', f'raw/{sample}_R2.fastq.gz']
        barcode_path = output_dir / 'raw' / 'barcodes.tsv'
        with barcode_path.open('w', encoding='ascii') as handle:
            for barcode in barcodes:
                handle.write(f'{barcode}\n')
        output_files.append(str(barcode_path.relative_to(output_dir)))
    elif protocol == 'smartseq2':
        raw_dir = output_dir / 'raw'
        output_files: list[str] = []
        for cell_index, cell_id in enumerate(cell_ids):
            current_weights = cell_expression_weights(weights, cell_types[cell_index])
            r1_path = raw_dir / f'{cell_id}_R1.fastq.gz'
            r2_path = raw_dir / f'{cell_id}_R2.fastq.gz'
            output_files.extend([str(r1_path.relative_to(output_dir)), str(r2_path.relative_to(output_dir))])
            with fastq_writer(r1_path) as r1, fastq_writer(r2_path) as r2:
                for molecule in range(1, reads_per_cell + 1):
                    transcript = choose_transcript(
                        rng, transcripts, current_weights, smartseq_read_length * 2
                    )
                    transcript_id = transcript.identifier
                    gene_id = transcript.feature_id
                    fragment = sample_transcript_fragment(
                        rng,
                        transcript,
                        smartseq_read_length * 2,
                        fragment_mean,
                        fragment_sd,
                    )
                    read1, s1, i1, d1 = introduce_errors(
                        fragment.sequence[:smartseq_read_length], rng, illumina
                    )
                    read2, s2, i2, d2 = introduce_errors(
                        reverse_complement(fragment.sequence[-smartseq_read_length:]),
                        rng,
                        illumina,
                    )
                    read_id = f'{cell_id}:{molecule:08d}'
                    write_fastq_record(r1, f'{read_id}/1', read1, 'I' * len(read1))
                    write_fastq_record(r2, f'{read_id}/2', read2, 'I' * len(read2))
                    counts[(gene_index[gene_id], cell_index)] += 1
                    read_truth.append(
                        [
                            read_id,
                            cell_id,
                            '.',
                            '.',
                            transcript_id,
                            gene_id,
                            fragment.start,
                            s1 + s2,
                            i1 + i2,
                            d1 + d2,
                        ]
                    )
    else:
        raise ValueError(f'unsupported protocol: {protocol}')

    _write_matrix(
        output_dir,
        cell_ids,
        [(gene, gene) for gene in genes],
        counts,
        feature_type,
    )
    write_tsv(
        output_dir / 'truth' / 'reads.tsv',
        [
            'read_id',
            'cell_id',
            'barcode',
            'umi',
            'transcript_id',
            'gene_id',
            'transcript_start_0based',
            'substitutions',
            'insertions',
            'deletions',
        ],
        read_truth,
    )
    write_tsv(
        output_dir / 'truth' / 'cells.tsv',
        ['cell_id', 'barcode', 'cell_type'],
        (
            [cell_id, barcode, f'type_{cell_type}']
            for cell_id, barcode, cell_type in zip(cell_ids, barcodes, cell_types, strict=True)
        ),
    )
    metadata: dict[str, object] = {
        'schema_version': 2,
        'assay': 'single-cell-rna-sequencing',
        'sample': sample,
        'protocol': protocol,
        'transcript_reference': str(transcripts_path.resolve()),
        'annotation': str(annotation_path.resolve()) if annotation_path else None,
        'annotation_mapping_rate': mapping_rate,
        'matrix_feature_type': feature_type,
        'seed': seed,
        'cells': cells,
        'reads_per_cell': reads_per_cell,
        'output_files': output_files,
        'truth_matrix': {
            'matrix': 'truth/expression_matrix/matrix.mtx',
            'features': 'truth/expression_matrix/features.tsv',
            'cells': 'truth/cells.tsv',
            'barcodes': 'raw/barcodes.tsv' if protocol == '10x-3prime' else None,
        },
        'barcode_structure': '16 bp cell barcode + 12 bp UMI' if protocol == '10x-3prime' else None,
        'parameters': {
            'cdna_read_length': cdna_read_length,
            'smartseq_read_length': smartseq_read_length,
            'fragment_mean': fragment_mean,
            'fragment_sd': fragment_sd,
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success('Generated {} cells using {} in {}', cells, protocol, output_dir)
    return metadata
