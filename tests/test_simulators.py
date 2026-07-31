from __future__ import annotations

import gzip
import json
from pathlib import Path

from click.testing import CliRunner

from bioflow_sim.cli import cli


def write_reference(path: Path) -> None:
    path.write_text(
        '>chrA synthetic\n' + 'ACGT' * 1_000 + '\n'
        '>tx1 gene=GENE1\n' + 'ACGT' * 200 + '\n'
        '>tx2 gene=GENE2\n' + 'TGCA' * 180 + '\n',
        encoding='utf-8',
    )


def write_transcripts(path: Path) -> None:
    path.write_text(
        '>tx1 gene=GENE1\n' + 'ACGT' * 200 + '\n>tx2 gene=GENE2\n' + 'TGCA' * 180 + '\n',
        encoding='utf-8',
    )


def write_annotation(path: Path) -> None:
    path.write_text(
        'chrA\ttest\ttranscript\t1\t800\t.\t+\t.\t'
        'gene_id "GENE1"; transcript_id "tx1";\n'
        'chrA\ttest\ttranscript\t1\t720\t.\t+\t.\t'
        'gene_id "GENE2"; transcript_id "tx2";\n',
        encoding='utf-8',
    )


def write_assay_reference(path: Path) -> None:
    path.write_text(
        '>chr1\n' + 'AAGCTTGATCCGATCG' * 300 + '\n>chr2\n' + 'CGGATCAAGCTTATCG' * 250 + '\n',
        encoding='utf-8',
    )


def write_targets(path: Path) -> None:
    path.write_text(
        'chr1\t100\t700\tGENE1\nchr1\t800\t1400\tGENE2\nchr2\t100\t700\tGENE3\n',
        encoding='utf-8',
    )


def fastq_records(path: Path) -> list[list[str]]:
    with gzip.open(path, 'rt', encoding='ascii') as handle:
        lines = [line.rstrip('\n') for line in handle]
    assert len(lines) % 4 == 0
    records = [lines[index : index + 4] for index in range(0, len(lines), 4)]
    for name, sequence, plus, quality in records:
        assert name.startswith('@')
        assert plus == '+'
        assert len(sequence) == len(quality)
    return records


def test_illumina_pe_is_reproducible(tmp_path: Path) -> None:
    reference = tmp_path / 'reference.fa'
    write_reference(reference)
    runner = CliRunner()
    arguments = [
        'dna',
        '--reference',
        str(reference),
        '--technology',
        'illumina-pe',
        '--reads',
        '12',
        '--read-length',
        '75',
        '--seed',
        '42',
    ]
    first = tmp_path / 'first'
    second = tmp_path / 'second'
    result = runner.invoke(cli, [*arguments, '--output-dir', str(first)])
    assert result.exit_code == 0, result.output
    result = runner.invoke(cli, [*arguments, '--output-dir', str(second)])
    assert result.exit_code == 0, result.output
    assert (first / 'raw/S1_R1.fastq.gz').read_bytes() == (second / 'raw/S1_R1.fastq.gz').read_bytes()
    assert len(fastq_records(first / 'raw/S1_R1.fastq.gz')) == 12
    assert len(fastq_records(first / 'raw/S1_R2.fastq.gz')) == 12


def test_long_read_technologies(tmp_path: Path) -> None:
    reference = tmp_path / 'reference.fa'
    write_reference(reference)
    runner = CliRunner()
    for technology in ('pacbio-clr', 'pacbio-hifi', 'ont'):
        output = tmp_path / technology
        result = runner.invoke(
            cli,
            [
                'dna',
                '--reference',
                str(reference),
                '--output-dir',
                str(output),
                '--technology',
                technology,
                '--reads',
                '3',
                '--long-read-mean',
                '500',
                '--long-read-sd',
                '20',
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(fastq_records(output / 'raw/S1.fastq.gz')) == 3


def test_illumina_single_end_length(tmp_path: Path) -> None:
    reference = tmp_path / 'reference.fa'
    write_reference(reference)
    output = tmp_path / 'single'
    result = CliRunner().invoke(
        cli,
        [
            'dna',
            '--reference',
            str(reference),
            '--output-dir',
            str(output),
            '--technology',
            'illumina-se',
            '--reads',
            '4',
            '--read-length',
            '75',
        ],
    )
    assert result.exit_code == 0, result.output
    assert {len(record[1]) for record in fastq_records(output / 'raw/S1.fastq.gz')} == {75}


def test_10x_layout_and_truth_matrix(tmp_path: Path) -> None:
    reference = tmp_path / 'transcripts.fa'
    annotation = tmp_path / 'genes.gtf'
    write_transcripts(reference)
    write_annotation(annotation)
    output = tmp_path / 'tenx'
    result = CliRunner().invoke(
        cli,
        [
            'scrna',
            '--transcripts-path',
            str(reference),
            '--output-dir',
            str(output),
            '--annotation-path',
            str(annotation),
            '--protocol',
            '10x-3prime',
            '--cells',
            '3',
            '--reads-per-cell',
            '7',
            '--seed',
            '9',
        ],
    )
    assert result.exit_code == 0, result.output
    r1 = fastq_records(output / 'raw/SC1_R1.fastq.gz')
    r2 = fastq_records(output / 'raw/SC1_R2.fastq.gz')
    assert len(r1) == len(r2) == 21
    assert {len(record[1]) for record in r1} == {28}
    assert {len(record[1]) for record in r2} == {90}
    manifest = json.loads((output / 'manifest.json').read_text(encoding='utf-8'))
    assert manifest['barcode_structure'] == '16 bp cell barcode + 12 bp UMI'
    assert manifest['matrix_feature_type'] == 'Gene Expression'
    with (output / 'truth/expression_matrix/features.tsv').open('r', encoding='utf-8') as handle:
        assert len(handle.readlines()) == 2
    assert (output / 'truth/expression_matrix/matrix.mtx').exists()
    barcodes = (output / 'raw/barcodes.tsv').read_text(encoding='ascii').splitlines()
    assert len(barcodes) == 3
    assert {len(barcode) for barcode in barcodes} == {16}
    assert not list((output / 'truth').rglob('*.gz'))


def test_smartseq2_writes_one_pair_per_cell(tmp_path: Path) -> None:
    reference = tmp_path / 'transcripts.fa'
    write_transcripts(reference)
    output = tmp_path / 'smartseq'
    result = CliRunner().invoke(
        cli,
        [
            'scrna',
            '--transcripts-path',
            str(reference),
            '--output-dir',
            str(output),
            '--protocol',
            'smartseq2',
            '--cells',
            '2',
            '--reads-per-cell',
            '5',
            '--smartseq-read-length',
            '60',
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(list((output / 'raw').glob('*_R1.fastq.gz'))) == 2
    for path in (output / 'raw').glob('*_R1.fastq.gz'):
        assert len(fastq_records(path)) == 5


def test_unsafe_sample_name_is_rejected(tmp_path: Path) -> None:
    reference = tmp_path / 'reference.fa'
    write_reference(reference)
    result = CliRunner().invoke(
        cli,
        [
            'dna',
            '--reference',
            str(reference),
            '--output-dir',
            str(tmp_path / 'output'),
            '--sample',
            '../escape',
            '--technology',
            'illumina-se',
            '--reads',
            '1',
        ],
    )
    assert result.exit_code == 1
    assert 'sample must start' in result.output


def test_dna_snv_truth(tmp_path: Path) -> None:
    reference = tmp_path / 'reference.fa'
    write_reference(reference)
    output = tmp_path / 'variants'
    result = CliRunner().invoke(
        cli,
        [
            'dna',
            '--reference',
            str(reference),
            '--output-dir',
            str(output),
            '--technology',
            'illumina-se',
            '--reads',
            '5',
            '--snvs',
            '4',
            '--seed',
            '2',
        ],
    )
    assert result.exit_code == 0, result.output
    variants = (output / 'truth/variants.vcf').read_text(encoding='utf-8').splitlines()
    assert len([line for line in variants if not line.startswith('#')]) == 4


def test_tumor_normal_pair_has_expected_truth_and_no_sample_sheet(tmp_path: Path) -> None:
    reference = tmp_path / 'assay.fa'
    targets = tmp_path / 'genes.bed'
    write_assay_reference(reference)
    write_targets(targets)
    output = tmp_path / 'tumor_normal'
    result = CliRunner().invoke(
        cli,
        [
            'tumor-normal',
            '--reference',
            str(reference),
            '--candidate-targets',
            str(targets),
            '--output-dir',
            str(output),
            '--pair-name',
            'P01',
            '--normal-sample',
            'P01_N',
            '--tumor-sample',
            'P01_T',
            '--target-count',
            '2',
            '--target-width',
            '400',
            '--normal-depth',
            '2',
            '--tumor-depth',
            '4',
            '--tumor-purity',
            '0.6',
            '--germline-snvs',
            '2',
            '--clonal-snvs',
            '2',
            '--subclonal-snvs',
            '2',
            '--subclone-fraction',
            '0.25',
            '--read-length',
            '50',
            '--fragment-mean',
            '150',
            '--fragment-sd',
            '10',
            '--seed',
            '31',
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(fastq_records(output / 'raw/P01_N_R1.fastq.gz')) == 16
    assert len(fastq_records(output / 'raw/P01_T_R1.fastq.gz')) == 32
    assert (output / 'raw/targets.bed').exists()
    assert not (output / 'raw/samples.tsv').exists()

    germline = (output / 'truth/germline.vcf').read_text(encoding='utf-8').splitlines()
    somatic = (output / 'truth/somatic.vcf').read_text(encoding='utf-8').splitlines()
    assert len([line for line in germline if not line.startswith('#')]) == 2
    assert len([line for line in somatic if not line.startswith('#')]) == 4

    vaf_rows = [
        line.split('\t') for line in (output / 'truth/expected_vaf.tsv').read_text(encoding='utf-8').splitlines()[1:]
    ]
    observed = {(row[1], row[4], row[5], row[6]) for row in vaf_rows}
    assert ('germline', '0.500000', '0.500000', '1.000000') in observed
    assert ('somatic_clonal', '0.000000', '0.300000', '1.000000') in observed
    assert ('somatic_subclonal', '0.000000', '0.075000', '0.250000') in observed
    assert (output / 'truth/clones.tsv').exists()
    assert (output / 'truth/copy_number.tsv').exists()
    assert (output / 'truth/read_origins.tsv').exists()


def test_bulk_rna_command(tmp_path: Path) -> None:
    transcripts = tmp_path / 'transcripts.fa'
    annotation = tmp_path / 'genes.gtf'
    write_transcripts(transcripts)
    write_annotation(annotation)
    output = tmp_path / 'bulk'
    result = CliRunner().invoke(
        cli,
        [
            'bulk-rna',
            '--transcripts-path',
            str(transcripts),
            '--annotation-path',
            str(annotation),
            '--output-dir',
            str(output),
            '--samples-per-group',
            '1',
            '--reads-per-sample',
            '4',
            '--read-length',
            '50',
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(list((output / 'raw').glob('*_R1.fastq.gz'))) == 2
    assert (output / 'truth/expression.tsv').exists()


def test_chromatin_methylation_and_hic_commands(tmp_path: Path) -> None:
    reference = tmp_path / 'assay.fa'
    write_assay_reference(reference)
    runner = CliRunner()

    chromatin = tmp_path / 'chromatin'
    result = runner.invoke(
        cli,
        [
            'chromatin',
            '--reference',
            str(reference),
            '--output-dir',
            str(chromatin),
            '--assay',
            'atac',
            '--reads',
            '5',
            '--peaks',
            '2',
            '--peak-width',
            '50',
            '--read-length',
            '30',
            '--fragment-mean',
            '100',
        ],
    )
    assert result.exit_code == 0, result.output
    assert len((chromatin / 'truth/peaks.bed').read_text().splitlines()) == 2

    methylation = tmp_path / 'methylation'
    result = runner.invoke(
        cli,
        [
            'methylation',
            '--reference',
            str(reference),
            '--output-dir',
            str(methylation),
            '--protocol',
            'wgbs',
            '--reads',
            '5',
            '--sites',
            '4',
            '--read-length',
            '30',
            '--fragment-mean',
            '100',
        ],
    )
    assert result.exit_code == 0, result.output
    assert len((methylation / 'truth/methylation.tsv').read_text().splitlines()) == 5

    hic = tmp_path / 'hic'
    result = runner.invoke(
        cli,
        [
            'hic',
            '--reference',
            str(reference),
            '--output-dir',
            str(hic),
            '--enzyme',
            'mboi',
            '--reads',
            '5',
            '--read-length',
            '30',
        ],
    )
    assert result.exit_code == 0, result.output
    assert (hic / 'raw/restriction_sites.tsv').exists()
    assert (hic / 'truth/contacts.tsv').exists()
