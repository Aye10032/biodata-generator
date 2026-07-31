import math
import random
from contextlib import ExitStack
from pathlib import Path

from loguru import logger

from bioflow_sim.core.io import (
    fastq_writer,
    prepare_output_directory,
    read_bed,
    read_fasta,
    validate_sample_name,
    write_fastq_record,
    write_json,
    write_lines,
    write_rows,
    write_tsv,
)
from bioflow_sim.generators.random_values import sample_positive_normal
from bioflow_sim.generators.read_models import TECHNOLOGIES, introduce_errors
from bioflow_sim.generators.sequences import reverse_complement
from bioflow_sim.generators.tumor import (
    TumorNormalVariant,
    apply_variants_to_fragment,
    generate_tumor_normal_variants,
    sample_target_fragment,
    select_target_regions,
)


def simulate_tumor_normal(
    *,
    reference: Path,
    candidate_targets: Path,
    output_dir: Path,
    pair_name: str,
    normal_sample: str,
    tumor_sample: str,
    target_count: int,
    target_width: int,
    normal_depth: float,
    tumor_depth: float,
    tumor_purity: float,
    germline_snvs: int,
    clonal_snvs: int,
    subclonal_snvs: int,
    subclone_fraction: float,
    read_length: int,
    fragment_mean: int,
    fragment_sd: int,
    seed: int,
) -> dict[str, object]:
    if normal_depth <= 0 or tumor_depth <= 0:
        raise ValueError('normal and tumor depths must be positive')
    if not 0 <= tumor_purity <= 1:
        raise ValueError('tumor purity must be between 0 and 1')
    for name in (pair_name, normal_sample, tumor_sample):
        validate_sample_name(name)

    sequences = [(name, sequence) for name, _, sequence in read_fasta(reference)]
    rng = random.Random(seed)
    targets = select_target_regions(
        rng,
        read_bed(candidate_targets),
        sequences,
        target_count,
        target_width,
    )
    variants = generate_tumor_normal_variants(
        rng,
        sequences,
        targets,
        germline_snvs,
        clonal_snvs,
        subclonal_snvs,
        subclone_fraction,
    )
    prepare_output_directory(output_dir)
    raw_dir = output_dir / 'raw'
    callable_bases = sum(target.end - target.start for target in targets)
    illumina = TECHNOLOGIES['illumina-pe']
    truth_rows: list[list[object]] = []
    output_files = ['raw/targets.bed']

    samples = (
        (normal_sample, 'normal', normal_depth),
        (tumor_sample, 'tumor', tumor_depth),
    )
    for sample, role, depth in samples:
        read_pairs = max(1, math.ceil(callable_bases * depth / (2 * read_length)))
        r1_path = raw_dir / f'{sample}_R1.fastq.gz'
        r2_path = raw_dir / f'{sample}_R2.fastq.gz'
        output_files.extend(
            [
                str(r1_path.relative_to(output_dir)),
                str(r2_path.relative_to(output_dir)),
            ]
        )
        with ExitStack() as stack:
            r1 = stack.enter_context(fastq_writer(r1_path))
            r2 = stack.enter_context(fastq_writer(r2_path))
            for ordinal in range(1, read_pairs + 1):
                fragment_length = sample_positive_normal(
                    rng,
                    fragment_mean,
                    fragment_sd,
                    read_length * 2,
                    target_width,
                )
                fragment = sample_target_fragment(rng, sequences, targets, fragment_length)
                haplotype = rng.randrange(2)
                if role == 'normal' or rng.random() >= tumor_purity:
                    cellular_source = 'normal'
                    clone = 'normal'
                    subclone = False
                else:
                    cellular_source = 'tumor'
                    subclone = rng.random() < subclone_fraction
                    clone = 'tumor_subclone' if subclone else 'tumor_ancestral'
                template, applied = apply_variants_to_fragment(
                    fragment,
                    variants,
                    haplotype=haplotype,
                    cellular_source=cellular_source,
                    subclone=subclone,
                )
                if rng.random() < 0.5:
                    template = reverse_complement(template)
                    strand = '-'
                else:
                    strand = '+'
                read1, s1, i1, d1 = introduce_errors(template[:read_length], rng, illumina)
                read2, s2, i2, d2 = introduce_errors(reverse_complement(template[-read_length:]), rng, illumina)
                read_id = f'{sample}:{ordinal:09d}'
                write_fastq_record(r1, f'{read_id}/1', read1, 'I' * len(read1))
                write_fastq_record(r2, f'{read_id}/2', read2, 'I' * len(read2))
                truth_rows.append(
                    [
                        read_id,
                        sample,
                        role,
                        cellular_source,
                        clone,
                        haplotype,
                        fragment.target.contig,
                        fragment.start,
                        fragment.end,
                        strand,
                        ','.join(applied) if applied else '.',
                        s1 + s2,
                        i1 + i2,
                        d1 + d2,
                    ]
                )

    write_rows(
        raw_dir / 'targets.bed',
        ([target.contig, target.start, target.end, target.name] for target in targets),
    )
    _write_vcf(
        output_dir / 'truth' / 'germline.vcf',
        [variant for variant in variants if variant.category == 'germline'],
        normal_sample,
        tumor_sample,
        tumor_purity,
    )
    _write_vcf(
        output_dir / 'truth' / 'somatic.vcf',
        [variant for variant in variants if variant.category != 'germline'],
        normal_sample,
        tumor_sample,
        tumor_purity,
    )
    write_tsv(
        output_dir / 'truth' / 'expected_vaf.tsv',
        [
            'variant_id',
            'category',
            'contig',
            'position_1based',
            'normal_vaf',
            'tumor_vaf',
            'cancer_cell_fraction',
        ],
        (
            [
                variant.identifier,
                variant.category,
                variant.contig,
                variant.position + 1,
                f'{variant.normal_vaf():.6f}',
                f'{variant.tumor_vaf(tumor_purity):.6f}',
                f'{variant.cancer_cell_fraction:.6f}',
            ]
            for variant in variants
        ),
    )
    write_tsv(
        output_dir / 'truth' / 'clones.tsv',
        ['sample', 'clone', 'fraction'],
        [
            [normal_sample, 'normal', '1.000000'],
            [tumor_sample, 'normal_contamination', f'{1 - tumor_purity:.6f}'],
            [tumor_sample, 'tumor_ancestral', f'{tumor_purity * (1 - subclone_fraction):.6f}'],
            [tumor_sample, 'tumor_subclone', f'{tumor_purity * subclone_fraction:.6f}'],
        ],
    )
    write_tsv(
        output_dir / 'truth' / 'copy_number.tsv',
        [
            'contig',
            'start_0based',
            'end_0based_exclusive',
            'normal_total_cn',
            'tumor_total_cn',
            'tumor_minor_cn',
            'loh',
        ],
        ([target.contig, target.start, target.end, 2, 2, 1, 0] for target in targets),
    )
    write_tsv(
        output_dir / 'truth' / 'read_origins.tsv',
        [
            'read_id',
            'sample',
            'role',
            'cellular_source',
            'clone',
            'haplotype',
            'contig',
            'start_0based',
            'end_0based_exclusive',
            'strand',
            'variant_ids',
            'substitutions',
            'insertions',
            'deletions',
        ],
        truth_rows,
    )
    metadata: dict[str, object] = {
        'schema_version': 1,
        'assay': 'tumor-normal-targeted-dna',
        'pair_name': pair_name,
        'normal_sample': normal_sample,
        'tumor_sample': tumor_sample,
        'reference': str(reference.resolve()),
        'candidate_targets': str(candidate_targets.resolve()),
        'seed': seed,
        'output_files': output_files,
        'truth_files': [
            'truth/germline.vcf',
            'truth/somatic.vcf',
            'truth/expected_vaf.tsv',
            'truth/clones.tsv',
            'truth/copy_number.tsv',
            'truth/read_origins.tsv',
        ],
        'parameters': {
            'target_count': target_count,
            'target_width': target_width,
            'callable_bases': callable_bases,
            'normal_depth': normal_depth,
            'tumor_depth': tumor_depth,
            'tumor_purity': tumor_purity,
            'germline_snvs': germline_snvs,
            'clonal_snvs': clonal_snvs,
            'subclonal_snvs': subclonal_snvs,
            'subclone_fraction': subclone_fraction,
            'read_length': read_length,
            'fragment_mean': fragment_mean,
            'fragment_sd': fragment_sd,
            'copy_number_model': 'diploid copy-neutral',
        },
    }
    write_json(output_dir / 'manifest.json', metadata)
    logger.success(
        'Generated tumor-normal pair {} with purity {:.3f} in {}',
        pair_name,
        tumor_purity,
        output_dir,
    )
    return metadata


def _write_vcf(
    path: Path,
    variants: list[TumorNormalVariant],
    normal_sample: str,
    tumor_sample: str,
    tumor_purity: float,
) -> None:
    rows: list[list[object]] = [
        ['##fileformat=VCFv4.3'],
        ['##source=bioflow-sim'],
        ['##INFO=<ID=CATEGORY,Number=1,Type=String,Description="Simulation category">'],
        ['##INFO=<ID=CCF,Number=1,Type=Float,Description="Cancer cell fraction">'],
        ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'],
        ['##FORMAT=<ID=AF,Number=1,Type=Float,Description="Expected allele fraction">'],
        ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', normal_sample, tumor_sample],
    ]
    for variant in variants:
        normal_gt = '0/1' if variant.category == 'germline' else '0/0'
        tumor_gt = '0/1'
        rows.append(
            [
                variant.contig,
                variant.position + 1,
                variant.identifier,
                variant.reference,
                variant.alternate,
                '.',
                'PASS',
                f'CATEGORY={variant.category};CCF={variant.cancer_cell_fraction:.6f}',
                'GT:AF',
                f'{normal_gt}:{variant.normal_vaf():.6f}',
                f'{tumor_gt}:{variant.tumor_vaf(tumor_purity):.6f}',
            ]
        )
    write_lines(path, ('\t'.join(str(field) for field in row) for row in rows))
