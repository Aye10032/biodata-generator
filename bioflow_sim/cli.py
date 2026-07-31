import sys
from pathlib import Path

import click
from loguru import logger

from . import __version__
from .simulators.bulk_rna import (
    BULK_RNA_LAYOUTS,
    BULK_RNA_STRANDEDNESS,
    simulate_bulk_rna,
)
from .simulators.chromatin import CHROMATIN_ASSAYS, simulate_chromatin
from .simulators.dna import DNA_TECHNOLOGY_NAMES, simulate_dna
from .simulators.hic import HIC_ENZYMES, simulate_hic
from .simulators.methylation import METHYLATION_PROTOCOLS, simulate_methylation
from .simulators.scrna import SCRNA_PROTOCOLS, simulate_scrna


def configure_logging(verbose: bool) -> None:
    logger.remove()
    logger.add(
        sys.stderr,
        level='DEBUG' if verbose else 'INFO',
        format='<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | {message}',
        colorize=sys.stderr.isatty(),
    )


@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(__version__)
@click.option('-v', '--verbose', is_flag=True, help='Enable debug logging.')
def cli(verbose: bool) -> None:
    """Generate deterministic sequencing fixtures from reference sequences."""
    configure_logging(verbose)


@cli.command('dna')
@click.option(
    '--reference',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    required=True,
    help='Genome FASTA, optionally gzip-compressed.',
)
@click.option('--output-dir', type=click.Path(path_type=Path, file_okay=False), required=True)
@click.option('--sample', default='S1', show_default=True)
@click.option(
    '--technology',
    'technology_name',
    type=click.Choice(DNA_TECHNOLOGY_NAMES),
    required=True,
)
@click.option('--reads', type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option('--seed', type=int, default=114514, show_default=True)
@click.option('--read-length', type=click.IntRange(min=20), default=None)
@click.option('--fragment-mean', type=click.IntRange(min=40), default=350, show_default=True)
@click.option('--fragment-sd', type=click.IntRange(min=0), default=35, show_default=True)
@click.option('--long-read-mean', type=click.IntRange(min=100), default=None)
@click.option('--long-read-sd', type=click.IntRange(min=0), default=None)
@click.option('--snvs', type=click.IntRange(min=0), default=0, show_default=True)
def dna_command(**kwargs: object) -> None:
    """Simulate Illumina, PacBio, or Oxford Nanopore genomic reads."""
    try:
        simulate_dna(**kwargs)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error


@cli.command('bulk-rna')
@click.option(
    '--transcripts-path',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    required=True,
)
@click.option(
    '--annotation-path',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    default=None,
)
@click.option('--output-dir', type=click.Path(path_type=Path, file_okay=False), required=True)
@click.option('--sample-prefix', default='RNA', show_default=True)
@click.option('--samples-per-group', type=click.IntRange(min=1), default=2, show_default=True)
@click.option('--reads-per-sample', type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option('--layout', type=click.Choice(BULK_RNA_LAYOUTS), default='pe', show_default=True)
@click.option(
    '--strandedness',
    type=click.Choice(BULK_RNA_STRANDEDNESS),
    default='unstranded',
    show_default=True,
)
@click.option('--read-length', type=click.IntRange(min=20), default=150, show_default=True)
@click.option('--fragment-mean', type=click.IntRange(min=40), default=350, show_default=True)
@click.option('--fragment-sd', type=click.IntRange(min=0), default=35, show_default=True)
@click.option('--fold-change', type=click.FloatRange(min=0.01), default=4.0, show_default=True)
@click.option('--seed', type=int, default=114514, show_default=True)
def bulk_rna_command(**kwargs: object) -> None:
    """Simulate control/treatment bulk RNA-seq samples."""
    try:
        simulate_bulk_rna(**kwargs)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error


@cli.command('chromatin')
@click.option(
    '--reference',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    required=True,
)
@click.option('--output-dir', type=click.Path(path_type=Path, file_okay=False), required=True)
@click.option('--sample', default='EP1', show_default=True)
@click.option('--assay', type=click.Choice(CHROMATIN_ASSAYS), required=True)
@click.option('--reads', type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option('--peaks', type=click.IntRange(min=1), default=20, show_default=True)
@click.option('--peak-width', type=click.IntRange(min=10), default=500, show_default=True)
@click.option('--enrichment', type=click.FloatRange(min=0, max=1), default=0.8, show_default=True)
@click.option('--read-length', type=click.IntRange(min=20), default=75, show_default=True)
@click.option('--fragment-mean', type=click.IntRange(min=40), default=250, show_default=True)
@click.option('--fragment-sd', type=click.IntRange(min=0), default=50, show_default=True)
@click.option('--seed', type=int, default=114514, show_default=True)
def chromatin_command(**kwargs: object) -> None:
    """Simulate ATAC-seq, ChIP-seq, or CUT&Tag paired reads."""
    try:
        simulate_chromatin(**kwargs)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error


@cli.command('methylation')
@click.option(
    '--reference',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    required=True,
)
@click.option('--output-dir', type=click.Path(path_type=Path, file_okay=False), required=True)
@click.option('--sample', default='METH1', show_default=True)
@click.option('--protocol', type=click.Choice(METHYLATION_PROTOCOLS), required=True)
@click.option('--reads', type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option('--sites', type=click.IntRange(min=1), default=500, show_default=True)
@click.option('--methylation-rate', type=click.FloatRange(min=0, max=1), default=0.7, show_default=True)
@click.option('--conversion-rate', type=click.FloatRange(min=0, max=1), default=0.99, show_default=True)
@click.option('--read-length', type=click.IntRange(min=20), default=100, show_default=True)
@click.option('--fragment-mean', type=click.IntRange(min=40), default=300, show_default=True)
@click.option('--fragment-sd', type=click.IntRange(min=0), default=40, show_default=True)
@click.option('--seed', type=int, default=114514, show_default=True)
def methylation_command(**kwargs: object) -> None:
    """Simulate WGBS or EM-seq paired reads and CpG truth."""
    try:
        simulate_methylation(**kwargs)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error


@cli.command('hic')
@click.option(
    '--reference',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    required=True,
)
@click.option('--output-dir', type=click.Path(path_type=Path, file_okay=False), required=True)
@click.option('--sample', default='HIC1', show_default=True)
@click.option('--enzyme', 'enzyme_name', type=click.Choice(HIC_ENZYMES), default='mboi', show_default=True)
@click.option('--reads', type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option('--read-length', type=click.IntRange(min=20), default=100, show_default=True)
@click.option('--intra-rate', type=click.FloatRange(min=0, max=1), default=0.9, show_default=True)
@click.option('--mean-distance', type=click.IntRange(min=1), default=50_000, show_default=True)
@click.option('--seed', type=int, default=114514, show_default=True)
def hic_command(**kwargs: object) -> None:
    """Simulate restriction-enzyme Hi-C reads and contact truth."""
    try:
        simulate_hic(**kwargs)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error


@cli.command('scrna')
@click.option(
    '--transcripts-path',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    required=True,
    help='Transcript FASTA, optionally gzip-compressed.',
)
@click.option(
    '--annotation-path',
    type=click.Path(path_type=Path, exists=True, dir_okay=False, readable=True),
    default=None,
    help='Matching GTF used to produce a gene-by-cell matrix.',
)
@click.option('--output-dir', type=click.Path(path_type=Path, file_okay=False), required=True)
@click.option('--sample', default='SC1', show_default=True)
@click.option(
    '--protocol',
    type=click.Choice(SCRNA_PROTOCOLS),
    required=True,
)
@click.option('--cells', type=click.IntRange(min=1), default=8, show_default=True)
@click.option('--reads-per-cell', type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option('--seed', type=int, default=114514, show_default=True)
@click.option('--cdna-read-length', type=click.IntRange(min=20), default=90, show_default=True)
@click.option('--smartseq-read-length', type=click.IntRange(min=20), default=150, show_default=True)
@click.option('--fragment-mean', type=click.IntRange(min=40), default=400, show_default=True)
@click.option('--fragment-sd', type=click.IntRange(min=0), default=50, show_default=True)
def scrna_command(**kwargs: object) -> None:
    """Simulate Smart-seq2 paired FASTQ or 10x-style barcode/UMI reads."""
    try:
        simulate_scrna(**kwargs)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error


if __name__ == '__main__':
    cli()
