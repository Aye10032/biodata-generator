import sys
from pathlib import Path

import click
from loguru import logger

from . import __version__
from .simulators.dna import DNA_TECHNOLOGY_NAMES, simulate_dna
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
def dna_command(**kwargs: object) -> None:
    """Simulate Illumina, PacBio, or Oxford Nanopore genomic reads."""
    try:
        simulate_dna(**kwargs)
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
