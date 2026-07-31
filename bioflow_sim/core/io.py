import csv
import gzip
import io
import json
import re
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

SAFE_NAME = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')


def _open_text(path: Path, mode: str = 'rt') -> TextIO:
    if path.suffix == '.gz':
        return gzip.open(path, mode, encoding='utf-8', newline='')
    return path.open(mode, encoding='utf-8', newline='')


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    """Return (identifier, description, sequence) records from FASTA."""
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append(_fasta_record(header, chunks))
                header = line[1:]
                chunks = []
            elif header is None:
                raise ValueError(f'{path}: sequence data appears before a FASTA header')
            else:
                chunks.append(line)

    if header is not None:
        records.append(_fasta_record(header, chunks))
    if not records:
        raise ValueError(f'{path}: no FASTA records found')
    return records


def _fasta_record(header: str, chunks: list[str]) -> tuple[str, str, str]:
    identifier = header.split(maxsplit=1)[0]
    sequence = ''.join(chunks).upper().replace('U', 'T')
    if not sequence:
        raise ValueError(f'FASTA record {identifier!r} has no sequence')
    return identifier, header, sequence


def validate_sample_name(sample: str) -> None:
    if not SAFE_NAME.fullmatch(sample):
        raise ValueError(
            "sample must start with an alphanumeric character and contain only letters, numbers, '.', '_', or '-'"
        )


def read_gtf_transcript_gene_map(path: Path) -> dict[str, str]:
    """Read transcript_id -> gene_id assignments from a GTF annotation."""
    transcript_pattern = re.compile(r'(?:^|;\s*)transcript_id "([^"]+)"')
    gene_pattern = re.compile(r'(?:^|;\s*)gene_id "([^"]+)"')
    mapping: dict[str, str] = {}
    with _open_text(path) as handle:
        for raw_line in handle:
            if raw_line.startswith('#'):
                continue
            fields = raw_line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            transcript_match = transcript_pattern.search(fields[8])
            gene_match = gene_pattern.search(fields[8])
            if transcript_match and gene_match:
                mapping[transcript_match.group(1)] = gene_match.group(1)
    if not mapping:
        raise ValueError(f'{path}: no transcript_id/gene_id pairs found in GTF')
    return mapping


def read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    """Read the first four columns of a BED file."""
    records: list[tuple[str, int, int, str]] = []
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip('\n')
            if not line or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                raise ValueError(f'{path}:{line_number}: BED record has fewer than three columns')
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(f'{path}:{line_number}: invalid BED coordinates') from error
            if start < 0 or end <= start:
                raise ValueError(f'{path}:{line_number}: invalid BED interval {start}-{end}')
            name = fields[3] if len(fields) > 3 and fields[3] else f'interval_{line_number}'
            records.append((fields[0], start, end, name))
    if not records:
        raise ValueError(f'{path}: no BED records found')
    return records


@contextmanager
def fastq_writer(path: Path) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Blank embedded filename and mtime=0 make the compressed bytes reproducible.
    with (
        path.open('wb') as raw_handle,
        gzip.GzipFile(filename='', mode='wb', fileobj=raw_handle, compresslevel=6, mtime=0) as gzip_handle,
        io.TextIOWrapper(gzip_handle, encoding='ascii', newline='\n', write_through=True) as text_handle,
    ):
        yield text_handle


def write_fastq_record(handle: TextIO, name: str, sequence: str, quality: str) -> None:
    if len(sequence) != len(quality):
        raise ValueError(f'{name}: sequence and quality lengths differ')
    handle.write(f'@{name}\n{sequence}\n+\n{quality}\n')


def write_tsv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def write_rows(
    path: Path,
    rows: Iterable[Iterable[object]],
    delimiter: str = '\t',
) -> None:
    """Write delimiter-separated rows without a header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator='\n')
        writer.writerows(rows)


def write_lines(path: Path, lines: Iterable[str]) -> None:
    """Write text lines verbatim with normalized line endings."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        for line in lines:
            handle.write(line)
            handle.write('\n')


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')


def prepare_output_directory(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f'output directory is not empty: {path}')
    path.mkdir(parents=True, exist_ok=True)
