from __future__ import annotations

import gzip
import json
from pathlib import Path

from click.testing import CliRunner

from bioflow_sim.cli import cli


def write_reference(path: Path) -> None:
    path.write_text(
        ">chrA synthetic\n" + "ACGT" * 1_000 + "\n"
        ">tx1 gene=GENE1\n" + "ACGT" * 200 + "\n"
        ">tx2 gene=GENE2\n" + "TGCA" * 180 + "\n",
        encoding="utf-8",
    )


def write_transcripts(path: Path) -> None:
    path.write_text(
        ">tx1 gene=GENE1\n" + "ACGT" * 200 + "\n"
        ">tx2 gene=GENE2\n" + "TGCA" * 180 + "\n",
        encoding="utf-8",
    )


def write_annotation(path: Path) -> None:
    path.write_text(
        'chrA\ttest\ttranscript\t1\t800\t.\t+\t.\t'
        'gene_id "GENE1"; transcript_id "tx1";\n'
        'chrA\ttest\ttranscript\t1\t720\t.\t+\t.\t'
        'gene_id "GENE2"; transcript_id "tx2";\n',
        encoding="utf-8",
    )


def fastq_records(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        lines = [line.rstrip("\n") for line in handle]
    assert len(lines) % 4 == 0
    records = [lines[index : index + 4] for index in range(0, len(lines), 4)]
    for name, sequence, plus, quality in records:
        assert name.startswith("@")
        assert plus == "+"
        assert len(sequence) == len(quality)
    return records


def test_illumina_pe_is_reproducible(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    write_reference(reference)
    runner = CliRunner()
    arguments = [
        "dna",
        "--reference",
        str(reference),
        "--technology",
        "illumina-pe",
        "--reads",
        "12",
        "--read-length",
        "75",
        "--seed",
        "42",
    ]
    first = tmp_path / "first"
    second = tmp_path / "second"
    result = runner.invoke(cli, [*arguments, "--output-dir", str(first)])
    assert result.exit_code == 0, result.output
    result = runner.invoke(cli, [*arguments, "--output-dir", str(second)])
    assert result.exit_code == 0, result.output
    assert (first / "raw/S1_R1.fastq.gz").read_bytes() == (
        second / "raw/S1_R1.fastq.gz"
    ).read_bytes()
    assert len(fastq_records(first / "raw/S1_R1.fastq.gz")) == 12
    assert len(fastq_records(first / "raw/S1_R2.fastq.gz")) == 12


def test_long_read_technologies(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    write_reference(reference)
    runner = CliRunner()
    for technology in ("pacbio-clr", "pacbio-hifi", "ont"):
        output = tmp_path / technology
        result = runner.invoke(
            cli,
            [
                "dna",
                "--reference",
                str(reference),
                "--output-dir",
                str(output),
                "--technology",
                technology,
                "--reads",
                "3",
                "--long-read-mean",
                "500",
                "--long-read-sd",
                "20",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(fastq_records(output / "raw/S1.fastq.gz")) == 3


def test_illumina_single_end_length(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    write_reference(reference)
    output = tmp_path / "single"
    result = CliRunner().invoke(
        cli,
        [
            "dna",
            "--reference",
            str(reference),
            "--output-dir",
            str(output),
            "--technology",
            "illumina-se",
            "--reads",
            "4",
            "--read-length",
            "75",
        ],
    )
    assert result.exit_code == 0, result.output
    assert {len(record[1]) for record in fastq_records(output / "raw/S1.fastq.gz")} == {75}


def test_10x_layout_and_truth_matrix(tmp_path: Path) -> None:
    reference = tmp_path / "transcripts.fa"
    annotation = tmp_path / "genes.gtf"
    write_transcripts(reference)
    write_annotation(annotation)
    output = tmp_path / "tenx"
    result = CliRunner().invoke(
        cli,
        [
            "scrna",
            "--transcripts-path",
            str(reference),
            "--output-dir",
            str(output),
            "--annotation-path",
            str(annotation),
            "--protocol",
            "10x-3prime",
            "--cells",
            "3",
            "--reads-per-cell",
            "7",
            "--seed",
            "9",
        ],
    )
    assert result.exit_code == 0, result.output
    r1 = fastq_records(output / "raw/SC1_R1.fastq.gz")
    r2 = fastq_records(output / "raw/SC1_R2.fastq.gz")
    assert len(r1) == len(r2) == 21
    assert {len(record[1]) for record in r1} == {28}
    assert {len(record[1]) for record in r2} == {90}
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["barcode_structure"] == "16 bp cell barcode + 12 bp UMI"
    assert manifest["matrix_feature_type"] == "Gene Expression"
    with (
        output / "truth/expression_matrix/features.tsv"
    ).open("r", encoding="utf-8") as handle:
        assert len(handle.readlines()) == 2
    assert (output / "truth/expression_matrix/matrix.mtx").exists()
    barcodes = (output / "raw/barcodes.tsv").read_text(encoding="ascii").splitlines()
    assert len(barcodes) == 3
    assert {len(barcode) for barcode in barcodes} == {16}
    assert not list((output / "truth").rglob("*.gz"))


def test_smartseq2_writes_one_pair_per_cell(tmp_path: Path) -> None:
    reference = tmp_path / "transcripts.fa"
    write_transcripts(reference)
    output = tmp_path / "smartseq"
    result = CliRunner().invoke(
        cli,
        [
            "scrna",
            "--transcripts-path",
            str(reference),
            "--output-dir",
            str(output),
            "--protocol",
            "smartseq2",
            "--cells",
            "2",
            "--reads-per-cell",
            "5",
            "--smartseq-read-length",
            "60",
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(list((output / "raw").glob("*_R1.fastq.gz"))) == 2
    for path in (output / "raw").glob("*_R1.fastq.gz"):
        assert len(fastq_records(path)) == 5


def test_unsafe_sample_name_is_rejected(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    write_reference(reference)
    result = CliRunner().invoke(
        cli,
        [
            "dna",
            "--reference",
            str(reference),
            "--output-dir",
            str(tmp_path / "output"),
            "--sample",
            "../escape",
            "--technology",
            "illumina-se",
            "--reads",
            "1",
        ],
    )
    assert result.exit_code == 1
    assert "sample must start" in result.output
