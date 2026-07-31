# bioflow-sim

Reference-driven, deterministic sequencing fixtures for BioFlowAI. See
[`PLAN.md`](PLAN.md) for the technology matrix and delivery phases.

## Environment

The project requires Python 3.13 and uses `uv`:

```bash
uv sync --dev
uv run bioflow-sim --help
```

Runtime dependencies are deliberately small:

- Click for command-line parsing.
- Loguru for structured human-readable logs.

## Genomic sequencing

Supported technology names:

- `illumina-se`
- `illumina-pe`
- `pacbio-clr`
- `pacbio-hifi`
- `ont`

Example Illumina paired-end Lambda smoke fixture:

```bash
uv run bioflow-sim dna \
  --reference ../Lambda/reference/GCF_000840245.1_ViralProj14204_genomic.fna \
  --output-dir ../Lambda/wgs/smoke_illumina_pe \
  --sample LAMBDA_PE \
  --technology illumina-pe \
  --reads 1000 \
  --read-length 150 \
  --fragment-mean 350 \
  --fragment-sd 35 \
  --seed 20260731
```

Example yeast Oxford Nanopore fixture:

```bash
uv run bioflow-sim dna \
  --reference ../Yeast/reference/GCF_000146045.2_R64_genomic.fna \
  --output-dir ../Yeast/wgs/smoke_ont \
  --sample YEAST_ONT \
  --technology ont \
  --reads 100 \
  --long-read-mean 12000 \
  --long-read-sd 5000 \
  --seed 20260731
```

For paired-end technologies, `--reads` means read pairs. For single-end and
long-read technologies, it means individual reads.

## Single-cell RNA sequencing

The 10x-style mode writes one R1/R2 pair for the library:

- R1: 16 bp cell barcode followed by a 12 bp UMI.
- R2: transcript-derived cDNA sequence.

```bash
uv run bioflow-sim scrna \
  --transcripts-path ../Homo_chr21/reference/transcripts.fa \
  --annotation-path ../Homo_chr21/reference/genes.gtf \
  --output-dir ../Homo_chr21/sc-rna-seq/smoke_10x \
  --sample HUMAN21_10X \
  --protocol 10x-3prime \
  --cells 8 \
  --reads-per-cell 500 \
  --seed 20260731
```

Smart-seq2 mode writes one ordinary paired-end FASTQ pair per cell:

```bash
uv run bioflow-sim scrna \
  --transcripts-path ../Yeast/reference/rna.fna \
  --annotation-path ../Yeast/reference/genomic.gtf \
  --output-dir ../Yeast/sc-rna-seq/smoke_smartseq2 \
  --sample YEAST_SS2 \
  --protocol smartseq2 \
  --cells 4 \
  --reads-per-cell 500 \
  --smartseq-read-length 100 \
  --seed 20260731
```

## Output

```text
<case>/
  manifest.json
  raw/
    *.fastq.gz
    barcodes.tsv                       # 10x only
  truth/
    reads.tsv
    cells.tsv                         # scRNA only
    expression_matrix/                # scRNA only
      features.tsv
      matrix.mtx
```

With `--annotation-path`, the Matrix Market output is gene-by-cell and at least
95% of FASTA transcript identifiers must map to GTF `transcript_id` values.
Without an annotation, the command explicitly produces a transcript-by-cell
matrix. The per-read truth always retains transcript and feature identifiers.

`raw/` contains files presented to an analysis tool as experimental inputs:
reads, 10x barcode lists, and later assay metadata such as restriction-enzyme
sites. `truth/` contains simulator-only expected results such as read origins,
true variants, cell labels, and expression matrices.

Only sequencing reads (`*.fastq.gz`) are compressed. Manifests, truth TSVs,
barcodes, features, and Matrix Market files remain plain text because the
fixtures are intentionally small and should be easy to inspect.

## Package layout

```text
bioflow_sim/
  cli.py
  core/
    io.py
  sequencing/
    dna.py
    models.py
  single_cell/
    rna.py
```

## Tests

```bash
uv run pytest
```

Tests cover deterministic Illumina output, all three long-read profiles, the
10x barcode/UMI layout, expression-matrix output, and per-cell Smart-seq2 files.
