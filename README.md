# bioflow-sim

Reference-driven, deterministic sequencing fixtures for BioFlowAI. See
[`PLAN.md`](PLAN.md) for the technology matrix and delivery phases.

## Environment

The project requires Python 3.13 and uses `uv`:

```bash
uv sync --dev
uv run bioflow-sim --help
```

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

Add simulated SNVs and a truth VCF with:

```bash
uv run bioflow-sim dna \
  --reference ../Yeast/reference/GCF_000146045.2_R64_genomic.fna \
  --output-dir ../Yeast/wgs/smoke_illumina_variants \
  --technology illumina-pe \
  --reads 1000 \
  --snvs 20
```

## Bulk RNA sequencing

Generate control and treatment replicates with differential expression:

```bash
uv run bioflow-sim bulk-rna \
  --transcripts-path ../Homo_chr21/reference/transcripts.fa \
  --annotation-path ../Homo_chr21/reference/genes.gtf \
  --output-dir ../Homo_chr21/rna-seq/smoke_bulk \
  --samples-per-group 2 \
  --reads-per-sample 1000 \
  --layout pe \
  --strandedness forward \
  --read-length 100 \
  --fold-change 4
```

The sample sheet is written to `raw/samples.tsv`. True transcript weights,
observed read counts, and per-read origins are written under `truth/`.

## Chromatin assays

The `chromatin` command accepts `atac`, `chip`, or `cuttag`:

```bash
uv run bioflow-sim chromatin \
  --reference ../Yeast/reference/GCF_000146045.2_R64_genomic.fna \
  --output-dir ../Yeast/atac-seq/smoke \
  --assay atac \
  --reads 1000 \
  --peaks 20 \
  --peak-width 500 \
  --enrichment 0.8
```

The FASTQ and library metadata are under `raw/`; true enriched regions are
written as `truth/peaks.bed`.

## DNA methylation

```bash
uv run bioflow-sim methylation \
  --reference ../Homo_chr21/reference/genome.fa \
  --output-dir ../Homo_chr21/wgbs/smoke \
  --protocol wgbs \
  --reads 1000 \
  --sites 500 \
  --methylation-rate 0.7 \
  --conversion-rate 0.99
```

Use `--protocol emseq` for the EM-seq fixture. Selected CpG states are written
to `truth/methylation.tsv`; unlisted cytosines are treated as methylated by the
current small-fixture model.

## Hi-C

```bash
uv run bioflow-sim hic \
  --reference ../Yeast/reference/GCF_000146045.2_R64_genomic.fna \
  --output-dir ../Yeast/Hi-C/smoke \
  --enzyme hindiii \
  --reads 1000 \
  --read-length 100 \
  --intra-rate 0.9 \
  --mean-distance 50000
```

Supported enzymes are MboI, DpnII, and HindIII. The derived enzyme cut-site
list is tool-facing input at `raw/restriction_sites.tsv`; true ligation contacts
are in `truth/contacts.tsv`.

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
  generators/
    contacts.py
    expression.py
    methylation.py
    random_values.py
    read_models.py
    regions.py
    sequences.py
    variants.py
  simulators/
    bulk_rna.py
    chromatin.py
    dna.py
    hic.py
    methylation.py
    scrna.py
```

`generators/` contains data classes and pure value-generation functions, such
as platform error profiles, random lengths, sequence changes, barcodes, and
expression weights. It does not assemble output directories.

`simulators/` contains complete assay scenarios called by the CLI. These
functions combine reference input, generators, FASTQ writing, truth generation,
and manifests. `core/` is restricted to shared I/O, formats, and validation.

## Tests

```bash
uv run pytest
```

Tests cover deterministic Illumina output, long-read profiles, SNV truth,
bulk RNA, chromatin assays, methylation, Hi-C, the 10x barcode/UMI layout,
expression-matrix output, and per-cell Smart-seq2 files.
