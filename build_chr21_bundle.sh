#!/usr/bin/env bash

set -euo pipefail

GENOME_SOURCE="GRCh38.primary_assembly.genome.fa"
GTF_SOURCE="gencode.v50.basic.annotation.gtf.gz"
GFF3_SOURCE="gencode.v50.basic.annotation.gff3.gz"

OUTPUT_DIR="Homo_sapiens_GRCh38p14_chr21_GENCODE_v50"
REGION="chr21"

for command_name in samtools gffread gawk gzip sha256sum; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "Error: required command not found: ${command_name}" >&2
        exit 1
    fi
done

for input_file in \
    "${GENOME_SOURCE}" \
    "${GTF_SOURCE}" \
    "${GFF3_SOURCE}"
do
    if [[ ! -f "${input_file}" ]]; then
        echo "Error: input file not found: ${input_file}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

echo "[1/8] Indexing source genome..."

if [[ ! -f "${GENOME_SOURCE}.fai" ]]; then
    samtools faidx "${GENOME_SOURCE}"
fi

if ! cut -f1 "${GENOME_SOURCE}.fai" | grep -Fxq "${REGION}"; then
    echo "Error: ${REGION} is not present in ${GENOME_SOURCE}" >&2
    exit 1
fi

echo "[2/8] Extracting ${REGION} genome..."

samtools faidx "${GENOME_SOURCE}" "${REGION}" \
    > "${OUTPUT_DIR}/genome.fa"

samtools faidx "${OUTPUT_DIR}/genome.fa"

echo "[3/8] Filtering GTF..."

gzip -cd "${GTF_SOURCE}" |
    gawk -v region="${REGION}" '
        BEGIN {
            FS = OFS = "\t"
        }

        /^#/ {
            print
            next
        }

        $1 == region {
            print
        }
    ' > "${OUTPUT_DIR}/genes.gtf"

echo "[4/8] Filtering GFF3..."

gzip -cd "${GFF3_SOURCE}" |
    gawk -v region="${REGION}" '
        BEGIN {
            FS = OFS = "\t"
        }

        /^##FASTA/ {
            exit
        }

        /^##sequence-region/ {
            if ($2 == region) {
                print
            }
            next
        }

        /^#/ {
            print
            next
        }

        $1 == region {
            print
        }
    ' > "${OUTPUT_DIR}/genes.gff3"

echo "[5/8] Generating gene-level BED6..."

gawk '
    BEGIN {
        FS = OFS = "\t"
    }

    /^#/ {
        next
    }

    $3 != "gene" {
        next
    }

    {
        gene_id = ""
        gene_name = ""

        if (match($9, /gene_id "[^"]+"/)) {
            gene_id = substr($9, RSTART + 9, RLENGTH - 10)
        }

        if (match($9, /gene_name "[^"]+"/)) {
            gene_name = substr($9, RSTART + 11, RLENGTH - 12)
        }

        name = gene_name

        if (name == "") {
            name = gene_id
        }

        if (name == "") {
            name = "."
        }

        # GTF: 1-based closed
        # BED: 0-based half-open
        print $1, $4 - 1, $5, name, 0, $7
    }
' "${OUTPUT_DIR}/genes.gtf" |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    > "${OUTPUT_DIR}/genes.bed"

echo "[6/8] Generating transcript and protein FASTA files..."

gffread "${OUTPUT_DIR}/genes.gtf" \
    -g "${OUTPUT_DIR}/genome.fa" \
    -w "${OUTPUT_DIR}/transcripts.fa" \
    -y "${OUTPUT_DIR}/proteins.fa"

echo "[7/8] Running validation checks..."

GENOME_REGIONS="$(
    grep '^>' "${OUTPUT_DIR}/genome.fa" |
        sed 's/^>//; s/[[:space:]].*$//' |
        sort -u
)"

GTF_REGIONS="$(
    gawk -F '\t' '!/^#/ {print $1}' "${OUTPUT_DIR}/genes.gtf" |
        sort -u
)"

GFF3_REGIONS="$(
    gawk -F '\t' '!/^#/ {print $1}' "${OUTPUT_DIR}/genes.gff3" |
        sort -u
)"

if [[ "${GENOME_REGIONS}" != "${REGION}" ]]; then
    echo "Error: unexpected genome sequence names:" >&2
    printf '%s\n' "${GENOME_REGIONS}" >&2
    exit 1
fi

if [[ "${GTF_REGIONS}" != "${REGION}" ]]; then
    echo "Error: unexpected GTF sequence names:" >&2
    printf '%s\n' "${GTF_REGIONS}" >&2
    exit 1
fi

if [[ "${GFF3_REGIONS}" != "${REGION}" ]]; then
    echo "Error: unexpected GFF3 sequence names:" >&2
    printf '%s\n' "${GFF3_REGIONS}" >&2
    exit 1
fi

for output_file in \
    genome.fa \
    genes.gtf \
    genes.gff3 \
    genes.bed \
    transcripts.fa \
    proteins.fa
do
    if [[ ! -s "${OUTPUT_DIR}/${output_file}" ]]; then
        echo "Error: output is missing or empty: ${output_file}" >&2
        exit 1
    fi
done

GENOME_LENGTH="$(
    cut -f2 "${OUTPUT_DIR}/genome.fa.fai"
)"

GENE_COUNT="$(
    gawk -F '\t' '$3 == "gene" {count++} END {print count + 0}' \
        "${OUTPUT_DIR}/genes.gtf"
)"

TRANSCRIPT_COUNT="$(
    grep -c '^>' "${OUTPUT_DIR}/transcripts.fa"
)"

PROTEIN_COUNT="$(
    grep -c '^>' "${OUTPUT_DIR}/proteins.fa"
)"

echo "[8/8] Writing metadata and checksums..."

cat > "${OUTPUT_DIR}/about.md" <<ABOUT
# Homo sapiens GRCh38.p14 chr21 reference bundle

## Reference definition

- Species: Homo sapiens
- Assembly: GRCh38.p14
- Annotation release: GENCODE Human Release 50
- Retained sequence region: chr21
- Sequence naming convention: GENCODE/UCSC-style \`chr21\`
- Annotation set: GENCODE basic chromosome annotation

## Source files

- Genome: \`${GENOME_SOURCE}\`
- GTF annotation: \`${GTF_SOURCE}\`
- GFF3 annotation: \`${GFF3_SOURCE}\`

## Generated files

- \`genome.fa\`: chr21 extracted from the GENCODE primary-assembly genome
- \`genes.gtf\`: chr21 records filtered from the GENCODE basic CHR GTF
- \`genes.gff3\`: chr21 records filtered from the GENCODE basic CHR GFF3
- \`genes.bed\`: gene-level BED6 generated from GTF gene features
- \`transcripts.fa\`: spliced transcript sequences generated by gffread
- \`proteins.fa\`: CDS translation sequences generated by gffread

## Coordinate conventions

- GTF/GFF3: 1-based, closed intervals
- BED: 0-based, half-open intervals
- FASTA sequence name: chr21

## Summary

- chr21 length: ${GENOME_LENGTH} bp
- Annotated genes: ${GENE_COUNT}
- Transcript sequences: ${TRANSCRIPT_COUNT}
- Protein sequences: ${PROTEIN_COUNT}

## Build software

\`\`\`text
$(samtools --version | head -n 1)
$(gffread --version 2>&1 | head -n 1)
\`\`\`
ABOUT

(
    cd "${OUTPUT_DIR}"
    sha256sum \
        genome.fa \
        genes.gtf \
        genes.gff3 \
        genes.bed \
        transcripts.fa \
        proteins.fa \
        > SHA256SUMS
)

echo
echo "Build completed successfully:"
echo "  ${OUTPUT_DIR}"
echo
echo "Summary:"
echo "  chr21 length: ${GENOME_LENGTH} bp"
echo "  genes:        ${GENE_COUNT}"
echo "  transcripts:  ${TRANSCRIPT_COUNT}"
echo "  proteins:     ${PROTEIN_COUNT}"
