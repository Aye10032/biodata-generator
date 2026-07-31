#!/usr/bin/env bash
# Generate command-line upstream fixtures from the compact simulated FASTQs.
# This script deliberately does not publish expression matrices or downstream
# interpretation products. Run from any directory; paths are resolved from the
# script location.

set -Eeuo pipefail
IFS=$'\n\t'

# Keep pristine descriptors for orchestration logs. Third-party tools are sent
# through prefixed stdout/stderr streams later, so their output stays visually
# separate from the script's own phase and step messages.
exec 3>&1 4>&2

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
THREADS="${THREADS:-4}"
PHASES="all"
SPECIES_FILTER="all"
JUICER_TOOLS_JAR="${JUICER_TOOLS_JAR:-}"
WORK_ROOT=""
COLOR_MODE="auto"
STEP_NUMBER=0
CURRENT_STEP=""
CURRENT_STEP_STARTED=0
CURRENT_PHASE=""
CURRENT_PHASE_STARTED=0
RUN_STARTED="$SECONDS"
TOOL_STDOUT_PID=""
TOOL_STDERR_PID=""
TOOL_LOG_MARKER=""
TOOL_STDOUT_ACK=""
TOOL_STDERR_ACK=""
TOOL_LOG_FILE=""
INTERACTIVE_LOGGING=0
LIVE_TOOL_WIDTH=100

C_RESET=""
C_BOLD=""
C_DIM=""
C_BLUE=""
C_CYAN=""
C_GREEN=""
C_YELLOW=""
C_RED=""

usage() {
  cat <<'EOF'
Usage:
  generate_upstream.sh [--phase NAME[,NAME...]] [--species NAME[,NAME...]] [--no-color]

Phases:
  reference,wgs,rna,scrna,chromatin,hic,methylation,somatic,all

Species:
  Lambda,Ecoli,Yeast,Homo_chr21,all

Environment:
  THREADS             Number of worker threads (default: 4)
  TMPDIR              Parent directory for temporary work
  JUICER_TOOLS_JAR    Required by the hic phase to create .hic files

Safety:
  The script refuses to write into an existing non-empty upstream directory.
  Hi-C is resumable and skips complete pairs/Cooler stages after interruption.
  It never removes existing dataset files.

Logging:
  Color is enabled on terminals. Use --no-color or NO_COLOR=1 to disable it.
EOF
}

die() {
  printf '%bERROR%b %s\n' "$C_RED$C_BOLD" "$C_RESET" "$*" >&4
  exit 1
}

format_duration() {
  local seconds="$1"
  if (( seconds >= 3600 )); then
    printf '%dh%02dm%02ds' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
  elif (( seconds >= 60 )); then
    printf '%dm%02ds' "$((seconds / 60))" "$((seconds % 60))"
  else
    printf '%ds' "$seconds"
  fi
}

tool_log_formatter() {
  local ack_path="$1" line=""
  while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" == "$TOOL_LOG_MARKER" ]]; then
      printf 'ok\n' > "$ack_path"
      continue
    fi
    line="${line//$'\r'/}"
    printf '%s\n' "$line" >> "$TOOL_LOG_FILE"
    if (( INTERACTIVE_LOGGING )); then
      # Keep one transient tool line directly below the current STEP. It is
      # replaced in place instead of adding thousands of terminal lines.
      printf '\r\033[2K%b│ tool │%b %.*s' "$C_DIM" "$C_RESET" "$LIVE_TOOL_WIDTH" "$line" >&4
    fi
  done
}

setup_logging() {
  if [[ "$COLOR_MODE" != "never" && -z "${NO_COLOR:-}" && "${TERM:-}" != "dumb" && -t 4 ]]; then
    C_RESET=$'\033[0m'
    C_BOLD=$'\033[1m'
    C_DIM=$'\033[2m'
    C_BLUE=$'\033[34m'
    C_CYAN=$'\033[36m'
    C_GREEN=$'\033[32m'
    C_YELLOW=$'\033[33m'
    C_RED=$'\033[31m'
  fi
  if [[ -t 4 && "${TERM:-}" != "dumb" ]]; then
    local terminal_columns
    INTERACTIVE_LOGGING=1
    terminal_columns="$(tput cols 2>/dev/null || printf '120')"
    LIVE_TOOL_WIDTH=$((terminal_columns - 12))
    (( LIVE_TOOL_WIDTH < 20 )) && LIVE_TOOL_WIDTH=20
  fi

  TOOL_LOG_MARKER="__BIOFLOW_LOG_SYNC_${BASHPID}_${RANDOM}__"
  TOOL_STDOUT_ACK="$WORK_ROOT/.stdout-ack"
  TOOL_STDERR_ACK="$WORK_ROOT/.stderr-ack"
  TOOL_LOG_FILE="$WORK_ROOT/current-step.tool.log"
  : > "$TOOL_LOG_FILE"

  # Preserve data-producing pipelines: only uncaptured tool stdout/stderr is
  # decorated. A marker/ack handshake lets each step synchronously drain both
  # formatters before its DONE line is printed.
  exec > >(tool_log_formatter "$TOOL_STDOUT_ACK")
  TOOL_STDOUT_PID="$!"
  exec 2> >(tool_log_formatter "$TOOL_STDERR_ACK")
  TOOL_STDERR_PID="$!"
}

flush_tool_logs() {
  local attempt
  rm -f -- "$TOOL_STDOUT_ACK" "$TOOL_STDERR_ACK"
  printf '%s\n' "$TOOL_LOG_MARKER" >&1
  printf '%s\n' "$TOOL_LOG_MARKER" >&2
  for ((attempt = 0; attempt < 500; attempt++)); do
    if [[ -s "$TOOL_STDOUT_ACK" && -s "$TOOL_STDERR_ACK" ]]; then
      return 0
    fi
    sleep 0.01
  done
  printf '%b! WARNING%b tool log formatter did not acknowledge flush within 5 seconds\n' \
    "$C_YELLOW$C_BOLD" "$C_RESET" >&4
}

info() {
  printf '%b●%b %s\n' "$C_BLUE" "$C_RESET" "$*" >&4
}

clear_live_tool_line() {
  if (( INTERACTIVE_LOGGING )); then
    printf '\r\033[2K' >&4
  fi
}

begin_phase() {
  CURRENT_PHASE="$1"
  CURRENT_PHASE_STARTED="$SECONDS"
  printf '\n%b━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━%b\n' "$C_BLUE$C_BOLD" "$C_RESET" >&4
  printf '%b  PHASE  %s%b\n' "$C_BLUE$C_BOLD" "$CURRENT_PHASE" "$C_RESET" >&4
  printf '%b━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━%b\n' "$C_BLUE$C_BOLD" "$C_RESET" >&4
}

end_phase() {
  local elapsed=$((SECONDS - CURRENT_PHASE_STARTED))
  printf '%b✓ PHASE COMPLETE%b  %s  %b(%s)%b\n' \
    "$C_GREEN$C_BOLD" "$C_RESET" "$CURRENT_PHASE" "$C_DIM" "$(format_duration "$elapsed")" "$C_RESET" >&4
  CURRENT_PHASE=""
}

begin_step() {
  STEP_NUMBER=$((STEP_NUMBER + 1))
  CURRENT_STEP="$1"
  CURRENT_STEP_STARTED="$SECONDS"
  clear_live_tool_line
  : > "$TOOL_LOG_FILE"
  printf '\n%b┌─ STEP %02d%b  %s\n' "$C_CYAN$C_BOLD" "$STEP_NUMBER" "$C_RESET" "$CURRENT_STEP" >&4
}

end_step() {
  local elapsed
  flush_tool_logs
  clear_live_tool_line
  elapsed=$((SECONDS - CURRENT_STEP_STARTED))
  printf '%b└─ DONE%b     %s  %b(%s)%b\n' \
    "$C_GREEN$C_BOLD" "$C_RESET" "$CURRENT_STEP" "$C_DIM" "$(format_duration "$elapsed")" "$C_RESET" >&4
  : > "$TOOL_LOG_FILE"
  CURRENT_STEP=""
}

skip_step() {
  STEP_NUMBER=$((STEP_NUMBER + 1))
  printf '%b↷ STEP %02d  SKIP%b  %s\n' "$C_YELLOW" "$STEP_NUMBER" "$C_RESET" "$1" >&4
}

run_if_missing() {
  local label="$1" output="$2"
  shift 2
  if [[ -s "$output" ]]; then
    skip_step "$label — already exists"
    return 0
  fi
  begin_step "$label"
  "$@"
  end_step
}

on_error() {
  local status="$1" line="$2" command="$3" tool_line
  trap - ERR
  shutdown_tool_logging
  clear_live_tool_line
  printf '\n%b✗ FAILED%b' "$C_RED$C_BOLD" "$C_RESET" >&4
  [[ -n "$CURRENT_PHASE" ]] && printf '  phase=%s' "$CURRENT_PHASE" >&4
  [[ -n "$CURRENT_STEP" ]] && printf '  step=%s' "$CURRENT_STEP" >&4
  printf '\n%b  line %s · exit %s · %s%b\n' "$C_RED" "$line" "$status" "$command" "$C_RESET" >&4
  if [[ -s "$TOOL_LOG_FILE" ]]; then
    printf '%b  last 20 tool lines:%b\n' "$C_YELLOW$C_BOLD" "$C_RESET" >&4
    while IFS= read -r tool_line; do
      printf '%b  │ %b%s\n' "$C_DIM" "$C_RESET" "$tool_line" >&4
    done < <(tail -n 20 "$TOOL_LOG_FILE")
  fi
  exit "$status"
}

shutdown_tool_logging() {
  exec 1>&- 2>&-
  if [[ -n "$TOOL_STDOUT_PID" ]]; then
    wait "$TOOL_STDOUT_PID" || true
    TOOL_STDOUT_PID=""
  fi
  if [[ -n "$TOOL_STDERR_PID" ]]; then
    wait "$TOOL_STDERR_PID" || true
    TOOL_STDERR_PID=""
  fi
}

cleanup() {
  local status="$?"
  trap - EXIT
  shutdown_tool_logging
  [[ -n "$WORK_ROOT" && -d "$WORK_ROOT" ]] && rm -rf -- "$WORK_ROOT"
  exit "$status"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

phase_selected() {
  [[ "$PHASES" == "all" || ",$PHASES," == *",$1,"* ]]
}

species_selected() {
  [[ "$SPECIES_FILTER" == "all" || ",$SPECIES_FILTER," == *",$1,"* ]]
}

validate_selection() {
  local item
  IFS=',' read -r -a selected_phases <<< "$PHASES"
  for item in "${selected_phases[@]}"; do
    case "$item" in
      all|reference|wgs|rna|scrna|chromatin|hic|methylation|somatic) ;;
      *) die "unknown phase: $item" ;;
    esac
  done
  IFS=',' read -r -a selected_species <<< "$SPECIES_FILTER"
  for item in "${selected_species[@]}"; do
    case "$item" in
      all|Lambda|Ecoli|Yeast|Homo_chr21) ;;
      *) die "unknown species: $item" ;;
    esac
  done
}

prepare_upstream() {
  local directory="$1"
  if [[ -d "$directory" && -n "$(find "$directory" -mindepth 1 -print -quit)" ]]; then
    die "refusing to overwrite non-empty directory: $directory"
  fi
  mkdir -p "$directory"
}

all_files_nonempty() {
  local path
  for path in "$@"; do
    [[ -s "$path" ]] || return 1
  done
}

reference_path() {
  case "$1" in
    Lambda) printf '%s\n' "$REPO_ROOT/Lambda/reference/GCF_000840245.1_ViralProj14204_genomic.fna" ;;
    Ecoli) printf '%s\n' "$REPO_ROOT/Ecoli/reference/GCF_000005845.2_ASM584v2_genomic.fna" ;;
    Yeast) printf '%s\n' "$REPO_ROOT/Yeast/reference/GCF_000146045.2_R64_genomic.fna" ;;
    Homo_chr21) printf '%s\n' "$REPO_ROOT/Homo_chr21/reference/genome.fa" ;;
    *) die "no reference configured for $1" ;;
  esac
}

annotation_path() {
  case "$1" in
    Ecoli) printf '%s\n' "$REPO_ROOT/Ecoli/reference/genomic.gtf" ;;
    Yeast) printf '%s\n' "$REPO_ROOT/Yeast/reference/genomic.gtf" ;;
    Homo_chr21) printf '%s\n' "$REPO_ROOT/Homo_chr21/reference/genes.gtf" ;;
    *) die "no RNA annotation configured for $1" ;;
  esac
}

bwa_prefix() {
  local reference
  reference="$(reference_path "$1")"
  printf '%s/index/bwa/%s\n' "$(dirname "$reference")" "$(basename "$reference")"
}

dictionary_path() {
  local reference="$1"
  printf '%s.dict\n' "${reference%.*}"
}

sample_prefix() {
  case "$1" in
    Lambda) printf 'LAMBDA' ;;
    Ecoli) printf 'ECOLI' ;;
    Yeast) printf 'YEAST' ;;
    Homo_chr21) printf 'H21' ;;
    *) die "no sample prefix configured for $1" ;;
  esac
}

read_group() {
  local sample="$1"
  # BWA expects backslash-escaped tab sequences in -R, not literal tab bytes.
  printf '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA' "$sample" "$sample"
}

build_reference_indexes() {
  require_cmd samtools
  require_cmd gatk
  require_cmd bwa
  require_cmd minimap2

  local species reference ref_dir index_dir prefix dictionary mmi
  for species in Lambda Ecoli Yeast Homo_chr21; do
    species_selected "$species" || continue
    reference="$(reference_path "$species")"
    ref_dir="$(dirname "$reference")"
    index_dir="$ref_dir/index"
    prefix="$(bwa_prefix "$species")"
    dictionary="$(dictionary_path "$reference")"
    mmi="$index_dir/minimap2/$(basename "$reference").mmi"
    mkdir -p "$index_dir/bwa" "$index_dir/minimap2"

    info "Reference system: $species"
    run_if_missing "$species · samtools FASTA index" "$reference.fai" samtools faidx "$reference"
    run_if_missing "$species · GATK sequence dictionary" "$dictionary" \
      gatk CreateSequenceDictionary -R "$reference" -O "$dictionary"
    run_if_missing "$species · BWA index" "$prefix.bwt" bwa index -p "$prefix" "$reference"
    run_if_missing "$species · minimap2 index" "$mmi" minimap2 -t "$THREADS" -d "$mmi" "$reference"

    cat > "$index_dir/README.txt" <<EOF
Generated by Scripts/other_scripts/generate_upstream.sh
Reference: $(basename "$reference")
FASTA index: $(basename "$reference").fai
Sequence dictionary: $(basename "$dictionary")
BWA prefix: bwa/$(basename "$reference")
minimap2 index: minimap2/$(basename "$mmi")
EOF
  done
}

align_bwa_se() {
  local species="$1" fastq="$2" sample="$3" output="$4"
  bwa mem -t "$THREADS" -R "$(read_group "$sample")" "$(bwa_prefix "$species")" "$fastq" \
    | samtools sort -@ "$THREADS" -o "$output" -
  samtools index -@ "$THREADS" "$output"
}

align_bwa_pe() {
  local species="$1" r1="$2" r2="$3" sample="$4" output="$5"
  bwa mem -t "$THREADS" -R "$(read_group "$sample")" "$(bwa_prefix "$species")" "$r1" "$r2" \
    | samtools sort -@ "$THREADS" -o "$output" -
  samtools index -@ "$THREADS" "$output"
}

call_small_variants() {
  local reference="$1" bam="$2" sample="$3" output_dir="$4"
  local vcf="$output_dir/${sample}.calls.vcf.gz"
  local bcf="$output_dir/${sample}.calls.bcf"
  bcftools mpileup -Ou -f "$reference" "$bam" \
    | bcftools call -mv -Oz -o "$vcf"
  tabix -f -p vcf "$vcf"
  bcftools view -Ob -o "$bcf" "$vcf"
  bcftools index -f -c "$bcf"
}

align_long_reads() {
  local species="$1" fastq="$2" sample="$3" preset="$4" output_dir="$5"
  local reference bam paf
  reference="$(reference_path "$species")"
  bam="$output_dir/${sample}.sorted.bam"
  paf="$output_dir/${sample}.alignments.paf.gz"

  minimap2 -t "$THREADS" -x "$preset" "$reference" "$fastq" | bgzip -c > "$paf"
  minimap2 -t "$THREADS" -a --MD -x "$preset" "$reference" "$fastq" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index -@ "$THREADS" "$bam"
}

generate_wgs() {
  require_cmd bwa
  require_cmd minimap2
  require_cmd samtools
  require_cmd bcftools
  require_cmd bgzip
  require_cmd tabix

  local species base prefix case_dir upstream sample reference bam
  for species in Lambda Ecoli Yeast Homo_chr21; do
    species_selected "$species" || continue
    base="$REPO_ROOT/$species/wgs"
    prefix="$(sample_prefix "$species")"
    reference="$(reference_path "$species")"

    case_dir="$base/illumina_se"
    upstream="$case_dir/upstream"
    sample="${prefix}_SE"
    prepare_upstream "$upstream"
    begin_step "$species · WGS · Illumina single-end alignment"
    align_bwa_se "$species" "$case_dir/${sample}.fastq.gz" "$sample" "$upstream/${sample}.sorted.bam"
    end_step

    case_dir="$base/illumina_pe"
    upstream="$case_dir/upstream"
    sample="${prefix}_PE"
    prepare_upstream "$upstream"
    begin_step "$species · WGS · Illumina paired-end alignment"
    bam="$upstream/${sample}.sorted.bam"
    align_bwa_pe "$species" "$case_dir/${sample}_R1.fastq.gz" "$case_dir/${sample}_R2.fastq.gz" "$sample" "$bam"
    end_step
    begin_step "$species · WGS · small-variant calling"
    call_small_variants "$reference" "$bam" "$sample" "$upstream"
    end_step

    case_dir="$base/pacbio_clr"
    upstream="$case_dir/upstream"
    sample="${prefix}_CLR"
    prepare_upstream "$upstream"
    begin_step "$species · WGS · PacBio CLR alignment"
    align_long_reads "$species" "$case_dir/${sample}.fastq.gz" "$sample" map-pb "$upstream"
    end_step

    case_dir="$base/pacbio_hifi"
    upstream="$case_dir/upstream"
    sample="${prefix}_HIFI"
    prepare_upstream "$upstream"
    begin_step "$species · WGS · PacBio HiFi alignment"
    align_long_reads "$species" "$case_dir/${sample}.fastq.gz" "$sample" map-hifi "$upstream"
    end_step

    case_dir="$base/ont"
    upstream="$case_dir/upstream"
    sample="${prefix}_ONT"
    prepare_upstream "$upstream"
    begin_step "$species · WGS · Oxford Nanopore alignment"
    align_long_reads "$species" "$case_dir/${sample}.fastq.gz" "$sample" map-ont "$upstream"
    end_step
  done
}

star_sa_bases() {
  case "$1" in
    Ecoli) printf '10' ;;
    Yeast) printf '10' ;;
    Homo_chr21) printf '12' ;;
    *) die "no STAR index tuning configured for $1" ;;
  esac
}

build_star_index() {
  local species="$1" index_dir="$2"
  mkdir -p "$index_dir"
  STAR --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$index_dir" \
    --genomeFastaFiles "$(reference_path "$species")" \
    --sjdbGTFfile "$(annotation_path "$species")" \
    --sjdbOverhang 99 \
    --genomeSAindexNbases "$(star_sa_bases "$species")"
}

star_align_bulk_sample() {
  local index_dir="$1" r1="$2" r2="$3" sample="$4" output_dir="$5" work="$6"
  local prefix="$work/${sample}."
  STAR --runThreadN "$THREADS" \
    --genomeDir "$index_dir" \
    --readFilesIn "$r1" "$r2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$prefix" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrRGline "ID:$sample" "SM:$sample" "PL:ILLUMINA"
  cp "$prefix"Aligned.sortedByCoord.out.bam "$output_dir/${sample}.sorted.bam"
  cp "$prefix"SJ.out.tab "$output_dir/${sample}.SJ.out.tab"
  cp "$prefix"Log.final.out "$output_dir/${sample}.Log.final.out"
  samtools index -@ "$THREADS" "$output_dir/${sample}.sorted.bam"
}

generate_rna() {
  require_cmd STAR
  require_cmd samtools

  local species case_dir upstream index_dir work sample base
  for species in Ecoli Yeast Homo_chr21; do
    species_selected "$species" || continue
    case_dir="$REPO_ROOT/$species/rna-seq/illumina_pe"
    upstream="$case_dir/upstream"
    prepare_upstream "$upstream"
    work="$WORK_ROOT/rna-$species"
    index_dir="$work/star-index"
    mkdir -p "$work"
    begin_step "$species · build temporary STAR index"
    build_star_index "$species" "$index_dir"
    end_step
    base="$(sample_prefix "$species")_RNA"
    for sample in "${base}_control_1" "${base}_treatment_1"; do
      begin_step "$species · bulk RNA-seq alignment · $sample"
      star_align_bulk_sample "$index_dir" "$case_dir/${sample}_R1.fastq.gz" "$case_dir/${sample}_R2.fastq.gz" "$sample" "$upstream" "$work"
      end_step
    done
  done
}

star_align_10x() {
  local index_dir="$1" case_dir="$2" sample="$3" output_dir="$4" work="$5"
  local prefix="$work/${sample}."
  STAR --runThreadN "$THREADS" \
    --genomeDir "$index_dir" \
    --readFilesIn "$case_dir/${sample}_R2.fastq.gz" "$case_dir/${sample}_R1.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$prefix" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM CR CY UR UY CB UB GX GN \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$case_dir/barcodes.tsv" \
    --soloFeatures Gene
  cp "$prefix"Aligned.sortedByCoord.out.bam "$output_dir/${sample}.tagged.sorted.bam"
  cp "$prefix"Log.final.out "$output_dir/${sample}.Log.final.out"
  samtools index -@ "$THREADS" "$output_dir/${sample}.tagged.sorted.bam"
}

generate_scrna() {
  require_cmd STAR
  require_cmd samtools

  local species index_dir work case_dir upstream sample barcode ordinal cell
  for species in Yeast Homo_chr21; do
    species_selected "$species" || continue
    work="$WORK_ROOT/scrna-$species"
    index_dir="$work/star-index"
    mkdir -p "$work"
    begin_step "$species · build temporary STAR index for scRNA-seq"
    build_star_index "$species" "$index_dir"
    end_step

    case_dir="$REPO_ROOT/$species/sc-rna-seq/10x_3prime"
    upstream="$case_dir/upstream"
    prepare_upstream "$upstream"
    sample="$(sample_prefix "$species")_10X"
    begin_step "$species · 10x tagged alignment · $sample"
    star_align_10x "$index_dir" "$case_dir" "$sample" "$upstream" "$work"
    end_step

    case_dir="$REPO_ROOT/$species/sc-rna-seq/smartseq2"
    upstream="$case_dir/upstream"
    prepare_upstream "$upstream"
    ordinal=0
    while IFS= read -r barcode; do
      [[ -n "$barcode" ]] || continue
      ordinal=$((ordinal + 1))
      printf -v cell 'CELL%04d-%s' "$ordinal" "$barcode"
      begin_step "$species · Smart-seq2 alignment · $cell"
      star_align_bulk_sample "$index_dir" "$case_dir/${cell}_R1.fastq.gz" "$case_dir/${cell}_R2.fastq.gz" "$cell" "$upstream" "$work"
      end_step
    done < "$REPO_ROOT/$species/sc-rna-seq/10x_3prime/barcodes.tsv"
  done
}

align_markdup_pe() {
  local species="$1" r1="$2" r2="$3" sample="$4" output_dir="$5" work="$6"
  local qname="$work/${sample}.qname.bam"
  local fixmate="$work/${sample}.fixmate.bam"
  local sorted="$output_dir/${sample}.sorted.bam"
  local dedup="$output_dir/${sample}.dedup.bam"

  bwa mem -t "$THREADS" -R "$(read_group "$sample")" "$(bwa_prefix "$species")" "$r1" "$r2" \
    | samtools view -u - \
    | samtools sort -n -@ "$THREADS" -o "$qname" -
  samtools fixmate -m -@ "$THREADS" "$qname" "$fixmate"
  samtools sort -@ "$THREADS" -o "$sorted" "$fixmate"
  samtools index -@ "$THREADS" "$sorted"
  samtools markdup -s -@ "$THREADS" "$sorted" "$dedup" 2> "$output_dir/${sample}.dup_metrics.txt"
  samtools index -@ "$THREADS" "$dedup"
}

write_fragments() {
  local bam="$1" sample="$2" output_dir="$3" work="$4"
  local qname="$work/${sample}.dedup.qname.bam"
  samtools sort -n -@ "$THREADS" -o "$qname" "$bam"
  bedtools bamtobed -bedpe -i "$qname" \
    | awk -v sample="$sample" 'BEGIN {OFS="\t"} $1==$4 && $2>=0 && $5>=0 {start=($2<$5?$2:$5); end=($3>$6?$3:$6); print $1,start,end,sample,1}' \
    | sort -k1,1 -k2,2n \
    | bgzip -c > "$output_dir/${sample}.fragments.tsv.gz"
  tabix -f -p bed "$output_dir/${sample}.fragments.tsv.gz"
}

call_peaks() {
  local bam="$1" sample="$2" reference="$3" output_dir="$4" work="$5"
  local genome_size peak_dir
  genome_size="$(awk '{total += $2} END {print total}' "$reference.fai")"
  peak_dir="$work/${sample}.macs3"
  mkdir -p "$peak_dir"
  macs3 callpeak \
    -t "$bam" \
    -f BAMPE \
    -g "$genome_size" \
    --keep-dup all \
    -n "$sample" \
    --outdir "$peak_dir"
  cp "$peak_dir/${sample}_peaks.narrowPeak" "$output_dir/${sample}.peaks.narrowPeak"
}

generate_chromatin() {
  require_cmd bwa
  require_cmd samtools
  require_cmd bedtools
  require_cmd bgzip
  require_cmd tabix
  require_cmd macs3

  local species assay case_name case_dir upstream sample reference work
  for species in Yeast Homo_chr21; do
    species_selected "$species" || continue
    reference="$(reference_path "$species")"
    for assay in atac-seq chip-seq cuttag; do
      case "$assay" in
        atac-seq) case_name="atac"; sample="$(sample_prefix "$species")_ATAC" ;;
        chip-seq) case_name="chip"; sample="$(sample_prefix "$species")_CHIP" ;;
        cuttag) case_name="cuttag"; sample="$(sample_prefix "$species")_CUTTAG" ;;
      esac
      case_dir="$REPO_ROOT/$species/$assay/$case_name"
      upstream="$case_dir/upstream"
      work="$WORK_ROOT/chromatin-$species-$case_name"
      mkdir -p "$work"
      prepare_upstream "$upstream"
      begin_step "$species · $assay · align and mark duplicates"
      align_markdup_pe "$species" "$case_dir/${sample}_R1.fastq.gz" "$case_dir/${sample}_R2.fastq.gz" "$sample" "$upstream" "$work"
      end_step
      begin_step "$species · $assay · MACS3 peak calling"
      call_peaks "$upstream/${sample}.dedup.bam" "$sample" "$reference" "$upstream" "$work"
      end_step
      if [[ "$assay" == "atac-seq" || "$assay" == "cuttag" ]]; then
        begin_step "$species · $assay · fragment generation and indexing"
        write_fragments "$upstream/${sample}.dedup.bam" "$sample" "$upstream" "$work"
        end_step
      fi
    done
  done
}

hic_bin_size() {
  case "$1" in
    Ecoli) printf '1000' ;;
    Yeast) printf '5000' ;;
    Homo_chr21) printf '10000' ;;
    *) die "no Hi-C bin size configured for $1" ;;
  esac
}

hic_resolutions() {
  case "$1" in
    Ecoli) printf '1000,5000,10000' ;;
    Yeast) printf '5000,10000,50000' ;;
    Homo_chr21) printf '10000,50000,100000' ;;
    *) die "no Hi-C resolutions configured for $1" ;;
  esac
}

generate_hic() {
  require_cmd bwa
  require_cmd pairtools
  require_cmd pairix
  require_cmd bgzip
  require_cmd cooler
  require_cmd java
  [[ -n "$JUICER_TOOLS_JAR" && -s "$JUICER_TOOLS_JAR" ]] \
    || die "hic phase requires JUICER_TOOLS_JAR to point to a Juicer Tools jar"

  local species enzyme case_dir upstream sample reference chromsizes work raw_pairs valid_pairs
  local pairs_gz pairs_index pair_stats bin_size resolutions short_pairs cool mcool
  local staged_pairs staged_pairs_index staged_stats staged_cool staged_mcool staged_hic
  for species in Ecoli Yeast Homo_chr21; do
    species_selected "$species" || continue
    case "$species" in
      Ecoli) enzyme="mboi"; sample="ECOLI_HIC" ;;
      Yeast) enzyme="hindiii"; sample="YEAST_HIC" ;;
      Homo_chr21) enzyme="hindiii"; sample="H21_HIC" ;;
    esac
    case_dir="$REPO_ROOT/$species/Hi-C/$enzyme"
    upstream="$case_dir/upstream"
    mkdir -p "$upstream"
    reference="$(reference_path "$species")"
    work="$WORK_ROOT/hic-$species"
    mkdir -p "$work"
    chromsizes="$work/chrom.sizes"
    cut -f1,2 "$reference.fai" > "$chromsizes"
    raw_pairs="$work/${sample}.dedup.pairs"
    valid_pairs="$work/${sample}.valid.pairs"
    pairs_gz="$upstream/${sample}.pairs.gz"
    pairs_index="$pairs_gz.px2"
    pair_stats="$upstream/${sample}.pairtools_stats.txt"
    bin_size="$(hic_bin_size "$species")"
    resolutions="$(hic_resolutions "$species")"
    cool="$upstream/${sample}.${bin_size}.cool"
    mcool="$upstream/${sample}.mcool"

    if all_files_nonempty "$pairs_gz" "$pairs_index" "$pair_stats"; then
      skip_step "$species · Hi-C · pairs — complete outputs already exist"
    else
      staged_pairs="$work/${sample}.pairs.gz"
      staged_pairs_index="$staged_pairs.px2"
      staged_stats="$work/${sample}.pairtools_stats.txt"
      begin_step "$species · Hi-C · align, parse and deduplicate pairs"
      bwa mem -SP -t "$THREADS" -R "$(read_group "$sample")" "$(bwa_prefix "$species")" \
        "$case_dir/${sample}_R1.fastq.gz" "$case_dir/${sample}_R2.fastq.gz" \
        | pairtools parse --chroms-path "$chromsizes" --assembly "$species" --drop-sam \
        | pairtools sort --nproc "$THREADS" --tmpdir "$work" \
        | pairtools dedup --output-stats "$staged_stats" -o "$raw_pairs"
      pairtools select '(pair_type == "UU")' "$raw_pairs" -o "$valid_pairs"
      bgzip -c "$valid_pairs" > "$staged_pairs"
      pairix -f "$staged_pairs"
      mv -f "$staged_pairs" "$pairs_gz"
      mv -f "$staged_pairs_index" "$pairs_index"
      mv -f "$staged_stats" "$pair_stats"
      end_step
    fi

    if all_files_nonempty "$cool" "$mcool"; then
      skip_step "$species · Hi-C · Cooler matrices — complete outputs already exist"
    else
      staged_cool="$work/${sample}.${bin_size}.cool"
      staged_mcool="$work/${sample}.mcool"
      begin_step "$species · Hi-C · Cooler matrices"
      cooler cload pairix "$chromsizes:$bin_size" "$pairs_gz" "$staged_cool"
      cooler zoomify \
        --nproc "$THREADS" \
        --resolutions "$resolutions" \
        -o "$staged_mcool" \
        "$staged_cool"
      mv -f "$staged_cool" "$cool"
      mv -f "$staged_mcool" "$mcool"
      end_step
    fi

    begin_step "$species · Hi-C · Juicer .hic"
    short_pairs="$work/${sample}.juicer.short.txt"
    staged_hic="$work/${sample}.hic"
    gzip -cd "$pairs_gz" \
      | awk 'BEGIN {OFS="\t"} $1 !~ /^#/ {s1=($6=="+"?0:16); s2=($7=="+"?0:16); print s1,$2,$3,0,s2,$4,$5,1}' \
      > "$short_pairs"
    java -Xmx2g -jar "$JUICER_TOOLS_JAR" pre \
      -n \
      -r "$resolutions" \
      "$short_pairs" \
      "$staged_hic" \
      "$chromsizes"
    [[ -s "$staged_hic" ]] || die "Juicer did not create a non-empty .hic file for $sample"
    mv -f "$staged_hic" "$upstream/${sample}.hic"
    end_step
  done
}

copy_first_match() {
  local search_dir="$1" pattern="$2" destination="$3"
  local source
  source="$(find "$search_dir" -type f -name "$pattern" -print -quit)"
  [[ -n "$source" ]] || die "expected output not found: $search_dir/$pattern"
  cp "$source" "$destination"
}

generate_one_methylation_case() {
  local case_dir="$1" sample="$2" output_dir="$3" work="$4"
  local reference genome_dir alignment_dir dedup_dir extract_dir aligned dedup_bam
  reference="$(reference_path Homo_chr21)"
  genome_dir="$work/bismark-genome"
  alignment_dir="$work/alignment"
  dedup_dir="$work/dedup"
  extract_dir="$work/extract"
  mkdir -p "$genome_dir" "$alignment_dir" "$dedup_dir" "$extract_dir"
  cp "$reference" "$genome_dir/genome.fa"

  bismark_genome_preparation --bowtie2 --parallel "$THREADS" "$genome_dir"
  bismark \
    --genome "$genome_dir" \
    --parallel "$THREADS" \
    --basename "$sample" \
    --output_dir "$alignment_dir" \
    -1 "$case_dir/${sample}_R1.fastq.gz" \
    -2 "$case_dir/${sample}_R2.fastq.gz"
  aligned="$(find "$alignment_dir" -type f -name "${sample}*pe.bam" -print -quit)"
  [[ -n "$aligned" ]] || die "Bismark paired BAM not found for $sample"

  deduplicate_bismark --paired --bam --output_dir "$dedup_dir" "$aligned"
  dedup_bam="$(find "$dedup_dir" -type f -name '*.deduplicated.bam' -print -quit)"
  [[ -n "$dedup_bam" ]] || die "deduplicated Bismark BAM not found for $sample"

  samtools sort -@ "$THREADS" -o "$output_dir/${sample}.sorted.bam" "$dedup_bam"
  samtools index -@ "$THREADS" "$output_dir/${sample}.sorted.bam"

  bismark_methylation_extractor \
    --paired-end \
    --comprehensive \
    --merge_non_CpG \
    --bedGraph \
    --counts \
    --cytosine_report \
    --gzip \
    --genome_folder "$genome_dir" \
    --output_dir "$extract_dir" \
    "$dedup_bam"

  copy_first_match "$alignment_dir" '*_PE_report.txt' "$output_dir/${sample}.alignment_report.txt"
  copy_first_match "$dedup_dir" '*deduplication_report.txt' "$output_dir/${sample}.duplication_report.txt"
  copy_first_match "$extract_dir" '*.bismark.cov.gz' "$output_dir/${sample}.coverage.cov.gz"
  copy_first_match "$extract_dir" '*CpG_context*.txt.gz' "$output_dir/${sample}.CpG_context.txt.gz"
  copy_first_match "$extract_dir" '*CpG_report.txt.gz' "$output_dir/${sample}.cytosine_report.txt.gz"
  copy_first_match "$extract_dir" '*.bedGraph.gz' "$output_dir/${sample}.bedGraph.gz"
  copy_first_match "$extract_dir" '*splitting_report.txt' "$output_dir/${sample}.methylation_report.txt"
  copy_first_match "$extract_dir" '*M-bias.txt' "$output_dir/${sample}.M-bias.txt"
}

generate_methylation() {
  species_selected Homo_chr21 || return 0
  require_cmd bismark_genome_preparation
  require_cmd bismark
  require_cmd deduplicate_bismark
  require_cmd bismark_methylation_extractor
  require_cmd samtools

  local assay case_name sample case_dir upstream work publish_dir
  local -a publish_files
  for assay in wgbs em-seq; do
    case "$assay" in
      wgbs) case_name="wgbs"; sample="H21_WGBS" ;;
      em-seq) case_name="emseq"; sample="H21_EMSEQ" ;;
    esac
    case_dir="$REPO_ROOT/Homo_chr21/$assay/$case_name"
    upstream="$case_dir/upstream"
    work="$WORK_ROOT/methylation-$case_name"
    publish_dir="$work/publish"
    mkdir -p "$work" "$publish_dir" "$upstream"
    begin_step "Homo_chr21 · $assay · Bismark upstream workflow"
    generate_one_methylation_case "$case_dir" "$sample" "$publish_dir" "$work"
    shopt -s nullglob
    publish_files=("$publish_dir"/*)
    shopt -u nullglob
    (( ${#publish_files[@]} > 0 )) || die "no methylation outputs were staged for $sample"
    mv -f -- "${publish_files[@]}" "$upstream/"
    end_step
  done
}

generate_somatic() {
  species_selected Homo_chr21 || return 0
  require_cmd bwa
  require_cmd samtools
  require_cmd gatk
  require_cmd tabix

  local case_name pair normal tumor case_dir upstream normal_bam tumor_bam unfiltered output_vcf
  for case_name in tumor_normal_clonal tumor_normal_low_purity tumor_normal_subclonal; do
    case "$case_name" in
      tumor_normal_clonal) pair="H21_CLONAL" ;;
      tumor_normal_low_purity) pair="H21_LOW_PURITY" ;;
      tumor_normal_subclonal) pair="H21_SUBCLONAL" ;;
    esac
    normal="${pair}_N"
    tumor="${pair}_T"
    case_dir="$REPO_ROOT/Homo_chr21/somatic/$case_name"
    upstream="$case_dir/upstream"
    prepare_upstream "$upstream"
    normal_bam="$upstream/${normal}.sorted.bam"
    tumor_bam="$upstream/${tumor}.sorted.bam"
    unfiltered="$upstream/${pair}.unfiltered.vcf.gz"
    output_vcf="$upstream/${pair}.somatic.vcf.gz"

    begin_step "Homo_chr21 · somatic · matched-normal alignment · $case_name"
    align_bwa_pe Homo_chr21 "$case_dir/${normal}_R1.fastq.gz" "$case_dir/${normal}_R2.fastq.gz" "$normal" "$normal_bam"
    end_step
    begin_step "Homo_chr21 · somatic · tumor alignment · $case_name"
    align_bwa_pe Homo_chr21 "$case_dir/${tumor}_R1.fastq.gz" "$case_dir/${tumor}_R2.fastq.gz" "$tumor" "$tumor_bam"
    end_step

    begin_step "Homo_chr21 · somatic · Mutect2 and filtering · $case_name"
    gatk Mutect2 \
      -R "$(reference_path Homo_chr21)" \
      -I "$tumor_bam" -tumor "$tumor" \
      -I "$normal_bam" -normal "$normal" \
      -L "$case_dir/targets.bed" \
      --native-pair-hmm-threads "$THREADS" \
      -O "$unfiltered"
    gatk FilterMutectCalls \
      -R "$(reference_path Homo_chr21)" \
      -V "$unfiltered" \
      --filtering-stats "$upstream/${pair}.filtering_stats.tsv" \
      -O "$output_vcf"
    [[ -s "$output_vcf.tbi" ]] || tabix -f -p vcf "$output_vcf"
    if [[ -s "$unfiltered.stats" ]]; then
      cp "$unfiltered.stats" "$upstream/${pair}.caller_stats.txt"
    fi
    end_step
  done
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --phase)
      [[ $# -ge 2 ]] || die "--phase requires a value"
      PHASES="$2"
      shift 2
      ;;
    --species)
      [[ $# -ge 2 ]] || die "--species requires a value"
      SPECIES_FILTER="$2"
      shift 2
      ;;
    --no-color)
      COLOR_MODE="never"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

validate_selection
WORK_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/bioflow-upstream.XXXXXX")"
setup_logging
trap cleanup EXIT
trap 'on_error "$?" "$LINENO" "$BASH_COMMAND"' ERR

printf '%bBioFlowAI upstream fixture generation%b\n' "$C_BOLD" "$C_RESET" >&4
info "Phases: $PHASES"
info "Species: $SPECIES_FILTER"
info "Threads: $THREADS"
info "Temporary workspace: $WORK_ROOT"

if phase_selected reference; then
  begin_phase "REFERENCE INDEXES"
  build_reference_indexes
  end_phase
fi
if phase_selected wgs; then
  begin_phase "WHOLE-GENOME SEQUENCING"
  generate_wgs
  end_phase
fi
if phase_selected rna; then
  begin_phase "BULK RNA-SEQ"
  generate_rna
  end_phase
fi
if phase_selected scrna; then
  begin_phase "SINGLE-CELL RNA-SEQ"
  generate_scrna
  end_phase
fi
if phase_selected chromatin; then
  begin_phase "CHROMATIN ASSAYS"
  generate_chromatin
  end_phase
fi
if phase_selected hic; then
  begin_phase "HI-C"
  generate_hic
  end_phase
fi
if phase_selected methylation; then
  begin_phase "DNA METHYLATION"
  generate_methylation
  end_phase
fi
if phase_selected somatic; then
  begin_phase "SOMATIC TUMOR-NORMAL"
  generate_somatic
  end_phase
fi

printf '\n%b✓ ALL REQUESTED PHASES COMPLETE%b  %b(%s, %d steps)%b\n' \
  "$C_GREEN$C_BOLD" "$C_RESET" "$C_DIM" "$(format_duration "$((SECONDS - RUN_STARTED))")" "$STEP_NUMBER" "$C_RESET" >&4
