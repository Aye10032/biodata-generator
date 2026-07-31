# Upstream Fixture Generation

`generate_upstream.sh` generates reference indexes and command-line upstream
intermediates from the compact simulated FASTQs in this repository. It does not
publish expression matrices, differential-analysis results, annotations,
clusters, motifs, TADs, loops, or visualization products.

## Outputs

The script writes assay outputs below each case:

```text
<species>/<assay>/<case>/upstream/
```

Reference preparation writes `.fai` and `.dict` beside each FASTA and stores
BWA/minimap2 indexes below:

```text
<species>/reference/index/
  bwa/
  minimap2/
```

The generated upstream formats are:

- WGS: sorted BAM/BAI for every technology, PAF for long reads, and
  VCF/TBI plus BCF/CSI for each Illumina paired-end case.
- Bulk RNA-seq: coordinate-sorted BAM/BAI, STAR splice junctions, and STAR
  alignment logs. No count matrix is retained.
- 10x scRNA-seq: CB/UB-tagged BAM/BAI and an alignment log. STARsolo matrices
  are created only in the temporary work directory and are deleted afterward.
- Smart-seq2: one coordinate-sorted BAM/BAI and alignment log per cell.
- ATAC-seq, ChIP-seq, and CUT&Tag: sorted and deduplicated BAM/BAI,
  duplication metrics, and MACS3 narrowPeak. ATAC-seq and CUT&Tag also receive
  bgzip/Tabix-indexed fragment files.
- Hi-C: pairtools pairs/PX2, pairtools statistics, Cooler `.cool`/`.mcool`, and
  Juicer `.hic`. Juicer normalization is disabled because these smoke fixtures
  are intentionally too sparse for stable normalization vectors.
- WGBS and EM-seq: Bismark BAM/BAI, CpG context, cytosine report, coverage,
  bedGraph, and native command reports.
- Somatic: matched tumor/normal BAM/BAI, Mutect2 VCF/TBI, and caller stats.

Because the fixtures are deliberately shallow, a caller may produce a valid
but empty VCF or peak file. Simulator truth is not renamed or copied into
`upstream/`.

## 1. Install the Pixi Environment

From `Scripts/other_scripts/`:

```bash
pixi install --locked
```

The committed `pixi.toml` defines the direct dependencies and `pixi.lock`
pins the resolved environment. `pixi run` also installs the environment
automatically when needed, but the explicit install step makes environment
problems visible before a long generation phase starts.

The Hi-C phase additionally requires a Juicer Tools jar. Obtain the jar through
the normal Juicer distribution channel and provide its absolute path:

```bash
export JUICER_TOOLS_JAR=/absolute/path/to/juicer_tools.jar
```

The other phases do not use this variable.

## 2. Generate Reference Indexes First

Run the reference phase before phases that align reads:

```bash
pixi run reference
```

The reference phase is idempotent. Existing non-empty `.fai`, `.dict`, BWA,
and minimap2 indexes are reported as colored `SKIP` steps.

Limit work to one or several reference systems with `--species`:

```bash
pixi run upstream --phase reference --species Ecoli,Yeast
```

Accepted names are `Lambda`, `Ecoli`, `Yeast`, and `Homo_chr21`.

## 3. Run Upstream Phases

Run one phase at a time:

```bash
pixi run wgs
pixi run rna
pixi run scrna
pixi run chromatin
pixi run hic
pixi run methylation
pixi run somatic
```

The Hi-C phase is resumable. If pairs and Cooler outputs already exist after an
interrupted run, they are reported as `SKIP`; Juicer `.hic` is rebuilt through
a temporary file and published atomically. Run `pixi run hic` to rebuild all
three `.hic` files consistently without repeating alignment or Cooler work. To
limit the resume to human chr21:

```bash
pixi run upstream --phase hic --species Homo_chr21
```

The methylation phase is transactionally published: WGBS/EM-seq outputs are
completed under the temporary workspace first and moved into `upstream/` only
after the whole case succeeds. A failed run can therefore be rerun directly;
existing partial files do not need to be removed manually.

Several phases can be selected in one invocation:

```bash
pixi run upstream \
  --phase rna,scrna,chromatin \
  --species Yeast,Homo_chr21
```

After setting `JUICER_TOOLS_JAR`, all phases can be run in dependency order:

```bash
pixi run all
```

Use `THREADS` to control parallelism. Temporary STAR/Bismark indexes and other
working files use `TMPDIR` and are removed when the command exits:

```bash
THREADS=8 TMPDIR=/path/with/free/space \
  pixi run upstream --phase rna --species Homo_chr21
```

## Logging

Interactive terminal output uses three visually separate levels:

```text
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  PHASE  WHOLE-GENOME SEQUENCING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

┌─ STEP 01  Ecoli · WGS · Illumina paired-end alignment
│ tool │ [latest bwa/samtools message, updated in place]
└─ DONE     Ecoli · WGS · Illumina paired-end alignment  (1s)
```

- Blue banners identify phases.
- Cyan step headers and green completion lines identify orchestration work.
- During a running step, third-party stdout/stderr is collapsed into one
  transient `│ tool │` line that updates in place.
- A successful step clears its transient tool line and proceeds directly to
  `DONE`; full successful tool output is discarded.
- Failed commands show the phase, current step, line, exit status, command, and
  the final 20 combined stdout/stderr lines from that step.
- Existing reference indexes appear as yellow `SKIP` steps.

Colors are enabled only for an interactive terminal. Redirected logs remain
plain text. Disable colors explicitly with either form:

```bash
pixi run upstream --phase wgs --no-color
NO_COLOR=1 pixi run wgs
```

To keep a concise plain-text orchestration log:

```bash
NO_COLOR=1 pixi run wgs 2>&1 | tee wgs.log
```

## Safety and Restart Behavior

- Paths are resolved relative to the script, so the command may be launched
  from any working directory. Pixi tasks themselves should be launched from
  `Scripts/other_scripts/`, or with Pixi's manifest-path option.
- The script never deletes dataset files.
- It refuses to enter a case whose `upstream/` directory is already non-empty,
  except for resumable Hi-C and transactionally published methylation phases.
- Reference index creation is idempotent and skips existing non-empty indexes.
- After another interrupted assay phase, inspect and remove only that phase's
  incomplete `upstream/` directory before rerunning it. Hi-C detects completed
  pairs/Cooler outputs, while methylation replaces outputs only after successful
  staging; neither requires removing the directory.

This strict behavior is intentional: reruns cannot silently mix files produced
by different tool versions or parameter sets.
