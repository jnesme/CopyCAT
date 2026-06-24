# CopyCAT — Copy number from Coverage Analysis Tool

Estimate plasmid copy numbers from short-read coverage using anvi'o 9 under the hood. Contigs are split into ~20 kb windows, enabling split-level coverage inspection to verify that coverage is uniform along each contig before reporting a copy number. This guards against chimeric contigs or localized mapping artifacts inflating estimates. Results are cross-validated with geNomad sequence-based classification.

## How to use

### Prerequisites

- [anvi'o 9](https://anvio.org/) (`conda activate anvio-9`)
- [geNomad](https://github.com/apcamargo/genomad) (optional, for sequence-based cross-validation)
- Input genomes in `input_genomes/` and paired-end reads in `input_reads/<sample>/`
- `fasta.txt` and `samples.txt` configured (see [Data](#data) section)

### 1. Run the anvi'o metagenomics workflow

```bash
conda activate anvio-9

# Validate config (dry-run)
anvi-run-workflow -w metagenomics -c config.json --dry-run

# Submit to LSF (or run locally)
bsub < run_anvio_workflow.sh
```

This runs QC, read mapping, contig profiling with 20 kb splits, and HMMs for all samples defined in `samples.txt`.

### 2. Export coverage data

```bash
conda activate anvio-9
bash scripts/00_export_coverage.sh
```

Exports contig-level and split-level coverages to `07_COVERAGE/`. Skips samples already exported. Handles both merged profiles (multi-sample groups) and single-sample profiles automatically.

### 3. Annotate contigs

```bash
python3 scripts/01a_annotate_contigs.py
```

Collects per-contig features: split-smoothed coverage (arithmetic mean, geometric mean, median), coefficient of variation, HMM hit counts per source (Bacteria_71, Archaea_76, Protista_83, Ribosomal_RNA_16S), and assembler metadata when available.

### 4. Cross-validate with geNomad (optional)

```bash
bsub < run_genomad.sh
# After completion:
python3 scripts/01b_integrate_genomad.py
```

Adds geNomad plasmid/virus classification columns to the annotation table. Skips cleanly if geNomad hasn't been run.

### 5. Compute copy numbers

```bash
python3 scripts/02_compute_copy_numbers.py
```

Reads the annotated table, identifies the chromosomal baseline, and computes copy numbers. Classification logic lives here — edit this script to refine rules (e.g., require bacterial SCGs for baseline).

## Data

| Sample | Assembly | Contigs | Reads | Notes |
|--------|----------|---------|-------|-------|
| S2052 | `input_genomes/S2052.scaffolds.fa.gz` | 69 | PE Illumina, ~300 MB/file | Contig 9: 224 kb, circular, 1.59x assembler depth |
| S2753 | `input_genomes/S2753.scaffolds.fa.gz` | 79 | PE Illumina, ~300–320 MB/file | Very similar to S2754 |
| S2754 | `input_genomes/S2754.scaffolds.fa.gz` | 76 | PE Illumina, ~320–340 MB/file | Very similar to S2753 |

Raw paired-end reads are in `input_reads/<sample>/`.

## Method

### 1. Coverage-based copy number estimation (anvi'o)

Uses the **anvi'o metagenomics workflow** in `references_mode`:

1. Quality-filter raw reads (`iu-filter-quality-minoche`)
2. Reformat and simplify contig names (`anvi-script-reformat-fasta`)
3. Build contigs databases with gene calling and 20 kb splits (`anvi-gen-contigs-database`)
4. Map reads to their own assembly with bowtie2
5. Profile BAMs to compute per-split and per-contig coverage (`anvi-profile`, min contig length 1000 bp)
6. Run HMMs for single-copy core gene detection (`anvi-run-hmms`)
7. Extract coverage and compute plasmid copy number as `coverage_plasmid / median_coverage_chromosome`

### Copy number calculation

For each sample, the copy number of a contig is computed as:

```
copy_number = coverage_geom_mean(contig) / median_coverage_geom_mean(chromosomal contigs)
```

**Contig coverage from splits:**

Rather than using anvi'o's position-level mean directly, contig coverage is derived from its 20 kb splits using the **geometric mean** (exp(mean(log(split coverages)))). This naturally dampens outlier splits in log-space — a 2x spike and a 0.5x dip cancel out, which is the right behavior for coverage data. The arithmetic mean and median of split coverages are also reported as diagnostic columns for QC.

For contigs with a single split (< 20 kb), all estimators are identical.

**Assembler-agnostic chromosomal baseline:**

Chromosomal contigs are identified purely from coverage, without requiring assembler-reported depth or circularity:

1. Large contigs (>= 20 kb) form the candidate pool
2. The median of their geometric-mean coverages defines the center
3. Contigs within **±10%** of this median are classified as chromosomal

The percentage window avoids the over-sensitivity of MAD-based approaches on clean data, and remains stable across noisy datasets.

Assembler metadata (depth, circularity) is preserved as annotation columns when available (supports Unicycler, Flye, and similar assemblers) but does not drive the classification.

**Split-level verification:**

The coefficient of variation (CV) of coverage across a contig's 20 kb splits is reported. A low CV (< 0.15) confirms uniform coverage — ruling out chimeric assemblies, partial duplications, or localized mapping artifacts.

**Interpretation:**

- Copy number ~1.0: single-copy, chromosomal
- Copy number 1.5–2.0: low-copy plasmid or recent duplication
- Copy number >> 1: multi-copy plasmid (e.g., ~10x = ~10 copies per chromosome)

### 2. rRNA filtering

Elevated copy number does not imply plasmid — collapsed multi-copy chromosomal elements (rRNA operons, IS elements) inflate coverage the same way. The script queries anvi'o's HMM results to flag contigs containing rRNA genes. Note: anvi'o's default HMMs only detect 16S/18S rRNA (via barrnap); 23S and 5S are not included. BLAST validation against NCBI nt is recommended for small elevated contigs not flagged by HMMs.

### 3. Sequence-based classification (geNomad, optional)

Runs `genomad end-to-end` on profiled contigs (>=1000 bp) to independently classify contigs as plasmid, virus, or chromosome using marker genes and neural networks. Results are merged with copy numbers by `scripts/02_integrate_genomad.py`.

## Results

### S2052 (median chromosomal coverage: 216.7x)

| Contig | Length | Circular | Copy # | rRNA | geNomad | Notes |
|--------|--------|----------|--------|------|---------|-------|
| 9 | 223,858 bp | yes | **1.60** | | plasmid (0.985) | Confirmed plasmid — matches known V. coralliilyticus plasmid |
| 31 | 3,253 bp | no | **9.98** | 23S | chromosome | Collapsed rRNA operon (BLAST: 100% identity to V. coralliilyticus 23S) |

Also detected: provirus in contig 8 (236 kb).

### S2753 (median chromosomal coverage: 226.8x)

| Contig | Length | Copy # | rRNA | geNomad | Notes |
|--------|--------|--------|------|---------|-------|
| 26 | 2,786 bp | **11.23** | | chromosome | Collapsed rRNA operon (pending BLAST confirmation) |
| 30 | 1,077 bp | **10.37** | 16S | chromosome | Collapsed rRNA operon (flagged by HMM) |

Also detected: proviruses in contigs 1 (1.46 Mb), 2 (609 kb), 6 (316 kb).

### S2754 (median chromosomal coverage: 240.5x)

| Contig | Length | Copy # | rRNA | geNomad | Notes |
|--------|--------|--------|------|---------|-------|
| 26 | 2,599 bp | **11.21** | | chromosome | Collapsed rRNA operon (pending BLAST confirmation) |
| 30 | 1,106 bp | **2.74** | | chromosome | Pending BLAST confirmation |
| 31 | 1,077 bp | **10.48** | 16S | chromosome | Collapsed rRNA operon (flagged by HMM) |

Also detected: proviruses in contigs 1 (1.46 Mb), 2 (609 kb), 6 (316 kb).

### Validation against complete reference: S2052ref

The same S2052 reads were mapped to the complete reference assembly of V. coralliilyticus strain S2052 ([GCF_000967465.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000967465.2/)), which consists of 3 replicons:

| Replicon | Accession | Length | Coverage | CN (geom) | rRNA |
|----------|-----------|--------|----------|-----------|------|
| Chromosome 1 | NZ_CP063051.1 | 3,328,634 bp | 216.2x | **1.00** | yes |
| Chromosome 2 | NZ_CP063052.1 | 1,884,472 bp | 216.7x | **1.00** | |
| Plasmid | NZ_CP063053.1 | 223,859 bp | 347.2x | **1.60** | |

Key observations:
- Plasmid copy number (**1.60**) matches the draft assembly result exactly, validating the method
- Both chromosomes are at 1.00x — the two-chromosome Vibrio architecture is correctly resolved
- No elevated rRNA contigs — the ~10 rRNA operons are resolved within the chromosomes instead of collapsing into separate high-copy contigs
- rRNA is detected on chromosome 1 (16S via HMM); 23S is present but not flagged by anvi'o's default HMMs

This confirms that CopyCAT's coverage-based approach is robust: it produces the same copy number estimate regardless of whether the input is a fragmented draft assembly or a complete reference genome.

## Files

| File | Description |
|------|-------------|
| `fasta.txt` | Workflow input: genome names and FASTA paths |
| `samples.txt` | Workflow input: sample names, read paths, group assignments |
| `config.json` | Workflow config (references_mode, QC on, HMMs on, heavy annotations off) |
| `run_anvio_workflow.sh` | LSF script: anvi'o metagenomics workflow (12 cores, 8 GB/core, 24 h) |
| `run_genomad.sh` | LSF script: geNomad end-to-end on all 3 samples (12 cores, 8 GB/core, 4 h) |
| `workflow.pdf` | DAG of the Snakemake workflow |

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/00_export_coverage.sh` | Export contig and split-level coverages from anvi'o profiles to `07_COVERAGE/` |
| `scripts/01a_annotate_contigs.py` | Annotate contigs: split-smoothed coverage, HMM hit counts, assembler metadata |
| `scripts/01b_integrate_genomad.py` | Add geNomad plasmid/virus classification to annotations (optional) |
| `scripts/02_compute_copy_numbers.py` | Compute copy numbers from annotated contigs, classify chromosome/plasmid |

## Output directories

| Directory | Contents |
|-----------|----------|
| `00_LOGS/` | Snakemake rule logs |
| `01_QC/` | Quality-filtered reads and QC report |
| `02_FASTA/` | Reformatted reference FASTAs and reformat reports |
| `03_CONTIGS/` | Anvi'o contigs databases (with HMMs) |
| `04_MAPPING/` | BAM files (sorted, indexed) |
| `05_ANVIO_PROFILE/` | Per-sample anvi'o profile databases |
| `06_MERGED/` | README files (single sample per group, no merge needed) |
| `07_COVERAGE/` | Exported contig and split coverage files |
| `genomad_out/` | geNomad output per sample |
| `results/` | Final tables (see below) |

## Result files

| File | Description |
|------|-------------|
| `results/contig_annotations.tsv` | Per-contig features: coverage stats, HMM counts, assembler metadata (+geNomad if 01b run) |
| `results/copy_numbers.tsv` | Final copy numbers with classification |
| `results/copy_numbers_summary.txt` | Human-readable summary |

## Roadmap: multi-condition copy number comparison

The current implementation estimates copy numbers from a single read set mapped to a single assembly. For a single isolate with reads, the assembler already provides depth and circularity — CopyCAT adds split-level verification and rRNA filtering, but the core result is the same.

The real power of a coverage-based approach is **comparing copy numbers across conditions**, where assembly alone cannot answer the question. This is the planned next feature.

### Concept

One reference assembly, multiple read sets from different conditions. For each condition, map reads independently and compute per-contig copy numbers. Output a matrix of copy numbers (contigs x conditions) that reveals how plasmid dosage responds to environmental changes.

### Example use cases

- **Antibiotic stress**: Does plasmid copy number increase under sub-MIC antibiotic exposure? AMR gene dosage directly affects resistance level.
- **Growth phase**: Log vs stationary phase — plasmid replication may decouple from chromosomal replication.
- **Temperature / salinity shifts**: Relevant for marine Vibrio — environmental conditions may modulate plasmid maintenance.
- **Serial passaging**: Track plasmid stability over generations — does copy number drift or stay stable?

### Design sketch

The anvi'o workflow already supports this natively. In `samples.txt`, multiple read sets can map to the same group (assembly):

```
sample          r1              r2              group
condition_A     reads_A_R1.fq   reads_A_R2.fq   reference
condition_B     reads_B_R1.fq   reads_B_R2.fq   reference
condition_C     reads_C_R1.fq   reads_C_R2.fq   reference
```

With `references_mode: true`, anvi'o maps each sample to the reference, profiles each independently, then **merges** them — producing a single profile DB with per-sample coverage for every split. The merged profile already contains the multi-condition coverage matrix; the analysis script just needs to compute copy numbers per sample against a shared chromosomal baseline.

### Statistical considerations

- **Shared vs per-condition baseline**: Use a shared chromosomal baseline (median across all conditions) or compute one per condition? Per-condition is more robust if sequencing depth varies, but a shared baseline makes fold-changes directly comparable.
- **Normalization**: Total read count differs between conditions. Normalizing to chromosomal coverage (which CopyCAT already does) inherently accounts for this — copy number is a ratio, not an absolute value.
- **Replication gradient**: Fast-growing cells show a coverage gradient from origin to terminus on the chromosome. This inflates coverage near the origin and could bias the baseline. The split-level geometric mean helps dampen this, but for growth-rate comparisons, origin-proximal and terminus-proximal splits should be compared explicitly.
- **Statistical testing**: With biological replicates, copy number differences between conditions can be tested (e.g., Mann-Whitney on per-split copy numbers). Without replicates, only descriptive comparison is possible.
