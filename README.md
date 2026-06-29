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

| Sample | Species | Assembly | Replicons | Notes |
|--------|---------|----------|-----------|-------|
| S2052 | *V. coralliilyticus* | Draft (Unicycler) | 69 contigs | Contig 9: 224 kb plasmid |
| S2753 | *P. galatheae* | Draft (Unicycler) | 79 contigs | Very similar to S2754 |
| S2754 | *P. galatheae* | Draft (Unicycler) | 76 contigs | Very similar to S2753 |
| S2052ref | *V. coralliilyticus* | [GCF_000967465.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000967465.2/) | 2 chr + 1 plasmid | Complete reference for S2052 |
| S2753ref | *P. galatheae* | Complete | 2 chr + 1 plasmid (323 kb) | Complete reference for S2753 |
| S2754ref | *P. galatheae* | Complete | 2 chr + 1 plasmid (323 kb) | Complete reference for S2754 |
| S26ref | *P. piscinae* | [GCF_000826835.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000826835.2/) | 1 chr + 4 plasmids | No matching draft assembly |

Raw paired-end reads are in `input_reads/<sample>/`. For `*ref` samples, the same reads as the corresponding draft are mapped to the complete reference.

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

**Multi-criteria chromosomal baseline:**

The chromosomal baseline is identified through a decision tree that combines coverage, single-copy core genes (SCGs), and geNomad classification:

1. **Anchor**: The largest contig's geometric-mean coverage sets the center. In bacterial genomes, the largest contig is almost certainly chromosomal and has the most splits (lowest variance).
2. **Coverage window**: All contigs >= 20 kb within **±10%** of the anchor coverage form the candidate pool.
3. **SCG positive filter**: Candidates with `Bacteria_71 > 0` (bacterial single-copy core genes detected by anvi'o HMMs) are kept. This confirms chromosomal identity — plasmids rarely carry SCGs.
4. **geNomad negative filter**: Candidates classified as `plasmid` by geNomad's filtered summary are excluded.

Each filter is defensive: if it would empty the candidate set (e.g., no SCG data available, or geNomad not run), it is skipped, falling back to the previous step's result.

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

### V. coralliilyticus S2052 — draft vs complete reference

| | Draft (S2052) | Complete (S2052ref) |
|---|---|---|
| **Baseline** | 216.1x (13 contigs with SCGs) | 216.5x (2 chromosomes) |
| **Plasmid CN** | **1.61** (contig 9, 224 kb, circular) | **1.60** (NZ_CP063053.1, 224 kb) |
| **rRNA artifacts** | contig 31: 10.01x (collapsed 23S) | none — resolved within chromosomes |

Elevated contigs in draft:

| Contig | Length | CN | Notes |
|--------|--------|----|-------|
| 9 | 223,858 bp | **1.61** | Plasmid (geNomad: 0.985) |
| 31 | 3,253 bp | **10.01** | Collapsed rRNA operon (BLAST: 100% V. coralliilyticus 23S) |

Complete reference:

| Replicon | Accession | Length | Coverage | CN | Bacteria_71 |
|----------|-----------|--------|----------|----|-------------|
| Chromosome 1 | NZ_CP063051.1 | 3,328,634 bp | 216.2x | **1.00** | 69 |
| Chromosome 2 | NZ_CP063052.1 | 1,884,472 bp | 216.7x | **1.00** | 2 |
| Plasmid | NZ_CP063053.1 | 223,859 bp | 347.2x | **1.60** | 0 |

### P. galatheae S2753 — draft vs complete reference

| | Draft (S2753) | Complete (S2753ref) |
|---|---|---|
| **Baseline** | 226.6x (9 contigs with SCGs) | 227.5x (3 replicons) |
| **Plasmids** | none | none |
| **rRNA artifacts** | contig 26: 11.24x, contig 30: 10.38x | none |

Complete reference — 3 replicons, all at ~1.0x:

| Replicon | Length | Coverage | CN | Bacteria_71 |
|----------|--------|----------|----|-------------|
| Chromosome 1 | 3,603,140 bp | 225.9x | **0.99** | 70 |
| Chromosome 2 | 1,628,073 bp | 227.5x | **1.00** | 3 |
| Plasmid (323 kb) | 323,546 bp | 236.8x | **1.04** | 1 (SecY) |

The 323 kb replicon is at unit copy number and carries a duplicate SecY (also present on chr1). geNomad does not classify it as plasmid — no recognized plasmid replication markers detected.

### P. galatheae S2754 — draft vs complete reference

| | Draft (S2754) | Complete (S2754ref) |
|---|---|---|
| **Baseline** | 237.6x (9 contigs with SCGs) | 243.4x (3 replicons) |
| **Plasmids** | none | none |
| **rRNA artifacts** | contig 26: 11.34x, contig 30: 2.78x, contig 31: 10.60x | none |

Complete reference — same architecture as S2753ref:

| Replicon | Length | Coverage | CN | Bacteria_71 |
|----------|--------|----------|----|-------------|
| Chromosome 1 | 3,601,282 bp | 239.2x | **0.98** | 70 |
| Chromosome 2 | 1,628,064 bp | 243.4x | **1.00** | 3 |
| Plasmid (323 kb) | 323,532 bp | 260.1x | **1.07** | 1 (SecY) |

### P. piscinae S26ref — complete reference

*Phaeobacter piscinae* strain S26 ([GCF_000826835.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000826835.2/)) — 1 chromosome + 4 plasmids. The SCG filter is critical here: only the chromosome has bacterial SCGs (Bacteria_71 = 72), so the baseline is anchored on a single replicon.

| Replicon | Accession | Length | Coverage | CN | Bacteria_71 |
|----------|-----------|--------|----------|----|-------------|
| Chromosome | NZ_CP080275.1 | 3,723,659 bp | 101.6x | **1.00** | 72 |
| pS26_248 | NZ_CP080276.1 | 248,209 bp | 83.9x | **0.83** | 0 |
| pS26_80 | NZ_CP080278.1 | 80,357 bp | 93.7x | **0.92** | 0 |
| pS26_68 | NZ_CP080279.1 | 68,485 bp | 169.8x | **1.67** | 0 |
| pS26_106 | NZ_CP080277.1 | 106,329 bp | 190.4x | **1.87** | 0 |

Key observations:
- **pS26_68** and **pS26_106**: elevated copy number (1.67–1.87x), actively maintained above 1:1 ratio
- **pS26_248** and **pS26_80**: sub-1x copy number (0.83–0.92x), suggesting plasmid loss during growth or slower replication than the chromosome
- Without the SCG filter, the old median-of-large-contigs approach would have failed — 4 out of 5 large replicons are plasmids

### Plasmid typing (mob_typer)

Plasmid sequences were extracted from each reference genome (`scripts/extract_ref_plasmids.sh`) and typed with [MOB-suite](https://github.com/phac-nml/mob-suite) `mob_typer` (`scripts/run_mob_typer.sh`), one annotation per FASTA entry. Results in `results/mob_typer/`.

| Sample | Accession | Size (bp) | Replicon type | Relaxase | Predicted mobility | Mash nearest neighbor | Distance | Primary cluster |
|--------|-----------|-----------|---------------|----------|-------------------|----------------------|----------|-----------------|
| S2052ref | NZ_CP063053.1 | 223,859 | — | — | non-mobilizable | *V. coralliilyticus* | 0.031 | AC648 |
| S2753ref | 12165.1 | 323,546 | — | — | non-mobilizable | no match | 1.0 | AA379 |
| S2754ref | 31931.2 | 323,532 | — | — | non-mobilizable | no match | 1.0 | AA379 |
| S26ref — pS26_106 | NZ_CP080277.1 | 106,329 | — | — | non-mobilizable | no match | 1.0 | AA379 |
| S26ref — pS26_248 | NZ_CP080276.1 | 248,209 | rep_cluster_377 | — | non-mobilizable | *P. piscinae* | 0.015 | AC739 |
| S26ref — pS26_68 | NZ_CP080279.1 | 68,485 | rep_cluster_252, rep_cluster_735 | — | non-mobilizable | *P. piscinae* | 0.010 | AA429 |
| S26ref — pS26_80 | NZ_CP080278.1 | 80,357 | rep_cluster_288 | — | non-mobilizable | *P. piscinae* | 0.015 | AB852 |

Key observations:
- **All elements are predicted non-mobilizable** — no relaxase or MPF type detected in any replicon.
- **S2052ref plasmid**: no replicon type recognized, but Mash nearest neighbor is *V. coralliilyticus* (distance 0.031) — consistent with a host-specific plasmid.
- **S26ref pS26_248, pS26_68, pS26_80**: replicon types identified (rep_cluster_252/288/377/735), all with close Mash matches to *Phaeobacter piscinae* — confirmed plasmids with known replication origins.
- **S26ref pS26_106**: no replicon type, no Mash match (distance 1.0) — replication machinery unrecognized by the current MOB-suite database.
- **323 kb replicons (S2753ref / S2754ref)**: no replicon type, no relaxase, no orit, Mash distance 1.0 — no recognized plasmid replication or mobilization machinery. Combined with unit copy number and carriage of essential genes (SecY), this reinforces a chromid-like classification rather than a typical plasmid.

## Files

| File | Description |
|------|-------------|
| `fasta.txt` | Workflow input: genome names and FASTA paths |
| `samples.txt` | Workflow input: sample names, read paths, group assignments |
| `config.json` | Workflow config (references_mode, QC on, HMMs on, heavy annotations off) |
| `run_anvio_workflow.sh` | LSF script: anvi'o metagenomics workflow (12 cores, 8 GB/core, 24 h) |
| `run_genomad.sh` | LSF script: geNomad end-to-end on all 7 samples (12 cores, 8 GB/core, 4 h) |
| `workflow.pdf` | DAG of the Snakemake workflow |

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/00_export_coverage.sh` | Export contig and split-level coverages from anvi'o profiles to `07_COVERAGE/` |
| `scripts/01a_annotate_contigs.py` | Annotate contigs: split-smoothed coverage, HMM hit counts, assembler metadata |
| `scripts/01b_integrate_genomad.py` | Add geNomad plasmid/virus classification to annotations (optional) |
| `scripts/02_compute_copy_numbers.py` | Compute copy numbers from annotated contigs, classify chromosome/plasmid |
| `scripts/extract_ref_plasmids.sh` | Extract plasmid sequences from reference FASTAs by header keyword into `input_plasmidsRef/` |
| `scripts/run_mob_typer.sh` | LSF script: run mob_typer per plasmid sequence on all reference genomes (4 cores, 4 GB, 1 h) |

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
| `results/mob_typer/` | mob_typer typing results per reference sample (one TSV per sample, one row per plasmid) |
| `input_plasmidsRef/` | Plasmid sequences extracted from reference FASTAs, one FASTA per sample |

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
