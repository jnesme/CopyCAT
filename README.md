# CopyCAT — Copy number from Coverage Analysis Tool

Initial test case applying the CopyCAT approach to three Vibrio isolate genomes (S2052, S2753, S2754) from the Galathea3 expedition. Built on anvi'o 9 under the hood: contigs are split into ~20 kb windows, enabling split-level coverage inspection to verify that coverage is uniform along each contig before reporting a copy number. This guards against chimeric contigs or localized mapping artifacts inflating estimates. Results are cross-validated with geNomad sequence-based classification.

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
copy_number = coverage(contig) / median_coverage(chromosomal contigs)
```

**Establishing the chromosomal baseline:**

1. Contigs are classified as chromosomal if they satisfy all three criteria:
   - Length >= 20 kb (long enough to be reliable)
   - Assembler-reported depth between 0.8x and 1.3x (consistent with single-copy)
   - Not marked as circular by the assembler
2. The **median** coverage across all chromosomal contigs is used as the baseline. Median is preferred over mean because it is robust to outliers from rRNA operons or repeat regions.

**Split-level verification:**

Before trusting a contig's copy number, the coefficient of variation (CV) of coverage across its 20 kb splits is checked. A low CV (< 0.15) confirms that coverage is uniform along the contig — ruling out chimeric assemblies, partial duplications, or localized mapping artifacts that would inflate the mean coverage. All elevated contigs in this dataset passed this check.

**Interpretation:**

- Copy number ~1.0: single-copy, chromosomal
- Copy number 1.5–2.0: low-copy plasmid or recent duplication
- Copy number >> 1: multi-copy plasmid (e.g., ~10x = ~10 copies per chromosome)

### 2. Sequence-based classification (geNomad)

Runs `genomad end-to-end` on profiled contigs (>=1000 bp) to independently classify contigs as plasmid, virus, or chromosome using marker genes and neural networks.

### 3. Split-smoothed copy numbers (assembler-agnostic)

An alternative approach (`scripts/04_copy_numbers_split_smoothed.py`) that removes the dependency on assembler-specific metadata:

**Contig coverage from splits:**

Instead of using anvi'o's position-level mean directly, the contig coverage is derived from its 20 kb splits. Three estimators are computed side by side:

- **Arithmetic mean**: standard average of split coverages. Sensitive to outlier splits.
- **Geometric mean**: exp(mean(log(split coverages))). Naturally dampens outliers in log-space — a 2x spike and a 0.5x dip cancel out, which is the right behavior for coverage. Used as the primary estimator.
- **Median**: most robust to outliers but ignores information from non-central splits.

For contigs with a single split (< 20 kb), all three are identical.

**Assembler-agnostic chromosomal baseline:**

Chromosomal contigs are identified purely from coverage, without requiring assembler-reported depth or circularity:

1. Large contigs (>= 20 kb) form the candidate pool
2. The median of their geometric-mean coverages defines the center
3. Contigs within **±10%** of this median are classified as chromosomal

The percentage window avoids the over-sensitivity of MAD-based approaches on clean data, and remains stable across noisy datasets.

Assembler metadata (depth, circularity) is preserved as annotation columns when available (supports Unicycler, Flye, and similar assemblers) but does not drive the classification.

### 4. Integration

Results from both methods are merged. The small high-copy contigs (1–3 kb, ~10x copy number) are missed by geNomad, likely because they are too short for confident marker-based detection.

## Results

### S2052 (median chromosomal coverage: 216.7x)

| Contig | Length | Circular | Copy # | geNomad | Notes |
|--------|--------|----------|--------|---------|-------|
| 9 | 223,858 bp | yes | **1.60** | plasmid (0.985) | Large plasmid, 7 hallmark genes, uniform split coverage |
| 31 | 3,253 bp | no | **9.98** | chromosome | Small high-copy element |

Also detected: provirus in contig 8 (236 kb).

### S2753 (median chromosomal coverage: 226.8x)

| Contig | Length | Copy # | geNomad | Notes |
|--------|--------|--------|---------|-------|
| 26 | 2,786 bp | **11.23** | chromosome | Small high-copy element |
| 30 | 1,077 bp | **10.37** | chromosome | Small high-copy element |

Also detected: proviruses in contigs 1 (1.46 Mb), 2 (609 kb), 6 (316 kb).

### S2754 (median chromosomal coverage: 240.5x)

| Contig | Length | Copy # | geNomad | Notes |
|--------|--------|--------|---------|-------|
| 26 | 2,599 bp | **11.21** | chromosome | Homolog of S2753 contig 26 |
| 30 | 1,106 bp | **2.74** | chromosome | Different copy # from S2753 homolog |
| 31 | 1,077 bp | **10.48** | chromosome | Homolog of S2753 contig 30 |

Also detected: proviruses in contigs 1 (1.46 Mb), 2 (609 kb), 6 (316 kb).

### Split-smoothed results

Comparison of the three split-based coverage estimators for elevated contigs. The chromosomal baseline is computed assembler-agnostically using the ±10% window on geometric-mean coverages of large contigs.

#### S2052 (baseline: 216.7x, 26 chromosomal contigs)

| Contig | Length | Splits | Anvi'o mean | Arith mean | Geom mean | Median | CN (geom) | Circular |
|--------|--------|--------|-------------|------------|-----------|--------|-----------|----------|
| 9 | 223,858 bp | 11 | 347.3 | 347.0 | 346.8 | 345.5 | **1.60** | yes |
| 31 | 3,253 bp | 1 | 2163.7 | 2163.7 | 2163.7 | 2163.7 | **9.98** | |

#### S2753 (baseline: 226.9x, 19 chromosomal contigs)

| Contig | Length | Splits | Anvi'o mean | Arith mean | Geom mean | Median | CN (geom) |
|--------|--------|--------|-------------|------------|-----------|--------|-----------|
| 26 | 2,786 bp | 1 | 2547.0 | 2547.0 | 2547.0 | 2547.0 | **11.23** |
| 30 | 1,077 bp | 1 | 2351.7 | 2351.7 | 2351.7 | 2351.7 | **10.37** |

#### S2754 (baseline: 240.5x, 19 chromosomal contigs)

| Contig | Length | Splits | Anvi'o mean | Arith mean | Geom mean | Median | CN (geom) |
|--------|--------|--------|-------------|------------|-----------|--------|-----------|
| 26 | 2,599 bp | 1 | 2695.7 | 2695.7 | 2695.7 | 2695.7 | **11.21** |
| 30 | 1,106 bp | 1 | 660.1 | 660.1 | 660.1 | 660.1 | **2.74** |
| 31 | 1,077 bp | 1 | 2519.3 | 2519.3 | 2519.3 | 2519.3 | **10.48** |

In this dataset, the three estimators converge because coverage is highly uniform across splits (CV < 0.03 for all contigs). The small elevated contigs have only a single split, making the estimators identical. The difference between methods will matter more in noisier datasets (metagenomes, fragmented assemblies) where individual splits may be affected by repeat-driven multi-mapping or assembly artifacts.

## Running

### Anvi'o workflow

```bash
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

# Validate config
anvi-run-workflow -w metagenomics -c config.json --dry-run

# Submit
bsub < run_anvio_workflow.sh
```

### Extract coverage (after workflow completes)

```bash
conda activate anvio-9
for sample in S2052 S2753 S2754; do
    anvi-export-splits-and-coverages \
        -p 05_ANVIO_PROFILE/${sample}/${sample}/PROFILE.db \
        -c 03_CONTIGS/${sample}-contigs.db \
        --report-contigs -o coverage -O ${sample}_contigs
    anvi-export-splits-and-coverages \
        -p 05_ANVIO_PROFILE/${sample}/${sample}/PROFILE.db \
        -c 03_CONTIGS/${sample}-contigs.db \
        --splits-mode -o coverage -O ${sample}_splits
done
```

### Copy number and split-level analysis

```bash
python3 scripts/01_compute_copy_numbers.py      # assembler-metadata-based
python3 scripts/02_split_coverage_analysis.py    # split CV verification
python3 scripts/04_copy_numbers_split_smoothed.py  # split-smoothed, assembler-agnostic
```

### geNomad

```bash
bsub < run_genomad.sh
# After completion:
python3 scripts/03_integrate_genomad.py
```

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
| `scripts/01_compute_copy_numbers.py` | Compute copy numbers from contig coverage vs. median chromosomal baseline |
| `scripts/02_split_coverage_analysis.py` | Analyze split-level coverage for intra-contig variation |
| `scripts/03_integrate_genomad.py` | Merge geNomad classifications with copy number results |
| `scripts/04_copy_numbers_split_smoothed.py` | Split-smoothed copy numbers (arith/geom/median), assembler-agnostic baseline |

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
| `coverage/` | Exported contig and split coverage files |
| `genomad_out/` | geNomad output per sample |
| `results/` | Final tables (see below) |

## Result files

| File | Description |
|------|-------------|
| `results/copy_numbers.tsv` | Per-contig copy numbers for all samples |
| `results/copy_numbers_summary.txt` | Human-readable summary |
| `results/split_coverage_stats.tsv` | Per-contig split-level coverage statistics (mean, std, CV) |
| `results/elevated_contigs_splits.tsv` | Split-level coverage for elevated copy number contigs |
| `results/copy_numbers_with_genomad.tsv` | Integrated copy numbers + geNomad classification |
| `results/copy_numbers_split_smoothed.tsv` | Split-smoothed copy numbers (3 estimators) + assembler annotations |
| `results/copy_numbers_split_smoothed_summary.txt` | Human-readable summary of split-smoothed results |
