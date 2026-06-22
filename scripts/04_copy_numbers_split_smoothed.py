#!/usr/bin/env python3
"""
Compute plasmid copy numbers using split-smoothed coverage estimates.

Differences from 01_compute_copy_numbers.py:
  - Contig coverage is derived from splits (arithmetic mean, geometric mean, median)
    rather than using anvi'o's position-level mean directly
  - Chromosomal baseline is identified purely from coverage distribution
    (no assembler metadata required)
  - Assembler metadata (depth, circularity) is included as annotation columns
    but does not drive the classification

Chromosomal identification:
  Large contigs (>= 20 kb) form the candidate pool. The median of their
  geometric-mean coverages defines the center. Contigs within a percentage
  window (default ±10%) of this median are classified as chromosomal.
  This is assembler-agnostic and robust across tight and noisy distributions.

Usage:
    python scripts/04_copy_numbers_split_smoothed.py
"""

import os
import re
import csv
import math
import statistics
import sqlite3

SAMPLES = ["S2052", "S2753", "S2754", "S2052ref"]
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFORMAT_DIR = os.path.join(BASE_DIR, "02_FASTA")
COVERAGE_DIR = os.path.join(BASE_DIR, "07_COVERAGE")
CONTIGS_DIR = os.path.join(BASE_DIR, "03_CONTIGS")
OUTPUT_DIR = os.path.join(BASE_DIR, "results")

MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL = 20000
CHROMOSOMAL_PERCENT_WINDOW = 0.10
COPY_NUMBER_THRESHOLD = 1.5


def geometric_mean(values):
    if not values or any(v <= 0 for v in values):
        return 0.0
    log_sum = sum(math.log(v) for v in values)
    return math.exp(log_sum / len(values))


def get_rrna_contigs(sample):
    """
    Query the anvi'o contigs database to find contigs containing rRNA genes.
    Joins hmm_hits (source like 'Ribosomal_RNA%') -> genes_in_contigs -> contig name.
    Returns a set of contig names that contain rRNA hits.
    """
    db_path = os.path.join(CONTIGS_DIR, f"{sample}-contigs.db")
    if not os.path.exists(db_path):
        return set()

    conn = sqlite3.connect(db_path)
    cursor = conn.execute("""
        SELECT DISTINCT g.contig
        FROM hmm_hits h
        JOIN genes_in_contigs g ON h.gene_callers_id = g.gene_callers_id
        WHERE h.source LIKE 'Ribosomal_RNA%'
    """)
    rrna_contigs = {row[0] for row in cursor}
    conn.close()
    return rrna_contigs


def get_contig_lengths_from_fasta(sample):
    """Get contig lengths from the exported CONTIGS FASTA."""
    fasta_path = os.path.join(COVERAGE_DIR, f"{sample}_contigs-CONTIGS.fa")
    lengths = {}
    if not os.path.exists(fasta_path):
        return lengths
    name = None
    seq_len = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].strip().split()[0]
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if name is not None:
            lengths[name] = seq_len
    return lengths


def parse_reformat_report(sample):
    """Parse the reformat report to extract assembler metadata (optional annotations).
    Falls back to FASTA-derived lengths when headers lack length= tags."""
    path = os.path.join(REFORMAT_DIR, sample, f"{sample}-reformat-report.txt")
    contigs = {}
    if not os.path.exists(path):
        return contigs

    fasta_lengths = get_contig_lengths_from_fasta(sample)

    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            new_name = parts[0]
            original_header = parts[1]

            length_match = re.search(r"length=(\d+)", original_header)
            depth_match = re.search(r"depth=([\d.]+)x", original_header)
            circular = "circular=true" in original_header

            length = int(length_match.group(1)) if length_match else fasta_lengths.get(new_name, 0)

            contigs[new_name] = {
                "original_id": original_header.split()[0],
                "length": length,
                "assembler_depth": float(depth_match.group(1)) if depth_match else None,
                "circular": circular,
            }
    return contigs


def parse_split_coverage(sample):
    """Parse split-level coverage and group by parent contig."""
    path = os.path.join(COVERAGE_DIR, f"{sample}_splits-COVs.txt")
    splits_by_contig = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            split_name = parts[0]
            cov = float(parts[1])
            contig = re.sub(r"_split_\d+$", "", split_name)
            if contig not in splits_by_contig:
                splits_by_contig[contig] = []
            splits_by_contig[contig].append(cov)
    return splits_by_contig


def parse_contig_coverage(sample):
    """Parse anvi'o contig-level coverage (position-level mean) for comparison."""
    path = os.path.join(COVERAGE_DIR, f"{sample}_contigs-COVs.txt")
    coverages = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            coverages[parts[0]] = float(parts[1])
    return coverages


def compute_split_smoothed_coverages(splits_by_contig):
    """Compute three coverage estimates per contig from splits."""
    results = {}
    for contig, split_covs in splits_by_contig.items():
        results[contig] = {
            "n_splits": len(split_covs),
            "arith_mean": statistics.mean(split_covs),
            "geom_mean": geometric_mean(split_covs),
            "median": statistics.median(split_covs),
            "cv": (statistics.stdev(split_covs) / statistics.mean(split_covs)
                   if len(split_covs) > 1 and statistics.mean(split_covs) > 0
                   else 0.0),
        }
    return results


def identify_chromosomal_contigs(contigs_meta, smoothed_covs):
    """
    Identify chromosomal contigs from coverage distribution alone.
    Uses large contigs (>= 20kb) and selects those within a percentage window
    of the median coverage. This is robust across clean isolates and noisier
    datasets without being over-sensitive to tight distributions.
    """
    large_contigs = {}
    for name, meta in contigs_meta.items():
        if meta["length"] >= MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL and name in smoothed_covs:
            large_contigs[name] = smoothed_covs[name]["geom_mean"]

    if not large_contigs:
        return []

    cov_values = list(large_contigs.values())
    med = statistics.median(cov_values)
    lo = med * (1 - CHROMOSOMAL_PERCENT_WINDOW)
    hi = med * (1 + CHROMOSOMAL_PERCENT_WINDOW)

    chromosomal = [
        name for name, cov in large_contigs.items()
        if lo <= cov <= hi
    ]
    return chromosomal


def compute_copy_numbers(sample):
    """Compute copy numbers using split-smoothed coverage."""
    contigs_meta = parse_reformat_report(sample)
    splits_by_contig = parse_split_coverage(sample)
    anvio_contig_covs = parse_contig_coverage(sample)
    smoothed = compute_split_smoothed_coverages(splits_by_contig)
    rrna_contigs = get_rrna_contigs(sample)

    chromosomal = identify_chromosomal_contigs(contigs_meta, smoothed)
    chr_covs = [smoothed[c]["geom_mean"] for c in chromosomal]
    baseline = statistics.median(chr_covs)

    results = []
    for name in sorted(smoothed.keys()):
        meta = contigs_meta.get(name, {})
        sm = smoothed[name]
        anvio_cov = anvio_contig_covs.get(name, 0)

        cn_arith = sm["arith_mean"] / baseline
        cn_geom = sm["geom_mean"] / baseline
        cn_median = sm["median"] / baseline
        is_chromosomal = name in chromosomal
        has_rrna = name in rrna_contigs
        is_elevated = cn_geom >= COPY_NUMBER_THRESHOLD

        results.append({
            "sample": sample,
            "contig": name,
            "original_id": meta.get("original_id", ""),
            "length": meta.get("length", 0),
            "n_splits": sm["n_splits"],
            "split_cv": round(sm["cv"], 4),
            "anvio_mean_cov": round(anvio_cov, 2),
            "split_arith_mean_cov": round(sm["arith_mean"], 2),
            "split_geom_mean_cov": round(sm["geom_mean"], 2),
            "split_median_cov": round(sm["median"], 2),
            "chromosomal_baseline": round(baseline, 2),
            "cn_arith_mean": round(cn_arith, 2),
            "cn_geom_mean": round(cn_geom, 2),
            "cn_median": round(cn_median, 2),
            "has_rrna": has_rrna,
            "classification": "chromosomal" if is_chromosomal else (
                "putative_plasmid" if is_elevated else "uncertain"
            ),
            "assembler_depth": meta.get("assembler_depth", ""),
            "circular": meta.get("circular", ""),
        })

    return results, baseline, len(chromosomal)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_results = []
    summary_lines = []

    for sample in SAMPLES:
        results, baseline, n_chr = compute_copy_numbers(sample)
        all_results.extend(results)

        elevated = [r for r in results if r["classification"] == "putative_plasmid"]
        summary_lines.append(f"\n{'='*90}")
        summary_lines.append(f"Sample: {sample}")
        summary_lines.append(f"{'='*90}")
        summary_lines.append(f"  Total contigs profiled:        {len(results)}")
        summary_lines.append(f"  Chromosomal contigs (±{CHROMOSOMAL_PERCENT_WINDOW:.0%}):  {n_chr}")
        summary_lines.append(f"  Chromosomal baseline (geom):   {baseline:.2f}x")
        summary_lines.append(f"  Elevated copy number (>={COPY_NUMBER_THRESHOLD}x): {len(elevated)}")
        summary_lines.append("")
        summary_lines.append(
            f"  {'Contig':<25} {'Length':>8} {'Splits':>6} "
            f"{'Anvi.o':>8} {'Arith':>8} {'Geom':>8} {'Median':>8} "
            f"{'CN_geo':>7} {'rRNA':>5} {'Circ':>5} {'Class'}"
        )
        summary_lines.append(
            f"  {'-'*25} {'-'*8} {'-'*6} "
            f"{'-'*8} {'-'*8} {'-'*8} {'-'*8} "
            f"{'-'*7} {'-'*5} {'-'*5} {'-'*16}"
        )
        for r in results:
            circ_str = "yes" if r["circular"] is True else ""
            rrna_str = "yes" if r["has_rrna"] else ""
            summary_lines.append(
                f"  {r['contig']:<25} {r['length']:>8,} {r['n_splits']:>6} "
                f"{r['anvio_mean_cov']:>8.1f} {r['split_arith_mean_cov']:>8.1f} "
                f"{r['split_geom_mean_cov']:>8.1f} {r['split_median_cov']:>8.1f} "
                f"{r['cn_geom_mean']:>7.2f} {rrna_str:>5} {circ_str:>5} {r['classification']}"
            )

    out_tsv = os.path.join(OUTPUT_DIR, "copy_numbers_split_smoothed.tsv")
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_results)

    out_summary = os.path.join(OUTPUT_DIR, "copy_numbers_split_smoothed_summary.txt")
    summary_text = "\n".join(summary_lines)
    with open(out_summary, "w") as f:
        f.write(summary_text + "\n")

    print(summary_text)
    print(f"\nFull results written to: {out_tsv}")
    print(f"Summary written to:     {out_summary}")


if __name__ == "__main__":
    main()
