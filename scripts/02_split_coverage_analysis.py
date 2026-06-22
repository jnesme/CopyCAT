#!/usr/bin/env python3
"""
Analyze split-level coverage to detect intra-contig coverage variation.

For each contig, this script:
  1. Collects coverage of all splits (20kb windows)
  2. Computes mean, median, std, CV (coefficient of variation), min, max
  3. Flags contigs with high CV (>0.15) as having heterogeneous coverage
     (potential chimeras, partial duplications, or mosaic regions)
  4. Outputs per-split coverage for elevated-copy-number contigs for plotting

Usage:
    python scripts/02_split_coverage_analysis.py
"""

import os
import re
import csv
import statistics
import math

SAMPLES = ["S2052", "S2753", "S2754", "S2052ref"]
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COVERAGE_DIR = os.path.join(BASE_DIR, "07_COVERAGE")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

CV_THRESHOLD = 0.15
COPY_NUMBER_THRESHOLD = 1.5


def parse_contig_coverage(sample):
    """Parse contig-level coverage to get median chromosomal baseline."""
    path = os.path.join(COVERAGE_DIR, f"{sample}_contigs-COVs.txt")
    coverages = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            coverages[parts[0]] = float(parts[1])
    return coverages


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
            splits_by_contig[contig].append({"split": split_name, "coverage": cov})
    return splits_by_contig


def analyze_sample(sample):
    """Analyze split-level coverage for one sample."""
    contig_covs = parse_contig_coverage(sample)
    splits_by_contig = parse_split_coverage(sample)

    contig_stats = []
    for contig in sorted(splits_by_contig.keys()):
        splits = splits_by_contig[contig]
        covs = [s["coverage"] for s in splits]
        n_splits = len(covs)
        mean_cov = statistics.mean(covs)
        median_cov = statistics.median(covs)
        std_cov = statistics.stdev(covs) if n_splits > 1 else 0.0
        cv = std_cov / mean_cov if mean_cov > 0 else 0.0
        min_cov = min(covs)
        max_cov = max(covs)

        contig_stats.append({
            "sample": sample,
            "contig": contig,
            "n_splits": n_splits,
            "contig_mean_cov": round(mean_cov, 2),
            "contig_median_cov": round(median_cov, 2),
            "contig_std_cov": round(std_cov, 2),
            "contig_cv": round(cv, 4),
            "min_split_cov": round(min_cov, 2),
            "max_split_cov": round(max_cov, 2),
            "heterogeneous": cv > CV_THRESHOLD,
        })

    return contig_stats, splits_by_contig


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    all_stats = []
    all_elevated_splits = []

    for sample in SAMPLES:
        contig_covs = parse_contig_coverage(sample)
        cov_values = list(contig_covs.values())
        median_chr = statistics.median(cov_values)

        contig_stats, splits_by_contig = analyze_sample(sample)
        all_stats.extend(contig_stats)

        print(f"\n{'='*70}")
        print(f"Sample: {sample} — Split-level coverage analysis")
        print(f"{'='*70}")
        print(f"  Median contig coverage: {median_chr:.2f}x")
        print()

        heterogeneous = [s for s in contig_stats if s["heterogeneous"]]
        if heterogeneous:
            print(f"  Contigs with heterogeneous coverage (CV > {CV_THRESHOLD}):")
            print(f"  {'Contig':<25} {'Splits':>6} {'Mean':>8} {'Std':>8} {'CV':>6} {'Min':>8} {'Max':>8}")
            print(f"  {'-'*25} {'-'*6} {'-'*8} {'-'*8} {'-'*6} {'-'*8} {'-'*8}")
            for s in heterogeneous:
                print(f"  {s['contig']:<25} {s['n_splits']:>6} {s['contig_mean_cov']:>8.1f} "
                      f"{s['contig_std_cov']:>8.1f} {s['contig_cv']:>6.3f} "
                      f"{s['min_split_cov']:>8.1f} {s['max_split_cov']:>8.1f}")
        else:
            print("  No contigs with heterogeneous coverage detected.")

        for contig, cov in contig_covs.items():
            copy_num = cov / median_chr
            if copy_num >= COPY_NUMBER_THRESHOLD and contig in splits_by_contig:
                for split_info in splits_by_contig[contig]:
                    split_num = int(re.search(r"split_(\d+)$", split_info["split"]).group(1))
                    all_elevated_splits.append({
                        "sample": sample,
                        "contig": contig,
                        "copy_number": round(copy_num, 2),
                        "split": split_info["split"],
                        "split_num": split_num,
                        "split_coverage": round(split_info["coverage"], 2),
                        "median_chromosomal_cov": round(median_chr, 2),
                    })

    stats_path = os.path.join(RESULTS_DIR, "split_coverage_stats.tsv")
    with open(stats_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_stats[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_stats)
    print(f"\nSplit statistics written to: {stats_path}")

    if all_elevated_splits:
        splits_path = os.path.join(RESULTS_DIR, "elevated_contigs_splits.tsv")
        with open(splits_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_elevated_splits[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(all_elevated_splits)
        print(f"Elevated contig splits written to: {splits_path}")


if __name__ == "__main__":
    main()
