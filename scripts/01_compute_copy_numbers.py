#!/usr/bin/env python3
"""
Compute plasmid copy numbers from anvi'o coverage data.

For each sample:
  1. Parse the reformat report to get contig metadata (length, assembler depth, circularity)
  2. Parse the contig-level coverage from anvi'o profiling
  3. Identify chromosomal contigs (large, single-copy based on assembler depth ~1x)
  4. Compute median chromosomal coverage as the baseline
  5. Compute copy number = contig_coverage / median_chromosomal_coverage
  6. Flag contigs with elevated copy number or circularity as putative plasmids

Usage:
    python scripts/01_compute_copy_numbers.py
"""

import os
import re
import csv
import statistics

SAMPLES = ["S2052", "S2753", "S2754"]
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFORMAT_DIR = os.path.join(BASE_DIR, "02_FASTA")
COVERAGE_DIR = os.path.join(BASE_DIR, "07_COVERAGE")
OUTPUT_DIR = os.path.join(BASE_DIR, "results")

MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL = 20000
ASSEMBLER_DEPTH_RANGE_FOR_CHROMOSOMAL = (0.8, 1.3)
COPY_NUMBER_THRESHOLD = 1.5


def parse_reformat_report(sample):
    """Parse the reformat report to extract contig metadata."""
    path = os.path.join(REFORMAT_DIR, sample, f"{sample}-reformat-report.txt")
    contigs = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            new_name = parts[0]
            original_header = parts[1]

            length_match = re.search(r"length=(\d+)", original_header)
            depth_match = re.search(r"depth=([\d.]+)x", original_header)
            circular = "circular=true" in original_header

            contigs[new_name] = {
                "original_id": original_header.split()[0],
                "length": int(length_match.group(1)) if length_match else 0,
                "assembler_depth": float(depth_match.group(1)) if depth_match else 0,
                "circular": circular,
            }
    return contigs


def parse_coverage(sample):
    """Parse the contig-level coverage file from anvi'o."""
    path = os.path.join(COVERAGE_DIR, f"{sample}_contigs-COVs.txt")
    coverages = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            coverages[parts[0]] = float(parts[1])
    return coverages


def identify_chromosomal_contigs(contigs, coverages):
    """
    Identify chromosomal contigs based on:
    - Length >= MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL
    - Assembler depth within ASSEMBLER_DEPTH_RANGE_FOR_CHROMOSOMAL (suggesting single-copy)
    - Not marked as circular
    """
    chromosomal = []
    lo, hi = ASSEMBLER_DEPTH_RANGE_FOR_CHROMOSOMAL
    for name, meta in contigs.items():
        if name not in coverages:
            continue
        if (meta["length"] >= MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL
                and lo <= meta["assembler_depth"] <= hi
                and not meta["circular"]):
            chromosomal.append(name)
    return chromosomal


def compute_copy_numbers(sample):
    """Compute copy numbers for all contigs in a sample."""
    contigs = parse_reformat_report(sample)
    coverages = parse_coverage(sample)

    chromosomal = identify_chromosomal_contigs(contigs, coverages)
    chromosomal_coverages = [coverages[c] for c in chromosomal]
    median_chr_cov = statistics.median(chromosomal_coverages)

    results = []
    for name in sorted(coverages.keys()):
        meta = contigs.get(name, {})
        cov = coverages[name]
        copy_number = cov / median_chr_cov
        is_chromosomal = name in chromosomal
        is_elevated = copy_number >= COPY_NUMBER_THRESHOLD

        results.append({
            "sample": sample,
            "contig": name,
            "original_id": meta.get("original_id", ""),
            "length": meta.get("length", 0),
            "circular": meta.get("circular", False),
            "assembler_depth": meta.get("assembler_depth", 0),
            "mapping_coverage": round(cov, 2),
            "median_chromosomal_cov": round(median_chr_cov, 2),
            "copy_number": round(copy_number, 2),
            "classification": "chromosomal" if is_chromosomal else (
                "putative_plasmid" if is_elevated else "uncertain"
            ),
        })

    return results, median_chr_cov, len(chromosomal), len(chromosomal_coverages)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_results = []
    summary_lines = []

    for sample in SAMPLES:
        results, median_cov, n_chr, n_chr_with_cov = compute_copy_numbers(sample)
        all_results.extend(results)

        elevated = [r for r in results if r["classification"] == "putative_plasmid"]
        summary_lines.append(f"\n{'='*70}")
        summary_lines.append(f"Sample: {sample}")
        summary_lines.append(f"{'='*70}")
        summary_lines.append(f"  Total contigs profiled:     {len(results)}")
        summary_lines.append(f"  Chromosomal contigs:        {n_chr_with_cov}")
        summary_lines.append(f"  Median chromosomal coverage: {median_cov:.2f}x")
        summary_lines.append(f"  Elevated copy number (>={COPY_NUMBER_THRESHOLD}x): {len(elevated)}")
        summary_lines.append("")
        summary_lines.append(f"  {'Contig':<25} {'Length':>10} {'Circ':>5} {'Coverage':>10} {'CopyNum':>8} {'Class'}")
        summary_lines.append(f"  {'-'*25} {'-'*10} {'-'*5} {'-'*10} {'-'*8} {'-'*16}")
        for r in results:
            flag = "*" if r["classification"] == "putative_plasmid" else ""
            summary_lines.append(
                f"  {r['contig']:<25} {r['length']:>10,} {str(r['circular']):>5} "
                f"{r['mapping_coverage']:>10.1f} {r['copy_number']:>8.2f} {r['classification']}{flag}"
            )

    out_tsv = os.path.join(OUTPUT_DIR, "copy_numbers.tsv")
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_results)

    out_summary = os.path.join(OUTPUT_DIR, "copy_numbers_summary.txt")
    summary_text = "\n".join(summary_lines)
    with open(out_summary, "w") as f:
        f.write(summary_text + "\n")

    print(summary_text)
    print(f"\nFull results written to: {out_tsv}")
    print(f"Summary written to:     {out_summary}")


if __name__ == "__main__":
    main()
