#!/usr/bin/env python3
"""
CopyCAT — Compute copy numbers from annotated contigs.

Reads contig_annotations.tsv (from 01a, optionally enriched by 01b),
identifies the chromosomal baseline, and computes copy numbers.

Usage:
    python scripts/02_compute_copy_numbers.py
"""

import os
import csv
import statistics

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
ANNOTATIONS_TSV = os.path.join(RESULTS_DIR, "contig_annotations.tsv")

MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL = 20000
CHROMOSOMAL_PERCENT_WINDOW = 0.10
COPY_NUMBER_THRESHOLD = 1.5


def load_annotations():
    with open(ANNOTATIONS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    for row in rows:
        row["length"] = int(row["length"])
        row["split_geom_mean_cov"] = float(row["split_geom_mean_cov"])
        row["split_arith_mean_cov"] = float(row["split_arith_mean_cov"])
        row["split_median_cov"] = float(row["split_median_cov"])
        if "Bacteria_71" in row and row["Bacteria_71"]:
            row["Bacteria_71"] = int(row["Bacteria_71"])
        else:
            row["Bacteria_71"] = None
    return rows


def identify_chromosomal_contigs(sample_rows):
    large = [r for r in sample_rows
             if r["length"] >= MIN_CONTIG_LENGTH_FOR_CHROMOSOMAL]

    if not large:
        return set()

    # Step 1: anchor on largest contig
    anchor = max(large, key=lambda r: r["length"])
    anchor_cov = anchor["split_geom_mean_cov"]
    lo = anchor_cov * (1 - CHROMOSOMAL_PERCENT_WINDOW)
    hi = anchor_cov * (1 + CHROMOSOMAL_PERCENT_WINDOW)

    # Step 2: candidate pool — large contigs within ±10% of anchor
    candidates = [r for r in large if lo <= r["split_geom_mean_cov"] <= hi]

    # Step 3: positive filter — keep candidates with bacterial SCGs
    has_scg_data = any(r["Bacteria_71"] is not None for r in candidates)
    if has_scg_data:
        scg_filtered = [r for r in candidates if r["Bacteria_71"] and r["Bacteria_71"] > 0]
        if scg_filtered:
            candidates = scg_filtered

    # Step 4: negative filter — exclude geNomad-classified plasmids
    if any("genomad_class" in r for r in candidates):
        filtered = [r for r in candidates if r.get("genomad_class") != "plasmid"]
        if filtered:
            candidates = filtered

    return {r["contig"] for r in candidates}


def compute_for_sample(sample_rows):
    chromosomal = identify_chromosomal_contigs(sample_rows)
    chr_covs = [r["split_geom_mean_cov"] for r in sample_rows if r["contig"] in chromosomal]
    baseline = statistics.median(chr_covs)

    results = []
    for r in sample_rows:
        copy_number = r["split_geom_mean_cov"] / baseline
        is_chromosomal = r["contig"] in chromosomal
        is_elevated = copy_number >= COPY_NUMBER_THRESHOLD

        out = dict(r)
        out["chromosomal_baseline"] = round(baseline, 2)
        out["copy_number"] = round(copy_number, 2)
        out["classification"] = (
            "chromosomal" if is_chromosomal
            else "putative_plasmid" if is_elevated
            else "uncertain"
        )
        results.append(out)

    return results, baseline, len(chromosomal)


def main():
    all_rows = load_annotations()

    samples = []
    seen = set()
    for r in all_rows:
        if r["sample"] not in seen:
            samples.append(r["sample"])
            seen.add(r["sample"])

    all_results = []
    summary_lines = []

    for sample in samples:
        sample_rows = [r for r in all_rows if r["sample"] == sample]
        results, baseline, n_chr = compute_for_sample(sample_rows)
        all_results.extend(results)

        elevated = [r for r in results if r["classification"] == "putative_plasmid"]
        summary_lines.append(f"\n{'='*90}")
        summary_lines.append(f"Sample: {sample}")
        summary_lines.append(f"{'='*90}")
        summary_lines.append(f"  Total contigs:               {len(results)}")
        summary_lines.append(f"  Chromosomal contigs (±{CHROMOSOMAL_PERCENT_WINDOW:.0%}): {n_chr}")
        summary_lines.append(f"  Chromosomal baseline (geom): {baseline:.2f}x")
        summary_lines.append(f"  Elevated (>={COPY_NUMBER_THRESHOLD}x):          {len(elevated)}")
        summary_lines.append("")
        summary_lines.append(
            f"  {'Contig':<25} {'Length':>8} {'GeomCov':>8} "
            f"{'CN':>7} {'Circ':>5} {'Class'}"
        )
        summary_lines.append(
            f"  {'-'*25} {'-'*8} {'-'*8} "
            f"{'-'*7} {'-'*5} {'-'*16}"
        )
        for r in results:
            circ_str = "yes" if r.get("circular") == "True" else ""
            summary_lines.append(
                f"  {r['contig']:<25} {r['length']:>8,} "
                f"{r['split_geom_mean_cov']:>8.1f} "
                f"{r['copy_number']:>7.2f} {circ_str:>5} {r['classification']}"
            )

    out_tsv = os.path.join(RESULTS_DIR, "copy_numbers.tsv")
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_results)

    out_summary = os.path.join(RESULTS_DIR, "copy_numbers_summary.txt")
    summary_text = "\n".join(summary_lines)
    with open(out_summary, "w") as f:
        f.write(summary_text + "\n")

    print(summary_text)
    print(f"\nResults written to: {out_tsv}")
    print(f"Summary written to: {out_summary}")


if __name__ == "__main__":
    main()
