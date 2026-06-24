#!/usr/bin/env python3
"""
Test sensitivity of chromosomal baseline to window size.

Sweeps the ±window from 5% to 20% and reports how the baseline,
number of chromosomal contigs, and key copy numbers change.

Usage:
    python scripts/test_baseline_window.py
"""

import os
import csv
import statistics

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
ANNOTATIONS_TSV = os.path.join(RESULTS_DIR, "contig_annotations.tsv")

MIN_CONTIG_LENGTH = 20000
WINDOWS = [0.05, 0.08, 0.10, 0.12, 0.15, 0.20]
CN_THRESHOLD = 1.5


def load_annotations():
    with open(ANNOTATIONS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    for row in rows:
        row["length"] = int(row["length"])
        row["split_geom_mean_cov"] = float(row["split_geom_mean_cov"])
        if "Bacteria_71" in row and row["Bacteria_71"]:
            row["Bacteria_71"] = int(row["Bacteria_71"])
        else:
            row["Bacteria_71"] = None
    return rows


def identify_chromosomal(sample_rows, window):
    large = [r for r in sample_rows if r["length"] >= MIN_CONTIG_LENGTH]
    if not large:
        return set(), 0.0

    anchor = max(large, key=lambda r: r["length"])
    anchor_cov = anchor["split_geom_mean_cov"]
    lo = anchor_cov * (1 - window)
    hi = anchor_cov * (1 + window)

    candidates = [r for r in large if lo <= r["split_geom_mean_cov"] <= hi]

    has_scg_data = any(r["Bacteria_71"] is not None for r in candidates)
    if has_scg_data:
        scg_filtered = [r for r in candidates if r["Bacteria_71"] and r["Bacteria_71"] > 0]
        if scg_filtered:
            candidates = scg_filtered

    if any("genomad_class" in r for r in candidates):
        filtered = [r for r in candidates if r.get("genomad_class") != "plasmid"]
        if filtered:
            candidates = filtered

    chr_set = {r["contig"] for r in candidates}
    baseline = statistics.median([r["split_geom_mean_cov"] for r in candidates])
    return chr_set, baseline


def main():
    all_rows = load_annotations()

    samples = []
    seen = set()
    for r in all_rows:
        if r["sample"] not in seen:
            samples.append(r["sample"])
            seen.add(r["sample"])

    for sample in samples:
        sample_rows = [r for r in all_rows if r["sample"] == sample]

        # Find elevated contigs to track across windows
        large = [r for r in sample_rows if r["length"] >= MIN_CONTIG_LENGTH]
        if not large:
            continue
        anchor = max(large, key=lambda r: r["length"])

        print(f"\n{'='*90}")
        print(f"Sample: {sample}  |  Anchor: {anchor['contig']} ({anchor['length']:,} bp, {anchor['split_geom_mean_cov']:.1f}x)")
        print(f"{'='*90}")
        print(f"  {'Window':>8} {'N_chr':>6} {'Baseline':>10} ", end="")

        # Collect contigs with notable CN at any window
        notable = {}
        for r in sample_rows:
            for w in WINDOWS:
                _, bl = identify_chromosomal(sample_rows, w)
                cn = r["split_geom_mean_cov"] / bl if bl > 0 else 0
                if cn >= CN_THRESHOLD or r.get("genomad_class") == "plasmid":
                    notable[r["contig"]] = r
                    break

        for contig in sorted(notable):
            r = notable[contig]
            label = f"{contig}({r['length']//1000}k)"
            print(f" {label:>18}", end="")
        print()
        print(f"  {'-'*8} {'-'*6} {'-'*10} ", end="")
        for _ in notable:
            print(f" {'-'*18}", end="")
        print()

        for w in WINDOWS:
            chr_set, baseline = identify_chromosomal(sample_rows, w)
            print(f"  ±{w:>5.0%} {len(chr_set):>6} {baseline:>10.2f}x", end="")
            for contig in sorted(notable):
                r = notable[contig]
                cn = r["split_geom_mean_cov"] / baseline if baseline > 0 else 0
                in_baseline = "*" if contig in chr_set else " "
                print(f" {cn:>16.2f}{in_baseline}", end="")
            print()


if __name__ == "__main__":
    main()
