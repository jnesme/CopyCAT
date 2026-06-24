#!/usr/bin/env python3
"""
CopyCAT — Annotate contigs with coverage stats, HMM hits, and assembler metadata.

Usage:
    python scripts/01a_annotate_contigs.py
"""

import os
import re
import csv
import math
import statistics
import sqlite3

SAMPLES = ["S2052", "S2753", "S2754", "S2052ref", "S2753ref", "S2754ref", "S26ref"]
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFORMAT_DIR = os.path.join(BASE_DIR, "02_FASTA")
COVERAGE_DIR = os.path.join(BASE_DIR, "07_COVERAGE")
CONTIGS_DIR = os.path.join(BASE_DIR, "03_CONTIGS")
OUTPUT_DIR = os.path.join(BASE_DIR, "results")


def geometric_mean(values):
    if not values or any(v <= 0 for v in values):
        return 0.0
    log_sum = sum(math.log(v) for v in values)
    return math.exp(log_sum / len(values))


def get_hmm_hits_per_contig(sample):
    """
    Query contigs.db for per-contig HMM hit counts grouped by source.
    Returns {contig_name: {source: count}}.
    """
    db_path = os.path.join(CONTIGS_DIR, f"{sample}-contigs.db")
    if not os.path.exists(db_path):
        return {}

    conn = sqlite3.connect(db_path)
    cursor = conn.execute("""
        SELECT g.contig, h.source, COUNT(*) as n_hits
        FROM hmm_hits h
        JOIN genes_in_contigs g ON h.gene_callers_id = g.gene_callers_id
        GROUP BY g.contig, h.source
    """)
    hits = {}
    for contig, source, count in cursor:
        if contig not in hits:
            hits[contig] = {}
        hits[contig][source] = count
    conn.close()
    return hits


def get_hmm_sources(sample):
    """Get the list of HMM sources present in the contigs.db."""
    db_path = os.path.join(CONTIGS_DIR, f"{sample}-contigs.db")
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    sources = [row[0] for row in conn.execute("SELECT DISTINCT source FROM hmm_hits ORDER BY source")]
    conn.close()
    return sources


def get_contig_lengths_from_fasta(sample):
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
    path = os.path.join(COVERAGE_DIR, f"{sample}_contigs-COVs.txt")
    coverages = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            coverages[parts[0]] = float(parts[1])
    return coverages


def compute_split_smoothed_coverages(splits_by_contig):
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


def annotate_sample(sample, hmm_sources):
    contigs_meta = parse_reformat_report(sample)
    splits_by_contig = parse_split_coverage(sample)
    anvio_contig_covs = parse_contig_coverage(sample)
    smoothed = compute_split_smoothed_coverages(splits_by_contig)
    hmm_hits = get_hmm_hits_per_contig(sample)

    rows = []
    for name in sorted(smoothed.keys()):
        meta = contigs_meta.get(name, {})
        sm = smoothed[name]
        anvio_cov = anvio_contig_covs.get(name, 0)
        contig_hmm = hmm_hits.get(name, {})

        row = {
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
            "assembler_depth": meta.get("assembler_depth", ""),
            "circular": meta.get("circular", ""),
        }
        for source in hmm_sources:
            row[source] = contig_hmm.get(source, 0)

        rows.append(row)

    return rows


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_hmm_sources = set()
    for sample in SAMPLES:
        all_hmm_sources.update(get_hmm_sources(sample))
    hmm_sources = sorted(all_hmm_sources)

    all_rows = []
    for sample in SAMPLES:
        print(f"Annotating {sample}...")
        rows = annotate_sample(sample, hmm_sources)
        all_rows.extend(rows)
        print(f"  {len(rows)} contigs, HMM sources: {hmm_sources}")

    out_path = os.path.join(OUTPUT_DIR, "contig_annotations.tsv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nAnnotations written to: {out_path}")
    print(f"Total: {len(all_rows)} contigs across {len(SAMPLES)} samples")
    print(f"HMM sources: {', '.join(hmm_sources)}")


if __name__ == "__main__":
    main()
