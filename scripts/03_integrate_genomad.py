#!/usr/bin/env python3
"""
Integrate geNomad plasmid/virus predictions with coverage-based copy numbers.

Parses geNomad summary files for each sample, merges with the existing
copy_numbers.tsv, and outputs a combined table with both coverage-based
and sequence-based classification.

Usage:
    python scripts/03_integrate_genomad.py
"""

import os
import csv

SAMPLES = ["S2052", "S2753", "S2754", "S2052ref"]
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMAD_DIR = os.path.join(BASE_DIR, "genomad_out")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
COPY_NUMBERS_TSV = os.path.join(RESULTS_DIR, "copy_numbers.tsv")


def parse_genomad_summary(sample, summary_type):
    """
    Parse a geNomad plasmid or virus summary TSV.
    Returns dict of contig_name -> {score, n_hallmarks, ...}
    """
    prefix = f"{sample}_contigs-CONTIGS"
    summary_path = os.path.join(
        GENOMAD_DIR, sample, f"{prefix}_summary",
        f"{prefix}_{summary_type}_summary.tsv"
    )
    results = {}
    if not os.path.exists(summary_path):
        return results

    with open(summary_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq_name = row.get("seq_name", "")
            # geNomad appends |provirus to provirus names
            base_contig = seq_name.split("|")[0]

            score_key = "plasmid_score" if summary_type == "plasmid" else "virus_score"
            score = float(row.get(score_key, 0))
            n_hallmarks = int(row.get("n_hallmarks", 0))

            if base_contig not in results or score > results[base_contig]["score"]:
                results[base_contig] = {
                    "seq_name": seq_name,
                    "score": score,
                    "n_hallmarks": n_hallmarks,
                    "topology": row.get("topology", ""),
                    "length": int(row.get("length", 0)),
                }
    return results


def load_copy_numbers():
    """Load the existing copy_numbers.tsv."""
    rows = []
    with open(COPY_NUMBERS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def integrate():
    """Merge geNomad results with copy number data."""
    copy_numbers = load_copy_numbers()

    genomad_data = {}
    for sample in SAMPLES:
        plasmids = parse_genomad_summary(sample, "plasmid")
        viruses = parse_genomad_summary(sample, "virus")
        genomad_data[sample] = {"plasmid": plasmids, "virus": viruses}

    integrated = []
    for row in copy_numbers:
        sample = row["sample"]
        contig = row["contig"]

        plasmid_info = genomad_data.get(sample, {}).get("plasmid", {}).get(contig)
        virus_info = genomad_data.get(sample, {}).get("virus", {}).get(contig)

        if plasmid_info and virus_info:
            if plasmid_info["score"] >= virus_info["score"]:
                genomad_class = "plasmid"
                genomad_score = plasmid_info["score"]
                genomad_hallmarks = plasmid_info["n_hallmarks"]
                genomad_topology = plasmid_info["topology"]
            else:
                genomad_class = "virus"
                genomad_score = virus_info["score"]
                genomad_hallmarks = virus_info["n_hallmarks"]
                genomad_topology = virus_info["topology"]
        elif plasmid_info:
            genomad_class = "plasmid"
            genomad_score = plasmid_info["score"]
            genomad_hallmarks = plasmid_info["n_hallmarks"]
            genomad_topology = plasmid_info["topology"]
        elif virus_info:
            genomad_class = "virus"
            genomad_score = virus_info["score"]
            genomad_hallmarks = virus_info["n_hallmarks"]
            genomad_topology = virus_info["topology"]
        else:
            genomad_class = "chromosome"
            genomad_score = ""
            genomad_hallmarks = ""
            genomad_topology = ""

        row["genomad_class"] = genomad_class
        row["genomad_score"] = genomad_score
        row["genomad_n_hallmarks"] = genomad_hallmarks
        row["genomad_topology"] = genomad_topology
        integrated.append(row)

    return integrated


def print_summary(integrated):
    """Print a combined summary."""
    for sample in SAMPLES:
        rows = [r for r in integrated if r["sample"] == sample]
        print(f"\n{'='*80}")
        print(f"Sample: {sample}")
        print(f"{'='*80}")
        print(f"  {'Contig':<25} {'Length':>8} {'CopyNum':>8} {'Coverage-class':<18} {'geNomad-class':<12} {'Score':>6} {'Hallmarks':>9} {'Topology'}")
        print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*18} {'-'*12} {'-'*6} {'-'*9} {'-'*10}")
        for r in rows:
            score_str = f"{float(r['genomad_score']):.3f}" if r["genomad_score"] != "" else ""
            hallmark_str = str(r["genomad_n_hallmarks"]) if r["genomad_n_hallmarks"] != "" else ""
            topo_str = r.get("genomad_topology", "")
            print(
                f"  {r['contig']:<25} {int(r['length']):>8,} {float(r['copy_number']):>8.2f} "
                f"{r['classification']:<18} {r['genomad_class']:<12} {score_str:>6} {hallmark_str:>9} {topo_str}"
            )

        # Agreement check
        elevated_cov = [r for r in rows if r["classification"] == "putative_plasmid"]
        genomad_plasmids = [r for r in rows if r["genomad_class"] == "plasmid"]
        both = [r for r in rows if r["classification"] == "putative_plasmid" and r["genomad_class"] == "plasmid"]

        print(f"\n  Coverage-based putative plasmids: {len(elevated_cov)}")
        print(f"  geNomad-classified plasmids:      {len(genomad_plasmids)}")
        print(f"  Agreed by both methods:           {len(both)}")

        cov_only = [r for r in rows if r["classification"] == "putative_plasmid" and r["genomad_class"] != "plasmid"]
        genomad_only = [r for r in rows if r["classification"] != "putative_plasmid" and r["genomad_class"] == "plasmid"]
        if cov_only:
            print(f"  Coverage-only (not confirmed by geNomad): {', '.join(r['contig'] for r in cov_only)}")
        if genomad_only:
            print(f"  geNomad-only (no elevated copy number):   {', '.join(r['contig'] for r in genomad_only)}")


def main():
    integrated = integrate()

    out_path = os.path.join(RESULTS_DIR, "copy_numbers_with_genomad.tsv")
    fieldnames = list(integrated[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(integrated)

    print_summary(integrated)
    print(f"\nIntegrated results written to: {out_path}")


if __name__ == "__main__":
    main()
