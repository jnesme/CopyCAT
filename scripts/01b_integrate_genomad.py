#!/usr/bin/env python3
"""
CopyCAT — Integrate geNomad scores into contig annotations.

Reads contig_annotations.tsv (from 01a), joins geNomad aggregated
classification scores (plasmid_score, virus_score) for all contigs,
plus hallmarks and topology from the filtered summaries.
Optional step — if genomad_out/ doesn't exist, exits cleanly.

Usage:
    python scripts/01b_integrate_genomad.py
"""

import os
import sys
import csv

SAMPLES = ["S2052", "S2753", "S2754", "S2052ref", "S2753ref", "S2754ref", "S26ref"]
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMAD_DIR = os.path.join(BASE_DIR, "genomad_out")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
ANNOTATIONS_TSV = os.path.join(RESULTS_DIR, "contig_annotations.tsv")


def parse_aggregated_classification(sample):
    """Read aggregated scores (plasmid_score, virus_score) for all contigs."""
    prefix = f"{sample}_contigs-CONTIGS"
    path = os.path.join(
        GENOMAD_DIR, sample,
        f"{prefix}_aggregated_classification",
        f"{prefix}_aggregated_classification.tsv"
    )
    scores = {}
    if not os.path.exists(path):
        return scores

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            scores[row["seq_name"]] = {
                "plasmid_score": float(row["plasmid_score"]),
                "virus_score": float(row["virus_score"]),
            }
    return scores


def parse_summary(sample, summary_type):
    """Read hallmarks and topology from filtered plasmid/virus summaries."""
    prefix = f"{sample}_contigs-CONTIGS"
    path = os.path.join(
        GENOMAD_DIR, sample, f"{prefix}_summary",
        f"{prefix}_{summary_type}_summary.tsv"
    )
    results = {}
    if not os.path.exists(path):
        return results

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            base_contig = row["seq_name"].split("|")[0]
            n_hallmarks = int(row.get("n_hallmarks", 0))
            if base_contig not in results or n_hallmarks > results[base_contig]["n_hallmarks"]:
                results[base_contig] = {
                    "n_hallmarks": n_hallmarks,
                    "topology": row.get("topology", ""),
                }
    return results


def main():
    if not os.path.isdir(GENOMAD_DIR):
        print(f"geNomad output not found ({GENOMAD_DIR}), skipping.")
        sys.exit(0)

    if not os.path.exists(ANNOTATIONS_TSV):
        print(f"Error: {ANNOTATIONS_TSV} not found. Run 01a_annotate_contigs.py first.")
        sys.exit(1)

    with open(ANNOTATIONS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    genomad_data = {}
    n_found = 0
    for sample in SAMPLES:
        agg = parse_aggregated_classification(sample)
        plasmid_summary = parse_summary(sample, "plasmid")
        virus_summary = parse_summary(sample, "virus")
        if agg:
            n_found += 1
        genomad_data[sample] = {
            "scores": agg,
            "plasmid": plasmid_summary,
            "virus": virus_summary,
        }

    if n_found == 0:
        print("No geNomad results found for any sample, skipping.")
        sys.exit(0)

    new_cols = [
        "genomad_plasmid_score", "genomad_virus_score",
        "genomad_class", "genomad_n_hallmarks", "genomad_topology",
    ]
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    for row in rows:
        sample = row["sample"]
        contig = row["contig"]
        data = genomad_data.get(sample, {})

        scores = data.get("scores", {}).get(contig, {})
        row["genomad_plasmid_score"] = round(scores["plasmid_score"], 4) if scores else ""
        row["genomad_virus_score"] = round(scores["virus_score"], 4) if scores else ""

        plasmid_info = data.get("plasmid", {}).get(contig)
        virus_info = data.get("virus", {}).get(contig)
        if plasmid_info and virus_info:
            row["genomad_class"] = "plasmid" if scores.get("plasmid_score", 0) >= scores.get("virus_score", 0) else "virus"
            best = plasmid_info if row["genomad_class"] == "plasmid" else virus_info
            row["genomad_n_hallmarks"] = best["n_hallmarks"]
            row["genomad_topology"] = best["topology"]
        elif plasmid_info:
            row["genomad_class"] = "plasmid"
            row["genomad_n_hallmarks"] = plasmid_info["n_hallmarks"]
            row["genomad_topology"] = plasmid_info["topology"]
        elif virus_info:
            row["genomad_class"] = "virus"
            row["genomad_n_hallmarks"] = virus_info["n_hallmarks"]
            row["genomad_topology"] = virus_info["topology"]
        else:
            row["genomad_class"] = ""
            row["genomad_n_hallmarks"] = ""
            row["genomad_topology"] = ""

    with open(ANNOTATIONS_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"geNomad scores integrated into: {ANNOTATIONS_TSV}")
    print(f"Samples with geNomad data: {n_found}/{len(SAMPLES)}")


if __name__ == "__main__":
    main()
