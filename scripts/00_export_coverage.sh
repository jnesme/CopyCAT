#!/bin/bash
# Export contig and split-level coverages from anvi'o merged profiles.
# Requires: anvio-9 conda environment active.

set -euo pipefail

SAMPLES=("S2052" "S2753" "S2754" "S2052ref" "S2753ref" "S2754ref" "S26ref")

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MERGED_DIR="${BASE_DIR}/06_MERGED"
PROFILE_DIR="${BASE_DIR}/05_ANVIO_PROFILE"
CONTIGS_DIR="${BASE_DIR}/03_CONTIGS"
OUTPUT_DIR="${BASE_DIR}/07_COVERAGE"

mkdir -p "${OUTPUT_DIR}"

for sample in "${SAMPLES[@]}"; do
    if [[ -f "${OUTPUT_DIR}/${sample}_contigs-COVs.txt" ]]; then
        echo "Skipping ${sample} — already exported"
        continue
    fi

    # Single-sample groups have no merged profile; use the single profile instead
    if [[ -f "${MERGED_DIR}/${sample}/PROFILE.db" ]]; then
        PROFILE_DB="${MERGED_DIR}/${sample}/PROFILE.db"
    else
        PROFILE_DB="${PROFILE_DIR}/${sample}/${sample}/PROFILE.db"
    fi

    echo "=== Exporting ${sample} ==="
    anvi-export-splits-and-coverages \
        -p "${PROFILE_DB}" \
        -c "${CONTIGS_DIR}/${sample}-contigs.db" \
        -o "${OUTPUT_DIR}" \
        -O "${sample}_contigs" \
        --report-contigs
    anvi-export-splits-and-coverages \
        -p "${PROFILE_DB}" \
        -c "${CONTIGS_DIR}/${sample}-contigs.db" \
        -o "${OUTPUT_DIR}" \
        -O "${sample}_splits" \
        --splits-mode
    echo "Done: ${sample}"
done

echo "All exports complete."
