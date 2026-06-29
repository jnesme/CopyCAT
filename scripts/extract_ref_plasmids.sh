#!/bin/bash
# Extract plasmid sequences from reference genomes based on FASTA headers.
# Outputs one file per sample to input_plasmidsRef/.

set -euo pipefail

BASE_DIR="$(dirname "$(dirname "$(realpath "$0")")")"
IN_DIR="${BASE_DIR}/input_genomes"
OUT_DIR="${BASE_DIR}/input_plasmidsRef"

mkdir -p "${OUT_DIR}"

extract_plasmids() {
    local infile="$1"
    local outfile="$2"

    if [[ "${infile}" == *.gz ]]; then
        zcat "${infile}"
    else
        cat "${infile}"
    fi | awk '/^>/ { keep = /plasmid/ } keep { print }' > "${outfile}"

    local n
    n=$(grep -c "^>" "${outfile}" || true)
    echo "  ${n} plasmid sequence(s) → ${outfile}"
}

echo "Extracting plasmid sequences from reference genomes..."

echo "S2052ref:"
extract_plasmids \
    "${IN_DIR}/GCF_000967465.2_ASM96746v2_genomic.fna.gz" \
    "${OUT_DIR}/S2052ref_plasmids.fasta"

echo "S2753ref:"
extract_plasmids \
    "${IN_DIR}/Photobacterium_galatheae_S2753-2.fasta" \
    "${OUT_DIR}/S2753ref_plasmids.fasta"

echo "S2754ref:"
extract_plasmids \
    "${IN_DIR}/Photobacterium_galatheae_S2754.fasta" \
    "${OUT_DIR}/S2754ref_plasmids.fasta"

echo "S26ref:"
extract_plasmids \
    "${IN_DIR}/S26_GCF_000826835.2_ASM82683v2_genomic.fna.gz" \
    "${OUT_DIR}/S26ref_plasmids.fasta"

echo "Done."
