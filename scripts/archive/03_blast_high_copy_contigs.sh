#!/bin/bash
# Extract high-copy contigs and BLAST against NCBI nt
#
# Requires: anvio-9 conda env (blastn, biopython)
# Usage: bash scripts/03_blast_high_copy_contigs.sh

set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS="${BASEDIR}/results"
COVERAGE="${BASEDIR}/07_COVERAGE"

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

python3 - "${RESULTS}" "${COVERAGE}" <<'PYEOF'
import sys, csv
from Bio import SeqIO

results_dir = sys.argv[1]
coverage_dir = sys.argv[2]

elevated = set()
with open(f"{results_dir}/copy_numbers.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["classification"] == "putative_plasmid":
            elevated.add((row["sample"], row["contig"]))

with open(f"{results_dir}/high_copy_contigs.fa", "w") as out:
    for sample in sorted(set(s for s, _ in elevated)):
        wanted = {c for s, c in elevated if s == sample}
        for rec in SeqIO.parse(f"{coverage_dir}/{sample}_contigs-CONTIGS.fa", "fasta"):
            if rec.id in wanted:
                SeqIO.write(rec, out, "fasta")
PYEOF

echo "Extracted contigs:"
grep "^>" "${RESULTS}/high_copy_contigs.fa"

echo ""
echo "Running remote BLAST against NCBI nt..."
blastn \
    -query "${RESULTS}/high_copy_contigs.fa" \
    -db nt \
    -remote \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 3 \
    -evalue 1e-10 \
    -out "${RESULTS}/high_copy_blast_hits.tsv"

echo ""
echo "Top hits:"
column -t -s $'\t' "${RESULTS}/high_copy_blast_hits.tsv"
