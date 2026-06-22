#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J genomad_pcn
#BSUB -n 12
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8500MB
#BSUB -W 4:00
#BSUB -u josne@dtu.dk
#BSUB -B
#BSUB -N
#BSUB -o genomad_pcn_%J.out
#BSUB -e genomad_pcn_%J.err

GENOMAD_DB="/work3/josne/Databases/genomad_db"
THREADS=12
BASEDIR="/work3/josne/Projects/Vibrio_Galathea3/PlasdmidCopyNum_Shengda"
SAMPLES="S2052 S2753 S2754 S2052ref"

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate /work3/josne/miniconda3/envs/genomad

cd "${BASEDIR}"

echo "=========================================="
echo "geNomad — plasmid/virus detection"
echo "Job started:  $(date)"
echo "Job ID:       ${LSB_JOBID}"
echo "Host:         $(hostname)"
echo "Database:     ${GENOMAD_DB}"
echo "Threads:      ${THREADS}"
echo "=========================================="

for SAMPLE in ${SAMPLES}; do
    CONTIGS="${BASEDIR}/07_COVERAGE/${SAMPLE}_contigs-CONTIGS.fa"
    OUTDIR="${BASEDIR}/genomad_out/${SAMPLE}"

    if [ ! -f "${CONTIGS}" ]; then
        echo "ERROR: Contigs file not found: ${CONTIGS}" >&2
        exit 1
    fi

    echo ""
    echo "--- ${SAMPLE} ---"
    echo "  Input:  ${CONTIGS}"
    echo "  Output: ${OUTDIR}"

    mkdir -p "${OUTDIR}"

    genomad end-to-end \
        --threads "${THREADS}" \
        --cleanup \
        "${CONTIGS}" \
        "${OUTDIR}" \
        "${GENOMAD_DB}"

    EXIT_CODE=$?
    if [ ${EXIT_CODE} -ne 0 ]; then
        echo "ERROR: geNomad failed for ${SAMPLE} (exit ${EXIT_CODE})" >&2
        exit ${EXIT_CODE}
    fi

    PREFIX="${SAMPLE}_contigs-CONTIGS"
    PLASMID_SUMMARY="${OUTDIR}/${PREFIX}_summary/${PREFIX}_plasmid_summary.tsv"
    VIRUS_SUMMARY="${OUTDIR}/${PREFIX}_summary/${PREFIX}_virus_summary.tsv"

    N_PLASMIDS=0
    [ -f "${PLASMID_SUMMARY}" ] && N_PLASMIDS=$(( $(wc -l < "${PLASMID_SUMMARY}") - 1 ))
    N_VIRUSES=0
    [ -f "${VIRUS_SUMMARY}" ] && N_VIRUSES=$(( $(wc -l < "${VIRUS_SUMMARY}") - 1 ))

    echo "  Plasmid contigs: ${N_PLASMIDS}"
    echo "  Virus contigs:   ${N_VIRUSES}"
done

echo ""
echo "=========================================="
echo "Job finished: $(date)"
echo "Exit code:    0"
echo "=========================================="
