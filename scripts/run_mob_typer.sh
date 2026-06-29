#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -J mob_typer_pcn
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=4GB]"
#BSUB -M 4500MB
#BSUB -W 1:00
#BSUB -u josne@dtu.dk
#BSUB -B
#BSUB -N
#BSUB -o mob_typer_pcn_%J.out
#BSUB -e mob_typer_pcn_%J.err

THREADS=4
BASEDIR="/work3/josne/Projects/Vibrio_Galathea3/PlasdmidCopyNum_Shengda"
IN_DIR="${BASEDIR}/input_plasmidsRef"
OUT_DIR="${BASEDIR}/results/mob_typer"
SAMPLES="S2052ref S2753ref S2754ref S26ref"

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate /work3/josne/miniconda3/envs/mob_suite

mkdir -p "${OUT_DIR}"

echo "=========================================="
echo "mob_typer — plasmid typing (one entry per run)"
echo "Job started:  $(date)"
echo "Job ID:       ${LSB_JOBID}"
echo "Host:         $(hostname)"
echo "Threads:      ${THREADS}"
echo "=========================================="

# Run mob_typer on a single-sequence FASTA; append to combined output
run_entry() {
    local seq_id="$1"
    local header="$2"
    local seq="$3"
    local idx="$4"
    local tmpdir="$5"
    local outfile="$6"
    local header_written="$7"

    local entry_fasta="${tmpdir}/${seq_id}.fasta"
    local entry_out="${tmpdir}/${seq_id}_mob.tsv"

    printf "%s\n%s\n" "${header}" "${seq}" > "${entry_fasta}"

    mob_typer -i "${entry_fasta}" -o "${entry_out}" \
        -a "${tmpdir}/run_${idx}" -n "${THREADS}"

    if [[ "${header_written}" -eq 0 ]]; then
        cat "${entry_out}" >> "${outfile}"
    else
        tail -n +2 "${entry_out}" >> "${outfile}"
    fi
}

for SAMPLE in ${SAMPLES}; do
    INFILE="${IN_DIR}/${SAMPLE}_plasmids.fasta"
    OUTFILE="${OUT_DIR}/${SAMPLE}_mob_typer.tsv"
    TMPDIR="${OUT_DIR}/${SAMPLE}_tmp"

    echo ""
    echo "--- ${SAMPLE} ---"
    echo "  Input:  ${INFILE}"
    echo "  Output: ${OUTFILE}"

    if [ ! -f "${INFILE}" ]; then
        echo "ERROR: Input file not found: ${INFILE}" >&2
        exit 1
    fi

    mkdir -p "${TMPDIR}"
    rm -f "${OUTFILE}"

    HEADER_WRITTEN=0
    ENTRY_IDX=0
    CURRENT_HEADER=""
    CURRENT_SEQ_ID=""
    CURRENT_SEQ=""

    while IFS= read -r line || [[ -n "${line}" ]]; do
        if [[ "${line}" == ">"* ]]; then
            if [[ -n "${CURRENT_HEADER}" ]]; then
                run_entry "${CURRENT_SEQ_ID}" "${CURRENT_HEADER}" "${CURRENT_SEQ}" \
                    "${ENTRY_IDX}" "${TMPDIR}" "${OUTFILE}" "${HEADER_WRITTEN}"
                HEADER_WRITTEN=1
                ENTRY_IDX=$(( ENTRY_IDX + 1 ))
            fi
            CURRENT_HEADER="${line}"
            CURRENT_SEQ_ID="${line:1}"
            CURRENT_SEQ_ID="${CURRENT_SEQ_ID%% *}"
            CURRENT_SEQ=""
        else
            CURRENT_SEQ="${CURRENT_SEQ}${line}"$'\n'
        fi
    done < "${INFILE}"

    # Process last entry
    if [[ -n "${CURRENT_HEADER}" ]]; then
        run_entry "${CURRENT_SEQ_ID}" "${CURRENT_HEADER}" "${CURRENT_SEQ}" \
            "${ENTRY_IDX}" "${TMPDIR}" "${OUTFILE}" "${HEADER_WRITTEN}"
    fi

    rm -rf "${TMPDIR}"

    N=$(tail -n +2 "${OUTFILE}" | wc -l)
    echo "  Records typed: ${N}"
done

echo ""
echo "=========================================="
echo "Job finished: $(date)"
echo "Exit code:    0"
echo "=========================================="
