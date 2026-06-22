#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J anvio_pcn
#BSUB -n 12
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8500MB
#BSUB -W 24:00
#BSUB -u josne@dtu.dk
#BSUB -B
#BSUB -N
#BSUB -o anvio_pcn_%J.out
#BSUB -e anvio_pcn_%J.err

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

cd /work3/josne/Projects/Vibrio_Galathea3/PlasdmidCopyNum_Shengda

echo "=========================================="
echo "Anvi'o metagenomics workflow (references_mode)"
echo "Job started:  $(date)"
echo "Job ID:       ${LSB_JOBID}"
echo "Host:         $(hostname)"
echo "=========================================="

anvi-run-workflow -w metagenomics \
    -c config.json \
    -A --unlock

anvi-run-workflow -w metagenomics \
    -c config.json \
    -A --jobs 3 --resources threads=12

EXIT_CODE=$?

echo "=========================================="
echo "Job finished: $(date)"
echo "Exit code:    ${EXIT_CODE}"
echo "=========================================="

exit ${EXIT_CODE}
