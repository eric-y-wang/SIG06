#!/bin/bash

# The name of the job
#BSUB -J deseq
# The name of the queue you wan't to submit the job to
#BSUB -q cpuqueue
# The number of processors you want
#BSUB -n 24
# Sets memory requirements for the job (per processor)
#BSUB -R rusage[mem=12]
# How long the job will take (you job will be killed if it runs over the time limit)
#BSUB -W 24:00
# Output and error log files (optional but recommended)
#BSUB -o /lila/data/rudensky/EYW/R-out.%J
#BSUB -e /lila/data/rudensky/EYW/R-err.%J

# load bulkseq conda environment
source ~/.bashrc
mamba activate R-deseq2

# set directory (with fail safe in case it fails)
cd /lila/data/rudensky/EYW/git_projects/SIG06/analysis || { echo "Failure"; exit 1; }

Rscript job_deseq.r