#!/bin/bash
#SBATCH -c 1															# Request one core

#SBATCH -t 0-1:00														# Runtime in D-HH:MM format
#SBATCH -p short														# Partition to run in
#SBATCH -o rtsne_logs/e%j       										# File to which STDOUT will be written
#SBATCH -e rtsne_logs/o%j												# File to which STDERR will be written
#SBATCH --mem-per-cpu 4G
#SBATCH --array=1-1200               									# 1-1200 in total

cd "$DEPMAP_SRC_DIR"

Rscript rtsne_original.r ${SLURM_ARRAY_TASK_ID}
 