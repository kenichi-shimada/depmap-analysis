#!/bin/bash
#SBATCH -c 1							# Request one core
#SBATCH -t 0-12:00						# Runtime in D-HH:MM format
#SBATCH -p short						# Partition to run in
#SBATCH -o clue_logs/e%j			    # File to which STDOUT will be written
#SBATCH -e clue_logs/o%j				# File to which STDERR will be written
#SBATCH --mem-per-cpu 80G
#SBATCH --array=30-200		            

cd $DEPMAP_SRC_DIR

Rscript assemble_rtsne_dbscan.r ${SLURM_ARRAY_TASK_ID}
 