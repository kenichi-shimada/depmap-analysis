#!/bin/bash
#SBATCH -c 1							# Request one core
#SBATCH -t 2-00:00						# Runtime in D-HH:MM format
#SBATCH -p medium						# Partition to run in
#SBATCH -o clue_logs/e%j			    # File to which STDOUT will be written
#SBATCH -e clue_logs/o%j				# File to which STDERR will be written
#SBATCH --mem-per-cpu 80G
#SBATCH --array=819,990					# 290,305,322,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,480,488,490,493,494,495,498,500,502,504,505,506,507,508,510,511,512,513

cd "$DEPMAP_SRC_DIR"

Rscript assemble_rtsne_dbscan.r ${SLURM_ARRAY_TASK_ID}
