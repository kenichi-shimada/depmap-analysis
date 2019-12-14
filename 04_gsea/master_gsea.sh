#!/bin/bash
#SBATCH -c 1                            # Request one core

#SBATCH -t 0-12:00                      # Runtime in D-HH:MM format
#SBATCH -p short                        # Partition to run in
#SBATCH -e fgsea_logs/e7_%j	        	# File to which STDERR will be written
#SBATCH -o fgsea_logs/o7_%j	        	# File to which OUTPUT will be written
#SBATCH --mem-per-cpu 16G				# 1e7 => 16G, 1e6 => 8G, 1e5 => 4G
#SBATCH --array=1-2						# 1-2

cd "$DEPMAP_SRC_DIR"
Rscript gsea.r ${SLURM_ARRAY_TASK_ID}
