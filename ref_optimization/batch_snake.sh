#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=ycga_long
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=16           # number of cores
#SBATCH --mem-per-cpu=10G            # shared memory, scaling with CPU request

module load miniconda
snakemake --cores $SLURM_CPUS_PER_TASK --snakefile workflow/Snakefile
