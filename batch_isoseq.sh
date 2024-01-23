#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=1-00:00:00
#SBATCH --partition=ycga
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=32           # number of cores
#SBATCH --mem-per-cpu=5G            # shared memory, scaling with CPU request

module load miniconda
conda activate isoseq
snakemake --snakefile isoseq.smk --cores $SLURM_CPUS_PER_TASK --config species=Cyanea_sp transcriptome=W7.clustered.hq.fasta
