# Dunn Lab isoseq workflow

Adapted from the excellent workflow developed by Natasha Picciani at https://github.com/npicciani/podocoryna/tree/main/ref_optimization .

## Installation

To create a conda environment to run this workflow:

    conda create -n isoseq -c bioconda snakemake busco transdecoder orthofinder pandas biopython seaborn matplotlib -c conda-forge ete3