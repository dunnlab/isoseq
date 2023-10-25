# Dunn Lab isoseq workflow

Adapted from the excellent workflow developed by Natasha Picciani at https://github.com/npicciani/podocoryna/tree/main/ref_optimization .

## Installation

To create a conda environment to run this workflow:

    conda create -y -n isoseq -c bioconda -c conda-forge snakemake busco transdecoder orthofinder pandas biopython seaborn matplotlib ete3 graphviz eggnog-mapper

From the `ref_optimization` directory, do a dry run with:

    snakemake -n --snakefile workflow/Snakefile

To view the DAG:

   snakemake --dag | dot -Tpng > dag.png