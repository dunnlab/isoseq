# Dunn Lab isoseq workflow

Adapted from the excellent workflow developed by Natasha Picciani at https://github.com/npicciani/podocoryna/tree/main/ref_optimization .

## Installation

To create a conda environment to run this workflow:

    conda create -y -n isoseq -c bioconda -c conda-forge -c anaconda scipy snakemake bioconda busco=5.5.0 transdecoder orthofinder pandas biopython seaborn matplotlib ete3 graphviz eggnog-mapper

## Preparing for a run

You need two inputs for this workflow:

- An isoseq dataset. It should be the `.hq.fasta` output of the clustering workflow described at https://isoseq.how/clustering/cli-workflow.html.
- A set of complete proteomes from closely related species.

Clone this repo. Place the `.hq.fasta` isoseq file in `resources/isoseq`. Update the reference path and reference species in the `config.yaml`.

Place the proteomes in `resources/sequences`. Each file should be named `Genus_species.pep.fasta`, the species names that will appear in the trees will be parsed from these names.

## Running the workflow

From the `isoseq_optimization` directory, do a dry run with:

    snakemake -np

Do an actual run with the batch script, or with:

    snakemake --cores=4
    
To view the DAG:

   snakemake --dag | dot -Tpng > dag.png