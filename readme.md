# Dunn Lab isoseq workflow

Adapted from the excellent workflow developed by Natasha Picciani at https://github.com/npicciani/podocoryna/tree/main/ref_optimization .

## Installation

To create a conda environment to run this workflow:

    conda create -y -n isoseq -c bioconda -c conda-forge -c anaconda scipy snakemake bioconda busco=5.5.0 transdecoder orthofinder pandas biopython seaborn matplotlib ete3 graphviz eggnog-mapper

## Preparing for a run

You need two inputs for this workflow:

- An isoseq dataset. It should be the `.hq.fasta` output of the clustering workflow described at https://isoseq.how/clustering/cli-workflow.html.
- A directory with set of complete proteomes from closely related species. Each file should be named `Genus_species.pep.fasta`, the species names that will appear in the trees will be parsed from these names.

## Running the workflow

You can use this repo to run an analysis without modifing any of the files here, just specify the configurations at the command line
to override the configurations in `config.yaml``.

Do a run with the following, substitute each [] with the values specific to your analysis:

    conda activate isoseq
    snakemake --cores=[cores] --config species=[Genus_species] transcriptome=[path to transcriptome] proteomes=[path to proteome folder]

To do a dry run, replace `--cores=[cores]` with `-np`.

Or copy `batch_snake.sh` and modify it for your needs.