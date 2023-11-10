# Dunn Lab isoseq workflow

Adapted from the excellent workflow developed by Natasha Picciani at https://github.com/npicciani/podocoryna/tree/main/ref_optimization . This workflow is based on the phylogenetic refinement of transcriptomes described in [this manuscript](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244202). 

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

## Results

There are now two kinds of sequence output, collapsed and strict. If there is a clade with multiple sequences from the same species, it is collapsed to the single longest sequence if the branch lengths fall below a certain threshold. strict is a subset of collapsed . It includes only sequences for genes in the gene trees. The motivation for this is that if a gene isn't in the gene tree, it can't be collapsed. 

When the workflow completes, check `output/aggregate_stats.tsv` to see how the number of transcripts and BUSCO genes varies with threshold and by collapsed vs strict. Once you have identified the threshold that gives the best results, you can find the protein and transcript files in `output/treeinform/`.