rule collapse_with_treeinform:
    """
    Generate a collapsed protein file based on gene trees and threshold
    of branch lengths for collapsing sequences. If sequences are connected by branches
    with length less than threshold, keep the longest sequence variant.
    """
    input:
        gene_trees=get_orthofinder_outdir(),
        script="workflow/scripts/treeinform_collapse.py"
    output:
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"])
    params:
        peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"]),
        sp=config["species"],
        outdir="results/reference/treeinform/threshold_{threshold}"
    shell:
        """
        mkdir {params.outdir}
        python {input.script} -s {params.peptides} -gt {input.gene_trees} -t {wildcards.threshold} -sp {params.sp} -o {params.outdir}
        """

rule plot_subtree_histogram:
    """
    Plot distribution of subtree branch lengths from all gene trees in a histogram.
    """
    input:
        gene_trees=get_orthofinder_outdir(),
        script="workflow/scripts/plot_subtree_histogram.py",
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"])
    output:
        image="results/reference/treeinform/threshold_{threshold}/branch.length.hist_{threshold}.png",
    params:
        outdir="results/reference/treeinform/threshold_{threshold}"
    shell:
        """
        python {input.script} -gt {input.gene_trees} -t {wildcards.threshold} -b 20000 -o {params.outdir}
        """

rule match_transcripts:
    """
    Generate transcript list and fasta file with transcriptome that corresponds to protein file.
    """
    input:
        original_reference=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"]),
        collapsed_proteins=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta", species=config["species"]),
        script="workflow/scripts/pull_transcripts.py"
    output:
        transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"]),
        list_of_transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.list.txt", species=config["species"])
    shell:
        """
        python {input.script} {input.original_reference} {input.collapsed_proteins}
        """

rule select_from_gtf:
    """
    Generate a gtf file only with sequences included in the collapsed transcriptome.
    """
    input:
        gtf=expand("{original_gtf}", original_gtf=config["reference"]["gtf"]),
        list_of_transcripts=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.list.txt", species=config["species"])
    output:
        expand("results/reference/treeinform/threshold_{{threshold}}/{original_gtf_name}.selected.gtf", original_gtf_name=config["reference"]["gtfname"])
    params:
        script="workflow/scripts/select_from_gtf.py",
        outdir="results/reference/treeinform/threshold_{threshold}"
    shell:
        """
        python {params.script} {input.list_of_transcripts} {input.gtf} {params.outdir}
        """
