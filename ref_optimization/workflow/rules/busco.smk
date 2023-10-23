rule busco_scores:
    """
    Assess completeness of collapsed reference transcriptome
    with BUSCO based on dataset metazoa_obd10.
    """
    input:
        expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"])
    output:
        directory("results/reference/treeinform/threshold_{threshold}/busco")
    threads: 20
    conda:
        "../envs/busco.yaml" #busco=5.1.3
    params:
        download_path="resources/busco_downloads",
        mode="transcriptome",
        lineage="metazoa",
        filename=expand("{species}.collapsed.fasta.transcripts.fasta", species=config["species"])
    shell:
        "busco -i {input} -o {params.filename} --force --out_path {output} -l {params.lineage} -m {params.mode} --download_path {params.download_path} -c {threads}"

rule plot_QC:
    """
    Plot BUSCO scores, cellranger metrics, and number of transcripts in collapsed transcriptome files
    """
    input:
        expand("results/reference/treeinform/threshold_{threshold}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"], threshold=THRESHOLD_VALS),
        expand("results/reference/treeinform/threshold_{threshold}/busco", threshold=THRESHOLD_VALS),
        expand("results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/outs/metrics_summary.csv", sample=SAMPLE_IDS, threshold=THRESHOLD_VALS)
    output:
        "results/reports/busco_duplication_scores.png",
        "results/reports/busco_completeness_scores.png",
        "results/reports/number_of_transcript_sequences.png",
        "results/reports/reads_mapped_confidently_to_trancriptome.png",
        "results/reports/reads_mapped_confidently_to_genome.png"
    params:
        script="workflow/scripts/plot_QC.py",
        treeinform_dir="results/reference/treeinform",
        species=config["species"],
        outdir="results/reports"
    shell:
        """
        python {params.script} {params.treeinform_dir} {params.species} {params.outdir}
        """
