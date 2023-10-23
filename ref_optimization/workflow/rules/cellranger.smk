def get_batch(wildcards):
    """
    Return string with batch information from each sample
    """
    return config["batch"][wildcards.sample]

rule mkref_cellranger:
    """
    Make references using cellranger and the thresholded transcriptome files
    """
    input:
        transcript=expand("results/reference/treeinform/threshold_{{threshold}}/{species}.collapsed.fasta.transcripts.fasta", species=config["species"]),
        gtf=expand("results/reference/treeinform/threshold_{{threshold}}/{original_gtf_name}.selected.gtf", original_gtf_name=config["reference"]["gtfname"])
    output:
        directory("results/reference/treeinform/threshold_{threshold}/cellranger/reference")
    threads: 8
    params:
        outdir="results/reference/treeinform/threshold_{threshold}/cellranger"
    shell:
        """
        cd {params.outdir}
        cellranger mkref --genome=reference \
                        --fasta=../../../../../{input.transcript} \
                        --genes=../../../../../{input.gtf}
        """

rule count_cellranger:
    """
    Run cellranger count to align single cell reads to reference transcriptome.
    """
    input:
        input_dir="results/reference/treeinform/threshold_{threshold}/cellranger/reference"
    output:
        outfile="results/reference/treeinform/threshold_{threshold}/cellranger/{sample}/outs/metrics_summary.csv",
    threads: 8
    params:
        outdir="results/reference/treeinform/threshold_{threshold}/cellranger",
        reference_dir="./reference",
        fastqs_dir="../../../../../resources/rawdata/{sample}"
    shell:
        """
        rm -rf results/reference/treeinform/threshold_{wildcards.threshold}/cellranger/{wildcards.sample}
        cd {params.outdir}
        cellranger count --id={wildcards.sample} \
                         --transcriptome={params.reference_dir} \
                         --fastqs={params.fastqs_dir} \
                         --localcores={threads} \
                         --localmem=64
        """
