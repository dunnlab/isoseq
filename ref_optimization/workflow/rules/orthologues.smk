from datetime import date
import os

def get_orthofinder_outdir():
    """
    Generate path to orthofinder gene trees folder with current date as written by orthofinder
    """
    ortho_dir='results/orthofinder'
    if os.path.exists(ortho_dir):
        for results_folder in os.listdir(ortho_dir):
            results_path=f"results/orthofinder/{results_folder}/Gene_Trees"
            return results_path
    else:
        today = date.today()
        monthDay = today.strftime("%b%d")
        outdir= f"results/orthofinder/Results_{monthDay}/Gene_Trees"
        return outdir

rule generate_longest_ORFs:
    """
    Generate open reading frames from reference transcriptome.
    """
    input:
        transcriptome_path=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"])
    output:
        expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"])
    params:
        reference=expand("{transcriptome_stem}", transcriptome_stem=config["reference"]["filestem"])
    conda:
        "../../workflow/envs/transdecoder.yaml" #transdecoder v5.5.0
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptome_path} --output_dir results/reference/{params.reference}.transdecoder_dir"

rule gunzip:
    """
    Decompress (if gzipped) protein files downloaded from public databases.
    """
    input:
        expand("resources/sequences/ensembl/{species}.pep.fasta.gz",species=ensembl_targets.loc[:,"species"]),
        expand("resources/sequences/ensemblgenomes/{species}.pep.fasta.gz",species=ensemblgenomes_targets.loc[:,"species"]),
        expand("resources/sequences/other/{species}.pep.fasta",species=other_targets.loc[:,"species"]),
        expand("resources/sequences/other_gz/{species}.pep.fasta.gz",species=other_gz_targets.loc[:,"species"]),
        expand("resources/sequences/gdrive/{species}.pep.fasta",species=gdrive_targets.loc[:,"species"]),
        reference_peptides=expand("results/reference/{transcriptome_stem}.transdecoder_dir/longest_orfs.pep", transcriptome_stem=config["reference"]["filestem"])
    output:
        expand("resources/sequences/{species}.pep.fasta", species=targets.index)
    params:
        copyfile=expand("resources/sequences/{species}.pep.fasta", species=config["species"])
    shell:
        """
        dir="resources/sequences"
        gzfiles=`find $dir/* -name '*.gz'`
        for file in $gzfiles; do
            gunzip $file
        done
        fastafiles=`find $dir/* -name '*.fasta'`
        for file in $fastafiles; do
            mv $file $dir
        done
        subdirs=`ls -d $dir/*/`
        rm -R $subdirs

        cp {input.reference_peptides} {params.copyfile}
        cp resources/local/Nemopilema_nomurai.pep.fasta resources/sequences/Nemopilema_nomurai.pep.fasta
        """

rule orthofinder:
    """
    Infer gene trees from set of protein sequences downloaded from public databases.
    """
    input:
        expand("resources/sequences/{species}.pep.fasta", species=targets.index)
    output:
        directory(get_orthofinder_outdir())
    conda:
        "../../workflow/envs/orthofinder.yaml"
    threads: 20
    shell:
        """
        rm -rf results/orthofinder
        orthofinder -t {threads} -f resources/sequences -o results/orthofinder
        """
