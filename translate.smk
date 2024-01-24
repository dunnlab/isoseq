import sys
import os
from datetime import date
from Bio import SeqIO
import pandas as pd
import re

configfile: "config.yaml"

include: "rules/preliminaries.smk"

print(f"species: {config['species']}")

rule all:
    input:
        expand("resources/sequences/{species}.annotated.pep.fasta", species=config['species']) + 
	expand("output/busco_threshold_{species}/short_summary.specific.metazoa_odb10.busco_threshold_{species}.txt", species=config['species'])

rule busco_scores:
    input:
        fasta="resources/sequences/{species}.annotated.pep.fasta"
    output:
        busco="output/busco_threshold_{species}/short_summary.specific.metazoa_odb10.busco_threshold_{species}.txt"
    wildcard_constraints:
        threshold="\d+(\.\d+)?"
    threads: workflow.cores
    params:
        mode="protein",
        lineage="/gpfs/gibbs/data/db/busco/metazoa_odb10",
        filename="busco_threshold_{species}"
    shell:
        """
        # Create a sanitized version of the input file
        sanitized_fasta=$(mktemp)
        cat {input.fasta} | sed 's|/|_|g' > $sanitized_fasta

        # Run BUSCO using the sanitized fasta file
        busco -i $sanitized_fasta -o {params.filename} --force --out_path output/ -l {params.lineage} -m {params.mode} -c {threads}

        # Remove the temporary sanitized fasta file
        rm $sanitized_fasta
        """
