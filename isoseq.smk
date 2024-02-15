import sys
import os
from datetime import date
from Bio import SeqIO
import pandas as pd
import re

configfile: "config.yaml"

include: "rules/preliminaries.smk"

print(f"thresholds: {config['thresholds']}")
print(f"species: {config['species']}")

rule all:
    input:
        ["output/aggregate_stats.tsv"] + 
        expand("output/treeinform/{species}_{threshold}_transcripts.{kind}.gtf", 
                threshold=config['thresholds'], 
                species=config['species'],
                kind=["collapsed","strict"])

rule copy_proteomes:
    output:
        touch("resources/sequences/proteomes_copied.flag")
    params:
        src = config['proteomes'],
        dest = "resources/sequences/"
    shell:
        """
        mkdir -p {params.dest}
        cp {params.src}/*.fasta {params.dest}
        """

rule orthofinder:
    """
    Infer gene trees from set of protein sequences downloaded from public databases.
    """
    input:
        "resources/sequences/{species}.annotated.pep.fasta",
        "resources/sequences/proteomes_copied.flag"
    output:
        directory("output/orthofinder/{species}/Gene_Trees/")
    log:
        "logs/run_orthofinder_{species}.log"
    threads: workflow.cores
    shell:
        """
        mkdir -p output/orthofinder
        rm -rf output/orthofinder/{wildcards.species}
        
        orthofinder -t {threads} -f resources/sequences -o output/orthofinder/{wildcards.species} > {log} 2>&1
        
        # Copy the Gene_Trees to a location that does not depend on the date
        mkdir -p output/orthofinder/{wildcards.species}/Gene_Trees
        cp output/orthofinder/{wildcards.species}/Results_*/Gene_Trees/* output/orthofinder/{wildcards.species}/Gene_Trees/
        """

rule collapse_with_treeinform:
    input:
        gene_trees="output/orthofinder/{species}/Gene_Trees/",
        fasta="resources/sequences/{species}.annotated.pep.fasta"
    output:
        collapsed_proteins="output/treeinform/{species}_{threshold}_protein.collapsed.fasta",
        strict_proteins="output/treeinform/{species}_{threshold}_protein.strict.fasta"
    params:
        outroot="output/treeinform/{species}_{threshold}_protein"
    log:
        "logs/treeprune_{species}_{threshold}.log"
    shell:
        """
        mkdir -p output/treeinform/

        # Need to add _annotated to species since this is added to the filename in rule update_fasta_headers  
        python scripts/treeinform_collapse.py -s {input.fasta} -gt {input.gene_trees} -t {wildcards.threshold} -sp {wildcards.species}_annotated -o {params.outroot} > {log} 2>&1
        """

rule proteins_to_transcripts:
    input:
        # kind is either collapsed or strict
        proteins="output/treeinform/{species}_{threshold}_protein.{kind}.fasta",
        transcriptome = "output/{species}.filtered.fasta"
    output:
        transcripts="output/treeinform/{species}_{threshold}_transcripts.{kind}.fasta"
    run:
        # A protein header:
        # >transcript_10 type:complete gc:universal transcript_10:2803-5241(+) SPARC-related modular calcium-binding protein
        #
        # A transcript header:
        # >transcript_10
        def proteins_to_transcripts(protein_file, transcripts_file, out_file):
            transcripts = {}
            with open(transcripts_file) as input_transcript_seqs:
                for record in SeqIO.parse(input_transcript_seqs, "fasta"):
                    match = re.search(r'transcript_(\d+)', record.id)
                    if match:
                        index = match.group(1)
                        transcripts[index] = record
                    else:
                        raise ValueError(f"No match found for the transcript index in: {record.id}")

            proteins = {}
            annotations = {}
            n = 0
            with open(protein_file) as input_prot_seqs:
                for record in SeqIO.parse(input_prot_seqs, "fasta"):
                    match = re.search(r'transcript_(\d+)', record.id)
                    if match:
                        index = match.group(1)
                        proteins[index] = record
                        
                        match_annotation = re.search(r' (.+)', record.description)
                        annotations[index] = match_annotation.group(1) if match_annotation else ''

                        # Print first few for debugging:
                        if n < 0:
                            print(f"index: {index}")
                            print(f"  id: {record.id}")
                            print(f"  description: {record.description}")
                            print(f"  annotation: {annotations[index]}")
                        n = n + 1

                    else:
                        raise ValueError(f"No match found for the transcript index in the protein ID: {record.id}")

            print(f"Number of annotations for {protein_file}: {len(annotations)}")
            n = 0
            with open(out_file, 'w') as output_seqs:
                for index in proteins.keys():
                    if index not in transcripts:
                        raise ValueError(f"Protein index {index} not found in transcripts.")
                    
                    transcript_record = transcripts[index]
                    transcript_record.description = f"{transcript_record.description} {annotations[index]}"
                    #transcript_record.id = f"{transcript_record.id} {annotations[index]}"
                    SeqIO.write(transcript_record, output_seqs, "fasta")
                    # print first few records for debugging
                    if (n < 0):
                        print(f"index: {index}")
                        print(f"   annotation: {annotations[index]}")
                        print(f"   id: {transcript_record.id}")
                        print(f"   description: {transcript_record.description}")
                    n = n + 1

        proteins_to_transcripts(input.proteins, input.transcriptome, output.transcripts)


rule generate_gtf:
    input:
        transcripts="output/treeinform/{species}_{threshold}_transcripts.{kind}.fasta"
    output:
        gtf="output/treeinform/{species}_{threshold}_transcripts.{kind}.gtf"
    run:
        # Example fasta header:
        # >transcript_23 triglyceride mobilization
        # Example gtf line:
        # transcript_23	x	exon	1	9573	1	+	.	gene_id "transcript_23"; transcript_id "transcript_23"; gene_name "triglyceride mobilization";

        with open(output.gtf, "w") as gtf_file:
            with open(input.transcripts) as input_transcript_seqs:
                for record in SeqIO.parse(input_transcript_seqs, "fasta"):
                    fields = record.description.split()
                    id = fields[0]
                    name = "NA"
                    if len(fields)  > 1:
                        name = " ".join(fields[1:])
                    length = len(record.seq)
                    gtf_file.write(f'{id}\tx\texon\t1\t{length}\t1\t+\t.\tgene_id "{id}"; transcript_id "{id}"; gene_name "{name}";\n')


rule busco_scores:
    input:
        fasta="output/treeinform/{species}_{threshold}_protein.{kind}.fasta"
    output:
        busco="output/busco_threshold_{threshold}_{species}_{kind}/short_summary.specific.metazoa_odb10.busco_threshold_{threshold}_{species}_{kind}.txt"
    wildcard_constraints:
        threshold="\d+(\.\d+)?"
    threads: workflow.cores
    params:
        mode="protein",
        lineage="/gpfs/gibbs/data/db/busco/metazoa_odb10",
        filename="busco_threshold_{threshold}_{species}_{kind}"
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

rule aggregate_stats:
    input:
        transcripts=expand("output/treeinform/{species}_{threshold}_transcripts.{kind}.fasta",
                           threshold=config['thresholds'], 
                           species=config['species'],
                           kind=["collapsed", "strict"]),
        busco=expand("output/busco_threshold_{threshold}_{species}_{kind}/short_summary.specific.metazoa_odb10.busco_threshold_{threshold}_{species}_{kind}.txt",
                     threshold=config['thresholds'], 
                     species=config['species'],
                     kind=["collapsed", "strict"])
    output:
        "output/aggregate_stats.tsv"
    run:
        # Create the DataFrame with the desired columns
        df = pd.DataFrame(columns=['threshold', 'species', 'kind', 'num_seqs', 'busco_single', 'busco_duplicated', 'busco_fragmented', 'busco_missing', 'busco_total'])

        # Fill in the DataFrame for transcripts
        for transcript_file in input.transcripts:
            seq_count = sum(1 for _ in SeqIO.parse(transcript_file, "fasta"))
            # threshold, species, kind = re.findall(r"threshold_([^/]+)/([^/]+)/[^/]+_transcripts.([^/]+).fasta", transcript_file)[0]
            species, threshold, kind = re.findall(r"output/treeinform/([^/]+)_([\.\d]+)_transcripts.(\w+).fasta", transcript_file)[0]
            # Add seq_count to DataFrame
            df = df.append({'threshold': threshold, 'species': species, 'kind': kind, 'num_seqs': seq_count}, ignore_index=True)

        # Fill in the DataFrame for busco
        for busco_file in input.busco:
            # Example busco_file:
            # output/busco_threshold_20_Cyanea_sp_strict/short_summary.specific.metazoa_odb10.busco_threshold_20_Cyanea_sp_strict.txt
            match = re.search(r"odb10\.busco_threshold_([\d\.]+)_(.+)_([^/]+?)\.txt", busco_file)
            if match:
                threshold, species, kind = match.groups()
            else:
                raise ValueError(f"No match found in the file name: {busco_file}")

            complete_buscos, complete_single_buscos, complete_dup_buscos, fragmented_buscos, missing_buscos, total_buscos = [0] * 6  # Initialize all counts
            with open(busco_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    if 'Complete BUSCOs' in line:
                        complete_buscos = int(line.split()[0])
                    elif 'Complete and single-copy BUSCOs' in line:
                        complete_single_buscos = int(line.split()[0])
                    elif 'Complete and duplicated BUSCOs' in line:
                        complete_dup_buscos = int(line.split()[0])
                    elif 'Fragmented BUSCOs' in line:
                        fragmented_buscos = int(line.split()[0])
                    elif 'Missing BUSCOs' in line:
                        missing_buscos = int(line.split()[0])
                    elif 'Total BUSCO groups searched' in line:
                        total_buscos = int(line.split()[0])
            
            print(f"BUSCO file: {busco_file}")
            print(f"Extracted values - single: {complete_single_buscos}, duplicated: {complete_dup_buscos}, fragmented: {fragmented_buscos}, missing: {missing_buscos}, total: {total_buscos}")

            # Update the DataFrame with BUSCO stats
            mask = (df['threshold'] == threshold) & (df['species'] == species) & (df['kind'] == kind)
            if df[mask].empty:
                print(f"No matching rows found in DataFrame for {threshold}, {species}, {kind}")
            else:
                df.loc[mask, 'busco_single'] = complete_single_buscos
                df.loc[mask, 'busco_duplicated'] = complete_dup_buscos
                df.loc[mask, 'busco_fragmented'] = fragmented_buscos
                df.loc[mask, 'busco_missing'] = missing_buscos
                df.loc[mask, 'busco_total'] = total_buscos

        # Writing DataFrame to file
        df.to_csv(output[0], sep='\t', index=False)
