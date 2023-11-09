import sys
import os
from datetime import date
from Bio import SeqIO
import pandas as pd

configfile: "config.yaml"

print(f"thresholds: {config['thresholds']}")
print(f"species: {config['species']}")

rule all:
    input:
        expand("output/treeinform/{species}_{threshold}_protein.collapsed.fasta", 
                threshold=config['thresholds'], 
                species=config['species']) +
        expand("output/busco_threshold_{threshold}_{species}_{kind}/short_summary.specific.metazoa_odb10.busco_threshold_{threshold}_{species}_{kind}.txt",
                threshold=config['thresholds'], 
                species=config['species'],
                kind=["collapsed","strict"]) +
        expand("output/treeinform/threshold_{threshold}/{species}/{species}_transcripts.{kind}.fasta",
                threshold=config['thresholds'], 
                species=config['species'],
                kind=["collapsed","strict"]) +
        ["output/summary/aggregate_stats.txt"]


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
    

rule sanitize_headers:
    input:
        raw_fasta = config["transcriptome"]
    output:
        sanitized_fasta="output/{species}.sanitized.fasta"
    run:
        # So that we can ingest transcripts from other sources, constrain headers to isoseq
        # format if needed
        # An isoseq transcript header:
        # >transcript/10 full_length_coverage=2;length=5637

        records = []
        clean = False
        i = 0
        with open(input.raw_fasta) as input_transcript_seqs:
            for record in SeqIO.parse(input_transcript_seqs, "fasta"):
                i = i + 1
                match = re.search(r'transcript/\d+ full_length_coverage=\d+;length=\d+', record.id)
                if match:
                    if n == 1:
                        clean = True
                    else:
                        raise ValueError(f"Mixed header naming scheme")
                if not clean:
                    seq_length = len(record.seq)
                    record.description = f"transcript/{i} full_length_coverage=0;length={seq_length}"
                records.append(record)

        with open(output.sanitized_fasta, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

rule diamond_blastx:
    input:
        seqs="output/{species}.sanitized.fasta"
    output:
        table='output/blastx_{species}.tsv'
    params:
        evalue=1e-5,
        db=config['diamond']
    threads: workflow.cores
    shell:
        """
        diamond blastx \
            --threads {threads} \
            --db {params.db} \
            --query {input.seqs} \
            --out {output.table} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs 1 \
            --evalue {params.evalue}
        """

rule filter_fasta_based_on_diamond_hits:
    input:
        seqs="output/{species}.sanitized.fasta",
        hits="output/blastx_{species}.tsv"
    output:
        filtered_fasta="output/{species}.filtered.fasta"
    log:
        "logs/filter_{species}.log"
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq

        # Read the BLASTx results and store hits with reverse orientation
        reverse_hits = set()
        forward_hits = set()
        with open(input.hits, 'r') as blastx_f:
            for line in blastx_f:
                parts = line.strip().split('\t')
                if len(parts) > 11:
                    seq_id, evalue, sstart, send = parts[0], float(parts[10]), int(parts[8]), int(parts[9])
                    if evalue <= 1e-5:
                        if send < sstart:  # This implies the hit is in reverse orientation
                            reverse_hits.add(seq_id)
                        else:  # Forward orientation
                            forward_hits.add(seq_id)

        # Filter the FASTA file
        filtered_seqs = []
        total_records = 0
        records_with_hits = 0
        records_reverse_complemented = 0

        for record in SeqIO.parse(input.seqs, 'fasta'):
            total_records += 1
            if record.id in reverse_hits:
                # If in reverse orientation, add the reverse complement
                record.seq = record.seq.reverse_complement()
                filtered_seqs.append(record)
                records_with_hits += 1
                records_reverse_complemented += 1
            elif record.id in forward_hits:
                # If a forward hit, add the original sequence
                filtered_seqs.append(record)
                records_with_hits += 1

        # Write the filtered sequences to the output file
        SeqIO.write(filtered_seqs, output.filtered_fasta, 'fasta')

        # Log the counts
        with open(log[0], "w") as log_file:
            log_file.write(f"Total records read: {total_records}\n")
            log_file.write(f"Records with hits: {records_with_hits}\n")
            log_file.write(f"Records that were reverse complemented: {records_reverse_complemented}\n")


rule generate_longest_ORFs:
    message: "Generate open reading frames from reference transcriptome."
    input:
        transcriptome_path = "output/{species}.filtered.fasta"
    output:
        pep_fasta="output/{species}.pep.fasta"
    params:
        outdir = "output/{species}.transdecoder"
    log:
        "logs/transdecoder_{species}.log"
    shell:
        """
        # Use -S flag to be strand specific (top strand) since all transcripts already oriented
        TransDecoder.LongOrfs -S -t {input.transcriptome_path} --output_dir {params.outdir} > {log} 2>&1

		transcriptome_base=$(basename {input.transcriptome_path})
        cp {params.outdir}/${{transcriptome_base}}.transdecoder_dir/longest_orfs.pep {output.pep_fasta}
        """

rule run_emapper:
    input:
        pep_fasta="output/{species}.pep.fasta"
    output:
        emapper_out= "output/{species}.emapper.annotations"
    params:
        out_prefix="{species}",
        out_dir="output"
    threads: workflow.cores
    log:
        "logs/run_emapper_{species}.log"
    shell:
        """
        emapper.py -i {input.pep_fasta} -m diamond -o {params.out_prefix} --cpu {threads} --output_dir {params.out_dir} > {log} 2>&1
        """

rule update_fasta_headers:
    input:
        fasta = "output/{species}.pep.fasta",
        emapper = "output/{species}.emapper.annotations"
    output:
        updated_fasta = "resources/sequences/{species}.annotated.pep.fasta"
    run:
        import pandas as pd

        # Function to parse the FASTA file
        def read_fasta(file_path):
            with open(file_path, 'r') as f:
                sequences = {}
                current_seq = ""
                current_header = ""
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_seq:
                            sequences[current_header] = current_seq
                        current_header = line
                        current_seq = ""
                    else:
                        current_seq += line
                if current_seq:
                    sequences[current_header] = current_seq
            return sequences

        # Function to parse the eggNOG mapper output
        def read_emapper(file_path):
            annotations = pd.read_csv(file_path, sep='\t', comment='#', header=None)
            annotations = annotations[[0, 7]]
            annotations.columns = ['id', 'description']
            annotations_dict = annotations.set_index('id').to_dict()['description']
            return annotations_dict

        # Read protein sequences and eggNOG annotations
        sequences = read_fasta(input.fasta)
        annotations = read_emapper(input.emapper)

        # Update FASTA headers
        with open(output.updated_fasta, 'w') as f:
            for header, sequence in sequences.items():
                id = header.split()[0][1:]
                if id in annotations:
                    f.write(header + " " + annotations[id] + "\n")
                else:
                    f.write(header + "\n")
                f.write(sequence + "\n")

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
        outdir="output/treeinform/threshold_{threshold}/{species}"
    log:
        "logs/treeprune_{species}_{threshold}.log"
    shell:
        """
        mkdir -p {params.outdir}

        # Need to add _annotated to species since this is added to the filename in rule update_fasta_headers  
        python scripts/treeinform_collapse.py -s {input.fasta} -gt {input.gene_trees} -t {wildcards.threshold} -sp {wildcards.species}_annotated -o {params.outdir} > {log} 2>&1
        """

rule proteins_to_transcripts:
    input:
        # kind is either collapsed or strict
        proteins="output/treeinform/{species}_{threshold}_protein.{kind}.fasta",
        transcriptome = "output/{species}.filtered.fasta"
    output:
        transcripts="output/treeinform/threshold_{threshold}/{species}/{species}_transcripts.{kind}.fasta"
    run:
        # A protein header:
        # >transcript/10.p2 type:complete gc:universal transcript/10:2803-5241(+) SPARC-related modular calcium-binding protein
        #
        # A transcript header:
        # >transcript/10 full_length_coverage=2;length=5637
        def proteins_to_transcripts(protein_file, transcripts_file, out_file):
            transcripts = {}
            with open(transcripts_file) as input_transcript_seqs:
                for record in SeqIO.parse(input_transcript_seqs, "fasta"):
                    match = re.search(r'transcript/(\d+)', record.id)
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
                    match = re.search(r'transcript/(\d+)', record.id)
                    if match:
                        index = match.group(1)
                        proteins[index] = record
                        
                        match_annotation = re.search(r'\([+-]\) (.+)', record.description)
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
        transcripts=expand("output/treeinform/threshold_{threshold}/{species}/{species}_transcripts.{kind}.fasta",
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
            threshold, species, kind = re.findall(r"threshold_([^/]+)/([^/]+)/[^/]+_transcripts.([^/]+).fasta", transcript_file)[0]
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
