rule sanitize_headers:
    input:
        raw_fasta = config["transcriptome"]
    output:
        sanitized_fasta="output/{species}.sanitized.fasta"
    run:

        records = []
        has_transcript_index = True  # Initially assume the file has transcript indexes

        # Check each record's header format
        with open(input.raw_fasta) as input_transcript_seqs:
            for record in SeqIO.parse(input_transcript_seqs, "fasta"):
                if not re.search(r'transcript/\d+', record.id):
                    has_transcript_index = False
                    break

        # Re-read the file to process records
        with open(input.raw_fasta) as input_transcript_seqs:
            for i, record in enumerate(SeqIO.parse(input_transcript_seqs, "fasta"), start=1):
                if has_transcript_index:
                    match = re.search(r'transcript/(\d+)', record.id)
                    transcript_index = match.group(1)
                    record.id = f"transcript_{transcript_index}"
                else:
                    record.id = f"transcript_{i}"
                record.description = record.id
                records.append(record)

        # Write to output
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
                    seq_id, evalue, qstart, qend, sstart, send = parts[0], float(parts[10]), int(parts[6]), int(parts[7]), int(parts[8]), int(parts[9])
                    if evalue <= 1e-5:
                        # reverse hits seem to have reversed queries, but will check both to be sure
                        if (send < sstart) ^ (qend < qstart):  # Use exclusive or to see if just one of them is flipped
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
        pep_fasta="output/{species}.pep.all.fasta"
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

rule filter_translations:
    input:
        pep_fasta="output/{species}.pep.all.fasta"
    output:
        pep_fasta="output/{species}.pep.fasta"
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq

        # Filter the FASTA file
        filtered_seqs = []

        # Raw headers have the following form:
        # >transcript_24.p1 type:complete gc:universal transcript_24:66-4712(+)
        # >transcript_24.p2 type:complete gc:universal transcript_24:3934-4371(+)
        # >transcript_24.p3 type:complete gc:universal transcript_24:4438-4824(+)
        #
        # Retain only those with .p1, and trim the header to just the transcript index

        for record in SeqIO.parse(input.pep_fasta, 'fasta'):
            if record.id.endswith(".p1"):
                record.id = record.id.split(".")[0]
                record.description = ""
                filtered_seqs.append(record)

        # Write the filtered sequences to the output file
        SeqIO.write(filtered_seqs, output.filtered_pep_fasta, 'fasta')


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
