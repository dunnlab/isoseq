# -*- coding: utf-8 -*-
# Written by Natasha Picciani
# Oct 31, 2022
# usage: pull_transcripts.py original_transcriptome collapsed_protein_file

from Bio import SeqIO
import re
import sys

original_transcriptome = sys.argv[1]
collapsed_protein_file = sys.argv[2]

transcript_seqs = set()
with open(collapsed_protein_file) as proteinfile:
    for record in SeqIO.parse(proteinfile, "fasta"):
        # transcript_id = re.search(r"(PB.\d+.\d+|)", record.id)
        # transcript_seqs.add(transcript_id[0])
        transcript_id = re.search(r"(transcript_\d+).p\d+", record.id)
        transcript_seqs.add(transcript_id.group(1))

collapsed = []
with open(original_transcriptome) as transcriptfile:
    for record in SeqIO.parse(transcriptfile, "fasta"):
        for transcript in transcript_seqs:
            if transcript == record.id:
                # if transcript == record.id.split("|")[0]:
                collapsed.append(record)

outfile = f"{collapsed_protein_file}.transcripts.fasta"
outfile2 = f"{collapsed_protein_file}.transcripts.list.txt"
with open(outfile, "w") as out:
    with open(outfile2, "w") as out2:
        for record in collapsed:
            SeqIO.write(record, out, "fasta")
            out2.write(record.id + "\n")
