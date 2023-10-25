#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Jacob Musser
modified by Natasha Picciani
last modified on March 14, 2022

Generate annotations with eggNOG mapper, summarize emapper file,
generate gene names and go terms for each gene, and produce gtf file
Usage: makeGTF_emapper.py longORFS.fasta longORFS.pep outDir
Arguments:
longORFS.fasta -- nucleotide fasta file with open reading frames
longORFS.pep -- protein fasta file with open reading frames
outDir -- path to output directory
"""

import re
import subprocess
from pathlib import PurePosixPath
import sys
from Bio import SeqIO

# Path to programs and files
emapper = "/gpfs/ysm/project/dunn/nnp9/conda_envs/emapper/bin/emapper.py"  # path to eggNOG mapper 2.0.6
python = (
    "/gpfs/ysm/project/dunn/nnp9/conda_envs/emapper/bin/python"  # path to python 3.7
)
gonames_file = "/home/nnp9/local/datasets/go_terms_2019.txt"  # GO Terms IDs and their corresponding names from 2019 GO release


# User inputs
transcriptomeFile = sys.argv[1]  # longORFs nucleotide file
proteinFile = sys.argv[2]  # longORFs peptide file
outDir = sys.argv[3]  # output directory
filename = PurePosixPath(transcriptomeFile).name

# # Functionally annotate the ORF sequences with eggNOG-mapper
subprocess.call(
    [
        python,
        emapper,
        "-i",
        proteinFile,
        "-m",
        "diamond",
        "-o",
        filename,
        "--cpu",
        "15",
        "--output_dir",
        outDir,
    ]
)

# Generate a gene list and the gene ID/transcript ID map from the original transcriptome file
geneListFile = (
    outDir + "/" + filename + ".gene_list.txt"
)  # unique list of genes, one gene per line. SHOULD REPRESENT ALL TRANSCRIPTS (NOT JUST PROTEINS)
transcriptIDtogeneIDFile = outDir + "/" + filename + ".transcriptID_to_gene.txt"
geneIDtotranscriptIDFile = outDir + "/" + filename + ".geneID_to_transcript.txt"

searchStr = "(transcript_\d+)"
protein_searchStr = "((transcript_\d+).p\d+)"

with open(transcriptomeFile, "r") as infile:
    with open(geneListFile, "w") as outfile1:
        with open(transcriptIDtogeneIDFile, "w") as outfile2:
            with open(geneIDtotranscriptIDFile, "w") as outfile3:
                for line in infile:
                    if line.startswith(">"):
                        line = line.strip("\n")
                        p = re.search(searchStr, line)
                        geneID = "Gene." + p.group(1)
                        transcriptID = p.group(1)
                        outfile1.write(geneID + "\n")
                        outfile2.write(transcriptID + "\t" + geneID + "\n")
                        outfile3.write(geneID + "\t" + transcriptID + "\n")

# Generate the gene ID to protein ID map from the eggNOG annotation file
emapperFile = outDir + "/" + filename + ".emapper.annotations"
proteinIDtogeneIDFile = (
    outDir + "/" + filename + ".proteinID_to_gene.txt"
)  # no header, tab delimited. First column proteinID, second column geneID.

with open(emapperFile, "r") as infile:
    with open(proteinIDtogeneIDFile, "w") as outfile:
        for line in infile:
            if line[0] == "#":
                continue
            else:
                line = line.strip("\n")
                elementList = line.split("\t")
                p = re.search(protein_searchStr, elementList[0])
                geneID = "Gene." + p.group(2)
                proteinID = p.group(1)
                outfile.write(proteinID + "\t" + geneID + "\n")

emapper_proteinID_field = 0
emapper_gene_name_field = 8  # Preferred_name
emapper_gene_description_field = 7  # Description
emapper_goterms_field = 9  # GOs

# Read in proteinID-to-geneID dictionary, gene list, and tab delimited fasta/sequence file to calculate length

# list of all genes
gene_list = []
for line in open(geneListFile, "r"):
    line = line.strip("\n")
    gene_list.append(line)

# proteinID_geneID_dict
proteinID_geneID_dict = {}
for line in open(proteinIDtogeneIDFile, "r"):
    line = line.strip("\n").split("\t")
    proteinID_geneID_dict[line[0]] = line[1]

genes_with_proteins_list = list(proteinID_geneID_dict.values())

# protein_length_dict
protein_sequence_dict = {}
protein_length_dict = {}
protein_sequence_dict = SeqIO.to_dict(SeqIO.parse(proteinFile, "fasta"))
for protein in protein_sequence_dict:
    protein_length_dict[protein] = len(protein_sequence_dict[protein])

# gene_protein_dict
geneID_proteinID_dict = {}
for protein in proteinID_geneID_dict:
    if proteinID_geneID_dict[protein] in geneID_proteinID_dict:
        geneID_proteinID_dict[proteinID_geneID_dict[protein]].append(protein)
    else:
        geneID_proteinID_dict[proteinID_geneID_dict[protein]] = [protein]

# sort proteins for each gene by length
for gene in geneID_proteinID_dict:
    proteins = geneID_proteinID_dict[gene]
    if len(proteins) > 1:
        geneID_proteinID_dict[gene] = sorted(
            proteins, key=protein_length_dict.get, reverse=True
        )

# add genes without proteins to gene_protein_dict
for gene in gene_list:
    if gene not in geneID_proteinID_dict:
        geneID_proteinID_dict[gene] = []

# transcript sequences dict
transcript_sequence_dict = SeqIO.to_dict(SeqIO.parse(transcriptomeFile, "fasta"))

# transcript length dict
transcript_length_dict = {}
for transcript in transcript_sequence_dict:
    transcript_length_dict[transcript] = len(transcript_sequence_dict[transcript])

# transcriptID to geneID dict
transcriptID_geneID_dict = {}
for line in open(transcriptIDtogeneIDFile, "r"):
    line = line.strip("\n").split("\t")
    transcriptID_geneID_dict[line[0]] = line[1]

# Build dictionary of goterm names for each goID
gonames_dict = {}

for line in open(gonames_file, "r"):
    line = line.strip().split("\t")
    ID = line[0]
    term = line[1]
    gonames_dict[ID] = term

# Parse the emapper annotations file
prot_name_dict = {}
prot_goterm_dict = {}
prot_description_dict = {}
prot_list = []

for line in open(emapperFile, "r"):
    if line[0] == "#":
        continue
    else:
        line = line.strip("\n").split("\t")
        prot = line[emapper_proteinID_field]
        prot_list.append(prot)
        name = line[emapper_gene_name_field]
        goterms = line[emapper_goterms_field].split(",")
        description = line[emapper_gene_description_field]

        prot_name_dict[prot] = name
        prot_description_dict[prot] = description

        if len(goterms) > 0:
            prot_goterm_dict[prot] = goterms

# Create final gene name and goterm dictionaries
gene_shortname_dict = {}
gene_longname_dict = {}
gene_longest_protein_dict = {}
gene_goterm_dict = {}

for gene in geneID_proteinID_dict:
    if len(geneID_proteinID_dict[gene]) > 0:
        longest_protein = geneID_proteinID_dict[gene][0]
        gene_longest_protein_dict[gene] = longest_protein

        if longest_protein in prot_name_dict:
            if len(prot_name_dict[longest_protein]) > 0:
                gene_longname_dict[gene] = (
                    gene
                    + " "
                    + prot_name_dict[longest_protein]
                    + " ("
                    + prot_description_dict[longest_protein]
                    + ")"
                )
                gene_shortname_dict[gene] = gene + " " + prot_name_dict[longest_protein]
            elif len(prot_description_dict[longest_protein]) > 0:
                gene_longname_dict[gene] = (
                    gene + " (" + prot_description_dict[longest_protein] + ")"
                )
                gene_shortname_dict[gene] = gene
            else:
                gene_longname_dict[gene] = gene
                gene_shortname_dict[gene] = gene
        else:
            gene_longname_dict[gene] = gene
            gene_shortname_dict[gene] = gene

        if longest_protein in prot_goterm_dict:
            if len(prot_goterm_dict[longest_protein][0]) > 0:
                gene_goterm_dict[gene] = prot_goterm_dict[longest_protein]
            else:
                gene_goterm_dict[gene] = []
        else:
            gene_goterm_dict[gene] = []
    else:
        gene_longest_protein_dict[gene] = ""
        gene_shortname_dict[gene] = gene
        gene_longname_dict[gene] = gene

# Write output files

# goterm file
goterm_header = filename + "_geneID\tgo_id\tterm"
goterm_outlist = [goterm_header]
goterm_outfile_name = outDir + "/" + filename + "_goterms.txt"

for gene in gene_goterm_dict:
    goterms = gene_goterm_dict[gene]
    if len(goterms) > 0:
        for goterm in goterms:
            try:
                goterm_outlist.append(
                    gene + "\t" + goterm + "\t" + gonames_dict[goterm]
                )
            except KeyError:
                pass

goterm_outfile = open(goterm_outfile_name, "w")
for i in goterm_outlist:
    goterm_outfile.write(i + "\n")
goterm_outfile.close()


# gtf file
gtf_outfile_name = outDir + "/" + filename + ".eggnog.gtf"
gtf_outlist = []
for transcript in transcript_length_dict:
    gene = transcriptID_geneID_dict[transcript]
    gtf_outlist.append(
        transcript
        + "\ttranscriptome\texon\t1\t"
        + str(transcript_length_dict[transcript])
        + '\t.\t+\t.\tgene_id "'
        + gene
        + '"; transcript_id "'
        + transcript
        + '"; gene_name "'
        + gene_shortname_dict[gene]
        + '";'
    )

gtf_outfile = open(gtf_outfile_name, "w")
for i in gtf_outlist:
    gtf_outfile.write(i + "\n")
gtf_outfile.close()


# names files
names_outfile_name = outDir + "/" + filename + "_gene_final_names.txt"
names_header = "geneID\tshortname\tlongname\tlongest_protein\tprotein_list"
names_outlist = [names_header]
for gene in gene_longest_protein_dict:
    shortname = gene_shortname_dict[gene]
    longname = gene_longname_dict[gene]
    longest_protein = gene_longest_protein_dict[gene]
    protein_list = geneID_proteinID_dict[gene]

    names_outlist.append(
        gene
        + "\t"
        + shortname
        + "\t"
        + longname
        + "\t"
        + longest_protein
        + "\t"
        + ",".join(protein_list)
    )

names_outfile = open(names_outfile_name, "w")
for i in names_outlist:
    names_outfile.write(i + "\n")
names_outfile.close()
