# -*- coding: utf-8 -*-
# Written by Daniel Stadtmauer, modified by Natasha Picciani
# Feb 08, 2022

import os
from Bio import SeqIO
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

treeinform_dir = sys.argv[1]  # path to treeinform folder
species = sys.argv[2]  # species used for reference
outdir = sys.argv[3]  # output directory for plots


def get_busco_duplication_score(buscofile):
    with open(buscofile) as f:
        for line in f:
            if ",D:" in line:
                score = line.strip().split(",D:")[1].split("%]")[0]
                result = float(score)
                return result


def get_busco_completeness_score(buscofile):
    with open(buscofile) as f:
        for line in f:
            if "C:" in line:
                score = line.strip().split("C:")[1].split("%[")[0]
                result = float(score)
                return result


def get_number_of_seqs(collapsed_transcriptome):
    count = 0
    for record in SeqIO.parse(collapsed_transcriptome, "fasta"):
        count = count + 1
    return int(count)


# Modified from https://stackoverflow.com/questions/12432663/what-is-a-clean-way-to-convert-a-string-percent-to-a-float
def p2f(x):
    if "%" in x:
        return float(x.strip("%")) / 100
    else:
        return x


df = {}
busco_trans = {}
for ithresh in os.listdir(treeinform_dir):
    nthresh = ithresh.split("_")[1]
    crdir = os.path.join(treeinform_dir, ithresh, "cellranger")
    buscofile = os.path.join(
        treeinform_dir,
        ithresh,
        f"busco/{species}.collapsed.fasta.transcripts.fasta/short_summary.specific.metazoa_odb10.{species}.collapsed.fasta.transcripts.fasta.txt",
    )
    collapsed_transcriptome = os.path.join(
        treeinform_dir, ithresh, f"{species}.collapsed.fasta.transcripts.fasta"
    )
    sample_list = [
        sample
        for sample in os.listdir(crdir)
        if os.path.isdir(os.path.join(crdir, sample)) and sample != "reference"
    ]
    for isample in sample_list:
        busco_metrics = pd.DataFrame()
        cr_metrics = pd.read_csv(
            os.path.join(crdir, isample, "outs", "metrics_summary.csv"),
            thousands=",",
            converters={
                "Valid Barcodes": p2f,
                "Sequencing Saturation": p2f,
                "Reads Mapped Confidently to Genome": p2f,
                "Q30 Bases in Barcode": p2f,
                "Q30 Bases in RNA Read": p2f,
                "Q30 Bases in UMI": p2f,
                "Reads Mapped to Genome": p2f,
                "Reads Mapped Confidently to Intronic Regions": p2f,
                "Reads Mapped Confidently to Exonic Regions": p2f,
                "Reads Mapped Confidently to Intergenic Regions": p2f,
                "Reads Mapped Antisense to Gene": p2f,
                "Fraction Reads in Cells": p2f,
                "Reads Mapped Confidently to Transcriptome": p2f,
            },
        )
        cr_metrics.index = [isample + "_t" + str(nthresh)]
        cr_metrics["Threshold"] = float(nthresh)
        cr_metrics["Sample"] = str(isample)
        df[cr_metrics.index[0]] = cr_metrics.iloc[0]
        busco_metrics.index = [ithresh]
        busco_metrics["Threshold"] = float(nthresh)
        busco_metrics["Busco Duplication Score"] = get_busco_duplication_score(
            buscofile
        )
        busco_metrics["Busco Completeness Score"] = get_busco_completeness_score(
            buscofile
        )
        busco_metrics["Number of transcript sequences"] = get_number_of_seqs(
            collapsed_transcriptome
        )
        busco_trans[busco_metrics.index[0]] = busco_metrics.iloc[0]

res = pd.DataFrame(df).T
busco = pd.DataFrame(busco_trans).T

### Busco Duplication Score
plt.rcParams["figure.dpi"] = 300
ax = sns.stripplot(data=busco, x="Threshold", y="Busco Duplication Score")
xticks = busco["Threshold"].sort_values()
xticklabels = (f"{n}" for n in xticks)
ax.set_xticklabels(xticklabels)
ax.set_title("Busco Duplication Score (%)")
plt.savefig(f"{outdir}/busco_duplication_scores.png", format="png", facecolor="w")
plt.clf()

### Busco Completeness Score
plt.rcParams["figure.dpi"] = 300
ax = sns.stripplot(data=busco, x="Threshold", y="Busco Completeness Score")
xticks = busco["Threshold"].sort_values()
xticklabels = (f"{n}" for n in xticks)
ax.set_xticklabels(xticklabels)
ax.set_title("Busco Completeness Score (%)")
plt.savefig(f"{outdir}/busco_completeness_scores.png", format="png", facecolor="w")
plt.clf()

### Number of transcripts
plt.rcParams["figure.dpi"] = 300
ax = sns.stripplot(data=busco, x="Threshold", y="Number of transcript sequences")
xticks = busco["Threshold"].sort_values()
xticklabels = (f"{n}" for n in xticks)
ax.set_xticklabels(xticklabels)
ax.set_title("Number of transcript sequences")
plt.savefig(f"{outdir}/number_of_transcript_sequences.png", format="png", facecolor="w")
plt.clf()

### Cell Ranger metrics
#### Reads mapped confidently to transcriptome
plt.rcParams["figure.dpi"] = 300
ax = sns.stripplot(
    data=res, x="Threshold", y="Reads Mapped Confidently to Transcriptome", hue="Sample"
)
xticks = list(set(res["Threshold"].sort_values()))[:-1]
xticklabels = (f"{n}" for n in xticks)
ax.set_xticklabels(xticklabels)
ax.set_title("Reads Mapped Confidently to Transcriptome")
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.savefig(
    f"{outdir}/reads_mapped_confidently_to_trancriptome.png",
    format="png",
    facecolor="w",
    bbox_inches="tight",
)
plt.clf()

#### Reads mapped confidently to genome
plt.rcParams["figure.dpi"] = 300
ax = sns.stripplot(
    data=res, x="Threshold", y="Reads Mapped Confidently to Genome", hue="Sample"
)
xticks = list(set(res["Threshold"].sort_values()))[:-1]
xticklabels = (f"{n}" for n in xticks)
ax.set_xticklabels(xticklabels)
ax.set_title("Reads Mapped Confidently to Genome")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.savefig(
    f"{outdir}/reads_mapped_confidently_to_genome.png",
    format="png",
    facecolor="w",
    bbox_inches="tight",
)
plt.clf()

#### Reads mapped to genome
plt.rcParams["figure.dpi"] = 300
ax = sns.stripplot(data=res, x="Threshold", y="Reads Mapped to Genome", hue="Sample")
xticks = list(set(res["Threshold"].sort_values()))[:-1]
xticklabels = (f"{n}" for n in xticks)
ax.set_xticklabels(xticklabels)
ax.set_title("Reads Mapped to Genome")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.savefig(
    f"{outdir}/reads_mapped_to_genome.png",
    format="png",
    facecolor="w",
    bbox_inches="tight",
)
plt.clf()
