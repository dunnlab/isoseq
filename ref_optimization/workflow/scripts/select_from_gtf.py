#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Dec 2, 2021
by Natasha Picciani
Python 3.7
usage: select_from_gtf.py [transcriptlist] [gtf] [outdir]

"""

from pathlib import PurePosixPath
import sys
import pandas as pd
import time

if len(sys.argv) < 2:
    sys.stderr.write(__doc__)
    sys.stderr.write("Please provide a list of transcripts and a gtf file")

else:

    # start_time=time.perf_counter()     ## Uncomment to verify performance
    transcriptList = sys.argv[1]
    gtffile = sys.argv[2]
    outputDir = sys.argv[3]

    querylist = []
    with open(transcriptList, "r") as file:
        for line in file:
            line = line.strip("\n")
            querylist.extend(line.split(" "))

    gtfname = PurePosixPath(gtffile).name
    gtf = pd.read_table(gtffile, header=None)
    newgtfFile = f"{outputDir}/{gtfname}.selected.gtf"

    # this oneliner is 4,000 times faster than for loops below
    newgtf = gtf[gtf[0].isin(querylist)]
    newgtf.to_csv(
        newgtfFile, sep="\t", quoting=3, header=False, index=False
    )  # quoting=3 is the same as quoting=csv.QUOTE_NONE

    # for name in gtf[0]:  # for transcript name in the first column of gtf file
    #     if (
    #         name not in querylist
    #     ):  # if transcript name is not in the list of transcripts to keep, then drop it
    #         index_val = gtf[
    #             gtf[0] == name
    #         ].index.values  # get the row number where it is in
    #         gtf.drop(
    #             index=index_val, inplace=True
    #         )  # drop that row and modify the actual df
    # gtf.to_csv(
    #     newgtfFile, sep="\t", quoting=3, header=False, index=False
    # )  # quoting=3 is the same as quoting=csv.QUOTE_NONE

    ## Uncomment to verify performance
    # end_time=time.perf_counter()
    # print(f"Finished filtering gtf file in {end_time - start_time:0.4f} seconds")
