#!/bin/env python3

import math
import re
import sys
import pandas as pd
from pysam.libcfaidx import FastaFile
import pyhgvs
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory


def main():
    if len(sys.argv) != 1:
        sys.stderr.write(f"Usage {sys.argv[0]} hgvs_searches_combined.csv\n")
        sys.exit(1)

    filename = sys.argv[1]
    df = pd.read_csv(filename)

    non_resolve_mask = df["can_resolve"] is False
    hgvs_errors_df = df[non_resolve_mask]

    genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
    factory = JSONPyHGVSTranscriptFactory(["/home/dlawrence/Downloads/cdot-0.2.12.refseq.grch37_grch38.json.gz",
                                           "/home/dlawrence/Downloads/cdot-0.2.12.ensembl.grch37_grch38.json.gz"])


if __name__ == "__main__":
    main() 
