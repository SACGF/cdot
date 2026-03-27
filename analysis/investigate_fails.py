#!/bin/env python3

import math
import re
import sys
import pandas as pd
from pysam.libcfaidx import FastaFile
import pyhgvs
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory


def write_clinvar_facts(df: pd.DataFrame):
    """Write ClinVar resolution stats to output/facts/clinvar.csv."""
    import os
    n_variants = len(df)
    if n_variants == 0:
        return
    n_resolved_cdot = df["can_resolve"].sum()
    cdot_resolution_pct = round(100 * n_resolved_cdot / n_variants, 1)

    # UTA resolution requires a separate UTA run — placeholder 0 if not available
    uta_resolution_pct = 0.0
    if "can_resolve_uta" in df.columns:
        n_resolved_uta = df["can_resolve_uta"].sum()
        uta_resolution_pct = round(100 * n_resolved_uta / n_variants, 1)

    facts = pd.DataFrame({
        "n_variants":          [n_variants],
        "cdot_resolution_pct": [cdot_resolution_pct],
        "uta_resolution_pct":  [uta_resolution_pct],
    })
    out = "output/facts/clinvar.csv"
    os.makedirs("output/facts", exist_ok=True)
    facts.to_csv(out, index=False)
    print(f"Written: {out}")
    print(f"  N variants:    {n_variants:,}")
    print(f"  cdot resolved: {cdot_resolution_pct}%")
    print(f"  UTA resolved:  {uta_resolution_pct}%")


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
    write_clinvar_facts(df)


if __name__ == "__main__":
    main() 
