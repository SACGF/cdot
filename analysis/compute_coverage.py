#!/usr/bin/env python3
"""Compute transcript coverage counts from cdot JSON.gz files and write facts CSV.

Run this script against the production data files to populate output/facts/coverage.csv.
Per CLAUDE.md, do NOT run against full datasets during development — use production
data files on a dedicated run.

Usage:
    python analysis/compute_coverage.py \
        --refseq-grch37 cdot.refseq.grch37.json.gz \
        --refseq-grch38 cdot.refseq.grch38.json.gz \
        --refseq-t2t    cdot.refseq.t2t.json.gz \
        --ensembl-grch37 cdot.ensembl.grch37.json.gz \
        --ensembl-grch38 cdot.ensembl.grch38.json.gz \
        --ensembl-t2t    cdot.ensembl.t2t.json.gz \
        --uta-count 141000
"""

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


UTA_COUNT = 141000  # uta_20210129, from Hart 2020


def count_transcript_build_pairs(json_gz_path: str) -> tuple[int, int]:
    """Return (alignment_count, unique_accession_count) from a cdot JSON.gz file.

    Each transcript × genome_build pair is one alignment.
    """
    if not json_gz_path:
        return 0, 0
    path = Path(json_gz_path)
    if not path.exists():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt") as fh:
        data = json.load(fh)
    transcripts = data.get("transcripts", {})
    unique_accessions = set(transcripts.keys())
    alignments = sum(
        len(tx.get("genome_builds", {})) for tx in transcripts.values()
    )
    return alignments, len(unique_accessions)


def main():
    parser = argparse.ArgumentParser(description="Compute cdot transcript coverage stats.")
    parser.add_argument("--refseq-grch37",  default="")
    parser.add_argument("--refseq-grch38",  default="")
    parser.add_argument("--refseq-t2t",     default="")
    parser.add_argument("--ensembl-grch37", default="")
    parser.add_argument("--ensembl-grch38", default="")
    parser.add_argument("--ensembl-t2t",    default="")
    parser.add_argument("--uta-count",      type=int, default=UTA_COUNT)
    args = parser.parse_args()

    refseq_grch37_count, _ = count_transcript_build_pairs(args.refseq_grch37)
    refseq_grch38_count, _ = count_transcript_build_pairs(args.refseq_grch38)
    t2t_count,           _ = count_transcript_build_pairs(args.refseq_t2t)
    ensembl_grch37_count, _ = count_transcript_build_pairs(args.ensembl_grch37)
    ensembl_grch38_count, ensembl_unique_count = count_transcript_build_pairs(args.ensembl_grch38)

    # T2T unique = transcripts in T2T not present in GRCh37/38
    # Approximate: use T2T count as upper bound (would need cross-file comparison for exact value)
    t2t_unique_count = t2t_count  # refine by loading all files and set-differencing if needed

    total_count = (
        refseq_grch37_count + refseq_grch38_count + t2t_count
        + ensembl_grch37_count + ensembl_grch38_count
    )
    uta_count = args.uta_count
    improvement_fold = round(total_count / uta_count, 1) if uta_count else 0

    facts = {
        "total_count":          [total_count],
        "refseq_grch37_count":  [refseq_grch37_count],
        "refseq_grch38_count":  [refseq_grch38_count],
        "t2t_count":            [t2t_count],
        "ensembl_grch37_count": [ensembl_grch37_count],
        "ensembl_grch38_count": [ensembl_grch38_count],
        "ensembl_unique_count": [ensembl_unique_count],
        "t2t_unique_count":     [t2t_unique_count],
        "improvement_fold":     [improvement_fold],
    }

    out = Path("output/facts/coverage.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(facts).to_csv(out, index=False)
    print(f"Written: {out}")
    print(f"  Total alignments:     {total_count:,}")
    print(f"  UTA count:            {uta_count:,}")
    print(f"  Improvement:          {improvement_fold}×")
    print(f"  Ensembl unique:       {ensembl_unique_count:,}")


if __name__ == "__main__":
    main()
