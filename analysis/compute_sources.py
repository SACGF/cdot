#!/usr/bin/env python3
"""
Count the annotation releases ingested per consortium/build and write
output/facts/sources.csv (the `sources.*_releases` facts used in Methods).

Source of truth is generate_transcript_data/cdot_transcripts.yaml — the committed
list of every RefSeq GFF3 / Ensembl GTF release cdot merges, per genome build.
This is what Supplementary Tables S1/S2 enumerate, so the counts match the paper.

Usage:
    python analysis/compute_sources.py
"""
import csv
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
YAML = REPO / "generate_transcript_data" / "cdot_transcripts.yaml"
OUT = REPO / "output" / "facts" / "sources.csv"

BUILDS = {
    "refseq_grch37":  ("refseq", "GRCh37"),
    "refseq_grch38":  ("refseq", "GRCh38"),
    "refseq_t2t":     ("refseq", "T2T-CHM13v2.0"),
    "ensembl_grch37": ("ensembl", "GRCh37"),
    "ensembl_grch38": ("ensembl", "GRCh38"),
    "ensembl_t2t":    ("ensembl", "T2T-CHM13v2.0"),
}


def main():
    cfg = yaml.safe_load(YAML.read_text())["config"]
    row = {}
    for key, (consortium, build) in BUILDS.items():
        entries = cfg.get(consortium, {}).get(build, {}) or {}
        # each release is a mapping (name -> {url: ...}); ignore stray scalars/comments
        row[f"{key}_releases"] = sum(1 for v in entries.values() if isinstance(v, dict))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    print(f"Written: {OUT}")
    for k, v in row.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
