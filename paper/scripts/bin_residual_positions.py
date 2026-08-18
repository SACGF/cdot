#!/usr/bin/env python3
"""Bin the ClinVar VCF-pass projections by relative CDS position (R5 positional drift).

Companion to the positional-drift deciles in compute_version_stability.py: that
script asks where along the CDS a *version bump* moves coordinates; this one asks
where along the CDS the residual `incorrect` ClinVar projections sit. If
projection errors concentrated toward the 3' end, the incorrect rows' relative
CDS positions would skew high against the correct rows' baseline distribution.

Consumes the per-variant table from resolve_clinvar_pass.py (VCF mode; the run
documented in claude/clinvar_diff_coordinates.md) plus the same cdot release
JSON.gz the pass used (only start/stop codons are kept, then the JSON is freed).
Every `correct` and `incorrect` row's cited c. position is anchored to a relative
CDS position (intronic offsets anchor to their exon boundary; ranges to their
5'-most end) and binned into deciles (1 = 5'-most tenth, 10 = 3'-most); 5'UTR
(c.-N) and 3'UTR (c.*N) citations are counted separately, as are n. transcripts.

Writes a single-row facts CSV: per-decile *shares* (percent of that bucket's
binned coding rows) for the incorrect and correct buckets, plus the UTR /
non-coding / unbinnable counts. The committed snapshot lives at
paper/empirical_results/clinvar_residual_positions.csv (no Snakefile rule, like
the other clinvar_vcf facts: the per-variant table is a dedicated multi-hour run).

Per CLAUDE.md: never run against full datasets in development; verify on a small
table from tests/test_data/clinvar_hgvs/ first.

Usage:
    python paper/scripts/bin_residual_positions.py \
        output/clinvar_pass/refseq_full_vcf.csv \
        --json cdot-0.2.33.refseq.GRCh38.json.gz \
        --out output/facts/clinvar_residual_positions.csv
"""
import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

N_BINS = 10

# First cited position after the kind letter: "c.100A>G", "c.-49_12del",
# "c.*103del", "c.4072-1234G>T" -> the leading (possibly -/*) anchor number.
_POS = re.compile(r":[cn]\.(\*?)(-?\d+)")


def load_cds_lengths(path):
    """Return {accession.version: cds_len} for coding transcripts in a release."""
    with gzip.open(path, "rt") as fh:
        transcripts = json.load(fh)["transcripts"]
    out = {}
    for key, t in transcripts.items():
        sc, ec = t.get("start_codon"), t.get("stop_codon")
        if sc is not None and ec is not None and ec > sc:
            out[key] = ec - sc
    return out


def bin_row(c_hgvs, cds_lengths):
    """Classify one row: decile 1..N_BINS, '5putr', '3putr', 'noncoding_n',
    'no_cds_len' or 'unparsed'."""
    m = _POS.search(c_hgvs)
    if m is None:
        return "unparsed"
    if ":n." in c_hgvs:
        return "noncoding_n"
    star, pos = m.group(1), int(m.group(2))
    if star:
        return "3putr"
    if pos < 0:
        return "5putr"
    key = c_hgvs.split(":", 1)[0]
    cds_len = cds_lengths.get(key)
    if not cds_len:
        return "no_cds_len"
    pos = min(max(pos, 1), cds_len)
    return min((pos - 1) * N_BINS // cds_len, N_BINS - 1) + 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pass_csv", help="per-variant table from resolve_clinvar_pass.py")
    ap.add_argument("--json", required=True, help="cdot release JSON.gz the pass used")
    ap.add_argument("--out", default="output/facts/clinvar_residual_positions.csv")
    args = ap.parse_args()

    print(f"loading CDS lengths from {args.json} ...", flush=True)
    cds_lengths = load_cds_lengths(args.json)
    print(f"  {len(cds_lengths)} coding transcript versions")

    bins = {"correct": Counter(), "incorrect": Counter()}
    with open(args.pass_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            counter = bins.get(row["bucket"])
            if counter is not None:
                counter[bin_row(row["c_hgvs"], cds_lengths)] += 1

    facts = {}
    for bucket, counter in bins.items():
        binned = sum(counter[i] for i in range(1, N_BINS + 1))
        facts[f"{bucket}_binned"] = binned
        for i in range(1, N_BINS + 1):
            facts[f"{bucket}_decile{i}_pct"] = (
                round(100 * counter[i] / binned, 1) if binned else 0.0)
        for extra in ("5putr", "3putr", "noncoding_n", "no_cds_len", "unparsed"):
            facts[f"{bucket}_{extra}"] = counter[extra]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(facts))
        w.writeheader()
        w.writerow(facts)
    print(f"Written: {out}")
    for k, v in facts.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
