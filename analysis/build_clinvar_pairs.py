#!/usr/bin/env python3
"""
Build a large (g.HGVS, c.HGVS) benchmark set from ClinVar, by joining:
  * the ClinVar VCF        -> normalised genomic g.HGVS (INFO/CLNHGVS), keyed by ALLELEID
  * variant_summary.txt.gz -> transcript c.HGVS (the "Name" column),     keyed by AlleleID

variant_summary's ReferenceAllele/AlternateAllele are often "na", so we take the
ground-truth g.HGVS from the VCF's CLNHGVS (which is already normalised HGVS and
matches the biocommons c_to_g output), and the c.HGVS from variant_summary.

Used to scale analysis/benchmark_resolution.py beyond the small committed test
sets (paper throughput/coverage numbers). The output TSV is large and
data-derived, so it lives under a data dir and is NOT committed (only this
builder is).

variant_summary Name: ``NM_000059.4(BRCA2):c.9976A>T (p.Lys3326Ter)``
  -> c.HGVS ``NM_000059.4:c.9976A>T`` (gene paren + protein suffix stripped)

Usage:
    python analysis/build_clinvar_pairs.py variant_summary.txt.gz clinvar.vcf.gz out.tsv [--max N]
"""
import argparse
import gzip
import re
import sys

# 1-based columns in variant_summary.txt.gz
COL = {"AlleleID": 1, "Name": 3, "Assembly": 17}

_TX = re.compile(r"^((?:NM_|NR_|XM_|XR_)\d+\.\d+)\([^)]*\):([cn]\.[^ ]+)")
_ALLELEID = re.compile(r"(?:^|;)ALLELEID=(\d+)")
_CLNHGVS = re.compile(r"(?:^|;)CLNHGVS=([^;]+)")


def parse_c(name):
    m = _TX.match(name)
    return f"{m.group(1)}:{m.group(2)}" if m else None


def load_vcf_g(vcf_path):
    """alleleid(int) -> genomic g.HGVS (first NC_ CLNHGVS value)."""
    g = {}
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            info = line.split("\t", 8)[7] if line.count("\t") >= 7 else ""
            a = _ALLELEID.search(info)
            h = _CLNHGVS.search(info)
            if not a or not h:
                continue
            g_hgvs = h.group(1).split(",")[0]
            if g_hgvs.startswith("NC_") and ":g." in g_hgvs:
                g[a.group(1)] = g_hgvs
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant_summary")
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--assembly", default="GRCh38")
    ap.add_argument("--max", type=int, default=0)
    args = ap.parse_args()

    print("reading VCF CLNHGVS...", file=sys.stderr)
    g_by_allele = load_vcf_g(args.vcf)
    print(f"  {len(g_by_allele):,} alleles with genomic HGVS", file=sys.stderr)

    ci = {k: v - 1 for k, v in COL.items()}
    seen = set()
    kept = scanned = 0
    with gzip.open(args.variant_summary, "rt") as fh, open(args.out, "w") as out:
        fh.readline()
        for line in fh:
            scanned += 1
            f = line.rstrip("\n").split("\t")
            if len(f) <= ci["Assembly"] or f[ci["Assembly"]] != args.assembly:
                continue
            g_hgvs = g_by_allele.get(f[ci["AlleleID"]])
            if not g_hgvs:
                continue
            c_hgvs = parse_c(f[ci["Name"]])
            if not c_hgvs:
                continue
            key = (g_hgvs, c_hgvs)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{g_hgvs}\t{c_hgvs}\n")
            kept += 1
            if args.max and kept >= args.max:
                break

    print(f"scanned {scanned:,} variant_summary rows -> wrote {kept:,} unique "
          f"(g,c) pairs to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
