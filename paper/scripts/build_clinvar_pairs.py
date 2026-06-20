#!/usr/bin/env python3
"""
Build a large benchmark set from ClinVar, by joining:
  * the ClinVar VCF        -> ground-truth genomic coordinate, keyed by ALLELEID
  * variant_summary.txt.gz -> transcript c.HGVS (the "Name" column),     keyed by AlleleID

variant_summary's ReferenceAllele/AlternateAllele are often "na", so we take the
ground truth from the VCF and the c.HGVS from variant_summary.

The ground truth is captured two ways from each VCF line: the **VCF coordinate**
(CHROM/POS/REF/ALT, the primary-assembly columns) and the normalised genomic
**g.HGVS** (INFO/CLNHGVS). The benchmark prefers the VCF coordinate: comparing
``c_to_g`` output as a VCF tuple (via biocommons Babelfish) is representation-robust
(3'-shift, equivalent del/dup spellings) and, crucially, is defined even for the
ClinVar variants whose CLNHGVS uses tandem-repeat ``ref[N]`` or identity ``=``
notation that the HGVS parser cannot read. The g.HGVS column is kept for reference.

Output TSV has a header and columns ``chrom  pos  ref  alt  g_hgvs  c_hgvs``. It is
large and data-derived, so it lives under a data dir and is NOT committed (only this
builder is).

variant_summary Name: ``NM_000059.4(BRCA2):c.9976A>T (p.Lys3326Ter)``
  -> c.HGVS ``NM_000059.4:c.9976A>T`` (gene paren + protein suffix stripped)

Usage:
    python paper/scripts/build_clinvar_pairs.py variant_summary.txt.gz clinvar.vcf.gz out.tsv [--max N]
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
    """alleleid(str) -> (chrom, pos, ref, alt, g_hgvs).

    chrom/pos/ref/alt are the VCF primary-assembly columns; g_hgvs is the first
    NC_ CLNHGVS value (kept for reference). Only alleles that carry both an
    NC_ genomic CLNHGVS and a usable ALT are kept."""
    g = {}
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            info = f[7]
            a = _ALLELEID.search(info)
            h = _CLNHGVS.search(info)
            if not a or not h:
                continue
            g_hgvs = h.group(1).split(",")[0]
            if g_hgvs.startswith("NC_") and ":g." in g_hgvs:
                g[a.group(1)] = (chrom, pos, ref, alt, g_hgvs)
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
        out.write("chrom\tpos\tref\talt\tg_hgvs\tc_hgvs\n")
        for line in fh:
            scanned += 1
            f = line.rstrip("\n").split("\t")
            if len(f) <= ci["Assembly"] or f[ci["Assembly"]] != args.assembly:
                continue
            gt = g_by_allele.get(f[ci["AlleleID"]])
            if not gt:
                continue
            chrom, pos, ref, alt, g_hgvs = gt
            c_hgvs = parse_c(f[ci["Name"]])
            if not c_hgvs:
                continue
            key = (g_hgvs, c_hgvs)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{g_hgvs}\t{c_hgvs}\n")
            kept += 1
            if args.max and kept >= args.max:
                break

    print(f"scanned {scanned:,} variant_summary rows -> wrote {kept:,} unique "
          f"(g,c) pairs to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
