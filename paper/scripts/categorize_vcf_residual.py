#!/usr/bin/env python3
"""Break down the ``incorrect`` rows of a VCF-coordinate ClinVar pass by cause.

Input is the VCF-format pairs file (``build_clinvar_pairs.py``: chrom, pos, ref, alt,
g_hgvs, c_hgvs) and the per-variant table ``resolve_clinvar_pass.py`` wrote for it, whose
``converted_g`` column holds cdot's normalised ``chrom-pos-ref-alt``. Each incorrect row
gets the first matching category:

    identity_or_symbolic   ClinVar identity (``=``) or symbolic allele, no comparable ALT
    ambiguity_code_allele  same position, ClinVar or cdot ALT holds an IUPAC ambiguity code.
                           ClinVar's VCF writes these as N; on a minus-strand transcript
                           cdot's ALT is also left uncomplemented, because bioutils
                           ``reverse_complement`` passes ambiguity codes through unchanged
    babelfish_megaallele   cdot's allele is over 10 bases longer than ClinVar's (the VCF
                           normaliser shuffled an insertion through a repeat)
    snv_diff_pos           ClinVar SNV at a different position (also counts transcripts)
    indel_mismatch         any other indel representation / left-alignment difference

Writes ``clinvar_vcf_residual.csv``, plus, with ``--examples-out``, up to
``--n-examples`` rows per category for inspecting the residual by hand.

    python paper/scripts/categorize_vcf_residual.py clinvar_pairs.GRCh38.tsv \\
        vcf_pass_cdot.csv --out clinvar_vcf_residual.csv --examples-out vcf_residual_examples.csv
"""
import argparse
import csv
from collections import Counter

CATEGORIES = ["indel_mismatch", "snv_diff_pos", "babelfish_megaallele",
              "identity_or_symbolic", "ambiguity_code_allele"]
COLUMNS = ["n_incorrect", "indel_mismatch", "snv_diff_pos", "snv_diff_pos_transcripts",
           "babelfish_megaallele", "identity_or_symbolic", "ambiguity_code_allele"]
IUPAC_AMBIGUITY = set("RYSWKMBDHVN")
MEGAALLELE_EXTRA_BASES = 10


def categorize(g_hgvs, c_hgvs, gt, cdot):
    _chrom, pos, ref, alt = gt
    _c_chrom, c_pos, c_ref, c_alt = cdot
    if alt in (".", "=", "") or alt.startswith("<") or ref == alt or "=" in g_hgvs or "=" in c_hgvs:
        return "identity_or_symbolic"
    if c_pos == pos and (set(alt.upper()) | set(c_alt.upper())) & IUPAC_AMBIGUITY:
        return "ambiguity_code_allele"
    if max(len(c_ref), len(c_alt)) > max(len(ref), len(alt)) + MEGAALLELE_EXTRA_BASES:
        return "babelfish_megaallele"
    if len(ref) == 1 and len(alt) == 1:
        return "snv_diff_pos"
    return "indel_mismatch"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs_file", help="VCF-format pairs TSV")
    ap.add_argument("pass_csv", help="resolve_clinvar_pass.py per-variant table for it")
    ap.add_argument("--out", required=True)
    ap.add_argument("--examples-out", help="CSV of example rows per category")
    ap.add_argument("--n-examples", type=int, default=20)
    args = ap.parse_args()
    csv.field_size_limit(10**9)

    incorrect = {}
    with open(args.pass_csv) as f:
        for r in csv.DictReader(f):
            if r["bucket"] == "incorrect":
                incorrect[(r["g_hgvs"], r["c_hgvs"])] = r

    counts = Counter()
    snv_tx = set()
    examples = []
    with open(args.pairs_file) as f:
        for p in csv.DictReader(f, delimiter="\t"):
            r = incorrect.pop((p["g_hgvs"], p["c_hgvs"]), None)
            if r is None:
                continue
            gt = (p["chrom"], int(p["pos"]), p["ref"], p["alt"])
            c_chrom, c_pos, c_ref, c_alt = r["converted_g"].split("-", 3)
            category = categorize(p["g_hgvs"], p["c_hgvs"], gt, (c_chrom, int(c_pos), c_ref, c_alt))
            counts[category] += 1
            if counts[category] <= args.n_examples:
                examples.append({"category": category, **p, "tx": r["tx"],
                                 "cdot_vcf": r["converted_g"]})
            if category == "snv_diff_pos":
                snv_tx.add(r["tx"])
    if incorrect:
        raise SystemExit(f"{len(incorrect)} incorrect rows have no pair in {args.pairs_file}")

    counts["n_incorrect"] = sum(counts[c] for c in CATEGORIES)
    counts["snv_diff_pos_transcripts"] = len(snv_tx)
    row = {c: counts[c] for c in COLUMNS}
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    print(row)
    if args.examples_out:
        examples.sort(key=lambda e: CATEGORIES.index(e["category"]))
        with open(args.examples_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(examples[0]))
            w.writeheader()
            w.writerows(examples)


if __name__ == "__main__":
    main()
