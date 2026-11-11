#!/usr/bin/env python3
"""Break down the ``incorrect`` rows of a VCF-coordinate ClinVar pass by cause.

Input is the VCF-format pairs file (``build_clinvar_pairs.py``: chrom, pos, ref, alt,
g_hgvs, c_hgvs) and the per-variant table ``resolve_clinvar_pass.py`` wrote for it, whose
``converted_g`` column holds cdot's normalised ``chrom-pos-ref-alt``. Each incorrect row
gets the first matching category:

    unaligned_between_exons  transcript (with ``--json``) has transcript bases that align to
                           neither of two adjacent exons. Where those bases go is an
                           alignment choice, and ClinVar's numbering can place them
                           differently from cdot's (and UTA's)
    identity_or_symbolic   ClinVar identity (``=``) or symbolic allele, no comparable ALT
    ambiguity_code_allele  same position, ClinVar or cdot ALT holds an IUPAC ambiguity code.
                           ClinVar's VCF writes these as N; on a minus-strand transcript
                           cdot's ALT is also left uncomplemented, because bioutils
                           ``reverse_complement`` passes ambiguity codes through unchanged
    long_cdot_allele       cdot's allele is over 10 bases longer than ClinVar's. Mostly a
                           c. range whose two ends are in different exons (eg c.8331_8332del
                           with c.8331 the last base of an exon): biocommons HGVS maps it to
                           a genomic span that includes the intron
    snv_diff_pos           ClinVar SNV at a different position (also counts transcripts)
    indel_mismatch         any other indel representation / left-alignment difference

Writes ``clinvar_vcf_residual.csv``, plus, with ``--examples-out``, up to
``--n-examples`` rows per category for inspecting the residual by hand.

    python paper/scripts/categorize_vcf_residual.py clinvar_pairs.GRCh38.tsv \\
        vcf_pass_cdot.csv --json cdot.refseq.grch38.json.gz --out clinvar_vcf_residual.csv --examples-out vcf_residual_examples.csv
"""
import argparse
import csv
import gzip
import json
from collections import Counter
from facts_io import write_typed_facts  # noqa: E402

CATEGORIES = ["unaligned_between_exons", "indel_mismatch", "snv_diff_pos", "long_cdot_allele",
              "identity_or_symbolic", "ambiguity_code_allele"]
COLUMNS = ["n_incorrect", "unaligned_between_exons", "unaligned_between_exons_transcripts",
           "indel_mismatch", "snv_diff_pos", "snv_diff_pos_transcripts",
           "long_cdot_allele", "identity_or_symbolic", "ambiguity_code_allele"]
IUPAC_AMBIGUITY = set("RYSWKMBDHVN")
LONG_ALLELE_EXTRA_BASES = 10


def categorize(g_hgvs, c_hgvs, gt, cdot):
    _chrom, pos, ref, alt = gt
    _c_chrom, c_pos, c_ref, c_alt = cdot
    if alt in (".", "=", "") or alt.startswith("<") or ref == alt or "=" in g_hgvs or "=" in c_hgvs:
        return "identity_or_symbolic"
    if c_pos == pos and (set(alt.upper()) | set(c_alt.upper())) & IUPAC_AMBIGUITY:
        return "ambiguity_code_allele"
    if max(len(c_ref), len(c_alt)) > max(len(ref), len(alt)) + LONG_ALLELE_EXTRA_BASES:
        return "long_cdot_allele"
    if len(ref) == 1 and len(alt) == 1:
        return "snv_diff_pos"
    return "indel_mismatch"


def unaligned_between_exons_transcripts(json_filename, build):
    """ Accessions whose exons (cdot ``[alt_start, alt_end, exon_id, tx_start, tx_end, gap]``,
        1-based transcript coordinates) leave transcript bases between two exons """
    with gzip.open(json_filename) as f:
        transcripts = json.load(f)["transcripts"]
    accessions = set()
    for ac, tx in transcripts.items():
        if b := tx["genome_builds"].get(build):
            exons = sorted(b["exons"], key=lambda e: e[3])
            if any(e2[3] != e1[4] + 1 for e1, e2 in zip(exons, exons[1:])):
                accessions.add(ac)
    return accessions


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs_file", help="VCF-format pairs TSV")
    ap.add_argument("pass_csv", help="resolve_clinvar_pass.py per-variant table for it")
    ap.add_argument("--out", required=True)
    ap.add_argument("--examples-out", help="CSV of example rows per category")
    ap.add_argument("--n-examples", type=int, default=20)
    ap.add_argument("--json", help="cdot JSON.gz the pass used, for unaligned_between_exons")
    ap.add_argument("--build", default="GRCh38")
    args = ap.parse_args()
    csv.field_size_limit(10**9)
    gap_tx = unaligned_between_exons_transcripts(args.json, args.build) if args.json else set()

    incorrect = {}
    with open(args.pass_csv) as f:
        for r in csv.DictReader(f):
            if r["bucket"] == "incorrect":
                incorrect[(r["g_hgvs"], r["c_hgvs"])] = r

    counts = Counter()
    snv_tx = set()
    gap_cited_tx = set()
    examples = []
    with open(args.pairs_file) as f:
        for p in csv.DictReader(f, delimiter="\t"):
            r = incorrect.pop((p["g_hgvs"], p["c_hgvs"]), None)
            if r is None:
                continue
            gt = (p["chrom"], int(p["pos"]), p["ref"], p["alt"])
            c_chrom, c_pos, c_ref, c_alt = r["converted_g"].split("-", 3)
            accession = p["c_hgvs"].split(":", 1)[0]
            if accession in gap_tx:
                category = "unaligned_between_exons"
                gap_cited_tx.add(accession)
            else:
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
    counts["unaligned_between_exons_transcripts"] = len(gap_cited_tx)
    row = {c: counts[c] for c in COLUMNS}
    write_typed_facts(args.out, row)
    print(row)
    if args.examples_out:
        examples.sort(key=lambda e: CATEGORIES.index(e["category"]))
        with open(args.examples_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(examples[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(examples)


if __name__ == "__main__":
    main()
