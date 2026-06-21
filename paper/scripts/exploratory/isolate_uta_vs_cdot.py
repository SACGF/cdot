#!/usr/bin/env python3
"""Resolve the cdot-mismatched ClinVar variants through UTA to isolate the cause.

For every variant cdot resolved to a different VCF coordinate than ClinVar
(output/clinvar_pass/vcf_errors.csv, bucket=incorrect), resolve the SAME c.HGVS through
UTA (different transcript backend, identical biocommons engine + Babelfish) and compare:

  uta_matches_clinvar  UTA == ClinVar VCF, cdot differs  -> cdot-specific (cdot alignment)
  uta_matches_cdot     UTA == cdot,        both differ    -> not cdot-specific (shared
                                                             engine/alignment, or ClinVar
                                                             ground-truth/representation)
  uta_differs_both     three-way disagreement
  uta_no_data          UTA lacks the transcript/version (cannot compare)
  uta_error            UTA raised on the input

Needs the local UTA + SeqRepo (see the benchmark-data-locations note).
"""
import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("UTA_DB_URL", "postgresql://uta:uta@localhost:5432/uta/uta_20241220")
os.environ.setdefault("HGVS_SEQREPO_DIR", "/data/annotation/data/seqrepo/2024-12-20")

import hgvs.parser
from hgvs.assemblymapper import AssemblyMapper
from hgvs.extras.babelfish import Babelfish
import hgvs.dataproviders.uta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_resolution import classify_vcf, normalize_vcf

ERRORS = "output/clinvar_pass/vcf_errors.csv"
DEST = "output/clinvar_pass/uta_vs_cdot.csv"
BUILD = "GRCh38"

rows = []
with open(ERRORS) as fh:
    for r in csv.DictReader(fh):
        if r["bucket"] == "incorrect" and r["gt_chrom"]:
            rows.append(r)
print(f"cdot-mismatched variants to re-resolve through UTA: {len(rows)}")

provider = hgvs.dataproviders.uta.connect()
am = AssemblyMapper(provider, assembly_name=BUILD, alt_aln_method="splign",
                    replace_reference=False)
hp = hgvs.parser.Parser()
bf = Babelfish(provider, BUILD)

outcome = Counter()
by_cat = defaultdict(Counter)
out_rows = []
for r in rows:
    gt = (r["gt_chrom"], r["gt_pos"], r["gt_ref"], r["gt_alt"])
    bucket, uta_vcf = classify_vcf(am, hp, bf, gt, r["c_hgvs"])
    if bucket == "no_data":
        o = "uta_no_data"
    elif bucket == "error":
        o = "uta_error"
    elif bucket == "correct":
        o = "uta_matches_clinvar"          # UTA == ClinVar, cdot differs -> cdot-specific
    else:
        # UTA resolved but != ClinVar. Compare to cdot's call.
        cdot_norm = "-".join(str(x) for x in normalize_vcf(*r["cdot_vcf"].split("-", 3))) \
            if r["cdot_vcf"].count("-") >= 3 else r["cdot_vcf"]
        o = "uta_matches_cdot" if uta_vcf == cdot_norm else "uta_differs_both"
    outcome[o] += 1
    by_cat[r["category"]][o] += 1
    out_rows.append({**{k: r[k] for k in ("category","c_hgvs","tx","version")},
                     "gt": "-".join(gt), "cdot_vcf": r["cdot_vcf"],
                     "uta_vcf": uta_vcf or "", "verdict": o})

with open(DEST, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["category","c_hgvs","tx","version","gt","cdot_vcf","uta_vcf","verdict"])
    w.writeheader(); w.writerows(out_rows)

print("\n== UTA verdict on cdot-mismatched variants ==")
for o, n in outcome.most_common():
    print(f"  {n:>5}  {o}")
print("\n== by category ==")
for cat in sorted(by_cat):
    parts = ", ".join(f"{o}={n}" for o, n in by_cat[cat].most_common())
    print(f"  {cat:<24} {parts}")
print(f"\nWritten: {DEST}")
