#!/usr/bin/env python3
"""Measure the recall cost of adding a UTR-length equality gate to the version
safety check, on top of the (already applied) CDS alignment-gap gate.

For every (transcript base, cited version V) that appears resolved in the A1 table,
look at the alternate versions W present in the build. Under the current gap-fixed
provider, count which (V,W) are deemed coordinate-safe, and of those, which differ
in 5'UTR length (start_codon) or 3'UTR length (total_cDNA - stop_codon). Weight by
how many ClinVar variants cite each (base, V), and also report how many variants
would be left with NO safe alternate version if UTR-length equality were required
(the true recall loss for the version fallback).
"""
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cdot.hgvs.dataproviders import JSONDataProvider

RELEASE = "/data/cdot_data/refseq/cdot-0.2.33.refseq.GRCh38.json.gz"
A1 = "output/clinvar_pass/refseq_full.csv"
BUILD = "GRCh38"

p = JSONDataProvider([RELEASE])
with gzip.open(RELEASE, "rt") as fh:
    txs = json.load(fh)["transcripts"]

base_re = re.compile(r"^([A-Za-z]+_?\d+)\.(\d+)$")
versions_by_base = defaultdict(list)
for acc in txs:
    m = base_re.match(acc)
    if m:
        versions_by_base[m.group(1)].append(int(m.group(2)))
for b in versions_by_base:
    versions_by_base[b] = sorted(versions_by_base[b])


def utr_lengths(acc):
    """(five_utr, three_utr) in transcript bases for the GRCh38 placement, or None."""
    r = txs.get(acc)
    if not r:
        return None
    gb = r.get("genome_builds", {}).get(BUILD)
    if gb is None or r.get("start_codon") is None or r.get("stop_codon") is None:
        return None
    total_cdna = sum(e[1] - e[0] for e in gb["exons"])
    return (r["start_codon"], total_cdna - r["stop_codon"])


# A1 variant counts per (base, V), resolved rows with an integer version only.
df = pd.read_csv(A1, dtype=str, usecols=["tx", "version", "bucket"]).fillna("")
df = df[df["bucket"].isin(["correct", "incorrect"]) & (df["tx"] != "") & (df["version"] != "")]
counts = df.groupby(["tx", "version"]).size()

# Per (base, V): classify alternate versions under the gap-fixed provider.
safe_cache = {}


def classify_alts(base, v):
    key = (base, v)
    if key in safe_cache:
        return safe_cache[key]
    v_utr = utr_lengths(f"{base}.{v}")
    gap_safe = 0          # W deemed safe by the current (gap-fixed) check
    gap_safe_utr_eq = 0   # of those, UTR length also matches
    for w in versions_by_base.get(base, []):
        if w == v:
            continue
        safe, _ = p.is_version_substitution_safe(base, v, w, BUILD)
        if not safe:
            continue
        gap_safe += 1
        w_utr = utr_lengths(f"{base}.{w}")
        if v_utr is not None and w_utr is not None and v_utr == w_utr:
            gap_safe_utr_eq += 1
    res = (gap_safe, gap_safe_utr_eq)
    safe_cache[key] = res
    return res


tot_variants = 0
var_with_gap_safe = 0          # variants that have >=1 gap-safe alternate
var_with_utr_safe = 0         # variants that still have >=1 safe alternate under UTR gate
pair_gap_safe = 0             # variant-weighted (V,W) gap-safe count
pair_utr_eq = 0              # of those, UTR length equal

for (tx, version), n in counts.items():
    v = int(version)
    gap_safe, gap_safe_utr_eq = classify_alts(tx, v)
    tot_variants += n
    pair_gap_safe += n * gap_safe
    pair_utr_eq += n * gap_safe_utr_eq
    if gap_safe > 0:
        var_with_gap_safe += n
        if gap_safe_utr_eq > 0:
            var_with_utr_safe += n


def pct(a, b):
    return f"{100*a/b:.3f}%" if b else "n/a"


print("== UTR-length gate recall cost (variant-weighted over A1 resolved set) ==")
print(f"  total resolved variants (with version)     : {tot_variants}")
print(f"  variants with >=1 gap-safe alternate W     : {var_with_gap_safe}")
print(f"  variants still safe if UTR-equality required: {var_with_utr_safe}")
print(f"  variants that LOSE all safe alternates (recall loss): "
      f"{var_with_gap_safe - var_with_utr_safe}  "
      f"({pct(var_with_gap_safe - var_with_utr_safe, var_with_gap_safe)} of gap-safe variants)")
print(f"  (V,W) substitution options, gap-safe        : {pair_gap_safe}")
print(f"    of those with EQUAL UTR length            : {pair_utr_eq}  ({pct(pair_utr_eq, pair_gap_safe)})")
print(f"    of those that DIFFER in UTR length        : {pair_gap_safe - pair_utr_eq}  "
      f"({pct(pair_gap_safe - pair_utr_eq, pair_gap_safe)})")
