#!/usr/bin/env python3
"""Enumerate the mode-2 'genomic re-placement' version pairs for a blocklist.

These are version pairs the structural safety check (intrinsic CDS structure + CDS
alignment gaps + UTR length, ie is_version_substitution_safe) calls coordinate-safe,
yet whose actual CDS genomic placement differs in the build - so a coding variant
moves (NM_018263 is the ClinVar-observed example). The build-independent structural
check cannot see this; a precomputed blocklist can, and works even in the cross-build
fallback (the requested version absent from the target build) because it is built
offline where every version is present.

Output: distinct transcripts + (base, V, W) pairs, and a sample.
"""
import gzip
import json
import re
from collections import defaultdict

from cdot.hgvs.dataproviders import JSONDataProvider
from cdot.hgvs.version_safety import cds_genomic_map

RELEASE = "/data/cdot_data/refseq/cdot-0.2.33.refseq.GRCh38.json.gz"
BUILD = "GRCh38"

p = JSONDataProvider([RELEASE])
with gzip.open(RELEASE, "rt") as fh:
    txs = json.load(fh)["transcripts"]

base_re = re.compile(r"^([A-Za-z]+_?\d+)\.(\d+)$")
versions = defaultdict(list)
for acc in txs:
    m = base_re.match(acc)
    if m:
        versions[m.group(1)].append(int(m.group(2)))
for b in versions:
    versions[b].sort()


def placement(acc):
    r = txs.get(acc)
    if not r:
        return None
    m = cds_genomic_map(r, BUILD)
    if m is None:
        return None
    return (m["contig"], tuple(m["segs"]))


blocklist = []          # (base, V, W) deemed safe but placement differs
blocklist_tx = set()
safe_checked = 0
multi = 0
for base, vs in versions.items():
    if len(vs) < 2:
        continue
    multi += 1
    plc = {v: placement(f"{base}.{v}") for v in vs}
    for v in vs:
        for w in vs:
            if v == w:
                continue
            safe, _ = p.is_version_substitution_safe(base, v, w, BUILD)
            if not safe:
                continue
            safe_checked += 1
            pv, pw = plc[v], plc[w]
            if pv is None or pw is None:
                continue
            if pv != pw:
                blocklist.append((base, v, w))
                blocklist_tx.add(base)

print(f"transcripts with >=2 versions      : {multi}")
print(f"structurally-safe pairs checked     : {safe_checked}")
print(f"mode-2 pairs (safe but CDS re-placed): {len(blocklist)}")
print(f"distinct transcripts affected       : {len(blocklist_tx)}")
print("sample:")
for base, v, w in blocklist[:25]:
    print(f"   {base} .{v}->.{w}")
