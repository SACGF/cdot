#!/usr/bin/env python3
"""Generate the transcript-version re-placement blocklist shipped with cdot.

Most version bumps that keep the intrinsic CDS structure (and CDS alignment gaps)
also keep the genomic placement, so ``is_version_substitution_safe`` can decide them
build-independently. A tiny minority re-align the *same* CDS structure to a different
genomic locus (eg NM_018263, whose CDS jumps ~9.75 kb between versions). The
build-independent structure cannot see that. When both versions are loaded the
provider catches it by comparing genomic CDS maps directly, but when the requested
version is absent from all loaded data there is nothing to compare - so we ship a
small precomputed blocklist of these (accession, requested_version, substitute_version)
re-placements, found here by scanning a full release where every version is present.

Writes ``cdot/hgvs/_version_replacement_blocklist.py`` (a generated constant module).
Regenerate after each data release.

Usage:
    python generate_transcript_data/generate_version_replacement_blocklist.py \
        /path/to/cdot-X.Y.Z.refseq.GRCh38.json.gz [more releases...] --build GRCh38
"""
import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

from cdot.hgvs.version_safety import (
    intrinsic_cds_structure, cds_alignment_gaps, cds_genomic_map,
)

_BASE = re.compile(r"^([A-Za-z]+_?\d+)\.(\d+)$")
_OUT = Path(__file__).resolve().parents[1] / "cdot" / "hgvs" / "_version_replacement_blocklist.py"


def replacements_in_release(path, build):
    """Yield (accession, requested_version, substitute_version) re-placement triples:
    pairs with identical intrinsic CDS structure and CDS alignment gaps (so the
    structural check would call them safe) but a different genomic CDS placement."""
    with gzip.open(path, "rt") as fh:
        txs = json.load(fh)["transcripts"]
    versions = defaultdict(list)
    for acc in txs:
        m = _BASE.match(acc)
        if m:
            versions[m.group(1)].append(int(m.group(2)))

    for base, vs in versions.items():
        if len(vs) < 2:
            continue
        info = {}
        for v in vs:
            rec = txs[f"{base}.{v}"]
            struct = intrinsic_cds_structure(rec, build)
            gmap = cds_genomic_map(rec, build)
            if struct is None or gmap is None:
                continue
            info[v] = (struct, cds_alignment_gaps(rec, build),
                       (gmap["contig"], tuple(gmap["segs"])))
        for v in info:
            for w in info:
                if v == w:
                    continue
                (sv, gv, pv), (sw, gw, pw) = info[v], info[w]
                # Structurally safe (would pass the check) but placement differs.
                if sv == sw and gv == gw and pv != pw:
                    yield (base, v, w)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("releases", nargs="+", help="cdot JSON.gz release file(s) to scan")
    ap.add_argument("--build", default="GRCh38")
    args = ap.parse_args()

    triples = set()
    for path in args.releases:
        n0 = len(triples)
        for t in replacements_in_release(path, args.build):
            triples.add(t)
        print(f"{path}: +{len(triples) - n0} (total {len(triples)})")

    triples = sorted(triples)
    sources = ", ".join(Path(p).name for p in args.releases)
    lines = [
        '"""Transcript-version re-placement blocklist (GENERATED - do not edit).',
        "",
        "Regenerate with generate_transcript_data/generate_version_replacement_blocklist.py.",
        f"Source: {sources} ({args.build}).",
        "",
        "Each (accession, requested_version, substitute_version) is a substitution whose",
        "intrinsic CDS structure and alignment gaps match (so the structural safety check",
        "would call it coordinate-safe) but whose genomic CDS placement differs, so it is",
        "refused by is_version_substitution_safe even when the requested version is absent",
        "from the loaded data (the case the provider cannot recompute).",
        '"""',
        "KNOWN_REPLACEMENT_SUBSTITUTIONS = frozenset({",
    ]
    lines += [f"    ({acc!r}, {v}, {w})," for acc, v, w in triples]
    lines += ["})", ""]
    _OUT.write_text("\n".join(lines))
    print(f"wrote {len(triples)} triples to {_OUT}")


if __name__ == "__main__":
    main()
