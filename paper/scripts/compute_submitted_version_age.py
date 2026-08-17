#!/usr/bin/env python3
"""Version-age distribution of the submitted-HGVS ClinVar corpus (paper R2).

Of the SCV-submitted strings citing RefSeq transcripts (built by
build_clinvar_submitted_pairs.py), how many cite a transcript version that is no
longer the current version, and how many cite a version absent from the current
RefSeq annotation release? This is the number that makes the historical-depth
argument concrete: variant_summary's Name column is always at the current
version, but what laboratories actually submitted is not.

Reference for "current": the released cdot RefSeq GRCh38 JSON. cdot merges
annotation releases newest-wins and records each transcript's winning source
``url``, so a transcript entry whose URL is the newest RefSeq annotation release
(auto-detected as the highest RS_YYYY_MM tag among ``annotation_releases`` URLs,
e.g. RS_2025_08) is in the current annotation; every other entry exists only
because cdot carries historical releases (or UTA-derived alignments). Because a
RefSeq annotation release contains exactly one version per accession:

  * version_current      cited version == the current release's version
  * version_not_current  accession is in the current release at a different
                         (newer) version: the cited version is retired
  * base_retired         no version of the accession is in the current release
                         at all (transcript dropped or suppressed)
  * unversioned          the submitted string cites no version

version_not_current + base_retired together are the strings whose cited version
is absent from the current annotation, i.e. resolvable only via historical
transcript data. Of those, the script also reports how many the cdot GRCh38 file
actually carries, and how many of the remainder appear in the merged all-builds
file (GRCh37-era depth).

The cdot JSONs are scanned with a streaming regex (transcript key -> first
in-entry URL), not json.load, so peak memory stays low.

Counts are over unique (AlleleID, string) pairs; an SCV-weighted not-current
percentage (using the corpus scv_count column) is also emitted.

Usage:
    python paper/scripts/compute_submitted_version_age.py \
        clinvar_submitted_pairs.GRCh38.tsv \
        --refseq-grch38 cdot-0.2.34.refseq.GRCh38.json.gz \
        --refseq-allbuilds cdot-0.2.34.all-builds-refseq-....json.gz \
        --out output/facts/clinvar_submitted_version_age.csv
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

_KEY_URL = re.compile(
    r'"((?:NM|NR|XM|XR)_\d+\.\d+)": \{"biotype"|"url": "([^"]+)"'
)
_RS_TAG = re.compile(r"RS_\d{4}_\d{2}")
_CORPUS_TX = re.compile(r"^((?:NM|NR|XM|XR)_\d+)(?:\.(\d+))?:")


def scan_transcript_urls(json_gz, keys_only=False):
    """Stream a cdot JSON.gz, returning {accession.version: first url in entry}
    (url None with keys_only). Uses a chunked regex scan with overlap so the
    decompressed file is never held in memory."""
    tx_url = {}
    pending = None
    tail = ""
    with gzip.open(json_gz, "rt") as fh:
        while True:
            chunk = fh.read(16 << 20)
            if not chunk:
                break
            text = tail + chunk
            # Only trust matches that start before the retained tail region;
            # the last 4KB is rescanned with the next chunk.
            safe_end = len(text) - 4096 if len(chunk) == (16 << 20) else len(text)
            for m in _KEY_URL.finditer(text):
                if m.start() >= safe_end:
                    break
                key, url = m.group(1), m.group(2)
                if key:
                    pending = key
                    tx_url[key] = None
                elif pending and tx_url.get(pending) is None and not keys_only:
                    tx_url[pending] = url
                    pending = None
            tail = text[safe_end:]
    return tx_url


def current_release_tag(urls):
    """Highest RS_YYYY_MM tag among annotation-release URLs (the current release)."""
    tags = {t for u in urls if u and "annotation_releases" in u
            for t in _RS_TAG.findall(u)}
    if not tags:
        sys.exit("no RS_YYYY_MM annotation_releases URLs found; is this a RefSeq cdot JSON?")
    return max(tags)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs", help="submitted-pairs corpus TSV (build_clinvar_submitted_pairs.py)")
    ap.add_argument("--refseq-grch38", required=True, help="cdot RefSeq GRCh38 JSON.gz (reference for current/historical)")
    ap.add_argument("--refseq-allbuilds", help="cdot all-builds RefSeq JSON.gz (extra depth check)")
    ap.add_argument("--out", default="output/facts/clinvar_submitted_version_age.csv")
    args = ap.parse_args()

    print(f"scanning {args.refseq_grch38} ...", file=sys.stderr)
    tx_url = scan_transcript_urls(args.refseq_grch38)
    tag = current_release_tag(tx_url.values())
    print(f"  {len(tx_url):,} transcript versions; current release tag: {tag}", file=sys.stderr)

    versions_of = defaultdict(dict)  # base -> {version:int -> url}
    for acc, url in tx_url.items():
        base, ver = acc.rsplit(".", 1)
        versions_of[base][int(ver)] = url

    def is_current(base, ver):
        url = versions_of[base].get(ver)
        return bool(url and tag in url and "annotation_releases" in url)

    current_version_of = {
        base: next((v for v, u in vers.items() if u and tag in u and "annotation_releases" in u), None)
        for base, vers in versions_of.items()
    }

    n_pairs = n_refseq = n_unversioned = 0
    n_current = n_not_current = n_base_retired = 0
    n_not_current_in_cdot = n_absent_cdot = 0
    scv_refseq_versioned = scv_not_current = 0
    absent_accessions = set()
    with open(args.pairs) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_pairs += 1
            m = _CORPUS_TX.match(row["c_hgvs"])
            if not m:
                continue
            n_refseq += 1
            base, ver = m.group(1), m.group(2)
            scv = int(row.get("scv_count") or 1)
            if ver is None:
                n_unversioned += 1
                continue
            ver = int(ver)
            scv_refseq_versioned += scv
            in_cdot = ver in versions_of.get(base, {})
            if is_current(base, ver):
                n_current += 1
            else:
                n_not_current += 1
                scv_not_current += scv
                if current_version_of.get(base) is None:
                    n_base_retired += 1
                if in_cdot:
                    n_not_current_in_cdot += 1
            if not in_cdot:
                n_absent_cdot += 1
                absent_accessions.add(f"{base}.{ver}")

    n_absent_in_allbuilds = 0
    if args.refseq_allbuilds and absent_accessions:
        print(f"scanning {args.refseq_allbuilds} (keys only) ...", file=sys.stderr)
        allbuilds = scan_transcript_urls(args.refseq_allbuilds, keys_only=True)
        n_absent_in_allbuilds = sum(1 for a in absent_accessions if a in allbuilds)

    n_versioned = n_refseq - n_unversioned

    def pct(n, d):
        return round(100 * n / d, 2) if d else 0.0

    row = {
        "n_pairs": n_pairs,
        "n_refseq_pairs": n_refseq,
        "n_unversioned": n_unversioned,
        "n_versioned": n_versioned,
        "current_release_tag": tag,
        "n_version_current": n_current,
        "version_current_pct": pct(n_current, n_versioned),
        "n_version_not_current": n_not_current,
        "version_not_current_pct": pct(n_not_current, n_versioned),
        "n_base_retired": n_base_retired,
        "base_retired_pct": pct(n_base_retired, n_versioned),
        "n_not_current_in_cdot": n_not_current_in_cdot,
        "not_current_in_cdot_pct": pct(n_not_current_in_cdot, n_not_current),
        "n_absent_cdot_grch38": n_absent_cdot,
        "absent_cdot_grch38_pct": pct(n_absent_cdot, n_versioned),
        "n_absent_but_in_allbuilds": n_absent_in_allbuilds,
        "scv_weighted_not_current_pct": pct(scv_not_current, scv_refseq_versioned),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    for k, v in row.items():
        print(f"  {k:<32} {v}")
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
