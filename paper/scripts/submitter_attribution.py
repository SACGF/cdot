#!/usr/bin/env python3
"""Which submitters cite transcript versions cdot does not hold? (paper R2 / Discussion)

cdot ingests every RefSeq/Ensembl annotation release off the FTP sites, so in
principle it should hold every *published* transcript version. Yet a small fraction
of submitted ClinVar strings cite a version absent from cdot's data
(``absent_cdot_grch38`` in compute_submitted_version_age.py). Hypothesis (Dave): those
versions never appeared in an official annotation release, because a small number of
laboratories align transcripts to the genome themselves and cite the resulting
idiosyncratic version. This script tests that by attributing each absent-version
citation to its submitting laboratory.

Method: stream the ClinVar VCV XML; for each ClinicalAssertion (SCV) take its submitter
(``ClinVarAccession@SubmitterName``, else ``ClinVarSubmissionID@submitter``) and its
first RefSeq transcript c./n. expression. A cited ``base.version`` is *absent* if that
exact version is not among the versions cdot holds for the accession (same test as
compute_submitted_version_age.py, reusing ``scan_transcript_urls``). We then accumulate,
per submitter, total versioned RefSeq citations and absent-version citations, and per
absent version how many distinct submitters cite it.

Concentration is the question: if absent citations come from a few labs, and absent
versions are each cited by a single submitter, that is the self-alignment signature; if
absent versions are widely cited, they are a genuine ingest gap.

Usage:
    python paper/scripts/submitter_attribution.py \
        --xml ClinVarVCVRelease_2026-06.xml.gz \
        --refseq-grch38 cdot-0.2.34.refseq.GRCh38.json.gz \
        --out output/facts/submitter_attribution.csv [--limit N] [--top-out top.csv]
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compute_submitted_version_age import scan_transcript_urls  # noqa: E402

_TX = re.compile(r"^((?:NM|NR|XM|XR)_\d+)\.(\d+):[cn]\.")


def iter_scv_submitter_hgvs(xml_path, limit=0):
    """Yield (submitter, base, version) for each SCV's first RefSeq transcript c./n.
    expression, streaming the VCV XML one VariationArchive at a time."""
    try:
        from lxml import etree
    except ImportError:  # pragma: no cover
        import xml.etree.ElementTree as etree
    opener = gzip.open if str(xml_path).endswith(".gz") else open
    n_va = 0
    with opener(xml_path, "rb") as fh:
        for _ev, elem in etree.iterparse(fh, events=("end",), tag="VariationArchive"):
            n_va += 1
            record = elem.find("ClassifiedRecord")
            if record is None:
                record = elem.find("InterpretedRecord")
            if record is not None:
                for ca in record.iter("ClinicalAssertion"):
                    acc = ca.find("ClinVarAccession")
                    sub = (acc.get("SubmitterName") if acc is not None else None)
                    if not sub:
                        sid = ca.find("ClinVarSubmissionID")
                        sub = sid.get("submitter") if sid is not None else None
                    sub = (sub or "").strip() or "(unknown)"
                    ca_sa = ca.find("SimpleAllele")
                    if ca_sa is None:
                        continue
                    for attr in ca_sa.iter("Attribute"):
                        if attr.get("Type") == "HGVS" and attr.text:
                            m = _TX.match(attr.text.strip())
                            if m:
                                yield sub, m.group(1), int(m.group(2))
                                break
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
            if limit and n_va >= limit:
                return


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xml", required=True, help="ClinVar VCV XML release (.xml.gz)")
    ap.add_argument("--refseq-grch38", required=True, help="cdot RefSeq GRCh38 JSON.gz")
    ap.add_argument("--limit", type=int, default=0, help="stop after N VariationArchives (smoke test)")
    ap.add_argument("--out", default="output/facts/submitter_attribution.csv")
    ap.add_argument("--top-out", help="optional: per-submitter table CSV (all submitters with absent citations)")
    args = ap.parse_args()

    print(f"scanning {args.refseq_grch38} ...", file=sys.stderr)
    tx_url = scan_transcript_urls(args.refseq_grch38, keys_only=True)
    versions_of = defaultdict(set)
    for acc in tx_url:
        base, ver = acc.rsplit(".", 1)
        versions_of[base].add(int(ver))
    print(f"  {len(tx_url):,} versions across {len(versions_of):,} accessions", file=sys.stderr)

    total = defaultdict(int)          # submitter -> versioned RefSeq citations
    absent = defaultdict(int)         # submitter -> absent-version citations
    ver_submitters = defaultdict(set)  # base.ver -> {submitters citing it}
    n_scv = 0
    for sub, base, ver in iter_scv_submitter_hgvs(args.xml, args.limit):
        n_scv += 1
        total[sub] += 1
        if ver not in versions_of.get(base, ()):
            absent[sub] += 1
            ver_submitters[f"{base}.{ver}"].add(sub)

    n_absent = sum(absent.values())
    subs_with_absent = sorted(absent.items(), key=lambda kv: -kv[1])

    def share(k):
        return round(100 * sum(c for _, c in subs_with_absent[:k]) / n_absent, 1) if n_absent else 0.0

    # self-alignment signature: absent versions cited by exactly one submitter
    n_single_versions = sum(1 for subs in ver_submitters.values() if len(subs) == 1)

    row = {
        "n_versioned_citations": n_scv,
        "n_absent_citations": n_absent,
        "absent_pct": round(100 * n_absent / n_scv, 3) if n_scv else 0.0,
        "n_submitters_total": len(total),
        "n_submitters_with_absent": len(subs_with_absent),
        "top1_absent_share_pct": share(1),
        "top3_absent_share_pct": share(3),
        "top5_absent_share_pct": share(5),
        "top10_absent_share_pct": share(10),
        "n_absent_versions": len(ver_submitters),
        "n_absent_versions_single_submitter": n_single_versions,
        "single_submitter_version_pct": round(100 * n_single_versions / len(ver_submitters), 1) if ver_submitters else 0.0,
        "top_submitter_absent_n": subs_with_absent[0][1] if subs_with_absent else 0,
        "top_submitter_absent_rate_pct": round(100 * subs_with_absent[0][1] / total[subs_with_absent[0][0]], 1) if subs_with_absent else 0.0,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    for k, v in row.items():
        print(f"  {k:<36} {v}")
    print(f"Written: {out}")

    if args.top_out:
        with open(args.top_out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["submitter", "absent_citations", "total_citations", "absent_rate_pct"])
            for sub, a in subs_with_absent:
                w.writerow([sub, a, total[sub], round(100 * a / total[sub], 1)])
        print(f"Written per-submitter table: {args.top_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
