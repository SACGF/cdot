#!/usr/bin/env python3
"""
Build the SUBMITTED-HGVS ClinVar benchmark corpus: SCV-submitted transcript HGVS
strings (as the submitting laboratory wrote them) joined to the variant's genomic
ground truth via AlleleID.

Why this exists: build_clinvar_pairs.py takes c.HGVS from the ``Name`` column of
variant_summary.txt.gz, which is ClinVar's own recomputed preferred-transcript name
and is always at the CURRENT transcript version. Benchmarks on that corpus measure
the easy case. The per-SCV HGVS attributes in the VCV XML preserve what each lab
actually submitted, including transcript versions that were current when the
variant was classified but have since been superseded, so this corpus can exercise
historical transcript depth.

Sources (two input modes):

  * ``--xml ClinVarVCVRelease_YYYY-MM.xml.gz`` - stream the full VCV XML release
    (ncbi.nlm.nih.gov/clinvar; large, tens of GB uncompressed) with iterparse.
    Per VariationArchive: AlleleID from ClassifiedRecord/SimpleAllele (or
    InterpretedRecord in pre-2025 releases); per ClinicalAssertion: the first
    transcript c./n. expression among its ``AttributeSet/Attribute[@Type="HGVS"]``
    values, kept verbatim. Fully self-contained but slow (hours).
  * ``--scv-csv-dir DIR`` - consume per-SCV CSVs previously extracted from that
    XML (columns ``hgvs`` = first HGVS attribute value, ``allele_id``). Fast path
    when an extraction already exists. The two modes differ only in that the CSV
    keeps the first HGVS value rather than the first *transcript* value; on a
    50,000-row sample the two never disagreed (no SCV had a transcript expression
    preceded by a non-transcript one).

Ground truth: the ClinVar VCF (clinvar.GRCh38.vcf.gz), joined on INFO/ALLELEID,
exactly as in build_clinvar_pairs.py: CHROM/POS/REF/ALT plus the CLNHGVS g.HGVS
kept for reference. Scoring against the VCF coordinate (not the g.HGVS string) is
essential here: a submitted string may legitimately spell an indel differently
from ClinVar's normalised form.

Inclusion rule (documented, deliberately loose so submitted messiness survives):
a string is kept if it starts with a RefSeq/Ensembl transcript accession
(NM_/NR_/XM_/XR_/ENST + digits, version optional) and contains ``:c.`` or
``:n.``. Strings are kept verbatim apart from stripping outer whitespace and
internal tab/newline characters (TSV safety).

Dedup rule (documented): SCV rows are collapsed to unique
(AlleleID, submitted string) pairs; ``scv_count`` records how many SCVs
contributed each pair. Rationale: submitted HGVS is not deduplicated across SCVs
in ClinVar, and repeated identical submissions measure submitter volume, not
string diversity; distinct strings for the same variant (e.g. two labs citing
different transcript versions) are all kept.

Sampling rule (documented): ``--sample N --seed S`` draws a uniform random sample
of N pairs with ``random.Random(S).sample`` and writes them in original file
order; the committed test sample and the benchmark sample both use seed 42.

Output TSV: header ``chrom pos ref alt g_hgvs c_hgvs scv_count`` - the VCF-format
pairs layout consumed by resolve_clinvar_pass.py (extra columns are ignored by its
reader). The full corpus is large and data-derived so it lives under a data dir
and is NOT committed; a 500-pair sample is committed to
tests/test_data/clinvar_hgvs/.

Usage:
    # fast path: consume an existing per-SCV CSV extraction
    python paper/scripts/build_clinvar_submitted_pairs.py \
        --scv-csv-dir extracted/ClinVarVCVRelease_2026-06 \
        clinvar.GRCh38.vcf.gz clinvar_submitted_pairs.GRCh38.tsv

    # self-contained: stream the VCV XML release (slow; --limit for a smoke test)
    python paper/scripts/build_clinvar_submitted_pairs.py \
        --xml ClinVarVCVRelease_2026-06.xml.gz \
        clinvar.GRCh38.vcf.gz clinvar_submitted_pairs.GRCh38.tsv

    # sample an existing corpus (no rebuild)
    python paper/scripts/build_clinvar_submitted_pairs.py \
        --from-pairs clinvar_submitted_pairs.GRCh38.tsv \
        --sample 500 --seed 42 --sample-out clinvar_submitted_500.tsv
"""
import argparse
import csv
import glob
import gzip
import random
import re
import sys
from pathlib import Path

# Reuse the VCF ground-truth loader so the two corpora share one join convention.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_clinvar_pairs import load_vcf_g, source_of  # noqa: E402

# Loose inclusion filter: transcript accession at the start (version optional,
# wrappers/messiness allowed) and a coding/non-coding kind somewhere after it.
_TX_START = re.compile(r"^(?:NM_|NR_|XM_|XR_|ENST)\d+")
_HAS_KIND = re.compile(r":[cn]\.")


def is_submitted_tx_hgvs(s):
    return bool(_TX_START.match(s)) and bool(_HAS_KIND.search(s))


def sanitize(s):
    """Strip outer whitespace and TSV-breaking characters, keep the rest verbatim."""
    return s.strip().replace("\t", " ").replace("\n", " ").replace("\r", " ")


def iter_scv_csv(csv_dir):
    """Yield (allele_id, submitted_hgvs) from an extract_xml_to_csv.py CSV dir.

    The extraction's ``hgvs`` column holds the first HGVS attribute value of the
    SCV (which may be genomic or protein; the transcript filter is applied here)."""
    files = sorted(glob.glob(str(Path(csv_dir) / "*.csv")))
    if not files:
        sys.exit(f"no CSVs found in {csv_dir}")
    for fn in files:
        with open(fn, newline="") as fh:
            for row in csv.DictReader(fh):
                h = sanitize(row.get("hgvs") or "")
                allele_id = row.get("allele_id") or ""
                if h and allele_id:
                    yield allele_id, h


def iter_scv_xml(xml_path, limit=0):
    """Yield (allele_id, submitted_hgvs) by streaming a ClinVar VCV XML release.

    One VariationArchive at a time (iterparse + clear), so memory stays flat.
    Yields the FIRST transcript c./n. expression per ClinicalAssertion."""
    try:
        from lxml import etree
    except ImportError:  # pragma: no cover - stdlib fallback
        import xml.etree.ElementTree as etree
    opener = gzip.open if str(xml_path).endswith(".gz") else open
    n_va = 0
    with opener(xml_path, "rb") as fh:
        for _event, elem in etree.iterparse(fh, events=("end",), tag="VariationArchive"):
            n_va += 1
            record = elem.find("ClassifiedRecord")
            if record is None:  # pre-Aug-2025 releases
                record = elem.find("InterpretedRecord")
            if record is not None:
                sa = record.find("SimpleAllele")
                allele_id = sa.get("AlleleID") if sa is not None else None
                if allele_id:
                    for ca in record.iter("ClinicalAssertion"):
                        ca_sa = ca.find("SimpleAllele")
                        if ca_sa is None:
                            continue
                        for attr in ca_sa.iter("Attribute"):
                            if attr.get("Type") == "HGVS" and attr.text:
                                h = sanitize(attr.text)
                                if is_submitted_tx_hgvs(h):
                                    yield allele_id, h
                                    break
            elem.clear()
            while elem.getprevious() is not None:  # lxml only; frees siblings
                del elem.getparent()[0]
            if limit and n_va >= limit:
                return


def write_sample(pairs_path, n, seed, out_path):
    """Uniform random sample of n data lines, written in original order."""
    with open(pairs_path) as fh:
        header = fh.readline()
        lines = fh.readlines()
    if n > len(lines):
        sys.exit(f"--sample {n} > corpus size {len(lines)}")
    idx = sorted(random.Random(seed).sample(range(len(lines)), n))
    with open(out_path, "w") as out:
        out.write(header)
        for i in idx:
            out.write(lines[i])
    print(f"  sampled {n}/{len(lines):,} (seed {seed}) -> {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--scv-csv-dir", help="dir of per-SCV CSVs (extract_xml_to_csv.py output)")
    src.add_argument("--xml", help="ClinVar VCV XML release (.xml.gz); slow streaming parse")
    src.add_argument("--from-pairs", help="existing corpus TSV; skip the build and only sample")
    ap.add_argument("vcf", nargs="?", help="ClinVar VCF (ground truth, joined on ALLELEID)")
    ap.add_argument("out", nargs="?", help="output corpus TSV")
    ap.add_argument("--limit", type=int, default=0, help="(--xml) stop after N VariationArchives (smoke test)")
    ap.add_argument("--sample", type=int, default=0, help="also write a random sample of N pairs")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--sample-out", help="path for the sample TSV")
    args = ap.parse_args()

    if args.from_pairs:
        if not (args.sample and args.sample_out):
            sys.exit("--from-pairs needs --sample and --sample-out")
        write_sample(args.from_pairs, args.sample, args.seed, args.sample_out)
        return
    if not (args.vcf and args.out):
        sys.exit("vcf and out are required when building")

    print("reading VCF ground truth...", file=sys.stderr)
    g_by_allele = load_vcf_g(args.vcf)
    print(f"  {len(g_by_allele):,} alleles with genomic coordinate", file=sys.stderr)

    scv_iter = (iter_scv_xml(args.xml, args.limit) if args.xml
                else iter_scv_csv(args.scv_csv_dir))

    n_scv = n_tx = n_joined = 0
    comp = {"refseq": 0, "ensembl": 0}
    pairs = {}  # (allele_id, c_hgvs) -> scv_count
    for allele_id, h in scv_iter:
        n_scv += 1
        if not is_submitted_tx_hgvs(h):
            continue
        n_tx += 1
        if allele_id not in g_by_allele:
            continue
        n_joined += 1
        key = (allele_id, h)
        if key not in pairs:
            comp[source_of(h.split(":", 1)[0])] += 1
        pairs[key] = pairs.get(key, 0) + 1

    with open(args.out, "w") as out:
        out.write("chrom\tpos\tref\talt\tg_hgvs\tc_hgvs\tscv_count\n")
        for (allele_id, c_hgvs), count in pairs.items():
            chrom, pos, ref, alt, g_hgvs = g_by_allele[allele_id]
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{g_hgvs}\t{c_hgvs}\t{count}\n")

    n_pairs = len(pairs)
    print(f"scanned {n_scv:,} SCV HGVS values -> {n_tx:,} transcript c./n. strings "
          f"-> {n_joined:,} with VCF ground truth -> {n_pairs:,} unique "
          f"(AlleleID, string) pairs", file=sys.stderr)
    if n_pairs:
        print(f"  source mix: refseq {comp['refseq']:,} ({100*comp['refseq']/n_pairs:.2f}%) "
              f"ensembl {comp['ensembl']:,} ({100*comp['ensembl']/n_pairs:.2f}%)", file=sys.stderr)
    print(f"Written: {args.out}", file=sys.stderr)

    if args.sample:
        if not args.sample_out:
            sys.exit("--sample needs --sample-out")
        write_sample(args.out, args.sample, args.seed, args.sample_out)


if __name__ == "__main__":
    main()
