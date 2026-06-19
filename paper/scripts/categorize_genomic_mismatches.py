#!/usr/bin/env python3
"""Categorise the ClinVar genomic-match failures, the residual `incorrect` bucket (R6).

Most c->g conversions reproduce ClinVar's CLNHGVS exactly; a small residual does
not. A reviewer will ask what those are. The strong prior is that most are
*representation* differences (3'-shift vs ClinVar's normalisation, equivalent
del/dup spellings, reference-base disagreements), not cdot coordinate errors, so
categorising them likely shows the true accuracy is higher than the raw match rate.

Consumes the per-variant table from resolve_clinvar_pass.py (paper infra A1) and
re-examines only its `incorrect` rows, where cdot resolved the c.HGVS to some
genomic coordinate (`converted_g`) that differs from ClinVar's ground truth
(`g_hgvs`). Each is bucketed with DETERMINISTIC rules (no LLM):

    normalization_equivalent  the two g.HGVS normalise to the same string
                              (3'-shift / equivalent del/dup/inv spelling)
    ref_base                  same coordinate and alt, the reference base differs
                              (a reference-allele disagreement, not a cdot error)
    build                     different contig accession (assembly mismatch)
    gap_related               positions still differ after normalisation AND the
                              cited transcript has an alignment gap in this build
                              (a gap-representation difference)
    genuine_disagreement      a real coordinate disagreement after all of the above
    unparseable               one side would not parse / normalise

The headline is the corrected genuine-error rate: genuine_disagreement over all
resolved variants, i.e. what the match rate understates by counting representation
differences as errors. Written to output/facts/genomic_mismatch.csv.

This replaces the stale PyHGVS/GRCh37 investigate_fails.py stub.

Per CLAUDE.md: never run against full datasets in development. Verify the logic
with a small table built from tests/test_data/clinvar_hgvs/; the real breakdown is
a dedicated run over the full-ClinVar A1 table.

Usage:
    python paper/scripts/categorize_genomic_mismatches.py \
        output/clinvar_pass/refseq.csv --rest \
        --out output/facts/genomic_mismatch.csv
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import hgvs.parser
import hgvs.normalizer
from hgvs.exceptions import HGVSError

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_resolution import build_provider  # noqa: E402

CATEGORIES = ["normalization_equivalent", "ref_base", "build",
              "gap_related", "genuine_disagreement", "unparseable"]


def transcript_has_gap(provider, tx, version, build):
    """True if the cited transcript's alignment carries a gap in this build.

    A gapped alignment is where cdot's genomic projection can legitimately differ
    from ClinVar's representation of the same edit. Best-effort: unknown transcript
    or build -> False (falls through to genuine_disagreement)."""
    if not tx or version in (None, ""):
        return False
    try:
        record = provider._get_transcript(f"{tx}.{int(version)}")
    except Exception:  # noqa: BLE001
        return False
    if not record:
        return False
    gb = record.get("genome_builds", {}).get(build)
    if not gb:
        return False
    # exon = [alt_start, alt_end, exon_id, cds_start, cds_end, gap]
    return any(len(exon) > 5 and exon[5] for exon in gb.get("exons", []))


def categorize(hp, nz, provider, build, g_hgvs, converted_g, tx, version):
    """Return (category, normalized_clinvar, normalized_cdot) for one incorrect row."""
    try:
        v_clin = hp.parse_hgvs_variant(g_hgvs)
        v_cdot = hp.parse_hgvs_variant(converted_g)
    except HGVSError as e:
        logging.debug("parse failed %s / %s: %s", g_hgvs, converted_g, e)
        return "unparseable", "", ""

    if v_clin.ac != v_cdot.ac:
        return "build", str(v_clin), str(v_cdot)

    try:
        n_clin = nz.normalize(v_clin)
        n_cdot = nz.normalize(v_cdot)
    except HGVSError as e:
        logging.debug("normalize failed %s / %s: %s", g_hgvs, converted_g, e)
        return "unparseable", str(v_clin), str(v_cdot)

    s_clin, s_cdot = str(n_clin), str(n_cdot)
    if s_clin == s_cdot:
        return "normalization_equivalent", s_clin, s_cdot

    # Same coordinate span and same alternate, but the reference base differs ->
    # a reference-allele disagreement (ClinVar's submitted ref vs the assembly).
    pe_clin, pe_cdot = n_clin.posedit, n_cdot.posedit
    same_span = (str(pe_clin.pos) == str(pe_cdot.pos))
    edit_clin, edit_cdot = pe_clin.edit, pe_cdot.edit
    alt_clin = getattr(edit_clin, "alt", None)
    alt_cdot = getattr(edit_cdot, "alt", None)
    ref_clin = getattr(edit_clin, "ref", None)
    ref_cdot = getattr(edit_cdot, "ref", None)
    if same_span and alt_clin == alt_cdot and ref_clin != ref_cdot:
        return "ref_base", s_clin, s_cdot

    if transcript_has_gap(provider, tx, version, build):
        return "gap_related", s_clin, s_cdot

    return "genuine_disagreement", s_clin, s_cdot


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="A1 per-variant table (resolve_clinvar_pass.py output)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rest", action="store_true", help="cdot REST (default)")
    g.add_argument("--rest-insecure", action="store_true")
    g.add_argument("--json", help="cdot JSON.gz file")
    g.add_argument("--uta", action="store_true")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--fasta", nargs="+", help="genome FASTA(s) for local sequence fetching")
    ap.add_argument("--out", default="output/facts/genomic_mismatch.csv")
    ap.add_argument("--show", action="store_true", help="print each categorised mismatch")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.WARNING)
    # build_provider reads these even when the relevant mode is off.
    args.replace_reference = False

    df = pd.read_csv(args.results, dtype=str).fillna("")
    total_resolved = int((df["bucket"].isin(["correct", "incorrect"])).sum())
    n_correct = int((df["bucket"] == "correct").sum())
    incorrect = df[df["bucket"] == "incorrect"]
    print(f"== Categorising {len(incorrect)} incorrect of {total_resolved} resolved "
          f"({n_correct} correct) ==")

    provider, label = build_provider(args)
    hp = hgvs.parser.Parser()
    nz = hgvs.normalizer.Normalizer(provider)

    counts = {c: 0 for c in CATEGORIES}
    for _, row in incorrect.iterrows():
        cat, n_clin, n_cdot = categorize(
            hp, nz, provider, args.build,
            row["g_hgvs"], row["converted_g"], row["tx"], row["version"])
        counts[cat] += 1
        if args.show:
            print(f"  {cat:<24} clinvar={n_clin or row['g_hgvs']}  cdot={n_cdot or row['converted_g']}")

    n_incorrect = len(incorrect)
    genuine = counts["genuine_disagreement"]
    # Corrected error rate: genuine disagreements over all resolved variants. The
    # raw match rate counts every representation difference as an error too.
    raw_error_pct = round(100 * n_incorrect / total_resolved, 3) if total_resolved else 0.0
    genuine_error_pct = round(100 * genuine / total_resolved, 3) if total_resolved else 0.0

    facts = {
        "n_resolved":          [total_resolved],
        "n_incorrect":         [n_incorrect],
        "raw_error_pct":       [raw_error_pct],
        "genuine_error_pct":   [genuine_error_pct],
        **{f"n_{c}": [counts[c]] for c in CATEGORIES},
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(facts).to_csv(out, index=False)

    print(f"  provider          : {label}")
    for c in CATEGORIES:
        share = f"({100 * counts[c] / n_incorrect:.1f}% of incorrect)" if n_incorrect else ""
        print(f"  {c:<24}: {counts[c]:>5}  {share}")
    print(f"  raw error rate    : {raw_error_pct}%  ({n_incorrect}/{total_resolved})")
    print(f"  genuine error rate: {genuine_error_pct}%  ({genuine}/{total_resolved})")
    print(f"  Written: {out}")


if __name__ == "__main__":
    main()
