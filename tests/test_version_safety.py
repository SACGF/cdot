"""
Tests for cdot.hgvs.version_safety: the pure ``is_version_substitution_safe``
coordinate-safety check on two transcript records (issue #28) and
``is_version_substitution_safe_by_accession``, which fetches those records from a
data provider.

All transcript records here are *synthesised* from public accessions (BRCA2
NM_000059, RUNX1 NM_001754) with toy coordinates — no real release data is read,
so the structural arithmetic is checked exactly. The intrinsic CDS structure is
build-independent (CDS length + coding-exon-segment lengths in transcript
coordinates), so two versions on different genomic placements but the same
transcript layout must compare equal.
"""
import io
import json

import pytest

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider, LocalDataProvider
from cdot.hgvs.version_safety import (
    _genomic_at,
    cds_alignment_gaps,
    cds_genomic_map,
    describe_structure_change,
    exon_alignment_gaps,
    genomic_bracket_safe,
    genomic_maps_fully_agree,
    intrinsic_cds_structure,
    intrinsic_exon_structure,
    is_version_substitution_safe,
    is_version_substitution_safe_by_accession,
    utr_lengths,
    utr_regions_touched,
)


# ---------------------------------------------------------------------------
# Synthetic record helpers
# ---------------------------------------------------------------------------

def _build(exons, strand="+", contig="NC_000013.11"):
    return {"contig": contig, "strand": strand, "url": None,
            "cds_start": None, "cds_end": None, "exons": exons}


def _tx(acc, start_codon, stop_codon, builds, gene_name="BRCA2", strand="+"):
    """builds: {build_name: exons-list}. start_codon is the 0-based count of
    pre-CDS bases (cdot convention); CDS spans transcript [start_codon+1,
    stop_codon]."""
    return {
        "id": acc, "gene_name": gene_name,
        "start_codon": start_codon, "stop_codon": stop_codon,
        "genome_builds": {name: _build(exons, strand=strand) for name, exons in builds.items()},
    }


def _nc_tx(acc, builds, gene_name="RUNX1", strand="+"):
    """A non-coding (NR_) record: no start/stop codon."""
    return {
        "id": acc, "gene_name": gene_name,
        "genome_builds": {name: _build(exons, strand=strand) for name, exons in builds.items()},
    }


# Two coding exons, split at transcript position 11; CDS = transcript 1..30.
# Structure: strand '+', cds_len 30, segment lengths (10, 20).
EXONS = [[100, 110, 0, 1, 10, None], [200, 220, 1, 11, 30, None]]
# Same transcript layout, different genomic placement (a coordinate-preserving bump).
EXONS_MOVED = [[900, 910, 0, 1, 10, None], [2000, 2020, 1, 11, 30, None]]
# CDS extended by 3 bp (stop_codon 33): a CDS-length change → structure differs.
EXONS_LONGER = [[900, 910, 0, 1, 10, None], [2000, 2023, 1, 11, 33, None]]
# Same CDS length (30) but the exon boundary moved to transcript position 16:
# segment lengths (15, 15) instead of (10, 20) → structure differs.
EXONS_RESPLIT = [[100, 115, 0, 1, 15, None], [200, 215, 1, 16, 30, None]]
# Identical transcript/CDS structure to EXONS (cds_start/cds_end unchanged) but the
# second coding exon carries a 1-base alignment gap, so its genomic span is one
# longer and every coding base downstream of the gap shifts. Intrinsic structure
# compares equal; only the gap signature differs.
EXONS_GAP = [[100, 110, 0, 1, 10, None], [200, 221, 1, 11, 30, "M5 D1 M15"]]


# ---------------------------------------------------------------------------
# intrinsic_cds_structure
# ---------------------------------------------------------------------------

def test_intrinsic_structure_basic():
    rec = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    assert intrinsic_cds_structure(rec) == ("+", 30, (10, 20))


def test_intrinsic_structure_build_independent():
    # Same transcript layout on different genomic coordinates → identical structure.
    a = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    b = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED})
    assert intrinsic_cds_structure(a) == intrinsic_cds_structure(b)


def test_intrinsic_structure_cds_length_change_differs():
    a = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    b = _tx("NM_000059.5", 0, 33, {"GRCh38": EXONS_LONGER})
    assert intrinsic_cds_structure(a) != intrinsic_cds_structure(b)


def test_intrinsic_structure_resplit_same_length_differs():
    a = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    b = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_RESPLIT})
    sa, sb = intrinsic_cds_structure(a), intrinsic_cds_structure(b)
    assert sa[1] == sb[1] == 30          # same CDS length
    assert sa != sb                       # but different exon split


def test_intrinsic_structure_reads_from_other_build_when_target_absent():
    # Record only carries GRCh37; asking for GRCh38 still returns the (build-
    # independent) structure from GRCh37.
    rec = _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS})
    assert intrinsic_cds_structure(rec, "GRCh38") == ("+", 30, (10, 20))


def test_intrinsic_structure_none_for_noncoding():
    rec = {"id": "NR_1.1", "gene_name": "X",
           "genome_builds": {"GRCh38": _build(EXONS)}}  # no start/stop codon
    assert intrinsic_cds_structure(rec) is None


def test_cds_alignment_gaps_clean_is_empty():
    rec = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    assert cds_alignment_gaps(rec) == ()


def test_cds_alignment_gaps_reports_coding_exon_gap():
    rec = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_GAP})
    # gap on the second coding exon (CDS-segment index 1)
    assert cds_alignment_gaps(rec) == ((1, "M5 D1 M15"),)


def test_intrinsic_structure_blind_to_gap_but_signature_differs():
    # The whole point: intrinsic structure is identical, the gap signature is not.
    clean = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    gapped = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_GAP})
    assert intrinsic_cds_structure(clean) == intrinsic_cds_structure(gapped)
    assert cds_alignment_gaps(clean) != cds_alignment_gaps(gapped)


def test_utr_lengths():
    # start_codon=5 (5'UTR len 5), stop_codon=25, total cDNA = 30 → 3'UTR len 5.
    rec = _tx("NM_000059.4", 5, 25, {"GRCh38": EXONS})
    assert utr_lengths(rec) == (5, 5)


def test_utr_lengths_none_for_noncoding():
    rec = {"id": "NR_1.1", "gene_name": "X",
           "genome_builds": {"GRCh38": _build(EXONS)}}
    assert utr_lengths(rec) is None


@pytest.mark.parametrize("cited, expected", [
    ("-49=", {"five_prime"}),
    ("-49+10A>G", {"five_prime"}),       # 5'UTR intronic
    ("-3_5del", {"five_prime"}),         # range spanning into the CDS
    ("*100A>G", {"three_prime"}),
    ("*5_*10del", {"three_prime"}),
    ("123A>G", set()),                   # coding
    ("123+5G>A", set()),                 # CDS intronic (+ offset)
    ("2604-1G>C", set()),                # CDS intronic (- offset, not 5'UTR)
    ("123_124insAT", set()),
    ("", set()),
    # Uncertain / bracketed ranges (#158)
    ("(-30_-20)del", {"five_prime"}),
    ("(-30_5)del", {"five_prime"}),      # uncertain range starting in the 5'UTR
    ("(?_-30)_(5_?)del", {"five_prime"}),
    ("(-30_5)_(10_20)del", {"five_prime"}),
    ("(*5_*10)del", {"three_prime"}),
    ("(?_-30)_(*5_?)del", {"five_prime", "three_prime"}),
    ("(2604-10_2604-1)del", set()),      # uncertain CDS-intronic range
    ("(123_130)_(200_210)del", set()),
])
def test_utr_regions_touched(cited, expected):
    assert utr_regions_touched(cited) == expected


def test_describe_structure_change_cds_length():
    a = ("+", 30, (10, 20))
    b = ("+", 33, (10, 23))
    assert "CDS length changed 30→33" in describe_structure_change(a, b)


def test_uncertain_five_prime_range_checks_utr_length():
    # .4 adds 5 bases of 5'UTR but keeps the CDS and its genomic placement, so a
    # coding position is safe. An uncertain range starting in the 5'UTR
    # ("(-4_2)del") must still be caught as a 5'UTR position (#158).
    transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    substitute_transcript = _tx("NM_000059.4", 5, 35, {"GRCh38": [[95, 110, 0, 1, 15, None],
                                                       [200, 220, 1, 16, 35, None]]})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38",
                                                cited_position="5del")
    assert safe, reason
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38",
                                                cited_position="(-4_2)del")
    assert not safe
    assert "5'UTR" in reason


def test_describe_structure_change_strand():
    a = ("+", 30, (10, 20))
    b = ("-", 30, (10, 20))
    assert "strand changed +→-" in describe_structure_change(a, b)


def test_describe_structure_change_same_length():
    a = ("+", 30, (10, 20))
    b = ("+", 30, (15, 15))
    assert "coding-exon structure changed" in describe_structure_change(a, b)


# ---------------------------------------------------------------------------
# Malformed records (#158): None instead of a crash
# ---------------------------------------------------------------------------

# Exon transcript coordinates are None (as in the versionless records built in
# tests/test_gene_hgvs.py).
EXONS_MALFORMED = [[100, 110, 0, None, None, None], [200, 220, 1, 11, 30, None]]
# Exons with no trailing gap element: treated as gap-free, not malformed.
EXONS_NO_GAP_FIELD = [[100, 110, 0, 1, 10], [200, 220, 1, 11, 30]]


@pytest.mark.parametrize("func", [
    intrinsic_cds_structure, cds_alignment_gaps, utr_lengths,
    intrinsic_exon_structure, exon_alignment_gaps, cds_genomic_map,
])
def test_malformed_record_returns_none(func):
    rec = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MALFORMED})
    assert func(rec, "GRCh38") is None


def test_exons_without_gap_field_are_consistent():
    rec = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_NO_GAP_FIELD})
    assert intrinsic_cds_structure(rec) == ("+", 30, (10, 20))
    assert cds_alignment_gaps(rec) == ()
    assert cds_genomic_map(rec, "GRCh38")["has_gap"] is False


# ---------------------------------------------------------------------------
# Minus strand
# ---------------------------------------------------------------------------

# Same transcript layout as EXONS (segments 10, 20) on the minus strand. Exons are
# listed in genomic order, so transcript exon 1 (tx 1..10) is the higher one.
EXONS_MINUS = [[100, 120, 1, 11, 30, None], [200, 210, 0, 1, 10, None]]
EXONS_MINUS_MOVED = [[1100, 1120, 1, 11, 30, None], [1200, 1210, 0, 1, 10, None]]


def test_cds_genomic_map_minus_strand():
    rec = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    sig = cds_genomic_map(rec, "GRCh38")
    assert sig["strand"] == "-" and sig["sign"] == -1
    # c.1 (offset 0) is the last base of the higher exon; c.11 (offset 10) the last
    # base of the lower one.
    assert sig["segs"] == [(0, 209), (10, 119)]
    assert intrinsic_cds_structure(rec) == ("-", 30, (10, 20))


def test_cds_alignment_gaps_minus_strand_in_cds_order():
    exons = [[100, 121, 1, 11, 30, "M5 D1 M15"], [200, 210, 0, 1, 10, None]]
    rec = _tx("NM_000059.4", 0, 30, {"GRCh38": exons}, strand="-")
    # The gapped exon is the second coding segment in transcript order.
    assert cds_alignment_gaps(rec) == ((1, "M5 D1 M15"),)


def test_substitution_minus_strand_same_placement_safe():
    transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    substitute_transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert safe, reason


def test_substitution_minus_strand_moved_unsafe():
    transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    substitute_transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS_MOVED}, strand="-")
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_genomic_bracket_minus_strand():
    # Minus-strand flanks agree → safe, and a moved flank → unsafe.
    transcript_below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    transcript_above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    safe, reason = genomic_bracket_safe(transcript_below, transcript_above, "GRCh38")
    assert safe, reason
    transcript_above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS_MOVED}, strand="-")
    safe, reason = genomic_bracket_safe(transcript_below, transcript_above, "GRCh38")
    assert not safe


def test_substitution_strand_change_unsafe():
    # Requested only in GRCh37 (so no same-build placement check), same layout but
    # on the other strand.
    transcript = _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS})
    substitute_transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "strand changed" in reason


# ---------------------------------------------------------------------------
# Provider: is_version_substitution_safe


# ---------------------------------------------------------------------------
# is_version_substitution_safe: pure, on the two versions' records
# ---------------------------------------------------------------------------

def _provider(transcripts, builds=("GRCh37", "GRCh38")):
    data = {"cdot_version": "0.2.33", "genome_builds": list(builds),
            "transcripts": transcripts}
    return JSONDataProvider([io.StringIO(json.dumps(data))])


def test_substitution_safe_identical_structure():
    # Same structure AND same genomic placement → safe.
    transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    substitute_transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_substitution_unsafe_genomic_replacement_same_build():
    # Identical intrinsic CDS structure but the CDS is placed at a different genomic
    # locus in the same build (mode-2 re-placement, eg NM_018263). Recomputed from the
    # loaded annotation: both versions are present, so the placement difference is seen
    # directly and the substitution is refused.
    transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    substitute_transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_substitution_unsafe_gap_differs_same_structure():
    # Identical intrinsic CDS structure, but one version has an alignment gap the
    # other lacks: the coding coordinate would shift, so it must not be called safe.
    transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    substitute_transcript = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_GAP})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "gap" in reason


def test_substitution_safe_when_gaps_match():
    # Both versions carry the same alignment gap → coding coordinates agree → safe.
    transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_GAP})
    substitute_transcript = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_GAP})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_blocklist_refuses_known_replacement_from_record_ids(monkeypatch):
    # The pure check reads accession and versions from the records' ids, so a
    # blocklisted pair is refused even though the structures are identical.
    import cdot.hgvs._version_replacement_blocklist as bl
    monkeypatch.setattr(bl, "KNOWN_REPLACEMENT_SUBSTITUTIONS",
                        frozenset({("NM_000059", 2, 4)}))
    transcript = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    substitute_transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "blocklist" in reason


def test_substitution_unsafe_cds_length_change():
    transcript = _tx("NM_000059.5", 0, 33, {"GRCh38": EXONS_LONGER})
    substitute_transcript = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "CDS length changed" in reason


def test_substitution_decidable_from_other_build():
    # Requested .2 lives ONLY in GRCh37; we want to map in GRCh38 (where .3/.4
    # live). Its structure is read from GRCh37 and matches → safe.
    transcript = _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS})
    substitute_transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_substitution_other_build_structure_differs():
    transcript = _tx("NM_000059.2", 0, 33, {"GRCh37": EXONS_LONGER})
    substitute_transcript = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "CDS length changed" in reason


# ---- genomic_bracket_safe: pure, on the two flanking records ----------------

def test_genomic_bracket_safe_when_flanks_agree():
    # Flanking .2 and .4 share the same genomic CDS map in GRCh38 → bracket
    # agrees → safe (probabilistic).
    transcript_below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    transcript_above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = genomic_bracket_safe(transcript_below, transcript_above, "GRCh38")
    assert safe
    assert "bracket" in reason


def test_genomic_bracket_unsafe_when_flanks_disagree():
    # Flanking .2 and .4 sit on different genomic coordinates → bracket disagrees.
    transcript_below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    transcript_above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = genomic_bracket_safe(transcript_below, transcript_above, "GRCh38")
    assert not safe
    assert "may have moved" in reason


def test_genomic_bracket_unverified_without_bracket():
    # Nothing on one side → cannot bracket.
    transcript_below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    safe, reason = genomic_bracket_safe(transcript_below, None, "GRCh38")
    assert not safe
    assert "cannot bracket" in reason


# ---------------------------------------------------------------------------
# is_version_substitution_safe_by_accession: fetches the records from a provider
# and picks the structural check or, for a version held in no build, the bracket
# ---------------------------------------------------------------------------

def test_check_substitution_both_present_uses_structural_check():
    dp = _provider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4", "GRCh38")
    assert not safe
    assert "placement" in reason


def test_check_substitution_substitute_missing():
    dp = _provider({"NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})})
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4", "GRCh38")
    assert not safe
    assert ".4 not found" in reason


def test_check_substitution_blocklist_applies_when_requested_absent(monkeypatch):
    # The shipped blocklist exists for the case the provider cannot recompute: the
    # requested version is absent from the loaded data. Here only .4 is present; the
    # blocklist still refuses substituting it for the absent .2.
    import cdot.hgvs._version_replacement_blocklist as bl
    monkeypatch.setattr(bl, "KNOWN_REPLACEMENT_SUBSTITUTIONS",
                        frozenset({("NM_000059", 2, 4)}))
    dp = _provider({"NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})})
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.2", "NM_000059.4", "GRCh38")
    assert not safe
    assert "blocklist" in reason


def test_check_substitution_absent_version_brackets_nearest_placed_neighbours():
    # Requested .3 absent everywhere. Nearest held versions placed in GRCh38 are
    # .2 and .5 (.4 is GRCh37-only and is skipped); they agree → safe (probabilistic).
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh37": EXONS_MOVED}),
        "NM_000059.5": _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.5", "GRCh38")
    assert safe
    assert "bracket" in reason


def test_check_substitution_absent_version_cannot_bracket():
    # Requested .9 absent and nothing above it → cannot bracket.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.9", "NM_000059.4", "GRCh38")
    assert not safe
    assert "cannot bracket" in reason


def test_check_substitution_unverified_without_genome_build():
    # No genome_build → the structural check can still run when the record is
    # present, but the bracket fallback can't; an absent requested version is
    # therefore unverified.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4")
    assert not safe
    assert "cannot verify" in reason


def test_check_substitution_provider_cannot_enumerate_versions():
    class _NoVersions(JSONDataProvider):
        def get_tx_versions(self, accession, genome_build=None):
            raise NotImplementedError("nope")

    data = {"cdot_version": "0.2.33", "genome_builds": ["GRCh38"],
            "transcripts": {"NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})}}
    dp = _NoVersions([io.StringIO(json.dumps(data))])
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4", "GRCh38")
    assert not safe
    assert "cannot verify" in reason


def test_by_accession_rejects_unversioned_or_mismatched_accessions():
    dp = _provider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059", "NM_000059.4", "GRCh38")
    assert not safe
    assert "versioned" in reason
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_001754.5", "GRCh38")
    assert not safe
    assert "different transcripts" in reason


def test_provider_method_delegates_and_is_deprecated():
    # The provider method shipped in 0.2.30 (#120); it now delegates to
    # is_version_substitution_safe_by_accession and warns so callers migrate (#186).
    dp = _provider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    with pytest.warns(DeprecationWarning, match="is_version_substitution_safe_by_accession"):
        safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4, "GRCh38")
    assert (safe, reason) == is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4", "GRCh38")
    assert safe


# ---------------------------------------------------------------------------
# get_tx_versions build filter
# ---------------------------------------------------------------------------

def test_get_tx_versions_build_filter():
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS}),
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh37": EXONS, "GRCh38": EXONS}),
    })
    assert dp.get_tx_versions("NM_000059") == [2, 3, 4]
    assert dp.get_tx_versions("NM_000059", "GRCh38") == [3, 4]
    assert dp.get_tx_versions("NM_000059", "GRCh37") == [2, 4]


# ---------------------------------------------------------------------------
# Malformed records through the provider (#158)
# ---------------------------------------------------------------------------

class _RecordsProvider(LocalDataProvider):
    """A LocalDataProvider over raw dict records, skipping the msgspec validation
    JSONDataProvider does on load, so a malformed record (as a Redis or versionless
    store could hand back) reaches is_version_substitution_safe_by_accession."""
    def __init__(self, records):
        self._records = records

    def _get_transcript(self, tx_ac):
        return self._records.get(tx_ac)

    def get_tx_versions(self, accession, genome_build=None):
        return sorted(int(k.rsplit(".", 1)[1]) for k in self._records
                      if k.rsplit(".", 1)[0] == accession)

    def _get_transcript_ids_for_gene(self, gene):
        raise NotImplementedError

    def _get_gene(self, gene):
        raise NotImplementedError

    def _get_contig_interval_tree(self, alt_ac):
        raise NotImplementedError


def test_substitution_malformed_requested_does_not_crash():
    dp = _RecordsProvider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MALFORMED}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4", "GRCh38")
    assert not safe
    # The requested version is present, so the reason must not say it is absent.
    assert "absent" not in reason
    assert "no usable CDS structure" in reason
    # The pure check on the records alone says the same, without the bracket.
    safe, reason = is_version_substitution_safe(dp._get_transcript("NM_000059.3"),
                                                dp._get_transcript("NM_000059.4"), "GRCh38")
    assert not safe
    assert "requested version .3 has no usable CDS structure" in reason


def test_substitution_malformed_substitute_does_not_crash():
    dp = _RecordsProvider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MALFORMED}),
    })
    safe, reason = is_version_substitution_safe_by_accession(dp, "NM_000059.3", "NM_000059.4", "GRCh38")
    assert not safe
    assert "substitute version .4 has no usable CDS structure" in reason


# ---------------------------------------------------------------------------
# Non-coding (n.) path (#158)
# ---------------------------------------------------------------------------

# Three exons, transcript lengths (10, 20, 5).
NC_EXONS = [[100, 110, 0, 1, 10, None], [200, 220, 1, 11, 30, None],
            [300, 305, 2, 31, 35, None]]
# Same layout, different genomic placement.
NC_EXONS_MOVED = [[1100, 1110, 0, 1, 10, None], [1200, 1220, 1, 11, 30, None],
                  [1300, 1305, 2, 31, 35, None]]
# Same total length, boundary moved (15, 15, 5).
NC_EXONS_RESPLIT = [[100, 115, 0, 1, 15, None], [200, 215, 1, 16, 30, None],
                    [300, 305, 2, 31, 35, None]]
# Same layout with an alignment gap in exon 2.
NC_EXONS_GAP = [[100, 110, 0, 1, 10, None], [200, 221, 1, 11, 30, "M5 D1 M15"],
                [300, 305, 2, 31, 35, None]]
# Minus strand: listed in genomic order, transcript exon 1 is the highest.
NC_EXONS_MINUS = [[100, 105, 2, 31, 35, None], [200, 220, 1, 11, 30, None],
                  [300, 310, 0, 1, 10, None]]
NC_EXONS_MINUS_MOVED = [[1100, 1105, 2, 31, 35, None], [1200, 1220, 1, 11, 30, None],
                        [1300, 1310, 0, 1, 10, None]]


def test_intrinsic_exon_structure():
    rec = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    assert intrinsic_exon_structure(rec) == ("+", (10, 20, 5))


def test_intrinsic_exon_structure_minus_strand_in_transcript_order():
    rec = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    assert intrinsic_exon_structure(rec) == ("-", (10, 20, 5))


def test_exon_alignment_gaps():
    rec = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS_GAP})
    assert exon_alignment_gaps(rec) == ((1, "M5 D1 M15"),)


def test_noncoding_identical_safe():
    transcript = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38", cited_position="12A>G")
    assert safe, reason
    assert "identical exon structure" in reason


def test_noncoding_minus_strand_identical_safe():
    transcript = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert safe, reason


def test_noncoding_minus_strand_moved_unsafe():
    transcript = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MINUS_MOVED}, strand="-")
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_noncoding_exon_structure_changed_unsafe():
    transcript = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_RESPLIT})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "exon structure changed" in reason
    assert "absent" not in reason


def test_noncoding_gap_differs_unsafe():
    transcript = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_GAP})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "gap" in reason


def test_noncoding_cross_build_structure_safe():
    # Requested .3 only in GRCh37: structure (build-independent) and gaps match, so
    # it is safe even though the genomic placement can't be compared.
    transcript = _nc_tx("NR_003051.3", {"GRCh37": NC_EXONS})
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MOVED})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert safe, reason


def test_noncoding_strand_change_unsafe():
    transcript = _nc_tx("NR_003051.3", {"GRCh37": NC_EXONS})
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "strand changed" in reason


def test_coding_vs_noncoding_unsafe():
    transcript = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute_transcript = _tx("NR_003051.4", 0, 30, {"GRCh38": NC_EXONS})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "only .4 is coding" in reason


def test_noncoding_malformed_does_not_crash():
    malformed = [[100, 110, 0, None, None, None]]
    transcript = _nc_tx("NR_003051.3", {"GRCh38": malformed})
    substitute_transcript = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS})
    safe, reason = is_version_substitution_safe(transcript, substitute_transcript, "GRCh38")
    assert not safe
    assert "version .3 has a malformed exon record" in reason


# ---------------------------------------------------------------------------
# _genomic_at and genomic_maps_fully_agree: the genomic CDS map arithmetic
# behind the placement check and the bracket fallback
# ---------------------------------------------------------------------------

def test_genomic_at_plus_strand_walks_each_segment():
    sig = cds_genomic_map(_tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}), "GRCh38")
    assert sig["segs"] == [(0, 100), (10, 200)]
    # first segment: c.1..c.10 at genomic 100..109
    assert _genomic_at(sig, 0) == 100
    assert _genomic_at(sig, 9) == 109
    # second segment: c.11..c.30 at genomic 200..219
    assert _genomic_at(sig, 10) == 200
    assert _genomic_at(sig, 29) == 219


def test_genomic_at_minus_strand_counts_down():
    sig = cds_genomic_map(_tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-"), "GRCh38")
    assert sig["segs"] == [(0, 209), (10, 119)]
    # c.1 is the last base of the higher exon, then genomic coordinates decrease
    assert _genomic_at(sig, 0) == 209
    assert _genomic_at(sig, 9) == 200
    assert _genomic_at(sig, 10) == 119
    assert _genomic_at(sig, 29) == 100


def test_genomic_at_single_segment():
    sig = cds_genomic_map(_tx("NM_000059.4", 0, 30, {"GRCh38": [[100, 130, 0, 1, 30, None]]}), "GRCh38")
    assert sig["segs"] == [(0, 100)]
    assert [_genomic_at(sig, off) for off in (0, 15, 29)] == [100, 115, 129]


def test_genomic_at_many_segments_picks_the_right_one():
    # Five one-codon exons so the binary search has several breakpoints to choose from
    exons = [[100 + 50 * i, 103 + 50 * i, i, 1 + 3 * i, 3 + 3 * i, None] for i in range(5)]
    sig = cds_genomic_map(_tx("NM_000059.4", 0, 15, {"GRCh38": exons}), "GRCh38")
    assert sig["segs"] == [(0, 100), (3, 150), (6, 200), (9, 250), (12, 300)]
    for off in range(15):
        assert _genomic_at(sig, off) == 100 + 50 * (off // 3) + off % 3


def _map(exons, strand="+", stop_codon=30, start_codon=0, contig="NC_000013.11"):
    rec = _tx("NM_000059.4", start_codon, stop_codon, {"GRCh38": exons}, strand=strand)
    rec["genome_builds"]["GRCh38"]["contig"] = contig
    return cds_genomic_map(rec, "GRCh38")


def test_genomic_maps_fully_agree_identical():
    assert genomic_maps_fully_agree(_map(EXONS), _map(EXONS))
    assert genomic_maps_fully_agree(_map(EXONS_MINUS, strand="-"), _map(EXONS_MINUS, strand="-"))


def test_genomic_maps_fully_agree_moved_placement():
    assert not genomic_maps_fully_agree(_map(EXONS), _map(EXONS_MOVED))
    assert not genomic_maps_fully_agree(_map(EXONS_MINUS, strand="-"), _map(EXONS_MINUS_MOVED, strand="-"))


def test_genomic_maps_fully_agree_different_exon_split_same_genomic_bases():
    # The second exon is split in two in b, but every coding base maps to the same
    # genomic coordinate, so the maps agree (b's extra breakpoint is checked in a)
    split = [[100, 110, 0, 1, 10, None], [200, 210, 1, 11, 20, None], [210, 220, 2, 21, 30, None]]
    assert genomic_maps_fully_agree(_map(EXONS), _map(split))
    assert genomic_maps_fully_agree(_map(split), _map(EXONS))


def test_genomic_maps_fully_agree_longer_cds_agrees_on_shared_bases():
    # Same placement, the CDS of b runs 3 bases further: compared on the shared 30
    longer = [[100, 110, 0, 1, 10, None], [200, 223, 1, 11, 33, None]]
    assert genomic_maps_fully_agree(_map(EXONS), _map(longer, stop_codon=33))
    assert genomic_maps_fully_agree(_map(longer, stop_codon=33), _map(EXONS))


def test_genomic_maps_disagree_after_shared_breakpoint():
    # First exon identical, second exon shifted by one base
    shifted = [[100, 110, 0, 1, 10, None], [201, 221, 1, 11, 30, None]]
    assert not genomic_maps_fully_agree(_map(EXONS), _map(shifted))


def test_genomic_maps_disagree_within_a_resplit():
    # Same CDS length, exon boundary moved: c.11..c.15 sit at 110..114 in b, 200..204 in a
    assert not genomic_maps_fully_agree(_map(EXONS), _map(EXONS_RESPLIT))


def test_genomic_maps_fully_agree_gap_is_non_comparable():
    assert not genomic_maps_fully_agree(_map(EXONS), _map(EXONS_GAP))
    assert not genomic_maps_fully_agree(_map(EXONS_GAP), _map(EXONS_GAP))


def test_genomic_maps_fully_agree_none_or_different_contig_or_strand():
    assert not genomic_maps_fully_agree(None, _map(EXONS))
    assert not genomic_maps_fully_agree(_map(EXONS), None)
    assert not genomic_maps_fully_agree(_map(EXONS), _map(EXONS, contig="NC_000012.12"))
    # EXONS_MINUS covers the same genomic bases as EXONS, on the other strand
    assert not genomic_maps_fully_agree(_map(EXONS), _map(EXONS_MINUS, strand="-"))


def test_cds_genomic_map_minus_strand_with_utrs():
    # Minus strand with a 5-base 5'UTR and a 3'UTR: the CDS starts part way into the
    # higher exon and ends before the lower exon does
    exons = [[100, 120, 1, 16, 35, None], [200, 215, 0, 1, 15, None]]
    sig = cds_genomic_map(_tx("NM_000059.4", 5, 30, {"GRCh38": exons}, strand="-"), "GRCh38")
    assert sig["cds_len"] == 25
    # c.1 (tx 6) is five bases down from the top of the higher exon (214)
    assert sig["segs"] == [(0, 209), (10, 119)]
    assert _genomic_at(sig, 24) == 105  # last CDS base (tx 30), five above the exon start at 100


def test_cds_genomic_map_skips_utr_only_exons():
    exons = [[50, 60, 0, 1, 10, None], [100, 110, 1, 11, 20, None], [200, 220, 2, 21, 40, None],
             [300, 310, 3, 41, 50, None]]
    sig = cds_genomic_map(_tx("NM_000059.4", 10, 40, {"GRCh38": exons}), "GRCh38")
    assert sig["cds_len"] == 30
    assert sig["segs"] == [(0, 100), (10, 200)]

