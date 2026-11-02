"""
Tests for cdot.hgvs.version_safety (the pure ``is_version_substitution_safe``
coordinate-safety check on two transcript records, issue #28) and for
``cdot.hgvs.gene_hgvs.check_version_substitution``, which fetches those records
from a data provider.

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
from cdot.hgvs.gene_hgvs import check_version_substitution
from cdot.hgvs.version_safety import (
    cds_alignment_gaps,
    cds_genomic_map,
    describe_structure_change,
    exon_alignment_gaps,
    genomic_bracket_safe,
    intrinsic_cds_structure,
    intrinsic_exon_structure,
    is_version_substitution_safe,
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
    requested = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    substitute = _tx("NM_000059.4", 5, 35, {"GRCh38": [[95, 110, 0, 1, 15, None],
                                                       [200, 220, 1, 16, 35, None]]})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38",
                                                cited_position="5del")
    assert safe, reason
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38",
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
    requested = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    substitute = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert safe, reason


def test_substitution_minus_strand_moved_unsafe():
    requested = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    substitute = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS_MOVED}, strand="-")
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_genomic_bracket_minus_strand():
    # Minus-strand flanks agree → safe, and a moved flank → unsafe.
    below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    safe, reason = genomic_bracket_safe(below, above, "GRCh38")
    assert safe, reason
    above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MINUS_MOVED}, strand="-")
    safe, reason = genomic_bracket_safe(below, above, "GRCh38")
    assert not safe


def test_substitution_strand_change_unsafe():
    # Requested only in GRCh37 (so no same-build placement check), same layout but
    # on the other strand.
    requested = _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS})
    substitute = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
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
    requested = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    substitute = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_substitution_unsafe_genomic_replacement_same_build():
    # Identical intrinsic CDS structure but the CDS is placed at a different genomic
    # locus in the same build (mode-2 re-placement, eg NM_018263). Recomputed from the
    # loaded annotation: both versions are present, so the placement difference is seen
    # directly and the substitution is refused.
    requested = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})
    substitute = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_substitution_unsafe_gap_differs_same_structure():
    # Identical intrinsic CDS structure, but one version has an alignment gap the
    # other lacks: the coding coordinate would shift, so it must not be called safe.
    requested = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    substitute = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_GAP})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "gap" in reason


def test_substitution_safe_when_gaps_match():
    # Both versions carry the same alignment gap → coding coordinates agree → safe.
    requested = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_GAP})
    substitute = _tx("NM_000059.5", 0, 30, {"GRCh38": EXONS_GAP})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_blocklist_refuses_known_replacement_from_record_ids(monkeypatch):
    # The pure check reads accession and versions from the records' ids, so a
    # blocklisted pair is refused even though the structures are identical.
    import cdot.hgvs._version_replacement_blocklist as bl
    monkeypatch.setattr(bl, "KNOWN_REPLACEMENT_SUBSTITUTIONS",
                        frozenset({("NM_000059", 2, 4)}))
    requested = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    substitute = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "blocklist" in reason


def test_substitution_unsafe_cds_length_change():
    requested = _tx("NM_000059.5", 0, 33, {"GRCh38": EXONS_LONGER})
    substitute = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "CDS length changed" in reason


def test_substitution_decidable_from_other_build():
    # Requested .2 lives ONLY in GRCh37; we want to map in GRCh38 (where .3/.4
    # live). Its structure is read from GRCh37 and matches → safe.
    requested = _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS})
    substitute = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_substitution_other_build_structure_differs():
    requested = _tx("NM_000059.2", 0, 33, {"GRCh37": EXONS_LONGER})
    substitute = _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "CDS length changed" in reason


# ---- genomic_bracket_safe: pure, on the two flanking records ----------------

def test_genomic_bracket_safe_when_flanks_agree():
    # Flanking .2 and .4 share the same genomic CDS map in GRCh38 → bracket
    # agrees → safe (probabilistic).
    below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})
    safe, reason = genomic_bracket_safe(below, above, "GRCh38")
    assert safe
    assert "bracket" in reason


def test_genomic_bracket_unsafe_when_flanks_disagree():
    # Flanking .2 and .4 sit on different genomic coordinates → bracket disagrees.
    below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    above = _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED})
    safe, reason = genomic_bracket_safe(below, above, "GRCh38")
    assert not safe
    assert "may have moved" in reason


def test_genomic_bracket_unverified_without_bracket():
    # Nothing on one side → cannot bracket.
    below = _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS})
    safe, reason = genomic_bracket_safe(below, None, "GRCh38")
    assert not safe
    assert "cannot bracket" in reason


# ---------------------------------------------------------------------------
# check_version_substitution: fetches the records from a provider and picks
# the structural check or, for a version held in no build, the bracket
# ---------------------------------------------------------------------------

def test_check_substitution_both_present_uses_structural_check():
    dp = _provider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED}),
    })
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 4, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_check_substitution_substitute_missing():
    dp = _provider({"NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS})})
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 4, "GRCh38")
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
    safe, reason = check_version_substitution(dp, "NM_000059", 2, 4, "GRCh38")
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
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 5, "GRCh38")
    assert safe
    assert "bracket" in reason


def test_check_substitution_absent_version_cannot_bracket():
    # Requested .9 absent and nothing above it → cannot bracket.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = check_version_substitution(dp, "NM_000059", 9, 4, "GRCh38")
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
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 4)
    assert not safe
    assert "cannot verify" in reason


def test_check_substitution_provider_cannot_enumerate_versions():
    class _NoVersions(JSONDataProvider):
        def get_tx_versions(self, accession, genome_build=None):
            raise NotImplementedError("nope")

    data = {"cdot_version": "0.2.33", "genome_builds": ["GRCh38"],
            "transcripts": {"NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS})}}
    dp = _NoVersions([io.StringIO(json.dumps(data))])
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 4, "GRCh38")
    assert not safe
    assert "cannot verify" in reason


def test_provider_method_delegates_and_is_deprecated():
    # The provider method shipped in 0.2.30 (#120); it now delegates to
    # check_version_substitution and warns so callers migrate (#186).
    dp = _provider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    with pytest.warns(DeprecationWarning, match="check_version_substitution"):
        safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4, "GRCh38")
    assert (safe, reason) == check_version_substitution(dp, "NM_000059", 3, 4, "GRCh38")
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
    store could hand back) reaches check_version_substitution."""
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
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 4, "GRCh38")
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
    safe, reason = check_version_substitution(dp, "NM_000059", 3, 4, "GRCh38")
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
    requested = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38", cited_position="12A>G")
    assert safe, reason
    assert "identical exon structure" in reason


def test_noncoding_minus_strand_identical_safe():
    requested = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert safe, reason


def test_noncoding_minus_strand_moved_unsafe():
    requested = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MINUS_MOVED}, strand="-")
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "placement" in reason


def test_noncoding_exon_structure_changed_unsafe():
    requested = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_RESPLIT})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "exon structure changed" in reason
    assert "absent" not in reason


def test_noncoding_gap_differs_unsafe():
    requested = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_GAP})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "gap" in reason


def test_noncoding_cross_build_structure_safe():
    # Requested .3 only in GRCh37: structure (build-independent) and gaps match, so
    # it is safe even though the genomic placement can't be compared.
    requested = _nc_tx("NR_003051.3", {"GRCh37": NC_EXONS})
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MOVED})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert safe, reason


def test_noncoding_strand_change_unsafe():
    requested = _nc_tx("NR_003051.3", {"GRCh37": NC_EXONS})
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS_MINUS}, strand="-")
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "strand changed" in reason


def test_coding_vs_noncoding_unsafe():
    requested = _nc_tx("NR_003051.3", {"GRCh38": NC_EXONS})
    substitute = _tx("NR_003051.4", 0, 30, {"GRCh38": NC_EXONS})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "only .4 is coding" in reason


def test_noncoding_malformed_does_not_crash():
    malformed = [[100, 110, 0, None, None, None]]
    requested = _nc_tx("NR_003051.3", {"GRCh38": malformed})
    substitute = _nc_tx("NR_003051.4", {"GRCh38": NC_EXONS})
    safe, reason = is_version_substitution_safe(requested, substitute, "GRCh38")
    assert not safe
    assert "version .3 has a malformed exon record" in reason

