"""
Tests for cdot.hgvs.version_safety and the provider-level
``is_version_substitution_safe`` coordinate-safety check (issue #28).

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

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider
from cdot.hgvs.version_safety import (
    describe_structure_change,
    intrinsic_cds_structure,
)


# ---------------------------------------------------------------------------
# Synthetic record helpers
# ---------------------------------------------------------------------------

def _build(exons, strand="+", contig="NC_000013.11"):
    return {"contig": contig, "strand": strand, "url": None,
            "cds_start": None, "cds_end": None, "exons": exons}


def _tx(acc, start_codon, stop_codon, builds, gene_name="BRCA2"):
    """builds: {build_name: exons-list}. start_codon is the 0-based count of
    pre-CDS bases (cdot convention); CDS spans transcript [start_codon+1,
    stop_codon]."""
    return {
        "id": acc, "gene_name": gene_name,
        "start_codon": start_codon, "stop_codon": stop_codon,
        "genome_builds": {name: _build(exons) for name, exons in builds.items()},
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


def test_describe_structure_change_cds_length():
    a = ("+", 30, (10, 20))
    b = ("+", 33, (10, 23))
    assert "CDS length changed 30→33" in describe_structure_change(a, b)


def test_describe_structure_change_same_length():
    a = ("+", 30, (10, 20))
    b = ("+", 30, (15, 15))
    assert "coding-exon structure changed" in describe_structure_change(a, b)


# ---------------------------------------------------------------------------
# Provider: is_version_substitution_safe
# ---------------------------------------------------------------------------

def _provider(transcripts, builds=("GRCh37", "GRCh38")):
    data = {"cdot_version": "0.2.33", "genome_builds": list(builds),
            "transcripts": transcripts}
    return JSONDataProvider([io.StringIO(json.dumps(data))])


def test_substitution_safe_identical_structure():
    dp = _provider({
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_substitution_unsafe_cds_length_change():
    dp = _provider({
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.5": _tx("NM_000059.5", 0, 33, {"GRCh38": EXONS_LONGER}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 5, 4, "GRCh38")
    assert not safe
    assert "CDS length changed" in reason


def test_substitution_decidable_from_other_build():
    # Requested .2 lives ONLY in GRCh37; we want to map in GRCh38 (where .3/.4
    # live). Its structure is read from GRCh37 and matches → safe.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh37": EXONS}),
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MOVED}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 2, 3, "GRCh38")
    assert safe
    assert "identical CDS structure" in reason


def test_substitution_other_build_structure_differs():
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 33, {"GRCh37": EXONS_LONGER}),
        "NM_000059.3": _tx("NM_000059.3", 0, 30, {"GRCh38": EXONS_MOVED}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 2, 3, "GRCh38")
    assert not safe
    assert "CDS length changed" in reason


# ---- genomic-bracket fallback (requested version exists in no build) --------

def test_genomic_bracket_safe_when_flanks_agree():
    # Requested .3 absent everywhere. Flanking .2 and .4 share the same genomic
    # CDS map in GRCh38 → bracket agrees → safe (probabilistic).
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4, "GRCh38")
    assert safe
    assert "bracket" in reason


def test_genomic_bracket_unsafe_when_flanks_disagree():
    # Flanking .2 and .4 sit on different genomic coordinates → bracket disagrees.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS_MOVED}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4, "GRCh38")
    assert not safe
    assert "may have moved" in reason


def test_genomic_bracket_unverified_without_bracket():
    # Requested .9 absent and nothing above it → cannot bracket.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 9, 4, "GRCh38")
    assert not safe
    assert "cannot bracket" in reason


def test_unverified_without_genome_build():
    # No genome_build → the structural check can still run when the record is
    # present, but the bracket fallback can't; an absent requested version is
    # therefore unverified.
    dp = _provider({
        "NM_000059.2": _tx("NM_000059.2", 0, 30, {"GRCh38": EXONS}),
        "NM_000059.4": _tx("NM_000059.4", 0, 30, {"GRCh38": EXONS}),
    })
    safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4)
    assert not safe
    assert "cannot verify" in reason


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
