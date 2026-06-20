"""
Regression tests for ``is_version_substitution_safe`` against *real* transcript
records, exercising the alignment-gap and UTR-length coordinate-safety gates.

The records in ``tests/test_data/version_safety/clinvar_real_records.json`` are the
GRCh38 placements of a handful of RefSeq transcripts extracted from a cdot release.
Each case below corresponds to a real ClinVar variant whose genomic coordinate is
*not* preserved across the two versions, even though the versions share an identical
intrinsic CDS structure - exactly the class of substitution the structural check
used to wrongly call safe:

  * NM_000314 .8 vs .7  (PTEN, + strand): .7 carries a 1-base alignment gap in the
    first coding exon that .8 does not, so every coding base downstream shifts.
  * NM_000277 .1 vs .2  (GALT, - strand): same alignment-gap class.
  * NM_005577 .4 vs .2  (LPA, - strand): identical CDS map and no gaps, but a
    different 5'UTR length, so a 5'UTR (c.-N) variant moves while coding and
    CDS-intronic variants on the same pair do not.

ClinVar is public data, so these accessions/positions may appear here (unlike the
private search-log corpus).
"""
from pathlib import Path

import pytest

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider

RECORDS = Path(__file__).parent / "test_data" / "version_safety" / "clinvar_real_records.json"
BUILD = "GRCh38"


@pytest.fixture(scope="module")
def dp():
    return JSONDataProvider([str(RECORDS)])


# (accession, requested_version, substitute_version, cited_position, expected_safe, reason_substr)
CASES = [
    # Alignment-gap mode: the whole pair is unsafe (a coding-exon gap differs), so
    # even a plain coding position is correctly refused.
    ("NM_000314", 8, 7, "386G>A", False, "gap"),
    ("NM_000277", 1, 2, "100C>T", False, "gap"),

    # UTR mode: identical CDS structure, no gaps, but a different 5'UTR length.
    ("NM_005577", 4, 2, "-49=", False, "5'UTR"),           # 5'UTR variant moves
    ("NM_005577", 4, 2, "109C>T", True, "identical"),       # coding stays safe
    ("NM_005577", 4, 2, "2604-1G>C", True, "identical"),    # CDS-intronic stays safe

    # A common bump that changes BOTH UTRs but keeps the CDS: coding safe, UTR not.
    ("NM_000059", 3, 4, "36del", True, "identical"),
    ("NM_000059", 3, 4, "-15A>G", False, "5'UTR"),
    ("NM_000059", 3, 4, "*20del", False, "3'UTR"),
]


@pytest.mark.parametrize("acc, req, sub, pos, expected_safe, reason_substr", CASES)
def test_substitution_safety_real(dp, acc, req, sub, pos, expected_safe, reason_substr):
    safe, reason = dp.is_version_substitution_safe(acc, req, sub, BUILD, cited_position=pos)
    assert safe is expected_safe, f"{acc} .{req}->.{sub} c.{pos}: {reason}"
    assert reason_substr in reason


def test_coding_default_position_ignores_utr_difference(dp):
    # Without a cited position (the backward-compatible default) the check gives the
    # coding guarantee: a UTR-only difference does not make a pair unsafe.
    safe, reason = dp.is_version_substitution_safe("NM_000059", 3, 4, BUILD)
    assert safe
    assert "identical CDS structure" in reason


def test_genomic_replacement_is_refused(dp):
    # NM_018263 .6 vs .4 share an identical intrinsic CDS structure with no alignment
    # gaps, but .6 aligns the CDS exon at CDS-offset 140 ~9.75 kb away from where .4/.5
    # place it, so a coding/CDS-intronic variant in that exon (eg c.141-2A>G) moves.
    # The build-independent structure is blind to placement, but with both versions
    # present the check recomputes the genomic CDS map from the loaded annotation and
    # refuses the substitution.
    safe, reason = dp.is_version_substitution_safe("NM_018263", 6, 4, BUILD,
                                                   cited_position="141-2A>G")
    assert not safe, reason
    assert "placement" in reason


def test_genomic_replacement_refused_even_for_plain_coding(dp):
    # The re-placement moves the whole CDS, so it is refused regardless of position
    # (no cited_position needed): a plain coding variant would move too.
    safe, reason = dp.is_version_substitution_safe("NM_018263", 6, 4, BUILD)
    assert not safe
    assert "placement" in reason
