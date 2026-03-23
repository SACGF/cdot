"""
Tests for cdot.hgvs.gene_hgvs — gene-symbol HGVS resolution via MANE/canonical tags.

Uses the Ensembl GRCh38 test data fixture which has a single gene (AOAH) with
ENST00000617537.5 tagged as MANE_Select.
"""
import os
import pytest

from cdot.hgvs.clean import HGVSFixCode, HGVSFixSeverity, HGVSInputError
from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider
from cdot.hgvs.gene_hgvs import (
    DEFAULT_TAG_PRIORITY,
    _parse_gene_only_hgvs,
    _rank_transcripts_by_tags,
    fix_hgvs,
    resolve_gene_hgvs,
)

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
ENSEMBL_JSON = os.path.join(TEST_DATA_DIR, "cdot.ensembl.grch38.json")

C = HGVSFixCode
W = HGVSFixSeverity.WARNING
E = HGVSFixSeverity.ERROR


@pytest.fixture(scope="module")
def ensembl_provider():
    return JSONDataProvider([ENSEMBL_JSON])


# ---------------------------------------------------------------------------
# _parse_gene_only_hgvs
# ---------------------------------------------------------------------------

def test_parse_gene_only_hgvs_returns_symbol_and_allele():
    gene, allele = _parse_gene_only_hgvs("BRCA2:c.36del")
    assert gene == "BRCA2"
    assert allele == "c.36del"


def test_parse_gene_only_hgvs_already_has_transcript():
    gene, allele = _parse_gene_only_hgvs("NM_000059.4:c.36del")
    assert gene is None
    assert allele is None


def test_parse_gene_only_hgvs_transcript_with_gene_parens():
    gene, allele = _parse_gene_only_hgvs("NM_000059.4(BRCA2):c.36del")
    assert gene is None
    assert allele is None


def test_parse_gene_only_hgvs_no_colon():
    gene, allele = _parse_gene_only_hgvs("c.36del")
    assert gene is None
    assert allele is None


def test_parse_gene_only_hgvs_ensembl_transcript():
    gene, allele = _parse_gene_only_hgvs("ENST00000617537.5:c.36del")
    assert gene is None
    assert allele is None


# ---------------------------------------------------------------------------
# _rank_transcripts_by_tags
# ---------------------------------------------------------------------------

def test_rank_prefers_mane_select():
    tx_and_tags = [
        ("NM_002", ["Ensembl_canonical"]),
        ("NM_001", ["MANE_Select"]),
        ("NM_003", []),
    ]
    ranked = _rank_transcripts_by_tags(tx_and_tags, DEFAULT_TAG_PRIORITY)
    assert ranked[0][0] == "NM_001"
    assert ranked[0][2] == "MANE_Select"


def test_rank_mane_select_beats_mane_plus_clinical():
    tx_and_tags = [
        ("NM_002", ["MANE_Plus_Clinical"]),
        ("NM_001", ["MANE_Select"]),
    ]
    ranked = _rank_transcripts_by_tags(tx_and_tags, DEFAULT_TAG_PRIORITY)
    assert ranked[0][0] == "NM_001"


def test_rank_no_matching_tag_last():
    tx_and_tags = [
        ("NM_001", ["MANE_Select"]),
        ("NM_002", []),
    ]
    ranked = _rank_transcripts_by_tags(tx_and_tags, DEFAULT_TAG_PRIORITY)
    assert ranked[-1][0] == "NM_002"
    assert ranked[-1][2] is None


def test_rank_matched_tag_field():
    tx_and_tags = [("NM_001", ["Ensembl_canonical", "MANE_Select"])]
    ranked = _rank_transcripts_by_tags(tx_and_tags, DEFAULT_TAG_PRIORITY)
    assert ranked[0][2] == "MANE_Select"


# ---------------------------------------------------------------------------
# resolve_gene_hgvs — already has transcript (passthrough)
# ---------------------------------------------------------------------------

def test_resolve_passthrough_when_transcript_present(ensembl_provider):
    resolved, fixes = resolve_gene_hgvs(
        "ENST00000617537.5:c.36del", ensembl_provider, "GRCh38"
    )
    assert resolved == "ENST00000617537.5:c.36del"
    assert fixes == []


# ---------------------------------------------------------------------------
# resolve_gene_hgvs — gene symbol resolution
# ---------------------------------------------------------------------------

def test_resolve_gene_symbol_picks_mane_select(ensembl_provider):
    resolved, fixes = resolve_gene_hgvs("AOAH:c.36del", ensembl_provider, "GRCh38")
    assert resolved == "ENST00000617537.5:c.36del"
    assert len(fixes) == 1
    fix = fixes[0]
    assert fix.severity == W
    assert fix.code == C.RESOLVED_GENE_TO_TRANSCRIPT
    assert "AOAH" in fix.message
    assert "ENST00000617537.5" in fix.message
    assert "MANE_Select" in fix.message
    assert fix.original == "AOAH"
    assert fix.fixed == "ENST00000617537.5"


def test_resolve_unknown_gene_returns_error(ensembl_provider):
    resolved, fixes = resolve_gene_hgvs(
        "UNKNOWNGENE:c.36del", ensembl_provider, "GRCh38"
    )
    assert resolved == "UNKNOWNGENE:c.36del"
    assert len(fixes) == 1
    assert fixes[0].severity == E
    assert fixes[0].code == C.NO_TRANSCRIPT_FOR_GENE


def test_resolve_wrong_genome_build_returns_error(ensembl_provider):
    resolved, fixes = resolve_gene_hgvs(
        "AOAH:c.36del", ensembl_provider, "GRCh37"
    )
    assert resolved == "AOAH:c.36del"
    assert any(f.code == C.NO_TRANSCRIPT_FOR_GENE for f in fixes)


# ---------------------------------------------------------------------------
# fix_hgvs — cleaning + resolution combined
# ---------------------------------------------------------------------------

def test_fix_hgvs_cleaning_only_no_provider():
    result, fixes = fix_hgvs("NM_000059.4 c.316+5G>A")
    assert result == "NM_000059.4:c.316+5G>A"
    assert any(f.code == C.STRIPPED_WHITESPACE for f in fixes)


def test_fix_hgvs_gene_resolution_with_provider(ensembl_provider):
    result, fixes = fix_hgvs("AOAH:c.36del", data_provider=ensembl_provider, genome_build="GRCh38")
    assert result == "ENST00000617537.5:c.36del"
    assert any(f.code == C.RESOLVED_GENE_TO_TRANSCRIPT for f in fixes)


def test_fix_hgvs_cleaning_then_resolution(ensembl_provider):
    # uppercase mutation type should be cleaned, then gene resolved
    result, fixes = fix_hgvs("AOAH:c.36DEL", data_provider=ensembl_provider, genome_build="GRCh38")
    assert result == "ENST00000617537.5:c.36del"
    fix_codes = {f.code for f in fixes}
    assert C.LOWERCASED_MUTATION_TYPE in fix_codes
    assert C.RESOLVED_GENE_TO_TRANSCRIPT in fix_codes


def test_fix_hgvs_raise_on_errors_for_unknown_gene(ensembl_provider):
    with pytest.raises(HGVSInputError):
        fix_hgvs(
            "UNKNOWNGENE:c.36del",
            data_provider=ensembl_provider,
            genome_build="GRCh38",
            raise_on_errors=True,
        )


def test_fix_hgvs_no_raise_when_no_errors(ensembl_provider):
    result, fixes = fix_hgvs(
        "AOAH:c.36del",
        data_provider=ensembl_provider,
        genome_build="GRCh38",
        raise_on_errors=True,
    )
    assert result == "ENST00000617537.5:c.36del"
    assert all(f.severity == W for f in fixes)


def test_fix_hgvs_no_provider_no_genome_build():
    # Omitting data_provider entirely — only cleaning
    result, fixes = fix_hgvs("nm_000059.4:c.316+5G>A")
    assert result == "NM_000059.4:c.316+5G>A"
    assert any(f.code == C.UPPERCASED_TRANSCRIPT for f in fixes)
