"""
Tests for cdot.hgvs.gene_hgvs — gene-symbol HGVS resolution via MANE/canonical tags.

Uses the Ensembl GRCh38 test data fixture which has a single gene (AOAH) with
ENST00000617537.5 tagged as MANE_Select.
"""
import os
import pytest

from cdot.hgvs.clean import HGVSFixCode, HGVSFixSeverity, HGVSInputError, VersionStrategy
from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider, LocalDataProvider
from cdot.hgvs.gene_hgvs import (
    DEFAULT_TAG_PRIORITY,
    Consortium,
    _consortium_of,
    _filter_by_consortium,
    _parse_gene_only_hgvs,
    _parse_versioned_transcript,
    _rank_transcripts_by_tags,
    consortium_of,
    fix_hgvs,
    rank_transcripts_for_gene,
    resolve_gene_hgvs,
    resolve_transcript_version,
    UnsafeVersionPolicy,
)

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
ENSEMBL_JSON = os.path.join(TEST_DATA_DIR, "cdot.ensembl.grch38.json")

C = HGVSFixCode
W = HGVSFixSeverity.WARNING
E = HGVSFixSeverity.ERROR


@pytest.fixture(scope="module")
def ensembl_provider():
    return JSONDataProvider([ENSEMBL_JSON])


class _StubProvider:
    """Minimal data provider for consortium-filter tests.

    resolve_gene_hgvs only calls get_tx_ac_tags_for_gene(gene, build), so we
    stub that with a gene→[(tx_ac, tags)] map (case-sensitive, like Redis).
    """
    def __init__(self, gene_map):
        self._gene_map = gene_map

    def get_tx_ac_tags_for_gene(self, gene, genome_build):
        return self._gene_map.get(gene, [])


@pytest.fixture
def refseq_mixed_provider():
    return _StubProvider({
        # BRCA2: MANE Select pair present in both consortiums
        "BRCA2": [
            ("ENST00000380152.8", ["MANE_Select", "Ensembl_canonical"]),
            ("NM_000059.4", ["MANE Select"]),  # RefSeq space-form tag
        ],
        # Gene that only exists in Ensembl
        "ENSEMBLONLY": [
            ("ENST00000999999.1", ["MANE_Select"]),
        ],
    })


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


@pytest.mark.parametrize("hgvs_string", [
    "nm_000059.4:c.36del",
    "enst00000617537.5:c.36del",
])
def test_parse_gene_only_hgvs_lowercase_transcript_is_passthrough(hgvs_string):
    # A lowercase transcript accession must be recognised as a transcript (and
    # passed through), not mistaken for a gene symbol. Otherwise resolve_gene_hgvs
    # tries to look up e.g. "enst00000617537.5" as a gene and errors.
    gene, allele = _parse_gene_only_hgvs(hgvs_string)
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


def test_rank_refseq_space_form_tags_match():
    # RefSeq stores tags with spaces ("MANE Select"); must still match the
    # underscore-form entries in DEFAULT_TAG_PRIORITY.
    tx_and_tags = [
        ("NM_002", []),
        ("NM_001", ["MANE Select"]),
    ]
    ranked = _rank_transcripts_by_tags(tx_and_tags, DEFAULT_TAG_PRIORITY)
    assert ranked[0][0] == "NM_001"
    assert ranked[0][2] == "MANE_Select"          # matched tag is normalised
    assert ranked[0][1] == ["MANE Select"]        # original tag preserved


def test_rank_refseq_select_space_form_matches():
    tx_and_tags = [("NM_001", ["RefSeq Select"])]
    ranked = _rank_transcripts_by_tags(tx_and_tags, DEFAULT_TAG_PRIORITY)
    assert ranked[0][2] == "RefSeq_Select"


# ---------------------------------------------------------------------------
# _filter_by_consortium / _consortium_of
# ---------------------------------------------------------------------------

# The MANE Select pair: same biological transcript in both consortiums.
MANE_PAIR = [
    ("ENST00000380152.8", ["MANE_Select"]),
    ("NM_000059.4", ["MANE Select"]),
]


def test_consortium_of():
    assert _consortium_of("ENST00000380152.8") == Consortium.ENSEMBL
    assert _consortium_of("NM_000059.4") == Consortium.REFSEQ
    assert _consortium_of("XR_001.1") == Consortium.REFSEQ


def test_consortium_of_public_alias():
    # Promoted from _consortium_of (#114); the private name stays as an alias.
    assert consortium_of is _consortium_of
    assert consortium_of("ENST00000380152.8") == Consortium.ENSEMBL
    assert consortium_of("NM_000059.4") == Consortium.REFSEQ


# ---------------------------------------------------------------------------
# rank_transcripts_for_gene — the ranking core split out of resolution (#114)
# ---------------------------------------------------------------------------

def test_rank_transcripts_for_gene_returns_full_ordered_list(refseq_mixed_provider):
    # prefer_consortium=None keeps both consortiums so we get the full ranking
    ranked, fixes = rank_transcripts_for_gene(
        "BRCA2", refseq_mixed_provider, "GRCh38", prefer_consortium=None,
    )
    assert [tx_ac for tx_ac, _tags, _matched in ranked] == [
        "ENST00000380152.8", "NM_000059.4",
    ]
    assert all(matched == "MANE_Select" for _tx, _tags, matched in ranked)
    assert fixes == []


def test_rank_transcripts_for_gene_consortium_hard_filter(refseq_mixed_provider):
    ranked, fixes = rank_transcripts_for_gene(
        "BRCA2", refseq_mixed_provider, "GRCh38", prefer_consortium=Consortium.REFSEQ,
    )
    assert [tx_ac for tx_ac, _tags, _matched in ranked] == ["NM_000059.4"]
    assert fixes == []


def test_rank_transcripts_for_gene_no_transcript_errors(refseq_mixed_provider):
    ranked, fixes = rank_transcripts_for_gene(
        "NOSUCHGENE", refseq_mixed_provider, "GRCh38",
    )
    assert ranked == []
    assert fixes[-1].severity == E
    assert fixes[-1].code == C.NO_TRANSCRIPT_FOR_GENE


def test_rank_transcripts_for_gene_preferred_consortium_absent_errors(refseq_mixed_provider):
    # ENSEMBLONLY has only an Ensembl transcript; prefer=REFSEQ must error, not cross over
    ranked, fixes = rank_transcripts_for_gene(
        "ENSEMBLONLY", refseq_mixed_provider, "GRCh38", prefer_consortium=Consortium.REFSEQ,
    )
    assert ranked == []
    assert fixes[-1].code == C.NO_TRANSCRIPT_FOR_GENE
    assert "refseq" in fixes[-1].message.lower()


def test_filter_by_consortium_refseq():
    filtered = _filter_by_consortium(MANE_PAIR, Consortium.REFSEQ)
    assert [t[0] for t in filtered] == ["NM_000059.4"]


def test_filter_by_consortium_ensembl():
    filtered = _filter_by_consortium(MANE_PAIR, Consortium.ENSEMBL)
    assert [t[0] for t in filtered] == ["ENST00000380152.8"]


# ---------------------------------------------------------------------------
# resolve_gene_hgvs — consortium hard filter (uses RefSeq test fixture)
# ---------------------------------------------------------------------------

def test_resolve_prefers_refseq_from_mixed(refseq_mixed_provider):
    resolved, fixes = resolve_gene_hgvs(
        "BRCA2:c.36del", refseq_mixed_provider, "GRCh38",
        prefer_consortium=Consortium.REFSEQ,
    )
    assert resolved == "NM_000059.4:c.36del"


def test_resolve_prefers_ensembl_from_mixed(refseq_mixed_provider):
    resolved, fixes = resolve_gene_hgvs(
        "BRCA2:c.36del", refseq_mixed_provider, "GRCh38",
        prefer_consortium=Consortium.ENSEMBL,
    )
    assert resolved == "ENST00000380152.8:c.36del"


def test_resolve_refseq_filter_errors_when_only_ensembl(refseq_mixed_provider):
    # ENSG-only gene with prefer=REFSEQ must error, never cross to Ensembl
    resolved, fixes = resolve_gene_hgvs(
        "ENSEMBLONLY:c.36del", refseq_mixed_provider, "GRCh38",
        prefer_consortium=Consortium.REFSEQ,
    )
    assert resolved == "ENSEMBLONLY:c.36del"  # unchanged
    assert fixes[-1].severity == E
    assert fixes[-1].code == C.NO_TRANSCRIPT_FOR_GENE
    assert "refseq" in fixes[-1].message.lower()


def test_resolve_none_allows_either(refseq_mixed_provider):
    resolved, _fixes = resolve_gene_hgvs(
        "ENSEMBLONLY:c.36del", refseq_mixed_provider, "GRCh38",
        prefer_consortium=None,
    )
    assert resolved == "ENST00000999999.1:c.36del"


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
    # Fixture is Ensembl data, so allow Ensembl (default is RefSeq)
    resolved, fixes = resolve_gene_hgvs("AOAH:c.36del", ensembl_provider, "GRCh38",
                                        prefer_consortium=Consortium.ENSEMBL)
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


def test_resolve_lowercase_gene_symbol_retries_uppercased(ensembl_provider):
    # Gene lookups are case-sensitive; a lowercase symbol should retry uppercased
    resolved, fixes = resolve_gene_hgvs("aoah:c.36del", ensembl_provider, "GRCh38",
                                        prefer_consortium=Consortium.ENSEMBL)
    assert resolved == "ENST00000617537.5:c.36del"
    assert any(f.code == C.UPPERCASED_GENE_SYMBOL for f in fixes)
    uc_fix = next(f for f in fixes if f.code == C.UPPERCASED_GENE_SYMBOL)
    assert uc_fix.original == "aoah"
    assert uc_fix.fixed == "AOAH"
    assert any(f.code == C.RESOLVED_GENE_TO_TRANSCRIPT for f in fixes)


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
    result, fixes = fix_hgvs("AOAH:c.36del", data_provider=ensembl_provider,
                             genome_build="GRCh38", prefer_consortium=Consortium.ENSEMBL)
    assert result == "ENST00000617537.5:c.36del"
    assert any(f.code == C.RESOLVED_GENE_TO_TRANSCRIPT for f in fixes)


def test_fix_hgvs_cleaning_then_resolution(ensembl_provider):
    # uppercase mutation type should be cleaned, then gene resolved
    result, fixes = fix_hgvs("AOAH:c.36DEL", data_provider=ensembl_provider,
                             genome_build="GRCh38", prefer_consortium=Consortium.ENSEMBL)
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
        prefer_consortium=Consortium.ENSEMBL,
    )
    assert result == "ENST00000617537.5:c.36del"
    assert all(f.severity == W for f in fixes)


def test_fix_hgvs_no_provider_no_genome_build():
    # Omitting data_provider entirely — only cleaning
    result, fixes = fix_hgvs("nm_000059.4:c.316+5G>A")
    assert result == "NM_000059.4:c.316+5G>A"
    assert any(f.code == C.UPPERCASED_TRANSCRIPT for f in fixes)


# ---------------------------------------------------------------------------
# _parse_versioned_transcript
# ---------------------------------------------------------------------------

def test_parse_versioned_transcript_plain():
    assert _parse_versioned_transcript("NM_000059.4:c.36del") == ("NM_000059", 4)


def test_parse_versioned_transcript_with_gene_parens():
    assert _parse_versioned_transcript("NM_000059.4(BRCA2):c.36del") == ("NM_000059", 4)


def test_parse_versioned_transcript_ensembl():
    assert _parse_versioned_transcript("ENST00000380152.8:c.36del") == ("ENST00000380152", 8)


def test_parse_versioned_transcript_no_version():
    assert _parse_versioned_transcript("NM_000059:c.36del") is None


def test_parse_versioned_transcript_gene_only():
    assert _parse_versioned_transcript("BRCA2:c.36del") is None


def test_parse_versioned_transcript_no_colon():
    assert _parse_versioned_transcript("NM_000059.4") is None


# ---------------------------------------------------------------------------
# resolve_transcript_version — adjacent-version fallback
# ---------------------------------------------------------------------------

class _VersionStubProvider:
    """Stub exposing only get_tx_versions(accession) -> list[int]."""
    def __init__(self, version_map):
        self._version_map = version_map

    def get_tx_versions(self, accession):
        return self._version_map.get(accession, [])


@pytest.fixture
def version_provider():
    return _VersionStubProvider({"NM_000059": [2, 3, 4]})


def test_resolve_version_exact_match_unchanged(version_provider):
    resolved, fixes = resolve_transcript_version("NM_000059.4:c.36del", version_provider)
    assert resolved == "NM_000059.4:c.36del"
    assert fixes == []


def test_resolve_version_falls_back_up(version_provider):
    resolved, fixes = resolve_transcript_version("NM_000059.1:c.36del", version_provider)
    assert resolved == "NM_000059.2:c.36del"
    assert len(fixes) == 1
    assert fixes[0].severity == W
    assert fixes[0].code == C.USED_ADJACENT_VERSION


def test_resolve_version_falls_back_down_when_no_higher(version_provider):
    # requested .9 — nothing higher, UP_THEN_DOWN drops to the highest available
    resolved, fixes = resolve_transcript_version("NM_000059.9:c.36del", version_provider)
    assert resolved == "NM_000059.4:c.36del"
    assert fixes[0].code == C.USED_ADJACENT_VERSION


def test_resolve_version_latest_strategy(version_provider):
    resolved, _fixes = resolve_transcript_version(
        "NM_000059.1:c.36del", version_provider, strategy=VersionStrategy.LATEST
    )
    assert resolved == "NM_000059.4:c.36del"


def test_resolve_version_preserves_gene_parens(version_provider):
    resolved, fixes = resolve_transcript_version("NM_000059.1(BRCA2):c.36del", version_provider)
    assert resolved == "NM_000059.2(BRCA2):c.36del"
    assert fixes[0].code == C.USED_ADJACENT_VERSION


def test_resolve_version_unknown_accession_errors(version_provider):
    resolved, fixes = resolve_transcript_version("NM_999999.1:c.36del", version_provider)
    assert resolved == "NM_999999.1:c.36del"  # unchanged
    assert fixes[0].severity == E
    assert fixes[0].code == C.NO_TRANSCRIPT_VERSIONS


def test_resolve_version_no_versioned_transcript_noop(version_provider):
    # gene-only / versionless input — nothing to do
    assert resolve_transcript_version("BRCA2:c.36del", version_provider) == ("BRCA2:c.36del", [])
    assert resolve_transcript_version("NM_000059:c.36del", version_provider) == ("NM_000059:c.36del", [])


def test_resolve_version_unsupported_provider_errors():
    class _NoVersions:
        def get_tx_versions(self, accession):
            raise NotImplementedError("nope")

    resolved, fixes = resolve_transcript_version("NM_000059.1:c.36del", _NoVersions())
    assert resolved == "NM_000059.1:c.36del"
    assert fixes[0].severity == E
    assert fixes[0].code == C.NO_TRANSCRIPT_VERSIONS


# ---------------------------------------------------------------------------
# resolve_transcript_version — coordinate-safety check (#28)
#
# A real JSONDataProvider (which implements is_version_substitution_safe) backs
# these, so the substitution is checked structurally before it is applied. The
# requested version sits only in GRCh37, the substitution candidates in GRCh38, so
# the build-independent structure decides the verdict (matching the cross-build
# result in docs/transcript_version_safety.md).
# ---------------------------------------------------------------------------

def _safety_build(exons, contig="NC_000013.11"):
    return {"contig": contig, "strand": "+", "url": None,
            "cds_start": None, "cds_end": None, "exons": exons}


def _safety_tx(acc, stop_codon, exons, build):
    return {"id": acc, "gene_name": "BRCA2", "start_codon": 0, "stop_codon": stop_codon,
            "genome_builds": {build: _safety_build(exons)}}


# Same transcript layout (CDS 1..30, split at 11) on two genomic placements; and a
# CDS-length-changed layout that is NOT coordinate-safe.
_SAFE_EXONS  = [[100, 110, 0, 1, 10, None], [200, 220, 1, 11, 30, None]]
_SAFE_MOVED  = [[900, 910, 0, 1, 10, None], [2000, 2020, 1, 11, 30, None]]
_UNSAFE_EXONS = [[100, 110, 0, 1, 10, None], [200, 223, 1, 11, 33, None]]  # cds_len 33


def _safety_provider(requested_exons, requested_stop):
    """Provider with requested .2 in GRCh37 only, and .3/.4 in GRCh38."""
    import io
    import json
    data = {"cdot_version": "0.2.33", "genome_builds": ["GRCh37", "GRCh38"],
            "transcripts": {
                "NM_000059.2": _safety_tx("NM_000059.2", requested_stop, requested_exons, "GRCh37"),
                "NM_000059.3": _safety_tx("NM_000059.3", 30, _SAFE_MOVED, "GRCh38"),
                "NM_000059.4": _safety_tx("NM_000059.4", 30, _SAFE_MOVED, "GRCh38"),
            }}
    return JSONDataProvider([io.StringIO(json.dumps(data))])


def test_resolve_version_coord_safe_substitutes():
    dp = _safety_provider(_SAFE_EXONS, 30)   # requested .2 same structure as .3
    resolved, fixes = resolve_transcript_version(
        "NM_000059.2:c.36del", dp, genome_build="GRCh38")
    assert resolved == "NM_000059.3:c.36del"
    assert fixes[0].severity == W
    assert fixes[0].code == C.USED_ADJACENT_VERSION_COORD_SAFE


def test_resolve_version_unsafe_refuses_by_default():
    dp = _safety_provider(_UNSAFE_EXONS, 33)  # requested .2 differs in CDS length
    resolved, fixes = resolve_transcript_version(
        "NM_000059.2:c.36del", dp, genome_build="GRCh38")
    assert resolved == "NM_000059.2:c.36del"   # unchanged — refused
    assert fixes[0].severity == E
    assert fixes[0].code == C.REFUSED_UNSAFE_VERSION


def test_resolve_version_unsafe_substitute_policy_warns():
    dp = _safety_provider(_UNSAFE_EXONS, 33)
    resolved, fixes = resolve_transcript_version(
        "NM_000059.2:c.36del", dp, genome_build="GRCh38",
        on_unsafe_version=UnsafeVersionPolicy.SUBSTITUTE)
    assert resolved == "NM_000059.3:c.36del"   # substituted under SUBSTITUTE policy
    assert fixes[0].severity == W
    assert fixes[0].code == C.USED_ADJACENT_VERSION_COORD_UNVERIFIED


def test_fix_hgvs_version_fallback_unsafe_refuses_by_default():
    dp = _safety_provider(_UNSAFE_EXONS, 33)
    result, fixes = fix_hgvs(
        "NM_000059.2:c.36del", data_provider=dp, genome_build="GRCh38",
        version_fallback=VersionStrategy.UP_THEN_DOWN)
    assert result == "NM_000059.2:c.36del"
    assert any(f.code == C.REFUSED_UNSAFE_VERSION for f in fixes)


def test_fix_hgvs_version_fallback_unsafe_substitute_opt_in():
    dp = _safety_provider(_UNSAFE_EXONS, 33)
    result, fixes = fix_hgvs(
        "NM_000059.2:c.36del", data_provider=dp, genome_build="GRCh38",
        version_fallback=VersionStrategy.UP_THEN_DOWN,
        on_unsafe_version=UnsafeVersionPolicy.SUBSTITUTE)
    assert result == "NM_000059.3:c.36del"
    assert any(f.code == C.USED_ADJACENT_VERSION_COORD_UNVERIFIED for f in fixes)


# ---------------------------------------------------------------------------
# fix_hgvs — opt-in version_fallback
# ---------------------------------------------------------------------------

def test_fix_hgvs_version_fallback_off_by_default(version_provider):
    # No version_fallback → version untouched even though .1 isn't available
    result, fixes = fix_hgvs("NM_000059.1:c.36del", data_provider=version_provider)
    assert result == "NM_000059.1:c.36del"
    assert not any(f.code == C.USED_ADJACENT_VERSION for f in fixes)


def test_fix_hgvs_version_fallback_opt_in(version_provider):
    result, fixes = fix_hgvs(
        "NM_000059.1:c.36del", data_provider=version_provider,
        version_fallback=VersionStrategy.UP_THEN_DOWN,
    )
    assert result == "NM_000059.2:c.36del"
    assert any(f.code == C.USED_ADJACENT_VERSION for f in fixes)


def test_fix_hgvs_version_fallback_with_cleaning(version_provider):
    # messy input is cleaned first, then the version is bumped
    result, fixes = fix_hgvs(
        "nm_000059.1:c.36DEL", data_provider=version_provider,
        version_fallback=VersionStrategy.UP_THEN_DOWN,
    )
    assert result == "NM_000059.2:c.36del"
    codes = {f.code for f in fixes}
    assert C.LOWERCASED_MUTATION_TYPE in codes
    assert C.USED_ADJACENT_VERSION in codes


def test_fix_hgvs_version_fallback_raise_on_errors(version_provider):
    with pytest.raises(HGVSInputError):
        fix_hgvs(
            "NM_999999.1:c.36del", data_provider=version_provider,
            version_fallback=VersionStrategy.UP_THEN_DOWN, raise_on_errors=True,
        )


# ---------------------------------------------------------------------------
# JSONDataProvider.get_tx_ac_tags_for_gene — ranks by spliced (exonic) length
# ---------------------------------------------------------------------------

def test_json_provider_tags_ranked_by_exonic_not_genomic_length():
    """Ranking must use spliced (exonic) length, not genomic span.

    WIDE_SPAN has a huge genomic footprint but little exonic sequence; COMPACT
    has a small footprint but more exonic sequence. The longer *transcript*
    (COMPACT) must rank first - ranking on genomic span would wrongly pick
    WIDE_SPAN.
    """
    import io
    import json

    def _build(exons):
        return {
            "contig": "NC_000001.11", "strand": "+", "url": None,
            "cds_start": None, "cds_end": None, "exons": exons,
        }

    data = {
        "cdot_version": "0.2.26",
        "genome_builds": ["GRCh38"],
        "transcripts": {
            # 2 tiny exons spread across a 10kb span -> exonic length 100
            "WIDE_SPAN.1": {
                "id": "WIDE_SPAN.1", "gene_name": "GENEX",
                "genome_builds": {"GRCh38": _build(
                    [[0, 50, 1, 1, 50, None], [9950, 10000, 2, 51, 100, None]])},
            },
            # 1 big exon in a 500bp span -> exonic length 500
            "COMPACT.1": {
                "id": "COMPACT.1", "gene_name": "GENEX",
                "genome_builds": {"GRCh38": _build([[0, 500, 1, 1, 500, None]])},
            },
        },
    }
    dp = JSONDataProvider([io.StringIO(json.dumps(data))])
    ranked = dp.get_tx_ac_tags_for_gene("GENEX", "GRCh38")
    assert [tx_ac for tx_ac, _tags in ranked] == ["COMPACT.1", "WIDE_SPAN.1"]


# ---------------------------------------------------------------------------
# LocalDataProvider.get_tx_ac_tags_for_gene — returns the canonical versioned
# accession (transcript_data["id"]) even when _get_transcript_ids_for_gene
# yields a versionless id (#114)
# ---------------------------------------------------------------------------

class _VersionlessLocalProvider(LocalDataProvider):
    """Minimal LocalDataProvider whose gene→tx map yields *versionless* ids
    (like a store keyed without version) while each record carries the canonical
    versioned accession in its ``"id"`` field.

    Exercises the real get_tx_ac_tags_for_gene / _get_tags_by_tx_ac /
    _get_transcript_tags path. _get_transcript resolves a lookup by either the
    versionless id (from _get_transcript_ids_for_gene) or the versioned id (what
    get_tx_ac_tags_for_gene now ranks by).
    """
    def __init__(self, gene_to_tx_ids, records):
        # Deliberately skip AbstractJSONDataProvider.__init__ (assembly loading);
        # the gene→tags path under test doesn't need it.
        self._gene_to_tx_ids = gene_to_tx_ids
        self._records = records

    def _get_transcript(self, tx_ac):
        if tx_ac in self._records:
            return self._records[tx_ac]
        for rec in self._records.values():
            if rec.get("id") == tx_ac:
                return rec
        return None

    def _get_transcript_ids_for_gene(self, gene):
        return self._gene_to_tx_ids.get(gene, [])

    def _get_gene(self, gene):
        raise NotImplementedError

    def _get_contig_interval_tree(self, alt_ac):
        raise NotImplementedError

    def get_tx_for_gene(self, gene):
        raise NotImplementedError

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        raise NotImplementedError


def _versionless_record(tx_id, gene, tag):
    build = {
        "contig": "NC_000013.11", "strand": "+", "url": None,
        "cds_start": None, "cds_end": None,
        "exons": [[0, 100, 1, None, None, None]],
    }
    if tag is not None:
        build["tag"] = tag
    rec = {"gene_name": gene, "genome_builds": {"GRCh38": build}}
    if tx_id is not None:
        rec["id"] = tx_id
    return rec


def test_get_tx_ac_tags_for_gene_returns_versioned_id_from_record():
    # gene map yields versionless "NM_000059", record carries id "NM_000059.4"
    dp = _VersionlessLocalProvider(
        {"BRCA2": ["NM_000059"]},
        {"NM_000059": _versionless_record("NM_000059.4", "BRCA2", "MANE Select")},
    )
    result = dp.get_tx_ac_tags_for_gene("BRCA2", "GRCh38")
    # Returns the canonical versioned accession, with tags resolved against it
    assert result == [("NM_000059.4", ["MANE Select"])]


def test_rank_transcripts_for_gene_versionless_provider_yields_versioned():
    dp = _VersionlessLocalProvider(
        {"BRCA2": ["NM_000059"]},
        {"NM_000059": _versionless_record("NM_000059.4", "BRCA2", "MANE Select")},
    )
    ranked, fixes = rank_transcripts_for_gene("BRCA2", dp, "GRCh38")
    assert fixes == []
    assert len(ranked) == 1
    tx_ac, tags, matched = ranked[0]
    assert tx_ac == "NM_000059.4"            # versioned, not "NM_000059"
    assert tags == ["MANE Select"]           # tags intact
    assert matched == "MANE_Select"          # matched_tag intact


def test_resolve_gene_hgvs_versionless_provider_rewrites_to_versioned():
    dp = _VersionlessLocalProvider(
        {"BRCA2": ["NM_000059"]},
        {"NM_000059": _versionless_record("NM_000059.4", "BRCA2", "MANE Select")},
    )
    resolved, fixes = resolve_gene_hgvs("BRCA2:c.36del", dp, "GRCh38")
    assert resolved == "NM_000059.4:c.36del"
    assert any(f.code == C.RESOLVED_GENE_TO_TRANSCRIPT for f in fixes)
    resolve_fix = next(f for f in fixes if f.code == C.RESOLVED_GENE_TO_TRANSCRIPT)
    assert resolve_fix.fixed == "NM_000059.4"


def test_get_tx_ac_tags_for_gene_falls_back_to_id_in_hand_when_no_id():
    # Record has no "id" → fall back to the id from _get_transcript_ids_for_gene
    dp = _VersionlessLocalProvider(
        {"BRCA2": ["NM_FALLBACK"]},
        {"NM_FALLBACK": _versionless_record(None, "BRCA2", "MANE Select")},
    )
    result = dp.get_tx_ac_tags_for_gene("BRCA2", "GRCh38")
    assert result == [("NM_FALLBACK", ["MANE Select"])]


def test_get_tx_ac_tags_for_gene_skips_transcript_with_no_data_for_build():
    # A transcript with no genome_builds entry for the build is skipped, while
    # the one with data is still returned by its versioned accession.
    dp = _VersionlessLocalProvider(
        {"BRCA2": ["NM_000059", "NM_OTHER"]},
        {
            "NM_000059": _versionless_record("NM_000059.4", "BRCA2", "MANE Select"),
            # No GRCh38 build for this one
            "NM_OTHER": {"id": "NM_OTHER.1", "gene_name": "BRCA2",
                         "genome_builds": {"GRCh37": {"exons": [[0, 100, 1, None, None, None]]}}},
        },
    )
    result = dp.get_tx_ac_tags_for_gene("BRCA2", "GRCh38")
    assert result == [("NM_000059.4", ["MANE Select"])]


# ---------------------------------------------------------------------------
# JSONDataProvider.get_tx_versions — enumerate versions from in-memory data
# ---------------------------------------------------------------------------

def test_json_provider_get_tx_versions(ensembl_provider):
    # Fixture has ENST00000617537.5
    assert ensembl_provider.get_tx_versions("ENST00000617537") == [5]


def test_json_provider_get_tx_versions_unknown(ensembl_provider):
    assert ensembl_provider.get_tx_versions("NM_999999") == []


def test_json_provider_get_tx_versions_no_prefix_overmatch(ensembl_provider):
    # "ENST0000061753" must not match "ENST00000617537.5" (dot guard)
    assert ensembl_provider.get_tx_versions("ENST0000061753") == []


# ---------------------------------------------------------------------------
# RESTDataProvider.get_tx_versions — drives the versionless /transcript endpoint
# ---------------------------------------------------------------------------

def test_rest_provider_get_tx_versions(monkeypatch):
    from cdot.hgvs.dataproviders.json_data_provider import RESTDataProvider

    dp = RESTDataProvider(url="https://example.org")
    captured = {}

    def fake_get(url):
        captured["url"] = url
        # versionless lookup returns every version keyed by full accession
        return {"NM_000059.3": {"id": "NM_000059.3"}, "NM_000059.4": {"id": "NM_000059.4"}}

    monkeypatch.setattr(dp, "_get_from_url", fake_get)

    assert dp.get_tx_versions("NM_000059") == [3, 4]
    assert captured["url"] == "https://example.org/transcript/NM_000059"
    # versioned transcripts are warmed into the cache as a bonus
    assert dp.transcripts["NM_000059.4"] == {"id": "NM_000059.4"}


def test_rest_provider_get_tx_versions_unknown(monkeypatch):
    from cdot.hgvs.dataproviders.json_data_provider import RESTDataProvider

    dp = RESTDataProvider(url="https://example.org")
    monkeypatch.setattr(dp, "_get_from_url", lambda url: None)  # 404 → None
    assert dp.get_tx_versions("NM_999999") == []


# ---------------------------------------------------------------------------
# RESTDataProvider.get_tx_ac_tags_for_gene — drives the cdot_rest tags endpoint
# ---------------------------------------------------------------------------

def test_rest_provider_get_tx_ac_tags_for_gene(monkeypatch):
    from cdot.hgvs.dataproviders.json_data_provider import RESTDataProvider

    dp = RESTDataProvider(url="https://example.org")

    captured = {}

    def fake_get(url):
        captured["url"] = url
        # Server JSON: tuples come back as lists nested under "results"
        return {"results": [["NM_000059.4", ["MANE_Select", "basic"]],
                            ["NM_000059.3", []]]}

    monkeypatch.setattr(dp, "_get_from_url", fake_get)

    result = dp.get_tx_ac_tags_for_gene("BRCA2", "GRCh38")

    assert captured["url"] == "https://example.org/transcripts/gene/BRCA2/tags/GRCh38"
    assert result == [("NM_000059.4", ["MANE_Select", "basic"]),
                      ("NM_000059.3", [])]


def test_rest_provider_get_tx_ac_tags_for_gene_empty(monkeypatch):
    from cdot.hgvs.dataproviders.json_data_provider import RESTDataProvider

    dp = RESTDataProvider(url="https://example.org")
    monkeypatch.setattr(dp, "_get_from_url", lambda url: None)  # 404 → None
    assert dp.get_tx_ac_tags_for_gene("NOPE", "GRCh38") == []


# ---------------------------------------------------------------------------
# LocalDataProvider tag hooks: _get_transcript_tags signature + batch hook (#114)
# ---------------------------------------------------------------------------

def test_get_transcript_tags_signature_takes_tx_ac(ensembl_provider):
    # #114 item 4: _get_transcript_tags now receives tx_ac explicitly so overrides
    # don't have to re-derive it from transcript_data.
    tx_ac = "ENST00000617537.5"
    transcript_data = ensembl_provider._get_transcript(tx_ac)
    tags = ensembl_provider._get_transcript_tags(tx_ac, transcript_data, "GRCh38")
    assert "MANE_Select" in tags


def test_get_tx_ac_tags_for_gene_uses_batch_hook(ensembl_provider, monkeypatch):
    # #114 item 5: get_tx_ac_tags_for_gene must source tags via the overridable
    # batch hook _get_tags_by_tx_ac, so a subclass can answer in one query.
    calls = {}

    def fake_batch(tx_acs, genome_build):
        calls["tx_acs"] = list(tx_acs)
        calls["genome_build"] = genome_build
        return {tx_ac: ["MANE_Select"] for tx_ac in tx_acs}

    monkeypatch.setattr(ensembl_provider, "_get_tags_by_tx_ac", fake_batch)
    result = ensembl_provider.get_tx_ac_tags_for_gene("AOAH", "GRCh38")

    assert calls["genome_build"] == "GRCh38"
    assert calls["tx_acs"]  # got the gene's accessions in one call
    assert all(tags == ["MANE_Select"] for _tx_ac, tags in result)


def test_get_tags_by_tx_ac_default_loops_per_transcript_hook(ensembl_provider):
    # Default batch implementation defers to the per-transcript hook.
    tx_acs = ensembl_provider._get_transcript_ids_for_gene("AOAH")
    by_ac = ensembl_provider._get_tags_by_tx_ac(list(tx_acs), "GRCh38")
    assert "MANE_Select" in by_ac["ENST00000617537.5"]
