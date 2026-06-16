"""
Tests for cdot.hgvs.gene_hgvs — gene-symbol HGVS resolution via MANE/canonical tags.

Uses the Ensembl GRCh38 test data fixture which has a single gene (AOAH) with
ENST00000617537.5 tagged as MANE_Select.
"""
import os
import pytest

from cdot.hgvs.clean import HGVSFixCode, HGVSFixSeverity, HGVSInputError, VersionStrategy
from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider
from cdot.hgvs.gene_hgvs import (
    DEFAULT_TAG_PRIORITY,
    Consortium,
    _consortium_of,
    _filter_by_consortium,
    _parse_gene_only_hgvs,
    _parse_versioned_transcript,
    _rank_transcripts_by_tags,
    fix_hgvs,
    resolve_gene_hgvs,
    resolve_transcript_version,
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
