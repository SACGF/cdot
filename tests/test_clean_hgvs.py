"""
Tests for cdot.hgvs.clean — pure string operations, no data provider required.

Test cases ported from vg_code/test_hgvs.py (test_clean_hgvs and
test_fix_gene_transcript) and extended for the version-fallback helpers.
"""
import pytest

from cdot.hgvs.clean import (
    ALL_CLEAN_OPS,
    HGVSCleanOp,
    HGVSFix,
    HGVSFixCode,
    HGVSFixSeverity,
    HGVSInputError,
    VersionStrategy,
    _rank_versions,
    clean_hgvs,
    error_messages,
    get_best_transcript_version,
    messages,
    rank_transcript_versions,
    warning_messages,
)

C = HGVSFixCode
W = HGVSFixSeverity.WARNING
E = HGVSFixSeverity.ERROR


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def codes(fixes: list[HGVSFix]) -> set[HGVSFixCode]:
    return {f.code for f in fixes}


def severities(fixes: list[HGVSFix]) -> set[HGVSFixSeverity]:
    return {f.severity for f in fixes}


# ---------------------------------------------------------------------------
# clean_hgvs — op selection (subset of cleaning operations)
# ---------------------------------------------------------------------------

def test_ops_default_runs_everything():
    # whitespace + lowercase transcript both fixed by default
    result, fixes = clean_hgvs("nm_000059.4 c.316+5G>A")
    assert result == "NM_000059.4:c.316+5G>A"
    assert C.STRIPPED_WHITESPACE in codes(fixes)


def test_ops_allowlist_runs_only_selected():
    # Only strip whitespace; do NOT uppercase/reconstruct the transcript
    result, fixes = clean_hgvs(
        "nm_000059.4 c.316+5G>A",
        ops={HGVSCleanOp.STRIP_WHITESPACE},
        validate=False,
    )
    assert result == "nm_000059.4c.316+5G>A"  # whitespace gone, casing untouched
    assert codes(fixes) == {C.STRIPPED_WHITESPACE}


def test_ops_blocklist_via_set_algebra():
    # Run everything except whitespace stripping
    result, fixes = clean_hgvs(
        "NM_000059.4:c.316+5G>A p.Arg106*",
        ops=ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_PROTEIN_SUFFIX},
    )
    # protein suffix NOT stripped because that op was excluded
    assert "p.Arg106*" in result
    assert C.STRIPPED_PROTEIN_SUFFIX not in codes(fixes)


def test_ops_empty_set_is_noop():
    s = "nm_000059.4 c.316+5G>A"
    result, fixes = clean_hgvs(s, ops=set(), validate=False)
    assert result == s
    assert fixes == []


def test_validate_false_skips_error_fixes():
    # "ins" with no sequence is an ERROR from validation
    _result, fixes = clean_hgvs("NM_000059.4:c.123ins", validate=False)
    assert E not in severities(fixes)
    _result, fixes = clean_hgvs("NM_000059.4:c.123ins", validate=True)
    assert E in severities(fixes)


def test_all_clean_ops_covers_every_op():
    assert ALL_CLEAN_OPS == frozenset(HGVSCleanOp)


# ---------------------------------------------------------------------------
# clean_hgvs — ported from vg_code/test_hgvs.py::test_clean_hgvs
# ---------------------------------------------------------------------------

CLEAN_CASES = [
    # (input,                                    expected_output,                              required_codes)
    # Missing colon (no version)
    ("NM_205768 c.44A>G",                        "NM_205768:c.44A>G",                          {C.STRIPPED_WHITESPACE}),
    # Missing dot after kind
    ("NM_005629.3:c1403A>C",                     "NM_005629.3:c.1403A>C",                      set()),
    # Missing colon
    ("NM_001101.4 c.95C>G",                      "NM_001101.4:c.95C>G",                        {C.STRIPPED_WHITESPACE}),
    # Space after colon
    ("NM_00380.3: c.648_649delGA",               "NM_00380.3:c.648_649delGA",                  {C.STRIPPED_WHITESPACE}),
    # Space after g.
    ("NC_000023.10:g. 31496384G>A",              "NC_000023.10:g.31496384G>A",                 {C.STRIPPED_WHITESPACE}),
    # Double colon
    ("NM_004245: :c.337G>T",                     "NM_004245:c.337G>T",                         {C.FIXED_DOUBLE_COLON}),
    # Space after numbers
    ("NC_000017.10:g.21085664 G>C",              "NC_000017.10:g.21085664G>C",                 {C.STRIPPED_WHITESPACE}),
    # Space after g. (variant 2)
    ("NC_000023.10:g. 133547943G>A",             "NC_000023.10:g.133547943G>A",                {C.STRIPPED_WHITESPACE}),
    # Missing underscore + missing colon + missing "g." + space
    ("NC000002.10g39139341 C>T",                 "NC_000002.10:g.39139341C>T",                 {C.ADDED_TRANSCRIPT_UNDERSCORE, C.STRIPPED_WHITESPACE}),
    # Unbalanced bracket (trailing)
    ("NM_001754.5):c.557T>A",                    "NM_001754.5:c.557T>A",                       {C.STRIPPED_UNBALANCED_BRACKETS}),
    # Unbalanced bracket (leading)
    ("(NM_004991.4:c.2577+4A>T",                 "NM_004991.4:c.2577+4A>T",                    {C.STRIPPED_UNBALANCED_BRACKETS}),
    # Good brackets — gene symbol should be preserved
    ("NM_001754.5(RUNX1):c.1415T>C",             "NM_001754.5(RUNX1):c.1415T>C",              set()),
    # Uppercase mutation type
    ("NM_032638:c.1126_1133DUP",                 "NM_032638:c.1126_1133dup",                   {C.LOWERCASED_MUTATION_TYPE}),
    # Missing "c."
    ("NM_001754.5:557T>A",                       "NM_001754.5:c.557T>A",                       {C.ADDED_MISSING_KIND}),
    # Missing "." after c
    ("NM_001754.5:c557T>A",                      "NM_001754.5:c.557T>A",                       set()),
    # Has gene, missing "c."
    ("NM_001754.5(RUNX1):557T>A",                "NM_001754.5(RUNX1):c.557T>A",               {C.ADDED_MISSING_KIND}),
    # Has gene, missing "." after c
    ("NM_001754.5(RUNX1):c557T>A",               "NM_001754.5(RUNX1):c.557T>A",              set()),
    # Missing "g."
    ("NC_000007.13:117199563G>T",                "NC_000007.13:g.117199563G>T",                {C.ADDED_MISSING_KIND}),
    # Missing "." after g
    ("NC_000007.13:g117199563G>T",               "NC_000007.13:g.117199563G>T",               set()),
    # Missing "." with "-" in position (c- should become c.-)
    ("NM_000350.2(ABCA4):c-52delC",              "NM_000350.2(ABCA4):c.-52delC",              set()),
    # Gene symbol contains "G" which looks like g. — should not be confused
    ("NM_003560.2(PLA2G6):c.2221C>T",            "NM_003560.2(PLA2G6):c.2221C>T",             set()),
    # Gene symbol run together with transcript (no parens)
    ("NM_003560.2PLA2G6:c.2221C>T",              "NM_003560.2(PLA2G6):c.2221C>T",             set()),
    # #112 — trailing stray quote (search-box / copy-paste artifact)
    ("NM_000059.4:c.316+5G>A'",                  "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    # #112 — wrapping quotes both ends
    ("'NM_001754.5(RUNX1):c.1415T>C'",           "NM_001754.5(RUNX1):c.1415T>C",              {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    # #112 — trailing separator
    ("NM_000059.4:c.316+5G>A;",                  "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    # #112 — doubled dot in transcript version
    ("NM_000059..4:c.316+5G>A",                  "NM_000059.4:c.316+5G>A",                     {C.FIXED_DOUBLE_DOT}),
    # #112 — doubled dot after kind
    ("NM_001754.5(RUNX1):c..1415T>C",            "NM_001754.5(RUNX1):c.1415T>C",              {C.FIXED_DOUBLE_DOT}),
    # #112 — misplaced colon in transcript prefix
    ("NM:_000059.4:c.316+5G>A",                  "NM_000059.4:c.316+5G>A",                     {C.FIXED_PREFIX_COLON}),
    # #112 — leading junk before the transcript accession
    ("#NM_000059.4:c.316+5G>A",                  "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_LEADING_JUNK}),
    (":NM_000059.4:c.316+5G>A",                  "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_LEADING_JUNK}),
    ("GRCh38.p2 NM_000059.4:c.316+5G>A",         "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_LEADING_JUNK}),
    # #112 — gene symbol wedged between extra colons / stray parens
    ("NM_001754.5:(RUNX1):c.1415T>C",            "NM_001754.5(RUNX1):c.1415T>C",              {C.FIXED_GENE_WRAPPER}),
    ("NM_001754.5:RUNX1:c.1415T>C",              "NM_001754.5(RUNX1):c.1415T>C",              {C.FIXED_GENE_WRAPPER}),
    ("NM_001754.5(RUNX1:):c.1415T>C",            "NM_001754.5(RUNX1):c.1415T>C",              {C.FIXED_GENE_WRAPPER}),
    # #112 — comma for the kind dot, period for the substitution '>'
    ("NM_000059.4:c,123del",                     "NM_000059.4:c.123del",                       {C.FIXED_SEPARATOR_TYPO}),
    ("NM_000059.4:c.123C.T",                     "NM_000059.4:c.123C>T",                       {C.FIXED_SEPARATOR_TYPO}),
    # #112 — redundant base count after a ranged del/dup, normalised to plain del/dup
    ("NM_000059.4:c.1315_1337dup23",             "NM_000059.4:c.1315_1337dup",                 {C.DROPPED_DEL_DUP_COUNT}),
    ("NM_000059.4:c.123_127del5",                "NM_000059.4:c.123_127del",                   {C.DROPPED_DEL_DUP_COUNT}),
    # #112 — trailing comma separator
    ("NM_000059.4:c.316+5G>A,",                  "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    # Gene symbol contains a mutation-type substring (INS) AND the allele has an
    # uppercase mutation type — only the allele must be lowercased, not the gene.
    ("NM_000208.4(INSR):c.215_216DEL",           "NM_000208.4(INSR):c.215_216del",            {C.LOWERCASED_MUTATION_TYPE}),
]


@pytest.mark.parametrize("bad,expected,required_codes", CLEAN_CASES)
def test_clean_hgvs(bad, expected, required_codes):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected, f"Input: {bad!r}"
    fix_codes = codes(fixes)
    for code in required_codes:
        assert code in fix_codes, (
            f"Expected fix code {code.value!r} for input {bad!r}, got {[c.value for c in fix_codes]}"
        )


# ---------------------------------------------------------------------------
# clean_hgvs — gene/transcript swap and case
# Ported from vg_code/test_hgvs.py::test_fix_gene_transcript
# ---------------------------------------------------------------------------

GENE_TRANSCRIPT_CASES = [
    ("nm_000059.4:c.316+5G>A",           "NM_000059.4:c.316+5G>A",          {C.UPPERCASED_TRANSCRIPT}),
    ("nm_000059.4(BRCA1):c.316+5G>A",    "NM_000059.4(BRCA1):c.316+5G>A",   {C.UPPERCASED_TRANSCRIPT}),
    ("BRCA1(NM_000059.4):c.316+5G>A",    "NM_000059.4(BRCA1):c.316+5G>A",   {C.SWAPPED_GENE_TRANSCRIPT}),
    ("BRCA1(nm_000059.4):c.316+5G>A",    "NM_000059.4(BRCA1):c.316+5G>A",   {C.SWAPPED_GENE_TRANSCRIPT, C.UPPERCASED_TRANSCRIPT}),
]


@pytest.mark.parametrize("bad,expected,required_codes", GENE_TRANSCRIPT_CASES)
def test_fix_gene_transcript(bad, expected, required_codes):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected, f"Input: {bad!r}"
    fix_codes = codes(fixes)
    for code in required_codes:
        assert code in fix_codes, (
            f"Expected fix code {code.value!r} for input {bad!r}, got {[c.value for c in fix_codes]}"
        )


# ---------------------------------------------------------------------------
# clean_hgvs — protein suffix removal
# ---------------------------------------------------------------------------

def test_drop_del_dup_count_requires_range():
    # With a range the count is redundant -> dropped to plain dup
    cleaned, fixes = clean_hgvs("NM_000059.4:c.1315_1337dup23", validate=False)
    assert cleaned == "NM_000059.4:c.1315_1337dup"
    assert C.DROPPED_DEL_DUP_COUNT in codes(fixes)
    # Single position: "del5" means delete 5 bases -> must NOT be altered
    cleaned, fixes = clean_hgvs("NM_000059.4:c.123del5", validate=False)
    assert cleaned == "NM_000059.4:c.123del5"
    assert C.DROPPED_DEL_DUP_COUNT not in codes(fixes)


def test_clean_hgvs_protein_suffix():
    cleaned, fixes = clean_hgvs("NM_000059.4:c.316+5G>A p.Arg106*")
    assert cleaned == "NM_000059.4:c.316+5G>A"
    assert C.STRIPPED_PROTEIN_SUFFIX in codes(fixes)


def test_clean_hgvs_protein_suffix_in_parens():
    cleaned, fixes = clean_hgvs("NM_000059.4:c.316+5G>A (p.Arg106*)")
    assert cleaned == "NM_000059.4:c.316+5G>A"
    assert C.STRIPPED_PROTEIN_SUFFIX in codes(fixes)


# ---------------------------------------------------------------------------
# clean_hgvs — validation errors
# ---------------------------------------------------------------------------

def test_validate_ins_integer_length():
    _cleaned, fixes = clean_hgvs("NM_000441.2(SLC26A4):c.1246_2341ins23")
    error_fixes = [f for f in fixes if f.severity == E]
    assert any(f.code == C.INS_WITH_INTEGER_LENGTH for f in error_fixes)


def test_validate_ins_missing_sequence():
    _cleaned, fixes = clean_hgvs("NM_000441.2:c.100ins")
    error_fixes = [f for f in fixes if f.severity == E]
    assert any(f.code == C.INS_MISSING_SEQUENCE for f in error_fixes)


def test_raise_on_errors():
    with pytest.raises(HGVSInputError):
        clean_hgvs("NM_000441.2(SLC26A4):c.1246_2341ins23", raise_on_errors=True)


def test_no_raise_on_warnings_only():
    # Should not raise even though fixes are made
    cleaned, fixes = clean_hgvs("nm_000059.4:c.316+5G>A", raise_on_errors=True)
    assert cleaned == "NM_000059.4:c.316+5G>A"
    assert all(f.severity == W for f in fixes)


# ---------------------------------------------------------------------------
# clean_hgvs — already-valid strings pass through unchanged
# ---------------------------------------------------------------------------

ALREADY_VALID = [
    "NM_000059.4:c.316+5G>A",
    "NM_001145661.2(GATA2):c.1121G>A",
    "NM_001145661.2(GATA2):c.1142del",
    "NM_001145661.2(GATA2):c.1035_1036insTCTGGCC",
    "NM_001145661.2(GATA2):c.890_903dup",
    "NC_000017.10:g.21085664G>C",
    "GATA2:c.1113dup",
    # Gene symbols containing a mutation-type substring must not be corrupted
    # (INSR -> insR / INVS -> invS were the bug); both transcript-qualified and
    # bare gene-symbol forms must pass through untouched.
    "NM_000208.4(INSR):c.215A>G",
    "NM_014425.5(INVS):c.100A>G",
    "INSR:c.215A>G",
    # Multiple *balanced* brackets are valid HGVS (gene symbol in parens plus
    # uncertain-range notation) and must NOT be stripped as "unbalanced".
    "NM_004006.2(DMD):c.(4071+1_4072-1)_(5154+1_5155-1)del",
    "NC_000023.11:g.(31180435_31200854)_(33274278_33357726)del",
]


@pytest.mark.parametrize("hgvs_string", ALREADY_VALID)
def test_already_valid_unchanged(hgvs_string):
    cleaned, fixes = clean_hgvs(hgvs_string)
    assert cleaned == hgvs_string, f"Valid HGVS was modified: {hgvs_string!r} → {cleaned!r}"
    warning_fixes = [f for f in fixes if f.severity == W]
    assert not warning_fixes, (
        f"Unexpected warnings for valid HGVS {hgvs_string!r}: "
        f"{[f.code.value for f in warning_fixes]}"
    )


# ---------------------------------------------------------------------------
# HGVSFix structure
# ---------------------------------------------------------------------------

def test_fix_has_original_and_fixed():
    """Fixes should carry before/after values where applicable."""
    _cleaned, fixes = clean_hgvs("nm_000059.4:c.316+5G>A")
    transcript_fix = next((f for f in fixes if f.code == C.UPPERCASED_TRANSCRIPT), None)
    assert transcript_fix is not None
    assert transcript_fix.original is not None
    assert transcript_fix.fixed is not None


# ---------------------------------------------------------------------------
# _rank_versions
# ---------------------------------------------------------------------------

def test_rank_versions_up_then_down():
    result = _rank_versions(4, [1, 2, 3, 5, 6], VersionStrategy.UP_THEN_DOWN)
    assert result == [5, 6, 3, 2, 1]


def test_rank_versions_closest():
    result = _rank_versions(4, [1, 2, 3, 5, 6], VersionStrategy.CLOSEST)
    assert result == [5, 3, 6, 2, 1]


def test_rank_versions_latest():
    result = _rank_versions(4, [1, 2, 3, 5, 6], VersionStrategy.LATEST)
    assert result == [6, 5, 3, 2, 1]


def test_rank_versions_exact_match_not_special():
    # _rank_versions doesn't exclude the requested version — it just sorts.
    # get_best_transcript_version handles the exact-match short-circuit.
    # For UP_THEN_DOWN with requested=3, v=3 has distance=0 and prefer_later=False,
    # so it sorts before v=5 (distance=2, prefer_later=False).
    result = _rank_versions(3, [1, 2, 3, 5], VersionStrategy.UP_THEN_DOWN)
    assert result == [3, 5, 2, 1]


def test_rank_transcript_versions_public_name():
    # Promoted from _rank_versions (#114); the private name stays as an alias.
    assert rank_transcript_versions is _rank_versions
    # Strategy defaults to UP_THEN_DOWN
    assert rank_transcript_versions(4, [1, 2, 3, 5, 6]) == [5, 6, 3, 2, 1]


# ---------------------------------------------------------------------------
# HGVSFix.__str__ and message helpers (#114)
# ---------------------------------------------------------------------------

def _sample_fixes() -> list[HGVSFix]:
    return [
        HGVSFix(W, C.STRIPPED_WHITESPACE, "Removed whitespace"),
        HGVSFix(E, C.NO_COLON, "No colon ':' found in HGVS string"),
        HGVSFix(W, C.UPPERCASED_TRANSCRIPT, "Upper cased transcript"),
    ]


def test_hgvs_fix_str_is_message():
    fix = HGVSFix(W, C.STRIPPED_WHITESPACE, "Removed whitespace")
    assert str(fix) == "Removed whitespace"


def test_warning_messages():
    assert warning_messages(_sample_fixes()) == [
        "Removed whitespace", "Upper cased transcript",
    ]


def test_error_messages():
    assert error_messages(_sample_fixes()) == ["No colon ':' found in HGVS string"]


def test_messages_all_and_filtered():
    fixes = _sample_fixes()
    assert messages(fixes) == [f.message for f in fixes]  # no severity → all
    assert messages(fixes, W) == warning_messages(fixes)
    assert messages(fixes, E) == error_messages(fixes)


# ---------------------------------------------------------------------------
# get_best_transcript_version
# ---------------------------------------------------------------------------

def test_get_best_version_exact_match():
    version, fix = get_best_transcript_version("NM_000059", 4, [1, 2, 4, 5])
    assert version == 4
    assert fix is None


def test_get_best_version_fallback_up_then_down():
    version, fix = get_best_transcript_version("NM_000059", 4, [1, 2, 3, 5, 6])
    assert version == 5
    assert fix is not None
    assert fix.severity == W
    assert fix.code == C.USED_ADJACENT_VERSION
    assert fix.original == "NM_000059.4"
    assert fix.fixed == "NM_000059.5"


def test_get_best_version_fallback_closest():
    version, fix = get_best_transcript_version(
        "NM_000059", 4, [1, 2, 3, 5, 6], strategy=VersionStrategy.CLOSEST
    )
    assert version == 5
    assert fix is not None


def test_get_best_version_fallback_latest():
    version, fix = get_best_transcript_version(
        "NM_000059", 4, [1, 2, 3], strategy=VersionStrategy.LATEST
    )
    assert version == 3


def test_get_best_version_empty_raises():
    with pytest.raises(HGVSInputError):
        get_best_transcript_version("NM_000059", 4, [])


def test_get_best_version_message():
    _version, fix = get_best_transcript_version("NM_000059", 4, [1, 2, 3, 5, 6])
    assert "NM_000059.4" in fix.message
    assert "NM_000059.5" in fix.message or "5" in fix.message
