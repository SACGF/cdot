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
    ("NM_004245: :c.337G>T",                     "NM_004245:c.337G>T",                         {C.FIXED_MULTIPLE_COLON}),
    # Run of 3+ colons collapses to one
    ("NM_004245:::c.337G>T",                      "NM_004245:c.337G>T",                         {C.FIXED_MULTIPLE_COLON}),
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
    # Missing kind on non-coding RNA (NR_/XR_) → "n.", not "c."
    ("NR_003051.3:601A>G",                        "NR_003051.3:n.601A>G",                       {C.ADDED_MISSING_KIND}),
    ("XR_001737837.2:601A>G",                     "XR_001737837.2:n.601A>G",                    {C.ADDED_MISSING_KIND}),
    # Missing kind on predicted mRNA (XM_) and Ensembl (ENST) → "c."
    ("XM_005260000.3:557T>A",                     "XM_005260000.3:c.557T>A",                    {C.ADDED_MISSING_KIND}),
    ("ENST00000357654.9:557T>A",                  "ENST00000357654.9:c.557T>A",                 {C.ADDED_MISSING_KIND}),
    # Missing kind on other genomic prefixes (NG_/NW_/NT_) → "g."
    ("NG_012337.3:5000A>G",                       "NG_012337.3:g.5000A>G",                      {C.ADDED_MISSING_KIND}),
    ("NW_009646201.1:100G>T",                     "NW_009646201.1:g.100G>T",                    {C.ADDED_MISSING_KIND}),
    # Missing kind with a UTR coordinate ("-" 5' UTR, "*" 3' UTR)
    ("NM_000350.2:-52delC",                       "NM_000350.2:c.-52delC",                      {C.ADDED_MISSING_KIND}),
    ("NM_000350.2:*52A>G",                        "NM_000350.2:c.*52A>G",                       {C.ADDED_MISSING_KIND}),
    # Missing leading "N" on transcript/genomic accessions (XM_/XR_ have no N, so
    # the dropped letter is always N: NM/NC/NR/NG/NT/NW/NS/NZ/NP)
    ("M_001754.5:c.557T>A",                       "NM_001754.5:c.557T>A",                       {C.ADDED_N_PREFIX}),
    ("R_003051.3:n.601A>G",                       "NR_003051.3:n.601A>G",                       {C.ADDED_N_PREFIX}),
    ("G_012337.3:g.5000A>G",                      "NG_012337.3:g.5000A>G",                      {C.ADDED_N_PREFIX}),
    # Missing underscore on non-NM/NC prefixes (needs 6+ digits to fire)
    ("NR003051.3:n.601A>G",                       "NR_003051.3:n.601A>G",                       {C.ADDED_TRANSCRIPT_UNDERSCORE}),
    ("XM005260000.3:c.557T>A",                    "XM_005260000.3:c.557T>A",                    {C.ADDED_TRANSCRIPT_UNDERSCORE}),
    ("XR001737837.2:n.601A>G",                    "XR_001737837.2:n.601A>G",                    {C.ADDED_TRANSCRIPT_UNDERSCORE}),
    ("NG012337.3:g.5000A>G",                      "NG_012337.3:g.5000A>G",                      {C.ADDED_TRANSCRIPT_UNDERSCORE}),
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
    ("NM_000059..4:c.316+5G>A",                  "NM_000059.4:c.316+5G>A",                     {C.FIXED_MULTIPLE_DOT}),
    # #112 — doubled dot after kind
    ("NM_001754.5(RUNX1):c..1415T>C",            "NM_001754.5(RUNX1):c.1415T>C",              {C.FIXED_MULTIPLE_DOT}),
    # Run of 3+ dots collapses to one
    ("NM_001754.5(RUNX1):c....1415T>C",          "NM_001754.5(RUNX1):c.1415T>C",              {C.FIXED_MULTIPLE_DOT}),
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
    # #112 — colon used in place of the kind dot ("c:" -> "c.")
    ("NM_000059.4:c:123A>G",                     "NM_000059.4:c.123A>G",                       {C.FIXED_SEPARATOR_TYPO}),
    # Repeated kind prefix collapses to one ("c.c." and "c.c.c." both -> "c.")
    ("NM_000059.4:c.c.123A>G",                   "NM_000059.4:c.123A>G",                       {C.FIXED_MULTIPLE_KIND}),
    ("NM_000059.4:c.c.c.123A>G",                 "NM_000059.4:c.123A>G",                       {C.FIXED_MULTIPLE_KIND}),
    # #112 — parens wrapping the whole accession
    ("(NM_001754.5):c.1415T>C",                  "NM_001754.5:c.1415T>C",                      {C.FIXED_GENE_WRAPPER}),
    # Kind letter glued to the version reconstructs the same way whether or not a
    # stray colon sits after it ("2c1188" and "2c:1188" both -> "2:c.1188").
    ("NM_001017995.2c1188+1773_2733+6592del",    "NM_001017995.2:c.1188+1773_2733+6592del",    {C.RECONSTRUCTED_STRUCTURE}),
    ("NM_001017995.2c:1188+1773_2733+6592del",   "NM_001017995.2:c.1188+1773_2733+6592del",    {C.RECONSTRUCTED_STRUCTURE}),
    # #112 — genomic accession wedged into the gene-symbol parenthetical slot
    ("NM_000059.4(NC_000013.11):c.68del",        "NM_000059.4:c.68del",                        {C.DROPPED_GENOMIC_REF_IN_PARENS}),
    ("NM_001754.5(NG_042763.1):c.1415T>C",       "NM_001754.5:c.1415T>C",                      {C.DROPPED_GENOMIC_REF_IN_PARENS}),
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
    # #112 — gene-first wrapper missing the colon before the kind ")c." -> "):c."
    ("BRCA1(NM_000059.4)c.316+5G>A",     "NM_000059.4(BRCA1):c.316+5G>A",   {C.SWAPPED_GENE_TRANSCRIPT}),
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


def test_drop_genomic_ref_in_parens():
    # #112 — a genomic accession (NC_/NG_/NW_) sitting in the gene-symbol
    # parenthetical slot is not a gene symbol; dropping it lets the string parse.
    for genomic in ("NC_000013.11", "NG_012772.3", "NW_009646201.1"):
        bad = f"NM_000059.4({genomic}):c.68del"
        cleaned, fixes = clean_hgvs(bad, validate=False)
        assert cleaned == "NM_000059.4:c.68del", f"Input: {bad!r}"
        assert C.DROPPED_GENOMIC_REF_IN_PARENS in codes(fixes)


def test_drop_genomic_ref_in_parens_keeps_real_gene_symbol():
    # A real gene symbol in the parens must be kept, not dropped.
    cleaned, fixes = clean_hgvs("NM_000059.4(BRCA2):c.68del", validate=False)
    assert cleaned == "NM_000059.4(BRCA2):c.68del"
    assert C.DROPPED_GENOMIC_REF_IN_PARENS not in codes(fixes)


def test_drop_genomic_ref_in_parens_collapses_genomic_selector_form():
    # The "genomic reference, transcript selector" form NC_(NM_) is first
    # swapped to NM_(NC_) by the gene/transcript step, then the redundant
    # genomic parenthetical is dropped — collapsing to the same parseable
    # transcript reference (a semantically equivalent simplification).
    cleaned, fixes = clean_hgvs("NC_000013.11(NM_000059.4):c.68del", validate=False)
    assert cleaned == "NM_000059.4:c.68del"
    assert C.DROPPED_GENOMIC_REF_IN_PARENS in codes(fixes)


def test_clean_hgvs_protein_suffix():
    cleaned, fixes = clean_hgvs("NM_000059.4:c.316+5G>A p.Arg106*")
    assert cleaned == "NM_000059.4:c.316+5G>A"
    assert C.STRIPPED_PROTEIN_SUFFIX in codes(fixes)


def test_clean_hgvs_protein_suffix_in_parens():
    cleaned, fixes = clean_hgvs("NM_000059.4:c.316+5G>A (p.Arg106*)")
    assert cleaned == "NM_000059.4:c.316+5G>A"
    assert C.STRIPPED_PROTEIN_SUFFIX in codes(fixes)


def test_structure_reconstruction_not_reported_as_uppercased_bases():
    # Missing colon between an already-uppercase transcript and its c. kind:
    # the only change is inserting the ':', so the fix must be reported as a
    # structural reconstruction, NOT as uppercased_bases (no bases changed).
    cleaned, fixes = clean_hgvs("NM_000059.4c.100+5_200-3del")
    assert cleaned == "NM_000059.4:c.100+5_200-3del"
    assert C.RECONSTRUCTED_STRUCTURE in codes(fixes)
    assert C.UPPERCASED_BASES not in codes(fixes)


def test_multiple_underscore_not_reported_as_added_underscore():
    # Collapsing repeated underscores removes them, it does not add a missing
    # one, so it must be reported as FIXED_MULTIPLE_UNDERSCORE, not the
    # ADDED_TRANSCRIPT_UNDERSCORE code used when a separator is actually added.
    cleaned, fixes = clean_hgvs("NM__000059.4:c.123del")
    assert cleaned == "NM_000059.4:c.123del"
    assert C.FIXED_MULTIPLE_UNDERSCORE in codes(fixes)
    assert C.ADDED_TRANSCRIPT_UNDERSCORE not in codes(fixes)
    # A run of 3+ collapses to a single underscore in one pass
    cleaned, fixes = clean_hgvs("NM___000059.4:c.123del")
    assert cleaned == "NM_000059.4:c.123del"
    assert C.FIXED_MULTIPLE_UNDERSCORE in codes(fixes)


# ---------------------------------------------------------------------------
# clean_hgvs — add missing kind, picking the kind from the accession type
# ---------------------------------------------------------------------------

# (accession_prefix_example, expected_kind_letter)
ADD_MISSING_KIND_CASES = [
    # mRNA / Ensembl transcript → c.
    ("NM_001754.5",     "c"),
    ("XM_005260000.3",  "c"),
    ("ENST00000357654.9", "c"),
    # non-coding RNA → n.
    ("NR_003051.3",     "n"),
    ("XR_001737837.2",  "n"),
    # genomic → g.
    ("NC_000007.13",    "g"),
    ("NG_012337.3",     "g"),
    ("NW_009646201.1",  "g"),
    ("NT_187633.1",     "g"),
    ("AC_000156.1",     "g"),
]


@pytest.mark.parametrize("accession,kind", ADD_MISSING_KIND_CASES)
def test_add_missing_kind_by_reference_type(accession, kind):
    cleaned, fixes = clean_hgvs(f"{accession}:601A>G", validate=False)
    assert cleaned == f"{accession}:{kind}.601A>G"
    fix = next(f for f in fixes if f.code == C.ADDED_MISSING_KIND)
    assert fix.message == f"Added missing '{kind}.' kind prefix"


@pytest.mark.parametrize("accession,kind", ADD_MISSING_KIND_CASES)
def test_add_missing_kind_with_gene_symbol(accession, kind):
    # The gene-symbol parenthetical must not stop the kind being inserted.
    cleaned, _ = clean_hgvs(f"{accession}(GENE):601A>G", validate=False)
    assert cleaned == f"{accession}(GENE):{kind}.601A>G"


@pytest.mark.parametrize("coord", ["-52", "*52", "123"])
def test_add_missing_kind_utr_and_normal_coords(coord):
    # The coordinate can start with a digit, "-" (5' UTR) or "*" (3' UTR).
    cleaned, fixes = clean_hgvs(f"NM_000350.2:{coord}delC", validate=False)
    assert cleaned == f"NM_000350.2:c.{coord}delC"
    assert C.ADDED_MISSING_KIND in codes(fixes)


def test_add_missing_kind_lowercase_prefix():
    # Lowercase accession with no kind: kind is added and the prefix uppercased.
    cleaned, fixes = clean_hgvs("nr_003051.3:601A>G", validate=False)
    assert cleaned == "NR_003051.3:n.601A>G"
    assert C.ADDED_MISSING_KIND in codes(fixes)


# ---------------------------------------------------------------------------
# clean_hgvs — restore a dropped leading "N"
# ---------------------------------------------------------------------------

# (broken, expected) — every N-prefixed RefSeq type with the leading N dropped
ADD_N_PREFIX_CASES = [
    ("M_001754.5:c.557T>A",   "NM_001754.5:c.557T>A"),
    ("R_003051.3:n.601A>G",   "NR_003051.3:n.601A>G"),
    ("C_000007.13:g.1A>G",    "NC_000007.13:g.1A>G"),
    ("G_012337.3:g.5000A>G",  "NG_012337.3:g.5000A>G"),
    ("T_187633.1:g.1A>G",     "NT_187633.1:g.1A>G"),
    ("W_009646201.1:g.1A>G",  "NW_009646201.1:g.1A>G"),
    ("P_000050.2:p.Arg1Gly",  "NP_000050.2:p.Arg1Gly"),
]


@pytest.mark.parametrize("broken,expected", ADD_N_PREFIX_CASES)
def test_add_n_prefix_all_refseq_types(broken, expected):
    cleaned, fixes = clean_hgvs(broken, validate=False)
    assert cleaned == expected
    assert C.ADDED_N_PREFIX in codes(fixes)


@pytest.mark.parametrize("not_an_accession", [
    "X_001234.5:c.1A>G",   # NX_ is not a real prefix
    "A_001234.5:c.1A>G",   # NA_ is not a real prefix
    "B_001234.5:c.1A>G",   # NB_ is not a real prefix
])
def test_add_n_prefix_skips_unknown_prefixes(not_an_accession):
    cleaned, _ = clean_hgvs(not_an_accession, validate=False)
    assert not cleaned.startswith("N")


# ---------------------------------------------------------------------------
# clean_hgvs — add a missing underscore to any RefSeq prefix
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("broken,expected", [
    ("NM000059.4:c.1A>G",      "NM_000059.4:c.1A>G"),
    ("NC000007.13:g.1A>G",     "NC_000007.13:g.1A>G"),
    ("NR003051.3:n.1A>G",      "NR_003051.3:n.1A>G"),
    ("XM005260000.3:c.1A>G",   "XM_005260000.3:c.1A>G"),
    ("XR001737837.2:n.1A>G",   "XR_001737837.2:n.1A>G"),
    ("NG012337.3:g.1A>G",      "NG_012337.3:g.1A>G"),
    ("NW009646201.1:g.1A>G",   "NW_009646201.1:g.1A>G"),
])
def test_add_transcript_underscore_all_refseq_types(broken, expected):
    cleaned, fixes = clean_hgvs(broken, validate=False)
    assert cleaned == expected
    assert C.ADDED_TRANSCRIPT_UNDERSCORE in codes(fixes)


# Gene symbols that begin with a RefSeq accession prefix plus a digit. None has
# 6+ consecutive digits, so neither the underscore nor the reconstruct step may
# touch them. (NR3C1 is a real, important gene — the regression this guards.)
GENE_SYMBOL_LOOKALIKES = [
    "NR3C1:c.1411A>G",
    "NR1H3:c.100A>G",
    "NS1:c.100A>G",
    "ZP3:c.100A>G",
    "NC2:c.5A>G",
    "XRCC1:c.100A>G",
    "NM2:c.5A>G",
    "AP3B1:c.100A>G",
]


@pytest.mark.parametrize("hgvs_string", GENE_SYMBOL_LOOKALIKES)
def test_gene_symbol_lookalikes_not_mangled(hgvs_string):
    cleaned, fixes = clean_hgvs(hgvs_string, validate=False)
    assert cleaned == hgvs_string, f"Gene symbol mangled: {hgvs_string!r} → {cleaned!r}"
    assert C.ADDED_TRANSCRIPT_UNDERSCORE not in codes(fixes)
    assert C.RECONSTRUCTED_STRUCTURE not in codes(fixes)


def test_six_digits_required_for_underscore():
    # 5 digits is not a real accession (a gene-like token) → left alone;
    # 6 digits is the RefSeq minimum → underscore added.
    cleaned, _ = clean_hgvs("NR12345:c.1A>G", validate=False)
    assert cleaned == "NR12345:c.1A>G"
    cleaned, fixes = clean_hgvs("NR123456:n.1A>G", validate=False)
    assert cleaned == "NR_123456:n.1A>G"
    assert C.ADDED_TRANSCRIPT_UNDERSCORE in codes(fixes)


# ---------------------------------------------------------------------------
# clean_hgvs — validation errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("leading", ["c", "g", "n", "m", "p", "r", ".", ":"])
def test_validate_missing_reference_for_all_kind_letters(leading):
    # A string that starts with a bare kind letter (or ':'/'.') has no reference
    # sequence and must be flagged, for every kind letter cdot recognises.
    _cleaned, fixes = clean_hgvs(f"{leading}123A>G")
    assert C.MISSING_REFERENCE_SEQUENCE in codes(fixes)

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
    # #112 — LRG transcript references (the "t1"/"p1" suffix must survive
    # structure reconstruction, not be mistaken for the gene symbol).
    "LRG_199t1(RUNX1):c.1415T>C",
    "LRG_199t2(BRCA2):c.316+5G>A",
    # Gene symbols that start with a RefSeq accession prefix plus a digit must NOT
    # be mistaken for an underscore-less accession (NR3 + "C1", NS1, etc.) and
    # mangled. Real accessions always have an underscore and 6+ digits.
    "NR3C1:c.1411A>G",
    "NS1:c.100A>G",
    "ZP3:c.100A>G",
    "XRCC1:c.100A>G",
    "NM_001145661.2(NR3C1):c.1121G>A",
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
