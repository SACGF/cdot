"""
Tests for cdot.hgvs.clean — pure string operations, no data provider required.

Test cases ported from vg_code/test_hgvs.py (test_clean_hgvs and
test_fix_gene_transcript) and extended for the version-fallback helpers.
"""
import warnings

import pytest

from cdot.hgvs.clean import (
    ALL_CLEAN_OPS,
    HGVSCleanOpWarning,
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
    ("GRCh38.p2 NM_000059.4:c.316+5G>A",         "NM_000059.4:c.316+5G>A",                     {C.STRIPPED_BUILD_NAME}),
    ("GRCh37.p12#NM_000059.4:c.68del",           "NM_000059.4:c.68del",                        {C.STRIPPED_BUILD_NAME}),
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
    # #112: semicolon/underscore in place of the accession:allele colon
    ("NM_000059.4;c.68del",                      "NM_000059.4:c.68del",                        {C.FIXED_SEPARATOR_TYPO}),
    ("NM_001754.5_c.749G>A",                     "NM_001754.5:c.749G>A",                       {C.FIXED_SEPARATOR_TYPO}),
    ("NM_000059.4; c.68del",                     "NM_000059.4:c.68del",                        {C.FIXED_SEPARATOR_TYPO}),
    # #112: dangling dot of an empty (dropped) transcript version
    ("NM_001754.:c.749G>A",                      "NM_001754:c.749G>A",                         {C.DROPPED_EMPTY_VERSION}),
    ("NM_001754.::c..749G>A",                    "NM_001754:c.749G>A",                         {C.DROPPED_EMPTY_VERSION}),
    # #112: leading colon/dash junk separated from the accession by whitespace
    (": NM_000059.4(BRCA2):c.68del",             "NM_000059.4(BRCA2):c.68del",                 {C.STRIPPED_LEADING_JUNK}),
    ("- NM_001754.5:c.749G>A",                   "NM_001754.5:c.749G>A",                       {C.STRIPPED_LEADING_JUNK}),
    # #112: space in place of the reference:kind colon after a gene symbol.
    # Must also work for symbols that resemble an accession (2-letter RefSeq
    # prefix + digit, eg NR2E3), which reconstruction deliberately refuses.
    ("BRCA2 c.68del",                            "BRCA2:c.68del",                              {C.FIXED_SEPARATOR_TYPO}),
    ("NR2E3 c.349G>A",                           "NR2E3:c.349G>A",                             {C.FIXED_SEPARATOR_TYPO}),
    # #148: gene symbol, space, then transcript accession -> accession(symbol)
    ("APC NM_001407446.1:c.-40G>A",              "NM_001407446.1(APC):c.-40G>A",               {C.WRAPPED_GENE_SYMBOL}),
    ("BRCA2 NM_000059.4 c.68del",                "NM_000059.4(BRCA2):c.68del",                 {C.WRAPPED_GENE_SYMBOL}),
    ("BRCA2 NM_000059:c.68del",                  "NM_000059(BRCA2):c.68del",                   {C.WRAPPED_GENE_SYMBOL}),
    ("RUNX1  NM_001754.5:c.1415T>C",             "NM_001754.5(RUNX1):c.1415T>C",               {C.WRAPPED_GENE_SYMBOL}),
    # #165: a genome build name before the reference is dropped, not taken for a
    # gene symbol (was "NM_000059.4(hg19):c.68del" / "GRCh38BRCA2NM_000059.4:c.68del")
    ("hg19 NM_000059.4:c.68del",                 "NM_000059.4:c.68del",                        {C.STRIPPED_BUILD_NAME}),
    ("GRCh38 BRCA2 NM_000059.4:c.68del",         "NM_000059.4(BRCA2):c.68del",                 {C.STRIPPED_BUILD_NAME, C.WRAPPED_GENE_SYMBOL}),
    # #225: the accession already has the symbol, so the matching one before it is dropped
    ("BRCA2 NM_000059.4(BRCA2):c.68del",         "NM_000059.4(BRCA2):c.68del",                 {C.WRAPPED_GENE_SYMBOL}),
    # #225: protein annotation after a comma/semicolon goes with its separator
    ("NM_000059.4:c.68del, p.(Glu23fs)",         "NM_000059.4:c.68del",                        {C.STRIPPED_PROTEIN_SUFFIX}),
    ("NM_001754.5(RUNX1):c.1415T>C; p.?",        "NM_001754.5(RUNX1):c.1415T>C",               {C.STRIPPED_PROTEIN_SUFFIX}),
    # #225: bracket typed for the accession underscore
    ("NM)000059.4:c.68del",                      "NM_000059.4:c.68del",                        {C.ADDED_TRANSCRIPT_UNDERSCORE}),
    # #225: trailing stray "(" opens nothing
    ("NM_000059.4:c.68del (",                    "NM_000059.4:c.68del",                        {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    # #227: gene symbol before the accession separated by ':' or '.'
    ("BRCA2:NM_000059.4:c.68del",                "NM_000059.4(BRCA2):c.68del",                 {C.WRAPPED_GENE_SYMBOL}),
    ("BRCA2.NM_000059.4:c.68del",                "NM_000059.4(BRCA2):c.68del",                 {C.WRAPPED_GENE_SYMBOL}),
    ("BRCA2: NM_000059.4:c.68del",               "NM_000059.4(BRCA2):c.68del",                 {C.WRAPPED_GENE_SYMBOL}),
    # #227: dot or colon in place of the accession prefix underscore
    ("NM.000059.4:c.68del",                      "NM_000059.4:c.68del",                        {C.ADDED_TRANSCRIPT_UNDERSCORE}),
    ("NM:000059.4:c.68del",                      "NM_000059.4:c.68del",                        {C.ADDED_TRANSCRIPT_UNDERSCORE}),
    # #227: stray leading bullet/asterisk, trailing dot/footnote mark/backslash
    ("\u2022NM_000059.4:c.68del",                "NM_000059.4:c.68del",                        {C.STRIPPED_LEADING_JUNK}),
    ("*NM_000059.4:c.68del",                     "NM_000059.4:c.68del",                        {C.STRIPPED_LEADING_JUNK}),
    ("NM_000059.4:c.68del.",                     "NM_000059.4:c.68del",                        {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    ("NM_000059.4:c.68del\u2020",                "NM_000059.4:c.68del",                        {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    ("NM_000059.4:c.68del\\",                   "NM_000059.4:c.68del",                        {C.STRIPPED_SURROUNDING_PUNCTUATION}),
    # #227: typographic dash in place of '-'
    ("NM_000059.4:c.67\u20102A>G",               "NM_000059.4:c.67-2A>G",                      {C.REPLACED_UNICODE_DASH}),
    ("NM_000059.4:c.\u221240G>A",                "NM_000059.4:c.-40G>A",                       {C.REPLACED_UNICODE_DASH}),
    # #227: kind repeated on the second position of a range
    ("NM_000059.4:c.67+3_c.67+6del",             "NM_000059.4:c.67+3_67+6del",                 {C.FIXED_MULTIPLE_KIND}),
    # #227: gene symbol repeated after its parens
    ("NM_000059.4(BRCA2):BRCA2:c.68del",         "NM_000059.4(BRCA2):c.68del",                 {C.FIXED_GENE_WRAPPER}),
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


# ---------------------------------------------------------------------------
# #165: a del/dup count is checked against the range before it is dropped
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad,expected", [
    # Plain integer positions on every kind: the count agrees with the range
    ("NM_000059.4:c.123_127del5",                "NM_000059.4:c.123_127del"),
    ("NM_000059.4:c.1315_1337dup23",             "NM_000059.4:c.1315_1337dup"),
    ("NC_000013.11:g.32316528_32316530del3",     "NC_000013.11:g.32316528_32316530del"),
    ("NR_003051.3:n.601_602dup2",                "NR_003051.3:n.601_602dup"),
    ("NC_012920.1:m.3243_3245del3",              "NC_012920.1:m.3243_3245del"),
    # Kind dot missing: the count is still checked before the structure is rebuilt
    ("NM_000059.4:c123_127del5",                 "NM_000059.4:c.123_127del"),
])
def test_del_dup_count_matching_range_is_dropped(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.DROPPED_DEL_DUP_COUNT in codes(fixes)
    assert E not in severities(fixes)


@pytest.mark.parametrize("bad", [
    "NM_000059.4:c.123_127del3",                 # 5 bp range, count says 3
    "NM_000059.4:c.1315_1337dup3",               # 23 bp range, count says 3
    "NC_000013.11:g.32316528_32316530del5",      # 3 bp range, count says 5
    "NM_000059.4(BRCA2):c.123_127dup7",
])
def test_del_dup_count_contradicting_range_is_error(bad):
    # The count and the range disagree, so which one is right is a guess:
    # neither is changed and the contradiction is reported as an ERROR.
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == bad
    assert C.DEL_DUP_COUNT_MISMATCH in codes(fixes)
    assert C.DROPPED_DEL_DUP_COUNT not in codes(fixes)
    error = next(f for f in fixes if f.code == C.DEL_DUP_COUNT_MISMATCH)
    assert error.severity == E
    assert error.original == bad
    with pytest.raises(HGVSInputError):
        clean_hgvs(bad, raise_on_errors=True)


def test_del_dup_count_mismatch_message_gives_both_lengths():
    _cleaned, fixes = clean_hgvs("NM_000059.4:c.123_127del3")
    error = next(f for f in fixes if f.code == C.DEL_DUP_COUNT_MISMATCH)
    assert error.message == "Base count 3 after del contradicts the 5 bp range 123_127"


@pytest.mark.parametrize("bad,expected", [
    # Intronic offsets: the range's length needs the transcript, so the count
    # cannot be checked and is dropped as before
    ("NM_000059.4:c.100+5_200-3del5",            "NM_000059.4:c.100+5_200-3del"),
    ("NM_000059.4:c.100+5_200-3dup2",            "NM_000059.4:c.100+5_200-3dup"),
    # A range starting in the 5' UTR (there is no position 0, so the length is
    # not end - start + 1) is not checked either
    ("NM_000059.4:c.-30_5dup3",                  "NM_000059.4:c.-30_5dup"),
])
def test_del_dup_count_dropped_unchecked_when_range_length_unknown(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.DROPPED_DEL_DUP_COUNT in codes(fixes)
    assert C.DEL_DUP_COUNT_MISMATCH not in codes(fixes)


@pytest.mark.parametrize("unchanged", [
    "NM_000059.4:c.*10_*20del3",                 # 3' UTR positions
    "NM_000059.4:c.(100_200)del3",               # uncertain range
])
def test_del_dup_count_left_alone_on_utr_or_uncertain_range(unchanged):
    # The count rule does not match these at all, so they pass through untouched
    # and without a count error, as before.
    cleaned, fixes = clean_hgvs(unchanged)
    assert cleaned == unchanged
    assert C.DROPPED_DEL_DUP_COUNT not in codes(fixes)
    assert C.DEL_DUP_COUNT_MISMATCH not in codes(fixes)


# ---------------------------------------------------------------------------
# #165: an unbalanced bracket in the allele is kept, not stripped
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [
    "NM_000059.4:c.(100_200del",                 # closing bracket of an uncertain range lost
    "NM_000059.4:c.100_200)del",                 # opening bracket lost
    "NM_000059.4(BRCA2):c.(100_200del",          # gene wrapper balanced, allele not
    "NM_000059.4:c.(100_200)_(300_400del",
    "NP_000050.2:p.(Glu23fs",                    # predicted protein consequence
])
def test_unbalanced_bracket_in_allele_is_error(bad):
    # Stripping would turn an uncertain range into a certain one (was
    # "NM_000059.4:c.100_200del"), and where the missing bracket belongs is a
    # guess, so the string is left alone and reported as an ERROR.
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == bad
    assert C.UNBALANCED_BRACKETS in codes(fixes)
    assert C.STRIPPED_UNBALANCED_BRACKETS not in codes(fixes)
    error = next(f for f in fixes if f.code == C.UNBALANCED_BRACKETS)
    assert error.severity == E
    assert error.original == bad


@pytest.mark.parametrize("bad,expected", [
    # A stray bracket around the reference is still stripped, and a balanced
    # uncertain range in the allele survives it
    ("(NM_000059.4:c.(100_200)del",              "NM_000059.4:c.(100_200)del"),
    ("NM_000059.4(BRCA2:c.68del",                "NM_000059.4(BRCA2):c.68del"),
    ("NM_000059.4(BRCA2)):c.68del",              "NM_000059.4(BRCA2):c.68del"),
    # A wrapping bracket that closes/opens nothing in the allele is unambiguous
    ("NM_000059.4:c.68del)",                     "NM_000059.4:c.68del"),
    ("NM_000059.4:c.(100_200)del)",              "NM_000059.4:c.(100_200)del"),
    ("NM_000059.4:(c.68del",                     "NM_000059.4:c.68del"),
    # No colon yet: the allele still starts at the kind letter
    ("(NM_000059.4 c.(100_200)del",              "NM_000059.4:c.(100_200)del"),
])
def test_stray_bracket_stripped_keeping_uncertain_range(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.STRIPPED_UNBALANCED_BRACKETS in codes(fixes)
    assert E not in severities(fixes)


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


@pytest.mark.parametrize("bad", [
    "NC_000013.11(NM_000059.4):g.32316528del",
    "NM_000059.4(NC_000013.11):g.32316528del",
])
def test_genomic_kind_keeps_genomic_accession(bad):
    # A g. position is on the genomic accession, so the transcript is the one dropped
    cleaned, fixes = clean_hgvs(bad, validate=False)
    assert cleaned == "NC_000013.11:g.32316528del", f"Input: {bad!r}"
    assert codes(fixes) == {C.DROPPED_TRANSCRIPT_FOR_GENOMIC_KIND}


def test_genomic_kind_keeps_mitochondrial_accession():
    cleaned, _fixes = clean_hgvs("NC_012920.1(MT-ND1):m.3460G>A", validate=False)
    assert cleaned == "NC_012920.1(MT-ND1):m.3460G>A"  # a gene symbol, not a transcript: untouched


def test_accession_separator_semicolon_inside_allele_untouched():
    # A ';' inside a (balanced-bracket) multi-variant allele is real HGVS
    # punctuation, not a separator typo, so the anchored rule must not fire.
    s = "NM_000059.4:c.[68del;70A>G]"
    cleaned, fixes = clean_hgvs(s, validate=False)
    assert cleaned == s
    assert not fixes


def test_empty_version_requires_colon_after_dot():
    # A real version must never be dropped; the rule only fires on a dangling
    # dot immediately before the colon.
    s = "NM_001754.5:c.749G>A"
    cleaned, fixes = clean_hgvs(s, validate=False)
    assert cleaned == s
    assert not fixes


def test_space_before_kind_leaves_accessions_to_whitespace_strip():
    # An accession head (contains '_'/'.') is glued by the whitespace strip and
    # then repaired by reconstruction. The space-separator rule must not fire,
    # so the reported fixes stay the same as before the rule existed.
    cleaned, fixes = clean_hgvs("nm_000059.4 c.316+5G>A")
    assert cleaned == "NM_000059.4:c.316+5G>A"
    assert C.STRIPPED_WHITESPACE in codes(fixes)
    assert C.FIXED_SEPARATOR_TYPO not in codes(fixes)


CONTIG_NAME_CASES = [
    # #148: a chromosome name before its own accession is redundant and dropped,
    # for each supported build, unversioned, upper case, and chrM/chrMT.
    ("chr13 NC_000013.11:g.32316528del",         "NC_000013.11:g.32316528del"),
    ("chr13 NC_000013.10:g.32890665del",         "NC_000013.10:g.32890665del"),
    ("CHR13 NC_000013:g.32316528del",            "NC_000013:g.32316528del"),
    ("chrX NC_000023.11:g.73820651A>G",          "NC_000023.11:g.73820651A>G"),
    ("chr1 NC_060925.1:g.1000A>G",               "NC_060925.1:g.1000A>G"),
    ("chrM NC_012920.1:m.3243A>G",               "NC_012920.1:m.3243A>G"),
    ("chrMT NC_012920.1:m.3243A>G",              "NC_012920.1:m.3243A>G"),
    ("GRCh38 chr13 NC_000013.11:g.32316528del",  "NC_000013.11:g.32316528del"),
]


@pytest.mark.parametrize("bad,expected", CONTIG_NAME_CASES)
def test_strip_matching_contig_name(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.STRIPPED_CONTIG_NAME in codes(fixes)
    assert C.WRAPPED_GENE_SYMBOL not in codes(fixes)


@pytest.mark.parametrize("bad", [
    "chr7 NC_000013.11:g.32316528del",  # contig disagrees with the accession
    "chr13 NM_000059.4:c.68del",        # can't check a transcript's contig
])
def test_contig_name_not_stripped_or_wrapped_unless_it_matches(bad):
    # Strict: never guess which of two disagreeing references was meant, and
    # never treat a chromosome name as a gene symbol.
    cleaned, fixes = clean_hgvs(bad)
    assert C.STRIPPED_CONTIG_NAME not in codes(fixes)
    assert C.WRAPPED_GENE_SYMBOL not in codes(fixes)
    assert not cleaned.startswith("NC_000013.11:")
    assert "(chr" not in cleaned


@pytest.mark.parametrize("bad,expected", [
    # Many genes start with CHR; they must still be treated as gene symbols
    ("CHRNA7 NM_000746.6:c.1A>G",                "NM_000746.6(CHRNA7):c.1A>G"),
    ("CHRM1 NM_000738.3:c.1A>G",                 "NM_000738.3(CHRM1):c.1A>G"),
])
def test_chr_gene_symbol_before_accession(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.WRAPPED_GENE_SYMBOL in codes(fixes)
    assert C.STRIPPED_CONTIG_NAME not in codes(fixes)


@pytest.mark.parametrize("valid", [
    "NC_000013.11:g.32316528del",
    "NM_000059.4:c.68del",
    "NM_000059.4(BRCA2):c.68del",
    "NM_000746.6(CHRNA7):c.1A>G",
    "NC_012920.1:m.3243A>G",
])
def test_new_ops_leave_valid_hgvs_unchanged(valid):
    cleaned, fixes = clean_hgvs(valid)
    assert cleaned == valid
    assert not fixes


# ---------------------------------------------------------------------------
# #165: a genome build name before the reference is dropped
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("build", [
    "hg19", "hg38", "hg37", "HG19",
    "GRCh37", "GRCh38", "grch38", "GRCh38.p14",
    "T2T", "CHM13", "T2T-CHM13v2.0", "CHM13v2.0",
])
def test_build_name_before_accession_stripped(build):
    cleaned, fixes = clean_hgvs(f"{build} NM_000059.4:c.68del")
    assert cleaned == "NM_000059.4:c.68del"
    assert C.STRIPPED_BUILD_NAME in codes(fixes)
    assert C.WRAPPED_GENE_SYMBOL not in codes(fixes)
    assert C.STRIPPED_LEADING_JUNK not in codes(fixes)


@pytest.mark.parametrize("bad,expected", [
    # Build name, then a gene symbol, which is wrapped as usual
    ("GRCh38 BRCA2 NM_000059.4:c.68del",         "NM_000059.4(BRCA2):c.68del"),
    ("hg19 BRCA2 c.68del",                       "BRCA2:c.68del"),
    # Colon or '#' as the separator
    ("hg38:NM_000059.4:c.68del",                 "NM_000059.4:c.68del"),
    ("GRCh37.p12#NM_000059.4:c.68del",           "NM_000059.4:c.68del"),
    # Build name, then a chromosome name, then the chromosome's accession
    ("T2T chr1 NC_060925.1:g.1000A>G",           "NC_060925.1:g.1000A>G"),
    ("GRCh38 chr13 NC_000013.11:g.32316528del",  "NC_000013.11:g.32316528del"),
])
def test_build_name_before_symbol_or_contig(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.STRIPPED_BUILD_NAME in codes(fixes)
    assert E not in severities(fixes)


@pytest.mark.parametrize("symbol", [
    "CHM",     # choroideremia: a real gene symbol that starts like CHM13
    "HGS",     # a real gene symbol that starts like hg19/hg38
    "HG19X",   # a build name is only recognised as a whole token
])
def test_build_name_lookalike_gene_symbols_kept(symbol):
    cleaned, fixes = clean_hgvs(f"{symbol} NM_000059.4:c.68del")
    assert cleaned == f"NM_000059.4({symbol}):c.68del"
    assert C.WRAPPED_GENE_SYMBOL in codes(fixes)
    assert C.STRIPPED_BUILD_NAME not in codes(fixes)


def test_build_name_alone_is_not_stripped():
    # The op never empties the string: a bare build name is left for validation
    cleaned, fixes = clean_hgvs("hg19")
    assert cleaned == "hg19"
    assert C.STRIPPED_BUILD_NAME not in codes(fixes)
    assert C.NO_COLON in codes(fixes)


def test_symbol_only_input_skips_symbol_before_accession():
    # #148: with no accession after the symbol, the space-before-kind rule fixes it
    cleaned, fixes = clean_hgvs("BRCA2 c.68del")
    assert cleaned == "BRCA2:c.68del"
    assert C.FIXED_SEPARATOR_TYPO in codes(fixes)
    assert C.WRAPPED_GENE_SYMBOL not in codes(fixes)


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
    # #152: every kind letter and the genomic/protein accession types
    "NC_012920.1:m.3243A>G",
    "NM_000059.4:r.68a>g",
    "NM_000059.4:r.68del",
    "NG_012337.3:g.5000A>G",
    "NP_000050.2:p.Arg1Gly",
    "NP_000050.2:p.(Glu23fs)",
    "NR_003051.3:n.68del",
    # Gene symbols that start with a two-letter RefSeq prefix (AC, NS, ...)
    "ACTB:c.100A>G",
    "NSD1:c.100A>G",
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
    assert not fixes, f"Unexpected errors for valid HGVS {hgvs_string!r}: {[f.code.value for f in fixes]}"


@pytest.mark.filterwarnings("ignore::cdot.hgvs.clean.HGVSCleanOpWarning")
@pytest.mark.parametrize("op", list(HGVSCleanOp))
@pytest.mark.parametrize("hgvs_string", ALREADY_VALID)
def test_already_valid_is_fixed_point_of_every_op(hgvs_string, op):
    # #150: every valid string must be left alone by each op on its own, not only
    # by the whole pipeline (where a later op could undo an earlier one).
    cleaned, fixes = clean_hgvs(hgvs_string, ops={op}, validate=False)
    assert cleaned == hgvs_string, f"{op.name} modified valid HGVS {hgvs_string!r} → {cleaned!r}"
    assert not fixes


def test_clean_op_enum_matches_pipeline_order():
    # HGVSCleanOp documents that its members are listed in pipeline order
    from cdot.hgvs.clean import _PIPELINE
    assert [step.op for step in _PIPELINE] == list(HGVSCleanOp)


def test_pipeline_prerequisites_run_earlier():
    # An op can only rely on output from an op that runs before it
    from cdot.hgvs.clean import _PIPELINE
    earlier = set()
    for step in _PIPELINE:
        assert step.requires <= earlier, f"{step.op.name} requires a later op"
        earlier.add(step.op)


def test_missing_prerequisite_warns():
    # #166: without the reconstruction that inserts the colon, the mutation type
    # in a colon-less string is not lowercased
    s = "NM_000059.4 c.68_69DELINSAT"
    with pytest.warns(HGVSCleanOpWarning, match="LOWERCASE_MUTATION_TYPE.*RECONSTRUCT_STRUCTURE"):
        cleaned, _fixes = clean_hgvs(s, ops={HGVSCleanOp.STRIP_WHITESPACE,
                                             HGVSCleanOp.LOWERCASE_MUTATION_TYPE}, validate=False)
    assert cleaned == "NM_000059.4c.68_69DELINSAT"


@pytest.mark.parametrize("ops", [
    None,
    ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_PROTEIN_SUFFIX},
    {HGVSCleanOp.STRIP_WHITESPACE},
    {HGVSCleanOp.RECONSTRUCT_STRUCTURE, HGVSCleanOp.LOWERCASE_MUTATION_TYPE,
     HGVSCleanOp.ADD_N_PREFIX, HGVSCleanOp.ADD_TRANSCRIPT_UNDERSCORE},
])
def test_prerequisites_selected_no_warning(ops):
    with warnings.catch_warnings():
        warnings.simplefilter("error", HGVSCleanOpWarning)
        clean_hgvs("NM_000059.4 c.68_69DELINSAT", ops=ops)


# ---------------------------------------------------------------------------
# #150: reconstruction must not build garbage when the colon is missing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [
    # Missing colon and kind: the version digits ran into the coordinate.
    "NM_000059.4 123A>G",
    # Missing colon and kind after a gene symbol in parens
    "NM_000059.4(BRCA2) 100G>A",
    # Last letter of a gene symbol (FANCG) followed by a coordinate is not a kind
    "NM_000059.4 FANCG 100A>G",
])
def test_reconstruct_refuses_ambiguous_missing_colon(bad):
    # These can't be repaired unambiguously: must be reported as an ERROR, never
    # rewritten into a "valid looking" string.
    cleaned, fixes = clean_hgvs(bad)
    assert C.RECONSTRUCTED_STRUCTURE not in codes(fixes)
    assert C.NO_COLON in codes(fixes)
    assert E in severities(fixes)
    assert ":" not in cleaned
    assert cleaned == bad.replace(" ", "")


@pytest.mark.parametrize("bad,expected", [
    # Kind letter glued to the version with no colon and no dot
    ("NM_005629.3c1403A>C",          "NM_005629.3:c.1403A>C"),
    ("NM_000059.4c100G>A",           "NM_000059.4:c.100G>A"),
    # Gene symbol glued between accession and kind
    ("NM_000059.4 PALB2 c.100A>G",   "NM_000059.4(PALB2):c.100A>G"),
    # Gene symbols starting with RefSeq prefix letters reconstruct as symbols
    ("ACTBc.100A>G",                 "ACTB:c.100A>G"),
    # Uppercase kind letter with its dot
    ("NM_000059.4:C.100A>G",         "NM_000059.4:c.100A>G"),
])
def test_reconstruct_unambiguous_missing_colon(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.RECONSTRUCTED_STRUCTURE in codes(fixes)
    assert E not in severities(fixes)


# ---------------------------------------------------------------------------
# #151: ops that need input a later op used to create
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad,expected,required_codes", [
    # Mutation type lowercased after reconstruction has inserted the colon
    ("NM_000059.4 c.68_69DELINSAT",       "NM_000059.4:c.68_69delinsAT",
     {C.RECONSTRUCTED_STRUCTURE, C.LOWERCASED_MUTATION_TYPE}),
    # Quotes stripped before the start-anchored symbol-before-accession op
    ('"BRCA2 NM_000059.4:c.68del"',       "NM_000059.4(BRCA2):c.68del",
     {C.STRIPPED_SURROUNDING_PUNCTUATION, C.WRAPPED_GENE_SYMBOL}),
    ("'#NM_000059.4:c.68del'",            "NM_000059.4:c.68del",
     {C.STRIPPED_SURROUNDING_PUNCTUATION, C.STRIPPED_LEADING_JUNK}),
    # Underscore restored before the separator-typo op looks for "NM_"
    ("NM000059.4;c.68del",                "NM_000059.4:c.68del",
     {C.ADDED_TRANSCRIPT_UNDERSCORE, C.FIXED_SEPARATOR_TYPO}),
    # Leading "N" restored before the empty-version op looks for "NM_"
    ("M_001754.:c.749G>A",                "NM_001754:c.749G>A",
     {C.ADDED_N_PREFIX, C.DROPPED_EMPTY_VERSION}),
    # A tab (or any whitespace) before the protein annotation
    ("NM_000059.4:c.68del\tp.(Glu23fs)",  "NM_000059.4:c.68del",
     {C.STRIPPED_PROTEIN_SUFFIX}),
])
def test_pipeline_ordering(bad, expected, required_codes):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert required_codes <= codes(fixes), codes(fixes)
    assert E not in severities(fixes)


def test_protein_allele_not_stripped_as_suffix():
    # "p." is the allele itself here, not a trailing annotation to drop
    cleaned, fixes = clean_hgvs("NP_000050.2 p.Arg1Gly")
    assert cleaned == "NP_000050.2:p.Arg1Gly"
    assert C.STRIPPED_PROTEIN_SUFFIX not in codes(fixes)


def test_surrounding_whitespace_reported_as_whitespace():
    # The punctuation op runs first but leaves plain whitespace to its own step
    cleaned, fixes = clean_hgvs("NM_000059.4:c.68del ")
    assert cleaned == "NM_000059.4:c.68del"
    assert codes(fixes) == {C.STRIPPED_WHITESPACE}


# ---------------------------------------------------------------------------
# #152: lowercase genomic/protein accessions and m./r. kinds
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad,expected", [
    ("ng_012337.3:g.5000A>G",   "NG_012337.3:g.5000A>G"),
    ("np_000050.2:p.Arg1Gly",   "NP_000050.2:p.Arg1Gly"),
    ("nc_012920.1:m.3243A>G",   "NC_012920.1:m.3243A>G"),
    # Not matched by reconstruction (dots in the allele), so the gene/transcript
    # step has to uppercase the genomic prefix itself
    ("nc_000023.11:g.(31180435_31200854)_(33274278_33357726)del",
     "NC_000023.11:g.(31180435_31200854)_(33274278_33357726)del"),
    # Only the LRG prefix is uppercased, the "t1" suffix stays lowercase
    ("lrg_199t1:c.1415T>C",     "LRG_199t1:c.1415T>C"),
])
def test_lowercase_accession_any_type(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert C.MISSING_REFERENCE_SEQUENCE not in codes(fixes)
    assert E not in severities(fixes)


def test_fix_gene_transcript_uppercases_genomic_prefix():
    # The op on its own (reconstruction not run) uppercases a genomic prefix
    cleaned, fixes = clean_hgvs("nc_000023.11:g.(1_2)_(3_4)del",
                                ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT}, validate=False)
    assert cleaned == "NC_000023.11:g.(1_2)_(3_4)del"
    assert C.UPPERCASED_HGVS_PREFIX in codes(fixes)


@pytest.mark.parametrize("bad,expected", [
    ("NC_012920.1 m.3243A>G",   "NC_012920.1:m.3243A>G"),
    ("NC_012920.1m.3243A>G",    "NC_012920.1:m.3243A>G"),
    ("NM_000059.4 r.68del",     "NM_000059.4:r.68del"),
    # r. bases are lowercase in HGVS and must stay lowercase
    ("nm_000059.4 r.68a>g",     "NM_000059.4:r.68a>g"),
    ("NP_000050.2 p.Arg1Gly",   "NP_000050.2:p.Arg1Gly"),
    ("NM_000059.4:r,68del",     "NM_000059.4:r.68del"),
    ("NM_000059.4:r.r.68del",   "NM_000059.4:r.68del"),
])
def test_m_r_p_kinds(bad, expected):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == expected
    assert E not in severities(fixes)


@pytest.mark.parametrize("not_missing", [
    "ng_012337.3:g.5000A>G",
    "np_000050.2:p.Arg1Gly",
    "nc_012920.1:m.3243A>G",
    "gata2:c.1A>G",
])
def test_validate_lowercase_reference_is_not_missing(not_missing):
    # Validation alone: starting with a kind letter is not the same as having no
    # reference sequence
    _cleaned, fixes = clean_hgvs(not_missing, ops=set())
    assert C.MISSING_REFERENCE_SEQUENCE not in codes(fixes)


def test_prefix_and_kind_tables_drive_patterns():
    # Every RefSeq prefix in the table is handled by the derived patterns
    from cdot.hgvs.clean import _REFSEQ_PREFIX_TYPE
    for prefix in _REFSEQ_PREFIX_TYPE:
        lower = prefix.lower() + "000123.1"
        cleaned, _fixes = clean_hgvs(f"{lower}:g.1A>G", ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT},
                                     validate=False)
        assert cleaned.startswith(prefix), prefix


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


# ---------------------------------------------------------------------------
# UPPERCASED_HGVS_PREFIX: a lowercase non-transcript accession (eg NG_, NP_) is
# uppercased in either position of the prefix, and is reported distinctly from a
# lowercase transcript (UPPERCASED_TRANSCRIPT)
# ---------------------------------------------------------------------------

def test_uppercased_hgvs_prefix_in_gene_parens():
    # Lowercase genomic accession where the gene symbol goes: not a transcript, so
    # not a gene/transcript swap, just uppercased in place
    cleaned, fixes = clean_hgvs("BRCA2(ng_012337.3):g.5000A>G",
                                ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT}, validate=False)
    assert cleaned == "BRCA2(NG_012337.3):g.5000A>G"
    assert fixes == [HGVSFix(W, C.UPPERCASED_HGVS_PREFIX, "Upper cased HGVS prefix",
                             original="BRCA2(ng_012337.3)", fixed="BRCA2(NG_012337.3)")]


def test_uppercased_hgvs_prefix_in_gene_parens_full_pipeline():
    cleaned, fixes = clean_hgvs("BRCA2(ng_012337.3):g.5000A>G")
    assert cleaned == "BRCA2(NG_012337.3):g.5000A>G"
    assert codes(fixes) == {C.UPPERCASED_HGVS_PREFIX}


@pytest.mark.parametrize("bad,expected,original,fixed", [
    ("ng_012337.3:g.5000A>G", "NG_012337.3:g.5000A>G", "ng_012337.3", "NG_012337.3"),
    ("np_000050.2:p.Arg1Gly", "NP_000050.2:p.Arg1Gly", "np_000050.2", "NP_000050.2"),
    ("nc_012920.1:m.3243A>G", "NC_012920.1:m.3243A>G", "nc_012920.1", "NC_012920.1"),
])
def test_uppercased_hgvs_prefix_accession_position(bad, expected, original, fixed):
    cleaned, fixes = clean_hgvs(bad, ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT}, validate=False)
    assert cleaned == expected
    assert fixes == [HGVSFix(W, C.UPPERCASED_HGVS_PREFIX, "Upper cased HGVS prefix",
                             original=original, fixed=fixed)]


def test_lowercase_transcript_is_not_reported_as_hgvs_prefix():
    # A lowercase mRNA accession gets the more specific UPPERCASED_TRANSCRIPT code
    cleaned, fixes = clean_hgvs("nm_000059.4:c.68del",
                                ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT}, validate=False)
    assert cleaned == "NM_000059.4:c.68del"
    assert codes(fixes) == {C.UPPERCASED_TRANSCRIPT}


def test_lowercase_transcript_in_gene_parens_is_swap_not_hgvs_prefix():
    cleaned, fixes = clean_hgvs("BRCA2(nm_000059.4):c.68del",
                                ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT}, validate=False)
    assert cleaned == "NM_000059.4(BRCA2):c.68del"
    assert codes(fixes) == {C.SWAPPED_GENE_TRANSCRIPT, C.UPPERCASED_TRANSCRIPT}


def test_uppercase_prefix_already_uppercase_no_fix():
    cleaned, fixes = clean_hgvs("NG_012337.3:g.5000A>G",
                                ops={HGVSCleanOp.FIX_GENE_TRANSCRIPT}, validate=False)
    assert cleaned == "NG_012337.3:g.5000A>G"
    assert fixes == []


# ---------------------------------------------------------------------------
# NO_COLON: validation of a string with no accession/allele separator
# ---------------------------------------------------------------------------

def test_no_colon_is_the_only_error_for_bare_accession():
    cleaned, fixes = clean_hgvs("NM_000059.4")
    assert cleaned == "NM_000059.4"
    assert fixes == [HGVSFix(E, C.NO_COLON, "No colon ':' found in HGVS string")]
    assert error_messages(fixes) == ["No colon ':' found in HGVS string"]


def test_no_colon_not_reported_when_validation_off():
    cleaned, fixes = clean_hgvs("NM_000059.4", validate=False)
    assert cleaned == "NM_000059.4"
    assert fixes == []


def test_no_colon_absent_for_well_formed_input():
    _cleaned, fixes = clean_hgvs("NM_000059.4:c.68del")
    assert C.NO_COLON not in codes(fixes)
    # and after a repair that inserts the colon
    _cleaned, fixes = clean_hgvs("NM_000059.4 c.68del")
    assert C.NO_COLON not in codes(fixes)
    assert E not in severities(fixes)


def test_no_colon_reported_alongside_other_warnings():
    # Whitespace is cleaned (WARNING) but the missing colon is still an ERROR
    cleaned, fixes = clean_hgvs("NM_000059.4 ")
    assert cleaned == "NM_000059.4"
    assert codes(fixes) == {C.STRIPPED_WHITESPACE, C.NO_COLON}
    assert severities(fixes) == {W, E}


# ---------------------------------------------------------------------------
# STRIPPED_LEADING_JUNK: a tab character, and the two-character literal "\t"
# left by a copy from an escaped export, before the accession
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [
    "\tNM_000059.4:c.68del",            # a real tab
    "\\tNM_000059.4:c.68del",           # the two characters backslash, t
    "\\t\\tNM_000059.4:c.68del",        # repeated
    "GRCh38\\tNM_000059.4:c.68del",     # build name joined to the accession by a literal tab
    "#\\tNM_000059.4:c.68del",          # mixed with other junk
])
def test_literal_and_real_tab_before_accession_is_leading_junk(bad):
    cleaned, fixes = clean_hgvs(bad)
    assert cleaned == "NM_000059.4:c.68del"
    assert C.STRIPPED_LEADING_JUNK in codes(fixes)
    assert E not in severities(fixes)
    junk = [f for f in fixes if f.code == C.STRIPPED_LEADING_JUNK][0]
    assert junk.original == bad
    assert junk.fixed == "NM_000059.4:c.68del"


def test_literal_tab_needs_the_leading_junk_op():
    # Without the op the literal "\t" is not whitespace, so nothing else removes it
    cleaned, fixes = clean_hgvs("\\tNM_000059.4:c.68del",
                                ops=ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_LEADING_JUNK})
    assert cleaned == "\\tNM_000059.4:c.68del"
    assert C.STRIPPED_LEADING_JUNK not in codes(fixes)


def test_literal_tab_not_followed_by_accession_is_kept():
    # The op only fires when an accession follows the junk
    cleaned, fixes = clean_hgvs("\\tBRCA2:c.68del", validate=False)
    assert C.STRIPPED_LEADING_JUNK not in codes(fixes)
    assert cleaned.endswith("BRCA2:c.68del")


@pytest.mark.parametrize("s", [
    "NM_000059.4:c.68del",
    "NM_000059.4(BRCA2):c.68del",
    "NM_000059.4(BRCA2)\n:c.68del",
    "BRCA2(NM_000059.4):68del",
    "NM_000059.4(BRCA2)(X):C.68DEL",
    "(BRCA2):c.(100_200)del",
    "NM_000059.4():c.68del",
])
def test_parse_round_trips(s):
    # #166: the structural ops rely on parsing being lossless, so one that
    # changes nothing returns its input
    from cdot.hgvs.clean import _parse
    assert str(_parse(s)) == s


@pytest.mark.parametrize("hgvs_string, message", [
    ("NM_000059.4:g.32316528del", "NM_000059.4 is a coding transcript accession, so it can't take 'g.' (use c. or r.)"),
    ("NM_000059.4:p.Arg10Gly", "NM_000059.4 is a coding transcript accession, so it can't take 'p.' (use c. or r.)"),
    ("NM_000059.4:n.68del", "NM_000059.4 is a coding transcript accession, so it can't take 'n.' (use c. or r.): "
                            "n. counts from the transcript start, not the start codon"),
    ("NR_003051.3:c.10del", "NR_003051.3 is a non-coding transcript accession, so it can't take 'c.' (use n. or r.)"),
    ("NC_000013.11:c.68del", "NC_000013.11 is a genomic accession, so it can't take 'c.' (use g. or m.)"),
    ("NP_000050.2:c.68del", "NP_000050.2 is a protein accession, so it can't take 'c.' (use p.)"),
    ("ENST00000380152.8:g.32316528del", "ENST00000380152.8 is a transcript accession, so it can't take 'g.' (use c., n. or r.)"),
    ("NM_000059.4(BRCA2):g.32316528del", "NM_000059.4 is a coding transcript accession, so it can't take 'g.' (use c. or r.)"),
])
def test_kind_reference_mismatch(hgvs_string, message):
    cleaned, fixes = clean_hgvs(hgvs_string)
    assert cleaned == hgvs_string
    assert [(f.severity, f.code, f.message) for f in fixes] == [(E, C.KIND_REFERENCE_MISMATCH, message)]


@pytest.mark.parametrize("hgvs_string", [
    "NM_000059.4:c.68del", "NM_000059.4:r.68del", "NR_003051.3:n.10del",
    "NC_000013.11:g.32316528del", "NC_012920.1:m.3243A>G", "NP_000050.2:p.Arg10Gly",
    "ENST00000380152.8:n.68del", "BRCA2:c.68del", "LRG_199t1:c.100del",
])
def test_kind_reference_match(hgvs_string):
    _cleaned, fixes = clean_hgvs(hgvs_string)
    assert C.KIND_REFERENCE_MISMATCH not in codes(fixes)


@pytest.mark.parametrize("hgvs_string", ["BRCA2:g.32316528del", "brca2:g.32316528del", "MT-ND1:m.3460G>A"])
def test_gene_symbol_with_genomic_kind(hgvs_string):
    gene, kind = hgvs_string.split(":")[0], hgvs_string.split(":")[1][0]
    _cleaned, fixes = clean_hgvs(hgvs_string)
    assert [(f.severity, f.code, f.message) for f in fixes] == [(
        E, C.KIND_REFERENCE_MISMATCH,
        f"{gene} looks like a gene symbol, so it can't take '{kind}.': {kind}. positions are on a "
        f"genomic accession (eg NC_000013.11)")]


@pytest.mark.parametrize("hgvs_string", ["chr13:g.32316528del", "13:g.32316528del", "X:g.100del", "chrM:m.3460G>A",
                                         "MT:m.3460G>A"])
def test_chromosome_name_with_genomic_kind_not_a_gene_symbol(hgvs_string):
    _cleaned, fixes = clean_hgvs(hgvs_string)
    assert C.KIND_REFERENCE_MISMATCH not in codes(fixes)


@pytest.mark.parametrize("hgvs_string, message", [
    ("chr13:c.68del", "chr13 is a chromosome name, so it can't take 'c.': c. positions are on a transcript "
                      "accession (eg NM_000059.4)"),
    ("13:n.68del", "13 is a chromosome name, so it can't take 'n.': n. positions are on a transcript "
                   "accession (eg NM_000059.4)"),
    ("X:p.Arg10Gly", "X is a chromosome name, so it can't take 'p.': p. positions are on a protein "
                     "accession (eg NP_000050.3)"),
])
def test_chromosome_name_with_transcript_kind(hgvs_string, message):
    _cleaned, fixes = clean_hgvs(hgvs_string)
    assert [(f.severity, f.code, f.message) for f in fixes] == [(E, C.KIND_REFERENCE_MISMATCH, message)]
