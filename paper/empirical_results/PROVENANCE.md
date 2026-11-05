# Stored results: provenance

This directory is the checked-in snapshot of measured results for the paper. Each
CSV is a one-row "facts" file (one column per fact) that the render reads via
`output/facts/` (see `paper/README.md` for the two build modes).

Two kinds of file live here:

* **Computed** results have a Snakefile rule that regenerates them from data or a
  benchmark run. When the required inputs are absent the rule copies the committed
  file from here instead, so a build always renders. See the rule docstrings in
  `paper/Snakefile`.
* **Frozen** results have no compute path: they were measured once (or transcribed
  from published papers) and committed. The Snakefile's `frozen_fact` rule just
  copies them into the facts dir. Their provenance is recorded below, because there
  is no script or rule to read it off.

Refreshing a frozen result means re-running the command recorded below by hand and
editing the CSV here, then updating the note in this file (date, data version, what
changed).

---

## literature.csv (frozen)

Static constants transcribed from published papers, not measured here.

* `mutalyzer_correct_pct` / `mutalyzer_error_pct` / `mutalyzer_autocorrect_pct`:
  Lefter 2021 (Mutalyzer), the breakdown of the ~26 million unique descriptions
  from five years of production logs (per unique description, not per submission).
  The error figure is 15.8% syntactic plus 25.4% semantic failures of the Name
  Checker.
* `clinvar_variants`, `lof_agreement_pct`, `brca2_accuracy_pct`, `uta_count`,
  `uta_remote_tps`, `mane_coverage_pct`, `seqrepo_speedup_fold`: see the citations
  in `paper/references.bib` and `paper/literature/literature_review.md`.

## clinvar.csv (frozen)

R2 cdot-vs-UTA ClinVar resolution facts, with a RefSeq-vs-Ensembl breakdown. Two
things are recorded here, deliberately from two samples:

1. Headline rates over the seeded 2,818-variant ClinVar sample (cdot 99.0% vs UTA
   70.3%). Frozen constants from the original cdot+UTA run; the seeded sample is not
   committed (built from a ClinVar download), so these stay frozen.
2. Per-source split, measured on the *committed* `tests/test_data/clinvar_hgvs/*.tsv`
   (500 RefSeq NM_ + 818 Ensembl ENST) so a referee can reproduce it. It shows WHERE
   the gap is: RefSeq is at parity (both hold every cited current version), the whole
   gap is Ensembl (UTA holds no ENST at all).

cdot sequences are served by cdot's own `FastaSeqFetcher` (transcript sequence spliced
from a local genome FASTA), so no pair is dropped for a missing sequence: this covers
the Ensembl transcripts SeqRepo does not hold, and reproduces Ensembl transcript
sequence exactly. The whole committed set is counted (500 RefSeq + 818 Ensembl).

How the per-source columns were measured (GRCh38, cdot 0.2.34 / UTA uta_20241220; re-run
2026-11-04 on PyPI cdot 0.2.36 with data_v0.2.36, identical rates):

```bash
# cdot end-to-end resolution over the committed pairs, --fasta = FastaSeqFetcher:
F=GCF_000001405.39_GRCh38.p13_genomic.fna.gz   # NCBI GRCh38, NC_ contig names
python paper/scripts/benchmark_resolution.py tests/test_data/clinvar_hgvs/clinvar_hgvs_500.tsv \
    --json cdot-0.2.34.refseq.GRCh38.json.gz --fasta $F --out cdot_rs.json
python paper/scripts/benchmark_resolution.py tests/test_data/clinvar_hgvs/clinvar_hgvs_ensembl.tsv \
    --json cdot-0.2.34.ensembl.GRCh38.json.gz --fasta $F --out cdot_en.json
```

cdot: RefSeq 500/500 = 100.0%; Ensembl 816/818 = 99.8% (the 2 misses are Ensembl
patch-scaffold-only transcripts absent from the data, issue #113, not sequence).
Rates identical on 0.2.33 and 0.2.34.

UTA coverage is membership of the exact cited accession.version in UTA's GRCh38
alignment set (biocommons has no version fallback, so holding the exact version is
what UTA resolution requires). Against the UTA GRCh38 alignment dump (`ac` column of
`Homo_sapiens_GRCh38_UTA_20241220.csv`): RefSeq 500/500 = 100.0%, Ensembl 0 (no ENST
rows in UTA at all) = 0.0%. The live UTA DB gives the same rates; the offline
cdot-JSON conversion of UTA is NOT used because its splign projection under-counts.

## clinvar_submitted.csv, clinvar_submitted_residual.csv (frozen)

R2 submitted-HGVS ClinVar corpus facts: version-age distribution plus the cdot-vs-UTA
comparison on strings as the submitting labs wrote them.

The variant_summary-based corpora take c.HGVS from ClinVar's recomputed `Name` column,
always at the current transcript version, so they are a current-version ceiling. This
corpus instead uses the per-SCV HGVS attributes of the VCV XML (submitted strings
verbatim), joined to the ClinVar VCF coordinate via AlleleID, and tags each pair with
the earliest SCV submission date so it can be sampled by era. The corpus (built from a
ClinVar XML download) is not committed, but a 500-pair seed-42 random sample is
(`tests/test_data/clinvar_hgvs/clinvar_submitted_500.tsv`). Nothing reads it for a
reported number, so it was drawn once (2026-08, ClinVarVCVRelease_2026-06) and is not redrawn with each ClinVar
release.

The cdot-vs-UTA head-to-head is run over the **whole** corpus (no sampling) and broken
down by submission-year era (`era_*` columns); a committed 500-pair random sample
supports quick checks. The full run supersedes the earlier two-sample approach (a
time-bucketed and a random draw); `build_clinvar_submitted_pairs.py` still emits those
draws (`--sample-out` / `--time-bucketed-out`) but R2 now reports the full corpus.

Measured on ClinVarVCVRelease_2026-10 with the data_v0.2.36 RefSeq GRCh38 JSON and local
uta_20241220 (FastaSeqFetcher for cdot, SeqRepo for UTA sequence). The cdot pass is cdot
0.2.36 from PyPI, 2026-11-04. The UTA pass, version age and submitter attribution do not
run cdot resolution code and come from the 2026-11-02 run:

```bash
# corpus: 4,082,030 SCV transcript expressions -> 3,245,277 unique (AlleleID, string)
# pairs; 100% RefSeq, 0 ENST; all dated (2008-2026). The XML path scans to the first
# *transcript* HGVS per assertion, capturing SCVs whose first HGVS attribute is
# genomic/protein (the --scv-csv-dir path drops those).
python paper/scripts/build_clinvar_submitted_pairs.py \
    --xml ClinVarVCVRelease_2026-10.xml.gz clinvar_20261025.GRCh38.vcf.gz \
    clinvar_submitted_pairs.GRCh38.tsv

# version age vs the current annotation release (RS_2025_08): 74.65% not-current
python paper/scripts/compute_submitted_version_age.py clinvar_submitted_pairs.GRCh38.tsv \
    --refseq-grch38 cdot-0.2.36.refseq.GRCh38.json.gz --refseq-allbuilds ...

# full-corpus resolution (VCF-coordinate scoring), cdot (1 h 25 min wall, 635 HGVS/s)
# then UTA (3 h 19 min, 272 HGVS/s):
F=GCF_000001405.39_GRCh38.p13_genomic.fna.gz
python paper/scripts/resolve_clinvar_pass.py clinvar_submitted_pairs.GRCh38.tsv \
    --json cdot-0.2.36.refseq.GRCh38.json.gz --fasta $F --with-fixes --out submitted_full_cdot.csv
UTA_DB_URL=postgresql://postgres@127.0.0.1:5433/uta/uta_20241220 HGVS_SEQREPO_DIR=... \
python paper/scripts/resolve_clinvar_pass.py clinvar_submitted_pairs.GRCh38.tsv \
    --uta --out submitted_full_uta.csv
```

Full head-to-head over 3,245,277 pairs: cdot 99.2% vs UTA 81.9% matched; 560,821 resolve
through cdot alone, 971 through UTA alone, and 0 pairs that resolve as written fail
after `fix_hgvs`. Per-era: 2008-2015 cdot 97.9% / UTA 76.2%; 2016-2020 99.3% / 88.9%;
2021-2026 99.2% / 81.1%.

Against the cdot 0.2.35 pass, 6,859 pairs move to `correct`. 5,808 are `incorrect` rows
citing one of the 4,788 GRCh38 versions whose first exon starts after transcript
position 0, which 0.2.35 resolved shifted by the unaligned length (#230). The other 1,051
(mostly `reference_mismatch` errors) are on minus-strand transcripts with alignment gaps,
whose genome-built sequence `FastaSeqFetcher` assembled wrongly in 0.2.35 (fixed in the
same release). 10 further pairs, all on minus-strand gapped transcripts, change between
non-correct buckets or leave `correct` (1).

The one-row CSV is aggregated from the two pass tables, joined row by row on the pairs
file. A pair counts for cdot when its `fixed_bucket` (after `fix_hgvs`) is `correct` and
for UTA when its `bucket` is `correct`; `rescued_by_fix` is cdot `bucket != correct`
with `fixed_bucket == correct`, `regressions` the reverse. Eras use the earliest SCV
submission year (2008-2015, 2016-2020, 2021-2026). The version-age columns come from
`compute_submitted_version_age.py` and `n_scv_tx_strings` from
`submitter_attribution.py` (`n_versioned_citations`). Every percentage is rounded once,
from the counts.

`clinvar_submitted_residual.csv` is derived from the full cdot pass rows with
`fixed_bucket != correct` (25,923 of 3,245,277); the `no_data` bucket is split into
`version_refused` (accession holds other versions, via `get_tx_versions`) vs
`unknown_accession`, `incorrect` maps to `coordinate_drift`, and the `error` bucket is
re-resolved to recover the exception class (`HGVSInvalidVariantError` reference_mismatch,
`HGVSInvalidIntervalError` position_out_of_bounds, `HGVSParseError` grammar_unsupported).

## submitter_attribution.csv (frozen)

R2 submitter attribution of absent-from-cdot transcript versions (why version
substitution is needed at all). `paper/scripts/submitter_attribution.py` streams the
VCV XML, attributes each SCV's first RefSeq transcript citation to its submitting
laboratory, and flags versions cdot does not hold (same test as
`compute_submitted_version_age.py`). Measured 2026-11-02 (ClinVarVCVRelease_2026-10,
data_v0.2.36 RefSeq GRCh38):

```bash
python paper/scripts/submitter_attribution.py --xml ClinVarVCVRelease_2026-10.xml.gz \
    --refseq-grch38 cdot-0.2.36.refseq.GRCh38.json.gz --out submitter_attribution.csv
```

18,155 of 4,082,030 versioned RefSeq citations (0.445%) cite a version absent from cdot.
Of 3,028 submitters only 88 ever do; a single laboratory contributes 94.8% of absent
citations, and 82.9% of the 287 distinct absent versions are cited by one submitter
alone (the self-alignment signature). Only the aggregate CSV is committed; the
per-submitter table (`--top-out`, which names laboratories) is NOT checked in.

## benchmark_fullscale.csv (frozen)

R3 full-scale throughput of a single local-JSON process over the whole current ClinVar
(g.HGVS, c.HGVS) pair set (`clinvar_pairs.GRCh38.tsv` from ClinVar 2026-10, cdot 0.2.36
from PyPI, data_v0.2.36 RefSeq GRCh38, local SeqRepo), measured 2026-11-04 with the
hot-cache protocol: a first pass warms the OS sequence cache and is discarded, the next
is timed. Nothing else ran on the machine during either pass.

```bash
HGVS_SEQREPO_DIR=... SEQREPO_FD_CACHE_MAXSIZE=128 \
python paper/scripts/resolve_clinvar_pass.py clinvar_pairs.GRCh38.tsv \
    --json cdot-0.2.36.refseq.GRCh38.json.gz --out /dev/null   # x2 (discard 1st)
```

Single timed passes vary by several percent between days on this VM (below), so the next
measurement should use `paper/scripts/measure_fullscale.py`, which runs the same pass in
fresh processes, one warm-up plus `--runs` timed (default 3), and writes this CSV with the
mean `hot_tps` and its range (`hot_tps_min`, `hot_tps_max`, `n_runs`). The current values
are from one warm-up and one timed pass.

`n_pairs` 4,583,798; `resolved_pct` 99.4; `cold_tps` 606 HGVS/s (first pass, wall
7570 s); `hot_tps` 603 HGVS/s (timed pass, wall 7603 s), under 1% apart, so full-scale
local-JSON throughput is not sequence-cache-limited. The same protocol on cdot 0.2.35
(2026-11-02) gave 646 / 640 HGVS/s. That gap is run-to-run variation on this VM, not the
client: 0.2.35 and 0.2.36 run alternately over the first 100,000 pairs gave 793 / 790 and
795 / 782 HGVS/s. The controlled backend comparison
(local JSON vs prefetched REST) stays in Table 1 / `benchmark.csv`.

## injection_benchmark.csv, residual_taxonomy.csv, cleaning_fix_categories.csv (table CSVs, not facts)

Multi-row table CSVs rendered inline by the `<!-- include-csv: ... -->` directive in
`supplementary.md` (Tables S5 and S6) and `results.md` (Table 2, frozen with
`cleaning_corpus.csv`), so the numbers live in a CSV rather than hardcoded
in the markdown. `vibepaper`'s fact loader skips multi-row CSVs, so these are not `{{ }}`
facts; edit the CSV to change the table.

* `injection_benchmark.csv` (Table S5): per-category recovery of `clean_hgvs()` vs LOVD /
  VariantValidator / Mutalyzer on the injection benchmark. Body rows from
  `inject_and_clean.py` + the LOVD / VV / Mutalyzer head-to-heads; the **Total** row is
  the frozen aggregate (matches `lovd_comparison.csv` / `vv_mutalyzer_comparison.csv`).
  Refresh by re-running those scripts and re-transcribing.
* `residual_taxonomy.csv` (Table S6): the seven repair-relevant residual classes of the
  private cleaning corpus (frozen LLM classification; see `cleaning_corpus.csv`). Percentages
  are of the 887 genuine-HGVS residual queries.

## cleaning_corpus.csv (frozen, Tier 2)

R4 production cleaning-corpus headline numbers, transcribed from a deterministic run of
`cdot_private/analyze_cleaning.py` (`clean_hgvs()` plus provider-verified accession
restoration) over the private search-box corpus (issue #112). Not reproducible here (the
corpus is private); these are the literal constants R4 / Table 2 / Table S6 render from.

`corpus_n` (32,671) is the loosely-HGVS search-box corpus after removing a small residue
of non-HGVS input (81 queries: pasted URLs, report templates, prose) that slipped the
collection regex and is a data-collection artifact, not something cleaning could repair
(issue #112 feedback). Removing it from the denominator: as-submitted parseable
`as_submitted_pct` 91.7%, after-cleaning `after_pct` 97.3% (`gain_pct` +5.6%, `rescued`
1,828 = `rescued_share_pct` ~67% of the `failed_pp` 8.3 points that failed
as-submitted), residual `residual_n` 887 (`residual_pct` 2.7%); `nonhgvs_n` 81.

Measured 2026-11-04 with PyPI cdot 0.2.36 (`combined_search_hgvs.csv`, 32,752 rows):
`clean_hgvs()` takes 29,956 parseable rows to 31,784 (1,828 rescued, 0 regressed, 968
still failing), and provider-verified accession-prefix restoration
(`resolve_missing_accession_prefix` against the 0.2.36 RefSeq GRCh38 JSON) rescues 1 more.
Three rows of the excluded 81 are non-HGVS strings that cleaning changes: a pasted report
row holding two variants (twice), which still fails to parse, and a sentence quoting one
variant, from which 0.2.36 extracts it. That one is not counted as a rescue. So
`rescued` = 1,828 - 1 + 1 and `residual_n` = 968 - 1 - (81 - 1). cdot 0.2.35 parsed the
report row too, keeping only its first variant, and gives the same `rescued`,
`residual_n` and Table 2 counts once all three are excluded. The same procedure on PyPI cdot 0.2.33 reproduces the
previous snapshot (1,721 rescued, 994 residual) and the previous Table 2 row counts exactly.

Table 2 groups the rescued rows by `HGVSFixCode` (a row counts once per category):
whitespace `STRIPPED_WHITESPACE`; structural punctuation `FIXED_SEPARATOR_TYPO`,
`STRIPPED_SURROUNDING_PUNCTUATION`, `FIXED_MULTIPLE_COLON`, `STRIPPED_UNBALANCED_BRACKETS`,
`FIXED_MULTIPLE_DOT`, `FIXED_MULTIPLE_UNDERSCORE`, `STRIPPED_LEADING_JUNK`,
`STRIPPED_BUILD_NAME`, `REPLACED_UNICODE_DASH`, `DROPPED_EMPTY_VERSION`,
`FIXED_PREFIX_COLON`, `FIXED_MULTIPLE_KIND`; wrapper `SWAPPED_GENE_TRANSCRIPT`,
`FIXED_GENE_WRAPPER`, `WRAPPED_GENE_SYMBOL`, `UPPERCASED_TRANSCRIPT`; reconstruction
`RECONSTRUCTED_STRUCTURE`; protein suffix `STRIPPED_PROTEIN_SUFFIX`; genomic ref in parens
`DROPPED_GENOMIC_REF_IN_PARENS`; base re-casing `UPPERCASED_BASES`; other
`DROPPED_DEL_DUP_COUNT`, `LOWERCASED_MUTATION_TYPE`; prefix/kind `ADDED_N_PREFIX`,
`ADDED_TRANSCRIPT_UNDERSCORE`; plus the one provider-verified restoration. The Table S6
taxonomy is an LLM classification, refreshed by delta re-rating only the residual
queries a new release rescues (Supplementary Table S6). Refresh by re-running the
procedure above and editing this CSV, `cleaning_fix_categories.csv` (Table 2) and
`residual_taxonomy.csv`.

## historical.csv (frozen, Tier 2)

Historical clinical-lab resolution benchmark (cdot vs UTA). The source corpus is the
set of HGVS descriptions imported into the Shariant platform: real, patient-derived
clinical-lab classifications that are NOT shareable, so a referee cannot regenerate
these numbers (hence Tier 2). They were produced by resolving the full corpus against
a local cdot RefSeq+Ensembl JSON and a locally loaded UTA (uta_20241220), identical
biocommons engine, pure coordinate projection. Reproduction (on a machine with the
private corpus):

```bash
UTA_DB_URL=... HGVS_SEQREPO_DIR=... python paper/scripts/benchmark_historical.py \
    shariant_imported_alleles_uniq.txt \
    --json cdot.refseq.grch38.json.gz cdot.ensembl.grch38.json.gz \
    --uta-uri postgresql://uta:uta@localhost:5432/uta/uta_20241220
```

Measured 2026-08-14 (cdot data 0.2.34, refseq+ensembl GRCh38 JSON; re-run 2026-11-04 on
PyPI cdot 0.2.36 with data_v0.2.36, identical): cdot 99.6% vs UTA
93.7% over 36,010 strings; of the 2,141 (5.9%) cdot resolves and UTA does not, 98.6%
are RefSeq historical versions UTA holds no GRCh38 alignment for, 1.3% Ensembl, 0.1%
other. The jump from the previous measurement (97.5% on 0.2.33-era data, 2026-06-18)
is the #51 historical GRCh38 RefSeq alignments plus the #95 UTA-source restore; the
corpus itself also grew slightly (35,695 -> 36,010 lines). The corpus source mix,
per-source rates and the cdot-only cause split come from `benchmark_historical.py`'s
`--out` JSON (`source_mix`, `cdot.by_source`, `uta.by_source`,
`head_to_head.cdot_only_by_cause`); `n_unique_tx` / `n_multi_version` are counted from
the corpus accessions.

## clinvar_vcf.csv, clinvar_vcf_residual.csv (frozen)

Full-ClinVar cdot-only resolution accuracy scored against the ClinVar VCF coordinate
(R1 / Supplementary S4). Built by `paper/scripts/build_clinvar_pairs.py` (joining
variant_summary c.HGVS to the VCF ground truth via AlleleID) and resolved by
`paper/scripts/resolve_clinvar_pass.py` over the full pair set: a dedicated multi-hour
run against a production RefSeq GRCh38 JSON.gz, which is why there is no rule.
Measured 2026-11-04: ClinVar `variant_summary` and `clinvar_20261025.GRCh38.vcf.gz`,
cdot 0.2.36 from PyPI, data_v0.2.36 RefSeq GRCh38, FastaSeqFetcher (1 h 51 min wall).
Comparing `c_to_g` output as a VCF tuple (via biocommons Babelfish) is
representation-robust (3'-shift, equivalent del/dup spellings) and is defined even for
the ClinVar variants whose CLNHGVS uses tandem-repeat `ref[N]` or identity `=`
notation the HGVS parser cannot read.

4,583,798 pairs: 99.40% resolved, and of those 99.959% match the ClinVar VCF
coordinate exactly. `clinvar_vcf_residual.csv` breaks the 1,850 incorrect rows down by
cause with `paper/scripts/categorize_vcf_residual.py` (pairs TSV plus the pass table;
first matching rule wins: identity or symbolic ALT, IUPAC ambiguity code at the same
position, cdot allele over 10 bases longer than ClinVar's, ClinVar SNV, other indel). Analysis writeup: `paper/notes/clinvar_diff_coordinates.md`.
`--examples-out` writes the first 20 rows of each category to
`clinvar_vcf_residual_examples.csv` (public ClinVar data) for inspecting by hand. Of the 68
`ambiguity_code_allele` rows, ClinVar's VCF writes the ALT with N in all of them. 27 cite an
R/Y/K/M/B/V/D/H code in the c.HGVS; the 16 of those on minus-strand transcripts are also
wrong against ClinVar's g.HGVS, because bioutils `reverse_complement` leaves ambiguity
codes uncomplemented (c.315G>Y maps to g.C>Y, not g.C>R; biocommons/bioutils#121). The other 41 are long
insertions holding N runs. Of the 71 `long_cdot_allele` rows, 67 are c. ranges whose two
ends are in different exons (eg `NM_000059.4:c.8331_8332del`, c.8331 being the last base of
an exon), which biocommons HGVS maps to a genomic deletion that includes the intron
(`g.32363533_32370402del`); 2 have intronic offsets, and 2 are ClinVar joins of a
multi-base deletion to an SNV.

## clinvar_residual_positions.csv (frozen)

Where along the CDS the residual `incorrect` projections of that same full-ClinVar VCF
pass sit, binned into deciles (1 = 5'-most tenth) against the `correct` rows' baseline
distribution, with 5'UTR / 3'UTR / non-coding citations counted separately. Companion
to the version-bump positional drift in `positional_drift.csv`. Produced by
`paper/scripts/bin_residual_positions.py` from the per-variant table of the run above
plus the same release JSON.gz; frozen for the same reason (the table is a dedicated
full run).

## clinvar_source_breakdown.csv (frozen)

Per-source resolution outcome breakdown of the full-ClinVar pass (Supplementary S4),
produced by `paper/scripts/summarize_clinvar_pass.py` from the same per-variant table.
Multi-row (one row per transcript source plus a total), unlike the one-row facts files.

## genomic_mismatch.csv (frozen)

Cause breakdown of the residual `incorrect` bucket under the *older* g.HGVS string
comparison, produced by `paper/scripts/categorize_genomic_mismatches.py`. It shows
that most string mismatches are representation differences rather than cdot
coordinate errors (raw 0.513% vs genuine 0.025%). Superseded for the headline numbers
by the VCF-coordinate scoring in `clinvar_vcf.csv`, kept as the supporting analysis.

## version_safety_validation.csv (frozen)

Empirical ClinVar validation of safe transcript-version substitution (R5b): take every
ClinVar variant where another version of its transcript exists and
`is_version_substitution_safe` calls the V->W substitution safe, resolve under W, and
confirm the genomic coordinate is identical to resolving under V. Produced by
`paper/scripts/validate_safe_versions.py` from the full-ClinVar per-variant table
above. A VCF-format pass table is compared as a VCF tuple (Babelfish). Measured
2026-11-04 on the 2026-10 VCF pass (PyPI cdot 0.2.36, data_v0.2.36 RefSeq GRCh38,
FastaSeqFetcher; 3 h 39 min wall):

```bash
python paper/scripts/validate_safe_versions.py vcf_pass_cdot.csv \
    --json cdot-0.2.36.refseq.GRCh38.json.gz --fasta $F \
    --out version_safety_validation.csv --show-changes
```

4,342,420 variants have another version of their transcript; 8,409,757 substitutions are
called safe and tested, 8,403,420 resolve under the substitute (6,337 do not, which is not
a coordinate change), and 0 move the coordinate. The structural Tier-1 claim of
`version_stability.csv` becomes an observed zero error rate. Frozen because it consumes
the same dedicated full-ClinVar run.

---

## Computed results

These have a Snakefile rule; the rule docstring records what regenerates them and what
inputs it needs. The measurement conditions behind the committed snapshot are noted
here.

### coverage.csv

`paper/scripts/compute_coverage.py` over the release JSON.gz files (needs
`--config data_dir=...`). Committed snapshot: data 0.2.36, measured 2026-10-13 (counts
identical to 0.2.34, as no transcript sources were added between them); re-run
2026-11-04 on PyPI cdot 0.2.36, identical.

### benchmark.csv

`paper/scripts/compute_benchmark.py`. Committed snapshot measured 2026-11-04 (cdot
0.2.36 from PyPI, data_v0.2.36 GRCh38 RefSeq, uta_20241220 local, public
uta.biocommons.org remote, cdotlib.org REST, biocommons engine, one shared local SeqRepo
with fd caching), run alone on the machine with `UTA_DB_URL` unset, since
`uta.connect()` with no URL reads that variable and would turn the remote row local:
median (q1, q3) HGVS/s over 5 passes of the committed 500-pair set (remote UTA: first
25 pairs). See R3 / Table 1. The remote UTA row does not depend on the cdot client and
is kept from the 2026-11-03 run (`--skip-uta-remote` on 2026-11-04).

### version_stability.csv, positional_drift.csv

`paper/scripts/compute_version_stability.py`. Committed snapshot: GRCh38,
`--sample 12000 --seed 1`, data 0.2.36 (re-measured 2026-10-13; re-run 2026-11-04 on
PyPI cdot 0.2.36, identical). The only change from
0.2.34 is `refseq_struct_portable_pct` 99.3 -> 99.4; everything else, including
positional_drift.csv, is identical. Relative to the pre-0.2.33 snapshot, the
RefSeq drift increase comes from the #51 historical alignments (old versions drift
more), and the large Ensembl stability improvement from Ensembl 116 entering the data
(verified by re-running on the 0.2.33 files: RefSeq matched the old snapshot exactly,
Ensembl already matched the new numbers).

### cleaning.csv

`paper/scripts/inject_and_clean.py` over public ClinVar c.HGVS committed to this repo.
Fully reproducible, no inputs needed.

### lovd_comparison.csv

`paper/scripts/lovd_head_to_head.py`. Committed snapshot measured 2026-11-04 (PyPI cdot
0.2.36), LOVD HGVS syntax checker v1.2.2 (commit 4cd074ba9cbd), PHP 8.5.

### vv_mutalyzer_comparison.csv

`paper/scripts/vv_mutalyzer_head_to_head.py` against the public VariantValidator and
Mutalyzer services, with responses cached in `output/rest_cache/`. Most responses date
from the 2026-08-18 run (`service_date`); the 218 VariantValidator and 21 Mutalyzer
inputs that were new on 2026-11-02 (the redundant del/dup count cases, which now carry
the true range length) were fetched then and report the same VariantValidator version.
The cdot side was re-scored from that cache on PyPI cdot 0.2.36 (2026-11-04), identical.

### sources.csv

`paper/scripts/compute_sources.py` over `generate_transcript_data/cdot_transcripts.yaml`.
Fully reproducible, no inputs needed.
