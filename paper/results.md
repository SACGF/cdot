# Results

*Full Original Paper (Bioinformatics): Results section. New for the full-paper switch.*
*Numbers come from the two-tier fact model (plan §3): Tier 1 = reproducible templated
facts; Tier 2 =
frozen aggregate constants from the private corpus, flagged non-reproducible.
Synthesised examples
only; never corpus strings.*

> **Provenance flags.** Every reported number is tagged **[Tier 1]** (reproducible from
public data
> committed to this repo) or **[Tier 2]** (aggregate statistics from a private
production corpus,
> published as frozen constants, not reproducible by a referee; see Methods /
data-availability).

---

## R1: Transcript coverage and ClinVar resolution

**[Tier 1]** The unit of coverage is the *transcript-version alignment*: a particular
transcript version aligned to a particular genome build. The same transcript version is
counted separately per build, because each alignment is what a resolution against that
build actually needs. cdot covers {{ coverage.total_count | commas }} such alignments
across all builds and sources, compared with ~{{ literature.uta_count | commas }} in UTA,
a {{ coverage.improvement_fold | fmt('.1f') }}× increase. The gain has two
distinct origins. First, historical depth: UTA does hold several versions per transcript,
but cdot ingests the complete run of RefSeq and Ensembl releases straight from the FTP
archives, so it retains many more historical versions. An NM_ version cited in an older
clinical report or ClinVar submission still resolves even after NCBI has retired it from
the current annotation. Second, Ensembl is covered in
full: {{ coverage.ensembl_unique_count | commas }} Ensembl transcript accessions are
present in cdot but absent from UTA entirely. T2T-CHM13v2.0 adds a further
{{ coverage.t2t_unique_count | commas }} alignments; many transcript versions also exist
on GRCh37/38, but their T2T alignments are distinct and are counted as additional
transcript-version alignments.

**[Tier 1]** To measure practical impact rather than raw counts, a seeded sample of
{{ clinvar.n_variants | commas }} ClinVar [@Landrum2025] variant descriptions,
spanning both RefSeq (NM_) and Ensembl (ENST) transcripts, was resolved against
cdot and a locally loaded UTA (release `uta_20241220`) through the identical
biocommons/hgvs code path (genomic→transcript projection, with sequences served from a
shared local SeqRepo so only the transcript data layer differs). A sample is used here
because this comparison is gated by UTA's throughput; the full ClinVar set is resolved
through cdot alone in R6. cdot resolved
{{ clinvar.cdot_resolution_pct | dp(1) }}% versus
{{ clinvar.uta_resolution_pct | dp(1) }}% for UTA. On the RefSeq subset the two are at
parity (~99% each); the entire difference is Ensembl, which UTA cannot resolve at all
(0% of ENST inputs, as it stores no Ensembl alignments) whereas cdot resolves it
natively (Supplementary Table S4). The benchmark is reproducible: the ClinVar
build script, random seed, and resolution harness are committed.

**[Tier 2: production validation, not reproducible].** The same gap holds on genuinely
historical clinical data, the data that motivated cdot. We resolved the complete set of
{{ historical.n_lines | commas }} unique HGVS descriptions imported into the Australian
Genomics Shariant variant-sharing platform [@Tudini2022]: real classifications submitted
by clinical genetic-testing laboratories over many years, each written against whatever
transcript version was current when the variant was first classified
({{ historical.n_multi_version | commas }} of the {{ historical.n_unique_tx | commas }}
distinct transcripts are cited at more than one version). Run through the identical
biocommons engine with only the transcript-data layer swapped (pure coordinate
projection, `replace_reference=False`, so the sequence layer never differs), cdot
produced a genomic coordinate for {{ historical.cdot_resolved_pct | dp(1) }}% of them
versus {{ historical.uta_resolved_pct | dp(1) }}% for the same locally loaded UTA release
(`uta_20241220`). Of the strings cdot resolved but UTA could not
({{ historical.cdot_only_pct | dp(1) }}% of the corpus),
{{ historical.cdot_only_historical_pct | dp(0) }}% were RefSeq transcript versions for
which UTA holds no GRCh38 alignment and {{ historical.cdot_only_ensembl_pct | dp(0) }}%
were Ensembl transcripts (absent from UTA entirely), confirming, on the real corpus cdot
was built to serve, that the resolution gap is driven by historical-version depth plus
Ensembl coverage. Where the ClinVar comparison above is a sanity check that the
GTF→JSON→biocommons pipeline works at scale (ClinVar cites current transcript versions,
which any up-to-date backend holds), this corpus is the demonstration that the historical
depth actually matters: it is exactly the older-version traffic a working clinical lab
generates. There is no ground-truth genomic coordinate for the private corpus, so the
metric is resolution rate rather than correctness, and because the classifications are
patient-derived the corpus cannot be shared (Methods, data availability).

## R2: String cleaning recovers malformed real-world HGVS

The headline test of cleaning is its effect on a real production query stream; a
reproducible injection benchmark provides a supporting safety check.

**[Tier 2: production validation, not reproducible] Headline result.** The corpus is
**N = 32,752** real queries typed into the HGVS search box of production clinical and
research variant-curation platforms based on VariantGrid [@VariantGrid]. Users treat that box as a shortcut to jump straight
to a variant or its classification, so the strings are whatever a clinician or curator
happened to paste or type, not curated HGVS. They arrive carrying the damage of their
route to the box: stray whitespace and non-printable characters from copying out of Word
documents and report PDFs and pasting between systems, lost casing, transposed
punctuation, and trailing protein annotations. Run over this corpus, `clean_hgvs()` raised
the fraction parseable by biocommons/hgvs from **91.5%** as-submitted to **96.6%** after
cleaning: a **+5.1%** absolute gain (≈1,700 additional strings rescued) with **0
regressions** (no already-valid string was broken). This is the result that matters
operationally: it measures what cleaning recovers from genuinely messy, human-entered
input rather than from synthetic errors. The rescues break down by fix type as shown in
Table 1; they are dominated by whitespace removal and base re-casing, followed by
protein-suffix stripping and gene/transcript-wrapper repair. Because the underlying
queries contain real patient-derived data they cannot be shared or regenerated by a
referee.

**Table 1. Fixes applied across the production corpus (N = 32,752).** Each row is a
`clean_hgvs()` fix category, with the number of rescued queries in which it fired and that
number as a share of the 1,678 rescued queries. Categories overlap (a single query may
need several fixes), so the counts sum to more than the total. *(Tier 2; counts are
frozen constants from a deterministic run of `clean_hgvs()` over the production corpus.)*

| Fix category | Example (→ repaired) | Rescued queries | % of rescued |
|---|---|---|---|
| Whitespace / non-printable removal | `NM_000059.4: c.1A>G` → `NM_000059.4:c.1A>G` | 940 | 56.0% |
| Base re-casing | `NM_000059.4:c.1delg` → `…delG` | 553 | 33.0% |
| Structural-punctuation repair | `NM_000059..4:c.1A>G` → `NM_000059.4:c.1A>G` | 310 | 18.5% |
| Gene/transcript-wrapper repair | `BRCA2(NM_000059.4):c.1A>G` → `NM_000059.4(BRCA2):c.1A>G` | 209 | 12.5% |
| Protein-suffix stripping | `NM_000059.4:c.1A>G p.(Met1?)` → `NM_000059.4:c.1A>G` | 130 | 7.7% |
| Genomic-ref-in-parens removal | `NM_000059.4(NC_000013.11):c.68del` → `NM_000059.4:c.68del` | 57 | 3.4% |
| Other (del/dup count, mutation-type case, …) | `NM_000059.4:c.1_2del2` → `…del` | 21 | 1.3% |
| Prefix / kind restoration | `NM_000059.4:1A>G` → `NM_000059.4:c.1A>G` | 11 | 0.7% |
| **Total unique queries rescued** | | **1,678** | **100%** |

**[Tier 1] Injection safety check.** As a reproducible control,
`paper/scripts/inject_and_clean.py` injects each `clean_hgvs()` fix category into a seeded
sample of clean, parseable ClinVar c.HGVS strings committed to this repo and confirms the
cleaner recovers the canonical target with **{{ cleaning.inject_regressions | int }}
regressions**: no already-valid string is ever broken (Supplementary Table S5). Because
this benchmark injects the very errors it then repairs, its recovery rate is not an
independent measure of real-world performance and is reported in the supplement only; its
purpose here is to demonstrate the no-regression guarantee on which the production result
depends.

## R3: The ceiling of cleaning, a taxonomy of residual errors  *(Table S6)*

**[Tier 2: not reproducible].** The 3.4% of the production corpus (1,118 queries; 860
unique strings) that still fail to parse after cleaning define the ceiling of what
pure string repair can achieve. Each residual string was assigned to one single-label
error class under a fixed decision-tree taxonomy. Of the eight classes, the seven
repair-relevant ones are shown below with synthesised examples; the eighth was non-HGVS
input (81 queries, 7.2%: pasted URLs, report templates, or prose) and is excluded from
the table, as there is nothing in it for cleaning to repair.

**Residual error classes after cleaning** (counts and % of the 1,118 residual
queries; examples synthesised from public BRCA2 `NM_000059.4`). *(Tier 2; frozen
constants from a deterministic run over the production corpus.)*

| Class | Queries | What it is (*example*) |
|---|---|---|
| Truncated | 284 (25.4%) | cut off before a complete variant — `NM_000059.4:c.68_69` (range, no edit) |
| No reference | 277 (24.8%) | a bare variant body, no transcript/gene/accession — `c.68_69delAG` |
| Bad accession | 167 (14.9%) | missing prefix, or misplaced/truncated version — `000059.4:c.68del` (`NM_` prefix dropped) |
| Edit syntax | 113 (10.1%) | malformed or non-standard edit operation — `NM_000059.4:c.68AG>T` (multi-base reference in a substitution) |
| Trailing / concatenated | 85 (7.6%) | extra characters after a complete variant, or several run together — `NM_000059.4:c.68delAG;c.70A>G` |
| Grammar gap | 81 (7.2%) | legitimate HGVS the biocommons grammar rejects — `NM_000059.4:c.(67+1_68-1)_(70+1_71-1)del` (uncertain-range deletion) |
| Insertion (length only) | 30 (2.7%) | an insertion given as a base count, not a sequence — `NM_000059.4:c.68_69ins5` (position and length recoverable; inserted bases not) |

The residual falls into three regimes. Just over half (~50%: Truncated + No reference) is
incomplete or reference-less user input — information the user never supplied, which no
string-level repair can invent; the integer-length insertions (2.7%) belong with these,
since the inserted bases were never given (only the position and length are recoverable, so
the variant is not properly resolvable). About 33% (Bad accession + Edit syntax + Trailing /
concatenated) is in principle fixable and marks the frontier for future cleaning rules.
The remaining ~7% is a grammar gap (valid HGVS the parser rejects rather than the input;
the Grammar gap class), and a further ~7% (excluded above) was non-HGVS strings that should
not be parsed at all. Most of what remains is therefore either out of scope for any repair
or a downstream grammar limitation.

*Method and limitation:* classification was performed by a large language model (Claude
Opus 4, Anthropic; 2026-06-17) applying the shared decision tree to each unique string,
single-label and single-rater; no second-rater adjudication was done, so no inter-rater
agreement (κ) is reported. The taxonomy is version `v1`. Synthesised examples (from public
NM_000059.4 / NM_001754.5) illustrate each class; no corpus string is reproduced.

## R4: Transcript version fallback

**[Tier 1]** The opt-in adjacent-version fallback (Methods) is evaluated by an
end-to-end ablation (`paper/scripts/benchmark_resolution.py`): removing the requested
version from each test variant and resolving through the fallback recovered the correct
genomic coordinate with no false rescues (a false rescue being a substitution that
resolves to a different coordinate). This is a client-layer feature rather than a backend one: biocommons/hgvs
has no adjacent-version fallback regardless of its data provider, so a UTA-backed pipeline
gains nothing here even though UTA itself stores several versions per transcript. cdot's
multi-release depth gives the fallback more versions to choose from, and it is never
applied automatically, preserving exact-version semantics by default. Whether a given
substitution is *coordinate-safe*, and how cdot decides, is the subject of R5.

## R5: Transcript-version stability and safe version fallback

**[Tier 1]** The version fallback (R4) is only safe if substituting an adjacent version
does not move the variant: a coding `c.` position must still map to the same genomic
coordinate. We measured this directly on the released cdot data
(`paper/scripts/compute_version_stability.py`, GRCh38, a seeded
{{ version_stability.sample_n | commas }}-accession sample drawn from accessions cdot
holds at two or more versions, the only accessions where a version bump can be assessed),
exploiting the fact that
cdot stores the genome *alignment* rather than the sequence, so a version bump can only
move a coordinate if it changes that alignment. Across consecutive RefSeq version bumps
({{ version_stability.refseq_pairs | commas }} pairs),
{{ version_stability.refseq_preserving_pct | dp(1) }}% preserved every coding coordinate;
for Ensembl ({{ version_stability.ensembl_pairs | commas }} pairs)
{{ version_stability.ensembl_preserving_pct | dp(1) }}%. Weighting by coding base (the
chance a *random* variant is unaffected), safety is higher still,
{{ version_stability.refseq_pervariant_safety | dp(3) }} (RefSeq) and
{{ version_stability.ensembl_pervariant_safety | dp(3) }} (Ensembl); and the genuinely
dangerous case, a *partial* bump that mis-places some variants but not others, is rare
({{ version_stability.refseq_partial_drift_pct | dp(1) }}% of RefSeq bumps,
{{ version_stability.ensembl_partial_drift_pct | dp(1) }}% Ensembl); when a coordinate
does move it is almost always the whole CDS, driven by a re-annotation of the coding
region. Drift is also highly concentrated: only
{{ version_stability.refseq_accessions_drift_pct | dp(1) }}% of RefSeq accessions ever
drift, and the worst-affected 5% hold {{ version_stability.refseq_top5_drift_share_pct | dp(0) }}%
of it.

This is predictable from a single version, with no
flanking pair or genomic context. A transcript version's **intrinsic CDS structure** (its
CDS length plus the lengths of its coding-exon segments in transcript coordinates)
is **build-independent** (identical across GRCh37/GRCh38/T2T for
{{ version_stability.refseq_struct_portable_pct | dp(1) }}% of RefSeq and
{{ version_stability.ensembl_struct_portable_pct | dp(1) }}% of Ensembl versions), and a
change in that structure is almost exactly equivalent to a genomic-coordinate drift:
every drifting RefSeq bump ({{ version_stability.refseq_drift_struct_flagged_pct | dp(0) }}%)
and {{ version_stability.ensembl_drift_struct_flagged_pct | dp(1) }}% of drifting Ensembl
bumps carry an intrinsic-structure change, while conversely
{{ version_stability.refseq_struct_unchanged_preserved_pct | dp(0) }}% (RefSeq) /
{{ version_stability.ensembl_struct_unchanged_preserved_pct | dp(1) }}% (Ensembl) of
structure-unchanged bumps are genomically preserved. This makes the safety decision
work even when the requested version was never aligned to the target build. To judge
substituting an available version for a requested one in, say, GRCh38, cdot reads the
requested version's intrinsic CDS structure from any build that carries it (the structure
is build-independent, so a GRCh37 or T2T record serves), reads the substitute version's
structure from GRCh38, and substitutes only when the two match; a match guarantees every
coding position maps to the same GRCh38 coordinate. Because the comparison is on the
structure itself rather than on flanking genomic coordinates, it also catches
transient-revert versions (a coordinate that goes A→B→A) that a genomic-bracket check
between neighbours cannot see.

cdot ships this as the safety gate behind the version fallback
(`intrinsic_cds_structure()` and the data-provider method
`is_version_substitution_safe()`): when the requested version is absent, the substitute
is applied only if its intrinsic CDS structure matches the requested version's, reported
as a coordinate-safe `HGVSFix`; a structural mismatch (or a version absent from every
build, where only a probabilistic genomic-bracket check remains) is refused by default
rather than silently substituted, preserving exact-variant semantics.

## R6: Throughput

**[Tier 1]** To compare transcript backends fairly we held the sequence layer constant
(every configuration was served by the *same local SeqRepo*), so the only thing that
varies across rows of Table 2 is the transcript-data layer itself.

**Table 2. End-to-end HGVS resolution throughput by transcript backend**, sequence layer
held constant (shared local SeqRepo), identical biocommons/hgvs engine. *(Tier 1.)*

| Configuration | Throughput (HGVS/s) |
|---|---|
| UTA: public remote database | ~{{ benchmark.uta_remote_tps | dp(1) }} |
| UTA: local PostgreSQL | ~{{ benchmark.uta_local_tps | int }} |
| cdot REST (one request per transcript) | ~{{ benchmark.cdot_rest_tps | int }} |
| cdot REST (after one batch `prefetch()`) | ~{{ benchmark.cdot_rest_prefetch_tps | int }} |
| cdot local JSON | {{ benchmark.cdot_local_min_tps | int }}–{{ benchmark.cdot_local_max_tps | int }} |

A GRCh38 RefSeq JSON file loads in ~{{ benchmark.grch38_load_time_s | dp(0) }} s and then
resolves at {{ benchmark.cdot_local_min_tps | int }}–{{ benchmark.cdot_local_max_tps |
int }} HGVS/s through the biocommons engine. The REST provider serves
~{{ benchmark.cdot_rest_tps | int }} HGVS/s when each transcript version is fetched in its
own request; batching those lookups into one `prefetch()` request amortises the
per-request network and processing overhead across the whole set and warms the transcript
cache, so later lookups are in-memory hits and REST throughput rises to match local. A locally
loaded UTA running the identical engine reached only
~{{ benchmark.uta_local_tps | int }} HGVS/s, and the public remote UTA database only
~{{ benchmark.uta_remote_tps | dp(1) }} HGVS/s. cdot's local data layer is thus roughly
30× faster than a local UTA running the identical engine; the full set of remote and
local timings is in Table 2. Once the data layer is local, the remaining bottleneck is
sequence fetching rather than transcript lookup.

**[Tier 1]** At scale, a single local-JSON process resolved the entire set of 3,660,452
unique ClinVar (g.HGVS, c.HGVS) pairs in ~92 minutes (665 HGVS/s; 99.3% produced a
genomic coordinate, 98.8% matched the ClinVar genomic HGVS exactly). The REST provider
matched this over the network: after one batch cache-warming pass (21,277 distinct
transcripts warmed in ~6 s) it resolved the same set in ~83 minutes (731 HGVS/s),
effectively matching local JSON: once the cache is warm both providers resolve from an
in-memory dict with no further network I/O, so throughput is bounded by the shared
sequence layer. The marginal difference is within run-to-run variance, plausibly because
the warmed REST cache holds only the few thousand transcripts the set actually touches
whereas local JSON holds the entire dataset in memory. The same exhaustive pass is impractical against the
public remote UTA database (extrapolated at hundreds of days from its ~0.1 HGVS/s).
(`paper/scripts/build_clinvar_pairs.py` builds the pair set by joining ClinVar's
variant_summary with the ClinVar VCF; the residual ~1% non-exact/error rate is dominated
by genomic-HGVS normalisation differences between the ClinVar source string and the
biocommons output, not resolution failures.)

## R7: T2T-CHM13v2.0 coverage

**[Tier 1]** cdot is the first transcript data source to bring the T2T-CHM13v2.0 assembly
to the Python HGVS libraries (biocommons/hgvs and PyHGVS), contributing
{{ coverage.t2t_unique_count | commas }} transcript-version alignments that are absent
from UTA and from every other backend these libraries can use. (Ensembl VEP can already
generate HGVS against T2T, but it is a standalone annotation tool, not a transcript
backend for these libraries.) cdot's JSON format also stores
per-exon alignment-gap information (indels of the transcript relative to the genome) so
that downstream libraries can apply it during coordinate conversion; this is a property
of the stored data rather than a separate cdot result.

---

*Notes: keep Tier-1/Tier-2 provenance explicit in every reported number. Audit
"incorrect" buckets
before reporting (alignment diffs #95). Normalise before string-compare.*
