# Supplementary Material

## Supplementary Methods

### Data sources and generation (detail)

cdot merges three transcript-annotation sources into a single dataset per build
(Figure 1A), downloading the complete run of historical releases from the RefSeq and
Ensembl FTP sites:

1. **RefSeq GFF3** [@OLeary2016]: {{ sources.refseq_grch37_releases | int }} GRCh37,
   {{ sources.refseq_grch38_releases | int }} GRCh38, and
   {{ sources.refseq_t2t_releases | int }} T2T-CHM13v2.0 NCBI annotation releases
   through RS_2025_08 (Supplementary Table S1), plus NCBI's historical
   transcript-alignment file (RS_2024_08), covering GRCh38 versions replaced or
   suppressed before any archived release.

2. **Ensembl GTF** [@Martin2023]: {{ sources.ensembl_grch37_releases | int }} GRCh37,
   {{ sources.ensembl_grch38_releases | int }} GRCh38, and
   {{ sources.ensembl_t2t_releases | int }} T2T-CHM13v2.0 releases, spanning Ensembl
   releases 81 to 116 (Supplementary Table S2); GTF is used because Ensembl's GFF3
   files omit transcript protein versions.

3. **UTA**: the Universal Transcript Archive computes its own alignments, so it holds
   some accessions never published in any annotation file; these fill gaps predating
   the oldest releases and are overridden by official annotation where available.

Releases are merged chronologically, newer entries overwriting older ones for the same
transcript version; this also backfills fields missing from early files, such as the
alignment-gap CIGAR strings absent from early RefSeq GFF3. A Snakemake pipeline parses
each source with HTSeq, normalises contig names via bioutils [@bioutils], extracts
coding-sequence (CDS) boundaries from start/stop codon features, attaches gene metadata
(symbols, HGNC IDs, biotypes) from HGNC and NCBI gene-info files, and serialises to
gzip-compressed JSON.

### JSON format (detail)

JSON parses quickly in every major language and serialises cleanly over HTTP, so the
same files drive both the in-memory provider and the REST API. Each transcript entry
stores shared metadata plus a `genome_builds` dict keyed by assembly name, so one entry
carries a transcript's GRCh37, GRCh38, and T2T-CHM13v2.0 alignments and a single-build
file is a drop-in GTF/GFF replacement. Per-build fields include a per-exon `gap` string
recording transcript-vs-genome indels, converted to HGVS CIGAR format at query time
(transcript bases that align to neither of two adjacent exons become a transcript-only
run at the end of the first), and
canonical tags (MANE Select and Plus Clinical, RefSeq Select, Ensembl canonical) where
present; a root-level `schema_version` lets clients reject incompatible files on load.
Supplementary Table S3 summarises the schema; the full JSON Schema ships in the
repository (`docs/cdot-json-schema.json`).

### String cleaning (`clean_hgvs()`) (detail)

Many of the descriptions that reach a curation platform do not parse: copy/paste
whitespace, lost casing, transposed punctuation, a swapped gene symbol and accession, a
trailing protein annotation. The intended variant is usually unmistakable, but a strict
grammar rejects the string; cdot repairs it instead.

`clean_hgvs()` is a pure string operation: it takes an HGVS string and returns the
repaired string with structured `HGVSFix` records describing every change. It needs no
genome build, sequence, parser, or data provider, so it can run as a pre-pass anywhere,
whichever downstream library consumes the result. Cleaning is an ordered pipeline of
single-purpose operations that strip extraneous characters, repair structural
punctuation, normalise casing and prefixes, and rebuild the canonical
`transcript(GENE):c.` shape; the per-operation catalogue, with examples, is
Supplementary Table S7. Each `HGVSFix` carries a severity, a stable code, a message,
and the before/after values, so a caller can audit exactly what changed. What cleaning
cannot fix (no colon at all, a bare variant body, an insertion given as a length instead
of a sequence) is flagged as an `ERROR`-level fix; cleaning never raises and returns its
best attempt.

One repair needs transcript data: a bare-number accession whose RefSeq prefix was
dropped entirely (`000059.4:c.68del`). `resolve_missing_accession_prefix()`, applied by
`fix_hgvs()` when a data provider is supplied, generates the candidate accessions the
kind letter allows and restores the prefix only when exactly one exists in the loaded
data, reported as an `HGVSFix` like every other repair.

### Version substitution (detail)

A separate, opt-in helper, `get_best_transcript_version()`, addresses transcript-version
drift. Substituting a different version is a heuristic that can be wrong, so it is
never applied automatically: when the requested version is absent from a caller's data,
the helper returns the best available adjacent version (up-then-down, closest, or
latest strategies) as a reported `HGVSFix`, and the caller decides whether to accept
it. Before a substitution is offered as safe, cdot applies a build-independent check: a
transcript version's intrinsic CDS structure (CDS length plus coding-exon segment
lengths in transcript coordinates) is the same in every build, and a substitution
passes only when the two versions' structures match, their CDS alignment gaps match,
and, for a UTR variant, the relevant UTR length is preserved
(`intrinsic_cds_structure()`, `is_version_substitution_safe()`; Supplementary Methods).
The check is validated against ClinVar in Results R5.

### Access and client libraries (detail)

cdot's client stack (Figure 1B) offers three data providers behind the same interface.

**Local JSON**: `JSONDataProvider` loads a JSON.gz file into memory with lazy interval
trees and lookup dictionaries; retrieval is O(1).

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves the same JSON
data at cdotlib.org. `RESTDataProvider` fetches one transcript version per request,
suited to occasional lookups without downloading the full file, and can warm its cache
with a single batched `prefetch()` request (Results R3).

**Ensembl TARK**: `EnsemblTarkDataProvider` exposes the Ensembl Transcript Archive
(TARK) REST service through the same interface.

**biocommons/hgvs integration**: cdot implements the full
`biocommons.hgvs.dataproviders.interface.Interface`, so it is a drop-in replacement for
UTA in any biocommons/hgvs pipeline:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
hdp = JSONDataProvider(["cdot-0.2.36.refseq.GRCh38.json.gz"])
```

Sequence data comes from SeqRepo [@Hart2020] or cdot's `FastaSeqFetcher`, which
reconstructs a transcript sequence by splicing its exon ranges out of a local genome
FASTA, enabling fully offline operation. That reproduces an Ensembl transcript exactly
but is not guaranteed to match a curated RefSeq transcript, whose sequence can differ
from the genome at a few positions (Supplementary Methods). `ChainedSeqFetcher` tries
several sequence sources in a caller-defined order.

**PyHGVS integration**: `JSONPyHGVSTranscriptFactory` exposes the same data to the
Counsyl PyHGVS library. PyHGVS is no longer maintained, so new development targets the
biocommons/hgvs path.

### Canonical transcript selection (detail)

`resolve_gene_hgvs()` maps a gene symbol to a MANE Select (or MANE Plus Clinical)
transcript accession [@Morales2022; @Wright2023] via the data provider's
`get_tx_ac_tags_for_gene()`, which ranks a gene's transcripts by tag priority, so a
pipeline given `BRCA1:c.68_69del` can resolve it through `NM_007294.4`. The selection
is a sensible default for search and gene-name lookup, not an authoritative clinical
choice.

### Benchmarking (detail)

Resolution accuracy and throughput are measured with scripts committed to the
repository (`paper/scripts/`); protocol detail is in Supplementary Methods.
`benchmark_resolution.py` resolves real (g.HGVS, c.HGVS) pairs through a pluggable
provider (local JSON, REST, or UTA) and reports resolution rate, recovery from cleaning
and version substitution, and speed; the ClinVar pair set is built by
`build_clinvar_pairs.py`. Cleaning is evaluated on a production query corpus and, as a
reproducible control, with `inject_and_clean.py`, which injects each fix category into
clean ClinVar strings. `lovd_head_to_head.py` runs the same injected cases through both
`clean_hgvs()` and the LOVD HGVS syntax checker [@LovdHgvsChecker]
({{ lovd_comparison.lovd_version }}, run locally as a PHP CLI), scoring a case as
recovered when the tool's output (for LOVD, its top-ranked correction) exactly matches
the canonical target; the same false-correction check runs over the uncorrupted
originals.

`vv_mutalyzer_head_to_head.py` extends the same cases and scoring rule to the two
sequence-aware validation services, VariantValidator [@Freeman2018] and Mutalyzer
[@Lefter2021], over their public REST APIs
({{ vv_mutalyzer_comparison.service_date }}; remote services cannot be version-pinned,
so the service date and each API's reported version are recorded instead). On the
uncorrupted originals, a valid input the service *alters* is a false correction while one
it *rejects* is counted separately as a validity or coverage position. Version-substitution
safety is measured by `compute_version_stability.py` on GRCh38, over a seeded
{{ version_stability.sample_n | commas }}-accession sample of accessions cdot holds at
two or more versions; the same run bins preserved coding bases by relative CDS position
(Supplementary Figure S1).

The submitted-string corpus (Results R2) is built by
`build_clinvar_submitted_pairs.py` from the per-submission HGVS attributes of a ClinVar
VCV XML release (ClinVarVCVRelease_2026-10): each SCV's first transcript c./n.
expression is joined verbatim to the variant's VCF coordinate via its AlleleID and
collapsed to unique (AlleleID, string) pairs, each tagged with the earliest SCV
submission date (Supplementary Methods). Scoring uses the VCF coordinate rather than
the g.HGVS string, since a submitted string may legitimately spell an indel differently
from ClinVar's normalised form. Transcript version age is computed by
`compute_submitted_version_age.py` against the released cdot RefSeq JSON, whose
per-transcript source URL identifies whether an entry survives only through cdot's
historical releases.

Throughput (Table 1) is measured by `compute_benchmark.py`: every configuration
resolves the identical committed set of {{ benchmark.n_pairs | commas }} ClinVar pairs
through the same biocommons/hgvs engine with the sequence layer held constant (one
shared local SeqRepo, file-descriptor cache enabled), so only the transcript-data layer
varies between rows. Each configuration is timed over
{{ benchmark.n_repeats | int }} passes and reported as median (IQR), steady-state
resolution only (exclusions in Supplementary Methods). The REST rows run against the
production server at cdotlib.org over the public internet; the public remote UTA row
uses the first {{ benchmark.uta_remote_n | int }} pairs, a full-size pass being
impractical at ~{{ benchmark.uta_remote_tps | dp(1) }} HGVS/s.

Results marked **[private data]** are aggregate statistics from private production
corpora, reported as frozen constants that a referee cannot regenerate; every other
result is reproducible from public data with the scripts in the repository.

### Submitted-string corpus construction

`build_clinvar_submitted_pairs.py` streams the ClinVar VCV XML (tens of gigabytes
uncompressed; a pre-extracted per-SCV table is also accepted). For each
ClinicalAssertion (SCV), the first RefSeq/Ensembl transcript c./n. expression among its
`Attribute[@Type="HGVS"]` values is kept verbatim; this first-attribute selection is
what defines "the" submitted string when an SCV carries several HGVS expressions.
Submitted HGVS is not deduplicated across SCVs in ClinVar, so the builder collapses SCV
rows to unique (AlleleID, submitted string) pairs, retaining the collapsed submission
count; distinct strings for the same variant (for example, two laboratories citing
different transcript versions) are all kept. The cdot-versus-UTA comparison (Results R2)
is run over the whole corpus; a committed 500-pair random sample (seed 42) is for quick
reproduction.

### LOVD head-to-head scoring detail

The two ecosystems canonicalise a parenthesised gene symbol in opposite directions
(biocommons keeps `NM_x.y(GENE):c.`, LOVD removes the symbol), so the comparison
ignores that annotation on both sides. A secondary metric accepting the target anywhere
in LOVD's ranked correction list gave identical results, so top-1 ranking cost LOVD
nothing. Every injected category is string-repairable by construction
(`inject_and_clean.py` does not inject errors whose repair needs transcript data, such
as a missing accession prefix), so both tools are eligible on every case.

### Shariant corpus resolution protocol

The Shariant corpus (Results R2, private data) is resolved as pure coordinate projection with
`replace_reference=False`, so the sequence layer never differs between the two
backends.

### RefSeq transcript sequence and `FastaSeqFetcher`

`FastaSeqFetcher` reconstructs transcript sequence from the genome, which is not
guaranteed to match a curated RefSeq transcript at every position. In practice variants
are usually called against a genome-mapped read alignment, so any such discrepancy lies
in the reference the variant was already described against. Where an exact RefSeq
transcript sequence is needed, SeqRepo supplies it when it holds the accession (it does
not hold every version cdot covers), and `ChainedSeqFetcher` can prefer SeqRepo with a
FASTA fallback to maximise the fraction of transcripts served.

### Throughput measurement detail

The timed portion is steady-state resolution throughput (HGVS parse plus coordinate
projection per string); provider setup (JSON load, database connection), the REST
cache-warming `prefetch()`, and one untimed warm-up pass per configuration are
excluded. The absolute numbers are sensitive to the shared sequence layer: with
SeqRepo's file-descriptor cache disabled (its default), SeqRepo re-opens a compressed
FASTA file on every fetch, every configuration is bounded by sequence fetching at well
under 200 HGVS/s, and the differences between transcript backends are masked. All
configurations therefore run with the cache enabled.

### Version-substitution safety gate detail

To judge substituting an available version for a requested one in, say, GRCh38, cdot
reads the requested version's intrinsic CDS structure from any build that carries it
(the structure is build-independent, so a GRCh37 or T2T record serves), reads the
substitute version's structure from GRCh38, and substitutes only when the two match.
Because the comparison is on the structure rather than on flanking genomic coordinates,
it also catches transient-revert versions (a coordinate that goes A→B→A) that a
genomic-bracket check between neighbours cannot see.

Identical transcript structure does not on its own guarantee safety, so the gate adds
two refinements. An alignment gap (a transcript-vs-genome indel) present in one version
and not the other shifts every coding base downstream of it, so the gate also requires
the two versions' CDS alignment gaps to match. And because the structure is
build-independent it cannot see the variant's position relative to the UTRs: a 5' or 3'
UTR variant can move when only the UTR length changes (UTR annotations change between
versions far more often than the CDS), so for a cited UTR position the matching UTR
length must be preserved too, while a coding variant on a UTR-only change stays safe.

As a predictor of genomic drift, an intrinsic-structure change flags
{{ version_stability.refseq_drift_struct_flagged_pct | pct(1) }} of drifting RefSeq
bumps and {{ version_stability.ensembl_drift_struct_flagged_pct | pct(1) }} of drifting
Ensembl bumps, while conversely
{{ version_stability.refseq_struct_unchanged_preserved_pct | pct(0) }} (RefSeq) /
{{ version_stability.ensembl_struct_unchanged_preserved_pct | pct(0) }} (Ensembl) of
structure-unchanged bumps are genomically preserved. The build-independent structure
cannot see a re-alignment of the same CDS structure to a different genomic locus, so
cdot additionally compares the two versions' genomic CDS maps when both are loaded, and
ships a small precomputed blocklist of these re-placements for the case where the
requested version is absent from the loaded data. A structural mismatch, or a version
absent from every build (where only a probabilistic genomic-bracket check remains), is
refused by default.

## Supplementary Results

These results support the main-text Results but are not about cdot itself: what the
ClinVar and Shariant corpora show about how laboratories cite transcripts, what string
repair cannot fix, and how transcript versions drift.

### Current-version ClinVar check

Resolution rate is compared against a locally loaded UTA (release `{{ clinvar_submitted._prov.tools.uta }}`),
first on current-version descriptions and then on the historical submitted-string
corpus. For the current-version check, the committed ClinVar [@Landrum2025] pair sets
({{ clinvar.n_refseq | commas }} RefSeq NM_ and {{ clinvar.n_ensembl | commas }} Ensembl
ENST preferred-transcript descriptions) were resolved through cdot, with transcript
sequence spliced from a local genome FASTA so that no Ensembl pair is dropped for a
sequence SeqRepo lacks.
UTA's rate is membership of the exact cited accession.version
in its GRCh38 alignment set, which is what resolution through UTA requires; the live
database gives the same rates. On RefSeq the two backends are at parity (cdot
{{ clinvar.cdot_refseq_pct | pct(1) }}, UTA {{ clinvar.uta_refseq_pct | pct(1) }}).
Because the pair sets are Ensembl-enriched, their aggregate rate is not a like-for-like
accuracy score and we do not report one.

On the full variant_summary corpus, the residual also includes
{{ clinvar_vcf_residual.unaligned_between_exons_transcripts | number_word }} transcripts
with bases that align to neither of two adjacent exons, where cdot mostly gives the same
coordinate as UTA and ClinVar places the unaligned bases differently (Supplementary
Table S4).

### Transcript versions cited in ClinVar submissions

The submitted-string corpus reflects the natural submission mix, which is essentially
all RefSeq (not one submitted string cites an Ensembl transcript; Ensembl accessions
enter ClinVar through its recomputed names, not submitter descriptions). Weighted by
submission count, {{ clinvar_submitted.scv_weighted_not_current_pct | pct(1) }} of
submitted strings cite a superseded version, and
{{ clinvar_submitted.base_retired_pct | pct(1) }} cite a transcript no longer annotated
at any version. Most submitted traffic is historical, and
only {{ clinvar_submitted.absent_cdot_pct | pct(1) }} of cited versions are absent from
cdot's GRCh38 data.

Because cdot ingests every RefSeq and Ensembl annotation release, a version absent from
its data should be one that never appeared in an official release, and the absent
citations bear that out: of {{ submitter_attribution.n_submitters_total | commas }}
submitters only {{ submitter_attribution.n_submitters_with_absent | int }} ever cite a
version cdot lacks, a single large laboratory accounts for
{{ submitter_attribution.top1_absent_share_pct | pct(1) }} of all absent-version
citations, and {{ submitter_attribution.single_submitter_version_pct | pct(1) }} of the
{{ submitter_attribution.n_absent_versions | int }} distinct absent versions are cited by
one submitter alone.
These are the signature of transcripts a laboratory aligned to the
genome itself rather than drawing from a published release, so they are addressed by
version substitution, not by deeper ingest.

The RefSeq gap invisible at the current-version ceiling (both backends
{{ clinvar.cdot_refseq_pct | pct(1) }} above) opens to
{{ clinvar_submitted.cdot_only_pct | dp(1) }} points on the submitted strings:
{{ clinvar_submitted.cdot_only | commas }} pairs resolve through cdot alone and
{{ clinvar_submitted.uta_only | commas }} through UTA alone, because UTA holds no GRCh38
alignment for {{ clinvar_submitted.uta_no_data_pct | pct(1) }} of the cited versions.
Splitting by submission era, cdot reproduces the coordinate for
{{ clinvar_submitted.era_old_cdot_pct | pct(1) }} of the oldest (2008-2015)
submissions, {{ clinvar_submitted.era_mid_cdot_pct | pct(1) }} of 2016-2020 and
{{ clinvar_submitted.era_recent_cdot_pct | pct(1) }} of 2021-2026, while UTA reaches
{{ clinvar_submitted.era_old_uta_pct | pct(1) }},
{{ clinvar_submitted.era_mid_uta_pct | pct(1) }} and
{{ clinvar_submitted.era_recent_uta_pct | pct(1) }}: lowest on the oldest submissions,
which cite the most superseded versions, and lower on the newest than on 2016-2020.

The residual {{ clinvar_submitted.residual_pct | pct(1) }} after cleaning and version
substitution ({{ clinvar_submitted.residual_n | commas }} of
{{ clinvar_submitted.full_n | commas }}) is dominated by historical-version effects:
{{ clinvar_submitted_residual.version_refused | commas }} cite a version absent from the
data where substitution declines to act because coordinate safety cannot be verified,
and {{ clinvar_submitted_residual.coordinate_drift | commas }} resolve through the cited
version to a coordinate that differs from ClinVar's current interpretation; the remainder
cite a position or base that does not exist on the cited version, notation the
biocommons grammar rejects, or an accession absent entirely (Supplementary Table S8).

### Shariant: what cdot resolves that UTA does not

**[private data].** Of the Shariant strings cdot resolved but UTA could not
({{ historical.cdot_only_pct | pct(1) }} of the corpus),
{{ historical.cdot_only_historical_pct | pct(0) }} were RefSeq transcript versions for
which UTA holds no GRCh38 alignment and
{{ historical.cdot_only_ensembl_pct | pct(0) }} were Ensembl transcripts.

### Residual errors: the ceiling of cleaning

**[private data].** A small residue of non-HGVS input ({{ cleaning_corpus.nonhgvs_n | int }}
queries: pasted URLs, report templates, prose) that slipped the corpus regex is a
data-collection artifact and is excluded from the corpus above. The remaining
{{ cleaning_corpus.residual_pct | pct(1) }} ({{ cleaning_corpus.residual_n | commas }}
queries) of genuine HGVS-shaped input that still fails to parse after cleaning defines
the ceiling of pure string repair. Classified under a fixed decision-tree taxonomy
(Supplementary Table S6), well over half is incomplete or reference-less input that no
string-level repair can invent, about a quarter is in principle fixable and marks the
frontier for future cleaning rules, and the remaining
~{{ cleaning_corpus.residual_grammar_gap_pct | pct(0) }} is valid HGVS the biocommons
grammar rejects.

### Injection control

As a control, `paper/scripts/inject_and_clean.py` injects each `clean_hgvs()` fix
category into a seeded sample of clean, parseable ClinVar c.HGVS strings and confirms
the cleaner recovers the canonical target with
{{ cleaning.inject_regressions | int }} regressions (Supplementary Table S5).
Because this benchmark injects the very errors it then repairs, its recovery rate is
reported in the supplement only; its purpose is the no-regression guarantee on which
the production result depends.

### How transcript versions drift

Weighted by coding base (the chance a *random* variant is unaffected by substituting
the next version), safety is
{{ version_stability.refseq_pervariant_safety | pct(1) }} (RefSeq) and
{{ version_stability.ensembl_pervariant_safety | pct(1) }} (Ensembl).
When a coordinate does move it is usually the whole CDS
({{ version_stability.refseq_full_drift_pct | pct(1) }} of RefSeq pairs,
{{ version_stability.ensembl_full_drift_pct | pct(1) }} Ensembl), driven by a
re-annotation of the coding region; the most dangerous case, a *partial* substitution
that mis-places some variants but not others, is rare
({{ version_stability.refseq_partial_drift_pct | pct(1) }} of RefSeq pairs,
{{ version_stability.ensembl_partial_drift_pct | pct(1) }} Ensembl). Within that tail
the risk is positional: a partial drift keeps a 5' prefix intact up to its first
alignment change and moves what lies downstream, so preservation falls from
{{ positional_drift.refseq_partial_decile1_pct | pct(0) }} of coding bases in the
5'-most decile to {{ positional_drift.refseq_partial_decile10_pct | pct(0) }} in the
3'-most for RefSeq ({{ positional_drift.ensembl_partial_decile1_pct | pct(0) }} to
{{ positional_drift.ensembl_partial_decile10_pct | pct(0) }} for Ensembl; Supplementary
Figure S1). Drift is also concentrated: only
{{ version_stability.refseq_accessions_drift_pct | pct(1) }} of RefSeq accessions
({{ version_stability.ensembl_accessions_drift_pct | pct(1) }} of Ensembl) ever drift,
and the 5% of accessions with the most drifting pairs hold
{{ version_stability.refseq_top5_drift_share_pct | pct(1) }} of them
({{ version_stability.ensembl_top5_drift_share_pct | pct(1) }} for Ensembl), so the risky
transcripts form a short, identifiable list that can be flagged or blocklisted.

## Supplementary Tables

### Table S1: RefSeq GFF3 annotation releases

| Genome build | Annotation release | Transcripts added |
|-------------|-------------------|------------------|
| GRCh37 | ... | ... |
| GRCh38 | ... | ... |
| T2T-CHM13v2.0 | ... | ... |

*Generated from Snakemake pipeline summary stats.*

### Table S2: Ensembl GTF releases

| Genome build | Ensembl release | Transcript count |
|-------------|-----------------|-----------------|
| GRCh37 | 82, 85, 87 | ... |
| GRCh38 | 81–116 | ... |
| T2T-CHM13v2.0 | rapid release 2022_06, 2022_07 | ... |

### Table S3: JSON schema fields (data release 0.2.36)

A summary of the fields biocommons/hgvs resolution consumes; the full JSON Schema is
`docs/cdot-json-schema.json` in the repository.

| Field | Level | Description |
|-------|-------|-------------|
| `schema_version` | root | Schema compatibility version (major.minor checked on load) |
| `transcripts` | root | Dict keyed by transcript accession.version |
| `gene_name` | transcript | Gene symbol |
| `hgnc` | transcript | HGNC ID |
| `biotype` | transcript | Sorted list of transcript types, e.g. `["mRNA", "protein_coding"]` |
| `protein` | transcript | Protein accession.version (coding only) |
| `start_codon` / `stop_codon` | transcript | 0-based transcript (cDNA) CDS start and exclusive end (coding only) |
| `genome_builds` | transcript | Dict keyed by build name (`GRCh37`, `GRCh38`, `T2T-CHM13v2.0`) |
| `contig` | build | Chromosome/contig accession, e.g. `NC_000007.14` |
| `strand` | build | `"+"` or `"-"` |
| `exons` | build | List of `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]`; `gap` records transcript-vs-genome indels, `null` when the exon aligns without one |
| `cds_start` / `cds_end` | build | 0-based genomic CDS coordinates (coding only) |
| `tag` | build | Comma-separated canonical tags from the source, e.g. `MANE_Select,Ensembl_canonical` (MANE Select, MANE Plus Clinical, RefSeq Select, Ensembl canonical) |
| `url` | build | Source annotation file the alignment was taken from |

### Table S4: ClinVar benchmark details

Full-scale resolution of every c.HGVS in ClinVar's variant_summary through cdot alone
(GRCh38, cdot {{ clinvar_vcf._prov.tools.cdot }} with data release {{ clinvar_vcf._prov.tools.cdot_data }}; ClinVar {{ clinvar_vcf._prov.tools.clinvar }}), scored as a VCF
coordinate (CHROM/POS/REF/ALT) against ClinVar's own VCF, not as a g.HGVS string, so
equivalent representations and ClinVar's tandem-repeat / identity notations are not
miscounted.

**Caveat.** ClinVar submissions are dominated by a handful of large (largely US) clinical
laboratories citing mostly current RefSeq versions, so this is a reproducible scale
check, not an unbiased sample of the transcripts clinical labs use. ClinVar's
`variant_summary` Name column is RefSeq-centric, so the Ensembl share at this scale is
near zero; ClinVar's comprehensive Ensembl HGVS lives in `hgvs4variation.txt.gz`, which
this pass does not ingest. The Ensembl resolution evidence is therefore the R2
per-source split, not this table.

Total pairs: {{ clinvar_vcf.n_pairs | commas }}.

| Outcome | Count | % |
|---|---|---|
| Resolved | {{ clinvar_vcf.n_resolved | commas }} | {{ clinvar_vcf.resolved_pct | dp(2) }} |
| Matched (cdot coordinate = ClinVar VCF) | {{ clinvar_vcf.n_correct | commas }} | {{ clinvar_vcf.matched_pct | dp(2) }} |
| Coordinate differs | {{ clinvar_vcf.incorrect | commas }} | {{ clinvar_vcf.incorrect_pct | dp(2) }} |
| Unresolved (no data) | {{ clinvar_vcf.no_data | commas }} | |
| Unresolved (input not parseable / convertible) | {{ clinvar_vcf.error | commas }} | |

Of the {{ clinvar_vcf.matched_of_resolved_pct | pct(2) }} match rate among resolved
variants, the {{ clinvar_vcf_residual.n_incorrect | commas }} coordinate differences are
overwhelmingly representation, multi-mapping or alignment choice, not cdot errors:
{{ clinvar_vcf_residual.unaligned_between_exons | commas }} on
{{ clinvar_vcf_residual.unaligned_between_exons_transcripts | int }} transcripts with
transcript bases that align to neither of two adjacent exons (SHANK3 `NM_001372044.2` and
others), most of them the same coordinate UTA gives, where ClinVar's numbering
distributes the unaligned bases differently;
{{ clinvar_vcf_residual.indel_mismatch | commas }} indel / left-alignment representation
differences; {{ clinvar_vcf_residual.snv_diff_pos | commas }} SNVs concentrated in
{{ clinvar_vcf_residual.snv_diff_pos_transcripts }} transcripts (paralog and copy-number
genes that map to more than one genomic locus, each with a constant per-transcript offset);
{{ clinvar_vcf_residual.long_cdot_allele | commas }} alleles much longer than ClinVar's, almost
all c. ranges whose two ends lie in different exons, which biocommons/hgvs maps to a
genomic span that includes the intron; {{ clinvar_vcf_residual.ambiguity_code_allele }} IUPAC ambiguity-code
alleles (position matches, only the degenerate base differs); and
{{ clinvar_vcf_residual.identity_or_symbolic }} identity/symbolic alleles. On the
transcripts with unaligned bases, an intronic offset anchored on one of those bases stays
unresolved, because biocommons/hgvs drops the offset when its anchor base has no
genomic position.

Binning these coordinate differences by relative CDS position (the decile scheme of
Figure S1, cited position over CDS length) shows no 3'-end concentration: correct projections
are nearly uniform across deciles
({{ clinvar_residual_positions.correct_decile10_pct | pct(1) }} in the 3'-most decile)
and the {{ clinvar_residual_positions.incorrect_binned | commas }} binned differences
track them with a mild mid-CDS excess and a depleted 3'-most decile
({{ clinvar_residual_positions.incorrect_decile10_pct | pct(1) }}), consistent with
their representation, multi-mapping and alignment-choice origins rather than positional alignment drift
(`bin_residual_positions.py`).

### Table S5: Injection cleaning benchmark, per category

Reproducible injection benchmark (`paper/scripts/inject_and_clean.py`): each
`clean_hgvs()` fix category is injected into a seeded sample of {{ cleaning.inject_sample_size | commas }}
clean, parseable public ClinVar c.HGVS strings (committed to the repository, seed 112,
at most 200 cases per category), and recovery is an exact string match to the known
canonical target. The LOVD columns come from `paper/scripts/lovd_head_to_head.py`
(LOVD HGVS syntax checker {{ lovd_comparison.lovd_version }}, local PHP CLI; "top-1" scores the checker's
highest-confidence suggested correction) and the VariantValidator (VV) and Mutalyzer
columns from `paper/scripts/vv_mutalyzer_head_to_head.py` (public REST APIs,
{{ vv_mutalyzer_comparison.service_date }}), all under the scoring rule in Methods.
Production ops with no string-level injector (structure reconstruction, empty-version
dropping, provider-verified accession-prefix restoration) are absent by design.
Examples are synthesised from public `NM_000059.4` (BRCA2). Totals also appear in
`paper/empirical_results/cleaning.csv`, `lovd_comparison.csv` and
`vv_mutalyzer_comparison.csv`.

<!-- include-csv: empirical_results/injection_benchmark.csv
  align: left
-->

Weighted by the production rescue-op distribution (Results Table 2) the totals are
{{ lovd_comparison.cdot_weighted_pct | pct(1) }} for `clean_hgvs()`,
{{ lovd_comparison.lovd_top1_weighted_pct | pct(1) }} for LOVD top-1,
{{ vv_mutalyzer_comparison.vv_weighted_pct | pct(1) }} for VariantValidator and
{{ vv_mutalyzer_comparison.mut_weighted_pct | pct(1) }} for Mutalyzer; accepting the
target anywhere in LOVD's ranked correction list changes no case. On the
{{ lovd_comparison.originals_n | commas }} uncorrupted originals no tool falsely
corrects any input ({{ lovd_comparison.cdot_false_corrections | int }} for `clean_hgvs()`,
{{ lovd_comparison.lovd_false_corrections | int }} for LOVD,
{{ vv_mutalyzer_comparison.vv_false_corrections | int }} for VariantValidator,
{{ vv_mutalyzer_comparison.mut_false_corrections | int }} for Mutalyzer); LOVD flags
{{ lovd_comparison.lovd_flagged_invalid_pct | pct(1) }} of them (intronic positions on
a transcript reference) as requiring a genomic reference while leaving the string
unchanged. The LOVD partial rates are one-sided accession-family support: the
gene/transcript swap is repaired for RefSeq but not Ensembl accessions, accession
re-casing for Ensembl but not RefSeq, and the two rates track the corpus's roughly
even RefSeq/Ensembl split. The protein-suffix misses absorb the stray `p` into the
edit (`c.68del p.` becomes `c.68delP`).

The services' sub-100% ceilings in the categories they otherwise repair are
service-level rejections independent of the injected error, measured on the same
uncorrupted originals: VariantValidator rejects
{{ vv_mutalyzer_comparison.vv_rejected_valid | int }} of the
{{ vv_mutalyzer_comparison.originals_n | commas }} valid originals (all recent Ensembl
transcripts absent from its vvta_2025_02 database), and Mutalyzer rejects
{{ vv_mutalyzer_comparison.mut_rejected_valid | int }} (`EINTRONIC`: intronic positions
on a RefSeq transcript reference; it resolves intronic Ensembl descriptions).
Mutalyzer's lowercased-accession rate is one-sided in the opposite direction to LOVD's
(it re-cases `enst...` but not `nm_...`). VariantValidator's responses for input it
could not parse embedded the LOVD syntax checker's ranked suggestions in
{{ vv_mutalyzer_comparison.vv_lovd_suggestion_cases | commas }} of the
{{ vv_mutalyzer_comparison.n_cases | commas }} cases without applying them (Discussion).

### Table S6: Residual error classes after cleaning

**[private data].** Single-label classification of the {{ cleaning_corpus.residual_n | commas }}
genuine-HGVS production queries that still fail to parse after cleaning (Supplementary
Results, "Residual errors"), under a fixed decision-tree taxonomy. The
{{ cleaning_corpus.nonhgvs_n | int }} non-HGVS queries that slipped the corpus regex
are excluded (Results). Counts and % below are of the
{{ cleaning_corpus.residual_n | commas }} residual queries; examples are synthesised
from public BRCA2 `NM_000059.4`.

<!-- include-csv: empirical_results/residual_taxonomy.csv
  align: left
-->

*Method and limitation:* classification was performed by a large language model (Claude
Opus 4, Anthropic; 2026-06-17) applying the shared decision tree to each unique string,
single-label and single-rater; no second-rater adjudication was done, so no inter-rater
agreement (κ) is reported. The taxonomy is version `v1`. When a later cdot release
rescued further residual queries, only those queries were removed from the class the
earlier rating had given them (delta re-rating by Claude Fable 5 and Claude Opus 5.5,
Anthropic, on 2026-08-17 and 2026-11-02); no other class changed. Counts are for cdot
0.2.37. Synthesised examples (from public `NM_000059.4` / `NM_001754.5`) illustrate
each class; no corpus string is reproduced.

### Figure S1: Positional drift along the CDS across version bumps

![](paper/figures/figure_s1_positional_drift.svg)

**Figure S1. Coordinate preservation by relative CDS position across consecutive
transcript version bumps** (GRCh38, seeded
{{ version_stability.sample_n | commas }}-accession sample per consortium; deciles of the
shared CDS, 1 = 5'-most tenth). **(A)** Conditioned on the partial-drift bumps
({{ positional_drift.refseq_partial_pairs }} RefSeq and
{{ positional_drift.ensembl_partial_pairs }} Ensembl pairs, the
{{ version_stability.refseq_partial_drift_pct | pct(1) }} and
{{ version_stability.ensembl_partial_drift_pct | pct(1) }} of bumps where some coding
bases move and others do not): preservation declines monotonically toward the 3' end,
from {{ positional_drift.refseq_partial_decile1_pct | pct(1) }} to
{{ positional_drift.refseq_partial_decile10_pct | pct(1) }} of coding bases for RefSeq
and {{ positional_drift.ensembl_partial_decile1_pct | pct(1) }} to
{{ positional_drift.ensembl_partial_decile10_pct | pct(1) }} for Ensembl, because a
partial drift keeps a 5' prefix intact up to its first alignment change and shifts the
bases downstream of it. **(B)** Unconditioned over every compared bump the curve is
essentially flat (RefSeq {{ positional_drift.refseq_all_decile1_pct | pct(1) }} to
{{ positional_drift.refseq_all_decile10_pct | pct(1) }}; Ensembl
{{ positional_drift.ensembl_all_decile1_pct | pct(1) }} to
{{ positional_drift.ensembl_all_decile10_pct | pct(1) }}): most drift is a whole-CDS
relocation, position independent and reliably flagged by the intrinsic-structure check
(Results R5). Produced by `compute_version_stability.py` (facts) and
`make_positional_figure.py` (rendering).

### Table S7: `clean_hgvs()` operation catalogue

The cleaning pipeline (Methods) applies these operation groups in a fixed canonical order.
Each operation inspects the string, makes at most one class of change, and records an
`HGVSFix` if it fired.

- **Stripping**: stray leading characters before the accession (e.g. a `GRCh38.p2` build tag or a
  stray `#`/`:`), all internal and surrounding whitespace and non-printable characters,
  wrapping quotes/backticks, trailing separators, and brackets only when unbalanced
  (balanced parentheses are preserved, since they are valid in uncertain-range and
  gene-symbol notation).
- **Structural punctuation**: collapsing a doubled colon, dot, or underscore
  (`NM_000059..4` → `NM_000059.4`); repairing a misplaced colon in the accession prefix;
  normalising a gene symbol wedged between extra colons or stray parentheses so that
  `NM_000059.4:(BRCA2):c.…` and `BRCA1(NM_000059.4)c.…` become the canonical
  `transcript(GENE):c.…`; collapsing a doubled kind token (`c.c.` → `c.`); fixing a
  comma or colon used in place of the kind dot, a period used in place of a
  substitution `>`, a semicolon, underscore, or space used in place of the
  reference:allele colon (`NM_000059.4;c.68del` and `BRCA2 c.68del` →
  `NM_000059.4:c.68del` / `BRCA2:c.68del`); and dropping the dangling dot of an
  empty transcript version (`NM_000059.:c.68del` → `NM_000059:c.68del`).
- **Casing and prefixes**: uppercasing nucleotides in substitutions and del/ins/dup
  edits (`c.123delg` → `c.123delG`), lowercasing an uppercased mutation type while
  protecting gene symbols that contain those letters (so `NM_000059.4(INSR):c.…` is
  never corrupted to `insR`), restoring a missing `N` prefix
  (`M_000059.4` → `NM_000059.4`) or transcript underscore,
  and adding a missing `c.`/`g.` kind where the accession type makes it unambiguous.
- **Reconstruction and gene/transcript repair**: a lenient pattern that rebuilds the
  canonical `transcript(gene):kind.variant` shape from a mangled one (inserting a
  missing `:` or `.`, uppercasing the accession prefix and lowercasing the kind letter),
  and a final step that detects and repairs the common clinical mistake of swapping the
  gene symbol and transcript accession (`BRCA2(NM_000059.4):c.…` →
  `NM_000059.4(BRCA2):c.…`).

One further repair sits outside `clean_hgvs()` because it needs transcript data:
`resolve_missing_accession_prefix()` (applied by `fix_hgvs()` when a data provider is
supplied) restores a fully dropped RefSeq prefix (`000059.4:c.68del` →
`NM_000059.4:c.68del`) by generating the candidates the kind letter allows (`c.` →
`NM_`/`XM_`, `n.` → `NR_`/`XR_`; 7-8 digit fields zero-padded to the 9-digit form) and
applying the fix only when exactly one candidate accession exists in the loaded data
(Methods).

### Table S8: Residual outcomes on the submitted-string corpus

The {{ clinvar_submitted.residual_n | commas }} of
{{ clinvar_submitted.full_n | commas }} submitted-string pairs (Results R2)
that remain unresolved after cleaning and version substitution, by cause.

| Cause | n |
|---|---|
| Cited version absent from the data; substitution declines to act because coordinate safety cannot be verified (no false rescues, by design) | {{ clinvar_submitted_residual.version_refused | commas }} |
| Resolves through the cited historical version to a coordinate that differs from ClinVar's current interpretation | {{ clinvar_submitted_residual.coordinate_drift | commas }} |
| Cited position does not exist on the cited version | {{ clinvar_submitted_residual.position_out_of_bounds | commas }} |
| Cited reference base does not exist on the cited version | {{ clinvar_submitted_residual.reference_mismatch | commas }} |
| Repeat or allele notation the biocommons grammar rejects | {{ clinvar_submitted_residual.grammar_unsupported | commas }} |
| Accession absent entirely | {{ clinvar_submitted_residual.unknown_accession | commas }} |
