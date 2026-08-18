# Supplementary Material

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
| GRCh37 | 82–87 | ... |
| GRCh38 | 76–116 | ... |
| T2T-CHM13v2.0 | ... | ... |

### Table S3: JSON schema fields (v0.2.34)

| Field | Level | Description |
|-------|-------|-------------|
| `schema_version` | root | Schema compatibility version |
| `transcripts` | root | Dict keyed by transcript accession |
| `gene_name` | transcript | HGNC gene symbol |
| `biotype` | transcript | e.g. `"protein_coding"` |
| `genome_builds` | transcript | Dict keyed by build name |
| `contig` | build | Chromosome/contig accession |
| `strand` | build | `+1` or `-1` |
| `exons` | build | List of `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]` |
| `cds_start` / `cds_end` | build | CDS coordinates |
| `mane_select` | build | MANE Select accession if applicable |
| `mane_plus_clinical` | build | MANE Plus Clinical accession if applicable |
| `refseq_select` | build | RefSeq Select flag |
| `ensembl_canonical` | build | Ensembl canonical flag |

### Table S4: ClinVar benchmark details

Full-scale resolution of every RefSeq and Ensembl c.HGVS in ClinVar through cdot alone
(GRCh38, cdot 0.2.33). Unlike the cdot-vs-UTA comparison in Results R2 (gated to a sample
by UTA's throughput), every variant is summarised here because only cdot is in the loop. The
projection is scored as a VCF coordinate (CHROM/POS/REF/ALT) against ClinVar's own VCF,
not as a g.HGVS string, so equivalent representations and ClinVar's tandem-repeat / identity
notations are not miscounted.

**Caveat.** ClinVar submissions are dominated by a handful of large (largely US) clinical
laboratories citing mostly current RefSeq versions, so this is a clean, public,
reproducible scale check that cdot resolves real variants at scale, not an unbiased sample
of the transcripts clinical labs use. The unbiased real-world complement is the Shariant
historical corpus (Results R2, Tier 2). The pair builder extracts both RefSeq and Ensembl
c.HGVS and reports the measured source mix, but ClinVar's `variant_summary` Name column is
RefSeq-centric, so the Ensembl share at this scale is near zero. ClinVar's comprehensive
Ensembl HGVS lives in `hgvs4variation.txt.gz`, which this pass does not ingest; the
Ensembl resolution evidence is therefore the R2 sample (per-source split) and the Shariant
corpus, not this table.

Total pairs: {{ clinvar_vcf.n_pairs | commas }}.

| Outcome | Count | % |
|---|---|---|
| Resolved | {{ clinvar_vcf.n_resolved | commas }} | {{ clinvar_vcf.resolved_pct | dp(1) }} |
| Matched (cdot coordinate = ClinVar VCF) | {{ clinvar_vcf.n_correct | commas }} | {{ clinvar_vcf.matched_pct | dp(2) }} |
| Coordinate differs | {{ clinvar_vcf.incorrect | commas }} | {{ clinvar_vcf.incorrect_pct | dp(2) }} |
| Unresolved (no data) | {{ clinvar_vcf.no_data | commas }} | |
| Unresolved (input not parseable / convertible) | {{ clinvar_vcf.error | commas }} | |

Of the {{ clinvar_vcf.matched_of_resolved_pct | dp(2) }}% match rate among resolved
variants, the {{ clinvar_vcf_residual.n_incorrect | commas }} coordinate differences are
overwhelmingly representation or multi-mapping, not cdot errors:
{{ clinvar_vcf_residual.indel_mismatch | commas }} indel / left-alignment representation
differences; {{ clinvar_vcf_residual.snv_diff_pos | commas }} SNVs concentrated in just
{{ clinvar_vcf_residual.snv_diff_pos_transcripts }} transcripts (paralog and copy-number
genes that map to more than one genomic locus, each with a constant per-transcript offset);
{{ clinvar_vcf_residual.babelfish_megaallele | commas }} insertions over-shuffled by the
VCF normaliser; {{ clinvar_vcf_residual.ambiguity_code_allele }} IUPAC ambiguity-code
alleles (position matches, only the degenerate base differs); and
{{ clinvar_vcf_residual.identity_or_symbolic }} identity/symbolic alleles.

Binning these coordinate differences by relative CDS position (the decile scheme of
Figure S1, cited position over CDS length; a re-run of the pass, which reproduced the
committed totals within six variants) shows no 3'-end concentration: correct projections
are nearly uniform across deciles
({{ clinvar_residual_positions.correct_decile10_pct | dp(1) }}% in the 3'-most decile)
and the {{ clinvar_residual_positions.incorrect_binned | commas }} binned differences
track them with a mild mid-CDS excess and a depleted 3'-most decile
({{ clinvar_residual_positions.incorrect_decile10_pct | dp(1) }}%), consistent with
their representation and multi-mapping origins rather than positional alignment drift
(`bin_residual_positions.py`).

### Table S5: Injection cleaning benchmark, per category

Reproducible injection benchmark (`paper/scripts/inject_and_clean.py`): each
`clean_hgvs()` fix category is injected into a seeded sample of {{ cleaning.inject_sample_size | commas }}
clean, parseable public ClinVar c.HGVS strings (committed to the repository, seed 112,
at most 200 cases per category), and recovery is an exact string match to the known
canonical target. The LOVD columns come from `paper/scripts/lovd_head_to_head.py`, which
runs the same cases through the LOVD HGVS syntax checker (v1.2.2, local PHP CLI) under
the same scoring rule (Methods); "top-1" scores the checker's highest-confidence
suggested correction. The VariantValidator (VV) and Mutalyzer columns come from
`paper/scripts/vv_mutalyzer_head_to_head.py`, which runs the same cases through the two
services' public REST APIs ({{ vv_mutalyzer_comparison.service_date }}; Methods) under
the same rule, scoring the validated (VV) or corrected/normalized (Mutalyzer)
description. Production ops with no string-level injector (structure reconstruction,
empty-version dropping, provider-verified accession-prefix restoration) are absent by
design. Examples are synthesised from public `NM_000059.4` (BRCA2). Regenerate with the
commands in each script's docstring; totals also appear in
`paper/empirical_results/cleaning.csv`, `lovd_comparison.csv` and
`vv_mutalyzer_comparison.csv`.

| Injected error (example → target) | n | `clean_hgvs()` | LOVD top-1 | VV | Mutalyzer |
|---|---|---|---|---|---|
| Whitespace (` NM_000059.4: c.68del`) | 200 | 100% | 100% | 96.5% | 89.5% |
| Lowercased bases (`c.316g>a`) | 200 | 100% | 100% | 96.0% | 89.5% |
| Trailing protein suffix (`c.68del p.Arg100Ter`) | 200 | 100% | 91.5% | 0% | 0% |
| Gene wrapper with colon (`NM_000059.4:(BRCA2):c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Gene/transcript swapped (`BRCA2(NM_000059.4):c.68del`) | 200 | 100% | 49.5% | 94.5% | 0% |
| Surrounding quotes (`"NM_000059.4:c.68del"`) | 200 | 100% | 0% | 97.5% | 0% |
| Doubled colon (`NM_000059.4::c.68del`) | 200 | 100% | 0% | 99.0% | 0% |
| Unbalanced bracket (`(NM_000059.4:c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Separator typo (`NM_000059.4:c,68del`) | 200 | 100% | 100% | 0% | 0% |
| Doubled version dot (`NM_000059..4:c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Leading junk (`GRCh38.p2 NM_000059.4:c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Doubled kind (`NM_000059.4:c.c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Lowercased accession (`nm_000059.4:c.68del`) | 200 | 100% | 52.0% | 94.0% | 35.5% |
| Redundant del/dup count (`c.68_69del23`) | 77 | 100% | 0% | 0% | 0% |
| Missing accession underscore (`NM000059.4:c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Colon in accession prefix (`NM:_000059.4:c.68del`) | 200 | 100% | 0% | 0% | 0% |
| Uppercased mutation type (`c.68DEL`) | 142 | 100% | 100% | 94.4% | 0% |
| Dropped accession letter (`M_000059.4:c.68del`) | 200 | 100% | 0% | 0% | 0% |
| **Total** | **{{ lovd_comparison.n_cases | commas }}** | **{{ lovd_comparison.cdot_pct | dp(1) }}%** | **{{ lovd_comparison.lovd_top1_pct | dp(1) }}%** | **{{ vv_mutalyzer_comparison.vv_pct | dp(1) }}%** | **{{ vv_mutalyzer_comparison.mut_pct | dp(1) }}%** |

Weighted by the production rescue-op distribution (Results Table 2) the totals are
{{ lovd_comparison.cdot_weighted_pct | dp(1) }}% for `clean_hgvs()`,
{{ lovd_comparison.lovd_top1_weighted_pct | dp(1) }}% for LOVD top-1,
{{ vv_mutalyzer_comparison.vv_weighted_pct | dp(1) }}% for VariantValidator and
{{ vv_mutalyzer_comparison.mut_weighted_pct | dp(1) }}% for Mutalyzer; accepting the
target anywhere in LOVD's ranked correction list changes no case. On the
{{ lovd_comparison.originals_n | commas }} uncorrupted originals no tool falsely
corrects any input ({{ lovd_comparison.cdot_false_corrections | int }} for `clean_hgvs()`,
{{ lovd_comparison.lovd_false_corrections | int }} for LOVD,
{{ vv_mutalyzer_comparison.vv_false_corrections | int }} for VariantValidator,
{{ vv_mutalyzer_comparison.mut_false_corrections | int }} for Mutalyzer); LOVD flags
{{ lovd_comparison.lovd_flagged_invalid_pct | dp(1) }}% of them (intronic positions on
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

**[Tier 2].** Single-label classification of the 1,075 production queries (826 unique
strings) that still fail to parse after cleaning (Results, "Residual errors"), under a
fixed decision-tree taxonomy. Of the eight classes, the seven repair-relevant ones are
shown; the eighth was non-HGVS input (81 queries, 7.5%: pasted URLs, report templates,
or prose), excluded here as there is nothing in it for cleaning to repair. Counts and %
are of the 1,075 residual queries; examples are synthesised from public BRCA2
`NM_000059.4`. *(Tier 2; frozen constants from a deterministic run over the production
corpus.)*

| Class | Queries | What it is (*example*) |
|---|---|---|
| Truncated | 284 (26.4%) | cut off before a complete variant: `NM_000059.4:c.68_69` (range, no edit) |
| No reference | 277 (25.8%) | a bare variant body, no transcript/gene/accession: `c.68_69delAG` |
| Bad accession | 124 (11.5%) | misplaced or truncated version, or a missing prefix with no unique data match: `NM_000059/4:c.68del` (slash in place of the version dot) |
| Edit syntax | 113 (10.5%) | malformed or non-standard edit operation: `NM_000059.4:c.68AG>T` (multi-base reference in a substitution) |
| Trailing / concatenated | 85 (7.9%) | extra characters after a complete variant, or several run together: `NM_000059.4:c.68delAG;c.70A>G` |
| Grammar gap | 81 (7.5%) | legitimate HGVS the biocommons grammar rejects: `NM_000059.4:c.(67+1_68-1)_(70+1_71-1)del` (uncertain-range deletion) |
| Insertion (length only) | 30 (2.8%) | an insertion given as a base count, not a sequence: `NM_000059.4:c.68_69ins5` (position and length recoverable; inserted bases not) |

*Method and limitation:* classification was performed by a large language model (Claude
Opus 4, Anthropic; 2026-06-17) applying the shared decision tree to each unique string,
single-label and single-rater; no second-rater adjudication was done, so no inter-rater
agreement (κ) is reported. The taxonomy is version `v1`. Counts were refreshed on
2026-08-17 after the accession-repair cleaning rules landed: the 43 residual queries the
new rules rescue all sat in the Bad accession class (delta re-rated by Claude Fable 5,
Anthropic), which drops from 167 (14.9% of the previous 1,118-query residual) to 124;
the other classes are unchanged. Synthesised examples (from public `NM_000059.4` /
`NM_001754.5`) illustrate each class; no corpus string is reproduced.

### Figure S1: Positional drift along the CDS across version bumps

![](paper/figures/figure_s1_positional_drift.svg)

**Figure S1. Coordinate preservation by relative CDS position across consecutive
transcript version bumps** (GRCh38, seeded
{{ version_stability.sample_n | commas }}-accession sample per consortium; deciles of the
shared CDS, 1 = 5'-most tenth). **(A)** Conditioned on the partial-drift bumps
({{ positional_drift.refseq_partial_pairs }} RefSeq and
{{ positional_drift.ensembl_partial_pairs }} Ensembl pairs, the
{{ version_stability.refseq_partial_drift_pct | dp(1) }}% and
{{ version_stability.ensembl_partial_drift_pct | dp(1) }}% of bumps where some coding
bases move and others do not): preservation declines monotonically toward the 3' end,
from {{ positional_drift.refseq_partial_decile1_pct | dp(1) }}% to
{{ positional_drift.refseq_partial_decile10_pct | dp(1) }}% of coding bases for RefSeq
and {{ positional_drift.ensembl_partial_decile1_pct | dp(1) }}% to
{{ positional_drift.ensembl_partial_decile10_pct | dp(1) }}% for Ensembl, because a
partial drift keeps a 5' prefix intact up to its first alignment change and shifts the
bases downstream of it. **(B)** Unconditioned over every compared bump the curve is
essentially flat (RefSeq {{ positional_drift.refseq_all_decile1_pct | dp(1) }}% to
{{ positional_drift.refseq_all_decile10_pct | dp(1) }}%; Ensembl
{{ positional_drift.ensembl_all_decile1_pct | dp(1) }}% to
{{ positional_drift.ensembl_all_decile10_pct | dp(1) }}%): most drift is a whole-CDS
relocation, position independent and reliably flagged by the intrinsic-structure check
(Results R5). The positional effect is confined to the rare partial-drift tail, which is
what makes a 5' coding variant safer to substitute than a 3' one when a version must be
swapped. Produced by `compute_version_stability.py` (facts) and
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
`NM_`/`XM_`, `n.` → `NR_`/`XR_`) and applying the fix only when exactly one candidate
accession exists in the loaded data (Methods).
