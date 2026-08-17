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
