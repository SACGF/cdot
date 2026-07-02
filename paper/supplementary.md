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
| GRCh37 | 75–81 | ... |
| GRCh38 | 76–115 | ... |
| T2T-CHM13v2.0 | ... | ... |

### Table S3: JSON schema fields (v0.2.33)

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
  `transcript(GENE):c.…`; collapsing a doubled kind token (`c.c.` → `c.`); and fixing a
  comma or colon used in place of the kind dot, or a period used in place of a
  substitution `>`.
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
