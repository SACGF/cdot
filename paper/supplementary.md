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

Full-scale resolution of every RefSeq c.HGVS in ClinVar through cdot alone (GRCh38, cdot
0.2.33). Unlike the cdot-vs-UTA comparison in Results R1 (gated to a sample by UTA's
throughput), every variant is summarised here because only cdot is in the loop. The
projection is scored as a VCF coordinate (CHROM/POS/REF/ALT) against ClinVar's own VCF,
not as a g.HGVS string, so equivalent representations and ClinVar's tandem-repeat / identity
notations are not miscounted.

**Caveat.** ClinVar submissions are dominated by a handful of large (largely US) clinical
laboratories citing mostly current RefSeq versions, so this is a clean, public,
reproducible scale check that cdot resolves real variants at scale, not an unbiased sample
of the transcripts clinical labs use. The unbiased real-world complement is the Shariant
historical corpus (Results R1, Tier 2). The input is RefSeq-only (ClinVar's Name column is
RefSeq-centric), so there is no Ensembl row here.

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
