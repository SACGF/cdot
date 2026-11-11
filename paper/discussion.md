# Discussion

cdot is, first, a general-purpose source of versioned transcript-to-genome alignments:
RefSeq and Ensembl together, across GRCh37, GRCh38, and T2T-CHM13v2.0, with the full
release history, in a JSON format that loads in seconds and can be read from any
language, not only Python. The per-release JSON on
the GitHub releases page is a faster-loading drop-in for the corresponding GTF/GFF, and
the REST API returns transcript coordinates on demand and in batches, so thin clients
need not bundle large annotation downloads. For biocommons/hgvs users specifically, it
removes UTA's constraints (a PostgreSQL database, RefSeq only, limited history) with a
drop-in provider; for pipelines that require Ensembl's own authoritative source,
`EnsemblTarkDataProvider` exposes the Ensembl Transcript Archive (TARK) through the
same interface.

cdot is not the first tool to repair broken HGVS. VariantValidator [@Freeman2018],
built on biocommons/hgvs with a self-hosted copy of UTA, corrects the mistakes it can
interpret; Mutalyzer's Name Checker returned a corrected description for
~{{ literature.mutalyzer_autocorrect_pct | pct(0) }} of the unique descriptions in its
production logs [@Lefter2021]; the LOVD HGVS syntax checker [@LovdHgvsChecker] suggests
ranked corrections without needing a reference sequence; and the ClinGen Allele
Registry [@Pawliczek2018] normalises descriptions on registration. For these tools
repair happens inside validation or registration, usually behind a web service.
`clean_hgvs()` differs in placement: an offline function that repairs the string before
parsing, returns every change as an auditable `HGVSFix`, and never breaks a description
that already parsed. cdot is complementary to these services, supplying the
transcript-coordinate layer, not validation.

The LOVD checker is the closest comparator (offline, no reference sequence), so we ran {# lint-ok: 'no reference sequence' describes the checker #}
it head-to-head with `clean_hgvs()` on the injection corpus
({{ lovd_comparison.n_cases | commas }} corrupted strings, checker
{{ lovd_comparison.lovd_version }}; scoring in Methods, per-category breakdown in
Supplementary Table S5). The corpus injects the error classes `clean_hgvs()` targets,
{% assert cleaning.inject_total_recovered == cleaning.inject_total_attempted, "'clean_hgvs() recovers every case by construction'" %}
so `clean_hgvs()` recovers every case by construction and the comparison measures how
much of that territory a general syntax checker also covers. Weighted by the production
error mix, LOVD's top-ranked correction restored
{{ lovd_comparison.lovd_top1_weighted_pct | pct(0) }} of cases
({{ lovd_comparison.lovd_top1_pct | pct(0) }} unweighted), matching `clean_hgvs()` on
common single-token errors but not on free-text damage or cross-accession-family
{% assert lovd_comparison.cdot_false_corrections == 0 and lovd_comparison.lovd_false_corrections == 0, "'neither tool altered any valid input'" %}
repairs; neither tool altered any valid input.

The same corpus puts the sequence-aware services in context. Run over their public REST {# lint-ok: 'over' is a preposition #}
APIs ({{ vv_mutalyzer_comparison.service_date }}; Methods), VariantValidator's validated
description recovered {{ vv_mutalyzer_comparison.vv_weighted_pct | pct(0) }} of cases
weighted by the production error mix, and Mutalyzer
{{ vv_mutalyzer_comparison.mut_weighted_pct | pct(0) }} (Supplementary Table S5). Each
{% assert vv_mutalyzer_comparison.vv_false_corrections == 0 and vv_mutalyzer_comparison.mut_false_corrections == 0, "'neither falsely corrected a valid input' (VariantValidator, Mutalyzer)" %}
repairs the damage it can interpret and rejects the rest as invalid; neither falsely
corrected a valid input. These services answer a different question from
`clean_hgvs()`: they validate a description against the reference sequence, so a string
that cannot be parsed is an error to report, not text to repair. The roles are
converging, though. Since release 4.0.0 (June 2026) VariantValidator embeds the LOVD
syntax checker [@VariantValidator400], and its responses for input it cannot parse
carry the checker's ranked suggestions alongside the validation error, although the
service does not apply them. The LOVD measurement above therefore also characterises
the syntax-repair layer surfaced by a current VariantValidator pipeline.

The general-purpose value of the format shows in reuse beyond Python and beyond
biocommons/hgvs. Where Ensembl's own REST service serves only Ensembl transcripts, cdot
serves both consortia with history; and ferro-hgvs [@FerroHgvs],
an independent HGVS parser and normaliser written in Rust, loads cdot JSON directly as
its transcript-to-genome alignment source, needing only a JSON parser.

cdot separates unambiguous string cleaning, safe to apply automatically, from
heuristics that can be wrong (an adjacent transcript version, a canonical transcript
for a gene symbol), which are opt-in and always reported as an `HGVSFix` the caller can
inspect or reject. The main limitation is currency: a transcript version published
after the most recent ingested release is not covered until the dataset is regenerated.
Automating that regeneration so each new RefSeq and Ensembl release is ingested and
published without manual intervention is a planned improvement. Two further limits are
visible in the results. A share of the residue in every corpus is valid HGVS the
biocommons/hgvs grammar does not parse (repeat and allele notation; Supplementary
Tables S6 and S8), a parser limitation cdot does not address. And the full-ClinVar
scale check inherits ClinVar's submitter bias toward current RefSeq versions
(Supplementary Table S4), so it is a reproducible check at scale rather than an
unbiased sample of the transcripts clinical laboratories use.
