# Discussion

For biocommons/hgvs users, cdot removes UTA's constraints (a PostgreSQL database,
RefSeq only, limited history) with a drop-in provider that adds Ensembl,
T2T-CHM13v2.0, and the full release history. For the wider HGVS ecosystem, it recovers
the malformed and version-drifted descriptions real pipelines receive, with every
change reported for audit and a version substitution refused unless verified
coordinate-safe.

cdot also integrates Ensembl TARK: `EnsemblTarkDataProvider` is, to our knowledge, the
only client that exposes the Ensembl Transcript Archive through the biocommons/hgvs
data-provider interface, for pipelines that require Ensembl's own authoritative source.
Tools such as VariantValidator [@Freeman2018], built on biocommons/hgvs with a
self-hosted copy of UTA, and Mutalyzer [@Lefter2021], which uses its own normalisation
stack, are widely used to check and correct HGVS descriptions; cdot is complementary to
them, supplying the transcript-coordinate layer, not a validation service.

Nor is cdot the first tool to repair broken HGVS. VariantValidator corrects the
mistakes it can interpret [@Freeman2018]; Mutalyzer's Name Checker returned a corrected
description for ~{{ literature.mutalyzer_autocorrect_pct | dp(0) }}% of the unique
descriptions in its production logs [@Lefter2021]; the LOVD HGVS syntax checker
[@LovdHgvsChecker] suggests ranked corrections without needing a reference sequence;
and the ClinGen Allele Registry [@Pawliczek2018] normalises descriptions on
registration. For these tools repair happens inside validation or registration, usually
behind a web service. `clean_hgvs()` differs in placement: an offline function that
repairs the string before parsing, returns every change as an auditable `HGVSFix`, and
never breaks a description that already parsed.

The LOVD checker is the closest comparator (offline, no reference sequence), so we ran
it head-to-head with `clean_hgvs()` on the injection corpus
({{ lovd_comparison.n_cases | commas }} corrupted strings, checker
{{ lovd_comparison.lovd_version }}; scoring in Methods, per-category breakdown in
Supplementary Table S5). The corpus injects the error classes `clean_hgvs()` targets,
so `clean_hgvs()` recovers every case by construction and the comparison measures how
much of that territory a general syntax checker also covers. Weighted by the production
error mix, LOVD's top-ranked correction restored
{{ lovd_comparison.lovd_top1_weighted_pct | dp(0) }}% of cases
({{ lovd_comparison.lovd_top1_pct | dp(0) }}% unweighted), matching `clean_hgvs()` on
common single-token errors but not on free-text damage or cross-accession-family
repairs; neither tool altered any valid input.

The same corpus puts the sequence-aware services in context. Run over their public REST
APIs ({{ vv_mutalyzer_comparison.service_date }}; Methods), VariantValidator's validated
description recovered {{ vv_mutalyzer_comparison.vv_weighted_pct | dp(0) }}% of cases
weighted by the production error mix, and Mutalyzer
{{ vv_mutalyzer_comparison.mut_weighted_pct | dp(0) }}% (Supplementary Table S5). Each
repairs the damage it can interpret and rejects the rest as invalid: whitespace and
letter case for both, with VariantValidator also handling quotes, doubled colons,
accession re-casing and swapped gene and transcript forms at 94 to 99% per category,
capped by transcripts absent from its database rather than by the injected error.
Neither falsely corrected a valid input. These services answer a different question from
`clean_hgvs()`: they validate a description against the reference sequence, so a string
that cannot be parsed is an error to report, not text to repair. The roles are
converging, though. Since release 4.0.0 (June 2026) VariantValidator embeds the LOVD
syntax checker [@VariantValidator400], and its API responses for input it cannot parse
carry the checker's ranked suggestions alongside the validation error, although the
service does not apply them. The LOVD measurement above therefore also characterises
the syntax-repair layer surfaced by a current VariantValidator pipeline.

Beyond HGVS resolution, the JSON representation is useful in its own right. It parses
far faster than the GTF/GFF files it is built from, and the per-release JSON published
on the GitHub releases page is a faster-loading drop-in for the corresponding GTF/GFF.
The REST API returns transcript coordinates on demand and in batches, so thin clients
need not bundle large annotation downloads; Ensembl's public REST service covers only
Ensembl transcripts at the latest version, while cdot serves both consortia with
history. The format is also read outside Python: ferro-hgvs [@FerroHgvs], an
independent HGVS parser and normaliser written in Rust, loads cdot JSON directly as its
transcript-to-genome alignment source. A consumer in another language needs only a JSON
parser, not an HGVS library or the Python data-provider interface.

cdot separates unambiguous string cleaning, safe to apply automatically, from
heuristics that can be wrong (an adjacent transcript version, a canonical transcript
for a gene symbol), which are opt-in, never applied silently, and always reported as an
`HGVSFix` the caller can inspect or reject. The main limitation is currency: a
transcript version published after the most recent ingested release is not covered
until the dataset is regenerated. Automating that regeneration so each new RefSeq and
Ensembl release is ingested and published without manual intervention is a planned
improvement.
