# Discussion

cdot fills two gaps in the HGVS resolution stack. The first is specific to biocommons/hgvs:
the library needs a transcript data provider, and the standard one (UTA) requires a
PostgreSQL database, omits Ensembl, and keeps a limited set of historical versions. cdot is
a drop-in replacement that needs no database, covers both RefSeq and Ensembl, adds
T2T-CHM13v2.0, and retains the full release history. The second gap is general to the HGVS
ecosystem, independent of any one library: the descriptions that clinical and research
systems actually receive are often malformed or cite a transcript version that has since
been retired, and most tools reject them. cdot's string cleaning and safe version
fallback recover many of these.

The practical benefit shows up in pipelines that handle messy input. A
clinical laboratory or research group collecting HGVS from report PDFs, spreadsheets, and
free-text search boxes can run `clean_hgvs()` as a pre-pass, so that whitespace, casing,
and copy/paste damage no longer break otherwise-valid descriptions, with every change
reported for audit rather than applied silently. A variant cited against a transcript version that
has since been retired can still resolve through the opt-in version fallback, which
substitutes the nearest available version only when a build-independent check confirms the
substitution does not move the variant. These two steps turn input that would otherwise
fail into resolved coordinates, without asking the user to reformat anything.

cdot also integrates Ensembl TARK, the Ensembl transcript archive: its
`EnsemblTarkDataProvider` is, to our knowledge, the only client that exposes TARK through
the biocommons/hgvs data-provider interface. This lets a pipeline draw transcript data
straight from Ensembl's own official REST source when that authoritative provenance is
required, without bespoke code. Beyond what TARK's REST service offers, cdot's own data
adds RefSeq coverage, fully offline operation, and support for T2T-CHM13v2.0. Tools such
as VariantValidator [@Freeman2018], built on the biocommons/hgvs library with a
self-hosted copy of UTA, and Mutalyzer [@Lefter2021], which uses its own independent
normalisation stack and retrieves transcripts directly from NCBI and Ensembl, are widely
used to check and correct HGVS descriptions; cdot is complementary to them, supplying
the transcript-coordinate layer rather than a validation service.

cdot is also not the first tool to repair broken HGVS rather than merely reject it.
VariantValidator automatically corrects the mistakes it can interpret, including
intronic positions described relative to the wrong exon boundary [@Freeman2018].
Mutalyzer's Name Checker returned a corrected description for
~{{ literature.mutalyzer_autocorrect_pct | dp(0) }}% of the unique descriptions in its
production logs [@Lefter2021]. The LOVD HGVS syntax checker [@LovdHgvsChecker] checks
syntax without needing a reference sequence and suggests corrections for invalid
descriptions, ranked by likelihood; the ClinGen Allele Registry [@Pawliczek2018]
normalises the descriptions it registers to canonical allele identifiers; and the
biocommons/hgvs parser itself tolerates a few common deviations, such as a gene symbol
in parentheses after the accession. For these tools repair happens as part of validation
or registration, usually behind a web service. `clean_hgvs()` differs in role and
placement: it is an offline Python function that repairs the string before parsing, so
anything that consumes the string downstream benefits; every change is returned as an
`HGVSFix` the caller can audit; and it guarantees no regressions, never breaking a
description that already parsed.

The LOVD checker is the closest comparator, being the one tool in this group that runs
offline without a reference sequence, so we ran it head-to-head with `clean_hgvs()` on
the reproducible injection corpus ({{ lovd_comparison.n_cases | commas }} corrupted
strings, checker {{ lovd_comparison.lovd_version }}; scoring in Methods, per-category
breakdown in Supplementary Table S5). The corpus
injects the error classes `clean_hgvs()` targets, so `clean_hgvs()` recovers every case
by construction and the comparison measures how much of that territory a general syntax
checker also covers, not the overall quality of either tool. Weighted by the production
error mix, LOVD's top-ranked correction restored
{{ lovd_comparison.lovd_top1_weighted_pct | dp(0) }}% of cases
({{ lovd_comparison.lovd_top1_pct | dp(0) }}% unweighted across categories), matching
`clean_hgvs()` on the common single-token errors: whitespace, letter case, separator
typos and trailing protein annotations. The rest splits into inputs that must be treated
as free text rather than a malformed variant description (surrounding quotes, unbalanced
brackets, leading assembly text), which sit outside the checker's intended input, and
repairs applied to only one accession family (it restores a swapped gene and transcript
for RefSeq but not Ensembl accessions, and re-cases a lowercase Ensembl accession but
not a RefSeq one). Neither tool altered any valid input, though LOVD flags intronic
descriptions on a transcript reference
({{ lovd_comparison.lovd_flagged_invalid_pct | dp(0) }}% of the valid originals) as
requiring a genomic reference, a deliberate design position rather than a defect.

Beyond HGVS resolution, the JSON representation is useful in its own right. It parses far
faster than the GTF/GFF files it is built from and loads trivially over HTTP, so cdot
doubles as a lightweight, queryable gene/transcript reference. We publish the per-release
JSON for each annotation version on the GitHub releases page, where a single file is a
much faster-loading drop-in for the corresponding GTF/GFF. Because the REST API returns
JSON quickly and in batches, downstream software can query transcript coordinates on
demand instead of bundling large annotation downloads, which is convenient for thin
clients, and for AI agents that call the API directly. Ensembl offers a public REST service, but only
for Ensembl transcripts and only at the latest version of each; cdot serves both RefSeq
and Ensembl and retains historical versions. cdot's JSON is also read outside Python:
ferro-hgvs [@FerroHgvs], an independent HGVS parser and normaliser written in Rust, loads
cdot JSON directly as its transcript-to-genome alignment source, and uses it to provide
Ensembl support. Because the format is plain JSON with a documented schema, a consumer in
another language needs only a JSON parser, not an HGVS library or the Python data-provider
interface.

cdot separates unambiguous string cleaning, which is safe to apply
automatically, from heuristics that can be wrong (choosing an adjacent transcript
version, or mapping a gene symbol to a canonical transcript), which are opt-in, never
applied silently, and always reported as an `HGVSFix` the caller can inspect or reject.
A general limitation is that resolution is only as current as the ingested annotation
releases: a transcript version published after the most recent ingested release is not
covered until the dataset is regenerated. Automating that regeneration so each new RefSeq
and Ensembl release is ingested and published as a data release without manual
intervention is a planned improvement.

