# Discussion

cdot fills two gaps in the HGVS resolution stack. The first is specific to biocommons/hgvs:
the library needs a transcript data provider, and the standard one (UTA) requires a
PostgreSQL database, omits Ensembl, and keeps a limited set of historical versions. cdot is
a drop-in replacement that needs no database, covers both RefSeq and Ensembl, adds
T2T-CHM13v2.0, and retains the full release history. The second gap is general to the HGVS
ecosystem, independent of any one library: the descriptions that clinical and research
systems actually receive are often malformed or cite a transcript version that has since
been retired, and most tools simply reject them. cdot's string cleaning and safe version
fallback recover many of these.

The practical payoff is robustness where it is otherwise hard to achieve. A clinical
laboratory or research group collecting HGVS from report PDFs, spreadsheets, and free-text
search boxes can run `clean_hgvs()` as a pre-pass, so that whitespace, casing, and
copy/paste damage no longer break otherwise-valid descriptions, with every change reported
for audit rather than applied silently. A variant cited against a transcript version that
has since been retired can still resolve through the opt-in version fallback, which
substitutes the nearest available version only when a build-independent check confirms the
substitution does not move the variant. Together these turn input that would otherwise fail
into resolved coordinates without asking the user to reformat anything, which is what makes
cdot worth dropping into a working clinical or research pipeline.

cdot also integrates Ensembl TARK, the Ensembl transcript archive: its
`EnsemblTarkDataProvider` is, to our knowledge, the only client that exposes TARK through
the biocommons/hgvs data-provider interface. This lets a pipeline draw transcript data
straight from Ensembl's own official REST source when that authoritative provenance is
required, without bespoke code. Beyond what TARK's REST service offers, cdot's own data
adds RefSeq coverage, fully offline operation, and support for T2T-CHM13v2.0. Tools such
as VariantValidator [@Freeman2018], built on the biocommons/hgvs library with a
self-hosted copy of UTA, and Mutalyzer [@Lefter2021], which uses its own independent
normalisation stack and retrieves transcripts directly from NCBI and Ensembl, are widely
used to check HGVS correctness; cdot is complementary to them, supplying the
transcript-coordinate layer rather than validating descriptions.

Beyond HGVS resolution, the JSON representation is useful in its own right. It parses far
faster than the GTF/GFF files it is built from and loads trivially over HTTP, so cdot
doubles as a lightweight, queryable gene/transcript reference. We publish the per-release
JSON for each annotation version on the GitHub releases page, where a single file is a
much faster-loading drop-in for the corresponding GTF/GFF. Because the REST API returns
JSON quickly and in batches, downstream software can query transcript coordinates on
demand instead of bundling large annotation downloads, which is convenient for thin
clients, and for AI agents that call the API directly. Ensembl offers a public REST service, but only
for Ensembl transcripts and only at the latest version of each; cdot serves both RefSeq
and Ensembl and retains historical versions.

By design, cdot separates unambiguous string cleaning, which is safe to apply
automatically, from heuristics that can be wrong (choosing an adjacent transcript
version, or mapping a gene symbol to a canonical transcript), which are opt-in, never
applied silently, and always reported as an `HGVSFix` the caller can inspect or reject.
A general limitation is that resolution is only as current as the ingested annotation
releases: a transcript version published after the most recent ingested release is not
covered until the dataset is regenerated. Automating that regeneration so each new RefSeq
and Ensembl release is ingested and published as a data release without manual
intervention is a planned improvement.

