# Discussion

*Target: ~300 words. Positioning, limitations, future. Do not repeat implementation
detail.*

---

cdot fills a specific gap in the biocommons HGVS stack: comprehensive, versioned,
multi-source transcript coordinate data accessible without database infrastructure.

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
covered until the dataset is regenerated. In practice cdot tracks new releases closely
and is regenerated regularly (the current data already incorporates RefSeq RS_2025_08 and
Ensembl 115), so this is normally a short window after each release rather than a lasting
gap.

