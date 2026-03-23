# Methods

*Title in Bioinformatics Original Papers: "System and Methods" or "Implementation"*
*Target: ~1200–1500 words.*

---

## 2.1 Data sources

cdot integrates transcript annotation from three sources, applied in priority order:

1. **RefSeq GFF3** — NCBI RefSeq annotation releases for GRCh37 (historical releases spanning [year range]) and GRCh38 (releases [range]) are processed to extract transcript-level alignment coordinates, CDS boundaries, and alignment gap strings. Multiple historical releases are included to maximise coverage of versioned transcript accessions that appear in historical clinical HGVS strings.

2. **Ensembl GTF** — Ensembl gene annotation releases 75–81 (GRCh37) and 76–[115] (GRCh38), plus T2T-CHM13v2.0 releases, are parsed from GTF format. Ensembl GTF files provide transcript coordinates but express alignments without gap strings; gap information for Ensembl transcripts is derived from the RefSeq GFF3 where accession overlap exists.

3. **UTA CSV** — The Universal Transcript Archive (UTA) CSV dump is ingested as a baseline source providing transcript alignments for accessions that predate the oldest available GFF3/GTF files. UTA data is overridden by official GFF3/GTF data where available.

For T2T-CHM13v2.0, both Ensembl and RefSeq annotation files are processed using the same pipeline.

Gene-symbol-to-HGNC-ID mappings are derived from GENCODE HGNC files.

## 2.2 JSON format

[Describe the JSON schema: top-level structure, transcript dict, genome-build-specific fields, 6-tuple exon representation, gap string, MANE/canonical tags, schema versioning. Reference the schema version (v0.2.33 at time of writing).]

Key design decisions:
- **Exon 6-tuple**: `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]` — gap field stores the GFF3 gap string (`"M196 I1 M61"`) which is converted on read to HGVS CIGAR format with inverted I/D operators
- **Shared vs build-specific fields**: transcript metadata (gene, biotype, sequence accession) stored once; coordinate data (exon positions, strand, contig) stored per genome build
- **Schema versioning**: JSON files carry a `schema_version` field; clients check this on load
- **Compression**: files are gzip-compressed JSON; the Python client transparently handles both `.json` and `.json.gz`

## 2.3 Data generation pipeline

[Describe the Snakemake pipeline: input GFF3/GTF files, parsing steps, merge strategy, output JSON.gz files, how source priority is applied.]

Classes:
- `GFFParser` (base) → `GTFParser` (Ensembl) / `GFF3Parser` (RefSeq) — parse annotation files into a canonical internal transcript representation
- Merge step applies source priority: RefSeq GFF3 > Ensembl GTF > UTA for overlapping accessions
- Output is split by genome build; files are named `cdot.[version].refseq.[build].json.gz` and `cdot.[version].ensembl.[build].json.gz`

## 2.4 Access methods

**Local JSON**: The `LocalDataProvider` / `JSONDataProvider` classes load a JSON.gz file into memory on initialisation. Transcript lookup is O(1) dictionary access. Typical load time: [X] seconds for GRCh38 RefSeq (~[N] MB uncompressed).

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves JSON data over HTTP at cdotlib.org. The `RESTDataProvider` client fetches individual transcripts on demand. This is analogous to SeqRepo REST [Hart 2020] for sequence data — suitable for occasional lookups where the full file download is not desired.

**Ensembl TARK**: An `EnsemblTarkDataProvider` wraps the Ensembl Transcript Archive (TARK) REST API as an alternative Ensembl data source.

## 2.5 biocommons/hgvs integration

cdot implements the `biocommons.hgvs.dataproviders.interface.Interface` abstract class. The methods `get_tx_info`, `get_tx_exons`, `get_tx_seq`, `get_tx_mapping_options`, and related accessors are implemented using data from the JSON files. Alignment gap strings are converted from GFF3 format (`"M196 I1 M61"`) to HGVS CIGAR format on read, with insertion/deletion operators inverted to match the HGVS convention (gap relative to transcript rather than genome).

Sequence data is provided via `SeqRepo` (biocommons) or cdot's `FastaSeqFetcher` (using local genome FASTA files), allowing fully offline operation.

## 2.6 PyHGVS integration

cdot provides a `PyHGVSTranscriptFactory` that reads from the same JSON data and returns `pygr`-compatible transcript objects. This adds alignment gap support to PyHGVS, which lacks it natively.

## 2.7 Canonical transcript selection

[Describe the `CanonicalTranscriptSelector` API once implemented (#36). Fields stored: `mane_select`, `mane_plus_clinical`, `refseq_select`, `ensembl_canonical`. Priority order for lookup. Gene-name → transcript → HGVS conversion use case.]

---

*Notes:*
- *Expand 2.2 with exact JSON schema field names once finalised*
- *Add Snakemake pipeline detail (number of rules, input/output)*
- *Section 2.7 requires #36 to be implemented before submission*
