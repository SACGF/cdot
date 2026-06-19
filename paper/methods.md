# Methods

*Full Original Paper (Bioinformatics): Methods section. See `claude/paper_plan.md`.*
*Target: ~1,200–1,600 words. Subsections: Data sources & generation · JSON format ·
String cleaning (`clean_hgvs()`) · Access & client libraries · Canonical transcript
selection.*

---

## Data sources and generation

cdot draws on three transcript-annotation sources, merged into a single dataset per
build. We download the complete run of historical releases directly from the RefSeq and
Ensembl FTP sites rather than only the current annotation, because resolving historical
HGVS strings depends on it: a clinical report from 2015 may reference an NM_ version that
was retired in a later RefSeq release but is still cited in patient records and ClinVar
submissions. The releases are combined in chronological order, with each newer release
overwriting older entries for the same transcript version. This ordering matters beyond
simple recency: older files are sometimes less complete than later ones (for example,
early RefSeq GFF3 releases omitted the transcript/genome alignment-gap CIGAR strings that
later releases include), so letting newer data win improves the merged result rather than
merely keeping the latest version.

1. **RefSeq GFF3**: {{ sources.refseq_grch37_releases | int }} GRCh37,
   {{ sources.refseq_grch38_releases | int }} GRCh38, and
   {{ sources.refseq_t2t_releases | int }} T2T-CHM13v2.0 historical NCBI annotation
   releases (Supplementary Table S1), running from the earliest archived annotation
   releases through to the dated RS_2025_08 release (August 2025).

2. **Ensembl GTF**: {{ sources.ensembl_grch37_releases | int }} GRCh37,
   {{ sources.ensembl_grch38_releases | int }} GRCh38, and
   {{ sources.ensembl_t2t_releases | int }} T2T-CHM13v2.0 releases, spanning Ensembl
   releases 81 to 115 (Supplementary Table S2). We use Ensembl GTF rather than GFF3
   because the GFF3 files omit transcript
   protein versions. Ensembl coverage adds
   {{ coverage.ensembl_unique_count | commas }} transcript accessions absent from
   RefSeq.

3. **UTA**: the Universal Transcript Archive is included because UTA computes its own
   transcript-to-genome alignments, so it holds alignments for some accessions that were
   never published in any RefSeq GFF3 or Ensembl GTF release; these fill gaps for
   accessions predating the oldest available annotation files and are overridden by
   official annotation where available.

A Snakemake pipeline parses each source using `GFF3Parser` or `GTFParser` (HTSeq-based),
normalises contig names via bioutils, extracts CDS boundaries from start/stop codon
features, and serialises to gzip-compressed JSON. Gene-level metadata (symbols, HGNC
IDs, and biotypes) is drawn from the HGNC dataset and NCBI RefSeq gene-info files and
attached to each transcript.

## JSON format

Each transcript entry stores shared metadata (`gene_name`, `hgnc`, `biotype`, transcript
accession and version) plus a `genome_builds` dict keyed by assembly name (e.g.
`"GRCh38"`). The build is always a dict key, even in a single-build file, so the format
is identical whether a file stands in for one GTF/GFF for a single assembly or carries
several builds at once: the same structure lets a single-build file act as a drop-in
GTF/GFF replacement and lets the REST API return a transcript's GRCh37, GRCh38, and
T2T-CHM13v2.0 alignments together in one response, with no reshaping when builds are
combined. Build-specific fields include `contig`, `strand`, `cds_start`, `cds_end`,
and `exons`, a list of 6-tuples `[alt_start, alt_end, exon_id, cds_start, cds_end,
gap]`. The `gap` field stores the GFF3 gap string (e.g. `"M196 I1 M61"`) for transcripts
with indels relative to the genome; the Python client converts this to HGVS CIGAR format
(with I/D operators inverted per HGVS convention) at query time. Schema versioning
(`schema_version` at root level) allows clients to reject incompatible files on load.

Canonical transcript tags (`mane_select`, `mane_plus_clinical`, `refseq_select`, and
`ensembl_canonical`) are stored as build-specific fields where present. MANE Select
covers >{{ literature.mane_coverage_pct | dp(0) }}% of protein-coding genes
[@Morales2022] and is available for GRCh38; Ensembl canonical tags are available for
GRCh37 and GRCh38.

## String cleaning (`clean_hgvs()`)

Resolving an HGVS string is only possible once it parses, yet a large fraction of the
descriptions that reach a variant-curation platform do not. They arrive from clinical
report PDFs, spreadsheets, literature, and free-text search boxes, and accumulate
formatting errors along the way: copy/paste whitespace and quotes, lost casing,
transposed punctuation, a gene symbol and transcript accession swapped, or a trailing
protein annotation. The intended variant is usually unmistakable to a human reader, but
such errors are fatal to a strict grammar. The conventional response is to validate the
string: report it as invalid and stop. cdot instead attempts to repair it.

`clean_hgvs()` is a pure string operation. It takes an HGVS string and returns the
repaired string together with a list of structured `HGVSFix` records describing every
change. It requires no genome build, no sequence, no HGVS parser, and no data provider,
so it can run as a pre-pass anywhere: in a search box, a batch importer, or ahead
of any parser. It is independent of which downstream library (biocommons/hgvs or
PyHGVS) ultimately consumes the result.

Cleaning is an ordered pipeline of single-purpose operations applied in a fixed
canonical order. Each operation inspects the string, makes at most one class of change,
and records an `HGVSFix` if it fired; operations are pure and order-dependent (for
example, stray leading characters and surrounding whitespace are stripped before structural
punctuation is examined, and casing is normalised before the gene/transcript
relationship is resolved). Callers may restrict cleaning to a subset of operations, such
as an allowlist for conservative use, but selection only filters the pipeline and does
not reorder it, so the canonical order is preserved regardless of the subset chosen. The
operations group into:

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
  `NM_000059.4(BRCA2):c.…`). The standard HGVS gene-in-parentheses form
  (`NM_000059.4(BRCA2):c.…`) is itself valid and is parsed by biocommons/hgvs unchanged;
  these operations target only the non-standard variants of it: accession and gene
  transposed, or the separating colon or dot missing.

Each repair is reported as an `HGVSFix` carrying a severity (`WARNING` for an
unambiguous correction), a stable machine-readable code, a human-readable message, and
the before/after values, so a caller can surface, log, or audit exactly what was
changed instead of silently mutating user input. After the repair pipeline, a
validation pass flags problems that cleaning cannot fix as `ERROR`-level fixes (no
colon at all, a bare variant body with no reference sequence, or an insertion given an
integer length instead of the inserted sequence); these mark incomplete input rather
than formatting noise. By default cleaning never raises and always returns its best
attempt, though callers may opt into raising on the first error.

A separate, opt-in helper, `get_best_transcript_version()`, addresses
transcript-version drift. Unlike the cleaning operations above, which are unambiguous
formatting corrections, substituting a different transcript version is a heuristic that
can be wrong (a coordinate may shift, or the variant may not exist, in the substituted
version), so it is never applied automatically. When the requested accession version is
absent from a caller's data, it returns the best available adjacent version (under a
configurable up-then-down, closest, or latest strategy) as a reported `HGVSFix`, and the
caller decides whether to accept the change and rewrite the string, which keeps control
of HGVS compatibility with the user. Before a substitution is accepted, cdot checks
whether it is coordinate-safe using a build-independent test of the two versions' CDS
structure (Results R5); the strategy chooses the candidate version, the safety check
decides whether it can be used.

## Access and client libraries

**Local JSON**: `JSONDataProvider` loads a JSON.gz file into memory on initialisation
(typically ~{{ benchmark.grch38_load_time_s | dp(0) }} seconds for GRCh38 RefSeq),
lazily building interval trees for region queries on first use and dictionaries for
transcript and gene lookup. Transcript retrieval is then O(1), giving throughput of
{{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }}
transcripts/second (Results, Table 2).

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves the same JSON data
at cdotlib.org. `RESTDataProvider` fetches transcripts one request per transcript
version, suitable for occasional lookups without downloading the full file, and can
warm its cache with a single batched `prefetch()` request (Results R6).

**biocommons/hgvs integration**: cdot implements
`biocommons.hgvs.dataproviders.interface.Interface`, providing `get_tx_info`,
`get_tx_exons`, `get_tx_seq`, `get_tx_mapping_options`, and related methods. This is the
same data-provider interface implemented by other tools such as hgvs-weaver
[@HgvsWeaver], so cdot serves as a drop-in replacement for UTA in any biocommons/hgvs
pipeline:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
hdp = JSONDataProvider(["cdot.0.2.33.refseq.grch38.json.gz"])
```

Sequence data is supplied by SeqRepo [@Hart2020] or cdot's `FastaSeqFetcher` (local genome FASTA),
enabling fully offline operation without SeqRepo's installation overhead.
`FastaSeqFetcher` reconstructs transcript sequence by splicing the transcript's exon
ranges out of the genome FASTA. This reproduces an Ensembl transcript exactly, but is not
guaranteed to match a RefSeq transcript, whose curated sequence can differ from the
genome at a small number of positions. In practice variants are usually called against a
genome-mapped read alignment, so any such discrepancy lies in the reference the variant
was already described against. Where an exact RefSeq transcript sequence is needed, SeqRepo
supplies it when it has the accession, but it does not hold every version cdot covers.
`ChainedSeqFetcher` addresses this by trying several sequence sources in a caller-defined
order, falling through to the next when one cannot supply a sequence, so a pipeline can,
for example, prefer SeqRepo and fall back to `FastaSeqFetcher` (or the reverse) to
maximise the fraction of transcripts it can serve.

**PyHGVS integration**: `JSONPyHGVSTranscriptFactory` provides a transcript factory for
the Counsyl PyHGVS library, exposing the same cdot transcript data to PyHGVS-based
pipelines. PyHGVS is no longer actively maintained, so new development targets the
biocommons/hgvs path; the PyHGVS factory is retained for legacy compatibility.

## Canonical transcript selection

cdot's `CanonicalTranscriptSelector` maps gene symbols to MANE Select (or MANE Plus
Clinical) transcript accessions, enabling gene-name HGVS lookup, a common requirement
in clinical reporting pipelines that receive a gene name rather than a transcript
accession. For example, a pipeline given `BRCA1:c.68_69del` with no transcript can look
up the MANE Select accession for BRCA1 (`NM_007294.4`) and resolve
`NM_007294.4(BRCA1):c.68_69del` against the genome. To our knowledge this is the first
HGVS data provider to expose programmatic canonical transcript selection aligned to the
MANE standard [@Morales2022; @Wright2023]. Picking a single canonical transcript for a
gene is not unambiguous: the clinically preferred transcript can differ from MANE Select,
and a gene may have more than one transcript of interest. The selection is intended as a
sensible default for search and gene-name lookup, not an authoritative clinical choice,
and callers should apply their own discretion where it matters.

---

*Notes:*
- *Fill all bracketed numbers before submission (run compute_coverage.py,
  compute_benchmark.py, and the ClinVar benchmark). Coverage/ClinVar prose now lives in
  `results.md` R1; Methods describes how, Results reports the numbers.*
- *String-cleaning subsection: examples synthesised from public NM_000059.4 (BRCA2) /
  NM_001754.5 (RUNX1) only, never corpus strings (CLAUDE.md / plan §3).*
- *Gap-correctness framing removed per feedback: cdot stores per-exon gaps (JSON-format
  detail) but applying them is a downstream-library concern, not promoted as a cdot
  result. PyHGVS paragraph trimmed to legacy-compatibility only.*
- *Canonical transcript section requires #36 to be implemented*
- *If over word count, the UTA CSV baseline source paragraph can be moved to
  supplementary*
