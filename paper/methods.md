# Methods

*Full Original Paper (Bioinformatics) — Methods section. See `claude/paper_plan.md`.*
*Target: ~1,200–1,600 words. Subsections: Data sources & generation · JSON format ·
String cleaning (`clean_hgvs()`) · Access & client libraries · Canonical transcript
selection.*

---

## Data sources and generation

cdot processes transcript annotation from three sources, applied in priority order:

1. **RefSeq GFF3** — multiple historical NCBI annotation releases for GRCh37 and GRCh38
   (see Supplementary Table S1). Ingesting multiple releases is essential for resolving
   historical HGVS strings: a clinical report from 2015 may reference an NM_ version
   that was retired in subsequent RefSeq releases but is still cited in patient records
   and ClinVar submissions.

2. **Ensembl GTF** — releases 75–81 (GRCh37) and 76–115 (GRCh38), plus T2T-CHM13v2.0
   releases (Supplementary Table S2). Ensembl coverage adds
   {{ coverage.ensembl_unique_count | commas }} transcript accessions absent from
   RefSeq.

3. **UTA CSV** — the Universal Transcript Archive CSV dump provides alignments for
   accessions predating the oldest available GFF3/GTF files; overridden by official
   annotation where available.

A Snakemake pipeline parses each source using `GFF3Parser` or `GTFParser` (HTSeq-based),
normalises contig names via bioutils, extracts CDS boundaries from start/stop codon
features, and serialises to gzip-compressed JSON. Gene-symbol-to-HGNC-ID mappings come
from GENCODE HGNC files.

## JSON format

Each transcript entry stores shared metadata (`gene_name`, `hgnc`, `biotype`, transcript
accession and version) plus a `genome_builds` dict keyed by assembly name (e.g.
`"GRCh38"`). Build-specific fields include `contig`, `strand`, `cds_start`, `cds_end`,
and `exons` — a list of 6-tuples `[alt_start, alt_end, exon_id, cds_start, cds_end,
gap]`. The `gap` field stores the GFF3 gap string (e.g. `"M196 I1 M61"`) for transcripts
with indels relative to the genome; the Python client converts this to HGVS CIGAR format
(with I/D operators inverted per HGVS convention) at query time. Schema versioning
(`schema_version` at root level) allows clients to reject incompatible files on load.

Canonical transcript tags — `mane_select`, `mane_plus_clinical`, `refseq_select`, and
`ensembl_canonical` — are stored as build-specific fields where present. MANE Select
covers >{{ literature.mane_coverage_pct | dp(0) }}% of protein-coding genes
[@Morales2022] and is available for GRCh38; Ensembl canonical tags are available for
GRCh37 and GRCh38.

## String cleaning (`clean_hgvs()`)

Resolving an HGVS string is only possible once it parses, yet a large fraction of the
descriptions that reach a variant-curation platform do not. They arrive from clinical
report PDFs, spreadsheets, literature, and free-text search boxes, and carry the
formatting damage of that journey: copy/paste whitespace and quotes, lost casing,
transposed punctuation, a gene symbol and transcript accession swapped, or a trailing
protein annotation. These are not biological ambiguities — the intended variant is
usually unmistakable to a human reader — but they are fatal to a strict grammar. The
conventional response is to *validate*: report the string as invalid and stop. cdot
instead attempts to *repair* it.

`clean_hgvs()` is a pure string operation. It takes an HGVS string and returns the
repaired string together with a list of structured `HGVSFix` records describing every
change. It requires no genome build, no sequence, no HGVS parser, and no data provider,
so it can run as a cheap pre-pass anywhere — in a search box, a batch importer, or ahead
of any parser — and is independent of which downstream library (biocommons/hgvs or
PyHGVS) ultimately consumes the result.

Cleaning is an **ordered pipeline of single-purpose operations**, applied in a fixed
canonical order. Each operation inspects the string, makes at most one class of change,
and records an `HGVSFix` if it fired; operations are pure and order-dependent (for
example, leading junk and surrounding whitespace are stripped before structural
punctuation is examined, and casing is normalised before the gene/transcript
relationship is resolved). Callers may restrict cleaning to a subset of operations — an
allowlist for conservative use, or the full set minus a few — but selection only
*filters* the pipeline, never reorders it, so the canonical order is preserved
regardless of the subset chosen. The operations group into:

- **Stripping** — leading junk before the accession (e.g. a `GRCh38.p2` build tag or a
  stray `#`/`:`), all internal and surrounding whitespace and non-printable characters,
  wrapping quotes/backticks, trailing separators, and brackets only when genuinely
  unbalanced (balanced parentheses are preserved, since they are valid in
  uncertain-range and gene-symbol notation).
- **Structural punctuation** — collapsing a doubled colon, dot, or underscore
  (`NM_000059..4` → `NM_000059.4`); repairing a misplaced colon in the accession prefix;
  normalising a gene symbol wedged between extra colons or stray parentheses so that
  `NM_000059.4:(BRCA2):c.…` and `BRCA1(NM_000059.4)c.…` become the canonical
  `transcript(GENE):c.…`; collapsing a doubled kind token (`c.c.` → `c.`); and fixing a
  comma or colon used in place of the kind dot, or a period used in place of a
  substitution `>`.
- **Casing and prefixes** — uppercasing nucleotides in substitutions and del/ins/dup
  edits (`c.123delg` → `c.123delG`), lowercasing an uppercased mutation type while
  protecting gene symbols that contain those letters (so `NM_000059.4(INSR):c.…` is
  never corrupted to `insR`), restoring a missing `N` prefix or transcript underscore,
  and adding a missing `c.`/`g.` kind where the accession type makes it unambiguous.
- **Reconstruction and gene/transcript repair** — a lenient pattern that rebuilds the
  canonical `transcript(gene):kind.variant` shape from a mangled one (inserting a
  missing `:` or `.`, uppercasing the accession prefix and lowercasing the kind letter),
  and a final step that detects and repairs the common clinical mistake of swapping the
  gene symbol and transcript accession (`BRCA2(NM_000059.4):c.…` →
  `NM_000059.4(BRCA2):c.…`).

Each repair is reported as an `HGVSFix` carrying a severity (`WARNING` for an
unambiguous correction), a stable machine-readable code, a human-readable message, and
the before/after values — so a caller can surface, log, or audit exactly what was
changed rather than silently mutating user input. After the repair pipeline, a
validation pass flags problems that cleaning *cannot* fix as `ERROR`-level fixes (no
colon at all, a bare variant body with no reference sequence, or an insertion given an
integer length instead of the inserted sequence); these mark genuinely incomplete input
rather than formatting noise. Cleaning never raises by default — the returned string is
always the best attempt — but callers may opt into raising on the first error.

A separate, deliberately opt-in helper, `get_best_transcript_version()`, addresses
transcript-version drift: when the requested accession version is absent from a caller's
data, it returns the best available adjacent version (under a configurable up-then-down,
closest, or latest strategy) as a reported `HGVSFix`. It never substitutes a version
automatically — the caller decides whether to accept the change and rewrite the string —
keeping the user in control of HGVS compatibility. Examples in this section are
synthesised from the public ClinVar transcripts NM_000059.4 (BRCA2) and NM_001754.5
(RUNX1).

## Access and client libraries

**Local JSON**: `JSONDataProvider` loads a JSON.gz file into memory on initialisation
(typically {{ benchmark.grch38_load_time_s | dp(1) }} seconds for GRCh38 RefSeq),
building interval trees for region queries and dictionaries for transcript and gene
lookup. Transcript retrieval is O(1). Throughput of
{{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }}
transcripts/second is
{{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }}×
faster than UTA remote access (~{{ literature.uta_remote_tps }} transcript/second),
comparable to SeqRepo's {{ literature.seqrepo_speedup_fold }}× speedup over remote
sequence retrieval [@Hart2020].

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves the same JSON data
at cdotlib.org. `RESTDataProvider` fetches transcripts on demand — suitable for
occasional lookups without downloading the full file.

**biocommons/hgvs integration**: cdot implements
`biocommons.hgvs.dataproviders.interface.Interface`, providing `get_tx_info`,
`get_tx_exons`, `get_tx_seq`, `get_tx_mapping_options`, and related methods. This makes
cdot a drop-in replacement for UTA in any biocommons/hgvs pipeline:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
hdp = JSONDataProvider(["cdot.0.2.33.refseq.grch38.json.gz"])
```

Sequence data is supplied by SeqRepo or cdot's `FastaSeqFetcher` (local genome FASTA),
enabling fully offline operation without SeqRepo's installation overhead.

**PyHGVS integration**: `JSONPyHGVSTranscriptFactory` provides a transcript factory for
the Counsyl PyHGVS library, adding alignment gap support that PyHGVS lacks natively.
This fixes incorrect coordinate conversion for the fraction of RefSeq transcripts with
documented indels relative to the genome — the class of error responsible for the
{{ literature.brca2_accuracy_pct | dp(0) }}% BRCA2 accuracy reported by @Munz2015 in
tools without gap support.

## Canonical transcript selection

cdot's `CanonicalTranscriptSelector` maps gene symbols to MANE Select (or MANE Plus
Clinical) transcript accessions, enabling gene-name HGVS lookup — a common requirement
in clinical reporting pipelines that are given a gene name rather than a transcript
accession. This is to our knowledge the first HGVS data provider to expose programmatic
canonical transcript selection aligned to the MANE standard [@Morales2022; @Wright2023].

---

*Notes:*
- *Fill all bracketed numbers before submission (run compute_coverage.py,
  compute_benchmark.py, and the ClinVar benchmark). Coverage/ClinVar prose now lives in
  `results.md` R1 — Methods describes how, Results reports the numbers.*
- *String-cleaning subsection: examples synthesised from public NM_000059.4 (BRCA2) /
  NM_001754.5 (RUNX1) only — never corpus strings (CLAUDE.md / plan §3).*
- *The PyHGVS/gap support paragraph is important — don't cut it*
- *Canonical transcript section requires #36 to be implemented*
- *If over word count, the UTA CSV baseline source paragraph can be moved to
  supplementary*
