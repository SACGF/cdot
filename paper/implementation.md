# Implementation

*This is the main body section for an Application Note. Weave benchmarks in inline rather than a separate Results section.*
*Target: ~1,200 words across 4–5 subsections.*

---

## Data sources and generation

cdot processes transcript annotation from three sources, applied in priority order:

1. **RefSeq GFF3** — multiple historical NCBI annotation releases for GRCh37 and GRCh38 (see Supplementary Table S1). Ingesting multiple releases is essential for resolving historical HGVS strings: a clinical report from 2015 may reference an NM_ version that was retired in subsequent RefSeq releases but is still cited in patient records and ClinVar submissions.

2. **Ensembl GTF** — releases 75–81 (GRCh37) and 76–115 (GRCh38), plus T2T-CHM13v2.0 releases (Supplementary Table S2). Ensembl coverage adds [X] transcript accessions absent from RefSeq.

3. **UTA CSV** — the Universal Transcript Archive CSV dump provides alignments for accessions predating the oldest available GFF3/GTF files; overridden by official annotation where available.

A Snakemake pipeline parses each source using `GFF3Parser` or `GTFParser` (HTSeq-based), normalises contig names via bioutils, extracts CDS boundaries from start/stop codon features, and serialises to gzip-compressed JSON. Gene-symbol-to-HGNC-ID mappings come from GENCODE HGNC files.

## JSON format

Each transcript entry stores shared metadata (`gene_name`, `hgnc`, `biotype`, transcript accession and version) plus a `genome_builds` dict keyed by assembly name (e.g. `"GRCh38"`). Build-specific fields include `contig`, `strand`, `cds_start`, `cds_end`, and `exons` — a list of 6-tuples `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]`. The `gap` field stores the GFF3 gap string (e.g. `"M196 I1 M61"`) for transcripts with indels relative to the genome; the Python client converts this to HGVS CIGAR format (with I/D operators inverted per HGVS convention) at query time. Schema versioning (`schema_version` at root level) allows clients to reject incompatible files on load.

Canonical transcript tags — `mane_select`, `mane_plus_clinical`, `refseq_select`, and `ensembl_canonical` — are stored as build-specific fields where present. MANE Select covers >97% of protein-coding genes [Morales 2022] and is available for GRCh38; Ensembl canonical tags are available for GRCh37 and GRCh38.

## Access and client libraries

**Local JSON**: `JSONDataProvider` loads a JSON.gz file into memory on initialisation (typically [X] seconds for GRCh38 RefSeq), building interval trees for region queries and dictionaries for transcript and gene lookup. Transcript retrieval is O(1). Throughput of [500–1000] transcripts/second is [X]× faster than UTA remote access (~1 transcript/second), comparable to SeqRepo's 1300× speedup over remote sequence retrieval [Hart 2020].

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves the same JSON data at cdotlib.org. `RESTDataProvider` fetches transcripts on demand — suitable for occasional lookups without downloading the full file.

**biocommons/hgvs integration**: cdot implements `biocommons.hgvs.dataproviders.interface.Interface`, providing `get_tx_info`, `get_tx_exons`, `get_tx_seq`, `get_tx_mapping_options`, and related methods. This makes cdot a drop-in replacement for UTA in any biocommons/hgvs pipeline:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
hdp = JSONDataProvider(["cdot.0.2.33.refseq.grch38.json.gz"])
```

Sequence data is supplied by SeqRepo or cdot's `FastaSeqFetcher` (local genome FASTA), enabling fully offline operation without SeqRepo's installation overhead.

**PyHGVS integration**: `JSONPyHGVSTranscriptFactory` provides a transcript factory for the Counsyl PyHGVS library, adding alignment gap support that PyHGVS lacks natively. This fixes incorrect coordinate conversion for the fraction of RefSeq transcripts with documented indels relative to the genome — the class of error responsible for the 32% BRCA2 accuracy reported by Münz et al. [2015] in tools without gap support.

## Coverage and ClinVar benchmark

cdot covers [X] versioned transcript alignments across all builds and sources, compared to ~141,000 in UTA — a [9]×  increase (Figure 1). Ensembl transcripts ([X] unique accessions) are absent from UTA entirely. T2T-CHM13v2.0 contributes [Z] transcripts with no GRCh37/38 equivalent, in genomic regions previously inaccessible to sequencing.

To assess practical impact, [N] HGVS variant descriptions were retrieved from ClinVar [Landrum 2025] and resolved against cdot and UTA. cdot resolved [X]% versus [Y]% for UTA, with the improvement driven by Ensembl coverage and historical RefSeq versions absent from UTA's single snapshot (Figure 1, Supplementary Table S4).

## Canonical transcript selection

cdot's `CanonicalTranscriptSelector` maps gene symbols to MANE Select (or MANE Plus Clinical) transcript accessions, enabling gene-name HGVS lookup — a common requirement in clinical reporting pipelines that are given a gene name rather than a transcript accession. This is to our knowledge the first HGVS data provider to expose programmatic canonical transcript selection aligned to the MANE standard [Morales 2022; Wright 2023].

---

*Notes:*
- *Fill all bracketed numbers before submission*
- *The PyHGVS/gap support paragraph is important — don't cut it*
- *Canonical transcript section requires #36 to be implemented*
- *If over word count, the UTA CSV baseline source paragraph can be moved to supplementary*
