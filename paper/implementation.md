# Implementation

*This is the main body section for an Application Note. Weave benchmarks in inline rather than a separate Results section.*
*Target: ~1,200 words across 4–5 subsections.*

---

## Data sources and generation

cdot processes transcript annotation from three sources, applied in priority order:

1. **RefSeq GFF3** — multiple historical NCBI annotation releases for GRCh37 and GRCh38 (see Supplementary Table S1). Ingesting multiple releases is essential for resolving historical HGVS strings: a clinical report from 2015 may reference an NM_ version that was retired in subsequent RefSeq releases but is still cited in patient records and ClinVar submissions.

2. **Ensembl GTF** — releases 75–81 (GRCh37) and 76–115 (GRCh38), plus T2T-CHM13v2.0 releases (Supplementary Table S2). Ensembl coverage adds {{ coverage.ensembl_unique_count | commas }} transcript accessions absent from RefSeq.

3. **UTA CSV** — the Universal Transcript Archive CSV dump provides alignments for accessions predating the oldest available GFF3/GTF files; overridden by official annotation where available.

A Snakemake pipeline parses each source using `GFF3Parser` or `GTFParser` (HTSeq-based), normalises contig names via bioutils, extracts CDS boundaries from start/stop codon features, and serialises to gzip-compressed JSON. Gene-symbol-to-HGNC-ID mappings come from GENCODE HGNC files.

## JSON format

Each transcript entry stores shared metadata (`gene_name`, `hgnc`, `biotype`, transcript accession and version) plus a `genome_builds` dict keyed by assembly name (e.g. `"GRCh38"`). Build-specific fields include `contig`, `strand`, `cds_start`, `cds_end`, and `exons` — a list of 6-tuples `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]`. The `gap` field stores the GFF3 gap string (e.g. `"M196 I1 M61"`) for transcripts with indels relative to the genome; the Python client converts this to HGVS CIGAR format (with I/D operators inverted per HGVS convention) at query time. Schema versioning (`schema_version` at root level) allows clients to reject incompatible files on load.

Canonical transcript tags — `mane_select`, `mane_plus_clinical`, `refseq_select`, and `ensembl_canonical` — are stored as build-specific fields where present. MANE Select covers >{{ literature.mane_coverage_pct | pct(0) }} of protein-coding genes [@Morales2022] and is available for GRCh38; Ensembl canonical tags are available for GRCh37 and GRCh38.

## Access and client libraries

**Local JSON**: `JSONDataProvider` loads a JSON.gz file into memory on initialisation (typically {{ benchmark.grch38_load_time_s | dp(1) }} seconds for GRCh38 RefSeq), building interval trees for region queries and dictionaries for transcript and gene lookup. Transcript retrieval is O(1). Throughput of {{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }} transcripts/second is {{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }}× faster than UTA remote access (~{{ literature.uta_remote_tps }} transcript/second), comparable to SeqRepo's {{ literature.seqrepo_speedup_fold }}× speedup over remote sequence retrieval [@Hart2020].

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves the same JSON data at cdotlib.org. `RESTDataProvider` fetches transcripts on demand — suitable for occasional lookups without downloading the full file.

**biocommons/hgvs integration**: cdot implements `biocommons.hgvs.dataproviders.interface.Interface`, providing `get_tx_info`, `get_tx_exons`, `get_tx_seq`, `get_tx_mapping_options`, and related methods. This makes cdot a drop-in replacement for UTA in any biocommons/hgvs pipeline:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
hdp = JSONDataProvider(["cdot.0.2.33.refseq.grch38.json.gz"])
```

Sequence data is supplied by SeqRepo or cdot's `FastaSeqFetcher` (local genome FASTA), enabling fully offline operation without SeqRepo's installation overhead.

**PyHGVS integration**: `JSONPyHGVSTranscriptFactory` provides a transcript factory for the Counsyl PyHGVS library, adding alignment gap support that PyHGVS lacks natively. This fixes incorrect coordinate conversion for the fraction of RefSeq transcripts with documented indels relative to the genome — the class of error responsible for the {{ literature.brca2_accuracy_pct | pct(0) }} BRCA2 accuracy reported by @Munz2015 in tools without gap support.

## Coverage and ClinVar benchmark

cdot covers {{ coverage.total_count | commas }} versioned transcript alignments across all builds and sources, compared to ~{{ literature.uta_count | commas }} in UTA — a {{ coverage.improvement_fold | fmt('.0f') }}× increase (Figure 1). Ensembl transcripts ({{ coverage.ensembl_unique_count | commas }} unique accessions) are absent from UTA entirely. T2T-CHM13v2.0 contributes {{ coverage.t2t_unique_count | commas }} transcripts with no GRCh37/38 equivalent, in genomic regions previously inaccessible to sequencing.

To assess practical impact, {{ clinvar.n_variants | commas }} HGVS variant descriptions were retrieved from ClinVar [@Landrum2025] and resolved against cdot and UTA. cdot resolved {{ clinvar.cdot_resolution_pct | pct(1) }} versus {{ clinvar.uta_resolution_pct | pct(1) }} for UTA, with the improvement driven by Ensembl coverage and historical RefSeq versions absent from UTA's single snapshot (Figure 1, Supplementary Table S4).

## Canonical transcript selection

cdot's `CanonicalTranscriptSelector` maps gene symbols to MANE Select (or MANE Plus Clinical) transcript accessions, enabling gene-name HGVS lookup — a common requirement in clinical reporting pipelines that are given a gene name rather than a transcript accession. This is to our knowledge the first HGVS data provider to expose programmatic canonical transcript selection aligned to the MANE standard [@Morales2022; @Wright2023].

---

*Notes:*
- *Fill all bracketed numbers before submission (run compute_coverage.py, compute_benchmark.py, and ClinVar benchmark)*
- *The PyHGVS/gap support paragraph is important — don't cut it*
- *Canonical transcript section requires #36 to be implemented*
- *If over word count, the UTA CSV baseline source paragraph can be moved to supplementary*
