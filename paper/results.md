# Results

*Target: ~1000–1200 words. Lead with the numbers. Every claim needs a table or figure.*

---

## 3.1 Transcript coverage

*[Table 1 — transcript counts by source and genome build]*

| Source | GRCh37 | GRCh38 | T2T-CHM13v2.0 | Total unique accessions |
|--------|--------|--------|----------------|------------------------|
| RefSeq (cdot) | [X] | [X] | [X] | [X] |
| Ensembl (cdot) | [X] | [X] | [X] | [X] |
| cdot combined | [X] | [X] | [X] | **[X]** |
| UTA (uta_20210129) | ~86,000 | ~55,000 | — | ~141,000 |

*[All figures to be filled from Snakemake summary stats — see #pipeline stats issue]*

cdot provides [X] versioned transcript alignments across all builds and sources, compared to ~141,000 in UTA uta_20210129, a [9]× increase. Ensembl transcripts ([X]) are absent from UTA entirely.

## 3.2 Speed comparison

*[Figure 1 — throughput comparison bar chart: cdot local vs cdot REST vs UTA remote vs UTA local]*

| Method | Throughput (transcripts/s) |
|--------|--------------------------|
| cdot local JSON | ~500–1000 |
| cdot REST (cdotlib.org) | ~[X] |
| UTA remote (public server) | ~1 |
| UTA local (PostgreSQL) | ~[X] |

*[Numbers to be confirmed with scripted benchmark — see pre-submission checklist]*

cdot local JSON access achieves approximately [500–1000]× higher throughput than UTA remote access. This speedup is comparable to the 1300× speedup reported for SeqRepo over remote sequence retrieval [Hart 2020].

## 3.3 ClinVar HGVS resolution benchmark

*[Table 2 — ClinVar benchmark: number of variants, resolution rate by tool]*

[N] HGVS variant descriptions were retrieved from ClinVar [Landrum 2025]. Each was attempted with: (i) cdot (GRCh38 RefSeq + Ensembl combined), (ii) UTA (uta_20210129), (iii) cdot (GRCh37 for older accessions).

| Tool | Variants resolved | Resolution rate |
|------|-------------------|----------------|
| cdot (GRCh38) | [X] | [X]% |
| cdot (GRCh37) | [X] | [X]% |
| cdot (combined, any build) | [X] | [X]% |
| UTA (uta_20210129) | [X] | [X]% |

*[Benchmark to be run — see issue #5 and pre-submission checklist]*

The improvement in resolution rate is driven primarily by: (a) Ensembl transcripts absent from UTA; (b) historical RefSeq versions covered by multiple GFF3 annotation releases but absent from UTA's single snapshot.

## 3.4 Alignment gap correctness

[Describe the BRCA2 / gapped alignment benchmark if run. Compare with/without gap support. Reference Münz 2015 (32% without gaps). Show cdot achieves [X]% on the same test set.]

*[To be filled — requires running the gap support benchmark described in pre-submission checklist]*

## 3.5 T2T-CHM13v2.0 support

*[Table 3 or supplementary — T2T unique transcripts]*

The T2T-CHM13v2.0 assembly [Nurk 2022] contains [X] annotated transcripts from RefSeq and [Y] from Ensembl. Of these, [Z] have no equivalent in GRCh37 or GRCh38, representing genes in regions previously inaccessible to sequencing. To our knowledge, cdot is the first HGVS coordinate conversion resource to support T2T-CHM13v2.0.

*[T2T unique count to be confirmed — see pre-submission checklist]*

## 3.6 Canonical transcript selection

[Once #36 is implemented: demonstrate gene-name HGVS lookup via MANE Select. Show coverage: MANE Select covers [X]% of protein-coding genes in cdot GRCh38. Show a worked example: BRCA2 → MANE Select transcript → HGVS string.]

*[Section to be written after #36 is implemented]*

---

*Notes:*
- *All bracketed numbers must be filled before submission*
- *Speed benchmark should be scripted and reproducible (commit the benchmark script)*
- *ClinVar benchmark (#5) is the single highest-impact result — prioritise this*
- *Consider whether gap correctness benchmark is feasible before submission deadline*
