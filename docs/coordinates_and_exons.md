# Coordinates & exon alignments

This page explains how cdot stores exon coordinates and the alignment `gap` strings — the part of
the [JSON data format](json_data_format.md) that most often trips people up. See that page for the
full field-by-field reference; this one is the conceptual walk-through (issue
[#102](https://github.com/SACGF/cdot/issues/102)).

## The exon array

Each transcript has a `genome_builds` dict; each build has an `exons` array. To keep the files small,
every exon is a **positional array**, not an object:

```
[alt_start, alt_end, exon_id, cds_start, cds_end, gap]
```

| # | Field       | Coordinate system            | Meaning |
|---|-------------|------------------------------|---------|
| 0 | `alt_start` | genomic, **0-based**         | Start of the exon on the contig. |
| 1 | `alt_end`   | genomic, **0-based, exclusive** | End of the exon on the contig (half-open). |
| 2 | `exon_id`   | —                            | Exon ordinal in **stranded (transcript) order**, starting at 0. |
| 3 | `cds_start` | transcript/cDNA, **1-based** | Start of the exon in cDNA coordinates. |
| 4 | `cds_end`   | transcript/cDNA, **1-based, inclusive** | End of the exon in cDNA coordinates. |
| 5 | `gap`       | —                            | Alignment gap string (see below), or `null` when the exon aligns cleanly. |

> **Two coordinate systems in one array.** `alt_*` are 0-based half-open genomic coordinates (the
> contig is named by the build's `contig` field, e.g. `NC_000007.14`). `cds_*` are 1-based inclusive
> transcript (cDNA) coordinates. The names `cds_start`/`cds_end` are historical — they hold the
> exon's position along the **whole transcript**, not just the coding region.

The genomic span of an exon is `alt_end - alt_start`. The cDNA span is `cds_end - cds_start + 1`.
**When the exon aligns cleanly these are equal and `gap` is `null`.** When they differ, the `gap`
string explains the difference.

## Strand and exon order

`exon_id` is in transcript order: exon `0` is always the first exon transcribed.

- On the **`+`** strand, transcript order matches genomic order — `alt_start` increases with `exon_id`.
- On the **`-`** strand, transcript order is the reverse of genomic order — exon `0` has the
  **largest** genomic coordinates, and `cds_start`/`cds_end` increase as `alt_start` decreases.

`cds_start`/`cds_end` always increase with `exon_id` regardless of strand, because they are measured
along the transcript.

## The `gap` string

Most exons align base-for-base against the genome and store `gap = null`. But RefSeq/Ensembl
alignments sometimes contain small indels between the transcript sequence and the reference genome
(sequencing/assembly differences, edits, etc.). cdot records these the way the source GTF/GFF3 does —
a space-separated, CIGAR-like run-length string:

| Op    | Meaning |
|-------|---------|
| `M<n>`| `n` bases **align** (present in both transcript and genome). |
| `I<n>`| `n` bases are in the **transcript but not the genome** (insertion relative to the genome). |
| `D<n>`| `n` bases are in the **genome but not the transcript** (deletion relative to the genome). |

From this it follows that:

```
genomic span (alt_end - alt_start)   = sum of M + sum of D
cDNA span    (cds_end - cds_start+1) = sum of M + sum of I
```

These two identities are a good sanity check when reading the data.

### Worked example

```
[36552548, 36552986, 20, 2001, 2440, "M196 I1 M61 I1 M181"]
```

- Genomic span: `36552986 - 36552548 = 438`
- cDNA span: `2440 - 2001 + 1 = 440`
- gap `M196 I1 M61 I1 M181`: `M` total `196 + 61 + 181 = 438` (= genomic span ✓);
  `I` total `1 + 1 = 2`, so `M + I = 440` (= cDNA span ✓).

So this exon has **two single-base insertions in the transcript** relative to the genome. Reading the
gap left-to-right along the transcript:

```
transcript (cDNA, 1-based)   genome (contig, 0-based)
2001 ─ M196 ─ 2196           36552548 ─ 196 bp ─ 36552744   (196 bases align)
2197 ─ I1                     (no genomic base)              (1 base in transcript only)
2198 ─ M61  ─ 2258           36552744 ─ 61 bp  ─ 36552805   (61 bases align)
2259 ─ I1                     (no genomic base)              (1 base in transcript only)
2260 ─ M181 ─ 2440           36552805 ─ 181 bp ─ 36552986   (181 bases align)
```

A `D<n>` op would be the mirror image: the genomic cursor advances by `n` while the transcript cursor
stays put.

## How cdot uses the gap

When serving the biocommons HGVS interface, cdot converts each gap string into a biocommons CIGAR in
`get_tx_exons` (`AbstractJSONDataProvider._convert_gap_to_cigar`). Because biocommons describes the
alignment from the *genome → transcript* direction (the opposite of the GTF gap), the op meanings are
swapped on the way out: `M→=`, transcript-insertion `I→D`, genome-deletion `D→I`. You normally never
see this — it happens internally — but it explains why an `I` in the JSON shows up as a `D` in a
biocommons CIGAR.

## See also

- [JSON data format reference](json_data_format.md) — every field, auto-generated from the models.
- [Advanced usage](advanced_usage.md) — fixing HGVS input and bulk read-ahead retrieval.
