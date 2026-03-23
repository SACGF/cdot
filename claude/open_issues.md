# cdot — Open GitHub Issues

Fetched 2026-03-23. All issues are open. Repo: https://github.com/SACGF/cdot/issues

---

## Infrastructure & Ops

### #106 — Release notes are too long
Release notes on GitHub take up too much space. Options: collapsible sections, or minimal notes that link to a README in source.

### #100 — Move domain off .cc
`cdot.cc` is blocked by SA Path (and presumably others). New domain `cdotlib.org` has been obtained. Plan: redirect old to new.

### #85 — Make a proper release
Depends on #82 (c_to_p) and #84 (FastaSeqFetcher improvements). Goal after those are done: get cdot officially recognised by biocommons HGVS — either listed as an alternative data provider in their docs, or at minimum get linked from [biocommons/hgvs#634](https://github.com/biocommons/hgvs/issues/634).

### #63 — Create VEP-compatible releases in CI
Build and release VEP-compatible JSON files from CI. Memory usage may be a concern; `orjson` suggested as a more memory-efficient alternative. Snakemake GitHub Actions could trigger the workflow automatically.

### #62 — Add tests to CI
No automated CI test suite yet. Snakemake also supports unit tests — could hook in with small fake data.

---

## Data Generation & Sources

### #104 — Bump minimum Python version in setup
`generate_transcript_data/` uses type annotations that require Python >3.8, but `setup.cfg` still declares `python_requires = >=3.8`. Also needs a separate optional dependency group (e.g. `cdot-build` extras) for packages only needed by data builders: `ijson`, `pyhgvs`, `HTSeq`, `biopython`. **Note:** HTSeq 2.1.2 has been released fixing the duplicated-tags issue.

### #103 — Snakemake support for Mouse
Old Bash scripts under `generate_transcript_data/Mus_musculus/` supported mouse. It's probably possible to pass a different YAML file to Snakemake to support mouse without major changes. If it works, add basic mouse versions to release data uploads. Known issues from a contributor:
- Scripts must be run from repo root, not from their directory.
- Mouse gene info file location: `https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/Mus_musculus.gene_info.gz`

### #101 — Refactor GTF/GFF code: split RefSeq vs Ensembl
Currently split on GTF vs GFF3 format, but the natural split is RefSeq vs Ensembl, as there is growing source-specific logic. Step 1: write comprehensive unit tests for GTF/GFF → JSON before refactoring. Also: Ensembl GFF3 handling is broken for protein versions (only GTF works); add explicit tests and a clear error message.

### #90 — Investigate GENCODE as a transcript source
Spun out of #89. GENCODE is Ensembl's curated gene annotation; Ensembl itself is supposed to use GENCODE transcripts. cdot already downloads GENCODE HGNC gz files (for name/HGNC mapping, added in #97). Question is whether GENCODE GTFs should also be used as a primary transcript source.

### #84 — FastaSeqFetcher improvements
FastaSeqFetcher is needed whenever transcripts aren't in seqrepo (which is common with cdot's broader coverage). Improvements needed:
- Add to docs / wiki.
- Detect when it's appropriate: Ensembl is fine; RefSeq GFFs sometimes flag mismatches in the `Note` field (e.g. `gap_count=1`). Should warn or refuse to map based on a configurable setting.
- Storing mismatches in the JSON may require a schema version bump, or a separate field alongside the existing CIGAR.

### #69 — Transcripts from NCBI/Ensembl API
NCBI GFF3 is available per-transcript via API (e.g. `https://eutils.ncbi.nlm.nih.gov/sviewer/viewer.cgi?report=gff3;id=NM_001001890.3`). Ensembl only exposes latest versions via REST. Idea: make a DataProvider that chains to these APIs as a last resort for transcripts not in the local JSON. (Ensembl Tark was spun off into #86 and is now implemented.)

### #54 — Diff coding coordinates from latest
When merging historical transcripts, some are discarded because their coding coordinates differ from the latest version. Needs investigation: is this a real coordinate change, a data quality issue, or a cdot bug?

### #51 — Historical GRCh38 RefSeq
NCBI released access to historical human transcript alignments ([blog post](https://ncbiinsights.ncbi.nlm.nih.gov/2023/06/29/access-to-historical-human-transcript-alignments/)). Could add ~40k new transcript versions. Some historical GFFs had exons being read twice; investigation found ~80% pass validation, ~80% of failures are `_validate_cdna_match()` errors (e.g. "cDNA match starts at 3 not 1").

### #16 — Add build patch version
Store genome build patch version (e.g. GRCh38.p14 vs GRCh38.p13). Could be extracted from the source URL for RefSeq; Ensembl builds would need a lookup. Useful for the paper as a known limitation section — patch versions can affect transcript alignments.

### #5 — Do some benchmarks
Benchmark plan: generate 100/1k/10k random HGVS from ClinVar, install SeqRepo locally, time JSON vs REST client. Existing informal results (500 random ClinVar HGVS, cdot.cc hosted in Australia):
- cdot REST: median 0.1 s/HGVS, 100% resolved
- UTA: median 1.84 s/HGVS, 17% resolved (missing data for the rest)

UTA missing ~80% of cdot-held transcripts skews the UTA timing — should restrict benchmark to transcripts present in both. A [truth set from genome medicine](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0396-7#Sec14) exists. Should make a proper script rather than ad-hoc runs.

---

## JSON Format & Schema

### #102 — Update docs: how exons work
Wiki page https://github.com/SACGF/cdot/wiki/Transcript-JSON-format needs expansion to document the exon tuple fields (chrom, start, end, exon ID, cDNA match gap string format `"M21 D1 M20"`, etc.).

### #78 — JSON breaking schema changes
Desired changes if/when a breaking schema change is ever made:
- **Exons:** change from fixed-length arrays to key/value dicts so new fields can be added without breaking length-dependent parsers.
- **Phasing:** add phasing array (see #76) — only needed if exons become dicts first.
- **Rename:** `cds_start`/`cds_end` to `tx_start`/`tx_end` was suggested by a contributor, but Dave notes these *are* coding coordinates; transcript extents are already recoverable as `exons[0][0]` and `exons[-1][1]`.

### #77 — Document JSON.gz file format
Wiki documentation needed for:
- Top-level JSON keys (`cdot_version`, `genome_builds`, `transcripts`, `genes`, `metadata`).
- Transcript format (shared fields vs build-specific fields).
- The gene index quirk: RefSeq uses numeric GeneID as the key; UTA-sourced transcripts use a fake symbol (prefixed `_`) because GeneID is not available from UTA.

### #76 — CDS phase / ribosomal frameshift support
Some transcripts use −1 ribosomal frameshifting (e.g. `NM_015068.3`, gene PEG10), where one base is "read twice" by overlapping CDS exons in the GFF3 (`join(480..1436,1436..2605)`). cdot currently doesn't model this — it just stores the exon boundaries without the frame-shift overlap. The GFF3 `phase` column (column 8) would be needed for correct handling. Dave's proposed approach: add an optional `phasing` array to the JSON (omit if all zeros, to avoid breaking existing consumers). Blocked until biocommons HGVS adds phase support, since there's no point storing it otherwise. GRCh38 GFF3 for this transcript curiously shows phase=0 even though the slippage is well-established.

### #37 — Provide python-attrs / dataclass typed classes
Request from external contributor to provide typed Python classes (via `attrs` or `dataclasses`) for the transcript data structures, which would serve as both documentation and enable typed import via `cattrs`. Dave is open to it as an addition alongside JSON, not a replacement.

---

## Features / API

### #93 — Move HGVSConverter abstraction from VariantGrid to cdot
VariantGrid has an interface that abstracts differences between PyHGVS and Biocommons HGVS, enabling gradual migration between them. Could be useful to others if extracted into cdot.

### #36 — Canonical transcripts
Relates to [biocommons/hgvs#747](https://github.com/biocommons/hgvs/issues/747). Goal: resolve `BRCA2:c.36del` by picking a canonical transcript automatically. The challenge is that canonical is stored differently per source/build:

| Annotation | Build | Tag |
|---|---|---|
| RefSeq | GRCh37 | RefSeq Select |
| RefSeq | GRCh38 | RefSeq Select, MANE |
| Ensembl | GRCh37 | n/a (basic etc., not definitive) — **now added via external file in v0.2.30** |
| Ensembl | GRCh38 | Ensembl_canonical, MANE |

Partial progress: GRCh37 Ensembl canonical transcripts now deployed in data release 0.2.30. Still TODO: cdot REST API should expose `canonical: true` in gene transcript listings; design a `CanonicalTranscriptSelector` class; think about Biocommons integration.

### #28 — Move transcript version up/down logic from VariantGrid to cdot
VariantGrid has logic to try adjacent transcript versions when the exact one isn't found. Could be moved to cdot as an opt-in feature, paired with warnings similar to the HGVS fixing framework in #27.

### #27 — cdot responsible for fixing bad HGVS / surfacing warnings
VariantGrid has substantial HGVS string cleaning code (handling spaces, missing colons, wrong case, unbalanced brackets, etc.) that should be moved to cdot for general use. Also needs a warnings/errors framework so callers can see what was auto-corrected. Known categories of bad input from Shariant:
- Incorrect case (e.g. `nm_002342.3:C.5094-11G>A`)
- Trailing quotes or leading tabs

---

## Data Quality

### #95 — Transcript version choice changed with UTA 20241220 (CDS off-by-1)
`NM_001017995.2` in cdot 0.2.28 has `cds_start`/`cds_end` shifted by 1 base vs cdot 0.2.24. Investigation confirmed this is a difference in UTA's alignment (not a cdot conversion bug) — actual HGVS c_to_g results still match VariantGrid. **Follow-up (from EmmaTudini):** document transcript source priority/order (currently: UTA first, then official GTFs override), and run a systematic benchmark of how often UTA vs RefSeq GTF alignments differ and why (genome patch versions? aligner parameters?). Interesting material for the cdot paper.

### #55 — FastaSeqFetcher: differences between RefSeq and FASTA-reconstructed transcripts
`NM_000399.3` (GRCh37/38) has a CIGAR gap (`M2190 D1 M282`). The RefSeq sequence has "T" inserted at position 790, but the FASTA-reconstructed sequence has "A" inserted at position 2697. Both have the same total length. The GFF flags `gap_count=1` and a `Note` about a non-frameshifting indel. Needs deeper investigation to determine whether cdot's reconstruction is wrong or RefSeq's stated sequence differs from the genomic FASTA.

### #89 — BAC gene transcripts returning `gene_name: None` (external user question)
A user processing gene-fusion data found that Ensembl transcripts for BAC-derived gene symbols (e.g. `RP5-899B16.3`) returned `gene_name: None`. Explanation: the GFF3 for that gene (`ENSG00000287820`) has no gene name in the `Name` attribute — the BAC symbol only appears in older GENCODE annotations, not in Ensembl's GFF3. The transcript exists in cdot but without a symbol. Workaround: iterate the JSON directly and build a `transcripts_per_gene` dict. Spun off #90 to investigate GENCODE as a source.
