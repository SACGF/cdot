## [unreleased]

### Fixed

- #235 - Transcripts with unaligned bases between two exons (eg SHANK3 `NM_001372044.2`) failed with "Exons X and Y are not adjacent"

## [0.2.36] - 2026-11-04

### Fixed

- #230 - Positions on alignments whose first transcript bases don't align (eg `NM_000274.1`) resolved shifted by that many bases
- #230 - `FastaSeqFetcher` placed alignment gaps mirrored within minus-strand exons
- #232 - HGVS cleaning: a second variant after a protein annotation was dropped with it (`NM_001754.5:c.1A>G p.(Met1?) NM_000059.4:c.68del`)

## [0.2.35] - 2026-11-02

### Added

- #227 - HGVS cleaning: gene symbol before the accession with `:`/`.` (`BRCA2:NM_000059.4:c.68del`), `.`/`:` for the prefix underscore (`NM.000059.4`), typographic dashes, kind repeated in a range (`c.67+3_c.67+6del`), gene repeated after its parens, any stray leading/trailing punctuation or symbol (`•`, `†`, `.`, `:`...)

### Fixed

- #225 - HGVS cleaning regressions from 0.2.34: protein suffix after `,`/`;` left the separator, gene symbol doubled when the accession already had one, `NM)` typo, trailing `(`

## [0.2.34] - 2026-11-02

### Added

- #148 - HGVS cleaning: gene symbol before the accession becomes `ACCESSION(SYMBOL)` (`APC NM_001407446.1:c.-40G>A` → `NM_001407446.1(APC):c.-40G>A`)
- #148 - HGVS cleaning: drop a chromosome name before its own accession (`chr13 NC_000013.11:g.32316528del`)
- #145 - `models.Transcript.fake` for RefSeq mitochondrial fake transcripts (needs data schema >= 0.2.37)
- #166 - `clean_hgvs` warns (`HGVSCleanOpWarning`) when `ops` selects an op without one it depends on
- #171 - `get_transcript(tx_ac)` on data providers, the public way to get a transcript record
- #181 - Export `EnsemblTarkDataProvider` / `EnsemblTarkSeqFetcher` from `cdot.hgvs.dataproviders`, and `__all__` on public modules
- #181 - `FastaSeqFetcher` without pysam gives a helpful `ImportError` on construction
- #206 - Version fallback gives a versionless transcript a version (`NM_000059:c.36del` → `NM_000059.4:c.36del`) when all held versions agree at that position (`RESOLVED_VERSIONLESS_TRANSCRIPT`), otherwise refuses
- #215 - HGVS cleaning reports a kind that doesn't fit the accession (`NC_000013.11:c.68del`, `NM_000059.4:n.…`) as an ERROR, `KIND_REFERENCE_MISMATCH`. With a data provider, `fix_hgvs` also checks `c.`/`n.` against whether the transcript is coding (Ensembl, and gene resolution)
- #217 - Gene resolution handles `p.` (`BRCA2:p.Arg10Gly` → the MANE transcript's protein, `RESOLVED_GENE_TO_PROTEIN`) and `r.`. A gene symbol with `g.`/`m.`, or a chromosome name (`chr13`, `13`, `X`) with `c.`/`n.`/`r.`/`p.`, is a `KIND_REFERENCE_MISMATCH` ERROR

### Changed

- #170 - `get_tx_for_region` returns versioned accessions (matters for versionless-keyed providers like cdot_rest's)
- #170 - Unknown contig in `get_tx_exons` raises `HGVSDataNotAvailableError` not `ValueError`
- #171 - `prefer_consortium` renamed `annotation_consortium` (it was always a hard filter)
- #171 - `fallback_to_longest` compares transcript lengths itself rather than trusting provider order
- #171 - `fix_c_to_p` returns the input variant (not `None`) on error
- #171 - `fix_hgvs` warns `SKIPPED_GENE_RESOLUTION` for gene-only input with no `genome_build`
- #168 - `JSONDataProvider.get_tx_versions()` uses an index instead of scanning every transcript
- #179 - Python 3.10-3.14 classifiers
- #180 - Removed the no-op `cdot[hgvs]` install target (`hgvs` is a core dependency)

### Deprecated

- #171 - `prefer_consortium` (use `annotation_consortium`)
- #186 - `AbstractJSONDataProvider.is_version_substitution_safe()`, use `cdot.hgvs.version_safety.is_version_substitution_safe` (two transcript records) or `is_version_substitution_safe_by_accession`

### Fixed

- #221 - `FastaSeqFetcher` raises `GenomeMismatchError` for a slice containing unknown (`N`) bases, which broke normalization
- #150 - HGVS cleaning: missing colon no longer produces wrong but valid-looking HGVS; ambiguous cases are a `NO_COLON` ERROR
- #151 - HGVS cleaning: fix op ordering (quotes, separators, case, tabs)
- #152 - HGVS cleaning: handle all accession prefixes and kinds (lowercase `ng_`/`np_`, `m.`, `r.`, protein alleles)
- #165 - HGVS cleaning: del/dup count is checked against the range (`DEL_DUP_COUNT_MISMATCH` ERROR if it disagrees)
- #165 - HGVS cleaning: unbalanced bracket inside the allele is an ERROR (`UNBALANCED_BRACKETS`), not silently removed
- #165 - HGVS cleaning: drop a genome build name (`hg19`, `GRCh38` etc) before the reference
- #158 - Version safety check: malformed records give "cannot verify" not `TypeError`; non-coding transcripts are checked on exon structure
- #156 - `get_latest_combo_file_urls` finds the all-builds files (`cdot.data_release.ALL_BUILDS`)
- #154 - `RESTDataProvider`: only cache 404s as missing (other HTTP errors raise), use a session, 30s default timeout
- #155 - `RESTDataProvider` supports T2T-CHM13v2.0 (new `assemblies` argument)
- #169 - Close JSON file handles (`ResourceWarning`)
- #219 - Data providers match biocommons' UTA provider: `get_acs_for_protein_seq` returns the `MD5_` accession (fixes `IndexError` in `c_to_p` for a transcript with no protein ID), `get_tx_identity_info` raises `HGVSDataNotAvailableError` for an unknown transcript, region queries return `[]` for a non-splign method
- #209 - `fix_hgvs` no longer gives a `NO_TRANSCRIPT_FOR_GENE` ERROR for genomic or protein HGVS (`NC_000013.11:g.32316528del`)
- #211 - `get_tx_exons` raises `HGVSDataNotAvailableError` for an unknown transcript or a contig it isn't on, as biocommons' UTA provider does (was `None`, or exons labelled with the wrong contig)
- #213 - HGVS cleaning: `NC_000013.11(NM_000059.4):g.32316528del` keeps the genomic accession (`NC_000013.11:g.32316528del`, `DROPPED_TRANSCRIPT_FOR_GENOMIC_KIND`) instead of becoming `NM_000059.4:g.…`

## [0.2.33] - 2026-10-06

### Added

- #115 - `FastaSeqFetcher` warns when the transcript differs from the genome (`genome_mismatch` arg: warn/raise/off)
- #115 - `models.GenomeMismatch`

### Changed

- #113 - GRCh38 uses the `GRCh38.p14` assembly map, so patch scaffold (`NW_`) contigs work
- #113 - `FastaSeqFetcher` fetches `NW_` and `NT_` contigs

## [0.2.32] - 2026-09-27

### Added

- #131, #76 - `models`: `Transcript.translation`, `GenomeBuild.genome_mismatch` and `GenomeBuild.warnings`
- #131, #76 - `cdot.hgvs.fix_c_to_p` wraps biocommons `c_to_p()` using translation data (mito/selenocysteine codes, ribosomal slippage, translation exceptions), returning `(var_p, fixes)`

## [0.2.31] - 2026-08-19

### Added

- #112 - `resolve_missing_accession_prefix` restores a dropped RefSeq prefix (`000059.4:c.68del` → `NM_000059.4:c.68del`), run by `fix_hgvs`
- #112 - HGVS cleaning: semicolon/underscore for colon, space before kind, dangling version dot, more leading junk
- #129 - `_get_transcripts(tx_acs)` batch hook on data providers, so a DB/REST provider can override one method

### Changed

- #5 - Performance: LRU cache for `get_tx_exons()`

## [0.2.30] - 2026-06-25

### Fixed

- Packaging: wheel/sdist only contain the `cdot` package

## [0.2.29] - 2026-06-25

### Added

- #120 - Coordinate-safety check for the transcript version fallback (`USED_ADJACENT_VERSION_COORD_SAFE`), `cdot.hgvs.version_safety`
- #28 - `get_tx_versions(accession, genome_build=None)`
- #112 - HGVS cleaning: drop a genomic accession in the gene symbol slot (`NM_000059.4(NC_000013.11):c.68del`)
- #119 - HGVS cleaning works on all RefSeq accession types, not just `NM_`/`NC_`

### Changed

- #119 - Renamed `FIXED_DOUBLE_*` fix codes / `FIX_DOUBLE_*` ops to `*_MULTIPLE_*`
- #28, #120 - Version fallback refuses unsafe substitutions by default (`REFUSED_UNSAFE_VERSION`), see `on_unsafe_version`

### Fixed

- #119 - HGVS cleaning fixes (gene symbols that look like accessions, wrong fix codes)

## [0.2.28] - 2026-06-18

### Added

- #114 - Public ranking functions `rank_transcript_versions`, `rank_transcripts_for_gene`, `consortium_of`, `HGVSFix` message helpers, `_get_tags_by_tx_ac` batch hook
- #112 - HGVS cleaning: `:c:123`, parens around the accession, `BRCA1(NM_…)c.…`

### Changed

- #114 - `_get_transcript_tags` takes `tx_ac`
- #114 - Gene resolution and `get_tx_for_gene` return versioned accessions

### Fixed

- #112 - HGVS cleaning preserves LRG `t1`/`p1` suffixes
- `JSONDataProvider` loads 0.2.33+ data (`source` as list)
- #108, #109 - `RESTDataProvider.prefetch()` sends batches (and halves on 413/414)

## [0.2.27] - 2026-06-16

### Added

- #86 - Ensembl Tark Data Provider
- #27, #112 - HGVS cleaning: `cdot.hgvs.clean` (`clean_hgvs`), selectable ops
- #27, #36 - Gene HGVS: `resolve_gene_hgvs`, `fix_hgvs` (eg `BRCA2:c.36del` via MANE/canonical tags)
- #28 - Transcript version fallback (`resolve_transcript_version`, `fix_hgvs(version_fallback=...)`), off by default
- #36 - `get_tx_ac_tags_for_gene` returns transcripts ranked by tag
- #36 - `Ensembl_canonical` tag for GRCh37 Ensembl (data only)
- #27 - `RESTDataProvider.get_tx_ac_tags_for_gene`
- #37 - msgspec typed data models in `cdot.models`
- #108, #109 - `RESTDataProvider.prefetch()` / `prefetch_from_hgvs()` to warm the cache in one request

### Changed

- #27 - Gene resolution: match `MANE Select` tags with spaces, retry lowercase gene symbols uppercased
- #27 - `prefer_consortium` arg (RefSeq default) on gene resolution
- #62 - Dropped Python 3.9
- #111 - `cdot.models` accepts gene `biotype` as str or list
- #88 get_acs_for_protein_seq should return list not None
- #83 Ensembl files missing protein - breaking c_to_p (only affects data not client code)
- #17 RefSeq missing MT transcripts (only affects data not client code)
- #96 - All builds files now only contain 1 annotation consortium (RefSeq OR Ensembl not both)  
- #97 - Ensembl now has HGNC codes (used Gencode lookup) to match RefSeq
- #100 - Swapped cdot.cc domain references to cdotlib.org
- #106 - Collapse release nodes (only affects data not client code)
- Adjusted requirements - depend on hgvs core, bring in pysam for FASTA SeqFetcher
- Bug fixes + optimisations

### Fixed

- #27 - `clean_hgvs` no longer corrupts gene symbols like `INSR`, or strips balanced brackets
- #86 - `EnsemblTarkDataProvider`: coding transcripts without a UTR, `get_tx_exons` on an unaligned contig
- #36 - Longest transcript uses spliced length, not genomic span
- #36 - `resolve_gene_hgvs` accepts lowercase accessions
- #37 - `models.GenomeBuild` keeps `source`, `ccds`, `transcript_support_level`
- #39 - `ExonsFromGenomeFastaSeqFetcher(cache=False)` disables caching

## [0.2.26] 2024-08-15

Bumped version to 0.2.26 to catch up with data release. Only new client functionality is #81 'data_release' helper functions

All other changes in this release were for data (and contained in data_v0.2.26)

### Added

- #81 New 'data_release' code eg 'get_latest_combo_file_urls' that looks on GitHub to find latest data
- New GFFs: RefSeq RS_2023_10, Ensembl 111, 112
- #79 - RefSeq MT transcripts
- #66 - We now store 'Note' field (thanks holtgrewe for suggestion)
- client / JSON data schema version compatability check

### Changed

- #56 - Fix occasional UTA duplicated exons
- #57 - Correctly handle retrieving genomic position and dealing w/indels in GFF (thanks ltnetcase for reporting)
- #60 - Fix for missing protein IDs due to Genbank / GenBank (thanks holtgrewe)
- #64 - Split code/data versions. json.gz are now labelled according to data schema version (thanks holtgrewe)
- Renamed 'CHM13v2.0' to 'T2T-CHM13v2.0' so it could work with biocommons bioutils
- #72 - Correctly handle ncRNA_gene genes (thanks holtgrewe for reporting)
- #73 - HGNC ID was missing for some chrMT genes in Ensembl

## [0.2.21] - 2023-08-14

### Changed

- #45 - FastaSeqFetcher - fix alignment gaps properly
- #52 - Added transcripts from Ensembl 110 GRCh38 release
- #53 - UTA to cdot transcript start/end conversion issue  

## [0.2.20] - 2023-07-10

### Changed

- #50 - Biotype was missing in Ensembl transcripts 

## [0.2.19] - 2023-07-06

### Changed

- #49 - MT not converted to contigs correctly (GRCh37/Ensembl only) #49 
- Removed accidental logging

## [0.2.18] - 2023-07-05

### Added

- #44 - Support for mouse transcripts (Mus Musculus GRCm38 and GRCm39). Bash scripts only, never released as data, and removed in a later version (see #103)
- #47 - Implement HGVS DataProvider get_alignments_for_region

### Changed

- #45 - FastaSeqFetcher - handle deletions correctly (had swapped HGVS cigar projections around)
- #46 - HGVS DataProvider get_tx_info should properly handle alt_ac and alt_aln_method

## [0.2.17] - 2023-05-08

### Added

- #42 - Ensembl T2T CHM13v2.0

### Changed

- #43 - Contigs not converted to accession numbers properly (this was breaking local Biocommons HGVS conversion using 0.2.16 data)  

## [0.2.16] - 2023-04-12

### Added

- Added historical release 110 (2022-04-12) for T2T CHM13v2.0
- Added latest GRCh38.p14 release (2023-03-21)

## [0.2.15] - 2023-04-03

### Added

- Support for T2T CHM13v2.0

## [0.2.14] - 2023-03-21

### Added

- #39 - Fasta file SeqFetcher implementation
- Add Ensembl 109 GTF

### Changed

- #38 - Differing implementation of get_tx_for_region to hgvs one (reported by Manuel Holtgrewe) 
- #35 - Tags (ie MANE Select / RefSeq select etc) should be genome build specific 
- #34 - Stick to PyHGVS conventions, throw ValueError: transcript is required on missing transcript

## [0.2.13] - 2023-02-23

### Changed

- Fix for #25 - Pyhgvs data conversion - non-coding transcripts have bad cds start/end conversion
- Fix for #32 - Signature of get_pyhgvs_data consistent for all return statements

## [0.2.12] - 2022-12-08

### Added

- #30 - We now store "tag" attributes (eg "MANE Select", "RefSeq Select")
- Switch to using Ensembl GFF3 (so we can get tags out)
- Add Ensembl 108 GFF3

### Changed

- Fix for #25 -  GeneInfo currently fails for some records
- Fix for #27 -  Change URL for missing RefSeq GFFs

## [0.2.11] - 2022-09-27

### Added

- Now support all methods (get_gene_info, get_tx_for_gene, get_tx_for_region) for REST
- Add Ensembl 107 GTF

### Changed

- Ensembl gene info was missing "description"

## [0.2.10] - 2022-09-19

### Added

- [Implement get_gene_info](https://github.com/SACGF/cdot/issues/20) - For local JSON data only

### Changed

- Fixed issue [#23 UTA transcripts for PyHGVS](https://github.com/SACGF/cdot/issues/23)

## [0.2.9] - 2022-09-01

### Changed

- [BugFix for get_tx_for_region](https://github.com/SACGF/cdot/issues/22)


## [0.2.8] - 2022-08-29

### Added

- [Implemented get_pro_ac_for_tx_ac](https://github.com/SACGF/cdot/issues/14) (c_to_p can now generate p.HGVS)
- [Implemented get_tx_for_region](https://github.com/SACGF/cdot/issues/18) for local JSON data only

## [0.2.7] - 2022-05-19

### Added

- Add transcripts from latest RefSeq GRCh37 (105) and RefSeq GRCh38 (110)

### Changed

- Fixed default arguments bug where PyHGVS only worked on SACGF fork
- gtf_to_json now goes straight to cdot format (without intermediary PyReference format)
- UTA is not included in generation scripts by default, to enable, set environment variable UTA_TRANSCRIPTS=True
- Handle mismatches in UTA CIGAR alignments (convert to match (no indels) as GFF format has no support for mismatch)

## [0.2.6] - 2022-05-19

### Changed

- Fixed issue [Ensembl contigs g_to_c](https://github.com/SACGF/cdot/issues/9) - Ensembl JSON was using chrom names ie "17" instead of "NC_000017.11" for contig 

## [0.2.5] - 2022-04-14

### Changed

- PyHGVS conversion fix - non-coding cds_start/cds_end is set to start/end (not None)

## [0.2.4] - 2022-04-13

### Added

- Latest RefSeq (110) and Ensembl (106) transcripts

### Changed

- Fixed bug where all UTA transcripts were '-' strand
- Add "other_chroms" to combined historical file

## [0.2.3] - 2022-03-29

### Changed

- Fixed bug where HGNC not extracted properly from Ensembl GTFs
- Gene information is now included by default (only adds 5%)
- Clean artifacts from UTA data
- Support for [SACGF PyHGVS fork](https://github.com/SACGF/hgvs) (which adds alignment gap support)

## [0.2.2] - 2022-03-03

### Added

- Support for HTTPS (bought SSL certificate for REST server)

## [0.2.1] - 2022-03-03

### Added

- [Download/Convert UTA transcripts](https://github.com/SACGF/cdot/issues/1)
- [REST client](https://github.com/SACGF/cdot/issues/4) for [REST Service](https://github.com/SACGF/cdot_rest/)

### Changed

- [JSON format changed](https://github.com/SACGF/cdot/issues/2), separating common/build specific coordinates. This is so a transcript can contain data for multiple builds.
- [Use ijson to reduce RAM usage](https://github.com/SACGF/cdot/issues/7) - uses iterator vs loading all JSON into RAM

## [0.1.1] - 2022-01-19

### Added

- Initial commit

[unreleased]: https://github.com/SACGF/cdot/compare/v0.2.33...HEAD
[0.2.33]: https://github.com/SACGF/cdot/compare/v0.2.32...v0.2.33
[0.2.32]: https://github.com/SACGF/cdot/compare/v0.2.31...v0.2.32
[0.2.31]: https://github.com/SACGF/cdot/compare/v0.2.30...v0.2.31
[0.2.30]: https://github.com/SACGF/cdot/compare/v0.2.28...v0.2.30
[0.2.28]: https://github.com/SACGF/cdot/compare/v0.2.27...v0.2.28
[0.2.27]: https://github.com/SACGF/cdot/compare/v0.2.26...v0.2.27
[0.2.26]: https://github.com/SACGF/cdot/compare/v0.2.21...v0.2.26
[0.2.21]: https://github.com/SACGF/cdot/compare/v0.2.20...v0.2.21
[0.2.20]: https://github.com/SACGF/cdot/compare/v0.2.19...v0.2.20
[0.2.19]: https://github.com/SACGF/cdot/compare/v0.2.18...v0.2.19
[0.2.18]: https://github.com/SACGF/cdot/compare/v0.2.17...v0.2.18
[0.2.17]: https://github.com/SACGF/cdot/compare/v0.2.16...v0.2.17
[0.2.16]: https://github.com/SACGF/cdot/compare/v0.2.15...v0.2.16
[0.2.15]: https://github.com/SACGF/cdot/compare/v0.2.14...v0.2.15
[0.2.14]: https://github.com/SACGF/cdot/compare/v0.2.13...v0.2.14
[0.2.13]: https://github.com/SACGF/cdot/compare/v0.2.12...v0.2.13
[0.2.12]: https://github.com/SACGF/cdot/compare/v0.2.11...v0.2.12
[0.2.11]: https://github.com/SACGF/cdot/compare/v0.2.10...v0.2.11
[0.2.10]: https://github.com/SACGF/cdot/compare/v0.2.8...v0.2.10
[0.2.8]: https://github.com/SACGF/cdot/compare/v0.2.7...v0.2.8
[0.2.7]: https://github.com/SACGF/cdot/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/SACGF/cdot/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/SACGF/cdot/compare/v0.2.4...v0.2.5
[0.2.4]: https://github.com/SACGF/cdot/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/SACGF/cdot/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/SACGF/cdot/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/SACGF/cdot/compare/v0.1.1...v0.2.1
[0.1.1]: https://github.com/SACGF/cdot/releases/tag/v0.1.1
